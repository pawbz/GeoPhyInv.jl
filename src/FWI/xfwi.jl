

"""
Update inversion variable with the initial model.
At the same time, update the lower and upper bounds on the inversion variable.
"""
function initialize!(pa::Param)

	randn!(pa.mx.last_x)  # reset last_x

	# use mod_initial as starting model according to parameterization
	Models.Seismic_get!(pa.mx.x, pa.mod_initial, pa.parameterization)

	# bounds
	Seismic_xbound!(pa.mx.lower_x, pa.mx.upper_x, pa)

end

# finalize xfwi or wfwi on pa
function finalize!(pa::Param, xminimizer, attrib_mod)
	# update modm
	x_to_modm!(pa, xminimizer)

	# update calculated data using the minimizer
	F!(pa, nothing, attrib_mod)

	# put minimizer into modi
	Models.Seismic_reparameterize!(pa.modi, xminimizer, pa.parameterization) 

	# update initial model using minimizer, so that running xfwi! next time will start with updated model
	copyto!(pa.mod_initial, pa.modi)
end

"""
Full Waveform Inversion using `Optim` and `Ipopt` packages.
This method updates `pa.modm` and `pa.dcal`. More details about the
optional parameters can be found in the documentation of the `Optim` 
package.
The gradiet is computed using the adjoint state method.
`pa.modi` is used as initial model if non-zero.

# Arguments

* `pa::Param` : inversion object (updated inside method)
* `obj::Union{LS,LS_prior}` : which objective function?
  * `=LS()`
  * `=LS_prior([1.0, 0.5])`

# Optional Arguments

* `optim_options` : see Optim.jl package for optimization options
* `optim_scheme=LBFGS()` : see Optim.jl for other options
* `bounded_flag=true` : use box constraints, see Optim.jl

# Outputs

* updated `modm` in the input `Param`
* returns the result of optimization as an Optim.jl object
  * `=:migr_finite_difference` same as above but *not* using adjoint state method; time consuming; only for testing, TODO: implement autodiff here
"""
function fit!(pa::Param, obj::Union{LS,LS_prior}, attrib_mod; 
	   optim_scheme=LBFGS(),
	   optim_options=Optim.Options(show_trace=true,
	   store_trace=true, 
	   extended_trace=false, 
	   outer_iterations=2,
	   iterations=5,
	   f_tol=1e-5, 
	   g_tol=1e-8, 
	   x_tol=1e-5, ),
	   bounded_flag=false, solver=nothing, 
	   ipopt_options=[["max_iter", 5]])

	global to
	reset_timer!(to)

	println("updating modm and modi...")
	println("> xfwi: number of inversion variables:\t", xfwi_ninv(pa)) 

	initialize!(pa)

	f(x) = ζfunc(x, pa.mx.last_x,  pa, obj, attrib_mod)
	g!(storage, x) = ζgrad!(storage, x, pa.mx.last_x, pa, obj, attrib_mod)
	if(!bounded_flag)
		"""
		Unbounded LBFGS inversion, only for testing
		"""
		@timeit to "xfwi!" begin
			res = optimize(f, g!, pa.mx.x, optim_scheme, optim_options)
		end
	else
		"""
		Bounded LBFGS inversion
		"""
		if (solver == :ipopt)
			function eval_f(x)
			    return ζfunc(x, pa.mx.last_x,  pa, obj, attrib_mod)
			end


			function eval_grad_f(x, grad_f)
			    ζgrad!(grad_f, x, pa.mx.last_x, pa, obj, attrib_mod)
			end


			function void_g(x, g)
			    
			end

			function void_g_jac(x, mode, rows, cols, values)
			    
			end

			global prob = createProblem(size(pa.mx.x)[1], pa.mx.lower_x, pa.mx.upper_x, 0, 
				Array{Float64}(undef,0), Array{Float64}(undef,0), 0, 0,
					eval_f, void_g, eval_grad_f, void_g_jac, nothing)

			addOption(prob, "hessian_approximation", "limited-memory")

			if !(ipopt_options === nothing)
			    for i in 1:size(ipopt_options)[1]
				addOption(prob, ipopt_options[i][1], ipopt_options[i][2])
			    end
			end

			prob.x = pa.mx.x
			    
			@timeit to "xfwi!" begin
				res = solveProblem(prob)
			end
		else
			@timeit to "xfwi!" begin
				res = optimize(f, g!, pa.mx.lower_x, pa.mx.upper_x,pa.mx.x, Fminbox(optim_scheme), optim_options)
			end
                end
	end
	println(to)

	# print some stuff after optimization...
	if (solver == :ipopt)
		pa.verbose && println(ApplicationReturnStatus[res])
	else
		pa.verbose && println(res)
	end
	

	# finalize optimization, using the optimizer... 
	if solver == :ipopt
		finalize!(pa, prob.x, attrib_mod)
	else
		finalize!(pa, Optim.minimizer(res), attrib_mod)
	end


	return res
end


"""
Return gradient at the first iteration, i.e., a migration image
"""
function migr!(pa::Param, attrib_mod)

	global to
	reset_timer!(to)

	initialize!(pa)

	g1=pa.mx.gm[1]
	@timeit to "xfwi!" begin
		grad!(g1, pa.mx.x, pa.mx.last_x,  pa, attrib_mod)
	end
	println(to)

	# convert gradient vector to model
	visualize_gx!(pa.gmodm, pa.modm, pa.gmodi, pa.modi, g1, pa)

	return pa.gmodi, g1
end


"""
Return gradient at the first iteration, i.e., a migration image, without using
adjoint state method.
"""
function migr_fd!(pa::Param, attrib_mod)
	global to
	reset_timer!(to)

	initialize!(pa)
	gx=pa.mx.gm[1]

	println("> number of functional evaluations:\t", xfwi_ninv(pa)) 

	f(x) = ζfunc(x, pa.mx.last_x,  pa, LS(), attrib_mod)

	gx = Calculus.gradient(f,pa.mx.x) # allocating new gradient vector here
	
	println(to)

	# convert gradient vector to model
	visualize_gx!(pa.gmodm, pa.modm, pa.gmodi, pa.modi, gx, pa)

	return pa.gmodi, gx
end



		#=
		Harvest the Optim result to plot functional and gradient --- to be updated later
		if(extended_trace)
			# convert gradient vector to model
			gmodi = [Models.Seismic_zeros(pa.modi.mgrid) for itr=1:Optim.iterations(res)]
			gmodm = [Models.Seismic_zeros(pa.modm.mgrid) for itr=1:Optim.iterations(res)] 
			modi = [Models.Seismic_zeros(pa.modi.mgrid) for itr=1:Optim.iterations(res)]
			modm = [Models.Seismic_zeros(pa.modm.mgrid) for itr=1:Optim.iterations(res)]
			for itr=1:Optim.iterations(res)
				# update modm and modi
				Seismic_x!(modm[itr], modi[itr], Optim.x_trace(res)[itr], pa, -1)

				# update gmodm and gmodi
				Seismic_gx!(gmodm[itr],modm[itr],gmodi[itr],modi[itr],Optim.trace(res)[itr].metadata["g(x)"],pa,-1)
			end
			
		else
			f = Optim.minimum(res)
		end
		=#

