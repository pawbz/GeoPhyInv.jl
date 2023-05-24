

"""
Update inversion variable with the initial model.
At the same time, update the lower and upper bounds on the inversion variable.
"""
function initialize!(pa::PFWI)

	randn!(pa.mx.last_x)  # reset last_x

	# use mod_initial as starting model according to parameterization
	copyto!(pa.mx.x, pa.mod_initial, pa.parameterization)

	# bounds
	Seismic_xbound!(pa.mx.lower_x, pa.mx.upper_x, pa)

end

# finalize xfwi or wfwi on pa
function finalize!(pa::PFWI, xminimizer)
	# update mediumm
	x_to_mediumm!(pa, xminimizer)

	# update calculated data using the minimizer
	F!(pa, nothing)

	# put minimizer into mediumi
	copyto!(pa.mediumi, xminimizer, pa.parameterization) 

	# update initial model using minimizer, so that running xfwi! next time will start with updated model
	copyto!(pa.mod_initial, pa.mediumi)
end

"""
```julia
update!(pa)
```
Full Waveform Inversion using `Optim` and `Ipopt` packages.
This method updates `pa.mediumm` and `pa.dcal`. More details about the
optional parameters can be found in the documentation of the `Optim` 
package.
The gradient is computed using the adjoint state method.
`pa.mediumi` is used as initial model if non-zero.

# Arguments

* `pa::PFWI` : inversion object (updated inside method)

# Keywords Arguments
* `solver` Choose `Symbol`
  * `=:ipopt` for `Ipopt.jl` 
  * `=:optim` for `Optim.jl`

# Optional Arguments

* `optim_options` : see `Optim.jl` package for optimization options
* `optim_scheme=LBFGS()` : see `Optim.jl` for other options
* `bounded_flag=true` : constrained inversion using medium parameter bounds (`Ipopt.jl` only has bounded inversion)
* `ipopt_options` : see `Ipopt.jl` 
"""
function update!(pa::PFWI{T1,T2,T3};
	   optim_scheme=LBFGS(),
	   optim_options=Optim.Options(show_trace=true,
	   store_trace=true, 
	   extended_trace=false, 
	   outer_iterations=2,
	   iterations=5,
	   f_tol=1e-5, 
	   g_tol=1e-8, 
	   x_tol=1e-5, ),
	   bounded_flag=true, solver=nothing, 
	   ipopt_options=[["max_iter", 5]]) where {T1,T2,T3<:Union{LS,LS_prior}}

	global fwi_to
	reset_timer!(fwi_to)

	println("updating mediumm and mediumi...")
	println("> xfwi: number of inversion variables:\t", xfwi_ninv(pa)) 

	initialize!(pa)

	f(x) = ζfunc(x, pa.mx.last_x,  pa)
	g!(storage, x) = ζgrad!(storage, x, pa.mx.last_x, pa)
	if(!bounded_flag)
		(solver==:ipopt) && @info("bounded_flag not functional for Ipopt, switching to Optim")
		"""
		Unbounded LBFGS inversion, only for testing
		"""
		@timeit_debug fwi_to "xfwi!" begin
			res = optimize(f, g!, pa.mx.x, optim_scheme, optim_options)
		end
	else
		"""
		Bounded LBFGS inversion
		"""
		if (solver == :ipopt)
			function eval_f(x)
				return ζfunc(x, pa.mx.last_x,  pa)
			end


			function eval_grad_f(x, grad_f)
				ζgrad!(grad_f, x, pa.mx.last_x, pa)
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
			    
			@timeit_debug fwi_to "xfwi!" begin
				res = solveProblem(prob)
			end
		else
			@timeit_debug fwi_to "xfwi!" begin
				res = optimize(f, g!, pa.mx.lower_x, pa.mx.upper_x,pa.mx.x, Fminbox(optim_scheme), optim_options)
			end
                end
	end
	println(fwi_to)

	# print some stuff after optimization...
	if (solver == :ipopt)
		pa.verbose && println(ApplicationReturnStatus[res])
	else
		pa.verbose && println(res)
	end
	

	# finalize optimization, using the optimizer... 
	if solver == :ipopt
		finalize!(pa, prob.x)
	else
		finalize!(pa, Optim.minimizer(res))
	end


	return res
end


"""
```julia
update!(pa)
```
Return gradient at the first iteration, i.e., a migration image, when `attrib_inv==Migr()`.
"""
function update!(pa::PFWI{T1,T2,Migr}) where {T1,T2}

	global fwi_to
	reset_timer!(fwi_to)

	initialize!(pa)

	g1=pa.mx.gm[1]
	@timeit_debug fwi_to "xfwi!" begin
		grad!(g1, pa.mx.x, pa.mx.last_x,  pa)
	end
	println(fwi_to)

	# convert gradient vector to model
	visualize_gx!(pa.gmediumm, pa.mediumm, pa.gmediumi, pa.mediumi, g1, pa)

	return pa.gmediumi, g1
end


"""
```julia
update!(pa)
```
Return gradient at the first iteration, i.e., a migration image, without using
adjoint state method, when `attrib_inv==Migr_FD()`.
"""
function update!(pa::PFWI{T1,T2,Migr_FD}) where {T1,T2}   
	global fwi_to
	reset_timer!(fwi_to)

	initialize!(pa)
	gx=pa.mx.gm[1]

	println("> number of functional evaluations:\t", xfwi_ninv(pa)) 

	f(x) = ζfunc(x, pa.mx.last_x,  pa)

	gx = Calculus.gradient(f,pa.mx.x) # allocating new gradient vector here
	
	println(fwi_to)

	# convert gradient vector to model
	visualize_gx!(pa.gmediumm, pa.mediumm, pa.gmediumi, pa.mediumi, gx, pa)

	return pa.gmediumi, gx
end
                

"""
Return posterior covariance matrix and estimate as given by Best Linear Unbiased Analysis.
    # Arguments
    - H_obs : observation operator
    - C_obs : observation error covariance matrix
    - C_prior : prior covariance matrix
    - est_0 : prior estimate
    - y_obs : observation vector
    - compute_update : if true returns posterior mean and covariance, if false only covariance
    - diag : if true, assumes that observation covariance matrix is diag and uses it to siply inverse calculation
                
    # Output
    Posterior covariance and mean (if compute_update = true)
"""
function blue(H_obs, C_obs, C_prior, est_0, y_obs; compute_update = false, diag = false)
    C_prior_inv = inv(C_prior)
    C_prior_inv = 0.5 * (C_prior_inv + C_prior_inv')
    
    if diag
        C_obs_inv = 1 ./ C_obs
    else
        C_obs_inv = inv(C_obs)
        C_obs_inv = 0.5 * (C_obs_inv + C_obs_inv')
    end
    
    C_post = inv((H_obs' * C_obs_inv * H_obs + C_prior_inv))
    
    if compute_update
        est = C_post * (C_prior_inv * est_0 + H_obs' * C_obs_inv * y_obs)
        return C_post, est;
    else
        return C_post;    
    end
end



		#=
		Harvest the Optim result to plot functional and gradient --- to be updated later
		if(extended_trace)
			# convert gradient vector to model
			gmediumi = [Medium_zeros(pa.mediumi.grid) for itr=1:Optim.iterations(res)]
			gmediumm = [Medium_zeros(pa.mediumm.grid) for itr=1:Optim.iterations(res)] 
			mediumi = [Medium_zeros(pa.mediumi.grid) for itr=1:Optim.iterations(res)]
			mediumm = [Medium_zeros(pa.mediumm.grid) for itr=1:Optim.iterations(res)]
			for itr=1:Optim.iterations(res)
				# update mediumm and mediumi
				Seismic_x!(mediumm[itr], mediumi[itr], Optim.x_trace(res)[itr], pa, -1)

				# update gmediumm and gmediumi
				Seismic_gx!(gmediumm[itr],mediumm[itr],gmediumi[itr],mediumi[itr],Optim.trace(res)[itr].metadata["g(x)"],pa,-1)
			end
			
		else
			f = Optim.minimum(res)
		end
		=#

