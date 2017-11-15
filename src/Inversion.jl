__precompile__()

"""
This module has different inversion schemes.
* Alternate Optimization
* Multi-objective Optimization
"""

module Inversion
using Optim, LineSearches

"""
Parameters for alternating minimization of a single objective function, while updating different model parameters
"""
type ParamAM
	name::String			# prints name of the optimization
	max_roundtrips::Int64		# limit number of roundtrips
	max_reroundtrips::Int64		# there will be reroundtrips, if roundtrips don't converge 
	noptim::Int64 			# number of optimizations in every roundtrip
	optim_func::Vector{Function}	# optimization functions, call them to perform optimizations
	reinit_func::Function		# re-initialize function that is executed if roundtrips fail to converge
	fvec::Matrix{Float64}		# functionals
	fvec_init::Vector{Float64}	# 
	optim_tols::Vector{Float64}	# tolerance for each optimization
	roundtrip_tol::Float64		# tolerance for stopping roundtrips
	verbose::Bool
end


function ParamAM(optim_func;
	       name="",
	       noptim=length(optim_func),
	       optim_tols=[1e-3 for iop in 1:noptim],
	       roundtrip_tol=1e-3,
	       verbose=true,
	       reinit_func=x->randn(),
	       max_reroundtrips=1, re_init_flag=true, max_roundtrips=1)

	fvec=zeros(noptim, 2)
	fvec_init=zeros(noptim)


	pa=ParamAM(name, max_roundtrips, max_reroundtrips, noptim, optim_func, 
	  reinit_func,
	  fvec, fvec_init, optim_tols, roundtrip_tol, verbose)

	return pa



end

"""
Perform
alternating optimizations, updating different model parameters, computing a
same objective functional
"""
function go(pa::ParamAM)


	reroundtrip_converge=false
	itrr=0

	pa.verbose && println(pa.name, "\t alternate optimization")  

	while ((!reroundtrip_converge && itrr < pa.max_reroundtrips))
		itrr += 1
		pa.verbose && (itrr > 1) && println("failed to converge.. reintializing (",itrr,"/",pa.max_reroundtrips,")")
		pa.verbose && println("=================================================================================")  

	
		pa.fvec[:]=0.0
		pa.fvec_init[:]=0.0

		# execute re-initialization function
		pa.reinit_func(nothing)

		itr=0
		roundtrip_converge=false

		# print
		if(pa.verbose)
			@printf("trip\t|")
			for iop in 1:pa.noptim
				@printf("\top %d (%0.1e)\t|",iop, pa.optim_tols[iop])
			end
			@printf("\tvar(op) (%0.1e)\t",pa.roundtrip_tol)
			@printf("\n")
		end
		while !roundtrip_converge && itr < pa.max_roundtrips


			itr += 1

			# optimizations in each roundtrip
			for iop in 1:pa.noptim
				pa.fvec[iop,2]=pa.fvec[iop,1]
				pa.fvec[iop,1]=pa.optim_func[iop](nothing)
			end

			# store functionals at the first roundtrip
			if(iszero(pa.fvec_init))
				for iop in 1:pa.noptim
					pa.fvec_init[iop]=pa.fvec[iop,1]
				end
			end

			# normalize func for each optim
			for iop in 1:pa.noptim
				pa.fvec[iop,1] /= pa.fvec_init[iop]
			end

			if(pa.verbose)
				@printf("%d\t|",itr)
				for iop in 1:pa.noptim
					@printf("\t%0.6e\t|",pa.fvec[iop,1])
				end
			end
			# variance b/w objectives
			rf=vecnorm(pa.fvec[:,1])/vecnorm(fill(pa.fvec[1,1],pa.noptim))
			if(pa.verbose)
				@printf("\t%0.6e\t|",rf)
				@printf("\n")
			end
			(itr > 3) && (roundtrip_converge = (rf < pa.roundtrip_tol)) # do atleast 3 round trips before quitting
		end
		reroundtrip_converge = all(pa.fvec[:,1] .< pa.optim_tols[:])
		pa.verbose && println("=================================================================================")  
	end
end

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Multi-objective Inversion Framework

# Credits:
#	Pawan Bharadwaj
#	November 2017

"""
Parameters for multi-objective inversion
"""
type ParamMO
	noptim				# number of optimizations

	optim_func::Vector{Function}	# optimization functions
	optim_grad::Vector{Function} 	# gradient computing functions

	αvec::Vector{Float64}		# weights for each parameters 
	fvec_init::Vector{Float64}	# store initial functionals for scaling
	fvec::Vector{Float64}		# store functional in each call
	storage_temp::Vector{Float64}	# temp storage of gradient

	func::Function			# call this multi-objective functional
	grad!::Function			# call this multi-objective gradient 

end


"""
Scalarize a vector of multiple objective functions
* `x_init::Vector{Float64}` : model to compute initial functionals
* `optim_func::Vector{Function}` : to compute functionals
* `optim_grad::Vector{Function}` : to compute gradients

Only allocation is creating an additional vector for temporary gradient storage
"""
function ParamMO(;
		 noptim=1,
		 αvec=ones(noptim),
		 x_init=nothing,
		 optim_func=nothing,
		 optim_grad=nothing,
		 )
	(length(optim_func) != noptim) && error("length func")
	(length(optim_grad) != noptim) && error("length grad")
	(length(αvec) != noptim) && error("length αvec")

	fvec=zeros(noptim)
	fvec_init=zeros(noptim)

	# compute functional using x_init and store them 
	for iop in 1:noptim
		fvec_init[iop]=optim_func[iop](x_init)
	end
	any(iszero.(fvec_init)) && error("fvec_init cannot be zero")


	# func scalarizes fvec
	func=function func(x, pa)
		f=0.0
		for iop in 1:pa.noptim
			pa.fvec[iop]=pa.optim_func[iop](x)/pa.fvec_init[iop]
			f += (pa.fvec[iop]*pa.αvec[iop])
		end
		return f
	end

	# sum gradients with proper weights
	grad! =function grad!(storage, x, pa)
		f=0.0
		storage[:]=0.0
		for iop in 1:pa.noptim
			pa.fvec[iop]=pa.optim_grad[iop](pa.storage_temp, x)/pa.fvec_init[iop]
			f += (pa.fvec[iop]*pa.αvec[iop])
			for i in eachindex(storage)
				storage[i] += pa.storage_temp[i]*pa.αvec[iop]*inv(pa.fvec_init[iop])
			end
		end
		return f
	end

	# allocate storage_temp
	storage_temp=similar(x_init)
	pa=ParamMO(noptim,optim_func,optim_grad,αvec,fvec_init,fvec,storage_temp,func,grad!)
end


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# finite-differencing from gradient
macro forwardrule(x, e)
	x, e = esc(x), esc(e)
	quote
	$e = sqrt(eps(eltype($x))) * max(one(eltype($x)), abs($x))
	end
end

macro centralrule(x, e)
	x, e = esc(x), esc(e)
	quote
	$e = cbrt(eps(eltype($x))) * max(one(eltype($x)), abs($x))
	end
end

function finite_difference{T<:Number}(f::Function,x::T, dtype::Symbol = :central)
	if dtype == :forward
		@forwardrule x epsilon
		xplusdx = x + epsilon
		return (f(xplusdx) - f(x)) / epsilon
	elseif dtype == :central
		@centralrule x epsilon
		xplusdx, xminusdx = x + epsilon, x - epsilon
		return (f(xplusdx) - f(xminusdx)) / (epsilon + epsilon)
	else
		error("dtype must be :forward, :central")
	end
end



function finite_difference!{S <: Number, T <: Number}(f::Function,
	x::Vector{S},
	g::Vector{T},
	dtype::Symbol)
	# What is the dimension of x?
	n = length(x)

	gpar = SharedArray{eltype(g)}(size(g))
	xpar = SharedArray{eltype(x)}(size(x)); xpar = x;
	# Iterate over each dimension of the gradient separately.
	# Use xplusdx to store x + dx instead of creating a new vector on each pass.
	# Use xminusdx to store x - dx instead of creating a new vector on each pass.
	if dtype == :forward
		# Establish a baseline value of f(x).
		f_x = f(x)
		for i = 1:n
			@forwardrule x[i] epsilon
			oldx = x[i]
			x[i] = oldx + epsilon
			f_xplusdx = f(x)
			x[i] = oldx
			g[i] = (f_xplusdx - f_x) / epsilon
		end
	elseif dtype == :central
		@sync @parallel for i = 1:n
			@centralrule xpar[i] epsilon
			oldx = xpar[i]
			xpar[i] = oldx + epsilon
			f_xplusdx = f(xpar)
			xpar[i] = oldx - epsilon
			f_xminusdx = f(xpar)
			xpar[i] = oldx
			gpar[i] = (f_xplusdx - f_xminusdx) / (epsilon + epsilon)
		end
	else
		error("dtype must be :forward or :central")
	end

	g[:] = copy(gpar);
	return g
end



end # module
