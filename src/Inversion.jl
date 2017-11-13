__precompile__()

"""
This module has different inversion schemes.
* Alternate Optimization
* Multi-objective Optimization
"""

module Inversion

"""
Parameters for alternating minimization of a single objective function, while updating different model parameters
"""
type ParamAM
	name::String			# prints name of the optimization
	max_roundtrips::Int64		# limit number of roundtrips
	max_reroundtrips::Int64		# there will be reroundtrips, if roundtrips don't converge 
	noptim::Int64 			# number of optimizations in every roundtrip
	optim_func::Vector{Function}	# optimization functions
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
			rf=var(pa.fvec[:,1])
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


"""
Parameters for multi-objective inversion
"""
type ParamMO

	optim_func::Vector{Function}	# optimization functions

	αvec::Vector{Float64}		# 




end



function ParamMO(optim_func)




end

end # module
