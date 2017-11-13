# blind deconvolution
__precompile__()

module Decon
import JuMIT.DSP
import JuMIT.Misfits
using Optim, LineSearches

type Param
	gfobs::Array{Float64,2}
	gf::Array{Float64,2}
	gfprecon::Array{Float64,2}
	dgf::Array{Float64,2}
	ntgf::Int64
	dobs::Array{Float64,2}
	dcal::Array{Float64,2}
	ddcal::Array{Float64,2}
	wavobs::Vector{Float64}
	wav::Vector{Float64}
	wavprecon::Vector{Float64}
	wavmat::Array{Float64,2}
	dwavmat::Array{Float64,2}
	spow2::Array{Complex{Float64},2}
	rpow2::Array{Complex{Float64},2}
	wpow2::Array{Complex{Float64},2}
	nt::Int64
	nr::Int64
	attrib_inv::Symbol
	verbose::Bool
	fftplan::Base.DFT.FFTW.cFFTWPlan
	ifftplan::Base.DFT.ScaledPlan
	xgf::Vector{Float64}
	last_xgf::Vector{Float64}
	dfgf::Optim.UninitializedOnceDifferentiable{Void}
	xwav::Vector{Float64}
	last_xwav::Vector{Float64}
	dfwav::Optim.UninitializedOnceDifferentiable{Void}
	func_max_save::Vector{Float64}
	α::Vector{Float64}
end

"""
`gfprecon` : a preconditioner applied to each Greens functions [ntgf]
"""
function Param(ntgf, nt, nr; 
	       gfprecon=nothing,
	       wavprecon=nothing,
	       dobs=nothing, gfobs=nothing, wavobs=nothing, verbose=false, attrib_inv=:gf,) 
	(dobs===nothing) && (dobs=zeros(nt, nr))
	dcal=zeros(nt, nr)
	ddcal=zeros(nt, nr)
	
	# initial values are random
	wav=randn(nt)
	wavmat=repeat(wav,outer=(1,nr))
	dwavmat=similar(wavmat)
	gf=randn(ntgf,nr)
	dgf=similar(gf)

	# use maximum threads for fft
	# FFTW.set_num_threads(Sys.CPU_CORES)
	# fft dimension for plan
	np2=nextpow2(maximum([2*nt, 2*ntgf]))
	# create plan
	fftplan = plan_fft!(complex.(zeros(np2,nr)),[1],flags=FFTW.PATIENT)
	ifftplan = plan_ifft!(complex.(zeros(np2,nr)),[1],flags=FFTW.PATIENT)

	# preallocate pow2 vectors
	spow2=complex.(zeros(np2,nr))
	rpow2=complex.(zeros(np2,nr))
	wpow2=complex.(zeros(np2,nr))

	
	# convert initial model to the inversion variable
	xgf = zeros(ntgf*nr);
	last_xgf = randn(size(xgf)) # reset last_x

	# create df for optimization
	dfgf = OnceDifferentiable(x -> func_grad!(nothing, x, last_xgf, pa),
	  (storage, x) -> func_grad!(storage, x, last_xgf, pa))

	xwav = zeros(nt);
	last_xwav = randn(size(xwav)) # reset last_x

	# create df for optimization
	dfwav = OnceDifferentiable(x -> func_grad!(nothing, x, last_xwav, pa),
	  (storage, x) -> func_grad!(storage, x, last_xwav, pa))

	# create dummy gfobs if necessary
	(gfobs===nothing) && (gfobs=zeros(gf))
	# create dummy wavobs if necessary
	(wavobs===nothing) && (wavobs=zeros(wav))

	# create gf precon
	(gfprecon===nothing) && (gfprecon=ones(gf))

	# create gf precon
	(wavprecon===nothing) && (wavprecon=ones(wav))

	pa = Param(gfobs, gf, gfprecon, dgf, ntgf, dobs, dcal, ddcal, wavobs, wav, wavprecon, wavmat, dwavmat,
	     spow2, rpow2, wpow2,
	     nt, nr, attrib_inv, verbose, fftplan, ifftplan,
	     xgf,last_xgf,dfgf,
	     xwav,last_xwav,dfwav, zeros(3), [1., 0.0, 0.0])

	if iszero(pa.dobs) 
		((gfobs===nothing) | (wavobs===nothing)) && error("need gfobs and wavobs")
		F!(pa, 
     		pa.xwav, pa.last_xwav, # dummy
     		mattrib=:gfwav, dattrib=:dobs, gf=gfobs, wav=wavobs)
	end

	return pa
end

function ninv(pa)
	if(pa.attrib_inv == :wav)
		return pa.nt
	else(pa.attrib_inv == :gf)
		return pa.ntgf*pa.nr
	end
end


function error(pa) 
	fwav, α = Misfits.error_after_autocorr_scaling(pa.wav, pa.wavobs)
	fgf, α = Misfits.error_after_autocorr_scaling(pa.gf, pa.gfobs)
	f = Misfits.fg_cls!(nothing, pa.dcal, pa.dobs)

	return fwav, fgf, f
end 


function model_to_x!(x, pa)
	if(pa.attrib_inv == :wav)
		for i in eachindex(x)
			x[i]=pa.wav[i]*pa.wavprecon[i]
		end
	else(pa.attrib_inv == :gf)
		for i in eachindex(x)
			x[i]=pa.gf[i]*pa.gfprecon[i] 		# multiply by gfprecon
		end
	end
	return x
end


function x_to_model!(x, pa)
	if(pa.attrib_inv == :wav)
		xn=vecnorm(x)
		for i in eachindex(pa.wav)
			if(iszero(pa.wavprecon[i]))
				pa.wav[i]=0.0
			else
				pa.wav[i]=x[i]/pa.wavprecon[i]
			end
		end
		scale!(pa.wav, inv(xn))
	else(pa.attrib_inv == :gf)
		for i in eachindex(pa.gf)
			if(iszero(pa.gfprecon[i]))
				pa.gf[i]=0.0
			else
				pa.gf[i]=x[i]/pa.gfprecon[i]
			end
		end
	end
	return pa
end

function F!(
			pa::Param,
			x::Vector{Float64},
			last_x::Vector{Float64};
			gf=nothing,
			wav=nothing,
			mattrib::Symbol=:x,
			dattrib::Symbol=:dcal,
			   )

	if(x!=last_x)

		if(mattrib == :x)
			x_to_model!(x, pa)
			gf=pa.gf
			wav=pa.wav
		elseif(mattrib == :gf)
			wav=pa.wav
			(gf===nothing) && error("need gf")
		elseif(mattrib == :wav)
			gf=pa.gf
			(wav===nothing) && error("need wav")
		elseif(mattrib == :gfwav)
			(gf===nothing) && error("need gf")
			(wav===nothing) && error("need wav")
		else
			error("invalid mattrib")
		end

		#pa.verbose && println("updating buffer")
		copy!(last_x, x)
		for ir in 1:pa.nr
			pa.wavmat[:,ir]=wav
		end
		DSP.fast_filt!(pa.ddcal, gf, pa.wavmat, :s, nwplags=pa.nt-1,
		 nwnlags=0,nsplags=pa.nt-1,nsnlags=0,nrplags=pa.ntgf-1,nrnlags=0,
		  np2=size(pa.fftplan,1), fftplan=pa.fftplan, ifftplan=pa.ifftplan,
		  spow2=pa.spow2, rpow2=pa.rpow2, wpow2=pa.wpow2,
		  			)
		copy!(getfield(pa, dattrib), pa.ddcal)
		return pa
	end
end

function func_grad!(storage, x::Vector{Float64},
			last_x::Vector{Float64}, 
			pa)

	α=1e-3

	#pa.verbose && println("computing gradient...")

	# x to w 
	x_to_model!(x, pa)

	F!(pa, x, last_x)

	if(storage === nothing)
		# compute misfit and δdcal
		f1 = Misfits.error_squared_euclidean!(nothing, pa.dcal, pa.dobs, nothing)
	else
		f1 = Misfits.error_squared_euclidean!(pa.ddcal, pa.dcal, pa.dobs, nothing)
		Fadj!(pa, x, storage, pa.ddcal)
	end

	if(pa.func_max_save[1] == 0.0) 
		pa.func_max_save[1]=f1 # store first value of f
	end
#
	fac1=pa.α[1]/pa.func_max_save[1]
#
#
	if((pa.attrib_inv == :gf) && !(storage === nothing))
		f2 = Misfits.error_weighted_norm!(pa.dgf,pa.gf, pa.gfprecon)
	else
		f2 = Misfits.error_weighted_norm!(nothing,pa.gf, pa.gfprecon)
	end
	if(pa.func_max_save[2] == 0.0) 
		pa.func_max_save[2]=f2 # store first value of f
	end

	fac2=pa.α[2]/pa.func_max_save[2]

	if(!(storage === nothing) && (pa.attrib_inv == :gf))
		for i in eachindex(storage)
			storage[i] = fac1*storage[i]+fac2*pa.dgf[i]
		end
	end
#
	J = f1*fac1 + f2*fac2


	return J

end


# add model based constraints here

# all the greens' functions have to be correlated

# exponential-weighted norm for the green functions

#  



"""
Apply Fadj to 
"""
function Fadj!(pa, x, storage, dcal)
	if(pa.attrib_inv == :wav)
		DSP.fast_filt!(dcal, pa.gf, pa.wavmat, :w, nwplags=pa.nt-1, 
		 nwnlags=0,nsplags=pa.nt-1,nsnlags=0,nrplags=pa.ntgf-1,nrnlags=0,
		  np2=size(pa.fftplan,1), fftplan=pa.fftplan, ifftplan=pa.ifftplan,
		  spow2=pa.spow2, rpow2=pa.rpow2, wpow2=pa.wpow2,
			)
		storage[:] = 0.
		for i in 1:size(pa.wavmat,2)
			for j in 1:size(pa.wavmat,1)
				storage[j] += pa.wavmat[j,i]
			end
		end
		# factor, because wav was divided by norm of x
		xn=vecnorm(x)
		for i in eachindex(storage)
			if(iszero(pa.wavprecon[i]))
				storage[i]=0.0
			else
				storage[i] = storage[i]*(xn-x[i]*x[i]/xn)/xn/xn/pa.wavprecon[i]
			end
		end

	else(pa.attrib_inv == :gf)
		for ir in 1:pa.nr
			for it in 1:size(pa.wavmat,1)
				pa.wavmat[it,ir]=pa.wav[it]
			end
		end
		DSP.fast_filt!(dcal, pa.dgf, pa.wavmat, :r, nwplags=pa.nt-1,
		 nwnlags=0,nsplags=pa.nt-1,nsnlags=0,nrplags=pa.ntgf-1,nrnlags=0,
		  np2=size(pa.fftplan,1), fftplan=pa.fftplan, ifftplan=pa.ifftplan,
		  spow2=pa.spow2, rpow2=pa.rpow2, wpow2=pa.wpow2,
			)
		copy!(storage, pa.dgf)
		for i in eachindex(storage)
			if(iszero(pa.gfprecon[i]))
				storage[i]=0.0
			else
				storage[i]=pa.dgf[i]/pa.gfprecon[i]
			end
		end

	end
	return storage
end

# core algorithm
function update!(pa::Param, x, last_x, df; store_trace::Bool=false, 
		 extended_trace::Bool=false, 
	     time_limit=Float64=2.0*60., 
	     f_tol::Float64=1e-8, g_tol::Float64=1e-8, x_tol::Float64=1e-8)

	# initial w to x
	model_to_x!(x, pa)

	"""
	Unbounded LBFGS inversion, only for testing
	"""
	res = optimize(df, x, 
		       LBFGS(),
		       Optim.Options(g_tol = g_tol, f_tol=f_tol, x_tol=x_tol,
		       iterations = 200, store_trace = store_trace,
		       extended_trace=extended_trace, show_trace = false))
	# pa.verbose && println(res)

	x_to_model!(Optim.minimizer(res), pa)

	return res
end

function update_gf!(pa, xgf, last_xgf, dfgf)
	pa.attrib_inv=:gf    
	resgf = update!(pa, xgf, last_xgf, dfgf)
	fgf = Optim.minimum(resgf)
	return fgf
end

function update_wav!(pa, xwav, last_xwav, dfwav)
	pa.attrib_inv=:wav    
	reswav = update!(pa, xwav, last_xwav, dfwav)
	fwav = Optim.minimum(reswav)
	return fwav
end

#function remove_gfprecon(pa)
#	for i in eachindex(pa.gfprecon)
#		if(pa.gfprecon[i]≠0.0)
#			pa.gfprecon[i]=1.0
#		end
#	end
#end

"""
* re_init_flag :: re-initialize inversions with random input or not?
"""
function update_all!(pa; rtol=1e-3, tol=1e-3, maxitr_all=10, re_init_flag=true, maxitr=1000)


	converged_all=false
	itr_all=0

	pa.verbose && println("Blind Decon")  

	while ((!converged_all && itr_all < maxitr_all))
		itr_all += 1
		pa.verbose && println("===============================================================================")  
		pa.verbose && (itr_all > 1) && println("\t", "failed to converge.. reintializing (",itr_all,"/",maxitr_all,")")

	
		re_init_flag && initialize!(pa)

		fgfmin=Inf
		fgfmax=0.0
		fwavmin=Inf
		fwavmax=0.0

		itr=0
		converged=false

		pa.verbose && println("trip\t|\tfgf(<",rtol,"?)\t|\tfwav(<",rtol,"?)\t|\tconvergence(<",tol,"?)",) 
		while !converged && itr < maxitr

			pa.func_max_save[:]=0.0 # reset functional scales

			itr += 1
			fgf = update_gf!(pa, pa.xgf, pa.last_xgf, pa.dfgf)

			fgfmin = min(fgfmin, fgf) 
			fgfmax = max(fgfmax, fgf) 

			fwav = update_wav!(pa, pa.xwav, pa.last_xwav, pa.dfwav)
			fwavmin = min(fwavmin, fwav) 
			fwavmax = max(fwavmax, fwav) 

			rf = abs(fwav-fgf)/fwav
			pa.verbose && @printf("%d\t|\t%0.6e\t|\t%0.6e\t|\t%0.6e\t\n",itr, fgfmin/fgfmax, fwavmin/fwavmax,rf)

			converged = (rf < tol)
		end
		converged_all = (fgfmin/fgfmax < rtol) & (fwavmin/fwavmax < rtol)
	end
end


function initialize!(pa)
	# starting models
	randn!(pa.wav)
	randn!(pa.gf)
end


end # module
