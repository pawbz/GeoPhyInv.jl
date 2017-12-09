# blind deconvolution
__precompile__()

module Decon
using DSP
import JuMIT.DSP
import JuMIT.Inversion
import JuMIT.Misfits
import JuMIT.Inversion
using Optim, LineSearches
using RecipesBase
using StatsBase

type Param
	gfobs::Array{Float64,2}
	gf::Array{Float64,2}
	gfprecon::Array{Float64,2}
	gfweights::Array{Float64,2}
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
	wavnorm_flag::Bool 			# restrict wav along a unit circle during optimization
	wavnormmat::Matrix{Float64}             # stored outer product of wav
	dwavnorm::Vector{Float64}		# gradient w.r.t. normalized wavelet
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
end

@userplot Plot


@recipe function f(p::Plot, rvec=nothing)
	pa=p.args[1]
	(rvec===nothing) && (rvec=1:pa.nr)

	# time vectors
	# autocorr wav
	awavobs=autocor(pa.wavobs, 1:pa.nt-1, demean=true)
	awav=autocor(pa.wav, 1:pa.nt-1, demean=true)
	wavli=max(maximum(abs,awavobs), maximum(abs,awav))
	# autocorr gf 
	agfobs=autocor(pa.gfobs,1:pa.ntgf-1, demean=true)
	agf=autocor(pa.gf,1:pa.ntgf-1, demean=true)
	gfli=max(maximum(abs,agfobs), maximum(abs,agf))

	# cut receivers
	dcal=pa.dcal[:,rvec]
	dobs=pa.dobs[:,rvec]

	layout := (5,3)

	@series begin        
		subplot := 1
#		aspect_ratio := :auto
		legend := false
		pa.wavobs
	end
	@series begin        
		subplot := 2
		legend := false
		pa.gfobs
	end
#
#	@series begin        
#		subplot := 3
#		legend := false
#		pa.dobs
#	end
#
#	@series begin        
#		subplot := 4
#		legend := false
#		pa.wav
#	end
#	@series begin        
#		subplot := 5
#		legend := false
#		pa.gf[:,rvec]
#	end
#
#	@series begin        
#		subplot := 6
#		legend := false
#		pa.dcal[:,rvec]
#	end
#
#	@series begin        
#		subplot := 7
#		legend := false
#		awavobs
#
#		
#	end
#	@series begin        
#		subplot := 8
#		legend := false
#		agfobs[:, rvec]
#	end
#
#	@series begin        
#		subplot := 9
#		legend := false
#		
#	end
#
#	@series begin        
#		subplot := 10
#		legend := false
#		awav
#
#		
#	end
#	@series begin        
#		subplot := 11
#		legend := false
#		agf[:,rvec]
#	end
#
#	@series begin        
#		subplot := 12
#		legend := false
#		pa.dcal[:,rvec]-pa.dobs[:,rvec]
#	end
#
#
#
#
#	@series begin        
#		subplot := 13
#		aspect_ratio := :equal
#		seriestype := :histogram2d
#		title := "Scatter Wav"
#		legend := false
#		awavobs, awav
#	end
#
#	@series begin        
#		subplot := 14
#		aspect_ratio := :equal
#		seriestype := :histogram2d
#		title := "Scatter Wav"
#		legend := false
#		agfobs, agf
#	end
#
#	@series begin        
#		subplot := 15
#		aspect_ratio := 1
#		seriestype := :histogram2d
#		title := "Scatter Wav"
#		legend := false
#		pa.dobs, pa.dcal
#	end
#


end

"""
`gfprecon` : a preconditioner applied to each Greens functions [ntgf]
"""
function Param(ntgf, nt, nr; 
	       gfprecon=nothing,
	       gfweights=nothing,
	       gfoptim=nothing,
	       gfαvec=nothing,
	       wavoptim=nothing,
	       wavαvec=nothing,
	       wavprecon=nothing,
	       wavnorm_flag=false,
	       fft_threads=false,
	       dobs=nothing, gfobs=nothing, wavobs=nothing, verbose=false, attrib_inv=:gf,) 
	(dobs===nothing) && (dobs=zeros(nt, nr))
	dcal=zeros(nt, nr)
	ddcal=zeros(nt, nr)
	
	# initial values are random
	wav=randn(nt)
	wavmat=repeat(wav,outer=(1,nr))
	dwavmat=similar(wavmat)
	gf=zeros(ntgf,nr)
	dgf=similar(gf)

	# use maximum threads for fft
	fft_threads &&  (FFTW.set_num_threads(Sys.CPU_CORES))
	# fft dimension for plan
	np2=nextfastfft(maximum([2*nt, 2*ntgf]))
	# create plan
	fftplan = plan_fft!(complex.(zeros(np2,nr)),[1],flags=FFTW.MEASURE)
	ifftplan = plan_ifft!(complex.(zeros(np2,nr)),[1],flags=FFTW.MEASURE)

	# preallocate pow2 vectors
	spow2=complex.(zeros(np2,nr))
	rpow2=complex.(zeros(np2,nr))
	wpow2=complex.(zeros(np2,nr))

	
	# convert initial model to the inversion variable
	xgf = zeros(ntgf*nr);
	last_xgf = randn(size(xgf)) # reset last_x

	# dummy dfgf for optimization, to be updated later
	dfgf = OnceDifferentiable(x -> randn(),   (storage, x) -> randn())

	xwav = zeros(nt);
	last_xwav = randn(size(xwav)) # reset last_x

	# dummy, to be updated later
	dfwav = OnceDifferentiable(x -> randn(),   (storage, x) -> randn())

	# create dummy gfobs if necessary
	(gfobs===nothing) && (gfobs=zeros(gf))
	# create dummy wavobs if necessary
	(wavobs===nothing) && (wavobs=zeros(wav))

	# create gf precon
	(gfprecon===nothing) && (gfprecon=ones(gf))
	(gfweights===nothing) && (gfweights=ones(gf))

	# create gf precon
	(wavprecon===nothing) && (wavprecon=ones(wav))

	wavnorm_flag ?	(wavnormmat=zeros(nt, nt)) : (wavnormmat=zeros(1,1))
	wavnorm_flag ?	(dwavnorm=zeros(nt)) : (dwavnorm=zeros(1))

	pa = Param(gfobs, gf, gfprecon, gfweights, dgf, ntgf, dobs, dcal, ddcal, wavobs, wav, wavprecon, wavmat, dwavmat,
	    wavnorm_flag, wavnormmat, dwavnorm,
	     spow2, rpow2, wpow2,
	     nt, nr, attrib_inv, verbose, fftplan, ifftplan,
	     xgf, last_xgf, dfgf,
	     xwav, last_xwav, dfwav)
 
	if iszero(pa.dobs) 
		((gfobs===nothing) | (wavobs===nothing)) && error("need gfobs and wavobs")
		pa.attrib_inv=:wav
		F!(pa, 
     		pa.xwav,  # dummy
     		mattrib=:gfwav, dattrib=:dobs, gf=gfobs, wav=wavobs)
	end

	initialize!(pa)
	update_func_grad!(pa,gfoptim=gfoptim,wavoptim=wavoptim,gfαvec=gfαvec,wavαvec=wavαvec)

	return pa
	
end

function update_func_grad!(pa; gfoptim=nothing, wavoptim=nothing, gfαvec=nothing, wavαvec=nothing)
	# they will be changed in this program, so make a copy 
	wavsave=copy(pa.wav);
	gfsave=copy(pa.gf);
	dcalsave=copy(pa.dcal);

	(gfoptim===nothing) && (gfoptim=[:ls])
	(gfαvec===nothing) && (gfαvec=ones(length(gfoptim)))

	(wavoptim===nothing) && (wavoptim=[:ls])
	(wavαvec===nothing) && (wavαvec=ones(length(wavoptim)))

	# dfgf for optimization functions
	optim_funcgf=Vector{Function}(length(gfoptim))
	optim_gradgf=Vector{Function}(length(gfoptim))
	for iop in 1:length(gfoptim)
		if (gfoptim[iop]==:ls)
			optim_funcgf[iop]= x->func_grad!(nothing, x,  pa) 
			optim_gradgf[iop]=(storage, x)->func_grad!(storage, x,  pa)
		elseif(gfoptim[iop]==:weights)
			optim_funcgf[iop]= x -> func_grad_gf_weights!(nothing, x, pa) 
			optim_gradgf[iop]= (storage, x) -> func_grad_gf_weights!(storage, x, pa)
		else
			error("invalid optim_funcgf")
		end
	end
	pa.attrib_inv=:gf
	# multi-objective framework
	paMOgf=Inversion.ParamMO(noptim=length(gfoptim), ninv=length(pa.xgf), αvec=gfαvec,
			    		optim_func=optim_funcgf,optim_grad=optim_gradgf,
					x_init=randn(length(pa.xgf),10))
	# create dfgf for optimization
	pa.dfgf = OnceDifferentiable(x -> paMOgf.func(x, paMOgf),       
			    (storage, x) -> paMOgf.grad!(storage, x, paMOgf))


	# dfwav for optimization functions
	optim_funcwav=Vector{Function}(length(wavoptim))
	optim_gradwav=Vector{Function}(length(wavoptim))
	for iop in 1:length(wavoptim)
		if (wavoptim[iop]==:ls)
			optim_funcwav[iop]=x->func_grad!(nothing, x,  pa) 
			optim_gradwav[iop]=(storage, x)->func_grad!(storage, x,  pa) 
		else
			error("invalid optim_funcwav")
		end
	end

	pa.attrib_inv=:wav
	# multi-objective framework
	paMOwav=Inversion.ParamMO(noptim=length(wavoptim), ninv=pa.nt, αvec=wavαvec,
			    		optim_func=optim_funcwav,optim_grad=optim_gradwav,
					x_init=randn(length(pa.xwav),10))
	pa.dfwav = OnceDifferentiable(x -> paMOwav.func(x, paMOwav),         
			    (storage, x) -> paMOwav.grad!(storage, x, paMOwav))


	copy!(pa.wav, wavsave)
	copy!(pa.gf, gfsave)
	copy!(pa.dcal,dcalsave)

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
	fwav = Misfits.error_after_normalized_autocor(pa.wav, pa.wavobs)
	fgf = Misfits.error_after_normalized_autocor(pa.gf, pa.gfobs)
	f = Misfits.error_squared_euclidean!(nothing, pa.dcal, pa.dobs, nothing, norm_flag=true)

	println("Blind Decon\t")
	println("===========")
	println("error in estimated wavelet:\t", fwav)
	println("error after autocor in estimated Green Functions:\t", fgf)
	println("normalized error in the data:\t", f)

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
		pa.wavnorm_flag && (scale!(pa.wav, inv(xn)))
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
			x::AbstractVector{Float64};
			gf=nothing,
			wav=nothing,
			mattrib::Symbol=:x,
			dattrib::Symbol=:dcal,
			   )
	if(pa.attrib_inv==:wav)
		compute=(x!=pa.last_xwav)
	elseif(pa.attrib_inv==:gf)
		compute=(x!=pa.last_xgf)
	else
		compute=false
	end

	if(compute)

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
		if(pa.attrib_inv==:wav)
			copy!(pa.last_xwav, x)
		elseif(pa.attrib_inv==:gf)
			copy!(pa.last_xgf, x)
		end


		for ir in 1:pa.nr
			for it in 1:pa.nt
				pa.wavmat[it,ir]=wav[it]
			end
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

function func_grad!(storage, x::AbstractVector{Float64},
			pa)


	#pa.verbose && println("computing gradient...")

	# x to w 
	x_to_model!(x, pa)

	F!(pa, x)

	if(storage === nothing)
		# compute misfit and δdcal
		f = Misfits.error_squared_euclidean!(nothing, pa.dcal, pa.dobs, nothing, norm_flag=true)
	else
		f = Misfits.error_squared_euclidean!(pa.ddcal, pa.dcal, pa.dobs, nothing, norm_flag=true)
		Fadj!(pa, x, storage, pa.ddcal)
	end
	return f

end


# add model based constraints here

# all the greens' functions have to be correlated

# exponential-weighted norm for the green functions
function func_grad_gf_weights!(storage, x, pa)
	x_to_model!(x, pa)
	!(pa.attrib_inv == :gf) && error("only for gf inversion")
	if(!(storage === nothing)) #
		f = Misfits.error_weighted_norm!(pa.dgf,pa.gf, pa.gfweights) #
		for i in eachindex(storage)
			storage[i]=pa.dgf[i]
		end
	else	
		f = Misfits.error_weighted_norm!(nothing,pa.gf, pa.gfweights)
	end
	return f
end

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
		for i in eachindex(storage)
			if(iszero(pa.wavprecon[i]))
				storage[i]=0.0
			else
				storage[i] = storage[i]/pa.wavprecon[i]
			end
		end
		if(pa.wavnorm_flag)
			copy!(pa.dwavnorm, storage)
			Misfits.derivative_vector_magnitude!(storage,pa.dwavnorm,x,pa.wavnormmat)
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
function update!(pa::Param, x, df; store_trace::Bool=false, 
		 extended_trace::Bool=false, 
	     f_tol::Float64=1e-8, g_tol::Float64=1e-30, x_tol::Float64=1e-30)

	# initial w to x
	model_to_x!(x, pa)

	"""
	Unbounded LBFGS inversion, only for testing
	"""
	res = optimize(df, x, 
		       LBFGS(),
		       Optim.Options(g_tol = g_tol, f_tol=f_tol, x_tol=x_tol,
		       iterations = 2000, store_trace = store_trace,
		       extended_trace=extended_trace, show_trace = false))
	pa.verbose && println(res)

	x_to_model!(Optim.minimizer(res), pa)

	return res
end

function update_gf!(pa, xgf,  dfgf)
	pa.attrib_inv=:gf    
	resgf = update!(pa, xgf,  dfgf)
	fgf = Optim.minimum(resgf)
	return fgf
end

function update_wav!(pa, xwav, dfwav)
	pa.attrib_inv=:wav    
	reswav = update!(pa, xwav, dfwav)
	fwav = Optim.minimum(reswav)
	return fwav
end

function remove_gfprecon!(pa)
	for i in eachindex(pa.gfprecon)
		if(pa.gfprecon[i]≠0.0)
			pa.gfprecon[i]=1.0
		end
	end
end

"""
* re_init_flag :: re-initialize inversions with random input or not?
"""
function update_all!(pa; max_roundtrips=100, max_reroundtrips=10, ParamAM_func=nothing, roundtrip_tol=1e-3)

	if(ParamAM_func===nothing)
		ParamAM_func=x->Inversion.ParamAM(x, optim_tols=[1e-5, 1e-5],name="Blind Decon",
				    roundtrip_tol=roundtrip_tol, max_roundtrips=max_roundtrips,
				    max_reroundtrips=max_reroundtrips,
				    min_roundtrips=10,
				    reinit_func=x->initialize!(pa))
	end

	
	# create alternating minimization parameters
	f1=x->update_wav!(pa, pa.xwav,  pa.dfwav)
	f2=x->update_gf!(pa, pa.xgf, pa.dfgf)
	paam=ParamAM_func([f1, f2])

	# do inversion
	Inversion.go(paam)

	# print errors
	error(pa)
end


function initialize!(pa)
	# starting random models
	randn!(pa.wav)
	randn!(pa.gf)
end





"""
Create preconditioners using the observed Green Functions.
* `cflag` : causaulity flag
"""
function create_weights(ntgf, nt, gfobs; αexp=0.0, cflag=true)

	nr=size(gfobs,2)
	wavprecon=ones(nt)
	gfprecon=ones(ntgf, nr); 
	gfweights=ones(ntgf, nr); 
	minindz=ntgf
	gfweights=ones(ntgf, nr)
	for ir in 1:nr
		gf=normalize(view(gfobs,:,ir))
		indz=findfirst(x->abs(x)>1e-6, gf)
		if(cflag)
			if(indz==0)
				gfprecon[:,ir]=0.0    
				gfweights[:,ir]=0.0
			else
				gfprecon[1:indz-1,ir]=0.0    
				gfweights[1:indz-1,ir]=0.0
#				wavprecon[end-indz+1:end]=0.0
			end
		end
		if(indz≠0)
			for i in indz+1:ntgf        
				gfweights[i,ir]=exp(αexp*(i-indz-1)/ntgf)  # exponential weights
				gfprecon[i,ir]=exp(αexp*(i-indz-1)/ntgf)  # exponential weights
			end
		end
	end

	return gfprecon, gfweights, wavprecon

end

end # module
