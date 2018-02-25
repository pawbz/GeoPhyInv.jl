# blind deconvolution
__precompile__()

module DeConv
using DSP
import JuMIT.DSP
import JuMIT.Inversion
import JuMIT.Misfits
import JuMIT.Conv
import JuMIT.Grid
import JuMIT.Inversion
using Optim, LineSearches
using RecipesBase
using DataFrames
using StatsBase
using JLD
using CSV

mutable struct Param
	ntgf::Int64
	nt::Int64
	nr::Int64
	nra::Int64 # save actual nr here, as for mode==2: nr=nr*nr
	gfobs::Array{Float64, 2} # save, before modifying it in obs
	wavobs::Vector{Float64} # save, before modifying it in obs
	dobs::Array{Float64,2} # save before modifying it
	obs::Conv.Param{Float64,2} # observed convolutional model
	cal::Conv.Param{Float64,2} # calculated convolutional model
	calsave::Conv.Param{Float64,2} # save the best result
	dgf::Array{Float64,2}
	dwav::Array{Float64,2}
	ddcal::Array{Float64,2}
	gfprecon::Array{Float64,2}
	gfpreconI::Array{Float64,2}
	gfweights::Array{Float64,2}
	wavprecon::Vector{Float64}
	wavpreconI::Vector{Float64}
	wavnorm_flag::Bool 			# restrict wav along a unit circle during optimization
	wavnormmat::Matrix{Float64}             # stored outer product of wav
	dwavnorm::Vector{Float64}		# gradient w.r.t. normalized wavelet
	attrib_inv::Symbol
	verbose::Bool
	xgf::Vector{Float64}
	last_xgf::Vector{Float64}
	gf_func::Function
	gf_grad!::Function
	xwav::Vector{Float64}
	last_xwav::Vector{Float64}
	wav_func::Function
	wav_grad!::Function
	err::DataFrames.DataFrame
	mode::Int32
	gf_acorr::Conv.Param{Float64,2}
	dgf_acorr::Array{Float64,2}
end



"""
`gfprecon` : a preconditioner applied to each Greens functions [ntgf]
"""
function Param(ntgf, nt, nr; 
	       mode=1,
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

	# use maximum threads for fft
	fft_threads &&  (FFTW.set_num_threads(Sys.CPU_CORES))


	# create models depending on mode
	nra=nr
	if(mode==1)
		obs=Conv.Param(ntwav=nt, ntd=nt, ntgf=ntgf, dims=(nr,), wavlags=[nt-1, 0])
		cal=deepcopy(obs)
	elseif(mode==2)
		nr=nr*nr
		obs=Conv.Param(ntwav=2*nt-1, ntd=2*nt-1, ntgf=2*ntgf-1, dims=(nr,), wavlags=[nt-1, nt-1], dlags=[nt-1, nt-1], gflags=[ntgf-1, ntgf-1])
		cal=deepcopy(obs)
	end
	calsave=deepcopy(cal);

	# initial values are random
	dwav=zeros(cal.wav)
	dgf=zeros(cal.gf)
	ddcal=zeros(cal.d)
	
	# inversion variables allocation
	xgf = zeros(length(cal.gf));
	if(mode==1)
		xwav = zeros(nt);
	elseif(mode==2)
		xwav = zeros(nt-1);
	end
	last_xgf = randn(size(xgf)) # reset last_x
	last_xwav = randn(size(xwav)) # reset last_x

	wavnorm_flag ?	(wavnormmat=zeros(nt, nt)) : (wavnormmat=zeros(1,1))
	wavnorm_flag ?	(dwavnorm=zeros(nt)) : (dwavnorm=zeros(1))

	err=DataFrame(gf=[], gf_nodecon=[], wav=[],d=[])

	gf_acorr=Conv.Param(ntgf=ntgf, ntd=ntgf, ntwav=2*ntgf-1, dims=(nra,), wavlags=[ntgf-1, ntgf-1])
	dgf_acorr=zeros(2*ntgf-1, nra)

	pa=Param(ntgf,nt,nr,
		 nra,
	  zeros(ntgf,nra), zeros(nt),
	  zeros(nt, nra),
	  obs,cal,calsave, dgf,dwav,ddcal,
	  zeros(cal.gf), zeros(cal.gf),
	  zeros(cal.gf),
	  zeros(xwav), zeros(xwav),
	  wavnorm_flag,wavnormmat,
	  dwavnorm,attrib_inv,verbose,xgf,last_xgf,x->randn(),x->randn(),xwav,last_xwav,x->randn(),x->randn(), err,
	  mode, gf_acorr, dgf_acorr)

	add_gfprecon!(pa, gfprecon)
	add_gfweights!(pa, gfweights)
	add_wavprecon!(pa, wavprecon)
 
	if(!(gfobs===nothing))
		for i in eachindex(pa.gfobs)
			pa.gfobs[i]=gfobs[i]
		end# save gfobs, before modifying
	end

	if(!(wavobs===nothing))
		for i in eachindex(pa.wavobs)
			pa.wavobs[i]=wavobs[i]
		end# save gfobs, before modifying
	end

	if(!(dobs===nothing))
		for i in eachindex(pa.dobs)
			pa.dobs[i]=dobs[i]
		end
	else # otherwise perform modelling
		(iszero(pa.gfobs) || iszero(pa.wavobs)) && error("need gfobs and wavobs")
		obstemp=Conv.Param(ntwav=nt, ntd=nt, ntgf=ntgf, dims=(nra,), wavlags=[nt-1, 0])
		copy!(obstemp.gf, pa.gfobs)
		copy!(obstemp.wav, repmat(pa.wavobs,1,nra))
		Conv.mod!(obstemp, :d) # model observed data
		copy!(pa.dobs, obstemp.d)
	end

	if(mode==1)
		gfobs=pa.gfobs
		wavobs=pa.wavobs
		dobs=pa.dobs
	elseif(mode==2)
		gfobs=hcat(Conv.xcorr(pa.gfobs)...)
		wavobs=hcat(Conv.xcorr(pa.wavobs)...)
		dobs=hcat(Conv.xcorr(pa.dobs)...) # do a cross-correlation 
	end

	# obs.gf <-- gfobs
	copy!(pa.obs.gf, gfobs)
	# obs.wav <-- wavobs
	replace_obswav!(pa, wavobs)
	# obs.d <-- dobs
	copy!(pa.obs.d, dobs) 

	initialize!(pa)
	update_func_grad!(pa,gfoptim=gfoptim,wavoptim=wavoptim,gfαvec=gfαvec,wavαvec=wavαvec)

	return pa
	
end

function add_gfprecon!(pa::Param, gfprecon=nothing)
	# create gf precon
	(gfprecon===nothing) && (gfprecon=ones(pa.cal.gf))

	copy!(pa.gfprecon, gfprecon)		
	copy!(pa.gfpreconI, pa.gfprecon)
	for i in eachindex(gfprecon)
		if(!(iszero(gfprecon[i])))
			pa.gfpreconI[i]=inv(pa.gfprecon[i])
		end
	end
end
function add_wavprecon!(pa::Param, wavprecon=nothing)
	# create gf precon
	(wavprecon===nothing) && (wavprecon=ones(pa.xwav))
	copy!(pa.wavprecon, wavprecon)
	copy!(pa.wavpreconI, pa.wavprecon)
	for i in eachindex(wavprecon)
		if(!(iszero(wavprecon[i])))
			pa.wavpreconI[i]=inv(pa.wavprecon[i])
		end
	end
end
function add_gfweights!(pa::Param, gfweights=nothing)
	(gfweights===nothing) && (gfweights=ones(pa.cal.gf))
	copy!(pa.gfweights, gfweights)
end



function update_func_grad!(pa; gfoptim=nothing, wavoptim=nothing, gfαvec=nothing, wavαvec=nothing)
	# they will be changed in this program, so make a copy 
	wavsave=copy(pa.cal.wav);
	gfsave=copy(pa.cal.gf);
	dcalsave=copy(pa.cal.d);

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
		elseif(gfoptim[iop]==:acorr_weights)
			optim_funcgf[iop]= x -> func_grad_gf_acorr_weights!(nothing, x, pa) 
			optim_gradgf[iop]= (storage, x) -> func_grad_gf_acorr_weights!(storage, x, pa)
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
	pa.gf_func = x -> paMOgf.func(x, paMOgf)
	pa.gf_grad! = (storage, x) -> paMOgf.grad!(storage, x, paMOgf)
#	pa.dfgf = OnceDifferentiable(x -> paMOgf.func(x, paMOgf),       
#			    (storage, x) -> paMOgf.grad!(storage, x, paMOgf), )


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
	paMOwav=Inversion.ParamMO(noptim=length(wavoptim), ninv=length(pa.xwav), αvec=wavαvec,
			    		optim_func=optim_funcwav,optim_grad=optim_gradwav,
					x_init=randn(length(pa.xwav),10))
#	pa.dfwav = OnceDifferentiable(x -> paMOwav.func(x, paMOwav),         
#			    (storage, x) -> paMOwav.grad!(storage, x, paMOwav))
	pa.wav_func = x -> paMOwav.func(x, paMOwav)
	pa.wav_grad! =  (storage, x) -> paMOwav.grad!(storage, x, paMOwav)


	copy!(pa.cal.wav, wavsave)
	copy!(pa.cal.gf, gfsave)
	copy!(pa.cal.d,dcalsave)

	return pa
	
end


function ninv(pa)
	if(pa.attrib_inv == :wav)
		return pa.nt
	else(pa.attrib_inv == :gf)
		return pa.ntgf*pa.nr
	end
end

"""
compute errors
update pa.err
print?
give either cal or calsave?
"""
function err!(pa::Param; cal=pa.cal) 
	xgf_nodecon=hcat(Conv.xcorr(pa.dobs, lags=[pa.ntgf-1, pa.ntgf-1])...)
	xgfobs=hcat(Conv.xcorr(pa.gfobs)...) # compute xcorr with reference gf
	if(pa.mode==1) 
		fwav = Misfits.error_after_normalized_autocor(cal.wav, pa.obs.wav)
		xgfcal=hcat(Conv.xcorr(cal.gf)...) # compute xcorr with reference gf
		fgf = Misfits.error_squared_euclidean!(nothing, xgfcal, xgfobs, nothing, norm_flag=true)
	elseif(pa.mode==2)
		fgf = Misfits.error_squared_euclidean!(nothing, cal.gf, pa.obs.gf, nothing, norm_flag=true)
		fwav = Misfits.error_squared_euclidean!(nothing, cal.wav, pa.obs.wav, nothing, norm_flag=true)
	end
	fgf_nodecon = Misfits.error_squared_euclidean!(nothing, xgf_nodecon, xgfobs, nothing, norm_flag=true)
	f = Misfits.error_squared_euclidean!(nothing, cal.d, pa.obs.d, nothing, norm_flag=true)

	push!(pa.err[:wav],fwav)
	push!(pa.err[:d],f)
	push!(pa.err[:gf],fgf)
	push!(pa.err[:gf_nodecon],fgf_nodecon)
	println("Blind Decon Errors\t")
	println("==================")
	show(pa.err)
end 


function model_to_x!(x, pa)
	if(pa.attrib_inv == :wav)
		if(pa.mode==1)
			for i in eachindex(x)
				x[i]=pa.cal.wav[i,1]*pa.wavprecon[i] # just take any one receiver
			end
		elseif(pa.mode==2)
			for i in eachindex(x)
				x[i]=pa.cal.wav[i+pa.nt,1]*pa.wavprecon[i] # just take any one receiver and only positive lags
			end
		end
	else(pa.attrib_inv == :gf)
		for i in eachindex(x)
			x[i]=pa.cal.gf[i]*pa.gfprecon[i] 		# multiply by gfprecon
		end
	end
	return x
end


function x_to_model!(x, pa)
	if(pa.attrib_inv == :wav)
		if(pa.mode==1)
			for j in 1:pa.nr
				for i in 1:pa.nt
					# put same in all receivers
					pa.cal.wav[i,j]=x[i]*pa.wavpreconI[i]
				end
			end
			if(pa.wavnorm_flag)
				xn=vecnorm(x)
				scale!(pa.cal.wav, inv(xn))
			end
		elseif(pa.mode==2)
			for j in 1:pa.nr
				pa.cal.wav[pa.nt,j]=1.0
				for i in 1:pa.nt-1
					# put same in all receivers
					pa.cal.wav[pa.nt+i,j]=x[i]*pa.wavpreconI[i]
					# put same in negative lags
					pa.cal.wav[pa.nt-i,j]=x[i]*pa.wavpreconI[i]
				end
			end
		end
	else(pa.attrib_inv == :gf)
		for i in eachindex(pa.cal.gf)
			pa.cal.gf[i]=x[i]*pa.gfpreconI[i]
		end
	end
	return pa
end

function F!(pa::Param,	x::AbstractVector{Float64}  )
	if(pa.attrib_inv==:wav)
		compute=(x!=pa.last_xwav)
	elseif(pa.attrib_inv==:gf)
		compute=(x!=pa.last_xgf)
	else
		compute=false
	end

	if(compute)

		x_to_model!(x, pa) # modify pa.cal.wav or pa.cal.gf

		#pa.verbose && println("updating buffer")
		if(pa.attrib_inv==:wav)
			copy!(pa.last_xwav, x)
		elseif(pa.attrib_inv==:gf)
			copy!(pa.last_xgf, x)
		end


		Conv.mod!(pa.cal, :d) # modify pa.cal.d
		return pa
	end
end

function func_grad!(storage, x::AbstractVector{Float64},pa)

	# x to pa.cal.wav or pa.cal.gf 
	x_to_model!(x, pa)

	F!(pa, x) # forward

	if(storage === nothing)
		# compute misfit and δdcal
		f = Misfits.error_squared_euclidean!(nothing, pa.cal.d, pa.obs.d, nothing, norm_flag=true)
	else
		f = Misfits.error_squared_euclidean!(pa.ddcal, pa.cal.d, pa.obs.d, nothing, norm_flag=true)
		Fadj!(pa, x, storage, pa.ddcal)
	end
	return f

end

"""
update calsave only when error in d is low
"""
function update_calsave!(pa)
	f1=Misfits.error_squared_euclidean!(nothing, pa.calsave.d, pa.obs.d, nothing, norm_flag=true)
	f2=Misfits.error_squared_euclidean!(nothing, pa.cal.d, pa.obs.d, nothing, norm_flag=true)
	if(f2<f1)
		copy!(pa.calsave.d, pa.cal.d)
		copy!(pa.calsave.gf, pa.cal.gf)
		copy!(pa.calsave.wav, pa.cal.wav)
	end
end

# add model based constraints here

# all the greens' functions have to be correlated

# exponential-weighted norm for the green functions
function func_grad_gf_weights!(storage, x, pa)
	x_to_model!(x, pa)
	!(pa.attrib_inv == :gf) && error("only for gf inversion")
	if(!(storage === nothing)) #
		f = Misfits.error_weighted_norm!(pa.dgf,pa.cal.gf, pa.gfweights) #
		for i in eachindex(storage)
			storage[i]=pa.dgf[i]
		end
	else	
		f = Misfits.error_weighted_norm!(nothing,pa.cal.gf, pa.gfweights)
	end
	return f
end

# exponential-weighted norm for the green functions
function func_grad_gf_acorr_weights!(storage, x, pa)
	x_to_model!(x, pa)
	!(pa.attrib_inv == :gf) && error("only for gf inversion")

	if(!(storage === nothing)) #
		f = Misfits.error_acorr_weighted_norm!(pa.dgf,pa.cal.gf, 
					 paconv=pa.gf_acorr,dfdwav=pa.dgf_acorr) #
		for i in eachindex(storage)
			storage[i]=pa.dgf[i]
		end
	else	
		f = Misfits.error_acorr_weighted_norm!(nothing,pa.cal.gf, 
					 paconv=pa.gf_acorr,dfdwav=pa.dgf_acorr)
	end

	return f
end
#  



"""
Apply Fadj to 
x is not used?
"""
function Fadj!(pa, x, storage, dcal)
	storage[:] = 0.
	if(pa.attrib_inv == :wav)
		Conv.mod!(pa.cal, :wav, d=dcal, wav=pa.dwav)
		if(pa.mode==1)
			# stack ∇wav along receivers
			for i in 1:size(pa.dwav,2)             
				for j in 1:size(pa.dwav,1)
					storage[j] += pa.dwav[j,i]
				end
			end
		elseif(pa.mode==2)
			for i in 1:size(pa.dwav,2)             
				for j in 1:pa.nt-1
					storage[j] += pa.dwav[pa.nt-j,i] # -ve lags
					storage[j] += pa.dwav[pa.nt+j,i] # -ve lags
				end
			end
		end

		# apply precon
		for i in eachindex(storage)
			if(iszero(pa.wavprecon[i]))
				storage[i]=0.0
			else
				storage[i] = storage[i]/pa.wavprecon[i]
			end
		end
		# factor, because wav was divided by norm of x
		if(pa.wavnorm_flag)
			copy!(pa.dwavnorm, storage)
			Misfits.derivative_vector_magnitude!(storage,pa.dwavnorm,x,pa.wavnormmat)
		end

	else(pa.attrib_inv == :gf)
		Conv.mod!(pa.cal, :gf, gf=pa.dgf, d=dcal)
		copy!(storage, pa.dgf) # remove?

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
function update!(pa::Param, x, f, g!; 
		 store_trace::Bool=false, 
		 extended_trace::Bool=false, 
	     f_tol::Float64=1e-8, g_tol::Float64=1e-30, x_tol::Float64=1e-30)

	# initial w to x
	model_to_x!(x, pa)

	"""
	Unbounded LBFGS inversion, only for testing
	"""
	res = optimize(f, g!, x, 
		       LBFGS(),
		       Optim.Options(g_tol = g_tol, f_tol=f_tol, x_tol=x_tol,
		       iterations = 2000, store_trace = store_trace,
		       extended_trace=extended_trace, show_trace = false))
	pa.verbose && println(res)

	x_to_model!(Optim.minimizer(res), pa)

	return res
end

function replace_obswav!(pa::Param, wav)
	if(length(wav)==length(pa.obs.wav))
		copy!(pa.obs.wav, wav)
	elseif(length(wav)==size(pa.obs.wav,1))
		for j in 1:pa.nr, i in 1:size(pa.obs.wav,1)
				pa.obs.wav[i,j]=wav[i]
		end
	else
		error("invalid input wav")
	end
	# observed data changes as well :)
	Conv.mod!(pa.obs,:d)
end

function update_gf!(pa, xgf)
	pa.attrib_inv=:gf    
	resgf = update!(pa, xgf,  pa.gf_func, pa.gf_grad!)
	fgf = Optim.minimum(resgf)
	return fgf
end

function update_wav!(pa, xwav)
	pa.attrib_inv=:wav    
	reswav = update!(pa, xwav, pa.wav_func, pa.wav_grad!)
	fwav = Optim.minimum(reswav)
	return fwav
end


"""
Remove preconditioners from pa
"""
function remove_gfprecon!(pa; all=false)
	for i in eachindex(pa.gfprecon)
		if((pa.gfprecon[i]≠0.0) || all)
			pa.gfprecon[i]=1.0
			pa.gfpreconI[i]=1.0
		end
	end
end

"""
Remove weights from pa
"""
function remove_gfweights!(pa; all=false)
	for i in eachindex(pa.gfweights)
		if((pa.gfweights[i]≠0.0) || all)
			pa.gfweights[i]=1.0
		end
	end
end


"""
* re_init_flag :: re-initialize inversions with random input or not?
"""
function update_all!(pa; max_roundtrips=100, max_reroundtrips=10, ParamAM_func=nothing, roundtrip_tol=1e-6,
		     optim_tols=[1e-6, 1e-6])

	if(ParamAM_func===nothing)
		ParamAM_func=x->Inversion.ParamAM(x, optim_tols=optim_tols,name="Blind Decon",
				    roundtrip_tol=roundtrip_tol, max_roundtrips=max_roundtrips,
				    max_reroundtrips=max_reroundtrips,
				    min_roundtrips=10,
				    reinit_func=x->initialize!(pa),
				    after_reroundtrip_func=x->(err!(pa); update_calsave!(pa);),
				    )
	end

	
	# create alternating minimization parameters
	f1=x->update_wav!(pa, pa.xwav)
	f2=x->update_gf!(pa, pa.xgf)
	paam=ParamAM_func([f1, f2])

	# do inversion
	Inversion.go(paam)

	# print errors
	err!(pa)
	println(" ")
end


function initialize!(pa)
	if(pa.mode==1)
		# starting random models
		for i in 1:pa.nt
			x=(pa.wavprecon[i]≠0.0) ? randn() : 0.0
			for j in 1:pa.nr
				pa.cal.wav[i,j]=x
			end
		end
	elseif(pa.mode==2)
		for i in 1:pa.nt-1
			x=(pa.wavprecon[i]≠0.0) ? randn() : 0.0
			for j in 1:pa.nr
				pa.cal.wav[pa.nt+i,j]=x*0.0
				pa.cal.wav[pa.nt-i,j]=x*0.0
			end
		end
		pa.cal.wav[pa.nt,:]=1.0

	end
	for i in eachindex(pa.cal.gf)
		x=(pa.gfprecon[i]≠0.0) ? randn() : 0.0
		pa.cal.gf[i]=x
	end
end





"""
Create preconditioners using the observed Green Functions.
* `cflag` : impose causaulity by creating gfprecon using gfobs
* `max_tfrac_gfprecon` : maximum length of precon windows on gf
"""
function create_weights(ntgf, nt, gfobs; αexp=0.0, cflag=true,
		       max_tfrac_gfprecon=1.0)
	
	ntgfprecon=round(Int,max_tfrac_gfprecon*ntgf);

	nr=size(gfobs,2)
	wavprecon=ones(nt)
	gfprecon=ones(ntgf, nr); 
	gfweights=ones(ntgf, nr); 
	minindz=ntgf
	gfweights=ones(ntgf, nr)
	for ir in 1:nr
		gf=normalize(view(gfobs,:,ir))
		indz=findfirst(x->abs(x)>1e-6, gf)
	#	if(indz > 1) 
	#		indz -= 1 # window one sample less than actual
	#	end
		if(!cflag && indz≠0)
			indz=1
		end
		if(indz≠0)
			for i in 1:indz-1
				gfprecon[i,ir]=0.0
				gfweights[i,ir]=0.0
			end
			for i in indz:indz+ntgfprecon
				if(i≤ntgf)
					gfweights[i,ir]=exp(αexp*(i-indz-1)/ntgf)  # exponential weights
					gfprecon[i,ir]=exp(αexp*(i-indz-1)/ntgf)  # exponential weights
				end
			end
			for i in indz+ntgfprecon+1:ntgf
				gfprecon[i,ir]=0.0
				gfweights[i,ir]=0.0
			end
		else
			gfprecon[:,ir]=0.0
			gfweights[:,ir]=0.0
		end
	end
	return gfprecon, gfweights, wavprecon
end

function create_white_weights(ntgf, nt, nr)
	wavprecon=ones(nt)
	gfprecon=ones(2*ntgf-1, nr*nr);
	gfweights=zeros(2*ntgf-1, nr*nr);
	for ir in 1:nr
		gfprecon[:,ir+(ir-1)*nr]=1.0
		for i in 1:ntgf-1    
			gfprecon[ntgf+i,ir+(ir-1)*nr]=0.0    
			gfprecon[ntgf-i,ir+(ir-1)*nr]=0.0
		end
		gfweights[ntgf,ir+(ir-1)*nr]=0.0
		for i in 1:ntgf-1
			gfweights[ntgf+i,ir+(ir-1)*nr]=i*2
			gfweights[ntgf-i,ir+(ir-1)*nr]=i*2
		end
	end
	return gfprecon, gfweights, wavprecon
end

@userplot Plot


@recipe function f(p::Plot; rvec=nothing, δt=1.0)
	pa=p.args[1]
	(rvec===nothing) && (rvec=1:pa.nr)

	gfxobs=zeros(2*size(pa.obs.gf,1)-1, size(pa.obs.gf,2))
	gfxcal=similar(gfxobs)
	scobs=vecnorm(pa.obs.gf[:,1])^2
	sccal=vecnorm(pa.cal.gf[:,1])^2
	for ir in 1:pa.nr
		gfxobs[:,ir] = xcorr(pa.obs.gf[:,1], pa.obs.gf[:, ir])/scobs
		gfxcal[:,ir] = xcorr(pa.cal.gf[:,1], pa.cal.gf[:, ir])/sccal
	end
	#
	gf=collect(1:pa.ntgf)*δt
	gfx=collect(-pa.ntgf+1:1:pa.ntgf-1)*δt

	# time vectors
	# autocorr wav
	awavobs=autocor(pa.obs.wav[:,1], 1:pa.nt-1, demean=true)
	awav=autocor(pa.cal.wav[:,1], 1:pa.nt-1, demean=true)
	wavli=max(maximum(abs,awavobs), maximum(abs,awav))
	# autocorr gf 
	agfobs=autocor(pa.obs.gf,1:pa.ntgf-1, demean=true)
	agf=autocor(pa.cal.gf,1:pa.ntgf-1, demean=true)
	gfli=max(maximum(abs,agfobs), maximum(abs,agf))

	nwav=length(awavobs)
	fact=(nwav>1000) ? round(Int,nwav/1000) : 1

	awavobs=awavobs[1:fact:nwav] # resample
	awav=awav[1:fact:nwav] # resample

	fact=(pa.nt*pa.nr>1000) ? round(Int,pa.nt*pa.nr/1000) : 1
	# cut receivers
#	dcal=pa.cal.d[1:fact:pa.nt*pa.nr]
#	dobs=pa.obs.d[1:fact:pa.nt*pa.nr]

	layout := (3,2)

	@series begin        
		subplot := 1
#		aspect_ratio := :auto
		legend := false
		l := :stem
		title := "Observed GF"
		w := 3
		gf, pa.obs.gf
	end
	@series begin        
		subplot := 2
#		aspect_ratio := :auto
		legend := false
		l := :stem
		title := "Modelled GF"
		w := 3
		gf, pa.cal.gf
	end
	@series begin        
		subplot := 3
#		aspect_ratio := :auto
		legend := false
		l := :stem
		title := "Observed GFX"
		w := 3
		gfx, gfxobs
	end
	@series begin        
		subplot := 4
#		aspect_ratio := :auto
		legend := false
		title := "Calculated GFX"
		l := :stem
		w := 3
		gfx, gfxcal
	end


	@series begin        
		subplot := 5
		aspect_ratio := :equal
		seriestype := :scatter
		title := "Scatter GF"
		legend := false
		gfxobs, gfxcal
	end

	@series begin        
		subplot := 6
		aspect_ratio := :equal
		seriestype := :scatter
		title := "Scatter Wav"
		legend := false
		awavobs, awav
	end


#
#	@series begin        
#		subplot := 3
#		legend := false
#		pa.obs.d
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
#		pa.cal.gf[:,rvec]
#	end
#
#	@series begin        
#		subplot := 6
#		legend := false
#		pa.cal.d[:,rvec]
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
#		pa.cal.d[:,rvec]-pa.obs.d[:,rvec]
#	end
#
#
#
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
#		pa.obs.d, pa.cal.d
#	end
#
end



"""
Save Param
"""
function save(pa::Param, folder; tgridgf=nothing, tgrid=nothing)
	!(isdir(folder)) && error("invalid directory")


	(tgridgf===nothing) && (tgridgf = Grid.M1D(0.0, (pa.ntgf-1)*1.0, pa.ntgf))

	# save original gf
	file=joinpath(folder, "gfobs.csv")
	CSV.write(file,DataFrame(hcat(tgridgf.x, pa.gfobs)))
	# save for imagesc
	file=joinpath(folder, "imgfobs.csv")
	CSV.write(file,DataFrame(hcat(repeat(tgridgf.x,outer=pa.nra),
				      repeat(1:pa.nra,inner=pa.ntgf),vec(pa.gfobs))),)
	# file=joinpath(folder, "gfcal.csv")
	# CSV.write(file,DataFrame(hcat(tgridgf.x, pa.calsave.gf)))

	# compute cross-correlations of gf
	xtgridgf=Grid.M1D_xcorr(tgridgf); # grid
	if(pa.mode==1)
		xgfobs=hcat(Conv.xcorr(pa.gfobs)...)
	elseif(pa.mode==2)
		xgfobs=pa.obs.gf
	end
	file=joinpath(folder, "xgfobs.csv")
	CSV.write(file,DataFrame(hcat(xtgridgf.x, xgfobs)))

	file=joinpath(folder, "imxgfobs.csv")
	CSV.write(file,DataFrame(hcat(repeat(xtgridgf.x,outer=pa.nra),
				      repeat(1:pa.nra,inner=xtgridgf.nx),vec(xgfobs[:,1:pa.nra]))),)
	if(pa.mode==1)
		xgfcal=hcat(Conv.xcorr(pa.cal.gf)...)
	elseif(pa.mode==2)
		xgfcal=pa.calsave.gf
	end
	file=joinpath(folder, "xgfcal.csv")
	CSV.write(file,DataFrame(hcat(xtgridgf.x,xgfcal)))
	file=joinpath(folder, "imxgfcal.csv")
	CSV.write(file,DataFrame(hcat(repeat(xtgridgf.x,outer=pa.nra),
				      repeat(1:pa.nra,inner=xtgridgf.nx),vec(xgfcal[:,1:pa.nra]))),)

	# compute cross-correlations without blind Decon
	xgf_nodecon=hcat(Conv.xcorr(pa.dobs, lags=[pa.ntgf-1, pa.ntgf-1])...)
	file=joinpath(folder, "xgf_nodecon.csv")
	CSV.write(file,DataFrame(hcat(xtgridgf.x,xgf_nodecon)))

	file=joinpath(folder, "imxgf_nodecon.csv")
	CSV.write(file,DataFrame(hcat(repeat(xtgridgf.x,outer=pa.nra),
				      repeat(1:pa.nra,inner=xtgridgf.nx),vec(xgf_nodecon[:,1:pa.nra]))),)

	# compute cross-correlation of source wavelet
	xwavobs=hcat(Conv.xcorr(pa.wavobs, lags=[pa.ntgf-1, pa.ntgf-1])...)
	file=joinpath(folder, "xwavobs.csv")
	CSV.write(file,DataFrame(hcat(xtgridgf.x,xwavobs)))

	# cross-plot of gf
	file=joinpath(folder, "gfcross.csv")
	CSV.write(file,DataFrame( hcat(vec(xgfobs), vec(xgfcal))))
	file=joinpath(folder, "gfcross_nodecon.csv")
	CSV.write(file,DataFrame( hcat(vec(xgfobs), vec(xgf_nodecon))))


	# compute autocorrelations of source
	if(pa.mode==1)
		xwavobs=autocor(pa.obs.wav[:,1], 1:pa.nt-1, demean=true)
		xwavcal=autocor(pa.calsave.wav[:,1], 1:pa.nt-1, demean=true)
	elseif(pa.mode==2)
		xwavobs=pa.obs.wav[:,1]
		xwavcal=pa.calsave.wav[:,1]
	end
	wavmat=hcat(vec(xwavobs), vec(xwavcal))
	scale!(wavmat, inv(maximum(abs, wavmat)))
	
	# resample wav, final sample should not exceed 1000
	nwav=length(xwavobs)
	fact=(nwav>1000) ? round(Int,nwav/1000) : 1
	xwavobs=xwavobs[1:fact:nwav] # resample
	xwavcal=xwavcal[1:fact:nwav] # resample

	# x plots of wav
	file=joinpath(folder, "wavcross.csv")
	CSV.write(file,DataFrame(wavmat))

	# x plots of data
	dcal=pa.calsave.d
	dobs=pa.obs.d
	nd=length(dcal);
	fact=(nd>1000) ? round(Int,nd/1000) : 1
	dcal=dcal[1:fact:nd]
	dobs=dobs[1:fact:nd]
	datmat=hcat(vec(dcal), vec(dobs))
	scale!(datmat, inv(maximum(abs, datmat)))

	file=joinpath(folder, "datcross.csv")
	CSV.write(file,DataFrame( datmat))

	# finally save err
	err!(pa, cal=pa.calsave) # compute err in calsave 
	file=joinpath(folder, "err.csv")
	CSV.write(file, pa.err)
end




function phase_retrievel(Ay,
		 store_trace::Bool=false, 
		 extended_trace::Bool=false, 
	     f_tol::Float64=1e-12, g_tol::Float64=1e-30, x_tol::Float64=1e-30)



	nt1=size(Ay,1)
	nr1=size(Ay,2)

	nr=Int(sqrt(nr1))
	nt=Int((nt1+1)/2)

	Ayy=[Ay[:,1+(ir-1)*nr:ir*nr]for ir in 1:nr]
	pacse=Misfits.Param_CSE(nt,nr,Ay=Ayy)
	f=function f(x)
		x=reshape(x,nt,nr)
		J=Misfits.error_corr_squared_euclidean!(nothing,x,pacse)
		return J
	end
	g! =function g!(storage, x) 
		x=reshape(x,nt,nr)
		gg=zeros(nt,nr)
		Misfits.error_corr_squared_euclidean!(gg, x, pacse)
		copy!(storage,gg)
	end

	# initial w to x
	x=randn(nt*nr)

	"""
	Unbounded LBFGS inversion, only for testing
	"""
	res = optimize(f, g!, x, 
		       LBFGS(),
		       Optim.Options(g_tol = g_tol, f_tol=f_tol, x_tol=x_tol,
		       iterations = 2000, store_trace = store_trace,
		       extended_trace=extended_trace, show_trace = true))
	println(res)

	x=reshape(Optim.minimizer(res), nt, nr)
	
	return x

end







"""
Deterministic Decon, where wavelet is known
"""
mutable struct ParamD{T<:Real,N}
	np2::Int
	ntd::Int
	ntwav::Int
	d::Array{T,N}
	gf::Array{T,N}
	wav::Vector{T}
	dpad::Array{T,N}
	gfpad::Array{T,N}
	wavpad::Vector{T}
	dfreq::Array{Complex{T},N}
	gffreq::Array{Complex{T},N}
	wavfreq::Array{Complex{T},1}
	fftplan::Base.DFT.FFTW.rFFTWPlan
	ifftplan::Base.DFT.ScaledPlan
	fftplanwav::Base.DFT.FFTW.rFFTWPlan
	ϵ::T
end

function ParamD(;ntd=1, ntwav=1, dims=(), np2=nextfastfft(maximum([2*ntwav, 2*ntd])), # fft dimension for plan
			d=zeros(ntd, dims...), wav=zeros(ntwav), gf=zeros(d), ϵ=1e-2)
	T=eltype(d)

	dims=size(d)[2:end]
	nrfft=div(np2,2)+1
	fftplan=plan_rfft(zeros(T, np2,dims...),[1])
	ifftplan=plan_irfft(complex.(zeros(T, nrfft,dims...)),np2,[1])
	fftplanwav=plan_rfft(zeros(T, np2,),[1])
	
	dfreq=complex.(zeros(T,nrfft,dims...))
	gffreq=complex.(zeros(T,nrfft,dims...))
	wavfreq=complex.(zeros(T,nrfft))

	# preallocate padded arrays
	dpad=(zeros(T,np2,dims...))
	gfpad=(zeros(T,np2,dims...))
	wavpad=(zeros(T,np2,))

	wavv=normalize(wav) # make a copy, don't edit wav

	return ParamD(np2,ntd,ntwav,d,gf,wavv,dpad,gfpad,wavpad,dfreq,gffreq,wavfreq,
		fftplan, ifftplan, fftplanwav, ϵ)

end

"""
Convolution that allocates `Param` internally.
"""
function mod!{T,N}(
	   d::AbstractArray{T,N}, 
	   wav::AbstractVector{T}, attrib::Symbol)
	ntd=size(d,1)
	ntgf=size(gf,1)
	ntwav=size(wav,1)

	# allocation of freq matrices
	pa=ParamD(ntd=ntd, ntwav=ntwav, wav=wav, d=d)

	# using pa, return d, gf, wav according to attrib
	mod!(pa)
end

"""
Convolution modelling with no allocations at all.
By default, the fields `gf`, `d` and `wav` in pa are modified accordingly.
Otherwise use keyword arguments to input them.
"""
function mod!(pa::ParamD; 
	      gf=pa.gf, d=pa.d, wav=pa.wav # external arrays to be modified
	     )
	T=eltype(pa.d)
	ntd=size(pa.d,1)
	ntwav=size(pa.wav,1)
	
	# initialize freq vectors
	pa.dfreq[:] = complex(T(0))
	pa.gffreq[:] = complex(T(0))
	pa.wavfreq[:] = complex(T(0))

	pa.gfpad[:]=T(0)
	pa.dpad[:]=T(0)
	pa.wavpad[:]=T(0)

	# necessary zero padding
	Conv.pad_truncate!(gf, pa.gfpad, ntd-1, 0, pa.np2, 1)
	Conv.pad_truncate!(d, pa.dpad, ntd-1, 0, pa.np2, 1)
	Conv.pad_truncate!(wav, pa.wavpad, ntwav-1, 0, pa.np2, 1)

	A_mul_B!(pa.wavfreq, pa.fftplanwav, pa.wavpad)
	A_mul_B!(pa.dfreq, pa.fftplan, pa.dpad)
	for i in eachindex(pa.gffreq)
		ii=ind2sub(pa.gffreq,i)[1]
		pa.gffreq[i]=pa.dfreq[i]*inv(pa.wavfreq[ii]*conj(pa.wavfreq[ii])+pa.ϵ)
	end
	A_mul_B!(pa.gfpad, pa.ifftplan, pa.gffreq)

	Conv.pad_truncate!(gf, pa.gfpad, ntd-1, 0, pa.np2, -1)
	

	return gf

end



end # module
