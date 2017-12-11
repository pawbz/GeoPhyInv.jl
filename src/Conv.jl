__precompile__()

"""
This module is used in DeConv, ConvMix and Coupling
"""
module Conv
using DSP

"""
Model : d = convolution(gf, wav)
d, gf and wav can have arbitary +ve and -ve lags
"""
type Param{T<:Real,N}
	ntwav::Int
	ntgf::Int
	ntd::Int
	"length after zero padding"
	np2::Int
	d::Array{T,N}
	gf::Array{T,N}
	wav::Array{T,N}
	dpad::Array{T,N}
	gfpad::Array{T,N}
	wavpad::Array{T,N}
	dfreq::Array{Complex{T},N}
	gffreq::Array{Complex{T},N}
	wavfreq::Array{Complex{T},N}
	fftplan::Base.DFT.FFTW.rFFTWPlan
	ifftplan::Base.DFT.ScaledPlan
	"+ve and -ve lags of gf"
	gflags::Vector{Int}
	"+ve and -ve lags of wav"
	wavlags::Vector{Int}
	"+ve and -ve lags of d"
	dlags::Vector{Int}
end

function Param(;ntwav=1, ntgf=1, ntd=1, dims=(), 
	       d=zeros(ntd, dims...), wav=zeros(ntwav, dims...), gf=zeros(ntgf, dims...),
	      wavlags=nothing, dlags=nothing, gflags=nothing)
	(eltype(d) ≠ eltype(wav) ≠ eltype(gf)) && error("type error")
	T=eltype(d)

	# fft dimension for plan
	np2=nextfastfft(maximum([2*ntwav, 2*ntgf, 2*ntd]))
	dims=size(d)[2:end]

	# default lags
	if(wavlags===nothing)
		# equal +ve and -ve lags for wav
		nwplags=div(ntwav-1,2)
		nwnlags=ntwav-1-nwplags
		wavlags=[nwplags, nwnlags]
	end
	if(dlags===nothing)
		# no negative lags
		nsplags=ntd-1
		nsnlags=ntd-1-nsplags
		dlags=[nsplags, nsnlags]
	end
	if(gflags===nothing)
		# no negative lags
		nrplags=ntgf-1
		nrnlags=ntgf-1-nrplags
		gflags=[nrplags, nrnlags]
	end

	(size(d)[2:end] ≠ size(wav)[2:end] ≠ size(gf)[2:end]) && error("dimension error")

	# check if nlags are consistent with the first dimension of inputs
	(sum(gflags)+1 ≠ size(gf,1)) && error("gflags")
	(sum(dlags)+1 ≠ size(d,1)) && error("dlags")
	(sum(wavlags)+1 ≠ size(wav,1)) && error("wavlags")

	#FFTW.set_num_threads(Sys.CPU_CORES)
	nrfft=div(np2,2)+1
	fftplan=plan_rfft(zeros(T, np2,dims...),[1])
	ifftplan=plan_irfft(complex.(zeros(T, nrfft,dims...)),np2,[1])

	# preallocate freq domain vectors after rfft
	dfreq=complex.(zeros(T,nrfft,dims...))
	gffreq=complex.(zeros(T,nrfft,dims...))
	wavfreq=complex.(zeros(T,nrfft,dims...))

	# preallocate padded arrays
	dpad=(zeros(T,np2,dims...))
	gfpad=(zeros(T,np2,dims...))
	wavpad=(zeros(T,np2,dims...))

	pa=Param(ntwav, ntgf, ntd, np2, d, gf, wav, 
	  dpad, gfpad, wavpad,
	  dfreq, gffreq, wavfreq, 
	  fftplan, ifftplan, gflags, wavlags, dlags)
end


"""
Convolution that allocates pa.
"""
function mod!{T,N}(
	   d::AbstractArray{T,N}, 
	   gf::AbstractArray{T,N},
	   wav::AbstractArray{T,N}, attrib::Symbol)
	ntd=size(d,1)
	ntgf=size(gf,1)
	ntwav=size(wav,1)

	# allocation of freq matrices
	pa=Param(ntgf=ntgf, ntd=ntd, ntwav=ntwav, gf=gf, wav=wav, d=d)

	# using pa
	mod!(pa, attrib)
end

"""
Convolution modelling with no allocations at all.
"""
function mod!(pa::Param, attrib::Symbol)
	T=eltype(pa.d)
	
	# initialize freq vectors
	pa.dfreq[:] = complex(T(0))
	pa.gffreq[:] = complex(T(0))
	pa.wavfreq[:] = complex(T(0))

	# necessary zero padding
	pad_truncate!(pa.gf, pa.gfpad, pa.gflags[1], pa.gflags[2], pa.np2, 1)
	pad_truncate!(pa.d, pa.dpad, pa.dlags[1], pa.dlags[2], pa.np2, 1)
	pad_truncate!(pa.wav, pa.wavpad, pa.wavlags[1], pa.wavlags[2], pa.np2, 1)

	if(attrib == :d)
		A_mul_B!(pa.wavfreq, pa.fftplan, pa.wavpad)
		A_mul_B!(pa.gffreq, pa.fftplan, pa.gfpad)
		@. pa.dfreq = pa.gffreq * pa.wavfreq
		A_mul_B!(pa.dpad, pa.ifftplan, pa.dfreq)
		pad_truncate!(pa.d, pa.dpad, pa.dlags[1], pa.dlags[2], pa.np2, -1)
		return pa.d
	elseif(attrib == :gf)
		A_mul_B!(pa.wavfreq, pa.fftplan, pa.wavpad)
		A_mul_B!(pa.dfreq, pa.fftplan, pa.dpad)
		conj!(pa.wavfreq)
		@. pa.gffreq = pa.dfreq * pa.wavfreq
		A_mul_B!(pa.gfpad, pa.ifftplan, pa.gffreq)
		pad_truncate!(pa.gf, pa.gfpad, pa.gflags[1], pa.gflags[2], pa.np2, -1)
		
		return pa.gf
	elseif(attrib == :wav)
		A_mul_B!(pa.gffreq, pa.fftplan, pa.gfpad)
		A_mul_B!(pa.dfreq, pa.fftplan, pa.dpad)
		conj!(pa.gffreq)
		@. pa.wavfreq = pa.dfreq * pa.gffreq
		A_mul_B!(pa.wavpad, pa.ifftplan, pa.wavfreq)
		pad_truncate!(pa.wav, pa.wavpad, pa.wavlags[1], pa.wavlags[2], pa.np2, -1)
		return pa.wav
	end
end


"""
Method to perform zero padding and truncation.

# Arguments

* `x` : real signal with dimension nplag + nnlag + 1
	first has decreasing negative lags, 
	then has zero lag at nnlags + 1,
	then has increasing positive nplags lags,
	signal contains only positive lags and zero lag if nnlag=0 and vice versa
* `xpow2` : npow2 real vector with dimension npow2
* `nplags` : number of positive lags
* `nnlags` : number of negative lags
* `npow2` : number of samples in xpow2
* `flag` : = 1 means xpow2 is returned using x
	   = -1 means x is returned using xpow2
"""
function pad_truncate!{T}(
				  x::AbstractArray{T}, 
				  xpow2::AbstractArray{T}, 
				  nplags::Integer, 
				  nnlags::Integer, 
				  npow2::Integer, 
				  flag::Integer
				  )
	(size(x,1) ≠ nplags + nnlags + 1) && error("size x")
	(size(xpow2,1) ≠ npow2) && error("size xpow2")

	for id in 1:size(x,2)
		if(flag == 1)
			xpow2[1,id] = (x[nnlags+1,id]) # zero lag
			# +ve lags
			if (nplags > 0) 
				for i=1:nplags
					xpow2[i+1,id]= (x[nnlags+1+i,id])
				end
			end
			# -ve lags
			if(nnlags != 0) 
				for i=1:nnlags
					xpow2[npow2-i+1,id] =(x[nnlags+1-i,id])
				end
			end
		elseif(flag == -1)
			x[nnlags+1,id] = (xpow2[1,id]); # zero lag
			if(nplags != 0) 
				for i=1:nplags
					x[nnlags+1+i,id] = (xpow2[1+i,id]);
				end
			end
			if(nnlags != 0)
				for i=1:nnlags
					x[nnlags+1-i,id] = (xpow2[npow2-i+1,id])
				end
			end
		else
			error("invalid flag")
		end
	end
	return nothing
end




end
