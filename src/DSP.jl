__precompile__()

module DSP

import JuMIT.Grid
using Distributions
using DistributedArrays
using DSP # from julia

"""
Random source signal models for Param.
Generate a band-limited 
random signal of a particular maximum time, `tmax`.
Sinc function in the frequency domain corresponds to 
a box car function is applied in the time domain.

This method has some extra allocations.
"""
function freq_rand!(x;
	tgrid=nothing,
	fmin=nothing, # minimum fraction of Nyquist sampling
	fmax=nothing, # maximum fraction of Nyquist
	rng=Normal(),
	verbose=false,
	)

	nt=size(x,1)
	(tgrid===nothing) && (tgrid=Grid.M1D(0.0, (nt-1)*1.0, 1.0))  # create dummy tgrid
	fgrid=Grid.M1D_rfft(tgrid) # corresponding fgrid
	(nt≠tgrid.nx) && error("dimension")

	tmax=abs(tgrid.x[end]-tgrid.x[1])

	xfreq=complex.(zeros(fgrid.nx)) # allocation in freq domain
	xx=zeros(tgrid.nx) # allocation in time domain

	fs=inv(tgrid.δx)
	freqmin=(fmin===nothing) ? 0.0 : fmin*fs*0.5
	freqmax=(fmax===nothing) ? fs*0.5 : fmax*fs*0.5
	(freqmax ≤ freqmin) && error("fmin and fmax")

	Δf = inv(tmax) # resolution in the frequency domain
	(Δf ≤ fgrid.δx) ? error("desired resolution smaller than grid sampling") :
	(Δf ≥ (freqmax-freqmin)) ? error("need to increase tmax") :
	fvec = [f for f in freqmin:Δf:freqmax]
	ifvec = fill(0, size(fvec))

	verbose && println("number of frequencies added to signal:\t", size(fvec,1))
	verbose && println("interval between random variable:\t", Δf)
	verbose && println("minimum frequency added to signal\t",minimum(fvec))
	verbose && println("maximum frequency added to signal\t",maximum(fvec))
	for iff in eachindex(fvec)
		ifvec[iff] =  indmin((fgrid.x - fvec[iff]).^2.0)
	end

	for iff in eachindex(fvec)
		# translated sinc function
		fshift=fgrid.x[ifvec[iff]]
		X=rand(rng) * exp(im*rand(Uniform(-pi, pi)))
		for ifff in 1:fgrid.nx
			xfreq[ifff] += (X * sinc(tmax*(abs(fgrid.x[ifff] - fshift))))
		end
	end
	xx=irfft(xfreq, nt)

	copy!(x, xx)
	return x

end


"""
Tapering is necessary to be able to input random signal into finite-difference code
Filtering tapering are applied only if the length of the time series is greater than 10
"""
function get_tapered_random_tmax_signal(tgrid::Grid.M1D; 
					fmin=nothing,
					fmax=nothing,
					tmaxfrac::Float64=1.0,
					dist=Uniform(-2.0, 2.0),
					sparsep=1.0,
					taperperc=20.
					)
	filt_flag=(tgrid.nx > 5) && (!(fmin === nothing)) && (!(fmax===nothing))

	fs = 1/ tgrid.δx;
	if(filt_flag)
		designmethod = Butterworth(6);
		filtsource = Bandpass(fmin, fmax; fs=fs);
	end

	itind = indmin(abs.(tgrid.x-abs(tmaxfrac)*tgrid.x[end]))
	if(tmaxfrac>0.0)
		its=1:itind
	elseif(tmaxfrac<0.0)
		its=itind:tgrid.nx
	end
	# 20% taper window
	twin = taper(ones(length(its)),taperperc) 
	X = zeros(length(its))
	wavsrc = zeros(tgrid.nx) 
	if(filt_flag) 
		X[:] = rand(dist, length(its)) .* twin
	else
		X[:] = rand(dist, length(its))
	end
	if(sparsep ≠ 1.0)
		Xs=sprandn(length(X), sparsep)
		X[findn(Xs.==0.0)]=0.0
	end
	# band limit
	(filt_flag) && (filt!(X, digitalfilter(filtsource, designmethod), X))
	
	(length(X) ≠ 1) && normalize!(X)
	wavsrc[its] = X
	return wavsrc
end


"""
* `x` : first dimension is time
* `p` : sparsity
* `rng` : Random Number Generator
* `btperc` : Taper Perc
* `etperc` : Taper Perc
* `bfrac` : zeros at the beginning of each coloumn
* `efrac` : zeros fraction at the end of each coloumn
"""
function tapered_rand!(x;p=1.0,rng=Normal(),btperc=0., etperc=0.,bfrac=0.0, efrac=0.0)
	x[:]=0.0
	nt=size(x,1)

	a=1+round(Int,bfrac*nt)
	b=nt-round(Int,efrac*nt)

	dd=size(x)[2:end]
	for i in CartesianRange(dd)
		xx=view(x,a:b,i)
		rand!(rng, xx)
		taper!(xx, bperc=btperc, eperc=etperc)
		nxx=size(xx,1)
		if(p ≠ 1.0)
			xs=sprandn(nxx, p)
			xx[findn(xs.==0.0)]=0.0
		end
	end
	scale!(x, inv(vecnorm(x)))
	return x
end

"""
randomly shift and add x to itself
using circshift; not memory efficient
"""
function cshift_add_rand!(x::AbstractVector; ntimes=1, tminfrac=0.0, tmaxfrac=1.0)
	nt=size(x,1)
	a=1+round(Int,tminfrac*nt)
	b=round(Int,tmaxfrac*nt)

	for it in 1:ntimes
		# a random shift
		its=rand(DiscreteUniform(a,b))
		xx = circshift(x,(its,)) # pre-allocate for performance?
		for i in eachindex(x)
			x[i] += xx[i]
		end
	end
end


"""
Construct Toy Green's functions
Decaying peaks, control number of events, and their positions, depending on bfrac and efrac.
* `afrac` : control amplitude of peaks
"""
function toy_green!(x;nevents=1,afrac=[1./2.^(ie-1) for ie in 1:nevents],bfrac=0.0,efrac=0.0)
	nt=size(x,1)
	nr=size(x,2)
	x[:]=0.0
	a=1+round(Int,bfrac*nt)
	b=nt-round(Int,efrac*nt)
	(b-a<nevents+2) && error("not enough samples")

	itvec=zeros(Int64,nevents+1)
	for ir in 1:nr
		itvec[:]=0
		for ie in 1:nevents
			if(ie==1)
				# first event, direct arrival
				itvec[ie]=rand(DiscreteUniform(a,b-nevents-1))
			else
				# all other events after first event
				itvec[ie]=rand(DiscreteUniform(itvec[1]+1,b-1))
			end
			sgn=1.
			# randomly choose sign for all other events
			if(ie≠1)
				sgn=(rand(Bernoulli())==1) ? 1. : -1.
			end
			x[itvec[ie],ir]=sgn*afrac[ie]
		end
	end
	return x
end


function findfreq{ND}(
		  x::Array{Float64, ND},
		  tgrid::Grid.M1D;
		  attrib::Symbol=:peak,
		  threshold::Float64=-50.
		  )

nfft = nextpow2(tgrid.nx);
# npow2 grid for time
tnpow2grid = Grid.M1D_fft(nfft, tgrid.δx);
# corresponding npow2 frequency grid 
fnpow2grid = Grid.M1D_fft(tnpow2grid);

cx = fill(complex(0.0,0.0),nfft);
cx[1:tgrid.nx] = complex.(x,0.0);

cx = fft(cx);
ax = (abs.(cx).^2); # power spectrum in dB
ax[fnpow2grid.x .< 0.] = 0. # remove negative frequencies

if(maximum(ax) == 0.0)
	warn("x is zero"); return 0.0
else 
	ax /= maximum(ax);
	ax = 10. .* log10.(ax)
end

if(attrib == :max)
	return maximum(fnpow2grid.x[ax .>= threshold])
elseif(attrib == :min)
	return minimum(fnpow2grid.x[ax .>= threshold])
elseif(attrib == :peak)
	return fnpow2grid.x[indmax(ax)]
end

end


"""
Cosine taper a N-dimensional array along its first dimension.

# Arguments
* `x::Array{Float64,N}` : 
* `perc::Float64` : taper percentage
"""
function taper(x::AbstractArray,perc::Float64)
	xout=copy(x);
	taper!(xout,perc)
	return xout
end

function taper!{N}(x::AbstractArray{Float64,N}, perc::Float64=0.0; bperc=perc,eperc=perc)

	nt=size(x,1)
	nttb=min(round(Int,nt*bperc/100.), nt)
	ntte=min(round(Int,nt*eperc/100.), nt)
	dd=size(x)[2:N]
	for i in CartesianRange(dd)
		kb=inv(2.*round(Int,nt*bperc/100.)-1)*pi
		for it in 1:nttb
			x[it,i] *= sin((it-1)*kb)
		end
		ke=inv(2.*round(Int,nt*eperc/100.)-1)*pi
		for it in nt-ntte+1:nt 
			x[it,i] *= sin((-it+nt)*ke)
		end
	end
end



"""
A bandpass butterworth filter using `fmin` and `fmax` in the frequency domain.
Return either a zero-phase or minimum-phase filter using `attrib`.
More info on minimum-phase filtering: http://www.katjaas.nl/minimumphase/minimumphase.html.
! then it is converted to minimum phase filter in the cepstral domain
! more info: http://www.katjaas.nl/minimumphase/minimumphase.html
! positive quefrencies correspond to minimum phase component of the signal
! negetive quefrencies correspond to maximum phase component of the signal
! signal = minimum phase [convolved with] maximum phase


# Arguments

* ``

# Keyword Arguments

* `order::Int64` : order of the butterworth filter
* `attrib::Symbol`
  * `=:zp` means zero phase 
  * `=:mp` means zero phase 
* `fmin::Float64`
  * `=0.` means a lowpass filter at fmax
* `fmax::Float64`
  * `=0.` means a highpass filter at fmin

"""
function butterworth_filter(fvec; order::Int64=2, attrib::Symbol=:zp, fmin=0.0, fmax=0.0)
# Author : Pawan Bharadwaj
#          p.b.pisupati@tudelft.nl
# August 2017, imported to Julia from FORTRAN90

(fmin < 0.0) && error("fmin cannot be .lt. zero")
(fmax < 0.0) && error("fmax cannot be .lt. zero")
(fmax ≠ 0.0) && (fmax <= fmin) && error("fmax .le. fmin") 

if(fmin == 0.0) 
        # lowpass filter
	F  = complex.((1.0 + (fvec./fmax).^(2.0*order)).^(-1.0), 0.0)
elseif(fmax == 0.0) 
        # highpass filter
	F  = complex.((1.0 + (fmin./fvec).^(2.0*order)).^(-1.0), 0.0)
elseif((fmin ≠ 0.0) & (fmax ≠ 0.0))
        # bandpass filter
        x1 = sqrt(fmin*fmax) / (fmax-fmin) * (fvec./sqrt(fmin*fmax) + sqrt(fmax*fmin)./fvec);
        F  = complex.((1.0 + (x1).^(2.0*order)).^(-1.0), 0.0)
end
# DC component is always forced zero
F[1] = complex(0.0, 0.0);

# conversion to minimum phase
if(attrib == :mp)
# to prevent log(0)
	damp = 1e-20 * maximum(abs.(F))
	# logarithm 
	X = log.(complex.(abs.(F) + damp, 0.0)) 
	# to cepstral domain - IFFT
	ifft!(X)
	# only real part output
	X = complex.(real.(X), 0.0);
	# scaling
	# X = X / complex(real(npow2), 0.0)

	# positive cepstrum x 2
	X[2 : npow2/2 + 1] *= complex(2.0, 0.0)
	# remove negative quefrencies
	X[npow2/2 + 2 : npow2] = complex(0.0, 0.0) 

	# FFT
	fft!(X)

	# exponential
	F = exp.(X)

	F /= complex(maximum(abs.(F)), 0.0)
end

return F

end # butterworth


end # module
