__precompile__()

module DSP

import JuMIT.Grid
using Distributions
using DSP # from julia


function get_tapered_random_tmax_signal(tgrid::Grid.M1D, fmin::Float64, fmax::Float64, tmax::Float64)

	fs = 1/ tgrid.δx;
	designmethod = Butterworth(6);
	filtsource = Bandpass(fmin, fmax; fs=fs);

	itmax = indmin(abs.(tgrid.x-tmax))
	# 20% taper window
	twin = taper(ones(itmax),20.) 
	X = zeros(itmax)
	wavsrc = zeros(tgrid.nx) 
	X[:] = rand(Uniform(-1.0, 1.0), itmax) .* twin
	# band limit
	filt!(X, digitalfilter(filtsource, designmethod), X);
	
	#wavsrc[1:itmax] = X.*twin
	wavsrc[1:itmax] = X
	return wavsrc
end



"""
Generate a band-limited 
random signal of a particular maximum time, `tmax`.



such that a box car function is applied in the time domain
"""
function get_random_tmax_signal(;
	fgrid::Grid.M1D=nothing, # frequency domain grid
	fmin::Float64=0.0, # minimum frequency
	fmax::Float64=nothing, # maximum frequency
	tmax::Float64=nothing, # frequency sampling is decided based on the length in time
	dist::Symbol=:gaussian # distribution type
	)

	# initialize outputs
	S = fill(complex(0.0,0.0),fgrid.nx);
	Sreal = fill(0.0,fgrid.nx);
	s = fill(complex(0.0,0.0),fgrid.nx);

	Δf = 1.0 / tmax
	Δf <= fgrid.δx ? error("sampling smaller than grid sampling") :
	Δf >= (fmax-fmin) ? error("need to increase tmax") :
	fvec = [f for f in fmin:Δf:fmax]
	ifvec = fill(0, size(fvec))

	println("number of frequencies added to signal:\t", size(fvec,1))
	println("interval between random variable:\t", Δf)
	println("minimum frequency added to signal\t",minimum(fvec))
	println("maximum frequency added to signal\t",maximum(fvec))
	for iff in eachindex(fvec)
		ifvec[iff] =  indmin((fgrid.x - fvec[iff]).^2.0)
	end

	if(dist == :gaussian)
		X = randn(size(fvec));
	elseif(dist == :uniform)
		X = rand(Uniform(-2.0, 2.0), size(fvec))
	else
		error("invalid dist")
	end

	for iff in eachindex(fvec)
		Sreal += X[iff] .* sinc(tmax.*(abs(fgrid.x - fgrid.x[ifvec[iff]])))
	end

	S = complex.(Sreal, 0.0);

	# remove mean 
	#S[1] = 0.0; 
	s = ifft(S); 

	return S, s

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
function taper{N}(x::Array{Float64,N},perc::Float64)

	n = size(x,1)
	nt = Int64(floor(n*perc/100.0))
	# cosine window
	tap = cosine(2*nt); 
	win = vcat(tap[1:nt], ones(n-2*nt), tap[nt+1:2*nt]);
	
	# apply
	(N==1) ? (return x.*win) : (return x.*repeat(win, 
					      outer=vcat([1],[size(x,s) for s in 2:N])))
end

"""

# Arguments
* `RW` :  the first time series which is causal
"""
function fast_filt!{T<:Real}(
		   s::AbstractArray{T}, 
		   r::AbstractArray{T},
		   w::AbstractArray{T},
		   attrib::Symbol
		  ) 

	isodd(length(w)) ? nl=div(length(w)-1,2) : error("filt length should be even")
	nx = (length(r) == length(s)) ? length(s) : error("x and y length")
	np2 = nextpow2(maximum([2*length(s), 2*length(r), length(w)]));	
		
	rpow2=complex.(zeros(T,np2), zeros(T,np2)); 
	spow2=complex.(zeros(T,np2), zeros(T,np2)); 
	wpow2=complex.(zeros(T,np2), zeros(T,np2));
	nlag_npow2_pad_truncate!(r, rpow2, nx-1, 0, np2, 1)
	nlag_npow2_pad_truncate!(s, spow2, nx-1, 0, np2, 1)

	if(attrib == :s)
		if(w===nothing)
			wpow2 = filt_func()nlag_npow2_pad_truncate!(w, wpow2, nl, nl, np2, 1)
		else
			nlag_npow2_pad_truncate!(w, wpow2, nl, nl, np2, 1)
			fft!(wpow2);
		end
		fft!(rpow2)
		spow2 = rpow2 .* wpow2;
		ifft!(spow2)
		nlag_npow2_pad_truncate!(s, spow2, nx-1, 0, np2, -1)
		return s
	elseif(attrib == :r)
		nlag_npow2_pad_truncate!(flipdim(w,1), wpow2, nl, nl, np2, 1)
		fft!(spow2); fft!(wpow2);
		rpow2 = spow2 .* wpow2
		ifft!(rpow2)
		nlag_npow2_pad_truncate!(r, rpow2, nx-1, 0, np2, -1)
		return r
	elseif(attrib == :w)
		nlag_npow2_pad_truncate!(w, wpow2, nl, nl, np2, 1)
		fft!(rpow2); fft!(spow2);
		wpow2 = spow2 .* conj(rpow2);
		ifft!(wpow2)
		nlag_npow2_pad_truncate!(w, wpow2, nl, nl, np2, -1)
		return w
	end
end

"""

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
function nlag_npow2_pad_truncate!{T<:Number}(
				  x::AbstractVector{T}, 
				  xpow2::AbstractVector{Complex{T}}, 
				  nplags::Integer, 
				  nnlags::Integer, 
				  npow2::Integer, 
				  flag::Integer
				  )
	(length(x) == nplags + nnlags + 1) ? nothing : error("size x")
	(length(xpow2) == npow2) ? nothing : error("size xpow2")

	if(flag == 1)
		xpow2[1] = copy(x[nnlags+1]); # zero lag
		# +ve lags
		if (nplags > 0) 
			xpow2[2:nplags+1] = copy(complex.(x[nnlags+2:nnlags+1+nplags],T(0)))
		end
		# -ve lags
		if(nnlags != 0) 
			for i=1:nnlags
				xpow2[npow2-i+1] = copy(complex.(x[nnlags+1-i],T(0)))
			end
		end
		return xpow2
	elseif(flag == -1)
		x[nnlags+1] =  copy(real(xpow2[1])); # zero lag
		if(nplags != 0) 
			for i=1:nplags
				x[nnlags+1+i] = copy(real(xpow2[1+i]));
			end
		end
		if(nnlags != 0)
			for i=1:nnlags
				x[nnlags+1-i] = copy(real(xpow2[npow2-i+1]))
			end
		end
		return x
	else
		error("invalid flag")
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
