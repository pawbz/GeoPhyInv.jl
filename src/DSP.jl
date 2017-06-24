module DSP

import SIT.Grid
using Distributions
using DSP # from julia

"""
generate random signals in the frequency domain 
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

S = complex(Sreal, 0.0);

# remove mean 
#S[1] = 0.0; 
s = ifft(S); 

return S, s

end



function findfreq{ND}(
		  x::Array{Float64, ND},
		  tgrid::Grid.M1D;
		  attrib::Symbol=:peak,
		  threshold::Float64=1e-6
		  )

nfft = nextpow2(tgrid.nx);
# npow2 grid for time
tnpow2grid = Grid.M1D_npow2(nfft, tgrid.δx);
# corresponding npow2 frequency grid 
fnpow2grid = Grid.M1D_npow2_tf(tnpow2grid);

cx = fill(complex(0.0,0.0),nfft);
cx[1:tgrid.nx] = complex(x,0.0);

cx = fft(cx);
ax = real(cx.*conj(cx));
ax[fnpow2grid.x .< 0.] = 0. # remove negative frequencies

if(maximum(ax) == 0.0)
	warn("x is zero"); return 0.0
else 
	ax /= maximum(ax);
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
		   s::Vector{T}, 
		   r::Vector{T},
		   w::Vector{T},
		   attrib::Symbol
		  ) 

	isodd(length(w)) ? nl=div(length(w)-1,2) : error("filt length should be even")
	nx = (length(r) == length(s)) ? length(s) : error("x and y length")
	np2 = nextpow2(maximum([2*length(s), 2*length(r), length(w)]));	
		
	rpow2=complex(zeros(T,np2), zeros(T,np2)); 
	spow2=complex(zeros(T,np2), zeros(T,np2)); 
	wpow2=complex(zeros(T,np2), zeros(T,np2));
	nlag_npow2_pad_truncate!(r, rpow2, nx-1, 0, np2, 1)
	nlag_npow2_pad_truncate!(s, spow2, nx-1, 0, np2, 1)

	if(attrib == :s)
		nlag_npow2_pad_truncate!(w, wpow2, nl, nl, np2, 1)
		fft!(rpow2); fft!(wpow2);
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
				  x::Vector{T}, 
				  xpow2::Vector{Complex{T}}, 
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
			xpow2[2:nplags+1] = copy(complex(x[nnlags+2:nnlags+1+nplags],T(0)))
		end
		# -ve lags
		if(nnlags != 0) 
			for i=1:nnlags
				xpow2[npow2-i+1] = copy(complex(x[nnlags+1-i],T(0)))
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

end # module
