module DSP

import SIT.Grid
using Distributions

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

maximum(ax) == 0.0 ? error("x is zero") : ax /= maximum(ax);

if(attrib == :max)
	return maximum(fnpow2grid.x[ax .>= threshold])
elseif(attrib == :min)
	return minimum(fnpow2grid.x[ax .>= threshold])
elseif(attrib == :peak)
	return fnpow2grid.x[indmax(ax)]
end

end





end # module
