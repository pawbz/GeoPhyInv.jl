
module Utils

using DSP
using FFTW
using AxisArrays

include("Wavelets.jl")

for file in ["freq", "taper"]
	fn=joinpath(@__DIR__, string(file,".jl"))
	include(fn)
end



end
