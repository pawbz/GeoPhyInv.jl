
module Utils

using DSP
using AxisArrays

include("Wavelets.jl")

for file in ["freq"]
	fn=joinpath(@__DIR__, string(file,".jl"))
	#include(fn)
end



end
