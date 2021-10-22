
module Utils

using DSP
using Distances
using LinearAlgebra
using Test
using FFTW
using AxisArrays


for file in ["freq", "taper", "adjtest"]
	fn=joinpath(@__DIR__, string(file,".jl"))
	include(fn)
end



end
