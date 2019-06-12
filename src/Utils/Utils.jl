
module Utils

using DSP
using Misfits
using LinearAlgebra
using Test
using FFTW
using AxisArrays

include("Wavelets.jl")

for file in ["freq", "taper", "adjtest"]
	fn=joinpath(@__DIR__, string(file,".jl"))
	include(fn)
end



end
