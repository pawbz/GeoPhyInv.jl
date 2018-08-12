using LinearAlgebra
using JuMIT
using Test
using Signals
using Misfits
using BenchmarkTools
using Grid
using Profile
using Distributed

addprocs(5)

folder="FWI"
for t in ["gradient_accuracy"]
	fp = joinpath(folder, string(t, ".jl"))
	println(" *********** $(fp) ***************")
	include(fp)
end






#tests = ["Fdtd_accuracy", "Fdtd_backprop", 
#	 # "Poisson_testscript_RandomEigenfns", 
#	"Models", "Data", "FWI_grad_test"]
#
#println("Running tests:")
#
#end
