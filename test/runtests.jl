using Distributed
#addprocs(5)

using LinearAlgebra
using JuMIT
using Test
using Signals
using Misfits
using BenchmarkTools
using Grid
using Profile

function initialize(fp)
	println(" *********************************")
	println(" *********** $(fp) ***************")
	println(" *********************************")
	include(fp)
end

folder="Models"
for t in ["param_adj"]
	fp = joinpath(folder, string(t, ".jl"))
	initialize(fp)
end



folder="Fdtd"
for t in ["rho_projection"]
	fp = joinpath(folder, string(t, ".jl"))
	initialize(fp)
end



folder="FWI"
for t in ["gradient_accuracy", "born_map"]
	fp = joinpath(folder, string(t, ".jl"))
	initialize(fp)
end






#tests = ["Fdtd_accuracy", "Fdtd_backprop", 
#	 # "Poisson_testscript_RandomEigenfns", 
#	"Models", "Data", "FWI_grad_test"]
#
#println("Running tests:")
#
#end
