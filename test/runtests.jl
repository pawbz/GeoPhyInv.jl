using Distributed
addprocs(5)

using LinearAlgebra
using GeoPhyInv
using Test
using Misfits
using BenchmarkTools
using Profile


function initialize(fp)
	println(" *********************************")
	println(" *********** $(fp) ***************")
	println(" *********************************")
	include(fp)
end


folder="Poisson"
for t in ["testscript_RandomEigenfns", "adj_state_expt", "adj_state", "testdAdx", "forw", "test_born"]
	fp = joinpath(folder, string(t, ".jl"))
	initialize(fp)
end


folder="Interpolation"
for t in ["interp_tests"]
	fp = joinpath(folder, string(t, ".jl"))
	initialize(fp)
end



folder="Models"
for t in ["Models", "param_adj"]
	fp = joinpath(folder, string(t, ".jl"))
	initialize(fp)
end



folder="Fdtd"
for t in ["accuracy", "backprop", "rho_projection", "create_snaps"]
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
