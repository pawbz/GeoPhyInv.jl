using Distributed
#addprocs(5)

using LinearAlgebra
using GeoPhyInv
using Test
using Misfits
using BenchmarkTools
using Profile
using LinearMaps


function initialize(fp)
	println(" *********************************")
	println(" *********** $(fp) ***************")
	println(" *********************************")
	include(fp)
end

folder="media"
for t in ["convert", "param_adj", "doc"]
	fp = joinpath(folder, string(t, ".jl"))
	initialize(fp)
end

folder="fdtd"
for t in ["accuracy", "backprop", "rho_projection", "create_snaps"]
	fp = joinpath(folder, string(t, ".jl"))
	initialize(fp)
end

folder="fwi"
for t in ["gradient_accuracy", "born_map"]
	fp = joinpath(folder, string(t, ".jl"))
	initialize(fp)
end


folder="Interpolation"
for t in ["interp_tests"]
	fp = joinpath(folder, string(t, ".jl"))
	initialize(fp)
end

folder="Poisson"
for t in ["testscript_RandomEigenfns", "adj_state_expt", "adj_state", "testdAdx", "forw", "test_born"]
	fp = joinpath(folder, string(t, ".jl"))
	initialize(fp)
end


