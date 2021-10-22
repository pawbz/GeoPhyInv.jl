using Distributed
#addprocs(5)

using LinearAlgebra
using GeoPhyInv
using Test
using Distances
using BenchmarkTools
using Profile
using LinearMaps


"""
* `files` are in folder
* `gfiles` are files in the generated folder
"""
function run_test(files,folder,gfiles=[])
	for file in files
		fp = joinpath(folder, string(file, ".jl"))
		println(" ***********************************************")
		println(" *********** $(fp) ***************")
		println(" ***********************************************")
		include(fp)
	end
	for file in gfiles
		fpg = joinpath("generated",folder,string(file, ".jl"))
		println(" ***********************************************")
		println(" ******** $(fpg) *********")
		println(" ***********************************************")
		include(fpg)
	end
end


run_test(["convert", "param_adj"],"media", ["doc"])
run_test([],"ageom", ["doc"])
run_test([],"srcwav", ["doc"])
run_test([],"data", ["doc"])
# run_test(["accuracy2D", "rho_projection"],"fdtd", ["doc","create_snaps","reuse_expt"])
#run_test(["accuracy2D","backprop", "rho_projection"],"fdtd", ["doc","create_snaps","reuse_expt"])
#run_test(["gradient_accuracy", "born_map"],"fwi", ["doc", "pizza", "born_tutorial"])
run_test(["interp_tests"],"Interpolation", [])
run_test(["testscript_RandomEigenfns", "adj_state_expt", "adj_state", "testdAdx"], "Poisson", ["doc","forw", "test_born"])



