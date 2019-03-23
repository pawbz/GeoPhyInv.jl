
# run this script to update the pages for Documenter.jl

# to generate doc pages
using Literate


function run_literate(names, folder)
	for t in names
		fp = joinpath(folder, string(t, ".jl"))
		output_folder=joinpath(@__DIR__, "..", "docs", "src", folder) 
		println(output_folder)
		Literate.markdown(fp, output_folder, documenter=true)
	end
end


run_literate(["adj_state_expt"], "Poisson")

run_literate(["gradient_accuracy","born_map"], "FWI")

run_literate(["create_snaps"], "Fdtd")
