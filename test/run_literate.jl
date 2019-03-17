
# run this script to update the pages for Documenter.jl

# to generate doc pages
using Literate

folder="Poisson"
for t in ["adj_state_expt"]
	fp = joinpath(folder, string(t, ".jl"))
	output_folder=joinpath(@__DIR__, "..", "docs", "src", folder) 
	println(output_folder)
	Literate.markdown(fp, output_folder, documenter=true)
end
