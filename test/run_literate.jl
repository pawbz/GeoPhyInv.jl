
# run this script to update the pages for Documenter.jl

# to generate doc pages
using Literate


function run_literate(names, folder)
	for t in names
		fp = joinpath(folder, string(t, ".jl"))
		output_folder=joinpath(@__DIR__, "..", "docs", "generated", folder) 
		output_test_folder=joinpath(@__DIR__, "generated", folder) 
		println(output_folder)
		Literate.markdown(fp, output_folder, documenter=true,credit=false)
		Literate.script(fp, output_test_folder, documenter=true,credit=false)
#		Literate.notebook(fp, output_folder, documenter=true,credit=false)
	end
end

run_literate(["doc"], "media")

run_literate(["gallery", "doc"], "ageom")

run_literate(["doc"], "srcwav")

run_literate(["doc"], "data")

run_literate(["forw", "test_born"], "Poisson")

run_literate(["gradient_accuracy","born_map"], "fwi")

run_literate(["create_snaps","reuse_expt"], "fdtd")
