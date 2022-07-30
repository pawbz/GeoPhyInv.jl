# run this script to update the pages for Documenter.jl, either using Literate.jl or PlutoStaticHTML.jl

# to generate doc pages
using Pkg; Pkg.add.(["Plots", "DSP", "PlutoStaticHTML"])
using Literate
using PlutoStaticHTML
using Weave
using Dates


function update_date(content)
    content = replace(content, "DATEOFTODAY" => Date(now()))
    return content
end


output_folder = joinpath(@__DIR__, "..", "docs", "src", "generated")
output_test_folder = joinpath(@__DIR__, "generated")
foreach(x -> rm(x, recursive = true), readdir(output_test_folder, join = true))
foreach(x -> rm(x, recursive = true), readdir(output_folder, join = true))

# clearing files in the output folders 

function run_literate(literate_names, folder, pluto_notebooks=[])
    for t in literate_names
        fp = joinpath(@__DIR__, folder, string(t, ".jl"))
        of = joinpath(output_folder, folder)
        otf = joinpath(output_test_folder, folder)
        # fpm = joinpath(of, string(t, ".md"))
        mkpath(of)
        mkpath(otf)
        println("writing files to: \n", of, "\n", otf)
        Literate.markdown(
            fp,
            of,
            documenter = true,
            credit = false,
            preprocess = update_date,
        )
        # weave(fpm; doctype = "md2html", out_path = of)
        Literate.script(
            fp,
            otf,
            documenter = true,
            credit = false,
            preprocess = update_date,
        )
        #		Literate.notebook(fp, output_folder, documenter=true,credit=false)
    end
    for t in pluto_notebooks
        of = joinpath(output_folder, folder)
        tutorial_src_folder=joinpath(@__DIR__, folder)
        rendered = build_notebooks( #though not really parallel here
        BuildOptions(
            tutorial_src_folder;
            output_format=documenter_output,
            write_files=false,
            use_distributed=false,
        ),
        ["$t.jl"],
    )
     write(
        joinpath(of, string(t,".md")),
        """
        ```@meta
        EditURL = "$(tutorial_src_folder)$(t).jl"
        ```
        $(rendered[1])
        """,
    )

    end
end

run_literate(["doc"], "media")

run_literate(["doc"], "ageom")

run_literate(["doc"], "srcwav")

run_literate(["doc"], "data")

# run_literate(["doc","forw", "test_born"], "Poisson")

#run_literate(["gradient_accuracy","born_map"], "fwi")

# run_literate(["doc", "pizza", "born_tutorial"], "fwi")
run_literate(["doc", "create_snaps",], "fdtd", ["reuse_expt", "grid_visualize"])
