ENV["PLOTS_TEST"] = "true"
ENV["GKSwstype"] = "100"

using Documenter, GeoPhyInv, SparseArrays, LinearAlgebra, Random, LinearMaps
using Plots, Test, Luxor
using ColorSchemes

# standard setting for generating doc pages
# @init_parallel_stencil(2, false, Float32, 2)

gr()
theme(:juno)

#makedocs(modules=[GeoPhyInv],
#	         doctest=true)
#makedocs()

makedocs(
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    modules = [GeoPhyInv],
    sitename = "GeoPhyInv",
    pages = [
        "Home" => "index.md",
        "Medium" => "generated/media/doc.md",
        "AGeom" => "generated/ageom/doc.md",
        "SrcWav" => "generated/srcwav/doc.md",
        "Records" => "generated/data/doc.md",
        "SeisForwExpt" => Any[
            "Description"=>"generated/fdtd/doc.md",
            "Staggered Fields"=>"generated/fdtd/grid_visualize.md",
            "Acoustic"=>"generated/fdtd/marmousi_acoustic.md",
            # "Tutorial I"=>"generated/fdtd/create_snaps.html",
            # "Reuse Expt"=>"generated/fdtd/reuse_expt.md",
            "Reuse Expt"=>"generated/fdtd/reuse_expt.md",
        ],
        # "SeisInvExpt" => Any[
        # 		  "Description" => "generated/fwi/doc.md",
        # 		  "Tutorial" => "generated/fwi/pizza.md",

        # 		  "Adjoint Tests" => "generated/fwi/born_tutorial.md",
        # 		 ],
        # "PoissonExpt" => Any[
        #   "Description" => "generated/Poisson/doc.md",
        #  "Tutorial" => "generated/Poisson/forw.md",
        #  "Adjoint Tests" => "generated/Poisson/test_born.md",
        #  ],
    ],
    #	    "Seismic Born Modeling" => "FWI/born_map.md",
    #	    "Seismic Full Waveform Inversion" => "FWI/gradient_accuracy.md",
    #    modules = []
)

deploydocs(repo = "github.com/pawbz/GeoPhyInv.jl.git")
#=
deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/pawbz/GeoPhyInv.jl.git",
    julia  = "1.0",
    osname = "linux",
    target = "build",
    make   = nothing
)
=#
