using Pkg; Pkg.add("Gadfly"); Pkg.add("Plots");

using Documenter, GeoPhyInv
using Gadfly, Plots
gr()

#makedocs(modules=[GeoPhyInv],
#	         doctest=true)
#makedocs()

makedocs(
	 format = Documenter.HTML(
		  prettyurls = get(ENV, "CI", nothing) == "true"),
	 modules=[GeoPhyInv],
   sitename = "GeoPhyInv",
   pages = ["Home" => "index.md",
	    "Medium" => "generated/media/doc.md",
	    "AGeom" => "generated/ageom/doc.md",
	    "SrcWav" => "generated/srcwav/doc.md",
	    "Data" => "generated/data/doc.md",
	    "SeisForwExpt" => Any[
				  "Introduction" => "generated/fdtd/doc.md",
				  "Basic usage" => "generated/fdtd/reuse_expt.md",
				  "Generate snaps" => "generated/fdtd/create_snaps.md",
				  ],
	    "PoissonExpt" => Any[
				  "Introduction" => "generated/Poisson/doc.md",
				 "Record data" => "generated/Poisson/forw.md",
				 "Born map" => "generated/Poisson/test_born.md",
				 ],
#	    "SeisInvExpt" => Any[]
	    ]
#	    "Seismic Born Modeling" => "FWI/born_map.md",
#	    "Seismic Full Waveform Inversion" => "FWI/gradient_accuracy.md",
#    modules = []
)
 
deploydocs(
    repo   = "github.com/pawbz/GeoPhyInv.jl.git",
)
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
