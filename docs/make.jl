using Documenter, GeoPhyInv

#makedocs(modules=[GeoPhyInv],
#	         doctest=true)
#makedocs()

makedocs(
	 format = Documenter.HTML(
		  prettyurls = get(ENV, "CI", nothing) == "true"),
   sitename = "Toolbox for Geophysical Inversion",
   pages = ["Home" => "index.md",
	    "Medium" => "media/doc.md",
	    "Geom" => "geom/doc.md",
	    "SrcWav" => "srcwav/doc.md",
	    "Data" => "data/doc.md",
	    "SeisForwExpt" => Any[
				  "Introduction" => "Fdtd/intro.md",
				  "Basic usage" => "Fdtd/reuse_expt.md",
				  "Generate snaps" => "Fdtd/create_snaps.md",
				  ],
	    "PoissonExpt" => Any[
				  "Introduction" => "Poisson/intro.md",
				 "Record data" => "Poisson/forw.md",
				 "Born map" => "Poisson/test_born.md",
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
