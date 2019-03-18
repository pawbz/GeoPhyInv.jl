using Documenter, GeoPhyInv

#makedocs(modules=[GeoPhyInv],
#	         doctest=true)
#makedocs()

makedocs(
	 format = Documenter.HTML(
		  prettyurls = get(ENV, "CI", nothing) == "true"),
   sitename = "Toolbox for Geophysical Inversion",
   pages = ["Home" => "index.md",
	    "Seismic Born Modeling" => "FWI/born_map.md",
	    "Seismic Full Waveform Inversion" => "FWI/gradient_accuracy.md",
	    "Poisson Solver" => "Poisson/adj_state_expt.md",]
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
