using Documenter, GeoPhyInv

#makedocs(modules=[GeoPhyInv],
#	         doctest=true)
#makedocs()

makedocs(
	 format = Documenter.HTML(
		  prettyurls = get(ENV, "CI", nothing) == "true"),
   sitename = "Toolbox for Geophysical Inversion",
   pages = ["Home" => "index.md",
	    "Tutorial" => "test/page1.md",]
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
