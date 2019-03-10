using Documenter, JuMIT

#makedocs(modules=[JuMIT],
#	         doctest=true)
#makedocs()

makedocs(
    format = Documenter.HTML(),
   sitename = "Toolbox for Geophysical Inversion",
   pages = Any["index.md",
	       "modeling/page1.md"]
#    modules = []
)
 
deploydocs(
    repo   = "github.com/pawbz/JuMIT.jl.git",
)
#=
deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/pawbz/JuMIT.jl.git",
    julia  = "1.0",
    osname = "linux",
    target = "build",
    make   = nothing
)
=#
