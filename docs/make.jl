using Documenter, JuMIT

#makedocs(modules=[JuMIT],
#	         doctest=true)
makedocs()

#makedocs(
#    format = :html,
#    sitename = "Seismic Inversion Toolbox",
#    modules = []
#)
 
deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/pawbz/JuMIT.jl.git",
    julia  = "0.6",
    osname = "linux",
    target = "build",
    make   = nothing
)
