using Documenter, SIT
makedocs()


#makedocs(
#    format = :html,
#    sitename = "Seismic Inversion Toolbox",
#    modules = []
#)
 
deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/pawbz/SIT.jl.git",
    julia  = "0.5",
    osname = "linux",
    target = "build",
    make   = nothing
)
