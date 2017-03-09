using Documenter, SeismicInversion
makedocs()


makedocs(
    format = :html,
    sitename = "Seismic Inversion",
    modules = [SeismicInversion.Fdtd]
)
 
deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/pawbz/SeismicInversion.jl.git",
    julia  = "0.5",
    osname = "linux",
    target = "build",
    make   = nothing
)
