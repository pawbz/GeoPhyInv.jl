# Documentation

makedocs(
    format = :html,
    sitename = "Seismic Inversion",
    modules = [SeismicInversion.fdtd]
)
 
deploydocs(
    repo   = "github.com/pawbz/SeismicInversion.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing
)
