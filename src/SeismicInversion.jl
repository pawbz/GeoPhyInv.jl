module SeismicInversion

# include rest of the modules (note: due to dependencies order is important!)
include("mesh.jl")
include("models.jl")
include("f90libs.jl")
include("fdtd.jl")

end # module
