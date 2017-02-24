module SeismicInversion

# include modules (note: due to dependencies, order is important!)
include("Grid.jl")
include("Wavelets.jl")
include("Models.jl")
include("Acquisition.jl")
include("F90libs.jl")
include("Fdtd.jl")

end # module
