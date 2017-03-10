module SIT

# include modules (note: due to dependencies, order is important!)
include("F90libs.jl")
include("IO.jl")
include("Grid.jl")
include("Wavelets.jl")
include("Models.jl")
include("Acquisition.jl")
include("Fdtd.jl")
include("DSP.jl")
include("Banded_ICA.jl")

end # module
