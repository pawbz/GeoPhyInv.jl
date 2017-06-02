module SIT

# include modules (note: due to dependencies, order is important!)
include("F90libs.jl")
include("IO.jl")
include("Grid.jl")
include("Wavelets.jl")
include("Interpolation.jl")
include("Models.jl")
include("Acquisition.jl")
include("Coupling.jl")
include("Data.jl")
include("DSP.jl")
include("Interferometry.jl")
include("Gallery.jl")
include("Analytic.jl")
include("Fdtd.jl")
include("Misfits.jl")
include("Inversion.jl")
include("CICA.jl")
include("Plots.jl")

end # module
