__precompile__()

module SIT

#const depsfile = joinpath(dirname(@__FILE__), "..", "deps", "deps.jl")
#if isfile(depsfile)
#	include(depsfile)
#else
#	error("SIT not properly installed. Please run Pkg.build(\"SIT\") then restart Julia.")
#end

# include modules (note: due to dependencies, order is important!)
include("F90libs.jl")
include("IO.jl")
include("Grid.jl")
include("Wavelets.jl")
include("Interpolation.jl")
include("DSP.jl")
include("Smooth.jl")
include("Models.jl")
include("Acquisition.jl")
include("Coupling.jl")
include("Data.jl")
include("Interferometry.jl")
include("Gallery.jl")
include("Analytic.jl")
include("Fdtd.jl")
include("Misfits.jl")
include("Inversion.jl")
include("CICA.jl")
include("Plots.jl")

end # module
