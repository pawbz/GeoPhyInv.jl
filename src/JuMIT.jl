__precompile__()

module JuMIT

const JT = JuMIT
const SIT = JuMIT
export JT, SIT

#const depsfile = joinpath(dirname(@__FILE__), "..", "deps", "deps.jl")
#if isfile(depsfile)
#	include(depsfile)
#else
#	error("JuMIT not properly installed. Please run Pkg.build(\"JuMIT\") then restart Julia.")
#end

# include modules (note: due to dependencies, order is important!)
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
