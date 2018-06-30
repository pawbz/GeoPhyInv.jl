__precompile__()

module JuMIT

#const depsfile = joinpath(dirname(@__FILE__), "..", "deps", "deps.jl")
#if isfile(depsfile)
#	include(depsfile)
#else
#	error("JuMIT not properly installed. Please run Pkg.build(\"JuMIT\") then restart Julia.")
#end

# include modules (note: due to dependencies, order is important!)
include("Poisson.jl")
include("Operators.jl")
include("Smooth.jl")
include("IO.jl")
include("Models.jl")
include("Acquisition.jl")
include("Coupling.jl")
include("Data.jl")
include("Gallery.jl")
include("Analytic.jl")
include("Fdtd.jl")
include("FWI.jl")
include("Interferometry.jl")
include("SAC.jl")
include("Plots.jl")

end # module
