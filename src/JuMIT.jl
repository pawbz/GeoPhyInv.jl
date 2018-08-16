module JuMIT

#const depsfile = joinpath(dirname(@__FILE__), "..", "deps", "deps.jl")
#if isfile(depsfile)
#	include(depsfile)
#else
#	error("JuMIT not properly installed. Please run Pkg.build(\"JuMIT\") then restart Julia.")
#end

# include modules (note: due to dependencies, order is important!)
include("Utils/Utils.jl")
include("Poisson.jl")
include("Operators.jl")
include("Smooth.jl")
include("IO.jl")
include("Models/Models.jl")
include("Acquisition/Acquisition.jl")
include("Coupling.jl")
include("Data/Data.jl")
include("Born/Born.jl")
include("Fdtd/Fdtd.jl")
include("FWI/FWI.jl")
include("Gallery/Gallery.jl")
include("Interferometry.jl")
include("Plots.jl")

const J=JuMIT
export J

const JP=JuMIT.Plots
export JP

const JF=JuMIT.FWI
export JF

const JD=JuMIT.Data
export JD

const JA=JuMIT.Acquisition
export JA

const JG=JuMIT.Gallery
export JG

end # module
