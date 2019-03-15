module GeoPhyInv

# include modules (note: due to dependencies, order is important!)
include("Interpolation/Interpolation.jl")
include("Utils/Utils.jl")
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
include("Poisson/Poisson.jl")
include("Gallery/Gallery.jl")
include("Plots.jl")


const JuMIT=GeoPhyInv
export JuMIT

const GPI=GeoPhyInv
export GPI

const J=GeoPhyInv
export J

const JP=GeoPhyInv.Plots
export JP

const JF=GeoPhyInv.FWI
export JF

const JD=GeoPhyInv.Data
export JD

const JA=GeoPhyInv.Acquisition
export JA

const JG=GeoPhyInv.Gallery
export JG

end # module
