module JuMIT

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
