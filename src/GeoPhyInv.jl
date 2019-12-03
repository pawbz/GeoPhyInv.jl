module GeoPhyInv

# include modules (note: due to dependencies, order is important!)
include("Interpolation/Interpolation.jl")
include("Utils/Utils.jl")
include("Operators.jl")
include("Smooth.jl")
include("IO.jl")
include("Media/Media.jl")
include("Acquisition/Acquisition.jl")
include("Coupling.jl")
include("Data/Data.jl")
include("Born/Born.jl")
include("Fdtd/Fdtd.jl")
include("FWI/FWI.jl")
include("Poisson/Poisson.jl")
include("Gallery/Gallery.jl")
include("Plots.jl")

# export the Expt for Seismic Forward Modelling
const SeisForwExpt=GeoPhyInv.Fdtd.Param
export SeisForwExpt
mod!(a::SeisForwExpt)=GeoPhyInv.Fdtd.mod!(a)

# export the Expt for Seismic Inversion
const SeisInvExpt=GeoPhyInv.FWI.Param
export SeisInvExpt
export fit!

# export the Expt for Poisson
const PoissonExpt=GeoPhyInv.Poisson.ParamExpt
export PoissonExpt
mod!(a::PoissonExpt,b,c)=GeoPhyInv.Poisson.mod!(a,b,c)
mod!(a::PoissonExpt,b)=GeoPhyInv.Poisson.mod!(a,b)
mod!(a::PoissonExpt)=GeoPhyInv.Poisson.mod!(a)
operator_Born(a::PoissonExpt,b)=GeoPhyInv.Poisson.operator_Born(a,b)

export mod!, operator_Born

const JuMIT=GeoPhyInv
export JuMIT

const GIPh=GeoPhyInv
export GIPh

const J=GeoPhyInv
export J

const GInterp=GeoPhyInv.Interpolation
export GInterp

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
