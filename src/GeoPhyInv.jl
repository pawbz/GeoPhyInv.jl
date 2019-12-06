module GeoPhyInv

# include modules (note: due to dependencies, order is important!)
include("Interpolation/Interpolation.jl")
include("Utils/Utils.jl")
include("Operators.jl")
include("Smooth.jl")
include("IO.jl")
include("Media/medium.jl")
include("SrcWav/wavelets.jl")
export ricker, ormsby 

# need to define supersource, source and receiver structs and export them
struct Srcs
	n::Int
end
Srcs()=Srcs(0)
struct SSrcs
	n::Int
end
SSrcs()=SSrcs(0)
struct Recs
	n::Int
end
Recs()=Recs(0)
export Srcs, Recs, SSrcs

include("SrcWav/srcwav.jl")
export SrcWav
include("Geom/geom.jl")
export Geom, Geomss
include("Coupling.jl")
include("Data/Data.jl")
include("Born/Born.jl")
include("Fdtd/fdtd.jl")
export SeisForwExpt 
include("FWI/fwi.jl")
include("Poisson/Poisson.jl")
include("Gallery/Gallery.jl")
include("Plots.jl")

# export the Expt for Seismic Forward Modelling

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

const JG=GeoPhyInv.Gallery
export JG

end # module
