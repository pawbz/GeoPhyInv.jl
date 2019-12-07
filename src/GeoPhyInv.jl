module GeoPhyInv


# load all necessary packages
using Misfits
using TimerOutputs
using LinearMaps
using LinearAlgebra
using Ipopt
using Optim, LineSearches
using DistributedArrays
using Calculus
using Random
using ProgressMeter
using Distributed
using DistributedArrays
using SharedArrays
using LinearAlgebra
using AxisArrays
using Printf





# create a timer object, used throughout this package, see TimerOutputs.jl
global const to = TimerOutput();

# include modules (note: due to dependencies, order is important!)
include("Interpolation/Interpolation.jl")
include("Utils/Utils.jl")
include("Operators.jl")
include("Smooth.jl")
include("IO.jl")
include("media/medium.jl")
include("srcwav/wavelets.jl")
export ricker, ormsby 

# need to define supersource, source and receiver structs and export them (necessary for multiple dispatch)
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


include("geom/geom.jl")
export Geom, Geomss

include("database/core.jl")

include("srcwav/srcwav.jl")
export SrcWav
SrcWav=Array{NamedD,1}

include("Coupling.jl")
Data=Array{NamedD,1}



# Pressure and velocity fields (used for multiple dispatch)
struct P end
struct Vx end
struct Vz end



#include("Data/Data.jl")
include("Born/Born.jl")
include("fdtd/fdtd.jl")
export SeisForwExpt 

include("fwi/fwi.jl")
export SeisInvExpt, fit!

include("Poisson/Poisson.jl")
#include("gallery/Gallery.jl")
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

end # module
