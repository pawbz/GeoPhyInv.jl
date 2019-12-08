module GeoPhyInv


# load all necessary packages
using Misfits
using TimerOutputs
using LinearMaps
using Ipopt
using Optim, LineSearches
using DistributedArrays
using Calculus
using ProgressMeter
using Distributed
using DistributedArrays
using SharedArrays
using Printf
using DataFrames
using SparseArrays
using Interpolations
using OrderedCollections
using CSV
using Statistics
using LinearAlgebra
using Random
using ImageFiltering
using NamedArrays
using DSP
using Test
using FFTW
using AxisArrays
using Distributions










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


include("geom/geom.jl")

include("database/core.jl")

include("srcwav/srcwav.jl")

include("Coupling.jl")
include("data/core.jl")


# Pressure and velocity fields (used for multiple dispatch)
struct P end
struct Vx end
struct Vz end



#include("Data/Data.jl")
include("Born/Born.jl")
include("fdtd/fdtd.jl")

include("fwi/fwi.jl")

include("Poisson/Poisson.jl")
#include("gallery/Gallery.jl")
include("Plots.jl")

# export stuff from GeoPhyInv
export Data
export SrcWav
export update!, Medium
export ricker, ormsby 
export Srcs, Recs, SSrcs
export Geom, Geomss
export update!, SeisForwExpt, SeisInvExpt, Fdtd, FdtdBorn, LS, LS_prior

# export the Expt for Poisson
const PoissonExpt=GeoPhyInv.Poisson.ParamExpt
export PoissonExpt
mod!(a::PoissonExpt,b,c)=GeoPhyInv.Poisson.mod!(a,b,c)
mod!(a::PoissonExpt,b)=GeoPhyInv.Poisson.mod!(a,b)
mod!(a::PoissonExpt)=GeoPhyInv.Poisson.mod!(a)
operator_Born(a::PoissonExpt,b)=GeoPhyInv.Poisson.operator_Born(a,b)


end # module
