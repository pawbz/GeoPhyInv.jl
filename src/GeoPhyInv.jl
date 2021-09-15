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
using ParallelStencil
using ParallelStencil.FiniteDifferences3D
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
using AxisArrays
using Distributions
using StatsBase
using InteractiveUtils
using RecipesBase
using FFTW
using CUDA
using HDF5



USE_GPU=true
@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 3);
else
    @init_parallel_stencil(Threads, Float64, 3);
end
# import ..USE_GPU
# import ..NDIMS
# check whether 2D or 3D, and initialize ParallelStencils accordingly
# ParallelStencil.@reset_parallel_stencil()
# USE_GPU ? @init_parallel_stencil(CUDA, Float64, NDIMS) : @init_parallel_stencil(Threads, Float64, NDIMS)
# @static @init_parallel_stencil(CUDA, Float64, 3)
# use_gpu ? @init_parallel_stencil(CUDA, Float64, 2) : @init_parallel_stencil(Threads, Float64, 2)

# import ParallelStencil.FiniteDifferences3D
# import ParallelStencil.Data
# import ParallelStencil.FiniteDifferences3D


# This is extensively used to group arrays

NamedStack{T}=NamedArray{T,1,Array{T,1},Tuple{OrderedCollections.OrderedDict{Symbol,Int64}}}


# create a timer object, used throughout this package, see TimerOutputs.jl
global const to = TimerOutput();


struct FdtdElastic end

"""
2-D acoustic forward modeling using a finite-difference simulation of the acoustic wave-equation.
"""
struct FdtdAcou end

"""
2-D ViscoAcoustic forward modeling using a finite-difference simulation of the acoustic wave-equation.
"""
struct FdtdAcouVisco  end

"""
2-D Linearized forward modeling using a finite-difference simulation of the acoustic wave-equation.
"""
struct FdtdAcouBorn  end


global const npml = 20

# define structs for wavefields in 2D/3D
include("fields.jl")

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


include("ageom/core.jl")

include("database/core.jl")

include("srcwav/core.jl")

include("Coupling.jl")
include("records/core.jl")





#include("Records/Records.jl")
include("Born/Born.jl")
include("fdtd/fdtd.jl")

include("fwi/fwi.jl")

include("Poisson/Poisson.jl")
include("plots.jl")

# export stuff from GeoPhyInv
export Records
export SrcWav
export update!, Medium
export ricker, ormsby 
export Srcs, Recs, SSrcs
export AGeom, AGeomss
export update, update!, padarray, padarray!, SeisForwExpt, SeisInvExpt, FdtdAcou, FdtdElastic, FdtdAcouBorn, FdtdAcouVisco, LS, LS_prior, Migr, Migr_FD

# export the Expt for Poisson
const PoissonExpt=GeoPhyInv.Poisson.ParamExpt
export PoissonExpt
mod!(a::PoissonExpt,b,c)=GeoPhyInv.Poisson.mod!(a,b,c)
mod!(a::PoissonExpt,b)=GeoPhyInv.Poisson.mod!(a,b)
mod!(a::PoissonExpt)=GeoPhyInv.Poisson.mod!(a)
operator_Born(a::PoissonExpt,b)=GeoPhyInv.Poisson.operator_Born(a,b)


end # module
