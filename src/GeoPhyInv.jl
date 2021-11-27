module GeoPhyInv

# required packages
using ParallelStencil
using Distances
using TimerOutputs
using LinearMaps
using Ipopt
using Optim, LineSearches
using DistributedArrays
using Calculus
using ProgressMeter
using Distributed
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
using Test
using AxisArrays
using Distributions
using StatsBase
using InteractiveUtils
using RecipesBase
using ColorSchemes
using FFTW
using HDF5
using SpecialFunctions
using DSP
using CUDA.CUSPARSE
CUDA.allowscalar(false)



include("params.jl")

# This type is extensively used to create named arrays with Symbols
NamedStack{T} =
    NamedArray{T,1,Array{T,1},Tuple{OrderedCollections.OrderedDict{Symbol,Int64}}}

# a module with some util methods
include("Utils/Utils.jl")
include("Interpolation/Interpolation.jl")

# some structs used for multiple dispatch throughout this package
include("types.jl")
export Srcs, Recs, SSrcs
export FdtdAcoustic, FdtdElastic, FdtdAcousticBorn, AcousticBorn, ElasticBorn, FdtdAcousticVisco

# mutable data type for seismic medium + related methods
include("media/core.jl")
export Medium

# source-receiver geometry
include("ageom/core.jl")
export AGeom, AGeomss

# 
include("database/core.jl")

# store records
include("records/core.jl")
export Records

# mutable data type for storing time-domain source wavelets
include("srcwav/wavelets.jl")
export ricker, ormsby
include("srcwav/core.jl")
export SrcWav



"""
Initialize the package with ParallelStencil, giving access to its main functionality. 
```julia
@init_parallel_stencil(ndims, use_gpu, datatype, order)
```

# Arguments

* `ndims::Int`: the number of dimensions used for the stencil computations 2D or 3D 
* `use_gpu::Bool` : use GPU for stencil computations or not
* `datatype`: the type of numbers used by field arrays (e.g. Float32 or Float64)
* `order::Int ∈ [2,4,6,8]` : order of the finite-difference stencil 
"""
macro init_parallel_stencil(ndims::Int, use_gpu::Bool, datatype, order)
    quote
        # using ParallelStencil
        # ParallelStencil.is_initialized() && ParallelStencil.@reset_parallel_stencil()
        # ParallelStencil.@init_parallel_stencil(CUDA, $datatype, $ndims)
        @eval GeoPhyInv begin
            _fd.ndims = $ndims
            @assert _fd.ndims ∈ [2, 3]
            _fd.order = $order
            @assert _fd.order ∈ [2, 4, 6, 8]
            _fd.use_gpu = $use_gpu
            _fd.datatype = $datatype
            @assert _fd.datatype ∈ [Float32, Float64]
            # extend by more points as we increase order, not sure if that is necessary!!!
            _fd.npml = 20 + ($order - 1)
            _fd.npextend = 20 + ($order - 1) # determines exmedium

            ParallelStencil.is_initialized() && ParallelStencil.@reset_parallel_stencil()
            # check whether 2D or 3D, and initialize ParallelStencils accordingly
            @static if $use_gpu
                # using CUDA
                ParallelStencil.@init_parallel_stencil(CUDA, $datatype, $ndims)
            else
                ParallelStencil.@init_parallel_stencil(Threads, $datatype, $ndims)
            end
            include(
                joinpath(
                    dirname(pathof(GeoPhyInv)),
                    "fdtd",
                    string("diff", $ndims, "D.jl"),
                ),
            )
            # define structs for wavefields in 2D/3D
            include(joinpath(dirname(pathof(GeoPhyInv)), "fields.jl"))
            export Fields
            include(joinpath(dirname(pathof(GeoPhyInv)), "fdtd", "core.jl"))
            include(joinpath(dirname(pathof(GeoPhyInv)), "born", "core.jl"))
            nothing
        end
    end
end
export @init_parallel_stencil
export SeisForwExpt

# export update from GeoPhyInv, as it is commonly used
export update, update!

# # include modules (note: due to dependencies, order is important!)
# include("Operators.jl")
# include("Smooth.jl")
# include("IO.jl")

# include("Coupling.jl")

# include("fwi/fwi.jl")

# include("Poisson/Poisson.jl")
include("plots.jl")




#     padarray,
#     padarray!,
#     SeisInvExpt,
#     LS,
#     LS_prior,
#     Migr,
#     Migr_FD

# # export the Expt for Poisson
# const PoissonExpt=GeoPhyInv.Poisson.ParamExpt
# export PoissonExpt
# mod!(a::PoissonExpt, b, c) = GeoPhyInv.Poisson.mod!(a, b, c)
# mod!(a::PoissonExpt, b) = GeoPhyInv.Poisson.mod!(a, b)
# mod!(a::PoissonExpt) = GeoPhyInv.Poisson.mod!(a)
# operator_Born(a::PoissonExpt, b) = GeoPhyInv.Poisson.operator_Born(a, b)


end # module
