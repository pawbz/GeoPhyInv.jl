module GeoPhyInv

# required packages
using PrecompileTools
using Preferences
using Distributed
using ParallelStencil
using Distances
using TimerOutputs
using LinearMaps
using Luxor
using Ipopt
using IterativeSolvers
using Optim, LineSearches
using DistributedArrays
using Calculus
using ProgressMeter
using SharedArrays
using Printf
using Markdown
using DataFrames
using Measures
using SparseArrays
using OrderedCollections
using Interpolations: interpolate, extrapolate, Gridded, Linear, Flat
using CSV
using Statistics
using LossFunctions
using LinearAlgebra
using Random
using ImageFiltering
using NamedArrays
using Test
using AxisArrays
using Distributions
using StatsBase
using InteractiveUtils
using LaTeXStrings
using RecipesBase
using ColorSchemes
using FFTW
using HDF5
using SpecialFunctions
using DSP
using FourierTools
using CUDA
using MLUtils
using CUDA.CUSPARSE
import Base.@doc
CUDA.allowscalar(false)

# module to do convolutions using FFTW
include("Conv.jl")

"""
Setup ParallelStencil preferences.
```julia
set_preferences(; ndims, use_gpu, datatype, order)
```

# Arguments

* `ndims::Int=2`: the number of dimensions used for the stencil computations 2D or 3D 
* `use_gpu::Bool=false`: use GPU for stencil computations or not
* `datatype="Float32"`: the type of numbers used by field arrays (e.g. "Float32" or "Float64")
* `order::Int=2`: order ∈ [2,4,6,8] of the finite-difference stencil 
"""
function set_preferences(; ndims=2, use_gpu=false, datatype="Float32", order=2, verbose=false)
    @assert ndims ∈ [2, 3]
    @assert order ∈ [2, 4, 6, 8]
    @assert isa(use_gpu, Bool)
    @assert isa(verbose, Bool)
    @assert datatype ∈ ["Float32", "Float64"]
    # Set it in our runtime values, as well as saving it to disk
    @set_preferences!("ndims" => ndims)
    @set_preferences!("verbose" => verbose)
    @set_preferences!("order" => order)
    @set_preferences!("use_gpu" => use_gpu)
    @set_preferences!("datatype" => datatype)

    @info("GeoPhyInv: new parallel stencil preferences set; restart your Julia session for this change to take effect!")
end

const verbose = @load_preference("verbose", false)
const _fd_order = @load_preference("order", 2)
const _fd_ndims = @load_preference("ndims", 2)
const _fd_datatype = eval(Symbol(@load_preference("datatype", "Float32")))
const _fd_use_gpu = @load_preference("use_gpu", false)
const _fd_npml = 20 + (_fd_order - 1)
const _fd_npextend = 20 + (_fd_order - 1) # determines exmedium
const _fd_nbound = 3 # number of points to store the boundary values


ParallelStencil.is_initialized() && ParallelStencil.@reset_parallel_stencil()
@static if _fd_use_gpu
    ParallelStencil.@init_parallel_stencil(CUDA, _fd_datatype, _fd_ndims)
else
    ParallelStencil.@init_parallel_stencil(Threads, _fd_datatype, _fd_ndims)
end

# create a timer object, used throughout this package, see TimerOutputs.jl
global const to = TimerOutput();

"""
Return axis names of 1D, 2D or 3D fields
"""
function dim_names(N, prefix = "", suffix = ""; order = 1)
    if (N == 3)
        names = (order == 2) ? [:zz, :yy, :xx, :xz, :xy, :yz] : [:z, :y, :x]
    elseif (N == 2)
        names = (order == 2) ? [:zz, :xx, :xz] : [:z, :x]
    elseif (N == 1)
        names = [:z]
    else
        error("invalid dim num")
    end
    return [Symbol(prefix, string(m), suffix) for m in names]
end

# This type is extensively used to create named arrays with Symbols
NamedStack{T} =
    NamedArray{T,1,Array{T,1},Tuple{OrderedCollections.OrderedDict{Symbol,Int64}}}

# a module with some util methods
include("Utils/Utils.jl")
include("proj_mat.jl")

# supersources used in both database and ageom
struct SSrcs
    n::Int
end
SSrcs() = SSrcs(0)
include("database/database.jl")
export SSrcs, Srcs, Recs

# some structs used for multiple dispatch throughout this package
include("physics_types.jl")
export FdtdAcoustic, FdtdElastic, Born, FullWave, AcousticBorn, ElasticBorn, FdtdAcousticVisco

# mutable data type for seismic medium + related methods

global const marmousi_folder=joinpath(dirname(pathof(GeoPhyInv)), "media", "marmousi2")
global const overthrust_folder=joinpath(dirname(pathof(GeoPhyInv)), "media", "overthrust")
include("media/media.jl")
export Medium, AcousticMedium, ElasticMedium

# source-receiver geometry
include("ageom/ageom.jl")
export AGeom, AGeomss


# store records
include("records/records.jl")

# mutable data type for storing time-domain source wavelets
include("srcwav/wavelets.jl")
export ricker, ormsby
include("srcwav/srcwav.jl")


@static if _fd_ndims == 3
    include("fdtd/diff3D.jl")
else
    include("fdtd/diff2D.jl")
end

# define structs for wavefields in 2D/3D
include("fields.jl")
export Fields
include("fdtd/fdtd.jl")
include("fdtd/func_grad.jl")
include("born/core.jl")

export SeisForwExpt

# export update from GeoPhyInv, as it is commonly used
export update, update!


include("fwi/fwi.jl")
include("fwi/func_grad.jl")
export SeisInvExpt

# ============= Gallery =================================
abstract type Gallery end
struct Homogeneous <: Gallery end
struct RandScatterer <: Gallery end
struct Marmousi2 <: Gallery end
struct Overthrust <: Gallery end
export Homogeneous, RandScatterer, Marmousi2, Overthrust
include("ageom/gallery.jl")
include("media/gallery.jl")
include("fdtd/gallery.jl")
# =======================================================


# include("Poisson/Poisson.jl")
include("plots.jl")





# ============= Precompile =================================
@setup_workload begin
    @show CUDA.functional()
    @show _fd_use_gpu
    @show Data.Array(randn(2,2))
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.
    pa_test = SeisForwExpt(FdtdAcoustic{FullWave}(:forward, 1), RandScatterer())
    @compile_workload begin
        update!(pa_test);
    end
end


end # module GeoPhyInv



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

# # include modules (note: due to dependencies, order is important!)
# include("Operators.jl")
# include("Smooth.jl")
# include("IO.jl")

# include("Coupling.jl")
