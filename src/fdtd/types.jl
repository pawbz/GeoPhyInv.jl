# Define high-level structs and their zero initializations 


"""
Struct for modelling parameters specific to each supersource. 
Each worker needs a vector of this struct, as it models multiple superspources. 
"""
mutable struct P_x_worker_x_pw_x_ss{N}
    iss::Int64
    wavelets::Vector{NamedStack{T6}} where {T6<:Data.Array{1}} # [ .. for it in 1:nt]
    ssprayw::NamedStack{T3} where {T3<:AbstractMatrix{Data.Number}}
    records::NamedStack{Vector{T5}} where {T5<:Data.Array{1}}
    rinterpolatew::NamedStack{T4} where {T4<:AbstractMatrix{Data.Number}}
    boundary::NamedStack{NamedStack{Vector{Data.Array{N}}}} # [[[ .. for it in 1:nt] for dim, snap] for field]
    snaps::NamedVector{
        Array{Data.Number,N},
        Vector{Array{Data.Number,N}},
        Tuple{OrderedDict{String,Int64}},
    }
    illum::Matrix{Float64}
    gradients::NamedStack{T7} where {T7<:Data.Array{N}}# gradients stored w.r.t. medium parameters
end

function iszero_boundary(pa)
    result = false
    @sync begin
        for (ip, p) in enumerate(procs(pa.p))
            @async remotecall_wait(p) do
                pap = localpart(pa.p)
                for issp = 1:length(pap.ss)
                    pass = pap.ss[issp]
                    for i = 1:length(pass.boundary)
                        result = result | iszero(pass.boundary[i])
                    end
                end
            end
        end
    end
    return result
end
function initialize_boundary!(pa)
    # zero out results stored per worker
    for ipw = 1:pa.c.ic[:npw]
        @sync begin
            for (ip, p) in enumerate(procs(pa.p))
                @async remotecall_wait(p) do
                    pap = localpart(pa.p)[ipw]
                    for issp = 1:length(pap.ss)
                        pass = pap.ss[issp]
                        map(pass.boundary) do b1
                            map(b1) do b2
                                map(b2) do b3
                                    fill!(b3, zero(Data.Number))
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

"""
Initialize between each simulation
"""
function initialize!(pa::P_x_worker_x_pw_x_ss)
    map(x -> fill!.(x, 0.0), pa.records)
    fill!.(pa.snaps, 0.0)
    fill!(pa.illum, 0.0)
    fill!.(pa.gradients, 0.0)
end


"""
Parameters specific to each worker, not necessarily for every supersource.
This struct contains velocity and stress fields simulated by each worker.
Again, note that a single worker can take care of multiple supersources.
* T1==Matrix{Float64} for 2D simulation
"""
mutable struct P_x_worker_x_pw{N,Q<:Data.Array{N}}
    ss::Vector{P_x_worker_x_pw_x_ss{N}} # supersource specific arrays
    w1::NamedStack{NamedStack{Q}} # fields stored on GPU or CPU, for both t (current) and tp (previous) time steps
    w2::NamedStack{NamedStack{Q}} # required for attenuation, where third dimension is nsls (only used for 2D simulation right now)
    memory_pml::NamedStack{Q} # memory variables for CPML, only derivative fields stored 
    velocity_buffer::NamedStack{Q} # used for born modeling # only used for 2-D simulation
end

# P_x_worker=Vector{P_x_worker_x_pw}

function initialize!(pap::Vector{P_x_worker_x_pw{N,B}}) where {N,B}
    reset_w2!(pap)
    map(pap) do pap_x_pw
        initialize!.(pap_x_pw.ss)
    end

end

# reset wavefields for every worker
function reset_w2!(pap::Vector{P_x_worker_x_pw{N,B}}) where {N,B}
    for pap_x_pw in pap
        fill!.(pap_x_pw.velocity_buffer, zero(Data.Number))
        map(pap_x_pw.w1) do w
            fill!.(w, zero(Data.Number))
        end
        map(pap_x_pw.w2) do w
            fill!.(w, zero(Data.Number))
        end
        fill!.(pap_x_pw.memory_pml, zero(Data.Number))
    end
end

"""
Modelling parameters common for all supersources
# Keyword Arguments that are modified by the method (some of them are returned as well)

* `Records::Vector{Records.TD}=[Records(rfields,tgridmod,ageom[ip]) for ip in 1:length(findn(rflags))]`
* `illum::Array{Float64,2}=zeros(length(medium.mgrid[1]), length(medium.mgrid[2]))` : source energy if `illum_flag`
* `boundary::Array{Array{Float64,4},1}` : stored boundary values for first propagating wavefield 
* `snaps::NamedArray{Array{Float64}}` :snapshots saved at `tsnaps`

# Return (in order)

* modelled data for each propagating wavefield as `Vector{TD}`
* stored boundary values of the first propagating wavefield as `Array{Array{Float64,4},1}` (use for backpropagation)
* final conditions of the first propagating wavefield as `Array{Float64,4}` (use for back propagation)
* gradient model as `Medium`
* stored snaps shots at tsnaps as Array{Float64,4} 

"""
mutable struct P_common{T,N,Q1<:Data.Array{1},Q2<:Data.Array{N}}
    jobname::Symbol
    attrib_mod::T
    exmedium::Medium
    medium::Medium
    ageom::Vector{AGeom}
    srcwav::Vector{SrcWav}
    pml_faces::Vector{Symbol}
    rigid_faces::Vector{Symbol}
    stressfree_faces::Vector{Symbol}
    rfields::Vector{Symbol}
    fc::NamedStack{Data.Number}
    ic::NamedStack{Int64}
    pml::NamedStack{NamedStack{Q1}} # e.g., pml[:dvxdx][:x][:a], pml[:dtauxxdx][:z][:b]
    mod::NamedStack{Q2} # e.g., mod[:invK], mod[:K]
    δmod::NamedStack{Q2} # e.g., mod[:invK], mod[:K]
    # attenuation related parameters
    mod3::NamedStack{Array{Float64,3}} # (only used for 2-D attenuation modelling, so fixed)
    gradients::NamedStack{SharedArrays.SharedArray{Data.Number,N}}
    ref_mod::NamedStack{Data.Number} # reference medium values
    illum_flag::Bool
    illum_stack::SharedArrays.SharedArray{Float64,N}
    backprop_flag::Symbol
    snaps_field::Symbol
    itsnaps::NamedVector{
        Int64,
        Vector{Int64},
        Tuple{OrderedCollections.OrderedDict{String,Int64}},
    }
    datamat::SharedArrays.SharedArray{Data.Number,3}
    data::Vector{Records}
    verbose::Bool
end


function initialize!(pac::P_common)
    fill!.(pac.gradients, 0.0)
    fill!(pac.illum_stack, 0.0)
    map(pac.data) do dat
        fill!(dat, 0.0)
    end
    fill!(pac.datamat, 0.0)
end


"""
A single struct combining parameters for all workers.
"""
mutable struct PFdtd{T,N,B}
    sschunks::Vector{UnitRange{Int64}} # how supersources are distributed among workers
    p::DistributedArrays.DArray{Vector{P_x_worker_x_pw{N,B}},1,Vector{P_x_worker_x_pw{N,B}}} # distributed parameters among workers
    c::P_common{T,N} # common parameters
end

