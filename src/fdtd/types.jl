# Define high-level structs and their zero initializations 


"""
Struct for modelling parameters specific to each supersource. 
Each worker needs a vector of this struct, as it models multiple superspources. 
"""
mutable struct P_x_worker_x_pw_x_ss{N}
    iss::Int64
    wavelets::Vector{NamedStack{Vector{Float64}}} # [ .. for it in 1:nt]
    ssprayw::NamedStack{Vector{Vector{Float64}}}
    records::NamedStack{Matrix{Float64}}
    rinterpolatew::NamedStack{Vector{Vector{Float64}}}
    sindices::NamedStack{Vector{T1}} where T1<:CartesianIndices{N} # contains indices for each source, each named field
    rindices::NamedStack{Vector{T2}} where T2<:CartesianIndices{N}# contains indices for each receiver, each named field
    boundary::Vector{Array{Float64,3}}
    snaps::NamedVector{
        Array{Data.Number,N},
        Vector{Array{Data.Number,N}},
        Tuple{OrderedDict{String,Int64}},
    }
    illum::Matrix{Float64}
    grad_mod::NamedStack{Matrix{Float64}} # w.r.t different coeffs
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
                        for i = 1:length(pass.boundary)
                            fill!(pass.boundary[i], 0.0)
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
    fill!.(pa.records, 0.0)
    fill!.(pa.snaps, 0.0)
    fill!(pa.illum, 0.0)
    fill!.(pa.grad_mod, 0.0)
end


"""
Parameters specific to each worker, not necessarily for every supersource.
This struct contains velocity and stress fields simulated by each worker.
Again, note that a single worker can take care of multiple supersources.
* T1==Matrix{Float64} for 2D simulation
"""
mutable struct P_x_worker_x_pw{N,Q<:Data.Array{N}}
    ss::Vector{P_x_worker_x_pw_x_ss{N}} # supersource specific arrays
    # 2D arrays of p, vx, vz, their previous snapshots, and 
    # x derivatives of p, vx, vz
    # z derivatives of p, vx, vz
    w1::NamedStack{NamedStack{Q}} # p, vx, vz on GPU or CPU
    wr::NamedStack{Array{Data.Number,N}} # same as above but only for receiver fields when GPU to facilitate faster scalar indexing!?
    w2::NamedStack{NamedStack{Q}} # required for attenuation, where third dimension is nsls (only used for 2D simulation right now)
    memory_pml::NamedStack{Q} # memory variables for CPML 
    born_svalue_stack::Array{Float64,N} # used for born modeling # only used for 2-D simulation
end

# P_x_worker=Vector{P_x_worker_x_pw}

function initialize!(pap::Vector{P_x_worker_x_pw{N,B}}) where {N,B}
    reset_w2!(pap)
    for pap_x_pw in pap
        initialize!.(pap_x_pw.ss)
    end

end

# reset wavefields for every worker
function reset_w2!(pap::Vector{P_x_worker_x_pw{N,B}}) where {N,B}
    for pap_x_pw in pap
        fill!(pap_x_pw.born_svalue_stack, 0.0)
        for w in pap_x_pw.w1
            fill!.(w, 0.0)
        end
        for w in pap_x_pw.w2
            fill!.(w, 0.0)
        end
        fill!.(pap_x_pw.memory_pml, 0.0)
    end
end

"""
Modelling parameters common for all supersources
# Keyword Arguments that are modified by the method (some of them are returned as well)

* `gradient::Vector{Float64}=Medium(medium.mgrid)` : gradient medium modified only if `gmodel_flag`
* `TDout::Vector{Records.TD}=[Records.TD_zeros(rfields,tgridmod,ageom[ip]) for ip in 1:length(findn(rflags))]`
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
    activepw::Vector{Int64}
    exmedium::Medium
    medium::Medium
    ageom::Vector{AGeom}
    srcwav::Vector{SrcWav}
    pml_edges::Vector{Symbol}
    sfields::Vector{Vector{Symbol}}
    sflags::Vector{Int64}
    rfields::Vector{Symbol}
    rflags::Vector{Int64}
    fc::NamedStack{Float64}
    ic::NamedStack{Int64}
    pml::NamedStack{NamedStack{Q1}} # e.g., pml[:x][:a], pml[:z][:b]
    mod::NamedStack{Q2} # e.g., mod[:KI], mod[:K]
    mod3::NamedStack{Array{Float64,3}} # (only used for 2-D attenuation modelling, so fixed)
    #=
    	modKI::Matrix{Float64}
    	modK::Matrix{Float64} # just storing inv(modKI) for speed
    	modrr::Matrix{Float64}
    	modrhovxI::Matrix{Float64}
    	modrhovzI::Matrix{Float64}
    	=#
    δmod::NamedStack{Array{Float64,N}}
    δmodall::Vector{Float64}
    #=
    	δmodKI::Matrix{Float64}
    	δmodrr::Matrix{Float64} 
    	δmodrhovxI::Matrix{Float64}
    	δmodrhovzI::Matrix{Float64}
    	δmod::Vector{Float64} # perturbation vector (KI, rhoI)
    	=#
    gradient::Vector{Float64}  # output gradient vector w.r.t (KI, rhoI)
    grad_mod::NamedStack{SharedArrays.SharedArray{Float64,N}}
    #=
    	grad_modKI_stack::SharedArrays.SharedArray{Float64,2} # contains gmodKI
    	grad_modrhovxI_stack::SharedArrays.SharedArray{Float64,2}
    	grad_modrhovzI_stack::SharedArrays.SharedArray{Float64,2}
    	grad_modrr_stack::Array{Float64,2}
    	=#
    illum_flag::Bool
    illum_stack::SharedArrays.SharedArray{Float64,N}
    backprop_flag::Int64
    snaps_field::Symbol
    itsnaps::NamedVector{
        Int64,
        Vector{Int64},
        Tuple{OrderedCollections.OrderedDict{String,Int64}},
    }
    gmodel_flag::Bool
    bindices::NamedStack{Int64}
    #=
    	ibx0::Int64
    	ibz0::Int64
    	ibx1::Int64
    	ibz1::Int64
    	isx0::Int64
    	isz0::Int64
    	=#
    datamat::SharedArrays.SharedArray{Float64,3}
    data::Vector{Records}
    # attenuation related parameters
    verbose::Bool
end


function initialize!(pac::P_common)
    fill!(pac.gradient, 0.0)
    fill!.(pac.grad_mod, 0.0)
    # fill!(pac.grad_mod[:rhovxI],0.0)
    # fill!(pac.grad_mod[:rhovzI],0.0)
    # fill!(pac.grad_mod[:rhoI],0.0)
    fill!(pac.illum_stack, 0.0)
    for dat in pac.data
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

