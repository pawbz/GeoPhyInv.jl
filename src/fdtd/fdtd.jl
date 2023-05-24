
include("types.jl")
include("attenuation.jl")


"""
```julia
pa = SeisForwExpt(attrib_mod; ageom, srcwav, medium, tgrid);
```
Method to create an instance of `SeisForwExpt`. 
The output of this method can be used as an input to the in-place method `update!`, to actually perform a
finite-difference modeling.

# Arguments
* `attrib_mod` : attribute to choose the type of modeling. Choose from 
  * `=FdtdAcoustic()` for acoustic modeling  
  * `=FdtdElastic()` for solving the isotropic elastic wave equation 
  * `=FdtdAcoustic{Born}()` for Born modeling 

# Keyword Arguments
* `medium` : medium parameters 
* `tgrid` : modeling time grid, maximum time in `tgrid` should be greater than or equal to maximum source time, same sampling interval as in `srcwav`
* `ageom` :  acquisition 
* `srcwav` : source wavelets  

# Optional Keyword Arguments 

* `rfields=[:p]` : multi-component receiver flag. Choose `Vector{Symbol}`
  * `=[:p]` record pressure 
  * `=[:vx]` record horizontal component of particle velocity  
  * `=[:vz]` record vertical component of particle velocity  
  * `=[:p, :vx]` record both pressure and velocity 
* `tsnaps` : store snapshots at these modeling times (defaults to half time)
  * `=[0.1,0.2,0.3]` record at these instances of tgrid
* `snaps_field::Symbol=:tauxx` : which field to record? 
* `pml_faces::Vector{Symbol}=[:zmin, :zmax, :ymax, :ymin, :xmax, :xmin]` : control absorbing PML boundary conditions, e.g., in 2D 
  * `=[:zmin, :zmax]` apply PML conditions only at the top and bottom of the medium 
  * `=[:zmax, :xmax, :xmin]` top surface `zmin` is stress free and reflecting
* `rigid_faces::Vector{Symbol}` : rigid faces 
* `stressfree_faces::Vector{Symbol}` : stress-free faces
* `verbose::Bool=false` : verbose flag
"""
SeisForwExpt(
    attrib_mod::Union{FdtdAcoustic,FdtdElastic,FdtdAcousticVisco},
    args1...;
    args2...
) = PFdtd(attrib_mod, args1...; args2...)


"""
Primary method to generate Expt variable when FdtdAcoustic() and FdtdAcoustic{Born}() are used.

# Some internal arguments
* `jobname::Symbol` : name
* `npw::Int64=1` : number of independently propagating wavefields in `medium`
* `backprop_flag::Symbol` : final state variables and the boundary conditions for later use
  * `=:save` save boundary and final values in `boundary` 
  * `=:force` use stored values in `boundary` for back propagation
* `illum_flag::Bool=false` : flag to output wavefield energy or source illumination; it can be used as preconditioner during inversion
"""
function PFdtd(
    attrib_mod;
    jobname::Symbol=:forward_propagation,
    medium::Medium=nothing,
    pml_faces::Vector{Symbol}=[:zmin, :zmax, :ymin, :ymax, :xmax, :xmin],
    rigid_faces=pml_faces,
    tgrid::StepRangeLen=nothing,
    ageom::Union{AGeom,Vector{AGeom}}=nothing,
    srcwav::Union{SrcWav,Vector{SrcWav}}=nothing,
    rfields::Vector{Symbol}=[:vz],
    stressfree_faces=[:dummy],
    backprop_flag::Symbol=:null,
    illum_flag::Bool=false,
    tsnaps::AbstractVector{Float64}=fill(0.5 * (tgrid[end] + tgrid[1]), 1),
    snaps_field=nothing,
    verbose::Bool=false,
    nworker=nothing
)
    N = ndims(medium)

    # get npw from attrib_mod
    npw = attrib_mod.npw

    # convert to vectors of length npw
    if (typeof(ageom) == AGeom) && (typeof(srcwav) == SrcWav)
        if (npw == 2)
            # acquisition geometry for the adjoint wavefield
            ageom = [ageom, get_adjoint_ageom(ageom)]
            # source wavelets for adjoint wavefield
            srcwav = [srcwav, SrcWav(tgrid, get_adjoint_ageom(ageom[1]), rfields)]
        else
            ageom = fill(ageom, npw)
            srcwav = fill(srcwav, npw)
        end
    else
        @assert length(ageom) == length(srcwav)
    end


    # check if fields input are meaningful for given attrib_mod and ndims
    @assert (N == _fd_ndims) "Cannot initiate SeisForwExpt due to ndims inconsistency with @init_parallel_stencil"
    !(snaps_field === nothing) && @assert snaps_field ∈ Fields(attrib_mod, ndims=N)
    foreach(rf -> @assert(rf ∈ Fields(attrib_mod; ndims=N)), rfields)
    foreach(
        s -> foreach(
            ss -> foreach(f -> @assert(f ∈ Fields(attrib_mod, ndims=N)), AxisArrays.names(ss.d)[1]),
            s,
        ),
        srcwav,
    )

    # check sizes and errors based on input
    @assert (length(ageom) == npw)
    @assert (length(srcwav) == npw)
    @assert (maximum(tgrid) >= maximum(srcwav[1][1].grid)) "modeling time is less than source time"
    #(any([getfield(TDout[ip],:tgrid).δx < tgridmod.δx for ip=1:length(TDout)])) && error("output time grid sampling finer than modeling")
    #any([maximum(getfield(TDout[ip],:tgrid).x) > maximum(tgridmod) for ip=1:length(TDout)]) && error("output time > modeling time")

    # all the propagating wavefields should have same supersources, check that
    nss = length(ageom[1])
    (fill(nss, npw) != [length(ageom[ip]) for ip = 1:npw]) && error("different supersources")

    # check if all sources are receivers are inside medium
    any(.![(ageom[ip] ∈ medium.grid) for ip = 1:npw]) && error("sources or receivers not inside medium")

    all([issimilar(ageom[ip], srcwav[ip]) for ip = 1:npw]) || error("ageom and srcwav mismatch")

    #(verbose) &&	println(string("\t> number of super sources:\t",nss))	

    # extend mediums in the PML layers
    exmedium = padarray(medium, _fd_npextend, pml_faces)
    # mod parameters are on the stress grid, so no need to worry
    mod = NamedArray(
        [
            Data.Array(zeros(length.(exmedium.grid)...)) for
            name in MediumParameters(attrib_mod)[1]
        ],
        MediumParameters(attrib_mod)[1],
    )
    # dmod parameters are NOT on stress grid, so need to use get_mgrid function to get sizes
    dmod = NamedArray(
        [
            Data.Array(zeros(length.(get_mgrid(field, attrib_mod, exmedium.grid...))...)) for
            field in MediumParameters(attrib_mod)[3]
        ],
        MediumParameters(attrib_mod)[2],
    )
    δmod = deepcopy(mod) # for perturbations in medium parameters

    # PML
    pml = get_pml(attrib_mod, exmedium.grid)

    # shared arrays required to reduce all the gradient from individual workers
    gradients = NamedArray(
        [
            SharedArray{_fd_datatype}(zeros(length.(exmedium.grid)...)) for
            name in MediumParameters(attrib_mod)[1]
        ],
        MediumParameters(attrib_mod)[1],
    )
    # need reference values for nondimensionalize and dimensionalize
    ref_mod = NamedArray([Data.Number(exmedium[name].ref) for name in MediumParameters(attrib_mod)[1]], MediumParameters(attrib_mod)[1])

    illum_stack = SharedArray{Float64}(zeros(length.(medium.grid)...))

    if (!(snaps_field === nothing))
        @assert eval(snaps_field) <: Fields "invalid snaps field"
        itsnaps = NamedArray(
            [argmin(abs.(tgrid .- tsnaps[i])) for i = 1:length(tsnaps)],
            string.(tsnaps),
        )
    else
        snaps_field = :Nothing
        itsnaps = NamedArray([0], ["0"])
    end

    nrmat = [ageom[ipw][iss].nr for ipw = 1:npw, iss = 1:nss]
    datamat = SharedArray{_fd_datatype}(length(tgrid), maximum(nrmat), nss)
    data = [Records(tgrid, ageom[ip], rfields) for ip = 1:npw]

    # dont need visco parameters, initialize with proper sizes later
    nsls = Int32(0)
    memcoeff1 = zeros(1, 1, 1)
    memcoeff2 = zeros(1, 1, 1)

    fc = get_fc(exmedium, tgrid)
    ic = get_ic(exmedium, tgrid, nsls, npw)

    pac = P_common(
        jobname,
        attrib_mod,
        exmedium,
        medium,
        ageom,
        srcwav,
        pml_faces,
        # pml_faces also need rigid_boundary conditions
        unique(vcat(rigid_faces, pml_faces)),
        stressfree_faces,
        rfields,
        fc,
        ic,
        pml,
        mod,
        dmod,
        δmod,
        NamedArray([memcoeff1, memcoeff2], ([:memcoeff1, :memcoeff2],)),
        gradients,
        ref_mod,
        illum_flag,
        illum_stack,
        backprop_flag,
        snaps_field,
        itsnaps,
        datamat,
        data,
        verbose,
    )


    # update mod arrays in pac
    update!(pac, medium)

    if (first(typeof(attrib_mod).parameters) == Born)
        update!(pac, medium)
    end

    # dividing the supersources to workers
    if (nworker === nothing)
        nworker = min(nss, Distributed.nworkers())
    end
    work = Distributed.workers()[1:nworker]
    ssi = [round(Int, s) for s in range(0, stop=nss, length=nworker + 1)]
    sschunks = Array{UnitRange{Int64}}(undef, nworker)
    for ib = 1:nworker
        sschunks[ib] = ssi[ib]+1:ssi[ib+1]
    end

    # a distributed array of P_x_worker --- note that the parameters for each super source are efficiently distributed here
    papa = ddata(
        T=_fd_use_gpu ?
          Vector{P_x_worker_x_pw{N,CUDA.CuArray{_fd_datatype,N,CUDA.Mem.DeviceBuffer}}} : Vector{P_x_worker_x_pw{N,Array{_fd_datatype,N}}},
        init=I -> Vector{P_x_worker_x_pw}(sschunks[I...][1], pac),
        pids=work[1:1], # disable distributed in order for Pluto to work (waiting for Pluto bug to enable distributed)
    )

    pa = PFdtd(sschunks, papa, pac)

    # update source wavelets with default source type
    update!(pa, pa.c.srcwav, fill(1, length(pa.c.srcwav)), verbose=true)

    # update pml coefficient, now that freqpeak is scanned from the source wavelets
    update_pml!(pa.c)

    check_stability(pa, verbose, H=div(20, _fd_order))


    # update viscoelastic/ viscoacoustic parameters here
    # if(typeof(attrib_mod)==FdtdAcousticVisco)
    # 	nsls=Int32(3)
    # 	exmedium.ic=vcat(exmedium.ic,NamedArray([nsls], ([:nsls],)))
    # 	exmedium.fc=vcat(exmedium.fc,NamedArray(2*pi .* [pa.c.fc[:freqmin], pa.c.fc[:freqmax]], ([:freqmin,:freqmax],)))
    # 	memcoeff1, memcoeff2=get_memcoeff(exmedium)
    # else

    return pa
end

"""
Get some integer constants to store them in pac, e.g., model sizes. ``
The idea is to use them later for a cleaner code.
"""
function get_ic(medium, tgrid, nsls, npw)
    N = ndims(medium)
    return NamedArray(
        vcat(length.(medium.grid), [length(tgrid), nsls, npw]),
        vcat(dim_names(N, "n"), [:nt, :nsls, :npw]),
    )
end

"""
Get some float constants, and store then in pac, e.g, spatial sampling.
The idea is to use them later inside the loops for faster modelling.
"""
function get_fc(medium, tgrid)
    N = ndims(medium)
    ds = step.(medium.grid)
    # denominator depending on _fd_order
    dsI = inv.(ds)
    if (_fd_order == 4)
        dsI = inv.(ds .* 24.0)
    elseif (_fd_order == 6)
        dsI = inv.(ds .* 1920.0)
    elseif (_fd_order == 8)
        dsI = inv.(ds .* 107520.0)
    end
    dt = step(tgrid)
    dtI = inv(dt)

    return map(Data.Number, NamedArray(
        vcat([dt, dtI], ds, dsI),
        vcat([:dt, :dtI], dim_names(N, "d"), dim_names(N, "d", "I")),
    ))
end


"""
Create field arrays for each worker.
Each worker performs the modeling of supersources in `sschunks`.
The parameters common to all workers are stored in `pac`.
"""
function P_x_worker_x_pw(ipw, sschunks::UnitRange{Int64}, pac::P_common{T,N}) where {T,N}
    n = length.(pac.exmedium.grid)


    fields = Fields(pac.attrib_mod) # for current time step
    velocity_fields = filter(x -> x ∉ Fields(pac.attrib_mod, "d"), Fields(pac.attrib_mod, "v"))
    fields_wo_derivatives = filter(x -> x ∉ Fields(pac.attrib_mod, "d"), Fields(pac.attrib_mod)) # remove derivatives, for previous time step
    w1 = NamedArray(
        [NamedArray([zeros(eval(f)(), T(), n...) for f in fields], Symbol.(fields)), # for all fields attached with attrib_mod
            NamedArray([zeros(eval(f)(), T(), n...) for f in fields_wo_derivatives], Symbol.(fields_wo_derivatives)) # for all fields, without derivatives
        ],
        [:t, :tp], # current and previous time step
    )
    velocity_buffer = NamedArray([zeros(eval(f)(), T(), n...) for f in velocity_fields], Symbol.(velocity_fields))

    # dummy (use for viscoelastic modeling later)
    w2 = NamedArray(
        [NamedArray([@zeros(fill(1, N)...) for i in [:r]], ([:r],)) for i = 1:2],
        ([:t, :tp],),
    )

    # memory fields for all derivatives
    dfields = Fields(pac.attrib_mod, "d")
    memory_pml = NamedArray(
        [zeros(eval(f)(), T(), n..., pml=true) for f in dfields],
        Symbol.(dfields),
    )

    ss = [P_x_worker_x_pw_x_ss(ipw, iss, pac) for (issp, iss) in enumerate(sschunks)]

    return P_x_worker_x_pw(ss, w1, w2, memory_pml, velocity_buffer)
end






function Vector{P_x_worker_x_pw}(sschunks::UnitRange{Int64}, pac::P_common)
    return [P_x_worker_x_pw(ipw, sschunks, pac) for ipw = 1:pac.ic[:npw]]
end

"""
Create modeling parameters for each supersource. 
Every worker models one or more supersources.
"""
function P_x_worker_x_pw_x_ss(ipw, iss::Int64, pac::P_common{T,N}) where {T,N}
    rfields = pac.rfields
    sfields = AxisArrays.names(pac.srcwav[ipw][iss].d)[1]
    nt = pac.ic[:nt]
    n = length.(pac.exmedium.grid)
    ageom = pac.ageom
    srcwav = pac.srcwav

    # records_output, distributed array among different procs
    records = NamedArray(
        [[Data.Array(zeros(ageom[ipw][iss].nr)) for it = 1:nt] for i = 1:length(rfields)],
        (rfields,),
    )

    snaps = NamedArray(
        [
            (eval(pac.snaps_field) <: Fields) ?
            Array(zeros(eval(pac.snaps_field)(), T(), n...)) :
            zeros(Data.Number, fill(1, N)...) for i = 1:length(pac.itsnaps)
        ],
        AxisArrays.names(pac.itsnaps)[1],
    )
    # source wavelets
    wavelets = [
        NamedArray(
            [Data.Array(zeros(ageom[ipw][iss].ns)) for i = 1:length(sfields)],
            (sfields,),
        ) for it = 1:nt
    ]

    # boundary fields stored (after removing derivatives)
    # these fields are stored for the boundary value problem (used for RTM and FWI)
    bfields = filter(x -> x ∉ Fields(pac.attrib_mod, "d"), Fields(pac.attrib_mod))
    boundary = NamedArray(
        [get_boundary_store(eval(f)(), T(), n..., nt) for f in bfields],
        Symbol.(bfields),
    )

    # initialize source_spray_weights per supersource and receiver interpolation weights per sequential source
    ssprayw = NamedArray(
        [
            spzeros(
                Data.Number,
                prod(length.(get_mgrid(eval(sf)(), T(), pac.exmedium.grid...))),
                ageom[ipw][iss].ns,
            ) for sf in sfields
        ],
        (sfields,),
    )
    rinterpolatew = NamedArray(
        [
            spzeros(
                Data.Number,
                prod(length.(get_mgrid(eval(rf)(), T(), pac.exmedium.grid...))),
                ageom[ipw][iss].nr,
            ) for rf in rfields
        ],
        rfields,
    )

    # transfer to GPUs if 
    if (_fd_use_gpu)
        ssprayw = map(x -> CuSparseMatrixCSC(x), ssprayw)
        rinterpolatew = map(x -> CuSparseMatrixCSC(x), rinterpolatew)
    end

    # named array to store gradients for each supersource
    gradients = NamedArray(
        [Data.Array(zeros(length.(pac.exmedium.grid)...)) for name in MediumParameters(pac.attrib_mod)[1]],
        MediumParameters(pac.attrib_mod)[1],
    )

    # saving illum (dummy)
    # illum =  (pac.illum_flag) ? zeros(nz, nx) : zeros(1,1)
    illum = zeros(1, 1)

    pass = P_x_worker_x_pw_x_ss(
        iss,
        wavelets,
        ssprayw,
        records,
        rinterpolatew,
        boundary,
        snaps,
        illum,
        gradients
    )

    # update acquisition
    update!(pass, ipw, iss, ageom[ipw][iss], pac)
    return pass
end

include("source.jl")
include("receiver.jl")
include("dirichlet.jl")
include("cpml.jl")
include("save_tp.jl")
include("advance_acou.jl")
include("advance_elastic.jl")
include("gradient.jl")
include("born.jl")
include("boundary.jl")


# update TDout after forming a vector and resampling
#	ipropout=0;
#	for iprop in 1:pac.ic[:npw]
#			ipropout += 1
##			Records.TD_resamp!(pac.data[ipropout], Records.TD_urpos((Array(records[:,:,iprop,:,:])), rfields, tgridmod, ageom[iprop],
##				ageom_urpos[1].nr[1],
##				(ageom_urpos[1].r[:z][1], ageom_urpos[1].r[:x][1])
##				)) 
#		end
#	end
# return without resampling for testing
#return [Records.TD(reshape(records[1+(iprop-1)*nd : iprop*nd],tgridmod.nx,recv_n,nss),
#		       tgridmod, ageom[1]) for iprop in 1:npw]

function stack_illums!(pac::P_common, pap::Vector{P_x_worker_x_pw{N}}) where {N}
    nx, nz = pac.ic[:nx], pac.ic[:nz]
    illums = pac.illum_stack
    pass = pap[1].ss
    for issp = 1:length(pass)
        gs = pass[issp].illum
        gss = view(gs, _fd_npml+1:nz-_fd_npml, _fd_npml+1:nx-_fd_npml)
        @. illums += gss
    end
end




# Need illumination to estimate the approximate diagonal of Hessian
@inbounds @fastmath function compute_illum!(
    issp::Int64,
    pap::Vector{P_x_worker_x_pw{N}},
) where {N}
    # saving illumination to be used as preconditioner 
    p = pap[1].w1[:t][:p]
    illum = pap[1].ss[issp].illum
    for i in eachindex(illum)
        illum[i] += abs2(p[i])
    end
end


"""
Save snapshots for ipw=1
"""
function snaps_save!(itsnap, issp::Int64, ipw::Int, pac::P_common, pap)
    p = pap[ipw].w1[:t][pac.snaps_field]
    snaps = pap[ipw].ss[issp].snaps[itsnap]
    copyto!(snaps, p)
end



include("stability.jl")
include("medium.jl")
include("ageom.jl")
include("propagate.jl")
include("getprop.jl")

