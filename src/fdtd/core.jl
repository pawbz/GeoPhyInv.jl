# 
# implementation greatly inspired from: https://github.com/geodynamics/seismic_cpml/blob/master

# # Credits 
# Author: Pawan Bharadwaj 
#         (bharadwaj.pawan@gmail.com)
# 
# * original code in FORTRAN90: March 2013
# * modified: 11 Sept 2013
# * major update: 25 July 2014
# * code optimization with help from Jan Thorbecke: Dec 2015
# * rewritten in Julia: June 2017
# * added parrallelization over supersources in Julia: July 2017
# * efficient parrallelization using distributed arrays: Sept 2017
# * optimized memory allocation: Oct 2017

# 

# global const nlayer_rand = 0

include("types.jl")
include("attenuation.jl")

# 
#As forward modeling method, the 
#finite-difference method is employed. 
#It uses a discrete version of the two-dimensional isotropic acoustic wave equation.
#
#```math
#\pp[\tzero] - \pp[\tmo] = \dt \mB \left({\partial_x\vx}[\tmh]
# + \partial_z \vz[\tmh]  + \dt\sum_{0}^{\tmo}\sfo\right)
# ```
# ```math
#\pp[\tpo] - \pp[\tzero] = \dt \mB \left(\partial_x \vx[\tph]
# + {\partial_z \vz}[\tph]  + \dt\sum_{0}^{\tzero}\sfo\right)
# ```

#attenuation related
#Rmemory::StackedArray2DVector{Array{Float64,3}} # memory variable for attenuation, see Carcione et al (1988), GJ
#Rmemoryp::StackedArray2DVector{Array{Float64,3}} # at previous time step

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

* `sflags=2` : source related flags 
  * `=0` inactive sources
  * `=1` sources with injection rate
  * `=2` volume injection sources
  * `=3` sources input after time reversal (use only during backpropagation)
* `rflags=1` : receiver related flags 
  * `=0` receivers do not record (or) inactive receivers
  * `=1` receivers are active only for the second propagating wavefield
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
    args2...,
) = PFdtd(attrib_mod, args1...; args2...)


"""
Primary method to generate Expt variable when FdtdAcoustic() and FdtdAcoustic{Born}() are used.

# Some internal arguments
* `jobname::String` : name
* `npw::Int64=1` : number of independently propagating wavefields in `medium`
* `backprop_flag::Bool=Int64` : save final state variables and the boundary conditions for later use
  * `=1` save boundary and final values in `boundary` 
  * `=-1` use stored values in `boundary` for back propagation
* `gmodel_flag=false` : flag that is used to  output gradient; there should be atleast two propagating wavefields in order to do so: 1) forward wavefield and 2) adjoint wavefield
* `illum_flag::Bool=false` : flag to output wavefield energy or source illumination; it can be used as preconditioner during inversion
"""
function PFdtd(
    attrib_mod;
    jobname::Symbol = :forward_propagation,
    npw::Int64 = 1,
    medium::Medium = nothing,
    pml_faces::Vector{Symbol} = [:zmin, :zmax, :ymin, :ymax, :xmax, :xmin],
    rigid_faces = pml_faces,
    tgrid::StepRangeLen = nothing,
    ageom::Union{AGeom,Vector{AGeom}} = nothing,
    srcwav::Union{SrcWav,Vector{SrcWav}} = nothing,
    sflags::Int = 2,
    rflags::Int = 1,
    rfields::Vector{Symbol} = [:vz],
    stressfree_faces = [:dummy],
    backprop_flag::Int64 = 0,
    gmodel_flag::Bool = false,
    illum_flag::Bool = false,
    tsnaps::AbstractVector{Float64} = fill(0.5 * (tgrid[end] + tgrid[1]), 1),
    snaps_field = nothing,
    verbose::Bool = false,
    nworker = nothing,
)
    N = ndims(medium)


    if (first(typeof(attrib_mod).parameters) == Born)
        # born modeling needs npw=2
        npw = 2
        sflags=[sflags, 0] # source is only active in the background medium 
        rflags = [0, 1] # record scatterred data only 
    else
    (typeof(sflags) == Int) &&        (sflags = fill(sflags, npw))
    (typeof(rflags) == Int) &&     (rflags = fill(rflags, npw))

    end

    
    # convert to vectors of length npw
    (typeof(ageom) == AGeom) && (ageom = fill(ageom, npw))
    (typeof(srcwav) == SrcWav) &&  (srcwav = fill(srcwav, npw))
    

    # check if fields input are meaningful for given attrib_mod and ndims
    @assert (N == _fd.ndims) "Cannot initiate SeisForwExpt due to ndims inconsistency with @init_parallel_stencil"
    !(snaps_field === nothing) && @assert snaps_field ∈ Fields(attrib_mod, ndims = N)
    foreach(rf -> @assert(rf ∈ Fields(attrib_mod; ndims = N)), rfields)
    foreach(
        s -> foreach(
            ss -> foreach(f -> @assert(f ∈ Fields(attrib_mod, ndims = N)), names(ss.d)[1]),
            s,
        ),
        srcwav,
    )

    # check sizes and errors based on input
    #(length(TDout) ≠ length(findn(rflags))) && error("TDout dimension")
    @assert (length(ageom) == npw) 
    @assert (length(srcwav) == npw) 
    @assert (length(sflags) == npw) 
    @assert (length(rflags) == npw) 
    @assert (maximum(tgrid) >= maximum(srcwav[1][1].grid)) "modeling time is less than source time"
    #(any([getfield(TDout[ip],:tgrid).δx < tgridmod.δx for ip=1:length(TDout)])) && error("output time grid sampling finer than modeling")
    #any([maximum(getfield(TDout[ip],:tgrid).x) > maximum(tgridmod) for ip=1:length(TDout)]) && error("output time > modeling time")

    # all the propagating wavefields should have same supersources, check that
    nss = length(ageom[1])
    (fill(nss, npw) != [length(ageom[ip]) for ip = 1:npw]) && error("different supersources")

    # check if all sources are receivers are inside medium
    any(.![(ageom[ip] ∈ medium.mgrid) for ip = 1:npw]) && error("sources or receivers not inside medium")

    all([issimilar(ageom[ip], srcwav[ip]) for ip = 1:npw]) || error("ageom and srcwav mismatch")

    #(verbose) &&	println(string("\t> number of super sources:\t",nss))	

    # extend mediums in the PML layers
    exmedium = padarray(medium, _fd.npextend, pml_faces)
    mod = NamedArray(
        [
            Data.Array(zeros(length.(exmedium.mgrid)...)) for
            name in get_medium_names(attrib_mod)
        ],
        get_medium_names(attrib_mod),
    )

    # PML
    pml = get_pml(attrib_mod, exmedium.mgrid)

    # initialize gradient arrays
    gradient = zeros(2 * prod(length.(medium.mgrid)))

    # shared arrays required to reduce all the gradient from individual workers
    grad_mod = NamedArray(
        [
            SharedArray{Float64}(zeros(length.(exmedium.mgrid)...)) for
            name in get_medium_names(attrib_mod)
        ],
        get_medium_names(attrib_mod),
    )

    illum_stack = SharedArray{Float64}(zeros(length.(medium.mgrid)...))

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
    datamat = SharedArray{Float64}(length(tgrid), maximum(nrmat), nss)
    data = [Records(tgrid, ageom[ip], rfields) for ip = 1:length(findall(!iszero, rflags))]

    # default is all prpagating wavefields are active
    activepw = [ipw for ipw = 1:npw]

    # dont need visco parameters, initialize with proper sizes later
    nsls = Int32(0)
    memcoeff1 = zeros(1, 1, 1)
    memcoeff2 = zeros(1, 1, 1)

    fc = get_fc(exmedium, tgrid)
    ic = get_ic(exmedium, tgrid, nsls, npw)

    pac = P_common(
        jobname,
        attrib_mod,
        activepw,
        exmedium,
        medium,
        ageom,
        srcwav,
        pml_faces,
        # pml_faces also need rigid_boundary conditions
        unique(vcat(rigid_faces, pml_faces)),
        stressfree_faces,
        sflags,
        rfields,
        rflags,
        fc,
        ic,
        pml,
        mod,
        NamedArray([memcoeff1, memcoeff2], ([:memcoeff1, :memcoeff2],)),
        gradient,
        grad_mod,
        illum_flag,
        illum_stack,
        backprop_flag,
        snaps_field,
        itsnaps,
        gmodel_flag,
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
    ssi = [round(Int, s) for s in range(0, stop = nss, length = nworker + 1)]
    sschunks = Array{UnitRange{Int64}}(undef, nworker)
    for ib = 1:nworker
        sschunks[ib] = ssi[ib]+1:ssi[ib+1]
    end

    # a distributed array of P_x_worker --- note that the parameters for each super source are efficiently distributed here
    papa = ddata(
        T = _fd.use_gpu ?
            Vector{P_x_worker_x_pw{N,CUDA.CuArray{_fd.datatype,N,CUDA.Mem.DeviceBuffer}}} : Vector{P_x_worker_x_pw{N,Array{_fd.datatype,N}}},
        init = I -> Vector{P_x_worker_x_pw}(sschunks[I...][1], pac),
        pids = work[1:1], # disable distributed in order for Pluto to work (waiting for Pluto bug to enable distributed)
    )

    pa = PFdtd(sschunks, papa, pac)

    # update source wavelets
    update!(pa, pa.c.srcwav, pa.c.sflags, verbose = true)

    # update pml coefficient, now that freqpeak is scanned from the source wavelets
    update_pml!(pa.c)

    check_stability(pa, verbose, H = div(20, _fd.order))


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
        vcat(length.(medium.mgrid), [length(tgrid), nsls, npw]),
        vcat(dim_names(N, "n"), [:nt, :nsls, :npw]),
    )
end

"""
Get some float constants, and store then in pac, e.g, spatial sampling.
The idea is to use them later inside the loops for faster modelling.
"""
function get_fc(medium, tgrid)
    N = ndims(medium)
    ds = step.(medium.mgrid)
    # denominator depending on _fd.order
    dsI = inv.(ds)
    if (_fd.order == 4)
        dsI = inv.(ds .* 24.0)
    elseif (_fd.order == 6)
        dsI = inv.(ds .* 1920.0)
    elseif (_fd.order == 8)
        dsI = inv.(ds .* 107520.0)
    end
    dt = step(tgrid)
    dtI = inv(dt)

    return NamedArray(
        vcat([dt, dtI], ds, dsI),
        vcat([:dt, :dtI], dim_names(N, "d"), dim_names(N, "d", "I")),
    )
end

function get_medium_names(::FdtdAcoustic{FWmod})
    # the following medium parameters are stored on XPU
    return [:K, :rhoI]
end
function get_medium_names(::FdtdAcoustic{Born})
    # the following medium parameters are stored on XPU
    return [:K, :rhoI, :δK, :δrhoI]
end

function get_medium_names(::FdtdElastic)
    # the following medium parameters are stored on XPU
    return [:lambda, :M, :mu, :rho]
end


"""
Create field arrays for each worker.
Each worker performs the modeling of supersources in `sschunks`.
The parameters common to all workers are stored in `pac`.
"""
function P_x_worker_x_pw(ipw, sschunks::UnitRange{Int64}, pac::P_common{T,N}) where {T,N}
    n = length.(pac.exmedium.mgrid)

    born_svalue_stack = zeros(n...)

    fields = Fields(pac.attrib_mod)
    w1 = NamedArray(
        [NamedArray([zeros(eval(f)(), T(), n...) for f in fields], Symbol.(fields))],
        [:t],
    )

    # temp field arrays copied on the CPU before performing scalar operations
    wr = NamedArray(
        [
            _fd.use_gpu ? Array(zeros(eval(f)(), T(), n...)) :
            zeros(Data.Number, fill(1, N)...) for f in pac.rfields
        ],
        Symbol.(pac.rfields),
    )

    # dummy (use for viscoelastic modeling later)
    w2 = NamedArray(
        [NamedArray([@zeros(fill(1, N)...) for i in [:r]], ([:r],)) for i = 1:2],
        ([:t, :tp],),
    )

    # memory fields for all derivatives
    dfields = Fields(pac.attrib_mod, "d")
    memory_pml = NamedArray(
        [zeros(eval(f)(), T(), n..., pml = true) for f in dfields],
        Symbol.(dfields),
    )

    ss = [P_x_worker_x_pw_x_ss(ipw, iss, pac) for (issp, iss) in enumerate(sschunks)]

    return P_x_worker_x_pw(ss, w1, wr, w2, memory_pml, born_svalue_stack)
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
    sfields = names(pac.srcwav[ipw][iss].d)[1]
    nt = pac.ic[:nt]
    n = length.(pac.exmedium.mgrid)
    ageom = pac.ageom
    srcwav = pac.srcwav
    sflags = pac.sflags

    # records_output, distributed array among different procs
    records = NamedArray(
        [[Data.Array(zeros(ageom[ipw][iss].nr)) for it = 1:nt] for i = 1:length(rfields)],
        (rfields,),
    )

    # gradient outputs
    grad_modlambda = zeros(1, 1)


    # saving illum
    # illum =  (pac.illum_flag) ? zeros(nz, nx) : zeros(1,1)
    illum = zeros(1, 1)

    snaps = NamedArray(
        [
            (eval(pac.snaps_field) <: Fields) ?
            Array(zeros(eval(pac.snaps_field)(), T(), n...)) :
            zeros(Data.Number, fill(1, N)...) for i = 1:length(pac.itsnaps)
        ],
        names(pac.itsnaps)[1],
    )
    # source wavelets
    wavelets = [
        NamedArray(
            [Data.Array(zeros(ageom[ipw][iss].ns)) for i = 1:length(sfields)],
            (sfields,),
        ) for it = 1:nt
    ]

    # boundary fields stored (remove derivatives)
    # fields are stored for backpropagation (used for RTM and FWI)
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
                prod(length.(get_mgrid(eval(sf)(), T(), pac.exmedium.mgrid...))),
                ageom[ipw][iss].ns,
            ) for sf in sfields
        ],
        (sfields,),
    )
    rinterpolatew = NamedArray(
        [
            spzeros(
                Data.Number,
                ageom[ipw][iss].nr,
                prod(length.(get_mgrid(eval(rf)(), T(), pac.exmedium.mgrid...))),
            ) for rf in rfields
        ],
        rfields,
    )

    # transfer to GPUs if 
    if (_fd.use_gpu)
        ssprayw = map(x -> CuSparseMatrixCSC(x), ssprayw)
        rinterpolatew = map(x -> CuSparseMatrixCSC(x), rinterpolatew)
    end


    pass = P_x_worker_x_pw_x_ss(
        iss,
        wavelets,
        ssprayw,
        records,
        rinterpolatew,
        boundary,
        snaps,
        illum,# grad_modKI,grad_modrhovxI,grad_modrhovzI)
        NamedArray([grad_modlambda], ([:lambda],)),
    )

    # update acquisition
    update!(pass, ipw, iss, ageom[ipw][iss], pac)
    return pass
end

include("proj_mat.jl")
include("source.jl")
include("receiver.jl")
include("dirichlet.jl")
include("cpml.jl")
include("advance_acou.jl")
include("advance_elastic.jl")
include("rho_projection.jl")
include("gradient.jl")
include("born_acoustic.jl")
include("boundary.jl")
include("gallery.jl")


# update TDout after forming a vector and resampling
#	ipropout=0;
#	for iprop in 1:pac.ic[:npw]
#		if(pac.rflags[iprop] ≠ 0)
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
        gss = view(gs, _fd.npml+1:nz-_fd.npml, _fd.npml+1:nx-_fd.npml)
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
include("updates.jl")
include("getprop.jl")

