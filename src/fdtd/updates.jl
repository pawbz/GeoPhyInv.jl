# Methods to update FDTD structs in several ways


"""
update the perturbation vector using the perturbed medium
in this case, medium will be treated as the background medium 
* `δmod` is [δKI, δrhoI]
"""
function update_δmods!(pac::P_common, δmod::Vector{Float64})
    nx = pac.ic[:nx]
    nz = pac.ic[:nz]
    nznxd = prod(length.(pac.medium.mgrid))
    copyto!(pac.δmodall, δmod)
    fill!(pac.δmod[:KI], 0.0)
    δmodKI = view(pac.δmod[:KI], npml+1:nz-npml, npml+1:nx-npml)
    for i = 1:nznxd
        # put perturbation due to KI as it is
        δmodKI[i] = δmod[i]
    end
    fill!(pac.δmod[:rhoI], 0.0)
    δmodrr = view(pac.δmod[:rhoI], npml+1:nz-npml, npml+1:nx-npml)
    for i = 1:nznxd
        # put perturbation due to rhoI here
        δmodrr[i] = δmod[nznxd+i]
    end
    # project δmodrr onto the vz and vx grids
    get_rhovxI!(pac.δmod[:rhovxI], pac.δmod[:rhoI])
    get_rhovzI!(pac.δmod[:rhovzI], pac.δmod[:rhoI])
    return nothing
end

"""
This method should be executed only after the updating the main medium.
Update the `δmods` when a perturbed `medium_pert` is input.
The medium through which the waves are propagating 
is assumed to be the background medium.
"""
function update_δmods!(pac::P_common, medium_pert::Medium)
    nznxd = prod(length.(pac.medium.mgrid))
    fill!(pac.δmodall, 0.0)
    copyto!(pac.δmodall, medium_pert, [:KI, :rhoI])

    for i = 1:nznxd
        pac.δmodall[i] -= pac.mod[:KI][i] # subtracting the background medium
        pac.δmodall[nznxd+i] -= pac.mod[:rhoI][i] # subtracting the background medium
    end
    update_δmods!(pac, pac.δmodall)
    return nothing
end

"""
```julia
update!(pa,medium_new)
```
Update `pa` with a new bundle of medium parameters `medium_new`, without additional memory allocation.
This routine is used during inversion, where medium parameters are iteratively updated. 
The ability to iteratively run the forward mediuming task (with no additional memory allocation) on  
various subsurface mediums is necessary while implementing inversion 
algorithms.
"""
function update!(pa::PFdtd, medium::Medium)
    return update!(pa.c, medium)
end
function update!(pac::T, medium::Medium) where {T<:P_common{FdtdAcou}}
    copyto!(pac.medium, medium)
    padarray!(pac.exmedium, pac.medium, npml, pac.pml_edges)
    copyto!(pac.mod[:K], pac.exmedium, [:K])
    copyto!(pac.mod[:KI], pac.exmedium, [:KI])
    copyto!(pac.mod[:rhoI], pac.exmedium, [:rhoI])
    # get_rhovxI!(pac.mod[:rhovxI], pac.mod[:rhoI])
    # get_rhovzI!(pac.mod[:rhovzI], pac.mod[:rhoI])
    return nothing
end
function update!(pac::T, medium::Medium) where {T<:P_common{FdtdElastic}}
    copyto!(pac.medium, medium)
    padarray!(pac.exmedium, pac.medium, npml, pac.pml_edges)
    copyto!(pac.mod[:mu], pac.exmedium, [:mu])
    copyto!(pac.mod[:lambda], pac.exmedium, [:lambda])
    copyto!(pac.mod[:M], pac.exmedium, [:M])
    copyto!(pac.mod[:rho], pac.exmedium, [:rho])
    return nothing
end


"""
```julia
update!(pa,srcwav_new,sflags)
```
Update `pa` with a new bundle of source wavelets `srcwav_new`, without additional memory allocation.
Optionally, `sflags` can be changed. 
"""
function update!(pa::PFdtd, srcwav::SrcWav, sflags::Any = nothing; verbose = false)
    update!(pa, [srcwav], sflags)
end
function update!(pa::PFdtd, srcwav::Vector{SrcWav}, sflags = nothing; verbose = false)
    # update srcwav in pa.c
    (length(srcwav) ≠ pa.c.ic[:npw]) && error("cannot update")
    for i = 1:length(srcwav)
        copyto!(pa.c.srcwav[i], srcwav[i])
    end
    if (!(sflags === nothing))
        copyto!(pa.c.sflags, sflags)
    end
    for ipw = 1:pa.c.ic[:npw]
        freqmin = 0.0
        freqmax = Inf
        freqpeaks = []
        # fill_wavelets for each supersource
        @sync begin
            for (ip, p) in enumerate(procs(pa.p))
                @async remotecall_wait(p) do
                    pap = localpart(pa.p)[ipw]
                    for is = 1:length(pap.ss)
                        iss = pap.ss[is].iss
                        wavelets = pap.ss[is].wavelets
                        broadcast(x -> fill!.(x, 0.0), wavelets)
                        fmin, fmax, fpeak =
                            fill_wavelets!(ipw, iss, wavelets, pa.c.srcwav, pa.c.sflags)
                        freqmin = max(freqmin, fmin)
                        freqmax = min(freqmax, fmax)
                        push!(freqpeaks, fpeak)
                    end
                end
            end
        end
        freqpeak = Statistics.mean(freqpeaks)
        # store frequency (in Hz) bounds for first propagating wavefield
        if (ipw == 1)
            for (f, nf) in
                zip([freqmin, freqmax, freqpeak], [:freqmin, :freqmax, :freqpeak])
                if (nf ∈ names(pa.c.fc)[1])
                    pa.c.fc[nf] = f
                else
                    pa.c.fc = vcat(pa.c.fc, NamedArray([f], [nf]))
                end
            end
        end
        if (verbose)
            freqmin = @sprintf("%0.2e", freqmin)
            freqmax = @sprintf("%0.2e", freqmax)
            freqpeak = @sprintf("%0.2e", freqpeak)
            @info "frequency bounds for propagating wavefield $ipw are: [$freqmin, $freqmax], with peak at: $freqpeak"
        end
    end
end

function update!(pass::P_x_worker_x_pw_x_ss, ipw, iss, ageomss::AGeomss, pac)
    @assert ageomss.ns == pac.ageom[ipw][iss].ns
    @assert ageomss.nr == pac.ageom[ipw][iss].nr

    ssprayw = pass.ssprayw
    rinterpolatew = pass.rinterpolatew
    
    for sfield in names(pass.sindices)[1]
        sindices=pass.sindices[sfield]
    for is = 1:ageomss.ns
        w = ssprayw[sfield][is]
        sindices[is], ww = get_indices_weights(eval(sfield)(),
            pac.exmedium.mgrid...,
            [s[is] for s in ageomss.s],
        )
        copyto!(w, ww)
    end
end
    for rfield in names(pass.rindices)[1]
    rindices = pass.rindices[rfield]
    for ir = 1:ageomss.nr
        w = rinterpolatew[rfield][ir]
        rindices[ir], ww = get_indices_weights(eval(rfield)(),
            pac.exmedium.mgrid...,
            [r[ir] for r in ageomss.r],
        )
        copyto!(w, ww)
    end
end

end

# if just one propagating field
update!(pa::PFdtd, ageom::AGeom) = update!(pa, [ageom])

function update!(pa::PFdtd, ageom::Vector{AGeom})
    for ipw = 1:pa.c.ic[:npw]
        copyto!(pa.c.ageom[ipw], ageom[ipw])
        @sync begin
            for (ip, p) in enumerate(procs(pa.p))
                @async remotecall_wait(p) do
                    pap = localpart(pa.p)[ipw]
                    for is = 1:length(pap.ss)
                        iss = pap.ss[is].iss
                        update!(localpart(pap).ss[is], ipw, iss, ageom[ipw][iss], pa.c)
                    end
                end
            end
        end
    end
end


"""
```julia
update!(pa)
```
In-place method to perform the experiment and update `pa` after wave propagation. After update, see
`pa[:data]` and `pa[:snaps]`.
"""
@fastmath function update!(pa::PFdtd)

    global to

    reset_timer!(to)

    # zero out all the results stored in pa.c
    initialize!(pa.c)

    # zero out results stored per worker
    @sync begin
        for (ip, p) in enumerate(procs(pa.p))
            @async remotecall_wait(p) do
                initialize!(localpart(pa.p))
            end
        end
    end


    # @timeit to "mod_x_proc!" begin
        # all localparts of DArray are input to this method
        # parallelization over shots
    @sync begin
        for (ip, p) in enumerate(procs(pa.p))
            @async remotecall_wait(p) do
                mod_x_proc!(pa.c, localpart(pa.p))
            end
        end
    end
    # end

    # stack gradients and illum over sources
    @sync begin
        for (ip, p) in enumerate(procs(pa.p))
            @sync remotecall_wait(p) do
                (pa.c.gmodel_flag) && stack_grads!(pa.c, localpart(pa.p))
                (pa.c.illum_flag) && stack_illums!(pa.c, localpart(pa.p))
            end
        end
    end

    # update gradient medium using grad_modKI_stack, grad_modrr_stack
    (pa.c.gmodel_flag) && update_gradient!(pa.c)


    for ipw in pa.c.activepw
        if (pa.c.rflags[ipw] ≠ 0) # record only if rflags is non-zero
            for rfield in pa.c.rfields

                fill!(pa.c.datamat, 0.0)
                @sync begin
                    for (ip, p) in enumerate(procs(pa.p))
                        @sync remotecall_wait(p) do
                            update_datamat!(rfield, ipw, pa.c, localpart(pa.p))
                        end
                    end
                end
                update_data!(rfield, ipw, pa.c)
            end
        end
    end
    if (pa.c.verbose)
        println("timer output of methods inside fdtd time loop")
        show(to)
        println("  ")
    end
    return nothing
end

# modelling for each worker
function mod_x_proc!(pac::P_common, pap::Vector{P_x_worker_x_pw{N,B}}) where {N,B}

    # source_loop
    for issp = 1:length(pap[1].ss) # note, all npw have same sources
        reset_w2!(pap)

        iss = pap[1].ss[issp].iss # same note as above

        # only for first propagating wavefield, i.e., pap[1]
        if (pac.backprop_flag == -1)
            "initial conditions from boundary for first propagating field only"
            boundary_force_snap_p!(issp, pac, pap)
            boundary_force_snap_vxvz!(issp, pac, pap)
        end

        prog = Progress(
            pac.ic[:nt],
            dt = 1,
            desc = "\tmodeling supershot $iss/$(length(pac.ageom[1])) ",
            color = :blue,
        )
        # time_loop
        """
        * don't use shared arrays inside this time loop, for speed when using multiple procs
        """
        for it = 1:pac.ic[:nt]

            pac.verbose && next!(prog, :white)

            @timeit to "advance time step" begin
                advance!(pac, pap)
            end

            # force p[1] on boundaries, only for ipw=1
            (pac.backprop_flag == -1) && boundary_force!(it, issp, pac, pap)

            @timeit to "add source" begin
                add_source!(it, issp, iss, pac, pap)
            end

            # no born flag for adjoint modelling
            if (!pac.gmodel_flag)
                (typeof(pac.attrib_mod) == FdtdAcouBorn) &&
                    add_born_sources!(issp, pac, pap)
            end

            # record boundaries after time reversal already
            (pac.backprop_flag == 1) && boundary_save!(pac.ic[:nt] - it + 1, issp, pac, pap)

            @timeit to "record at receivers" begin
                record!(it, issp, iss, pac, pap)
            end

            if(pac.gmodel_flag)
                @timeit to "compute gradient" begin
                    compute_gradient!(issp, pac, pap)
                end
            end

            (pac.illum_flag) && compute_illum!(issp, pap)

            if (eval(pac.snaps_field) <: Fields)
                iitsnaps = findall(x -> ==(x, it), pac.itsnaps)
                for itsnap in iitsnaps
                    snaps_save!(itsnap, issp, pac, pap)
                end
            end

        end # time_loop
        "now pressure is at [nt], velocities are at [nt-1/2]"

        "one more propagating step to save pressure at [nt+1] -- for time revarsal"
        advance!(pac, pap)

        "save last snap of pressure field"
        (pac.backprop_flag == 1) && boundary_save_snap_p!(issp, pac, pap)

        "one more propagating step to save velocities at [nt+3/2] -- for time reversal"
        advance!(pac, pap)

        "save last snap of velocity fields with opposite sign for adjoint propagation"
        (pac.backprop_flag == 1) && boundary_save_snap_vxvz!(issp, pac, pap)

        "scale gradients for each issp"
        (pac.gmodel_flag) && scale_gradient!(
            issp,
            pap,
            step(pac.medium.mgrid[2]) * step(pac.medium.mgrid[1]),
        )


    end # source_loop
end # mod_x_shot




function update_datamat!(
    rfield,
    ipw,
    pac::P_common,
    pap::Vector{P_x_worker_x_pw{N,B}},
) where {N,B}
    datamat = pac.datamat
    pass = pap[ipw].ss
    for issp = 1:length(pass)
        iss = pass[issp].iss
        records = pass[issp].records
        for ir = 1:pac.ageom[ipw][iss].nr
            for it = 1:pac.ic[:nt]
                datamat[it, ir, iss] = records[rfield][it, ir]
            end
        end
    end
end

function update_data!(rfield, ipw, pac::P_common)
    datamat = pac.datamat
    for iss = 1:length(pac.ageom[1])
        data = pac.data[ipw][iss].d[rfield]
        for ir = 1:pac.ageom[ipw][iss].nr
            for it = 1:pac.ic[:nt]
                data[it, ir] = datamat[it, ir, iss]
            end
        end
    end
end

