# Methods to update FDTD structs in several ways

#=
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
    δmodKI =
        view(pac.δmod[:KI], _fd.npextend+1:nz-_fd.npextend, _fd.npextend+1:nx-_fd.npextend)
    for i = 1:nznxd
        # put perturbation due to KI as it is
        δmodKI[i] = δmod[i]
    end
    fill!(pac.δmod[:rhoI], 0.0)
    δmodrr = view(
        pac.δmod[:rhoI],
        _fd.npextend+1:nz-_fd.npextend,
        _fd.npextend+1:nx-_fd.npextend,
    )
    for i = 1:nznxd
        # put perturbation due to rhoI here
        δmodrr[i] = δmod[nznxd+i]
    end
    # project δmodrr onto the vz and vx grids
    get_rhovxI!(pac.δmod[:rhovxI], pac.δmod[:rhoI])
    get_rhovzI!(pac.δmod[:rhovzI], pac.δmod[:rhoI])
    return nothing
end
=#

function get_update_parameters(pa::T) where {T<:Union{PFdtd{<:FdtdAcoustic{FullWave}},PFdtd{<:FdtdElastic{FullWave}}}}
    if (pa.c.ic[:npw] == 1)
        return (; activepw=[1], src_flags=[true], rec_flags=[true])
    elseif (pa.c.ic[:npw] == 2)
        if (pa.c.attrib_mod.mode == :adjoint)
            return (; activepw=[1, 2], src_flags=[true, true], rec_flags=[false, false])
        elseif (pa.c.attrib_mod.mode == :forward)
            return (; activepw=[1], src_flags=[true, false], rec_flags=[true, false])
        else
            error("unknown mode for $(pa.c.attrib_mod)")
        end
    end
end


function get_update_parameters(pa::T) where {T<:Union{PFdtd{<:FdtdAcoustic{Born}},PFdtd{<:FdtdElastic{Born}}}}
    @assert pa.c.ic[:npw] == 2
    if (pa.c.attrib_mod.mode == :adjoint)
        return (; activepw=[1, 2], src_flags=[true, true], rec_flags=[false, false])
    else
        return (; activepw=[1, 2], src_flags=[true, false], rec_flags=[false, true])
    end
end

"""
```julia
update!(pa)
```
In-place method to perform the experiment and update `pa` after wave propagation. After update, see
`pa[:data]` and `pa[:snaps]`.
Source added when src_flags is true, switch off the src_flag whenever necessary, useful when npw=2  

* `rec_flags=1` : Bool array to specify which wavefield is recorded 
  * `[true, false]` records the first pw
  * `[false, true]` records the second wavefield, used during Born
"""
@fastmath function update!(pa::PFdtd; upa=get_update_parameters(pa))
    (; activepw, src_flags, rec_flags) = upa

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
                mod_x_proc!(pa.c, localpart(pa.p), activepw, src_flags)
            end
        end
    end
    # end

    # stack gradients and illum over sources
    @sync begin
        for (ip, p) in enumerate(procs(pa.p))
            @sync remotecall_wait(p) do
                sum_grads!(Val(pa.c.attrib_mod.mode), pa.c, localpart(pa.p))
                # (pa.c.illum_flag) && stack_illums!(pa.c, localpart(pa.p))
            end
        end
    end

    for ipw in activepw
        if (rec_flags[ipw]) # record only if rec_flags is non-zero
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
    return to
end

# modelling for each worker
function mod_x_proc!(pac::P_common, pap::Vector{P_x_worker_x_pw{N,B}}, activepw, src_flags) where {N,B}

    # source_loop
    for issp = 1:length(pap[1].ss) # note, all npw have same sources
        reset_w2!(pap)

        iss = pap[1].ss[issp].iss # same note as above

        # only for first propagating wavefield, i.e., pap[1]
        "initial conditions from boundary for first propagating field only"
        boundary_force_snap_tau!(Val(pac.attrib_mod.mode), issp, pac, pap[1])
        boundary_force_snap_v!(Val(pac.attrib_mod.mode), issp, pac, pap[1])

        prog = Progress(
            pac.ic[:nt],
            dt=1,
            desc="$(pac.attrib_mod) supershot $iss/$(length(pac.ageom[1])) ",
            showspeed=true,
            color=:blue,
        )
        # time_loop
        #=
        don't use shared arrays inside this time loop, for speed when using multiple procs
        =#
        for it = 1:pac.ic[:nt]

            pac.verbose && next!(prog, :blue)

            # force p[1] on faces, only for ipw=1, after time reversal 
            boundary_force!(Val(pac.attrib_mod.mode), pac.ic[:nt] - it + 1, issp, pac, pap[1])

            @timeit to "compute stress derivatives" begin
                for ipw in activepw
                    update_dstress!(pap[ipw], pac)
                end
            end
            @timeit to "compute particle velocities" begin
                for ipw in activepw
                    update_v!(pap[ipw], pac)
                end
            end
            add_born_sources_velocity!(pap, pac)

            @timeit to "compute velocity derivatives" begin
                for ipw in activepw
                    update_dv!(pap[ipw], pac)
                end
            end
            @timeit to "compute stress" begin
                for ipw in activepw
                    update_stress!(pap[ipw], pac)
                end
            end
            add_born_sources_stress!(pap, pac)

            @timeit to "add source" begin
                add_source!(it, issp, iss, pac, pap, activepw, src_flags)
            end

            # record wavefield on all the faces for ipw=1
            boundary_save!(Val(pac.attrib_mod.mode), it, issp, pac, pap[1])


            @timeit to "record at receivers" begin
                record!(it, issp, iss, pac, pap, activepw)
            end

            @timeit to "compute gradient" begin
                compute_gradient!(Val(pac.attrib_mod.mode), issp, pac, pap)
            end

            # (pac.illum_flag) && compute_illum!(issp, pap)

            if (eval(pac.snaps_field) <: Fields)
                iitsnaps = findall(x -> ==(x, it), pac.itsnaps)
                for ipw in activepw
                    for itsnap in iitsnaps
                        snaps_save!(itsnap, issp, ipw, pac, pap)
                    end
                end
            end

        end # time_loop
        # now pressure is at [nt], velocities are at [nt-1/2]

        # save last snap of pressure field for solving initial-value problem later
        boundary_save_snap_tau!(Val(pac.attrib_mod.mode), issp, pac, pap[1])

        # one more advance step to move pressure to [nt+1] and 
        # velocities to [nt+3/2] for solving initial-value problem
        update_dstress!(pap[1], pac)
        update_v!(pap[1], pac)

        # save last snap of velocity fields
        boundary_save_snap_v!(Val(pac.attrib_mod.mode), issp, pac, pap[1])

        # scale gradients for each issp
        # (pac.gmodel_flag) && scale_gradient!(
        # issp,
        # pap,
        # step(pac.medium.mgrid[2]) * step(pac.medium.mgrid[1]),
        # )

    end # source_loop
end # mod_x_shot


