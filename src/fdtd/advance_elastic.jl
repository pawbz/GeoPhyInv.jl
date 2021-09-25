function advance!(pac::T, pap) where {T<:P_common{<:FdtdElastic}}
    for ipw in pac.activepw
        # store p for the last two steps
        # pppppp!(pap[ipw],pac.attrib_mod)
        advance_kernel!(pap[ipw], pac)
    end
    return nothing
end


# these relative indices of the arrays point to same location
# [ix,iy,iz]     --> tauxx, tauyy, and tauzz grid [1:nx,1:ny,1:nz]
# [ix+1/2,iy+1/2,iz+1]      --> tauxy
# [ix+1/2,iy+1,iz+1/2]      -->  tauxz
# [ix-1/2,iy,iz]       --> vx
# [ix,iy-1/2,iz]       --> vy
# [ix,iy,iz-1/2]       --> vz

@parallel function compute_dtau!(
    tauxx::Data.Array{3},
    tauyy,
    tauzz,
    tauxy,
    tauxz,
    tauyz,
    dtauxxdx,
    dtauxydx,
    dtauxzdx,
    dtauyydy,
    dtauxydy,
    dtauyzdy,
    dtauzzdz,
    dtauyzdz,
    dtauxzdz,
    dxI,
    dyI,
    dzI,
)

    @all(dtauxxdx) = @d_zi(tauxx) * dxI # at [ix+1/2,iy+1,iz+1] with indices []
    @all(dtauxydx) = @d_za(tauxy) * dxI # at [ix+1,iy+1/2,iz+1]
    @all(dtauxzdx) = @d_za(tauxz) * dxI # at [ix+1,iy+1,iz+1/2]


    @all(dtauyydy) = @d_yi(tauyy) * dyI # at [ix+1,iy+1/2,iz+1]
    @all(dtauxydy) = @d_ya(tauxy) * dyI # at [ix+1/2,iy+1,iz+1]
    @all(dtauyzdy) = @d_ya(tauyz) * dyI # at [ix+1,iy+1,iz+1/2]


    @all(dtauzzdz) = @d_xi(tauzz) * dzI # at [ix+1,iy+1,iz+1/2]
    @all(dtauxzdz) = @d_xa(tauxz) * dzI # at [ix+1/2,iy+1,iz+1]
    @all(dtauyzdz) = @d_xa(tauyz) * dzI # at [ix+1,iy+1/2,iz+1]

    return
end

@parallel function compute_dtau!(
    tauxx::Data.Array{2},
    tauzz,
    tauxz,
    dtauxxdx,
    dtauxzdx,
    dtauzzdz,
    dtauxzdz,
    dxI,
    dzI,
)

    @all(dtauxxdx) = @d_yi(tauxx) * dxI # at [ix+1/2,iz+1] with indices []
    @all(dtauxzdx) = @d_ya(tauxz) * dxI # at [ix+1,iz+1/2]

    @all(dtauzzdz) = @d_xi(tauzz) * dzI # at [ix+1,iz+1/2]
    @all(dtauxzdz) = @d_xa(tauxz) * dzI # at [ix+1/2,iz+1]

    return
end

@parallel function compute_v!(
    vx::Data.Array{3},
    vy,
    vz,
    dtauxxdx,
    dtauxydx,
    dtauxzdx,
    dtauyydy,
    dtauxydy,
    dtauyzdy,
    dtauzzdz,
    dtauyzdz,
    dtauxzdz,
    dt,
    rho,
)
    @inn(vx) =
        @inn(vx) - dt / @av_zi(rho) * (@all(dtauxxdx) + @all(dtauxydy) + @all(dtauxzdz))
    @inn(vy) =
        @inn(vy) - dt / @av_yi(rho) * (@all(dtauxydx) + @all(dtauyydy) + @all(dtauyzdz))
    @inn(vz) =
        @inn(vz) - dt / @av_xi(rho) * (@all(dtauxzdx) + @all(dtauyzdy) + @all(dtauzzdz))

    return
end
@parallel function compute_v!(
    vx::Data.Array{2},
    vz,
    dtauxxdx,
    dtauxzdx,
    dtauzzdz,
    dtauxzdz,
    dt,
    rho,
)
    @inn(vx) = @inn(vx) - dt / @av_yi(rho) * (@all(dtauxxdx) + @all(dtauxzdz))
    @inn(vz) = @inn(vz) - dt / @av_xi(rho) * (@all(dtauxzdx) + @all(dtauzzdz))

    return
end



@parallel function compute_dv!(
    vx::Data.Array{3},
    vy,
    vz,
    dvxdx,
    dvydy,
    dvzdz,
    dvxdy,
    dvxdz,
    dvydx,
    dvydz,
    dvzdx,
    dvzdy,
    dxI,
    dyI,
    dzI,
)
    @all(dvxdx) = @d_za(vx) * dxI # at [ix,iy,iz]
    @all(dvydy) = @d_ya(vy) * dyI # at      "
    @all(dvzdz) = @d_xa(vz) * dzI # at      "

    @all(dvxdy) = @d_yi(vx) * dyI # at [ix+1/2,iy+1/2,iz+1]
    @all(dvxdz) = @d_xi(vx) * dzI # at [ix+1/2,iy+1,iz+1/2]

    @all(dvydz) = @d_xi(vy) * dzI # at [ix+1,iy+1/2,iz+1/2]
    @all(dvydx) = @d_zi(vy) * dxI # at [ix+1/2,iy+1/2,iz+1]

    @all(dvzdx) = @d_zi(vz) * dxI # at [ix+1/2,iy+1,iz+1/2]
    @all(dvzdy) = @d_yi(vz) * dyI # at [ix+1,iy+1/2,iz+1/2]

    return
end


@parallel function compute_dv!(vx::Data.Array{2}, vz, dvxdx, dvzdz, dvxdz, dvzdx, dxI, dzI)
    @all(dvxdx) = @d_za(vx) * dxI # at [ix,iz]
    @all(dvzdz) = @d_xa(vz) * dzI # at      "

    @all(dvxdz) = @d_xi(vx) * dzI # at [ix+1/2,iz+1/2]
    @all(dvzdx) = @d_zi(vz) * dxI # at [ix+1/2,iy+1,iz+1/2]

    return
end





@parallel function compute_tauii!(
    tauxx::Data.Array{3},
    tauyy,
    tauzz,
    dvxdx,
    dvydy,
    dvzdz,
    dt,
    M,
    lambda,
)

    @all(tauxx) =
        @all(tauxx) -
        dt * ((@all(M) * @all(dvxdx)) + (@all(lambda) * (@all(dvydy) + @all(dvzdz))))
    @all(tauyy) =
        @all(tauyy) -
        dt * ((@all(M) * @all(dvydy)) + (@all(lambda) * (@all(dvxdx) + @all(dvydy))))
    @all(tauzz) =
        @all(tauzz) -
        dt * ((@all(M) * @all(dvzdz)) + (@all(lambda) * (@all(dvydy) + @all(dvxdx))))

    return
end

@parallel function compute_tauii!(tauxx::Data.Array{2}, tauzz, dvxdx, dvzdz, dt, M, lambda)
    @all(tauxx) =
        @all(tauxx) - dt * ((@all(M) * @all(dvxdx)) + (@all(lambda) * (@all(dvzdz))))
    @all(tauzz) =
        @all(tauzz) - dt * ((@all(M) * @all(dvzdz)) + (@all(lambda) * (@all(dvxdx))))

    return
end
@parallel function compute_tauij!(
    tauxy::Data.Array{3},
    tauxz,
    tauyz,
    dvxdy,
    dvxdz,
    dvydx,
    dvydz,
    dvzdx,
    dvzdy,
    dt,
    mu,
)
    @all(tauxz) = @all(tauxz) - dt * (@av_xzi(mu) * (@all(dvxdz) + @all(dvzdx)))
    @all(tauxy) = @all(tauxy) - dt * (@av_yzi(mu) * (@all(dvxdy) + @all(dvydx)))
    @all(tauyz) = @all(tauyz) - dt * (@av_xyi(mu) * (@all(dvydz) + @all(dvzdy)))

    return
end
@parallel function compute_tauij!(tauxz::Data.Array{2}, dvxdz, dvzdx, dt, mu)
    @all(tauxz) = @all(tauxz) - dt * (@av(mu) * (@all(dvxdz) + @all(dvzdx)))
    return
end




@parallel_indices (iy, ix) function free_surface!(tau::Data.Array{3}, iz_free)
    tau[iz_free, iy, ix] = 0.0
    return
end
@parallel_indices (iy, ix) function free_surface_mirror!(tau::Data.Array{3}, iz_free)
    tau[iz_free-1, iy, ix] = -1 * tau[iz_free, iy, ix]
    return
end
@parallel_indices (ix) function free_surface!(tau::Data.Array{2}, iz_free)
    tau[iz_free, ix] = 0.0
    return
end
@parallel_indices (ix) function free_surface_mirror!(tau::Data.Array{2}, iz_free)
    tau[iz_free-1, ix] = -1 * tau[iz_free, ix]
    return
end




function advance_kernel!(pap, pac::T) where {T<:P_common{FdtdElastic,3}}
    w1t = pap.w1[:t]
    mw = pap.memory_pml
    pml = pac.pml

    @parallel compute_dtau!(
        w1t[:tauxx],
        w1t[:tauyy],
        w1t[:tauzz],
        w1t[:tauxy],
        w1t[:tauxz],
        w1t[:tauyz],
        w1t[:dtauxxdx],
        w1t[:dtauxydx],
        w1t[:dtauxzdx],
        w1t[:dtauyydy],
        w1t[:dtauxydy],
        w1t[:dtauyzdy],
        w1t[:dtauzzdz],
        w1t[:dtauyzdz],
        w1t[:dtauxzdz],
        pac.fc[:dxI],
        pac.fc[:dyI],
        pac.fc[:dzI],
    )

    memoryx!(
        mw[:dtauxxdx],
        w1t[:dtauxxdx],
        pml[:x][:a_half],
        pml[:x][:b_half],
        pml[:x][:k_halfI],
    )
    memory1x!(mw[:dtauxydx], w1t[:dtauxydx], pml[:x][:a], pml[:x][:b], pml[:x][:kI])
    memory1x!(mw[:dtauxzdx], w1t[:dtauxzdx], pml[:x][:a], pml[:x][:b], pml[:x][:kI])

    memoryy!(
        mw[:dtauyydy],
        w1t[:dtauyydy],
        pml[:y][:a_half],
        pml[:y][:b_half],
        pml[:y][:k_halfI],
    )
    memory1y!(mw[:dtauxydy], w1t[:dtauxydy], pml[:y][:a], pml[:y][:b], pml[:y][:kI])
    memory1y!(mw[:dtauyzdy], w1t[:dtauyzdy], pml[:y][:a], pml[:y][:b], pml[:y][:kI])

    memoryz!(
        mw[:dtauzzdz],
        w1t[:dtauzzdz],
        pml[:z][:a_half],
        pml[:z][:b_half],
        pml[:z][:k_halfI],
    )
    memory1z!(mw[:dtauyzdz], w1t[:dtauyzdz], pml[:z][:a], pml[:z][:b], pml[:z][:kI])
    memory1z!(mw[:dtauxzdz], w1t[:dtauxzdz], pml[:z][:a], pml[:z][:b], pml[:z][:kI])

    @parallel compute_v!(
        w1t[:vx],
        w1t[:vy],
        w1t[:vz],
        w1t[:dtauxxdx],
        w1t[:dtauxydx],
        w1t[:dtauxzdx],
        w1t[:dtauyydy],
        w1t[:dtauxydy],
        w1t[:dtauyzdy],
        w1t[:dtauzzdz],
        w1t[:dtauyzdz],
        w1t[:dtauxzdz],
        pac.fc[:dt],
        pac.mod[:rho],
    )

    @parallel (1:pac.ic[:nz], 1:pac.ic[:ny]) dirichletx!(
        w1t[:vx],
        w1t[:vz],
        w1t[:vy],
        pac.ic[:nx],
    )
    @parallel (1:pac.ic[:nz], 1:pac.ic[:nx]) dirichlety!(
        w1t[:vy],
        w1t[:vz],
        w1t[:vx],
        pac.ic[:ny],
    )
    @parallel (1:pac.ic[:ny], 1:pac.ic[:nx]) dirichletz!(
        w1t[:vz],
        w1t[:vy],
        w1t[:vx],
        pac.ic[:nz],
    )
    @parallel compute_dv!(
        w1t[:vx],
        w1t[:vy],
        w1t[:vz],
        w1t[:dvxdx],
        w1t[:dvydy],
        w1t[:dvzdz],
        w1t[:dvxdy],
        w1t[:dvxdz],
        w1t[:dvydx],
        w1t[:dvydz],
        w1t[:dvzdx],
        w1t[:dvzdy],
        pac.fc[:dxI],
        pac.fc[:dyI],
        pac.fc[:dzI],
    )
    memoryx!(mw[:dvxdx], w1t[:dvxdx], pml[:x][:a], pml[:x][:b], pml[:x][:kI])
    memoryy!(mw[:dvydy], w1t[:dvydy], pml[:y][:a], pml[:y][:b], pml[:y][:kI])
    memoryz!(mw[:dvzdz], w1t[:dvzdz], pml[:z][:a], pml[:z][:b], pml[:z][:kI])

    memoryy!(mw[:dvxdy], w1t[:dvxdy], pml[:y][:a_half], pml[:y][:b_half], pml[:y][:k_halfI])
    memoryz!(mw[:dvxdz], w1t[:dvxdz], pml[:z][:a_half], pml[:z][:b_half], pml[:z][:k_halfI])
    memoryx!(mw[:dvydx], w1t[:dvydx], pml[:x][:a_half], pml[:x][:b_half], pml[:x][:k_halfI])

    memoryz!(mw[:dvydz], w1t[:dvydz], pml[:z][:a_half], pml[:z][:b_half], pml[:z][:k_halfI])
    memoryx!(mw[:dvzdx], w1t[:dvzdx], pml[:x][:a_half], pml[:x][:b_half], pml[:x][:k_halfI])
    memoryy!(mw[:dvzdy], w1t[:dvzdy], pml[:y][:a_half], pml[:y][:b_half], pml[:y][:k_halfI])

    @parallel compute_tauii!(
        w1t[:tauxx],
        w1t[:tauyy],
        w1t[:tauzz],
        w1t[:dvxdx],
        w1t[:dvydy],
        w1t[:dvzdz],
        pac.fc[:dt],
        pac.mod[:M],
        pac.mod[:lambda],
    )
    @parallel compute_tauij!(
        w1t[:tauxy],
        w1t[:tauxz],
        w1t[:tauyz],
        w1t[:dvxdy],
        w1t[:dvxdz],
        w1t[:dvydx],
        w1t[:dvydz],
        w1t[:dvzdx],
        w1t[:dvzdy],
        pac.fc[:dt],
        pac.mod[:mu],
    )

    if (pac.ic[:iz_free] ≠ 0)
        @parallel (1:size(w1t[:tauzz], 2), size(w1t[:tauzz], 3)) free_surface_mirror!(
            w1t[:tauzz],
            pac.ic[:iz_free] + 1,
        )
        @parallel (1:size(w1t[:tauxz], 2), size(w1t[:tauxz], 3)) free_surface!(w1t[:tauxz], pac.ic[:iz_free])
        @parallel (1:size(w1t[:tauyz], 2), size(w1t[:tauyz], 3)) free_surface!(w1t[:tauyz], pac.ic[:iz_free])
    end

end


function advance_kernel!(pap, pac::T) where {T<:P_common{FdtdElastic,2}}
    w1t = pap.w1[:t]
    mw = pap.memory_pml
    pml = pac.pml

    @parallel compute_dtau!(
        w1t[:tauxx],
        w1t[:tauzz],
        w1t[:tauxz],
        w1t[:dtauxxdx],
        w1t[:dtauxzdx],
        w1t[:dtauzzdz],
        w1t[:dtauxzdz],
        pac.fc[:dxI],
        pac.fc[:dzI],
    )

    memoryx!(
        mw[:dtauxxdx],
        w1t[:dtauxxdx],
        pml[:x][:a_half],
        pml[:x][:b_half],
        pml[:x][:k_halfI],
    )
    memory1x!(mw[:dtauxzdx], w1t[:dtauxzdx], pml[:x][:a], pml[:x][:b], pml[:x][:kI])

    memoryz!(
        mw[:dtauzzdz],
        w1t[:dtauzzdz],
        pml[:z][:a_half],
        pml[:z][:b_half],
        pml[:z][:k_halfI],
    )
    memory1z!(mw[:dtauxzdz], w1t[:dtauxzdz], pml[:z][:a], pml[:z][:b], pml[:z][:kI])

    @parallel compute_v!(
        w1t[:vx],
        w1t[:vz],
        w1t[:dtauxxdx],
        w1t[:dtauxzdx],
        w1t[:dtauzzdz],
        w1t[:dtauxzdz],
        pac.fc[:dt],
        pac.mod[:rho],
    )

    @parallel (1:pac.ic[:nz]) dirichletx!(w1t[:vx], w1t[:vz], pac.ic[:nx])
    @parallel (1:pac.ic[:nx]) dirichletz!(w1t[:vz], w1t[:vx], pac.ic[:nz])
    @parallel compute_dv!(
        w1t[:vx],
        w1t[:vz],
        w1t[:dvxdx],
        w1t[:dvzdz],
        w1t[:dvxdz],
        w1t[:dvzdx],
        pac.fc[:dxI],
        pac.fc[:dzI],
    )
    memoryx!(mw[:dvxdx], w1t[:dvxdx], pml[:x][:a], pml[:x][:b], pml[:x][:kI])
    memoryz!(mw[:dvzdz], w1t[:dvzdz], pml[:z][:a], pml[:z][:b], pml[:z][:kI])

    memoryz!(mw[:dvxdz], w1t[:dvxdz], pml[:z][:a_half], pml[:z][:b_half], pml[:z][:k_halfI])

    memoryx!(mw[:dvzdx], w1t[:dvzdx], pml[:x][:a_half], pml[:x][:b_half], pml[:x][:k_halfI])

    @parallel compute_tauii!(
        w1t[:tauxx],
        w1t[:tauzz],
        w1t[:dvxdx],
        w1t[:dvzdz],
        pac.fc[:dt],
        pac.mod[:M],
        pac.mod[:lambda],
    )
    @parallel compute_tauij!(
        w1t[:tauxz],
        w1t[:dvxdz],
        w1t[:dvzdx],
        pac.fc[:dt],
        pac.mod[:mu],
    )


    if (pac.ic[:izfree] ≠ 0)
        @parallel (1:size(w1t[:tauzz], 2)) free_surface_mirror!(w1t[:tauzz], pac.ic[:izfree] + 1)
        @parallel (1:size(w1t[:tauxz], 2)) free_surface!(w1t[:tauxz], pac.ic[:izfree])
    end

end
