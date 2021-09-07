@inbounds @fastmath function advance!(pac::T, pap) where {T<:P_common{<:FdtdElastic}}
    for ipw in pac.activepw
        # store p for the last two steps
        # pppppp!(pap[ipw],pac.attrib_mod)
        advance_kernel!(pap[ipw], pac)
    end
    return nothing
end


# these relative indices of the arrays point to same location
# [ix,iy,iz]     --> tauxx, tauyy, and tauzz grid
# [ix+1/2,iy+1/2,iz+1]      --> tauxy
# [ix+1/2,iy+1,iz+1/2]      -->  tauxz
# [ix-1/2,iy,iz]       --> vx
# [ix,iy-1/2,iz]       --> vy
# [ix,iy,iz-1/2]       --> vz

@parallel function compute_dtau!(
    tauxx::Data.Array,
    tauyy::Data.Array,
    tauzz::Data.Array,
    tauxy::Data.Array,
    tauxz::Data.Array,
    tauyz::Data.Array,
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

    @all(dtauxxdx) = @d_zi(tauxx) * dxI # at [ix+1/2,iy+1,iz+1]
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
@parallel function compute_v!(
    vx,
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
        @inn(vx) - dt / @av_xi(rho) * (@all(dtauxxdx) + @all(dtauxydy) + @all(dtauxzdz))
    @inn(vy) =
        @inn(vy) - dt / @av_yi(rho) * (@all(dtauxydx) + @all(dtauyydy) + @all(dtauyzdz))
    @inn(vz) =
        @inn(vz) - dt / @av_zi(rho) * (@all(dtauxzdx) + @all(dtauyzdy) + @all(dtauzzdz))

    return
end


@parallel function compute_dv!(
    vx,
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




@parallel function compute_tauii!(tauxx, tauyy, tauzz, dvxdx, dvydy, dvzdz, dt, M, lambda)

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
@parallel function compute_tauij!(
    tauxy,
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
    @all(tauxy) = @all(tauxy) - dt * (@av_xyi(mu) * (@all(dvxdy) + @all(dvydx)))
    @all(tauyz) = @all(tauyz) - dt * (@av_yzi(mu) * (@all(dvydz) + @all(dvzdy)))

    return
end


# for dim in [:x, :y, :z]
#     fname = Symbol("memory", string(dim), "!")
#     i = Symbol("i", string(dim))
#     @eval @parallel_indices (ix, iy, iz) function $fname(
#         memory::Data.Array,
#         d::Data.Array,
#         a,
#         b,
#         kI,
#     )
#         memory[ix, iy, iz] = b[$i] * memory[ix, iy, iz] + a[$i] * d[ix, iy, iz]
#         d[ix, iy, iz] = d[ix, iy, iz] * kI[$i] + memory[ix, iy, iz]
#         return
#     end
#     fname = Symbol("memory1", string(dim), "!")
#     @eval @parallel_indices (ix, iy, iz) function $fname(
#         memory::Data.Array,
#         d::Data.Array,
#         a,
#         b,
#         kI,
#     )
#         memory[ix, iy, iz] = b[$i+1] * memory[ix, iy, iz] + a[$i+1] * d[ix, iy, iz]
#         d[ix, iy, iz] = d[ix, iy, iz] * kI[$i+1] + memory[ix, iy, iz]
#         return
#     end
#     fname = Symbol("dirichlet", string(dim), "!")
#     is = [:ix, :iy, :iz]

#     is1, is2, is3 = [replace(is, i => ii) for ii = 1:3]
#     isn = replace(is, i => :n)
#     isnm1 = replace(is, i => :(n - 1))
#     isnp1 = replace(is, i => :(n + 1))
#     i1, i2 = filter(x -> x != i, is)
#     v = Symbol("v", string(dim))
#     vs = (:vx, :vy, :vz)
#     v1, v2 = filter(x -> x != v, vs)

#     println(string(is1,is2,is3,i1,i2,v,v1,v2,isnp1))

#     @eval @parallel_indices ($i1, $i2) function $fname(vx, vy, vz, n)
#         # relative to the tauii grid, velocity at i=2 should be zero, so making use of ghost cells 
#         # $v[$(is1...)] = -$v[$(is2...)]
#         # $v[$(isnp1...)] = -$v[$(isn...)]
#         # along other dimensions velocity grid matches tauii grid
#         $v1[$(is1...)] = 0
#         $v1[$(isn...)] = 0
#         $v2[$(is1...)] = 0
#         $v2[$(isn...)] = 0


#         $v[$(is1...)] = -$v[$(is2...)]
#         $v[$(isnp1...)] = -$v[$(isn...)]
#         return
#     end
# end




function advance_kernel!(pap, pac::T) where {T<:P_common{FdtdElastic}}
    w1t = pap.w1[:t]
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
    # @parallel (1:size(dvxdx,1),1:size(dvxdx,2),1:size(dvxdx,3)) memoryx!(memory_dvxdx,dvxdx,pmlx[:a],pmlx[:b],pmlx[:kI])
    # @parallel (1:size(dvydy,1),1:size(dvydy,2),1:size(dvydy,3)) memoryx!(memory_dvydy,dvydy,pmlx[:a],pmlx[:b],pmlx[:kI])
    # @parallel (1:size(dvydy,1),1:size(dvydy,2),1:size(dvydy,3)) memoryx!(memory_dvydy,dvydy,pmlx[:a],pmlx[:b],pmlx[:kI])

    # @parallel (1:size(dvxdy,1),1:size(dvxdy,2),1:size(dvxdy,3)) memoryy!(memory_dvxdy,dvxdy,pmly[:a_half],pmly[:b_half],pmly[:k_halfI])
    # @parallel (1:size(dvxdz,1),1:size(dvxdz,2),1:size(dvxdz,3)) memoryy!(memory_dvxdz,dvxdz,pmly[:a_half],pmly[:b_half],pmly[:k_halfI])
    # @parallel (1:size(dvydx,1),1:size(dvydx,2),1:size(dvydx,3)) memoryy!(memory_dvydx,dvydx,pmly[:a_half],pmly[:b_half],pmly[:k_halfI])

    # @parallel (1:size(dvydz,1),1:size(dvydz,2),1:size(dvydz,3)) memoryx!(memory_dvydz,dvydz,pmlz[:a_half],pmlz[:b_half],pmlz[:k_halfI])
    # @parallel (1:size(dvzdx,1),1:size(dvzdx,2),1:size(dvzdx,3)) memoryx!(memory_dvzdx,dvzdx,pmlz[:a_half],pmlz[:b_half],pmlz[:k_halfI])
    # @parallel (1:size(dvzdy,1),1:size(dvzdy,2),1:size(dvzdy,3)) memoryx!(memory_dvzdy,dvzdy,pmlz[:a_half],pmlz[:b_half],pmlz[:k_halfI])

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
    # if(it < nt)
    # tauxx[div(nx, 2),div(ny, 2),div(nz, 2)] += wav[it]*5e8
    # end
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

    # @parallel (1:size(dtauxxdx,1),1:size(dtauxxdx,2),1:size(dtauxxdx,3)) memoryx!(memory_dtauxxdx, dtauxxdx, pmlx[:a_half], pmlx[:b_half], pmlx[:k_halfI])
    # @parallel (1:size(dtauxydx,1),1:size(dtauxydx,2),1:size(dtauxydx,3)) memory1x!(memory_dtauxydx, dtauxydx, pmlx[:a], pmlx[:b], pmlx[:kI])
    # @parallel (1:size(dtauxzdx,1),1:size(dtauxzdx,2),1:size(dtauxzdx,3)) memory1x!(memory_dtauxzdx, dtauxzdx, pmlx[:a], pmlx[:b], pmlx[:kI])

    # @parallel (1:size(dtauyy_dy,1),1:size(dtauyy_dy,2),1:size(dtauyy_dy,3)) memoryy!(memory_dtauyy_dy, dtauyy_dy, pmly[:a_half], pmly[:b_half], pmly[:k_halfI])
    # @parallel (1:size(dtauxy_dy,1),1:size(dtauxy_dy,2),1:size(dtauxy_dy,3)) memory1y!(memory_dtauxy_dy, dtauxy_dy, pmly[:a], pmly[:b], pmly[:kI])
    # @parallel (1:size(dtauyz_dy,1),1:size(dtauyz_dy,2),1:size(dtauyz_dy,3)) memory1y!(memory_dtauyz_dy, dtauyz_dy, pmly[:a], pmly[:b], pmly[:kI])

    # @parallel (1:size(dtauzz_dz,1),1:size(dtauzz_dz,2),1:size(dtauzz_dz,3)) memoryz!(memory_dtauzz_dz, dtauzz_dz, pmlz[:a_half], pmlz[:b_half], pmlz[:k_halfI])
    # @parallel (1:size(dtauyz_dz,1),1:size(dtauyz_dz,2),1:size(dtauyz_dz,3)) memory1z!(memory_dtauyz_dz, dtauyz_dz, pmlz[:a], pmlz[:b], pmlz[:kI])
    # @parallel (1:size(dtauxz_dz,1),1:size(dtauxz_dz,2),1:size(dtauxz_dz,3)) memory1z!(memory_dtauxz_dz, dtauxz_dz, pmlz[:a], pmlz[:b], pmlz[:kI])

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

    # @parallel (1:ny,1:nz) dirichletx!(vx,vy,vz,nx)
    # @parallel (1:nx,1:nz) dirichlety!(vx,vy,vz,ny)
    # @parallel (1:nx,1:ny) dirichletz!(vx,vy,vz,nz)
end
