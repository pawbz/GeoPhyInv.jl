@inbounds @fastmath function advance!(pac::T, pap) where T<:P_common{<:FdtdElastic}
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

@parallel function compute_dtau!(tauxx::Data.Array, tauyy::Data.Array, tauzz::Data.Array, tauxy::Data.Array, tauxz::Data.Array, tauyz::Data.Array, 
    dtauxx_dx, dtauxy_dx, dtauxz_dx, dtauyy_dy, dtauxy_dy, dtauyz_dy, dtauzz_dz, dtauyz_dz, dtauxz_dz, dx::Data.Number, dy::Data.Number, dz::Data.Number)

    @all(dtauxx_dx) = @d_xi(tauxx) / dx # at [ix+1/2,iy+1,iz+1]
    @all(dtauxy_dx) = @d_xa(tauxy) / dx # at [ix+1,iy+1/2,iz+1]
    @all(dtauxz_dx) = @d_xa(tauxz) / dx # at [ix+1,iy+1,iz+1/2]


    @all(dtauyy_dy) = @d_yi(tauyy) / dy # at [ix+1,iy+1/2,iz+1]
    @all(dtauxy_dy) = @d_ya(tauxy) / dy # at [ix+1/2,iy+1,iz+1]
    @all(dtauyz_dy) = @d_ya(tauyz) / dy # at [ix+1,iy+1,iz+1/2]


    @all(dtauzz_dz) = @d_zi(tauzz) / dz # at [ix+1,iy+1,iz+1/2]
    @all(dtauxz_dz) = @d_za(tauxz) / dz # at [ix+1/2,iy+1,iz+1]
    @all(dtauyz_dz) = @d_za(tauyz) / dz # at [ix+1,iy+1/2,iz+1]

    return
end
# @parallel function compute_V!(vx::Data.Array, vy::Data.Array, vz::Data.Array, 
#     dtauxx_dx, dtauxy_dx, dtauxz_dx, dtauyy_dy, dtauxy_dy, dtauyz_dy, dtauzz_dz, dtauyz_dz, dtauxz_dz,
#     dt::Data.Number, rho)
  
#     @inn(pt[:vx]) = @inn(pt[:vx]) - dt / @av_xi(rho) * (@all(dtauxx_dx) + @all(dtauxy_dy) + @all(dtauxz_dz))
#     @inn(pt[:vy]) = @inn(pt[:vy]) - dt / @av_yi(rho) * (@all(dtauxy_dx) + @all(dtauyy_dy) + @all(dtauyz_dz))
#     @inn(pt[:vz]) = @inn(pt[:vz]) - dt / @av_zi(rho) * (@all(dtauxz_dx) + @all(dtauyz_dy) + @all(dtauzz_dz))

#     return
# end


@parallel function compute_dv!(dvx_dx,vx)
#     # vx, vy, vz,p[:dx][:vx],p[:dy][:vy],p[:dy][:vy],p[:dy][:vx],p[:dz][:vx],p[:dx][:vy],p[:dz][:vy],p[:dx][:vz],p[:dy][:vz], dx,dy,dz)
    # p=pap.w1
    # pt=p[:t]

@all(dvx_dx) = @d_za(vx) #/ dx # at [ix,iy,iz]
    # @all(pap.w1[:dx][:vx]) = @d_xa(p[:t][:vx]) / dx # at [ix,iy,iz]
    # @all(pap.w1[:dy][:vy]) = @d_ya(p[:t][:vy]) / dy # at      "
    # @all(pap.w1[:dy][:vy]) = @d_za(p[:t][:vz]) / dz # at      "


#     @all(p[:dy][:vx]) = @d_yi(pt[:vx]) / dy # at [ix+1/2,iy+1/2,iz+1]
#     @all(p[:dz][:vx]) = @d_zi(pt[:vx]) / dz # at [ix+1/2,iy+1,iz+1/2]

#     @all(p[:dz][:vy]) = @d_zi(pt[:vy]) / dz # at [ix+1,iy+1/2,iz+1/2]
#     @all(p[:dx][:vy]) = @d_xi(pt[:vy]) / dx # at [ix+1/2,iy+1/2,iz+1]

#     @all(p[:dx][:vz]) = @d_xi(pt[:vz]) / dx # at [ix+1/2,iy+1,iz+1/2]
#     @all(p[:dy][:vz]) = @d_yi(pt[:vz]) / dy # at [ix+1,iy+1/2,iz+1/2]


    return
end




# @parallel function compute_tauii!(tauxx::Data.Array, tauyy::Data.Array, tauzz::Data.Array,p[:dx][:vx]::Data.Array, p[:dy][:vy]::Data.Array, p[:dy][:vy]::Data.Array, dt::Data.Number, lambda2mu::Data.Array, lambda::Data.Array)

#     @all(tauxx) = @all(tauxx) - dt * ((@all(lambda2mu)*(@all(p[:dx][:vx])))  + (@all(lambda)*(@all(p[:dy][:vy]) + @all(p[:dy][:vy]))))
#     @all(tauyy) = @all(tauyy) - dt * ((@all(lambda2mu)*@all(p[:dy][:vy]))  + (@all(lambda)*(@all(p[:dx][:vx]) + @all(p[:dy][:vy]))))
#     @all(tauzz) = @all(tauzz) - dt * ((@all(lambda2mu)*@all(p[:dy][:vy]))  + (@all(lambda)*(@all(p[:dy][:vy]) + @all(p[:dx][:vx]))))
#     return
# end
# @parallel function compute_tauij!(tauxy::Data.Array, tauxz::Data.Array, tauyz::Data.Array,p[:dy][:vx],p[:dz][:vx],p[:dx][:vy],p[:dz][:vy],p[:dx][:vz],p[:dy][:vz], dt::Data.Number, mu::Data.Array)
#     @all(tauxz) = @all(tauxz) - dt * (@av_xzi(mu) * (@all(p[:dz][:vx]) + @all(p[:dx][:vz])))
#     @all(tauxy) = @all(tauxy) - dt * (@av_xyi(mu) * (@all(p[:dy][:vx]) + @all(p[:dx][:vy])))
#     @all(tauyz) = @all(tauyz) - dt * (@av_yzi(mu) * (@all(p[:dz][:vy]) + @all(p[:dy][:vz])))

#     return
# end


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




function advance_kernel!(pap, pac::T) where T<:P_common{FdtdElastic}
	@parallel compute_dv!(pap.w1[:t][:dvxdx],pap.w1[:t][:vx])
	# @parallel (1:size(p[:dx][:vx],1),1:size(p[:dx][:vx],2),1:size(p[:dx][:vx],3)) memoryx!(memory_p[:dx][:vx],p[:dx][:vx],pmlx[:a],pmlx[:b],pmlx[:kI])
	# @parallel (1:size(p[:dy][:vy],1),1:size(p[:dy][:vy],2),1:size(p[:dy][:vy],3)) memoryx!(memory_p[:dy][:vy],p[:dy][:vy],pmlx[:a],pmlx[:b],pmlx[:kI])
	# @parallel (1:size(p[:dy][:vy],1),1:size(p[:dy][:vy],2),1:size(p[:dy][:vy],3)) memoryx!(memory_p[:dy][:vy],p[:dy][:vy],pmlx[:a],pmlx[:b],pmlx[:kI])

	# @parallel (1:size(p[:dy][:vx],1),1:size(p[:dy][:vx],2),1:size(p[:dy][:vx],3)) memoryy!(memory_p[:dy][:vx],p[:dy][:vx],pmly[:a_half],pmly[:b_half],pmly[:k_halfI])
	# @parallel (1:size(p[:dz][:vx],1),1:size(p[:dz][:vx],2),1:size(p[:dz][:vx],3)) memoryy!(memory_p[:dz][:vx],p[:dz][:vx],pmly[:a_half],pmly[:b_half],pmly[:k_halfI])
	# @parallel (1:size(p[:dx][:vy],1),1:size(p[:dx][:vy],2),1:size(p[:dx][:vy],3)) memoryy!(memory_p[:dx][:vy],p[:dx][:vy],pmly[:a_half],pmly[:b_half],pmly[:k_halfI])

	# @parallel (1:size(p[:dz][:vy],1),1:size(p[:dz][:vy],2),1:size(p[:dz][:vy],3)) memoryx!(memory_p[:dz][:vy],p[:dz][:vy],pmlz[:a_half],pmlz[:b_half],pmlz[:k_halfI])
	# @parallel (1:size(p[:dx][:vz],1),1:size(p[:dx][:vz],2),1:size(p[:dx][:vz],3)) memoryx!(memory_p[:dx][:vz],p[:dx][:vz],pmlz[:a_half],pmlz[:b_half],pmlz[:k_halfI])
	# @parallel (1:size(p[:dy][:vz],1),1:size(p[:dy][:vz],2),1:size(p[:dy][:vz],3)) memoryx!(memory_p[:dy][:vz],p[:dy][:vz],pmlz[:a_half],pmlz[:b_half],pmlz[:k_halfI])

	# @parallel compute_tauii!(tauxx, tauyy, tauzz, p[:dx][:vx], p[:dy][:vy], p[:dy][:vy], step(tgrid), lambda2mu, lambda)
	# if(it < nt)
	# tauxx[div(nx, 2),div(ny, 2),div(nz, 2)] += wav[it]*5e8
	# end
	# @parallel compute_tauij!(tauxy, tauxz, tauyz, p[:dy][:vx],p[:dz][:vx],p[:dx][:vy],p[:dz][:vy],p[:dx][:vz],p[:dy][:vz], step(tgrid), mu)
	# @parallel compute_dtau!(tauxx, tauyy, tauzz, tauxy, tauxz, tauyz, dtauxx_dx, dtauxy_dx, dtauxz_dx, dtauyy_dy, dtauxy_dy, dtauyz_dy, dtauzz_dz, dtauyz_dz, dtauxz_dz, reverse(step.(mgrid))...)

	# @parallel (1:size(dtauxx_dx,1),1:size(dtauxx_dx,2),1:size(dtauxx_dx,3)) memoryx!(memory_dtauxx_dx, dtauxx_dx, pmlx[:a_half], pmlx[:b_half], pmlx[:k_halfI])
	# @parallel (1:size(dtauxy_dx,1),1:size(dtauxy_dx,2),1:size(dtauxy_dx,3)) memory1x!(memory_dtauxy_dx, dtauxy_dx, pmlx[:a], pmlx[:b], pmlx[:kI])
	# @parallel (1:size(dtauxz_dx,1),1:size(dtauxz_dx,2),1:size(dtauxz_dx,3)) memory1x!(memory_dtauxz_dx, dtauxz_dx, pmlx[:a], pmlx[:b], pmlx[:kI])

	# @parallel (1:size(dtauyy_dy,1),1:size(dtauyy_dy,2),1:size(dtauyy_dy,3)) memoryy!(memory_dtauyy_dy, dtauyy_dy, pmly[:a_half], pmly[:b_half], pmly[:k_halfI])
	# @parallel (1:size(dtauxy_dy,1),1:size(dtauxy_dy,2),1:size(dtauxy_dy,3)) memory1y!(memory_dtauxy_dy, dtauxy_dy, pmly[:a], pmly[:b], pmly[:kI])
	# @parallel (1:size(dtauyz_dy,1),1:size(dtauyz_dy,2),1:size(dtauyz_dy,3)) memory1y!(memory_dtauyz_dy, dtauyz_dy, pmly[:a], pmly[:b], pmly[:kI])

	# @parallel (1:size(dtauzz_dz,1),1:size(dtauzz_dz,2),1:size(dtauzz_dz,3)) memoryz!(memory_dtauzz_dz, dtauzz_dz, pmlz[:a_half], pmlz[:b_half], pmlz[:k_halfI])
	# @parallel (1:size(dtauyz_dz,1),1:size(dtauyz_dz,2),1:size(dtauyz_dz,3)) memory1z!(memory_dtauyz_dz, dtauyz_dz, pmlz[:a], pmlz[:b], pmlz[:kI])
	# @parallel (1:size(dtauxz_dz,1),1:size(dtauxz_dz,2),1:size(dtauxz_dz,3)) memory1z!(memory_dtauxz_dz, dtauxz_dz, pmlz[:a], pmlz[:b], pmlz[:kI])

	# @parallel compute_V!(vx, vy, vz, dtauxx_dx, dtauxy_dx, dtauxz_dx, dtauyy_dy, dtauxy_dy, dtauyz_dy, dtauzz_dz, dtauyz_dz, dtauxz_dz, step(tgrid), rho)

	# @parallel (1:ny,1:nz) dirichletx!(vx,vy,vz,nx)
	# @parallel (1:nx,1:nz) dirichlety!(vx,vy,vz,ny)
	# @parallel (1:nx,1:ny) dirichletz!(vx,vy,vz,nz)
end
