function advance!(pac::T, pap) where {T<:P_common{FdtdAcoustic}}
    for ipw in pac.activepw
        # store p for the last two steps
        # pppppp!(pap[ipw],pac.attrib_mod)
        advance_kernel!(pap[ipw], pac)
    end
    return nothing
end

function advance_kernel!(pap, pac::T) where {T<:P_common{FdtdAcoustic,2}}

    w1t = pap.w1[:t]
    mw = pap.memory_pml
    pml = pac.pml
    pml_faces = pac.pml_faces

    #compute dpdx and dpdz at [it-1] for all propagating fields
    @parallel compute_dp!(w1t[:p], w1t[:dpdx], w1t[:dpdz], pac.fc[:dzI], pac.fc[:dxI])
    memoryx!(
        mw[:dpdx],
        w1t[:dpdx],
        pml[:dpdx][:a],
        pml[:dpdx][:b],
        pml[:dpdx][:kI],
        pml_faces,
    )
    memoryz!(
        mw[:dpdz],
        w1t[:dpdz],
        pml[:dpdz][:a],
        pml[:dpdz][:b],
        pml[:dpdz][:kI],
        pml_faces,
    )

    #update velocity at [it-1/2] using 
    #velocity at [it-3/2] and dpdx and dpdz at [it-1] 
    @parallel compute_v!(
        w1t[:vx],
        w1t[:vz],
        pac.mod[:rhoI],
        w1t[:dpdx],
        w1t[:dpdz],
        pac.fc[:dt],
    )#

    #rigid boundary conditions 
    (:xmin ∈ pac.rigid_faces) &&
        @parallel (1:pac.ic[:nz]) dirichletxmin!(w1t[:vx], w1t[:vz], pac.ic[:nx])
    (:xmax ∈ pac.rigid_faces) &&
        @parallel (1:pac.ic[:nz]) dirichletxmax!(w1t[:vx], w1t[:vz], pac.ic[:nx])
    (:zmin ∈ pac.rigid_faces) &&
        @parallel (1:pac.ic[:nx]) dirichletzmin!(w1t[:vz], w1t[:vx], pac.ic[:nz])
    (:zmax ∈ pac.rigid_faces) &&
        @parallel (1:pac.ic[:nx]) dirichletzmax!(w1t[:vz], w1t[:vx], pac.ic[:nz])

    @parallel compute_dv!(
        w1t[:vx],
        w1t[:vz],
        w1t[:dvxdx],
        w1t[:dvzdz],
        pac.fc[:dxI],
        pac.fc[:dzI],
    )#
    memoryx!(
        mw[:dvxdx],
        w1t[:dvxdx],
        pml[:dvxdx][:a],
        pml[:dvxdx][:b],
        pml[:dvxdx][:kI],
        pml_faces,
    )
    memoryz!(
        mw[:dvzdz],
        w1t[:dvzdz],
        pml[:dvzdz][:a],
        pml[:dvzdz][:b],
        pml[:dvzdz][:kI],
        pml_faces,
    )

    #compute pressure at [it] using p at [it-1] and dvxdx
    #and dvzdz at [it-1/2]
    @parallel compute_p!(w1t[:p], w1t[:dvxdx], w1t[:dvzdz], pac.mod[:K], pac.fc[:dt])

end

function advance_kernel!(pap, pac::T) where {T<:P_common{FdtdAcoustic,3}}

    w1t = pap.w1[:t]
    mw = pap.memory_pml
    pml = pac.pml
    pml_faces = pac.pml_faces

    #compute dpdx and dpdz at [it-1] for all propagating fields
    @parallel compute_dp!(
        w1t[:p],
        w1t[:dpdx],
        w1t[:dpdy],
        w1t[:dpdz],
        pac.fc[:dzI],
        pac.fc[:dyI],
        pac.fc[:dxI],
    )
    memoryx!(
        mw[:dpdx],
        w1t[:dpdx],
        pml[:dpdx][:a],
        pml[:dpdx][:b],
        pml[:dpdx][:kI],
        pml_faces,
    )
    memoryy!(
        mw[:dpdy],
        w1t[:dpdy],
        pml[:dpdy][:a],
        pml[:dpdy][:b],
        pml[:dpdy][:kI],
        pml_faces,
    )
    memoryz!(
        mw[:dpdz],
        w1t[:dpdz],
        pml[:dpdz][:a],
        pml[:dpdz][:b],
        pml[:dpdz][:kI],
        pml_faces,
    )

    #update velocity at [it-1/2] using 
    #velocity at [it-3/2] and dpdx and dpdz at [it-1] 
    @parallel compute_v!(
        w1t[:vx],
        w1t[:vy],
        w1t[:vz],
        pac.mod[:rhoI],
        w1t[:dpdx],
        w1t[:dpdy],
        w1t[:dpdz],
        pac.fc[:dt],
    )#

    #rigid boundary conditions 
    (:xmin ∈ pac.rigid_faces) && @parallel (1:pac.ic[:nz], 1:pac.ic[:ny]) dirichletxmin!(
        w1t[:vx],
        w1t[:vz],
        w1t[:vy],
        pac.ic[:nx],
    )
    (:xmax ∈ pac.rigid_faces) && @parallel (1:pac.ic[:nz], 1:pac.ic[:ny]) dirichletxmax!(
        w1t[:vx],
        w1t[:vz],
        w1t[:vy],
        pac.ic[:nx],
    )
    (:ymin ∈ pac.rigid_faces) && @parallel (1:pac.ic[:nz], 1:pac.ic[:nx]) dirichletymin!(
        w1t[:vy],
        w1t[:vz],
        w1t[:vx],
        pac.ic[:ny],
    )
    (:ymax ∈ pac.rigid_faces) && @parallel (1:pac.ic[:nz], 1:pac.ic[:nx]) dirichletymax!(
        w1t[:vy],
        w1t[:vz],
        w1t[:vx],
        pac.ic[:ny],
    )

    (:zmin ∈ pac.rigid_faces) && @parallel (1:pac.ic[:ny], 1:pac.ic[:nx]) dirichletzmin!(
        w1t[:vz],
        w1t[:vy],
        w1t[:vx],
        pac.ic[:nz],
    )

    (:zmax ∈ pac.rigid_faces) && @parallel (1:pac.ic[:ny], 1:pac.ic[:nx]) dirichletzmax!(
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
        pac.fc[:dxI],
        pac.fc[:dyI],
        pac.fc[:dzI],
    )#
    memoryx!(
        mw[:dvxdx],
        w1t[:dvxdx],
        pml[:dvxdx][:a],
        pml[:dvxdx][:b],
        pml[:dvxdx][:kI],
        pml_faces,
    )
    memoryy!(
        mw[:dvydy],
        w1t[:dvydy],
        pml[:dvydy][:a],
        pml[:dvydy][:b],
        pml[:dvydy][:kI],
        pml_faces,
    )
    memoryz!(
        mw[:dvzdz],
        w1t[:dvzdz],
        pml[:dvzdz][:a],
        pml[:dvzdz][:b],
        pml[:dvzdz][:kI],
        pml_faces,
    )

    #compute pressure at [it] using p at [it-1] and dvxdx
    #and dvzdz at [it-1/2]
    @parallel compute_p!(
        w1t[:p],
        w1t[:dvxdx],
        w1t[:dvydy],
        w1t[:dvzdz],
        pac.mod[:K],
        pac.fc[:dt],
    )

end

"""
Exchange pointers (i.e., set names of NamedArray) instead of copying arrays around
"""
function pppppp!(pap, attrib_mod)
    w2 = pap.w1
    names_old = names(w2)[1]
    names_new = vcat(circshift(names_old[1:3], -1), names_old[4:end])
    setnames!(w2, names_new, 1)
    # (old method), use for debugging
    #copyto!.(w2[:tpp],w2[:tp])
    #copyto!.(w2[:tp],w2[:t])
end


# these relative indices of the arrays point to same location
# [ix,iy,iz]     --> pressure
# [ix-1/2,iy,iz]       --> vx
# [ix,iy-1/2,iz]       --> vy
# [ix,iy,iz-1/2]       --> vz

@parallel function compute_dp!(p, dpdx, dpdz, dzI, dxI)
    @all(dpdx) = @d_xi(p) * dxI #
    @all(dpdz) = @d_zi(p) * dzI #
    return
end


@parallel function compute_dp!(p, dpdx, dpdy, dpdz, dzI, dyI, dxI)
    @all(dpdx) = @d_xi(p) * dxI #
    @all(dpdy) = @d_yi(p) * dyI #
    @all(dpdz) = @d_zi(p) * dzI #
    return
end


@parallel function compute_v!(vx, vz, rhoI, dpdx, dpdz, dt)#
    @inn(vx) = @inn(vx) + dt * @av_xi(rhoI) * @all(dpdx)
    @inn(vz) = @inn(vz) + dt * @av_zi(rhoI) * @all(dpdz)
    return
end

@parallel function compute_v!(vx, vy, vz, rhoI, dpdx, dpdy, dpdz, dt)#
    @inn(vx) = @inn(vx) + dt * @av_xi(rhoI) * @all(dpdx)
    @inn(vy) = @inn(vy) + dt * @av_yi(rhoI) * @all(dpdy)
    @inn(vz) = @inn(vz) + dt * @av_zi(rhoI) * @all(dpdz)
    return
end


# 
@parallel function compute_dv!(vx, vz, dvxdx, dvzdz, dxI, dzI)#
    @all(dvxdx) = @d_xa(vx) * dxI # 
    @all(dvzdz) = @d_za(vz) * dzI # 
    return
end

@parallel function compute_dv!(vx, vy, vz, dvxdx, dvydy, dvzdz, dxI, dyI, dzI)#
    @all(dvxdx) = @d_xa(vx) * dxI # 
    @all(dvydy) = @d_ya(vy) * dyI # 
    @all(dvzdz) = @d_za(vz) * dzI # 
    return
end


# no attenuation (no memory in stress-strain relation)
@parallel function compute_p!(p, dvxdx, dvzdz, K, dt)
    @all(p) = @all(p) + @all(K) * (@all(dvxdx) + @all(dvzdz)) * dt
    return
end


# no attenuation (no memory in stress-strain relation)
@parallel function compute_p!(p, dvxdx, dvydy, dvzdz, K, dt)
    @all(p) = @all(p) + @all(K) * (@all(dvxdx) + @all(dvzdz) + @all(dvydy)) * dt
    return
end




# # 
#  function compute_p!(pap,pac,::T) where {T<:Union{FdtdAcousticVisco}}
# 	dvxdx=pap.w1[:dx][:vx]
# 	dvzdz=pap.w1[:dz][:vz]
# 	p=pap.w1[:t][:p]
# 	pp=pap.w1[:tp][:p]
# 	r=pap.w2[:t][:r]
# 	rp=pap.w2[:tp][:r]
# 	prvzvx!(p,pp,r,rp,dvzdz,dvxdx,pac.mod[:K],pac.mod3[:memcoeff1],pac.mod3[:memcoeff2],pac.fc[:dt],pac.ic[:nz],pac.ic[:nx],pac.ic[:nsls])
# end
# # viscoacoustic modeling (memory in stress-strain relation)
#  function prvzvx!(p,pp,r,rp,dvxdx,dvzdz,mod,mod31,mod32,fc1,nz,nx,nsls)#pw,dpdxw,dpdzw,modK,nz,nx,dt, ::FdtdAcousticVisco)
# 	for ix=1:nx  # see limits above
# 	@simd for iz=1:nz

# 		 # use the Auxiliary Differential Equation form, 
# 		 # which is second-order accurate in time if implemented 
# 		 #following eq (14) of Robertsson, Blanch and Symes, Geophysics, vol.
# 		 #	  59(9), pp 1444-1456 (1994), which is what we do here        

# 		 # loop over standard linear solids
# 		 sumr=0.0
# 		 for isls = 1:nsls
# 			# this average of the two terms comes from eq (14) of 
# 			#Robertsson, Blanch and Symes, Geophysics, vol. 59(9), pp 1444-1456 (1994) 
# 			# central finite-difference around the time of vx (not the staggerred time grid!); therefore the factor 0.5         
# 			r[isls,iz,ix] = r[isls,iz,ix] + ((dvxdx[iz,ix] + dvzdz[iz,ix]) * mod[iz,ix] * mod32[isls,iz,ix] - (mod31[isls,iz,ix] * r[isls,iz,ix])) * fc1 * 0.5
# 			sumr += r[isls,iz,ix] 
# 		 end

# 		# this average of the two terms comes from eq (13) of 
# 		#Robertsson, Blanch and Symes, Geophysics, vol. 59(9), pp 1444-1456 (1994)  
# 		#pressure(i,j) = pressure(i,j) + (- kappa_half_x * (value_dvx_dx + value_dvy_dy) + 
# 		#  0.5d0 * sum_of_memory_variables_kappa) * DELTAT
# 		# Here sumr should be at the time same as dvxdx
# 		@inbounds p[iz,ix] = pp[iz,ix] + ((mod[iz,ix] * (dvxdx[iz,ix] + dvzdz[iz,ix])) + (sumr)) * fc1  
# 	end
# 	end
# end





# ====================================
# not used at the moment
# ====================================


# function dtsample(p)
#     p[1] .= 0.0  # it is necessary to have these boundary conditions
#     nt = length(p)
#     dpdt = zero(p)
#     pout = zero(p)
#     for it = 1:nt-1 # use these bounds later
#         dpdt[it] = p[it+1] - p[it]
#     end

#     for it = 1:nt-1 # use these exact bounds later
#         pout[it+1] = dpdt[it+1] - dpdt[it]
#     end
#     return pout
# end

# """
# Simple test to see staggerred grid spatial derivatives that are used satisfy adjoint test.
# """
# function Δxsample(p)
#     nx = length(p)
#     # boundary conditions
#     p[1:2] .= 0.0
#     p[nx-1:nx] .= 0.0
#     dpdx = zero(p)
#     pout = zero(p)
#     for i = 1:100 # applying operators a 100 times, say
#         for ix = 2:nx-2 # use these bounds later
#             dpdx[ix] = (27.e0 * p[ix+1] - 27.e0 * p[ix] - p[ix+2] + p[ix-1]) * 1.0
#         end

#         for ix = 3:nx-1 # use these exact bounds later
#             pout[ix] =
#                 (27.e0 * dpdx[ix] - 27.e0 * dpdx[ix-1] - dpdx[ix+1] + dpdx[ix-2]) * 1.0
#         end
#     end

#     return pout
# end


# """
# advance by one time step for adjoint test
# """
# function advance_sample(ntimes, p, modrhovxI, modrhovzI, modK)
#     nz = size(p, 1)
#     nx = size(p, 2)
#     # boundary conditions
#     for i = 1:3
#         p[1:2, :, i] .= 0.0
#         p[:, 1:2, i] .= 0.0
#         p[nz-1:nz, :, i] .= 0.0
#         p[:, nx-1:nx, i] .= 0.0
#     end


#     dx24I = 1.0
#     dz24I = 1.0
#     dt = 1.0

#     a_x = zeros(nx)
#     a_x_half = zeros(nx)
#     b_x = zeros(nx)
#     b_x_half = zeros(nx)
#     k_x = ones(nx)
#     k_xI = ones(nx)
#     k_x_half = ones(nx)
#     k_x_halfI = ones(nx)

#     k_z = ones(nz)
#     k_zI = ones(nz)
#     k_z_half = ones(nz)
#     k_z_halfI = ones(nz)
#     a_z = zeros(nz)
#     a_z_half = zeros(nz)
#     b_z = zeros(nz)
#     b_z_half = zeros(nz)

#     modrhovxI = ones(nz, nx)
#     modrhovzI = ones(nz, nx)
#     modK = ones(nz, nx)

#     pw = zeros(nz, nx, 3)
#     #input
#     copyto!(pw, p)
#     #=	
#     	for iz in 1:nz
#     		for ix in 1:nx
#     			pw[iz,ix,3]=p[iz,ix,3]
#     		end
#     	end
#     =#
#     dpdxw = zeros(nz, nx, 3)
#     dpdzw = zeros(nz, nx, 3)
#     ppw = zeros(nz, nx, 3)
#     pppw = zeros(nz, nx, 3)

#     memory_dvx_dxw = zeros(nz, nx)
#     memory_dvx_dzw = zeros(nz, nx)
#     memory_dvz_dxw = zeros(nz, nx)
#     memory_dvz_dzw = zeros(nz, nx)
#     memory_dp_dxw = zeros(nz, nx)
#     memory_dp_dzw = zeros(nz, nx)

#     for i = 1:ntimes

#         advance_kernel!(
#             pw,
#             dpdxw,
#             dpdzw,
#             dx24I,
#             dz24I,
#             memory_dp_dxw,
#             memory_dp_dzw,
#             b_x_half,
#             b_z_half,
#             a_x_half,
#             a_z_half,
#             k_x_halfI,
#             k_z_halfI,
#             nx,
#             nz,
#             dt,
#             modrhovxI,
#             modrhovzI,
#             memory_dvx_dxw,
#             memory_dvz_dzw,
#             b_x,
#             b_z,
#             a_x,
#             a_z,
#             k_xI,
#             k_zI,
#             modK,
#         )
#     end


#     return pw
# end


