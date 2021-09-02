@inbounds @fastmath function advance!(pac, pap)
	for ipw in pac.activepw
		# store p for the last two steps
		pppppp!(pap[ipw],pac.attrib_mod)
		advance_kernel!(pap[ipw], pac, pac.attrib_mod)
	end
	return nothing
end

function advance_kernel!(pap,pac,attrib_mod)
	#compute dpdx and dpdz at [it-1] for all propagating fields
	update_dpdx!(pap,pac,attrib_mod)
	update_dpdz!(pap,pac,attrib_mod)

	#update velocity at [it-1/2] using 
	#velocity at [it-3/2] and dpdx and dpdz at [it-1] 
	update_vx!(pap,pac,attrib_mod)
	update_vz!(pap,pac,attrib_mod)

	#rigid boundary conditions (not currently implemented, need to do adjoint test)
	#such that:
	#p[3,.] is where velocities are zero and
	#p[.,n-1] is where velocities are zero
#	rigid_boundary_left!(p,nx,nz)
#	rigid_boundary_right!(p,nx,nz)
#	rigid_boundary_top!(p,nx,nz)
#	rigid_boundary_bottom!(p,nx,nz)


	update_dvdx!(pap,pac,attrib_mod)
	update_dvdz!(pap,pac,attrib_mod)

	#compute pressure at [it] using p at [it-1] and dvxdx
	#and dvzdz at [it-1/2]
	update_p!(pap,pac,attrib_mod)

end

"""
Exchange pointers (i.e., set names of NamedArray) instead of copying arrays around
"""
function pppppp!(pap, attrib_mod)
	w2=pap.w2
	names_old=names(w2)[1]
	names_new=vcat(circshift(names_old[1:3],-1), names_old[4:end])
	setnames!(w2, names_new, 1)
	# (old method), use for debugging
	#copyto!.(w2[:tpp],w2[:tp])
	#copyto!.(w2[:tp],w2[:t])
end



# 
@inbounds @fastmath function update_dvdx!(pap,pac,attrib_mod)#
	dvdx=pap.w2[:dx][:vx]
	vx=pap.w2[:t][:vx]
	mp=pap.memory_pml[:dvxdx]
	dvdx!(dvdx,vx,mp,pac.pml[:x][:b],pac.pml[:x][:a],pac.pml[:x][:kI],pac.fc[:δx24I],pac.ic[:nz],pac.ic[:nx])
end
@inbounds @fastmath function dvdx!(dvdx,vx,mp,pml1,pml2,pml3,fc1,nz,nx)
	for ix=3:nx-1 # ix=1,2,nx are useless, rigid boundary conditions apply on ix=3,nx-1
	@simd for iz=1:nz
		@inbounds dvdx[iz,ix] = (27.e0*vx[iz,ix]-27.e0*vx[iz,ix-1]-vx[iz,ix+1]+vx[iz,ix-2]) * fc1
		@inbounds mp[iz,ix] = pml1[ix] * mp[iz,ix] + pml2[ix] * dvdx[iz,ix] # pml 
		@inbounds dvdx[iz,ix] = muladd(dvdx[iz,ix] , pml3[ix] , mp[iz,ix]) # pml
	end
	end
end


@inbounds @fastmath function update_dvdz!(pap,pac,attrib_mod)#
	dvdz=pap.w2[:dz][:vz]
	vz=pap.w2[:t][:vz]
	mp=pap.memory_pml[:dvzdz]
	dvdz!(dvdz,vz,mp,pac.pml[:z][:b],pac.pml[:z][:a],pac.pml[:z][:kI],pac.fc[:δz24I],pac.ic[:nz],pac.ic[:nx])
end
@inbounds @fastmath function dvdz!(dvdz,vz,mp,pml1,pml2,pml3,fc1,nz,nx)
	for ix=1:nx
	@simd for iz=3:nz-1 # iz=1,2 are usless, boundary conditions on iz=3,nz-1
		@inbounds dvdz[iz,ix] = (27.e0*vz[iz,ix]-27.e0*vz[iz-1,ix]-vz[iz+1,ix]+vz[iz-2,ix]) * fc1
		@inbounds mp[iz,ix] = pml1[iz] * mp[iz,ix] + pml2[iz] * dvdz[iz,ix] # pml
		@inbounds dvdz[iz,ix] = muladd(dvdz[iz,ix] , pml3[iz] , mp[iz,ix]) # pml
	end
	end
end

# no attenuation (no memory in stress-strain relation)
@inbounds @fastmath function update_p!(pap,pac,::T) where {T<:Union{FdtdAcou,FdtdAcouBorn}}
	dvdx=pap.w2[:dx][:vx]
	dvdz=pap.w2[:dz][:vz]
	p=pap.w2[:t][:p]
	pp=pap.w2[:tp][:p]
	pvzvx!(p,pp,dvdz,dvdx,pac.mod[:K],pac.fc[:δt],pac.ic[:nz],pac.ic[:nx])
end
@inbounds @fastmath function pvzvx!(p,pp,dvdz,dvdx,mod,fc1,nz,nx)
	for ix=1:nx  # see limits above
	@simd for iz=1:nz
		@inbounds p[iz,ix] = pp[iz,ix] + (mod[iz,ix] * (dvdx[iz,ix] + dvdz[iz,ix])) * fc1 #* boundary_p(iz,ix)
	end
	end
end

@inbounds @fastmath function update_p!(pap,pac,::T) where {T<:Union{FdtdAcouVisco}}
	dvdx=pap.w2[:dx][:vx]
	dvdz=pap.w2[:dz][:vz]
	p=pap.w2[:t][:p]
	pp=pap.w2[:tp][:p]
	r=pap.w3[:t][:r]
	rp=pap.w3[:tp][:r]
	prvzvx!(p,pp,r,rp,dvdz,dvdx,pac.mod[:K],pac.mod3[:memcoeff1],pac.mod3[:memcoeff2],pac.fc[:δt],pac.ic[:nz],pac.ic[:nx],pac.ic[:nsls])
end
# viscoacoustic modeling (memory in stress-strain relation)
@inbounds @fastmath function prvzvx!(p,pp,r,rp,dvdx,dvdz,mod,mod31,mod32,fc1,nz,nx,nsls)#pw,dpdxw,dpdzw,modK,nz,nx,δt, ::FdtdAcouVisco)
	for ix=1:nx  # see limits above
	@simd for iz=1:nz

		 # use the Auxiliary Differential Equation form, 
		 # which is second-order accurate in time if implemented 
		 #following eq (14) of Robertsson, Blanch and Symes, Geophysics, vol.
		 #	  59(9), pp 1444-1456 (1994), which is what we do here        
		 
		 # loop over standard linear solids
		 sumr=0.0
		 for isls = 1:nsls
			# this average of the two terms comes from eq (14) of 
			#Robertsson, Blanch and Symes, Geophysics, vol. 59(9), pp 1444-1456 (1994) 
			# central finite-difference around the time of vx (not the staggerred time grid!); therefore the factor 0.5         
			r[isls,iz,ix] = r[isls,iz,ix] + ((dvdx[iz,ix] + dvdz[iz,ix]) * mod[iz,ix] * mod32[isls,iz,ix] - (mod31[isls,iz,ix] * r[isls,iz,ix])) * fc1 * 0.5
			sumr += r[isls,iz,ix] 
		 end

		# this average of the two terms comes from eq (13) of 
		#Robertsson, Blanch and Symes, Geophysics, vol. 59(9), pp 1444-1456 (1994)  
		#pressure(i,j) = pressure(i,j) + (- kappa_half_x * (value_dvx_dx + value_dvy_dy) + 
		#  0.5d0 * sum_of_memory_variables_kappa) * DELTAT
		# Here sumr should be at the time same as dvdx
		@inbounds p[iz,ix] = pp[iz,ix] + ((mod[iz,ix] * (dvdx[iz,ix] + dvdz[iz,ix])) + (sumr)) * fc1  
	end
	end
end


# update for ix=[3,...,nx-2]
# ix=1,2 are virtual nodes (apply rigid boundary condition later)
# ix=nx-1,nx are virtual nodes (apply rigid boundary condition later)
@inbounds @fastmath function update_dpdx!(pap,pac,attrib_mod)#
	dpdx=pap.w2[:dx][:p]
	p=pap.w2[:tp][:p]
	mp=pap.memory_pml[:dpdx]
	dpdx!(dpdx,p,mp,pac.pml[:x][:b_half],pac.pml[:x][:a_half],pac.pml[:x][:k_halfI],pac.fc[:δx24I],pac.ic[:nz],pac.ic[:nx])
end
@inbounds @fastmath function dpdx!(dpdx,p,mp,pml1,pml2,pml3,fc1,nz,nx)
	for ix=2:nx-2
	@simd for iz=1:nz
		@inbounds dpdx[iz,ix] = (27.e0*p[iz,ix+1]-27.e0*p[iz,ix]-p[iz,ix+2]+p[iz,ix-1]) * fc1
		@inbounds mp[iz,ix] = pml1[ix] * mp[iz,ix] + pml2[ix] * dpdx[iz,ix] # pml
		@inbounds dpdx[iz,ix] = muladd(dpdx[iz,ix] , pml3[ix] , mp[iz,ix]) # pml
	end
	end
end

# update for iz=[2,...,nz-2]
# iz=1,2 are virtual node
# iz=nz-1,nz are virtual nodes
@inbounds @fastmath function update_dpdz!(pap,pac,attrib_mod)#
	dpdz=pap.w2[:dz][:p]
	p=pap.w2[:tp][:p]
	mp=pap.memory_pml[:dpdz]
	dpdz!(dpdz,p,mp,pac.pml[:z][:b_half],pac.pml[:z][:a_half],pac.pml[:z][:k_halfI],pac.fc[:δz24I],pac.ic[:nz],pac.ic[:nx])
end
@inbounds @fastmath function dpdz!(dpdz,p,mp,pml1,pml2,pml3,fc1,nz,nx)
	for ix=1:nx
	@simd for iz=2:nz-2
		@inbounds dpdz[iz,ix] = (27.e0*p[iz+1,ix]-27.e0*p[iz,ix]-p[iz+2,ix]+p[iz-1,ix]) * fc1
		@inbounds mp[iz,ix] = pml1[iz] * mp[iz,ix] + pml2[iz] * dpdz[iz,ix] # pml
		@inbounds dpdz[iz,ix] = muladd(dpdz[iz,ix] , pml3[iz] , mp[iz,ix]) # pml
	end
	end
end

# dpdx is previously computed from [2,...,nx-2]
@inbounds @fastmath function update_vx!(pap,pac,attrib_mod)#
	dpdx=pap.w2[:dx][:p]
	vx=pap.w2[:t][:vx]
	vxp=pap.w2[:tp][:vx]
	vx!(vx,vxp,dpdx,pac.mod[:rhovxI],pac.fc[:δt],pac.ic[:nz],pac.ic[:nx])
end
@inbounds @fastmath function vx!(vx,vxp,dpdx,mod,fc1,nz,nx)
	for ix=1:nx # see dpdxw computation above
	@simd for iz=1:nz
		@inbounds vx[iz,ix] = vxp[iz,ix] + (dpdx[iz,ix]) * fc1 * mod[iz,ix] #* boundary_vx(iz,ix)
	end
	end
end

# dpdz is previously computed from [2,...,nz-2]
@inbounds @fastmath function update_vz!(pap,pac,attrib_mod)#
	dpdz=pap.w2[:dz][:p]
	vz=pap.w2[:t][:vz]
	vzp=pap.w2[:tp][:vz]
	vz!(vz,vzp,dpdz,pac.mod[:rhovzI],pac.fc[:δt],pac.ic[:nz],pac.ic[:nx])
end
@inbounds @fastmath function vz!(vz,vzp,dpdz,mod,fc1,nz,nx)
	for ix=1:nx
	@simd for iz=1:nz # see dpdzw computation above
		@inbounds vz[iz,ix] = vzp[iz,ix] + (dpdz[iz,ix]) * fc1 * mod[iz,ix] #* boundary_vz(iz,ix)
	end
	end
end





# ====================================
# not used at the moment
# ====================================
# rigid boundary conditions on the left side of the model
function rigid_boundary_left!(pw,nx,nz)
	@simd for iz=1:nz
		pw[iz,1,2] = -1.0*pw[iz,4,2]  # ix=1 is the virtual node for vx
		pw[iz,2,2] = -1.0*pw[iz,3,2]  # ix=2 is the virtual node for vx
		pw[iz,3,3] = 0.0 # ix=3, is where the vz should be zero
	end
end

function rigid_boundary_right!(pw,nx,nz)
	@simd for iz=1:nz
		pw[iz,nx-1,2] = -1.0*pw[iz,nx-2,2]  # ix=nx-1 is the virtual node for vx
		pw[iz,nx,2] = -1.0*pw[iz,nx-3,2]  # ix=nx is the virtual node for vx
		pw[iz,nx-1,3] = 0.0 # ix=nx-1, is where the vz should be zero
	end
end

# rigid boundary conditions on the top side of the model
function rigid_boundary_top!(pw,nx,nz)
	@simd for ix=1:nx
		pw[1,ix,3] = -1.0*pw[4,ix,3]  # iz=1 is the virtual node for vz
		pw[2,ix,3] = -1.0*pw[3,ix,3]  # iz=2 is the virtual node for vz
		pw[3,ix,2] = 0.0 # iz=3 is where the vx should be zero
	end
end

function rigid_boundary_bottom!(pw,nx,nz)
	@simd for ix=1:nx
		pw[nz-1,ix,3] = -1.0*pw[nz-2,ix,3]  # iz=nz-1 is the virtual node for vz
		pw[nz,ix,3] = -1.0*pw[nz-3,ix,3]  # iz=nz-1 is the virtual node for vz
		pw[nz-1,ix,2] = 0.0 # iz=nz-1, is where the vx should be zero
	end
end



function Δtsample(p)
	p[1] .= 0.0  # it is necessary to have these boundary conditions
	nt=length(p)
	dpdt=zero(p)
	pout=zero(p)
		for it=1:nt-1 # use these bounds later
			dpdt[it] = p[it+1]-p[it]
		end

		for it=1:nt-1 # use these exact bounds later
			pout[it+1] = dpdt[it+1]-dpdt[it]
		end
	return pout
end

"""
Simple test to see staggerred grid spatial derivatives that are used satisfy adjoint test.
"""
function Δxsample(p)	
	nx=length(p)
	# boundary conditions
	p[1:2].=0.0
	p[nx-1:nx].=0.0
	dpdx=zero(p)
	pout=zero(p)
	for i in 1:100 # applying operators a 100 times, say
		for ix=2:nx-2 # use these bounds later
			dpdx[ix] = (27.e0*p[ix+1]-27.e0*p[ix]-p[ix+2]+p[ix-1]) * 1.0
		end

		for ix=3:nx-1 # use these exact bounds later
			pout[ix] = (27.e0*dpdx[ix]-27.e0*dpdx[ix-1]-dpdx[ix+1]+dpdx[ix-2]) * 1.0
		end
	end
	
	return pout
end


"""
advance by one time step for adjoint test
"""
function advance_sample(ntimes, p, modrhovxI, modrhovzI, modK)
	nz=size(p,1)
	nx=size(p,2)
	# boundary conditions
	for i in 1:3
		p[1:2,:,i] .= 0.0
		p[:,1:2,i] .= 0.0
		p[nz-1:nz,:,i] .= 0.0
		p[:,nx-1:nx,i] .= 0.0
	end


	δx24I=1.; δz24I=1.; δt=1.

	a_x=zeros(nx)
	a_x_half=zeros(nx)
	b_x=zeros(nx);	b_x_half=zeros(nx)
	k_x=ones(nx);	k_xI=ones(nx)
	k_x_half=ones(nx);	k_x_halfI=ones(nx)

	k_z=ones(nz);	k_zI=ones(nz)
	k_z_half=ones(nz);	k_z_halfI=ones(nz)
	a_z=zeros(nz);	a_z_half=zeros(nz)
	b_z=zeros(nz);	b_z_half=zeros(nz)

	modrhovxI=ones(nz,nx)
	modrhovzI=ones(nz,nx)
	modK=ones(nz,nx)

	pw=zeros(nz,nx,3)
	#input
	copyto!(pw, p)
#=	
	for iz in 1:nz
		for ix in 1:nx
			pw[iz,ix,3]=p[iz,ix,3]
		end
	end
=#
	dpdxw=zeros(nz,nx,3); 	dpdzw=zeros(nz,nx,3) 
	ppw=zeros(nz,nx,3); 	pppw=zeros(nz,nx,3) 

	memory_dvx_dxw=zeros(nz,nx)
	memory_dvx_dzw=zeros(nz,nx)
	memory_dvz_dxw=zeros(nz,nx)
	memory_dvz_dzw=zeros(nz,nx)
	memory_dp_dxw=zeros(nz,nx)
	memory_dp_dzw=zeros(nz,nx)

	for i in 1:ntimes

		advance_kernel!(pw, dpdxw, dpdzw, δx24I, δz24I,
			 memory_dp_dxw,memory_dp_dzw,
			 b_x_half,b_z_half,a_x_half,a_z_half,k_x_halfI,k_z_halfI,
			 nx,nz,δt,modrhovxI,modrhovzI,memory_dvx_dxw,memory_dvz_dzw,
			 b_x,b_z,a_x,a_z,k_xI,k_zI,modK)
	end


	return pw
end


