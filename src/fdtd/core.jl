@inbounds @fastmath function advance!(pac, pap)
	# aliases
	#=
	p=pap.p; pp=pap.pp; ppp=pap.ppp;
	dpdx=pap.dpdx; dpdz=pap.dpdz;
	memory_dp_dx=pap.memory_dp_dx; memory_dp_dz=pap.memory_dp_dz; 
	memory_dvx_dx=pap.memory_dvx_dx; memory_dvz_dz=pap.memory_dvz_dz
	modttI=pac.modttI; modrrvx=pac.modrrvx; modrrvz=pac.modrrvz
	δx24I=pac.δx24I; δz24I=pac.δz24I; δt=pac.δt
	nx=pac.nx; nz=nz
	a_x=pac.a_x; b_x=pac.b_x; k_xI=pac.k_xI; a_x_half=pac.a_x_half; b_x_half=pac.b_x_half; k_x_halfI=pac.k_x_halfI 
	a_z=pac.a_z; b_z=pac.b_z; k_zI=pac.k_zI; a_z_half=pac.a_z_half; b_z_half=pac.b_z_half; k_z_halfI=pac.k_z_halfI
	attrib_mod=pac.attrib_mod
	=#

	for ipw in pac.activepw
		w2=pap[ipw].w2
		mpml=pap[ipw].memory_pml
		#=
		pw=p[ipw]
		ppw=pp[ipw]
		pppw=ppp[ipw]
		dpdxw=dpdx[ipw]
		dpdzw=dpdz[ipw]
		memory_dp_dxw=memory_dp_dx[ipw]
		memory_dp_dzw=memory_dp_dz[ipw]
		memory_dvx_dxw=memory_dvx_dx[ipw]
		memory_dvz_dzw=memory_dvz_dz[ipw]
		=#

		# store p for the last two steps
		pppppp!(w2)

		advance_kernel!(w2, pac.ic, pac.fc, mpml, pac.pml,pac.mod)
		  #dpdxw, dpdzw, δx24I, δz24I,
		#	 memory_dp_dxw, memory_dp_dzw,
		#	 b_x_half,b_z_half,a_x_half,a_z_half,k_x_halfI,k_z_halfI,
		##	 nx,nz,δt,modrrvx,modrrvz,memory_dvx_dxw,memory_dvz_dzw,
	#		 b_x,b_z,a_x,a_z,k_xI,k_zI,modttI,attrib_mod)

	end
	
	return nothing
end

function advance_kernel!(w2,ic,fc,mpml,pml,mod)
	#compute dpdx and dpdz at [it-1] for all propagating fields
	update_dpdx!(w2,ic,fc,mpml,pml)
	update_dpdz!(w2,ic,fc,mpml,pml)

	#update velocity at [it-1/2] using 
	#velocity at [it-3/2] and dpdx and dpdz at [it-1] 
	update_vx!(w2,ic,fc,mod)
	update_vz!(w2,ic,fc,mod)

	#rigid boundary conditions (not currently implemented, need to do adjoint test)
	#such that:
	#p[3,.] is where velocities are zero and
	#p[.,n-1] is where velocities are zero
#	rigid_boundary_left!(p,nx,nz)
#	rigid_boundary_right!(p,nx,nz)
#	rigid_boundary_top!(p,nx,nz)
#	rigid_boundary_bottom!(p,nx,nz)


	update_dvdx!(w2,ic,fc,mpml,pml)
	update_dvdz!(w2,ic,fc,mpml,pml)

	#compute pressure at [it] using p at [it-1] and dvxdx
	#and dvzdz at [it-1/2]
	update_p!(w2,ic,fc,mod)

end


function pppppp!(w2)
	copyto!.(w2[:tpp],w2[:tp])
	copyto!.(w2[:tp],w2[:t])
end



# 
@inbounds @fastmath function update_dvdx!(w2,ic,fc,mpml,pml)#dpdxw,pw,memory_dvx_dxw,b_x,a_x,k_xI,nz,nx,δx24I)
	dvdx=w2[:dx][:vx]
	vx=w2[:t][:vx]
	mp=mpml[:dvxdx]
	dvdx!(dvdx,vx,mp,pml[:b_x],pml[:a_x],pml[:k_xI],fc[:δx24I],ic[:nz],ic[:nx])
end
@inbounds @fastmath function dvdx!(dvdx,vx,mp,pml1,pml2,pml3,fc1,nz,nx)
	for ix=3:nx-1 # ix=1,2,nx are useless, rigid boundary conditions apply on ix=3,nx-1
	@simd for iz=1:nz
		@inbounds dvdx[iz,ix] = (27.e0*vx[iz,ix]-27.e0*vx[iz,ix-1]-vx[iz,ix+1]+vx[iz,ix-2]) * fc1
		@inbounds mp[iz,ix] = pml1[ix] * mp[iz,ix] + pml2[ix] * dvdx[iz,ix] # pml 
		@inbounds dvdx[iz,ix] = dvdx[iz,ix] * pml3[ix] + mp[iz,ix] # pml
	end
	end
end


@inbounds @fastmath function update_dvdz!(w2,ic,fc,mpml,pml)#dpdzw,pw,memory_dvz_dzw,b_z,a_z,k_zI,nz,nx,δz24I)
	dvdz=w2[:dz][:vz]
	vz=w2[:t][:vz]
	mp=mpml[:dvzdz]
	dvdz!(dvdz,vz,mp,pml[:b_z],pml[:a_z],pml[:k_zI],fc[:δz24I],ic[:nz],ic[:nx])
end
@inbounds @fastmath function dvdz!(dvdz,vz,mp,pml1,pml2,pml3,fc1,nz,nx)
	for ix=1:nx
	@simd for iz=3:nz-1 # iz=1,2 are usless, boundary conditions on iz=3,nz-1
		@inbounds dvdz[iz,ix] = (27.e0*vz[iz,ix]-27.e0*vz[iz-1,ix]-vz[iz+1,ix]+vz[iz-2,ix]) * fc1
		@inbounds mp[iz,ix] = pml1[iz] * mp[iz,ix] + pml2[iz] * dvdz[iz,ix] # pml
		@inbounds dvdz[iz,ix] = dvdz[iz,ix] * pml3[iz] + mp[iz,ix] # pml
	end
	end
end

# no attenuation (no memory in stress-strain relation)
@inbounds @fastmath function update_p!(w2,ic,fc,mod)#pw,dpdxw,dpdzw,modttI,nz,nx,δt, ::T) where {T<:Union{Fdtd,FdtdBorn}}
	dvdx=w2[:dx][:vx]
	dvdz=w2[:dz][:vz]
	p=w2[:t][:p]
	pvzvx!(p,dvdz,dvdx,mod[:ttI],fc[:δt],ic[:nz],ic[:nx])
end
@inbounds @fastmath function pvzvx!(p,dvdz,dvdx,mod,fc1,nz,nx)
	for ix=1:nx  # see limits above
	@simd for iz=1:nz
		@inbounds p[iz,ix] += (mod[iz,ix] * (dvdx[iz,ix] + dvdz[iz,ix])) * fc1 #* boundary_p(iz,ix)
	end
	end
end


# update for ix=[3,...,nx-2]
# ix=1,2 are virtual nodes (apply rigid boundary condition later)
# ix=nx-1,nx are virtual nodes (apply rigid boundary condition later)
@inbounds @fastmath function update_dpdx!(w2,ic,fc,mpml,pml)#pw, dpdxw, δx24I, memory_dp_dxw, b_x_half, a_x_half, k_x_halfI, nx, nz)
	dpdx=w2[:dx][:p]
	p=w2[:t][:p]
	mp=mpml[:dpdx]
	dpdx!(dpdx,p,mp,pml[:b_x_half],pml[:a_x_half],pml[:k_x_halfI],fc[:δx24I],ic[:nz],ic[:nx])
end
@inbounds @fastmath function dpdx!(dpdx,p,mp,pml1,pml2,pml3,fc1,nz,nx)
	for ix=2:nx-2
	@simd for iz=1:nz
		@inbounds dpdx[iz,ix] = (27.e0*p[iz,ix+1]-27.e0*p[iz,ix]-p[iz,ix+2]+p[iz,ix-1]) * fc1
		@inbounds mp[iz,ix] = pml1[ix] * mp[iz,ix] + pml2[ix] * dpdx[iz,ix] # pml
		@inbounds dpdx[iz,ix] = dpdx[iz,ix] * pml3[ix] + mp[iz,ix] # pml
	end
	end
end

# update for iz=[2,...,nz-2]
# iz=1,2 are virtual node
# iz=nz-1,nz are virtual nodes
@inbounds @fastmath function update_dpdz!(w2,ic,fc,mpml,pml)#pw, dpdzw, δz24I, memory_dp_dzw, b_z_half, a_z_half, k_z_halfI, nx, nz)
	dpdz=w2[:dz][:p]
	p=w2[:t][:p]
	mp=mpml[:dpdz]
	dpdz!(dpdz,p,mp,pml[:b_z_half],pml[:a_z_half],pml[:k_z_halfI],fc[:δz24I],ic[:nz],ic[:nx])
end
@inbounds @fastmath function dpdz!(dpdz,p,mp,pml1,pml2,pml3,fc1,nz,nx)
	for ix=1:nx
	@simd for iz=2:nz-2
		@inbounds dpdz[iz,ix] = (27.e0*p[iz+1,ix]-27.e0*p[iz,ix]-p[iz+2,ix]+p[iz-1,ix]) * fc1
		@inbounds mp[iz,ix] = pml1[iz] * mp[iz,ix] + pml2[iz] * dpdz[iz,ix] # pml
		@inbounds dpdz[iz,ix] = dpdz[iz,ix] * pml3[iz] + mp[iz,ix] # pml
	end
	end
end

# dpdx is previously computed from [2,...,nx-2]
@inbounds @fastmath function update_vx!(w2,ic,fc,mod)#pw, dpdxw, δt, modrrvx,  nx, nz)
	dpdx=w2[:dx][:p]
	vx=w2[:t][:vx]
	vx!(vx,dpdx,mod[:rrvx],fc[:δt],ic[:nz],ic[:nx])
end
@inbounds @fastmath function vx!(vx,dpdx,mod,fc1,nz,nx)
	for ix=1:nx # see dpdxw computation above
	@simd for iz=1:nz
		@inbounds vx[iz,ix] += (dpdx[iz,ix]) * fc1 * mod[iz,ix] #* boundary_vx(iz,ix)
	end
	end
end

# dpdz is previously computed from [2,...,nz-2]
@inbounds @fastmath function update_vz!(w2,ic,fc,mod)#pw, dpdzw, δt, modrrvz, nx, nz)
	dpdz=w2[:dz][:p]
	vz=w2[:t][:vz]
	vz!(vz,dpdz,mod[:rrvz],fc[:δt],ic[:nz],ic[:nx])
end
@inbounds @fastmath function vz!(vz,dpdz,mod,fc1,nz,nx)
	for ix=1:nx
	@simd for iz=1:nz # see dpdzw computation above
		@inbounds vz[iz,ix] +=  (dpdz[iz,ix]) * fc1 * mod[iz,ix] #* boundary_vz(iz,ix)
	end
	end
end



#=
# viscoacoustic modeling (memory in stress-strain relation)
@inbounds @fastmath function pvzvx!(pa,ic,fc,pml)#pw,dpdxw,dpdzw,modttI,nz,nx,δt, ::FdtdVisco)
	for ix=1:nx  # see limits above
	@simd for iz=1:nz
		#=

		 # use the Auxiliary Differential Equation form, 
		 # which is second-order accurate in time if implemented 
		 #following eq (14) of Robertsson, Blanch and Symes, Geophysics, vol.
		 #	  59(9), pp 1444-1456 (1994), which is what we do here        
		 sum_of_memory_variables_kappa = 0.d0        
		 
		 # loop over standard linear solids
		 for i_sls = 1,N_SLS
			 # this average of the two terms comes from eq (14) of 
			 #Robertsson, Blanch and Symes, Geophysics, vol. 59(9), pp 1444-1456 (1994)          
			 memory_variable_R_dot(i,j,i_sls) = 
			 	(memory_variable_R_dot_old(i,j,i_sls) + (value_dvx_dx + value_dvy_dy) * 
      kappa_unrelaxed(i,j) * 
      DELTAT_delta_relaxed_over_tau_sigma_without_Kappa(i_sls) - memory_variable_R_dot_old(i,j,i_sls) * 
      HALF_DELTAT_over_tau_sigma_kappa(i_sls))  * multiplication_factor_tau_sigma_kappa(i_sls)          
				
	sum_of_memory_variables_kappa = sum_of_memory_variables_kappa + &       
				    memory_variable_R_dot(i,j,i_sls) + memory_variable_R_dot_old(i,j,i_sls)        
		 end
		 =#

		@inbounds pw[iz,ix,1] += (modttI[iz,ix] * (dpdxw[iz,ix,2] + dpdzw[iz,ix,3])) * δt #* boundary_p(iz,ix)

		  # this average of the two terms comes from eq (13) of 
		  #Robertsson, Blanch and Symes, Geophysics, vol. 59(9), pp 1444-1456 (1994)  
	        #pressure(i,j) = pressure(i,j) + (- kappa_half_x * (value_dvx_dx + value_dvy_dy) + 
		#  0.5d0 * sum_of_memory_variables_kappa) * DELTAT
	end
	end
end
=#




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
function advance_sample(ntimes, p, modrrvx, modrrvz, modttI)
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

	modrrvx=ones(nz,nx)
	modrrvz=ones(nz,nx)
	modttI=ones(nz,nx)

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
			 nx,nz,δt,modrrvx,modrrvz,memory_dvx_dxw,memory_dvz_dzw,
			 b_x,b_z,a_x,a_z,k_xI,k_zI,modttI)
	end


	return pw
end


