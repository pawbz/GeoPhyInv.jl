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
	p[1:2,:].=0.0
	p[:,1:2].=0.0
	p[nz-1:nz,:].=0.0
	p[:,nx-1:nx].=0.0
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
function advance_sample(ntimes, p, modrrvx, modrrvz, a_x=nothing)
	nz=size(p,1)
	nx=size(p,2)
	# boundary conditions
	p[1:2,:].=0.0
	p[:,1:2].=0.0
	p[nz-1:nz,:].=0.0
	p[:,nx-1:nx].=0.0


	δx24I=1.; δz24I=1.; δt=1.

	(a_x===nothing) && (a_x=zeros(nx))

	a_x_half=zeros(nx)
	b_x=zeros(nx);	b_x_half=zeros(nx)
	k_x=ones(nx);	k_xI=ones(nx)
	k_x_half=ones(nx);	k_x_halfI=ones(nx)

	k_z=ones(nz);	k_zI=ones(nz)
	k_z_half=ones(nz);	k_z_halfI=ones(nz)
	a_z=zeros(nz);	a_z_half=zeros(nz)
	b_z=zeros(nz);	b_z_half=zeros(nz)

	#modrrvx=ones(nz,nx)
	#modrrvz=ones(nz,nx)
	modttI=ones(nz,nx)

	pw=zeros(nz,nx,3)
	#input
	for iz in 1:nz
		for ix in 1:nx
			pw[iz,ix,1]=p[iz,ix]
		end
	end

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


	return pw[:,:,1]
end

@inbounds @fastmath function advance!(pac, pap)
	# aliases
	p=pap.p; pp=pap.pp; ppp=pap.ppp;
	dpdx=pap.dpdx; dpdz=pap.dpdz;
	memory_dp_dx=pap.memory_dp_dx; memory_dp_dz=pap.memory_dp_dz; 
	memory_dvx_dx=pap.memory_dvx_dx; memory_dvz_dz=pap.memory_dvz_dz
	modttI=pac.modttI; modrrvx=pac.modrrvx; modrrvz=pac.modrrvz
	δx24I=pac.δx24I; δz24I=pac.δz24I; δt=pac.δt
	nx=pac.nx; nz=pac.nz
	a_x=pac.a_x; b_x=pac.b_x; k_xI=pac.k_xI; a_x_half=pac.a_x_half; b_x_half=pac.b_x_half; k_x_halfI=pac.k_x_halfI 
	a_z=pac.a_z; b_z=pac.b_z; k_zI=pac.k_zI; a_z_half=pac.a_z_half; b_z_half=pac.b_z_half; k_z_halfI=pac.k_z_halfI

	for ipw in pac.activepw
		pw=p[ipw]
		ppw=pp[ipw]
		pppw=ppp[ipw]
		dpdxw=dpdx[ipw]
		dpdzw=dpdz[ipw]
		memory_dp_dxw=memory_dp_dx[ipw]
		memory_dp_dzw=memory_dp_dz[ipw]
		memory_dvx_dxw=memory_dvx_dx[ipw]
		memory_dvz_dzw=memory_dvz_dz[ipw]

		# store p for the last two steps
		pppppp!(pw,ppw,pppw)

		advance_kernel!(pw, dpdxw, dpdzw, δx24I, δz24I,
			 memory_dp_dxw, memory_dp_dzw,
			 b_x_half,b_z_half,a_x_half,a_z_half,k_x_halfI,k_z_halfI,
			 nx,nz,δt,modrrvx,modrrvz,memory_dvx_dxw,memory_dvz_dzw,
			 b_x,b_z,a_x,a_z,k_xI,k_zI,modttI)

	end
	
	return nothing
end

function advance_kernel!(pw, dpdxw, dpdzw, δx24I, δz24I,
			 memory_dp_dxw,memory_dp_dzw,
			 b_x_half,b_z_half,a_x_half,a_z_half,k_x_halfI,k_z_halfI,
			 nx,nz,δt,modrrvx,modrrvz,memory_dvx_dxw,memory_dvz_dzw,
			 b_x,b_z,a_x,a_z,k_xI,k_zI,modttI)
	#compute dpdx and dpdz at [it-1] for all propagating fields
	update_dpdx!(pw, dpdxw, δx24I, memory_dp_dxw, b_x_half, a_x_half, k_x_halfI, nx, nz)
	update_dpdz!(pw, dpdzw, δz24I, memory_dp_dzw, b_z_half, a_z_half, k_z_halfI, nx, nz)

	#update velocity at [it-1/2] using 
	#velocity at [it-3/2] and dpdx and dpdz at [it-1] 
	update_vx!(pw, dpdxw, δt, modrrvx, nx, nz)
	update_vz!(pw, dpdzw, δt, modrrvz, nx, nz)

	#rigid boundary conditions (not currently implemented, need to do adjoint test)
	#such that:
	#p[3,.] is where velocities are zero and
	#p[.,n-1] is where velocities are zero
#	rigid_boundary_left!(p,nx,nz)
#	rigid_boundary_right!(p,nx,nz)
#	rigid_boundary_top!(p,nx,nz)
#	rigid_boundary_bottom!(p,nx,nz)


	dvdx!(dpdxw,pw,memory_dvx_dxw,b_x,a_x,k_xI,nz,nx,δx24I)
	dvdz!(dpdzw,pw,memory_dvz_dzw,b_z,a_z,k_zI,nz,nx,δz24I)

	#compute pressure at [it] using p at [it-1] and dvxdx
	#and dvzdz at [it-1/2]
	pvzvx!(pw,dpdxw,dpdzw,modttI,nz,nx,δt)

end


function pppppp!(pw,ppw,pppw)
	copyto!(pppw,ppw)
	copyto!(ppw,pw)
end


# 
@inbounds @fastmath function dvdx!(dpdxw,pw,memory_dvx_dxw,b_x,a_x,k_xI,nz,nx,δx24I)
	for ix=3:nx-1 # ix=1,2,nx are useless, rigid boundary conditions apply on ix=3,nx-1
	@simd for iz=1:nz
		@inbounds dpdxw[iz,ix,2] = (27.e0*pw[iz,ix,2]-27.e0*pw[iz,ix-1,2]-pw[iz,ix+1,2]+pw[iz,ix-2,2]) * (δx24I)
		@inbounds memory_dvx_dxw[iz,ix] = b_x[ix] * memory_dvx_dxw[iz,ix] + a_x[ix] * dpdxw[iz,ix,2] # pml 
		@inbounds dpdxw[iz,ix,2] = dpdxw[iz,ix,2] * k_xI[ix] + memory_dvx_dxw[iz,ix] # pml
	end
	end
end


@inbounds @fastmath function dvdz!(dpdzw,pw,memory_dvz_dzw,b_z,a_z,k_zI,nz,nx,δz24I)
	for ix=1:nx
	@simd for iz=3:nz-1 # iz=1,2 are usless, boundary conditions on iz=3,nz-1
		@inbounds dpdzw[iz,ix,3] = (27.e0*pw[iz,ix,3]-27.e0*pw[iz-1,ix,3]-pw[iz+1,ix,3]+pw[iz-2,ix,3]) * (δz24I)
		@inbounds memory_dvz_dzw[iz,ix] = b_z[iz] * memory_dvz_dzw[iz,ix] + a_z[iz] * dpdzw[iz,ix,3] # pml
		@inbounds dpdzw[iz,ix,3] = dpdzw[iz,ix,3] * k_zI[iz] + memory_dvz_dzw[iz,ix] # pml
	end
	end
end

@inbounds @fastmath function pvzvx!(pw,dpdxw,dpdzw,modttI,nz,nx,δt)
	for ix=1:nx  # see limits above
	@simd for iz=1:nz
		@inbounds pw[iz,ix,1] += (modttI[iz,ix] * (dpdxw[iz,ix,2] + dpdzw[iz,ix,3])) * δt #* boundary_p(iz,ix)
	end
	end
end

# update for ix=[3,...,nx-2]
# ix=1,2 are virtual nodes (apply rigid boundary condition later)
# ix=nx-1,nx are virtual nodes (apply rigid boundary condition later)
@inbounds @fastmath function update_dpdx!(pw, dpdxw, δx24I, memory_dp_dxw, b_x_half, a_x_half, k_x_halfI, nx, nz)
	for ix=2:nx-2
	@simd for iz=1:nz
		@inbounds dpdxw[iz,ix,1] = (27.e0*pw[iz,ix+1,1]-27.e0*pw[iz,ix,1]-pw[iz,ix+2,1]+pw[iz,ix-1,1]) * (δx24I)
		@inbounds memory_dp_dxw[iz,ix] = b_x_half[ix] * memory_dp_dxw[iz,ix] + a_x_half[ix] * dpdxw[iz,ix,1] # pml
		@inbounds dpdxw[iz,ix,1] = dpdxw[iz,ix,1] * k_x_halfI[ix] + memory_dp_dxw[iz,ix] # pml
	end
	end
end

# update for iz=[2,...,nz-2]
# iz=1,2 are virtual node
# iz=nz-1,nz are virtual nodes
@inbounds @fastmath function update_dpdz!(pw, dpdzw, δz24I, memory_dp_dzw, b_z_half, a_z_half, k_z_halfI, nx, nz)
	for ix=1:nx
	@simd for iz=2:nz-2
		@inbounds dpdzw[iz,ix,1] = (27.e0*pw[iz+1,ix,1]-27.e0*pw[iz,ix,1]-pw[iz+2,ix,1]+pw[iz-1,ix,1]) * (δz24I)
		@inbounds memory_dp_dzw[iz,ix] = b_z_half[iz] * memory_dp_dzw[iz,ix] + a_z_half[iz] * dpdzw[iz,ix,1] # pml
		@inbounds dpdzw[iz,ix,1] = dpdzw[iz,ix,1] * k_z_halfI[iz] + memory_dp_dzw[iz,ix] # pml
	end
	end
end

# dpdx is previously computed from [2,...,nx-2]
@inbounds @fastmath function update_vx!(pw, dpdxw, δt, modrrvx,  nx, nz)
	for ix=1:nx # see dpdxw computation above
	@simd for iz=1:nz
		@inbounds pw[iz,ix,2] += (dpdxw[iz,ix,1]) * δt * modrrvx[iz,ix] #* boundary_vx(iz,ix)
	end
	end
end

# dpdz is previously computed from [2,...,nz-2]
@inbounds @fastmath function update_vz!(pw, dpdzw, δt, modrrvz, nx, nz)
	for ix=1:nx
	@simd for iz=1:nz # see dpdzw computation above
		@inbounds pw[iz,ix,3] +=  (dpdzw[iz,ix,1]) * δt * modrrvz[iz,ix] #* boundary_vz(iz,ix)
	end
	end
end

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




