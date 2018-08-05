


@inbounds @fastmath function advance!(pac::Paramc, pap::Paramp)
	# aliases
	p=pap.p; pp=pap.pp; ppp=pap.ppp;
	dpdx=pap.dpdx; dpdz=pap.dpdz;
	memory_dp_dx=pap.memory_dp_dx; memory_dp_dz=pap.memory_dp_dz; memory_dvx_dx=pap.memory_dvx_dx; memory_dvz_dz=pap.memory_dvz_dz
	modttI=pac.modttI; modrrvx=pac.modrrvx; modrrvz=pac.modrrvz
	δx24I=pac.δx24I; δz24I=pac.δz24I; δt=pac.δt
	nx=pac.nx; nz=pac.nz
	a_x=pac.a_x; b_x=pac.b_x; k_xI=pac.k_xI; a_x_half=pac.a_x_half; b_x_half=pac.b_x_half; k_x_halfI=pac.k_x_halfI 
	a_z=pac.a_z; b_z=pac.b_z; k_zI=pac.k_zI; a_z_half=pac.a_z_half; b_z_half=pac.b_z_half; k_z_halfI=pac.k_z_halfI

	pppppp!(p,pp,ppp,pac.activepw)

	"""
	compute dpdx and dpdz at [it-1] for all propagating fields
	"""
	update_dpdx!(p, dpdx, δx24I, memory_dp_dx, b_x_half, a_x_half, k_x_halfI, nx, nz,pac.activepw)
	update_dpdz!(p, dpdz, δz24I, memory_dp_dz, b_z_half, a_z_half, k_z_halfI, nx, nz,pac.activepw)

	"""
	update velocity at [it-1/2] using 
	velocity at [it-3/2] and dpdx and dpdz at [it-1] 
	"""
	update_vx!(p, dpdx, δt, modrrvx, nx, nz,pac.activepw)
	update_vz!(p, dpdz, δt, modrrvz, nx, nz,pac.activepw)


	dvdx!(dpdx,p,memory_dvx_dx,b_x,a_x,k_xI,nz,nx,δx24I,pac.activepw)
	dvdz!(dpdz,p,memory_dvz_dz,b_z,a_z,k_zI,nz,nx,δz24I,pac.activepw)

	"""
	compute pressure at [it] using p at [it-1] and dvxdx
	and dvzdz at [it-1/2]
	"""
	pvzvx!(p,dpdx,dpdz,modttI,nz,nx,δt,pac.activepw)

end
function pppppp!(p,pp,ppp,activepw)
	for ipw in activepw 
		pw=p[ipw]
		ppw=pp[ipw]
		pppw=ppp[ipw]
		copy!(pppw,ppw)
		copy!(ppw,pw)
	end
end

@inbounds @fastmath function dvdx!(dpdx,p,memory_dvx_dx,b_x,a_x,k_xI,nz,nx,δx24I,activepw)
	for ipw in activepw
		pw=p[ipw]
		dpdxw=dpdx[ipw]
		memory_dvx_dxw=memory_dvx_dx[ipw]
		for ix=3:nx-2
		@simd for iz=3:nz-2
			@inbounds dpdxw[iz,ix,2] = (27.e0*pw[iz,ix,2]-27.e0*pw[iz,ix-1,2]-pw[iz,ix+1,2]+pw[iz,ix-2,2]) * (δx24I)
			@inbounds memory_dvx_dxw[iz,ix] = b_x[ix] * memory_dvx_dxw[iz,ix] + a_x[ix] * dpdxw[iz,ix,2] # pml 
			@inbounds dpdxw[iz,ix,2] = dpdxw[iz,ix,2] * k_xI[ix] + memory_dvx_dxw[iz,ix] # pml
		end
		end
	end
end

@inbounds @fastmath function dvdz!(dpdz,p,memory_dvz_dz,b_z,a_z,k_zI,nz,nx,δz24I,activepw)
	for ipw in activepw
		pw=p[ipw]
		dpdzw=dpdz[ipw]
		memory_dvz_dzw=memory_dvz_dz[ipw]
		for ix=3:nx-2
		@simd for iz=3:nz-2
			@inbounds dpdzw[iz,ix,3] = (27.e0*pw[iz,ix,3]-27.e0*pw[iz-1,ix,3]-pw[iz+1,ix,3]+pw[iz-2,ix,3]) * (δz24I)
			@inbounds memory_dvz_dzw[iz,ix] = b_z[iz] * memory_dvz_dzw[iz,ix] + a_z[iz] * dpdzw[iz,ix,3] # pml
			@inbounds dpdzw[iz,ix,3] = dpdzw[iz,ix,3] * k_zI[iz] + memory_dvz_dzw[iz,ix] # pml
		end
		end
	end
end

@inbounds @fastmath function pvzvx!(p,dpdx,dpdz,modttI,nz,nx,δt,activepw)
	for ipw in activepw
		pw=p[ipw]
		dpdxw=dpdx[ipw]
		dpdzw=dpdz[ipw]
		for ix=3:nx-2
		@simd for iz=3:nz-2
			@inbounds pw[iz,ix,1] += (modttI[iz,ix] * (dpdxw[iz,ix,2] + dpdzw[iz,ix,3])) * δt #* boundary_p(iz,ix)
		end
		end
	end
end

@inbounds @fastmath function update_dpdx!(p, dpdx, δx24I, memory_dp_dx, b_x_half, a_x_half, k_x_halfI, nx, nz,activepw)
	for ipw in activepw
		pw=p[ipw]
		dpdxw=dpdx[ipw]
		memory_dp_dxw=memory_dp_dx[ipw]
		for ix = 3:nx-2
		@simd for iz = 3:nz-2
			@inbounds dpdxw[iz,ix,1] = (27.e0*pw[iz,ix+1,1]-27.e0*pw[iz,ix,1]-pw[iz,ix+2,1]+pw[iz,ix-1,1]) * (δx24I)
			@inbounds memory_dp_dxw[iz,ix] = b_x_half[ix] * memory_dp_dxw[iz,ix] + a_x_half[ix] * dpdxw[iz,ix,1] # pml
			@inbounds dpdxw[iz,ix,1] = dpdxw[iz,ix,1] * k_x_halfI[ix] + memory_dp_dxw[iz,ix] # pml
		end
		end
	end
end

@inbounds @fastmath function update_dpdz!(p, dpdz, δz24I, memory_dp_dz, b_z_half, a_z_half, k_z_halfI, nx, nz,activepw)
	for ipw in activepw
		pw=p[ipw]
		dpdzw=dpdz[ipw]
		memory_dp_dzw=memory_dp_dz[ipw]
		for ix = 3:nx-2
		@simd for iz = 3:nz-2
			@inbounds dpdzw[iz,ix,1] = (27.e0*pw[iz+1,ix,1]-27.e0*pw[iz,ix,1]-pw[iz+2,ix,1]+pw[iz-1,ix,1]) * (δz24I)
			@inbounds memory_dp_dzw[iz,ix] = b_z_half[iz] * memory_dp_dzw[iz,ix] + a_z_half[iz] * dpdzw[iz,ix,1] # pml
			@inbounds dpdzw[iz,ix,1] = dpdzw[iz,ix,1] * k_z_halfI[iz] + memory_dp_dzw[iz,ix] # pml
		end
		end
	end
end

@inbounds @fastmath function update_vx!(p, dpdx, δt, modrrvx,  nx, nz,activepw)
	for ipw in activepw
		pw=p[ipw]
		dpdxw=dpdx[ipw]
		for ix=3:nx-2
		@simd for iz=3:nz-2
			@inbounds pw[iz,ix,2] += (dpdxw[iz,ix,1]) * δt * modrrvx[iz,ix] #* boundary_vx(iz,ix)
		end
		end
	end
end

@inbounds @fastmath function update_vz!(p, dpdz, δt,  modrrvz, nx, nz,activepw)
	for ipw in activepw
		pw=p[ipw]
		dpdzw=dpdz[ipw]
		for ix=3:nx-2
		@simd for iz=3:nz-2
			@inbounds pw[iz,ix,3] +=  (dpdzw[iz,ix,1]) * δt * modrrvz[iz,ix] #* boundary_vz(iz,ix)
		end
		end
	end
end


