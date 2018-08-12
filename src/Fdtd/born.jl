

function add_born_sources!(issp::Int64, pac::Paramc, pass::Vector{Paramss}, pap::Paramp)

	born_svalue_stack=pass[issp].born_svalue_stack
	δx24I=pac.δx24I; δz24I=pac.δz24I; 
	δxI=pac.δxI; δzI=pac.δzI; 
	δt=pac.δt
	δtI=pac.δtI
	δmodtt=pac.δmodtt; modttI=pac.modttI;
	δmodrrvx=pac.δmodrrvx; δmodrrvz=pac.δmodrrvz

	nx=pac.model.mgrid.nx; nz=pac.model.mgrid.nz
	np=pac.model.mgrid.npml # born sources are not added in the PML region

	p=pap.p; pp=pap.pp; ppp=pap.ppp;
	dpdx=pap.dpdx; dpdz=pap.dpdz;

	# secondary sources for Born modeling
	# adding born sources from pressure(:,:,1) to pressure(:,:,2)
	# upto until [it-2]
	# lambdaI scatterrer source term at [it-1]
	# p is at [it], pp is at [it-1], ppp is at [it-2]
	# dpdx is at [it-1] and dpdz is at [it-1]
	# modrrvx scatterrer source term at [it-1]
	# modrrvz scatterrer source term at [it-1]

	born_stacktt!(born_svalue_stack,p,pp,ppp,δmodtt,nx,nz,np,δtI,δt)
	born_stackrrvx!(born_svalue_stack,dpdx,δmodrrvx,nx,nz,np,δx24I,δt)
	born_stackrrvz!(born_svalue_stack,dpdz,δmodrrvz,nx,nz,np,δz24I,δt)
	born_stack!(p,born_svalue_stack,modttI,nx,nz,np,δt)
end
@inbounds @fastmath function born_stacktt!(born_svalue_stack,p,pp,ppp,δmodtt,nx,nz,np,δtI,δt)
	pppw=ppp[1]
	ppw=pp[1]
	pw=p[1]
	for ix=np+1:np+nx
		@simd for iz=np+1:np+nz
			born_svalue_stack[iz,ix] += 
			δt * ((-1.0 * (pppw[iz, ix, 1] + pw[iz, ix,  1] - 2.0 * ppw[iz, ix,  1]) * δmodtt[iz, ix] * δtI * δtI)) 
		end
	end
end
@inbounds @fastmath function born_stackrrvx!(born_svalue_stack,dpdx,δmodrrvx,nx,nz,np,δx24I,δt)
	dpdxw=dpdx[1]
	for ix=np+1:np+nx
		@simd for iz=np+1:np+nz
			born_svalue_stack[iz,ix] += 
				δt * ((27.e0*dpdxw[iz,ix,1] * δmodrrvx[iz,ix] -27.e0*dpdxw[iz,ix-1,1] * δmodrrvx[iz,ix-1] -dpdxw[iz,ix+1,1] * δmodrrvx[iz,ix+1] +dpdxw[iz,ix-2,1] * δmodrrvx[iz,ix-2] ) * (δx24I)) 
		end
	end

end
@inbounds @fastmath function born_stackrrvz!(born_svalue_stack,dpdz,δmodrrvz,nx,nz,np,δz24I,δt)
	dpdzw=dpdz[1]
	for ix=np+1:np+nx
		@simd for iz=np+1:np+nz
			born_svalue_stack[iz,ix] += 
				δt * ((27.e0*dpdzw[iz,ix,1] * δmodrrvz[iz,ix] -27.e0*dpdzw[iz-1,ix,1] * δmodrrvz[iz-1,ix] -dpdzw[iz+1,ix,1] * δmodrrvz[iz+1,ix] +dpdzw[iz-2,ix,1] * δmodrrvz[iz-2,ix] ) * (δz24I))  

		end
	end
end
@inbounds @fastmath function born_stack!(p,born_svalue_stack,modttI,nx,nz,np,δt)
	pw2=p[2]
	for ix=np+1:np+nx
		@simd for iz=np+1:np+nz
			pw2[iz,ix,1] += born_svalue_stack[iz,ix] * modttI[iz,ix] * δt  #* δxI * δzI 
		end
	end
end

