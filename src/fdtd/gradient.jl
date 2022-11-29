# This routine ABSOLUTELY should not allocate any memory, called inside time loop.
@inbounds @fastmath function compute_gradient!(issp::Int64, pac, pap)
	# aliases
	p1=pap[1].w1[:t][:p]
	p1p=pap[1].w1[:tp][:p]
	p1pp=pap[1].w1[:tpp][:p]
	p2p=pap[2].w1[:tp][:p]
	fc1=pac.fc[:δtI]

	dpdx1=pap[1].w1[:dx][:p]
	dpdx2=pap[2].w1[:dx][:p]
	dpdz1=pap[1].w1[:dz][:p]
	dpdz2=pap[2].w1[:dz][:p]

	gKI=pap[1].ss[issp].grad_mod[:KI]
	grhovxI=pap[1].ss[issp].grad_mod[:rhovxI]
	grhovzI=pap[1].ss[issp].grad_mod[:rhovzI]

	gmodKI!(gKI,p1,p1p,p1pp,p2p,pac.ic[:nx],pac.ic[:nz],fc1)
	gmodrhovxI!(grhovxI,dpdx1,dpdx2,pac.ic[:nx],pac.ic[:nz])
	gmodrhovzI!(grhovzI,dpdz1,dpdz2,pac.ic[:nx],pac.ic[:nz])
end

@inbounds @fastmath function gmodKI!(gKI,p1,p1p,p1pp,p2p,nx,nz,δtI)
	for ix=1:nx
		@simd for iz=1:nz
			# p at [it], pp at [it-1]	# dpdx and dpdz at [it]	# gradients w.r.t inverse of modKI, i.e., 1/rho/c^2 
			@inbounds gKI[iz,ix] += ((-1.0 * (p1pp[iz, ix] + p1[iz, ix] - 2.0 * p1p[iz, ix]) * δtI * δtI) *  p2p[iz,ix])
		end
	end
end


@inbounds @fastmath function gmodrhovxI!(grhovxI,dpdx1,dpdx2,nx,nz)
	for ix=1:nx
		@simd for iz=1:nz
			# gradient w.r.t. inverse of rho on vx and vz grids
			@inbounds grhovxI[iz,ix] += (- dpdx2[iz,ix,1]*dpdx1[iz,ix,1])
		end
	end

end
@inbounds @fastmath function gmodrhovzI!(grhovzI,dpdz1,dpdz2,nx,nz)
	for ix=1:nx
		@simd for iz=1:nz
			# gradient w.r.t. inverse of rho on vx and vz grids
			@inbounds grhovzI[iz,ix] += (- dpdz2[iz,ix,1]*dpdz1[iz,ix,1])
		end
	end
end

@inbounds @fastmath function scale_gradient!(issp::Int64,pap,δ)
	gKI=pap[1].ss[issp].grad_mod[:KI]
	grhovxI=pap[1].ss[issp].grad_mod[:rhovxI]
	grhovzI=pap[1].ss[issp].grad_mod[:rhovzI]
	"gradient is formed by intergration over time, hence multiply with δt, but why not?"
	"I don't completely understand where the factors δx and δz are coming from..."
	"probably the source term should not be multiplied by δxI and δzI during adjoint propagation"
	rmul!(gKI,δ)
	rmul!(grhovxI,δ)
	rmul!(grhovzI,δ)
end




function grad_modrr!(pac)
	grad_modrr_sprayrr!(pac.grad_mod[:rhoI],pac.grad_mod[:rhovxI],pac.grad_mod[:rhovzI])
	grad_modrr_sprayrhovxI!(pac.grad_mod[:rhoI],pac.grad_mod[:rhovxI])
	grad_modrr_sprayrhovzI!(pac.grad_mod[:rhoI],pac.grad_mod[:rhovzI])
end

function stack_grads!(pac, pap)
	# theses are SharedArrays
	gmodKI=pac.grad_mod[:KI]
	gmodrhovxI=pac.grad_mod[:rhovxI]
	gmodrhovzI=pac.grad_mod[:rhovzI]
	for issp in 1:length(pap[1].ss)
		gs=pap[1].ss[issp].grad_mod[:KI]
		for i in eachindex(gmodKI)
			gmodKI[i]+=gs[i]  # only update gmodKI
		end
		gs=pap[1].ss[issp].grad_mod[:rhovxI]
		for i in eachindex(gmodrhovxI)
			gmodrhovxI[i]+=gs[i]
		end
		gs=pap[1].ss[issp].grad_mod[:rhovzI]
		for i in eachindex(gmodrhovzI)
			gmodrhovzI[i]+=gs[i]
		end
	end
end

function update_gradient!(pac)
	nx, nz=pac.ic[:nx], pac.ic[:nz]
	nznxd = prod(length.(pac.model.mgrid))

	# combine rhovxI and rhovzI
	grad_modrr!(pac)

	gradient=pac.gradient
	# truncate
	gmodKI=view(pac.grad_mod[:KI],_fd.npml+1:nz-_fd.npml,_fd.npml+1:nx-_fd.npml)
	gmodrr=view(pac.grad_mod[:rhoI],_fd.npml+1:nz-_fd.npml,_fd.npml+1:nx-_fd.npml)
	for i in 1:nznxd
		# parameterization is  [:KI, :ρI, :null]
		gradient[i]=gmodKI[i]  # update gmodKI
		gradient[nznxd+i]=gmodrr[i]  # update gmodrr
	end
end

