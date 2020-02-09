# This routine ABSOLUTELY should not allocate any memory, called inside time loop.
@inbounds @fastmath function compute_gradient!(issp::Int64, pac, pap)
	# aliases
	p1=pap[1].w2[:t][:p]
	p1p=pap[1].w2[:tp][:p]
	p1pp=pap[1].w2[:tpp][:p]
	p2p=pap[2].w2[:tp][:p]
	fc1=pac.fc[:δtI]

	dpdx1=pap[1].w2[:dx][:p]
	dpdx2=pap[2].w2[:dx][:p]
	dpdz1=pap[1].w2[:dz][:p]
	dpdz2=pap[2].w2[:dz][:p]

	gtt=pap[1].ss[issp].grad_mod[:tt]
	grrvx=pap[1].ss[issp].grad_mod[:rrvx]
	grrvz=pap[1].ss[issp].grad_mod[:rrvz]

	gmodtt!(gtt,p1,p1p,p1pp,p2p,pac.ic[:nx],pac.ic[:nz],fc1)
	gmodrrvx!(grrvx,dpdx1,dpdx2,pac.ic[:nx],pac.ic[:nz])
	gmodrrvz!(grrvz,dpdz1,dpdz2,pac.ic[:nx],pac.ic[:nz])
end

@inbounds @fastmath function gmodtt!(gtt,p1,p1p,p1pp,p2p,nx,nz,δtI)
	for ix=1:nx
		@simd for iz=1:nz
			# p at [it], pp at [it-1]	# dpdx and dpdz at [it]	# gradients w.r.t inverse of modtt, i.e., 1/rho/c^2 
			@inbounds gtt[iz,ix] += ((-1.0 * (p1pp[iz, ix] + p1[iz, ix] - 2.0 * p1p[iz, ix]) * δtI * δtI) *  p2p[iz,ix])
		end
	end
end


@inbounds @fastmath function gmodrrvx!(grrvx,dpdx1,dpdx2,nx,nz)
	for ix=1:nx
		@simd for iz=1:nz
			# gradient w.r.t. inverse of rho on vx and vz grids
			@inbounds grrvx[iz,ix] += (- dpdx2[iz,ix,1]*dpdx1[iz,ix,1])
		end
	end

end
@inbounds @fastmath function gmodrrvz!(grrvz,dpdz1,dpdz2,nx,nz)
	for ix=1:nx
		@simd for iz=1:nz
			# gradient w.r.t. inverse of rho on vx and vz grids
			@inbounds grrvz[iz,ix] += (- dpdz2[iz,ix,1]*dpdz1[iz,ix,1])
		end
	end
end

@inbounds @fastmath function scale_gradient!(issp::Int64,pap,δ)
	gtt=pap[1].ss[issp].grad_mod[:tt]
	grrvx=pap[1].ss[issp].grad_mod[:rrvx]
	grrvz=pap[1].ss[issp].grad_mod[:rrvz]
	"gradient is formed by intergration over time, hence multiply with δt, but why not?"
	"I don't completely understand where the factors δx and δz are coming from..."
	"probably the source term should not be multiplied by δxI and δzI during adjoint propagation"
	rmul!(gtt,δ)
	rmul!(grrvx,δ)
	rmul!(grrvz,δ)
end




function grad_modrr!(pac)
	grad_modrr_sprayrr!(pac.grad_mod[:rr],pac.grad_mod[:rrvx],pac.grad_mod[:rrvz])
	grad_modrr_sprayrrvx!(pac.grad_mod[:rr],pac.grad_mod[:rrvx])
	grad_modrr_sprayrrvz!(pac.grad_mod[:rr],pac.grad_mod[:rrvz])
end

function stack_grads!(pac, pap)
	# theses are SharedArrays
	gmodtt=pac.grad_mod[:tt]
	gmodrrvx=pac.grad_mod[:rrvx]
	gmodrrvz=pac.grad_mod[:rrvz]
	for issp in 1:length(pap[1].ss)
		gs=pap[1].ss[issp].grad_mod[:tt]
		for i in eachindex(gmodtt)
			gmodtt[i]+=gs[i]  # only update gmodtt
		end
		gs=pap[1].ss[issp].grad_mod[:rrvx]
		for i in eachindex(gmodrrvx)
			gmodrrvx[i]+=gs[i]
		end
		gs=pap[1].ss[issp].grad_mod[:rrvz]
		for i in eachindex(gmodrrvz)
			gmodrrvz[i]+=gs[i]
		end
	end
end

function update_gradient!(pac)
	nx, nz=pac.ic[:nx], pac.ic[:nz]
	nznxd = prod(length.(pac.model.mgrid))

	# combine rrvx and rrvz
	grad_modrr!(pac)

	gradient=pac.gradient
	# truncate
	gmodtt=view(pac.grad_mod[:tt],npml+1:nz-npml,npml+1:nx-npml)
	gmodrr=view(pac.grad_mod[:rr],npml+1:nz-npml,npml+1:nx-npml)
	for i in 1:nznxd
		# parameterization is  [:KI, :ρI, :null]
		gradient[i]=gmodtt[i]  # update gmodtt
		gradient[nznxd+i]=gmodrr[i]  # update gmodrr
	end
end

