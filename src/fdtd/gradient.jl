# This routine ABSOLUTELY should not allocate any memory, called inside time loop.
@inbounds @fastmath function compute_gradient!(issp::Int64, pac, pap)
	# aliases
	p=pap.p
	pp=pap.pp
	ppp=pap.ppp
	dpdx=pap.dpdx
	dpdz=pap.dpdz
	δtI=pac.δtI
	grad_modtt=pass[issp].grad_modtt
	grad_modrrvx=pass[issp].grad_modrrvx
	grad_modrrvz=pass[issp].grad_modrrvz

	gmodtt!(grad_modtt,p,pp,ppp,pac.nx,pac.nz,δtI)
	gmodrrvx!(grad_modrrvx,dpdx,pac.nx,pac.nz)
	gmodrrvz!(grad_modrrvz,dpdz,pac.nx,pac.nz)
end

@inbounds @fastmath function gmodtt!(grad_modtt,p,pp,ppp,nx,nz,δtI,)
	pppw=ppp[1]
	pw=p[1]
	ppw=pp[1]
	ppw2=pp[2]
	for ix=1:nx
		@simd for iz=1:nz
			# p at [it], pp at [it-1]	# dpdx and dpdz at [it]	# gradients w.r.t inverse of modtt, i.e., 1/rho/c^2 
			@inbounds grad_modtt[iz,ix] += ((-1.0 * (pppw[iz, ix, 1] + pw[iz, ix,  1] - 2.0 * ppw[iz, ix,  1]) * δtI * δtI) *  ppw2[iz,ix,1])
		end
	end
end


@inbounds @fastmath function gmodrrvx!(grad_modrrvx,dpdx,nx,nz)
	dpdxw1=dpdx[1]
	dpdxw2=dpdx[2]
	for ix=1:nx
		@simd for iz=1:nz
			# gradient w.r.t. inverse of rho on vx and vz grids
			@inbounds grad_modrrvx[iz,ix] += (- dpdxw2[iz,ix,1]*dpdxw1[iz,ix,1])
		end
	end

end
@inbounds @fastmath function gmodrrvz!(grad_modrrvz,dpdz,nx,nz)
	dpdzw1=dpdz[1]
	dpdzw2=dpdz[2]
	for ix=1:nx
		@simd for iz=1:nz
			# gradient w.r.t. inverse of rho on vx and vz grids
			@inbounds grad_modrrvz[iz,ix] += (- dpdzw2[iz,ix,1]*dpdzw1[iz,ix,1])
		end
	end
end

@inbounds @fastmath function scale_gradient!(issp::Int64,pass,δ)
	grad_modtt=pass[issp].grad_modtt
	grad_modrrvx=pass[issp].grad_modrrvx
	grad_modrrvz=pass[issp].grad_modrrvz
	"gradient is formed by intergration over time, hence multiply with δt, but why not?"
	"I don't completely understand where the factors δx and δz are coming from..."
	"probably the source term should not be multiplied by δxI and δzI during adjoint propagation"
	rmul!(grad_modtt,δ)
	rmul!(grad_modrrvx,δ)
	rmul!(grad_modrrvz,δ)
end




function grad_modrr!(pac)
	grad_modrr_sprayrr!(pac.grad_modrr_stack,pac.grad_modrrvx_stack,pac.grad_modrrvz_stack)
	grad_modrr_sprayrrvx!(pac.grad_modrr_stack,pac.grad_modrrvx_stack)
	grad_modrr_sprayrrvz!(pac.grad_modrr_stack,pac.grad_modrrvz_stack)
end

function stack_grads!(pac, pap)
	# theses are SharedArrays
	gmodtt=pac.grad_modtt_stack
	gmodrrvx=pac.grad_modrrvx_stack
	gmodrrvz=pac.grad_modrrvz_stack
	pass=pap.ss
	for issp in 1:length(pass)
		gs=pass[issp].grad_modtt
		for i in eachindex(gmodtt)
			gmodtt[i]+=gs[i]  # only update gmodtt
		end
		gs=pass[issp].grad_modrrvx
		for i in eachindex(gmodrrvx)
			gmodrrvx[i]+=gs[i]
		end
		gs=pass[issp].grad_modrrvz
		for i in eachindex(gmodrrvz)
			gmodrrvz[i]+=gs[i]
		end
	end
end

function update_gradient!(pac)
	nx, nz=pac.nx, pac.nz
	nznxd = prod(length.(pac.model.mgrid))

	# combine rrvx and rrvz
	grad_modrr!(pac)

	gradient=pac.gradient
	# truncate
	gmodtt=view(pac.grad_modtt_stack,npml+1:nz-npml,npml+1:nx-npml)
	gmodrr=view(pac.grad_modrr_stack,npml+1:nz-npml,npml+1:nx-npml)
	for i in 1:nznxd
		# parameterization is  [:KI, :ρI, :null]
		gradient[i]=gmodtt[i]  # update gmodtt
		gradient[nznxd+i]=gmodrr[i]  # update gmodrr
	end
end

