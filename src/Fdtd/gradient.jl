


function grad_modrr!(pac::Paramc)
	@simd for i in eachindex(pac.grad_modrr_stack)
		@inbounds pac.grad_modrr_stack[i] = (pac.grad_modrrvx_stack[i] + pac.grad_modrrvz_stack[i]) * (0.5)
	end
	grad_modrr_sprayrrvx!(pac.grad_modrr_stack,pac.grad_modrrvx_stack)
	grad_modrr_sprayrrvz!(pac.grad_modrr_stack,pac.grad_modrrvz_stack)
end
function grad_modrr_sprayrrvx!(grad_modrr_stack,grad_modrrvx_stack)
	for ix=2:size(grad_modrr_stack,2)-1
		for iz=2:size(grad_modrr_stack,1)-1
			@inbounds grad_modrr_stack[iz,ix+1] +=  0.5e0 * grad_modrrvx_stack[iz,ix]
		end
	end
end
function grad_modrr_sprayrrvz!(grad_modrr_stack,grad_modrrvz_stack)
	for ix=2:size(grad_modrr_stack,2)-1
		for iz=2:size(grad_modrr_stack,1)-1
			@inbounds grad_modrr_stack[iz+1,ix] +=  0.5e0 * grad_modrrvz_stack[iz,ix]
		end
	end
end

function stack_grads!(pac::Paramc, pap::Paramp)
	np=pac.model.mgrid.npml # to truncate the gradients in PML region
	nx, nz=pac.nx, pac.nz
	nznxd = pac.model.mgrid.nz*pac.model.mgrid.nx

	# theses are SharedArrays
	grad_stack=pac.grad_stack
	gmodrrvx=pac.grad_modrrvx_stack
	gmodrrvz=pac.grad_modrrvz_stack
	pass=pap.ss
	for issp in 1:length(pass)
		gs=pass[issp].grad_modtt
		gss=view(gs,np+1:nz-np,np+1:nx-np)
		for i in 1:nznxd
			grad_stack[i]+=gss[i]  # only update gmodtt
		end
		gs=pass[issp].grad_modrrvx
		gss=view(gs,np+1:nz-np,np+1:nx-np)
		for i in 1:nznxd
			gmodrrvx[i] += gss[i]
		end
		gs=pass[issp].grad_modrrvz
		gss=view(gs,np+1:nz-np,np+1:nx-np)
		for i in 1:nznxd
			gmodrrvz[i] += gss[i]
		end
	end
	# combine rrvx and rrvz
	grad_modrr!(pac::Paramc)
	for i in 1:nznxd
		grad_stack[nznxd+i]+=pac.grad_modrr_stack[i]  # update gmodrr
	end

end

function update_gmodel!(pac::Paramc)
	nznxd = pac.model.mgrid.nz*pac.model.mgrid.nx
	
	for i in 1:nznxd
		pac.grad_stack[i]=Models.χg(pac.grad_stack[i],pac.model.ref.KI,1)
		pac.grad_stack[nznxd+i]=Models.χg(pac.grad_stack[nznxd+i],pac.model.ref.ρI,1)
	end

	Models.Seismic_chainrule!(pac.gmodel, pac.model, pac.grad_stack, [:χKI, :χρI, :null], 1)

	for i in 1:nznxd
		pac.grad_stack[i]=Models.χg(pac.grad_stack[i],pac.model.ref.KI,-1)
		pac.grad_stack[nznxd+i]=Models.χg(pac.grad_stack[nznxd+i],pac.model.ref.ρI,-1)
	end
end


# This routine ABSOLUTELY should not allocate any memory, called inside time loop.
@inbounds @fastmath function compute_gradient!(issp::Int64, pac::Paramc, pass::Vector{Paramss}, pap::Paramp)
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
	for ix=1:nx
		@simd for iz=1:nz
			# p at [it], pp at [it-1]	# dpdx and dpdz at [it]	# gradients w.r.t inverse of modtt, i.e., 1/rho/c^2 
			@inbounds grad_modtt[iz,ix] += ((-1.0 * (ppp[iz, ix, 1,1] + p[iz, ix,  1,1] - 2.0 * pp[iz, ix,  1,1]) * δtI * δtI) *  pp[iz,ix,1,2])
		end
	end
end
@inbounds @fastmath function gmodrrvx!(grad_modrrvx,dpdx,nx,nz)
	for ix=1:nx
		@simd for iz=1:nz
			# gradient w.r.t. inverse of rho on vx and vz grids
			@inbounds grad_modrrvx[iz,ix] += (- dpdx[iz,ix,1,2]*dpdx[iz,ix,1,1])
		end
	end

end
@inbounds @fastmath function gmodrrvz!(grad_modrrvz,dpdz,nx,nz)
	for ix=1:nx
		@simd for iz=1:nz
			# gradient w.r.t. inverse of rho on vx and vz grids
			@inbounds grad_modrrvz[iz,ix] += (- dpdz[iz,ix,1,2]*dpdz[iz,ix,1,1])
		end
	end
end

@inbounds @fastmath function scale_gradient!(issp::Int64,pass::Vector{Paramss},δ)
	grad_modtt=pass[issp].grad_modtt
	grad_modrrvx=pass[issp].grad_modrrvx
	grad_modrrvz=pass[issp].grad_modrrvz
	"gradient is formed by intergration over time, hence multiply with δt, but why not?"
	"I don't completely understand where the factors δx and δz are coming from..."
	"probably the source term should not be multiplied by δxI and δzI during adjoint propagation"
	scale!(grad_modtt,δ)
	scale!(grad_modrrvx,δ)
	scale!(grad_modrrvz,δ)
end


