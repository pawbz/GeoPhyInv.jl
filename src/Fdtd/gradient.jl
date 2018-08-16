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
	#gmodtt!(grad_modtt,p,pp,ppp,dpdx,dpdz,pac.nx,pac.nz)
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


#=
@inbounds @fastmath function gmodtt!(grad_modtt,p,pp,ppp,dpdx,dpdz,nx,nz)
	pw=p[1]
	pw2=p[2]
	ppw=pp[1]
	ppw2=pp[2]
	pppw2=ppp[2]
	pppw=ppp[1]
	dpdxw2=dpdx[2]
	dpdzw2=dpdz[2]
	for ix=1:nx
		@simd for iz=1:nz
			# p at [it], pp at [it-1]	# dpdx and dpdz at [it]	# gradients w.r.t inverse of modtt, i.e., 1/rho/c^2 
			@inbounds grad_modtt[iz,ix] += -1.0*(pw2[iz,ix,1]-ppw2[iz,ix,1])*pw[iz,ix,1]
		end
	end
end
=#
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

@inbounds @fastmath function scale_gradient!(issp::Int64,pass::Vector{Paramss},δ)
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
	gmodtt=pac.grad_modtt_stack
	gmodrrvx=pac.grad_modrrvx_stack
	gmodrrvz=pac.grad_modrrvz_stack
	pass=pap.ss
	for issp in 1:length(pass)
		gs=pass[issp].grad_modtt
		gss=view(gs,np+1:nz-np,np+1:nx-np)
		for i in 1:nznxd
			gmodtt[i]+=gss[i]  # only update gmodtt
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
end

function update_gradient!(pac::Paramc)
	nznxd = pac.model.mgrid.nz*pac.model.mgrid.nx

	# combine rrvx and rrvz
	grad_modrr!(pac)

	gradient=pac.gradient
	for i in 1:nznxd
		# parameterization is  [:KI, :ρI, :null]
		gradient[i]=pac.grad_modtt_stack[i]  # update gmodtt
		gradient[nznxd+i]=pac.grad_modrr_stack[i]  # update gmodrr
	end

end

