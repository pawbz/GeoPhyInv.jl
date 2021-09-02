

"""
Project density on to a staggerred grid using simple interpolation
"""
function get_rhovxI(rhoI::Array{Float64})
	rhovxI = zero(rhoI);
	get_rhovxI!(rhovxI, rhoI)
	return rhovxI
end
function get_rhovxI!(rhovxI, rhoI::Array{Float64})
	for ix = 1:size(rhoI, 2)-1
		for iz = 1:size(rhoI,1)
			rhovxI[iz, ix] = 0.5e0 *(rhoI[iz,ix+1] + rhoI[iz,ix])
		end
	end
	# last point, circular boundary conditions
	nx=size(rhoI, 2)
	for iz = 1:size(rhoI,1)
		rhovxI[iz, nx] = 0.5e0 *(rhoI[iz,1] + rhoI[iz,nx])
	end
	return nothing
end # get_rhovxI

"""
Project density on to a staggerred grid using simple interpolation
"""
function get_rhovzI(rhoI::Array{Float64})
	rhovzI = zero(rhoI);
	get_rhovzI!(rhovzI, rhoI)
	return rhovzI
end
function get_rhovzI!(rhovzI, rhoI::Array{Float64})
	nz=size(rhoI, 1)
	for ix = 1: size(rhoI, 2)
		for iz = 1: size(rhoI, 1)-1
			rhovzI[iz, ix] =  0.5e0 * (rhoI[iz+1,ix] + rhoI[iz,ix])
		end
		# last point, circular boundary conditions
		rhovzI[nz, ix] =0.5e0 * (rhoI[1,ix] + rhoI[nz,ix])
	end
	return nothing
end # get_rhovzI

function grad_modrr_sprayrr!(grad_modrr_stack,grad_modrhovxI_stack,grad_modrhovzI_stack)
	@simd for i in eachindex(grad_modrr_stack)
		@inbounds grad_modrr_stack[i] = (grad_modrhovxI_stack[i] + grad_modrhovzI_stack[i]) * (0.5)
	end
end

function grad_modrr_sprayrhovxI!(grad_modrr_stack,grad_modrhovxI_stack)
	for ix=1:size(grad_modrr_stack,2)-1
		for iz=1:size(grad_modrr_stack,1)
			@inbounds grad_modrr_stack[iz,ix+1] +=  0.5e0 * grad_modrhovxI_stack[iz,ix]
		end
	end
	nx=size(grad_modrr_stack,2)
	# circular boundary 
	for iz=1:size(grad_modrr_stack,1)
		@inbounds grad_modrr_stack[iz,1] +=  0.5e0 * grad_modrhovxI_stack[iz,nx]
	end
end
function grad_modrr_sprayrhovzI!(grad_modrr_stack,grad_modrhovzI_stack)
	nz=size(grad_modrr_stack,1)
	for ix=1:size(grad_modrr_stack,2)
		for iz=1:size(grad_modrr_stack,1)-1
			@inbounds grad_modrr_stack[iz+1,ix] +=  0.5e0 * grad_modrhovzI_stack[iz,ix]
		end
		# circular boundary
		grad_modrr_stack[1,ix] +=  0.5e0 * grad_modrhovzI_stack[nz,ix]
	end
end

