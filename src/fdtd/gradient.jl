
function sum_grads!(::Val{:adjoint}, pac, pap)
    g1 = pac.gradients
    for issp = 1:length(pap[1].ss)
        g = pap[1].ss[issp].gradients
        for name in names(g1)[1]
            for i in eachindex(g1[name])
                g1[name][i] = g1[name][i] + g[name][i]
            end
        end
    end
end
function sum_grads!(::Val{:forward}, pac, pap)
end


function gradlame!(issp, pap, pac::T) where {T<:P_common{<:FdtdAcoustic}}
    g = pap[1].ss[issp].gradients[:K]
    pforward = pap[1].w1[:t][:p]
    padjoint = pap[2].w1[:t][:p]
    padjoint_previous = pap[2].w1[:tp][:p]
	
    @parallel compute_gmodKI!(g, pforward, padjoint, padjoint_previous)
end


@parallel function compute_gmodKI!(g, pforward, padjoint, padjoint_previous)
    @all(g) = @all(g) + @all(pforward) * (@all(padjoint) - @all(padjoint_previous)) 
    return
end
function gradrho!(issp, pap, pac::T) where {T<:P_common{<:FdtdAcoustic,2}}
end



function compute_gradient!(::Val{:adjoint}, issp::Int64, pac, pap)
    gradlame!(issp, pap, pac)
    gradrho!(issp, pap, pac)
end

function compute_gradient!(::Val{:forward}, issp::Int64, pac, pap)
end
#=
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


=#