

function add_born_sources!(issp::Int64, pac, pap)

	δx24I=pac.fc[:δx24I]; δz24I=pac.fc[:δz24I]; 
	δxI=pac.fc[:δxI]; δzI=pac.fc[:δzI]; 
	δt=pac.fc[:δt]
	δtI=pac.fc[:δtI]
	δmodKI=pac.δmod[:KI]; modK=pac.mod[:K];
	δmodrhovxI=pac.δmod[:rhovxI]; δmodrhovzI=pac.δmod[:rhovzI]
	nx=pac.ic[:nx]; nz=pac.ic[:nz]

	born_svalue_stack=pap[1].born_svalue_stack

	p1=pap[1].w2[:t][:p]
	p2=pap[2].w2[:t][:p]
	p1p=pap[1].w2[:tp][:p]
	p1pp=pap[1].w2[:tpp][:p]
	dpdx1=pap[1].w2[:dx][:p]
	dpdx2=pap[2].w2[:dx][:p]
	dpdz1=pap[1].w2[:dz][:p]
	dpdz2=pap[2].w2[:dz][:p]

	# secondary sources for Born modeling
	# adding born sources from pressure(:,:,1) to pressure(:,:,2)
	# upto until [it-2]
	# lambdaI scatterrer source term at [it-1]
	# p is at [it], pp is at [it-1], ppp is at [it-2]
	# dpdx is at [it-1] and dpdz is at [it-1]
	# modrhovxI scatterrer source term at [it-1]
	# modrhovzI scatterrer source term at [it-1]

	# compute strength of secondary sources to δmod
	born_stackKI!(born_svalue_stack,p1,p1p,p1pp,δmodKI,nx,nz,δtI,δt)
	born_stackrhovxI!(born_svalue_stack,dpdx1,δmodrhovxI,nx,nz,δx24I,δt)
	born_stackrhovzI!(born_svalue_stack,dpdz1,δmodrhovzI,nx,nz,δz24I,δt)
	born_stack!(p2,born_svalue_stack,modK,nx,nz,δt)
end
@inbounds @fastmath function born_stackKI!(born_svalue_stack,p1,p1p,p1pp,δmodKI,nx,nz,δtI,δt)
	for ix=1:nx
		@simd for iz=1:nz
			born_svalue_stack[iz,ix] += δt * ((-1.0 * (p1pp[iz,ix] + p1[iz,ix] - 2.0 * p1p[iz,ix]) * δmodKI[iz,ix] * δtI * δtI)) 
		end
	end
end
@inbounds @fastmath function born_stackrhovxI!(born_svalue_stack,dpdx1,δmodrhovxI,nx,nz,δx24I,δt)
	for ix=3:nx-1
		@simd for iz=1:nz
			born_svalue_stack[iz,ix] += 
				δt * ((27.e0*dpdx1[iz,ix] * δmodrhovxI[iz,ix] -27.e0*dpdx1[iz,ix-1] * δmodrhovxI[iz,ix-1] -dpdx1[iz,ix+1] * δmodrhovxI[iz,ix+1] +dpdx1[iz,ix-2] * δmodrhovxI[iz,ix-2] ) * (δx24I)) 
		end
	end

end
@inbounds @fastmath function born_stackrhovzI!(born_svalue_stack,dpdz1,δmodrhovzI,nx,nz,δz24I,δt)
	for ix=1:nx
		@simd for iz=3:nz-1
			born_svalue_stack[iz,ix] += 
				δt * ((27.e0*dpdz1[iz,ix] * δmodrhovzI[iz,ix] -27.e0*dpdz1[iz-1,ix] * δmodrhovzI[iz-1,ix] -dpdz1[iz+1,ix] * δmodrhovzI[iz+1,ix] +dpdz1[iz-2,ix] * δmodrhovzI[iz-2,ix] ) * (δz24I))  
		end
	end
end
# add secondary sources to second pw
@inbounds @fastmath function born_stack!(p2,born_svalue_stack,modK,nx,nz,δt)
	for ix=1:nx
		@simd for iz=1:nz
			p2[iz,ix] += born_svalue_stack[iz,ix] * modK[iz,ix] * δt  #* δxI * δzI 
		end
	end
end

