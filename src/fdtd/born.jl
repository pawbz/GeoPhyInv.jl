function add_born_sources_stress!(::Val{:forward}, pap, pac::T) where {T<:P_common{FdtdAcoustic{Born},2}}
    @parallel compute_p!(pap[2].w1[:t][:p], pap[1].w1[:t][:dvxdx], pap[1].w1[:t][:dvzdz], pac.δmod[:K], pac.fc[:dt])
end
function add_born_sources_velocity!(::Val{:forward}, pap, pac::T) where {T<:P_common{FdtdAcoustic{Born},2}}
    @parallel compute_v!(
        pap[2].w1[:t][:vx],
        pap[2].w1[:t][:vz],
        pac.δmod[:rhoI],
        pap[1].w1[:t][:dpdx],
        pap[1].w1[:t][:dpdz],
        pac.fc[:dt],
    )#
end
# 

# function gradlame!(issp, pap, pac::T) where {T<:P_common{<:FdtdAcoustic}}
function add_born_sources_stress!(::Val{:forward}, pap, pac::T) where {T<:P_common{FdtdElastic{FullWave}}}
end
function add_born_sources_velocity!(::Val{:forward}, pap, pac::T) where {T<:P_common{FdtdElastic{FullWave}}}
end
function add_born_sources_stress!(::Val{:forward}, pap, pac::T) where {T<:P_common{FdtdAcoustic{FullWave}}}
end
function add_born_sources_velocity!(::Val{:forward}, pap, pac::T) where {T<:P_common{FdtdAcoustic{FullWave}}}
end

function add_born_sources_stress!(::Val{:adjoint}, ::Any, ::Any) 
end
function add_born_sources_velocity!(::Val{:adjoint}, ::Any, ::Any) 
end



"""
```
F=LinearMap(pa)
```
If `pa` is an instance of `SeisInvExpt`, then 
return the linearized forward modeling operator `F`, such that
`F*x` can be computed without explicitly storing the operator matrix (see `LinearMaps.jl`).
The imaging/migration operator is given by `transpose(F)`. 
These operators are the building blocks of iterative optimization schemes.
"""
function LinearMaps.LinearMap(pa::T) where {T<:PFdtd{FdtdAcoustic{Born}}}
    fw = (y, x) -> forward_map!(y, x, pa)
    bk = (y, x) -> adjoint_map!(y, x, pa)

    # data (output) length
    nd = mapreduce(+, pa.c.data[1]) do d
        return sum(length.(d.d))
    end
    # length of medium (input)
    nm = sum(length.(pa.c.δmod))

    return LinearMap(fw, bk,
        nd, nm, 
        ismutating=true)
end

# m is nondimensionalized model vector
function forward_map!(d, m, pa::PFdtd)
    # copy input m to pac.δmod
    broadcast(pa.c.δmod, chunk(m, size=length(first(pa.c.δmod)))) do m1, m2
        CUDA.@allowscalar copyto!(m1, m2)
    end

    update!(pa, pa.c.srcwav, [1, 1], verbose=false)
    
    mode_save = pa.c.attrib_mod.mode
    pa.c.attrib_mod.mode = :forward
    update!(pa)
    # copy pac.data to d (only implemented for first supersource for now
    copyto!(d, Iterators.flatten(pa.c.data[2][1].d))
    pa.c.attrib_mod.mode = mode_save
end

function adjoint_map!(gm, d, pa::PFdtd)

    # copy input d to pac.srcwav (only implemented first source for now
    broadcast(pa.c.srcwav[2][1].d, chunk(d, size=length(first(pa.c.srcwav[2][1].d)))) do d1, d2
        CUDA.@allowscalar copyto!(d1, d2)
    end

    # time reversal
    foreach(pa.c.srcwav[2]) do S # each supersource
        foreach(S.d) do s # each field
            reverse!(s, dims=1)
        end
    end

    update!(pa, pa.c.srcwav, [-1, 1], verbose=false)

    mode_save = pa.c.attrib_mod.mode
    pa.c.attrib_mod.mode = :adjoint
    update!(pa)

    # copy pac.gradients to m
    copyto!(gm, Iterators.flatten(pa.c.gradients))
    pa.c.attrib_mod.mode = mode_save
end



# function add_born_sources!(issp::Int64, pac, pap)

# 	δx24I=pac.fc[:δx24I]; δz24I=pac.fc[:δz24I]; 
# 	δxI=pac.fc[:δxI]; δzI=pac.fc[:δzI]; 
# 	δt=pac.fc[:δt]
# 	δtI=pac.fc[:δtI]
# 	δmodKI=pac.δmod[:KI]; modK=pac.mod[:K];
# 	δmodrhovxI=pac.δmod[:rhovxI]; δmodrhovzI=pac.δmod[:rhovzI]
# 	nx=pac.ic[:nx]; nz=pac.ic[:nz]

# 	born_svalue_stack=pap[1].born_svalue_stack

# 	p1=pap[1].w1[:t][:p]
# 	p2=pap[2].w1[:t][:p]
# 	p1p=pap[1].w1[:tp][:p]
# 	p1pp=pap[1].w1[:tpp][:p]
# 	dpdx1=pap[1].w1[:dx][:p]
# 	dpdx2=pap[2].w1[:dx][:p]
# 	dpdz1=pap[1].w1[:dz][:p]
# 	dpdz2=pap[2].w1[:dz][:p]

# 	# secondary sources for Born modeling
# 	# adding born sources from pressure(:,:,1) to pressure(:,:,2)
# 	# upto until [it-2]
# 	# lambdaI scatterrer source term at [it-1]
# 	# p is at [it], pp is at [it-1], ppp is at [it-2]
# 	# dpdx is at [it-1] and dpdz is at [it-1]
# 	# modrhovxI scatterrer source term at [it-1]
# 	# modrhovzI scatterrer source term at [it-1]

# 	# compute strength of secondary sources to δmod
# 	born_stackKI!(born_svalue_stack,p1,p1p,p1pp,δmodKI,nx,nz,δtI,δt)
# 	born_stackrhovxI!(born_svalue_stack,dpdx1,δmodrhovxI,nx,nz,δx24I,δt)
# 	born_stackrhovzI!(born_svalue_stack,dpdz1,δmodrhovzI,nx,nz,δz24I,δt)
# 	born_stack!(p2,born_svalue_stack,modK,nx,nz,δt)
# end
# @inbounds @fastmath function born_stackKI!(born_svalue_stack,p1,p1p,p1pp,δmodKI,nx,nz,δtI,δt)
# 	for ix=1:nx
# 		@simd for iz=1:nz
# 			born_svalue_stack[iz,ix] += δt * ((-1.0 * (p1pp[iz,ix] + p1[iz,ix] - 2.0 * p1p[iz,ix]) * δmodKI[iz,ix] * δtI * δtI)) 
# 		end
# 	end
# end
# @inbounds @fastmath function born_stackrhovxI!(born_svalue_stack,dpdx1,δmodrhovxI,nx,nz,δx24I,δt)
# 	for ix=3:nx-1
# 		@simd for iz=1:nz
# 			born_svalue_stack[iz,ix] += 
# 				δt * ((27.e0*dpdx1[iz,ix] * δmodrhovxI[iz,ix] -27.e0*dpdx1[iz,ix-1] * δmodrhovxI[iz,ix-1] -dpdx1[iz,ix+1] * δmodrhovxI[iz,ix+1] +dpdx1[iz,ix-2] * δmodrhovxI[iz,ix-2] ) * (δx24I)) 
# 		end
# 	end

# end
# @inbounds @fastmath function born_stackrhovzI!(born_svalue_stack,dpdz1,δmodrhovzI,nx,nz,δz24I,δt)
# 	for ix=1:nx
# 		@simd for iz=3:nz-1
# 			born_svalue_stack[iz,ix] += 
# 				δt * ((27.e0*dpdz1[iz,ix] * δmodrhovzI[iz,ix] -27.e0*dpdz1[iz-1,ix] * δmodrhovzI[iz-1,ix] -dpdz1[iz+1,ix] * δmodrhovzI[iz+1,ix] +dpdz1[iz-2,ix] * δmodrhovzI[iz-2,ix] ) * (δz24I))  
# 		end
# 	end
# end
# # add secondary sources to second pw
# @inbounds @fastmath function born_stack!(p2,born_svalue_stack,modK,nx,nz,δt)
# 	for ix=1:nx
# 		@simd for iz=1:nz
# 			p2[iz,ix] += born_svalue_stack[iz,ix] * modK[iz,ix] * δt  #* δxI * δzI 
# 		end
# 	end
# end

