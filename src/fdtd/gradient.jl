
function sum_grads!(::Val{:adjoint}, ::Val{2}, pac, pap)
    g1 = pac.gradients
    for issp = 1:length(pap[1].ss)
        g = pap[1].ss[issp].gradients
        for name in names(g1)[1]
            for i in eachindex(g1[name])
                CUDA.@allowscalar g1[name][i] = g1[name][i] + g[name][i]
            end
        end
    end
end
function sum_grads!(::Any, ::Any, pac, pap)
end


function gradlame!(issp, pap, pac::T) where {T<:P_common{<:FdtdAcoustic}}
    g = pap[1].ss[issp].gradients[:invK]
    pforward_previous = pap[1].w1[:tp][:p]
    pforward = pap[1].w1[:t][:p]
    padjoint_previous = pap[2].w1[:tp][:p]

    @parallel compute_gmodKI!(g, pforward, pforward_previous, padjoint_previous, pac.fc[:dtI])
end


@parallel function compute_gmodKI!(g, pforward, pforward_previous, padjoint_previous, dtI)
    @all(g) = @all(g) + @all(padjoint_previous) * (@all(pforward_previous) - @all(pforward)) * dtI
    return
end
function gradrho!(issp, pap, pac::T) where {T<:P_common{<:FdtdAcoustic,2}}
    vxforward_previous = pap[1].w1[:tp][:vx]
    vxforward = pap[1].w1[:t][:vx]
    vxadjoint_previous = pap[2].w1[:tp][:vx]
    vxbuffer = pap[1].velocity_buffer[:vx]
    @parallel compute_gmodrho!(vxbuffer, vxforward, vxforward_previous, vxadjoint_previous, pac.fc[:dtI])

    vzforward_previous = pap[1].w1[:tp][:vz]
    vzforward = pap[1].w1[:t][:vz]
    vzadjoint_previous = pap[2].w1[:tp][:vz]
    vzbuffer = pap[1].velocity_buffer[:vz]
    @parallel compute_gmodrho!(vzbuffer, vzforward, vzforward_previous, vzadjoint_previous, pac.fc[:dtI])

    g = pap[1].ss[issp].gradients[:rho]
	@parallel combine_gmodrho!(g, vxbuffer, vzbuffer)
end

@parallel function compute_gmodrho!(buffer, vforward, vforward_previous, vadjoint_previous, dtI)
    @all(buffer) = @all(vadjoint_previous) * (@all(vforward) - @all(vforward_previous)) * dtI
    return
end

@parallel function combine_gmodrho!(g, vxbuffer, vzbuffer)
    @inn(g) = @inn(g) - @av_xi(vxbuffer) - @av_zi(vzbuffer) 
    return
end

function compute_gradient!(::Val{:adjoint}, ::Val{2}, issp::Int64, pac, pap)
    gradlame!(issp, pap, pac)
    gradrho!(issp, pap, pac)
end

function compute_gradient!(::Any, ::Any, issp::Int64, pac, pap)
end