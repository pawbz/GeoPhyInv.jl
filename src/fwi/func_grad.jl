

function get_modelvector(pa::T) where {T<:NamedTuple{<:Any,<:Tuple{<:PFdtd,<:Records,Vararg}}}

    m = Data.Array(zeros(length(pa.paf.c.mod)*prod(length.(pa.migrid))))
    update!(pa.mfull, pa.paf)

    mc = chunk(m, length(pa.paf.c.mod))
    mfullc = chunk(pa.mfull, length(pa.paf.c.mod))
    broadcast(mc, mfullc) do mc1, mfullc1
        mul!(mc1, adjoint(pa.P), mfullc1)
    end
    return m
end


function lossvalue(m, pa::T) where {T<:NamedTuple{<:Any,<:Tuple{<:PFdtd,<:Records,Vararg}}}
    mc = chunk(m, length(pa.paf.c.mod))
    mfullc = chunk(pa.mfull, length(pa.paf.c.mod))
    broadcast(mc, mfullc) do mc1, mfullc1
        mul!(mfullc1, pa.P, mc1)
    end

    return lossvalue(pa.mfull, pa.loss, pa.dobs, pa.paf)
end

function gradient!(g, m, pa::T) where {T<:NamedTuple{<:Any,<:Tuple{<:PFdtd,<:Records,Vararg}}}

    mc = chunk(m, length(pa.paf.c.mod))
    mfullc = chunk(pa.mfull, length(pa.paf.c.mod))
    broadcast(mc, mfullc) do mc1, mfullc1
        mul!(mfullc1, pa.P, mc1)
    end

    gradient!(pa.gmfull, pa.mfull, pa.loss, pa.dobs, pa.paf)

    gmc = chunk(g, length(pa.paf.c.mod))
    gmfullc = chunk(pa.gmfull, length(pa.paf.c.mod))
    broadcast(gmc, gmfullc) do gmc1, gmfullc1
        mul!(gmc1, adjoint(pa.P), gmfullc1)
    end
end

