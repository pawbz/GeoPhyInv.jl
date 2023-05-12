

function get_modelvector(pa::T) where {T<:NamedTuple{<:Any,<:Tuple{<:PFdtd,<:Records,Vararg}}}

    m = Data.Array(zeros(length(pa.mparams)*prod(length.(pa.migrid))))
    update!(pa.mfull, pa.paf, pa.mparams)

    mc = chunk(m, length(pa.mparams))
    mfullc = chunk(pa.mfull, length(pa.mparams))
    broadcast(mc, mfullc) do mc1, mfullc1
        mul!(mc1, adjoint(pa.P), mfullc1)
    end
    return m
end


function lossvalue(m, pa::T) where {T<:NamedTuple{<:Any,<:Tuple{<:PFdtd,<:Records,Vararg}}}
    mc = chunk(m, length(pa.mparams))
    mfullc = chunk(pa.mfull, length(pa.mparams))
    broadcast(mc, mfullc) do mc1, mfullc1
        mul!(mfullc1, pa.P, mc1)
    end

    return lossvalue(pa.mfull, pa.loss, pa.dobs, pa.paf, pa.mparams)
end

function gradient!(g, m, pa::T) where {T<:NamedTuple{<:Any,<:Tuple{<:PFdtd,<:Records,Vararg}}}

    mc = chunk(m, length(pa.mparams))
    mfullc = chunk(pa.mfull, length(pa.mparams))
    broadcast(mc, mfullc) do mc1, mfullc1
        mul!(mfullc1, pa.P, mc1)
    end

    gradient!(pa.gmfull, pa.mfull, pa.loss, pa.dobs, pa.paf, pa.mparams)

    gmc = chunk(g, length(pa.mparams))
    gmfullc = chunk(pa.gmfull, length(pa.mparams))
    broadcast(gmc, gmfullc) do gmc1, gmfullc1
        mul!(gmc1, adjoint(pa.P), gmfullc1)
    end
end

