

function get_modelvector(pa::T) where {T<:NamedTuple{<:Any,<:Tuple{<:PFdtd,<:Records,Vararg}}}

    mmgrid = pa.paf.c.medium.grid

    P = get_proj_matrix(pa.migrid, mmgrid, use_gpu=_fd_use_gpu, number=Data.Number)

    m = Data.Array(zeros(length(pa.mparams) * prod(length.(pa.migrid))))
    update!(pa.mfull, pa.paf, pa.mparams)

    mc = chunk(m, length(pa.mparams))
    mfullc = chunk(pa.mfull, length(pa.mparams))
    broadcast(mc, mfullc) do mc1, mfullc1
        mul!(mc1, transpose(P), mfullc1)
    end
    return m
end

function lossvalue(m, pa::T) where {T<:NamedTuple{<:Any,<:Tuple{<:PFdtd,<:Records,Vararg}}}
    mc = chunk(m, length(pa.mparams))
    mfullc = chunk(pa.mfull, length(pa.mparams))

    # model parameters on modeling grid = P x model paramters on inversion grid
    broadcast(mc, mfullc) do mc1, mfullc1
        mul!(mfullc1, transpose(pa.P), mc1)
    end

    # simply do a forward solve using mfull, and return 
    return lossvalue(pa.mfull, pa.loss, pa.dobs, pa.paf, pa.mparams)
end

function gradient!(g, m, pa::T) where {T<:NamedTuple{<:Any,<:Tuple{<:PFdtd,<:Records,Vararg}}}

    mc = chunk(m, length(pa.mparams))
    mfullc = chunk(pa.mfull, length(pa.mparams))
    # model parameters on modeling grid = P x model paramters on inversion grid
    broadcast(mc, mfullc) do mc1, mfullc1
        mul!(mfullc1, transpose(pa.P), mc1)
    end

    gradient!(pa.gmfull, pa.mfull, pa.loss, pa.dobs, pa.paf, pa.mparams)

    gmc = chunk(g, length(pa.mparams))
    gmfullc = chunk(pa.gmfull, length(pa.mparams))
    # gradient w.r.t. model parameters on inversion grid = P' x gradient w.r.t. model paramters on modeling grid
    broadcast(gmc, gmfullc) do gmc1, gmfullc1
        mul!(gmc1, pa.P, gmfullc1)
    end
end

