

function get_modelvector(pa::T) where {T<:NamedTuple{<:Any,<:Tuple{<:PFdtd,<:Records,Vararg}}}

    mmgrid = pa.paf.c.medium.grid
    migrid = pa.migrid


    P = broadcast(migrid, mmgrid) do mi, mm
        Data.Array(Array(get_proj_matrix([mi], [mm], use_gpu=_fd_use_gpu, number=Data.Number)))
    end

    m = Data.Array(zeros(length(pa.mparams) * prod(length.(pa.migrid))))
    update!(pa.mfull, pa.paf, pa.mparams)

    mc = chunk(m, length(pa.mparams))
    mfullc = chunk(pa.mfull, length(pa.mparams))
    broadcast(mc, mfullc) do mc1, mfullc1
        mc11 = reshape(view(mc1, :), length.(migrid)...)
        mfullc11 = reshape(view(mfullc1, :), length.(mmgrid)...)
        apply_proj_matrix!(mc11, mfullc11, P...)
    end
    return m
end

function lossvalue(m, pa::T) where {T<:NamedTuple{<:Any,<:Tuple{<:PFdtd,<:Records,Vararg}}}

    mmgrid = pa.paf.c.medium.grid
    migrid = pa.migrid

    mc = chunk(m, length(pa.mparams))
    mfullc = chunk(pa.mfull, length(pa.mparams))

    # model parameters on modeling grid = P x model paramters on inversion grid
    broadcast(mc, mfullc) do mc1, mfullc1
        mc11 = reshape(view(mc1, :), length.(migrid)...)
        mfullc11 = reshape(view(mfullc1, :), length.(mmgrid)...)
        apply_proj_matrix!(mfullc11, mc11, pa.P...)
    end

    # simply do a forward solve using mfull, and return 
    return lossvalue(pa.mfull, pa.loss, pa.dobs, pa.paf, pa.mparams)
end

function gradient!(g, m, pa::T) where {T<:NamedTuple{<:Any,<:Tuple{<:PFdtd,<:Records,Vararg}}}

    mmgrid = pa.paf.c.medium.grid
    migrid = pa.migrid
    mc = chunk(m, length(pa.mparams))
    mfullc = chunk(pa.mfull, length(pa.mparams))
    # model parameters on modeling grid = P x model paramters on inversion grid
    broadcast(mc, mfullc) do mc1, mfullc1
        mc11 = reshape(view(mc1, :), length.(migrid)...)
        mfullc11 = reshape(view(mfullc1, :), length.(mmgrid)...)
        apply_proj_matrix!(mfullc11, mc11, pa.P...)
    end

    gradient!(pa.gmfull, pa.mfull, pa.loss, pa.dobs, pa.paf, pa.mparams)

    gmc = chunk(g, length(pa.mparams))
    gmfullc = chunk(pa.gmfull, length(pa.mparams))
    # gradient w.r.t. model parameters on inversion grid = P' x gradient w.r.t. model paramters on modeling grid
    broadcast(gmc, gmfullc) do gmc1, gmfullc1
        gmc11 = reshape(view(gmc1, :), length.(migrid)...)
        gmfullc11 = reshape(view(gmfullc1, :), length.(mmgrid)...)
        apply_proj_matrix!(gmc11, gmfullc11, map(transpose, pa.P)...)
        mul!(gmc1, pa.P, gmfullc1)
    end
end

