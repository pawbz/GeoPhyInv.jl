"""
```julia
pa=SeisInvExpt(paf, dobs, migrid, mparams; loss)
```
* `paf` struct for finite-difference modeling
* `dobs` observed data which is being compared to the modeled data (dobs = conv(w, dobs0))
* `migrid` inversion grid for medium parameters
* `mparams` medium parameterization
"""
function SeisInvExpt(paf::PFdtd, dobs::Records, migrid::Vector{<:AbstractRange}=paf.c.medium.grid, mparams=Medium(paf.c.attrib); loss=L2DistLoss())
    @assert length(migrid) == length(paf.c.medium.grid)
    # modeling mesh
    mmgrid = reverse(paf.c.exmedium.grid) # to x, y, z

    P = get_proj_matrix(migrid, mmgrid, use_gpu=_fd_use_gpu, number=Data.Number)

    # don't change the order of first two entries here
    # P will be repeated based on the number of medium parameters, e.g., KI and rho
    mfull=get_modelvector(paf, mparams)

    dobs0 = deepcopy(dobs)
    paconv = [Conv.Pconv(Data.Number, dsize=size(d), ssize=(length(d.grid), ), gsize=size(d)) for d in dobs0]
    return (; paf, dobs, dobs0, P, migrid, mparams, loss, mfull, gmfull=similar(mfull)), paconv
end

# same as above, but make migrid using the number of cells Nigrid
function SeisInvExpt(paf::PFdtd, dobs::Records, Nigrid::Vector{Int}, mparams=Medium(paf.c.attrib); loss=L2DistLoss())
    @assert length(Nigrid) == length(paf.c.medium.grid)
    migrid = broadcast(paf.c.medium.grid, Nigrid) do m, N
        range(first(m)+4*step(m), stop=last(m)-4*step(m), length=N)
    end
    return SeisInvExpt(paf, dobs, migrid, mparams, loss=loss)
end


# core algorithm using IterativeSolvers.jl to update the source wavelets so that dobs is modified
function update!(paw::T1, pac::T2) where {T1<:NamedTuple{<:Any,<:Tuple{<:PFdtd,<:Records,Vararg}}, T2<:Vector{<:Conv.Pconv}}

    # need to generate data first, paw.paf.c.data, will be used later
    m = get_modelvector(paw)
    loss_initial = lossvalue(m, paw)

    # broadcast over supersources
    broadcast(pac, paw.dobs0, paw.dobs, paw.paf.c.data[1]) do paconv, dobs0, dobs, dcal
        # the G is same as observed data (g=dobs0)
        copyto!(paconv.g, dobs0)
        # can we avoid these allocations?
        dvec = vec(dcal)
        w = zero(paconv.s)
        # matrix-free operator
        A=Conv.operator(paconv, Conv.G());
        IterativeSolvers.lsmr!(w, A, dvec)
        copyto!(paconv.s, w)
        # do conv with estimated w
        Conv.mod!(paconv, Conv.D())
        # update dobs in paw
        copyto!(dobs, paconv.d)
    end

    # estimate loss again with new dobs
    loss_final = lossvalue(m, paw)

    return (; loss_initial, loss_final)
end