
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

    paconv = [Conv.Pconv(Data.Number, dsize=size(d), ssize=(div(length(d.grid), 2), ), gsize=size(d)) for d in dobs]
    return (; paf, dobs, P, migrid, mparams, loss, mfull, gmfull=similar(mfull)), (; paconv, dobs, dcal=paf.c.data, dobs0=deepcopy(dobs))
end

# same as above, but make migrid using the number of cells Nigrid
function SeisInvExpt(paf::PFdtd, dobs::Records, Nigrid::Vector{Int}, mparams=Medium(paf.c.attrib); loss=L2DistLoss())
    @assert length(Nigrid) == length(paf.c.medium.grid)
    migrid = broadcast(paf.c.medium.grid, Nigrid) do m, N
        range(first(m)+4*step(m), stop=last(m)-4*step(m), length=N)
    end
    return SeisInvExpt(paf, dobs, migrid, mparams, loss=loss)
end


# core algorithm using IterativeSolvers.jl

function update!(pa::T) where {T<:NamedTuple{<:Any,<:Tuple{<:Vector{<:Conv.Pconv},<:Records,Vararg}}}
    # broadcast over supersources
    broadcast(pa.paconv, pa.dobs0, pa.dobs, pa.dcal) do paconv, dobs0, dobs, dcal
        copyto!(paconv.d, dcal)
        copyto!(paconv.g, dobs0)
   
        dvec = view(paconv.d, :)
        w = view(paconv.s, :)

        # create operator
        A=Conv.operator(paconv[1], Conv.G());

        IterativeSolvers.lsmr!(w, A, dvec)

        copyto!(paconv.s, w)
        Conv.mod!(paconv, Conv.D())
        copyto!(dobs, paconv.d)
    end
end