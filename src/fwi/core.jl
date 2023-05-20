
"""
```julia
pa=SeisInvExpt(attrib_mod,attrib_inv, attrib; rfields, parameterization)
```
Predefined gallery of `SeisInvExpt`. Choose `attrib::Symbol`
* `=:pizza` an experiment fast enough to be run on a laptop
* `=:downhole` sources and receivers on a drill-string 
"""
function SeisInvExpt(paf::PFdtd, dobs::Records, migrid::Vector{<:AbstractRange}=paf.c.medium.mgrid, mparams=Medium(paf.c.attrib); loss=L2DistLoss())
    @assert length(migrid) == length(paf.c.medium.mgrid)
    # modeling mesh
    mmgrid = reverse(paf.c.exmedium.mgrid) # to x, y, z

    P = get_proj_matrix(migrid, mmgrid, use_gpu=_fd_use_gpu, number=Data.Number)

    # don't change the order of first two entries here
    # P will be repeated based on the number of medium parameters, e.g., KI and rho
    mfull=get_modelvector(paf, mparams)
    return (; paf, dobs, P, migrid, mparams, loss, mfull, gmfull=similar(mfull))
end

# same as above, but make migrid using the number of cells Nigrid
function SeisInvExpt(paf::PFdtd, dobs::Records, Nigrid::Vector{Int}, mparams=Medium(paf.c.attrib); loss=L2DistLoss())
    @assert length(Nigrid) == length(paf.c.medium.mgrid)
    migrid = broadcast(paf.c.medium.mgrid, Nigrid) do m, N
        range(first(m)+4*step(m), stop=last(m)-4*step(m), length=N)
    end
    return SeisInvExpt(paf, dobs, migrid, mparams, loss=loss)
end