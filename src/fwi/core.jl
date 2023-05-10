
"""
```julia
pa=SeisInvExpt(attrib_mod,attrib_inv, attrib; rfields, parameterization)
```
Predefined gallery of `SeisInvExpt`. Choose `attrib::Symbol`
* `=:pizza` an experiment fast enough to be run on a laptop
* `=:downhole` sources and receivers on a drill-string 
"""
function SeisInvExpt(paf::PFdtd, dobs::Records, migrid=paf.c.medium.mgrid; loss=L2DistLoss())

    # modeling mesh
    mmgrid = reverse(paf.c.exmedium.mgrid) # to x, y, z

    P = get_proj_matrix(migrid, mmgrid, use_gpu=_fd_use_gpu, weight_fn=(w, cc) -> Data.Number(w))

    # don't change the order of first two entries here
    # P will be repeated based on the number of medium parameters, e.g., KI and rho
    return (; paf, dobs, P, migrid, loss, mfull=get_modelvector(paf), gmfull=get_modelvector(paf))
end