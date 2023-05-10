# returns non-dimensionalized model vector corresponding to pa.c.mod
function get_modelvector(pa::PFdtd)
    m = Data.Array(zeros(sum(length.(pa.c.mod))))
    update!(m, pa)
    return m
end

# updated pa.c.mod using a (non-dimensional) medium vector
# this is make pa.c.exmedium inconsistent with pa.c.mod, which is used in propagate functions
function update!(pa::PFdtd, m::AbstractVector)
    CUDA.allowscalar(true)
    mchunks = chunk(m, size=length(first(pa.c.mod)))
    # dimensionalize
    broadcast(pa.c.mod, mchunks) do mod, m
        copyto!(mod, m)
    end
    broadcast!(pa.c.mod, pa.c.mod, pa.c.ref_mod) do mod, r
        exp.(mod) .* r
    end
    CUDA.allowscalar(false)
end

# update m using pa.c.mod
function update!(m::AbstractVector, pa::PFdtd)
    CUDA.allowscalar(true)
    map!(pa.c.mod, pa.c.mod, pa.c.ref_mod) do mod, r
        log.(mod ./ r)
    end
    copyto!(m, Iterators.flatten(pa.c.mod))
    map!(pa.c.mod, pa.c.mod, pa.c.ref_mod) do mod, r
        exp.(mod) .* r
    end
    CUDA.allowscalar(false)
end

function Medium(::FdtdAcoustic)
    # the following medium parameters are stored on XPU
    return [:KI, :rho]
end
function Medium(::FdtdElastic)
    # the following medium parameters are stored on XPU
    return [:lambda, :M, :mu, :rho]
end


# function get_medium_grad_names(::FdtdElastic)
#     return [:KI, :rho]
# end

"""
This method should be executed only after the updating the main medium.
Update the `δmods` when a perturbed `medium_pert` is input.
The medium through which the waves are propagating 
is assumed to be the background medium.
"""
function update!(pa::PFdtd, medium::Medium, medium_pert::Medium)
    return update!(pa.c, medium, medium_pert)
end
function update!(pac::T, medium::Medium, medium_pert::Medium) where {T<:P_common{FdtdAcoustic{Born}}}
    update!(pac, medium)
    exmedium_pert = padarray(medium_pert, _fd_npextend, pac.pml_faces)
    broadcast(names(pac.δmod)[1]) do name
        copyto!(pac.δmod[name], exmedium_pert, [name])
    end
    # subtract reference values
    broadcast(pac.δmod, pac.mod) do δm, m
        @. δm = δm - m
    end
end

"""
```julia
update!(pa,medium_new)
```
Update `pa` with a new bundle of medium parameters `medium_new`, without additional memory allocation.
This routine is used during inversion, where medium parameters are iteratively updated. 
The ability to iteratively run the forward mediuming task (with no additional memory allocation) on  
various subsurface mediums is necessary while implementing inversion 
algorithms.
"""
function update!(pa::PFdtd, medium::Medium)
    return update!(pa.c, medium)
end
function update!(pac::T, medium::Medium) where {T<:P_common}
    copyto!(pac.medium, medium)
    padarray!(pac.exmedium, pac.medium, _fd_npextend, pac.pml_faces)
    broadcast(names(pac.mod)[1]) do name
        copyto!(pac.mod[name], pac.exmedium, [name])
    end
    return nothing
end
