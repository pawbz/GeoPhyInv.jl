# returns non-dimensionalized model vector corresponding to pa.c.mod
function get_modelvector(pa::PFdtd, mparams)
    m = Data.Array(zeros(length(first(pa.c.mod)) * length(mparams)))
    update!(m, pa, mparams)
    return m
end

# updated pa.c.mod using a (non-dimensional) medium vector
# NOTE: this makes pa.c.exmedium inconsistent with pa.c.mod, which is used in propagate functions
function update!(pa::PFdtd, m::AbstractVector, mparams)
    CUDA.allowscalar(true)
    mchunks = chunk(m, length(mparams))
    # dimensionalize
    broadcast(enumerate(mparams)) do (i, mname)
        mod = pa.c.mod[mname]
        x = mchunks[i]
        copyto!(mod, x)
    end
    broadcast(mparams) do mname
        mod = pa.c.mod[mname]
        r = pa.c.ref_mod[mname]
        @. mod = exp(mod) * r
        # @. mod = mod * r + r
    end
    CUDA.allowscalar(false)

    # update dependencies
    update_dmod!(pa.c)

end

# update m using pa.c.mod
function update!(m::AbstractVector, pa::PFdtd, mparams)
    CUDA.allowscalar(true)
    broadcast(mparams) do mname
        mod = pa.c.mod[mname]
        r = pa.c.ref_mod[mname]
        @. mod = log(mod * inv(r))
        # @. mod = (mod - r) * inv(r)
    end
    copyto!(m, Iterators.flatten(pa.c.mod[mparams]))
    broadcast(mparams) do mname
        mod = pa.c.mod[mname]
        r = pa.c.ref_mod[mname]
        @. mod = exp(mod) * r
        # @. mod = mod * r + r
    end
    CUDA.allowscalar(false)
end

# the following medium parameters are stored on XPU
# averaged values of rho on velocity grid
# these methods output the fields on which the medium parameters lie
function MediumParameters(::FdtdAcoustic)
    @static if _fd_ndims == 2
        return [:invK, :rho], [:dtinvavxirho, :dtinvavzirho, :dtK], [dpdx(), dpdz(), p()] # independent and dependent
    else
        return [:invK, :rho], [:dtinvavxirho, :dtinvavyirho, :dtinvavzirho, :dtK], [dpdx(), dpdy(), dpdz(), p()]
    end

end
function MediumParameters(::FdtdElastic)
    @static if _fd_ndims == 2
        return [:invlambda, :invmu, :rho], [:dtinvavxirho, :dtinvavzirho, :dtavmu, :dtM, :dtlambda], [dtauxxdx(), dtauzzdz(), dvxdz(), tauxx(), tauxx()]
    else
        return [:invlambda, :invmu, :rho], [:dtinvavxirho, :dtinvavyirho, :dtinvavzirho, :dtavxzimu, :dtavxyimu, :dtavyzimu, :dtM, :dtlambda], [dtauxxdx(), dtauyydy(), dtauzzdz(), dvxdz(), dvxdy(), dvydz(), tauxx(), tauxx()]
    end
end

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
    broadcast(AxisArrays.names(pac.mod)[1]) do name
        copyto!(pac.mod[name], pac.exmedium, name)
    end
    # update dependencies
    update_dmod!(pac)
    return nothing
end


function update_dmod!(pac::T) where {T<:P_common{<:FdtdAcoustic,2}}
    mod=pac.mod; dmod=pac.dmod
    @parallel store_invavxi!(dmod[:dtinvavxirho], mod[:rho], pac.fc[:dt])
    @parallel store_invavzi!(dmod[:dtinvavzirho], mod[:rho], pac.fc[:dt])
    broadcast!(inv, dmod[:dtK], mod[:invK])
    rmul!(dmod[:dtK], pac.fc[:dt])
end
function update_dmod!(pac::T) where {T<:P_common{<:FdtdAcoustic,3}}
    mod=pac.mod; dmod=pac.dmod
    @parallel store_invavxi!(dmod[:dtinvavxirho], mod[:rho], pac.fc[:dt])
    @parallel store_invavyi!(dmod[:dtinvavyirho], mod[:rho], pac.fc[:dt])
    @parallel store_invavzi!(dmod[:dtinvavzirho], mod[:rho], pac.fc[:dt])
    broadcast!(inv, dmod[:dtK], mod[:invK])
    rmul!(dmod[:dtK], pac.fc[:dt])
end

function update_dmod!(pac::T) where {T<:P_common{<:FdtdElastic,2}}
    mod=pac.mod; dmod=pac.dmod
    @parallel store_invavxi!(dmod[:dtinvavxirho], mod[:rho], pac.fc[:dt])
    @parallel store_invavzi!(dmod[:dtinvavzirho], mod[:rho], pac.fc[:dt])
    broadcast!(inv, dmod[:dtlambda], mod[:invlambda])
    rmul!(dmod[:dtlambda], pac.fc[:dt])
    broadcast!(dmod[:dtM], mod[:invlambda], mod[:invmu]) do invλ, invμ
        inv(invλ) + 2.0 * inv(invμ)
    end
    rmul!(dmod[:dtM], pac.fc[:dt])

    @parallel store_invav!(dmod[:dtavmu], mod[:invmu], pac.fc[:dt])
end
function update_dmod!(pac::T) where {T<:P_common{<:FdtdElastic,3}}
    mod=pac.mod; dmod=pac.dmod
    @parallel store_invavxi!(dmod[:dtinvavxirho], mod[:rho], pac.fc[:dt])
    @parallel store_invavyi!(dmod[:dtinvavyirho], mod[:rho], pac.fc[:dt])
    @parallel store_invavzi!(dmod[:dtinvavzirho], mod[:rho], pac.fc[:dt])
    broadcast!(inv, dmod[:dtlambda], mod[:invlambda])
    rmul!(dmod[:dtlambda], pac.fc[:dt])
    broadcast!(dmod[:dtM], mod[:invlambda], mod[:invmu]) do invλ, invμ
        inv(invλ) + 2.0 * inv(invμ)
    end
    rmul!(dmod[:dtM], pac.fc[:dt])
    @parallel store_invavxzi!(dmod[:dtavxzimu], mod[:invmu], pac.fc[:dt])
    @parallel store_invavxyi!(dmod[:dtavxyimu], mod[:invmu], pac.fc[:dt])
    @parallel store_invavyzi!(dmod[:dtavyzimu], mod[:invmu], pac.fc[:dt])
end




@parallel function store_invavxi!(a, b, dt)
    @all(a) = dt / @av_xi(b)
    return
end

@parallel function store_invavyi!(a, b, dt)
    @all(a) = dt / @av_yi(b) 
    return
end

@parallel function store_invavzi!(a, b, dt)
    @all(a) = dt / @av_zi(b)
    return
end

@parallel function store_invav!(a, b, dt)
    @all(a) = dt / @av(b)
    return
end
@parallel function store_invavxzi!(a, b, dt)
    @all(a) = dt / @av_xzi(b)
    return
end
@parallel function store_invavxyi!(a, b, dt)
    @all(a) = dt / @av_xyi(b)
    return
end
@parallel function store_invavyzi!(a, b, dt)
    @all(a) = dt / @av_yzi(b)
    return
end