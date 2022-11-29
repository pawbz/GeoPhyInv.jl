
"""
Initialize an model with zeros, given mgrid and names.
"""
function Base.zero(::Type{Medium}, mgrid, names = [:vp, :rho])
    nvec = length.(mgrid)
    return Medium(
        map(StepRangeLen{Float64}, mgrid),
        NamedArray([zeros(Float64, nvec...) for i in names], (names,)),
        NamedArray([fill(0.0, 2) for i in names], (names,)),
        NamedArray([0.0 for i in names], (names,)),
    )
end
Medium(mgrid, names = [:vp, :rho]) = zero(Medium, mgrid, names)
Medium(mgrid, m::Medium) = zero(Medium, mgrid, names(m.m)[1])



"""
Return true if a `Medium` model is just allocated with zeros.
"""
function Base.iszero(mod::Medium)
    if (any(broadcast(iszero, mod.bounds)) && any(broadcast(iszero, mod.ref)))
        return true
    else
        # vp and rho cannot be zero
        for field in [:vp, :rho]
            if (field ∈ names(mod.m)[1])
                if (any(iszero.(mod[field])))
                    return true
                end
            end
        end
    end
    return false
end


"""
Return true if a `Medium` model is just allocated with zeros.
"""
function Base.isassigned(mod::Medium, fields = [:vp, :rho, :vs])
    for field in fields
        @assert field ∈ [:vp, :rho, :vs]
    end

    return any(broadcast(iszero, mod.bounds)) && any(broadcast(iszero, mod.ref))
end

"""
Return a similar model (size, bounds and reference values) to the input model, used for allocation.
"""
function Base.similar(mod::Medium)
    mod_new = deepcopy(mod)
    fill!(mod_new)
    return mod_new
end


"Compare if two `Medium` models are equal"
function Base.isequal(mod1::T, mod2::T) where {T<:Union{Medium}}
    fnames = fieldnames(typeof(mod1))
    vec = [(isequal(getfield(mod1, name), getfield(mod2, name))) for name in fnames]
    return all(vec)
end

"""
Return if two `Medium` models have same dimensions and bounds.
"""
function Base.size(mod::Medium)
    return length.(mod.mgrid)
end

"""
Return the number of dimensions of `Medium` 
"""
function Base.ndims(mod::Medium)
    return length(mod.mgrid)
end

"""
Copy for `Medium`. The media should have same bounds and sizes.
Doesn't allocate any memory. Note that only [:vp, :rho, :vs] fields are stored in structs.
"""
function Base.copyto!(modo::Medium, mod::Medium)
    for name in [:vp, :vs, :rho]
        ((name ∈ names(modo.m)[1]) & (name ∈ names(mod.m)[1])) &&
            (copyto!(modo.m[name], mod.m[name]))
    end
    return modo
end

"""
Fill medium with reference medium parameters to make it homogeneous.
```julia
fill!(medium)
```
"""
function Base.fill!(mod::Medium)
    for s in names(mod.m)[1]
        fill!(mod.m[s], mod.ref[s])
    end
    return mod
end


"""
Print information about `Medium`
"""
function Base.show(io::Base.IO, mod::Medium)
    println(
        io,
        "> ",
        ndims(mod),
        "-D Medium with ",
        dim_names(ndims(mod)),
        " in ",
        [[m[1], m[end]] for m in mod.mgrid],
    )
    println(io, "> number of samples:\t", length.(mod.mgrid))
    println(io, "> sampling intervals:\t", step.(mod.mgrid))
    try
        println(io, "> vp:\t", "min\t", minimum(mod[:vp]), "\tmax\t", maximum(mod[:vp]))
    catch
    end
    try
        println(io, "> vs:\t", "min\t", minimum(mod[:vs]), "\tmax\t", maximum(mod[:vs]))
    catch
    end
    try
        println(io, "> rho:\t", "min\t", minimum(mod[:rho]), "\tmax\t", maximum(mod[:rho]))
    catch
    end
    println(io, "> Bounds:")
    println(io, mod.bounds)
end


"""
Re-parameterization routine 
that modifies the fields 
`χvp` and `χrho` of an input seismic model
using two input vectors.

# Arguments

* `mod::Medium` : to be updated
* `x1::Array{Float64,2}` : contrast of inverse bulk modulus
* `x2::Array{Float64,2}` : contrast of inverse density
* `fields:::Vector{Symbol}` : [:χKI, :χrhoI]
"""
function Base.copyto!(
    mod::Medium,
    x::AbstractVector,
    fields::Vector{Symbol} = [:χKI, :χrhoI, :null],
)
    iszero(mod) ? error("mod cannot be zero") : nothing
    nznx = prod(length.(mod.mgrid))
    (length(x) ≠ (count(fields .≠ :null) * nznx)) && error("size x")
    if (fields == [:χKI, :χrhoI, :null])
        @inbounds for i = 1:nznx
            K = inv(χ(x[i], mod.ref[:KI], -1))
            rho = inv(χ(x[nznx+i], mod.ref[:rhoI], -1))
            mod.m[:vp][i] = sqrt(K * inv(rho))
            mod.m[:rho][i] = rho
        end
    elseif (fields == [:χKI, :null, :null])
        @inbounds for i = 1:nznx
            K = inv(χ(x[i], mod.ref[:KI], -1))
            rho = mod.m[:rho][i]
            mod.m[:vp][i] = sqrt(K * inv(rho))
            mod.m[:rho][i] = rho
        end
    elseif (fields == [:χvp, :χrho, :null])
        @inbounds for i = 1:nznx
            mod.m[:vp][i] = χ(x[i], mod.ref[:vp], -1)
        end
        @inbounds for i = 1:nznx
            mod.m[:rho][i] = χ(x[nznx+i], mod.ref[:rho], -1)
        end
    elseif (fields == [:null, :χrho, :null])
        @inbounds for i = 1:nznx
            mod.m[:rho][i] = χ(x[i], mod.ref[:rho], -1)
        end
    elseif (fields == [:χvp, :null, :null])
        @inbounds for i = 1:nznx
            mod.m[:vp][i] = χ(x[i], mod.ref[:vp], -1)
        end
    else
        error("invalid fields")
    end
    return mod
end




"""
Input perturbations corresponding to `parameterization`.
Output perturbations of KI and rhoI, i.e., without the contrast.
`mod` is just used for reference values.
"""
function Base.copyto!(
    δxout::AbstractVector,
    δx::AbstractVector,
    mod::Medium,
    parameterization,
)
    nznx = prod(length.(mod.mgrid))
    fill!(δxout, 0.0)
    if (parameterization == [:χKI, :χrhoI, :null])
        @inbounds for i = 1:nznx
            δxout[i] = δx[i] * mod.ref[:KI]
            δxout[nznx+i] = δx[nznx+i] * mod.ref[:rhoI]
        end
    elseif (parameterization == [:χKI, :null, :null])
        @inbounds for i = 1:nznx
            δxout[i] = δx[i] * mod.ref[:KI]
        end
    elseif (parameterization == [:χvp, :χrho, :null])
        @inbounds for i = 1:nznx
            δxout[i] = (-1.0 * δx[nznx+i] - 2.0 * δx[i]) * mod.ref[:KI]
            δxout[nznx+i] = (-1.0 * δx[nznx+i]) * mod.ref[:rhoI]
        end
    elseif (parameterization == [:χvp, :null, :null])
        @inbounds for i = 1:nznx
            δxout[i] = -2.0 * δx[i] * mod.ref[:KI]
        end
    elseif (parameterization == [:null, :χrho, :null])
        @inbounds for i = 1:nznx
            δxout[i] = -1.0 * δx[i] * mod.ref[:rhoI]
        end
    end
end


#pert_of_KI(dvp, drho, vp0, rho0) = - 2. / (rho0 * vp0^3) * dvp - 1 / (rho0^2 * vp0^2) * drho
#pert_of_rhoI(drho, vp0, rho0) = -1. / (rho0^2) * drho



"""
Get other dependent model parameters of a seismic model
that are not present in `Medium`.

* `:rhoI` : inverse of density
* `:Zp` : P-wave impedance
* `:lambda` : first Lame parameter
* `:mu` : shear modulus
* `:K` : bulk modulus
* `:M` : P-wave modulus
"""
function Base.getindex(mod::Medium, s::Symbol)
    nznx = length(mod.m[:vp])
    # check if vp, vs, and rho are present 
    @assert :vp ∈ names(mod.m)[1]
    @assert :rho ∈ names(mod.m)[1]
    vp = mod.m[:vp]
    rho = mod.m[:rho]
    if (:vs ∈ names(mod.m)[1])
        vs = mod.m[:vs]
    end
    if (s == :vp)
        return vp
    elseif (s == :rho)
        return rho
    elseif (s in [:vs, :χmu])
        @assert :vs ∈ names(mod.m)[1]
        if (s == :vs)
            return mod.m[:vs]
        end
    elseif (s == :Q)
        @assert :Q ∈ names(mod.m)[1]
        return mod.m[:Q]
    end

    # these are 3D arrays
    if (s in [:tau_epsilon, :tau_sigma])
        @assert :Q ∈ names(mod.m)[1] # need quality factor
        @assert :freqmin ∈ names(mod.fc)[1] # need minimum freq that is being modelled
        @assert :freqmax ∈ names(mod.fc)[1] # need max freq
        @assert :nsls ∈ names(mod.ic)[1] # need number of standard linear solids

        nz, nx = size(mod.m[:vp])
        # add derived fields, if not present
        (s ∉ names(mod.m3)[1]) &&
            (mod.m3 = vcat(mod.m3, NamedArray([zeros(mod.ic[:nsls], nz, nx)], [s])))
        x = mod.m3[s] # going to return this 
    else
        # add derived fields, if not present
        (s ∉ names(mod.m)[1]) && (mod.m = vcat(mod.m, NamedArray([zero(vp)], [s])))
        x = mod.m[s] # going to return this
    end




    # output derived fields
    if (s == :rhoI)
        @inbounds for i = 1:nznx
            x[i] = inv(rho[i])
        end
    elseif (s == :lambda)
        if (:vs ∈ names(mod.m)[1])
            @inbounds for i = 1:nznx
                x[i] = (vp[i] * vp[i] - 2 * vs[i] * vs[i]) * rho[i]
            end
        else
            @inbounds for i = 1:nznx
                x[i] = vp[i] * vp[i] * rho[i]
            end
        end
    elseif (s == :M)
        @inbounds for i = 1:nznx
            x[i] = vp[i] * vp[i] * rho[i]
        end
    elseif (s == :χrhoI)
        @inbounds for i = 1:nznx
            x[i] = χ(inv(rho[i]), mod.ref[:rhoI], 1)
        end
    elseif (s == :χrho)
        @inbounds for i = 1:nznx
            x[i] = χ(rho[i], mod.ref[:rho], 1)
        end
    elseif (s == :χvp)
        @inbounds for i = 1:nznx
            x[i] = χ(vp[i], mod.ref[:vp], 1)
        end
    elseif (s == :χK)
        if (:vs ∈ names(mod.m)[1])
            @inbounds for i = 1:nznx
                x[i] =
                    χ((vp[i] * vp[i] - 4.0 / 3.0 * vs[i] * vs[i]) * rho[i], mod.ref[:K], 1)
            end
        else
            @inbounds for i = 1:nznx
                x[i] = χ(vp[i] * vp[i] * rho[i], mod.ref[:K], 1)
            end
        end
    elseif (s == :χM)
        @inbounds for i = 1:nznx
            x[i] = χ(vp[i] * vp[i] * rho[i], mod.ref[:K], 1)
        end
    elseif (s == :mu)
        @inbounds for i = 1:nznx
            x[i] = mod.m[:vs][i] * mod.m[:vs][i] * rho[i]
        end
    elseif (s == :χmu)
        @inbounds for i = 1:nznx
            x[i] = χ(mod.m[:vs][i] * mod.m[:vs][i] * rho[i], mod.ref[:mu], 1)
        end
    elseif (s == :χKI)
        if (:vs ∈ names(mod.m)[1])
            @inbounds for i = 1:nznx
                x[i] = χ(
                    inv((vp[i] * vp[i] - 4.0 / 3.0 * vs[i] * vs[i]) * rho[i]),
                    mod.ref[:KI],
                    1,
                )
            end
        else
            @inbounds for i = 1:nznx
                x[i] = χ(inv(vp[i] * vp[i] * rho[i]), mod.ref[:KI], 1)
            end
        end
    elseif (s == :χMI)
        @inbounds for i = 1:nznx
            x[i] = χ(inv(vp[i] * vp[i] * rho[i]), mod.ref[:KI], 1)
        end
    elseif (s == :KI)
        if (:vs ∈ names(mod.m)[1])
            @inbounds for i = 1:nznx
                x[i] = inv((vp[i] * vp[i] - 4.0 / 3.0 * vs[i] * vs[i]) * rho[i])
            end
        else
            @inbounds for i = 1:nznx
                x[i] = inv(vp[i] * vp[i] * rho[i])
            end
        end
    elseif (s == :K)
        if (:vs ∈ names(mod.m)[1])
            @inbounds for i = 1:nznx
                x[i] = (vp[i] * vp[i] - 4.0 / 3.0 * vs[i] * vs[i]) * rho[i]
            end
        else
            @inbounds for i = 1:nznx
                x[i] = vp[i] * vp[i] * rho[i]
            end
        end
    elseif (s == :Zp)
        @inbounds for i = 1:nznx
            x[i] = vp[i] * rho[i]
        end
    elseif (s == :tau_epsilon)
        # memory optimize later
        tau_sigma, tau_epsilon = get_relaxation_times(
            mod.m[:Q],
            mod.ic[:nsls],
            mod.fc[:freqmin],
            mod.fc[:freqmax],
        )
        copyto!(x, tau_epsilon)
    elseif (s == :tau_sigma)
        # memory optimize later
        tau_sigma, tau_epsilon = get_relaxation_times(
            mod.m[:Q],
            mod.ic[:nsls],
            mod.fc[:freqmin],
            mod.fc[:freqmax],
        )
        copyto!(x, tau_sigma)
    else
        error(string("invalid attrib\t", s))
    end
    return x
end





"""
Copy medium parameters from Medium `mod`, to vector `x`. Return `x`.
```julia
copyto!(x,mod,[:vp, :vs])
copyto!(x,mod,[:vp, :rho])
```
"""
function Base.copyto!(x::AbstractArray{T}, mod::Medium, fields::Vector{Symbol}) where {T}
    nznx = prod(length.(mod.mgrid))
    (length(x) ≠ (count(fields .≠ :null) * nznx)) && error("size x")
    i0 = 0
    for field in fields
        if (field ≠ :null)
            xx = view(x, i0+1:i0+nznx)
            m = Base.getindex(mod, field)
            copyto!(xx, m)
            i0 += nznx
        end
    end
    return x
end


"""
Return a vector of medium parameters `x`, from Medium `mod`.
```julia
x=vec(mod, field)
```
"""
function Base.vec(mod::Medium, field::Symbol)
    # allocate
    rout = zeros(length(mod.mgrid[1]), length(mod.mgrid[2]))
    copyto!(rout, mod, [field])
    return rout
end

