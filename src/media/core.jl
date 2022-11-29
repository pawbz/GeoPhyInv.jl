"""
Mutable type for storing medium parameters.
```julia
mod=Medium(mgrid, names=[:vp,:vs,:rho]) # initiate elastic medium
mod=Medium(mgrid, names=[:vp,:rho]) # initiate acoustic medium
```
This initializes a subsurface model with `:vp`, `:rho` and `:vs` parameters. Print the
names of the medium parameters stored in `mod`.
```julia
names(mod)
```

## Indexing
* `mod.mgrid` : returns the spatial-grid bundle
* `mod[:vp]` : P-wave velocity
* `mod[:vs]` : S-wave velocity
* `mod[:rho]` : mass density
* `mod[:Zp]` : P-wave impedance 
* `mod[:K]` : bulk modulus (unrelaxed when considering attenuation) 
* `mod[:M]` : P-wave modulus
* `mod[:mu]` : shear modulus
* `mod[:Q]` : quality factor (relaxation times are optimized to be constant over all frequencies; see Robertsson, et. al, 1994)
* `mod.ref` : reference medium parameters 
* `mod.bounds` : bounds of medium parameters

"""
mutable struct Medium{N} # N is 2 for 2D, and 3 for 3D
    mgrid::Vector{T} where {T<:StepRangeLen{Float64}}
    m::NamedStack{Array{Float64,N}}
    bounds::NamedArrays.NamedArray{
        Array{Float64,1},
        1,
        Array{Array{Float64,1},1},
        Tuple{OrderedCollections.OrderedDict{Symbol,Int64}},
    }
    ref::NamedArray{
        Float64,
        1,
        Array{Float64,1},
        Tuple{OrderedCollections.OrderedDict{Symbol,Int64}},
    }
    # store some floats specific to the medium
    fc::NamedArrays.NamedArray{
        Float64,
        1,
        Array{Float64,1},
        Tuple{OrderedCollections.OrderedDict{Symbol,Int64}},
    }
    # store some integer constants
    ic::NamedArrays.NamedArray{
        Int64,
        1,
        Array{Int64,1},
        Tuple{OrderedCollections.OrderedDict{Symbol,Int64}},
    }
    # store 3D arrays, if any (currently used for attenuation, when N=2)
    m3::NamedStack{Array{Float64,3}}
end

# default construction without fcs and ics, and m3
function Medium(mgrid, m::NamedStack{Array{Float64,N}}, bounds, ref) where {N}
    return Medium(
        map(StepRangeLen{Float64},(mgrid)),
        Array{Float64}.(m),
        Array{Float64}.(bounds),
        Float64.(ref),
        NamedArray(Float64.([0.0]), ([:o],)),
        NamedArray([0], ([:o],)),
        NamedArray([zeros(Float64,1, 1, 1)], ([:o],)),
    )
end

function NamedArrays.names(mod::Medium)
    return names(mod.m)
end


include("attenuation.jl")
include("base.jl")
include("gallery.jl")



#=


* `vp0::Vector{Float64}` : [vpmin, vpmax]
* `vs0::Vector{Float64}` : [vsmin, vsmax]
* `rho0::Vector{Float64}` : [rhomin, rhomax]
* `χvp::Array{Float64,2}` : two-dimensional contrast model (χ) for vp, for e.g., zeros(length(mgrid[1]), length(mgrid[2]))
* `χvs::Array{Float64}` : two-dimensional contrast model (χ) for vs, for e.g., zeros(length(mgrid[1]), length(mgrid[2]))
* `χrho::Array{Float64}` : two-dimensional contrast model (χ) for density, for e.g., zeros(length(mgrid[1]), length(mgrid[2]))
* `mgrid` : array of ranges to determine the dimensions of models
"""
#
#	vp0::Vector{Float64}
#	vs0::Vector{Float64}
#	rho0::Vector{Float64}
#	χvp::Array{Float64,2}
#	χvs::Array{Float64,2}
#	χrho::Array{Float64,2}
#	# all derived fields from previous
#	K0::Vector{Float64}
#	mu0::Vector{Float64}
#	KI0::Vector{Float64}
#	muI0::Vector{Float64}
#	rhoI0::Vector{Float64}
#	"adding conditions that are to be false while construction"
#	Medium(vp0,vs0,rho0,χvp,χvs,χrho,mgrid,K0,mu0,KI0,muI0,rhoI0,ref) = 
#		any([
#		     any(vp0.<0.0), 
#		     any(vs0.<0.0),
#		     any(rho0.<0.0),
#		     #all([all(vp0 .≠ 0.0), any(χ(χvp,vp0,-1) .< vp0[1])]), # check vp bounds
#		     #all([all(vp0 .≠ 0.0), any(χ(χvp,vp0,-1) .> vp0[2])]), # check vp bounds
#		     #all([all(vs0 .≠ 0.0), any(χ(χvs,vs0,-1) .< vs0[1])]), # check vs bounds
#		     #all([all(vs0 .≠ 0.0), any(χ(χvs,vs0,-1) .> vs0[2])]), # check vs bounds
#		     #all([all(rho0 .≠ 0.0), any(χ(χrho,rho0,-1) .< rho0[1])]), # check rho bounds
#		     #all([all(rho0 .≠ 0.0), any(χ(χrho,rho0,-1) .> rho0[2])]), # check rho bounds
#		     size(χvp) != (length(mgrid[1]), length(mgrid[2])), # dimension check
#		     size(χvs) != (length(mgrid[1]), length(mgrid[2])), # dimension check
#		     size(χrho) != (length(mgrid[1]), length(mgrid[2])) # dimension check
#		    ]) ? 
#		       error("error in Medium construction") : new(vp0,vs0,rho0,χvp,χvs,χrho,mgrid,K0,mu0,KI0,muI0,rhoI0,ref)
#end
#
# check whether a seismic model is bounded
function isbounded(mod::Medium)
	@warn "need to be implemented"
	result=true
end

# method to create Medium object without using derived fields
function Medium(vp0, vs0, rho0, χvp, χvs, χrho, 
		 mgrid::Vector{StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}})
	mod=Medium(vp0, vs0, rho0, χvp, χvs, χrho, mgrid,
	    zeros(2), zeros(2), zeros(2), zeros(2), zeros(2),
	    Medium_ref(zeros(8)...))
	update_derived!(mod)
	return mod
end

=#

# update the derived fields of the Seisimc model
function update!(mod::Medium)
    "note that: mean(K0) ≠ 1/mean(KI0)"
    mb = deepcopy(mod.bounds)

    rho = mb[:rho]
    rhoI = reverse(inv.(rho), dims = 1)

    if (:vs ∈ names(mb)[1])
        K = (mb[:vp] .* mb[:vp] .- 4.0 / 3.0 .* mb[:vs] .* mb[:vs]) .* mb[:rho]
        KI = reverse(inv.(K), dims = 1)
        newnames = [:K, :KI, :rhoI, :mu, :muI]
        mu = mb[:vs] .* mb[:vs] .* mb[:rho]
        muI = reverse(inv.(mu), dims = 1)
        mbnew = NamedArray([K, KI, rhoI, mu, muI], (newnames,))
    else
        K = (mb[:vp] .* mb[:vp]) .* mb[:rho]
        KI = reverse(inv.(K), dims = 1)
        newnames = [:K, :KI, :rhoI]
        mbnew = NamedArray([K, KI, rhoI], (newnames,))
    end


    for s in newnames
        # println(typeof(mod.bounds), '\'NamedArray([[0.0,0.0]],[[s],]))
        (s ∉ names(mod.bounds)[1]) &&
            (mod.bounds = Array{Float64}.(vcat(mod.bounds, NamedArray([[0.0, 0.0]], [s]))))
    end

    for s in newnames
        copyto!(mod.bounds[s], mbnew[s])
    end

    for s in names(mod.bounds)[1]
        meanvalue = Statistics.mean(mod.bounds[s])
        if (s ∈ names(mod.ref)[1])
            mod.ref[s] = meanvalue
        else
            mod.ref = vcat(mod.ref, NamedArray([meanvalue], [s]))
        end
    end

    return nothing
end


"""
Return medium property bounds based on maximum and minimum values of the array and frac.
The bounds cannot be less than zero
"""
function bounds(mod::AbstractArray, frac::Real = 0.1)
    any(mod .< 0.0) && error("model values less than zero")
    bounds = zeros(2)
    bound = frac * Statistics.mean(mod)
    bounds[1] = ((minimum(mod) - bound) < 0.0) ? 0.0 : (minimum(mod) - bound)
    bounds[2] = maximum(mod) + bound
    return bounds
end



function update!(mod::Medium, fields::AbstractVector{Symbol}, bounds)
    @assert !(bounds === nothing)
    for name in names(mod.bounds)[1]
        if (name ∈ fields)
            i = findall(x -> x == name, fields)[1]
            copyto!(mod.bounds[name], Float64.(bounds[i]))
        end
    end
    update!(mod)
end

"""
Adjust the bounds and hence the reference values.
Since the reference values are adjust the χ fields should also be changed
"""
function update!(mod::Medium, frac::Real)# 0.1
    for name in names(mod.m)[1]
        b = bounds(mod.m[name], frac)
        (name ∈ names(mod.bounds)[1]) && copyto!(mod.bounds[name], Float64.(b))
    end
    update!(mod)
    return mod
end


"""
Return dimensionless contrast model parameter using the
reference value.

# Arguments
* `mod::Array{Float64}` : subsurface parameter
* `mod0::Vector{Float64}` : reference value is mean of this vector
* `flag::Int64=1` : 
"""
function χ(mod::AbstractArray{T}, mod0::T, flag::Int64 = 1) where {T<:Real}
    mod_out = copy(mod)
    χ!(mod_out, mod0, flag)
    return mod_out
end
function χ!(mod::AbstractArray{T}, mod0::T, flag::Int64 = 1) where {T<:Real}
    m0 = mod0
    if (flag == 1)
        if (iszero(m0))
            return nothing
        else
            @inbounds for i in eachindex(mod)
                mod[i] = ((mod[i] - m0) * inv(m0))
            end
            return mod
        end
    elseif (flag == -1)
        @inbounds for i in eachindex(mod)
            mod[i] = mod[i] * m0 + m0
        end
        return mod
    end
end

function χ(m::T, m0::T, flag::Int64 = 1) where {T<:Real}
    if (flag == 1)
        return (m - m0) * inv(m0)
    else
        (flag == -1)
        return (m * m0 + m0)
    end
end

"""
Gradients
Return contrast model parameter using the
reference value.
"""
function χg(mod::AbstractArray{T}, mod0::T, flag::Int64 = 1) where {T<:Real}
    mod_out = copy(mod)
    χg!(mod_out, mod0, flag)
    return mod_out
end
function χg!(mod::AbstractArray{T}, mod0::T, flag::Int64 = 1) where {T<:Real}
    if (flag == 1)
        rmul!(mod, mod0)
    elseif (flag == -1)
        m0 = inv(mod0)
        rmul!(mod, m0)
    end
    return mod
end

function χg(m::T, m0::T, flag::Int64 = 1) where {T<:Real}
    if (flag == 1)
        return m * m0
    else
        (flag == -1)
        return m * inv(m0)
    end
end





include("addons.jl")

include("chainrule.jl")


"""
Interpolate, extrapolate or truncate model to a new grid.
Check for memory allocations, specially when using Interpolations library.
Does this library contain inplace methods?
"""
function update!(modex::Medium{N}, mod::Medium{N}) where {N}
    for m in names(mod.m)[1]
        itp = extrapolate(
            interpolate(Tuple(mod.mgrid), mod[m], Gridded(Linear())),
            Interpolations.Flat(),
        )
        copyto!(modex[m], itp[modex.mgrid...])
    end
    update!(modex, 0.1)
    return modex
end

"""
Similar to previous method, but with initialization.
"""
function update(
    mod::Medium,
    mgrid::Vector{T} where {T<:StepRangeLen{Float64}}
)
    modex = Medium(mgrid, names(mod)[1]) # initialize a new medium
    update!(modex, mod)
    return modex
end




function padarray!(modex::Medium{N}, mod::Medium{N}, npml, edges) where {N}
    nex = length.(modex.mgrid)
    ist = [
        any(edges .== Symbol(string(dim), "min")) ? -npml + 1 : 1 for
        dim in dim_names(ndims(mod))
    ]
    vw = [ist[i]:ist[i]+nex[i]-1 for i = 1:ndims(mod)]
    pad = ImageFiltering.Pad(:replicate, fill(npml, ndims(mod))...)
    for m in names(mod.m)[1]
        m1 = ImageFiltering.BorderArray(mod[m], pad)
        m11 = view(m1, vw...)
        copyto!(modex[m], m11)
    end
    update!(modex, 0.1)
    return modex
end

# edges contains :xmin, :xmax, ... etc
function padarray(mod::Medium, npml, edges)
    minpads = [
        any(edges .== Symbol(string(dim), "min")) ? npml : 0 for
        dim in dim_names(ndims(mod))
    ]
    maxpads = [
        any(edges .== Symbol(string(dim), "max")) ? npml : 0 for
        dim in dim_names(ndims(mod))
    ]
    mgrid_new = [
        range(
            mg[1] - minpads[i] * step(mg),
            stop = mg[end] + maxpads[i] * step(mg),
            length = length(mg) + minpads[i] + maxpads[i],
        ) for (i, mg) in enumerate(mod.mgrid)
    ]
    modex = Medium(mgrid_new, names(mod)[1]) # initialize a new medium
    padarray!(modex, mod, npml, edges)
end



"""
Apply smoothing to `Medium` using a Gaussian filter of zwidth and xwidth

# Arguments

* `mod::Medium` : argument that is modified
* `zperc::Real` : smoothing percentage in z-direction
* `xperc::Real=zperc` : smoothing percentage in x-direction

# Keyword Arguments

* `zmin::Real=mod.mgrid[1][1]` : 
* `zmax::Real=mod.mgrid[1][end]` : 
* `xmin::Real=mod.mgrid[2][1]` : 
* `xmax::Real=mod.mgrid[2][end]` : 
* `fields` : fields of seismic model that are to be smooth
"""
function Medium_smooth(
    mod::Medium,
    zperc::Real,
    xperc::Real = zperc;
    zmin::Real = mod.mgrid[1][1],
    zmax::Real = mod.mgrid[1][end],
    xmin::Real = mod.mgrid[2][1],
    xmax::Real = mod.mgrid[2][end],
    fields = [:vp, :rho],
)
    xwidth = Float64(xperc) * 0.01 * abs(mod.mgrid[2][end] - mod.mgrid[2][1])
    zwidth = Float64(zperc) * 0.01 * abs(mod.mgrid[1][end] - mod.mgrid[1][1])
    xnwin = Int(div(xwidth, step(mod.mgrid[2]) * 2.0))
    znwin = Int(div(zwidth, step(mod.mgrid[1]) * 2.0))

    izmin = Interpolation.indminn(mod.mgrid[1], Float64(zmin), 1)[1]
    izmax = Interpolation.indminn(mod.mgrid[1], Float64(zmax), 1)[1]
    ixmin = Interpolation.indminn(mod.mgrid[2], Float64(xmin), 1)[1]
    ixmax = Interpolation.indminn(mod.mgrid[2], Float64(xmax), 1)[1]

    # @warn "check this routine, smooth contrast values instead?"

    modg = deepcopy(mod)
    for (i, iff) in enumerate(fields)
        m = view(mod[iff], izmin:izmax, ixmin:ixmax)
        mg = view(modg[iff], izmin:izmax, ixmin:ixmax)
        imfilter!(mg, m, Kernel.gaussian([znwin, xnwin]))
    end
    return modg
end



"""
function to resample in the model domain

# Arguments
* `mod::Medium` : model
* `modi::Medium` : model after interpolation
"""
function interp_spray!(
    mod::Medium,
    modi::Medium,
    attrib::Symbol,
    Battrib::Symbol = :B2,
    fields = [:vp, :rho];
    pa = nothing,
)
    if (pa === nothing)
        pa = Interpolation.Kernel(
            [mod.mgrid[2], mod.mgrid[1]],
            [modi.mgrid[2], modi.mgrid[1]],
            Battrib,
        )
    end

    # :vs is missing here....
    "loop over fields in `Medium`, add vs later"
    for field in fields
        Interpolation.interp_spray!(mod[field], modi[field], pa, attrib)
    end
end

#
#macro inter

function save(mod::Medium, folder; N = 100)
    !(isdir(folder)) && error("invalid directory")
    error("need to be updated")

    nx = length(mod.mgrid[2])
    nz = length(mod.mgrid[1])
    n = max(nx, nz)
    fact = (n > N) ? round(Int, n / N) : 1
    #mgrid=resamp(mod.mgrid, step(mod.mgrid[2])*fact, step(mod.mgrid[1])*fact)
    x = mgrid[2]
    z = mgrid[1]
    nx = length(x)
    nz = length(z)
    modo = Medium_zeros(mgrid)
    adjust_bounds!(modo, mod)
    interp_spray!(mod, modo, :interp)

    for m in [:vp, :rho, :Zp]
        # save original gf
        file = joinpath(folder, string("im", m, ".csv"))
        CSV.write(
            file,
            DataFrame(hcat(repeat(z, outer = nx), repeat(x, inner = nz), vec(modo[m]))),
        )
    end
end




