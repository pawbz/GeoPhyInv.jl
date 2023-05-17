### A Pluto.jl notebook ###
# v0.19.21

using Markdown
using InteractiveUtils

# ╔═╡ ae7c128e-f0d0-4373-9497-60e2caefc24c
begin
    using Parameters, Statistics, PlutoUI
    using Interpolations: interpolate, extrapolate, Gridded, Linear, Flat
end

# ╔═╡ 8bb8e5e4-8603-45f0-bbf8-2aa2011c5aa9
begin
    abstract type MediumParameters end
    """
    Return medium property bounds based on maximum and minimum values of the array and frac.
    The bounds cannot be less than zero
    """
    function bounds(m::Array{T}, ref::T, frac) where {T}
        return [max(zero(eltype(m)), minimum(m) - convert(eltype(m), frac) * ref), maximum(m) + convert(eltype(m), frac) * ref]
    end
    for f in [:vp, :vs, :rho, :K, :M, :mu, :lambda, :invlambda, :invmu, :invrho, :invK]
        @eval(@with_kw_noshow mutable struct $f{T,N} <: MediumParameters
            m::Array{T,N}
            @assert all(m .> 0) "negative medium parameters"
            ref::T = Statistics.mean(m)
            bounds::Vector{T} = bounds(m, ref, 0.1)
        end)
        # if only one positional argument is given
        @eval $f(m) = $f(m=m)
        @eval $f{T,N}(m) where {T,N} = $f{T,N}(m=m)
    end
    abstract type Medium end

    mutable struct AcousticMedium{T,N} <: Medium
        mgrid::Vector{T1} where {T1<:StepRangeLen{Float64}}
        vp::vp{T,N}
        rho::rho{T,N}
        buffer::Array{T,N}
        function AcousticMedium{T,N}(mgrid, vp, rho, buffer) where {T,N}
            @assert all(size(vp.m) .== length.(mgrid))
            @assert all(size(rho.m) .== length.(mgrid))
            @assert all(size(buffer) .== length.(mgrid))
            new(map(StepRangeLen{Float64}, (mgrid)), vp, rho, buffer)
        end
    end

    mutable struct ElasticMedium{T,N} <: Medium
        mgrid::Vector{T1} where {T1<:StepRangeLen{Float64}}
        vp::vp{T,N}
        vs::vs{T,N}
        rho::rho{T,N}
        buffer::Array{T,N}
        function ElasticMedium{T,N}(mgrid, vp, vs, rho, buffer) where {T,N}
            @assert all(size(vp.m) .== length.(mgrid))
            @assert all(size(vs.m) .== length.(mgrid))
            @assert all(size(rho.m) .== length.(mgrid))
            @assert all(size(buffer) .== length.(mgrid))
            new(map(StepRangeLen{Float64}, (mgrid)), vp, vs, rho, buffer)
        end
    end
    # default vp, vs and rho
    function AcousticMedium{T,N}(mgrid) where {T,N}
        @assert length(mgrid) == N
        return AcousticMedium{T,N}(mgrid, vp(fill(T(2500.0), length.(mgrid)...)), rho(fill(T(2500.0), length.(mgrid)...)), zeros(T, length.(mgrid)...))
    end
    function ElasticMedium{T,N}(mgrid) where {T,N}
        @assert length(mgrid) == N
        return ElasticMedium{T,N}(mgrid, vp(fill(T(2500.0), length.(mgrid)...)), vs(fill(T(2000.0), length.(mgrid)...)), rho(fill(T(2500.0), length.(mgrid)...)), zeros(T, length.(mgrid)...))
    end
    AcousticMedium(mgrid) = AcousticMedium{Data.Number,length(mgrid)}(mgrid)
    ElasticMedium(mgrid) = ElasticMedium{Data.Number,length(mgrid)}(mgrid)
    AcousticMedium(mgrid, vp1, rho1) = AcousticMedium{eltype(vp1.m),length(mgrid)}(mgrid, vp1, rho1, zero(vp1.m))
    ElasticMedium(mgrid, vp1, vs1, rho1) = ElasticMedium{eltype(vp1.m),length(mgrid)}(mgrid, vp1, vs1, rho1, zero(vp1.m))
    MediumParameters(::AcousticMedium) = [:vp, :rho]
    MediumParameters(::ElasticMedium) = [:vp, :vs, :rho]

    for f in [:vp, :vs, :rho, :K, :M, :mu, :invmu, :lambda, :invlambda, :invrho, :invK]
        f1 = Symbol(f, "!")
        @eval(function $f1(m, medium::ElasticMedium)
            broadcast!(m, medium.vp.m, medium.vs.m, medium.rho.m) do vp1, vs1, rho1
                $f(vp1, vs1, rho1)
            end
            return m
        end)
        @eval(function $f(medium::ElasticMedium)
            m = zero(medium.vp.m)
            $f1(m, medium)
            return $f(m)
        end)
    end
    for f in [:vp, :rho, :K, :M, :mu, :lambda, :invlambda, :invrho, :invK]
        f1 = Symbol(f, "!")
        @eval(function $f1(m, medium::AcousticMedium)
            broadcast!(m, medium.vp.m, medium.rho.m) do vp1, rho1
                $f(vp1, rho1)
            end
            return m
        end)
        @eval(function $f(medium::AcousticMedium)
            m = zero(medium.vp.m)
            $f1(m, medium)
            return $f(m)
        end)
    end
    vp(vp, rho) = vp
    vp(vp, vs, rho) = vp
    vs(vp, vs, rho) = vs
    rho(vp, rho) = rho
    rho(vp, vs, rho) = rho
    K(vp, vs, rho) = (abs2(vp) - 4 / 3 * abs2(vs)) * rho
    K(vp, rho) = abs2(vp) * rho
    lambda(vp, rho) = abs2(vp) * rho
    lambda(vp, vs, rho) = (abs2(vp) - 2 * abs2(vs)) * rho
	invlambda(vp, vs, rho) = inv((abs2(vp) - 2 * abs2(vs)) * rho)
	invlambda(vp, rho) = inv(abs2(vp) * rho)
    M(vp, vs, rho) = abs2(vp) * rho
    M(vp, rho) = abs2(vp) * rho
    mu(vp, vs, rho) = abs2(vs) * rho
    mu(vp, rho) = zero(vp)
	invmu(vp, vs, rho) = inv(abs2(vs) * rho)
    invK(vp, vs, rho) = inv(K(vp, vs, rho))
    invK(vp, rho) = inv(K(vp, rho))
    invrho(vp, rho) = inv(rho)
    invrho(vp, vs, rho) = inv(rho)

end;

# ╔═╡ 6a2fff8e-cb1c-4aaf-9973-81bb457f600e
TableOfContents()

# ╔═╡ 2f2a3317-5a99-4ebd-b1cb-601123e6bcaa
"""
Return the number of dimensions of `Medium` 
"""
function Base.ndims(medium::Medium)
    return length(medium.mgrid)
end;

# ╔═╡ 5fab3bef-cc47-43ac-af58-94851311f09c
function Base.getindex(medium::Medium, s::Symbol)
    if (s ∈ MediumParameters(medium))
        return getfield(medium, s) # for vp, vs, and rho
    else
        eval(s)(medium) # derived fields
    end
end

# ╔═╡ 2cffedda-ce62-4328-bf66-450b28822758
md"## Base"

# ╔═╡ 3c29ef17-2e20-4a72-9211-0f002375bab0
Base.eltype(medium::Medium) = eltype(medium.vp.m)

# ╔═╡ a2b5a298-aca4-4bd0-b351-9534a7ea5e4f
function Base.copyto!(m1::AbstractArray, m::Medium, name)
    eval(Symbol(name, "!"))(m.buffer, m) # K!(array, m)
    copyto!(m1, m.buffer) # m1 can be on GPU
    return m1
end

# ╔═╡ 97109e57-0ba5-418b-bf73-c02050a29310
function Base.copyto!(m1::T, m::T) where {T<:MediumParameters}
    copyto!(m1.m, m.m)
    copyto!(m1.bounds, m.bounds)
    m1.ref = m.ref
    return m1
end

# ╔═╡ 16cc20c9-454d-454a-9e10-b8c0f9c44742
"""
Copy for `Medium`. The media should have same bounds and sizes.
Doesn't allocate any memory. Note that only [:vp, :rho, :vs] fields are stored in structs.
"""
function Base.copyto!(modo::T, mod::T) where {T<:Medium}
    @assert length.(mod.mgrid) == length.(modo.mgrid)
    for name in MediumParameters(mod)
        m1 = getfield(modo, name)
        m2 = getfield(mod, name)
        copyto!(m1, m2)
    end
    return modo
end;

# ╔═╡ 64065a7b-5bbd-4168-9eb1-fb5402076e59
md"## Medium Parameters"

# ╔═╡ c5dd4cd1-ace9-4e8c-b603-522f68450ea4
bounds(randn(3, 3), 0.1, 0.1)

# ╔═╡ b99d5714-f51a-451d-b96b-8a08866c3079
vp1 = vp(1000.0 .+ randn(3, 3))

# ╔═╡ cb5ed850-fd31-4b1f-8bc9-dc373bcd4932
vs1 = vs(500.0 .+ randn(3, 3))

# ╔═╡ 4bb4ce2c-f630-44f6-8685-6c2e2d4d5ee1
vp2 = vp(1500.0 .+ 10.0 * randn(3, 3))

# ╔═╡ 4de56401-d88d-444b-b264-c9aea237bc09
rho1 = rho{Float64,2}(1000.0 .+ randn(3, 3))

# ╔═╡ 8ce71374-9377-428d-b0e6-e7a75270cf97
rho2 = rho{Float32,2}(1000.0 .+ randn(3, 3))

# ╔═╡ 5623fd01-6ec9-4fa2-b998-7db4c61d8f7e
md"## Medium"

# ╔═╡ c1b6569e-d6d8-4c77-8cc4-1507e3e70603
medium1 = AcousticMedium([range(1, 2, length=3), range(1, 2, length=3)], vp1, rho1)

# ╔═╡ e405f8f5-4cf7-43ec-b534-a7f9b71d802b
medium2 = ElasticMedium([range(1, 2, length=3), range(1, 2, length=3)], vp1, vs1, rho1)

# ╔═╡ 801cab26-9939-41e7-aa6f-ba1ba57163ab
medium1D = AcousticMedium{Float32,1}([range(0, 100, length=10)])

# ╔═╡ d35fca7c-5cdb-46be-876e-9263f3ab8e3b
md"Replacing `vp` and `rho` fields in the medium is simple."

# ╔═╡ 1bd81ae0-021e-4bfc-9158-dfada19e7c00
medium1.vp = vp2

# ╔═╡ a074cbb1-c077-4ec6-aa79-3eb78116a8ba
md"## Derived Parameters"

# ╔═╡ ce6ae840-21f9-45f8-a2c1-4242435b8e69
K!(zeros(3, 3), medium1)

# ╔═╡ aa816da1-4b2f-4ace-807b-93d90293b3fd
rho!(zeros(3, 3), medium1)

# ╔═╡ 918a0fee-77d7-4e92-80d6-83a814e0f3d3
K(medium1)

# ╔═╡ 854debee-b6e4-4b0e-acde-fc09211669d0
invK(medium2)

# ╔═╡ c626f418-9a06-43b3-8081-b3c25be8e3a8
medium1[:K]

# ╔═╡ 20bd3b9f-4913-44d5-b652-3e0144f1b81a
medium1[:vp]

# ╔═╡ 8d5a4eb0-7b06-4db7-b458-dfbf5d27d5da
md"## Interpolations"

# ╔═╡ a3291f9d-88ee-4722-a5bd-a02bfae36005
function padarray!(modex::Medium, mod::Medium, npml, edges)
    number = eltype(mod)
    nex = length.(modex.mgrid)
    ist = [
        any(edges .== Symbol(string(dim), "min")) ? -npml + 1 : 1 for
        dim in dim_names(ndims(mod))
    ]
    vw = [ist[i]:ist[i]+nex[i]-1 for i = 1:ndims(mod)]
    pad = ImageFiltering.Pad(:replicate, fill(npml, ndims(mod))...)
    for mname in MediumParameters(mod)
        m1 = ImageFiltering.BorderArray(getfield(mod, mname).m, pad)
        m11 = view(m1, vw...)
        setfield!(modex, mname, eval(mname)(number.(m11)))
    end
    return modex
end

# ╔═╡ 05049ae0-3c72-4f22-aeac-3b5eacf5c0d2
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
            stop=mg[end] + maxpads[i] * step(mg),
            length=length(mg) + minpads[i] + maxpads[i],
        ) for (i, mg) in enumerate(mod.mgrid)
    ]
    modex = typeof(mod)(mgrid_new) # initialize a new medium
    padarray!(modex, mod, npml, edges)
end

# ╔═╡ ac6cec01-ce3d-4a92-842e-e124118ba34f
"""
Interpolate, extrapolate or truncate model to a new grid.
Check for memory allocations, specially when using Interpolations library.
Does this library contain inplace methods?
"""
function update!(modex::T, mod::T) where {T<:Medium}
    number = eltype(mod)
    for mname in MediumParameters(mod)
        itp = extrapolate(
            interpolate(Tuple(mod.mgrid), getfield(mod, mname).m, Gridded(Linear())),
            Flat(),
        )
        setfield!(modex, mname, eval(mname)(number.(itp[modex.mgrid...])))
    end
    return modex
end

# ╔═╡ eae711d0-f690-4b95-baeb-667eb48659c7
md"## Addon"

# ╔═╡ 9a8c7ee7-ca02-48bb-86bf-bb2afc3e5edc
# perturb an array m by percentage of m0, when mask is true
function perturb!(m, m0, perc, mask::Function)
    for i in CartesianIndices(m)
        ii = Tuple(i)
        if (mask(ii))
            m[i] = m[i] + perc * 1e-2 * m0
        end
    end
end

# ╔═╡ 196ca0eb-f1d3-44ab-9f21-58b9ed338741
function perturb!(m, m0, perc, idx::CartesianIndices)
    for i in idx
        m[i] = m[i] + perc * 1e-2 * m0
    end
end

# ╔═╡ c9ac2faa-492f-454f-a9a0-62075c27a2d4
"""
In-place method to add a rectangular box or cuboid in the medium. Its position is specified by its corners  
and perc denotes the percentage of the perturbation.

```julia
update!(medium, field; rectangle=[[0,10],[10,10]], perc=10) # for 2D medium
update!(medium, field; rectangle=[[0,10,15],[10,10,15]], perc=10) # for 3D medium
```

In-place method to add random noise to medium.

```julia
update!(medium, fields; randn_perc=10)
```
"""
function update!(medium::Medium, fields::Vector{Symbol}; rectangle=nothing,  perc=0.0, randn_perc=0.0)
	
    number = eltype(medium)
    function χ(m::T, m0::T, flag::Int64=1) where {T<:Real}
        if (flag == 1)
            return (m - m0) * inv(m0)
        else
            (flag == -1)
            return (m * m0 + m0)
        end
    end


    mgrid = medium.mgrid
    if (!(rectangle === nothing))
        @assert length(rectangle) == 2
        @assert length.(rectangle) == fill(ndims(medium), 2)
        rectangle_indices = CartesianIndices(
            Tuple([
                argmin(abs.(mg .- rectangle[1][im])):argmin(abs.(mg .- rectangle[2][im])) for (im, mg) in enumerate(mgrid)
            ]),
        )
        if (!iszero(perc))
            for field in fields
                m = getfield(medium, field).m
                m0 = getfield(medium, field).ref
                perturb!(m, m0, perc, rectangle_indices)
            end
        end
    end
    if (!iszero(randn_perc))
        for field in fields
            m = medium[field].m
            m0 = medium[field].ref
            broadcast!(m, m) do mi
                χ((χ(mi, m0, 1) + randn(number) * number(randn_perc * 0.01)), m0, -1)
            end
            setfield!(medium, field, eval(field)(m))
        end
    end

    return medium
end

# ╔═╡ 67f40155-47f1-471a-84b1-79950d91f2bd
"""
Similar to previous method, but with initialization.
"""
function update(
    mod::Medium,
    mgrid::Vector{T} where {T<:StepRangeLen{Float64}}
)
    modex = typeof(mod)(mgrid) # initialize a new medium
    update!(modex, mod)
    return modex
end

# ╔═╡ c23b328b-003c-4cf6-81c4-28b0fa2dd2bf
update(medium1, [range(-10, 200, length=100), range(-10, 200, length=100)])

# ╔═╡ d06d5b8a-e142-4bde-b1a0-d8eaa0186ffb
update!(medium1, [:vp]; randn_perc=5)

# ╔═╡ 22a591c4-9bd4-4b13-83ca-ba1ae67252b7
update!(medium1D, [:vp]; randn_perc=5)

# ╔═╡ f4ef3bba-8063-4855-985e-157af7c0c67f
md"## Appendix"

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
Parameters = "d96e819e-fc66-5662-9728-84c9c7592b0a"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
Interpolations = "~0.14.7"
Parameters = "~0.12.3"
PlutoUI = "~0.7.51"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.0"
manifest_format = "2.0"
project_hash = "f244600bde8b993e933c646b055cb9186833a626"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "76289dc51920fdc6e0013c872ba9551d54961c24"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.6.2"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "e30f2f4e20f7f186dc36529910beaedc60cfa644"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.16.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "7a60c856b9fa189eb34f5f8a6f6b5529b7942957"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.2+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "721ec2cf720536ad005cb38f50dbba7b02419a15"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.14.7"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "82d7c9e310fe55aa54996e6f7f94674e2a38fcb4"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.9"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OrderedCollections]]
git-tree-sha1 = "d321bf2de576bf25ec4d3e4360faca399afca282"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7302075e5e06da7d000d9bfa055013e3e85578ca"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.9"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "b478a748be27bd2f2c73a7690da219d0844db305"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.51"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "259e206946c293698122f63e2b513a7c99a244e8"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.1.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "7eb1686b4f04b82f96ed7a4ea5890a4f0c7a09f1"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "6d7bb727e76147ba18eed998700998e17b8e4911"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.4"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "8982b3607a212b070a5e46eea83eb62b4744ae12"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.25"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[deps.URIs]]
git-tree-sha1 = "074f993b0ca030848b897beff716d93aca60f06a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.7.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╠═8bb8e5e4-8603-45f0-bbf8-2aa2011c5aa9
# ╠═6a2fff8e-cb1c-4aaf-9973-81bb457f600e
# ╟─2cffedda-ce62-4328-bf66-450b28822758
# ╠═2f2a3317-5a99-4ebd-b1cb-601123e6bcaa
# ╠═5fab3bef-cc47-43ac-af58-94851311f09c
# ╠═3c29ef17-2e20-4a72-9211-0f002375bab0
# ╠═a2b5a298-aca4-4bd0-b351-9534a7ea5e4f
# ╠═97109e57-0ba5-418b-bf73-c02050a29310
# ╠═16cc20c9-454d-454a-9e10-b8c0f9c44742
# ╟─64065a7b-5bbd-4168-9eb1-fb5402076e59
# ╠═c5dd4cd1-ace9-4e8c-b603-522f68450ea4
# ╠═b99d5714-f51a-451d-b96b-8a08866c3079
# ╠═cb5ed850-fd31-4b1f-8bc9-dc373bcd4932
# ╠═4bb4ce2c-f630-44f6-8685-6c2e2d4d5ee1
# ╠═4de56401-d88d-444b-b264-c9aea237bc09
# ╠═8ce71374-9377-428d-b0e6-e7a75270cf97
# ╟─5623fd01-6ec9-4fa2-b998-7db4c61d8f7e
# ╠═c1b6569e-d6d8-4c77-8cc4-1507e3e70603
# ╠═e405f8f5-4cf7-43ec-b534-a7f9b71d802b
# ╠═801cab26-9939-41e7-aa6f-ba1ba57163ab
# ╟─d35fca7c-5cdb-46be-876e-9263f3ab8e3b
# ╠═1bd81ae0-021e-4bfc-9158-dfada19e7c00
# ╟─a074cbb1-c077-4ec6-aa79-3eb78116a8ba
# ╠═ce6ae840-21f9-45f8-a2c1-4242435b8e69
# ╠═aa816da1-4b2f-4ace-807b-93d90293b3fd
# ╠═918a0fee-77d7-4e92-80d6-83a814e0f3d3
# ╠═854debee-b6e4-4b0e-acde-fc09211669d0
# ╠═c626f418-9a06-43b3-8081-b3c25be8e3a8
# ╠═20bd3b9f-4913-44d5-b652-3e0144f1b81a
# ╟─8d5a4eb0-7b06-4db7-b458-dfbf5d27d5da
# ╠═a3291f9d-88ee-4722-a5bd-a02bfae36005
# ╠═05049ae0-3c72-4f22-aeac-3b5eacf5c0d2
# ╠═ac6cec01-ce3d-4a92-842e-e124118ba34f
# ╠═67f40155-47f1-471a-84b1-79950d91f2bd
# ╠═c23b328b-003c-4cf6-81c4-28b0fa2dd2bf
# ╟─eae711d0-f690-4b95-baeb-667eb48659c7
# ╠═9a8c7ee7-ca02-48bb-86bf-bb2afc3e5edc
# ╠═196ca0eb-f1d3-44ab-9f21-58b9ed338741
# ╠═c9ac2faa-492f-454f-a9a0-62075c27a2d4
# ╠═d06d5b8a-e142-4bde-b1a0-d8eaa0186ffb
# ╠═22a591c4-9bd4-4b13-83ca-ba1ae67252b7
# ╟─f4ef3bba-8063-4855-985e-157af7c0c67f
# ╠═ae7c128e-f0d0-4373-9497-60e2caefc24c
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
