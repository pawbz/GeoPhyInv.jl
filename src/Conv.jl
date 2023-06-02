### A Pluto.jl notebook ###
# v0.19.21

using Markdown
using InteractiveUtils

# ╔═╡ a9dd629b-e19e-464f-a9e6-a9675f9ee4ea
# ╠═╡ skip_as_script = true
#=╠═╡
using IterativeSolvers, PlutoTest, PlutoUI
  ╠═╡ =#

# ╔═╡ f01b439a-5641-488c-b175-a37da23dd94c
module Conv
using FFTW
using DSP
using LinearAlgebra
using LinearMaps


using DSP: nextfastfft
struct D end
struct G end
struct S end

"""
Model : d = convolution(g, s)
d, g and s can have arbitrary +ve and -ve lags
"""
mutable struct Pconv{T<:Real,Nd,Ng,Ns}
    dsize::NTuple{Nd,Int}
    gsize::NTuple{Ng,Int}
    ssize::NTuple{Ns,Int}
    "length after zero padding"
    np2::Int # nfastfft length
    d::AbstractArray{T,Nd}
    g::AbstractArray{T,Ng}
    s::AbstractArray{T,Ns}
    dpad::Array{T,Nd}
    gpad::Array{T,Ng}
    spad::Array{T,Ns}
    dfreq::Array{Complex{T},Nd}
    gfreq::Array{Complex{T},Ng}
    sfreq::Array{Complex{T},Ns}
    dfftp::FFTW.rFFTWPlan
    difftp::FFTW.Plan
    gfftp::FFTW.rFFTWPlan
    gifftp::FFTW.Plan
    sfftp::FFTW.rFFTWPlan
    sifftp::FFTW.Plan
    "+ve and -ve lags of g"
    glags::NTuple{2,Int}
    "+ve and -ve lags of s"
    slags::NTuple{2,Int}
    "+ve and -ve lags of d"
    dlags::NTuple{2,Int}
end

function Pconv(T=Float64;
    dsize=(1,),
    ssize=(1,),
    gsize=(1,),
    d=zeros(T, dsize...),
    s=zeros(T, ssize...),
    g=zeros(T, gsize...),
    slags=nothing,
    dlags=nothing,
    glags=nothing,
    np2=nextfastfft(ssize[1] + gsize[1] - 1), # zero-padded fft dimension for the plan, such that circular convolution approximates linear convolution
    fftwflag=FFTW.ESTIMATE
)
    @assert size(s) == ssize
    @assert size(g) == gsize
    @assert size(d) == dsize

    Nd = length(dsize)
    Ns = length(ssize)
    Ng = length(gsize)
    nts = ssize[1]
    ntg = gsize[1]
    ntd = dsize[1]

    # default lags
    if (slags === nothing)
        # equal +ve and -ve lags for s
        nwplags = div(nts - 1, 2)
        nwnlags = nts - 1 - nwplags
        slags = (nwplags, nwnlags)
    end
    if (dlags === nothing)
        # no negative lags
        nsplags = ntd - 1
        nsnlags = ntd - 1 - nsplags
        dlags = (nsplags, nsnlags)
    end
    if (glags === nothing)
        # no negative lags
        nrplags = ntg - 1
        nrnlags = ntg - 1 - nrplags
        glags = (nrplags, nrnlags)
    end


    # check if nlags are consistent with the first dimension of inputs
    (sum(glags) + 1 ≠ size(g, 1)) && error("glags")
    (sum(dlags) + 1 ≠ size(d, 1)) && error("dlags")
    (sum(slags) + 1 ≠ size(s, 1)) && error("slags")

    FFTW.set_num_threads(Sys.CPU_THREADS)
    nrfft = div(np2, 2) + 1


    dfftp = plan_rfft(zeros(T, np2, dsize[2:end]...), [1], flags=fftwflag, timelimit=Inf)
    difftp = plan_irfft(complex.(zeros(T, nrfft, dsize[2:end]...)), np2, [1], flags=fftwflag, timelimit=Inf)

    if (gsize[2:end] ≠ dsize[2:end])
        gfftp = plan_rfft(zeros(T, np2, gsize[2:end]...), [1], flags=fftwflag, timelimit=Inf)
        gifftp = plan_irfft(complex.(zeros(T, nrfft, gsize[2:end]...)), np2, [1], flags=fftwflag, timelimit=Inf)
    else
        gfftp = dfftp
        gifftp = difftp
    end

    if (ssize[2:end] ≠ dsize[2:end])
        sfftp = plan_rfft(zeros(T, np2, ssize[2:end]...), [1], flags=fftwflag, timelimit=Inf)
        sifftp = plan_irfft(complex.(zeros(T, nrfft, ssize[2:end]...)), np2, [1], flags=fftwflag, timelimit=Inf)
    else
        sfftp = dfftp
        sifftp = difftp

    end


    # preallocate freq domain vectors after rfft
    dfreq = complex.(zeros(T, nrfft, dsize[2:end]...))
    gfreq = complex.(zeros(T, nrfft, gsize[2:end]...))
    sfreq = complex.(zeros(T, nrfft, ssize[2:end]...))

    # preallocate padded arrays
    dpad = (zeros(T, np2, dsize[2:end]...))
    gpad = (zeros(T, np2, gsize[2:end]...))
    spad = (zeros(T, np2, ssize[2:end]...))

    pa = Pconv(dsize, gsize, ssize, np2, d, g, s,
        dpad, gpad, spad,
        dfreq, gfreq, sfreq,
        dfftp, difftp,
        gfftp, gifftp,
        sfftp, sifftp,
        glags, slags, dlags)
end

"""
Methods to perform zero padding and truncation.

# Arguments

* `x` : real signal with dimension nplag + nnlag + 1
	first has decreasing negative lags, 
	then has zero lag at nnlags + 1,
	then has increasing positive nplags lags,
	signal contains only positive lags and zero lag if nnlag=0 and vice versa
* `xpow2` : npow2 real vector with dimension npow2
* `nplags` : number of positive lags
* `nnlags` : number of negative lags
* `npow2` : number of samples in xpow2
* `pad` means xpow2 is returned using x
* `truncate` means x is returned using xpow2
"""
function pad!(
    x::AbstractArray{T},
    xpow2::AbstractArray{T},
    nplags::Integer,
    nnlags::Integer,
    npow2::Integer,
) where {T}
    @assert size(x, 1) == nplags + nnlags + 1
    @assert size(xpow2, 1) == npow2

    for id in 1:size(x, 2)
        xpow2[1, id] = (x[nnlags+1, id]) # zero lag
        # +ve lags
        if (nplags > 0)
            for i = 1:nplags
                @inbounds xpow2[i+1, id] = (x[nnlags+1+i, id])
            end
        end
        # -ve lags
        if (nnlags != 0)
            for i = 1:nnlags
                @inbounds xpow2[npow2-i+1, id] = (x[nnlags+1-i, id])
            end
        end
    end
    return nothing
end


function truncate!(
    x::AbstractArray{T},
    xpow2::AbstractArray{T},
    nplags::Integer,
    nnlags::Integer,
    npow2::Integer,
) where {T}
    (size(x, 1) ≠ nplags + nnlags + 1) && error("size x")
    (size(xpow2, 1) ≠ npow2) && error("size xpow2")

    for id in 1:size(x, 2)
        x[nnlags+1, id] = (xpow2[1, id]) # zero lag
        if (nplags != 0)
            for i = 1:nplags
                @inbounds x[nnlags+1+i, id] = (xpow2[1+i, id])
            end
        end
        if (nnlags != 0)
            for i = 1:nnlags
                @inbounds x[nnlags+1-i, id] = (xpow2[npow2-i+1, id])
            end
        end
    end
    return nothing
end


function initialize_d!(pa::Pconv, d=pa.d)
    T = eltype(pa.d)
    fill!(pa.dfreq, complex(T(0)))
    fill!(pa.dpad, T(0))
    pad!(d, pa.dpad, pa.dlags[1], pa.dlags[2], pa.np2)
end

function initialize_g!(pa::Pconv, g=pa.g)
    T = eltype(pa.g)
    fill!(pa.gfreq, complex(T(0)))
    fill!(pa.gpad, T(0))
    pad!(g, pa.gpad, pa.glags[1], pa.glags[2], pa.np2)
end


function initialize_s!(pa::Pconv, s=pa.s)
    T = eltype(pa.s)
    fill!(pa.sfreq, complex(T(0)))
    fill!(pa.spad, T(0))
    pad!(s, pa.spad, pa.slags[1], pa.slags[2], pa.np2)
end

# initialize freq vectors
function initialize_all!(pa::Pconv, d=pa.d, g=pa.g, s=pa.s)
    initialize_d!(pa, d)
    initialize_g!(pa, g)
    initialize_s!(pa, s)
end

"""
Convolution modelling with no allocations at all.
By default, the fields `g`, `d` and `s` in pa are modified accordingly.
Otherwise, use keyword arguments to input them.
"""
function mod!(pa::Pconv{T,N,N,N}, ::D; g=pa.g, d=pa.d, s=pa.s) where {N,T<:Real}
    initialize_all!(pa, d, g, s)
    mul!(pa.sfreq, pa.sfftp, pa.spad)
    mul!(pa.gfreq, pa.gfftp, pa.gpad)
    for i in eachindex(pa.dfreq)
        @inbounds pa.dfreq[i] = pa.gfreq[i] * pa.sfreq[i]
    end
    mul!(pa.dpad, pa.difftp, pa.dfreq)
    truncate!(d, pa.dpad, pa.dlags[1], pa.dlags[2], pa.np2)
    return pa
end

function mod!(pa::Pconv{T,N,N,N}, ::G; g=pa.g, d=pa.d, s=pa.s) where {N,T<:Real}
    initialize_all!(pa, d, g, s)
    mul!(pa.sfreq, pa.sfftp, pa.spad)
    mul!(pa.dfreq, pa.dfftp, pa.dpad)
    conj!(pa.sfreq)
    for i in eachindex(pa.dfreq)
        @inbounds pa.gfreq[i] = pa.dfreq[i] * pa.sfreq[i]
    end
    mul!(pa.gpad, pa.gifftp, pa.gfreq)
    truncate!(g, pa.gpad, pa.glags[1], pa.glags[2], pa.np2)
    return pa
end

function mod!(pa::Pconv{T,N,N,N}, ::S; g=pa.g, d=pa.d, s=pa.s) where {N,T<:Real}
    initialize_all!(pa, d, g, s)
    mul!(pa.gfreq, pa.gfftp, pa.gpad)
    mul!(pa.dfreq, pa.dfftp, pa.dpad)
    conj!(pa.gfreq)
    for i in eachindex(pa.dfreq)
        @inbounds pa.sfreq[i] = pa.dfreq[i] * pa.gfreq[i]
    end
    mul!(pa.spad, pa.sifftp, pa.sfreq)
    truncate!(s, pa.spad, pa.slags[1], pa.slags[2], pa.np2)
    return pa
end



# dsum=true
function mod!(pa::Pconv{T,1,2,2}, ::D; g=pa.g, d=pa.d, s=pa.s) where {T<:Real}
    initialize_all!(pa, d, g, s)
    mul!(pa.sfreq, pa.sfftp, pa.spad)
    mul!(pa.gfreq, pa.gfftp, pa.gpad)
    for i in CartesianIndices(size(pa.gfreq))
        @inbounds pa.dfreq[i[1]] += pa.gfreq[i] * pa.sfreq[i]
    end
    mul!(pa.dpad, pa.difftp, pa.dfreq)
    truncate!(d, pa.dpad, pa.dlags[1], pa.dlags[2], pa.np2)
    return pa
end

function mod!(pa::Pconv{T,1,2,2}, ::G; g=pa.g, d=pa.d, s=pa.s) where {T<:Real}
    initialize_all!(pa, d, g, s)
    mul!(pa.sfreq, pa.sfftp, pa.spad)
    mul!(pa.dfreq, pa.dfftp, pa.dpad)
    conj!(pa.sfreq)
    for i in CartesianIndices(size(pa.gfreq))
        @inbounds pa.gfreq[i] = pa.dfreq[i[1]] * pa.sfreq[i]
    end
    mul!(pa.gpad, pa.gifftp, pa.gfreq)
    truncate!(g, pa.gpad, pa.glags[1], pa.glags[2], pa.np2)
    return pa
end

function mod!(pa::Pconv{T,1,2,2}, ::S; g=pa.g, d=pa.d, s=pa.s) where {T<:Real}
    initialize_all!(pa, d, g, s)
    mul!(pa.gfreq, pa.gfftp, pa.gpad)
    mul!(pa.dfreq, pa.dfftp, pa.dpad)
    conj!(pa.gfreq)
    for i in CartesianIndices(size(pa.gfreq))
        @inbounds pa.sfreq[i] = pa.dfreq[i[1]] * pa.gfreq[i]
    end
    mul!(pa.spad, pa.sifftp, pa.sfreq)
    truncate!(s, pa.spad, pa.slags[1], pa.slags[2], pa.np2)
    return pa
end

# gsum=true
function mod!(pa::Pconv{T,2,1,2}, ::D; g=pa.g, d=pa.d, s=pa.s) where {T<:Real}
    initialize_all!(pa, d, g, s)
    mul!(pa.sfreq, pa.sfftp, pa.spad)
    mul!(pa.gfreq, pa.gfftp, pa.gpad)
    for i in CartesianIndices(size(pa.sfreq))
        @inbounds pa.dfreq[i] = pa.gfreq[i[1]] * pa.sfreq[i]
    end
    mul!(pa.dpad, pa.difftp, pa.dfreq)
    truncate!(d, pa.dpad, pa.dlags[1], pa.dlags[2], pa.np2)
    return pa
end
function mod!(pa::Pconv{T,2,1,2}, ::G; g=pa.g, d=pa.d, s=pa.s) where {T<:Real}
    initialize_all!(pa, d, g, s)
    mul!(pa.sfreq, pa.sfftp, pa.spad)
    mul!(pa.dfreq, pa.dfftp, pa.dpad)
    conj!(pa.sfreq)
    for i in CartesianIndices(size(pa.sfreq))
        @inbounds pa.gfreq[i[1]] += pa.dfreq[i] * pa.sfreq[i]
    end
    mul!(pa.gpad, pa.gifftp, pa.gfreq)
    truncate!(g, pa.gpad, pa.glags[1], pa.glags[2], pa.np2)
    return pa
end

function mod!(pa::Pconv{T,2,1,2}, ::S; g=pa.g, d=pa.d, s=pa.s) where {T<:Real}
    initialize_all!(pa, d, g, s)
    mul!(pa.gfreq, pa.gfftp, pa.gpad)
    mul!(pa.dfreq, pa.dfftp, pa.dpad)
    conj!(pa.gfreq)
    for i in CartesianIndices(size(pa.sfreq))
        @inbounds pa.sfreq[i] = pa.dfreq[i] * pa.gfreq[i[1]]
    end
    mul!(pa.spad, pa.sifftp, pa.sfreq)
    truncate!(s, pa.spad, pa.slags[1], pa.slags[2], pa.np2)
    return pa
end

# ssum=true
function mod!(pa::Pconv{T,2,2,1}, ::D; g=pa.g, d=pa.d, s=pa.s) where {T<:Real}
    initialize_all!(pa, d, g, s)
    mul!(pa.sfreq, pa.sfftp, pa.spad)
    mul!(pa.gfreq, pa.gfftp, pa.gpad)
    for i in CartesianIndices(size(pa.gfreq))
        @inbounds pa.dfreq[i] = pa.gfreq[i] * pa.sfreq[i[1]]
    end
    mul!(pa.dpad, pa.difftp, pa.dfreq)
    truncate!(d, pa.dpad, pa.dlags[1], pa.dlags[2], pa.np2)
    return pa
end
function mod!(pa::Pconv{T,2,2,1}, ::G; g=pa.g, d=pa.d, s=pa.s) where {T<:Real}
    initialize_all!(pa, d, g, s)
    mul!(pa.sfreq, pa.sfftp, pa.spad)
    mul!(pa.dfreq, pa.dfftp, pa.dpad)
    conj!(pa.sfreq)
    for i in CartesianIndices(size(pa.gfreq))
        @inbounds pa.gfreq[i] = pa.dfreq[i] * pa.sfreq[i[1]]
    end
    mul!(pa.gpad, pa.gifftp, pa.gfreq)
    truncate!(g, pa.gpad, pa.glags[1], pa.glags[2], pa.np2)
    return pa
end
function mod!(pa::Pconv{T,2,2,1}, ::S; g=pa.g, d=pa.d, s=pa.s) where {T<:Real}
    initialize_all!(pa, d, g, s)
    mul!(pa.gfreq, pa.gfftp, pa.gpad)
    mul!(pa.dfreq, pa.dfftp, pa.dpad)
    conj!(pa.gfreq)
    for i in CartesianIndices(size(pa.gfreq))
        @inbounds pa.sfreq[i[1]] += pa.dfreq[i] * pa.gfreq[i]
    end
    mul!(pa.spad, pa.sifftp, pa.sfreq)
    truncate!(s, pa.spad, pa.slags[1], pa.slags[2], pa.np2)
    return pa
end


"""
A linear operator corresponding to convolution with either s or g
"""
function operator(pa, attrib::S)

    fw = (y, x) -> F!(y, x, pa, attrib)
    bk = (y, x) -> Fadj!(y, x, pa, attrib)

    return LinearMap{collect(typeof(pa).parameters)[1]}(fw, bk,
        length(pa.d),  # length of output
        length(pa.g),  # length of input
        ismutating=true)
end

function operator(pa, attrib::G)

    fw = (y, x) -> F!(y, x, pa, attrib)
    bk = (y, x) -> Fadj!(y, x, pa, attrib)

    return LinearMap{collect(typeof(pa).parameters)[1]}(fw, bk,
        length(pa.d),  # length of output
        length(pa.s),  # length of input
        ismutating=true)
end


function F!(y, x, pa, ::S)
    copyto!(pa.g, x)
    mod!(pa, D())
    copyto!(y, pa.d)
    return nothing
end

function Fadj!(y, x, pa, ::S)
    copyto!(pa.d, x)
    mod!(pa, G())
    copyto!(y, pa.g)
    return nothing
end

function F!(y, x, pa, ::G)
    copyto!(pa.s, x)
    mod!(pa, D())
    copyto!(y, pa.d)
    return nothing
end

function Fadj!(y, x, pa, ::G)
    copyto!(pa.d, x)
    mod!(pa, S())
    copyto!(y, pa.s)
    return nothing
end

function F!(y, x, pa, ::D)
    copyto!(pa.g, x)
    mod!(pa, S())
    copyto!(y, pa.s)
    return nothing
end

function Fadj!(y, x, pa, ::D)
    copyto!(pa.s, x)
    mod!(pa, G())
    copyto!(y, pa.g)
    return nothing
end


end # module

# ╔═╡ 16250068-0132-45c6-9de8-971afcb7bf78
md"# Documentation"

# ╔═╡ c6d68e26-e95d-41d8-b860-6bc4d9c1514e
#=╠═╡
TableOfContents()
  ╠═╡ =#

# ╔═╡ 6c5788a6-a4e7-4edb-be5f-4203d9fe7a36
md"""
```math
D = Conv(G, S)
```
"""

# ╔═╡ 682f12c9-3753-4398-9135-ad62ea0995c0
s0=zeros(3); s0[3]=1.0;

# ╔═╡ 583445f3-4704-4465-ba48-fc1e78d370e4
begin
	pa=Conv.Pconv(dsize=(10, 2), gsize=(10,2), ssize=(3,), s=randn(3), g=randn(10,2), d=randn(10,2));
	Conv.mod!(pa, Conv.D())
end

# ╔═╡ d9020093-5f33-4671-aff3-498f39733011
A=Conv.operator(pa, Conv.G())

# ╔═╡ 5d76cc87-5134-4a60-9f7b-55481e0a2442
strue = randn(3)

# ╔═╡ 06ffb8ad-a456-4140-877f-778299d8455c
dtrue = A*strue

# ╔═╡ 4b9df61f-d5ff-41b9-a623-d7220bcdb513
#=╠═╡
w=zero(strue); IterativeSolvers.lsmr!(w, A, vec(dtrue))
  ╠═╡ =#

# ╔═╡ 4fa17014-2cd7-4ed2-8424-67bf42e971ca
#=╠═╡
w, strue
  ╠═╡ =#

# ╔═╡ 9bf0c426-93fa-4634-8e23-cb21d48a3d92
# adj_test(A)

# ╔═╡ 8645858f-eba8-4824-9ab2-1c3ac853a8b7
Array(A)

# ╔═╡ 013510e5-b81f-425b-8a38-1c36f2f7178e
Array(A')

# ╔═╡ 64bb10c8-c021-4167-bbc4-09f236fc8e0a
gtest=randn(10,2)

# ╔═╡ 9f3c1ea0-e590-42ef-bbc8-db4d13f6f5cf
# reshape(A*vec(gtest), 10, 2)

# ╔═╡ b78e3c5d-23ce-4bb7-bef1-ea41cd831526
Conv.operator(pa, Conv.G())

# ╔═╡ c190136a-4993-402c-8ba7-8e03e658b6c9
md"## Appendix"

# ╔═╡ a9a53efa-86d9-4fd4-adc7-672889b58e26
# ╠═╡ skip_as_script = true
#=╠═╡
function adj_test(F)
    x = randn(size(F, 2))
    y = randn(size(F, 1))
    a = LinearAlgebra.dot(y, F * x)
    b = LinearAlgebra.dot(x, adjoint(F) * y)
    c = LinearAlgebra.dot(x, transpose(F) * F * x)
    println("must be positive: ", c)
    println("adjoint test: ", a, "\t", b)
    @test isapprox(a, b, rtol=1e-5)
    @test c > 0.0
    return nothing
end
  ╠═╡ =#

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DSP = "717857b8-e6f2-59f4-9121-6e50c889abd2"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
IterativeSolvers = "42fd0dbc-a981-5370-80f2-aaf504508153"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
LinearMaps = "7a12625a-238d-50fd-b39a-03d52299707e"
PlutoTest = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
DSP = "~0.7.8"
FFTW = "~1.6.0"
IterativeSolvers = "~0.9.2"
LinearMaps = "~3.10.0"
PlutoTest = "~0.2.2"
PlutoUI = "~0.7.51"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.0"
manifest_format = "2.0"
project_hash = "42fd3d3f59a7c020ea869ea4d7be1f94f869fa57"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "16b6dbc4cf7caee4e1e75c49485ec67b667098a0"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.3.1"
weakdeps = ["ChainRulesCore"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

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

[[deps.DSP]]
deps = ["Compat", "FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "da8b06f89fce9996443010ef92572b193f8dca1f"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.7.8"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "f9818144ce7c8c41edf5c4c179c684d92aa4d9fe"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.6.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

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
git-tree-sha1 = "d75853a0bdbfb1ac815478bacd89cd27b550ace6"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.3"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0cb9352ef2e01574eeebdb102948a58740dcaf83"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2023.1.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "1169632f425f79429f245113b775a0e3d121457c"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.2"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

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

[[deps.LinearMaps]]
deps = ["ChainRulesCore", "LinearAlgebra", "SparseArrays", "Statistics"]
git-tree-sha1 = "4af48c3585177561e9f0d24eb9619ad3abf77cc7"
uuid = "7a12625a-238d-50fd-b39a-03d52299707e"
version = "3.10.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "0a1b7c2863e44523180fdb3146534e265a91870b"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.23"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "2ce8695e1e699b68702c03402672a69f54b8aca9"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2022.2.0+0"

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

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "a5aef8d4a6e8d81f171b2bd4be5265b01384c74c"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.10"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.0"

[[deps.PlutoTest]]
deps = ["HypertextLiteral", "InteractiveUtils", "Markdown", "Test"]
git-tree-sha1 = "17aa9b81106e661cffa1c4c36c17ee1c50a86eda"
uuid = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
version = "0.2.2"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "b478a748be27bd2f2c73a7690da219d0844db305"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.51"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "RecipesBase"]
git-tree-sha1 = "2857c96cdd343a13e8d78f77823bacc8266278ed"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "3.2.11"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsMakieCoreExt = "MakieCore"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    MakieCore = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"

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

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "ef28127915f4229c971eb43f3fc075dd3fe91880"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.2.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

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

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

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
# ╠═a9dd629b-e19e-464f-a9e6-a9675f9ee4ea
# ╠═f01b439a-5641-488c-b175-a37da23dd94c
# ╟─16250068-0132-45c6-9de8-971afcb7bf78
# ╠═c6d68e26-e95d-41d8-b860-6bc4d9c1514e
# ╟─6c5788a6-a4e7-4edb-be5f-4203d9fe7a36
# ╠═682f12c9-3753-4398-9135-ad62ea0995c0
# ╠═583445f3-4704-4465-ba48-fc1e78d370e4
# ╠═d9020093-5f33-4671-aff3-498f39733011
# ╠═5d76cc87-5134-4a60-9f7b-55481e0a2442
# ╠═06ffb8ad-a456-4140-877f-778299d8455c
# ╠═4fa17014-2cd7-4ed2-8424-67bf42e971ca
# ╠═4b9df61f-d5ff-41b9-a623-d7220bcdb513
# ╠═9bf0c426-93fa-4634-8e23-cb21d48a3d92
# ╠═8645858f-eba8-4824-9ab2-1c3ac853a8b7
# ╠═013510e5-b81f-425b-8a38-1c36f2f7178e
# ╠═64bb10c8-c021-4167-bbc4-09f236fc8e0a
# ╠═9f3c1ea0-e590-42ef-bbc8-db4d13f6f5cf
# ╠═b78e3c5d-23ce-4bb7-bef1-ea41cd831526
# ╟─c190136a-4993-402c-8ba7-8e03e658b6c9
# ╠═a9a53efa-86d9-4fd4-adc7-672889b58e26
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
