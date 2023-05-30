### A Pluto.jl notebook ###
# v0.19.21

using Markdown
using InteractiveUtils

# ╔═╡ a6ce608f-4e36-4626-afc3-64dc8d3a3acc
using AxisArrays, MLUtils, Parameters, Statistics, LinearAlgebra, Random, NamedArrays, OrderedCollections, PlutoUI, BenchmarkTools

# ╔═╡ 85e1393a-63df-4a76-a93f-fe33afe272cd
TableOfContents()

# ╔═╡ 6ff6935a-e8c6-4dbe-bf24-4389b5857fb6
begin
    for f in [:Srcs, :Recs]
        @eval(@with_kw_noshow mutable struct $f{T}
            n::Int
            grid::K where {K<:StepRangeLen{Float64}} = range(0.0, 1.0, length=2)
            fields::Vector{Symbol} = [:vz]
            d::NamedArrays.NamedArray{
                Array{T,2},
                1,
                Array{Array{T,2},1},
                Tuple{OrderedCollections.OrderedDict{Symbol,Int64}},
            } = NamedArray([zeros(T, length(grid), n) for i in fields], (fields,))
            function $f{T}(n, grid, fields, d) where {T}
                @assert all(map(x -> size(x) == (length(grid), n), d)) "inconsistent data size"
                new(n, StepRangeLen{Float64}(grid), fields, d)
            end
        end)
        # if only one positional argument is given
        @eval $f(n) = $f{Data.Number}(n=n)
        @eval $f{T}(n) where {T} = $f{T}(n=n)
        @eval $f(n, grid, fields) = $f{Data.Number}(n=n, grid=grid, fields=fields)
        @eval $f{T}(n, grid, fields) where {T} = $f{T}(n=n, grid=grid, fields=fields)
        @eval $f(n, grid) = $f{Data.Number}(n=n, grid=grid)
        @eval $f{T}(n, grid) where {T} = $f{T}(n=n, grid=grid)
    end
end;

# ╔═╡ cecf392d-fd64-4c1d-8088-73818f308453
"""
Simplify getting properties 
"""
function Base.getindex(obj::Union{Srcs,Recs}, sym::Symbol)
    if (sym ∈ [:nr, :ns, :n])
        return obj.n
    else
        return Base.getindex(obj.d, sym) # for fields
    end
end;

# ╔═╡ deb5c48c-1c6d-424b-9fe0-12707b7906af
md"## Structs"

# ╔═╡ 518e57fe-ff05-4f1c-a1c5-789e6d22ae66
md"## Base"

# ╔═╡ 4d63a3ef-5a6f-4bf2-b47b-6ab28f44fda8
"""
```julia
isequal(srcwav[1], srcwav[2])
```
Assert if wavelets due to first two supersources are equal.
"""
function Base.isequal(dat1::Union{Srcs,Recs}, dat2::Union{Srcs,Recs})
    return (dat1.grid == dat2.grid) && (dat1.d == dat2.d) && (dat1.n == dat2.n) && (dat1.fields == dat2.fields)
end;

# ╔═╡ 2e99149d-bbb4-466c-bcf6-d04381d83aec
begin
	"""
   	Flatten
    """
    function Iterators.flatten(data::Union{Srcs,Recs})
        return Iterators.flatten(data.d.array)
    end
	function Iterators.flatten(data::Vector{T}) where {T<:Union{Srcs,Recs}}
        return Iterators.flatten(map(Iterators.flatten, data))
    end
    """
    Return a vec of data object sorted in the order
    time, channels
    """
    function Base.vec(data::Union{Srcs,Recs})
        return mapreduce(vec, vcat, data.d)
    end
    """
    Return a vec of data object sorted in the order
    time, channels
    """
    function Base.vec(data::Vector{T}) where {T<:Union{Srcs,Recs}}
        return mapreduce(vec, vcat, data)
    end
end;

# ╔═╡ b0b83531-aeb5-49b9-8eec-4fb768b17c29
function Base.length(data::Union{Srcs,Recs})
    return length(data.grid) * data[:n] * length(data.fields)
end

# ╔═╡ b761d9d5-5023-4fbe-9f17-d186cabfae1d
# return (number of grid samples, number of channels)
function Base.size(data::Union{Srcs,Recs})
    return (length(data.grid), data.n * length(data.fields))
end

# ╔═╡ 87d8d3a9-46bf-4782-950e-6328e25a436d
function Base.hcat(data::Union{Srcs,Recs})
    return reduce(hcat, data.d)
end

# ╔═╡ 075a4d63-aa3a-4984-9048-8aa0c83184de
function Base.hcat(data::Vector{T}) where {T<:Union{Srcs,Recs}}
    return mapreduce(hcat, hcat, data)
end

# ╔═╡ aa4e23f9-a3a3-4672-b8fd-2e2dfec7c913
function Base.copyto!(d::Union{Srcs,Recs}, v::AbstractArray)
    @assert length(d) == length(v)
    vc = MLUtils.chunk(v, length(d.d)) # default dims is good!
    broadcast(d.d, vc) do d1, vc1
        copyto!(d1, vc1)
    end
end

# ╔═╡ 350ce42b-062d-4b60-8f23-265b04a97ca0
function Base.copyto!(v::AbstractArray, d::Union{Srcs,Recs})
    @assert length(d) == length(v)
    vc = MLUtils.chunk(v, length(d.d))
    broadcast(d.d, vc) do d1, vc1
        copyto!(vc1, d1)
    end
end

# ╔═╡ cc870915-be1f-470d-91f4-06f33ae916bc
begin
	function Base.copyto!(d1::T, d2::T) where {T<:Union{Srcs,Recs}}
		broadcast(Base.copyto!, d1.d, d2.d)
	end
	function Base.copyto!(d1::Vector{T}, d2::Vector{T}) where {T<:Union{Srcs,Recs}}
		broadcast(Base.copyto!, d1, d2)
	end
end

# ╔═╡ 6b3a2ede-54c3-4d70-a780-461b903a06cc
begin
    """
    ```julia
    v=vec(srcwav) 
    randn!(v)
    copyto!(srcwav, v)
    ```
    Copy `v` to `srcwav`. The sizes of the two objects much match.
    No memory allocations.
    """
    function Base.copyto!(d::Vector{T}, v::AbstractVector) where {T<:Union{Srcs,Recs}}
        @assert length(d) == length(v)
        vc = MLUtils.chunk(v, length(d))
        broadcast(d, v) do d1, v1
            copyto!(d1, v1)
        end
    end

    """
    ```julia
    copyto!(v, srcwav)
    ```
    Same as `vec(srcwav)`, but in-place operation.
    """
    function Base.copyto!(v::AbstractVector{Float64}, d::AbstractVector{T}) where {T<:Union{Srcs,Recs}}
        @assert length(d) == length(v)
        vc = MLUtils(v, length(d))
        broadcast(d, v) do d1, v1
            copyto!(v1, d1)
        end
    end
end;

# ╔═╡ 4c7edc3d-9456-4c07-b168-b1b5ef90173e
begin
    """
    ```julia
    iszero(srcwav)
    ```
    Returns bool if `srcwav` has all zeros.
    # """
    function Base.iszero(data::Vector{T}) where {T<:Union{Srcs,Recs}}
        return all([iszero(d) for d in data])
    end

    function Base.iszero(data::Union{Srcs,Recs})
        return all(broadcast(iszero, data.d))
    end
end

# ╔═╡ d51e190b-73e5-4d43-8ebf-697abc4e7cc9
function Base.fill!(data::Union{Srcs,Recs}, k::Float64)
    foreach(data.d) do dd
        fill!(dd, k)
    end
end

# ╔═╡ 19b03f1c-9666-498e-804c-7ded3f632151
"""
    ```julia
    fill!(srcwav, b)
    ```
    In-place method to fill `srcwav` with scalar `b`.
    """
function Base.fill!(data::Vector{T}, k::Float64) where {T<:Union{Srcs,Recs}}
    foreach(data) do d
        fill!(d, k)
    end
end;

# ╔═╡ 45c5cfe0-992a-4b3f-aea3-adb3e134d351
begin
    function Base.reverse!(data::Union{Srcs,Recs})
        # time reversal
        foreach(data.d) do d # each field
            reverse!(d, dims=1)
        end
    end
    """
 ```julia
 reverse!(srcwav)
 ```
 Perform in-place time-reversal operation for each wavelet in `srcwav`.
 """
    function Base.reverse!(data::Vector{T}) where {T<:Union{Srcs,Recs}}
        # time reversal
        foreach(data) do d # each supersource
            reverse!(d)
        end
    end
end;

# ╔═╡ e3e628cf-aa8d-46a9-8a97-b29fc6e1855b
md"## Print"

# ╔═╡ 35e0332b-c828-4344-9824-36c7bfff41bb
begin

    """
    Print information about `::Union{Srcs, Recs}`.
    """
    function Base.show(io::Base.IO, ::MIME"text/plain", src::Union{Srcs,Recs})
        s = isa(src, Srcs) ? "  ├─── # source(s): " : "  ├─ # receiver(s): "
        print(
            io,
            "Data of a supersource\n",
            "  ├───────  fields: ",
            AxisArrays.names(src.d)[1],
            '\n',
            s,
            src[:n],
            '\n',
            "  ├── time samples: ",
            length(src.grid),
            '\n',
        )
    end
    """
    Print information about `Vector{::Union{Srcs, Recs}}`.
    """
    function Base.show(io::Base.IO, ::MIME"text/plain", src::Vector{T}) where {T<:Union{Srcs,Recs}}
        s = (T == Srcs) ? "  ├─── # source(s): " : "  ├─ # receiver(s): "
        print(
            io,
            "Data with $(length(src)) supersource(s)\n",
            "  ├───────  fields: ",
            getindex.(AxisArrays.names.(getfield.(src, :d)), 1),
            '\n',
            s,
            getindex.(src, :n),
            '\n',
            "  ├── time samples: ",
            length.(getfield.(src, :grid)),
            '\n',
        )
    end

end;

# ╔═╡ 5ea2fffb-5f85-4d7d-9ace-854fd13f95a8
md"## Statistics"

# ╔═╡ e38e2c29-49e2-4a75-b85f-d8dd48845648
begin

    # """
    # Returns the variance of data
    # """
    # function Statistics.var(data1::V::Union{Srcs, Recs})
    #     σ = 0.0
    #     μ = Statistics.mean(data1)
    #     n = 0
    #     for iss in 1:length(data1)
    #         for ifield in names(data1[iss].d)[1]
    #             for ir = 1:data1[iss].n, it = 1:length(data1[iss].grid)
    #                 n += 1
    #                 σ += (data1[iss].d[ifield][it, ir] - μ)^2
    #             end
    #         end
    #     end
    #     return σ * inv(n)
    # end

    # function Statistics.mean(data1::V::Union{Srcs, Recs})
    #     n = 0
    #     μ = 0.0
    #     for iss in 1:length(data1)
    #         for ifield in names(data1[iss].d)[1]
    #             for ir = 1:data1[iss].n, it = 1:length(data1[iss].grid)
    #                 n += 1
    #                 μ += data1[iss].d[ifield][it, ir]
    #             end
    #         end
    #     end
    #     return μ * inv(n)
    # end


    # function addnoise!(dataN::V::Union{Srcs, Recs}, data::V::Union{Srcs, Recs}, SNR)

    #     σx = Statistics.var(data)

    #     σxN = sqrt(σx^2 * inv(10^(SNR / 10.0)))

    #     # factor to be multiplied to each scalar
    #     α = sqrt(σxN)
    #     for iss in 1:length(data)
    #         for ifield in names(data[iss].d)[1]
    #             for ir = 1:data[iss].n, it = 1:length(data[iss].grid)
    #                 dataN[iss].d[ifield][it, ir] = dataN[iss].d[ifield][it, ir] + α * Random.randn()
    #             end
    #         end
    #     end
    # end



end

# ╔═╡ 3e5e30a4-7d88-4a39-96c0-a15564ecef22
md"## Func Grad"

# ╔═╡ c97cfe9e-c8ed-4c30-a535-746b54103208
begin
    function lossvalue(loss, dataobs::Union{Srcs,Recs}, data::Union{Srcs,Recs})
        return mapreduce(+, dataobs.d, data.d) do d1, d2 # loop over fields
            mapreduce(+, d1, d2) do dd1, dd2 # loop over time, receiver
                loss(dd1, dd2)
            end
        end
    end
    function lossvalue(loss, dataobs::Vector{T}, data::Vector{T}) where {T<:Union{Srcs,Recs}}
        return mapreduce(+, dataobs, data) do d1, d2 # loop over supersources
            lossvalue(loss, d1, d2)
        end
    end


    function gradient!(buffer, loss, dataobs::Union{Srcs,Recs}, data::Union{Srcs,Recs})
        map(buffer.d, dataobs.d, data.d) do g, d1, d2 # loop over fields
            @. g = deriv(loss, d1, d2)
        end
        return buffer
    end
    function gradient!(buffer, loss, dataobs::Vector{T}, data::Vector{T}) where {T<:Union{Srcs,Recs}}
        map(buffer, dataobs, data) do g, d1, d2 # loop over supersources
            gradient!(g, loss, d1, d2)
        end
    end


end

# ╔═╡ c3455940-6f9b-47e5-ab71-ba41c6b5373d
md"## Methods"

# ╔═╡ 4618ade2-d189-44d1-87ae-856439027baa
begin
    function taper!(data::Vector{T}, perc=0.0; bperc=perc, eperc=perc) where {T<:Union{Srcs,Recs}}
        broadcast(data) do d
            taper!(d, perc, bperc=bperc, eperc=eperc)
        end
        return data
    end

    function taper!(data::Union{Srcs,Recs}, perc=0.0; bperc=perc, eperc=perc)
        broadcast(data.d) do dd
            Utils.taper!(dd, bperc=bperc, eperc=eperc)
        end
        return data
    end

end

# ╔═╡ 8f3561fe-158d-430a-8995-810f21254331
begin


    """
    ```
    rmul!(srcwav[1], b)
    ```
    Scale `srcwav` for first supersource by a scalar `b` overwriting in-place.
    """
    function LinearAlgebra.rmul!(dat::Union{Srcs,Recs}, x::Number)
        for dd in dat.d
            rmul!(dd, x)
        end
    end

    """
    ```
    rmul!(srcwav, b)
    ```
    Scale `srcwav` for all supersources by a scalar `b` overwriting in-place.
    """
    function LinearAlgebra.rmul!(dat::Vector{T}, x::Number) where {T<:Union{Srcs,Recs}}
        for d in dat
            rmul!(d, x)
        end
    end

    function Random.randn!(data::Union{Srcs,Recs})
        map(Random.randn!, data.d)
    end


    """
    Fill `srcwav` with Random values.
    """
    function Random.randn!(data::Vector{T}) where {T<:Union{Srcs,Recs}}
        map(Random.randn!, data)
    end



    function LinearAlgebra.dot(data1::Union{Srcs,Recs}, data2::Union{Srcs,Recs})
        return mapreduce(LinearAlgebra.dot, +, data1.d, data2.d)
    end


    """
    ```julia
    dot(srcwav1, srcwav2)
    ```
    Returns dot product. The sizes must match.
    """
    function LinearAlgebra.dot(data1::Vector{T}, data2::Vector{T}) where {T<:Union{Srcs,Recs}}
        return mapreduce(LinearAlgebra.dot, +, data1, data2)
    end
end;

# ╔═╡ 4beba22e-4cab-43e7-89a3-698ef3f8e5ea
begin
	    """
	    Return if two `::Union{Srcs, Recs}`'s have same dimensions and bounds.
	    """
	    function issimilar(dat1::Union{Srcs,Recs}, dat2::Union{Srcs,Recs})
	        return isequal(dat1.grid, dat2.grid) &&
	               isequal(AxisArrays.names(dat1.d), AxisArrays.names(dat2.d)) &&
	               (size(dat1.d) == size(dat2.d))
	    end
	
	    """
	    ```julia
	    issimilar(srcwav1, srcwav2)
	    ```
	    Assert if `srcwav1` and `srcwav2` have same dimensions and fields.
	    """
	    function issimilar(dat1::Vector{T}, dat2::Vector{T}) where {T<:Union{Srcs,Recs}}
	        @assert length(dat1) == length(dat2)
	        return all([issimilar(dat1[i], dat2[i]) for i = 1:length(dat1)])
	    end
	
end;

# ╔═╡ 494a485c-e2bb-42d6-a6f0-ef4281603f01
begin
    """
    ```julia
    srcwav_new=interp(srcwav, grid_new)
    ```
    Interpolates `srcwav` onto a new grid.
    """
    function interp(data::Vector{T}, grid::StepRangeLen, Battrib=:B1) where {T<:Union{Srcs,Recs}}
        nss = length(data)
        # dataout = [::Union{Srcs, Recs}(grid, data[iss].sr, AxisArrays.names(data[iss].d, 1)) for iss = 1:nss]
        interp_spray!(data, dataout, :interp, Battrib)
        return dataout
    end

    function interp(data::Union{Srcs,Recs}, grid::StepRangeLen, Battrib=:B1)
        # dataout = ::Union{Srcs, Recs}(grid, data.sr, AxisArrays.names(data.d, 1))
        interp_spray!(data, dataout, :interp, Battrib)
        return dataout
    end

    """
    ```julia
    interp_spray!(srcwav_new, srcwav, attrib, Battrib)
    ```
    Method to interpolate `srcwav` onto to a new grid.
    Can reduce allocations  ========
    """
    function interp_spray!(
        data::Union{Srcs,Recs},
        dataout::Union{Srcs,Recs},
        attrib=:interp,
        Battrib=:B1;
        pa=nothing
    )
        @assert length(data.d) == length(dataout.d)
        xin = data.grid
        xout = dataout.grid
        if (pa === nothing)
            pa = Interpolation.Kernel([xin], [xout], :B1)
        end
        for i = 1:length(data.d)
            dd = data.d[i]
            ddo = dataout.d[i]
            @assert size(dd, 2) == size(ddo, 2)
            for ir = 1:size(dd, 2)
                din = view(dd, :, ir)
                dout = view(ddo, :, ir)
                Interpolation.interp_spray!(din, dout, pa, attrib)
            end
        end
        return dataout
    end

    function interp_spray!(
        data::Vector{T},
        dataout::Vector{T},
        attrib=:interp,
        Battrib=:B1;
        pa=nothing
    ) where {T<:Union{Srcs,Recs}}
        @assert length(data) == length(dataout)
        for i = 1:length(data)
            d = data[i]
            dout = dataout[i]
            interp_spray!(d, dout, attrib, Battrib, pa=pa)
        end
    end

end

# ╔═╡ 10d5c2c1-1b74-4796-a30b-ac4dcb850d17
begin
    """
    ```julia
    update!(srcwav, [:p, :vx], w)
    ```
    Populate `srcwav` by a wavelet vector `w` for `:p` and `:vx` fields.
    """
    function update!(dat::Vector{T}, fields::Vector{Symbol}, w::AbstractArray) where {T<:Union{Srcs,Recs}}
        foreach(dat) do d
            update!(d, fields, w)
        end
    end


    """
    ```julia
    update!(srcwav[1], [:p, :vx], w)
    ```
    Populate first supersource of `srcwav` by a wavelet vector `w` for `:p` and `:vx` fields.
    """
    function update!(d::Union{Srcs,Recs}, fields::Vector{Symbol}, w::AbstractVector)
        @assert length(w) == length(d.grid)
        broadcast(fields) do f
            @inbounds for i = 1:d.n
                for it = 1:length(d.grid)
                    d.d[f][it, i] = w[it]
                end
            end
        end
        return d
    end

    """
    ```julia
    update!(srcwav[1], [:p, :vx], w)
    ```
    Populate first supersource of `srcwav` by a wavelets `w` for `:p` and `:vx` fields.
    """
    function update!(d::Union{Srcs,Recs}, fields::Vector{Symbol}, w::AbstractMatrix)
        @assert size(w, 1) == length(d.grid)
        @assert size(w, 2) == d.n
        @inbounds for f in fields
            @assert f ∈ AxisArrays.names(d.d)[1]
            @inbounds for i = 1:d.n
                for it = 1:length(d.grid)
                    d.d[f][it, i] = w[it, i]
                end
            end
        end
        return d
    end
end;

# ╔═╡ 0b5e45ae-0368-4d06-8eb0-4f2696c82868
md"# Documentation"

# ╔═╡ 480a4b78-d6b0-4b39-ada9-7347bf8c7817
grid = range(1.0, stop=100.0, step=1.0)

# ╔═╡ d36e5f9c-10cd-4581-989d-906502fcba2a
data = Srcs{Float32}(10, grid)

# ╔═╡ 00bc946c-a071-43be-9cc5-f7e153371d0a
length(data)

# ╔═╡ 7fbb9b5a-e921-4a30-b45e-7ba84a24829d
size(data)

# ╔═╡ 20bc27ea-1df4-4848-b854-1a428e65b28a
stack

# ╔═╡ 4bece999-623e-4695-a85b-8937f44e48ee
vcat(view(randn(10,10), :))

# ╔═╡ 81af3b1a-e61e-44a0-82f7-4cca66fc480a
hcat(data)

# ╔═╡ e6070bc0-27d0-453f-a1bc-662321321e45
dv = vec(data)

# ╔═╡ c661d5a8-8b4e-4d48-ac4b-10c9c5cddb7b
Iterators.flatten([data, data]) |> collect

# ╔═╡ 4bb3be76-7c4d-416f-bc2d-968fbcf3968e
@btime copyto!(dv, data)

# ╔═╡ bd32bfeb-d9f9-4eb5-b578-cccef76963d7
@btime copyto!(data, dv)

# ╔═╡ 20542ff1-c23e-41a4-b457-9adb3575bccf
md"## Appendix"

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
AxisArrays = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
MLUtils = "f1d291b0-491e-4a28-83b9-f70985020b54"
NamedArrays = "86f7a689-2022-50b4-a561-43c23ac3c673"
OrderedCollections = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
Parameters = "d96e819e-fc66-5662-9728-84c9c7592b0a"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
AxisArrays = "~0.4.6"
BenchmarkTools = "~1.3.2"
MLUtils = "~0.4.2"
NamedArrays = "~0.9.8"
OrderedCollections = "~1.6.0"
Parameters = "~0.12.3"
PlutoUI = "~0.7.51"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.0"
manifest_format = "2.0"
project_hash = "c39f102a07a88b86ef500ccc3407cdaa59a3ca33"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.Accessors]]
deps = ["Compat", "CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "LinearAlgebra", "MacroTools", "Requires", "Test"]
git-tree-sha1 = "2b301c2388067d655fe5e4ca6d4aa53b61f895b4"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.31"

    [deps.Accessors.extensions]
    AccessorsAxisKeysExt = "AxisKeys"
    AccessorsIntervalSetsExt = "IntervalSets"
    AccessorsStaticArraysExt = "StaticArrays"
    AccessorsStructArraysExt = "StructArrays"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "76289dc51920fdc6e0013c872ba9551d54961c24"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.6.2"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.ArgCheck]]
git-tree-sha1 = "a3a402a35a2f7e0b87828ccabbd5ebfbebe356b4"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.3.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Atomix]]
deps = ["UnsafeAtomics"]
git-tree-sha1 = "c06a868224ecba914baa6942988e2f2aade419be"
uuid = "a9b6321e-bd34-4604-b9c9-b65b8de01458"
version = "0.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "1dd4d9f5beebac0c03446918741b1a03dc5e5788"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.6"

[[deps.BangBang]]
deps = ["Compat", "ConstructionBase", "InitialValues", "LinearAlgebra", "Requires", "Setfield", "Tables"]
git-tree-sha1 = "54b00d1b93791f8e19e31584bd30f2cb6004614b"
uuid = "198e06fe-97b7-11e9-32a5-e1d131e6ad66"
version = "0.3.38"

    [deps.BangBang.extensions]
    BangBangChainRulesCoreExt = "ChainRulesCore"
    BangBangDataFramesExt = "DataFrames"
    BangBangStaticArraysExt = "StaticArrays"
    BangBangStructArraysExt = "StructArrays"
    BangBangTypedTablesExt = "TypedTables"

    [deps.BangBang.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    TypedTables = "9d95f2ec-7b3d-5a63-8d20-e2491e220bb9"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Baselet]]
git-tree-sha1 = "aebf55e6d7795e02ca500a689d326ac979aaf89e"
uuid = "9718e550-a3fa-408a-8086-8db961cd8217"
version = "0.1.1"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "d9a9701b899b30332bbcb3e1679c41cce81fb0e8"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.3.2"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

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

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

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

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "738fec4d684a9a6ee9598a8bfee305b26831f28c"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.2"
weakdeps = ["IntervalSets", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.ContextVariablesX]]
deps = ["Compat", "Logging", "UUIDs"]
git-tree-sha1 = "25cc3803f1030ab855e383129dcd3dc294e322cc"
uuid = "6add18c4-b38d-439d-96f6-d6bc489c04c5"
version = "0.1.3"

[[deps.DataAPI]]
git-tree-sha1 = "8da84edb865b0b5b0100c0666a9bc9a0b71c553c"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.15.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DefineSingletons]]
git-tree-sha1 = "0fba8b706d0178b4dc7fd44a96a92382c9065c2c"
uuid = "244e2a9f-e319-4986-a169-4d1fe445cd52"
version = "0.1.2"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FLoops]]
deps = ["BangBang", "Compat", "FLoopsBase", "InitialValues", "JuliaVariables", "MLStyle", "Serialization", "Setfield", "Transducers"]
git-tree-sha1 = "ffb97765602e3cbe59a0589d237bf07f245a8576"
uuid = "cc61a311-1640-44b5-9fba-1b764f453329"
version = "0.2.1"

[[deps.FLoopsBase]]
deps = ["ContextVariablesX"]
git-tree-sha1 = "656f7a6859be8673bf1f35da5670246b923964f7"
uuid = "b9860ae5-e623-471e-878b-f6a53c775ea6"
version = "0.1.1"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.FoldsThreads]]
deps = ["Accessors", "FunctionWrappers", "InitialValues", "SplittablesBase", "Transducers"]
git-tree-sha1 = "eb8e1989b9028f7e0985b4268dabe94682249025"
uuid = "9c68100b-dfe1-47cf-94c8-95104e173443"
version = "0.1.1"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "2d6ca471a6c7b536127afccfa7564b5b39227fe0"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.5"

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

[[deps.InitialValues]]
git-tree-sha1 = "4da0f88e9a39111c2fa3add390ab15f3a44f3ca3"
uuid = "22cec73e-a1b8-11e9-2c92-598750a2cf9c"
version = "0.3.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IntervalSets]]
deps = ["Dates", "Random", "Statistics"]
git-tree-sha1 = "16c0cc91853084cb5f58a78bd209513900206ce6"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.4"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "6667aadd1cdee2c6cd068128b3d226ebc4fb0c67"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.9"

[[deps.InvertedIndices]]
git-tree-sha1 = "0dc7b50b8d436461be01300fd8cd45aa0274b038"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

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

[[deps.JuliaVariables]]
deps = ["MLStyle", "NameResolution"]
git-tree-sha1 = "49fb3cb53362ddadb4415e9b73926d6b40709e70"
uuid = "b14d175d-62b4-44ba-8fb7-3064adc8c3ec"
version = "0.2.4"

[[deps.KernelAbstractions]]
deps = ["Adapt", "Atomix", "InteractiveUtils", "LinearAlgebra", "MacroTools", "PrecompileTools", "SparseArrays", "StaticArrays", "UUIDs", "UnsafeAtomics", "UnsafeAtomicsLLVM"]
git-tree-sha1 = "47be64f040a7ece575c2b5f53ca6da7b548d69f4"
uuid = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
version = "0.9.4"

[[deps.LLVM]]
deps = ["CEnum", "LLVMExtra_jll", "Libdl", "Printf", "Unicode"]
git-tree-sha1 = "26a31cdd9f1f4ea74f649a7bf249703c687a953d"
uuid = "929cbde3-209d-540e-8aea-75f648917ca0"
version = "5.1.0"

[[deps.LLVMExtra_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl", "TOML"]
git-tree-sha1 = "09b7505cc0b1cee87e5d4a26eea61d2e1b0dcd35"
uuid = "dad2f222-ce93-54a1-a47d-0025e8a3acab"
version = "0.0.21+0"

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

[[deps.MLStyle]]
git-tree-sha1 = "bc38dff0548128765760c79eb7388a4b37fae2c8"
uuid = "d8e11817-5142-5d16-987a-aa16d5891078"
version = "0.4.17"

[[deps.MLUtils]]
deps = ["ChainRulesCore", "Compat", "DataAPI", "DelimitedFiles", "FLoops", "FoldsThreads", "NNlib", "Random", "ShowCases", "SimpleTraits", "Statistics", "StatsBase", "Tables", "Transducers"]
git-tree-sha1 = "ca31739905ddb08c59758726e22b9e25d0d1521b"
uuid = "f1d291b0-491e-4a28-83b9-f70985020b54"
version = "0.4.2"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.MicroCollections]]
deps = ["BangBang", "InitialValues", "Setfield"]
git-tree-sha1 = "629afd7d10dbc6935ec59b32daeb33bc4460a42e"
uuid = "128add7d-3638-4c79-886c-908ea0c25c34"
version = "0.1.4"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.NNlib]]
deps = ["Adapt", "Atomix", "ChainRulesCore", "GPUArraysCore", "KernelAbstractions", "LinearAlgebra", "Pkg", "Random", "Requires", "Statistics"]
git-tree-sha1 = "99e6dbb50d8a96702dc60954569e9fe7291cc55d"
uuid = "872c559c-99b0-510c-b3b7-b6c96a88d5cd"
version = "0.8.20"

    [deps.NNlib.extensions]
    NNlibAMDGPUExt = "AMDGPU"

    [deps.NNlib.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"

[[deps.NameResolution]]
deps = ["PrettyPrint"]
git-tree-sha1 = "1a0fa0e9613f46c9b8c11eee38ebb4f590013c5e"
uuid = "71a1bf82-56d0-4bbc-8a3c-48b961074391"
version = "0.1.5"

[[deps.NamedArrays]]
deps = ["Combinatorics", "DataStructures", "DelimitedFiles", "InvertedIndices", "LinearAlgebra", "Random", "Requires", "SparseArrays", "Statistics"]
git-tree-sha1 = "b84e17976a40cb2bfe3ae7edb3673a8c630d4f95"
uuid = "86f7a689-2022-50b4-a561-43c23ac3c673"
version = "0.9.8"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

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
git-tree-sha1 = "a5aef8d4a6e8d81f171b2bd4be5265b01384c74c"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.10"

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

[[deps.PrettyPrint]]
git-tree-sha1 = "632eb4abab3449ab30c5e1afaa874f0b98b586e4"
uuid = "8162dcfd-2161-5ef2-ae6c-7681170c5f98"
version = "0.2.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

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

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.ShowCases]]
git-tree-sha1 = "7f534ad62ab2bd48591bdeac81994ea8c445e4a5"
uuid = "605ecd9f-84a6-4c9e-81e2-4798472b76a3"
version = "0.1.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SplittablesBase]]
deps = ["Setfield", "Test"]
git-tree-sha1 = "e08a62abc517eb79667d0a29dc08a3b589516bb5"
uuid = "171d559e-b47b-412a-8079-5efa626c420e"
version = "0.1.15"

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

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "45a7769a04a3cf80da1c1c7c60caf932e6f4c9f7"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.6.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "75ebe04c5bed70b91614d684259b661c9e6274a4"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "1544b926975372da01227b382066ab70e574a3ec"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Transducers]]
deps = ["Adapt", "ArgCheck", "BangBang", "Baselet", "CompositionsBase", "DefineSingletons", "Distributed", "InitialValues", "Logging", "Markdown", "MicroCollections", "Requires", "Setfield", "SplittablesBase", "Tables"]
git-tree-sha1 = "25358a5f2384c490e98abd565ed321ffae2cbb37"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.76"

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

[[deps.UnsafeAtomics]]
git-tree-sha1 = "6331ac3440856ea1988316b46045303bef658278"
uuid = "013be700-e6cd-48c3-b4a1-df204f14c38f"
version = "0.2.1"

[[deps.UnsafeAtomicsLLVM]]
deps = ["LLVM", "UnsafeAtomics"]
git-tree-sha1 = "ea37e6066bf194ab78f4e747f5245261f17a7175"
uuid = "d80eeb9a-aca5-4d75-85e5-170c8b632249"
version = "0.1.2"

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
# ╠═85e1393a-63df-4a76-a93f-fe33afe272cd
# ╟─deb5c48c-1c6d-424b-9fe0-12707b7906af
# ╠═6ff6935a-e8c6-4dbe-bf24-4389b5857fb6
# ╟─518e57fe-ff05-4f1c-a1c5-789e6d22ae66
# ╠═cecf392d-fd64-4c1d-8088-73818f308453
# ╠═4d63a3ef-5a6f-4bf2-b47b-6ab28f44fda8
# ╠═2e99149d-bbb4-466c-bcf6-d04381d83aec
# ╠═b761d9d5-5023-4fbe-9f17-d186cabfae1d
# ╠═b0b83531-aeb5-49b9-8eec-4fb768b17c29
# ╠═87d8d3a9-46bf-4782-950e-6328e25a436d
# ╠═075a4d63-aa3a-4984-9048-8aa0c83184de
# ╠═aa4e23f9-a3a3-4672-b8fd-2e2dfec7c913
# ╠═350ce42b-062d-4b60-8f23-265b04a97ca0
# ╠═cc870915-be1f-470d-91f4-06f33ae916bc
# ╠═6b3a2ede-54c3-4d70-a780-461b903a06cc
# ╠═4c7edc3d-9456-4c07-b168-b1b5ef90173e
# ╠═d51e190b-73e5-4d43-8ebf-697abc4e7cc9
# ╠═19b03f1c-9666-498e-804c-7ded3f632151
# ╠═45c5cfe0-992a-4b3f-aea3-adb3e134d351
# ╟─e3e628cf-aa8d-46a9-8a97-b29fc6e1855b
# ╠═35e0332b-c828-4344-9824-36c7bfff41bb
# ╟─5ea2fffb-5f85-4d7d-9ace-854fd13f95a8
# ╠═e38e2c29-49e2-4a75-b85f-d8dd48845648
# ╟─3e5e30a4-7d88-4a39-96c0-a15564ecef22
# ╠═c97cfe9e-c8ed-4c30-a535-746b54103208
# ╟─c3455940-6f9b-47e5-ab71-ba41c6b5373d
# ╠═4618ade2-d189-44d1-87ae-856439027baa
# ╠═8f3561fe-158d-430a-8995-810f21254331
# ╠═4beba22e-4cab-43e7-89a3-698ef3f8e5ea
# ╠═494a485c-e2bb-42d6-a6f0-ef4281603f01
# ╠═10d5c2c1-1b74-4796-a30b-ac4dcb850d17
# ╟─0b5e45ae-0368-4d06-8eb0-4f2696c82868
# ╠═480a4b78-d6b0-4b39-ada9-7347bf8c7817
# ╠═d36e5f9c-10cd-4581-989d-906502fcba2a
# ╠═00bc946c-a071-43be-9cc5-f7e153371d0a
# ╠═7fbb9b5a-e921-4a30-b45e-7ba84a24829d
# ╠═20bc27ea-1df4-4848-b854-1a428e65b28a
# ╠═4bece999-623e-4695-a85b-8937f44e48ee
# ╠═81af3b1a-e61e-44a0-82f7-4cca66fc480a
# ╠═e6070bc0-27d0-453f-a1bc-662321321e45
# ╠═c661d5a8-8b4e-4d48-ac4b-10c9c5cddb7b
# ╠═4bb3be76-7c4d-416f-bc2d-968fbcf3968e
# ╠═bd32bfeb-d9f9-4eb5-b578-cccef76963d7
# ╟─20542ff1-c23e-41a4-b457-9adb3575bccf
# ╠═a6ce608f-4e36-4626-afc3-64dc8d3a3acc
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
