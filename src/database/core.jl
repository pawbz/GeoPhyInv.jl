


"""
Type for storing NamedD 
"""
mutable struct NamedD{T<:Union{Srcs,Recs}}
    grid::StepRangeLen{Float64}
    d::NamedArrays.NamedArray{
        Array{Float64,2},
        1,
        Array{Array{Float64,2},1},
        Tuple{OrderedCollections.OrderedDict{Symbol,Int64}},
    }
    sr::T
    #	NamedD{T}(grid,d,sr)=any([(size(dd)≠(length(grid),sr.n)) for dd in d]) ? error("NamedD construction") : new(grid,d,sr)
end

"""
Simplify getting properties 
"""
function Base.getindex(obj::NamedD, sym::Symbol)
    if (sym ∈ [:nr, :ns, :n])
        return obj.sr.n
    else
        return Base.getindex(obj.d, sym)
    end
end



"""
Initialize with zeros, given grid and field names.
"""
function Base.zero(
    ::Type{NamedD},
    grid::StepRangeLen,
    s::Union{Srcs,Recs},
    fields::Vector{Symbol},
)
    for f in fields
        @assert eval(f) <: Fields
    end
    return NamedD(
        grid,
        NamedArray([zeros(length(grid), s.n) for i in fields], (fields,)),
        s,
    )
end
NamedD(grid, s::Union{Srcs,Recs}, fields::Vector{Symbol}) = zero(NamedD, grid, s, fields)

"""
Updates the fields.
"""
function update!(dat::NamedD, fields::Vector{Symbol})
    setnames!(dat.d, fields, 1)
    return dat
end

VNamedD = Union{Array{NamedD{Recs},1},Array{NamedD{Srcs},1}}


function Array{NamedD{Recs},1}(
    grid::StepRangeLen,
    ss::SSrcs,
    sr::Vector{Recs},
    fields::Vector{Symbol},
)
    return [NamedD(grid, sr[i], fields) for i = 1:ss.n]
end
function Array{NamedD{Srcs},1}(
    grid::StepRangeLen,
    ss::SSrcs,
    sr::Vector{Srcs},
    fields::Vector{Symbol},
)
    return [NamedD(grid, sr[i], fields) for i = 1:ss.n]
end

function Array{NamedD{Srcs},1}(
    grid::StepRangeLen,
    ss::SSrcs,
    sr::Srcs,
    fields::Vector{Symbol},
)
    return [NamedD(grid, sr, fields) for i = 1:ss.n]
end
function Array{NamedD{Recs},1}(
    grid::StepRangeLen,
    ss::SSrcs,
    sr::Recs,
    fields::Vector{Symbol},
)
    return [NamedD(grid, sr, fields) for i = 1:ss.n]
end


"""
Print information about `NamedD`.
"""
function Base.show(io::Base.IO, ::MIME"text/plain", src::NamedD)
    s = isa(src.sr, Srcs) ? "  ├─── # source(s): " : "  ├─ # receiver(s): "
    print(
        io,
        "Data of a supersource\n",
        "  ├───────  fields: ",
        names(src.d)[1],
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
Print information about `Vector{NamedD}`.
"""
function Base.show(io::Base.IO, ::MIME"text/plain", src::Vector{NamedD{T}}) where {T}
    s = (T == Srcs) ? "  ├─── # source(s): " : "  ├─ # receiver(s): "
    print(
        io,
        "Data with $(length(src)) supersource(s)\n",
        "  ├───────  fields: ",
        getindex.(names.(getfield.(src, :d)), 1),
        '\n',
        s,
        getindex.(src, :n),
        '\n',
        "  ├── time samples: ",
        length.(getfield.(src, :grid)),
        '\n',
    )
end





#function Base.isapprox(src1::NamedD, src2::NamedD)
#	vec=([(src1.nss==src2.nss), (src1.fields==src2.fields), (src1.ns==src2.ns), 
#       		(isequal(src1.grid, src2.grid)),
#		])
#	return all(vec)
#end


"""
```
rmul!(srcwav[1], b)
```
Scale `srcwav` for first supersource by a scalar `b` overwriting in-place.
"""
function LinearAlgebra.rmul!(dat::NamedD, x::Number)
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
function LinearAlgebra.rmul!(dat::VNamedD, x::Number)
    for d in dat
        rmul!(d, x)
    end
end



"""
```julia
update!(srcwav[1], [:p, :vx], w)
```
Populate first supersource of `srcwav` by a wavelet vector `w` for `:p` and `:vx` fields.
"""
function update!(d::NamedD, fields::Vector{Symbol}, w::AbstractVector)
    @assert length(w) == length(d.grid)
    @inbounds for f in fields
        @assert f ∈ names(d.d)[1]
        @inbounds for i = 1:d[:n]
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
function update!(d::NamedD, fields::Vector{Symbol}, w::AbstractMatrix)
    @assert size(w, 1) == length(d.grid)
    @assert size(w, 2) == d[:n]
    @inbounds for f in fields
        @assert f ∈ names(d.d)[1]
        @inbounds for i = 1:d[:n]
            for it = 1:length(d.grid)
                d.d[f][it, i] = w[it, i]
            end
        end
    end
    return d
end





"""
```julia
update!(srcwav, [:p, :vx], w)
```
Populate `srcwav` by a wavelet vector `w` for `:p` and `:vx` fields.
"""
function update!(dat::VNamedD, fields::Vector{Symbol}, w::AbstractArray)
    for d in dat
        update!(d, fields, w)
    end
end


"""
```julia
isequal(srcwav[1], srcwav[2])
```
Assert if wavelets due to first two supersources are equal.
"""
function Base.isequal(dat1::NamedD, dat2::NamedD)
    return (dat1.grid == dat2.grid) && (dat1.d == dat2.d)
end

"""
```julia
isequal(srcwav1, srcwav2)
```
Assert if `srcwav1` and `srcwav2` are equal.
"""
function Base.isequal(dat1::VNamedD, dat2::VNamedD)
    return all([isequal(dat1[i], dat2[i]) for i = 1:length(dat1)])
end

"""
Return if two `NamedD`'s have same dimensions and bounds.
"""
function issimilar(dat1::NamedD, dat2::NamedD)
    return isequal(dat1.grid, dat2.grid) &&
           isequal(names(dat1.d), names(dat2.d)) &&
           (size(dat1.d) == size(dat2.d))
end

"""
```julia
issimilar(srcwav1, srcwav2)
```
Assert if `srcwav1` and `srcwav2` have same dimensions and fields.
"""
function issimilar(dat1::VNamedD, dat2::VNamedD)
    @assert length(dat1) == length(dat2)
    return all([issimilar(dat1[i], dat2[i]) for i = 1:length(dat1)])
end

function Base.length(data::NamedD)
    return length(data.grid) * data[:n] * length(names(data.d)[1])
end

"""
Return a vec of data object sorted in the order
time, channels
"""
function Base.vec(data::NamedD)
    v = Vector{Float64}()
    for dd in data.d
        append!(v, dd)
    end
    return v
end


"""
```julia
v=vec(srcwav)
Reshape `srcwav` as a one-dimensional column vector.
```
"""
function Base.vec(data::VNamedD)
    v = Vector{Float64}()
    for d in data
        append!(v, vec(d))
    end
    return v
end

function Base.copyto!(d::NamedD, v::AbstractVector{Float64}, i0=1)
    #	@assert length(d)==length(v)
    for dd in d.d
        for i = 1:d[:n]
            for it = 1:length(d.grid)
                dd[it, i] = v[i0]
                i0 += 1
            end
        end
    end
    return i0
end


"""
```julia
v=vec(srcwav) 
randn!(v)
copyto!(srcwav, v)
```
Copy `v` to `srcwav`. The sizes of the two objects much match.
No memory allocations.
"""
function Base.copyto!(d::VNamedD, v::AbstractVector{Float64})
    i0 = 1
    for i = 1:length(d)
        dd = d[i]
        i0 = copyto!(dd, v, i0)
    end
    return d
end




function Base.copyto!(v::AbstractVector{Float64}, d::NamedD, i0=1)
    for dd in d.d
        for i = 1:d[:n]
            for it = 1:length(d.grid)
                v[i0] = dd[it, i]
                i0 += 1
            end
        end
    end
    return i0
end

"""
```julia
copyto!(v, srcwav)
```
Same as `vec(srcwav)`, but in-place operation.
"""
function Base.copyto!(v::AbstractVector{Float64}, d::VNamedD)
    i0 = 1
    for i = 1:length(d)
        dd = d[i]
        i0 = copyto!(v, dd, i0)
    end
    return d
end


function Random.randn!(data::NamedD)
    for dd in data.d
        Random.randn!(dd)
    end
    return data
end


"""
Fill `srcwav` with Random values.
"""
function Random.randn!(data::VNamedD)
    for d in data
        Random.randn!(d)
    end
end

function Base.copyto!(dataout::NamedD, data::NamedD)
    for i in names(dataout.d)[1]
        d1 = dataout.d[i]
        d2 = data.d[i]
        copyto!(d1, d2)
    end
end

"""
```julia
copyto!(srcwav1, srcwav2)
```
Copy `srcwav2` to `srcwav1`, which has same size.
"""
function Base.copyto!(dataout::VNamedD, data::VNamedD)
    for i = 1:length(dataout)
        d1 = dataout[i]
        d2 = data[i]
        copyto!(d1, d2)
    end
end


function Base.iszero(data::NamedD)
    return maximum(broadcast(maximum, abs, data.d)) == 0.0 ? true : false
end

"""
```julia
iszero(srcwav)
```
Returns bool if `srcwav` has all zeros.
"""
function Base.iszero(data::VNamedD)
    return all([iszero(d) for d in data])
end

function LinearAlgebra.dot(data1::NamedD, data2::NamedD)
    dotd = 0.0
    for iff in names(data1.d)[1]
        dd1 = data1.d[iff]
        dd2 = data2.d[iff]
        dotd += LinearAlgebra.dot(dd1, dd2)
    end
    return dotd
end


"""
```julia
dot(srcwav1, srcwav2)
```
Returns dot product. The sizes must match.
"""
function LinearAlgebra.dot(data1::VNamedD, data2::VNamedD)
    dotd = 0.0
    for i = 1:length(data1)
        dd1 = data1[i]
        dd2 = data2[i]
        dotd += LinearAlgebra.dot(dd1, dd2)
    end
    return dotd
end


function Base.fill!(data::NamedD, k::Float64)
    for dd in data.d
        fill!(dd, k)
    end
end

"""
```julia
fill!(srcwav, b)
```
In-place method to fill `srcwav` with scalar `b`.
"""
function Base.fill!(data::VNamedD, k::Float64)
    for d in data
        fill!(d, k)
    end
end

function Base.reverse!(data::NamedD)
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
function Base.reverse!(data::VNamedD)
    # time reversal
    foreach(data) do d # each supersource
        reverse!(d)
    end
end



include("misfits.jl")
include("statistics.jl")





"""
```julia
srcwav_new=interp(srcwav, grid_new)
```
Interpolates `srcwav` onto a new grid.
"""
function interp(data::VNamedD, grid::StepRangeLen, Battrib=:B1)
    nss = length(data)
    dataout = [NamedD(grid, data[iss].sr, names(data[iss].d, 1)) for iss = 1:nss]
    interp_spray!(data, dataout, :interp, Battrib)
    return dataout
end

function interp(data::NamedD, grid::StepRangeLen, Battrib=:B1)
    dataout = NamedD(grid, data.sr, names(data.d, 1))
    interp_spray!(data, dataout, :interp, Battrib)
    return dataout
end


function taper!(data::VNamedD, perc=0.0; bperc=perc, eperc=perc)
    for d in data
        taper!(d, perc, bperc=bperc, eperc=eperc)
    end
    return data
end

function taper!(data::NamedD, perc=0.0; bperc=perc, eperc=perc)
    for dd in data
        Utils.taper!(dd, bperc=bperc, eperc=eperc)
    end
    return data
end




"""
```julia
interp_spray!(srcwav_new, srcwav, attrib, Battrib)
```
Method to interpolate `srcwav` onto to a new grid.
Can reduce allocations =========


"""
function interp_spray!(
    data::NamedD,
    dataout::NamedD,
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
    data::VNamedD,
    dataout::VNamedD,
    attrib=:interp,
    Battrib=:B1;
    pa=nothing
)
    @assert length(data) == length(dataout)
    for i = 1:length(data)
        d = data[i]
        dout = dataout[i]
        interp_spray!(d, dout, attrib, Battrib, pa=pa)
    end
end



