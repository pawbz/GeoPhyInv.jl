


"""
Type for storing NamedD 
"""
mutable struct NamedD{T<:Union{Srcs,Recs}}
	grid::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
	d::NamedArrays.NamedArray{Array{Float64,2},1,Array{Array{Float64,2},1},Tuple{OrderedCollections.OrderedDict{Symbol,Int64}}}
	sr::T
#	NamedD{T}(grid,d,sr)=any([(size(dd)≠(length(grid),sr.n)) for dd in d]) ? error("NamedD construction") : new(grid,d,sr)
end

"""
Simplify getting properties 
"""
function Base.getindex(obj::NamedD, sym::Symbol)
	if(sym ∈ [:nr, :ns, :n])
		return obj.sr.n
	else
		return Base.getindex(obj.d, sym)
	end
end



"""
Initialize with zeros, given grid and names.
"""
function Base.zero(::Type{NamedD}, grid::StepRangeLen, s::Union{Srcs,Recs}, fields::Vector{Symbol})
	nt=length(grid); nf=length(fields)
	return NamedD(grid,NamedArray([zeros(nt,s.n) for i in fields], (fields,)),s)
end
NamedD(grid, s::Union{Srcs,Recs}, fields::Vector{Symbol})=zero(NamedD, grid, s, fields)

VNamedD=Union{Array{NamedD{Recs},1},Array{NamedD{Srcs},1}}


function Array{NamedD{Recs},1}(grid::StepRangeLen, ss::SSrcs, sr::Vector{Recs}, fields::Vector{Symbol}) 
	return [NamedD(grid,sr[i],fields) for i in 1:ss.n]
end
function Array{NamedD{Srcs},1}(grid::StepRangeLen, ss::SSrcs, sr::Vector{Srcs}, fields::Vector{Symbol})
	return [NamedD(grid,sr[i],fields) for i in 1:ss.n]
end

function Array{NamedD{Srcs},1}(grid::StepRangeLen, ss::SSrcs, sr::Srcs, fields::Vector{Symbol})
	return [NamedD(grid,sr,fields) for i in 1:ss.n]
end
function Array{NamedD{Recs},1}(grid::StepRangeLen, ss::SSrcs, sr::Recs, fields::Vector{Symbol})
	return [NamedD(grid,sr,fields) for i in 1:ss.n]
end


"""
Print information about `NamedD`.
To be implemented
"""
function Base.show(io::Base.IO, src::NamedD)
	println(io,"\tNamedD")
#	println(io,"\t> number of supersources:\t",src.nss)
#	println(io,"\t> sources per supersource:\t","min\t",minimum(src.ns[:]), "\tmax\t", maximum(src.ns[:]))
end



#function Base.isapprox(src1::NamedD, src2::NamedD)
#	vec=([(src1.nss==src2.nss), (src1.fields==src2.fields), (src1.ns==src2.ns), 
#       		(isequal(src1.grid, src2.grid)),
#		])
#	return all(vec)
#end




function LinearAlgebra.rmul!(dat::NamedD, x::Number)
	for dd in dat.d
		rmul!(dd, x)
	end
end
function LinearAlgebra.rmul!(dat::VNamedD, x::Number)
	for d in dat
		rmul!(d,x)
	end
end




function update!(d::NamedD, fields::Vector{Symbol}, w::Array{Float64},)
	@assert length(w)==length(d.grid)
	@inbounds for f in fields
		@assert f ∈ names(d.d)[1]
		@inbounds for i in 1:d[:n]
			for it in 1:length(d.grid)
				d.d[f][it,i]=w[it]
			end
		end
	end
	return d
end
function update!(dat::VNamedD, fields::Vector{Symbol}, w::Array{Float64},)
	for d in dat
		update!(d,fields,w)
	end
end


"""
Compare if two `NamedD`'s  are equal
"""
function Base.isequal(dat1::NamedD, dat2::NamedD)
	return (dat1.grid==dat2.grid) && (dat1.d==dat2.d)
end
function Base.isequal(dat1::VNamedD, dat2::VNamedD)
	return all([isequal(dat1[i],dat2[i]) for i in 1:length(dat1)])
end

"""
Return if two `NamedD`'s have same dimensions and bounds.
"""
function issimilar(dat1::NamedD, dat2::NamedD)
	return isequal(dat1.grid, dat2.grid) && isequal(names(dat1.d), names(dat2.d)) && (size(dat1.d)==size(dat2.d)) 
end
function issimilar(dat1::VNamedD, dat2::VNamedD)
	@assert length(dat1)==length(dat2)
	return all([issimilar(dat1[i],dat2[i]) for i in 1:length(dat1)])
end

function Base.length(data::NamedD)
	return length(data.grid)*data.n*length(names(data.d)[1])
end

"""
Return a vec of data object sorted in the order
time, channels
"""
function Base.vec(data::NamedD)
	v=Vector{Float64}()
	for dd in data.d
		append!(v,dd)
	end
	return v
end
function Base.vec(data::VNamedD)
	v=Vector{Float64}()
	for d in data
		append!(v,vec(d))
	end
	return v
end


"""
Method to Devectorize data!
No memory allocations
"""
function Base.copyto!(d::NamedD, v::AbstractVector{Float64}, i0=1)
#	@assert length(d)==length(v)
	for dd in d.d
		for i in 1:d[:n]
			for it in 1:length(d.grid)
				dd[it,i]=v[i0]
				i0+=1
			end
		end
	end
	return i0
end
function Base.copyto!(d::VNamedD, v::AbstractVector{Float64})
	i0=1
	for i in 1:length(d)
		dd=d[i]
		i0=copyto!(dd,v,i0)
	end
	return d
end


function Base.copyto!(v::AbstractVector{Float64}, d::NamedD,i0=1)
	for dd in d.d
		for i in 1:d[:n]
			for it in 1:length(d.grid)
				v[i0]=dd[it,i]
				i0+=1
			end
		end
	end
	return i0
end
function Base.copyto!(v::AbstractVector{Float64}, d::VNamedD)
	i0=1
	for i in 1:length(d)
		dd=d[i]
		i0=copyto!(v,dd,i0)
	end
	return d
end


"""
Fill with randn values
"""
function Random.randn!(data::NamedD)
	for dd in data.d
		Random.randn!(dd)
	end
	return data
end

function Random.randn!(data::VNamedD)
	for d in data
		Random.randn!(d)
	end
end

"""
Copy `NamedD`'s, which are similar.
"""
function Base.copyto!(dataout::NamedD, data::NamedD)
	for i in names(dataout.d)[1]
		d1=dataout.d[i]
		d2=data.d[i]
		copyto!(d1,d2)
	end
end
function Base.copyto!(dataout::VNamedD, data::VNamedD)
	for i in 1:length(dataout)
		d1=dataout[i]
		d2=data[i]
		copyto!(d1,d2)
	end
end


"""
Returns bool depending on if input `data::NamedD` has all zeros or not.
"""
function Base.iszero(data::NamedD)
	return maximum(broadcast(maximum,abs,data.d)) == 0.0 ? true : false
end
function Base.iszero(data::VNamedD)
	return all([iszero(d) for d in data])
end


"""
Returns dot product of data.

# Arguments 

* `data1::NamedD` : data 1
* `data2::NamedD` : data 2

# Return

* dot product as `Float64`
"""
function LinearAlgebra.dot(data1::NamedD, data2::NamedD)
	dotd = 0.0;
	for iff in names(data1.d)[1]
		dd1=data1.d[iff]
		dd2=data2.d[iff]
		dotd += LinearAlgebra.dot(dd1,dd2)
	end
	return dotd
end

function LinearAlgebra.dot(data1::VNamedD, data2::VNamedD)
	dotd = 0.0;
	for i in 1:length(data1)
		dd1=data1[i]
		dd2=data2[i]
		dotd += LinearAlgebra.dot(dd1,dd2)  
	end
	return dotd
end


function Base.fill!(data::NamedD, k::Float64)
	for dd in data.d
		fill!(dd,k)
	end
end
function Base.fill!(data::VNamedD, k::Float64)
	for d in data
		fill!(d,k)
	end
end

function Base.reverse!(data::NamedD)
	for dd in data.d
		for i in 1:data.n
			dv=view(dd,:,i)
			reverse!(dv)
		end
	end
end

function Base.reverse!(data::VNamedD)
	for d in data
		reverse!(d)
	end
end



include("misfits.jl")
include("statistics.jl")
