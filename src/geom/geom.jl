


"""
Acquisiton has supersources, sources and receivers.
Each supersource has `ns` multiple sources that are 
injected (or active) simultaneously.
For each supersource,
a set of `nr` receivers are 
recording waves.

# Fields

* `sx::Vector{Vector{Float64,1},1}` : ``x`` positions of sources
* `sz::Vector{Vector{Float64,1},1}` : ``z`` positions of sources
* `rx::Vector{Vector{Float64,1},1}` : ``x`` positions of receivers
* `rz::Vector{Vector{Float64,1},1}` : ``z`` positions of receivers
* `nss::Int64` : number of supersources
* `ns::Vector{Int64,1}` : number of sources for every supersource
* `nr::Vector{Int64,1}` : number of receivers for every supersource 
"""
mutable struct Geomss
	sx::Vector{Float64}
	sz::Vector{Float64}
	rx::Vector{Float64}
	rz::Vector{Float64}
	ns::Int64 # change to +ve only integer later
	nr::Int64 # change to +ve only integer later
	"adding conditions that are to be false while construction"
	Geomss(sx, sz, rx, rz, ns, nr) = 
		any([
                  length(sx) != ns,
                  length(sz) != ns,
                  length(rx) != nr,
                  length(rz) != nr,
		  ]) ? 
		error("Geomss construct") : new(sx, sz, rx, rz, ns, nr)
end # type

Geom=Array{Geomss,1}
 

"""
Randomly place a given number of source and receivers in a `mgrid`.
"""
function Geomss(mgrid::AbstractArray{T}, s::Srcs, r::Recs) where {T<:StepRangeLen}
	xmin=mgrid[2][1]
	xmax=mgrid[2][end]
	zmin=mgrid[1][1]
	zmax=mgrid[1][end]
	return Geomss(Random.rand(Uniform(xmin,xmax),s.n), Random.rand(Uniform(zmin,zmax),s.n),
		    Random.rand(Uniform(xmin,xmax),r.n), Random.rand(Uniform(zmin,zmax),r.n), s.n, r.n)
end

function Array{Geomss,1}(mgrid::AbstractVector{T}, ss::SSrcs, s::Vector{Srcs}, r::Vector{Recs}) where {T<:StepRangeLen}   
	return [Geomss(mgrid,s[i],r[i]) for i in 1:ss.n]
end
function Array{Geomss,1}(mgrid::AbstractVector{T}, ss::SSrcs, s::Srcs, r::Recs)  where {T<:StepRangeLen}   
	return Array{Geomss,1}(mgrid, ss, fill(s,ss.n), fill(r,ss.n))
end

"Compare if two `Geom` variables are equal"
function Base.isequal(geom1::Geomss, geom2::Geomss, attrib=:all)
	srcfnames = [:sx, :sz, :ns]
	recfnames = [:rx, :rz, :nr]
	srcvec=[(isequal(getfield(geom1, name),getfield(geom2, name))) for name in srcfnames]
	recvec=[(isequal(getfield(geom1, name),getfield(geom2, name))) for name in recfnames]
	if(attrib == :sources)
		return all(srcvec)
	elseif(attrib == :receivers)
		return all(recvec)
	else
		return all(vcat(srcvec, recvec))
	end
end
function Base.isequal(geom1::Geom, geom2::Geom, attrib=:all)
	return (length(geom1)==length(geom2)) ? all([isequal(geom1[iss],geom2[iss]) for iss in 1:length(geom1)]) : false
end

function Base.show(io::Base.IO, geomss::Geomss)
	print(io, "\t> acquisition with ",geomss.ns," source(s) and ",geomss.nr, " receiver(s)")
end

"""
Output `n` interpolated 2-D points between `p1` and `p2`. 
"""
function spread(n, p1::Vector{T1}, p2::Vector{T2}, rand_flag::Bool=false) where {T1<:Real,T2<:Real}
	@assert (length(p1)==2) && (length(p2)==2)

	if(n==1)
		return [p1[1]], [p1[2]]
	else
		knots = ([1,n],)
		itpz = Interpolations.interpolate(knots, [p1[1],p2[1]],  Gridded(Linear()));
		itpx = Interpolations.interpolate(knots, [p1[2],p2[2]],  Gridded(Linear()));
		return [itpz(i) for i in 1:n], [itpx(i) for i in 1:n]
	end
end

function spread(n, p0::Vector{T1}, rad::T2, angles::Vector{T3}, rand_flag::Bool=false)  where {T1<:Real,T2<:Real,T3<:Real}  
	th,ni=spread(n, [angles[1], 1],[angles[2], n], rand_flag)
	x=zeros(n); z=zeros(n)
	for i in 1:n
		x[i]=p0[2]+rad*cos(th[i])
		z[i]=p0[1]+rad*sin(th[i])
	end
	return z,x
end

"""
Update source positions
"""
function update!(geomss::Geomss, ::Srcs, args...) where {T<:Real}
	z,x=spread(geomss.ns,args...)
	copyto!(geomss.sx,x)
	copyto!(geomss.sz,z)
	return geomss
end

"""
Update receiver positions
"""
function update!(geomss::Geomss, ::Recs, args...)
	z,x=spread(geomss.nr,args...)
	copyto!(geomss.rx,x)
	copyto!(geomss.rz,z)
	return geomss
end

"""
Update supersource positions
"""
function update!(geom::Geom, ::SSrcs, args...)
	z,x=spread(length(geom), args...)
	for iss in 1:length(geom)
		update!(geom[iss], Srcs(), [z[iss],x[iss]], [z[iss],x[iss]])    
	end
	return geom
end

"""
Update receiver positions in all supersources
"""
function update!(geom::Geom, ::Recs, args...)
	for iss in 1:length(geom)
		update!(geom[iss], Recs(), args...)
	end
	return geom
end


"""
Check if all the sources and receivers in `Geom` are within the model 

# Return

* `true` if all the positions within the model, `false` otherwise
"""
function Base.in(geom::Geom, mgrid::AbstractVector{T}) where {T<:StepRangeLen}     
	xmin, zmin, xmax, zmax = mgrid[2][1], mgrid[1][1], mgrid[2][end], mgrid[1][end]
	checkvec=fill(false, length(geom))
	for iss=1:length(geom)
		checkvec[iss] = any(
		      vcat(
			   (((geom[iss].sx .- xmin).*(xmax .- geom[iss].sx)) .< 0.0),
			   (((geom[iss].sz .- zmin).*(zmax .- geom[iss].sz)) .< 0.0),
			   (((geom[iss].rx .- xmin).*(xmax .- geom[iss].rx)) .< 0.0),
			   (((geom[iss].rz .- zmin).*(zmax .- geom[iss].rz)) .< 0.0),
			   isnan.(geom[iss].sx),
			   isnan.(geom[iss].sz),
			   isnan.(geom[iss].rx),
			   isnan.(geom[iss].rz),
		     ) )
	end
	return !(any(checkvec))
end


#=
"""
Return a sparse ACQ matrix, for a given geom
data=ACQ*snapshot
"""
function ACQmat(geom::Geom,mgrid,iss=1)
	nz,nx=length.(mgrid)
	ACQ=spzeros(geom.nr[iss],prod(length.(mgrid)))
	for ir = 1:geom.nr[iss]
		irx=argmin(abs.(mgrid[2].-geom.rx[iss][ir]))
		irz=argmin(abs.(mgrid[1].-geom.rz[iss][ir]))
		ACQ[ir,irz+(irx-1)*nz]=1.0
	end
	return ACQ
end

=#
include("gallery.jl")

#=

"""
Advance either source or receiver array in an acquisition geometry in horizontal or vertical directions.

# Arguments
* `geom::Geom` : acquisition geometry that is updated
* `advances::Vector{Float64}=[[0.,0.], [0.,0.,]]` : source and receiver advancements
"""
function Geom_advance(geom::Geom, advances::Vector{Vector{Float64}}=[[0.,0.], [0.,0.]]) 
	geom_out=deepcopy(geom)
	for iss=1:geom.nss
		geom_out.sx[iss][:] .+= advances[1][2]
		geom_out.sz[iss][:] .+= advances[1][1]
		geom_out.rx[iss][:] .+= advances[2][2]
		geom_out.rz[iss][:] .+= advances[2][1]
	end
	return geom_out
end

"""
return a vector of the order 
"""
function Geom_getvec(geom::Vector{Geom}, field::Symbol)
	np = length(geom);
	vect = [getfield(geom[ip],field) for ip=1:np]
	return vec(vcat(vcat(vect...)...))
end


function save(geom, folder)
	!(isdir(folder)) && error("invalid directory")
	for iss in 1:geom.nss
		file=joinpath(folder, string("s", iss, ".csv"))
		CSV.write(file,DataFrame(hcat(geom.sx[iss], geom.sz[iss])))
		file=joinpath(folder, string("r", iss, ".csv"))
		CSV.write(file,DataFrame(hcat(geom.rx[iss], geom.rz[iss])))
	end
end


"""
Appends the input a vector of acquisition geometries.
"Merge" `Geom` objects in `geomv` to `geom1`.
Typically you have one acquisition, 
you change it and want to create the macro acquisition
where the supersources of both acquisitions are appended.
"""
function Geom_add!(geom1::Geom, geomv::Vector{Geom})
	for igeom=1:length(geomv)
		append!(geom1.sx, geomv[igeom].sx)
		append!(geom1.sz, geomv[igeom].sz)
		append!(geom1.rx, geomv[igeom].rx)
		append!(geom1.rz, geomv[igeom].rz)
		append!(geom1.ns, geomv[igeom].ns)
		append!(geom1.nr, geomv[igeom].nr)
		geom1.nss += geomv[igeom].nss
	end
	return geom1
end


"""
Adds input positions as either sources or receivers of every supershot.
"""
function Geom_add!(geom::Geom, zpos::Vector{Float64}, xpos::Vector{Float64}, attrib::Symbol)
	length(xpos) == length(zpos) ? nothing : error("same length vectors needed")
	for iss=1:geom.nss
		if(attrib == :receivers)
			geom.rx[iss] = vcat(geom.rx[iss][:], xpos)
			geom.rz[iss] = vcat(geom.rz[iss][:], zpos)
			geom.nr[iss] += length(xpos)
		elseif(attrib == :sources)
			geom.sx[iss] = vcat(geom.sx[iss][:], xpos)
			geom.sz[iss] = vcat(geom.sz[iss][:], zpos)
			geom.ns[iss] += length(xpos)
		end
	end
end




"""
A fixed spread acquisition has same set of sources and 
receivers for each supersource.
This method constructs a 
fixed spread acquisition geometry using either a
horizontal or vertical array of supersources/ receivers.
Current implementation has only one source for every supersource.

# Arguments 

* `ssmin::Float64` : minimum coordinate for sources
* `ssmax::Float64` : maximum coordinate for sources
* `ss0::Float64` : consant coordinate for sources
* `rmin::Float64` : minimum coordinate for receivers
* `rmax::Float64` : maximum coordinate for receivers
* `r0::Float64` : consant coordinate for receivers
* `nss::Int64` : number of supersources
* `nr::Int64` : number of receivers
* `ssattrib::Symbol=:horizontal` : supersource array kind
  `=:vertical` : vertical array of supersources
  `=:horizontal` horizontal array of supersources
* `rattrib::Symbol=:horizontal` : receiver array kind
  `=:vertical` : vertical array of receivers
  `=:horizontal` horizontal array of receivers
* `rand_flags::Vector{Bool}=[false, false]` : decide placement of supersources and receivers 
  `=[true, false]` : randomly place supersources for regularly spaced receivers
  `=[true, true]` : randomly place supersources and receivers
  `=[false, false]` : regularly spaced supersources and receivers
  `=[false, true]` : randomly place receivers for regularly spaced supersources 

# Return
* a fixed spread acquisition geometry `Geom`
"""
function Geom_fixed(
	      ssmin::Real,
	      ssmax::Real,
	      ss0::Real,
	      rmin::Real,
	      rmax::Real,
	      r0::Real,
	      nss::Int64,
	      nr::Int64,
	      ssattrib::Symbol=:horizontal,
	      rattrib::Symbol=:horizontal,
 	      rand_flags::Vector{Bool}=[false, false];
	      ssα::Real=0.0,
	      rα::Real=0.0,
	      ns::Vector{Int64}=ones(Int,nss),
	      srad::Real=0.0
	     )

	ssmin=Float64(ssmin); rmin=Float64(rmin)
	ssmax=Float64(ssmax); rmax=Float64(rmax)
	ss0=Float64(ss0); r0=Float64(r0)
	ssα=Float64(ssα)*pi/180.
	rα=Float64(rα)*pi/180.


	ssarray = isequal(nss,1) ? fill(ssmin,1) : (rand_flags[1] ? 
					     Random.rand(Uniform(minimum([ssmin,ssmax]),maximum([ssmin,ssmax])),nss) : range(ssmin,stop=ssmax,length=nss))
	if(ssattrib==:horizontal)
		ssz = ss0.+(ssarray.-minimum(ssarray)).*sin(ssα)/cos(ssα); ssx=ssarray
	elseif(ssattrib==:vertical)
		ssx = ss0.+(ssarray.-minimum(ssarray)).*sin(ssα)/cos(ssα); ssz=ssarray
	else
		error("invalid ssattrib")
	end

	rarray = isequal(nr,1) ? fill(rmin,1) : (rand_flags[2] ? 
					  Random.rand(Uniform(minimum([rmin,rmax]),maximum([rmin,rmax])),nr) : range(rmin,stop=rmax,length=nr))
	if(rattrib==:horizontal)
		rz = r0.+(rarray.-minimum(rarray)).*sin(rα)/cos(rα); rx = rarray
	elseif(rattrib==:vertical)
		rx = r0.+(rarray.-minimum(rarray)).*sin(rα)/cos(rα); rz = rarray
	else
		error("invalid rattrib")
	end
	rxall = [rx for iss=1:nss];
	rzall = [rz for iss=1:nss];
	ssxall=[zeros(ns[iss]) for iss=1:nss];
	sszall=[zeros(ns[iss]) for iss=1:nss];
	for iss in 1:nss
		for is in 1:ns[iss]
			θ=Random.rand(Uniform(-Float64(pi),Float64(pi)))
			r = iszero(srad) ? 0.0 : Random.rand(Uniform(0,srad))
			x=r*cos(θ)
			z=r*sin(θ)
			ssxall[iss][is]=ssx[iss]+x
			sszall[iss][is]=ssz[iss]+z
		end
	end
	return Geom(ssxall, sszall, rxall, rzall, nss, ns, fill(nr,nss))
end




=#
