"""
Acquisition geometry for each supersource.
"""
mutable struct AGeomss
    s::NamedStack{Vector{Float64}} # x, y, z coordinates
    r::NamedStack{Vector{Float64}} # x, y, z coordinates
    ns::Int64 # change to +ve only integer later
    nr::Int64 # change to +ve only integer later
    "adding conditions that are to be false while construction"
    AGeomss(s, r, ns, nr) =
        (
            any(AxisArrays.names(s)[1] != AxisArrays.names(r)[1]) |
            any([length(ss) != ns for ss in s]) |
            any([length(rr) != nr for rr in r])
        ) ? error("AGeomss construct") : new(s, r, ns, nr)
end # type

"""
Randomly place a given number of source and receivers in a `mgrid`.
"""
function AGeomss(mgrid::Vector{T}, s::Srcs, r::Recs) where {T<:StepRangeLen}
    nd = length(mgrid)
    names = dim_names(nd)
    mins = [m[1] for m in mgrid]
    maxs = [m[end] for m in mgrid]
    return AGeomss(
        NamedArray([Random.rand(Uniform(mins[i], maxs[i]), s.n) for i = 1:nd], (names,)),
        NamedArray([Random.rand(Uniform(mins[i], maxs[i]), r.n) for i = 1:nd], (names,)),
        s.n,
        r.n,
    )
end

"""
```julia
ageom=AGeom(mgrid, SSrcs(2), [Srcs(3), Srcs(2)], [Recs(10), Recs(20)])
```
Initialize an acquisition, with 2 supersources, by randomly placing a 
sources and receivers in `mgrid`.
"""
function Vector{AGeomss}(
    mgrid::Vector{T},
    ss::SSrcs,
    s::Vector{Srcs},
    r::Vector{Recs},
) where {T<:StepRangeLen}
    return [AGeomss(mgrid, s[i], r[i]) for i = 1:ss.n]
end

"""
```julia
ageom=AGeom(mgrid, SSrcs(10), Srcs(10), Recs(20))
```
Initialize an acquisition, with 10 supersources, by randomly placing a 
sources and receivers in `mgrid`. In this case, all the supersources have same source and receiver positions.
"""
function Vector{AGeomss}(
    mgrid::Vector{T},
    ss::SSrcs,
    s::Srcs,
    r::Recs,
) where {T<:StepRangeLen}
    return Array{AGeomss,1}(mgrid, ss, fill(s, ss.n), fill(r, ss.n))
end


AGeom = Array{AGeomss,1}


function Base.isequal(ageom1::AGeomss, ageom2::AGeomss, attrib = :all)
    @assert dims(ageom1) == dims(ageom2) "attempt to compare acquisitions of different dimensions"
    names = [:sx, :sz, :ns]
    recfnames = [:rx, :rz, :nr]
    src =
        all([(isequal(ageom1.s[name], ageom2.s[name])) for name in AxisArrays.names(ageom1.s)[1]]) & (ageom1.ns == ageom2.ns)
    rec =
        all([
            (isequal(getfield(ageom1, name), getfield(ageom2, name))) for name in recfnames
        ]) & (ageom1.nr == ageom2.nr)
    if (attrib == :sources)
        return src
    elseif (attrib == :receivers)
        return rec
    else
        return src & rec
    end
end

"""
```julia
isequal(ageom1, ageom2)
```
Assert if `ageom1` equals `ageom2`.
"""
function Base.isequal(ageom1::AGeom, ageom2::AGeom, attrib = :all)
    return (length(ageom1) == length(ageom2)) ?
           all([isequal(ageom1[iss], ageom2[iss], attrib) for iss = 1:length(ageom1)]) :
           false
end


function dim_names(ageomss::AGeomss)
    dN = dim_names(length(ageomss.s))
    @assert all(map(in(dN), AxisArrays.names(ageomss.s)[1]))
    @assert all(map(in(dN), AxisArrays.names(ageomss.r)[1]))
    return dN
end

function dim_names(ageom::Vector{AGeomss})
    dN = dim_names.(ageom)
    @assert all(y -> y == dN[1], dN)
    return dN[1]
end



function Base.show(io::Base.IO, ::MIME"text/plain", ageomss::AGeomss)
    print(
        io,
        "Acquisition geometry of a supersource\n",
        "  ├───────── ndims: ",
        length(dim_names(ageomss)),
        '\n',
        "  ├─── # source(s): ",
        ageomss.ns,
        '\n',
        "  ├─ # receiver(s): ",
        ageomss.nr,
        '\n',
    )
end


function Base.show(io::Base.IO, ::MIME"text/plain", ageom::AbstractVector{AGeomss})
    print(
        io,
        "Acquisition geometry with $(length(ageom)) supersource(s)\n",
        "  ├───────── ndims: ",
        length(dim_names(ageom)),
        '\n',
        "  ├─── # source(s): ",
        getfield.(ageom, :ns),
        '\n',
        "  ├─ # receiver(s): ",
        getfield.(ageom, :nr),
        '\n',
    )
end


"""
Output `n` interpolated sources/ receivers between positions `p1` and `p2`. 
"""
function spread(n, p1::Vector{T1}, p2::Vector{T2}) where {T1<:Real,T2<:Real}
    @assert (length(p1) == length(p2)) # check if 2D or 3D
    nd = length(p1)
    if (n == 1)
        return Tuple([[p1[i]] for i = 1:nd])
    else
        knots = ([1, n],)
        itp = [
            interpolate(
                knots,
                [Float64(p1[i]), Float64(p2[i])],
                Gridded(Linear()),
            ) for i = 1:nd
        ]
        return Tuple([[it(i) for i = 1:n] for it in itp])
    end
end

"""
Spread using difference sources/ receivers in each dimension given by `n_per_dim`.
"""
function spread(
    n,
    p1::Vector{T1},
    p2::Vector{T2},
    n_per_dim::Vector{T3},
) where {T1<:Real,T2<:Real,T3<:Int}
    @assert length(p1) == length(p2) == length(n_per_dim)
    nd = length(p1)
    @assert prod(n_per_dim) == n
    I = collect(
        Iterators.product(
            [spread(n_per_dim[i], [p1[i]], [p2[i]])[1] for i = 1:length(n_per_dim)]...,
        ),
    )
    return Tuple([[I[i][id] for i = 1:n] for id = 1:nd])
end

"""
theta=0 means in x-z plane
"""
function spread(
    n,
    p0::Vector{T1},
    rad::T2,
    phi::Vector{T4},
    theta::Vector{T3} = [0, 0],
) where {T1<:Real,T2<:Real,T3<:Real,T4<:Real}
    phi = (phi[1] ≈ phi[2]) ? fill(phi[1], n) : rand(Uniform(phi[1], phi[2]), n)
    theta = (theta[1] ≈ theta[2]) ? fill(theta[1], n) : rand(Uniform(theta[1], theta[2]), n)
    nd = length(p0)
    p = [zeros(n) for i = 1:nd]

    for i = 1:n
        p[1][i] = p0[1] + rad * cos(phi[i])
        if (nd > 1)
            p[2][i] = p0[2] + rad * sin(phi[i]) * sin(theta[i])
        end
        if (nd > 2)
            p[3][i] = p0[3] + rad * sin(phi[i]) * cos(theta[i])
        end
    end
    return Tuple(p)
end

"""
```julia
update!(ageom[1], Srcs(10), p1, p2)
```
Update source geometry for the first supersource, with 10 sources placed regularly between `p1=[z1,x1]` and `p2=[z2,x2]`

```julia
update!(ageom[2], Srcs(10), p0, rad, angles)
```
Update source geometry for the second supersource, with 10 sources placed equidistant from `p0=[z0,x0]` given radius `rad`
and angles `angles=[0,2pi]`.
"""
function update!(ageomss::AGeomss, ::Srcs, args...)
    p = spread(ageomss.ns, args...)
    for id = 1:length(dim_names(ageomss))
        copyto!(ageomss.s[id], p[id])
    end
    return ageomss
end

"""
```julia
update!(ageom[1], Recs(10), args...)
```
Similar to updating source positions, but for receivers.
"""
function update!(ageomss::AGeomss, ::Recs, args...)
    p = spread(ageomss.nr, args...)
    for id = 1:length(dim_names(ageomss))
        copyto!(ageomss.r[id], p[id])
    end
    return ageomss
end

"""
```julia
update!(ageom, SSrcs(), p1, p2)
update!(ageom, SSrcs(), p0, rad, angles)
```
Update all source positions in each supersources. 
"""
function update!(ageom::AGeom, ::SSrcs, args...)
    p = spread(length(ageom), args...)
    pz = collect(zip(p...))
    for iss = 1:length(ageom)
        update!(ageom[iss], Srcs(), collect(pz[iss]), collect(pz[iss]))
    end
    return ageom
end

"""
```julia
update!(ageom, Recs(), p1, p2)
update!(ageom, Recs(), p0, rad, angles)
```
Update receiver positions for all supersources.
"""
function update!(ageom::AGeom, ::Recs, args...)
    for iss = 1:length(ageom)
        update!(ageom[iss], Recs(), args...)
    end
    return ageom
end


"""
```julia
ageom ∈ mgrid
in(ageom, mgrid)
```
Assert if all the sources and receivers in `AGeom` are within bounds of `mgrid`.
"""
function Base.in(ageom::AGeom, mgrid::AbstractVector{T}) where {T<:StepRangeLen}
    # xmin, zmin, xmax, zmax = mgrid[2][1], mgrid[1][1], mgrid[2][end], mgrid[1][end]
    checkvec = fill(false, length(ageom))
    for iss = 1:length(ageom)
        nd = length(dim_names(ageom[iss]))
        @assert nd == length(mgrid)
        names = dim_names(nd)
        checkvec[iss] = any(
            hcat(
                [
                    vcat(
                        (
                            (
                                (ageom[iss].s[dn] .- mgrid[id][1]) .*
                                (mgrid[id][end] .- ageom[iss].s[dn])
                            ) .< 0.0
                        ),
                        (
                            (
                                (ageom[iss].r[dn] .- mgrid[id][1]) .*
                                (mgrid[id][end] .- ageom[iss].r[dn])
                            ) .< 0.0
                        ),
                        isnan.(ageom[iss].s[dn]),
                        isnan.(ageom[iss].r[dn]),
                    ) for (id, dn) in enumerate(names)
                ]...,
            ),
        )
    end
    return !(any(checkvec))
end


"""
```julia
A=SparseArrays.SparseMatrixCSC(ageom[1],mgrid)
D=A*P
```
Return a sparse matrix that restricts the field `P` on `mgrid` to receiver positions, to give `D`.
"""
function SparseArrays.SparseMatrixCSC(
    ageomss::AGeomss,
    mgrid::AbstractVector{T},
) where {T<:StepRangeLen}
    @assert length(mgrid) == 2 #  need to update this routine for 3D acq
    nz, nx = length.(mgrid)
    ACQ = spzeros(ageomss.nr, prod(length.(mgrid)))
    for ir = 1:ageomss.nr
        irx = argmin(abs.(mgrid[2] .- ageomss.r[:x][ir]))
        irz = argmin(abs.(mgrid[1] .- ageomss.r[:z][ir]))
        ACQ[ir, irz+(irx-1)*nz] = 1.0
    end
    return ACQ
end


"""
Return `[sz-rz, sy-ry, sx-rx]`
"""
offset(ageomss::AGeomss, is, ir) =
    [ageomss.s[d][is] - ageomss.r[d][ir] for d in dim_names(ageomss)]

#=

"""
Advance either source or receiver array in an acquisition ageometry in horizontal or vertical directions.

# Arguments
* `ageom::AGeom` : acquisition ageometry that is updated
* `advances::Vector{Float64}=[[0.,0.], [0.,0.,]]` : source and receiver advancements
"""
function AGeom_advance(ageom::AGeom, advances::Vector{Vector{Float64}}=[[0.,0.], [0.,0.]]) 
	ageom_out=deepcopy(ageom)
	for iss=1:ageom.nss
		ageom_out.s[:x][iss][:] .+= advances[1][2]
		ageom_out.s[:z][iss][:] .+= advances[1][1]
		ageom_out.r[:x][iss][:] .+= advances[2][2]
		ageom_out.r[:z][iss][:] .+= advances[2][1]
	end
	return ageom_out
end

"""
return a vector of the order 
"""
function AGeom_getvec(ageom::Vector{AGeom}, field::Symbol)
	np = length(ageom);
	vect = [getfield(ageom[ip],field) for ip=1:np]
	return vec(vcat(vcat(vect...)...))
end


function save(ageom, folder)
	!(isdir(folder)) && error("invalid directory")
	for iss in 1:ageom.nss
		file=joinpath(folder, string("s", iss, ".csv"))
		CSV.write(file,DataFrame(hcat(ageom.s[:x][iss], ageom.s[:z][iss])))
		file=joinpath(folder, string("r", iss, ".csv"))
		CSV.write(file,DataFrame(hcat(ageom.r[:x][iss], ageom.r[:z][iss])))
	end
end


"""
Appends the input a vector of acquisition ageometries.
"Merge" `AGeom` objects in `ageomv` to `ageom1`.
Typically you have one acquisition, 
you change it and want to create the macro acquisition
where the supersources of both acquisitions are appended.
"""
function AGeom_add!(ageom1::AGeom, ageomv::Vector{AGeom})
	for iageom=1:length(ageomv)
		append!(ageom1.s[:x], ageomv[iageom].s[:x])
		append!(ageom1.s[:z], ageomv[iageom].s[:z])
		append!(ageom1.r[:x], ageomv[iageom].r[:x])
		append!(ageom1.r[:z], ageomv[iageom].r[:z])
		append!(ageom1.ns, ageomv[iageom].ns)
		append!(ageom1.nr, ageomv[iageom].nr)
		ageom1.nss += ageomv[iageom].nss
	end
	return ageom1
end


"""
Adds input positions as either sources or receivers of every supershot.
"""
function AGeom_add!(ageom::AGeom, zpos::Vector{Float64}, xpos::Vector{Float64}, attrib::Symbol)
	length(xpos) == length(zpos) ? nothing : error("same length vectors needed")
	for iss=1:ageom.nss
		if(attrib == :receivers)
			ageom.r[:x][iss] = vcat(ageom.r[:x][iss][:], xpos)
			ageom.r[:z][iss] = vcat(ageom.r[:z][iss][:], zpos)
			ageom.nr[iss] += length(xpos)
		elseif(attrib == :sources)
			ageom.s[:x][iss] = vcat(ageom.s[:x][iss][:], xpos)
			ageom.s[:z][iss] = vcat(ageom.s[:z][iss][:], zpos)
			ageom.ns[iss] += length(xpos)
		end
	end
end




"""
A fixed spread acquisition has same set of sources and 
receivers for each supersource.
This method constructs a 
fixed spread acquisition ageometry using either a
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
* a fixed spread acquisition ageometry `AGeom`
"""
function AGeom_fixed(
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
	return AGeom(ssxall, sszall, rxall, rzall, nss, ns, fill(nr,nss))
end




=#
