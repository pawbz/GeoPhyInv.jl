module Acquisition

import SIT.Grid
import SIT.Wavelets

"""
"""
type Geom
	"``x`` positions of sources"
	sx::Array{Float64}
	"``z`` positions of sources"
	sz::Array{Float64}
	"``x`` positions of receivers"
	rx::Array{Float64}
	"``z`` positions of receivers"
	rz::Array{Float64}
	"number of supershources"
	nss::Int64
	"number of sources per each supersource"
	ns::Array{Int64}
	"number of receivers per each supersource"
	nr::Array{Int64}
end


"""
Modify input `Geom` such that the output has
either sources or receivers on the boundary of 
`mgrid`.

# Arguments
* `acqgeom::Geom` : input geometry
* `mgrid::Grid.M2D` : to determine the boundary
* `attrib::Symbol` : either `:srcborder` or `:recborder` 
"""
function Geom(acqgeom::Geom,
	      mgrid::Grid.M2D,
	      attrib::Symbol
	     )

	if((attrib == :recborder) | (attrib == :srcborder))
		bz, bx, nb = Grid.M2D_border(mgrid, 3, :inner)
	end
	"change the position of receivers to the boundary"
	if(attrib == :recborder)
		nr = nb; 		sx = acqgeom.sx;
		sz = acqgeom.sz;		ns = acqgeom.ns;
		nss = acqgeom.nss;
		return Geom(sx, sz, repmat(bx,1,nss), repmat(bz,1,nss), nss, ns, fill(nr,nss))
	"change the position of source (not supersources) to the boundary for back propagation"
	elseif(attrib == :srcborder)
		ns = nb;		rx = acqgeom.rx;
		rz = acqgeom.rz;		nr = acqgeom.nr;
		nss = acqgeom.nss;
		return Geom(repmat(bx,1,nss), repmat(bz,1,nss), rx, rz, nss, fill(ns,nss), nr)
	else
		error("invalid attrib")
	end
end


"""
return acquisition geometry depending 
on either horizontal or vertical array
It has only one source for every supersource
"""
function Geom(
	      smin::Float64,
	      smax::Float64,
	      s0::Float64,
	      rmin::Float64,
	      rmax::Float64,
	      r0::Float64,
	      nss::Int64,
	      nr::Int64,
	      sattrib::Symbol=:horizontal,
	      rattrib::Symbol=:horizontal
	     )
	if(sattrib==:horizontal)
		sz = fill(s0,nss)
		isequal(nss,1) ? sx = fill(smin,1) : sx = linspace(smin,smax,nss)
	elseif(sattrib==:vertical)
		sx = fill(s0,nss)
		isequal(nss,1) ? sz = fill(smin,1) : sz = linspace(smin,smax,nss)
	else
		error("invalid sattrib")
	end

	if(rattrib==:horizontal)
		rz = fill(r0,nr)
		isequal(nr,1) ? rx = fill(rmin,1) : rx = linspace(rmin,rmax,nr)
	elseif(rattrib==:vertical)
		rx = fill(r0,nr)
		isequal(nr,1) ? rz = fill(rmin,1) : rz = linspace(rmin,rmax,nr)
	else
		error("invalid rattrib")
	end
	rxall = repmat(rx, 1, nss);
	rzall = repmat(rz, 1, nss);
	sxall = reshape(sx, 1, nss);
	szall = reshape(sz, 1, nss);
	return Geom(sxall, szall, rxall, rzall, nss, fill(1,nss), fill(nr,nss))
end

"""
Circular acquisition. The sources and receivers can also be placed on a circle of radius
`rad`. The origin of the circle is at `loc`. 
This geometry is unrealistic, but useful for testing.
Receivers are placed such that the limits 
of the angular offset are given by `θlim`

# Arguments
`θlim::Vector{Float64}=[0.,2*pi]` : range of angular offset between source and receiver

"""
function Geom_circ(;
		   nss::Int64=10,
		   nr::Int64=10,
		   loc::Vector{Float64}=[0.,0.],
		   rad::Float64=100.,
		   θlim::Vector{Float64}=[0.,2*pi]
	       )
	rxall = zeros(nr,nss)
	rzall = zeros(nr,nss)
	sx = zeros(nss)
	sz = zeros(nss)
	δθr = abs((θlim[2] - θlim[1])) / nr;
	for iss =1:nss
		θs = (iss)*2.*pi / (nss);
		sx[iss] = rad*cos(θs) + loc[2]
		sz[iss] = rad*sin(θs) + loc[1]
		rxall[:,iss] = [rad * cos(θs + θlim[1] + (ir-1)*δθr) + loc[2] for ir in 1:nr]
		rzall[:,iss] = [rad * sin(θs + θlim[1] + (ir-1)*δθr) + loc[1] for ir in 1:nr]
	end
	sxall = reshape(sx, 1, nss);
	szall = reshape(sz, 1, nss);
	return Geom(sxall, szall, rxall, rzall, nss, fill(1,nss), fill(nr,nss))


end

"""
Modify the input acquisition geometry 
such that the adjoint source time functions can 
be propagated from the receiver positions.
The number of supersources will remain the same.
All the recievers will be fired as simultaneous sources.
"""
function Geom_adj(geomin::Geom)
	geomout = geomin;
	geomout.sx = geomin.rx; geomout.sz = geomin.rz;
	geomout.ns = geomin.nr; 

	return geomout
end


"""
return a vector of the order 
"""
function get_vecGeom(field::Symbol, geom::Array{Geom})

	np = length(geom);
	vecGeom = vec(getfield(geom[1],field));
	if(np > 1)
		for ip = 2:np
			vecGeom = vcat(vecGeom,vec(getfield(geom[ip],field)))
		end
	end
	return vecGeom
end


"""
Data type for the sources used.

# Fields
* `nss::Int64` : number of supersources
* `ns::Array{Int64}` : number of sources for each supersource
* `nfield::Int64` : number of fields
* `wav::Array{Float64}` : wavelets in time domain
* `tgrid::Grid.M1D` : time grid 
"""
type Src
	nss::Int64
	ns::Array{Int64}
	nfield::Int64
	wav::Array{Float64}
	tgrid::Grid.M1D
end


"""
Constructor for `Src` data type.
repeat same source wavelet for all sources and supersources

# Arguments
* `nss::Int64` : number of supersources
* `ns::Int64` : number of sources
* `wav::Array{Float64}`
"""
function Src(nss::Int64, 
	     ns::Int64, 
	     nfield::Int64,
	     wav::Array{Float64},
	     tgrid::Grid.M1D
	     )
	return Src(nss, fill(ns, nss), nfield, repeat(wav, inner=(1,ns,nss,nfield)), tgrid)
end


"""
Function that returns Src after time reversal
"""
function Src_tr(src::Src)
	return Src(src.nss,src.ns,src.nfield,src.wav[end:-1:1,:,:,:],src.tgrid)
end


"""
return a vector of the order 
"""
function get_vecSrc(field::Symbol, src::Array{Src})

	np = length(src);
	vecSrc = vec(getfield(src[1],field));
	if(np > 1)
		for ip = 2:np
			vecSrc = vcat(vecSrc,vec(getfield(src[ip],field)))
		end
	end
	return vecSrc
end

end # module
