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
modify exixting acquisitin geometries
same source but receivers at the border
"""
function Geom(acqgeom::Geom,
	      mgrid::Grid.M2D,
	      attrib::Symbol
	     )
	if(attrib == :recborder)
		rx = append!(mgrid.x, 
		     append!(mgrid.x[end]*ones(mgrid.nz-1),
		     append!(mgrid.x[end-1:-1:1],
		             mgrid.x[1]*ones(mgrid.nz-2))))
		rz = append!(mgrid.z[1]*ones(mgrid.nx), 
		     append!(mgrid.z[2:end],
		     append!(mgrid.z[end]*ones(mgrid.nx-1),
		             mgrid.z[end-1:-1:2])))
		nr = size(rx,1);
		sx = acqgeom.sx;
		sz = acqgeom.sz;
		ns = acqgeom.ns;
		nss = acqgeom.nss;
		return Geom(sx, sz, repmat(rx,1,nss), repmat(rz,1,nss), nss, ns, fill(nr,nss))
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

"swap the positions of sources and recievers"
function Geom_swapsr(geomin::Geom)
	geomout = geomin;
	geomout.rx = geomin.sx; geomout.rz = geomin.sz;
	geomout.sx = geomin.rx; geomout.sz = geomin.rz;
	geomout.ns = geomin.nr; geomout.nr = geomin.ns;

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
* `wav::Array{Float64}` : wavelets in time domain
* `tgrid::Grid.M1D` : time grid 
"""
type Src
	nss::Int64
	ns::Array{Int64}
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
	     wav::Array{Float64},
	     tgrid::Grid.M1D
	     )
	return Src(nss, fill(ns, nss), repeat(wav, inner=(1,ns,nss)), tgrid)
end


"""
Function that returns Src after time reversal
"""
function Src_tr(src::Src)
	return Src(src.nss,src.ns,src.wav[end:-1:1,:,:],src.tgrid)
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
