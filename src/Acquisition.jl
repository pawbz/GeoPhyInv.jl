__precompile__()

module Acquisition

import SIT.Grid
import SIT.Models
import SIT.Wavelets
using Distributions

"""
Acquisiton of supersources, sources and receivers.
Each supersource has multiple sources that are 
injected (active) simultaneously.
Each supersource 
has a set of receivers that 
record waves. 

# Fields

* `sx::Vector{Vector{Float64,1},1}` : ``x`` positions of sources
* `sz::Vector{Vector{Float64,1},1}` : ``z`` positions of sources
* `rx::Vector{Vector{Float64,1},1}` : ``x`` positions of receivers
* `rz::Vector{Vector{Float64,1},1}` : ``z`` positions of receivers
* `nss::Int64` : number of supersources
* `ns::Vector{Int64,1}` : number of sources for every supersource
* `nr::Vector{Int64,1}` : number of receivers for every supersource 
"""
type Geom
	sx::Vector{Vector{Float64}}
	sz::Vector{Vector{Float64}}
	rx::Vector{Vector{Float64}}
	rz::Vector{Vector{Float64}}
	nss::Int64 # change to +ve only integer later
	ns::Vector{Int64} # change to +ve only integer later
	nr::Vector{Int64} # change to +ve only integer later
	"adding conditions that are to be false while construction"
	Geom(sx, sz, rx, rz, nss, ns, nr) = 
		any([
                  broadcast(length, sx) != ns,
                  broadcast(length, sz) != ns,
                  broadcast(length, rx) != nr,
                  broadcast(length, rz) != nr,
       		  length(sx) != nss,
       		  length(sz) != nss,
       		  length(rz) != nss,
       		  length(rx) != nss 
		  ]) ? 
		error("Geom construct") : new(sx, sz, rx, rz, nss, ns, nr)
end # type

"""
Return some derived fields of `Geom`

# Outputs
`attrib=:uspos` : a tuple of x and z positions of unique sources
`attrib=:urpos` : a tuple of x and z position of unique receivers
`attrib=:nus` : number of unique source positions
`attrib=:nur` : number of unique source positions
"""
function Geom_get(acq::Vector{Geom}, attrib::Symbol)
	ngeom = length(acq)
	sposx = Geom_getvec(acq,:sx);	sposz = Geom_getvec(acq,:sz)
	isequal(length(sposx), length(sposz)) ? 
		spos = [[sposz[is],sposx[is]] for is in 1:length(sposx)] : error("input acq corrupt")

	rposx = Geom_getvec(acq,:rx);	rposz = Geom_getvec(acq,:rz)
	isequal(length(rposx), length(rposz)) ? 
		rpos = [[rposz[ir],rposx[ir]] for ir in 1:length(rposx)] : error("input acq corrupt")

	uspos = unique(spos); nus=length(uspos)
	urpos = unique(rpos); nur=length(urpos)
	uspost = ([uspos[iu][1] for iu in 1:nus], [uspos[iu][2] for iu in 1:nus]);
        urpost = ([urpos[iu][1] for iu in 1:nur], [urpos[iu][2] for iu in 1:nur]);


	if(attrib == :nus)
		return nus
	elseif(attrib == :nur)
		return nur
	elseif(attrib == :uspos)
		return uspost
	elseif(attrib == :urpos)
		return urpost
	elseif(attrib == :geomurpos)
		return [Geom(acq[igeom].sx, acq[igeom].sz, 
	       fill(urpost[2],acq[igeom].nss), 
	       fill(urpost[1],acq[igeom].nss), acq[igeom].nss, acq[igeom].ns, 
	       fill(nur,acq[igeom].nss)) for igeom=1:ngeom]
	elseif(attrib == :geomuspos)
		return [Geom(fill(uspost[2],acq[igeom].nss), 
	       fill(uspost[1],acq[igeom].nss),
	       acq[igeom].rx, acq[igeom].rz, 
	       acq[igeom].nss, fill(nus,acq[igeom].nss),
	       acq[igeom].nr) for igeom=1:ngeom]
	else
		error("invalid attrib")
	end
end

"""
Given receiver positions `rpos` and `rpos0`.
Returns an array Int indices of the dimension of number of supersources
with `true` at indices, if the waves due to that particular source are 
recorded.
"""
function Geom_find(acq::Geom; rpos::Array{Float64,1}=nothing, rpos0::Array{Float64,1}=nothing)
	rpos==nothing ? error("need rpos") : nothing
	sson = Array(Vector{Int64}, acq.nss);
	for iss = 1:acq.nss
		rvec = [[acq.rz[iss][ir],acq.rx[iss][ir]] for ir=1:acq.nr[iss]]
		ir=findfirst(rvec, rpos)
		if(rpos0==nothing)
			ir0=1;
		else
			ir0=findfirst(rvec, rpos0);
		end
		sson[iss] = ((ir != 0) && (ir0 != 0)) ? [ir,ir0] : [0]
	end
	return sson
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
function Geom_boundary(acqgeom::Geom,
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
		return Geom(sx, sz, fill(bx,nss), fill(bz,nss), nss, ns, fill(nr,nss))
	"change the position of source (not supersources) to the boundary for back propagation"
	elseif(attrib == :srcborder)
		ns = nb;		rx = acqgeom.rx;
		rz = acqgeom.rz;		nr = acqgeom.nr;
		nss = acqgeom.nss;
		return Geom(fill(bx,nss), fill(bz,nss), rx, rz, nss, fill(ns,nss), nr)
	else
		error("invalid attrib")
	end
end


"""
Return fixed spread acquisition geometry depending 
on either horizontal or vertical array
It has only one source for every supersource

* `rand_flags::Vector{Bool}=[false, false]` : randomly or equally spaced supersources and receivers?
"""
function Geom_fixed(
	      smin::Float64,
	      smax::Float64,
	      s0::Float64,
	      rmin::Float64,
	      rmax::Float64,
	      r0::Float64,
	      nss::Int64,
	      nr::Int64,
	      sattrib::Symbol=:horizontal,
	      rattrib::Symbol=:horizontal,
 	      rand_flags::Vector{Bool}=[false, false]
	     )


	sarray = isequal(nss,1) ? fill(smin,1) : (rand_flags[1] ? 
					   rand(Uniform(minimum([smin,smax]),maximum([smin,smax])),nss) : linspace(smin,smax,nss))
	if(sattrib==:horizontal)
		sz = fill(s0,nss); sx=sarray
	elseif(sattrib==:vertical)
		sx = fill(s0,nss); sz=sarray
	else
		error("invalid sattrib")
	end

	rarray = isequal(nr,1) ? fill(rmin,1) : (rand_flags[2] ? 
				        rand(Uniform(minimum([rmin,rmax]),maximum([rmin,rmax])),nr) : linspace(rmin,rmax,nr))
	if(rattrib==:horizontal)
		rz = fill(r0,nr); rx = rarray
	elseif(rattrib==:vertical)
		rx = fill(r0,nr); rz = rarray
	else
		error("invalid rattrib")
	end
	rxall = [rx for iss=1:nss];
	rzall = [rz for iss=1:nss];
	sxall = [[sx[iss]] for iss=1:nss];
	szall = [[sz[iss]] for iss=1:nss];
	return Geom(sxall, szall, rxall, rzall, nss, fill(1,nss), fill(nr,nss))
end

"""
Circular acquisition. The sources and receivers can also be placed on a circle of radius
`rad`. The origin of the circle is at `loc`. 
This geometry is unrealistic, but useful for testing.
Receivers are placed such that the limits 
of the angular offset are given by `θlim`

# Arguments
`nss::Int64=10` : number of supersources
`nr::Int64=10` : number receivers for each super source
`loc::Vector{Float64}=[0.,0.]` : location of origin
`rad::Float64=100.` : radius
`θlim::Vector{Float64}=[0.,2*pi]` : range of angular offset between source and receiver

"""
function Geom_circ(;
		   nss::Int64=10,
		   nr::Int64=10,
		   loc::Vector{Float64}=[0.,0.],
		   rad::Float64=100.,
		   θlim::Vector{Float64}=[0.,2*pi]
	       )
	rxall = fill(zeros(nr),nss)
	rzall = fill(zeros(nr),nss)
	sx = zeros(nss)
	sz = zeros(nss)
	δθr = abs((θlim[2] - θlim[1])) / nr;
	for iss =1:nss
		θs = (iss)*2.*pi / (nss);
		sx[iss] = rad*cos(θs) + loc[2]
		sz[iss] = rad*sin(θs) + loc[1]
		rxall[iss] = [rad * cos(θs + θlim[1] + (ir-1)*δθr) + loc[2] for ir in 1:nr]
		rzall[iss] = [rad * sin(θs + θlim[1] + (ir-1)*δθr) + loc[1] for ir in 1:nr]
	end
	sxall = [[sx[iss]] for iss=1:nss];
	szall = [[sz[iss]] for iss=1:nss];
	return Geom(sxall, szall, rxall, rzall, nss, fill(1,nss), fill(nr,nss))
end

"""
Appends the input a vector of acquisition geometries.
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
Advance either source or receiver array in an acquisition geometry in horizontal or vertical directions.

# Arguments
* `geom::Geom` : acquisition geometry that is updated
* `advances::Vector{Float64}=[[0.,0.], [0.,0.,]]` : source and receiver advancements
"""
function Geom_advance(geom::Geom, advances::Vector{Vector{Float64}}=[[0.,0.], [0.,0.]]) 
	geom_out=deepcopy(geom)
	for iss=1:geom.nss
		geom_out.sx[iss][:] += advances[1][2]
		geom_out.sz[iss][:] += advances[1][1]
		geom_out.rx[iss][:] += advances[2][2]
		geom_out.rz[iss][:] += advances[2][1]
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
	ns::Vector{Int64}
	nfield::Int64
	wav::Array{Array{Float64,2},2}
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
function Src_fixed(nss::Int64, 
	     ns::Int64, 
	     nfield::Int64,
	     wav::Array{Float64},
	     tgrid::Grid.M1D
	     )

	wavsrc = [repeat(wav,inner=(1,ns)) for iss=1:nss, ifield=1:nfield] 
	return Src(nss, fill(ns, nss), nfield, wavsrc, tgrid)
end

"""
Constructor of `Src`, which is typical for a input model such that 
the model has `nwav` wavelengths.

* `nwav::Int64=10` : number of wavelenghts in the model
nwav
"""
function Src_fixed_mod(nss::Int64, ns::Int64, nfield::Int64, model::Models.Seismic, nwav::Int64=10)

	x=model.mgrid.x; z=model.mgrid.z
	# dominant wavelength using model dimensions
	λdom=mean([(abs(x[end]-x[1])), (abs(z[end]-z[1]))])/real(nwav)
	# average P velocity
	vavg=Models.χ([mean(model.χvp)], model.vp0, -1)[1]

	fqdom = vavg/λdom
	println("\t choosing dominant frequency of Ricker as:\t",fqdom)

	# maximum distance the wave travels
	d = sqrt((x[1]-x[end]).^2+(z[1]-z[end]).^2)
	
	# use two-way maximum distance to get tmax
	tmax=2.*d/vavg
	δt=0.1*minimum([model.mgrid.δx, model.mgrid.δz])/Models.χ([(maximum(model.χvp))],model.vp0,-1)[1]
	# choose sampling interval to obey max freq of source wavelet
	tgrid=Grid.M1D(0.0, tmax, δt)
	wav = Wavelets.ricker(fqdom=fqdom, tgrid=tgrid);
	println("\t maximum modelling time is:\t",tmax," sec (",tgrid.nx," samples)")

	# choose ricker waveletes of fdom
	return Src_fixed(nss, ns, nfield, wav, tgrid)
end

"""
Function that returns Src after time reversal
"""
function Src_tr(src::Src)
	return Src(src.nss,src.ns,src.nfield,
	    [flipdim(src.wav[i,j],1) for i in 1:src.nss, j in 1:src.nfield],src.tgrid)
end



"""
Pad `Src` 
tgrids should be same in all Src
"""
function Src_uspos(src::Vector{Src}, acq::Vector{Geom})
	np = length(src) == length(acq) ? length(src) : error("unequal sizez")

	# unique source positions
	nus=Geom_get(acq,:nus) 
	uspos=Geom_get(acq,:uspos)
	# all zeros for all unique positions
	wavout = [[zeros(src[ip].tgrid.nx,nus) for iss=1:src[ip].nss, ifield=1:src[ip].nfield] for ip=1:np]
	# fill source wavelets when necessary
	for ip=1:np
		for ifield=1:src[ip].nfield, iss=1:src[ip].nss, is=1:acq[ip].ns[iss]
			is0=find([[uspos[1][i]-acq[ip].sz[iss][is],
		      uspos[2][i]-acq[ip].sx[iss][is]] == [0., 0.,] for i in 1:nus])

			wavout[ip][iss, ifield][:,is0] = src[ip].wav[iss, ifield][:,is] 
		end
	end
	# output src
	return [Src(src[ip].nss,fill(nus, src[ip].nss),
	     	src[ip].nfield,wavout[ip],src[ip].tgrid) for ip=1:np]
end


"""
return a vector of the order 

"""
function Src_getvec(src::Vector{Src}, field::Symbol)
	np = length(src);
	vect = [getfield(src[ip],field) for ip=1:np]
	return vec(hcat(hcat(vect...)...))
end

end # module
