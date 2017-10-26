__precompile__()

"""
This module defines the following data types:
* `Geom` : acquisition geometry, i.e., positions of supersources, sources and receivers
* `Src` : source related acquisition parameters, e.g., source wavelet
It also provides methods that either does operations on these data type or 
help their construction.
"""
module Acquisition

import JuMIT.Grid
import JuMIT.Models
import JuMIT.Wavelets
import JuMIT.DSP
using Distributions

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

"Compare if two `TD`'s  are equal"
function Base.isequal(geom1::Geom, geom2::Geom)
	fnames = fieldnames(Geom)
	vec=[(isequal(getfield(geom1, name),getfield(geom2, name))) for name in fnames]
	return all(vec)
end

"""
Print information about `Geom`
"""
function Base.print(geom::Geom, name::String="")
	println("\tAcquisition Geometry:\t",name)
	println("\t> number of supersources:\t",geom.nss)
	println("\t> sources per supersource:\t","min\t",minimum(geom.ns[:]), "\tmax\t", maximum(geom.ns[:]))
	println("\t> receivers per supersource:\t","min\t",minimum(geom.nr[:]), "\tmax\t", maximum(geom.nr[:]))
	println("\t> number of unique positions of:\t","sources\t",Geom_get([geom],:nus), "\treceivers\t", Geom_get([geom],:nur))
end

"""
Return some derived fields of `Geom`

# Arguments 

* `acq::Vector{Geom}` : a vector of `Geom`
* `attrib::Symbol` : attribute to determine the return object 
  * `=:nus` number of unique source positions in acquisition
  * `=:nur` number of unique receiver positions in acquisition
  * `=:uspos` a tuple of x and z positions of all the unique sources
  * `=:urpos` a tuple of x and z position of all the unique receivers
  * `=:geomurpos` a `Geom` vector as if all the unique receiver positions are used for each supersource
  * `=:geomuspos` a `Geom` vector as if all the unique source positions are used for each supersource
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
	sson = Array{Vector{Int64}}(acq.nss);
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
Modify input `Geom` such that the output `Geom` has
either sources or receivers on the boundary of 
`mgrid`.

# Arguments

* `acqgeom::Geom` : input geometry
* `mgrid::Grid.M2D` : grid to determine the boundary
* `attrib::Symbol` : decide return
  * `=:srcborder` sources on boundary (useful for back propagation)
  * `=:recborder` receivers on boundary
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
Check if all the sources and receivers in `Geom` are within the model 

# Return

* `true` if all the positions within the model, `false` otherwise
"""
function Geom_check(geom::Geom, mgrid::Grid.M2D)
	xmin, zmin, xmax, zmax = mgrid.x[1], mgrid.z[1], mgrid.x[end], mgrid.z[end]


	checkvec=fill(false, geom.nss)
	for iss=1:geom.nss
		checkvec[iss] = any(
		      vcat(
		     (((geom.sx[iss]-xmin).*(xmax-geom.sx[iss])) .< 0.0),
		     (((geom.sz[iss]-zmin).*(zmax-geom.sz[iss])) .< 0.0),
		     (((geom.rx[iss]-xmin).*(xmax-geom.rx[iss])) .< 0.0),
		     (((geom.rz[iss]-zmin).*(zmax-geom.rz[iss])) .< 0.0),
		     ) )
	end
	return !(any(checkvec))
end

"""
A fixed spread acquisition has same set of sources and 
receivers for each supersource.
This method constructs a 
fixed spread acquisition geometry using either a
horizontal or vertical array of supersources/ receivers.
Current implementation has only one source for every supersource.

# Arguments 

* `smin::Float64` : minimum coordinate for sources
* `smax::Float64` : maximum coordinate for sources
* `s0::Float64` : consant coordinate for sources
* `rmin::Float64` : minimum coordinate for receivers
* `rmax::Float64` : maximum coordinate for receivers
* `r0::Float64` : consant coordinate for receivers
* `nss::Int64` : number of supersources
* `nr::Int64` : number of receivers
* `sattrib::Symbol=:horizontal` : supersource array kind
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

* `nss::Int64=10` : number of supersources
* `nr::Int64=10` : number receivers for each super source
* `loc::Vector{Float64}=[0.,0.]` : location of origin
* `rad::Vector{Float64}=[100.,100.]` : radius for source and receiver circles, for example,
  * `=[0.,100.]` for sources at the center of circle
* `θlim::Vector{Float64}=[0.,π]` : acquisition is limited to these angular offsets between 0 and π

# Return

* a circular acquisition geometry `Geom`
"""
function Geom_circ(;
		   nss::Int64=10,
		   nr::Int64=10,
		   loc::Vector{Float64}=[0.,0.],
		   rad::Vector{Float64}=[100.,100.],
		   θlim::Vector{Float64}=[0.,π]
	       )

	# modify nr such that approximately nr receivers are present in θlim
	nra = nr * convert(Int64, floor(π/(abs(diff(θlim)[1]))))

	(nra==0) && error("θlim should be from 0 to π")

	sxx, szz = circ_coord(loc..., nss, rad[1])
	rxxa, rzza = circ_coord(loc..., nra, rad[2])
	sx=Array{Array{Float64}}(nss)
	sz=Array{Array{Float64}}(nss)
	rx=Array{Array{Float64}}(nss)
	rz=Array{Array{Float64}}(nss)
	ns=Array{Int64}(nss)
	nr=Array{Int64}(nss)

	for iss=1:nss
		rangles = [angle(complex(rxxa[ira]-loc[2], rzza[ira]-loc[1])) for ira in 1:nra]
		sangles = angle(complex(sxx[iss]-loc[2], szz[iss]-loc[1]))
		diff = abs.(rangles-sangles)
		diff = [min(diff[ira], (2π - diff[ira])) for ira in 1:nra]
		angles = find(minimum(θlim) .<= diff .<= maximum(θlim))
		sx[iss] = [sxx[iss]]
		sz[iss] = [szz[iss]]
		ns[iss] = 1 
		nr[iss] = length(angles)
		rx[iss] = rxxa[angles]
		rz[iss] = rzza[angles]
	end
	geom =  Geom(sx, sz, rx, rz, nss, ns, nr)

	print(geom, "circular acquisition geometry")
	return geom
end

function circ_coord(z, x, n, rad)
	θ = (n==1) ? [0.0] : collect(linspace(0, 2π, n+1))
	xx = (rad .* cos.(θ) + x)[1:n]
	zz = (rad .* sin.(θ) + z)[1:n]
	return xx, zz
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
Data type for the source related parameters during acquisiton.

# Fields
* `nss::Int64` : number of supersources
* `ns::Array{Int64}` : number of sources for each supersource
* `fields::Vector{Symbol}` : number of fields
* `wav::Array{Float64}` : wavelets in time domain
* `tgrid::Grid.M1D` : time grid 
"""
type Src
	nss::Int64
	ns::Vector{Int64}
	fields::Vector{Symbol}
	wav::Array{Array{Float64,2},2}
	tgrid::Grid.M1D
	Src(nss, ns, fields, wav, tgrid) = 
		any([
		  any([fields[ifield] ∉ [:P, :Vx, :Vz] for ifield in 1:length(fields)]),
		  length(fields) == 0,
		  broadcast(size,wav) ≠ [(tgrid.nx,ns[iss]) for iss=1:nss, ifield=1:length(fields)]
		  ]) ? 
		error("error in Src construction") : new(nss, ns, fields, wav, tgrid)
end

function Base.isapprox(src1::Src, src2::Src)
	vec=([(src1.nss==src2.nss), (src1.fields==src2.fields), (src1.ns==src2.ns), 
       		(isequal(src1.tgrid, src2.tgrid)),
		])
	return all(vec)
end

function Base.copy!(srco::Src, src::Src)
	if(isapprox(srco, src))
		for f in [:wav]
			srcowav=getfield(srco, f)
			srcwav=getfield(src, f)
			for iss=1:src.nss, ifield=1:length(src.fields) 
				srcowavv=srcowav[iss,ifield]
				srcwavv=srcwav[iss,ifield]
				copy!(srcowavv, srcwavv)
			end
		end
		return srco
	else
		error("attempt to copy dissimilar sources")
	end
end


"""
Allocate `Src` with zeros depending on the acquisition geometry.
"""
function Src_zeros(acqgeom::Geom,  fields::Vector{Symbol}, tgrid::Grid.M1D)
	wavsrc = [zeros(tgrid.nx,acqgeom.ns[iss]) for iss=1:acqgeom.nss, ifield=1:length(fields)] 
	return Src(acqgeom.nss, acqgeom.ns, fields, wavsrc, deepcopy(tgrid))
end

"""
Priiiint information about `Src`
"""
function Base.print(src::Src, name::String="")
	println("\tSource Acquisition:\t",name)
	println("\t> number of supersources:\t",src.nss)
	println("\t> sources per supersource:\t","min\t",minimum(src.ns[:]), "\tmax\t", maximum(src.ns[:]))
	
	freqmin, freqmax, freqpeak = freqs(src)
	println("\t> frequency:\t","min\t",freqmin, "\tmax\t",freqmax,"\tpeak\t",freqpeak)
	println("\t> time:\t","min\t",src.tgrid.x[1], "\tmax\t", src.tgrid.x[end])
	println("\t> samples:\t",src.tgrid.nx)
end


"""
Return minimum, maximum and peak frequencies of `Src`
"""
function freqs(src::Src)
	freqmin = minimum([DSP.findfreq(src.wav[i,j][:,:],src.tgrid,attrib=:min) for i in 1:src.nss, j in 1:length(src.fields)])
	freqmax = maximum([DSP.findfreq(src.wav[i,j][:,:],src.tgrid,attrib=:max) for i in 1:src.nss, j in 1:length(src.fields)])
	freqpeak = mean([DSP.findfreq(src.wav[i,j][:,:],src.tgrid,attrib=:peak) for i in 1:src.nss, j in 1:length(src.fields)])
	return freqmin, freqmax, freqpeak
end

"""
Constructor for `Src` data type.
Uses same source wavelet, i.e., `wav` for all sources and supersources

# Arguments

* `nss::Int64` : number of supersources
* `ns::Int64` : number of sources
* `fields::Vector{Symbol}` : number of fields the sources are exciting
* `wav::Array{Float64}` : a source wavelet that is used for all sources and supersources
* `tgrid::Grid.M1D` : time grid for the wavelet
"""
function Src_fixed(nss::Int64, 
	     ns::Int64, 
	     fields::Vector{Symbol},
	     wav::Array{Float64},
	     tgrid::Grid.M1D
	     )

	wavsrc = [repeat(wav,inner=(1,ns)) for iss=1:nss, ifield=1:length(fields)] 
	return Src(nss, fill(ns, nss), fields, wavsrc, deepcopy(tgrid))
end

"""
Constructor of `Src`, which is typical for a input model such that 
the model has `nλ` wavelengths.

# Arguments

* `nss::Int64` : number of supersources
* `ns::Int64` : number of source per each supersource
* `fields::Vector{Symbol}` :

# Keyword Arguments

* `mod::Models.Seismic` :
* `nλ::Int64=10` : number of wavelengths in the mod
* `wav_func::Function=(fqdom, tgrid)->Wavelets.ricker(fqdom,tgrid)` : which wavelet to generate, see Wavelets.jl
* `tmaxfrac::Float64=1.0` : by default the maximum modelling time is computed using the average velocity and the diagonal distance of the model, use this fraction to increase of reduce the maximum time
"""
function Src_fixed_mod(
		nss::Int64, ns::Int64, fields::Vector{Symbol}; 
		mod::Models.Seismic=nothing, 
		nλ::Int64=10,
		wav_func::Function=(fqdom, tgrid)->Wavelets.ricker(fqdom,tgrid),
		tmaxfrac::Float64=1.0
		)

	x=mod.mgrid.x; z=mod.mgrid.z
	# dominant wavelength using mod dimensions
	λdom=mean([(abs(x[end]-x[1])), (abs(z[end]-z[1]))])/real(nλ)
	# average P velocity
	vavg=Models.χ([mean(mod.χvp)], mod.vp0, -1)[1]

	fqdom = vavg/λdom

	# maximum distance the wave travels
	d = sqrt((x[1]-x[end]).^2+(z[1]-z[end]).^2)
	
	# use two-way maximum distance to get tmax
	tmax=2.*d/vavg .* tmaxfrac

	# choose sampling interval to obey max freq of source wavelet
	δmin = minimum([mod.mgrid.δx, mod.mgrid.δz])
	vpmax = mod.vp0[2]
	δt=0.5*δmin/vpmax

	# check if δt is reasonable
	#(δt > 0.1/fqdom) : error("decrease spatial sampling or nλ")

	tgrid=Grid.M1D(0.0, tmax, δt)
	wav = wav_func(fqdom, tgrid);

	src=Src_fixed(nss, ns, fields, wav, tgrid)
	print(src)
	# choose ricker waveletes of fdom
	return src
end

"""
Generate band-limited random source signals 
"""
function Src_fixed_random(nss::Int64, ns::Int64, fields::Vector{Symbol}, fmin::Float64, fmax::Float64, tgrid::Grid.M1D, tmax::Float64=tgrid.x[end] )
	wavsrc = [repeat(zeros(tgrid.nx),inner=(1,ns)) for iss=1:nss, ifield=1:length(fields)] 
	for ifield in 1:length(fields), iss in 1:nss, is in 1:ns
		wavsrc[iss, ifield][:,is] = DSP.get_tapered_random_tmax_signal(tgrid, fmin, fmax, tmax)
	end
	src= Src(nss, fill(ns, nss), fields, wavsrc, deepcopy(tgrid))
	print(src)
	return src
end

"""
Function that returns Src after time reversal
"""
function Src_tr(src::Src)
	return Src(src.nss,src.ns,src.fields,
	    [flipdim(src.wav[i,j],1) for i in 1:src.nss, j in 1:length(src.fields)],deepcopy(src.tgrid))
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
	wavout = [[zeros(src[ip].tgrid.nx,nus) for iss=1:src[ip].nss, ifield=1:length(src[ip].fields)] for ip=1:np]
	# fill source wavelets when necessary
	for ip=1:np
		for ifield=1:length(src[ip].fields), iss=1:src[ip].nss, is=1:acq[ip].ns[iss]
			is0=find([[uspos[1][i]-acq[ip].sz[iss][is],
		      uspos[2][i]-acq[ip].sx[iss][is]] == [0., 0.,] for i in 1:nus])

			wavout[ip][iss, ifield][:,is0] = src[ip].wav[iss, ifield][:,is] 
		end
	end
	# output src
	return [Src(src[ip].nss,fill(nus, src[ip].nss),
	     src[ip].fields,wavout[ip],src[ip].tgrid) for ip=1:np]
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
