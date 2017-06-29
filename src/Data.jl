__precompile__()

module Data

import SIT.Acquisition
import SIT.Grid
import SIT.Coupling
import SIT.DSP
using Interpolations

"""
Time domain representation of Seismic Data.
TODO: Also include acqsrc?

# Fields
* `d::Array{Array{Float64,2},2}` : data 
* `nfield::Int64` : number of components at each receiver
* `tgrid::Grid.M1D` : grid to represent time
* `acqgeom::Acquisition.Geom` : geometry used to generate the data
"""
type TD
	d::Array{Array{Float64,2},2}
	nfield::Int64
	tgrid::Grid.M1D
	acqgeom::Acquisition.Geom
	"adding conditions that are to be false while construction"
	TD(d, nfield, tgrid, acqgeom) = 
		any([
       		  nfield < 0.0,
		  broadcast(size,d) != [(tgrid.nx,acqgeom.nr[iss]) for iss=1:acqgeom.nss, ifield=1:nfield]
		  ]) ? 
		error("TD construct") : new(d, nfield, tgrid, acqgeom)

end


"""
function to resample data in time domain

# Arguments
* `data` : input data of type `TD`
* `tgrid` : resampling in time according to this time grid
"""
function TD_resamp(data::TD,
		tgrid::Grid.M1D
		)
	nss = data.acqgeom.nss
	nr = data.acqgeom.nr
	dataout = TD(
	      [zeros(tgrid.nx,data.acqgeom.nr[iss]) for iss=1:nss, ifield=1:data.nfield],
	      data.nfield,tgrid,data.acqgeom)
	TD_resamp!(data, dataout)
	for ifield = 1:data.nfield, iss = 1:nss, ir = 1:nr[iss]
		itp = interpolate((data.tgrid.x,),
		    data.d[iss, ifield][:, ir], 
			     Gridded(Linear()))
		dataout.d[iss, ifield][:,ir] = itp[tgrid.x]
	end
	return dataout
end

function TD_resamp!(data::TD, dataout::TD)
	nss = data.acqgeom.nss
	nr = data.acqgeom.nr
	for ifield = 1:data.nfield, iss = 1:nss, ir = 1:nr[iss]
		itp = interpolate((data.tgrid.x,),
		    data.d[iss, ifield][:, ir], 
			     Gridded(Linear()))
		dataout.d[iss, ifield][:,ir] = copy(itp[dataout.tgrid.x])
	end
	return dataout
end


"""
Return zeros
"""
function TD_zeros(nfield::Int64,
		 tgrid::Grid.M1D,
		 acqgeom::Acquisition.Geom
		 )
	return TD([zeros(tgrid.nx,acqgeom.nr[iss]) for iss=1:acqgeom.nss, ifield=1:nfield],
    				nfield,tgrid,acqgeom)
end

"""
Return bool depending on if `d` in `TD` is all zero. 
"""
function TD_iszero(data::TD)
	return maximum(broadcast(maximum,broadcast(abs,data.d))) == 0.0 ? true : false
end

"""
Time reverse the records of each receiver in `TD` 
"""
function TD_tr!(data::TD)
	data.d = copy([flipdim(data.d[i,j],1) for i in 1:data.acqgeom.nss, j in 1:data.nfield]);
end

"""
Returns dot product of the data
"""
function TD_dot(data1::TD, data2::TD)
	dotd = 0.0;
	for ifield = 1:data1.nfield, iss = 1:data1.acqgeom.nss, ir = 1:data1.acqgeom.nr[iss]
		dotd += dot(data1.d[iss, ifield][:, ir],data2.d[iss, ifield][:, ir])
	end
	return dotd
end


"""
normalize time-domain seismic data
"""
function TD_normalize(data::TD, attrib::Symbol)
	nr = data.acqgeom.nr;	nss = data.acqgeom.nss;	nt = data.tgrid.nx;
	datan = deepcopy(data);
	for ifield = 1:data.nfield, iss = 1:nss, ir = 1:nr[iss]
		if(attrib == :recrms)
			nval = sqrt(mean(datan.d[iss, ifield][:,ir].^2.))
		elseif(attrib == :recmax)
			nval = maximum(datan.d[iss, ifield][:,ir])
		else
			error("invalid attrib")
		end

		# normalize
		datan.d[iss, ifield][:, ir] = 
		isequal(nval, 0.0) ? zeros(nt) : datan.d[iss, ifield][:, ir]./nval  
	end
	return datan
end

"""
Construct TD using data at all the unique receiver positions
for all supersources.
"""
function TD_urpos(d::Array{Float64}, 
		   nfield::Int64, 
		   tgrid::Grid.M1D, 
		   acq::Acquisition.Geom,
		   nur::Int64,
		   urpos::Tuple{Array{Float64,1},Array{Float64,1}
		  }
		   )
	dout = [zeros(tgrid.nx,acq.nr[iss]) for iss=1:acq.nss, ifield=1:nfield] 

	for ifield=1:nfield, iss=1:acq.nss, ir=1:acq.nr[iss]
		# find index in urpos
		irr=find([[urpos[1][i]-acq.rz[iss][ir],
		       urpos[2][i]-acq.rx[iss][ir]] == [0., 0.,] for i in 1:nur])

		dout[iss, ifield][:,ir] = d[:,irr[1],ifield, iss] 
	end

	return TD(dout, nfield, tgrid, acq)

end

"""
Apply coupling functions to TD

* :s means s is returned in the place of rs
* :w means 
"""
function TDcoup!(
               s::TD,
	       r::TD,
	       w::Coupling.TD,
	       attrib::Symbol
	       )
	nr = r.acqgeom.nr;	nss = r.acqgeom.nss;	nt = r.tgrid.nx;
	nfield = (w.nfield == r.nfield == s.nfield) ? w.nfield : error("different nfields")
	sv=zeros(s.tgrid.nx)
	rv=zeros(r.tgrid.nx)
	wv=zeros(w.tgridssf.nx)
	for ifield = 1:nfield, iss = 1:nss, ir = 1:nr[iss]
		# receiver coupling
	#	DSP.fast_filt!(s.d[iss, ifield][:, ir],r.d[iss, ifield][:, ir],
#		 w.rf[iss, ifield][:,ir],:s)
		# source coupling
		sv=s.d[iss, ifield][:, ir]
		rv=r.d[iss, ifield][:, ir]
		wv=w.ssf[iss, ifield]
		DSP.fast_filt!(sv,rv,wv,attrib)
		s.d[iss, ifield][:, ir]=copy(sv)
		r.d[iss, ifield][:, ir]=copy(rv)
		w.ssf[iss, ifield][:]=copy(wv)
	end
end # TDcoup


end # module
