module Data

import SIT.Acquisition
import SIT.Grid
using Interpolations

"""
time domain representation of Seismic Data
TODO: Also include acqsrc?

# Fields
* `d` : data first sorted in time, then in receivers, and finally in sources
* `nfield` : number of components
* `tgrid` : `M1D` grid to represent time
* `acqgeom` : acquisition geometry
"""
type TD
	d::Array{Float64}
	nfield::Int64
	tgrid::Grid.M1D
	acqgeom::Acquisition.Geom
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
	nr = maximum(data.acqgeom.nr)
	nss = data.acqgeom.nss
	dataout = TD(zeros(tgrid.nx, nr, nss, data.nfield),data.nfield,tgrid,data.acqgeom)
	for ifield = 1:data.nfield
		for is = 1:nss
			for ig = 1:nr
				itp = interpolate((data.tgrid.x,),
					     data.d[:, ig, is, ifield], 
					     Gridded(Linear()))
				dataout.d[:,ig,is,ifield] = itp[tgrid.x]
			end
		end
	end
	return dataout
end

"""
normalize time-domain seismic data
"""
function TD_normalize(data::TD, attrib::Symbol)
	nr = maximum(data.acqgeom.nr);
	nss = data.acqgeom.nss;
	nt = data.tgrid.nx;
	datan = data;
	for ifield = 1:data.nfield, is = 1:nss, ig = 1:nr
		if(attrib == :recrms)
			nval = sqrt(mean(datan.d[:, ig, is, ifield].^2.))
		elseif(attrib == :recmax)
			nval = maximum(datan.d[:, ig, is, ifield])
		else
			error("invalid attrib")
		end

		# normalize
		datan.d[:, ig, is, ifield] = 
			isequal(nval, 0.0) ? zeros(nt) : datan.d[:, ig, is, ifield]./nval  
	end
	return datan
end
end # module
