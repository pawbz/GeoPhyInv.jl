module Data

import SIT.Acquisition
import SIT.Grid
using Interpolations

"""
time domain representation of Seismic Data

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
	ns = data.acqgeom.nss
	dataout = TD(zeros(tgrid.nx, nr, ns, data.nfield),data.nfield,tgrid,data.acqgeom)
	for ifield = 1:data.nfield
		for is = 1:ns
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

end # module
