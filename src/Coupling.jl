__precompile__()

"""
The data recorded at the receivers during a seismic experiment 
is a convolution of the Green's function,
the source signature and the receiver instrument
response.
The source filter `ssf` takes into account 
the effective wavelet present in the data. 
Similarly, the receiver filter `rf` takes the
receiver coupling into account.
The Inversion module can
can estimate both `ssf` and `rf` from the seismic data in order to achieve better data fitting.
A base source wavelet `Src` is always used to model the data 
before applying the time-domain filter data type  described in this module:
* `TD` : time domain filters for both source and receivers
"""
module Coupling


import JuMIT.Acquisition
import JuMIT.Grid

"""
Time-domain source and receiver filters.

# Fields

* `ssf::Array{Array{Float64,1},2}` : source filters for each supersource and recorded component
* `rf::Array{Array{Float64,2},2}` : receiver filters for each receiver, supersource and recorded component
* `nfield::Int64` :  number of recorded components at receivers
* `tgridssf::Grid.M1D` : a  time grid for source filters with both positive and negative lags
* `tgridrf::Grid.M1D` : a  time grid for receiver filters with both positive and negative lags
* `acqgeom::Acquisition.Geom` : acquisition geometry
"""
type TD
	ssf::Array{Array{Float64,1},2}
	rf::Array{Array{Float64,2},2}
	nfield::Int64
	tgridssf::Grid.M1D
	tgridrf::Grid.M1D
	acqgeom::Acquisition.Geom
end 

"""
Initialize coupling filters `TD` with  delta functions.

# Arguments
 
* `tlagssf<:Real` : maximum lag in seconds for source filter
* `tlagrf<:Real` : maximum lag in seconds  for receiver filter
* `δt:Float64` : sampling interval in time
* `nfield::Int64` : number of components
* `acqgeom::Acquisition.Geom` : acquisition geometry

# Return

* time-domain coupling filters as `TD`
"""
function TD_delta{T<:Real}(tlagssf::T,
		  tlagrf::T,
		  δt::Float64,
		  nfield::Int64,
		  acqgeom::Acquisition.Geom
		 )
	tgridssf = Grid.M1D_lag(tlagssf, δt)
	tgridrf = Grid.M1D_lag(tlagrf, δt)
	return TD_delta(tgridssf, tgridrf, nfield, acqgeom)
end

function TD_delta(tgridssf::Grid.M1D, tgridrf::Grid.M1D,
		  nfield::Int64, acqgeom::Acquisition.Geom )

	spiss = zeros(tgridssf.nx); spir = zeros(tgridrf.nx);
	# check where to put spikes
	nspikessf=(findn(tgridssf.x.==0.0) == []) ? error("no t=0") : findn(tgridssf.x.==0.0)
	nspikerf=(findn(tgridrf.x.==0.0) == []) ? error("no t=0") : findn(tgridrf.x.==0.0)
	# put spikes
	spiss[nspikessf] = 1.0; spir[nspikerf] = 1.0; 
	# number of unique receivers (implement in future)
	# one receiver filter for each unique receiver
	return TD([spiss for iss=1:acqgeom.nss, ifield=1:nfield],
		  [repeat(spir, outer=[1,acqgeom.nr[iss]]) for iss=1:acqgeom.nss, ifield=1:nfield],
    				nfield,tgridssf,tgridrf,acqgeom)
end

end # module
