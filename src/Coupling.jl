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


using Grid
import JuMIT.Acquisition

"""
Time-domain source and receiver filters.

# Fields

* `ssf::Array{Array{Float64,1},2}` : source filters for each supersource and recorded component
* `rf::Array{Array{Float64,2},2}` : receiver filters for each receiver, supersource and recorded component
* `fields::Vector{Symbol}` :  number of recorded components at receivers
* `tgridssf::Grid.M1D` : a  time grid for source filters with both positive and negative lags
* `tgridrf::Grid.M1D` : a  time grid for receiver filters with both positive and negative lags
* `acqgeom::Acquisition.Geom` : acquisition geometry
"""
mutable struct TD
	ssf::Array{Array{Float64,1},2}
	rf::Array{Array{Float64,2},2}
	fields::Vector{Symbol}
	tgridssf::Grid.M1D
	tgridrf::Grid.M1D
	ssflags::Vector{Int}
	rflags::Vector{Int}
	acqgeom::Acquisition.Geom
	TD(ssf, rf, fields, tgridssf, tgridrf, ssflags, rflags, acqgeom) = 
		any([
		  any([fields[ifield] ∉ [:P, :Vx, :Vz] for ifield in 1:length(fields)]),
		  length(fields) == 0,
		  broadcast(size,ssf) != [(tgridssf.nx,) for iss=1:acqgeom.nss, ifield=1:length(fields)]
		  ]) ? 
		error("error in TD construction") : new(ssf, rf, fields, tgridssf, tgridrf, ssflags, rflags, acqgeom)

end 

"""
Initialize coupling filters `TD` with  delta functions.

# Arguments

* `tgriddata` : time grid of the data
* `tlagssf_fracs::Vector{Real}` : +ve and -ve fractions of source filter
* `tlagrf_fracs::Vector{Real}` : +ve and -ve fractions for receiver filter
* `δt:Float64` : sampling interval in time
* `fields::Vector{Symbol}` : number of components
* `acqgeom::Acquisition.Geom` : acquisition geometry

# Return

* time-domain coupling filters as `TD`
"""
function TD_delta(
			   tgriddata, 
			   tlagssf_fracs,
			   tlagrf_fracs,
		  fields::Vector{Symbol},
		  acqgeom::Acquisition.Geom
		 )
	δt=tgriddata.δx
	tot_t=abs(tgriddata.x[end]-tgriddata.x[1])

	tgridssf, ssflags = Grid.M1D_lag(tot_t.*tlagssf_fracs, δt)
	tgridrf, rflags = Grid.M1D_lag(tot_t.*tlagrf_fracs, δt)

	spiss = zeros(tgridssf.nx); spir = zeros(tgridrf.nx);
	# check where to put spikes
	nspikessf=(findn(tgridssf.x.==0.0) == []) ? error("no t=0") : findn(tgridssf.x.==0.0)
	nspikerf=(findn(tgridrf.x.==0.0) == []) ? error("no t=0") : findn(tgridrf.x.==0.0)
	# put spikes
	spiss[nspikessf] = 1.0; spir[nspikerf] = 1.0; 
	# number of unique receivers (implement in future)
	# one receiver filter for each unique receiver

	ssf=[spiss for iss=1:acqgeom.nss, ifield=1:length(fields)]

	return TD(ssf,
	   [repeat(spir, outer=[1,acqgeom.nr[iss]]) for iss=1:acqgeom.nss, ifield=1:length(fields)],
	   fields,tgridssf,tgridrf,ssflags,rflags,deepcopy(acqgeom))

end


end # module
