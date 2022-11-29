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

import GeoPhyInv: AGeom

"""
Time-domain source and receiver filters.

# Fields

* `ssf::Array{Array{Float64,1},2}` : source filters for each supersource and recorded component
* `rf::Array{Array{Float64,2},2}` : receiver filters for each receiver, supersource and recorded component
* `fields::Vector{Symbol}` :  number of recorded components at receivers
* `tgridssf` : a  time grid for source filters with both positive and negative lags
* `tgridrf` : a  time grid for receiver filters with both positive and negative lags
* `ageom::AGeom` : acquisition ageometry
"""
mutable struct TD
	ssf::Array{Array{Float64,1},2}
	rf::Array{Array{Float64,2},2}
	fields::Vector{Symbol}
	tgridssf::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
	tgridrf::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
	ssflags::Vector{Int}
	rflags::Vector{Int}
	ageom::AGeom
	TD(ssf, rf, fields, tgridssf, tgridrf, ssflags, rflags, ageom) = 
		any([
		  any([fields[ifield] ∉ [:P, :vx, :vz] for ifield in 1:length(fields)]),
		  length(fields) == 0,
		  broadcast(size,ssf) != [(length(tgridssf),) for iss=1:ageom.nss, ifield=1:length(fields)]
		  ]) ? 
		error("error in TD construction") : new(ssf, rf, fields, tgridssf, tgridrf, ssflags, rflags, ageom)

end 

"""
Initialize coupling filters `TD` with  delta functions.

# Arguments

* `tgriddata` : time grid of the data
* `tlagssf_fracs::Vector{Real}` : +ve and -ve fractions of source filter
* `tlagrf_fracs::Vector{Real}` : +ve and -ve fractions for receiver filter
* `δt:Float64` : sampling interval in time
* `fields::Vector{Symbol}` : number of components
* `ageom::AGeom` : acquisition ageometry

# Return

* time-domain coupling filters as `TD`
"""
function TD_delta(
			   tgriddata, 
			   tlagssf_fracs,
			   tlagrf_fracs,
		  fields::Vector{Symbol},
		  ageom::AGeom
		 )
	δt=tgriddata.δx
	tot_t=abs(tgriddata[end]-tgriddata[1])

	error("fix these")
	#tgridssf, ssflags = lag(tot_t.*tlagssf_fracs, δt)
	#tgridrf, rflags =lag(tot_t.*tlagrf_fracs, δt)

	spiss = zeros(length(tgridssf)); spir = zeros(length(tgridrf));
	# check where to put spikes
	nspikessf=(findall(tgridssf.==0.0) == []) ? error("no t=0") : findall(tgridssf.==0.0)
	nspikerf=(findall(tgridrf.==0.0) == []) ? error("no t=0") : findall(tgridrf.==0.0)
	# put spikes
	spiss[nspikessf] = 1.0; spir[nspikerf] = 1.0; 
	# number of unique receivers (implement in future)
	# one receiver filter for each unique receiver

	ssf=[spiss for iss=1:ageom.nss, ifield=1:length(fields)]

	return TD(ssf,
	   [repeat(spir, outer=[1,ageom.nr[iss]]) for iss=1:ageom.nss, ifield=1:length(fields)],
	   fields,tgridssf,tgridrf,ssflags,rflags,deepcopy(ageom))

end


end # module
