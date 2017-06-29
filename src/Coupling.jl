__precompile__()

module Coupling


import SIT.Acquisition
import SIT.Grid

"""
The data recorded at the receiver during a seismic experiment 
is a convolution of the Green's function,
the source signature and the receiver instrument
response as explained in the Section~\ref{sec:sdata}.
The source filter, \srcf, takes into account 
the effective wavelet present in the data. 
Similarly, the receiver filter, \recf, takes the
receiver coupling into account.
\xfwi can estimate both \srcf and \recf during the inversion.
The base source wavelet, \sfo, is always used to model the data 
before applying the convolutional filters. 
\xfwi considers the source filter 
only if the input to \kword{WTD} contains 
\kwordi{[SRCF]}, as explained in Sub-section~\ref{subsec:tasks}.
Similarly, for the receiver filter,
\kword{WTD} should contain \kwordi{[RECF]}.
subsection{Keywords}
The keywords to input paramters related 
to the source and receiver filters contain \kword{[SF]} and are explained
below.
If the input is \kwordi{0}, then the filter is a spike at \tim=0,
with unknown amplitude.  
Which means, the filter only applies a scalar to the modelled data.
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
Return delta functions
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
