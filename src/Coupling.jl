module Coupling


import SIT.Acquisition
import SIT.Grid

type TD
	ssf::Array{Array{Float64,1},2}
	rf::Array{Array{Float64,1},2}
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
		  acqgeom::Acquisition.Geom
		 )
	tgridssf = Grid.M1D(-1.0*tlagssf, tlagssf, δt)
	tgridrf = Grid.M1D(-1.0*tlagrf, tlagrf, δt)
	spiss = zeros(tgridssf.nx); spiss[findn(tgridssf.==0.0)] = 1.0; 
	spir = zeros(tgridrf.nx); spir[findn(tgridrf.==0.0)] = 1.0; 
	
	# number of unique receivers (implement in future)

	return TDcoup([spiss for iss=1:acqgeom.nss, ifield=1:nfield],
		  [repeat(spir, outer=[1,acqgeom.nr[iss]]) for iss=1:acqgeom.nss, ifield=1:nfield],
    				nfield,tgridsf,tgridrf,acqgeom)
end


end 
