
"""
Gallery of acquisition ageometries `AGeom` using an input mesh.
The sources and receivers are not placed anywhere on the edges of the mesh.

# Arguments 
* `mgrid` : a 2-D mesh
* `attrib::Symbol` : attribute decides output
  * `=:xwell` cross-well acquisition
  * `=:surf` surface acquisition
  * `=:vsp` vertical seismic profiling
  * `=:rvsp`  reverse vertical seismic profiling
  * `=:downhole` downhole sources and receivers 

# Keyword Arguments

* `nss=2` : number of supersources
* `nr=2` : number of receivers per supersource
* `rand_flags::Vector{Bool}=[false, false]` : randomly or equally spaced supersources and receivers.
"""
function Array{AGeomss,1}(mgrid::AbstractArray{T}, attrib::Symbol, ss::SSrcs=SSrcs(1), r::Recs=Recs(10)) where {T<:StepRangeLen}
	otx=(0.9*mgrid[2][1]+0.1*mgrid[2][end]); ntx=(0.1*mgrid[2][1]+0.9*mgrid[2][end]);
	otwx=(0.95*mgrid[2][1]+0.05*mgrid[2][end]); ntwx=(0.05*mgrid[2][1]+0.95*mgrid[2][end]);
	otz=(0.9*mgrid[1][1]+0.1*mgrid[1][end]); ntz=(0.1*mgrid[1][1]+0.9*mgrid[1][end]);
	otwz=(0.95*mgrid[1][1]+0.05*mgrid[1][end]); ntwz=(0.05*mgrid[1][1]+0.95*mgrid[1][end]);
	quatx = (0.75*mgrid[2][1]+0.25*mgrid[2][end]); quatz = (0.75*mgrid[1][1]+0.25*mgrid[1][end]) 
	tquatx = (0.25*mgrid[2][1]+0.75*mgrid[2][end]); tquatz = (0.25*mgrid[1][1]+0.75*mgrid[1][end]) 
	halfx = 0.5*(mgrid[2][1]+mgrid[2][end]);	halfz = 0.5*(mgrid[1][1]+mgrid[1][end]);
	ageom=AGeom(mgrid, ss, Srcs(1), r)
	if(attrib == :xwell)
		update!(ageom, SSrcs(), [otz,otx], [ntz,otx])
		update!(ageom, Recs(), [otz,ntx], [ntz,ntx])
	elseif(attrib == :surf)
		update!(ageom, SSrcs(), [otwz,otx], [otwz,ntx])
		update!(ageom, Recs(), [otwz,otx], [otwz,ntx])
	elseif(attrib == :vsp)
		update!(ageom, SSrcs(), [otwz,otx], [otwz,ntx])
		update!(ageom, Recs(), [quatz,otx], [ntz,ntx])
	elseif(attrib == :rvsp)
		update!(ageom, SSrcs(), [quatz,otx], [ntz,otx])
		update!(ageom, Recs(),  [otwz,otx], [otwzz,ntx])
	elseif(attrib == :downhole)
		update!(ageom, SSrc(), [quatz,quatx], [otz,quatx])
		update!(ageom, Recs(),  [quatz,quatx], [otz,quatx])
	else
		error("invalid attrib")
	end
	return ageom
end

