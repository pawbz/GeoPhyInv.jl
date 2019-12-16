
"""
```julia
AGeom=Vector{AGeomss,1}
```
An acquisition geometry, bundled into the mutable `AGeom` type,
has to be specified in order to
create any `Expt` variable.
Acquisition has supersources, sources and receivers.
For `AGeom` is a mutable type to define source-receiver geometry. 
Each supersource has `ns` sources that are
injected (or active) simultaneously, and 
a set of `nr` receivers to 
record. 
`AGeomss` is a subtype for each supersource, therefore:
This package provides tools to easily create commonly-used `AGeom` instances, 
after deciding on a spatial grid `mgrid` for the experiment.

## Indexing and Fields
Let `ageom` be an instance of this type, then fields 
can be accessed using:

* `ageom[i]` : acquisition for ith supersource
* `ageom[i].sx` : x positions of sources
* `ageom[i].sz` : z positions of sources
* `ageom[i].rx` : x positions of receivers
* `ageom[i].rz` : z positions of receivers
* `ageom[i].ns` : number of sources 
* `ageom[i].nr` : number of receivers


```julia
ageom=AGeom(mgrid, attrib, SSrcs(10), Recs(20))
```
Return pre-defined acquisition geometries (with 10 sources and 20 receivers), based on `attrib`, on an input mesh `mgrid`.
The sources and receivers positions are same for all supersources, and are not placed on the edges of the mesh.
Choose `attrib::Symbol`
* `=:xwell` cross-well acquisition;
* `=:surf` surface acquisition;
* `=:vsp` vertical seismic profiling;
* `=:rvsp`  reverse vertical seismic profiling;
* `=:downhole` downhole sources and receivers.
"""
AGeom=Array{AGeomss,1}



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

