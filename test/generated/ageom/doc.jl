using GeoPhyInv
using Test
using SparseArrays

mgrid = [range(0.0, stop=10.,step=0.1), range(0.0, stop=30.,step=0.2)];

ageom=AGeom(mgrid, SSrcs(10), Srcs(10), Recs(10));

ageom=AGeom(mgrid, :xwell, SSrcs(10), Recs(10));

@test (ageom ∈ mgrid)

update!(ageom[1], Srcs(), [0,1], [10,20],);
update!(ageom[1], Recs(), [0,0], [10,20],);
update!(ageom, SSrcs(), [0,1], [10,20], );

ageom2=vcat(ageom, ageom);

