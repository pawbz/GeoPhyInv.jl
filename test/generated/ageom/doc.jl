using GeoPhyInv
using Test
using SparseArrays

mgrid2 = [range(0.0, stop = 10.0, step = 0.1), range(0.0, stop = 30.0, step = 0.2)];

ageom2 = AGeom(mgrid2, SSrcs(10), Srcs(10), Recs(10));

mgrid3 = fill(range(-10, stop=10, step=0.01),3);
ageom3 = AGeom(mgrid3, SSrcs(10), Srcs(10), Recs(10));

ndims(ageom2), ndims(ageom3) == 2, 3

ageom2 = AGeom(mgrid2, :xwell, SSrcs(10), Recs(10));

@test (ageom2 ∈ mgrid2)
@test (ageom3 ∈ mgrid3)

update!(ageom2[1], Srcs(), [0, 1], [10, 20]);
update!(ageom2[1], Recs(), [0, 0], [10, 20]);
update!(ageom2, SSrcs(), [0, 1], [10, 20]);

ageom2_new = vcat(ageom2, ageom2);

