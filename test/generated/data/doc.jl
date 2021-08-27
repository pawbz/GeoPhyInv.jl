using GeoPhyInv
using Test
using Random

mgrid = [range(0.0, stop=10.,step=0.1), range(0.0, stop=30.,step=0.2)];
ageom=AGeom(mgrid, SSrcs(10), Srcs(10), Recs(10));

tgrid=range(0, stop=1.0, step=0.1);

data=Data(tgrid, ageom, [:P,:vz]);

Random.randn!(data[3][:P]);

