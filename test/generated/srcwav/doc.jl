using GeoPhyInv
using Random
using LinearAlgebra
@init_parallel_stencil(3, false, Float32, 4)

mgrid = [range(0.0, stop=10.,step=0.1), range(0.0, stop=30.,step=0.2)];
ageom=AGeom(mgrid, SSrcs(10), Srcs(10), Recs(10));

tgrid=range(0, stop=1.0, step=0.1);

srcwav=Srcs(tgrid, ageom, [:p]);

randn!(srcwav[3][:p]);

x=randn(length(tgrid));
update!(srcwav, [:p,], x);

x1=randn(length(tgrid));
x2=randn(length(tgrid));
update!(srcwav[1], [:p,], x1);
update!(srcwav[2], [:p,], x2);

rmul!(srcwav, 2.0);

