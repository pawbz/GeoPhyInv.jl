using BenchmarkTools
using GeoPhyInv
using Test

mod=Medium(:marmousi2);

names(mod)

mod.ref

mod.bounds

mgrid = [range(0.0, stop=10.,step=0.1), range(0.0, stop=30.,step=0.2)];

mod = Medium(mgrid, [:vp,:rho,:vs]);

vpb=[2100.,2200.]; vsb=[1500, 1700]; rhob=[2100., 2300.];
update!(mod, [:vp,:vs,:rho], [vpb, vsb, rhob]);

fill!(mod);

mod[:Zp];

mod[:vp].=3000.;
mod[:vs].=2000.;

update!(mod, [:vp,:rho], randn_perc=1.);

