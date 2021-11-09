using BenchmarkTools
using GeoPhyInv
using Test

medium = Medium(:elastic_homo3D);

names(medium)

medium.ref

medium.bounds

mgrid = [range(0.0, stop = 10.0, step = 0.1), range(0.0, stop = 30.0, step = 0.2)];

mgrid = fill(range(-10, stop = 10.0, step = 0.1), 3)

medium = Medium(mgrid, [:vp, :rho, :vs]);

vpb = [2100.0, 2200.0];
vsb = [1500, 1700];
rhob = [2100.0, 2300.0];
update!(medium, [:vp, :vs, :rho], [vpb, vsb, rhob]);

fill!(medium);

medium[:Zp];

medium[:vp] .= 3000.0;
medium[:vs] .= 2000.0;

println(medium)

update!(medium, [:vp, :rho], randn_perc = 1.0);

