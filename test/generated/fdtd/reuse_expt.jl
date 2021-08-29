using GeoPhyInv
using Statistics
using Test

zgrid=range(-1000.0,stop=1000.0,length=201)
xgrid=range(-1000.0,stop=1000.0,length=201)
mgrid = [zgrid, xgrid];

vpb = [1500., 3500.] # bounds for vp
rhob = [1500., 3500.] # density bounds
medium=Medium(mgrid);
update!(medium, [:vp,:rho], [vpb,rhob]);
fill!(medium);
update!(medium, [:vp,:rho], randn_perc=5); # add some random noise

medium_new=similar(medium)
update!(medium_new, [:vp,:rho], randn_perc=5.) # add some noise
update!(medium_new, [:vp], constant_pert=500.) # perturb velocity

ageom=AGeom(medium.mgrid,:surf, SSrcs(3), Recs(30)); # surface seismic

tgrid = range(0.0,stop=1.0,length=1000) # generate a temporal grid
wav = ricker(10.0, tgrid, tpeak=0.25,); # Choose a source wavelet
srcwav = SrcWav(tgrid, ageom, [:P]) # initialize
update!(srcwav, [:P], wav) # distribute to all supersources

pa=SeisForwExpt(FdtdAcou(),medium=medium, ageom=ageom, srcwav=srcwav, tgrid=tgrid, verbose=true);

@time update!(pa);
d1=copy(pa[:data][:P])

update!(pa, medium_new)

@time update!(pa);
d2=copy(pa[:data][:P])

@test d1≠d2

