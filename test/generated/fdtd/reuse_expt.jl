using GeoPhyInv
using Statistics
using Test

medium = Medium(:elastic_homo2D, 5)

medium_new = similar(medium)
update!(medium_new, [:vp, :rho, :vs], rectangle = [[-500, -500], [500, 500]], perc = 5.0) # perturb velocity

ageom = AGeom(medium.mgrid, :surf, SSrcs(3), Recs(100)); # surface seismic

tgrid = range(0.0, stop = 2.0, length = 2500) # generate a temporal grid
wav = ricker(15.0, tgrid, tpeak = 0.25); # Choose a source wavelet
srcwav = SrcWav(tgrid, ageom, [:vz]) # initialize
update!(srcwav, [:vz], wav) # distribute to all supersources

pa = SeisForwExpt(
    FdtdElastic(),
    medium = medium,
    ageom = ageom,
    srcwav = srcwav,
    tgrid = tgrid,
    rfields = [:vz],
    verbose = true,
);

@time update!(pa);
d1 = copy(pa[:data][:vz])

update!(pa, medium_new)

@time update!(pa);
d2 = copy(pa[:data][:vz])

@test d1 ≠ d2

plot(p1,p2, size=(500, 300))

