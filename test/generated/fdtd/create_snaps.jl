using GeoPhyInv
using Statistics
@init_parallel_stencil(2, false, Float32, 4)

medium = Medium(:acou_homo2D, 5); # load a simple homogeneous acoustic medium from the gallery
update!(medium, [:vp, :rho], randn_perc = 5); # add some random noise to the medium
println(medium)

ageom = AGeom(medium.mgrid, :xwell, SSrcs(2)); # load a simple acquisition using `mgrid` of the medium
println(ageom)

tgrid = range(0.0, stop = 2.0, length = 2500); # generate a time grid
wav = ricker(15.0, tgrid, tpeak = 0.25); # ricker wavelet
srcwav = SrcWav(tgrid, ageom, [:p]);
update!(srcwav, [:p], wav);

tsnaps=tgrid[1:div(length(tgrid),20):end] # store 20 snapshots
pa = SeisForwExpt(
    FdtdAcoustic(),
    medium = medium,
    ageom = ageom,
    srcwav = srcwav,
    snaps_field = :p,
    tsnaps = tsnaps,
    tgrid = tgrid,
    rfields = [:p],
    verbose = true,
);

@time update!(pa);

snaps1 = pa[:snaps, 1]; # extracting snaps of the first supersource
snaps2 = pa[:snaps, 2]; # second supersource

medium = Medium(:elastic_homo2D, 5); # load a simple homogeneous elastic medium from the gallery
update!(medium, [:vp, :vs, :rho], randn_perc = 5); # add some random noise to the medium
println(medium)
ageom = AGeom(medium.mgrid, :xwell, SSrcs(2)); # load a simple acquisition using `mgrid` of the medium
println(ageom)

srcwav = SrcWav(tgrid, ageom, [:vz]);
update!(srcwav, [:vz], wav);

pa = SeisForwExpt(
    FdtdElastic(),
    medium = medium,
    ageom = ageom,
    srcwav = srcwav,
    snaps_field = :vz,
    tsnaps = tsnaps,
    rfields = [:vz],
    tgrid = tgrid,
    verbose = true,
)
@time update!(pa);

snaps = pa[:snaps, 1];
snaps = pa[:snaps, 2];

