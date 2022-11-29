using GeoPhyInv
using Statistics
#jl @init_parallel_stencil(2, false, Float32, 4)
#md using Plots; gr();
#md using ColorSchemes

# ### Medium
medium = Medium(:acou_homo2D, 5); # load a simple homogeneous acoustic medium from the gallery
update!(medium, [:vp, :rho], randn_perc = 5); # add some random noise to the medium
println(medium)

# ### AGeom
ageom = AGeom(medium.mgrid, :xwell, SSrcs(2)); # load a simple acquisition using `mgrid` of the medium
println(ageom)

# ### SrcWav
tgrid = range(0.0, stop = 2.0, length = 2500); # generate a time grid
wav = ricker(15.0, tgrid, tpeak = 0.25); # ricker wavelet
srcwav = SrcWav(tgrid, ageom, [:p]);
update!(srcwav, [:p], wav);

# ### Plotting the `vp` model
#md p1=heatmap(medium, :vp, seriescolor=cgrad(:roma))
#md scatter!(ageom, SSrcs())
#md scatter!(ageom, Recs())
#md plot(p1)


# ## Acoustic
# Now we have all the required variables to create `SeisForwExpt` object and 
# prepare the forward modeling.
# While creating, we switched the `snaps_flag` on, and instruct recording field at
# `tsnaps`.
# Once the `Expt` object is created, do the modeling without *additional* any 
# memory allocations using `update!`

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

# ### Update `pa` to perform modelling; prints time and allocations
@time update!(pa);

# ### Extracting snaps
snaps1 = pa[:snaps, 1]; # extracting snaps of the first supersource
snaps2 = pa[:snaps, 2]; # second supersource

# ### Plotting pressure snapshots after acoustic simulation
#md function plot_snapshots(name)
#md     clip_perc = 90
#md     pmax = maximum([maximum(pa[:snaps][it]) for it = 1:length(tsnaps)])
#md     pmax -= pmax * 0.01 * clip_perc
#md 
#md     anim = @animate for it = 1:length(tsnaps)
#md         plot(
#md             [
#md                 (
#md                     f = pa[:snaps, is][it];
#md                     heatmap(
#md                         f,
#md                         aspect_ratio = 1,
#md                         xlims = (1, size(f, 2)),
#md                         ylims = (1, size(f, 1)),
#md                         yflip = true,
#md                         seriescolor = cgrad(:seismic),
#md                         clims = (-pmax, pmax),
#md                         legend = false,
#md                     )
#md                 ) for is = 1:2
#md             ]...,
#md             title = string("snapshot at: ", round(tsnaps[it], digits = 5), " s"),
#md         )
#md     end
#md     return gif(anim, string(name,".gif"), fps = 1)
#md end
#md plot_snapshots("acoustic_snaps")


# ## Elastic
# Similarly, we can do the elastic wave-equation modeling
medium = Medium(:elastic_homo2D, 5); # load a simple homogeneous elastic medium from the gallery
update!(medium, [:vp, :vs, :rho], randn_perc = 5); # add some random noise to the medium
println(medium)
ageom = AGeom(medium.mgrid, :xwell, SSrcs(2)); # load a simple acquisition using `mgrid` of the medium
println(ageom)
# Add source term on `:vz` grid
srcwav = SrcWav(tgrid, ageom, [:vz]);
update!(srcwav, [:vz], wav);
# Perform modelling 
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

# ### Plotting is vs model
#md p1=heatmap(medium, :vs, seriescolor=cgrad(:roma))
#md scatter!(ageom, SSrcs())
#md scatter!(ageom, Recs())
#md plot(p1)

# ### Plotting snapshots after elastic simulation
#md plot_snapshots("elastic_snaps") # P and S waves!
