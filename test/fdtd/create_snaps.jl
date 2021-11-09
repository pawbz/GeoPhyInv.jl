using GeoPhyInv
using Statistics
#jl @init_parallel_stencil(2, false, Float32, 4)
#md using Plots; gr();

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
#md p1=heatmap(medium, :vp)
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

pa = SeisForwExpt(
    FdtdAcou(),
    medium = medium,
    ageom = ageom,
    srcwav = srcwav,
    snaps_field = :p,
    tsnaps = [0.4, 0.5, 0.8],
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

#md function plot_snapshots()
#md     p=[]
#md     for j in 1:2, i in 1:3
#md     	push!(p,heatmap(pa[:snaps, j][i], xlabel="x", ylabel="z", legend=:none))
#md     end
#md     Plots.plot(p..., layout=(2,3), size=(1200,600))
#md end
#md plot_snapshots()


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
# Perform modelling and extract snaps
pa = SeisForwExpt(
    FdtdElastic(),
    medium = medium,
    ageom = ageom,
    srcwav = srcwav,
    snaps_field = :vz,
    tsnaps = [0.4, 0.5, 0.8],
    rfields = [:vz],
    tgrid = tgrid,
    verbose = true,
)
@time update!(pa);

snaps = pa[:snaps, 1];
snaps = pa[:snaps, 2];

# ### Plotting is vs model
#md p1=heatmap(medium, :vs)
#md scatter!(ageom, SSrcs())
#md scatter!(ageom, Recs())
#md plot(p1)

# ### Plotting snapshots after elastic simulation
#md plot_snapshots() # P and S waves!
