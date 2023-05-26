using Revise
using GeoPhyInv
using Statistics
#md using Plots; gr();
using Test

# This tutorial demonstrates how an instance of `SeisForwExpt` can be used iteratively 
# to perform forward modeling using `update!` for 
# various bundles of medium parameters.


# ### Grid
zgrid=range(-1000.0,stop=1000.0,length=201)
xgrid=range(-1000.0,stop=1000.0,length=201)
mgrid = [zgrid, xgrid];


# ### Medium #1
vpb = [1500., 3500.] # bounds for vp
rhob = [1500., 3500.] # density bounds
medium=Medium(mgrid);
update!(medium, [:vp,:rho], [vpb,rhob]);
fill!(medium);
update!(medium, [:vp,:rho], randn_perc=5); # add some random noise

# ### Medium #2
medium_new=similar(medium)
update!(medium_new, [:vp,:rho], randn_perc=5.) # add some noise
update!(medium_new, [:vp], constant_pert=500.) # perturb velocity

# ### AGeom
ageom=AGeom(medium.grid,:surf, SSrcs(3), Recs(30)); # surface seismic

# ### Plotting #1
#md p1=heatmap(medium, :vp) 
#md scatter!(ageom, SSrcs())
#md scatter!(ageom, Recs())
#md p2=heatmap(medium_new, :vp) 
#md scatter!(ageom, SSrcs())
#md scatter!(ageom, Recs())
#md plot(p1, p2, size=(800,300))

# ### Srcs
tgrid = range(0.0,stop=1.0,length=1000) # generate a temporal grid
wav = ricker(10.0, tgrid, tpeak=0.25,); # Choose a source wavelet
srcwav = Srcs(tgrid, ageom, [:p]) # initialize 
update!(srcwav, [:p], wav) # distribute to all supersources

# ### SeisForwExpt
pa=SeisForwExpt(FdtdAcoustic(),medium=medium, ageom=ageom, srcwav=srcwav, tgrid=tgrid, verbose=true);

# ### Modeling #1
@time update!(pa);
d1=copy(pa[:data][:p])
#md p1=heatmap(pa[:data], grid=true);

# ### Change medium in `pa` without memory allocation
update!(pa, medium_new)

# ### Modeling #2
@time update!(pa);
d2=copy(pa[:data][:p])
#md p2=heatmap(pa[:data], grid=true);

# Test
@test d1≠d2

# ### Plotting #2
#md plot(p1,p2, size=(500, 300))
