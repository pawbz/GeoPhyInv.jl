# This page was generated on DATEOFTODAY
using GeoPhyInv
using Statistics
#md using Plots; gr();
using Test

# This tutorial demonstrates how an instance of `SeisForwExpt` can be used iteratively 
# to perform forward modeling using `update!` for 
# various bundles of medium parameters.


# ### Medium #1
medium = Medium(:elastic_homo2D, 5)

# ### Medium #2
medium_new = similar(medium)
update!(medium_new, [:vp, :rho, :vs], rectangle = [[-500, -500], [500, 500]], perc = 5.0) # perturb velocity

# ### AGeom
ageom = AGeom(medium.mgrid, :surf, SSrcs(3), Recs(100)); # surface seismic

# ### Plotting #1
#md p1=heatmap(medium, :rho, seriescolor=cgrad(colorschemes[:roma])) 
#md scatter!(ageom, SSrcs())
#md scatter!(ageom, Recs())
#md p2=heatmap(medium_new, :rho, seriescolor=cgrad(colorschemes[:roma])) 
#md scatter!(ageom, SSrcs())
#md scatter!(ageom, Recs())
#md plot(p1, p2, size=(800,300))

# ### SrcWav
tgrid = range(0.0, stop = 2.0, length = 2500) # generate a temporal grid
wav = ricker(15.0, tgrid, tpeak = 0.25); # Choose a source wavelet
srcwav = SrcWav(tgrid, ageom, [:vz]) # initialize 
update!(srcwav, [:vz], wav) # distribute to all supersources

# ### SeisForwExpt
pa = SeisForwExpt(
    FdtdElastic(),
    medium = medium,
    ageom = ageom,
    srcwav = srcwav,
    tgrid = tgrid,
    rfields = [:vz],
    verbose = true,
);

# ### Modeling #1
@time update!(pa);
d1 = copy(pa[:data][:vz])
#md p1=heatmap(pa[:data], :vz, 99, 99, grid=true, legend=:none, seriescolor=cgrad(colorschemes[:seismic]));

# ### Change medium in `pa` without memory allocation
update!(pa, medium_new)

# ### Modeling #2
@time update!(pa);
d2 = copy(pa[:data][:vz])
#md p2=heatmap(pa[:data], :vz, 99, 99, grid=true, legend=:none, seriescolor=cgrad(colorschemes[:seismic]));

# Test
@test d1 ≠ d2

# ### Plotting #2
#md plot(p1,p2, size=(500, 300))
