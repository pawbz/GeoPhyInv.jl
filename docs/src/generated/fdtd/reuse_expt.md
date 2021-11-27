```@meta
EditURL = "<unknown>/test/fdtd/reuse_expt.jl"
```

This page was generated on 2021-11-27

````@example reuse_expt
using GeoPhyInv
using Statistics
using Plots; gr();
using ColorSchemes
using Test
````

This tutorial demonstrates how an instance of `SeisForwExpt` can be used iteratively
to perform forward modeling using `update!` for
various bundles of medium parameters.

### Medium #1

````@example reuse_expt
medium = Medium(:elastic_homo2D, 5)
````

### Medium #2

````@example reuse_expt
medium_new = similar(medium)
update!(medium_new, [:vp, :rho, :vs], rectangle = [[-500, -500], [500, 500]], perc = 5.0) # perturb velocity
````

### AGeom

````@example reuse_expt
ageom = AGeom(medium.mgrid, :surf, SSrcs(3), Recs(100)); # surface seismic
nothing #hide
````

### Plotting #1

````@example reuse_expt
p1=heatmap(medium, :rho, seriescolor=cgrad(colorschemes[:roma]))
scatter!(ageom, SSrcs())
scatter!(ageom, Recs())
p2=heatmap(medium_new, :rho, seriescolor=cgrad(colorschemes[:roma]))
scatter!(ageom, SSrcs())
scatter!(ageom, Recs())
plot(p1, p2, size=(800,300))
````

### SrcWav

````@example reuse_expt
tgrid = range(0.0, stop = 2.0, length = 2500) # generate a temporal grid
wav = ricker(15.0, tgrid, tpeak = 0.25); # Choose a source wavelet
srcwav = SrcWav(tgrid, ageom, [:vz]) # initialize
update!(srcwav, [:vz], wav) # distribute to all supersources
````

### SeisForwExpt

````@example reuse_expt
pa = SeisForwExpt(
    FdtdElastic(),
    medium = medium,
    ageom = ageom,
    srcwav = srcwav,
    tgrid = tgrid,
    rfields = [:vz],
    verbose = true,
);
nothing #hide
````

### Modeling #1

````@example reuse_expt
@time update!(pa);
d1 = copy(pa[:data][:vz])
p1=heatmap(pa[:data], :vz, 99, 99, grid=true, legend=:none, seriescolor=cgrad(colorschemes[:seismic]));
nothing #hide
````

### Change medium in `pa` without memory allocation

````@example reuse_expt
update!(pa, medium_new)
````

### Modeling #2

````@example reuse_expt
@time update!(pa);
d2 = copy(pa[:data][:vz])
p2=heatmap(pa[:data], :vz, 99, 99, grid=true, legend=:none, seriescolor=cgrad(colorschemes[:seismic]));
nothing #hide
````

Test

````@example reuse_expt
@test d1 ≠ d2
````

### Plotting #2

````@example reuse_expt
plot(p1,p2, size=(500, 300))
````

