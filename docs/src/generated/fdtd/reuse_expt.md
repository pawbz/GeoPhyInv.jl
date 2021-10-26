```@meta
EditURL = "<unknown>/test/fdtd/reuse_expt.jl"
```

````@example reuse_expt
using Revise
using GeoPhyInv
using Statistics
using Plots; gr();
using Test
````

This tutorial demonstrates how an instance of `SeisForwExpt` can be used iteratively
to perform forward modeling using `update!` for
various bundles of medium parameters.

### Grid

````@example reuse_expt
zgrid=range(-1000.0,stop=1000.0,length=201)
xgrid=range(-1000.0,stop=1000.0,length=201)
mgrid = [zgrid, xgrid];
nothing #hide
````

### Medium #1

````@example reuse_expt
vpb = [1500., 3500.] # bounds for vp
rhob = [1500., 3500.] # density bounds
medium=Medium(mgrid);
update!(medium, [:vp,:rho], [vpb,rhob]);
fill!(medium);
update!(medium, [:vp,:rho], randn_perc=5); # add some random noise
nothing #hide
````

### Medium #2

````@example reuse_expt
medium_new=similar(medium)
update!(medium_new, [:vp,:rho], randn_perc=5.) # add some noise
update!(medium_new, [:vp], constant_pert=500.) # perturb velocity
````

### AGeom

````@example reuse_expt
ageom=AGeom(medium.mgrid,:surf, SSrcs(3), Recs(30)); # surface seismic
nothing #hide
````

### Plotting #1

````@example reuse_expt
p1=heatmap(medium, :vp)
scatter!(ageom, SSrcs())
scatter!(ageom, Recs())
p2=heatmap(medium_new, :vp)
scatter!(ageom, SSrcs())
scatter!(ageom, Recs())
plot(p1, p2, size=(800,300))
````

### SrcWav

````@example reuse_expt
tgrid = range(0.0,stop=1.0,length=1000) # generate a temporal grid
wav = ricker(10.0, tgrid, tpeak=0.25,); # Choose a source wavelet
srcwav = SrcWav(tgrid, ageom, [:p]) # initialize
update!(srcwav, [:p], wav) # distribute to all supersources
````

### SeisForwExpt

````@example reuse_expt
pa=SeisForwExpt(FdtdAcou(),medium=medium, ageom=ageom, srcwav=srcwav, tgrid=tgrid, verbose=true);
nothing #hide
````

### Modeling #1

````@example reuse_expt
@time update!(pa);
d1=copy(pa[:data][:p])
p1=heatmap(pa[:data], grid=true);
nothing #hide
````

### Change medium in `pa` without memory allocation

````@example reuse_expt
update!(pa, medium_new)
````

### Modeling #2

````@example reuse_expt
@time update!(pa);
d2=copy(pa[:data][:p])
p2=heatmap(pa[:data], grid=true);
nothing #hide
````

Test

````@example reuse_expt
@test d1≠d2
````

### Plotting #2

````@example reuse_expt
plot(p1,p2, size=(500, 300))
````

