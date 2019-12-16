```@meta
EditURL = "<unknown>/fdtd/reuse_expt.jl"
```

This tutorial first creates a variable `pa`, i.e. allocates
necessary memory to perform `SeisForwExpt`.
Then, we perform forward modeling using an in-place
function `mod!`. Finally, we will update `pa` with a different
subsurface model and re-run the modeling task with no additional memory
allocation.
The ability to iteratively run the forward modeling task on
various subsurface models is necessary while implementing inversion
algorithms.

### Load packages

```@example reuse_expt
using GeoPhyInv
using Statistics
using Plots
using Test
```

### Setup a spatial grid

```@example reuse_expt
zgrid=range(-1000.0,stop=1000.0,length=201)
xgrid=range(-1000.0,stop=1000.0,length=201)
mgrid = [zgrid, xgrid]
@info string("spatial sampling intervals (dz,dx)=", step.(mgrid))
```

### Allocate a `Seismic` model, and adjust bounds

```@example reuse_expt
vpb = [1500., 3500.] # bounds for vp
rhob = [1500., 3500.] # density bounds

model=Medium(mgrid)
update!(model, [:vp,:rho], [vpb,rhob])
fill!(model)
@info "a seismic model is created"
```

### Add some noise to the model (optional)

```@example reuse_expt
update!(model, [:vp,:rho], randn_perc=0.01); # add some random noise
nothing #hide
```

### A surface acquisition ageometry

```@example reuse_expt
ageom=AGeom(model.mgrid,:surf, SSrcs(3), Recs(30));
nothing #hide
```

### Plot the model and source, receivers

```@example reuse_expt
p1=heatmap(model, :vp)
scatter!(ageom, SSrcs())
scatter!(ageom, Recs())
plot(p1)
```

### Generate a temporal grid

```@example reuse_expt
tgrid = range(0.0,stop=2.0,length=1000)
```

### Choose a source wavelet

```@example reuse_expt
wav = ricker(10.0, tgrid, tpeak=0.25,);
nothing #hide
```

### Distribute the same source wavelet to all the supsersources

```@example reuse_expt
srcwav = SrcWav(tgrid, ageom, [:P])
update!(srcwav, [:P], wav)
```

create `Fdtd.Param` object to prepare forward modelling
* npw corresponds to the number of independently propagating wavefields (1 in most cases)
Once the `Param` object is created, do the modelling "without any memory allocations" using `mod!`

```@example reuse_expt
pa=SeisForwExpt(Fdtd(),npw=1,model=model,
	ageom=[ageom], srcwav=[srcwav],
	sflags=[2], rflags=[1],
	tgridmod=tgrid, verbose=true);

@time update!(pa);
nothing #hide
```

plot a record after modelling

```@example reuse_expt
pdata=heatmap(pa[:data])
plot(pdata)
```

create new seismic model

```@example reuse_expt
model_new=Medium(:acou_homo1) # prepare another model
update!(model_new, [:vp,:rho], randn_perc=0.01)
update!(model_new, [:vp,:rho], constant_pert=0.03) # perturb the model
p2=heatmap(model_new, :vp) # plot new model
scatter!(ageom, SSrcs())
scatter!(ageom, Recs())
plot(p2)
```

Now, we the change the model in the `Param` object without memory allocation
This routine can be used during FWI,
where medium parameters are itertively updated in the same `Fdtd.Param` object

```@example reuse_expt
update!(pa, model_new)
```

run modelling now and plot data again

```@example reuse_expt
@time update!(pa);
nothing #hide
```

plot a record after modelling

```@example reuse_expt
heatmap!(pdata, pa[:data])
plot(pdata)
```

