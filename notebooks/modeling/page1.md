```@meta
EditURL = "https://github.com/TRAVIS_REPO_SLUG/blob/master/"
```

load packages

```@example page1
using JuMIT
using Statistics
using Plots
pyplot()
```

create simple (almost) homogeneous acoustic model

```@example page1
model=J.Gallery.Seismic(:acou_homo1)
J.Models.Seismic_addon!(model, randn_perc=0.01)
```

a simple acquisition geometry

```@example page1
acqgeom = JuMIT.Gallery.Geom(model.mgrid,:xwell);
```

plot the model and source, receivers

```@example page1
p1=JP.seismic(model)
JP.geom!(acqgeom)
plot(p1)
```

generate time grid

```@example page1
tgrid = range(0.0,stop=2.0,length=1000)
```

Ricker wavelet

```@example page1
wav = JuMIT.Utils.Wavelets.ricker(10.0, tgrid, tpeak=0.25,);
```

distribute the same source wavelet to all the supsersources

```@example page1
acqsrc=JuMIT.Acquisition.Src_fixed(acqgeom.nss,1,[:P],wav,tgrid);
```

create `Fdtd.Param` object to prepare forward modelling
* npw corresponds to the number of independently propagating wavefields (1 in most cases)

```@example page1
pa=JuMIT.Fdtd.Param(npw=1,model=model,
	acqgeom=[acqgeom], acqsrc=[acqsrc],
	sflags=[2], rflags=[1],
	tgridmod=tgrid, verbose=true);
```

Once the `Param` object is created, do the modelling "without any memory allocations" using `mod!`

```@example page1
@time JuMIT.Fdtd.mod!(pa);
```

plot a record after modelling

```@example page1
pdata=plot(pa.c.data[1].d[1,1])
plot(pdata)
```

create new seismic model

```@example page1
model_new=J.Gallery.Seismic(:acou_homo1) # prepare another model
J.Models.Seismic_addon!(model_new, randn_perc=0.01)
J.Models.Seismic_addon!(model_new, constant_pert=0.03) # perturb the model
p2=JP.seismic(model_new) # plot new model
JP.geom!(acqgeom)
plot(p2)
```

Now, we the change the model in the `Param` object without memory allocation
This routine can be used during FWI,
where medium parameters are itertively updated in the same `Fdtd.Param` object

```@example page1
J.Fdtd.update_model!(pa.c, model_new)
```

run modelling now and plot data again

```@example page1
@time JuMIT.Fdtd.mod!(pa);
```

plot a record after modelling

```@example page1
plot!(pdata, pa.c.data[1].d[1,1])
plot(pdata)
```

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

