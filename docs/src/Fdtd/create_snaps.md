```@meta
EditURL = "https://github.com/TRAVIS_REPO_SLUG/blob/master/"
```

Forward problem, where the seismic data are generated
using synthetic Earth models and the acquisition parameters

```@example create_snaps
#corresponding to a seismic experiment.
```

Forward modeling consists of a finite-difference simulation, followed

\$a\otimes b\$

the aim is to save snapshots, at given `t` in Fdtd
we start with loading some packages

```@example create_snaps
using GeoPhyInv
using Statistics
```

Load a simple (almost) homogeneous acoustic model from the gallery.

```@example create_snaps
model=GIPh.Gallery.Seismic(:acou_homo1)


J.Models.Seismic_addon!(model, randn_perc=0.01)
```

a simple acquisition geometry

```@example create_snaps
acqgeom = GeoPhyInv.Gallery.Geom(model.mgrid,:xwell);
```

plot the model and source, receivers

```@example create_snaps
p1=JP.seismic(model)
JP.geom!(acqgeom)
#plot(p1)
```

generate time grid

```@example create_snaps
tgrid = range(0.0,stop=2.0,length=1000)
```

Ricker wavelet

```@example create_snaps
wav = GeoPhyInv.Utils.Wavelets.ricker(10.0, tgrid, tpeak=0.25,);
```

distribute the same source wavelet to all the supsersources

```@example create_snaps
acqsrc=GeoPhyInv.Acquisition.Src_fixed(acqgeom.nss,1,[:P],wav,tgrid);
```

create `Fdtd.Param` object to prepare forward modelling
* npw corresponds to the number of independently propagating wavefields (1 in most cases)

```@example create_snaps
pa=GeoPhyInv.Fdtd.Param(npw=1,model=model,
	acqgeom=[acqgeom], acqsrc=[acqsrc],
	sflags=[2], rflags=[1],
	snaps_flag=true,
	tsnaps=[0.3, 0.4, 0.5],
	tgridmod=tgrid, verbose=true);
```

Once the `Param` object is created, do the modelling "without any memory allocations" using `mod!`

```@example create_snaps
@time GeoPhyInv.Fdtd.mod!(pa);
```

plotting snapshots of the first supersource

```@example create_snaps
p1=[heatmap(pa.p.localpart.ss[1].snaps[:,:,ii]) for ii in 1:3]
#plot(p1..., layout=(1,3), aspect_ratio=:equal)
```

plotting snapshots of the second supersource

```@example create_snaps
p1=[heatmap(pa.p.localpart.ss[2].snaps[:,:,ii]) for ii in 1:3]
#plot(p1..., layout=(1,3), aspect_ratio=:equal)
```

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*
