```@meta
EditURL = "https://github.com/TRAVIS_REPO_SLUG/blob/master/"
```

A forward experiment, where the seismic data are generated
using some models and acquisition parameters from our gallery.
Forward modeling consists of a finite-difference simulation of the acoustic wave-equation.
We specifically aim to save the snapshots, at given time steps in `SeisForwExpt`.

### Loading some packages

```@example create_snaps
using GeoPhyInv
using Statistics
```

### Setting up the variables necessary to create the `Expt`

```@example create_snaps
model=GIPh.Gallery.Seismic(:acou_homo1); # load a simple homogeneous acoustic model from the gallery
GIPh.Models.Seismic_addon!(model, randn_perc=0.01); # add some random noise to the model
acqgeom=GIPh.Gallery.Geom(model.mgrid,:xwell); # load a simple acquisition geometry using `mgrid` of the seismic model
tgrid = range(0.0,stop=2.0,length=2000) # generate a time grid
wav = GIPh.Utils.Wavelets.ricker(10.0, tgrid, tpeak=0.25,); # ricker wavelet
acqsrc=GIPh.Acquisition.Src_fixed(acqgeom.nss,1,[:P],wav,tgrid); # distribute the same source wavelet to all the supsersources
@info "We are ready for the modeling."
```

### Final step

One can plot the model, source and receivers using these commands:
`using Plots;`
`p1=JP.seismic(model);`
`JP.geom!(acqgeom);`
`plot(p1);`
Now we have all the required variables to create `SeisForwExpt` object and
prepare the forward modelling.
While creating, we switched the `snaps_flag` on, and instructed recording field at
`tsnaps`.
Once the `Expt` object is created, do the modelling "without approximately any
memory allocations" using `mod!`

```@example create_snaps
paE=SeisForwExpt(model=model,
	acqgeom=[acqgeom], acqsrc=[acqsrc],
	snaps_flag=true,
	tsnaps=[0.3, 0.4, 0.5],
	tgridmod=tgrid, verbose=true);

@time mod!(paE);
```

### Extracting snaps from Expt

```@example create_snaps
snaps=paE[:snaps,1]; # extracting snaps of the first supersource
snaps=paE[:snaps,2]; # second supersource
@info "The dimensions of the snaps is [nz,nx,nt]."
```

We can now plot snapshots using these commands:
`p1=[heatmap(snaps[:,:,ii]) for ii in 1:3];`
`plot(p1..., layout=(1,3), aspect_ratio=:equal)`

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

