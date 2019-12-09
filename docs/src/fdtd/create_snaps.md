```@meta
EditURL = "<unknown>/fdtd/create_snaps.jl"
```

We aim to save the snapshots, at given time steps in `SeisForwExpt`.

### Loading some packages

```@example create_snaps
using GeoPhyInv
using Statistics
using Gadfly
```

### Setting up the variables necessary to create the `Expt`

```@example create_snaps
model=Medium(:acou_homo1); # load a simple homogeneous acoustic model from the gallery
update!(model, [:vp,:rho], randn_perc=0.1); # add some random noise to the model
ageom=AGeom(model.mgrid,:xwell, SSrcs(2)); # load a simple acquisition ageometry using `mgrid` of the seismic model
tgrid = range(0.0,stop=2.0,length=2000) # generate a time grid
wav = ricker(10.0, tgrid, tpeak=0.25,); # ricker wavelet

srcwav = SrcWav(tgrid, ageom, [:P])
update!(srcwav, [:P], wav)
@info "We are ready for the modeling."
```

### Final step

One can plot the model, source and receivers using these commands:
`using Plots;`
`p1=JP.seismic(model);`
`JP.ageom!(ageom);`
`plot(p1);`
Now we have all the required variables to create `SeisForwExpt` object and
prepare the forward modelling.
While creating, we switched the `snaps_flag` on, and instructed recording field at
`tsnaps`.
Once the `Expt` object is created, do the modelling "without approximately any
memory allocations" using `mod!`

```@example create_snaps
paE=SeisForwExpt(Fdtd(),model=model,
	ageom=[ageom], srcwav=[srcwav],
	snaps_flag=true,
	tsnaps=[0.3, 0.4, 0.5],
	tgridmod=tgrid, verbose=true);

@time update!(paE);
nothing #hide
```

### Extracting snaps from `Expt`

```@example create_snaps
snaps=paE[:snaps,1]; # extracting snaps of the first supersource
snaps=paE[:snaps,2]; # second supersource
@info string("The dimensions of the snaps are (nz,nx,nt)=", size(snaps))
```

We can now plot snapshots using these commands:

```@example create_snaps
p=[]
Gadfly.set_default_plot_size(20cm, 5cm)
for ii in 1:3
	push!(p, spy(snaps[:,:,ii], Guide.xlabel("x"), Guide.ylabel("z")))
end
hstack(p...)
```

