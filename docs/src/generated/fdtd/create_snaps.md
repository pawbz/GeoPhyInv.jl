```@meta
EditURL = "<unknown>/test/fdtd/create_snaps.jl"
```

````@example create_snaps
using GeoPhyInv
using Statistics
using Plots; gr();
nothing #hide
````

### Medium

````@example create_snaps
medium = Medium(:acou_homo1, 5); # load a simple homogeneous acoustic medium from the gallery
update!(medium, [:vp, :rho], randn_perc = 5); # add some random noise to the medium
nothing #hide
````

### AGeom

````@example create_snaps
ageom = AGeom(medium.mgrid, :xwell, SSrcs(2)); # load a simple acquisition using `mgrid` of the medium
nothing #hide
````

### SrcWav

````@example create_snaps
tgrid = range(0.0, stop = 2.0, length = 2000); # generate a time grid
wav = ricker(10.0, tgrid, tpeak = 0.25); # ricker wavelet
srcwav = SrcWav(tgrid, ageom, [:p]);
update!(srcwav, [:p], wav);
nothing #hide
````

### Plotting #1

````@example create_snaps
p1=heatmap(medium, :vp)
scatter!(ageom, SSrcs())
scatter!(ageom, Recs())
plot(p1)
````

### SeisForwExpt
Now we have all the required variables to create `SeisForwExpt` object and
prepare the forward modeling.
While creating, we switched the `snaps_flag` on, and instruct recording field at
`tsnaps`.
Once the `Expt` object is created, do the modeling "without approximately any
memory allocations" using `update!`

````@example create_snaps
pa = SeisForwExpt(
    FdtdAcou(),
    medium = medium,
    ageom = ageom,
    srcwav = srcwav,
    snaps_field = :p,
    tsnaps = [0.3, 0.4, 0.5],
    tgrid = tgrid,
    verbose = true,
);
nothing #hide
````

### Modeling

````@example create_snaps
@time update!(pa);
nothing #hide
````

### Extracting snaps

````@example create_snaps
snaps = pa[:snaps, 1]; # extracting snaps of the first supersource
snaps = pa[:snaps, 2]; # second supersource
nothing #hide
````

### Plotting #2

````@example create_snaps
p=[]
for ii in 1:3
	push!(p,heatmap(medium.mgrid[2], medium.mgrid[1], snaps[ii], xlabel="x", ylabel="z"))
end
plot(p..., layout=(1,3), size=(1100,250))
````

