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
medium = Medium(:acou_homo2D, 5); # load a simple homogeneous acoustic medium from the gallery
update!(medium, [:vp, :rho], randn_perc = 5); # add some random noise to the medium
println(medium)
````

### AGeom

````@example create_snaps
ageom = AGeom(medium.mgrid, :xwell, SSrcs(2)); # load a simple acquisition using `mgrid` of the medium
println(ageom)
````

### SrcWav

````@example create_snaps
tgrid = range(0.0, stop = 2.0, length = 2500); # generate a time grid
wav = ricker(15.0, tgrid, tpeak = 0.25); # ricker wavelet
srcwav = SrcWav(tgrid, ageom, [:p]);
update!(srcwav, [:p], wav);
nothing #hide
````

### Plotting the `vp` model

````@example create_snaps
p1=heatmap(medium, :vp)
scatter!(ageom, SSrcs())
scatter!(ageom, Recs())
plot(p1)
````

## Acoustic
Now we have all the required variables to create `SeisForwExpt` object and
prepare the forward modeling.
While creating, we switched the `snaps_flag` on, and instruct recording field at
`tsnaps`.
Once the `Expt` object is created, do the modeling without *additional* any
memory allocations using `update!`

````@example create_snaps
pa = SeisForwExpt(
    FdtdAcou(),
    medium = medium,
    ageom = ageom,
    srcwav = srcwav,
    snaps_field = :p,
    tsnaps = [0.4, 0.5, 0.8],
    tgrid = tgrid,
    rfields = [:p],
    verbose = true,
);
nothing #hide
````

### Update `pa` to perform modelling; prints time and allocations

````@example create_snaps
@time update!(pa);
nothing #hide
````

### Extracting snaps

````@example create_snaps
snaps1 = pa[:snaps, 1]; # extracting snaps of the first supersource
snaps2 = pa[:snaps, 2]; # second supersource
nothing #hide
````

### Plotting pressure snapshots after acoustic simulation

````@example create_snaps
function plot_snapshots()
    p=[]
    for j in 1:2, i in 1:3
    	push!(p,heatmap(pa[:snaps, j][i], xlabel="x", ylabel="z", legend=:none))
    end
    Plots.plot(p..., layout=(2,3), size=(1200,600))
end
plot_snapshots()
````

## Elastic
Similarly, we can do the elastic wave-equation modeling

````@example create_snaps
medium = Medium(:elastic_homo2D, 5); # load a simple homogeneous elastic medium from the gallery
update!(medium, [:vp, :vs, :rho], randn_perc = 5); # add some random noise to the medium
println(medium)
ageom = AGeom(medium.mgrid, :xwell, SSrcs(2)); # load a simple acquisition using `mgrid` of the medium
println(ageom)
````

Add source term on `:vz` grid

````@example create_snaps
srcwav = SrcWav(tgrid, ageom, [:vz]);
update!(srcwav, [:vz], wav);
nothing #hide
````

Perform modelling and extract snaps

````@example create_snaps
pa = SeisForwExpt(
    FdtdElastic(),
    medium = medium,
    ageom = ageom,
    srcwav = srcwav,
    snaps_field = :vz,
    tsnaps = [0.4, 0.5, 0.8],
    rfields = [:vz],
    tgrid = tgrid,
    verbose = true,
)
@time update!(pa);

snaps = pa[:snaps, 1];
snaps = pa[:snaps, 2];
nothing #hide
````

### Plotting is vs model

````@example create_snaps
p1=heatmap(medium, :vs)
scatter!(ageom, SSrcs())
scatter!(ageom, Recs())
plot(p1)
````

### Plotting snapshots after elastic simulation

````@example create_snaps
plot_snapshots() # P and S waves!
````

