```@meta
EditURL = "<unknown>/test/data/doc.jl"
```

This page was generated on 2022-03-18

````@example doc
using GeoPhyInv
using Test
using Random
````

# Intro

```@docs
Records
```

# Examples
Define an acquisition geometry.

````@example doc
mgrid = [range(0.0, stop=10.,step=0.1), range(0.0, stop=30.,step=0.2)];
ageom=AGeom(mgrid, SSrcs(10), Srcs(10), Recs(10));
nothing #hide
````

Need a time grid.

````@example doc
tgrid=range(0, stop=1.0, step=0.1);
nothing #hide
````

Lets initialize records for `:P` and `:vz` fields.

````@example doc
data=Records(tgrid, ageom, [:p,:vz]);
nothing #hide
````

Fill the `:P` field of 3rd supersource with random numbers.

````@example doc
Random.randn!(data[3][:p]);
nothing #hide
````

# Methods
Methods listed for `SrcWav` can used on instances of `Records` too.

