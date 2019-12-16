```@meta
EditURL = "<unknown>/data/doc.jl"
```

```@example doc
using GeoPhyInv
using Test
using Random
```

# Intro

```@docs
Data
```

# Examples
Define an acquisition geometry.

```@example doc
mgrid = [range(0.0, stop=10.,step=0.1), range(0.0, stop=30.,step=0.2)];
ageom=AGeom(mgrid, SSrcs(10), Srcs(10), Recs(10));
nothing #hide
```

Need a time grid.

```@example doc
tgrid=range(0, stop=1.0, step=0.1);
nothing #hide
```

Lets initialize records for `:P` and `:Vz` fields.

```@example doc
data=Data(tgrid, ageom, [:P,:Vz]);
nothing #hide
```

Fill the `:P` field of 3rd supersource with random numbers.

```@example doc
Random.randn!(data[3][:P]);
nothing #hide
```

# Methods
Methods listed for `SrcWav` can used on instances of `Data` too.

