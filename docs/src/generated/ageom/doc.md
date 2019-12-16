```@meta
EditURL = "<unknown>/ageom/doc.jl"
```

```@example doc
using GeoPhyInv
using Test
using SparseArrays
```

# Intro
```@docs
AGeom
```

# Examples

Lets create a `mgrid` for the experiment.

```@example doc
mgrid = [range(0.0, stop=10.,step=0.1), range(0.0, stop=30.,step=0.2)];
nothing #hide
```

Then initialize an acquisition on `mgrid`.

```@example doc
ageom=AGeom(mgrid, SSrcs(10), Srcs(10), Recs(10));
nothing #hide
```

Otherwise, initialize with one of the predefined acquisitions.

```@example doc
ageom=AGeom(mgrid, :xwell, SSrcs(10), Recs(10));
nothing #hide
```

By default, the sources and receivers are randomly placed
on the grid.
If unsure, test it.

```@example doc
@test (ageom ∈ mgrid)
```

The source and receiver positions can be updated as desired.

```@example doc
update!(ageom[1], Srcs(), [0,1], [10,20],);
update!(ageom[1], Recs(), [0,0], [10,20],);
update!(ageom, SSrcs(), [0,1], [10,20], );
nothing #hide
```

It is easy to combine supersources. Now `ageom2` has 20 supersources.

```@example doc
ageom2=vcat(ageom, ageom);
nothing #hide
```

# Methods
```@docs
update!(::AGeomss, ::Srcs)
update!(::AGeomss, ::Recs)
update!(::AGeom, ::SSrcs)
update!(::AGeom, ::Recs)
Base.in(::AGeom, ::AbstractVector{StepRangeLen})
Base.isequal(::AGeom, ::AGeom)
SparseArrays.SparseMatrixCSC(::AGeomss,::AbstractVector{StepRangeLen})
```

