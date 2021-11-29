```@meta
EditURL = "<unknown>/test/ageom/doc.jl"
```

This page was generated on 2021-11-29

````@example doc
using GeoPhyInv
using Test
````

# Intro
```@docs
AGeom
```

# Examples

Lets create a 2-D `mgrid` for the experiment.

````@example doc
mgrid2 = [range(0.0, stop = 10.0, step = 0.1), range(0.0, stop = 30.0, step = 0.2)];
nothing #hide
````

Then initialize an acquisition on `mgrid`, where the positions will be randomly chosen.

````@example doc
ageom2 = AGeom(mgrid2, SSrcs(10), Srcs(10), Recs(10));
nothing #hide
````

Similarly, we can do it for 3D too.

````@example doc
mgrid3 = fill(range(-10, stop=10, step=0.01),3);
ageom3 = AGeom(mgrid3, SSrcs(10), Srcs(10), Recs(10));
nothing #hide
````

You can check the number of dimensions.

````@example doc
ndims(ageom2), ndims(ageom3) == 2, 3
````

For 2D, we can also initialize with one of the predefined acquisitions.

````@example doc
ageom2 = AGeom(mgrid2, :xwell, SSrcs(10), Recs(10));
nothing #hide
````

By default, the sources and receivers are randomly placed
on the grid.
If unsure, test it.

````@example doc
@test (ageom2 ∈ mgrid2)
@test (ageom3 ∈ mgrid3)
````

The source and receiver positions can be updated as desired.

````@example doc
update!(ageom2[1], Srcs(), [0, 1], [10, 20]);
update!(ageom2[1], Recs(), [0, 0], [10, 20]);
update!(ageom2, SSrcs(), [0, 1], [10, 20]);
nothing #hide
````

It is easy to combine supersources. Now `ageom2` has 20 supersources.

````@example doc
ageom2_new = vcat(ageom2, ageom2);
nothing #hide
````

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

