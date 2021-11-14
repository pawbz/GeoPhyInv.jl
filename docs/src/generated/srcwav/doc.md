```@meta
EditURL = "<unknown>/test/srcwav/doc.jl"
```

This page was generated on 2021-11-14

````@example doc
using GeoPhyInv
using Random
using LinearAlgebra
````

# Intro

```@docs
SrcWav
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

Lets initialize records for `:p` field.

````@example doc
srcwav=SrcWav(tgrid, ageom, [:p]);
nothing #hide
````

Fill the `:P` field of 3rd supersource with random numbers.

````@example doc
randn!(srcwav[3][:p]);
nothing #hide
````

Often we want to populate the same source wavelet to all
the supersources and sources.

````@example doc
x=randn(length(tgrid));
update!(srcwav, [:p,], x);
nothing #hide
````

Populate two different wavelets for first and second supersources.

````@example doc
x1=randn(length(tgrid));
x2=randn(length(tgrid));
update!(srcwav[1], [:p,], x1);
update!(srcwav[2], [:p,], x2);
nothing #hide
````

Scale `srcwav` by a scalar overwriting it in-place.

````@example doc
rmul!(srcwav, 2.0);
nothing #hide
````

# Methods
Most of the methods listed below are also applicable to individual elements of `srcwav`.

```@docs
update!(GeoPhyInv.VNamedD, ::Vector{Symbol}, ::AbstractArray)
Base.reverse!(::GeoPhyInv.VNamedD)
Base.iszero(::GeoPhyInv.VNamedD)
Base.isequal(::GeoPhyInv.VNamedD, ::GeoPhyInv.VNamedD)
GeoPhyInv.issimilar(::GeoPhyInv.VNamedD, ::GeoPhyInv.VNamedD)
Base.vec(::GeoPhyInv.VNamedD)
Random.randn!(::GeoPhyInv.VNamedD)
Base.fill!(::GeoPhyInv.VNamedD, ::Float64)
LinearAlgebra.dot(::GeoPhyInv.VNamedD, ::GeoPhyInv.VNamedD)
LinearAlgebra.rmul!(::GeoPhyInv.VNamedD, ::Number)
Base.copyto!(::GeoPhyInv.VNamedD, ::GeoPhyInv.VNamedD)
Base.copyto!(::AbstractVector{Float64}, ::GeoPhyInv.VNamedD)
Base.copyto!(::GeoPhyInv.VNamedD, ::AbstractVector{Float64})
```

