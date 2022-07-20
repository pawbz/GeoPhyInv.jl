```@meta
EditURL = "<unknown>/test/fdtd/doc.jl"
```

This page was generated on 2022-03-18

````@example doc
using BenchmarkTools
using GeoPhyInv
using Test
using Luxor
````

## Staggered Grid
Luxor graphics are available to visualize 2-D grids.
3-D visualization is not available.
### Acoustic

````@example doc
@draw GeoPhyInv.luxor_mgrid(FdtdAcoustic())
````

We can avoid some clutter.

````@example doc
@draw GeoPhyInv.luxor_mgrid(FdtdAcoustic(), [:p, :vx, :vz])
````

### Elastic
Now lets visualize the grids of elastic simulation.

````@example doc
@draw GeoPhyInv.luxor_mgrid(FdtdElastic())
````

Grids for only stress and velocity fields.

````@example doc
@draw GeoPhyInv.luxor_mgrid(FdtdElastic(), [:tauxx, :tauxz, :tauzz, :vx, :vz])
````

## Methods

```@docs
SeisForwExpt
Base.getindex(::GeoPhyInv.PFdtd, ::Symbol, ::Int)
```

```@docs
update!(::GeoPhyInv.PFdtd)
update!(::GeoPhyInv.PFdtd, ::Medium)
update!(::GeoPhyInv.PFdtd, ::SrcWav, ::Any)
```

