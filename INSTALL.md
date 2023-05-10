In order to install `GeoPhyInv` enter these package manager commands in the REPL.
```julia
using Pkg
Pkg.add(PackageSpec(name="GeoPhyInv",url="https://github.com/pawbz/GeoPhyInv.jl.git"))
```
It is necessary to configure GeoPhyInv with a macro `@init_parallel_stencil` before using it. If you need to change this configuration, the julia kernel must be restarted.
```julia
using GeoPhyInv; @init_parallel_stencil(⋯)
```