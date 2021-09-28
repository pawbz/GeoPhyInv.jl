```@meta
EditURL = "<unknown>/GeoPhyInv/test/fdtd/doc.jl"
```

````@example doc
using BenchmarkTools
using GeoPhyInv
using Test
````

# Intro

```@docs
SeisForwExpt
Base.getindex(::GeoPhyInv.PFdtd, ::Symbol, ::Int)
```

# Methods
```@docs
update!(::GeoPhyInv.PFdtd)
update!(::GeoPhyInv.PFdtd, ::Medium)
update!(::GeoPhyInv.PFdtd, ::SrcWav, ::Any)
```

