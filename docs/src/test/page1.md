```@meta
EditURL = "https://github.com/TRAVIS_REPO_SLUG/blob/master/"
```

load packages

```@example page1
using GeoPhyInv
using Statistics
```

create simple (almost) homogeneous acoustic model

```@example page1
model=J.Gallery.Seismic(:acou_homo1)
J.Models.Seismic_addon!(model, randn_perc=0.01)
```

a simple acquisition geometry

```@example page1
acqgeom = GeoPhyInv.Gallery.Geom(model.mgrid,:xwell);
```

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

