```@meta
EditURL = "<unknown>/ageom/gallery.jl"
```

### Load packages

```@example gallery
using GeoPhyInv
```

### A surface acquisition

```@example gallery
mgrid=repeat(range(-1000.0,stop=1000.0,length=201), 2) # test spatial grid
#ageom_surf=GeoPhyInv.Gallery.AGeom(mgrid,:surf,nss=5,nr=30);
nothing #hide
```

