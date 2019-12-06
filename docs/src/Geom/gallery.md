```@meta
EditURL = "@__REPO_ROOT_URL__/"
```

# The Geom Datatype

An acquisition geometry, bundled into the mutable `Geom` type,
has to be specified in order to
create an `Expt` variable.
Acquisiton has supersources, sources and receivers.
Each supersource has `ns` multiple sources that are
injected (or active) simultaneously.
For each supersource,
a set of `nr` receivers are
recording waves.
This package provides tools to easily create commonly-used `Geom` variables,
after deciding on a spatial grid `mgrid` for the experiment.

### Load packages

```@example gallery
using GeoPhyInv
using Plots
```

### A surface acquisition

```@example gallery
mgrid=repeat(range(-1000.0,stop=1000.0,length=201), 2) # test spatial grid
geom_surf=GeoPhyInv.Gallery.Geom(mgrid,:surf,nss=5,nr=30);
GIPh.Plots.geom!(geom_surf)
```

