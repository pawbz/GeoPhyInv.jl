# load packages

using GeoPhyInv
using Statistics

# create simple (almost) homogeneous acoustic model

model=J.Gallery.Seismic(:acou_homo1)
J.Models.Seismic_addon!(model, randn_perc=0.01)

# a simple acquisition geometry

acqgeom = GeoPhyInv.Gallery.Geom(model.mgrid,:xwell);


