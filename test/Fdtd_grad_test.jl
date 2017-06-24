using SIT
using Base.Test

model_pert = SIT.Gallery.Seismic(:acou_homo2);
SIT.Models.Seismic_addon!(model_pert, circ_rad=30., circ_loc=[0.,0.],circ_pert=0.2, fields=[:χvp])
model0 = SIT.Gallery.Seismic(:acou_homo2);

tgrid = SIT.Gallery.M1D(:acou_homo2);
acqgeom = SIT.Gallery.Geom(model0.mgrid,:oneonev);
acqsrc = SIT.Gallery.Src(:acou_homo2);


