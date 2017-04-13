using SIT
model = SIT.Gallery.Seismic(:acou_homo1);
model = SIT.Models.Seismic_addon(model, circ_rad=30., circ_loc=[0.,0.],circ_pert=0.2)
model0 = SIT.Gallery.Seismic(:acou_homo1);

tgrid = SIT.Gallery.M1D(:acou_homo1);
acqgeom = SIT.Gallery.Geom(model.mgrid,:oneonev);
acqsrc = SIT.Gallery.Src(:acou_homo1);

# Store boundaries
rec1, b1, g1= SIT.Fdtd.fdtd_mod(model=model0, model0=model0, acqgeom=[acqgeom], acqsrc=[acqsrc], 
    tgridmod=tgrid, verbose=true, boundary_save_flag=true);

# Generate Born Data
recbh, bbb, gtemp = SIT.Fdtd.fdtd_mod(boundary_in=b1, npropwav=2, model=model, model0=model0, 
    acqgeom=[acqgeom,acqgeom], acqsrc=[acqsrc,acqsrc], src_flags=[0.0, 0.0], 
    tgridmod=tgrid, verbose=true, boundary_save_flag=false, born_flag=true);


#rec1 = SIT.Fdtd.fdtd_born_mod();

# migrate Born data
adjgeom = SIT.Acquisition.Geom_adj(acqgeom);
adjsrc = SIT.Acquisition.Src(adjgeom.nss,adjgeom.ns,1,recbh[2].d[end:-1:1,:,:,:], recbh[2].tgrid);
recb, bb, gb = SIT.Fdtd.fdtd_mod(npropwav=2, model=model, model0=model0, 
    acqgeom=[adjgeom,adjgeom], acqsrc=[adjsrc, adjsrc], src_flags=[0.0, -2.0], 
    tgridmod=tgrid, verbose=true, grad_out_flag=true, boundary_in=b1);


