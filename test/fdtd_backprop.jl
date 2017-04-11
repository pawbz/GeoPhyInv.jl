using SIT
model = SIT.Gallery.Seismic(:acou_homo1);
model0 = SIT.Gallery.Seismic(:acou_homo1);
tgrid = SIT.Gallery.M1D(:acou_homo1);

#acqgeom = SIT.Gallery.Geom(:acou_homo1);
acqgeom = SIT.Gallery.Geom(model.mgrid,:onefiftyv);

#src_n =2;
acqsrc = SIT.Gallery.Src(:acou_homo1);

rec1, b1 = SIT.Fdtd.fdtd_mod(model=model, model0=model0, acqgeom=[acqgeom], acqsrc=[acqsrc], 
    tgridmod=tgrid, verbose=true, boundary_save_flag=true);

# just to make sure 
acqsrc.wav = acqsrc.wav * 2.;
rec2, b2 = SIT.Fdtd.fdtd_mod(boundary_in = b1, model=model, model0=model0, 
			     acqgeom=[acqgeom], acqsrc=[acqsrc], src_flags=[0.],
			     tgridmod=tgrid, verbose=true);

# compare results
err = norm(rec1[1].d[:,1,1]-rec2[1].d[end:-1:1,1,1]) / norm(rec1[1].d[:,1,1]);

println(string("error\t", err))




