using SIT
using Base.Test


model = SIT.Gallery.Seismic(:acou_homo1);
tgrid = SIT.Gallery.M1D(:acou_homo1);

#acqgeom = SIT.Gallery.Geom(:acou_homo1);
acqgeom = SIT.Gallery.Geom(model.mgrid,:onefiftyv);

#src_n =2;
acqsrc = SIT.Gallery.Src(:acou_homo1);

rec1 = SIT.Fdtd.mod(model=model,  acqgeom=[acqgeom], acqsrc=[acqsrc], 
			tgridmod=tgrid, verbose=true, boundary_save_flag=true, src_flags=[2.]);

# sink source for propagation from boundary 
acqsrcsink = deepcopy(acqsrc); 
acqsrcsink.wav = [-1.0.*flipdim(acqsrc.wav[i,j],1) for i in 1:acqsrc.nss, j in 1:acqsrc.nfield];

rec2 = SIT.Fdtd.mod(boundary_in = rec1[2], model=model, 
			     acqgeom=[acqgeom], acqsrc=[acqsrcsink], src_flags=[2.],
			     tgridmod=tgrid, verbose=true);

# compare results

# time reverse
SIT.Data.TD_tr!(rec2[1][1]);
# least-squares misfit
err = SIT.Misfits.TD(rec1[1][1], rec2[1][1])

# normalization
error = err[1]/SIT.Data.TD_dot(rec1[1][1], rec1[1][1])

# desired accuracy? 
@test error<1e-5




