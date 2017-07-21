#addprocs(10)
using SIT
using Base.Test
reload("SIT")


"""
for single source
"""
model = SIT.Gallery.Seismic(:acou_homo1);
tgrid = SIT.Gallery.M1D(:acou_homo1_long);

#acqgeom = SIT.Gallery.Geom(:acou_homo1);
acqgeom = SIT.Gallery.Geom(model.mgrid,:surf);

#src_n =2;
acqsrc=SIT.Acquisition.Src_fixed_mod(acqgeom.nss,1,1,model,3)

@time rec1, b1,  gg = SIT.Fdtd.mod!(model=model,  acqgeom=[acqgeom], acqsrc=[acqsrc], 
		     backprop_flag=1,
		    tgridmod=tgrid, verbose=true, src_flags=[2], recv_flags=[2]);

@time rec2, b2,  gg = SIT.Fdtd.mod!( model=model, boundary=b1,
		     backprop_flag=-1,
		    acqgeom=[acqgeom], acqsrc=[acqsrc], src_flags=[3], recv_flags=[1],
			     tgridmod=tgrid, verbose=true);

# compare results

# time reverse
SIT.Data.TD_tr!(rec2[1]);
# least-squares misfit
err = SIT.Misfits.TD!(rec1[1], rec2[1])

# normalization
error = err[1]/SIT.Data.TD_dot(rec1[1], rec1[1])

println("error:\t", error)

# desired accuracy? 
@test error<1e-10

exit()

"""
for multiple source in parallel
"""
model = SIT.Gallery.Seismic(:acou_homo1);
tgrid = SIT.Gallery.M1D(:acou_homo1_long);

#acqgeom = SIT.Gallery.Geom(:acou_homo1);
acqgeom = SIT.Gallery.Geom(model.mgrid,:tentenv);

#src_n =2;
acqsrc = SIT.Gallery.Src(:acou_homo1, 10);

alloc = SIT.Fdtd.mod_alloc(model=model, acqgeom=[acqgeom], npropwav=1, tgridmod=tgrid)

rec1 = SIT.Fdtd.mod(model=model,  acqgeom=[acqgeom], acqsrc=[acqsrc], 
		    tgridmod=tgrid, verbose=true, boundary_save_flag=true, src_flags=[1], recv_flags=[1]);

rec2 = SIT.Fdtd.mod(boundary_in = rec1[2], model=model, 
		    acqgeom=[acqgeom], acqsrc=[acqsrc], src_flags=[3], recv_flags=[1],
			     tgridmod=tgrid, verbose=true);

# compare results

# time reverse
SIT.Data.TD_tr!(rec2[1][1]);
# least-squares misfit
err = SIT.Misfits.TD(rec1[1][1], rec2[1][1])

# normalization
error = err[1]/SIT.Data.TD_dot(rec1[1][1], rec1[1][1])

# desired accuracy? 
@test error<1e-10

