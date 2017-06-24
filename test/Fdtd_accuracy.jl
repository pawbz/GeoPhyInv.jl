using SIT
using Base.Test


model = SIT.Gallery.Seismic(:acou_homo1);
acqgeom = SIT.Gallery.Geom(model.mgrid,:oneonev);
tgrid = SIT.Gallery.M1D(:acou_homo1);
wav = SIT.Wavelets.ricker(fqdom=10.0, tgrid=tgrid, tpeak=0.25, );
# source wavelet for born modelling
acqsrc = SIT.Acquisition.Src_fixed(acqgeom.nss,1,1,wav,tgrid);


rec1 = SIT.Analytic.mod(vp0=model.vp0,
		             ρ0=model.ρ0,
			         acqgeom=acqgeom, acqsrc=acqsrc, tgridmod=tgrid, src_flag=2)



rec2, bb, gg= SIT.Fdtd.mod(npropwav=1,model=model, 
    acqgeom=[acqgeom], acqsrc=[acqsrc], 
        src_flags=[2], recv_flags=[1],
	    tgridmod=tgrid, verbose=true );


# least-squares misfit
err = SIT.Misfits.TD(rec1, rec2[1])

# normalization
error = err[1]/SIT.Data.TD_dot(rec1, rec1)

# desired accuracy? 
@test error<1e-2


