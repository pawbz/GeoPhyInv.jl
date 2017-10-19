using JuMIT
using Base.Test


model = JuMIT.Gallery.Seismic(:acou_homo1);
acqgeom = JuMIT.Gallery.Geom(model.mgrid,:xwell);
tgrid = JuMIT.Gallery.M1D(:acou_homo1);
wav = JuMIT.Wavelets.ricker(fqdom=10.0, tgrid=tgrid, tpeak=0.25, );
# source wavelet for born modelling
acqsrc = JuMIT.Acquisition.Src_fixed(acqgeom.nss,1,1,wav,tgrid);


vp0=mean(JuMIT.Models.χ(model.χvp,model.vp0,-1))
ρ0=mean(JuMIT.Models.χ(model.χρ,model.ρ0,-1))
rec1 = JuMIT.Analytic.mod(vp0=vp0,
		             ρ0=ρ0,
			         acqgeom=acqgeom, acqsrc=acqsrc, tgridmod=tgrid, src_flag=2)



rec2, bb, gg= JuMIT.Fdtd.mod(npropwav=1,model=model,
    acqgeom=[acqgeom], acqsrc=[acqsrc],
        src_flags=[2], recv_flags=[1],
	    tgridmod=tgrid, verbose=true );


# least-squares misfit
err = JuMIT.Misfits.TD(rec1, rec2[1])

# normalization
error = err[1]/JuMIT.Data.TD_dot(rec1, rec1)

# desired accuracy?
@test error<1e-2
