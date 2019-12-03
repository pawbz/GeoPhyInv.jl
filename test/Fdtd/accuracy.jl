
using Statistics

model = GeoPhyInv.Gallery.Seismic(:acou_homo1);
acqgeom = GeoPhyInv.Gallery.Geom(model.mgrid,:xwell);
tgrid = range(0.0,stop=2.0,length=1000)
wav = GeoPhyInv.Utils.Wavelets.ricker(10.0, tgrid, tpeak=0.25, );
# source wavelet for modelling
acqsrc = GeoPhyInv.Acquisition.Src_fixed(acqgeom.nss,1,[:P],wav,tgrid);


vp0=mean(model[:vp])
rho0=mean(model[:rho])
rec1 = GeoPhyInv.Born.mod(vp0=vp0,
            model_pert=model,
		     rho0=rho0,
			 acqgeom=acqgeom, acqsrc=acqsrc, tgridmod=tgrid, src_flag=2)



pa=GeoPhyInv.Fdtd.Param(npw=1,model=model,
    acqgeom=[acqgeom], acqsrc=[acqsrc],
        sflags=[2], rflags=[1],
	    tgridmod=tgrid, verbose=true );

@time GeoPhyInv.Fdtd.mod!(pa);


# least-squares misfit
paerr=GeoPhyInv.Data.P_misfit(rec1, pa.c.data[1])
err=GeoPhyInv.Data.func_grad!(paerr)

# normalize error
error = err[1]/paerr.ynorm

# desired accuracy?
@test error<1e-2
