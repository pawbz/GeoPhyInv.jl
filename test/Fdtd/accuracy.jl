
using Statistics

model = Medium(:acou_homo1);
acqgeom = Geom(model.mgrid,:xwell);
tgrid = range(0.0,stop=2.0,length=1000)
wav = ricker(10.0, tgrid, tpeak=0.25, );
acqsrc = SrcWav(tgrid, acqgeom, [:P])
update!(acqsrc, [:P], wav)


vp0=mean(model[:vp])
rho0=mean(model[:rho])
rec1 = GeoPhyInv.Born.mod(vp0=vp0,
            model_pert=model,
		     rho0=rho0,
			 acqgeom=acqgeom, acqsrc=acqsrc, tgridmod=tgrid, src_flag=2)



pa=SeisForwExpt(npw=1,model=model,
    acqgeom=[acqgeom], acqsrc=[acqsrc],
        sflags=[2], rflags=[1],
	    tgridmod=tgrid, verbose=true );

@time update!(pa);


# least-squares misfit
paerr=GeoPhyInv.VNamedD_misfit(rec1, pa.c.data[1])
err=GeoPhyInv.func_grad!(paerr)

# normalize error
error = err[1]/paerr.ynorm

# desired accuracy?
@test error<1e-2
