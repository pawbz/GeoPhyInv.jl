
using Statistics

medium = Medium(:acou_homo1);
ageom = AGeom(medium.mgrid,:xwell);
tgrid = range(0.0,stop=2.0,length=1000)
wav = ricker(10.0, tgrid, tpeak=0.25, );
srcwav = SrcWav(tgrid, ageom, [:P])
update!(srcwav, [:P], wav)


vp0=mean(medium[:vp])
rho0=mean(medium[:rho])
rec1 = GeoPhyInv.Born.mod(vp0=vp0,
            medium_pert=medium,
		     rho0=rho0,
			 ageom=ageom, srcwav=srcwav, tgridmod=tgrid, src_flag=2)



pa=SeisForwExpt(Fdtd(),npw=1,medium=medium,
    ageom=[ageom], srcwav=[srcwav],
        sflags=[2], rflags=[1],
	    tgrid=tgrid, verbose=true );

@time update!(pa);


# least-squares misfit
paerr=GeoPhyInv.VNamedD_misfit(rec1, pa.c.data[1])
err=GeoPhyInv.func_grad!(paerr)

# normalize error
error = err[1]/paerr.ynorm

# desired accuracy?
@test error<1e-2
