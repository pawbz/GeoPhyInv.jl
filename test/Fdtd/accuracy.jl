model = J.Gallery.Seismic(:acou_homo2);
model0 = J.Gallery.Seismic(:acou_homo2);
J.Models.Seismic_addon!(model,randn_perc=1, fields=[:χvp,:χρ])
acqgeom=J.Acquisition.Geom_fixed(model,5,10)
acqsrc=J.Acquisition.Src_fixed_mod(acqgeom.nss,1,[:P],mod=model, nλ=3, tmaxfrac=1.0)
tgrid=acqsrc.tgrid

attrib_mod=JF.ModFdtdBorn()

#for born_flag in [true, false]
pa=JF.Param(acqsrc, acqgeom, tgrid, attrib_mod, model0, 
modm_obs=model,  
modm0=model0,
igrid_interp_scheme=:B2,    
igrid=Grid.M2D_resamp(model.mgrid, 300.,300.,),     
parameterization=[:χvp, :χρ, :null],   verbose=false,
nworker=1)


pa_parallel=JF.Param(acqsrc, acqgeom, 
	tgrid, attrib_mod, model0, 
	modm_obs=model,  
	modm0=model0,
	igrid_interp_scheme=:B2,    
	igrid=Grid.M2D_resamp(model.mgrid, 300.,300.,),     
	parameterization=[:χvp, :χρ, :null],   verbose=false,
	nworker=nothing)


paerr=JuMIT.Data.P_misfit(pa_parallel.paTD.y, pa.paTD.y)
err=JuMIT.Data.func_grad!(paerr)

# normalize error
error = err[1]/paerr.ynorm


	result=JF.xfwi!(pa, JF.Migr())
	result_parallel=JF.xfwi!(pa_parallel, JF.Migr())


f=Misfits.error_squared_euclidean!(nothing, result[2], result_parallel[2], nothing, norm_flag=true)


ffff

model = JuMIT.Gallery.Seismic(:acou_homo1);
acqgeom = JuMIT.Gallery.Geom(model.mgrid,:xwell);
tgrid = JuMIT.Gallery.M1D(:acou_homo1);
wav = Signals.Wavelets.ricker(10.0, tgrid, tpeak=0.25, );
# source wavelet for modelling
acqsrc = JuMIT.Acquisition.Src_fixed(acqgeom.nss,1,[:P],wav,tgrid);


vp0=mean(JuMIT.Models.χ(model.χvp,model.ref.vp,-1))
ρ0=mean(JuMIT.Models.χ(model.χρ,model.ref.ρ,-1))
rec1 = JuMIT.Analytic.mod(vp0=vp0,
            model_pert=model,
		     ρ0=ρ0,
			 acqgeom=acqgeom, acqsrc=acqsrc, tgridmod=tgrid, src_flag=2)



pa=JuMIT.Fdtd.Param(npw=1,model=model,
    acqgeom=[acqgeom], acqsrc=[acqsrc],
        sflags=[2], rflags=[1],
	    tgridmod=tgrid, verbose=true );

@btime JuMIT.Fdtd.mod!(pa);


# least-squares misfit
paerr=JuMIT.Data.P_misfit(rec1, pa.c.data[1])
err=JuMIT.Data.func_grad!(paerr)

# normalize error
error = err[1]/paerr.ynorm

# desired accuracy?
@test error<1e-2
