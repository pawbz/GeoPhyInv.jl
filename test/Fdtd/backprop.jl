
model=GeoPhyInv.Gallery.Seismic(:acou_homo1)
update!(model, [:vp,:rho], randn_perc=0.1)
acqgeom =GeoPhyInv.Acquisition.Geom_circ(nss=4,nr=100,rad=[990.,990.]);
acqgeom =GeoPhyInv.Acquisition.Geom_circ(nss=4,nr=100,rad=[0.,200.]);
acqsrc=GeoPhyInv.Acquisition.Src_fixed_mod(acqgeom.nss,1,[:P],mod=model, nλ=3, tmaxfrac=0.4)

for sflags in [[1,-1],[2,-2]]
	pa=GeoPhyInv.Fdtd.Param(born_flag=false,npw=1, tgridmod=acqsrc.tgrid,
	#	abs_trbl=[:null],
		gmodel_flag=false,
		sflags=[sflags[1]],
		snaps_flag=true,
		verbose=true,
		backprop_flag=1,
		illum_flag=true,acqgeom=[acqgeom], acqsrc=[acqsrc],
		model=model);

	GeoPhyInv.Fdtd.mod!(pa);
	rec1=deepcopy(pa.c.data[1])

	# change source flag and update wavelets in pa
	pa.c.sflags=[sflags[2]];
	GeoPhyInv.Fdtd.update_acqsrc!(pa,[acqsrc])
	pa.c.backprop_flag=-1 # do backpropagation

	GeoPhyInv.Fdtd.mod!(pa);
	rec2=deepcopy(pa.c.data[1])

	# time reverse
	GeoPhyInv.Data.TD_tr!(rec2);

	# compare results
	# least-squares misfit
	paerr=GeoPhyInv.Data.P_misfit(rec1, rec2)
	err = GeoPhyInv.Data.func_grad!(paerr)

	# normalized error
	error = err[1]./paerr.ynorm

	# desired accuracy?
	@test error<1e-20
end
