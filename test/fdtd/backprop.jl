
medium=Medium(:acou_homo2D,5)
update!(medium, [:vp,:rho], randn_perc=0.1)
ageom=AGeom(medium.mgrid, SSrcs(4), Srcs(1), Recs(100))
update!(ageom, SSrcs(),[0,0],990.0,[0, 2π])
update!(ageom, Recs(),[0,0],990.0,[0, 2π])

wav, tgrid=ricker(medium, 3, 0.8)

for field in [:p, :vx,:vz]
	println("############ Testing Backprop for source type ", field)
	srcwav = SrcWav(tgrid, ageom, [field])
	update!(srcwav, [field], wav)

	for sflags in [[1,-1],[2,-2]]
		pa=SeisForwExpt(FdtdAcou(),npw=1, tgrid=tgrid,
		#	pml_edges=[:null],
			gmodel_flag=false,
			sflags=[sflags[1]],
			# snaps_field=true,
			verbose=true,
			backprop_flag=1,
			illum_flag=true,ageom=[ageom], srcwav=[srcwav],
			medium=medium);

		update!(pa);
		rec1=deepcopy(pa.c.data[1])

		# change source flag and update wavelets in pa
		pa.c.sflags=[sflags[2]];
		GeoPhyInv.update!(pa,[srcwav])
		pa.c.backprop_flag=-1 # do backpropagation

		update!(pa)
		rec2=deepcopy(pa.c.data[1])

		# time reverse
		reverse!(rec2);

		# compare results
		# least-squares misfit
		paerr=GeoPhyInv.VNamedD_misfit(rec1, rec2)
		err = GeoPhyInv.func_grad!(paerr)

		# normalized error
		error = err[1]./paerr.ynorm

		# desired accuracy?
		@test error<1e-20
	end
end
