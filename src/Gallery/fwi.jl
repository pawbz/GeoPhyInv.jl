
# create pizza problem

function xfwi_problem(attrib::Symbol)

	if(attrib==:pizza)

		model = Gallery.Seismic(:acou_homo2)
		print(model)

		# add perturbations
		Models.Seismic_addon!(model, ellip_rad=50., ellip_loc=[500.,0.],ellip_pert=0.1,randn_perc=0.1, fields=[:χvp,:χρ])
		Models.Seismic_addon!(model, ellip_rad=50., ellip_loc=[0.,500.],ellip_pert=0.1,randn_perc=0.1, fields=[:χvp,:χρ])
		Models.Seismic_addon!(model, ellip_rad=50., ellip_loc=[-500.,0.],ellip_pert=0.1,randn_perc=0.1, fields=[:χvp,:χρ])
		Models.Seismic_addon!(model, ellip_rad=50., ellip_loc=[0.,-500.],ellip_pert=0.1,randn_perc=0.1, fields=[:χvp,:χρ])

		# starting model
		model0 = Gallery.Seismic(:acou_homo2);

		# sources, receivers
		acqgeom=Acquisition.Geom_circ(nss=5,nr=20,rad=[900.,900.])

		acqsrc=Acquisition.Src_fixed_mod(acqgeom.nss,1,[:P],mod=model,nλ=3)
		tgrid=acqsrc.tgrid 

		pa = FWI.Param(acqsrc, acqgeom, 
		       tgrid, FWI.ModFdtd(),
		       model0,
		       modm_obs=model, 
		       igrid_interp_scheme=:B2, 
		       #igrid=model.mgrid, 
		       #igrid=Grid.M2D_resamp(model.mgrid, 50.,50.,),    
		       igrid=Grid.M2D_resamp(model.mgrid, 300.,300.,),    
		       parameterization=[:χvp, :χρ, :null],  verbose=false);
		println("ffffaf")

		return pa

	else
		error("invalid attrib")
	end
end


