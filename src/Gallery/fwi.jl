
# create pizza problem

function xfwi_problem(attrib::Symbol; born_flag=false)

	if(attrib==:pizza)
		# starting model
		model0 = Gallery.Seismic(:acou_homo2);
		model = deepcopy(model0)
		# add some noise to starting model
		Models.Seismic_addon!(model0, randn_perc=0.1, fields=[:χvp,:χρ])

		print(model)

		# add perturbations
		for ellip_loc in [[500.,0.], [0.,500.], [-500.,0.], [0.,-500.]]
			Models.Seismic_addon!(model, ellip_rad=50., ellip_loc=ellip_loc, 
				ellip_pert=0.1, fields=[:χvp,:χρ])
		end
		Models.Seismic_addon!(model, randn_perc=0.1, fields=[:χvp,:χρ])

		# sources, receivers
		acqgeom=Acquisition.Geom_circ(nss=5,nr=20,rad=[900.,900.])

		acqsrc=Acquisition.Src_fixed_mod(acqgeom.nss,1,[:P],mod=model,nλ=3)
		tgrid=acqsrc.tgrid 
		igrid=Grid.M2D_resamp(model.mgrid, 50.,50.,)
		parameterization=[:χvp, :χρ, :null]
		igrid_interp_scheme=:B2
	else
		error("invalid attrib")
	end

	if(born_flag)
		pa = FWI.Param(acqsrc, acqgeom, tgrid, FWI.ModFdtdBorn(),
		       model0, modm0=model0, modm_obs=model, 
		       igrid_interp_scheme=igrid_interp_scheme, 
		       igrid=igrid, parameterization=parameterization, verbose=false);


	else
		pa = FWI.Param(acqsrc, acqgeom, tgrid, FWI.ModFdtd(),
		       model0, modm_obs=model, modm0=model0, 
		       igrid_interp_scheme=igrid_interp_scheme, 
		       igrid=igrid, parameterization=parameterization, verbose=false);
	end

	return pa, model
end


