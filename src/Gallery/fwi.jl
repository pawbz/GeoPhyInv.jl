
# create pizza problem

function SeisInvExpt(attrib::Symbol; born_flag=false, rfields=[:Vx,:Vz,:P], α=0.0, parameterization=[:χvp, :χρ, :null])

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
		igrid=broadcast(x->range(x[1],stop=x[end],step=50.),model.mgrid)
		igrid_interp_scheme=:B2
	elseif(attrib==:downhole)
		gri = [range(-120.,stop=30.,step=1.0), range(-30.,stop=30.,step=1.0)]
		model = Models.Seismic_zeros(gri)
		Models.adjust_bounds!(model, [2500., 3500.], [2500., 3500.], [2500.,3500.])
		model0=deepcopy(model)
		print(model)
		Models.Seismic_addon!(model, randn_perc=0.1, fields=[:χvp,:χρ])
		Models.Seismic_addon!(model0, randn_perc=0.1, fields=[:χvp,:χρ])
		Models.Seismic_addon!(model, ellip_rad=[2000., 5.], ellip_loc=[0.,0.],ellip_pert=0.1, ellip_α=α, fields=[:χvp],)

		acqgeom=Acquisition.Geom_fixed(-75., -50., 0.0, -100., -50., 0.0, 2, 20, :vertical, :vertical)
		acqsrc=Acquisition.Src_fixed_mod(acqgeom.nss,1,[:P],mod=model,nλ=8, tmaxfrac=0.8)
		tgrid=acqsrc.tgrid 
		igrid=[range(-40.,stop=30.,length=100), range(-30.,length=1,step=60.)]
		igrid_interp_scheme=:B2
	else
		error("invalid attrib")
	end

	if(born_flag)
		pa = FWI.Param(acqsrc, acqgeom, tgrid, FWI.ModFdtdBorn(),
		       rfields=rfields,
		       model0, modm0=model0, modm_obs=model, 
		       igrid_interp_scheme=igrid_interp_scheme, 
		       igrid=igrid, parameterization=parameterization, verbose=false);


	else
		pa = FWI.Param(acqsrc, acqgeom, tgrid, FWI.ModFdtd(),
		       model0, modm_obs=model, modm0=model0, 
		       rfields=rfields,
		       igrid_interp_scheme=igrid_interp_scheme, 
		       igrid=igrid, parameterization=parameterization, verbose=false);
	end

	return pa
end


