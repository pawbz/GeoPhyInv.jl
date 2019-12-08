
"""
Gallery of `SeisInvExpt`.
# Arguments
* `attrib::Symbol` : name 
  * `=:pizza` all around sources and receivers
  * `=:downhole` sources and receivers along a drill-bit
"""
function SeisInvExpt(attrib::Symbol; born_flag=false, rfields=[:Vx,:Vz,:P], α=0.0, parameterization=[:χvp, :χrho, :null])

	if(attrib==:pizza)
		# starting model
		model0 = Gallery.Seismic(:acou_homo2);
		model = deepcopy(model0)
		# add some noise to starting model
		update!(model0, [:vp,:rho], randn_perc=0.1)

		print(model)

		# add perturbations
		for ellip_loc in [[500.,0.], [0.,500.], [-500.,0.], [0.,-500.]]
			update!(model, [:vp,:rho], ellip_rad=50., ellip_loc=ellip_loc, 
				ellip_pert=100.)
		end
		update!(model, [:vp,:rho], randn_perc=0.1)

		# sources, receivers
		geom=Geom_circ(nss=5,nr=20,rad=[900.,900.])

		srcwav=Src_fixed_mod(geom.nss,1,[:P],mod=model,nλ=3)
		tgrid=srcwav.tgrid 
		igrid=broadcast(x->range(x[1],stop=x[end],step=50.),model.mgrid)
		igrid_interp_scheme=:B2
	elseif(attrib==:downhole)
		gri = [range(-120.,stop=30.,step=1.0), range(-30.,stop=30.,step=1.0)]
		model = Medium(gri)
		update!(model, [:vp,:rho], [[2500., 3500.], [2500., 3500.]])
		fill!(model)
		model0=deepcopy(model)
		print(model)
		update!(model, [:vp,:rho], randn_perc=0.1)
		update!(model0, [:vp,:rho], randn_perc=0.1)
		update!(model, [:vp], ellip_rad=[2000., 5.], ellip_loc=[0.,0.],ellip_pert=0.1, ellip_α=α)

		geom=Geom_fixed(-75., -50., 0.0, -100., -50., 0.0, 2, 20, :vertical, :vertical)
		srcwav=Src_fixed_mod(geom.nss,1,[:P],mod=model,nλ=8, tmaxfrac=0.8)
		tgrid=srcwav.tgrid 
		igrid=[range(-40.,stop=30.,length=100), range(-30.,length=1,step=60.)]
		igrid_interp_scheme=:B2
	else
		error("invalid attrib")
	end

	if(born_flag)
		pa = FWI.Param(srcwav, geom, tgrid, FWI.ModFdtdBorn(),
		       rfields=rfields,
		       model0, modm0=model0, modm_obs=model, 
		       igrid_interp_scheme=igrid_interp_scheme, 
		       igrid=igrid, parameterization=parameterization, verbose=false);


	else
		pa = FWI.Param(srcwav, geom, tgrid, FWI.ModFdtd(),
		       model0, modm_obs=model, modm0=model0, 
		       rfields=rfields,
		       igrid_interp_scheme=igrid_interp_scheme, 
		       igrid=igrid, parameterization=parameterization, verbose=false);
	end

	return pa
end


