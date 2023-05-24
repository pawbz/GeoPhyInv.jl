



# function SeisInvExpt(attrib_mod::Union{FdtdAcoustic, FdtdAcoustic{Born}},attrib_inv::Union{LS,LS_prior,Migr,Migr_FD}, attrib::Symbol; 
# 		     rfields=[:vx,:vz,:p], α=0.0, parameterization=[:χvp, :χrho, :null])

# 	if(attrib==:pizza)
# 		# starting model
# 		model = Medium(:pizza);
# 		model0 = similar(model)

# 		ageom=AGeom(model.grid, SSrcs(5), Srcs(1), Recs(20))
# 		update!(ageom, SSrcs(),[0,0],900.0,[0, 2π])
# 		update!(ageom, Recs(),[0,0],900.0,[0, 2π])


# 		wav, tgrid=ricker(model, 3, 1.0)
# 		srcwav = SrcWav(tgrid, ageom, [:p])
# 		update!(srcwav, [:p], wav)

# 		igrid=broadcast(x->range(x[1],stop=x[end],step=50.),model.grid)
# 		igrid_interp_scheme=:B2
# 	elseif(attrib==:downhole)
# 		gri = [range(-120.,stop=30.,step=1.0), range(-30.,stop=30.,step=1.0)]
# 		model = Medium(gri)
# 		update!(model, [:vp,:rho], [[2500., 3500.], [2500., 3500.]])
# 		fill!(model)
# 		model0=deepcopy(model)
# 		print(model)
# 		update!(model, [:vp,:rho], randn_perc=0.1)
# 		update!(model0, [:vp,:rho], randn_perc=0.1)
# 		update!(model, [:vp], ellip_rad=[2000., 5.], ellip_loc=[0.,0.],ellip_pert=0.1, ellip_α=α)


# 		ageom=AGeom(model.grid, SSrcs(2), Srcs(1), Recs(20))

# 		update!(ageom, SSrcs(),[-75,0],[-50,0])
# 		update!(ageom, Recs(),[-100,0],[-50,0])

# 		wav, tgrid=ricker(model, 8, 0.8)
# 		srcwav = SrcWav(tgrid, ageom, [:p])
# 		update!(srcwav, [:p], wav)

# 		igrid=[range(-40.,stop=30.,length=100), range(-30.,length=1,step=60.)]
# 		igrid_interp_scheme=:B2
# 	else
# 		error("invalid attrib")
# 	end

# 	return  SeisInvExpt(attrib_mod, attrib_inv, srcwav=srcwav, ageom=ageom, tgrid=tgrid,
# 		       rfields=rfields,
# 		       mediumm=model0, mediumm0=model0, mediumm_obs=model, 
# 		       igrid_interp_scheme=igrid_interp_scheme, 
# 		       igrid=igrid, parameterization=parameterization, verbose=false);
# end


