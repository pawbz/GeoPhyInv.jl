
model = J.Gallery.Seismic(:acou_homo2);
J.Models.Seismic_addon!(model,randn_perc=1, fields=[:χvp,:χρ])

model0 = J.Gallery.Seismic(:acou_homo2);
J.Models.Seismic_addon!(model0, randn_perc=1, fields=[:χvp,:χρ])

acqgeom=J.Acquisition.Geom_fixed(model,1,10)
acqsrc=J.Acquisition.Src_fixed_mod(acqgeom.nss,1,[:P],mod=model, nλ=3, tmaxfrac=1.0)
tgrid=acqsrc.tgrid

parameterization=[:χvp, :χρ, :null]

@testset "test parallel implementation during gradient" begin
	for attrib_mod in [JF.ModFdtd(), JF.ModFdtdBorn()]
		pa=JF.Param(acqsrc, acqgeom, tgrid, attrib_mod, model0, 
				     modm_obs=model,  
				     modm0=model0,
				     igrid_interp_scheme=:B2,    
				     igrid=Grid.M2D_resamp(model.mgrid, 300.,300.,),     
				     parameterization=parameterization,   verbose=false,
				     nworker=1)


		pa_parallel=JF.Param(acqsrc, acqgeom, tgrid, attrib_mod, model0, 
				     modm_obs=model,  
				     modm0=model0,
				     igrid_interp_scheme=:B2,    
				     igrid=Grid.M2D_resamp(model.mgrid, 300.,300.,),     
				     parameterization=parameterization,   verbose=false,
				     nworker=nothing)

		result=JF.xfwi!(pa, JF.Migr(), attrib_mod)

		result_parallel=JF.xfwi!(pa_parallel, JF.Migr(), attrib_mod)

		@test result[2] ≈ result_parallel[2]
	end
end

@testset "Testing Born Modeling and its gradient" begin

	pa=JF.Param(acqsrc, acqgeom, tgrid, JF.ModFdtdBorn(), model0, 
	     		     modm0=model0,
			     modm_obs=model,  
			     igrid_interp_scheme=:B2,    
			     igrid=Grid.M2D_resamp(model.mgrid, 350.,350.,),     
			     parameterization=parameterization,   verbose=false)


	JF.xfwi!(pa, JF.LS(), JF.ModFdtdBorn(),  bounded_flag=true, solver=:ipopt,
			ipopt_options=[["max_iter", 0],["derivative_test", "first-order"]])

	result=JF.xfwi!(pa, JF.Migr(), JF.ModFdtdBorn())

	pa_fd=deepcopy(pa);
	result_fd=JF.xfwi!(pa_fd, JF.Migr_fd(), JF.ModFdtdBorn())

	f=Misfits.error_squared_euclidean!(nothing, result[2], result_fd[2], nothing, norm_flag=true)

	@test f<1e-15
end

@testset "Testing gradient LS FWI" begin
	pa=JF.Param(acqsrc, acqgeom, tgrid, JF.ModFdtd(), model0, 
			     modm_obs=model,  
			     igrid_interp_scheme=:B2,    
			     igrid=Grid.M2D_resamp(model.mgrid, 350.,350.,),     
			     parameterization=parameterization,   verbose=false)

	JF.xfwi!(pa, JF.LS(), JF.ModFdtd(),  bounded_flag=true, solver=:ipopt,
			ipopt_options=[["max_iter", 0],["derivative_test", "first-order"]])


	result=JF.xfwi!(pa, JF.Migr(), JF.ModFdtd())

	pa_fd=deepcopy(pa);
	result_fd=JF.xfwi!(pa_fd, JF.Migr_fd(), JF.ModFdtd())

	f=Misfits.error_squared_euclidean!(nothing, result[2], result_fd[2], nothing, norm_flag=true)

	@test f<1e-15
end


