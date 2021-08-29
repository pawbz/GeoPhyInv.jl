using GeoPhyInv
using Test

model = Medium(:acou_homo2);
update!(model, [:vp,:rho], randn_perc=1)

model0 = Medium(:acou_homo2);
update!(model0, [:vp,:rho], randn_perc=1)

ageom=AGeom(model.mgrid, :surf)
wav, tgrid=ricker(model, 3, 1.0)
srcwav=SrcWav(tgrid,ageom,[:P])
update!(srcwav,[:P], wav)

parameterization=[:χvp, :χrho, :null]

mgrid=model.mgrid

@testset "test parallel implementation during gradient" begin
	for attrib_mod in [FdtdAcou(), FdtdAcouBorn()]
		global pa=SeisInvExpt(attrib_mod, Migr(), srcwav=srcwav, ageom=ageom, tgrid=tgrid, modm=model0,
				     modm_obs=model,
				     modm0=model0,
				     igrid_interp_scheme=:B2,
				     igrid=broadcast(x->range(x[1],stop=x[end],step=300.),mgrid),
				     parameterization=parameterization,   verbose=false,
				     nworker=1)


		global pa_parallel=SeisInvExpt(attrib_mod, Migr(), srcwav=srcwav, ageom=ageom, tgrid=tgrid, modm=model0,
				     modm_obs=model,
				     modm0=model0,
				     igrid_interp_scheme=:B2,
				     igrid=broadcast(x->range(x[1],stop=x[end],step=300.),mgrid),
				     parameterization=parameterization,   verbose=false,
				     nworker=nothing)

		result=update!(pa)

		result_parallel=update!(pa_parallel)

		@test result[2] ≈ result_parallel[2]
	end
end



@testset "Testing Born Modeling and its gradient" begin
	global expt=x->SeisInvExpt(FdtdAcouBorn(), x,
			     srcwav=srcwav, ageom=ageom, tgrid=tgrid, modm=model0,
	     		     modm0=model0,
			     modm_obs=model,
			     igrid_interp_scheme=:B2,
			     igrid=broadcast(x->range(x[1],stop=x[end],step=350.),mgrid),
			     parameterization=parameterization,   verbose=false)

	global pa=expt(LS())


	update!(pa, bounded_flag=true, solver=:ipopt,
			ipopt_options=[["max_iter", 0],["derivative_test", "first-order"]])
	pa=expt(Migr())

	result=update!(pa)

	pa_fd=expt(Migr_FD());
	result_fd=update!(pa_fd)

	f=Misfits.error_squared_euclidean!(nothing, result[2], result_fd[2], nothing, norm_flag=true)

	@test f<1e-15
end


@testset "Testing gradient LS FWI" begin
	expt=x->SeisInvExpt(FdtdAcou(), x,
			     srcwav=srcwav, ageom=ageom, tgrid=tgrid, modm=model0,
	     		     modm0=model0,
			     modm_obs=model,
			     igrid_interp_scheme=:B2,
			     igrid=broadcast(x->range(x[1],stop=x[end],step=350.),mgrid),
			     parameterization=parameterization,   verbose=false)
	pa=expt(LS())


	update!(pa, bounded_flag=true, solver=:ipopt,
			ipopt_options=[["max_iter", 0],["derivative_test", "first-order"]])
	pa=expt(Migr())

	result=update!(pa)

	pa_fd=expt(Migr_FD());
	result_fd=update!(pa_fd)



	f=Misfits.error_squared_euclidean!(nothing, result[2], result_fd[2], nothing, norm_flag=true)

	@test f<1e-15
end

