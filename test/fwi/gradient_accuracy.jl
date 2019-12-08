using GeoPhyInv
using Test

# Can perform inversion of synthetic scenarios.
# First, the seismic data are modeled as in the forward problem. Then the
# data are used to perform full waveform inversion (FWI). The inverse
# problem estimates
# the Earth models and the source and receiver filters
# that resulted from the data.
# This task is necessary to test the performance of the inversion algorithm
# in various geological scenarios using different acquisition parameters.

model = Medium(:acou_homo2);
update!(model, [:vp,:rho], randn_perc=1)

model0 = Medium(:acou_homo2);
update!(model0, [:vp,:rho], randn_perc=1)

geom=Geom(model.mgrid, :surf)
wav, tgrid=ricker(model, 3, 1.0)
srcwav=SrcWav(tgrid,geom,[:P])
update!(srcwav,[:P], wav)

parameterization=[:χvp, :χrho, :null]

mgrid=model.mgrid

@testset "test parallel implementation during gradient" begin
	for attrib_mod in [GeoPhyInv.ModFdtd(), GeoPhyInv.ModFdtdBorn()]
		pa=SeisInvExpt(srcwav, geom, tgrid, attrib_mod, model0, 
				     modm_obs=model,  
				     modm0=model0,
				     igrid_interp_scheme=:B2,    
				     igrid=broadcast(x->range(x[1],stop=x[end],step=300.),mgrid),
				     parameterization=parameterization,   verbose=false,
				     nworker=1)


		pa_parallel=SeisInvExpt(srcwav, geom, tgrid, attrib_mod, model0, 
				     modm_obs=model,  
				     modm0=model0,
				     igrid_interp_scheme=:B2,    
				     igrid=broadcast(x->range(x[1],stop=x[end],step=300.),mgrid),
				     parameterization=parameterization,   verbose=false,
				     nworker=nothing)

		result=GeoPhyInv.migr!(pa, attrib_mod)

		result_parallel=GeoPhyInv.migr!(pa_parallel, attrib_mod)

		@test result[2] ≈ result_parallel[2]
	end
end


wefwefwef

@testset "Testing Born Modeling and its gradient" begin

	pa=JF.Param(srcwav, geom, tgrid, JF.ModFdtdBorn(), model0, 
	     		     modm0=model0,
			     modm_obs=model,  
			     igrid_interp_scheme=:B2,    
			     igrid=broadcast(x->range(x[1],stop=x[end],step=350.),mgrid),
			     parameterization=parameterization,   verbose=false)


	JF.fit!(pa, JF.LS(), JF.ModFdtdBorn(),  bounded_flag=true, solver=:ipopt,
			ipopt_options=[["max_iter", 0],["derivative_test", "first-order"]])

	result=JF.migr!(pa, JF.ModFdtdBorn())

	pa_fd=deepcopy(pa);
	result_fd=JF.migr_fd!(pa_fd,  JF.ModFdtdBorn())

	f=Misfits.error_squared_euclidean!(nothing, result[2], result_fd[2], nothing, norm_flag=true)

	@test f<1e-15
end

@testset "Testing gradient LS FWI" begin
	pa=JF.Param(srcwav, geom, tgrid, JF.ModFdtd(), model0, 
			     modm_obs=model,  
			     igrid_interp_scheme=:B2,    
			     igrid=broadcast(x->range(x[1],stop=x[end],step=350.),mgrid),
			     parameterization=parameterization,   verbose=false)

	JF.fit!(pa, JF.LS(), JF.ModFdtd(),  bounded_flag=true, solver=:ipopt,
			ipopt_options=[["max_iter", 0],["derivative_test", "first-order"]])


	result=JF.migr!(pa, JF.ModFdtd())

	pa_fd=deepcopy(pa);
	result_fd=JF.migr_fd!(pa_fd, JF.ModFdtd())

	f=Misfits.error_squared_euclidean!(nothing, result[2], result_fd[2], nothing, norm_flag=true)

	@test f<1e-15
end


