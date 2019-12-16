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

medium = Medium(:acou_homo2);
update!(medium, [:vp,:rho], randn_perc=1)

medium0 = Medium(:acou_homo2);
update!(medium0, [:vp,:rho], randn_perc=1)

ageom=AGeom(medium.mgrid, :surf)
wav, tgrid=ricker(medium, 3, 1.0)
srcwav=SrcWav(tgrid,ageom,[:P])
update!(srcwav,[:P], wav)

parameterization=[:χvp, :χrho, :null]

mgrid=medium.mgrid

@testset "test parallel implementation during gradient" begin
	for attrib_mod in [Fdtd(), FdtdBorn()]
		global pa=SeisInvExpt(attrib_mod, Migr(), srcwav=srcwav, ageom=ageom, tgrid=tgrid, mediumm=medium0, 
				     mediumm_obs=medium,  
				     mediumm0=medium0,
				     igrid_interp_scheme=:B2,    
				     igrid=broadcast(x->range(x[1],stop=x[end],step=300.),mgrid),
				     parameterization=parameterization,   verbose=false,
				     nworker=1)


		global pa_parallel=SeisInvExpt(attrib_mod, Migr(), srcwav=srcwav, ageom=ageom, tgrid=tgrid, mediumm=medium0, 
				     mediumm_obs=medium,  
				     mediumm0=medium0,
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
	global expt=x->SeisInvExpt(FdtdBorn(), x, 
			     srcwav=srcwav, ageom=ageom, tgrid=tgrid, mediumm=medium0, 
	     		     mediumm0=medium0,
			     mediumm_obs=medium,  
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
	expt=x->SeisInvExpt(Fdtd(), x, 
			     srcwav=srcwav, ageom=ageom, tgrid=tgrid, mediumm=medium0, 
	     		     mediumm0=medium0,
			     mediumm_obs=medium,  
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


