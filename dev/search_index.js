var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#Toolbox-for-GeoPhysical-Inversion-(GPI)-1",
    "page": "Home",
    "title": "Toolbox for GeoPhysical Inversion (GPI)",
    "category": "section",
    "text": "using GeoPhyInvThe methods in this software numerically solve some  differential equations of geophysics in Julia. Currently, the equations within the realm of this package  include:The methods within the realm of seismic modeling and inversion return:the linearized forward modeling operator F, such that Fx can be computed without explicitly storing the operator matrix (  see LinearMaps.jl);the imaging/migration operator F*;These maps are the building blocks of iterative optimization schemes.the acoustic $a\\otimes b$Forward problem, where the seismic data are generated using synthetic Earth models and the acquisition parameters  corresponding to a seismic experiment. Forward modeling consists of a finite-difference simulation, followed by convolutions in the time domain using the source and receiver filters. The details about our finite-difference  scheme are given in XXX.  And the filters corresponding to  sources and receivers are described in XXX.Can perform inversion of synthetic scenarios.First, the seismic data are modeled as in the forward problem. Then the  data are used to perform full waveform inversion (FWI). The inverse  problem estimates the Earth models and the source and receiver filters  that resulted from the data. This task is necessary to test the performance of the inversion algorithm  in various geological scenarios using different acquisition parameters.Read the measured seismic field data and parameters from a seismic experiment  to perform inversion like in the previous task.  The data measured in the field are not in a suitable format  yet for this software.  Pre-processing is necessary before it can be used as described. Also, the acquisition parameters from the field should be  converted to suitable 2-D coordinates as described."
},

{
    "location": "FWI/born_map/#",
    "page": "Seismic Born Modeling",
    "title": "Seismic Born Modeling",
    "category": "page",
    "text": "EditURL = \"https://github.com/TRAVIS_REPO_SLUG/blob/master/\"for scenario in [:downhole, :pizza]\n	println(\"@@@@@@@@@@@@TESTING \", scenario)\n	for rfields in [[:P], [:Vx], [:Vz]]\n		pa, model=JG.xfwi_problem(scenario, born_flag=true, rfields=rfields)\n\n\n		F=JF.operator_Born(pa);\n\n		x1=randn(size(F,2))\n		x2=randn(size(F,2))\n		x12=x1.+x2\n\n\n		d12=F*x12\n		d1=F*x1\n		δmodtt1=copy(pa.paf.c.δmodtt)\n		d2=F*x2\n		δmodtt2=copy(pa.paf.c.δmodtt)\n\n\n		d12new=d1.+d2\n\n		f=Misfits.error_squared_euclidean!(nothing, d12, d12new, nothing, norm_flag=true)\n\n		@test f<1e-25\n\n\n		function adjtest()\n			x=randn(size(F,2))\n			y=randn(size(F,1))\n			a=LinearAlgebra.dot(y,F*x)\n			b=LinearAlgebra.dot(x,adjoint(F)*y)\n			c=LinearAlgebra.dot(x, transpose(F)*F*x)\n			println(\"adjoint test: \", a, \"\\t\", b)\n			@test isapprox(a,b,rtol=1e-5)\n			println(\"must be positive: \", c)\n			@test c>0.0\n		end\n\n\n		for i in 1:3\n			adjtest()\n		end\n	end\nendThis page was generated using Literate.jl."
},

{
    "location": "FWI/gradient_accuracy/#",
    "page": "Seismic Full Waveform Inversion",
    "title": "Seismic Full Waveform Inversion",
    "category": "page",
    "text": "EditURL = \"https://github.com/TRAVIS_REPO_SLUG/blob/master/\"model = J.Gallery.Seismic(:acou_homo2);\nJ.Models.Seismic_addon!(model,randn_perc=1, fields=[:χvp,:χρ])\n\nmodel0 = J.Gallery.Seismic(:acou_homo2);\nJ.Models.Seismic_addon!(model0, randn_perc=1, fields=[:χvp,:χρ])\n\nacqgeom=J.Acquisition.Geom_fixed(model,1,10)\nacqsrc=J.Acquisition.Src_fixed_mod(acqgeom.nss,1,[:P],mod=model, nλ=3, tmaxfrac=1.0)\ntgrid=acqsrc.tgrid\n\nparameterization=[:χvp, :χρ, :null]\n\nmgrid=model.mgrid\n\n@testset \"test parallel implementation during gradient\" begin\n	for attrib_mod in [JF.ModFdtd(), JF.ModFdtdBorn()]\n		pa=JF.Param(acqsrc, acqgeom, tgrid, attrib_mod, model0,\n				     modm_obs=model,\n				     modm0=model0,\n				     igrid_interp_scheme=:B2,\n				     igrid=broadcast(x->range(x[1],stop=x[end],step=300.),mgrid),\n				     parameterization=parameterization,   verbose=false,\n				     nworker=1)\n\n\n		pa_parallel=JF.Param(acqsrc, acqgeom, tgrid, attrib_mod, model0,\n				     modm_obs=model,\n				     modm0=model0,\n				     igrid_interp_scheme=:B2,\n				     igrid=broadcast(x->range(x[1],stop=x[end],step=300.),mgrid),\n				     parameterization=parameterization,   verbose=false,\n				     nworker=nothing)\n\n		result=JF.xfwi!(pa, JF.Migr(), attrib_mod)\n\n		result_parallel=JF.xfwi!(pa_parallel, JF.Migr(), attrib_mod)\n\n		@test result[2] ≈ result_parallel[2]\n	end\nend\n\n@testset \"Testing Born Modeling and its gradient\" begin\n\n	pa=JF.Param(acqsrc, acqgeom, tgrid, JF.ModFdtdBorn(), model0,\n	     		     modm0=model0,\n			     modm_obs=model,\n			     igrid_interp_scheme=:B2,\n			     igrid=broadcast(x->range(x[1],stop=x[end],step=350.),mgrid),\n			     parameterization=parameterization,   verbose=false)\n\n\n	JF.xfwi!(pa, JF.LS(), JF.ModFdtdBorn(),  bounded_flag=true, solver=:ipopt,\n			ipopt_options=[[\"max_iter\", 0],[\"derivative_test\", \"first-order\"]])\n\n	result=JF.xfwi!(pa, JF.Migr(), JF.ModFdtdBorn())\n\n	pa_fd=deepcopy(pa);\n	result_fd=JF.xfwi!(pa_fd, JF.Migr_fd(), JF.ModFdtdBorn())\n\n	f=Misfits.error_squared_euclidean!(nothing, result[2], result_fd[2], nothing, norm_flag=true)\n\n	@test f<1e-15\nend\n\n@testset \"Testing gradient LS FWI\" begin\n	pa=JF.Param(acqsrc, acqgeom, tgrid, JF.ModFdtd(), model0,\n			     modm_obs=model,\n			     igrid_interp_scheme=:B2,\n			     igrid=broadcast(x->range(x[1],stop=x[end],step=350.),mgrid),\n			     parameterization=parameterization,   verbose=false)\n\n	JF.xfwi!(pa, JF.LS(), JF.ModFdtd(),  bounded_flag=true, solver=:ipopt,\n			ipopt_options=[[\"max_iter\", 0],[\"derivative_test\", \"first-order\"]])\n\n\n	result=JF.xfwi!(pa, JF.Migr(), JF.ModFdtd())\n\n	pa_fd=deepcopy(pa);\n	result_fd=JF.xfwi!(pa_fd, JF.Migr_fd(), JF.ModFdtd())\n\n	f=Misfits.error_squared_euclidean!(nothing, result[2], result_fd[2], nothing, norm_flag=true)\n\n	@test f<1e-15\nendThis page was generated using Literate.jl."
},

{
    "location": "Poisson/adj_state_expt/#",
    "page": "Poisson Solver",
    "title": "Poisson Solver",
    "category": "page",
    "text": "EditURL = \"https://github.com/TRAVIS_REPO_SLUG/blob/master/\"This module represents an explicit, direct sparse 2D finite-difference Poisson solver for heterogeneous media, i.e. media having spatially varying (space-dependent) medium parameters. The following functionality is currently available in this module:Apply operator A=(σ(xz)) on a field ψ.\nApply A^-1 in order to solve for ψ in Aψ=p, given p.Current implementation assumes Neumann boundary conditions at all the boundaries.As a demo, start with loading some packages.using SparseArrays\nusing StatsBase\nusing LinearAlgebra\nusing Random\nusing ProgressMeter\nusing LinearAlgebra\nusing Test\nusing ForwardDiff\nusing CalculusDimensions and spatial grids are allocated as follows.nx=20\nnz=20\nnt=4\nnznx=nz*nx\nmgrid=[range(0.0,step=0.5, length=nz), range(0.0,step=0.5, length=nz)]\ntgrid=range(0.0,step=0.5, length=nt)\n\np=randn(nznx+1)\np[end]=0.0\nQv=abs.(randn(nz,nx))\nη=abs.(randn(nz,nx))\nk=abs.(randn(nz,nx))\nσ=abs.(randn(nz,nx))\nσobs=abs.(randn(nz,nx))\nQobs=abs.(randn(nz,nx))\npsi0=randn(nznx)\n\nsnaps=randn(nz,nx,nt)Generate the parameters for the Poisson experiment.paE=J.Poisson.ParamExpt(snaps, tgrid, mgrid, Qv, k, η, σ, σobs=σobs, Qobs=Qobs)J.Poisson.updateLP!(paE, paE.Q)Calculate the data misfit.@time f=J.Poisson.func(σ,paE)Compute the gradient with finite-difference approximation.f = x->J.Poisson.func(x, paE);\ngcalc = Calculus.gradient(f);\ng_fd=reshape(gcalc(vec(σ)),nz,nx);Compute the gradient using adjoint state method.@time J.Poisson.mod!(paE,σ,J.Poisson.FGσ())\ng=paE.g;Check gradient accuracy.@test ≈(g,g_fd, rtol=1e-5)This page was generated using Literate.jl."
},

]}
