var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#The-Expt-Datatype-1",
    "page": "Home",
    "title": "The Expt Datatype",
    "category": "section",
    "text": "The methods in this package numerically solve some differential equations commonly faced in geophysical inverse problems. The functionality of this package revolves around the mutable Expt types. Firstly, most of the memory necessary to perform a given experiment is allocated while creating the Expt variables. Then these variables are input to in-place functions (e.g., mod!)  which as per Julia convention ends with an exclamation mark, to actually perform the experiment task. For example, the current Expt types within the realm of this package include:SeisForwExpt is the seismic (acoustic) forward modeling experiment  ;\nSeisInvExpt is the type for seismic inversion experiment, including migration;\nPoissonExpt is type for the solving the Poisson experiment.To get started, as an example, simply load a seismic inversion experiment already defined in our package gallery into REPL:using GeoPhyInv # load GIPh (after installation)\npaE=GIPh.Gallery.SeisInvExpt(:pizza); # \"pizza\" is the name of the experiment"
},

{
    "location": "Fdtd/create_snaps/#",
    "page": "SeismicExpt: generate snaps",
    "title": "SeismicExpt: generate snaps",
    "category": "page",
    "text": "EditURL = \"https://github.com/TRAVIS_REPO_SLUG/blob/master/\"A forward experiment, where the seismic data are generated using some models and acquisition parameters from our gallery. Forward modeling consists of a finite-difference simulation of the acoustic wave-equation. We specifically aim to save the snapshots, at given time steps in SeisForwExpt."
},

{
    "location": "Fdtd/create_snaps/#Loading-some-packages...-1",
    "page": "SeismicExpt: generate snaps",
    "title": "Loading some packages...",
    "category": "section",
    "text": "using GeoPhyInv\nusing StatisticsThen, load a simple homogeneous acoustic model from the gallery.model=GIPh.Gallery.Seismic(:acou_homo1);Add some random noise to the model.GIPh.Models.Seismic_addon!(model, randn_perc=0.01);Load a simple acquisition geometry from the gallery using the field mgrid of the seismic model.acqgeom=GIPh.Gallery.Geom(model.mgrid,:xwell);Plot the model, source and receivers using these commands: using Plots p1=JP.seismic(model); JP.geom!(acqgeom); plot(p1);Generate a time grid.tgrid = range(0.0,stop=2.0,length=1000)Ricker wavelet.wav = GIPh.Utils.Wavelets.ricker(10.0, tgrid, tpeak=0.25,);Distribute the same source wavelet to all the supsersources.acqsrc=GIPh.Acquisition.Src_fixed(acqgeom.nss,1,[:P],wav,tgrid);Now we have all the required variables to create SeisForwExpt object and prepare the forward modelling. While creating, we switched the snaps_flag on, and instructed recording field at tsnaps. Once the Expt object is created, do the modelling \"without approximately any memory allocations\" using mod!paE=SeisForwExpt(model=model,\n	acqgeom=[acqgeom], acqsrc=[acqsrc],\n	snaps_flag=true,\n	tsnaps=[0.3, 0.4, 0.5],\n	tgridmod=tgrid, verbose=true);\n\n@time mod!(paE);Extracting snaps of the first supersource. The dimensions of the snaps is [nz,nx,nt].snaps=paE[:snaps,1];Extracting snaps of the second supersource.snaps=paE[:snaps,2];We can plot snapshots using these commands: p1=[heatmap(snaps[:,:,ii]) for ii in 1:3]; plot(p1..., layout=(1,3), aspect_ratio=:equal)This page was generated using Literate.jl."
},

{
    "location": "Poisson/adj_state_expt/#",
    "page": "PoissonExpt",
    "title": "PoissonExpt",
    "category": "page",
    "text": "EditURL = \"https://github.com/TRAVIS_REPO_SLUG/blob/master/\"This module represents an explicit, direct sparse 2D finite-difference Poisson solver for heterogeneous media, i.e. media having spatially varying (space-dependent) medium parameters. The following functionality is currently available in this module:Apply operator A=(σ(xz)) on a field ψ.\nApply A^-1 in order to solve for ψ in Aψ=p, given p.Current implementation assumes Neumann boundary conditions at all the boundaries.As a demo, start with loading some packages.using SparseArrays\nusing StatsBase\nusing LinearAlgebra\nusing Random\nusing ProgressMeter\nusing LinearAlgebra\nusing Test\nusing ForwardDiff\nusing CalculusDimensions and spatial grids are allocated as follows.nx=20\nnz=20\nnt=4\nnznx=nz*nx\nmgrid=[range(0.0,step=0.5, length=nz), range(0.0,step=0.5, length=nz)]\ntgrid=range(0.0,step=0.5, length=nt)\n\np=randn(nznx+1)\np[end]=0.0\nQv=abs.(randn(nz,nx))\nη=abs.(randn(nz,nx))\nk=abs.(randn(nz,nx))\nσ=abs.(randn(nz,nx))\nσobs=abs.(randn(nz,nx))\nQobs=abs.(randn(nz,nx))\npsi0=randn(nznx)\n\nsnaps=randn(nz,nx,nt)Generate the parameters for the Poisson experiment.paE=J.Poisson.ParamExpt(snaps, tgrid, mgrid, Qv, k, η, σ, σobs=σobs, Qobs=Qobs)J.Poisson.updateLP!(paE, paE.Q)Calculate the data misfit.@time f=J.Poisson.func(σ,paE)Compute the gradient with finite-difference approximation.f = x->J.Poisson.func(x, paE);\ngcalc = Calculus.gradient(f);\ng_fd=reshape(gcalc(vec(σ)),nz,nx);Compute the gradient using adjoint state method.@time J.Poisson.mod!(paE,σ,J.Poisson.FGσ())\ng=paE.g;Check gradient accuracy.@test ≈(g,g_fd, rtol=1e-5)This page was generated using Literate.jl."
},

]}
