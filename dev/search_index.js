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
    "location": "Fdtd/create_snaps/#Loading-some-packages-1",
    "page": "SeismicExpt: generate snaps",
    "title": "Loading some packages",
    "category": "section",
    "text": "using GeoPhyInv\nusing Statistics"
},

{
    "location": "Fdtd/create_snaps/#Setting-up-the-variables-necessary-to-create-the-Expt-1",
    "page": "SeismicExpt: generate snaps",
    "title": "Setting up the variables necessary to create the Expt",
    "category": "section",
    "text": "model=GIPh.Gallery.Seismic(:acou_homo1); # load a simple homogeneous acoustic model from the gallery\nGIPh.Models.Seismic_addon!(model, randn_perc=0.01); # add some random noise to the model\nacqgeom=GIPh.Gallery.Geom(model.mgrid,:xwell); # load a simple acquisition geometry using `mgrid` of the seismic model\ntgrid = range(0.0,stop=2.0,length=2000) # generate a time grid\nwav = GIPh.Utils.Wavelets.ricker(10.0, tgrid, tpeak=0.25,); # ricker wavelet\nacqsrc=GIPh.Acquisition.Src_fixed(acqgeom.nss,1,[:P],wav,tgrid); # distribute the same source wavelet to all the supsersources\n@info \"We are ready for the modeling.\""
},

{
    "location": "Fdtd/create_snaps/#Final-step-1",
    "page": "SeismicExpt: generate snaps",
    "title": "Final step",
    "category": "section",
    "text": "One can plot the model, source and receivers using these commands: using Plots; p1=JP.seismic(model); JP.geom!(acqgeom); plot(p1); Now we have all the required variables to create SeisForwExpt object and prepare the forward modelling. While creating, we switched the snaps_flag on, and instructed recording field at tsnaps. Once the Expt object is created, do the modelling \"without approximately any memory allocations\" using mod!paE=SeisForwExpt(model=model,\n	acqgeom=[acqgeom], acqsrc=[acqsrc],\n	snaps_flag=true,\n	tsnaps=[0.3, 0.4, 0.5],\n	tgridmod=tgrid, verbose=true);\n\n@time mod!(paE);"
},

{
    "location": "Fdtd/create_snaps/#Extracting-snaps-from-Expt-1",
    "page": "SeismicExpt: generate snaps",
    "title": "Extracting snaps from Expt",
    "category": "section",
    "text": "snaps=paE[:snaps,1]; # extracting snaps of the first supersource\nsnaps=paE[:snaps,2]; # second supersource\n@info \"The dimensions of the snaps is [nz,nx,nt].\"We can now plot snapshots using these commands: p1=[heatmap(snaps[:,:,ii]) for ii in 1:3]; plot(p1..., layout=(1,3), aspect_ratio=:equal)This page was generated using Literate.jl."
},

{
    "location": "Poisson/adj_state_expt/#",
    "page": "PoissonExpt",
    "title": "PoissonExpt",
    "category": "page",
    "text": "EditURL = \"https://github.com/TRAVIS_REPO_SLUG/blob/master/\"This module represents an explicit, direct sparse 2D finite-difference Poisson solver for heterogeneous media, i.e. media having spatially varying (space-dependent) medium parameters. The following functionality is currently available in this module:Apply operator A=(σ(xz)) on a field ψ to get p.\nApply A^-1 in order to solve for ψ in Aψ=p, given p.Current implementation assumes Neumann boundary conditions at all the boundaries.As a demo, start with loading some packages.using SparseArrays\nusing StatsBase\nusing LinearAlgebra\nusing Random\nusing ProgressMeter\nusing LinearAlgebra\nusing Test\nusing ForwardDiff\nusing CalculusFrom here on, consider the following Poisson experiment:(σ(xz)) ψ(t) = (Q(xz)) p(t)Q = k * Qv  ηDimensions and spatial grids are allocated as follows.nx=21\nnz=21\nnt=4\nnznx=nz*nx\nmgrid=[range(-div(nz,2), step=1.0, length=nz), range(-div(nx,2), step=1.0, length=nx)]\ntgrid=range(0.0,step=0.5, length=nt)Now lets allocate the inputs for a toy experiment.Qv=abs.(randn(nz,nx))\nη=abs.(randn(nz,nx))\nk=abs.(randn(nz,nx))\nσ=abs.(randn(nz,nx))\np=randn(nz,nx,nt)\n\nfor it in 1:nt\n	for ix in 1:nx\n		for iz in 1:nz\n			p[iz,ix,it]=inv(abs(iz-11)+1)*inv(abs(ix-11)+1)\n		end\n	end\nendThese medium parameters are used to generate the observed field ψ.σobs=abs.(randn(nz,nx))\nQobs=abs.(randn(nz,nx))\n\nfill!(σobs, 1.0)\nfill!(Qobs, 1.0)\nfill!(Qv, 1.0)\nfill!(η, 1.0)\nfill!(k, 1.0)\nσobs[11,11]=20.0\nfill!(σ, 1.0)Generate the configuration for the Poisson experiment, using the arrays above.acqgeom=J.Acquisition.Geom_circ(nss=1,nr=30,rad=[5.,5.])\n#ACQ=sprandn(200,441,0.9)\nACQ=J.Acquisition.ACQmat(acqgeom,mgrid)\npaE=J.Poisson.ParamExpt(p, tgrid, mgrid, Qv, k, η, σ, ACQ, σobs=σobs, Qobs=Qobs,)\n#paE=J.Poisson.ParamExpt(p, tgrid, mgrid, Qv, k, η, σ, σobs=σobs, Qobs=Qobs,)Calculate the data misfit.@time f=J.Poisson.func(σ,paE)Compute the gradient with finite-difference approximation.f = x->J.Poisson.func(x, paE);\ngcalc = Calculus.gradient(f);\ng_fd=reshape(gcalc(vec(σ)),nz,nx);Compute the gradient using adjoint state method.@time J.Poisson.mod!(paE,σ,J.Poisson.FGσ())\ng=paE.g;Check gradient accuracy.@test ≈(g,g_fd, rtol=1e-5)This page was generated using Literate.jl."
},

]}
