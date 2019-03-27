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
    "page": "Generate snaps",
    "title": "Generate snaps",
    "category": "page",
    "text": "EditURL = \"https://github.com/TRAVIS_REPO_SLUG/blob/master/\"A forward experiment, where the seismic data are generated using some models and acquisition parameters from our gallery. Forward modeling consists of a finite-difference simulation of the acoustic wave-equation. We specifically aim to save the snapshots, at given time steps in SeisForwExpt."
},

{
    "location": "Fdtd/create_snaps/#Loading-some-packages-1",
    "page": "Generate snaps",
    "title": "Loading some packages",
    "category": "section",
    "text": "using GeoPhyInv\nusing Statistics"
},

{
    "location": "Fdtd/create_snaps/#Setting-up-the-variables-necessary-to-create-the-Expt-1",
    "page": "Generate snaps",
    "title": "Setting up the variables necessary to create the Expt",
    "category": "section",
    "text": "model=GIPh.Gallery.Seismic(:acou_homo1); # load a simple homogeneous acoustic model from the gallery\nGIPh.Models.Seismic_addon!(model, randn_perc=0.01); # add some random noise to the model\nacqgeom=GIPh.Gallery.Geom(model.mgrid,:xwell); # load a simple acquisition geometry using `mgrid` of the seismic model\ntgrid = range(0.0,stop=2.0,length=2000) # generate a time grid\nwav = GIPh.Utils.Wavelets.ricker(10.0, tgrid, tpeak=0.25,); # ricker wavelet\nacqsrc=GIPh.Acquisition.Src_fixed(acqgeom.nss,1,[:P],wav,tgrid); # distribute the same source wavelet to all the supsersources\n@info \"We are ready for the modeling.\""
},

{
    "location": "Fdtd/create_snaps/#Final-step-1",
    "page": "Generate snaps",
    "title": "Final step",
    "category": "section",
    "text": "One can plot the model, source and receivers using these commands: using Plots; p1=JP.seismic(model); JP.geom!(acqgeom); plot(p1); Now we have all the required variables to create SeisForwExpt object and prepare the forward modelling. While creating, we switched the snaps_flag on, and instructed recording field at tsnaps. Once the Expt object is created, do the modelling \"without approximately any memory allocations\" using mod!paE=SeisForwExpt(model=model,\n	acqgeom=[acqgeom], acqsrc=[acqsrc],\n	snaps_flag=true,\n	tsnaps=[0.3, 0.4, 0.5],\n	tgridmod=tgrid, verbose=true);\n\n@time mod!(paE);"
},

{
    "location": "Fdtd/create_snaps/#Extracting-snaps-from-Expt-1",
    "page": "Generate snaps",
    "title": "Extracting snaps from Expt",
    "category": "section",
    "text": "snaps=paE[:snaps,1]; # extracting snaps of the first supersource\nsnaps=paE[:snaps,2]; # second supersource\n@info string(\"The dimensions of the snaps are (nz,nx,nt)=\", size(snaps))We can now plot snapshots using these commands: p1=[heatmap(snaps[:,:,ii]) for ii in 1:3]; plot(p1..., layout=(1,3), aspect_ratio=:equal)"
},

{
    "location": "Poisson/forw/#",
    "page": "Record data",
    "title": "Record data",
    "category": "page",
    "text": "EditURL = \"https://github.com/TRAVIS_REPO_SLUG/blob/master/\""
},

{
    "location": "Poisson/forw/#Loading-some-packages-1",
    "page": "Record data",
    "title": "Loading some packages",
    "category": "section",
    "text": "using GeoPhyInv\nusing SparseArrays\nusing StatsBase\nusing LinearAlgebra\nusing Random\nusing ProgressMeter\nusing LinearAlgebra\nusing Test\nusing ForwardDiff\nusing Calculus"
},

{
    "location": "Poisson/forw/#Solve-for-ψ-in-a-PoissonExpt-1",
    "page": "Record data",
    "title": "Solve for ψ in a PoissonExpt",
    "category": "section",
    "text": "This module represents an explicit, direct sparse 2D finite-difference Poisson solver for heterogeneous media, i.e. media having spatially varying (space-dependent) medium parameters. Current implementation assumes Neumann boundary conditions at all the boundaries.Consider the following Poisson experiment:(σ(xz)) ψ(t) = (Q(xz)) p(t)Q = k * Q_v  ηWe start with the dimensions and spatial grids are allocated as follows.nx=21\nnz=21\nnt=4\nnznx=nz*nx\nmgrid=[range(-div(nz,2), step=1.0, length=nz), range(-div(nx,2), step=1.0, length=nx)]\ntgrid=range(0.0,step=0.5, length=nt)\n@info \"Grids are all set.\"Now lets allocate the inputs for a toy experiment. These medium parameters are used to generate the observed field ψ.Qv=abs.(randn(nz,nx))\nη=abs.(randn(nz,nx))\nk=abs.(randn(nz,nx))\nσ=abs.(randn(nz,nx))\np=randn(nz,nx,nt)\n@info \"Medium parameters allocated.\""
},

{
    "location": "Poisson/forw/#Acquisition-1",
    "page": "Record data",
    "title": "Acquisition",
    "category": "section",
    "text": "Now, we will generate an acquisition geometry and allocate a projection matrix ACQ.acqgeom=GIPh.Acquisition.Geom_circ(nss=1,nr=30,rad=[5.,5.]);\nACQ=GIPh.Acquisition.ACQmat(acqgeom,mgrid);\n@info \"ACQ will be used to project ψ onto receivers.\""
},

{
    "location": "Poisson/forw/#Generate-PoissonExpt-and-then-applying-mod!-1",
    "page": "Record data",
    "title": "Generate PoissonExpt and then applying mod!",
    "category": "section",
    "text": "This will firstapply operator A=(Q(xz)) on a field p;\nthen apply ((σ(xz)))^-1 in order to solve for ψ;\nfinally, records ψ at the receiver locations to generate data.paE=PoissonExpt(p, tgrid, mgrid, Qv, k, η, σ, ACQ)\nmod!(paE)"
},

{
    "location": "Poisson/forw/#Extracting-data-from-Expt-1",
    "page": "Record data",
    "title": "Extracting data from Expt",
    "category": "section",
    "text": "data=paE[:data]\n@info string(\"The dimensions of data are (nt,nr)=\",size(data))"
},

{
    "location": "Poisson/test_born/#",
    "page": "Born map",
    "title": "Born map",
    "category": "page",
    "text": "EditURL = \"https://github.com/TRAVIS_REPO_SLUG/blob/master/\"The linearized forward modeling operator F and its adjoint (aka Migration operator) are the building blocks of iterative optimization schemes. For the PoissonExpt, we have the functionality to compute F*x without explicitly storing the operator matrix (see LinearMaps.jl). The perturbed field δψ due to a perturbation in σ is given byδψ=-A¹(σ₀)A(δσ)ψ₀where(σ₀(xz)) ψ₀(t)=A(σ₀)ψ₀(t)=f(t)Lets start a tutorial."
},

{
    "location": "Poisson/test_born/#Load-some-packages-1",
    "page": "Born map",
    "title": "Load some packages",
    "category": "section",
    "text": "using GeoPhyInv\nusing SparseArrays\nusing StatsBase\nusing LinearAlgebra\nusing Random\nusing ProgressMeter\nusing LinearAlgebra\nusing Test\nusing ForwardDiff\nusing Calculus\nusing LinearMaps"
},

{
    "location": "Poisson/test_born/#Setting-up-the-spatial-and-temporal-grids-1",
    "page": "Born map",
    "title": "Setting up the spatial and temporal grids",
    "category": "section",
    "text": "nx=21\nnz=21\nnt=4\nnznx=nz*nx\nmgrid=[range(-div(nz,2), step=1.0, length=nz), range(-div(nx,2), step=1.0, length=nx)]\ntgrid=range(0.0,step=0.5, length=nt)\n@info \"Grids are all set.\""
},

{
    "location": "Poisson/test_born/#Allocating-medium-parameters-1",
    "page": "Born map",
    "title": "Allocating medium parameters",
    "category": "section",
    "text": "Qv=abs.(randn(nz,nx))\nη=abs.(randn(nz,nx))\nk=abs.(randn(nz,nx))\nσ=abs.(randn(nz,nx))\np=randn(nz,nx,nt)\n\nσobs=abs.(randn(nz,nx))\nQobs=abs.(randn(nz,nx))\nnr=10 # number of abstract receivers\nACQ=sprandn(nr,nz*nx,0.6); # choose a random acquisition operator\n@info \"We are ready for the PoissonExpt.\""
},

{
    "location": "Poisson/test_born/#Create-an-Expt,-and-then-extract-a-linear-forward-map-out-of-it-1",
    "page": "Born map",
    "title": "Create an Expt, and then extract a linear forward map out of it",
    "category": "section",
    "text": "paE=PoissonExpt(p, tgrid, mgrid, Qv, k, η, σ, ACQ, σobs=σobs, Qobs=Qobs,)\nF=operator_Born(paE, σ); # extract the linearized forward operator from `Expt`\nGIPh.Utils.test_linearmap(F) # finally do some tests on the linearmap"
},

{
    "location": "Poisson/test_born/#Usage-1",
    "page": "Born map",
    "title": "Usage",
    "category": "section",
    "text": "δx=randn(size(F,2)) # random model pertubation\nδd=F*δx # corresponding pertubation in data\n@info string(\"Length of data: (nt*nr)=\",length(δd))"
},

]}
