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
    "location": "Poisson/adj_state_expt/#",
    "page": "Poisson Solver",
    "title": "Poisson Solver",
    "category": "page",
    "text": "EditURL = \"https://github.com/TRAVIS_REPO_SLUG/blob/master/\"This module represents an explicit, direct sparse 2D finite-difference Poisson solver for heterogeneous media, i.e. media having spatially varying (space-dependent) medium parameters. The following functionality is currently available in this module:Apply operator A=(σ(xz)) on a field ψ.\nApply A^-1 in order to solve for ψ in Aψ=p, given p.Current implementation assumes Neumann boundary conditions at all boundaries.As a demo, start with loading some packages.using SparseArrays\nusing StatsBase\nusing LinearAlgebra\nusing Random\nusing ProgressMeter\nusing LinearAlgebra\nusing Test\nusing ForwardDiff\nusing Calculus\n#include(\"core.jl\")\n#include(\"expt.jl\")Dimensions and spatial grids are allocated as follows.nx=20\nnz=20\nnt=4\nnznx=nz*nx\nmgrid=[range(0.0,step=0.5, length=nz), range(0.0,step=0.5, length=nz)]\ntgrid=range(0.0,step=0.5, length=nt)\n\np=randn(nznx+1)\np[end]=0.0\nQv=abs.(randn(nz,nx))\nη=abs.(randn(nz,nx))\nk=abs.(randn(nz,nx))\nσ=abs.(randn(nz,nx))\nσobs=abs.(randn(nz,nx))\nQobs=abs.(randn(nz,nx))\npsi0=randn(nznx)\n\nsnaps=randn(nz,nx,nt)Generate the parameters for the Poisson experiment.paE=J.Poisson.ParamExpt(snaps, tgrid, mgrid, Qv, k, η, σ, σobs=σobs, Qobs=Qobs)J.Poisson.updateLP!(paE, paE.Q)Calculate the data misfit.@time f=J.Poisson.func(σ,paE)Compute the gradient with finite-difference approximation.f = x->J.Poisson.func(x, paE);\ngcalc = Calculus.gradient(f);\ng_fd=reshape(gcalc(vec(σ)),nz,nx);Compute the gradient using adjoint state method.@time J.Poisson.mod!(paE,σ,J.Poisson.FGσ())\ng=paE.g;Check gradient accuracy.@test ≈(g,g_fd, rtol=1e-5)This page was generated using Literate.jl."
},

]}
