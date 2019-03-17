```@meta
EditURL = "https://github.com/TRAVIS_REPO_SLUG/blob/master/"
```

This module represents an explicit, direct sparse 2D finite-difference Poisson solver for heterogeneous media,
i.e. media having spatially varying (space-dependent) medium parameters.
The following functionality is currently available in this module:
* Apply operator ``A=∇⋅(σ(x,z)∇)`` on a field ``ψ``.
* Apply ``A^{-1}`` in order to solve for ``ψ`` in ``Aψ=p``, given ``p``.
Current implementation assumes Neumann boundary conditions at all boundaries.

As a demo, start with loading some packages.

```@example adj_state_expt
using SparseArrays
using StatsBase
using LinearAlgebra
using Random
using ProgressMeter
using LinearAlgebra
using Test
using ForwardDiff
using Calculus
#include("core.jl")
#include("expt.jl")
```

Dimensions and spatial grids are allocated as follows.

```@example adj_state_expt
nx=20
nz=20
nt=4
nznx=nz*nx
mgrid=[range(0.0,step=0.5, length=nz), range(0.0,step=0.5, length=nz)]
tgrid=range(0.0,step=0.5, length=nt)

p=randn(nznx+1)
p[end]=0.0
Qv=abs.(randn(nz,nx))
η=abs.(randn(nz,nx))
k=abs.(randn(nz,nx))
σ=abs.(randn(nz,nx))
σobs=abs.(randn(nz,nx))
Qobs=abs.(randn(nz,nx))
psi0=randn(nznx)

snaps=randn(nz,nx,nt)
```

Generate the parameters for the Poisson experiment.

```@example adj_state_expt
paE=J.Poisson.ParamExpt(snaps, tgrid, mgrid, Qv, k, η, σ, σobs=σobs, Qobs=Qobs)
```

```@example adj_state_expt
J.Poisson.updateLP!(paE, paE.Q)
```

Calculate the data misfit.

```@example adj_state_expt
@time f=J.Poisson.func(σ,paE)
```

Compute the gradient with finite-difference approximation.

```@example adj_state_expt
f = x->J.Poisson.func(x, paE);
gcalc = Calculus.gradient(f);
g_fd=reshape(gcalc(vec(σ)),nz,nx);
```

Compute the gradient using adjoint state method.

```@example adj_state_expt
@time J.Poisson.mod!(paE,σ,J.Poisson.FGσ())
g=paE.g;
```

Check gradient accuracy.

```@example adj_state_expt
@test ≈(g,g_fd, rtol=1e-5)
```

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

