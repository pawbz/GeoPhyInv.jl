```@meta
EditURL = "@__REPO_ROOT_URL__/"
```

### Loading some packages

```@example forw
using GeoPhyInv
using SparseArrays
using StatsBase
using LinearAlgebra
using Random
using ProgressMeter
using LinearAlgebra
using Test
using ForwardDiff
using Calculus
```

### Solve for ``ψ`` in a `PoissonExpt`
This module represents an explicit, direct sparse 2D finite-difference Poisson solver for heterogeneous media,
i.e. media having spatially varying (space-dependent) medium parameters.
Current implementation assumes Neumann boundary conditions at all the boundaries.

Consider the following Poisson experiment:
```math
∇⋅(σ(x,z)∇) ψ(t) = ∇⋅(Q(x,z)∇) p(t),
```
```math
Q = k * Q_v / η.
```
We start with the dimensions and spatial grids are allocated as follows.

```@example forw
nx=21
nz=21
nt=4
nznx=nz*nx
mgrid=[range(-div(nz,2), step=1.0, length=nz), range(-div(nx,2), step=1.0, length=nx)]
tgrid=range(0.0,step=0.5, length=nt)
@info "Grids are all set."
```

Now lets allocate the inputs for a toy experiment.
These medium parameters are used to generate the *observed* field ``ψ``.

```@example forw
Qv=abs.(randn(nz,nx))
η=abs.(randn(nz,nx))
k=abs.(randn(nz,nx))
σ=abs.(randn(nz,nx))
p=randn(nz,nx,nt)
@info "Medium parameters allocated."
```

### Acquisition
Now, we will generate an acquisition geometry and allocate a projection matrix `ACQ`.

```@example forw
acqgeom=GIPh.Acquisition.Geom_circ(nss=1,nr=30,rad=[5.,5.]);
ACQ=GIPh.Acquisition.ACQmat(acqgeom,mgrid);
@info "ACQ will be used to project ψ onto receivers."
```

### Generate `PoissonExpt` and then applying `mod!`
This will first
* apply operator ``A=∇⋅(Q(x,z)∇)`` on a field ``p``;
* then apply ``(∇⋅(σ(x,z)∇))^{-1}`` in order to solve for ``ψ``;
* finally, records ``ψ`` at the receiver locations to generate data.

```@example forw
paE=PoissonExpt(p, tgrid, mgrid, Qv, k, η, σ, ACQ)
mod!(paE)
```

### Extracting data from `Expt`

```@example forw
data=paE[:data]
@info string("The dimensions of data are (nt,nr)=",size(data))
```

