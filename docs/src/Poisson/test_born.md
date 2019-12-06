```@meta
EditURL = "@__REPO_ROOT_URL__/"
```

The linearized forward modeling operator `F` and its adjoint (aka Migration operator)
are the building blocks of iterative optimization schemes.
For the `PoissonExpt`, we have the
functionality to compute `F*x` without explicitly storing the operator matrix
(see `LinearMaps.jl`).
The perturbed field `δψ` due to a perturbation in `σ` is given by
```math
δψ=-A⁻¹(σ₀)A(δσ)ψ₀,
```
where
```math
∇⋅(σ₀(x,z)∇) ψ₀(t)=A(σ₀)ψ₀(t)=f(t).
```
Lets start a tutorial.

### Load some packages

```@example test_born
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
using LinearMaps
```

### Setting up the spatial and temporal grids

```@example test_born
nx=21
nz=21
nt=4
nznx=nz*nx
mgrid=[range(-div(nz,2), step=1.0, length=nz), range(-div(nx,2), step=1.0, length=nx)]
tgrid=range(0.0,step=0.5, length=nt)
@info "Grids are all set."
```

### Allocating medium parameters

```@example test_born
Qv=abs.(randn(nz,nx))
η=abs.(randn(nz,nx))
k=abs.(randn(nz,nx))
σ=abs.(randn(nz,nx))
p=randn(nz,nx,nt)

σobs=abs.(randn(nz,nx))
Qobs=abs.(randn(nz,nx))
nr=10 # number of abstract receivers
ACQ=sprandn(nr,nz*nx,0.6); # choose a random acquisition operator
@info "We are ready for the PoissonExpt."
```

### Create an `Expt`, and then extract a linear forward map out of it

```@example test_born
paE=PoissonExpt(p, tgrid, mgrid, Qv, k, η, σ, ACQ, σobs=σobs, Qobs=Qobs,)
F=operator_Born(paE, σ); # extract the linearized forward operator from `Expt`
GIPh.Utils.test_linearmap(F) # finally do some tests on the linearmap
```

### Usage

```@example test_born
δx=randn(size(F,2)) # random model pertubation
δd=F*δx # corresponding pertubation in data
@info string("Length of data: (nt*nr)=",length(δd))
```

