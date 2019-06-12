# ### Solve for ``ψ`` in `PoissonExpt`.
# This module represents an explicit, direct sparse 2D finite-difference Poisson solver for heterogeneous media,
# i.e. media having spatially varying (space-dependent) medium parameters.
# Current implementation assumes Neumann boundary conditions at all the boundaries.

# Consider the following Poisson experiment: 
# ```math
# ∇⋅(σ(x,z)∇) ψ(t) = ∇⋅(Q(x,z)∇) p(t),
# ```
# ```math
# Q = k * Q_v / η.
# ```
# We start with the dimensions and spatial grids are allocated as follows.

nx=21
nz=21
nt=4
nznx=nz*nx
mgrid=[range(-div(nz,2), step=1.0, length=nz), range(-div(nx,2), step=1.0, length=nx)]
tgrid=range(0.0,step=0.5, length=nt)
@info "Grids are all set."


# The following functionality is currently available in this module:
# * Apply operator ``A=∇⋅(Q(x,z)∇)`` on a field ``p``.
# * Then apply ``A^{-1}`` in order to solve for ``ψ``.

#-

# Now lets allocate the inputs for a toy experiment.

Qv=abs.(randn(nz,nx))
η=abs.(randn(nz,nx))
k=abs.(randn(nz,nx))
σ=abs.(randn(nz,nx))
p=randn(nz,nx,nt)

# These medium parameters are used to generate the *observed* field ``ψ``.

σobs=abs.(randn(nz,nx))
Qobs=abs.(randn(nz,nx))

#-

# ### Loading some packages.

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
#src # include("core.jl")
#src # include("expt.jl")




# Generate the configuration for the Poisson experiment, using the arrays above.

acqgeom=J.Acquisition.Geom_circ(nss=1,nr=30,rad=[5.,5.])
ACQ=J.Acquisition.ACQmat(acqgeom,mgrid)
paE=J.Poisson.ParamExpt(p, tgrid, mgrid, Qv, k, η, σ, ACQ, σobs=σobs, Qobs=Qobs,)

#src # J.Poisson.updateLP!(paE, paE.Q)

# Calculate the data misfit.

@time f=J.Poisson.func(σ,paE)

# Compute the gradient with finite-difference approximation.

f = x->J.Poisson.func(x, paE);
gcalc = Calculus.gradient(f);
g_fd=reshape(gcalc(vec(σ)),nz,nx);

# Compute the gradient using adjoint state method.

@time mod!(paE,σ,J.Poisson.FGσ())
g=paE.g;

# Check gradient accuracy.

@test ≈(g,g_fd, rtol=1e-5)

