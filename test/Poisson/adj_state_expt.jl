# This module represents an explicit, direct sparse 2D finite-difference Poisson solver for heterogeneous media,
# i.e. media having spatially varying (space-dependent) medium parameters.
# The following functionality is currently available in this module:
# * Apply operator ``A=∇⋅(σ(x,z)∇)`` on a field ``ψ`` to get ``p``.
# * Apply ``A^{-1}`` in order to solve for ``ψ`` in ``Aψ=p``, given ``p``.
# Current implementation assumes Neumann boundary conditions at all the boundaries.

#-

# As a demo, start with loading some packages.

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



# From here on, consider the following Poisson experiment: 
# ```math
# ∇⋅(σ(x,z)∇) ψ(t) = ∇⋅(Q(x,z)∇) p(t),
# ```
# ```math
# Q = k * Qv / η.
# ```
# Dimensions and spatial grids are allocated as follows.

nx=20
nz=20
nt=4
nznx=nz*nx
mgrid=[range(0.0,step=0.5, length=nz), range(0.0,step=0.5, length=nz)]
tgrid=range(0.0,step=0.5, length=nt)

# Now lets allocate the inputs for a toy experiment.

Qv=abs.(randn(nz,nx))
η=abs.(randn(nz,nx))
k=abs.(randn(nz,nx))
σ=abs.(randn(nz,nx))
p=randn(nz,nx,nt)

# These medium parameters are used to generate the *observed* field ``ψ``.

σobs=abs.(randn(nz,nx))
Qobs=abs.(randn(nz,nx))

# Generate the configuration for the Poisson experiment, using the arrays above.

paE=J.Poisson.ParamExpt(p, tgrid, mgrid, Qv, k, η, σ, σobs=σobs, Qobs=Qobs)

#src # J.Poisson.updateLP!(paE, paE.Q)

# Calculate the data misfit.

@time f=J.Poisson.func(σ,paE)

# Compute the gradient with finite-difference approximation.

f = x->J.Poisson.func(x, paE);
gcalc = Calculus.gradient(f);
g_fd=reshape(gcalc(vec(σ)),nz,nx);

# Compute the gradient using adjoint state method.

@time J.Poisson.mod!(paE,σ,J.Poisson.FGσ())
g=paE.g;

# Check gradient accuracy.

@test ≈(g,g_fd, rtol=1e-5)

