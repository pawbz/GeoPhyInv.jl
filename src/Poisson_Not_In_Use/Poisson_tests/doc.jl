
# 
# This module represents an explicit, direct sparse 2D finite-difference Poisson solver for heterogeneous media,
# i.e. media having spatially varying (space-dependent) medium parameters.
# Current implementation assumes Neumann boundary conditions at all the boundaries.
# 
# Consider the following Poisson experiment:
# ```math
# ∇⋅(σ(x,z)∇) ψ(t) = ∇⋅(Q(x,z)∇) p(t),
# ```
# ```math
# Q = k * Q_v / η.
# ```
