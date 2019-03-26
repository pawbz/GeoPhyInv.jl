"""
This module represents an explicit, direct sparse 2D finite-difference Poisson solver for heterogeneous media,
i.e. media having spatially varying (space-dependent) medium parameters.
The following functionality is currently available in this module:
* 'Forward' problem: given the source, and the conductivity distribution, solve for the electrical potential ψ:
	 Example problem: solve Aψ=∇⋅j; A=∇⋅(σ(x,z)∇), for ψ.
* 'Inverse' problem: given wave propagation related pore pressure, and mechanical medium properties, calculate the source term:  
	 Example problem: solve Aψ=-∇⋅j; A=∇⋅([Q*k/η](x,z)∇B*P), for ∇⋅j;
* The Boundary Value Problem that is currently implemented assumes Neumann boundary conditions at all boundaries.
Developed by:\\\n
Niels Grobbe, Massachusetts Institute of Technology, USA.\\\n
In collaboration with: Aimé Fournier & Laurent Demanet, Massachusetts Institute of Technology, USA.\\\n
Date: October, 2017 \\\n
Contact: ngrobbe@gmail.com
"""
module Poisson

include("core.jl")
include("expt.jl")
include("expt_born.jl")
include("getprop.jl")

end # end module Poisson
