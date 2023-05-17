[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://pawbz.github.io/GeoPhyInv.jl/)
# Toolbox for Geophysical Modeling and Inverse Problems

GeoPhyInv provides architecture-agnostic elastic and acoustic wave equation solvers using either [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl) for [Base.Threads](https://docs.julialang.org/en/v1/base/multi-threading/) for high-performance computations on GPUs and CPUs, respectively.  For GPU computations, a performance similar to CUDA C is achieved, thanks to 
[ParallelStencil.jl](https://github.com/omlins/ParallelStencil.jl).
The finite-difference simulations are performed in both 2-D and 3-D
using a staggered-grid velocity-stress formulation.
Finally, [distributed computing](https://docs.julialang.org/en/v1/manual/distributed-computing/) shipped with Julia ensures that the modelling of the super-sources can be parallelized.

## Installation
For complete installation, enter these package manager commands in the REPL:
```julia
using Pkg
Pkg.add(PackageSpec(name="GeoPhyInv",url="https://github.com/pawbz/GeoPhyInv.jl.git"))
```

## Credits and References
* Some implementation ideas are borrowed from Jan Thorbecke's [fdelmodc](https://janth.home.xs4all.nl/Software/Software.html) software.
* [[paper]](https://library.seg.org/doi/abs/10.1190/1.1442147) P-SV wave propagation in heterogeneous media: Velocity‐stress finite‐difference method.
* Work of [Komatitsch and Martin (2007)](https://www.researchgate.net/publication/47503800_An_unsplit_convolutional_Perfectly_Matched_Layer_improved_at_grazing_incidence_for_the_seismic_wave_equation) on convolutional perfectly matched layers for seismic wave equation.
* Charles Clerget [@chclerget](https://github.com/chclerget) tested some methods of this package.
* The Poisson solver (`GeoPhyInv.Poisson`) was developed by Niels Grobbe, after adapting scripts from [Aime Fournier](https://erlweb.mit.edu/users/aimemitedu).
* Thanks to [Earth Resources Laboratory](https://erlweb.mit.edu), MIT. A few developments of this project were supported by Aime Fournier via research funds from Equinor.
