# Toolbox for GeoPhysical Inversion (GPI) 

```julia
using GeoPhyInv
```
The methods in this software numerically solve some 
differential equations of geophysics in Julia.
Currently, the equations within the realm of this package 
include:

The methods within the realm of seismic modeling and inversion return:

* the linearized forward modeling operator `F`, such that 
`Fx` can be computed without explicitly storing the operator matrix (
 see `LinearMaps.jl`);
* the imaging/migration operator `F*`;

These maps are the building blocks of iterative optimization schemes.
* the acoustic 
\$a\otimes b\$

* Forward problem, where the seismic data are generated 
using synthetic Earth models and the acquisition parameters 
corresponding to a seismic experiment.
Forward modeling consists of a finite-difference simulation, followed
by convolutions in the time domain using the source and
receiver filters. The details about our finite-difference 
scheme are given in XXX. 
And the filters corresponding to 
sources and receivers are described in XXX.

* Can perform inversion of synthetic scenarios.
First, the seismic data are modeled as in the forward problem. Then the 
data are used to perform full waveform inversion (FWI). The inverse 
problem estimates
the Earth models and the source and receiver filters 
that resulted from the data.
This task is necessary to test the performance of the inversion algorithm 
in various geological scenarios using different acquisition parameters.

* Read 
the measured seismic field data and parameters from a seismic experiment 
to perform inversion like in the previous task. 
The data measured in the field are not in a suitable format 
yet for this software. 
Pre-processing is necessary before it can be used as described.
Also, the acquisition parameters from the field
should be 
converted to suitable 2-D coordinates as described.


## Coding Conventions

This software is organised into various modules.
Each module has various type definitions and methods declared.
Commenting is 
done inside each module file to describe the purpose of the module and its usage.
Most of the comments in the code are inline with the text in this documentation. 
Within each module, variable naming is done 
to reduce the effort needed to understand the source code.
For example, 
`distance = velocity * time` is prefered over using 
`a = b * c` in most parts of the software.
The code inside each method is properly intended using spaces to facilitate 
redability.
We followed this documentation [![guide](http://www.stochasticlifestyle.com/finalizing-julia-package-documentation-testing-coverage-publishing/)](http://www.stochasticlifestyle.com/finalizing-julia-package-documentation-testing-coverage-publishing/)


* The methods ending with `!` ideally should not allocate and memory. They are supposed to be fast and iteratively called inside loops.

## Installation




[![Build Status](https://travis-ci.org/pawbz/JuMIT.jl.svg?branch=master)](https://travis-ci.org/pawbz/JuMIT.jl)

[![codecov.io](http://codecov.io/github/pawbz/JuMIT.jl/coverage.svg?branch=master)](http://codecov.io/github/pawbz/JuMIT.jl?branch=master)
