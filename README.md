# Julia Modeling and Inversion Toolbox

## Installation
Install these packages after pressing `]` in Julia:
```julia
<<<<<<< HEAD
add https://github.com/pawbz/Misfits.jl.git
add https://github.com/pawbz/Inversion.jl.git
add https://github.com/pawbz/Interpolation.jl.git
add https://github.com/pawbz/Conv.jl.git
add https://github.com/pawbz/JuMIT.jl.git
=======
using Pkg
Pkg.add("https://github.com/pawbz/JuMIT.jl.git")
>>>>>>> d4892e6afe7e0968dd37cef4629b25563dfb2775
```

## Tutorials

### Forward Modeling
* generate seismic data using 2D finite-difference simulation
solvers
[notebook] (https://github.com/pawbz/JuMIT/notebooks/modeling)


Various tutorials that demonstrate the use of this software are available 
[here](https://github.com/pawbz/JuMITtutorials). Installation if `IJulia` is necessary.


## Documentation
Online documentation of various modules of this package can be found 
[here](https://pawbz.github.io/JuMIT.jl/).


## Credits
* Thanks to Charles Clerget [@chclerget](https://github.com/chclerget) for testing a few methods of this package.



[![Build Status](https://travis-ci.org/pawbz/JuMIT.jl.svg?branch=master)](https://travis-ci.org/pawbz/JuMIT.jl)
[![codecov](https://codecov.io/gh/pawbz/JuMIT.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/pawbz/JuMIT.jl)

