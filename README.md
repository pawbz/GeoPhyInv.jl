[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://pawbz.github.io/GeoPhyInv.jl/dev)
[![Build Status](https://travis-ci.org/pawbz/GeoPhyInv.jl.svg?branch=master)](https://travis-ci.org/pawbz/GeoPhyInv.jl)
[![codecov](https://codecov.io/gh/pawbz/GeoPhyInv.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/pawbz/GeoPhyInv.jl)
# Toolbox for Geophysical Modeling and Inverse Problems


## Installation
At the moment, `GeoPhyInv` depends on another unregistered package `Misfits`. For complete installation,
enter these package manager commands in the REPL:
```julia
using Pkg
Pkg.add(PackageSpec(name="Misfits",url="git@github.com:pawbz/Misfits.jl.git"))
Pkg.add(PackageSpec(name="GeoPhyInv",url="git@github.com:pawbz/GeoPhyInv.jl.git"))
```
## Overview
### Seismic

### Poisson
applies to electrical, gravity and magnetic fields

## Tutorial Notebooks


### Forward Modeling
* generate seismic data using a 2D finite-difference solver
  * [notebook] (https://github.com/pawbz/GeoPhyInv.jl/tree/master/notebooks/modeling/page1.ipynb)


* saving time snapshots during modelling
  * [notebook] (https://github.com/pawbz/GeoPhyInv.jl/tree/master/notebooks/modeling/page2.ipynb)

## Credits
* Charles Clerget [@chclerget](https://github.com/chclerget) tested most of the methods of this package.
* The Poisson solver (`GeoPhyInv.Poisson`) was developed by Niels Grobbe, after adapting scripts from [Aime Fournier](https://erlweb.mit.edu/users/aimemitedu).
* Finally, thanks of [Earth Resources Laboratory](https://erlweb.mit.edu), MIT. The project was supported by Aime Fournier via research funds from Equinor.
