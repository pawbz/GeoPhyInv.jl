[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://pawbz.github.io/GeoPhyInv.jl/dev)
[![Build Status](https://travis-ci.org/pawbz/GeoPhyInv.jl.svg?branch=master)](https://travis-ci.org/pawbz/GeoPhyInv.jl)
[![codecov](https://codecov.io/gh/pawbz/GeoPhyInv.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/pawbz/GeoPhyInv.jl)
# Toolbox for Geophysical Modeling and Inverse Problems


## Installation
At the moment, `GeoPhyInv` depends on another unregistered package `Misfits`. For complete installation,
enter these package manager commands in the REPL:
```julia
using Pkg
Pkg.add(PackageSpec(name="Misfits",url="https://github.com/pawbz/Misfits.jl.git"))
Pkg.add(PackageSpec(name="GeoPhyInv",url="https://github.com/pawbz/GeoPhyInv.jl.git"))
```

## Credits
* Charles Clerget [@chclerget](https://github.com/chclerget) tested most of the methods of this package.
* The Poisson solver (`GeoPhyInv.Poisson`) was developed by Niels Grobbe, after adapting scripts from [Aime Fournier](https://erlweb.mit.edu/users/aimemitedu).
* Finally, thanks of [Earth Resources Laboratory](https://erlweb.mit.edu), MIT. The project was supported by Aime Fournier via research funds from Equinor.
