# The Expt Datatypes

The methods in this package numerically solve some
differential equations commonly faced in geophysical inverse problems.
The functionality of this package revolves around the mutable `Expt` types.
Julia's multiple dispatch is used to overload `Base` methods whenever possible. 
Using all of a function's arguments to choose which method should be invoked, rather than just the first, is known as multiple dispatch.
Firstly, most of the memory necessary to perform a given experiment
is allocated while creating the `Expt` variables.
Then these variables are input to in-place functions 
(e.g., `mod!`)  which as per Julia convention ends with an exclamation mark, to actually perform the experiment task.
For example, the
current `Expt` types within the realm of this package include:

* `SeisForwExpt` is the seismic (acoustic) forward modeling experiment  ;
* `SeisInvExpt` is the type for seismic inversion experiment, including migration;
* `PoissonExpt` is type for the solving the Poisson experiment.

To get started, as an example, simply load a seismic inversion experiment already defined in our package gallery into REPL:

```julia
using GeoPhyInv # load GIPh (after installation)
paE=GIPh.Gallery.SeisInvExpt(:pizza); # "pizza" is the name of the experiment
```


### Grids

It is necessary to input the evenly-spaced spatial and temporal grids while creating the `Expt` variables.
These grids can be simply created using `Base.range` in Julia, as shown below.

```julia
zgrid=range(0,stop=1000.0,length=201) # create vertical grid from 0 to 1000 m
xgrid=range(0,stop=1000.0,length=201) # create horizontal grid
mgrid=[zgrid, xgrid] # spatial-grid bundle
@info string("spatial sampling intervals (dz,dx)=", step.(mgrid))
tgrid=range(0,stop=1.0,step=0.001) # a temporal grid from 0 to 1.0 s
```
