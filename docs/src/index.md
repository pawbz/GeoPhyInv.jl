# The Expt Types

The methods in this package numerically solve some
differential equations commonly faced in geophysical inverse problems.
The functionality of this package revolves around the mutable `Expt` types.
Julia's multiple dispatch is used to overload `Base` methods whenever possible. 
Which means, if one is familiar with the `Base` methods of Julia, then (almost) no additional syntax is required 
to use this package. Easy!


While performing a given experiment,
firstly, most of the memory necessary
is allocated while creating the `Expt` variables.
Then these variables are input to in-place functions 
(e.g., `update!`)  which as per Julia convention ends with an exclamation mark, to actually perform the experiment task.
For example, the
current `Expt` types within the realm of this package include:

* `SeisForwExpt` is the seismic (both acoustic and elastic) forward modeling experiment;
* `SeisInvExpt` is the type for seismic inversion experiment, including migration;
* `PoissonExpt` is type for the solving the Poisson experiment.

Some of the commonly used (and exported) mutable types to create the `Expt` variables are:

* `Medium` for bundling medium parameters together;
* `AGeom` stores acquisition geometry related parameters;
* `SrcWav` allocates source signals input to an experiment;
* `Records` allocated the output records that are fitted during inversion.

### Loading the package
It is necessary to configure GeoPhyInv with a macro `@init_parallel_stencil` before using it.
If you need to change this configuration, julia kernel must be restarted.
```julia
using GeoPhyInv # load package (after installation)
```
```@docs
GeoPhyInv.@init_parallel_stencil
```
After this, to get started with modelling or inversion, as an example, simply load a seismic inversion experiment already defined in our package gallery into REPL:
```julia
pa=SeisInvExpt(FdtdAcoustic(), LS(), :pizza); # "pizza" is the name of the experiment
```
Then, simply use `update!` to perform least-squares inversion.
```julia
update!(pa, solver=:ipopt)
```

### Fields

For a given `attrib_mod` and `ndims`, you can always print a list of scalar fields that are active. 
The idea is that these symbols will be used 
* to generate mutable types like `SrcWav` or `Records`;
* as input to keywords like `rfields` or `snaps_field`.
First, lets look at the list of all the scalar fields defined in this package.
````@example
Fields()
````
We can always filter this list depending on the type of the `SeisForwExpt` i.e., `attrib_mod` and the number of dimensions.
For example, lets check out the fields for 3-D acoustic simulation.
````@example
Fields(FdtdAcoustic(), ndims=3)
````
Similarly, for 2-D elastic simulation.
````@example
Fields(FdtdElastic(), ndims=2)
````

### Grids

It is necessary to input the evenly-spaced spatial and temporal grids while creating the `Expt` variables.
These grids can be simply created using `Base.range` in Julia, as shown below.

```julia
zgrid=range(0,stop=1000.0,length=201) # create vertical grid from 0 to 1000 m
xgrid=range(0,stop=1000.0,length=201) # similarly, create horizontal grid
```
While constructing spatial-grid bundle, watch out for order.
```julia
mgrid=[zgrid, xgrid] # 2D
mgrid=[zgrid, ygrid, xgrid] # 3D
@info string("spatial sampling intervals (dz, dy, dx)=", step.(mgrid))
```
Similarly, a temporal grid.
```julia
tgrid=range(0,stop=1.0,step=0.001) # a temporal grid from 0 to 1.0 s
```
