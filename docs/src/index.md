# The Expt Datatype

The methods in this package numerically solve some
differential equations commonly faced in geophysical inverse problems.
The functionality of this package revolves around the mutable `Expt` types.
Firstly, most of the memory necessary to perform a given experiment
is allocated while creating the `Expt` variables.
Then these variables are input to in-place functions (e.g., `mod!`)  which as per Julia convention ends with an exclamation mark, to actually perform the experiment task.
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


