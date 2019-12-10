```@meta
EditURL = "<unknown>/media/doc.jl"
```

```@example doc
using BenchmarkTools
using GeoPhyInv
using Test
using Gadfly
import Cairo, Fontconfig
```

To construct a variable to type `Medium`, the first step is to create a 2-D grid

```@example doc
mgrid = [range(0.0, stop=10.,step=0.1), range(0.0, stop=30.,step=0.2)];
nothing #hide
```

Initiate the storage of medium parameters on the grid using...

```@example doc
mod = Medium(mgrid);
nothing #hide
```

By default certain parameters are populated, see

```@example doc
names(mod)
```

To define parameters yourself

```@example doc
mod = Medium(mgrid, [:vp,:rho,:vs])
```

Bounds for these parameters should be input for modeling or inversion. Use `update!`

```@example doc
vpb=[2100.,2200.]; vsb=[1500, 1700]; rhob=[2100., 2300.]
update!(mod, [:vp,:vs,:rho], [vpb, vsb, rhob]);
nothing #hide
```

To just fill the Medium with average values, just do

```@example doc
fill!(mod)
```

Otherwise, to manually fill in different parameters

```@example doc
mod[:vp].=3000.
mod[:vs].=2000.
```

In order to add random noise to the models

```@example doc
update!(mod, [:vp,:rho], randn_perc=1.)
```

Some plotting #1

```@example doc
spy(mod[:vp], Guide.xlabel("x"), Guide.ylabel("z"))  |> SVG("foo.svg")
```

![](foo.svg)

Some plotting #2

```@example doc
draw(SVG("foo.svg", 6inch, 3inch), spy(mod[:vs], Guide.xlabel("x"), Guide.ylabel("z")));
nothing #hide
```

![](foo.svg)

