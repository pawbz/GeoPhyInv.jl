```@meta
EditURL = "<unknown>/test/media/doc.jl"
```

````@example doc
using BenchmarkTools
using GeoPhyInv
using Test
using Plots
gr()
````

# Intro

```@docs
GeoPhyInv.Medium
Medium(::Symbol, ::Real)
```

# Examples

Load a predefined model.

````@example doc
medium=Medium(:elastic_homo1);
nothing #hide
````

Get the medium parameters that are stored.

````@example doc
names(medium)
````

Look at the reference values.

````@example doc
medium.ref
````

Inspect the bounds.

````@example doc
medium.bounds
````

Plotting is easy

````@example doc
p1=heatmap(medium, :vp)
plot(p1,size=(800,400))
````

To construct an instance of 2-D `Medium`, we need ranges for `z` and `x`.

````@example doc
mgrid = [range(0.0, stop=10.,step=0.1), range(0.0, stop=30.,step=0.2)];
nothing #hide
````

For 3-D, it is `z`, `y` and `x`.

````@example doc
mgrid = fill(range(-10, stop=10.,step=0.1), 3)
````

Allocate basic medium parameters on the grid.

````@example doc
medium = Medium(mgrid, [:vp,:rho,:vs]);
nothing #hide
````

Bounds for these parameters should be input for modeling or inversion. Use `update!`

````@example doc
vpb=[2100.,2200.]; vsb=[1500, 1700]; rhob=[2100., 2300.];
update!(medium, [:vp,:vs,:rho], [vpb, vsb, rhob]);
nothing #hide
````

Just fill `medium` with average (reference) values.

````@example doc
fill!(medium);
nothing #hide
````

Once the basic medium parameters are input, we can access some other derived parameters.

````@example doc
medium[:Zp];
nothing #hide
````

Otherwise, we can manually update parameters of `medium`.

````@example doc
medium[:vp].=3000.;
medium[:vs].=2000.;

println(medium)
````

A model can be also be updated by adding random noise.

````@example doc
update!(medium, [:vp,:rho], randn_perc=1.);
nothing #hide
````

# Methods
```@docs

````@example doc
md # Base.getindex(::Medium, ::Symbol)
````

GeoPhyInv.update!(::GeoPhyInv.Medium, ::Vector{Symbol},)
Base.copyto!(x::AbstractArray, medium::Medium, fields::Vector{Symbol})
Base.vec(medium::Medium, ::Symbol)
```

