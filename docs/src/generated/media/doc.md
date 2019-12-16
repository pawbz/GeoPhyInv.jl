```@meta
EditURL = "<unknown>/media/doc.jl"
```

```@example doc
using BenchmarkTools
using GeoPhyInv
using Test
using Plots
gr()
```

# Intro

```@docs
GeoPhyInv.Medium
Medium(::Symbol, ::Real)
```

# Examples

Load a predefined model.

```@example doc
mod=Medium(:marmousi2);
nothing #hide
```

Get the medium parameters that are stored.

```@example doc
names(mod)
```

Look at the reference values.

```@example doc
mod.ref
```

Inspect the bounds.

```@example doc
mod.bounds
```

Plotting is easy

```@example doc
p1=heatmap(mod, :vp)
plot(p1,size=(800,400))
```

To construct an instance of `Medium`, we need a 2-D grid.

```@example doc
mgrid = [range(0.0, stop=10.,step=0.1), range(0.0, stop=30.,step=0.2)];
nothing #hide
```

Allocate basic medium parameters on the grid.

```@example doc
mod = Medium(mgrid, [:vp,:rho,:vs]);
nothing #hide
```

Bounds for these parameters should be input for modeling or inversion. Use `update!`

```@example doc
vpb=[2100.,2200.]; vsb=[1500, 1700]; rhob=[2100., 2300.];
update!(mod, [:vp,:vs,:rho], [vpb, vsb, rhob]);
nothing #hide
```

Just fill `mod` with average (reference) values.

```@example doc
fill!(mod);
nothing #hide
```

Once the basic medium parameters are input, we can access some other derived parameters.

```@example doc
mod[:Zp];
nothing #hide
```

Otherwise, we can manually update parameters of `mod`.

```@example doc
mod[:vp].=3000.;
mod[:vs].=2000.;
nothing #hide
```

A model can be also be updated by adding random noise.

```@example doc
update!(mod, [:vp,:rho], randn_perc=1.);
nothing #hide
```

# Methods
```@docs
Base.getindex(::Medium, ::Symbol)
GeoPhyInv.update!(::GeoPhyInv.Medium, ::Vector{Symbol},)
Base.copyto!(x::AbstractArray, mod::Medium, fields::Vector{Symbol})
Base.vec(mod::Medium, ::Symbol)
```

