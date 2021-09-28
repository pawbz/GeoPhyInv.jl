using BenchmarkTools
using GeoPhyInv
using Test
#md using Plots
#md gr()


# # Intro

#md # ```@docs
#md # GeoPhyInv.Medium
#md # Medium(::Symbol, ::Real) 
#md # ```

# # Examples

# Load a predefined model.
mod=Medium(:marmousi2);

# Get the medium parameters that are stored.
names(mod)

# Look at the reference values.
mod.ref

# Inspect the bounds.
mod.bounds

# Plotting is easy
#md p1=heatmap(mod, :vp)
#md plot(p1,size=(800,400))


# To construct an instance of 2-D `Medium`, we need ranges for `z` and `x`.
mgrid = [range(0.0, stop=10.,step=0.1), range(0.0, stop=30.,step=0.2)];


# For 3-D, it is `z`, `y` and `x`.
mgrid = fill(range(-10, stop=10.,step=0.1), 3)

# Allocate basic medium parameters on the grid.
mod = Medium(mgrid, [:vp,:rho,:vs]);

# Bounds for these parameters should be input for modeling or inversion. Use `update!`
vpb=[2100.,2200.]; vsb=[1500, 1700]; rhob=[2100., 2300.];
update!(mod, [:vp,:vs,:rho], [vpb, vsb, rhob]);

# Just fill `mod` with average (reference) values.
fill!(mod);

# Once the basic medium parameters are input, we can access some other derived parameters.
mod[:Zp];

# Otherwise, we can manually update parameters of `mod`.
mod[:vp].=3000.;
mod[:vs].=2000.;

println(mod)

# A model can be also be updated by adding random noise.
update!(mod, [:vp,:rho], randn_perc=1.);

# # Methods 
#md # ```@docs
#md md # Base.getindex(::Medium, ::Symbol)
#md # GeoPhyInv.update!(::GeoPhyInv.Medium, ::Vector{Symbol},)
#md # Base.copyto!(x::AbstractArray, mod::Medium, fields::Vector{Symbol})  
#md # Base.vec(mod::Medium, ::Symbol)  
#md # ```

