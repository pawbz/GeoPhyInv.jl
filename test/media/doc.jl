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
medium = Medium(:elastic_homo1);

# Get the medium parameters that are stored.
names(medium)

# Look at the reference values.
medium.ref

# Inspect the bounds.
medium.bounds

# Plotting is easy
#md p1=heatmap(medium, :vp)
#md plot(p1,size=(800,400))


# To construct an instance of 2-D `Medium`, we need ranges for `z` and `x`.
mgrid = [range(0.0, stop = 10.0, step = 0.1), range(0.0, stop = 30.0, step = 0.2)];


# For 3-D, it is `z`, `y` and `x`.
mgrid = fill(range(-10, stop = 10.0, step = 0.1), 3)

# Allocate basic medium parameters on the grid.
medium = Medium(mgrid, [:vp, :rho, :vs]);

# Bounds for these parameters should be input for modeling or inversion. Use `update!`
vpb = [2100.0, 2200.0];
vsb = [1500, 1700];
rhob = [2100.0, 2300.0];
update!(medium, [:vp, :vs, :rho], [vpb, vsb, rhob]);

# Just fill `medium` with average (reference) values.
fill!(medium);

# Once the basic medium parameters are input, we can access some other derived parameters.
medium[:Zp];

# Otherwise, we can manually update parameters of `medium`.
medium[:vp] .= 3000.0;
medium[:vs] .= 2000.0;

println(medium)

# A model can be also be updated by adding random noise.
update!(medium, [:vp, :rho], randn_perc = 1.0);

# # Methods 
#md # ```@docs
#md md # Base.getindex(::Medium, ::Symbol)
#md # GeoPhyInv.update!(::GeoPhyInv.Medium, ::Vector{Symbol},)
#md # Base.copyto!(x::AbstractArray, medium::Medium, fields::Vector{Symbol})  
#md # Base.vec(medium::Medium, ::Symbol)  
#md # ```

