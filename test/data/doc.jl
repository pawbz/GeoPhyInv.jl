# This page was generated on DATEOFTODAY

using GeoPhyInv
using Test
using Random
#jl @init_parallel_stencil(3, false, Float32, 4)

# # Intro

#md # ```@docs
#md # Records
#md # ```


# # Examples
# Define an acquisition geometry.
mgrid = [range(0.0, stop=10.,step=0.1), range(0.0, stop=30.,step=0.2)];
ageom=AGeom(mgrid, SSrcs(10), Srcs(10), Recs(10));


# Need a time grid.
tgrid=range(0, stop=1.0, step=0.1);

# Lets initialize records for `:P` and `:vz` fields.
data=Records(tgrid, ageom, [:p,:vz]);

# Fill the `:P` field of 3rd supersource with random numbers.
Random.randn!(data[3][:p]);


# # Methods 
# Methods listed for `SrcWav` can used on instances of `Records` too.
