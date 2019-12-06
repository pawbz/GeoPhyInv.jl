using BenchmarkTools
using GeoPhyInv
using Test
using Gadfly


# To construct a variable to type `Medium`, the first step is to create a 2-D grid
mgrid = [range(0.0, stop=10.,step=0.1), range(0.0, stop=30.,step=0.2)];

# Initiate the storage of medium parameters on the grid using...
mod = Medium(mgrid);

# By default certain parameters are populated, see
names(mod)

# To define parameters yourself
mod = Medium(mgrid, [:vp,:rho,:vs])

# Bounds for these parameters should be input for modeling or inversion. Use `update!`
vpb=[2100.,2200.]; vsb=[1500, 1700]; rhob=[2100., 2300.]
update!(mod, [:vp,:vs,:rho], [vpb, vsb, rhob]);

# To just fill the Medium with average values, just do
fill!(mod)

# Otherwise, to manually fill in different parameters
mod[:vp].=3000. 
mod[:vs].=2000.

# In order to add random noise to the models 
update!(mod, [:vp,:rho], randn_perc=1.)

# Some plotting #1
spy(mod[:vp], Guide.xlabel("x"), Guide.ylabel("z"))


# Some plotting #2
spy(mod[:vs], Guide.xlabel("x"), Guide.ylabel("z"))

