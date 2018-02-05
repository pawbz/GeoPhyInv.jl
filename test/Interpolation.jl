using Revise
using JuMIT
using Base.Test
using BenchmarkTools

# also testing behaviour when out of bounds!!
nx=Array(linspace(1,20,100));
mx=Array(linspace(3,25,200));
nz=Array(linspace(1,20,30));
mz=Array(linspace(5,27,40));


Battrib=:B1

pa=JuMIT.Interpolation.Param([nx], [mx], Battrib)
println("=====================================")
## 1D
ny=randn(length(nx));
# interpolate (my=f(ny))
my=zeros(length(mx));
@time JuMIT.Interpolation.interp_spray!(ny, my, pa, :interp)

myp = randn(length(my));

# spray (nyp=fˣ(myp))
nyp=zeros(length(nx));
@time JuMIT.Interpolation.interp_spray!(nyp, myp, pa, :spray)

# dot product test
@test dot(my, myp) ≈ dot(ny, nyp)

"""
Interpolation on the same grid should not change
"""
pa=JuMIT.Interpolation.Param([nx, nz], [nx, nz], Battrib)
println("=====================================")
ny=randn(length(nz), length(nx));
my=zeros(length(nz), length(nx));
@time JuMIT.Interpolation.interp_spray!(ny, my, pa, :interp)
@test ny≈my
myp = randn(size(my));
nyp=zeros(length(nz), length(nx));
@time JuMIT.Interpolation.interp_spray!(nyp, myp, pa, :spray)
@test dot(my, myp) ≈ dot(ny, nyp)

pa=JuMIT.Interpolation.Param([nx, nz], [mx, mz], Battrib)
println("=====================================")

## 2D
ny=randn(length(nz), length(nx));
# interpolate
my=zeros(length(mz), length(mx));
@time JuMIT.Interpolation.interp_spray!(ny, my, pa, :interp)


myp = randn(size(my));

# spray
nyp=zeros(length(nz), length(nx));
@time JuMIT.Interpolation.interp_spray!(nyp, myp, pa, :spray)

# dot product test
@test dot(my, myp) ≈ dot(ny, nyp)


