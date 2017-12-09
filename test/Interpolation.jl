using Revise
using JuMIT
using Base.Test
using BenchmarkTools

# also testing behaviour when out of bounds!!
nx=Array(linspace(1,20,10));
mx=Array(linspace(3,25,20));
nz=Array(linspace(1,20,30));
mz=Array(linspace(5,27,40));

for Battrib in [:B1, :B2]

	println("=====================================")
	## 1D
	ny=randn(length(nx));
	# interpolate (my=f(ny))
	my=zeros(length(mx));
	@time JuMIT.Interpolation.interp_spray!(nx, ny, mx, my, :interp,
			Battrib)

	myp = randn(length(my));

	# spray (nyp=fˣ(myp))
	nyp=zeros(length(nx));
	@btime JuMIT.Interpolation.interp_spray!(mx, myp, nx, nyp, :spray, Battrib)

	# dot product test
	@test dot(my, myp) ≈ dot(ny, nyp)



	## 2D
	ny=randn(length(nz), length(nx));
	# interpolate
	my=zeros(length(mz), length(mx));
	@time JuMIT.Interpolation.interp_spray!(nx,nz, ny, mx,mz, my, :interp, Battrib)

	myp = randn(size(my));

	# spray
	nyp=zeros(length(nz), length(nx));
	@time JuMIT.Interpolation.interp_spray!(mx,mz, myp, nx,nz, nyp, :spray, Battrib)

	# dot product test
	@test dot(my, myp) ≈ dot(ny, nyp)
end



yout=randn(1000);
x=randn(1000);
@btime JuMIT.Interpolation.interp_B1_1D(4, yout, 3., 4., 5., 6., 3.5)
