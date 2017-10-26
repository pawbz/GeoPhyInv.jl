using JuMIT
using Base.Test

# also testing behaviour when out of bounds!!
nx=Array(linspace(1,20,10));
mx=Array(linspace(3,25,20)); 
nz=Array(linspace(1,20,30));
mz=Array(linspace(5,27,40)); 

for i in 1:3
	println("=====================================")
	for Battrib in [:B1, :B2]

		## 1D
		ny=randn(length(nx));
		# interpolate (my=f(ny))
		my=zeros(length(mx));
		@time JuMIT.Interpolation.interp_spray!(nx, ny, mx, my, :interp, Battrib)

		myp = randn(length(my));

		# spray (nyp=fˣ(myp))
		nyp=zeros(length(nx));
		@time JuMIT.Interpolation.interp_spray!(mx, myp, nx, nyp, :spray, Battrib)

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
end


