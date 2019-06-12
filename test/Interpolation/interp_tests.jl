using GeoPhyInv
using BenchmarkTools
using LinearAlgebra
using Test


function dptest(nx,nz,mx,mz,Battrib)
	pa=GInterp.Kernel([nx, nz], [mx, mz],Battrib)
	ny=randn(length(nz), length(nx));
	my=zeros(length(mz), length(mx));
	# interp
	@time GInterp.interp_spray!(ny, my, pa, :interp)
	myp = randn(size(my));
	# spray
	nyp=zeros(length(nz), length(nx));
	@time GInterp.interp_spray!(nyp, myp, pa, :spray)
	# dot product test
	@test LinearAlgebra.dot(my, myp) ≈ LinearAlgebra.dot(ny, nyp)
end

@testset "GInterp to a larger grid using  nx==1 or nz==1" begin
	nx=range(1.,step=20.,length=1);
	mx=range(1.,stop=40.,length=10);
	nz=range(1.,stop=20.,length=30);
	mz=range(4.,stop=16.,length=30);

	dptest(nx,nz,mx,mz,:B1)

	nx=range(4.,stop=16.,length=20);
	mx=range(1.,stop=30.,length=30);
	nz=range(1.,step=5.,length=1);
	mz=range(1.,stop=10.,length=66);

	dptest(nx,nz,mx,mz,:B1)
end


@testset "GInterp on same grid for nx==1 or nz==1" begin
	nx=range(1.,stop=1.,length=1);
	mx=range(1.,stop=1.,length=1);
	nz=range(1.,stop=20.,length=30);
	mz=range(4.,stop=16.,length=30);

	dptest(nx,nz,mx,mz,:B1)

	nx=range(4.,stop=16.,length=20);
	mx=range(1.,stop=30.,length=30);
	nz=range(1.,stop=1.,length=1);
	mz=range(1.,stop=1.,length=1);

	dptest(nx,nz,mx,mz,:B1)
end

@testset "no extrapolations" begin
	nx=range(21,stop=100,length=80);
	mx=range(1,stop=100,length=100);
	ny=ones(length(nx));
	c=randn()
	my=c.*ones(length(mx))

	pa=GInterp.P_core([nx], [mx], :B1)

	GInterp.interp_spray!(ny, my, pa, :interp)
	@test all(my[1:20] .== c)

	nx=range(3,stop=7,length=64);
	nz=range(3,stop=7,length=95);
	mx=range(1,stop=10,length=10);
	mz=range(1,stop=10,length=10);
	ny=ones(length(nz), length(nx));
	c=randn()
	my=c.*ones(length(mz), length(mx))

	pa=GInterp.P_core([nx, nz], [mx, mz], :B1)

	GInterp.interp_spray!(ny, my, pa, :interp)

	@test all(my[my .≠1] .== c)

end


# also testing behaviour when out of bounds!!
nx=range(1,stop=20,length=100);
mx=range(3,stop=25,length=200);
nz=range(1,stop=20,length=30);
mz=range(5,stop=27,length=40);

println("=====================================")
for Battrib in [:B1, :B2]
	global nx, mx

	pa=GInterp.Kernel([nx], [mx], Battrib)
	println("=====================================")
	## 1D
	ny=randn(length(nx));
	# interpolate (my=f(ny))
	my=zeros(length(mx));
	@time GInterp.interp_spray!(ny, my, pa, :interp)

	myp = randn(length(my));

	# spray (nyp=fˣ(myp))
	nyp=zeros(length(nx));
	@time GInterp.interp_spray!(nyp, myp, pa, :spray)

	# dot product test
	@test LinearAlgebra.dot(my, myp) ≈ LinearAlgebra.dot(ny, nyp)
end


@testset "GInterp on the same grid should not change for B1" begin
	Battrib=:B1
	pa=GInterp.Kernel([nx, nz], [nx, nz], Battrib)
	println("=====================================")
	ny=randn(length(nz), length(nx));
	my=zeros(length(nz), length(nx));
	@time GInterp.interp_spray!(ny, my, pa, :interp)
	@test ny≈my


	nyvec=vcat(vec(ny),vec(ny))
	myvec=zeros(length(my)*2)
	@time GInterp.interp_spray!(nyvec, myvec, pa, :interp, 2)
	@test nyvec≈myvec



	myp = randn(size(my));
	nyp=zeros(length(nz), length(nx));
	@time GInterp.interp_spray!(nyp, myp, pa, :spray)
	@test LinearAlgebra.dot(my, myp) ≈ LinearAlgebra.dot(ny, nyp)
end


for Battrib in [:B1, :B2]
	global nx, nz, mx, mz
	pa=GInterp.Kernel([nx, nz], [mx, mz], Battrib)
	## 2D
	ny=randn(length(nz), length(nx));
	# interpolate
	my=zeros(length(mz), length(mx));
	@time GInterp.interp_spray!(ny, my, pa, :interp)


	myp = randn(size(my));

	# spray
	nyp=zeros(length(nz), length(nx));
	@time GInterp.interp_spray!(nyp, myp, pa, :spray)

	# dot product test
	@test LinearAlgebra.dot(my, myp) ≈ LinearAlgebra.dot(ny, nyp)


	nmod=3
	ny=randn(length(nz)*length(nx)*nmod);
	# interpolate
	my=zeros(length(mz)*length(mx)*nmod);
	@time GInterp.interp_spray!(ny, my, pa, :interp,nmod)


	myp = randn(size(my));

	# spray
	nyp=zeros(length(nz)*length(nx)*nmod);
	@time GInterp.interp_spray!(nyp, myp, pa, :spray,nmod)

	# dot product test
	@test LinearAlgebra.dot(my, myp) ≈ LinearAlgebra.dot(ny, nyp)

end

