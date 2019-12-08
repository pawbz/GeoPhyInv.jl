
using BenchmarkTools, Test

grid=range(rand(), stop=randn(), length=10)


d=NamedD(grid,Srcs(10), [:a,:b, :c]);
d1=deepcopy(d);
	v=vec(d)

	@btime copyto!(v,d)

	randn!(v)

	@btime copyto!(d,v)
	@btime copyto!(d1,v)

	@test isequal(d,d1)
	@btime reverse!(d1)
	@btime reverse!(d1)
	@test isequal(d,d1)

	@btime randn!(d)
	@btime randn!(d1)

	@btime dot(d,d1)
	@test dot(d,d1)==dot(d1,d)

	@btime fill!(d, 0.0)
	@test iszero(d)


	println("##################################")

d=VNamedD(grid,SSrcs(10), Srcs(10), [:a,:b, :c]);
d1=deepcopy(d);
	v=vec(d)

	@btime copyto!(v,d)

	randn!(v)

	@btime copyto!(d,v)
	@btime copyto!(d1,v)

	@test isequal(d,d1)
	@btime reverse!(d1)
	@btime reverse!(d1)
	@test isequal(d,d1)

	@btime randn!(d)
	@btime randn!(d1)

	@btime dot(d,d1)
	@test dot(d,d1)==dot(d1,d)

	@btime fill!(d, 0.0)
	@test iszero(d)



