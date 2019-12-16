using Test
using GeoPhyInv
using LinearMaps
using Misfits
using LinearAlgebra

pa=SeisInvExpt(FdtdBorn(), LS(), :pizza, rfields=[:P, :Vz]); # use linearized modeling attribute
F=LinearMap(pa); # generate LinearMap

x1=randn(size(F,2)); # random input 1
x2=randn(size(F,2)); # random input 2
x12=x1.+x2; # sum

d12=F*x12;
d1=F*x1; #
d2=F*x2;
d12new=d1.+d2;

@test Misfits.error_squared_euclidean!(nothing, d12, d12new, nothing, norm_flag=true)<1e-25

x=randn(size(F,2)); # random x
y=randn(size(F,1)); # random y
a=LinearAlgebra.dot(y,F*x);
b=LinearAlgebra.dot(x,adjoint(F)*y);
c=LinearAlgebra.dot(x, transpose(F)*F*x);

println("adjoint test: ", a, "\t", b)
@test isapprox(a,b,rtol=1e-5)

println("x*transpose(F)*F*x be positive: ", c)
@test c>0.0

