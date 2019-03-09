using JuMIT
using Revise
using Calculus
using Test
#include("core.jl")

nx=5
nz=7
nznx=nz*nx

p=randn(nznx)
σ=randn(nz,nx)
Q=randn(nz,nx)

function func_test(x)
	x=reshape(x,nz,nx)
	paA=J.Poisson.Param(nx,nz,x)	

	global p
	psi=(paA.A*p)

	# some dummy functional
	f=sum(abs2,psi)
	return f
end


function grad_test(x)
	global nznx, nz, nx,p
	x=reshape(x,nz,nx)
	paA=J.Poisson.Param(nx,nz,x)	

	psi=paA.A*p

	# gradient, that involves dAdx for testing
	g1= vec(2.0*(paA.A)*p*transpose(p))'*hcat(vec.(paA.dAdx)...)
	return reshape(g1,nz,nx)
end

f = x->func_test(x);

gcalc = Calculus.gradient(f);
g1=reshape(gcalc(vec(Q)),nz,nx);

g2=grad_test(Q)

@test g1≈g2
