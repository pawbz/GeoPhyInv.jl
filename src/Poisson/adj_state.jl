
using SparseArrays
using StatsBase
using LinearAlgebra
using Random
using ProgressMeter
using LinearAlgebra
using Test
using ForwardDiff
using Calculus
include("core.jl")


nx=5
nz=5
nznx=nz*nx
ddx=22.5
ddz=10.3

p=randn(nznx+1)
p[end]=0.0
Q=abs.(randn(nz,nx))
σ=abs.(randn(nz,nx))
psi0=randn(nznx)

#paA=Param(nx,nz,σ)	
#paB=Param(nx,nz,Q)
#A=paA.A
#B=paB.A

function func_test(x,p)
	global psi0
	x=reshape(x,nz,nx)
	paA=Param(nx,nz,x)	

	psi=vec(applyinvA(p,paA))
	f= sum(abs2,psi.-psi0)
	return f
end


function grad_test(x,p)
	global psi0
	paA=Param(nx,nz,x)	

	psi=vec(applyinvA(p,paA))

	adjsrc=-2.0*(psi.-psi0)

	lambda=vec(applyinvAt(adjsrc,paA))
	#println(size(lambda))
	
	return vec((lambda*psi'))'*(paA.dAdx)
end

f = x->func_test(x, p);

gcalc = Calculus.gradient(f);
g3=reshape(gcalc(vec(σ)),nz,nx);

g2=reshape(grad_test(σ,p),nz,nx)


@test ≈(g2,g3, rtol=1e-5)
