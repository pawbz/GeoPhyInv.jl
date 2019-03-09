


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
include("expt.jl")




nx=20
nz=20
nt=4
nznx=nz*nx
mgrid=[range(0.0,step=0.5, length=nz), range(0.0,step=0.5, length=nz)]
tgrid=range(0.0,step=0.5, length=nt)

p=randn(nznx+1)
p[end]=0.0
Qv=abs.(randn(nz,nx))
η=abs.(randn(nz,nx))
k=abs.(randn(nz,nx))
σ=abs.(randn(nz,nx))
σobs=abs.(randn(nz,nx))
Qobs=abs.(randn(nz,nx))
psi0=randn(nznx)

snaps=randn(nz,nx,nt)

paE=ParamExpt(snaps, tgrid, mgrid, Qv, k, η, σ, σobs=σobs, Qobs=Qobs)

updateLP!(paE, paE.Q)

@time f=func(σ,paE)


f = x->func(x, paE);
gcalc = Calculus.gradient(f);
g3=reshape(gcalc(vec(σ)),nz,nx);


@time mod!(paE,σ,FGσ())
g=paE.g;

@test ≈(g,g3, rtol=1e-5)

