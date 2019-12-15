using GeoPhyInv
using SparseArrays
using StatsBase
using LinearAlgebra
using Random
using ProgressMeter
using LinearAlgebra
using Test
using ForwardDiff
using Calculus
using LinearMaps

nx=21
nz=21
nt=4
nznx=nz*nx
mgrid=[range(-div(nz,2), step=1.0, length=nz), range(-div(nx,2), step=1.0, length=nx)]
tgrid=range(0.0,step=0.5, length=nt)
@info "Grids are all set."

Qv=abs.(randn(nz,nx))
η=abs.(randn(nz,nx))
k=abs.(randn(nz,nx))
σ=abs.(randn(nz,nx))
p=randn(nz,nx,nt)

σobs=abs.(randn(nz,nx))
Qobs=abs.(randn(nz,nx))
nr=10 # number of abstract receivers
ACQ=sprandn(nr,nz*nx,0.6); # choose a random acquisition operator
@info "We are ready for the PoissonExpt."

paE=PoissonExpt(p, tgrid, mgrid, Qv, k, η, σ, ACQ, σobs=σobs, Qobs=Qobs,)
F=LinearMap(paE, σ); # extract the linearized forward operator from `Expt`
GeoPhyInv.Utils.test_linearmap(F) # finally do some tests on the linearmap

δx=randn(size(F,2)) # random model pertubation
δd=F*δx # corresponding pertubation in data
@info string("Length of data: (nt*nr)=",length(δd))

