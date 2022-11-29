
using Random
using Test
using GeoPhyInv
using LinearAlgebra

const JF=GeoPhyInv.Fdtd

function Δxtest(nx)
	(nx<6) && error("test requires nx>6")

	p=Random.randn(nx)
	p[1:2] .= 0.0  # it is necessary that starting models have these boundary conditions
	p[nx-1:nx] .= 0.0

	pout=JF.Δxsample(p)

	pout2=randn(nx)
	pout2=pout 

	p2=JF.Δxsample(pout2)


	@test dot(pout2,pout) ≈ dot(p,p2)
end

Δxtest(11)
