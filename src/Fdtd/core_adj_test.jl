
function Δxtest(nx)
	(nx<6) && error("test requires nx>6")

	p=Random.randn(nx)
	p[1:2] .= 0.0  # it is necessary that starting models have these boundary conditions
	p[nx-1:nx] .= 0.0

	pout=Δxsample(p)

	pout2=randn(nx)
	pout2[1:2] .= 0.0  # it is necessary that starting models have these boundary conditions
	pout2[nx-1:nx] .= 0.0

	p2=Δxsample(pout2)


	@test dot(pout2,pout) ≈ dot(p,p2)
end


function Δttest(nx)
	(nx<6) && error("test requires nx>6")

	p=Random.randn(nx)

	pout=Δtsample(p)
	println(pout)

	pout2=randn(nx)

	p2=Δtsample(pout2)


	@test dot(pout2,pout) ≈ dot(p,p2)
end





function testing2(nx,nz)
	ntimes=3

	p=zeros(nz,nx)
	p=Random.randn(nz,nx)
	modrrvz=Random.randn(nz,nx)
	modrrvx=Random.randn(nz,nx)
	a_x=Random.randn(nz,nx)


	pout=advance_sample(ntimes, p, modrrvx, modrrvz)


	pout2=randn(nz,nx)


	p2=advance_sample(ntimes, pout2, modrrvx, modrrvz)

	@test dot(pout2,pout) ≈ dot(p,p2)
end
	






