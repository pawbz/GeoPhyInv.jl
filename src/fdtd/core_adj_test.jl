
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
	ntimes=1

#	p=zeros(nz,nx,3)
#	p[div(nz,2),div(nx,2),1]=1.0
	p=Random.randn(nz,nx,3)
	modrhovzI=Random.randn(nz,nx)
	modrhovxI=Random.randn(nz,nx)
	modK=Random.randn(nz,nx)
	a_x=Random.randn(nz,nx)


	pout=advance_sample(ntimes, p, modrhovxI, modrhovzI, modK)
#	println(pout[:,:,1])


	pout2=randn(nz,nx,3)


	p2=advance_sample(ntimes, pout2, modrhovxI, modrhovzI, modK)

	i=3
	@test dot(pout2[:,:,i],pout[:,:,i]) ≈ dot(p[:,:,i],p2[:,:,i])
	#@test dot(pout2,pout-p) ≈ dot(p,p2-pout2)
end
	






