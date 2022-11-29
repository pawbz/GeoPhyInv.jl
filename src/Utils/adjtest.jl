
# testing a linear operator

function test_linearmap(F)
	for FF in [F, transpose(F)]
		# testing if (F) is linear
		x1=randn(size(FF,2)) 
		x2=randn(size(FF,2)) 
		x12=x1.+x2
		d12=FF*x12
		d1=FF*x1
		d2=FF*x2
		d12new=d1.+d2
		@test squared_euclidean!(nothing, d12, d12new, nothing, norm_flag=true)<1e-25
	end

	adjtest(F)

	return nothing
end


function adjtest(F)
	x=randn(size(F,2))
	y=randn(size(F,1))
	a=LinearAlgebra.dot(y,F*x)
	b=LinearAlgebra.dot(x,adjoint(F)*y)
	c=LinearAlgebra.dot(x, transpose(F)*F*x)
	println("must be positive: ", c)
	println("adjoint test: ", a, "\t", b)       
	@test isapprox(a,b,rtol=1e-5)
	@test c>0.0
	return nothing
end





