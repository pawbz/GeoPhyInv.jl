
# testing a linear operator
function linearmap_test(F)
    for FF in [F]
        # testing if (F) is linear
        x1 = randn(size(FF, 2))
        x2 = randn(size(FF, 2))


        x1 = zeros(size(FF, 2))
        x2 = zeros(size(FF, 2))

		p1, q1= chunk(x1, 2)
		p2, q2= chunk(x2, 2)

		N = length(p1)
  N2 = div(N, 2)

  p1[N2:N2+10-1] .= randn(10)
  p2[N2:N2+10-1] .= randn(10)
  q1[N2:N2+10-1] .= randn(10)
  q2[N2:N2+10-1] .= randn(10)


        x12 = x1 .+ x2
        d12 = FF * x12
        d1 = FF * x1
        d2 = FF * x2
        d12new = d1 .+ d2
        @test isapprox(d12, d12new)
    end
end


function adj_test(F)
    x = randn(size(F, 2))
    y = randn(size(F, 1))
    a = LinearAlgebra.dot(y, F * x)
    b = LinearAlgebra.dot(x, adjoint(F) * y)
    c = LinearAlgebra.dot(x, transpose(F) * F * x)
    println("must be positive: ", c)
    println("adjoint test: ", a, "\t", b)
    @test isapprox(a, b, rtol=1e-5)
    @test c > 0.0
    return nothing
end





