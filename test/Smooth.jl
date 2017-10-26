using JuMIT
using Base.Test


x1=randn(20)
x2=randn(20)
x3=copy(x1+x2)

JuMIT.Smooth.gaussian!(x1,[3])
JuMIT.Smooth.gaussian!(x2,[3])

JuMIT.Smooth.gaussian!(x3,[3])


@test isapprox((x1+x2),x3)



x1=randn(20,20)
x2=randn(20,20)
x3=copy(x1+x2)

JuMIT.Smooth.gaussian!(x1,[3,5])
JuMIT.Smooth.gaussian!(x2,[3,5])

JuMIT.Smooth.gaussian!(x3,[3,5])


@test isapprox((x1+x2),x3)
