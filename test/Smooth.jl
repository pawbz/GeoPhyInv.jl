using SIT
using Base.Test


x1=randn(20)
x2=randn(20)
x3=copy(x1+x2)

SIT.Smooth.gaussian!(x1,[3])
SIT.Smooth.gaussian!(x2,[3])

SIT.Smooth.gaussian!(x3,[3])


@test_approx_eq (x1+x2) x3



x1=randn(20,20)
x2=randn(20,20)
x3=copy(x1+x2)

SIT.Smooth.gaussian!(x1,[3,5])
SIT.Smooth.gaussian!(x2,[3,5])

SIT.Smooth.gaussian!(x3,[3,5])


@test_approx_eq (x1+x2) x3




