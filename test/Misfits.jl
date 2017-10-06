using JuMIT
using Base.Test

x = randn(10,1); α = randn(1); y = α[1] .* x;
J, α1 = JuMIT.Misfits.error_after_scaling(x,y)
@test_approx_eq α1.*x y
#@test_approx_eq J 0.0

x = randn(10,10); α = randn(1); y = α[1] .* x;
J, α1 = JuMIT.Misfits.error_after_scaling(x,y)
@test_approx_eq α1.*x y


x = complex.(randn(10,10), randn(10,10));
α = complex.(randn(1), randn(1)); y = α[1] .* x;
J, α1 = JuMIT.Misfits.error_after_scaling(x,y)
@test_approx_eq α1.*x y
