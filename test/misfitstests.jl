using SIT
using Base.Test

x = randn(10,1); α = randn(1); y = α[1] .* x;
J, α1 = SIT.Misfits.error_after_scaling(x,y) 
@test_approx_eq α1.*x y
#@test_approx_eq J 0.0
