using Revise
using JuMIT
using Base.Test

# real 1D
x = randn(10,1); α = randn(1); y = α[1] .* x;
J, α1 = JuMIT.Misfits.error_after_scaling(x,y)
@test isapprox(α1.*x, y)
#@test_approx_eq J 0.0

# real 2D
x = randn(10,10); α = randn(1); y = α[1] .* x;
J, α1 = JuMIT.Misfits.error_after_scaling(x,y)
@test isapprox(α1.*x, y)

# complex 2D
x = complex.(randn(10,10), randn(10,10));
α = complex.(randn(1), randn(1)); y = α[1] .* x;
J, α1 = JuMIT.Misfits.error_after_scaling(x,y)
@test isapprox(α1.*x, y)


# error invariant of translation or global phase
x=randn(100); y=randn() .* circshift(x,20)
J, α = JuMIT.Misfits.error_after_autocorr_scaling(x,y)
@test J < 1e-15


x=randn(100,10);
dfdx1=similar(x);
JuMIT.Misfits.error_pairwise_corr_dist(dfdx1,x)
xvec=vec(x)
dfdx2=similar(xvec);
JuMIT.Inversion.finite_difference!(x -> JuMIT.Misfits.error_pairwise_corr_dist(nothing, reshape(x,100,10)), xvec, dfdx2, :central)

@test dfdx1 ≈ reshape(dfdx2,100,10)


x=randn(10,10);
dfdx1=similar(x);
JuMIT.Misfits.error_autocorr_pairwise_corr_dist(dfdx1,x)
xvec=vec(x)
dfdx2=similar(xvec);
JuMIT.Inversion.finite_difference!(x -> JuMIT.Misfits.error_autocorr_pairwise_corr_dist(nothing, reshape(x,10,10)), xvec, dfdx2, :central)

@test dfdx1 ≈ reshape(dfdx2,10,10)
