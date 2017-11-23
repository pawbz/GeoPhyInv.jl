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
@time JuMIT.Misfits.error_pairwise_corr_dist(dfdx1,x)
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

x=randn(100,10);
w=randn(100,10);
dfdx1=similar(x);
@time JuMIT.Misfits.error_weighted_norm!(dfdx1,x,w)
xvec=vec(x)
dfdx2=similar(xvec);
JuMIT.Inversion.finite_difference!(x -> JuMIT.Misfits.error_weighted_norm!(nothing, reshape(x,100,10), w), xvec, dfdx2, :central)

@test dfdx1 ≈ reshape(dfdx2,100,10)

x=randn(100,10);
y=randn(100,10);
w=randn(100,10);
dfdx1=similar(x);
dfdx2=similar(x);

@time JuMIT.Misfits.error_cls!(dfdx1,x,y,w)
@time f(x)=JuMIT.Misfits.error_cls!(nothing,x,y,w)
@time ForwardDiff.gradient!(dfdx2,f, x);

@test dfdx1 ≈ reshape(dfdx2,100,10)


@time J=JuMIT.Misfits.error_corr_dist!(dfdx1,x,y)
@time f(x)=JuMIT.Misfits.error_corr_dist!(nothing,x,y)
@time ForwardDiff.gradient!(dfdx2, f, x);

@test dfdx1 ≈ reshape(dfdx2,100,10)


# some func
function f(x, z)
    y=x./vecnorm(x)
    J=sum((y-z).^2)
    return J
end

# test derivative_vector_magnitude
function g!(g, x, z)
    xn=vecnorm(x)
    scale!(x, inv(xn))
    g1=similar(g)
    for i in eachindex(g1)
        g1[i]=2.*(x[i]-z[i])
    end
         scale!(x, xn)
         nx=length(x)
         X=zeros(nx,nx)
         @time JuMIT.Misfits.derivative_vector_magnitude!(g,g1,x,X)
    return g
end

x=randn(10)
z=randn(10)
g1=zeros(x)
f1(x)=f(x,z)
ForwardDiff.gradient!(g1,f1, x)
g2=zeros(x)
@time g!(g2,x,z)
@test g1 ≈ g2
