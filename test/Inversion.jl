using Revise
using JuMIT
using Base.Test
using Optim


x=randn(10);
y1=randn(10);
y2=randn(10).*100.;

f1(x)=JuMIT.Misfits.error_squared_euclidean!(nothing, x, y1, nothing)
f2(x)=JuMIT.Misfits.error_squared_euclidean!(nothing, x, y2, nothing)

g1(st,x)=JuMIT.Misfits.error_squared_euclidean!(st, x, y1, nothing)
g2(st,x)=JuMIT.Misfits.error_squared_euclidean!(st, x, y2, nothing)

optim_func=[f1,f2]
optim_grad=[g1,g2]

pa=JuMIT.Inversion.ParamMO(noptim=2,x_init=randn(10),
    optim_func=optim_func,optim_grad=optim_grad)

pa.func(x, pa)

st=similar(x)
pa.grad!(st, x, pa)

using Optim
df = OnceDifferentiable(x -> pa.func(x, pa),
            (storage, x) -> pa.grad!(storage, x, pa))

# only first objective
pa.αvec=[1., 0.]
res=optimize(df, x, )
y11=Optim.minimizer(res)
@test y11 ≈  y1

# only second objective
pa.αvec=[0., 1.]
res=optimize(df, x, )
y22=Optim.minimizer(res)
@test y22 ≈ y2


# test finite difference gradient
for i in 1:3
    randn!(pa.αvec) # test for different weighting paramters
    g=similar(x)
    JuMIT.Inversion.finite_difference!(x -> pa.func(x,pa), x, g, :central)
    pa.grad!(st, x, pa)
    @test st ≈ g
end
