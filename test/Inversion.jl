using Revise
using JuMIT
using Base.Test
using Optim
using BenchmarkTools


x=randn(100)
y1=randn(100);
y2=randn(100).*100.;

f1(x)=JuMIT.Misfits.error_squared_euclidean!(nothing, x, y1, nothing)
f2(x)=JuMIT.Misfits.error_squared_euclidean!(nothing, x, y2, nothing)

g1(st,x)=JuMIT.Misfits.error_squared_euclidean!(st, x, y1, nothing)
g2(st,x)=JuMIT.Misfits.error_squared_euclidean!(st, x, y2, nothing)

optim_func=[f1,f2]
optim_grad=[g1,g2]

pa=JuMIT.Inversion.ParamMO(noptim=2,x_init=randn(size(x)),
    optim_func=optim_func,optim_grad=optim_grad)

@btime pa.func(x, pa)

st=similar(x)
@btime pa.grad!(st, x, pa)

f=x -> pa.func(x, pa)
g=(storage, x) -> pa.grad!(storage, x, pa)

# only first objective
pa.αvec=[1., 0.]
res=optimize(f,g, x, )
y11=Optim.minimizer(res)
@test y11 ≈  y1

# only second objective
pa.αvec=[0., 1.]
res=optimize(f,g, x, )
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
