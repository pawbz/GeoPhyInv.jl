
using GeoPhyInv
using Test
using LinearAlgebra

# test for a arbit acoustic model
mgrid = [range(0.0, stop = 10.0, step = 0.05), range(0.0, stop = 10.0, step = 0.05)];

model = Medium(mgrid);
vpb = [2100.0, 2200.0];
rhob = [2100.0, 2300.0];
update!(model, [:vp, :rho], [vpb, rhob])
fill!(model)

nznx = prod(length.(model.mgrid))

for parameterization in [
    [:χKI, :χrhoI, :null],
    [:χKI, :null, :null],
    [:χvp, :χrho, :null],
    [:χvp, :null, :null],
]
    δx1 = randn(count(parameterization .≠ :null) * nznx)
    δxout1 = zero(δx1)

    copyto!(δxout1, δx1, model, parameterization)

    δxout2 = randn(size(δxout1))
    δx2 = zero(δx1)
    GeoPhyInv.pert_chainrule!(δx2, δxout2, model, parameterization)

    @test LinearAlgebra.dot(δxout2, δxout1) ≈ LinearAlgebra.dot(δx1, δx2)
end
