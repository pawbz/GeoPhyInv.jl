
# some model
mgrid = Grid.M2D(0.0, 10., 0., 10., .05,.05,2);
model = JuMIT.Models.Seismic_zeros(mgrid);
vp0=[2100.,2200.];vs0=[-1., -1.]; ρ0=[2100., 2300.]
JuMIT.Models.adjust_bounds!(model, vp0,vs0,ρ0);


for parameterization in [[:χKI, :χρI, :null],[:χKI, :null, :null],[:χvp, :χρ, :null],[:χvp, :null, :null]]
	δx1=randn(count(parameterization .≠ :null)*mgrid.nz*mgrid.nx)
	δxout1=zero(δx1)

	J.Models.pert_reparameterize!(δxout1, δx1, model, parameterization)


	δxout2=randn(size(δxout1));
	δx2=zero(δx1)
	J.Models.pert_gradient_chainrule!(δx2, δxout2, model, parameterization)

	@test LinearAlgebra.dot(δxout2,δxout1) ≈ LinearAlgebra.dot(δx1,δx2)
end
