using JuMIT
reload("JuMIT")
using Base.Test

# some mem tests
mod=randn(1000); mod0=[2.0]
for i in 1:3
	@time JuMIT.Models.χ!(mod,mod0,1)
	@time JuMIT.Models.χ!(mod,mod0,-1)
	@time JuMIT.Models.χg!(mod,mod0,1)
	@time JuMIT.Models.χg!(mod,mod0,-1)
end
# some model
mgrid = JuMIT.Grid.M2D(0.0, 10., 0., 10., .05,.05,2);
model = JuMIT.Models.Seismic_zeros(mgrid);
model.vp0=[2100.,2200.]; model.vs0=[-1., -1.]; model.ρ0=[2100., 2300.]
JuMIT.Models.Seismic_addon!(model, randn_perc=1e-3)

# allocation of model0
model0 = JuMIT.Models.Seismic_zeros(model.mgrid);
model0.vp0 = model.vp0; model0.ρ0 = model.ρ0; model0.mgrid = deepcopy(model.mgrid);
model0.vs0 = [-1.0, -1.0];

# reparameterize twice to see conversion formulae are correct
JuMIT.Models.Seismic_reparameterize!(model0, JuMIT.Models.Seismic_get(model,:χKI),
				   JuMIT.Models.Seismic_get(model,:χρI),[:χKI,:χρI]);

@test_approx_eq model.χvp model0.χvp
@test_approx_eq model.χρ model0.χρ


nznx = model.mgrid.nz * model.mgrid.nx;
gχKI = randn(nznx); gχρI = randn(nznx);
g1 = similar(gχKI); g2 = similar(gχρI);

gmodel = JuMIT.Models.Seismic_zeros(model.mgrid);


for i in 1:2
	# chain rule 1
	@time JuMIT.Models.Seismic_chainrule!(gmodel, model, gχKI, gχρI, [:χKI, :χρI],1);

	# chain rule 2
	@time JuMIT.Models.Seismic_chainrule!(gmodel, model, g1, g2, [:χKI, :χρI], -1);
end


@test gχKI ≈ g1
@test gχρI ≈ g2


# testing pad_trun
a=randn(200,300)
np=50
b=zeros(200+2*np,300+2*np)
@time JuMIT.Models.pml_pad_trun!(b, a,1);
aa=similar(a);
@time JuMIT.Models.pml_pad_trun!(b, aa,-1);

@test aa==a
