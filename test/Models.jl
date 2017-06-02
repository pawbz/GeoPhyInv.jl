using SIT
using Base.Test

# some model
mgrid = SIT.Grid.M2D(0.0, 10., 0., 10., 5,5,2);
model = SIT.Models.Seismic_zeros(mgrid);
model.vp0=2000.; model.vs0=-1.; model.ρ0=2000.
SIT.Models.Seismic_addon!(model, randn_perc=1e-3)

# allocation of model0
model0 = SIT.Models.Seismic_zeros(model.mgrid);
model0.vp0 = model.vp0; model0.ρ0 = model.ρ0; model0.mgrid = deepcopy(model.mgrid);
model0.vs0 = -1.0;

# reparameterize twice to see conversion formulae are correct
SIT.Models.Seismic_reparameterize!(model0, SIT.Models.Seismic_get(model,:χKI), 
				   SIT.Models.Seismic_get(model,:χρI),[:χKI,:χρI]);

@test_approx_eq model.χvp model0.χvp
@test_approx_eq model.χρ model0.χρ


nznx = model.mgrid.nz * model.mgrid.nx;
gχKI = randn(nznx); gχρI = randn(nznx);
g1 = similar(gχKI); g2 = similar(gχρI);

gmodel = SIT.Models.Seismic_zeros(model.mgrid);

# chain rule 1
SIT.Models.Seismic_chainrule!(gmodel, model, gχKI, gχρI, [:χKI, :χρI],1);

# chain rule 2
SIT.Models.Seismic_chainrule!(gmodel, model, g1, g2, [:χKI, :χρI], -1);

@test_approx_eq gχKI g1
@test_approx_eq gχρI g2


