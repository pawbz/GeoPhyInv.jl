using SIT
using Base.Test

# some model
model = SIT.Gallery.Seismic(:acou_homo1);
# allocation of model0
model0 = SIT.Models.Seismic_zeros(model.mgrid)
model0.vp0 = model.vp0; model0.ρ0 = model.ρ0; model0.mgrid = deepcopy(model.mgrid);

# reparameterize twice to see conversion formulae are correct
SIT.Models.Seismic_χKI_χρI!(model0, SIT.Models.Seismic_get(model,:χKI), SIT.Models.Seismic_get(model,:χρI))

@test_approx_eq model.χvp model0.χvp
@test_approx_eq model.χρ model0.χρ


nznx = model.mgrid.nz * model.mgrid.nx;
gχKI = randn(nznx); gχρI = randn(nznx)
g1 = similar(gχKI); g2 = similar(gχρI);

gmodel = SIT.Models.Seismic_zeros(model.mgrid)

# chain rule 1
SIT.Models.Seismic_gχKI_gχρI!(gmodel, model, gχKI, gχρI)

# chain rule 2
SIT.Models.Seismic_getgrad!(gmodel, model, [:χKI, :χρI], g1, g2)

@test_approx_eq gχKI g1
@test_approx_eq gχρI g2


