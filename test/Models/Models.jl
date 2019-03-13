
# some mem tests
mod=randn(1000); mod0=2.0
@btime JuMIT.Models.χ!(mod,mod0,1)
@btime JuMIT.Models.χ!(mod,mod0,-1)
@btime JuMIT.Models.χg!(mod,mod0,1)
@btime JuMIT.Models.χg!(mod,mod0,-1)

# some model
mgrid = [range(0.0, stop=10.,step=0.05), range(0.0, stop=10.,step=0.05)];
model = JuMIT.Models.Seismic_zeros(mgrid);
vp0=[2100.,2200.];vs0=[-1., -1.]; ρ0=[2100., 2300.]
JuMIT.Models.adjust_bounds!(model, vp0,vs0,ρ0);

JuMIT.Models.Seismic_addon!(model, randn_perc=1e-3)
nznx=prod(length.(model.mgrid))

# test copyto!
model_new=deepcopy(model)
JuMIT.Models.Seismic_addon!(model_new, randn_perc=1e-3)

@btime copyto!(model_new, model)
@test isequal(model_new, model) 

# allocation of model0
model0 = JuMIT.Models.Seismic_zeros(model.mgrid);
JuMIT.Models.adjust_bounds!(model0, model);

for attribvec in [[:χKI,:χρI,:null], [:χKI,:null,:null], [:null,:χρ,:null], [:null,:χρ,:null], [:χvp,:χρ,:null]]
	x=zeros(nznx*count(attribvec.≠:null));
	# reparameterize twice to see conversion formulae are correct
        JuMIT.Models.Seismic_get!(x, model, attribvec);
	@time JuMIT.Models.Seismic_reparameterize!(model0, x, attribvec);
	@test model.χvp ≈ model0.χvp
	#@test model.χvs ≈ model0.χvs
	@test model.χρ ≈ model0.χρ
end

nznx=prod(length.(model.mgrid))
G1 = randn(nznx); G2 = randn(nznx);
g1 = similar(G1); g2 = similar(G2);

gmodel = JuMIT.Models.Seismic_zeros(model.mgrid);


for attribvec in [[:χKI,:χρI,:null], [:χKI,:null,:null], [:χvp,:χρ,:null] , [:null,:χρ,:null], [:χvp,:null,:null]]
	G=randn(nznx*count(attribvec.≠:null));
	g=zeros(nznx*count(attribvec.≠:null));
	# chain rule 1
	@time JuMIT.Models.Seismic_chainrule!(gmodel, model, G, attribvec,1);

	# chain rule 2
	@time JuMIT.Models.Seismic_chainrule!(gmodel, model, g, attribvec, -1);
	@test G ≈ g
end



