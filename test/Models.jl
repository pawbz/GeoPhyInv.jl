
# some mem tests
mod=randn(1000); mod0=2.0
@btime JuMIT.Models.χ!(mod,mod0,1)
@btime JuMIT.Models.χ!(mod,mod0,-1)
@btime JuMIT.Models.χg!(mod,mod0,1)
@btime JuMIT.Models.χg!(mod,mod0,-1)

# some model
mgrid = Grid.M2D(0.0, 10., 0., 10., .05,.05,2);
model = JuMIT.Models.Seismic_zeros(mgrid);
vp0=[2100.,2200.];vs0=[-1., -1.]; ρ0=[2100., 2300.]
JuMIT.Models.adjust_bounds!(model, vp0,vs0,ρ0);

JuMIT.Models.Seismic_addon!(model, randn_perc=1e-3)
nznx = model.mgrid.nz * model.mgrid.nx;

# test copy!
model_new=deepcopy(model)
JuMIT.Models.Seismic_addon!(model_new, randn_perc=1e-3)

@btime copy!(model_new, model)
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

nznx = model.mgrid.nz * model.mgrid.nx;
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




# testing pad_trun
a=randn(200,300)
np=50
b=zeros(200+2*np,300+2*np)
@time JuMIT.Models.pml_pad_trun!(b, a,1);
aa=similar(a);
@time JuMIT.Models.pml_pad_trun!(b, aa,-1);

@test aa==a
