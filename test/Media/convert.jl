

using BenchmarkTools
using GeoPhyInv
using Test

# some mem tests
mod=randn(1000); mod0=2.0 #src
@btime GIPh.χ!(mod,mod0,1) #src
@btime GIPh.χ!(mod,mod0,-1) #src
@btime GIPh.χg!(mod,mod0,1) #src
@btime GIPh.χg!(mod,mod0,-1) #src



# some model
mgrid = [range(0.0, stop=10.,step=0.05), range(0.0, stop=10.,step=0.05)];
model = Medium(mgrid);


vpb=[2100.,2200.];vsb=[1500, 1700]; rhob=[2100., 2300.]
update!(model, [:vp,:vs,:rho], [vpb, vsb, rhob]);
fill!(model)


update!(model, [:vp,:rho], randn_perc=1e-3)

nznx=prod(length.(model.mgrid))

# test copyto!
model0=similar(model);
update!(model0, model)
fill!(model0)
update!(model0, [:vp,:rho], randn_perc=1e-3)


@btime copyto!(model0, model)
@test isequal(model0, model) 

# allocation of model0

# reparameterize twice to see conversion formulae are correct
for fields in [[:χKI,:χrhoI,:null], [:χKI,:null,:null], [:null,:χrho,:null], [:null,:χrho,:null], [:χvp,:χrho,:null]]
	x=zeros(nznx*count(fields.≠:null));
        copyto!(x, model, fields);
	@time copyto!(model0, x, fields);
	@test model[:vp] ≈ model0[:vp]
	@test model[:rho] ≈ model0[:rho]
end


nznx=prod(length.(model.mgrid))
G1 = randn(nznx); G2 = randn(nznx);
g1 = similar(G1); g2 = similar(G2);

gmodel = Medium(model.mgrid, [:χvp,:χrho]);


for fields in [[:χKI,:χrhoI,:null], [:χKI,:null,:null], [:χvp,:χrho,:null] , [:null,:χrho,:null], [:χvp,:null,:null]]
	G=randn(nznx*count(fields.≠:null));
	g=zeros(nznx*count(fields.≠:null));
	@time GIPh.chainrule!(gmodel, model, G, fields,1);
	@time GIPh.chainrule!(gmodel, model, g, fields, -1);
	@test G ≈ g
end



