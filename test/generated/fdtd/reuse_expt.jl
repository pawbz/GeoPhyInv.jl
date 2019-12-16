using GeoPhyInv
using Statistics
using Test

zgrid=range(-1000.0,stop=1000.0,length=201)
xgrid=range(-1000.0,stop=1000.0,length=201)
mgrid = [zgrid, xgrid]
@info string("spatial sampling intervals (dz,dx)=", step.(mgrid))

vpb = [1500., 3500.] # bounds for vp
rhob = [1500., 3500.] # density bounds

model=Medium(mgrid)
update!(model, [:vp,:rho], [vpb,rhob])
fill!(model)
@info "a seismic model is created"

update!(model, [:vp,:rho], randn_perc=0.01); # add some random noise

ageom=AGeom(model.mgrid,:surf, SSrcs(3), Recs(30));

tgrid = range(0.0,stop=2.0,length=1000)

wav = ricker(10.0, tgrid, tpeak=0.25,);

srcwav = SrcWav(tgrid, ageom, [:P])
update!(srcwav, [:P], wav)

pa=SeisForwExpt(Fdtd(),npw=1,model=model,
	ageom=[ageom], srcwav=[srcwav],
	sflags=[2], rflags=[1],
	tgridmod=tgrid, verbose=true);

@time update!(pa);

model_new=Medium(:acou_homo1) # prepare another model
update!(model_new, [:vp,:rho], randn_perc=0.01)
update!(model_new, [:vp,:rho], constant_pert=0.03) # perturb the model

update!(pa, model_new)

@time update!(pa);

