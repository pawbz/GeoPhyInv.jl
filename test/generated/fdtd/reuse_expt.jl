using GeoPhyInv
using Statistics
using Plots
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

ageom=GeoPhyInv.Gallery.AGeom(model.mgrid,:surf, nss=3, nr=30);

p1=JP.seismic(model)
JP.ageom!(ageom)
plot(p1)

tgrid = range(0.0,stop=2.0,length=1000)

wav = GeoPhyInv.Utils.Wavelets.ricker(10.0, tgrid, tpeak=0.25,);

srcwav=GeoPhyInv.Acquisition.Src_fixed(ageom.nss,1,[:P],wav,tgrid);

pa=SeisForwExpt(npw=1,model=model,
	ageom=[ageom], srcwav=[srcwav],
	sflags=[2], rflags=[1],
	tgridmod=tgrid, verbose=true);

@time GeoPhyInv.Fdtd.mod!(pa);

pdata=plot(pa.c.data[1].d[1,1])
plot(pdata)

model_new=J.Gallery.Seismic(:acou_homo1) # prepare another model
update!(model_new, [:vp,:rho], randn_perc=0.01)
update!(model_new, [:vp,:rho], constant_pert=0.03) # perturb the model
p2=JP.seismic(model_new) # plot new model
JP.ageom!(ageom)
plot(p2)

update!(pa, model_new)

@time GeoPhyInv.Fdtd.mod!(pa);

plot!(pdata, pa.c.data[1].d[1,1])
plot(pdata)

