# This tutorial first creates a variable `pa`, i.e. allocates
# necessary memory to perform `SeisForwExpt`. 
# Then, we perform forward modeling using an in-place 
# function `mod!`. Finally, we will update `pa` with a different 
# subsurface model and re-run the modeling task with no additional memory 
# allocation. 
# The ability to iteratively run the forward modeling task on  
# various subsurface models is necessary while implementing inversion 
# algorithms.

# ### Load packages

using GeoPhyInv
using Statistics
using Plots


# ### Setup a spatial grid
zgrid=range(-1000.0,stop=1000.0,length=201)
xgrid=range(-1000.0,stop=1000.0,length=201)
mgrid = [zgrid, xgrid]
@info string("spatial sampling intervals (dz,dx)=", step.(mgrid))


# ### Allocate a `Seismic` model, and adjust bounds
vpb = [1500., 3500.] # bounds for vp
rhob = [1500., 3500.] # density bounds

model=Medium(mgrid)
update!(model, [:vp,:rho], [vpb,rhob])
fill!(model)
@info "a seismic model is created" 

# ### Add some noise to the model (optional) 
update!(model, [:vp,:rho], randn_perc=0.01); # add some random noise

# ### A surface acquisition ageometry
ageom=GeoPhyInv.Gallery.AGeom(model.mgrid,:surf, nss=3, nr=30);

# ### Plot the model and source, receivers
p1=JP.seismic(model) 
JP.ageom!(ageom)
plot(p1)

# ### Generate a temporal grid
tgrid = range(0.0,stop=2.0,length=1000)

# ### Choose a source wavelet
wav = GeoPhyInv.Utils.Wavelets.ricker(10.0, tgrid, tpeak=0.25,);

# ### Distribute the same source wavelet to all the supsersources 
srcwav=GeoPhyInv.Acquisition.Src_fixed(ageom.nss,1,[:P],wav,tgrid);

# create `Fdtd.Param` object to prepare forward modelling
# * npw corresponds to the number of independently propagating wavefields (1 in most cases)
# Once the `Param` object is created, do the modelling "without any memory allocations" using `mod!`

pa=SeisForwExpt(npw=1,model=model,
	ageom=[ageom], srcwav=[srcwav],
	sflags=[2], rflags=[1],
	tgridmod=tgrid, verbose=true);

@time GeoPhyInv.Fdtd.mod!(pa);

# plot a record after modelling
pdata=plot(pa.c.data[1].d[1,1])
plot(pdata)

# create new seismic model

model_new=J.Gallery.Seismic(:acou_homo1) # prepare another model
update!(model_new, [:vp,:rho], randn_perc=0.01)
update!(model_new, [:vp,:rho], constant_pert=0.03) # perturb the model
p2=JP.seismic(model_new) # plot new model 
JP.ageom!(ageom)
plot(p2)


# Now, we the change the model in the `Param` object without memory allocation
# This routine can be used during FWI, 
# where medium parameters are itertively updated in the same `Fdtd.Param` object

update!(pa, model_new)




# run modelling now and plot data again
@time GeoPhyInv.Fdtd.mod!(pa);

# plot a record after modelling
plot!(pdata, pa.c.data[1].d[1,1])
plot(pdata)

