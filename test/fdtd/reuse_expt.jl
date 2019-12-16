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
#md using Plots
using Test


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
ageom=AGeom(model.mgrid,:surf, SSrcs(3), Recs(30));

# ### Plot the model and source, receivers
#md p1=heatmap(model, :vp) 
#md scatter!(ageom, SSrcs())
#md scatter!(ageom, Recs())
#md plot(p1)

# ### Generate a temporal grid
tgrid = range(0.0,stop=2.0,length=1000)

# ### Choose a source wavelet
wav = ricker(10.0, tgrid, tpeak=0.25,);

# ### Distribute the same source wavelet to all the supsersources 
srcwav = SrcWav(tgrid, ageom, [:P])
update!(srcwav, [:P], wav)

# create `Fdtd.Param` object to prepare forward modelling
# * npw corresponds to the number of independently propagating wavefields (1 in most cases)
# Once the `Param` object is created, do the modelling "without any memory allocations" using `mod!`

pa=SeisForwExpt(Fdtd(),npw=1,model=model,
	ageom=[ageom], srcwav=[srcwav],
	sflags=[2], rflags=[1],
	tgridmod=tgrid, verbose=true);

@time update!(pa);

# plot a record after modelling
#md pdata=heatmap(pa[:data])
#md plot(pdata)

# create new seismic model

model_new=Medium(:acou_homo1) # prepare another model
update!(model_new, [:vp,:rho], randn_perc=0.01)
update!(model_new, [:vp,:rho], constant_pert=0.03) # perturb the model
#md p2=heatmap(model_new, :vp) # plot new model 
#md scatter!(ageom, SSrcs())
#md scatter!(ageom, Recs())
#md plot(p2)


# Now, we the change the model in the `Param` object without memory allocation
# This routine can be used during FWI, 
# where medium parameters are itertively updated in the same `Fdtd.Param` object

update!(pa, model_new)




# run modelling now and plot data again
@time update!(pa);

# plot a record after modelling
#md heatmap!(pdata, pa[:data])
#md plot(pdata)

