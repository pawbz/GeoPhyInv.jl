# A forward experiment, where the seismic data are generated
# using some models and acquisition parameters from our gallery.
# Forward modeling consists of a finite-difference simulation of the acoustic wave-equation.
# We specifically aim to save the snapshots, at given time steps in `SeisForwExpt`.

# ### Loading some packages...

using GeoPhyInv
using Statistics

# Then, load a simple homogeneous acoustic model from the gallery.

model=GIPh.Gallery.Seismic(:acou_homo1);

# Add some random noise to the model.

GIPh.Models.Seismic_addon!(model, randn_perc=0.01);

# Load a simple acquisition geometry from the gallery using the field `mgrid` of the seismic model.

acqgeom=GIPh.Gallery.Geom(model.mgrid,:xwell);

# Plot the model, source and receivers using these commands:
# `using Plots`
# `p1=JP.seismic(model);`
# `JP.geom!(acqgeom);`
# `plot(p1);`

# Generate a time grid.

tgrid = range(0.0,stop=2.0,length=1000)

# Ricker wavelet.

wav = GIPh.Utils.Wavelets.ricker(10.0, tgrid, tpeak=0.25,);

# Distribute the same source wavelet to all the supsersources.

acqsrc=GIPh.Acquisition.Src_fixed(acqgeom.nss,1,[:P],wav,tgrid);

# Now we have all the required variables to create `SeisForwExpt` object and 
# prepare the forward modelling.
# While creating, we switched the `snaps_flag` on, and instructed recording field at
# `tsnaps`.
# Once the `Expt` object is created, do the modelling "without approximately any 
# memory allocations" using `mod!`

paE=SeisForwExpt(model=model,
	acqgeom=[acqgeom], acqsrc=[acqsrc],
	snaps_flag=true,
	tsnaps=[0.3, 0.4, 0.5],
	tgridmod=tgrid, verbose=true);

@time mod!(paE);

# Extracting snaps of the first supersource. The dimensions of the snaps is `[nz,nx,nt]`.
snaps=paE[:snaps,1];

# Extracting snaps of the second supersource.
snaps=paE[:snaps,2];

# We can plot snapshots using these commands:
# `p1=[heatmap(snaps[:,:,ii]) for ii in 1:3];`
# `plot(p1..., layout=(1,3), aspect_ratio=:equal)`



