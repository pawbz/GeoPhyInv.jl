using GeoPhyInv
using Statistics
#md using Plots; gr();

# ### Medium
medium=Medium(:acou_homo1); # load a simple homogeneous acoustic medium from the gallery
update!(medium, [:vp,:rho], randn_perc=5); # add some random noise to the medium

# ### AGeom
ageom=AGeom(medium.mgrid,:xwell, SSrcs(2)); # load a simple acquisition using `mgrid` of the medium

# ### SrcWav
tgrid = range(0.0,stop=2.0,length=2000); # generate a time grid
wav = ricker(10.0, tgrid, tpeak=0.25,); # ricker wavelet
srcwav = SrcWav(tgrid, ageom, [:p]);
update!(srcwav, [:p], wav);

# ### Plotting #1
#md p1=heatmap(medium, :vp)
#md scatter!(ageom, SSrcs())
#md scatter!(ageom, Recs())
#md plot(p1)


# ### SeisForwExpt
# Now we have all the required variables to create `SeisForwExpt` object and 
# prepare the forward modeling.
# While creating, we switched the `snaps_flag` on, and instruct recording field at
# `tsnaps`.
# Once the `Expt` object is created, do the modeling "without approximately any 
# memory allocations" using `update!`

pa=SeisForwExpt(Fdtd(),medium=medium, ageom=ageom, srcwav=srcwav,
	snaps_flag=true, tsnaps=[0.3, 0.4, 0.5], tgrid=tgrid, verbose=true);

# ### Modeling
@time update!(pa);

# ### Extracting snaps
snaps=pa[:snaps,1]; # extracting snaps of the first supersource
snaps=pa[:snaps,2]; # second supersource

# ### Plotting #2

#md p=[]
#md for ii in 1:3
#md 	push!(p,heatmap(medium.mgrid[2], medium.mgrid[1], snaps[:,:,ii], xlabel="x", ylabel="z"))
#md end
#md plot(p..., layout=(1,3), size=(1100,250))




