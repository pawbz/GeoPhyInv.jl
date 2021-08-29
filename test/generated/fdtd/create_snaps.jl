using GeoPhyInv
using Statistics

medium=Medium(:acou_homo1); # load a simple homogeneous acoustic medium from the gallery
update!(medium, [:vp,:rho], randn_perc=5); # add some random noise to the medium

ageom=AGeom(medium.mgrid,:xwell, SSrcs(2)); # load a simple acquisition using `mgrid` of the medium

tgrid = range(0.0,stop=2.0,length=2000); # generate a time grid
wav = ricker(10.0, tgrid, tpeak=0.25,); # ricker wavelet
srcwav = SrcWav(tgrid, ageom, [:P]);
update!(srcwav, [:P], wav);

pa=SeisForwExpt(FdtdAcou(),medium=medium, ageom=ageom, srcwav=srcwav,
	snaps_flag=true, tsnaps=[0.3, 0.4, 0.5], tgrid=tgrid, verbose=true);

@time update!(pa);

snaps=pa[:snaps,1]; # extracting snaps of the first supersource
snaps=pa[:snaps,2]; # second supersource

