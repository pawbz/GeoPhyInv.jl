
function SeisForwExpt(attrib::Symbol)
	@assert attrib in [:acou_homo2D]

	model=Seismic(:acou_homo2D);
	ageom=AGeom(model.mgrid,:xwell);
	tgrid=range(0.0,stop=2.0,length=1000)

	wav = ricker(10.0, tgrid, tpeak=0.25, );
	srcwav = SrcWav(tgrid, SSrcs(length(ageom)), [Srcs(1)], [:P]);
	update!(srcwav, [:P], wav)

	vp0=mean(GeoPhyInv.Models.χ(model.χvp,model.ref.vp,-1))
	ρ0=mean(GeoPhyInv.Models.χ(model.χρ,model.ref.ρ,-1))

	pa=SeisForwExpt(npw=1,model=model,
	    ageom=[ageom], srcwav=[srcwav],
		sflags=[2], rflags=[1],
		    tgridmod=tgrid, verbose=true);
	return pa


end

