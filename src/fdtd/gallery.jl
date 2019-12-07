
"""
Gallery of `SeisForwExpt`.
# Arguments
* `attrib::Symbol` : name 
  * `=:acou_homo1` 
"""
function SeisForwExpt(attrib::Symbol)
	@assert attrib in [:acou_homo1]

	model=Seismic(:acou_homo1);
	acqgeom=Geom(model.mgrid,:xwell);
	tgrid=range(0.0,stop=2.0,length=1000)

	wav = ricker(10.0, tgrid, tpeak=0.25, );
	acqsrc = SrcWav(tgrid, SSrcs(length(acqgeom)), [Srcs(1)], [:P]);
	update!(acqsrc, [:P], wav)

	vp0=mean(GeoPhyInv.Models.χ(model.χvp,model.ref.vp,-1))
	ρ0=mean(GeoPhyInv.Models.χ(model.χρ,model.ref.ρ,-1))

	pa=SeisForwExpt(npw=1,model=model,
	    acqgeom=[acqgeom], acqsrc=[acqsrc],
		sflags=[2], rflags=[1],
		    tgridmod=tgrid, verbose=true);
	return pa


end

