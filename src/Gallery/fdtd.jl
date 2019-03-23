
"""
Gallery of `SeisForwExpt`.
# Arguments
* `attrib::Symbol` : name 
  * `=:acou_homo1` 
"""
function SeisForwExpt(attrib::Symbol)
	@assert attrib in [:acou_homo1]

	model=Seismic(:acou_homo1);
	acqgeom =Geom(model.mgrid,:xwell);
	tgrid=range(0.0,stop=2.0,length=1000)
	wav=GeoPhyInv.Utils.Wavelets.ricker(10.0, tgrid, tpeak=0.25, );
	# source wavelet for modelling
	acqsrc=GeoPhyInv.Acquisition.Src_fixed(acqgeom.nss,1,[:P],wav,tgrid);

	vp0=mean(GeoPhyInv.Models.χ(model.χvp,model.ref.vp,-1))
	ρ0=mean(GeoPhyInv.Models.χ(model.χρ,model.ref.ρ,-1))

	pa=GeoPhyInv.Fdtd.Param(npw=1,model=model,
	    acqgeom=[acqgeom], acqsrc=[acqsrc],
		sflags=[2], rflags=[1],
		    tgridmod=tgrid, verbose=true);
	return pa


end

