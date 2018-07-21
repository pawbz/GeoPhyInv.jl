using JuMIT
using Base.Test
using Signals
using Misfits
using BenchmarkTools



model = JuMIT.Gallery.Seismic(:acou_homo2);
JuMIT.Models.Seismic_addon!(model, ellip_rad=50., ellip_loc=[500.,0.],ellip_pert=0.1,randn_perc=0.0, fields=[:χvp,:χρ])

model0 = JuMIT.Gallery.Seismic(:acou_homo2);
acqgeom = JuMIT.Acquisition.Geom_circ(nss=1,nr=20,rad=[900.,900.])


acqsrc=JuMIT.Acquisition.Src_fixed_mod(acqgeom.nss,1,[:P],mod=model,nλ=3)
tgrid=acqsrc.tgrid

pa=JuMIT.FWI.Param(acqsrc, acqgeom, tgrid, :fdtd, model0, 
		     modm_obs=model,  
		     igrid_interp_scheme=:B2,    
		     igrid=Grid.M2D_resamp(model.mgrid,150.,150.,),     
	             parameterization=[:χvp, :χρ, :null],   verbose=false)

pa_fd=deepcopy(pa);

migr=JuMIT.FWI.xfwi!(pa, JuMIT.FWI.Migr())

Profile.clear_malloc_data()

migr_fd=JuMIT.FWI.xfwi!(pa_fd, JuMIT.FWI.Migr_fd())


f=Misfits.error_squared_euclidean!(nothing, migr.χvp, migr_fd.χvp, nothing, norm_flag=true)

@test f<1e-15

f=Misfits.error_squared_euclidean!(nothing, migr.χρ, migr_fd.χρ, nothing, norm_flag=true)

@test f<1e-15

