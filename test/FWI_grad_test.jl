using JuMIT
using Test
using Signals
using Misfits
using BenchmarkTools
using Grid
using Profile



const J=JuMIT
const JF=JuMIT.FWI

model = JuMIT.Gallery.Seismic(:acou_homo2);
model.mgrid.npml=5;
#JuMIT.Models.Seismic_addon!(model, ellip_rad=50., ellip_loc=[500.,0.],ellip_pert=0.1,randn_perc=0.01, fields=[:χvp,:χρ])

model0 = JuMIT.Gallery.Seismic(:acou_homo2);
model0.mgrid.npml=5;
acqgeom = JuMIT.Acquisition.Geom_circ(nss=1,nr=2,rad=[0.,200.])
JuMIT.Models.Seismic_addon!(model0, randn_perc=0.0001, fields=[:χvp,:χρ])


acqsrc=JuMIT.Acquisition.Src_fixed_mod(acqgeom.nss,1,[:P],mod=model,nλ=3, tmaxfrac=0.4)
tgrid=acqsrc.tgrid
println("Maximum modeling time", tgrid.x[end])

pa=JuMIT.FWI.Param(acqsrc, acqgeom, tgrid, JF.ModFdtd(), model0, 
		     modm_obs=model,  
		     igrid_interp_scheme=:B2,    
		     igrid=Grid.M2D_resamp(model.mgrid,300.,300.,),     
	             parameterization=[:χvp, :χρ, :null],   verbose=false)


JuMIT.FWI.xfwi!(pa, JuMIT.FWI.LS(),  bounded_flag=true, solver=:ipopt,
		ipopt_options=[["max_iter", 0],["derivative_test", "first-order"]])

pa_fd=deepcopy(pa);

migr=JuMIT.FWI.xfwi!(pa, JuMIT.FWI.Migr())

llll
Profile.clear_malloc_data()

migr_fd=JuMIT.FWI.xfwi!(pa_fd, JuMIT.FWI.Migr_fd())


f=Misfits.error_squared_euclidean!(nothing, migr.χvp, migr_fd.χvp, nothing, norm_flag=true)

@test f<1e-15

f=Misfits.error_squared_euclidean!(nothing, migr.χρ, migr_fd.χρ, nothing, norm_flag=true)

@test f<1e-15

