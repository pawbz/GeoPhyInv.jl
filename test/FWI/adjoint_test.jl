using JuMIT
using Base.Test
using BenchmarkTools


const J=JuMIT
const JF=JuMIT.FWI

model = J.Gallery.Seismic(:acou_homo2);
J.Models.Seismic_addon!(model, ellip_rad=50., ellip_loc=[500.,0.],ellip_pert=0.1,randn_perc=0.0, fields=[:χvp,:χρ])

model0 = J.Gallery.Seismic(:acou_homo2);
acqgeom = J.Acquisition.Geom_circ(nss=1,nr=20,rad=[900.,900.])


acqsrc=J.Acquisition.Src_fixed_mod(acqgeom.nss,1,[:P],mod=model,nλ=3)
tgrid=acqsrc.tgrid

pa=JF.Param(acqsrc, acqgeom, tgrid, JF.ModFdtdBorn(), model0, 
		     modm_obs=model,  
		     igrid_interp_scheme=:B1,    
		     igrid=Grid.M2D_resamp(model.mgrid,150.,150.,),     
	             parameterization=[:χKI, :χρI, :null],   verbose=false)


H=JF.Fadj_Fborn(pa)

x=randn(size(H,2))
x[:]=0.0; x[100]=1.
y=randn(size(H,1))
y[:]=0.0; y[3000]=1.


y1=H*x
x1=transpose(H)*y1

a=dot(x1,x)
b=dot(y1,y)
println(a, " ", b)
println(a*inv(b), " ",  b*inv(a))

#=
#println("######")
#
#JF.F!(pa, nothing, pa.attrib_mod)

println(dot(pa.paf.c.δmodtt, pa.paf.c.δmodtt))
println(dot(pa.paf.c.data[2], pa.paf.c.data[2]))

gmod=pa.gmodm
g=zeros(2*gmod.mgrid.nx*gmod.mgrid.nz)


J.Models.Seismic_chainrule!(gmod,pa.modm0,g);

println(dot(x1,x1))
=#

