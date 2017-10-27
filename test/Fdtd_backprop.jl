addprocs(4)
using JuMIT
using Base.Test
reload("JuMIT")

model=JuMIT.Gallery.Seismic(:acou_homo1)
acqgeom =JuMIT.Acquisition.Geom_circ(nss=4,nr=100,rad=[700.,700.]);
acqsrc=JuMIT.Acquisition.Src_fixed_mod(acqgeom.nss,1,[:P],mod=model, nλ=3, tmaxfrac=0.7)
@time pa=JuMIT.Fdtd.Param(born_flag=false,npw=1, tgridmod=acqsrc.tgrid,
        gmodel_flag=true,
		sflags=[2],
        snaps_flag=true,
        backprop_flag=1,
        illum_flag=true,acqgeom=[acqgeom], acqsrc=[acqsrc],
        model=model);

JuMIT.Fdtd.mod!(pa);
rec1=deepcopy(pa.c.data[1])

pa.c.sflags=[3];
pa.c.backprop_flag=-1
JuMIT.Fdtd.mod!(pa);
rec2=deepcopy(pa.c.data[1])

# compare results

# time reverse
JuMIT.Data.TD_tr!(rec2);
# least-squares misfit
err = JuMIT.Misfits.TD!(nothing, rec1, rec2)

# normalization
error = err/dot(rec1, rec1)

println("error:\t", error)

# desired accuracy?
@test error<1e-9
