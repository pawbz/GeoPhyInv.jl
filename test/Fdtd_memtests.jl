#addprocs(4)
using JuMIT
reload("JuMIT")
using DistributedArrays

model=JuMIT.Gallery.Seismic(:acou_homo1)
acqgeom =JuMIT.Acquisition.Geom_circ(nss=1,nr=100,rad=[700.,700.]);
acqsrc=JuMIT.Acquisition.Src_fixed_mod(acqgeom.nss,1,[:P],mod=model, nλ=3, tmaxfrac=0.7)
@time pa=JuMIT.Fdtd.Param(born_flag=true,npw=2, tgridmod=acqsrc.tgrid,
        gmodel_flag=true,
        snaps_flag=true,
        backprop_flag=1,
        illum_flag=true,acqgeom=[acqgeom,acqgeom], acqsrc=[acqsrc,acqsrc],
        model=model);
pac=pa.c;
pass=JuMIT.Fdtd.Paramss(1, pac)

# check memory of all these methods
for i in 1:3
        println("=============================")
        print("add source\t")
        @time JuMIT.Fdtd.add_source!(1,1,1,pac,[pass],localpart(pa.p))
        print("record\t")
        @time JuMIT.Fdtd.record!(1,1,1,pac,[pass],localpart(pa.p))
        print("compute gradient\t")
        @time JuMIT.Fdtd.compute_gradient!(1,pac,[pass],localpart(pa.p))
        print("advance\t")
        @time JuMIT.Fdtd.advance!(pac,localpart(pa.p))
        print("compute illum\t")
        @time JuMIT.Fdtd.compute_illum!(1, [pass], localpart(pa.p))
        print("add born sources\t")
        @time JuMIT.Fdtd.add_born_sources!(1,pac,[pass],localpart(pa.p))
        print("boundary force\t")
        @time JuMIT.Fdtd.boundary_force!(1,1,pac,[pass],localpart(pa.p))
        print("boundary save\t")
        @time JuMIT.Fdtd.boundary_save!(1,1,pac,[pass],localpart(pa.p))
        print("boundary save snap p\t")
        @time JuMIT.Fdtd.boundary_save_snap_p!(1,pac,[pass],localpart(pa.p))
        print("boundary save snap vzvx\t")
        @time JuMIT.Fdtd.boundary_save_snap_vxvz!(1,pac,[pass],localpart(pa.p))
        print("boundary force snap p\t")
        @time JuMIT.Fdtd.boundary_force_snap_p!(1,pac,[pass],localpart(pa.p))
        print("boundary force snap vxvz\t")
        @time JuMIT.Fdtd.boundary_force_snap_vxvz!(1,pac,[pass],localpart(pa.p))
        print("mod per proc\t")
        @time JuMIT.Fdtd.mod_per_proc!(pac,localpart(pa.p))
        print("stack grads\t")
        @time JuMIT.Fdtd.stack_grads!(pac, localpart(pa.p))
        print("stack illums\t")
        @time JuMIT.Fdtd.stack_illums!(pac, localpart(pa.p))
        print("update gmodel\t")
        @time JuMIT.Fdtd.update_gmodel!(pac)
        print("update data\t")
        @time JuMIT.Fdtd.update_datamat!(pac,localpart(pa.p))
        print("mod\t")
        @time JuMIT.Fdtd.mod!(pa);
        println("=============================")
end


#for nloop in [4,8,16]
#        print("looping mod! for:\t",nloop,"\t")
#        @time for iloop in 1:nloop
#                JuMIT.Fdtd.mod!(pa)
#        end
#end
