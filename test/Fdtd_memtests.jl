using JuMIT
reload("JuMIT")
using DistributedArrays

@time pa=JuMIT.Fdtd.Param(born_flag=true,npropwav=2);
pac=pa.c;
pass=JuMIT.Fdtd.Paramss(1, pac)

# check memory of all these methods
for i in 1:3
        println("=============================")
        @time JuMIT.Fdtd.add_source!(1,1,pac,[pass],localpart(pa.p))
        @time JuMIT.Fdtd.record!(1,1,pac,[pass],localpart(pa.p))
        @time JuMIT.Fdtd.compute_gradient!(1,pac,[pass],localpart(pa.p))
        @time JuMIT.Fdtd.advance!(pac,localpart(pa.p))
        @time JuMIT.Fdtd.compute_illum!(1, [pass], localpart(pa.p))
        @time JuMIT.Fdtd.add_born_sources!(1,pac,[pass],localpart(pa.p))
        @time JuMIT.Fdtd.boundary_force!(1,1,pac,[pass],localpart(pa.p))
        @time JuMIT.Fdtd.boundary_save!(1,1,pac,[pass],localpart(pa.p))
        @time JuMIT.Fdtd.boundary_save_snap_p!(1,pac,[pass],localpart(pa.p))
        @time JuMIT.Fdtd.boundary_save_snap_vxvz!(1,pac,[pass],localpart(pa.p))
        @time JuMIT.Fdtd.boundary_force_snap_p!(1,pac,[pass],localpart(pa.p))
        @time JuMIT.Fdtd.boundary_force_snap_vxvz!(1,pac,[pass],localpart(pa.p))
        @time JuMIT.Fdtd.mod_per_proc!(pac,localpart(pa.p))
        @time JuMIT.Fdtd.stack_grads!(pac, localpart(pa.p))
        @time JuMIT.Fdtd.stack_illums!(pac, localpart(pa.p))
        @time JuMIT.Fdtd.update_gmodel!(pac)
        @time JuMIT.Fdtd.mod!(pa);
        println("=============================")
end


addprocs(3)

@everywhere using DistributedArrays

a=dones(1000)

procs(a)

aa=ones(3)
@time b=sum(a)
@time sum(aa)
b=0.;
@sync begin
        for (ip,p) in enumerate(procs(a))
                @sync remotecall_wait(p) do
                        addd!(b, localpart(a))
                end
        end
end

for aaa in a
        b += aaa
end

@everywhere function addd!(a,b)
        a += sum(b)
end
