const USE_GPU = true  # Use GPU? If this is set false, then no GPU needs to be available
using ParallelStencil
using GeoPhyInv
using ParallelStencil.FiniteDifferences3D
using ProgressMeter
@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 3)
else
    @init_parallel_stencil(Threads, Float64, 3)
end
using Plots, Printf, Statistics

@parallel function compute_V!(Vx::Data.Array, Vy::Data.Array, Vz::Data.Array, tauxx::Data.Array, tauyy::Data.Array, tauzz::Data.Array, tauxy::Data.Array, tauxz::Data.Array, tauyz::Data.Array, dt::Data.Number, rho::Data.Number, dx::Data.Number, dy::Data.Number, dz::Data.Number)
    @inn(Vx) = @inn(Vx) - dt/rho*(@d_xi(tauxx)/dx+@d_yi(tauxy)/dy+@d_zi(tauxz)/dz)
    @inn(Vy) = @inn(Vy) - dt/rho*(@d_xi(tauxy)/dx+@d_yi(tauyy)/dy+@d_zi(tauyz)/dz)
    @inn(Vz) = @inn(Vz) - dt/rho*(@d_xi(tauxz)/dx+@d_yi(tauyz)/dy+@d_zi(tauzz)/dz)
    return
end

@parallel function compute_tauii!(tauxx::Data.Array, tauyy::Data.Array, tauzz::Data.Array,Vx::Data.Array, Vy::Data.Array, Vz::Data.Array, dt::Data.Number, lambda2mu::Data.Array, lambda::Data.Array, dx::Data.Number, dy::Data.Number, dz::Data.Number)
    @all(tauxx) = @all(tauxx) - dt*(@all(lambda2mu)*@d_xi(Vx)/dx + @all(lambda)*(@d_yi(Vy)/dy + @d_zi(Vz)/dz))
    @all(tauyy) = @all(tauyy) - dt*(@all(lambda2mu)*@d_yi(Vy)/dy + @all(lambda)*(@d_xi(Vx)/dx + @d_zi(Vz)/dz))
    @all(tauzz) = @all(tauzz) - dt*(@all(lambda2mu)*@d_zi(Vz)/dz + @all(lambda)*(@d_xi(Vx)/dx + @d_yi(Vy)/dy))
    return
end

@parallel function compute_tauij!(tauxy::Data.Array, tauxz::Data.Array, tauyz::Data.Array,Vx::Data.Array, Vy::Data.Array, Vz::Data.Array, dt::Data.Number, mu::Data.Array, dx::Data.Number, dy::Data.Number, dz::Data.Number)
    @all(tauxz) = @all(tauxz) - dt*(@all(mu)*(@d_zi(Vx)/dz + @d_xi(Vz)/dx))
    @all(tauxy) = @all(tauxy) - dt*(@all(mu)*(@d_yi(Vx)/dy + @d_xi(Vy)/dx))
    @all(tauyz) = @all(tauyz) - dt*(@all(mu)*(@d_zi(Vy)/dz + @d_yi(Vz)/dy))
    return
end

##################################################
# @views function acoustic3D()
    nx, ny, nz = 255, 255, 255     # numerical grid resolution; should be a mulitple of 32-1 for optimal GPU perf
    # Physics
    lx, ly, lz = 1000., 1000.0, 1000.0  # domain extends
    vp0=8000
    vs0=5000
    rho          = 1.0               # density
    lambda     = @zeros(nx,ny,nz)               # lame parameter
    mu     = @zeros(nx,ny,nz)  
    mu .=  vs0*vs0*rho            # mu parameter
    lambda .= vp0*vp0*rho .- 2 .* mu
    lambda2mu = lambda .+ 2 .* mu
    tmax          = 1
    nt         = 1000              # number of timesteps
    tgrid=range(0,tmax,length=nt)              # physical time
    wav=GeoPhyInv.ricker(5.,tgrid)
    # Numerics
    nout       = 100                # plotting frequency
    # Derived numerics
    dx, dy, dz = lx/(nx-1), ly/(ny-1), lz/(nz-1) # cell sizes
    # Array allocations
    tauxx          = @zeros(nx  ,ny  ,nz  )
    tauyy          = @zeros(nx  ,ny  ,nz  )
    tauzz          = @zeros(nx  ,ny  ,nz  )
    tauxy          = @zeros(nx  ,ny  ,nz  )
    tauxz          = @zeros(nx  ,ny  ,nz  )
    tauyz          = @zeros(nx  ,ny  ,nz  )

    Vx         = @zeros(nx+1,ny+1,nz+1)
    Vy         = @zeros(nx+1,ny+1,nz+1)
    Vz         = @zeros(nx+1,ny+1,nz+1)
    # Initial conditions
    # tauxx         .= Data.Array([exp(-((ix-1)*dx-0.5*lx)^2 -((iy-1)*dy-0.5*ly)^2 -((iz-1)*dz-0.5*lz)^2) for ix=1:size(tauxx,1), iy=1:size(tauxx,2), iz=1:size(tauxx,3)])
    # tauyy         .= Data.Array([exp(-((ix-1)*dx-0.5*lx)^2 -((iy-1)*dy-0.5*ly)^2 -((iz-1)*dz-0.5*lz)^2) for ix=1:size(tauxx,1), iy=1:size(tauxx,2), iz=1:size(tauxx,3)])
    # tauzz         .= Data.Array([exp(-((ix-1)*dx-0.5*lx)^2 -((iy-1)*dy-0.5*ly)^2 -((iz-1)*dz-0.5*lz)^2) for ix=1:size(tauxx,1), iy=1:size(tauxx,2), iz=1:size(tauxx,3)])
    dt         = min(dx,dy,dz)/sqrt(vp0*vp0+vs0*vs0)/10
    println(dt)
    # dt=1e-5
    # wefwrg
    # Preparation of visualisation
    ENV["GKSwstype"]="nul"; if isdir("viz3D_out")==false mkdir("viz3D_out") end; loadpath = "./viz3D_out/"; anim = Animation(loadpath,String[])
    println("Animation directory: $(anim.dir)")
    y_sl       = Int(ceil(ny/2))
    X, Y, Z    = -lx/2:dx:lx/2, -ly/2:dy:ly/2, -lz/2:dz:lz/2
    # Time loop
    @showprogress for it = 1:nt
        if (it==11)  global wtime0 = Base.time()  end
        @parallel compute_V!(Vx, Vy, Vz, tauxx,tauyy,tauzz,tauxy,tauxz,tauyz, dt, rho, dx, dy, dz)
        Vx[div(nx,2),div(ny,2),div(nz,2)]=wav[it]
        @parallel compute_tauii!(tauxx,tauyy,tauzz, Vx, Vy, Vz, dt, lambda2mu, lambda, dx, dy, dz)
        @parallel compute_tauij!(tauxy,tauxz,tauyz, Vx, Vy, Vz, dt, mu, dx, dy, dz)
        # t = t + dt
        # Visualisation
        if mod(it,nout)==0
            heatmap(X, Z, Array(tauxx)[:,y_sl,:]', aspect_ratio=1, xlims=(X[1],X[end]), ylims=(Z[1],Z[end]), c=:viridis, title="Pressure"); frame(anim)
        end
    end
    # Performance
    wtime    = Base.time() - wtime0
    A_eff    = (4*2)/1e9*nx*ny*nz*sizeof(Data.Number)  # Effective main memory access per iteration [GB] (Lower bound of required memory access: H and dHdτ have to be read and written (dHdτ for damping): 4 whole-array memaccess; B has to be read: 1 whole-array memaccess)
    wtime_it = wtime/(nt-10)                           # Execution time per iteration [s]
    T_eff    = A_eff/wtime_it                          # Effective memory throughput [GB/s]
    @printf("Total steps=%d, time=%1.3e sec (@ T_eff = %1.2f GB/s) \n", nt, wtime, round(T_eff, sigdigits=2))
    gif(anim, "acoustic3D.gif", fps = 15)
    # return
# end


# acoustic3D()
