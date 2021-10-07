
global const marmousi_folder=joinpath(@__DIR__, "marmousi2")
global const overthrust_folder=joinpath(@__DIR__, "overthrust")

"""
```julia
mod=Medium(attrib)
```
As of now, only seismic models are predefined in this package. Choose `attrib::Symbol`

* `=:acou_homo1` : a test homogeneous acoustic model
* `=:acou_homo2` : a test homogeneous acoustic model, but with coarser spatial sampling (faster testing)
* `=:marmousi2` : marmousi2 model with lower resolution; useful for surface seismic experiments
* `=:marmousi2_small` : a smaller section of marmousi2 
* `=:pizza` : `:acou_homo2` with some perturbations 
"""
function Medium(attrib::Symbol, δ::Real = 0.0; verbose = false)
    bfrac = 0.1
    δ = Float64(δ)
    if ((attrib == :acou_homo1))
        vp0 = [1500.0, 3500.0] # bounds for vp
        rho0 = [1500.0, 3500.0] # density bounds
        mgrid = repeat([range(-1000.0, stop = 1000.0, length = 201)], 2)
        nz, nx = length.(mgrid)
        model = Medium(mgrid, [:vp, :rho])
        update!(model, [:vp, :rho], [vp0, rho0])
        fill!(model)
    elseif ((attrib == :elastic_homo1))
        vp0 = [3000.0, 3500.0] # bounds for vp
        vs0 = [1900, 2000] # bounds for vs
        rho0 = [1000.0, 2000.0] # density bounds
        mgrid = fill(range(-500, stop = 500, length = 255), 3)
        vp0 = [1500.0, 3500.0] # bounds for vp
        rho0 = [1500.0, 3500.0] # density bounds
        model = Medium(mgrid, [:vp, :vs, :rho])
        update!(model, [:vp, :vs, :rho], [vp0, vs0, rho0])
        fill!(model)
    elseif ((attrib == :pizza))
        vp0 = [1700.0, 2300.0] # bounds for vp
        rho0 = [1700.0, 2300.0] # density bounds
        mgrid = repeat([range(-1000.0, stop = 1000.0, length = 51)], 2)
        nz, nx = length.(mgrid)
        model = Medium(mgrid, [:vp, :rho])
        update!(model, [:vp, :rho], [vp0, rho0])
        fill!(model)

        # add some noise to starting model
        update!(model, [:vp, :rho], randn_perc = 0.5)

        # add perturbations
        for ellip_loc in [[500.0, 0.0], [0.0, 500.0], [-500.0, 0.0], [0.0, -500.0]]
            update!(
                model,
                [:vp, :rho],
                ellip_rad = 50.0,
                ellip_loc = ellip_loc,
                ellip_pert = 100.0,
            )
        end

    elseif ((attrib == :acou_homo2))
        vp0 = [1700.0, 2300.0] # bounds for vp
        rho0 = [1700.0, 2300.0] # density bounds
        mgrid = repeat([range(-1000.0, stop = 1000.0, length = 51)], 2)
        nz, nx = length.(mgrid)
        model = Medium(mgrid, [:vp, :rho])
        update!(model, [:vp, :rho], [vp0, rho0])
        fill!(model)

    elseif (attrib == :marmousi2)
        c = h5open(joinpath(marmousi_folder, "marmousi2.h5"), "r") do file
            vp = read(file, "vp")
            vs = read(file, "vs")
            rho = read(file, "rho")
            xgrid = read(file, "xgrid")
            zgrid = read(file, "zgrid")
            mgrid = [
                range(zgrid[1], stop = zgrid[end], length = size(vp, 1)),
                range(xgrid[1], stop = xgrid[end], length = size(vp, 2)),
            ]
            model = Medium(mgrid, [:vp, :rho, :vs])
            copyto!(model[:vp], vp)
            copyto!(model[:rho], rho)
            copyto!(model[:vs], vs)
            update!(model, bfrac)
        end
	elseif (attrib == :marmousi2_small)
		fmodel=Medium(:marmousi2)
		fmgrid=fmodel.mgrid
		mgrid=[range(500,stop=3500, step=step(fmgrid[1])), range(4000,step=step(fmgrid[2]),stop=13000)]
		model=update(fmodel,mgrid)

    elseif (attrib == :overthrust)
        c = h5open(joinpath(overthrust_folder, "overthrust.h5"), "r") do file
            vp = read(file, "vp")
            xgrid = read(file, "xgrid")
            ygrid = read(file, "ygrid")
            zgrid = read(file, "zgrid")
            mgrid = [
                range(zgrid[1], stop = zgrid[end], length = size(vp, 1)),
                range(xgrid[1], stop = xgrid[end], length = size(vp, 2)),
                range(xgrid[1], stop = xgrid[end], length = size(vp, 3)),
            ]
            model = Medium(mgrid, [:vp])
            copyto!(model[:vp], vp)
            update!(model, bfrac)
        end

        #
        #		vs0=Models.bounds(vs,bfrac); 
        #		rho0=Models.bounds(rho, bfrac);
        #		update!(model,[:vp,:rho,:vs],[vp0,rho0,vs0])
        #=
        	elseif(attrib == :seismic_marmousi2_high_res)
        		vp, h= IO.readsu(joinpath(marmousi_folder,"vp_marmousi-ii.su"))
        		vs, h= IO.readsu(joinpath(marmousi_folder,"vs_marmousi-ii.su"))
        		rho,  h= IO.readsu(joinpath(marmousi_folder,"density_marmousi-ii.su"))
        		vp .*= 1000.; vs .*= 1000.; #rho .*=1000
        		vp0=Models.bounds(vp,bfrac); 
        		vs0=Models.bounds(vs,bfrac); 
        		rho0=Models.bounds(rho, bfrac);
        		mgrid=[range(0.,stop=3500.,length=size(vp,1)),range(0., stop=17000., length=size(vp,2))]
        		model=Medium(mgrid,[:vp,:rho,:vs])
        		update!(model,[:vp,:rho,:vs],[vp0,rho0,vs0])
        		copyto!(model[:vp],vp)
        		copyto!(model[:rho],rho)
        		copyto!(model[:vs],vs)

        	elseif(attrib == :seismic_marmousi2_xwell)
        		model=Medium_trun(Seismic(:seismic_marmousi2_high_res), 
        				     zmin=1000., zmax=2000., xmin=8500., xmax=9500.,)
        		update!(model, bfrac) # adjuts bounds just inside the bounds 
        	elseif(attrib == :seismic_marmousi2_surf)
        		model=Medium_trun(Seismic(:seismic_marmousi2_high_res), 
        				     xmin=6000., xmax=12000.,)
        		update!(model, bfrac) # adjust bounds just inside the bounds 
        	elseif(attrib == :seismic_marmousi2_downhole)
        		model=Medium_trun(Seismic(:seismic_marmousi2_high_res), 
        				     xmin=9025., xmax=9125., zmin=1400., zmax=1600.,)
        		update!(model, bfrac) # adjust bounds just inside the bounds 
        	elseif(attrib == :seismic_marmousi2_rvsp)
        		model=Medium_trun(Seismic(:seismic_marmousi2_high_res), 
        				     xmin=8000., xmax=10000., zmax=1700.,zmin=500.)
        		update!(model, bfrac) # adjust bounds just inside the bounds 
        		=#

    else
        error("invalid attrib")
    end
    if (δ == 0.0)
        verbose && Models.print(model, string(attrib))
        return model
    elseif (δ > 0.0)
        mgrid_out = broadcast(x -> range(x[1], stop = x[end], step = δ), model.mgrid)
        model_out = update(model, mgrid_out)
        verbose && Models.print(model_out, string(attrib))
        return model_out
    else
        error("invalid δ")
    end


end


