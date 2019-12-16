
global const marmousi_folder=joinpath(@__DIR__, "marmousi2")

"""
```julia
mod=Medium(attrib)
```
As of now, only seismic models are predefined in this package. Choose `attrib::Symbol`

* `=:acou_homo1` : a test homogeneous acoustic model
* `=:acou_homo2` : a test homogeneous acoustic model, but with coarser spatial sampling (faster testing)
* `=:seismic_marmousi2` : marmousi2 model with lower resolution; useful for surface seismic experiments
"""
function Medium(attrib::Symbol, δ::Real=0.0; verbose=false)
	bfrac=0.1; 
	δ=Float64(δ)
	if((attrib == :acou_homo1))
		vp0 = [1500., 3500.] # bounds for vp
		rho0 = [1500., 3500.] # density bounds
		mgrid = repeat([range(-1000.0,stop=1000.0,length=201)],2)
		nz,nx=length.(mgrid)
		model=Medium(mgrid,[:vp,:rho])
		update!(model,[:vp,:rho],[vp0,rho0])
		fill!(model)

	elseif((attrib == :acou_homo2))
		vp0 = [1700., 2300.] # bounds for vp
		rho0 = [1700., 2300.] # density bounds
		mgrid = repeat([range(-1000.0,stop=1000.0,length=51)],2)
		nz,nx=length.(mgrid)
		model=Medium(mgrid,[:vp,:rho])
		update!(model,[:vp,:rho],[vp0,rho0])
		fill!(model)

	elseif(attrib == :seismic_marmousi2)
		vp, h= IO.readsu(joinpath(marmousi_folder,"vp_marmousi-ii_0.1.su"))
		vs, h= IO.readsu(joinpath(marmousi_folder,"vs_marmousi-ii_0.1.su"))
		rho,  h= IO.readsu(joinpath(marmousi_folder,"density_marmousi-ii_0.1.su"))
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

	else
		error("invalid attrib")
	end
	if(δ==0.0)
		verbose && Models.print(model,string(attrib))
		return model
	elseif(δ > 0.0)
		mgrid_out=broadcast(x->range(x[1],stop=x[end],step=δ),model.mgrid)
		model_out=Medium_zeros(mgrid_out)
		update!(model_out,model)
		Models.interp_spray!(model, model_out, :interp, :B1)
		verbose && Models.print(model_out,string(attrib))
		return model_out
	else
		error("invalid δ")
	end


end


