module Gallery

using Signals
using Pkg
import JuMIT.IO
import JuMIT.Models
import JuMIT.Acquisition
import JuMIT.FWI
import JuMIT.Utils

global marmousi_folder=joinpath(@__DIR__, "marmousi2")


"""
Gallery of `Seismic` models.

# Arguments 
* `attrib::Symbol` : 

# Outputs
* `attrib=:acou_homo1` : an homogeneous acoustic model with `vp0=2000` and `ρ0=2000`
  * `attrib=:acou_homo2` : same as above, but with spatial sampling as 40 m (faster testing)
  * `attrib=:seismic_marmousi2` : marmousi model with lower resolution; ideal for surface seismic experiments
  * `attrib=:seismic_marmousi2_high_res` : marmousi model high resolution; slower to load
  * `attrib=:seismic_marmousi2_xwell` : boxed marmousi model ideal for crosswell experiments 
  * `attrib=:seismic_marmousi2_surf` : boxed marmousi2 for surface seismic experiments
  * `attrib=:seismic_marmousi2_downhole` : boxed marmousi2 for downhole seismic experiments  

"""
function Seismic(attrib::Symbol, δ::Float64=0.0; verbose=false)
	bfrac=0.1; 
	if((attrib == :acou_homo1))
		vp0 = [1500., 3500.] # bounds for vp
		vs0 = [1.0, 1.0] # dummy
		ρ0 = [1500., 3500.] # density bounds
		mgrid = repeat([range(-1000.0,stop=1000.0,length=201)],2)
		nz,nx=length.(mgrid)
		model= Models.Seismic(vp0, vs0, ρ0,
		      fill(0.0, (nz, nx)),
		      fill(0.0, (nz, nx)),
		      fill(0.0, (nz, nx)),
		      mgrid)

	elseif((attrib == :acou_homo2))
		vp0 = [1700., 2300.] # bounds for vp
		vs0 = [1.0, 1.0] # dummy
		ρ0 = [1700., 2300.] # density bounds
		mgrid = repeat([range(-1000.0,stop=1000.0,length=51)],2)
		nz,nx=length.(mgrid)
		model= Models.Seismic(vp0, vs0, ρ0,
		      fill(0.0, (nz, nx)),
		      fill(0.0, (nz, nx)),
		      fill(0.0, (nz, nx)),
		      mgrid)

	elseif(attrib == :seismic_marmousi2)
		vp, h= IO.readsu(joinpath(marmousi_folder,"vp_marmousi-ii_0.1.su"))
		vs, h= IO.readsu(joinpath(marmousi_folder,"vs_marmousi-ii_0.1.su"))
		ρ,  h= IO.readsu(joinpath(marmousi_folder,"density_marmousi-ii_0.1.su"))
		vp .*= 1000.; vs .*= 1000.; #ρ .*=1000
		vp0=Models.bounds(vp,bfrac); 
		vs0=Models.bounds(vs,bfrac); 
		ρ0=Models.bounds(ρ, bfrac);
		mgrid=[range(0.,stop=3500.,length=size(vp,1)),range(0., stop=17000., length=size(vp,2))]
		model= Models.Seismic(vp0, vs0, ρ0, Models.χ(vp,vp0,1), Models.χ(vs,vs0,1), Models.χ(ρ,ρ0,1), mgrid)
	elseif(attrib == :seismic_marmousi2_high_res)
		vp, h= IO.readsu(joinpath(marmousi_folder,"vp_marmousi-ii.su"))
		vs, h= IO.readsu(joinpath(marmousi_folder,"vs_marmousi-ii.su"))
		ρ,  h= IO.readsu(joinpath(marmousi_folder,"density_marmousi-ii.su"))
		vp .*= 1000.; vs .*= 1000.; #ρ .*=1000
		vp0=Models.bounds(vp,bfrac); 
		vs0=Models.bounds(vs,bfrac); 
		ρ0=Models.bounds(ρ, bfrac);
		mgrid=[range(0.,stop=3500.,length=size(vp,1)),range(0., stop=17000., length=size(vp,2))]
		model= Models.Seismic(vp0, vs0, ρ0, Models.χ(vp,vp0,1), Models.χ(vs,vs0,1), Models.χ(ρ,ρ0,1), mgrid)

	elseif(attrib == :seismic_marmousi2_xwell)
		model=Models.Seismic_trun(Seismic(:seismic_marmousi2_high_res), 
				     zmin=1000., zmax=2000., xmin=8500., xmax=9500.,)
		Models.adjust_bounds!(model, bfrac) # adjuts bounds just inside the bounds 
	elseif(attrib == :seismic_marmousi2_surf)
		model=Models.Seismic_trun(Seismic(:seismic_marmousi2_high_res), 
				     xmin=6000., xmax=12000.,)
		Models.adjust_bounds!(model, bfrac) # adjust bounds just inside the bounds 
	elseif(attrib == :seismic_marmousi2_downhole)
		model=Models.Seismic_trun(Seismic(:seismic_marmousi2_high_res), 
				     xmin=9025., xmax=9125., zmin=1400., zmax=1600.,)
		Models.adjust_bounds!(model, bfrac) # adjust bounds just inside the bounds 
	elseif(attrib == :seismic_marmousi2_rvsp)
		model=Models.Seismic_trun(Seismic(:seismic_marmousi2_high_res), 
				     xmin=8000., xmax=10000., zmax=1700.,zmin=500.)
		Models.adjust_bounds!(model, bfrac) # adjust bounds just inside the bounds 

	else
		error("invalid attrib")
	end
	if(δ==0.0)
		verbose && Models.print(model,string(attrib))
		return model
	elseif(δ > 0.0)
		mgrid_out=broadcast(x->range(x[1],stop=x[end],step=δ),model.mgrid)
		model_out=Models.Seismic_zeros(mgrid_out)
		Models.adjust_bounds!(model_out,model)
		Models.interp_spray!(model, model_out, :interp, :B1)
		verbose && Models.print(model_out,string(attrib))
		return model_out
	else
		error("invalid δ")
	end


end

"""
Gallery of acquisition geometries `Geom` using an input mesh.
The sources and receivers are not placed anywhere on the edges of the mesh.

# Arguments 
* `mgrid` : a 2-D mesh
* `attrib::Symbol` : attribute decides output
  * `=:xwell` cross-well acquisition
  * `=:surf` surface acquisition
  * `=:vsp` vertical seismic profiling
  * `=:rvsp`  reverse vertical seismic profiling
  * `=:downhole` downhole sources and receivers 

# Keyword Arguments

* `nss=2` : number of supersources
* `nr=2` : number of receivers per supersource
* `rand_flags::Vector{Bool}=[false, false]` : randomly or equally spaced supersources and receivers.
"""
function Geom(mgrid, attrib::Symbol; nss=2, nr=2, rand_flags=[false, false])
	otx=(0.9*mgrid[2][1]+0.1*mgrid[2][end]); ntx=(0.1*mgrid[2][1]+0.9*mgrid[2][end]);
	otwx=(0.95*mgrid[2][1]+0.05*mgrid[2][end]); ntwx=(0.05*mgrid[2][1]+0.95*mgrid[2][end]);
	otz=(0.9*mgrid[1][1]+0.1*mgrid[1][end]); ntz=(0.1*mgrid[1][1]+0.9*mgrid[1][end]);
	otwz=(0.95*mgrid[1][1]+0.05*mgrid[1][end]); ntwz=(0.05*mgrid[1][1]+0.95*mgrid[1][end]);
	quatx = (0.75*mgrid[2][1]+0.25*mgrid[2][end]); quatz = (0.75*mgrid[1][1]+0.25*mgrid[1][end]) 
	tquatx = (0.25*mgrid[2][1]+0.75*mgrid[2][end]); tquatz = (0.25*mgrid[1][1]+0.75*mgrid[1][end]) 
	halfx = 0.5*(mgrid[2][1]+mgrid[2][end]);	halfz = 0.5*(mgrid[1][1]+mgrid[1][end]);
	if(attrib == :xwell)
		geom=Acquisition.Geom_fixed(otz, ntz, otx, otz, ntz, ntx, nss, nr, :vertical, :vertical, rand_flags)
	elseif(attrib == :surf)
		geom=Acquisition.Geom_fixed(otx, ntx, otwz, otx, ntx, otwz, nss, nr, :horizontal, :horizontal, rand_flags)
	elseif(attrib == :vsp)
		geom=Acquisition.Geom_fixed(otx, ntx, otwz, quatz, ntz, otx, nss, nr, :horizontal, :vertical, rand_flags)
	elseif(attrib == :rvsp)
		geom=Acquisition.Geom_fixed(quatz, ntz, otx, otx, ntx, otwz, nss, nr, :vertical, :horizontal, rand_flags)
	elseif(attrib == :downhole)
		geom=Acquisition.Geom_fixed(quatz, otz, quatx, quatz, otz, quatx, nss, nr, :vertical, :vertical, rand_flags)
	else
		error("invalid attrib")
	end
	Acquisition.print(geom, string(attrib))
	return geom
end

"""
Gallery of source signals `Src`.

# Arguments 
* `attrib::Symbol` : 
* `nss::Int64=1` : number of supersources

# Outputs
* `attrib=:acou_homo1` : 
"""
function Src(attrib::Symbol, nss::Int64=1)
	if(attrib == :acou_homo1)
		tgrid=range(0.0,stop=2.0,length=1000)
		wav = Utils.Wavelets.ricker(10.0, tgrid, tpeak=0.25, )
		return Acquisition.Src_fixed(nss, 1, [:P], wav, tgrid)
	elseif(attrib == :acou_homo2)
		tgrid=range(0.0,stop=2.0,length=250)
		wav = Utils.Wavelets.ricker(3.0, tgrid, tpeak=0.3, )
		return Acquisition.Src_fixed(nss, 1, [:P], wav, tgrid)
	elseif(attrib == :vecacou_homo1)
		tgrid=range(0.0,stop=2.0,length=1000)
		wav = Utils.Wavelets.ricker(10.0, tgrid, tpeak=0.25, )
		return Acquisition.Src_fixed(nss, 1, [:P, :Vx, :Vz], wav, tgrid)
	end
end


for file in ["fwi.jl", "fdtd.jl"]   
	fn=joinpath(@__DIR__, file)
	include(fn)
	#Revise.track(@__MODULE__,fn)
end

		 	
end # module
