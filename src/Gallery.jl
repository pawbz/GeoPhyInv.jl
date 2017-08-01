__precompile__()

module Gallery

import JuMIT.IO
import JuMIT.Grid
import JuMIT.Models
import JuMIT.Acquisition
import JuMIT.Wavelets

global marmousi_folder=joinpath(Pkg.dir("JuMIT"), "marmousi2")

"""
Gallery of `M2D` grids.

# Arguments 
* `attrib::Symbol` : 

# Outputs
* `attrib=:acou_homo1` : a square grid for with 201 samples in each dimension, with 50 PML 
		points; both X and Z vary from -1000 to 1000.
* `attrib=:acou_homo2` : a square grid for with 51 samples in each dimension, with 50 PML 
		points; both X and Z vary from -1000 to 1000. 
"""

function M2D(attrib::Symbol)
	if(attrib == :acou_homo1)
		return Grid.M2D(-1000.0,1000.0,-1000.0,1000.0,201,201,50)
	elseif(attrib == :acou_homo2)
		return Grid.M2D(-1000.0,1000.0,-1000.0,1000.0,51,51,50)
	else
		error("invalid attrib")
	end
end


"""
Gallery of `M1D` grids.

# Arguments 
* `attrib::Symbol` : 

# Outputs
* `attrib=:acou_homo1` : a time grid for with 1000 samples; maximum time is 2 s
* `attrib=:acou_homo1_long` : a time grid for with 1000 samples; maximum time is 4 s
* `attrib=:npow2samp1` : a sample npow2 grid with 16 samples
"""
function M1D(attrib::Symbol)
	if(attrib == :acou_homo1)
		return Grid.M1D(0.0,2.0,1000)
	elseif(attrib == :acou_homo1_long)
		return Grid.M1D(0.0,4.0,2000)
	elseif(attrib == :acou_homo2)
		return Grid.M1D(0.0,2.0,250)
	elseif(attrib == :npow2samp)
		return Grid.M1D(npow2=16,δ=0.0001)
	else
		error("invalid attrib")
	end
end


"""
Gallery of `Seismic` models.

# Arguments 
* `attrib::Symbol` : 

# Outputs
* `attrib=:acou_homo1` : an homogeneous acoustic model with `vp0=2000` and `ρ0=2000`
  * `attrib=:acou_homo2` : same as above, but with spatial sampling as 40 m (faster testing)
  * `attrib=:seismic_marmousi2` : marmousi model with lower resolution; ideal for surface seismic experiments
  * `attrib=:seismic_marmousi2_high_res` : marmousi model high resolution; slower to load
  * `attrib=:seismic_marmousi2_box1` : 1x1 kilometer box of marmousi model; ideal for crosswell, borehole seismic studies
  * `attrib=:seismic_marmousi2_box2` : box for surface seismic studies
  * `attrib=:seismic_marmousi2_box3` : 100x100 meter box of marmousi model at the center

"""
function Seismic(attrib::Symbol, δ::Float64=0.0)
	bfrac=0.1; 
	if((attrib == :acou_homo1))
		vp0 = [1500., 3500.] # bounds for vp
		vs0 = [1.0, 1.0] # dummy
		ρ0 = [1500., 3500.] # density bounds
		mgrid = M2D(attrib)
		model= Models.Seismic(vp0, vs0, ρ0,
		      fill(0.0, (mgrid.nz, mgrid.nx)),
		      fill(0.0, (mgrid.nz, mgrid.nx)),
		      fill(0.0, (mgrid.nz, mgrid.nx)),
		      mgrid)

	elseif((attrib == :acou_homo2))
		vp0 = [1700., 2300.] # bounds for vp
		vs0 = [1.0, 1.0] # dummy
		ρ0 = [1700., 2300.] # density bounds
		mgrid = M2D(attrib)
		model= Models.Seismic(vp0, vs0, ρ0,
		      fill(0.0, (mgrid.nz, mgrid.nx)),
		      fill(0.0, (mgrid.nz, mgrid.nx)),
		      fill(0.0, (mgrid.nz, mgrid.nx)),
		      mgrid)

	elseif(attrib == :seismic_marmousi2)
		vp, h= IO.readsu(joinpath(marmousi_folder,"vp_marmousi-ii_0.1.su"))
		vs, h= IO.readsu(joinpath(marmousi_folder,"vs_marmousi-ii_0.1.su"))
		ρ,  h= IO.readsu(joinpath(marmousi_folder,"density_marmousi-ii_0.1.su"))
		vp .*= 1000.; vs .*= 1000.; #ρ .*=1000
		vp0=Models.bounds(vp,bfrac); 
		vs0=Models.bounds(vs,bfrac); 
		ρ0=Models.bounds(ρ, bfrac);
		mgrid = Grid.M2D(0., 17000., 0., 3500.,size(vp,2),size(vp,1),40)
		model= Models.Seismic(vp0, vs0, ρ0, Models.χ(vp,vp0,1), Models.χ(vs,vs0,1), Models.χ(ρ,ρ0,1), mgrid)
	elseif(attrib == :seismic_marmousi2_high_res)
		vp, h= IO.readsu(joinpath(marmousi_folder,"vp_marmousi-ii.su"))
		vs, h= IO.readsu(joinpath(marmousi_folder,"vs_marmousi-ii.su"))
		ρ,  h= IO.readsu(joinpath(marmousi_folder,"density_marmousi-ii.su"))
		vp .*= 1000.; vs .*= 1000.; #ρ .*=1000
		vp0=Models.bounds(vp,bfrac); 
		vs0=Models.bounds(vs,bfrac); 
		ρ0=Models.bounds(ρ, bfrac);
		mgrid = Grid.M2D(0., 17000., 0., 3500.,size(vp,2),size(vp,1),40)
		model= Models.Seismic(vp0, vs0, ρ0, Models.χ(vp,vp0,1), Models.χ(vs,vs0,1), Models.χ(ρ,ρ0,1), mgrid)

	elseif(attrib == :seismic_marmousi2_box1)
		mgrid=Grid.M2D(8500.,9500., 1000., 2000.,5.,5.,40)
		marm_box1=Models.Seismic_zeros(mgrid)
		Models.Seismic_interp_spray!(Seismic(:seismic_marmousi2), marm_box1, :interp, :B1)
		model= Models.adjust_bounds!(marm_box1, bfrac)
	elseif(attrib == :seismic_marmousi2_box2)
		marm_full = Seismic(:seismic_marmousi2)
		mgrid=Grid.M2D(6000.,12000., marm_full.mgrid.z[1],
		 	marm_full.mgrid.z[end] ,5.,5.,40)
		marm_box2=Models.Seismic_zeros(mgrid)
		Models.Seismic_interp_spray!(marm_full, marm_box2, :interp, :B1)
		model= Models.adjust_bounds!(marm_box2, bfrac)
	elseif(attrib == :seismic_marmousi2_box3)
		mgrid=Grid.M2D(9000.,9100., 1500., 1700.,0.1,0.1,40)
		marm_box3=Models.Seismic_zeros(mgrid)
		Models.Seismic_interp_spray!(Seismic(:seismic_marmousi2_high_res), marm_box3, :interp, :B1)
		model= Models.adjust_bounds!(marm_box3, bfrac)
	else
		error("invalid attrib")
	end
	if(δ==0.0)
		Models.print(model,string(attrib))
		return model
	elseif(δ > 0.0)
		mgrid_out=Grid.M2D_resamp(model.mgrid,δ,δ,)
		model_out=Models.Seismic_zeros(mgrid_out)
		Models.Seismic_interp_spray!(model, model_out, :interp, :B1)
		Models.print(model_out,string(attrib))
		return model_out
	else
		error("invalid δ")
	end


end

"""
Gallery of acquisition geometries `Geom` using an input mesh `M2D`.
The sources and receivers are not placed anywhere on the edges of the mesh.

# Arguments 
* `mgrid::Grid.M2D` : a 2-D mesh
* `attrib::Symbol` : attribute decides output
  * `=:xwell` cross-well acquisition
  * `=:surf` surface acquisition
  * `=:vsp` vertical seismic profiling
  * `=:rvsp`  reverse vertical seismic profiling
  * `=:downhole` downhole sources and receivers 

* `rand_flags::Vector{Bool}=[false, false]` : randomly or equally spaced?
"""
function Geom(mgrid::Grid.M2D, attrib::Symbol; nss=2, nr=2, rand_flags=[false, false])
	otx=(0.9*mgrid.x[1]+0.1*mgrid.x[end]); ntx=(0.1*mgrid.x[1]+0.9*mgrid.x[end]);
	otwx=(0.95*mgrid.x[1]+0.05*mgrid.x[end]); ntwx=(0.05*mgrid.x[1]+0.95*mgrid.x[end]);
	otz=(0.9*mgrid.z[1]+0.1*mgrid.z[end]); ntz=(0.1*mgrid.z[1]+0.9*mgrid.z[end]);
	otwz=(0.95*mgrid.z[1]+0.05*mgrid.z[end]); ntwz=(0.05*mgrid.z[1]+0.95*mgrid.z[end]);
	quatx = (0.75*mgrid.x[1]+0.25*mgrid.x[end]); quatz = (0.75*mgrid.z[1]+0.25*mgrid.z[end]) 
	tquatx = (0.25*mgrid.x[1]+0.75*mgrid.x[end]); tquatz = (0.25*mgrid.z[1]+0.75*mgrid.z[end]) 
	halfx = 0.5*(mgrid.x[1]+mgrid.x[end]);	halfz = 0.5*(mgrid.z[1]+mgrid.z[end]);
	if(attrib == :xwell)
		geom=Acquisition.Geom_fixed(otz, ntz, otx, otz, ntz, ntx, nss, nr, :vertical, :vertical, rand_flags)
	elseif(attrib == :surf)
		geom=Acquisition.Geom_fixed(otx, ntx, otwz, otx, ntx, otwz, nss, nr, :horizontal, :horizontal, rand_flags)
	elseif(attrib == :vsp)
		geom=Acquisition.Geom_fixed(otx, ntx, otwz, otz, ntz, otx, nss, nr, :horizontal, :vertical, rand_flags)
	elseif(attrib == :rvsp)
		geom=Acquisition.Geom_fixed(otz, ntz, otx, otx, ntx, otwz, nss, nr, :vertical, :horizontal, rand_flags)
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
		tgrid = M1D(attrib)
		wav = Wavelets.ricker(10.0, tgrid, tpeak=0.25, )
		return Acquisition.Src_fixed(nss, 1, 1, wav, tgrid)
	elseif(attrib == :acou_homo2)
		tgrid = M1D(attrib)
		wav = Wavelets.ricker(3.0, tgrid, tpeak=0.3, )
		return Acquisition.Src_fixed(nss, 1, 1, wav, tgrid)
	elseif(attrib == :vecacou_homo1)
		tgrid = M1D(:acou_homo1)
		wav = Wavelets.ricker(10.0, tgrid, tpeak=0.25, )
		return Acquisition.Src_fixed(nss, 1, 3, wav, tgrid)
	end
end


end # module
