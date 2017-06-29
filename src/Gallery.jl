__precompile__()

module Gallery

import SIT.IO
import SIT.Grid
import SIT.Models
import SIT.Acquisition
import SIT.Wavelets

global marmousi_folder="/math/home/pawbz/marmousi2/"

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
"""
function Seismic(attrib::Symbol, δ::Float64=0.0)
	if((attrib == :acou_homo1) | (attrib == :acou_homo2))
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
		vp, nz, nx = IO.readsu_data(fname=string(marmousi_folder,"vp_marmousi-ii_0.1.su"))
		vs, nz, nx = IO.readsu_data(fname=string(marmousi_folder,"vs_marmousi-ii_0.1.su"))
		ρ, nz, nx = IO.readsu_data(fname=string(marmousi_folder,"density_marmousi-ii_0.1.su"))
		vp .*= 1000.; vs .*= 1000.; #ρ .*=1000
		bound=0.1; vp0=zeros(2); vs0=zeros(2); ρ0=zeros(2);
		boundvp=bound*mean(vp); boundvs=bound*mean(vs); boundρ=bound*mean(ρ);
		vp0[1] = ((minimum(vp) - boundvp)<0.0) ? 0.0 : (minimum(vp) - boundvp)
		vp0[2] = maximum(vp)+boundvp
		vs0[1] = ((minimum(vs) - boundvs)<0.0) ? 0.0 : (minimum(vs) - boundvs)
		vs0[2] = maximum(vs)+boundvs
		ρ0[1] = ((minimum(ρ) - boundρ)<0.0) ? 0.0 : (minimum(ρ) - boundρ)
		ρ0[2] = maximum(ρ)+boundρ
		mgrid = Grid.M2D(0., 17000., 0., 3500.,nx,nz,40)
		model= Models.Seismic(vp0, vs0, ρ0, Models.χ(vp,vp0,1), Models.χ(vs,vs0,1), Models.χ(ρ,ρ0,1), mgrid)
	elseif(attrib == :seismic_marmousi2_high_res)
		vp, nz, nx = IO.readsu_data(fname=string(marmousi_folder,"vp_marmousi-ii.su"))
		vs, nz, nx = IO.readsu_data(fname=string(marmousi_folder,"vs_marmousi-ii.su"))
		ρ, nz, nx = IO.readsu_data(fname=string(marmousi_folder,"density_marmousi-ii.su"))
		vp .*= 1000.; vs .*= 1000.; #ρ .*=1000
		bound=0.1; vp0=zeros(2); vs0=zeros(2); ρ0=zeros(2);
		boundvp=bound*mean(vp); boundvs=bound*mean(vs); boundρ=bound*mean(ρ);
		vp0[1] = ((minimum(vp) - boundvp)<0.0) ? 0.0 : (minimum(vp) - boundvp)
		vp0[2] = maximum(vp)+boundvp
		vs0[1] = ((minimum(vs) - boundvs)<0.0) ? 0.0 : (minimum(vs) - boundvs)
		vs0[2] = maximum(vs)+boundvs
		ρ0[1] = ((minimum(ρ) - boundρ)<0.0) ? 0.0 : (minimum(ρ) - boundρ)
		ρ0[2] = maximum(ρ)+boundρ
		mgrid = Grid.M2D(0., 17000., 0., 3500.,nx,nz,40)
		model= Models.Seismic(vp0, vs0, ρ0, Models.χ(vp,vp0,1), Models.χ(vs,vs0,1), Models.χ(ρ,ρ0,1), mgrid)

	elseif(attrib == :seismic_marmousi2_box1)
		mgrid=Grid.M2D(8500.,9500., 1000., 2000.,5.,5.,40)
		marm_box1=Models.Seismic_zeros(mgrid)
		Models.Seismic_interp_spray!(Seismic(:seismic_marmousi2), marm_box1, :interp)
		model= marm_box1
	else
		error("invalid attrib")
	end
	if(δ==0.0)
		return model
	elseif(δ > 0.0)
		mgrid_out=Grid.M2D_resamp(model.mgrid,δ,δ,)
		model_out=Models.Seismic_zeros(mgrid_out)
		Models.Seismic_interp_spray!(model, model_out, :interp)
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
  * `=:surf` cross-well; but sources and receivers at unequal depths
  * `=:vsp` vertical seismic profiling
  * `=:rvsp`  reverse vertical seismic profiling
  * `=:downhole` downhole sources and receivers 

* `rand_flags::Vector{Bool}=[false, false]` : randomly or equally spaced?
"""
function Geom(mgrid::Grid.M2D, attrib::Symbol; nss=2, nr=2, rand_flags=[false, false])
	otx=(0.9*mgrid.x[1]+0.1*mgrid.x[end]); ntx=(0.1*mgrid.x[1]+0.9*mgrid.x[end]);
	otz=(0.9*mgrid.z[1]+0.1*mgrid.z[end]); ntz=(0.1*mgrid.z[1]+0.9*mgrid.z[end]);
	quatx = (0.75*mgrid.x[1]+0.25*mgrid.x[end]); quatz = (0.75*mgrid.z[1]+0.25*mgrid.z[end]) 
	tquatx = (0.25*mgrid.x[1]+0.75*mgrid.x[end]); tquatz = (0.25*mgrid.z[1]+0.75*mgrid.z[end]) 
	halfx = 0.5*(mgrid.x[1]+mgrid.x[end]);	halfz = 0.5*(mgrid.z[1]+mgrid.z[end]);
	if(attrib == :xwell)
		return Acquisition.Geom_fixed(otz, ntz, otx, otz, ntz, ntx, nss, nr, :vertical, :vertical, rand_flags)
	elseif(attrib == :surf)
		return Acquisition.Geom_fixed(otx, ntx, otz, otx, ntx, otz, nss, nr, :horizontal, :horizontal, rand_flags)
	elseif(attrib == :vsp)
		return Acquisition.Geom_fixed(otx, ntx, otz, otz, ntz, otx, nss, nr, :horizontal, :vertical, rand_flags)
	elseif(attrib == :rvsp)
		return Acquisition.Geom_fixed(otz, ntz, otx, otx, ntx, otz, nss, nr, :vertical, :horizontal, rand_flags)
	elseif(attrib == :downhole)
		return Acquisition.Geom_fixed(quatz, otz, quatx, quatz, otz, quatx, nss, nr, :vertical, :vertical, rand_flags)
	else
		error("invalid attrib")
	end
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
		wav = Wavelets.ricker(fqdom=10.0, tgrid=tgrid, tpeak=0.25, )
		return Acquisition.Src_fixed(nss, 1, 1, wav, tgrid)
	elseif(attrib == :acou_homo2)
		tgrid = M1D(attrib)
		wav = Wavelets.ricker(fqdom=3.0, tgrid=tgrid, tpeak=0.3, )
		return Acquisition.Src_fixed(nss, 1, 1, wav, tgrid)
	elseif(attrib == :vecacou_homo1)
		tgrid = M1D(:acou_homo1)
		wav = Wavelets.ricker(fqdom=10.0, tgrid=tgrid, tpeak=0.25, )
		return Acquisition.Src_fixed(nss, 1, 3, wav, tgrid)
	end
end


end # module
