module Gallery

import SIT.Grid
import SIT.Models
import SIT.Acquisition
import SIT.Wavelets

"""
Gallery of `M2D` grids.

# Arguments 
* `attrib::Symbol` : 

# Outputs
* `attrib=:acou_homo1` : a square grid for with 201 samples in each dimension, with 50 PML 
		points; both X and Z vary from -1000 to 1000.
"""

function M2D(attrib::Symbol)
	if(attrib == :acou_homo1)
		return Grid.M2D(-1000.0,1000.0,-1000.0,1000.0,201,201,50)
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
* `attrib=:npow2samp1` : a sample npow2 grid with 16 samples
"""
function M1D(attrib::Symbol)
	if(attrib == :acou_homo1)
		return Grid.M1D(0.0,2.0,1000)
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
* `attrib=:seismic_marmousi2` : marmousi model with lower resolution; ideal for surface seismic experiments
* `attrib=:seismic_marmousi2_high_res` : marmousi model high resolution; slower to load
* `attrib=:seismic_marmousi2_box1` : 1x1 kilometer box of marmousi model; ideal for crosswell, borehole seismic studies
"""

function Seismic(attrib::Symbol)
	if(attrib == :acou_homo1)
		vp0 = 2000.0;
		vs0 = 0.0;
		ρ0 = 2000.0;
		mgrid = M2D(:acou_homo1)
		return Models.Seismic(vp0, vs0, ρ0,
		      fill(0.0, (mgrid.nz, mgrid.nx)),
		      fill(0.0, (mgrid.nz, mgrid.nx)),
		      fill(0.0, (mgrid.nz, mgrid.nx)),
		      mgrid)
	elseif(attrib == :seismic_marmousi2)
		vp, nz, nx = IO.readsu_data(fname="/home/pawbz/marmousi2/vp_marmousi-ii_0.1.su")
		vs, nz, nx = IO.readsu_data(fname="/home/pawbz/marmousi2/vs_marmousi-ii_0.1.su")
		ρ, nz, nx = IO.readsu_data(fname="/home/pawbz/marmousi2/density_marmousi-ii_0.1.su")
		vp0 = 1000. * mean(vp); vs0 = 1000. * mean(vs); ρ0 = mean(ρ)
		mgrid = Grid.M2D(0., 17000., 0., 3500.,nx,nz,40)
		return Models.Seismic(vp0, vs0, ρ0, 1000.*vp, 1000.*vs, ρ,
		      mgrid)
	elseif(attrib == :seismic_marmousi2_high_res)
		vp, nz, nx = IO.readsu_data(fname="/home/pawbz/marmousi2/vp_marmousi-ii.su")
		vs, nz, nx = IO.readsu_data(fname="/home/pawbz/marmousi2/vs_marmousi-ii.su")
		ρ, nz, nx = IO.readsu_data(fname="/home/pawbz/marmousi2/density_marmousi-ii.su")
		vp0 = 1000. * mean(vp); vs0 = 1000. * mean(vs); ρ0 = mean(ρ)
		mgrid = Grid.M2D(0., 17000., 0., 3500.,nx,nz,40)
		return Models.Seismic(vp0, vs0, ρ0, 1000.*vp, 1000.*vs, ρ,
		      mgrid)

	elseif(attrib == :seismic_marmousi2_box1)
		mgrid=Grid.M2D(8500.,9500., 1000., 2000.,5.,5.,40)
		return Models.Seismic_resamp(Seismic(:gallery_marmousi2), mgrid)
	else
		error("invalid attrib")
	end
end


"""
Gallery of acquisition geometries `Geom`.

# Arguments 
* `attrib::Symbol` : 

# Outputs
* `attrib=:acou_homo1` : a simple one source and one receiver configuration
"""
function Geom(attrib::Symbol)
	if(attrib == :acou_homo1)
		return Acquisition.Geom_fixed(-300.0,-300.0,-300.0,300.0,300.0,300.0,1,1)
	else
		error("invalid attrib")
	end
end

"""
Gallery of acquisition geometries `Geom` based on input `M2D`.

# Arguments 
* `attrib::Symbol` : 

# Outputs
* `attrib=:oneonev` : one source at (xmin, mean(z)) and one receiver at (xmax, mean(z))
* `attrib=:twotwov` : two vertical wells, two sources at xmin and two receivers at xmax
* `attrib=:tentenv` : two vertical wells, two sources at xmin and two receivers at xmax
"""
function Geom(mgrid::Grid.M2D,
	      attrib::Symbol
	     )
	quatx = (0.75*mgrid.x[1]+0.25*mgrid.x[end]); quatz = (0.75*mgrid.z[1]+0.25*mgrid.z[end]) 
	tquatx = (0.25*mgrid.x[1]+0.75*mgrid.x[end]); tquatz = (0.25*mgrid.z[1]+0.75*mgrid.z[end]) 
	halfx = 0.5*(mgrid.x[1]+mgrid.x[end]);	halfz = 0.5*(mgrid.z[1]+mgrid.z[end]);
	if(attrib == :oneonev)
		return Acquisition.Geom_fixed(
		      maximum(mgrid.z), maximum(mgrid.z), mean(mgrid.x),
		      minimum(mgrid.z), minimum(mgrid.z), mean(mgrid.x),
		      1,1,:vertical,:vertical
				)
	elseif(attrib == :twotwov)
		return Acquisition.Geom_fixed(
		quatz, tquatz, quatx, quatz, tquatz, tquatx,
		      2,2,:vertical,:vertical
				)
	elseif(attrib == :twotwodv)
		return Acquisition.Geom_fixed(
		quatz, halfz, quatx, halfz, tquatz, tquatx,
		      2,2,:vertical,:vertical
				)
	elseif(attrib == :twotenv)
		return Acquisition.Geom_fixed(
		quatz, tquatz, quatx, quatz, tquatz, tquatx,
		      2,10,:vertical,:vertical
				)
	elseif(attrib == :tentenv)
		return Acquisition.Geom_fixed(
		quatz, tquatz, mgrid.x[1], quatz, tquatz, mgrid.x[end],
		      10,10,:vertical,:vertical
				)
	elseif(attrib == :onefiftyv)
		return Acquisition.Geom_fixed(
	      mgrid.z[round(Int,0.5*mgrid.nz)], mgrid.z[round(Int,0.5*mgrid.nz)], mgrid.x[1],
	      mgrid.z[round(Int,0.25*mgrid.nz)], mgrid.z[round(Int,0.75*mgrid.nz)], mgrid.x[end],
		      1,50,:vertical,:vertical
				)
	elseif(attrib == :onefiftys)
		return Acquisition.Geom_fixed(
	      mgrid.x[round(Int,0.5*mgrid.nx)], mgrid.x[round(Int,0.5*mgrid.nx)], mgrid.z[1],
	      mgrid.x[round(Int,0.25*mgrid.nx)], mgrid.x[round(Int,0.75*mgrid.nx)], mgrid.z[1],
		      1,50,:horizontal,:horizontal
				)
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
		tgrid = M1D(:acou_homo1)
		wav = Wavelets.ricker(fqdom=10.0, tgrid=tgrid, tpeak=0.25, )
		return Acquisition.Src_fixed(nss, 1, 1, wav, tgrid)
	elseif(attrib == :vecacou_homo1)
		tgrid = M1D(:acou_homo1)
		wav = Wavelets.ricker(fqdom=10.0, tgrid=tgrid, tpeak=0.25, )
		return Acquisition.Src_fixed(nss, 1, 3, wav, tgrid)
	end
end


end # module
