module Inversion

import SIT.Models
import SIT.Grid
import SIT.Acquisition
import SIT.Data
import SIT.Misfits
using Optim
import SIT.Fdtd


"""
Inversion Parameters, i.e., factors that 
are fixed throughout the inversion
"""
type Param
	"refe"
	vp0::Float64
	vs0::Float64
	ρ0::Float64
	mgrid::Grid.M2D
	igrid::Grid.M2D
	acqsrc::Acquisition.Src
	acqgeom::Acquisition.Geom
	tgrid::Grid.M1D
end


function Param(mod::Models.Seismic,acqsrc::Acquisition.Src,acqgeom::Acquisition.Geom,tgrid::Grid.M1D)
	return Param(mod.vp0, mod.vs0, mod.ρ0, mod.mgrid, mod.mgrid, acqsrc, acqgeom,tgrid)
end

"""
FWI using Optim

# Arguments


"""
function xfwi(
	      x::Vector{Float64}, 
	      pa::Param;
	      model::Models.Seismic=nothing,
	      model_init::Models.Seismic=nothing)

		# generate "observed" data using same forward modelling scheme
		buffer = Fdtd.fdtd_mod(model=model, acqgeom=[pa.acqgeom], 
					  acqsrc=[pa.acqsrc], 
					  tgridmod=pa.tgrid, 
					  boundary_save_flag=true);
		dd = deepcopy(buffer[1][1]);
		#buffer[1] = Data.TD_zeros(buffer[1]);

		#dd = dobs;
		#buffer()
	# convert initial model to the inversion variable
	x = zeros(2*model_init.mgrid.nz*model_init.mgrid.nx);
	Seismic_x!(model_init, x, pa, -1)

	last_x = similar(x)
	df = OnceDifferentiable(x -> func(x, buffer, x, dd, pa),
			  (x, storage) -> grad!(x, storage, buffer, last_x, dd, pa))
	res = optimize(df, 
		       x,
		       ConjugateGradient(),
		     Optim.Options(g_tol = 1e-12,
		     iterations = 10,
		  store_trace = true,
	       show_trace = true))
end

function forward_modelling!(x::Vector{Float64},
			    last_x::Vector{Float64}, buffer, pa::Param
			   )
	if(x!=last_x)
		println("forwaerd mide")
		copy!(last_x, x)

		# preallocate model
		mod = Models.Seismic_zeros(pa.mgrid);
		# x to model
		Seismic_x!(mod, x, pa, -1)		

		acq = pa.acqgeom;
		# generate modelled data
		buffer = Fdtd.fdtd_mod(model=mod, acqgeom=[acq], 
					  acqsrc=[pa.acqsrc], 
					  tgridmod=pa.tgrid, 
					  verbose=true,
					  boundary_save_flag=true);
	end
end

"""
Return functional and gradient of the CLS objective 
"""
function func(x::Vector{Float64}, buffer, last_x::Vector{Float64},
		       dobs::Data.TD, pa::Param)
	println("func")
	forward_modelling!(x, last_x, buffer, pa)
	# compute misfit and δdcal
	f, δdcal = Misfits.TD(buffer[1][1], dobs)
	return f
end

function grad!(x::Vector{Float64}, storage::Vector{Float64}, buffer, 
	       last_x::Vector{Float64}, dobs::Data.TD, pa::Param)
	println("grad")

	forward_modelling!(x, last_x, buffer, pa)

	# compute misfit and δdcal
	f, δdcal = Misfits.TD(buffer[1][1], dobs)

	# only compute if gradient necessary (saves time)
	# adjoint sources
	adjacq = AdjGeom(δdat.acqgeom);
	adjsrc = AdjSrc(δdcal)

	# source sinks
	acqsrcsink = deepcopy(pa.acqsrc); 
	acqsrcsink.wav = [-1.0.*flipdim(pa.acqsrc.wav[i,j],1) for i in 1:pa.acqsrc.nss, j in 1:pa.acqsrc.nfield];
	
	# adjoint simulation
	rectemp, btemp, gmod = Fdtd.fdtd_mod(npropwav=2, model=mod, model0=mod, 
		     acqgeom=[acq,adjacq], acqsrc=[acqsrcsink, adjsrc], src_flags=[2.0, -2.0], 
		     tgridmod=dobs.tgrid, grad_out_flag=true, boundary_in=buffer[2], verbose=true)

	Seismic_storage!(gmod,storage,pa,-1) 
	
end


"""
Convert `Seismic` model to x and vice versa
"""
function Seismic_x!(mod::Models.Seismic, 
		       x::Vector{Float64}, 
		       pa::Param, 
		       flag::Int64)
	if(flag ==1) # convert mod to x
		nznx = pa.mgrid.nz*pa.mgrid.nx
		# use igrid later
		x = zeros(2*nznx);
		x[1:nznx] = vec(mod.χvp);
		x[nznx+1:2*nznx] = vec(mod.χρ);
	elseif(flag == -1) # convert x to mod
		nznx = pa.mgrid.nz*pa.mgrid.nx # change to igrid later
		mod.χvp = reshape(x[1:nznx], pa.mgrid.nz, pa.mgrid.nx)	
		mod.χρ = reshape(x[nznx+1:2*nznx], pa.mgrid.nz, pa.mgrid.nx)	
		mod.χvs = zeros(mod.χvp);
		mod.vp0=pa.vp0; mod.vs0=pa.vs0; mod.ρ0=pa.ρ0;
		mod.mgrid=pa.mgrid
	else
		error("invalid flag")
	end
end

function Seismic_gx!(gmod::Models.Seismic,
		       gx::Vector{Float64},
		       pa::Param,
		       flag::Int64,
		      )
	if(flag ==1) # convert gmod to gx
		nznx = pa.mgrid.nz*pa.mgrid.nx
		gx = zeros(2*nznx);
		gx[1:nznx] = vec(gmod.χvp);
		gx[nznx+1:2*nznx] = vec(gmod.χρ);
	elseif(flag == -1) # convert gx to mod
		nznx = pa.mgrid.nz*pa.mgrid.nx
		gχvp = reshape(gx[1:nznx], pa.mgrid.nz, pa.mgrid.nx)	
		gχρ = reshape(gx[nznx+1:2*nznx], pa.mgrid.nz, pa.mgrid.nx)	
		gχvs = zeros(gχvp)
		gmod = Models.Seismic(pa.vp0, pa.vs0, pa.ρ0, gχvp, gχvs, gχρ, pa.mgrid)
	else
		error("invalid flag")
	end
end

function AdjSrc(δdat::Data.TD,
	       )
	adjacq = AdjGeom(δdat.acqgeom);
	return Acquisition.Src(adjacq.nss,adjacq.ns,δdat.nfield,
		  [flipdim(δdat.d[i,j],1) for i in 1:adjacq.nss, j in 1:δdat.nfield], δdat.tgrid);
end

"""

Modify the input acquisition geometry 
such that the adjoint source time functions can 
be propagated from the receiver positions.
The number of supersources will remain the same.
All the recievers will be fired as simultaneous sources.
"""
function AdjGeom(geomin::Acquisition.Geom)
	geomout = deepcopy(geomin);
	geomout.sx = geomin.rx; geomout.sz = geomin.rz;
	geomout.ns = geomin.nr; 

	return geomout
end


end # module
