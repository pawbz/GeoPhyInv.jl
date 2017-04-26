module Inversion

import SIT.Models
import SIT.Grid
import SIT.Acquisition
import SIT.Data
import SIT.Coupling
import SIT.Misfits
using Optim
import SIT.Fdtd

buffer = (Array(Data.TD,1),(Array(Float64,3),Array(Float64,5)),Models.Seismic)

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
	attrib_mod::Symbol
	attrib_inv::Symbol
	model0::Models.Seismic
end


function Param(mod::Models.Seismic,acqsrc::Acquisition.Src,acqgeom::Acquisition.Geom,tgrid::Grid.M1D,
	      attrib_mod::Symbol, attrib_inv::Symbol, model0::Models.Seismic)
	return Param(mod.vp0, mod.vs0, mod.ρ0, mod.mgrid, mod.mgrid, acqsrc, acqgeom,tgrid,
	      attrib_mod, attrib_inv, model0)
end

"""
FWI using Optim

# Arguments
"""
function xfwi(
	      pa::Param;
	      model::Models.Seismic=nothing,
	      model_init::Models.Seismic=nothing,
	      dobs::Data.TD=Data.TD_zeros(1,pa.tgrid,pa.acqgeom)
	      )

	# convert initial model to the inversion variable
	x = zeros(2*model_init.mgrid.nz*model_init.mgrid.nx);
	last_x = similar(x)

	# replace x with actual model
	Seismic_x!(model, x, pa, 1)

	# generate observed data
	if(Data.TD_iszero(dobs))
		# using same forward modelling scheme
		dobs = update_buffer!(x, last_x, pa)
	end

	last_x = rand(size(x)) # reset last_x
	# replace x with initial model
	Seismic_x!(model_init, x, pa, 1)

	if(pa.attrib_inv == :cls)
		df = OnceDifferentiable(x -> func(x, last_x, dobs, pa),
			  (x, storage) -> grad!(x, storage, last_x, dobs, pa))
		res = optimize(df, 
			       x,
			       ConjugateGradient(),
			     Optim.Options(g_tol = 1e-12,
			     iterations = 10,
			  store_trace = true,
		       show_trace = true))
		return res, dobs
	elseif(pa.attrib_inv == :migr)
		gx = similar(x);
		grad!(x, gx, last_x, dobs, pa)

		# convert gradient vector to model
		gmod = deepcopy(model); Seismic_gx!(gmod,gx,pa,-1)

		return gmod, dobs
	else
		error("invalid pa.attrib_inv")
	end
end

"""
SRCF

# Arguments
"""
function srcf(
	      pa::Param;
	      model::Models.Seismic=nothing,
	      model_init::Models.Seismic=nothing,
	      dobs::Data.TD=Data.TD_zeros(1,pa.tgrid,pa.acqgeom)
	      )

	# convert initial model to the inversion variable
	x = zeros(2*model_init.mgrid.nz*model_init.mgrid.nx);
	last_x = similar(x)

	# replace x with actual model
	Seismic_x!(model, x, pa, 1)

	# generate observed data
	if(Data.TD_iszero(dobs))
		# using same forward modelling scheme
		dobs = update_buffer!(x, last_x, pa)
	end

	last_x = rand(size(x)) # reset last_x
	# replace x with initial model
	Seismic_x!(model_init, x, pa, 1)

	if(pa.attrib_inv == :cls)
		df = OnceDifferentiable(x -> func(x, last_x, dobs, pa),
			  (x, storage) -> grad!(x, storage, last_x, dobs, pa))
		res = optimize(df, 
			       x,
			       ConjugateGradient(),
			     Optim.Options(g_tol = 1e-12,
			     iterations = 10,
			  store_trace = true,
		       show_trace = true))
		return res, dobs
	elseif(pa.attrib_inv == :migr)
		gx = similar(x);
		grad!(x, gx, last_x, dobs, pa)

		# convert gradient vector to model
		gmod = deepcopy(model); Seismic_gx!(gmod,gx,pa,-1)

		return gmod, dobs
	else
		error("invalid pa.attrib_inv")
	end
end


"""
Update the buffer, which consists of the modelled data
and boundary values for adjoint calculation.

# Arguments

* `x::Vector{Float64}` : inversion variable
* `last_x::Vector{Float64}` : buffer is only updated when x!=last_x, and modified such that last_x=x
* `pa::Param` : parameters that are constant during the inversion 

"""
function update_buffer!(x::Vector{Float64},
			    last_x::Vector{Float64}, pa::Param
			   )
	if(x!=last_x)
		println("updating buffer")
		println(maximum(x))
		copy!(last_x, x)

		# preallocate model and x --> mod;
		model = Models.Seismic_zeros(pa.mgrid); Seismic_x!(model, x, pa, -1)		

		# generate modelled data, border values, etc.
		if(pa.attrib_mod == :fdtd)
			global buffer = Fdtd.mod(model=model, acqgeom=[pa.acqgeom], 
						  acqsrc=[pa.acqsrc], 
						  tgridmod=pa.tgrid, 
						  boundary_save_flag=true);
		elseif(pa.attrib_mod == :fdtd_born)
			global buffer = Fdtd.mod(npropwav=2, model=pa.model0, model_pert=model, 
			      acqgeom=[pa.acqgeom,pa.acqgeom], acqsrc=[pa.acqsrc,pa.acqsrc], 
			      src_flags=[2.0, 0.0], recv_flags = [0.0, 2.0], 
			  tgridmod=pa.tgrid, verbose=true, boundary_save_flag=true, born_flag=true);
		elseif(pa.attrib_mod == :fdtd_hborn)
			# storing boundary values 
			hbuffer = Fdtd.mod(model=model, acqgeom=[pa.acqgeom], 
						  acqsrc=[pa.acqsrc], 
						  tgridmod=pa.tgrid, 
						  boundary_save_flag=true);

			# source sinks
			acqsrcsink = deepcopy(pa.acqsrc); 
			acqsrcsink.wav = [-1.0.*flipdim(pa.acqsrc.wav[i,j],1) for 
				   i in 1:pa.acqsrc.nss, j in 1:pa.acqsrc.nfield];


			global buffer = Fdtd.mod(boundary_in=hbuffer[2],
			      npropwav=2,model=pa.model0, model_pert=model, 
			       acqgeom=[pa.acqgeom,pa.acqgeom], 
			       acqsrc=[acqsrcsink,pa.acqsrc], src_flags=[2.0, 0.0], 
			       recv_flags=[0.0, 2.0],
			   tgridmod=pa.tgrid, verbose=true, boundary_save_flag=true, born_flag=true);
		else
			error("invalid pa.attrib_mod")
		end
		dat = deepcopy(buffer[1][1])

		return dat
	end
end

"""
Return functional and gradient of the CLS objective 
"""
function func(x::Vector{Float64},  last_x::Vector{Float64},
		       dobs::Data.TD, pa::Param)
	println("func+++++++")
	update_buffer!(x, last_x,  pa)
	# compute misfit and δdcal
	f, δdcal = Misfits.TD(buffer[1][1], dobs)
	return f
end

"""
Return functional and gradient of the CLS objective 
"""
function func_Coupling(x::Vector{Float64}, dcal, dobs::Data.TD, pa::Param)
	# x to w 
	w = Coupling.TD_delta(pa.tlagssf, 0.0, pa.tgrid.δx, pa.acqgeom)


	Coupling_x!(w, x, pa, -1)
	s = Data.TDcoup(dcal, w)
	# compute misfit and δdcal
	f, δdcal = Misfits.TD(s, dobs)
	return f
end

function grad!(x::Vector{Float64}, storage::Vector{Float64}, 
	       last_x::Vector{Float64}, dobs::Data.TD, pa::Param)
	println("grad++++++++")
	println(maximum(x))

	update_buffer!(x, last_x, pa)

	# compute misfit and δdcal
	f, δdcal = Misfits.TD(buffer[1][1], dobs)

	# only compute if gradient necessary (saves time)
	# adjoint sources
	adjacq = AdjGeom(δdcal.acqgeom);
	adjsrc = AdjSrc(δdcal)


	# preallocate model
	mod = Models.Seismic_zeros(pa.mgrid);
	# x to model
	Seismic_x!(mod, x, pa, -1)		

	if((pa.attrib_mod == :fdtd) | (pa.attrib_mod == :fdtd_born))
		# source sinks
		acqsrcsink = deepcopy(pa.acqsrc); 
		acqsrcsink.wav = [-1.0.*flipdim(pa.acqsrc.wav[i,j],1) for 
			   i in 1:pa.acqsrc.nss, j in 1:pa.acqsrc.nfield];

		
		# adjoint simulation
		adj = Fdtd.mod(npropwav=2, model=mod, model_pert=mod, 
		     acqgeom=[pa.acqgeom,adjacq], acqsrc=[acqsrcsink, adjsrc], src_flags=[2.0, -2.0], 
		     tgridmod=dobs.tgrid, grad_out_flag=true, boundary_in=buffer[2], verbose=true)
	elseif(pa.attrib_mod == :fdtd_hborn)
		# adjoint simulation
		adj = Fdtd.mod(npropwav=2, model=mod, model_pert=mod, 
		     acqgeom=[pa.acqgeom,adjacq], acqsrc=[pa.acqsrc, adjsrc], src_flags=[2.0, -2.0], 
		     tgridmod=dobs.tgrid, grad_out_flag=true, verbose=true)
	else
		error("invalid pa.attrib_mod")
	end

	# convert gradient model to gradient vector
	Seismic_gx!(adj[3],storage,pa,1) 
	return storage
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
		x[1:nznx] = copy(vec(mod.χvp));
		x[nznx+1:2*nznx] = copy(vec(mod.χρ));
	elseif(flag == -1) # convert x to mod
		nznx = pa.mgrid.nz*pa.mgrid.nx # change to igrid later
		mod.χvp = copy(reshape(x[1:nznx], pa.mgrid.nz, pa.mgrid.nx))	
		mod.χρ = copy(reshape(x[nznx+1:2*nznx], pa.mgrid.nz, pa.mgrid.nx))	
		mod.χvs =copy(zeros(mod.χvp))
		mod.vp0 = copy(pa.vp0); mod.vs0 = copy(pa.vs0); 
		mod.ρ0 = copy(pa.ρ0);
		mod.mgrid = deepcopy(pa.mgrid)
	else
		error("invalid flag")
	end
end


"""
Convert gradient vector to `Seismic` type and vice versa
This will be different from the previous one, once 
the parameterizations come in
"""
function Seismic_gx!(gmod::Models.Seismic,
		       gx::Vector{Float64},
		       pa::Param,
		       flag::Int64,
		      )
	if(flag ==1) # convert gmod to gx
		nznx = pa.mgrid.nz*pa.mgrid.nx
		gx[1:nznx] = copy(vec(gmod.χvp))
		gx[nznx+1:2*nznx] = copy(vec(gmod.χρ))
	elseif(flag == -1) # convert gx to mod
		nznx = pa.mgrid.nz*pa.mgrid.nx
		gmod.χvp = copy(reshape(gx[1:nznx], pa.mgrid.nz, pa.mgrid.nx))
		gmod.χρ = copy(reshape(gx[nznx+1:2*nznx], pa.mgrid.nz, pa.mgrid.nx))
		gmod.χvs = copy(zeros(gmod.χvp))
		gmod.vp0 = copy(pa.vp0); gmod.vs0 = copy(pa.vs0); 
		gmod.ρ0 = copy(pa.ρ0);
		gmod.mgrid = deepcopy(pa.mgrid)
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

"""
Convert `Seismic` model to x and vice versa
"""
function Coupling_x!(w::Coupling.TD, 
		   x::Vector{Float64}, 
		   pa::Param, 
		   flag::Int64)
	if(flag ==1) # convert w to x
		x[1:pa.tgridssf.nx] = copy(vec(w.ssf[1,1][:]));
	elseif(flag == -1) # convert x to mod
		w.ssf = copy([x for iss=1:acqgeom.nss, ifield=1:nfield])
	else
		error("invalid flag")
	end
end



end # module
