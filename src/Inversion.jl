__precompile__()

module Inversion

import SIT.Models
import SIT.Grid
import SIT.Acquisition
import SIT.Data
import SIT.Coupling
import SIT.Misfits
import SIT.Fdtd
using Optim


"""
Inversion Parameters, i.e., factors that 
are fixed throughout the inversion

# Fields

* `mgrid::Grid.M2D` : modelling grid
* `igrid::Grid.M2D` : inversion grid
* `acqsrc::Acquisition.Src` : base source wavelet for modelling data
* `acqgeom::Acquisition.Geom` : acquisition geometry
* `tgrid::Grid.M1D` : 
* `attrib_mod::Symbol`
* `model_obs` : model used for generating observed data
* `model0` : background velocity model (only used during Born modeling and inversion)
* `parameterization` : a vector of Symbols specifying parameterization of the inversion vector
* `verbose` : print level
* `attrib` : synthetic or real

TODO: add an extra attribute for coupling functions inversion and modelling
"""
type Param
	"base source wavelet"
	acqsrc::Acquisition.Src
	"acquisition geometry"
	acqgeom::Acquisition.Geom
	"modeling attribute"
	attrib_mod::Symbol
	"inversion attribute"
	attrib_inv::Symbol
	"Seismic model on modeling grid"
	modm::Models.Seismic
	"background Seismic model"
	modm0::Models.Seismic
	"Seismic on inversion grid"
	modi::Models.Seismic
	"parameterization"
	parameterization::Vector{Symbol}
	"calculated data"
	dcal::Data.TD
	"observed data"
	dobs::Data.TD
	"coupling functions"
	w::Coupling.TD
	verbose::Bool
	attrib::Symbol
end


"""
Constructor for `Param`

* `model` : initial model, the reference values are used from this

"Seismic model used to generate 'observed' synthetic data"
modm_obs::Models.Seismic

"source wavelet used to generate 'observed' synthetic data"
acqsrc_obs::Acquisition.Src
"""
function Param(
	       acqsrc::Acquisition.Src,
	       acqgeom::Acquisition.Geom,
	       tgrid::Grid.M1D,
	       attrib_mod::Symbol, 
	       attrib_inv::Symbol,
	       modm::Models.Seismic;
	       # other optional 
	       igrid::Grid.M2D=modm.mgrid,
	       dobs::Data.TD=Data.TD_zeros(1,tgrid,acqgeom),
	       tlagssf::Float64=0.0,
	       tlagrf::Float64=0.0,
	       acqsrc_obs::Acquisition.Src=acqsrc,
	       modm_obs::Models.Seismic=modm,
	       modm0::Models.Seismic=modm,
	       mod_inv_parameterization::Vector{Symbol}=[:χKI, :χρI],
	       verbose::Bool=false,
	       attrib::Symbol=:synthetic
	       )

	nfield=1

	# create modi according to igrid and  interpolation of modm
	modi = Models.Seismic_zeros(igrid);
	Models.Seismic_interp_spray!(modm, modi, :interp)


	pa = Param(acqsrc, acqgeom, attrib_mod, attrib_inv, 
	      modm, modm0, modi, mod_inv_parameterization,
	      Data.TD_zeros(nfield,tgrid,acqgeom), dobs, 
	      Coupling.TD_delta(tlagssf, tlagrf, tgrid.δx, nfield, acqgeom), verbose, attrib)

	# some allocations for modelling
	buffer1, buffer2, buffer3 = update_buffer_zeros(pa)
	x = zeros(xfwi_ninv(pa)); # dummy
	last_x = rand(size(x)) # dummy last_x

	# intial modelling dcal
	update_buffer!(buffer1, buffer2, buffer3, x, last_x, pa, pa.modm) # x is dummy here
	pa.dcal = deepcopy(buffer1[1])

	# generate observed data if attrib is synthetic
	if((attrib == :synthetic) & Data.TD_iszero(pa.dobs))
		# change source wavelets for modelling observed data
		pa.acqsrc = deepcopy(acqsrc_obs)
		last_x = rand(size(x)) # reset last_x
		update_buffer!(buffer1, buffer2, buffer3, x, last_x, pa, modm_obs) # x is dummy here
		pa.dobs = deepcopy(buffer1[1])
		# change back the source wavelets
		pa.acqsrc = deepcopy(acqsrc)
	elseif(Data.TD_iszero(pa.dobs) & (attrib == :real))
		error("input observed data for real data inversion")
	else
		error("invalid attrib")
	end

	return pa
end

"Return the number of inversion variables `Param`"
function xfwi_ninv(pa::Param)
	return 2*pa.modi.mgrid.nz*pa.modi.mgrid.nx
end

"Return the number of inversion variables `Param`"
function wfwi_ninv(pa::Param)
	return  pa.w.tgridssf.nx
end

"""
FWI using Optim, updates pa.modm and pa.dcal

# Arguments
* `model:Models.Seismic` : initial model and it is updated
* `dobs::Data.TD=Data.TD_zeros(1,pa.dcal.tgrid,pa.acqgeom)` : input observed data
"""
function xfwi!(pa::Param, extended_trace::Bool=true, time_limit=Float64=2.0*60.)

	# convert initial model to the inversion variable
	x = zeros(xfwi_ninv(pa));
	last_x = similar(x)
	pa.verbose ? println("xfwi: number of inversion variables:\t", length(x)) : nothing

	buffer1, buffer2, buffer3 = update_buffer_zeros(pa)

	last_x = rand(size(x)) # reset last_x
	modi = Models.Seismic_zeros(pa.modi.mgrid);
	# replace x and modi with initial model
	Seismic_x!(pa.modm, modi, x, pa, 1)

	# bounds
	lower_x = similar(x); upper_x = similar(x);
	Seismic_xbound!(modi, lower_x, upper_x, pa)

	if(pa.attrib_inv == :cls)
		df = OnceDifferentiable(x -> func_xfwi(x, last_x,  buffer1, buffer2, buffer3, pa),
			  (x, storage) -> grad_xfwi!(x, storage, last_x, buffer1, buffer2, buffer3, pa))
		"""
		Unbounded LBFGS inversion, only for testing
		"""
		#res = optimize(df, x, 
	        #		       LBFGS(),
	 	#	     Optim.Options(g_tol = 1e-12,
		# 			iterations = 10, store_trace = true,
		# 			extended_trace=true, show_trace = true))

		"""
		Bounded LBFGS inversion
		"""
		res = optimize(df, x, 
			       lower_x,
			       upper_x,
			       Fminbox(); optimizer=LBFGS, iterations=1, # barrier function iterations
			     optimizer_o=Optim.Options(g_tol = 1e-12,
		 			iterations = 2, store_trace = true, time_limit=time_limit,
					extended_trace=extended_trace, show_trace = true))
		pa.verbose ? println(res) : nothing
		# convert gradient vector to model
		if(extended_trace)
			gmodi = [Models.Seismic_zeros(pa.modi.mgrid) for itr=1:Optim.iterations(res)]
			gmodm = [Models.Seismic_zeros(pa.modm.mgrid) for itr=1:Optim.iterations(res)] 
		end
		modi = [Models.Seismic_zeros(pa.modi.mgrid) for itr=1:Optim.iterations(res)]
		modm = [Models.Seismic_zeros(pa.modm.mgrid) for itr=1:Optim.iterations(res)]
		for itr=1:Optim.iterations(res)
			# update modm and modi
			Seismic_x!(modm[itr], modi[itr], Optim.x_trace(res)[itr], pa, -1)
			# update gmodm and gmodi
			Seismic_gx!(gmodm[itr],modm[itr],gmodi[itr],modi[itr],Optim.trace(res)[itr].metadata["g(x)"],pa,-1)
		end
		# update pa.model
		pa.modm.χvp[:,:] = modm[Optim.iterations(res)].χvp
		pa.modm.χρ[:,:] = modm[Optim.iterations(res)].χρ
		pa.modi.χvp[:,:] = modi[Optim.iterations(res)].χvp
		pa.modi.χρ[:,:] = modi[Optim.iterations(res)].χρ
	
		# update dcal
		pa.dcal.d = buffer1[1].d

		if(extended_trace)
			return modm, modi, gmodm, gmodi, res
		else
			return modm, modi, res
		end
	elseif(pa.attrib_inv == :migr)
		df = OnceDifferentiable(x -> func_xfwi(x, last_x, buffer1, buffer2, buffer3, pa),
			  (x, storage) -> grad_xfwi!(x, storage, last_x, buffer1, buffer2, buffer3, pa))
		res = optimize(df, x, ConjugateGradient(), Optim.Options(iterations = 0,
			      extended_trace=true, store_trace = true,
		       show_trace = false))

		# convert gradient vector to model
		gmodi = Models.Seismic_zeros(pa.modi.mgrid); 
		gmodm = Models.Seismic_zeros(pa.modm.mgrid); 
		Seismic_gx!(gmodm,pa.modm,gmodi,modi,Optim.trace(res)[1].metadata["g(x)"],pa,-1)
		println("maximum value of g(x):\t",  maximum(Optim.trace(res)[1].metadata["g(x)"]))

		return gmodi
	elseif(pa.attrib_inv == :migr_finite_difference)

		gx = similar(x);
		gx = finite_difference!(x -> func_xfwi(x, last_x, buffer1, buffer2, buffer3, pa), x, gx, :central)

		# convert gradient vector to model
		gmodi = Models.Seismic_zeros(pa.modi.mgrid); 
		gmodm = Models.Seismic_zeros(pa.modm.mgrid); 
		Seismic_gx!(gmodm,pa.modm,gmodi,modi,gx,pa,-1)
		println("maximum value of g(x):\t",  maximum(gx))

		return gmodi

	else
		error("invalid pa.attrib_inv")
	end
end

"allocate buffer"
function update_buffer_zeros(pa::Param)
	nx = pa.modm.mgrid.nx+2*pa.modm.mgrid.npml
	nz = pa.modm.mgrid.nz+2*pa.modm.mgrid.npml

	# allocate buffer
	return ([Data.TD_zeros(1, pa.dcal.tgrid, pa.acqgeom)],
	 zeros(Grid.M2D_boundary(pa.modm.mgrid, 3, :outer, onlycount=true), pa.dcal.tgrid.nx,pa.acqgeom.nss),
		   zeros(nz,nx,3,pa.acqgeom.nss)) 
end


"""
Update the buffer, which consists of the modelled data
and boundary values for adjoint calculation.

# Arguments

* `x::Vector{Float64}` : inversion variable
* `last_x::Vector{Float64}` : buffer is only updated when x!=last_x, and modified such that last_x=x
* `pa::Param` : parameters that are constant during the inversion 
* `modm::Models.Seismic` : 
"""
function update_buffer!(
			buffer1::Array{Data.TD,1},
			buffer2::Array{Float64,3},
			buffer3::Array{Float64,4},
			x::Vector{Float64},
			last_x::Vector{Float64}, 
			pa::Param,
			modm::Models.Seismic=Models.Seismic_zeros(pa.modm.mgrid)
			   )
	if(x!=last_x)
		pa.verbose ? println("updating buffer") : nothing
		copy!(last_x, x)

		if(Models.Seismic_iszero(modm))
			# preallocate model and x --> mod;
			model = Models.Seismic_zeros(pa.modm.mgrid); 
			modi = Models.Seismic_zeros(pa.modi.mgrid); 
			Seismic_x!(model, modi, x, pa, -1)		
		else
			model = deepcopy(modm)
		end

		# generate modelled data, border values, etc.
		if(pa.attrib_mod == :fdtd)
			Fdtd.mod!(model=model, acqgeom=[pa.acqgeom], 
			    acqsrc=[pa.acqsrc],src_flags=[2], 
			    TDout=buffer1,
			    backprop_flag=1,
			    boundary=buffer2, initp=buffer3,
			    tgridmod=pa.dcal.tgrid) 
		elseif(pa.attrib_mod == :fdtd_born)
			buffer[1], buffer[2], buffer[3] = Fdtd.mod(npropwav=2, model=pa.modm0, model_pert=model, 
			      acqgeom=[pa.acqgeom,pa.acqgeom], acqsrc=[pa.acqsrc,pa.acqsrc], 
			      src_flags=[2, 0], recv_flags = [0, 1], 
			  tgridmod=pa.dcal.tgrid, verbose=false, boundary_save_flag=true, born_flag=true);
		elseif(pa.attrib_mod == :fdtd_hborn)
			# storing boundary values 
			hbuffer = Fdtd.mod(model=model, acqgeom=[pa.acqgeom], 
						  acqsrc=[pa.acqsrc], 
						  tgridmod=pa.dcal.tgrid, 
						  boundary_save_flag=true);


			buffer[1], buffer[2], buffer[3] = Fdtd.mod(boundary_in=hbuffer[2],
			      npropwav=2,model=pa.modm0, model_pert=model, 
			       acqgeom=[pa.acqgeom,pa.acqgeom], 
			       acqsrc=[acqsrcsink,pa.acqsrc], src_flags=[3.0, 0.0], 
			       recv_flags=[0, 1],
			   tgridmod=pa.dcal.tgrid, verbose=false, boundary_save_flag=true, born_flag=true);
		else
			error("invalid pa.attrib_mod")
		end

		return buffer1[1]
	end
end

"""
Return functional and gradient of the CLS objective 
"""
function func_xfwi(x::Vector{Float64},  last_x::Vector{Float64},
	      buffer1::Array{Data.TD,1},
	      buffer2::Array{Float64,3},
	      buffer3::Array{Float64,4},
	      pa::Param)
	pa.verbose ? println("computing functional...") : nothing
	update_buffer!(buffer1, buffer2, buffer3, x, last_x,  pa)

	# apply coupling to the modeled data
	dcalc = deepcopy(buffer1[1]);
	Data.TDcoup!(dcalc, buffer1[1], pa.w, :s)

	# compute misfit and δdcal
	f, δdcal = Misfits.TD(dcalc, pa.dobs)
	return f
end

function grad_xfwi!(x::Vector{Float64}, storage::Vector{Float64}, 
	       last_x::Vector{Float64}, 
	       buffer1::Array{Data.TD,1},
	       buffer2::Array{Float64,3},
	       buffer3::Array{Float64,4},
	       pa::Param)

	pa.verbose ? println("computing gradient...") : nothing
	update_buffer!(buffer1, buffer2, buffer3, x, last_x, pa)

	# apply coupling to the modeled data
	δdcal = deepcopy(buffer1[1]);
	Data.TDcoup!(δdcal, buffer1[1], pa.w, :s)

	# compute gradient w.r.t. coupled calculated data
	f, dcalc = Misfits.TD(δdcal, pa.dobs)

	# compute gradient w.r.t. calculated data
	Data.TDcoup!(dcalc, δdcal, pa.w, :r)


	# only compute if gradient necessary (saves time)
	# adjoint sources
	adjacq = AdjGeom(δdcal.acqgeom);
	adjsrc = AdjSrc(δdcal)

	# preallocate model
	modm = Models.Seismic_zeros(pa.modm.mgrid)
	modi = Models.Seismic_zeros(pa.modi.mgrid)
	# x to model
	Seismic_x!(modm, modi, x, pa, -1)		

	gmodel=Models.Seismic_zeros(pa.modm.mgrid)
	
	if(pa.attrib_mod == :fdtd)
		# adjoint simulation
		Fdtd.mod!(npropwav=2, model=modm,
	        	acqgeom=[pa.acqgeom, adjacq], acqsrc=[pa.acqsrc, adjsrc], src_flags=[3, 2],
			backprop_flag=-1, boundary=buffer2, 
			initp=buffer3, gmodel=gmodel, 
			tgridmod=pa.dobs.tgrid, gmodel_flag=true, verbose=false)

	elseif(pa.attrib_mod == :fdtd_born)
		
		# adjoint simulation
		adj = Fdtd.mod(npropwav=2, model=pa.modm0,  
		     acqgeom=[pa.acqgeom,adjacq], acqsrc=[acqsrcsink, adjsrc], src_flags=[2.0, -2.0], 
		     tgridmod=pa.dobs.tgrid, gmodel_flag=true, boundary_in=buffer[2], verbose=false)

	elseif(pa.attrib_mod == :fdtd_hborn)
		# adjoint simulation
		adj = Fdtd.mod(npropwav=2, model=pa.modm0, 
		     acqgeom=[pa.acqgeom,adjacq], acqsrc=[pa.acqsrc, adjsrc], src_flags=[-2.0, -2.0], 
		     tgridmod=pa.dobs.tgrid, grad_out_flag=true, verbose=false)
	else
		error("invalid pa.attrib_mod")
	end

	gmodi = Models.Seismic_zeros(pa.modi.mgrid);
	# convert gradient model to gradient vector
	Seismic_gx!(gmodel,modm,gmodi,modi,storage,pa,1) 

	return storage
end # grad_xfwi!

"""
Update pa.w
"""
function wfwi!(pa::Param)

	# convert initial model to the inversion variable
	x = zeros(wfwi_ninv(pa));
	pa.verbose ? println("wfwi: number of inversion variables:\t", length(x)) : nothing
	last_x = rand(size(x)) # reset last_x

	# initial w to x
	Coupling_x!(pa.w, x, pa, 1)

	(Data.TD_iszero(pa.dcal)) ? error("dcal zero wfwi") : buffer1=deepcopy(pa.dcal)

#	if(pa.attrib_inv == :cls)
	df = OnceDifferentiable(x -> func_Coupling(x, last_x, buffer1, pa),
		  (x, storage) -> grad_Coupling!(x, storage, last_x, buffer1, pa))
	"""
	Unbounded LBFGS inversion, only for testing
	"""
	res = optimize(df, x, 
		       LBFGS(),
		       Optim.Options(g_tol = 1e-12,
		       iterations = 10, store_trace = true,
		       extended_trace=true, show_trace = true))
	"testing gradient using auto-differentiation"
#	res = optimize(x -> func_Coupling(x, last_x, buffer1, pa), x, 
#		       LBFGS(),
#		       Optim.Options(g_tol = 1e-12,
#	 	       iterations = 10, store_trace = true,
#		       extended_trace=true, show_trace = true))

	pa.verbose ? println(res) : nothing
	# update w in pa
	Coupling_x!(pa.w, Optim.x_trace(res)[Optim.iterations(res)], pa, -1)

end # wfwi

"""
Convert `Seismic` model to x and vice versa

* `modm::Models.Seismic` : seismic model on modelling grid (input zeros to not use it)
* `modi::Models.Seismic` : seismic model on inversion grid
* `x::Vector{Float64}` : inversion vector
* `pa::Param` : fwi parameters
* `flag::Int64` : 
  * `=1` converts either modm or modi to x
  * `=-1` updates both modm and modi using x
"""
function Seismic_x!(
		    modm::Models.Seismic, 
		    modi::Models.Seismic,
		    x::Vector{Float64}, 
		    pa::Param, 
		    flag::Int64)
	if(flag ==1) # convert modm or modi to x
		all([Models.Seismic_iszero(modm), Models.Seismic_iszero(modi)]) ? 
					error("input modi or modm") : nothing

		# 1. get modi using interpolation and modm
		Models.Seismic_iszero(modm) ? nothing : Models.Seismic_interp_spray!(modm, modi, :interp)

		nznx = pa.modi.mgrid.nz*pa.modi.mgrid.nx;
		# 2. re-parameterization on the inversion grid (include other parameterizations here)
		x[1:nznx] = copy(vec(Models.Seismic_get(modi,pa.parameterization[1])));
		x[nznx+1:2*nznx] = copy(vec(Models.Seismic_get(modi,pa.parameterization[2])));

	elseif(flag == -1) # convert x to mod
		nznx = pa.modi.mgrid.nz*pa.modi.mgrid.nx
		modi.mgrid = deepcopy(pa.modi.mgrid); 
		modi.vp0 = pa.modm.vp0; modi.vs0 = pa.modm.vs0; modi.ρ0 = pa.modm.ρ0;
		# 1. re-parameterization on the inversion grid
		Models.Seismic_reparameterize!(modi, reshape(x[1:nznx],pa.modi.mgrid.nz, pa.modi.mgrid.nx), 
			     reshape(x[nznx+1:2*nznx],pa.modi.mgrid.nz,pa.modi.mgrid.nx), pa.parameterization)

		# 2. interpolation from the inversion grid to the modelling grid
		modm.mgrid = deepcopy(pa.modm.mgrid); 
		Models.Seismic_interp_spray!(modi, modm, :interp)
	else
		error("invalid flag")
	end
end

"""
Return bound vectors for the `Seismic` model, 
depeding on paramaterization
"""
function Seismic_xbound!(modi, lower_x, upper_x, pa)
	nznx = pa.modi.mgrid.nz*pa.modi.mgrid.nx;

	x1=Models.Seismic_get(modi, Symbol((replace("$(pa.parameterization[1])", "χ", "")),0))
	χx1 = Models.χ(x1,x1,1)
	x2=Models.Seismic_get(modi, Symbol((replace("$(pa.parameterization[2])", "χ", "")),0))
	χx2 = Models.χ(x2,x2,1)

	lower_x[1:nznx] = copy(fill(χx1[1], nznx))
	lower_x[nznx+1:2*nznx] = copy(fill(χx2[1], nznx))
	upper_x[1:nznx] = copy(fill(χx1[2], nznx))
	upper_x[nznx+1:2*nznx] = copy(fill(χx2[2], nznx))
end

"""
Convert gradient vector to `Seismic` type and vice versa
This will be different from the previous one, once 
the parameterizations come in

* `gmodm::Models.Seismic` : gradient model on the modelling grid
* `modm::Models.Seismic` : model on the modelling grid
* `gmodi::Models.Seismic` : gradient model on the inversion grid
* `modi::Models.Seismic` : model on the inversion grid
* `gx::Vector{Float64}` : gradient vector
* `pa::Param` :
* `flag::Int64` :
  * `=1` update the vector `gx` using `gmod`
  * `=-1` update gmod 
"""
function Seismic_gx!(gmodm::Models.Seismic,
		     modm:: Models.Seismic,
		     gmodi::Models.Seismic,
		     modi::Models.Seismic,
	             gx::Vector{Float64},
	             pa::Param,
	             flag::Int64 
		      )
	nznx=modi.mgrid.nz*modi.mgrid.nx
	g1=zeros(nznx); g2 =zeros(nznx)
	Models.Seismic_issimilar(gmodm, modm) ? nothing : error("gmodm, modm dimensions")
	Models.Seismic_issimilar(gmodi, modi) ? nothing : error("gmodi, modi dimensions")
	length(gx) == xfwi_ninv(pa) ? nothing : error("gx dimensions")
	if(flag ==1) # convert gmod to gx

		# either gmodi or gmodm should be nonzero
		all([Models.Seismic_iszero(gmodm), Models.Seismic_iszero(gmodi)]) ? 
					error("input gmodi or gmodm") : nothing

		# 1. update gmodi by spraying gmodm
		Models.Seismic_iszero(gmodm) ? nothing : Models.Seismic_interp_spray!(gmodm, gmodi, :spray)

		# 2. chain rule on gmodi 
		Models.Seismic_chainrule!(gmodi, modi, g1, g2, pa.parameterization,-1)

		# 3. copy 
		gx[1:nznx] = copy(g1);
		gx[nznx+1:2*nznx] = copy(g2)
	elseif(flag == -1) # convert gx to mod
		# 1. chain rule depending on re-parameterization
		#    and 
		g1 = gx[1:nznx]; g2 = gx[nznx+1:2*nznx];
		Models.Seismic_chainrule!(gmodi, modi, g1, g2, pa.parameterization, 1)

		# 2. get gradient after interpolation (just for visualization, not exact)
		Models.Seismic_interp_spray!(gmodi, gmodm, :interp)
	else
		error("invalid flag")
	end
end

function AdjSrc(δdat::Data.TD,
	       )
	adjacq = AdjGeom(δdat.acqgeom);
	return Acquisition.Src(adjacq.nss,adjacq.ns,δdat.nfield,
			[(flipdim(δdat.d[i,j],1)) for i in 1:adjacq.nss, j in 1:δdat.nfield], δdat.tgrid);
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
Return functional and gradient of the CLS objective 
"""
function func_Coupling(x::Vector{Float64}, last_x::Vector{Float64}, buffer1::Data.TD, pa::Param)
	# allocate for w
	w = Coupling.TD_delta(pa.w.tgridssf, pa.w.tgridrf, pa.w.nfield, pa.w.acqgeom)

	# x to w 
	Coupling_x!(w, x, pa, -1)

	if(x!=last_x)
		pa.verbose ? println("updating buffer") : nothing
		copy!(last_x, x)
		Data.TDcoup!(buffer1, pa.dcal, w, :s)
	end

	# compute misfit and δdcal
	f, δdcal = Misfits.TD(buffer1, pa.dobs)
	return f
end

"There is a bug in gradient computation"
function grad_Coupling!(x::Vector{Float64}, storage::Vector{Float64},
			last_x::Vector{Float64}, 
			buffer1::Data.TD, pa::Param)

	pa.verbose ? println("computing gradient...") : nothing

	# allocate for w
	w = Coupling.TD_delta(pa.w.tgridssf, pa.w.tgridrf, pa.w.nfield, pa.w.acqgeom)

	# x to w 
	Coupling_x!(w, x, pa, -1)

	if(x!=last_x)
		pa.verbose ? println("updating buffer") : nothing
		copy!(last_x, x)
		Data.TDcoup!(buffer1, pa.dcal, w, :s)
	end

	f, δdcal = Misfits.TD(buffer1, pa.dobs)

	# allocate for w
	gw = Coupling.TD_delta(pa.w.tgridssf, pa.w.tgridrf, pa.w.nfield, pa.w.acqgeom)

	# get gradient w.r.t. w
	Data.TDcoup!(δdcal, pa.dcal, gw, :w)

	# convert w to x
	Coupling_x!(gw, storage, pa, 1)

	return storage
end

"""
Convert coupling functions to x and vice versa
"""
function Coupling_x!(w::Coupling.TD, 
		   x::Vector{Float64}, 
		   pa::Param, 
		   flag::Int64)
	if(flag ==1) # convert w to x
		x[1:pa.w.tgridssf.nx] = copy(vec(w.ssf[1,1][:]));
	elseif(flag == -1) # convert x to w
		# update source functions
		for iss=1:pa.w.acqgeom.nss, ifield=1:pa.w.nfield
			w.ssf[iss, ifield][:] = copy(x)
		end
	else
		error("invalid flag")
	end
end


macro forwardrule(x, e)
	x, e = esc(x), esc(e)
	quote
	$e = sqrt(eps(eltype($x))) * max(one(eltype($x)), abs($x))
	end
end

macro centralrule(x, e)
	x, e = esc(x), esc(e)
	quote
	$e = cbrt(eps(eltype($x))) * max(one(eltype($x)), abs($x))
	end
end

function finite_difference{T<:Number}(f::Function,x::T, dtype::Symbol = :central)
	if dtype == :forward
		@forwardrule x epsilon
		xplusdx = x + epsilon
		return (f(xplusdx) - f(x)) / epsilon
	elseif dtype == :central
		@centralrule x epsilon
		xplusdx, xminusdx = x + epsilon, x - epsilon
		return (f(xplusdx) - f(xminusdx)) / (epsilon + epsilon)
	else
		error("dtype must be :forward, :central")
	end
end



function finite_difference!{S <: Number, T <: Number}(f::Function,
	x::Vector{S},
	g::Vector{T},
	dtype::Symbol)
	# What is the dimension of x?
	n = length(x)

	gpar = SharedArray(eltype(g), size(g))
	xpar = SharedArray(eltype(x), size(x)); xpar = x;
	# Iterate over each dimension of the gradient separately.
	# Use xplusdx to store x + dx instead of creating a new vector on each pass.
	# Use xminusdx to store x - dx instead of creating a new vector on each pass.
	if dtype == :forward
		# Establish a baseline value of f(x).
		f_x = f(x)
		for i = 1:n
			@forwardrule x[i] epsilon
			oldx = x[i]
			x[i] = oldx + epsilon
			f_xplusdx = f(x)
			x[i] = oldx
			g[i] = (f_xplusdx - f_x) / epsilon
		end
	elseif dtype == :central
		@sync @parallel for i = 1:n
			@centralrule xpar[i] epsilon
			oldx = xpar[i]
			xpar[i] = oldx + epsilon
			f_xplusdx = f(xpar)
			xpar[i] = oldx - epsilon
			f_xminusdx = f(xpar)
			xpar[i] = oldx
			gpar[i] = (f_xplusdx - f_xminusdx) / (epsilon + epsilon)
		end
	else
		error("dtype must be :forward or :central")
	end

	g[:] = copy(gpar);
	return g
end


end # module
