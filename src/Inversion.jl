__precompile__()

"""
This module defines a type called `Param` that is to be constructed 
before performing 
and kind of inversion.
The following functionality is currently added in this module:
* a simple RTM
* full-waveform inversion
* source and receiver filter inversion
Once a `Param` object is constructed, the following routines update certain fields
of the `Param` after performing the inversion.
* `xfwi!` : updates the initial subsurface models in `Param` using its observed data
* `wfwi!` : updates the source and receiver filters in `Param` using its observed data
"""
module Inversion

import SIT.Interpolation
import SIT.Models
import SIT.Grid
import SIT.Acquisition
import SIT.Data
import SIT.Coupling
import SIT.Misfits
import SIT.Fdtd
using Optim
using DistributedArrays


"""
Inversion Parameters

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
* `verbose` : print level on STOUT
* `attrib` : synthetic or real

TODO: add an extra attribute for coupling functions inversion and modelling
"""
type Param
	"base source wavelet"
	acqsrc::Acquisition.Src
	"adjoint source functions"
	adjsrc::Acquisition.Src
	"acquisition geometry"
	acqgeom::Acquisition.Geom
	"acquisition geometry for adjoint propagation"
	adjacqgeom::Acquisition.Geom
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
	"gradient Seismic model on modeling grid"
	gmodm::Models.Seismic
	"gradient Seismic on inversion grid"
	gmodi::Models.Seismic
	"parameterization"
	parameterization::Vector{Symbol}
	"model preconditioning"
	mprecon::AbstractArray{Float64,2}
	"calculated data after applying w"
	wdcal::Data.TD
	"calculated data"
	dcal::Data.TD
	"observed data"
	dobs::Data.TD
	"data preconditioning"
	dprecon::Data.TD
	"store temporary data"
	dtemp::Vector{Data.TD}
	"coupling functions"
	w::Coupling.TD
	verbose::Bool
	attrib::Symbol
	buffer::Any
end


"""
Constructor for `Param`

# Arguments

* `acqsrc::Acquisition.Src` : source time functions
* `acqgeom::Acquisition.Geom` : acquisition geometry
* `tgrid::Grid.M1D` : modelling time grid
* `attrib_mod::Symbol` : modelling attribute
* `attrib_inv::Symbol` : inversion attribute
* `modm::Models.Seismic` : seismic model on modelling mesh 

# Optional Arguments

* `igrid::Grid.M2D=modm.mgrid` : inversion grid if different from the modelling grid, i.e., `modm.mgrid`
* `mprecon_flag::Bool=false` : flag to use a model preconditioner
* `dobs::Data.TD` : observed data
* `dprecon::Data.TD=Data.TD_ones(1,dobs.tgrid,dobs.acqgeom)` : data preconditioning, defaults to one 
* `tlagssf::Float64=0.0` : maximum lagtime of unknown source filter
* `tlagrf::Float64=0.0` : maximum lagtime of unknown receiver filter
* `acqsrc_obs::Acquisition.Src=acqsrc` : source wavelets to generate *observed data*; can be different from `acqsrc`
* `modm_obs::Models.Seismic=modm` : actual seismic model to generate *observed data*
* `modm0::Models.Seismic=modm` : background seismic model for Born modelling and inversion (still being tested)
* `mod_inv_parameterization::Vector{Symbol}=[:χKI, :χρI]` : subsurface parameterization
* `verbose::Bool=false` : print level on STDOUT during inversion 
* `attrib::Symbol=:synthetic` : an attribute to control class
  * `=:synthetic` synthetic data inversion
  * `=:field` field data inversion
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
	       mprecon_flag::Bool=false,
	       dobs::Data.TD=Data.TD_zeros(1,tgrid,acqgeom),
	       dprecon::Data.TD=Data.TD_ones(1,dobs.tgrid,dobs.acqgeom),
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

        # acqgeom geometry for adjoint propagation
	adjacqgeom = AdjGeom(acqgeom)

	pa = Param(acqsrc, 
	     Acquisition.Src_zeros(adjacqgeom, nfield, tgrid),
	     acqgeom, 
	     adjacqgeom,
	     attrib_mod, attrib_inv, 
	     modm, modm0, modi, 
	     Models.Seismic_zeros(modm.mgrid), 
	     Models.Seismic_zeros(modi.mgrid), 
	     mod_inv_parameterization,
	     zeros(2,2), # dummy, update mprecon later
	     Data.TD_zeros(nfield,tgrid,acqgeom), 
	     Data.TD_zeros(nfield,tgrid,acqgeom), 
	     dobs, dprecon, 
	     [Data.TD_zeros(nfield,tgrid,acqgeom)], 
	     Coupling.TD_delta(tlagssf, tlagrf, tgrid.δx, nfield, acqgeom), verbose, attrib,
	     nothing)

	# buffer allocation
	forward_simulation_allocate_buffer!(pa)

	# intial modelling dcal
	illum = forward_simulation!(pa, modm=pa.modm, mattrib=:modm, illum_out=mprecon_flag) # x is dummy here, since actually using pa.modm

	# build precon based on illum
	build_mprecon!(pa, illum)

	# generate observed data if attrib is synthetic
	if((attrib == :synthetic) & Data.TD_iszero(pa.dobs))
		# change source wavelets for modelling observed data
		forward_simulation!(pa, modm=modm_obs, mattrib=:modm,acqsrc=acqsrc_obs, dattrib=:dobs) # x is dummy here
	elseif(Data.TD_iszero(pa.dobs) & (attrib == :real))
		error("input observed data for real data inversion")
	else
		error("invalid attrib")
	end

	return pa
end

"""
Return the number of inversion variables for FWI corresponding to `Param`.
This number of inversion variables depend on the size of inversion mesh.
"""
function xfwi_ninv(pa::Param)
	return 2*pa.modi.mgrid.nz*pa.modi.mgrid.nx
end

"""
Return the number of inversion variables for source and receiver filter inversion 
corresponding to `Param`.
This number depends on the maximum lagtimes of the filters. 
"""
function wfwi_ninv(pa::Param)
	return  pa.w.tgridssf.nx
end

"""
Full Waveform Inversion using `Optim`.
This method updates `pa.modm` and `pa.dcal`. More details about the
optional parameters can be found in the documentation of the `Optim` 
package.

# Arguments that are modified

* `pa::Param` : inversion parameters

# Optional Arguments

* `extended_trace::Bool=true` : save extended trace
* `time_limit=Float64=2.0*60.` : time limit for inversion (testing)
* `iterations::Int64=5` : maximum number of iterations
* `f_tol::Float64=1e-3` : functional tolerance
* `g_tol::Float64=1e-3` : gradient tolerance
* `x_tol::Float64=1e-3` : model tolerance

# Outputs

* depending on  `pa.attrib_inv`
  * `=:cls` classic least-squares inversion using adjoint state method
  * `=:migr` return gradient at the first iteration, i.e., a migration image
  * `=:migr_finite_difference` same as above but *not* using adjoint state method; time consuming; only for testing, TODO: implement autodiff here
"""
function xfwi!(pa::Param; extended_trace::Bool=false, time_limit=Float64=2.0*60., iterations::Int64=5, f_tol::Float64=1e-3, g_tol::Float64=1e-3, x_tol::Float64=1e-3)

	# convert initial model to the inversion variable
	x = zeros(xfwi_ninv(pa));
	last_x = similar(x)
	pa.verbose ? println("xfwi: number of inversion variables:\t", length(x)) : nothing

	forward_simulation_allocate_buffer!(pa)

	last_x = rand(size(x)) # reset last_x
	# replace x and modi with initial model (pa.modm)
	Seismic_x!(pa.modm, pa.modi, x, pa, 1)

	# bounds
	lower_x = similar(x); upper_x = similar(x);
	Seismic_xbound!(lower_x, upper_x, pa)

	if(pa.attrib_inv == :cls)
		df = OnceDifferentiable(x -> func_xfwi(x, last_x,  pa),
			  (storage, x) -> grad_xfwi!(storage, x, last_x, pa))
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
			       Fminbox(); optimizer=LBFGS, iterations=iterations, # barrier function iterations
				  #linesearch = Optim.LineSearches.morethuente!,
			     optimizer_o=Optim.Options(g_tol = g_tol, f_tol=f_tol, x_tol=x_tol,
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
	
		if(extended_trace)
			return modm, modi, gmodm, gmodi, res
		else
			return modm, modi, res
		end
	elseif(pa.attrib_inv == :migr)
		storage = similar(x)
		grad_xfwi!(storage, x, last_x,  pa)

		# convert gradient vector to model
		Seismic_gx!(pa.gmodm,pa.modm,pa.gmodi,pa.modi,storage, pa,-1)
		println("maximum value of g(x):\t",  maximum(storage))

		return pa.gmodi
	elseif(pa.attrib_inv == :migr_finite_difference)

		gx = similar(x);
		gx = finite_difference!(x -> func_xfwi(x, last_x, pa), x, gx, :central)

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

"""
Allocations necessary for forward simulations
"""
function forward_simulation_allocate_buffer!(pa::Param)
	nx = pa.modm.mgrid.nx
	nz = pa.modm.mgrid.nz
	npml = pa.modm.mgrid.npml
	pa.buffer = (Fdtd.initialize_boundary(nx, nz, npml, pa.acqgeom.nss, pa.dcal.tgrid.nx))
end


"""
Perform a forward simulation.
This simulation is common for both functional and gradient calculation.
During the computation of the gradient, we need an adjoint simulation.
Update the buffer, which consists of the modelled data
and boundary values for adjoint calculation.

# Arguments

* `x::Vector{Float64}` : inversion variable
* `last_x::Vector{Float64}` : buffer is only updated when x!=last_x, and modified such that last_x=x
* `pa::Param` : parameters that are constant during the inversion 
* `modm::Models.Seismic` : 
"""
function forward_simulation!(
			pa::Param,
			x::Vector{Float64}=zeros(xfwi_ninv(pa)),
			last_x::Vector{Float64}=ones(xfwi_ninv(pa));
			modm::Models.Seismic=Models.Seismic_zeros(pa.modm.mgrid),
			acqsrc::Acquisition.Src=pa.acqsrc,
			mattrib::Symbol=:x,
			dattrib::Symbol=:dcal,
			illum_out::Bool=false
			   )
	if(x!=last_x)
		pa.verbose && println("updating buffer")
		copy!(last_x, x)

		if(mattrib == :x)
			Seismic_x!(pa.modm, pa.modi, x, pa, -1)		
			model = pa.modm
		elseif(mattrib == :modm)
			(Models.Seismic_iszero(modm)) && error("need modm")
			model = modm
		else
			error("invalid mattrib")
		end

		illum=zeros(pa.modm.mgrid.nz, pa.modm.mgrid.nx)
		# generate modelled data, border values, etc.
		if(pa.attrib_mod == :fdtd)
			Fdtd.mod!(model=model, acqgeom=[pa.acqgeom], 
			    acqsrc=[pa.acqsrc],src_flags=[2], 
			    TDout=pa.dtemp,illum=illum, illum_flag=illum_out,
			    backprop_flag=1,
			    boundary=pa.buffer, 
			    tgridmod=pa.dcal.tgrid) 
		elseif(pa.attrib_mod == :fdtd_born)
			Fdtd.mod(born_flag=true,npropwav=2, model=pa.modm0, model_pert=model, 
	    		    acqgeom=[pa.acqgeom,pa.acqgeom], acqsrc=[pa.acqsrc,pa.acqsrc], 
			    TDout=buffer1, illum=illum, illum_flag=mprecon_flag,
			    backprop_flag=1,
			    boundary=buffer2,
			    src_flags=[2, 0], recv_flags = [0, 1], 
			    tgridmod=pa.dcal.tgrid);
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
		setfield!(pa, dattrib, deepcopy(pa.dtemp[1]))
		return illum
	end
end

"build a model precon"
function build_mprecon!(pa,illum::Array{Float64})

	illumi = zeros(pa.modi.mgrid.nz, pa.modi.mgrid.nx)
	if(maximum(abs, illum) != 0.0)
		Interpolation.interp_spray!(pa.modm.mgrid.x, pa.modm.mgrid.z, illum,
		      pa.modi.mgrid.x, pa.modi.mgrid.z, illumi, :spray, :B2)

		# use log to fix scaling of illum 
		# (note that illum is not the exact diagonal of the Hessian matrix, so I am doing tricks here) 
		# TODO: insert a factor here to control the amount of preconditioning
		illumi = 1.0+log.(illumi./minimum(illumi))

		illumi = vcat(vec(illumi), vec(illumi))

		length(illumi) == xfwi_ninv(pa) ? nothing : error("something went wrong with creating mprecon")
		pa.mprecon = spdiagm(illumi,(0),xfwi_ninv(pa), xfwi_ninv(pa))
	else
		pa.mprecon = spdiagm(ones(xfwi_ninv(pa)),(0),xfwi_ninv(pa),xfwi_ninv(pa))
	end
end


"""
Return functional and gradient of the CLS objective 
"""
function func_xfwi(x::Vector{Float64},  last_x::Vector{Float64}, pa::Param)

	pa.verbose && println("computing functional...")

	forward_simulation!(pa, x, last_x)

	# apply coupling to the modeled data
	Data.TDcoup!(pa.wdcal, pa.dcal, pa.w, :s)

	# compute misfit and δdcal
	f = Misfits.TD!(pa.wdcal, pa.dobs, pa.dprecon, :func)
	return f
end

function grad_xfwi!(storage::Vector{Float64}, x::Vector{Float64}, last_x::Vector{Float64}, pa::Param)

	pa.verbose && println("computing gradient...")

	forward_simulation!(pa, x, last_x)

	# apply coupling to the modeled data
	Data.TDcoup!(pa.wdcal, pa.dcal, pa.w, :s)

	# compute gradient w.r.t. coupled calculated data
	f = Misfits.TD!(pa.wdcal, pa.dobs, pa.dprecon, :grad)

	# compute gradient w.r.t. calculated data
	Data.TDcoup!(pa.wdcal, pa.dcal, pa.w, :r)

	# adjoint sources
	update_adjsrc!(pa.dcal, pa)

	# x to model
	Seismic_x!(pa.modm, pa.modi, x, pa, -1)		

	
	if(pa.attrib_mod == :fdtd)
		# adjoint simulation
		Fdtd.mod!(npropwav=2, model=pa.modm,
	        	acqgeom=[pa.acqgeom, pa.adjacqgeom], acqsrc=[pa.acqsrc, pa.adjsrc], src_flags=[3, 2],
			backprop_flag=-1, boundary=pa.buffer, 
			gmodel=pa.gmodm, 
			tgridmod=pa.dcal.tgrid, gmodel_flag=true, verbose=pa.verbose)

	elseif(pa.attrib_mod == :fdtd_born)
		
		# adjoint simulation
		adj = Fdtd.mod(npropwav=2, model=pa.modm0,  
		     acqgeom=[pa.acqgeom, pa.adjacqgeom], acqsrc=[acqsrcsink, adjsrc], src_flags=[2.0, -2.0], 
		     tgridmod=pa.dobs.tgrid, gmodel_flag=true, boundary_in=buffer[2], verbose=false)

	elseif(pa.attrib_mod == :fdtd_hborn)
		# adjoint simulation
		adj = Fdtd.mod(npropwav=2, model=pa.modm0, 
		     acqgeom=[pa.acqgeom, pa.adjacqgeom], acqsrc=[pa.acqsrc, adjsrc], src_flags=[-2.0, -2.0], 
		     tgridmod=pa.dobs.tgrid, grad_out_flag=true, verbose=false)
	else
		error("invalid pa.attrib_mod")
	end

	Seismic_gx!(pa.gmodm,pa.modm,pa.gmodi,pa.modi,storage,pa,1) 

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

	(Data.TD_iszero(pa.dcal)) && error("dcal zero wfwi")

	df = OnceDifferentiable(x -> func_Coupling(x, last_x, pa),
		  (storage, x) -> grad_Coupling!(storage, x, last_x, pa))
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

		# apply preconditioner
		x[:] =  copy(pa.mprecon * x)

	elseif(flag == -1) # convert x to mod

		nznx = pa.modi.mgrid.nz*pa.modi.mgrid.nx
		modi.mgrid = deepcopy(pa.modi.mgrid); 
		modi.vp0 = pa.modm.vp0; modi.vs0 = pa.modm.vs0; modi.ρ0 = pa.modm.ρ0;

		# apply preconditioner 
		Px = copy(pa.mprecon \ x)

		# 1. re-parameterization on the inversion grid
		Models.Seismic_reparameterize!(modi, reshape(Px[1:nznx],pa.modi.mgrid.nz, pa.modi.mgrid.nx), 
			     reshape(Px[nznx+1:2*nznx],pa.modi.mgrid.nz,pa.modi.mgrid.nx), pa.parameterization)

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
function Seismic_xbound!(lower_x, upper_x, pa)
	nznx = pa.modi.mgrid.nz*pa.modi.mgrid.nx;

	x1=Models.Seismic_get(pa.modi, Symbol((replace("$(pa.parameterization[1])", "χ", "")),0))
	χx1 = Models.χ(x1,x1,1)
	x2=Models.Seismic_get(pa.modi, Symbol((replace("$(pa.parameterization[2])", "χ", "")),0))
	χx2 = Models.χ(x2,x2,1)

	lower_x[1:nznx] = copy(fill(χx1[1], nznx))
	lower_x[nznx+1:2*nznx] = copy(fill(χx2[1], nznx))
	upper_x[1:nznx] = copy(fill(χx1[2], nznx))
	upper_x[nznx+1:2*nznx] = copy(fill(χx2[2], nznx))

	# apply preconditioner to the bounds as well
	lower_x[:] =  copy(pa.mprecon * lower_x)
	upper_x[:] =  copy(pa.mprecon * upper_x)

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

		# modify gradient as mprecon is applied to the model vector 
		gx[:] = copy(pa.mprecon \ gx)

	elseif(flag == -1) # convert gx to mod
		# this flag is only used while visuavalizing the gradient vector,
		# the operations are not exact adjoint of flag==1

		# apply preconditioner, commented as the gx to mod is inexact

		#gx =  copy(pa.mprecon * gx)

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

"""
Convert the data `TD` to `Src` after time reversal.
"""
function update_adjsrc!(δdat::Data.TD, pa::Param)
	for i in 1:pa.adjacqgeom.nss, j in 1:δdat.nfield
		pa.adjsrc.wav[i,j] = (flipdim(δdat.d[i,j],1))
	end
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

	geomout.rx = geomin.sx; geomout.rz = geomin.sz;
	geomout.nr = geomin.ns; 

	return geomout
end


"""
Return functional and gradient of the CLS objective 
"""
function func_Coupling(x::Vector{Float64}, last_x::Vector{Float64}, pa::Param)
	# allocate for w
	w = Coupling.TD_delta(pa.w.tgridssf, pa.w.tgridrf, pa.w.nfield, pa.w.acqgeom)

	# x to w 
	Coupling_x!(w, x, pa, -1)

	if(x!=last_x)
		pa.verbose ? println("updating buffer") : nothing
		copy!(last_x, x)
		Data.TDcoup!(pa.wdcal, pa.dcal, w, :s)
	end

	# compute misfit and δdcal
	f = Misfits.TD(pa.wdcal, pa.dobs, pa.dprecon)
	return f
end

"There is a bug in gradient computation"
function grad_Coupling!(x::Vector{Float64}, storage::Vector{Float64},
			last_x::Vector{Float64}, 
			pa::Param)

	pa.verbose ? println("computing gradient...") : nothing

	# allocate for w
	w = Coupling.TD_delta(pa.w.tgridssf, pa.w.tgridrf, pa.w.nfield, pa.w.acqgeom)

	# x to w 
	Coupling_x!(w, x, pa, -1)

	if(x!=last_x)
		pa.verbose ? println("updating buffer") : nothing
		copy!(last_x, x)
		Data.TDcoup!(pa.wdcal, pa.dcal, w, :s)
	end

	f = Misfits.TD(pa.wdcal, pa.dobs, pa.dprecon, :grad)

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
