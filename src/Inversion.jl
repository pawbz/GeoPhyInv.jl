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

# Fields

* `mgrid::Grid.M2D` : modelling grid
* `igrid::Grid.M2D` : inversion grid
* `acqsrc::Acquisition.Src` : base source wavelet for modelling data
* `acqgeom::Acquisition.Geom` : acquisition geometry
* `tgrid::Grid.M1D` : 
* `attrib_mod::Symbol`
* `model0` : background velocity model (only used during Born modeling and inversion)
* `parameterization` : a vector of Symbols specifying parameterization of the inversion vector
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
	parameterization::Vector{Symbol}
end


function Param(mod::Models.Seismic,acqsrc::Acquisition.Src,acqgeom::Acquisition.Geom,tgrid::Grid.M1D,
	      attrib_mod::Symbol, attrib_inv::Symbol, model0::Models.Seismic)
	return Param(mod.vp0, mod.vs0, mod.ρ0, mod.mgrid, mod.mgrid, acqsrc, acqgeom,tgrid,
	      attrib_mod, attrib_inv, model0, [:χKI, :χρI] )
end

"Return the number of inversion variables `Param`"
function fwi_ninv(pa::Param)
	return 2*pa.igrid.nz*pa.igrid.nx
end

"""
FWI using Optim

# Arguments
* `dobs::Data.TD=Data.TD_zeros(1,pa.tgrid,pa.acqgeom)` : input observed data
"""
function xfwi(
	      pa::Param;
	      model::Models.Seismic=nothing,
	      model_init::Models.Seismic=nothing,
	      dobs::Data.TD=Data.TD_zeros(1,pa.tgrid,pa.acqgeom)
	      )

	# convert initial model to the inversion variable
	x = zeros(fwi_ninv(pa));
	last_x = similar(x)
	println("xfwi: number of inversion variables:\t", length(x))

	# generate observed data using the same modelling scheme
	if(Data.TD_iszero(dobs))
		dobs = update_buffer!(x, last_x, pa, model) # x is dummy here
	end

	last_x = rand(size(last_x)) # reset last_x
	modi = Models.Seismic_zeros(pa.igrid);
	# replace x and modi with initial model
	Seismic_x!(model_init, modi, x, pa, 1)

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
		df = OnceDifferentiable(x -> func(x, last_x, dobs, pa),
			  (x, storage) -> grad!(x, storage, last_x, dobs, pa))
		res = optimize(df, x, ConjugateGradient(), Optim.Options(iterations = 0,
			      extended_trace=true, store_trace = true,
		       show_trace = false))

		# convert gradient vector to model
		gmodi = Models.Seismic_zeros(pa.igrid); 
		gmodm = Models.Seismic_zeros(pa.mgrid); 
		Seismic_gx!(gmodm,model_init,gmodi,modi,Optim.trace(res)[1].metadata["g(x)"],pa,-1)
		println("maximfgjfnjf",  maximum(Optim.trace(res)[1].metadata["g(x)"]))

		return gmodi, dobs
	elseif(pa.attrib_inv == :migr_finite_difference)

		gx = similar(x);
		gx = finite_difference!(x -> func(x, last_x, dobs, pa), x, gx, :central)

		# convert gradient vector to model
		gmodi = Models.Seismic_zeros(pa.igrid); 
		gmodm = Models.Seismic_zeros(pa.mgrid); 
		Seismic_gx!(gmodm,model_init,gmodi,modi,gx,pa,-1)

		println("maximfgjigggggggfnjf",  maximum(gx))

		return gmodi, dobs

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
* `modm::Models.Seismic` : 
"""
function update_buffer!(
			x::Vector{Float64},
			last_x::Vector{Float64}, 
			pa::Param,
			modm::Models.Seismic=Models.Seismic_zeros(pa.mgrid)
			   )
	if(x!=last_x)
		println("updating buffer")
		copy!(last_x, x)

		if(Models.Seismic_iszero(modm))
			# preallocate model and x --> mod;
			model = Models.Seismic_zeros(pa.mgrid); 
			modi = Models.Seismic_zeros(pa.igrid); 
			Seismic_x!(model, modi, x, pa, -1)		
		else
			model = deepcopy(modm)
		end

		# generate modelled data, border values, etc.
		if(pa.attrib_mod == :fdtd)
			global buffer = Fdtd.mod(model=model, acqgeom=[pa.acqgeom], 
						  acqsrc=[pa.acqsrc], 
						  tgridmod=pa.tgrid, 
						  boundary_save_flag=true);
		elseif(pa.attrib_mod == :fdtd_born)
			global buffer = Fdtd.mod(npropwav=2, model=pa.model0, model_pert=model, 
			      acqgeom=[pa.acqgeom,pa.acqgeom], acqsrc=[pa.acqsrc,pa.acqsrc], 
			      src_flags=[-2.0, 0.0], recv_flags = [0.0, 2.0], 
			  tgridmod=pa.tgrid, verbose=false, boundary_save_flag=true, born_flag=true);
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
			       acqsrc=[acqsrcsink,pa.acqsrc], src_flags=[-2.0, 0.0], 
			       recv_flags=[0.0, 2.0],
			   tgridmod=pa.tgrid, verbose=false, boundary_save_flag=true, born_flag=true);
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

	update_buffer!(x, last_x, pa)

	# compute misfit and δdcal
	f, δdcal = Misfits.TD(buffer[1][1], dobs)

	# only compute if gradient necessary (saves time)
	# adjoint sources
	adjacq = AdjGeom(δdcal.acqgeom);
	adjsrc = AdjSrc(δdcal)

	# preallocate model
	modm = Models.Seismic_zeros(pa.mgrid);
	modi = Models.Seismic_zeros(pa.igrid);
	# x to model
	Seismic_x!(modm, modi, x, pa, -1)		

	if(pa.attrib_mod == :fdtd)
		# source sinks
		acqsrcsink = deepcopy(pa.acqsrc); 
		acqsrcsink.wav = [-1.0.*flipdim(pa.acqsrc.wav[i,j],1) for 
			   i in 1:pa.acqsrc.nss, j in 1:pa.acqsrc.nfield];

		
		# adjoint simulation
		adj = Fdtd.mod(npropwav=2, model=modm,
		     acqgeom=[pa.acqgeom,adjacq], acqsrc=[acqsrcsink, adjsrc], src_flags=[-2.0, -2.0], 
		     tgridmod=dobs.tgrid, grad_out_flag=true, boundary_in=buffer[2], verbose=false)

	elseif(pa.attrib_mod == :fdtd_born)
		# source sinks
		acqsrcsink = deepcopy(pa.acqsrc); 
		acqsrcsink.wav = [-1.0.*flipdim(pa.acqsrc.wav[i,j],1) for 
			   i in 1:pa.acqsrc.nss, j in 1:pa.acqsrc.nfield];

		
		# adjoint simulation
		adj = Fdtd.mod(npropwav=2, model=pa.model0,  
		     acqgeom=[pa.acqgeom,adjacq], acqsrc=[acqsrcsink, adjsrc], src_flags=[-2.0, -2.0], 
		     tgridmod=dobs.tgrid, grad_out_flag=true, boundary_in=buffer[2], verbose=false)

	elseif(pa.attrib_mod == :fdtd_hborn)
		# adjoint simulation
		adj = Fdtd.mod(npropwav=2, model=pa.model0, 
		     acqgeom=[pa.acqgeom,adjacq], acqsrc=[pa.acqsrc, adjsrc], src_flags=[-2.0, -2.0], 
		     tgridmod=dobs.tgrid, grad_out_flag=true, verbose=false)
	else
		error("invalid pa.attrib_mod")
	end

	gmodi = Models.Seismic_zeros(pa.igrid);
	# convert gradient model to gradient vector
	Seismic_gx!(adj[3],modm,gmodi,modi,storage,pa,1) 

	return storage
end # grad!


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

		nznx = pa.igrid.nz*pa.igrid.nx;
		# 2. re-parameterization on the inversion grid (include other parameterizations here)
		x[1:nznx] = copy(vec(Models.Seismic_get(modi,pa.parameterization[1])));
		x[nznx+1:2*nznx] = copy(vec(Models.Seismic_get(modi,pa.parameterization[2])));

	elseif(flag == -1) # convert x to mod
		nznx = pa.igrid.nz*pa.igrid.nx
		modi.mgrid = deepcopy(pa.igrid); 
		modi.vp0 = pa.vp0; modi.vs0 = pa.vs0; modi.ρ0 = pa.ρ0;
		# 1. re-parameterization on the inversion grid
		Models.Seismic_reparameterize!(modi, reshape(x[1:nznx],pa.igrid.nz, pa.igrid.nx), 
			     reshape(x[nznx+1:2*nznx],pa.igrid.nz,pa.igrid.nx), pa.parameterization)

		# 2. interpolation from the inversion grid to the modelling grid
		modm.mgrid = deepcopy(pa.mgrid); 
		Models.Seismic_interp_spray!(modi, modm, :interp)
	else
		error("invalid flag")
	end
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
	length(gx) == fwi_ninv(pa) ? nothing : error("gx dimensions")
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

end # module
