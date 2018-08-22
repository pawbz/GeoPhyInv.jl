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
module FWI

using Interpolation
using Conv
using Grid
using Misfits
using Inversion
using TimerOutputs
using LinearMaps
using LinearAlgebra
using Ipopt

import JuMIT.Models
import JuMIT.Acquisition
import JuMIT.Data
import JuMIT.Coupling
import JuMIT.Fdtd
using Optim, LineSearches
using DistributedArrays
using Calculus
using Random


global const to = TimerOutput(); # create a timer object

# modeling types
struct ModFdtd end

struct ModFdtdBorn end

struct ModFdtdHBorn end

# fields
struct P end
struct Vx end
struct Vz end


"""
FWI Parameters

# Fields

* `mgrid::Grid.M2D` : modelling grid
* `igrid::Grid.M2D` : inversion grid
* `acqsrc::Acquisition.Src` : base source wavelet for modelling data
* `acqgeom::Acquisition.Geom` : acquisition geometry
* `tgrid::Grid.M1D` : 
* `attrib_mod`
* `model_obs` : model used for generating observed data
* `model0` : background velocity model (only used during Born modeling and inversion)
* `parameterization` : a vector of Symbols specifying parameterization of the inversion vector
* `verbose` : print level on STOUT
* `attrib` : synthetic or real

TODO: add an extra attribute for coupling functions inversion and modelling
"""
mutable struct Param{Tmod, Tdatamisfit}
	"model inversion variable"
	mx::Inversion.X{Float64,1}
	mxm::Inversion.X{Float64,1}
	"forward modelling parameters for"
	paf::Tmod
	"base source wavelet"
	acqsrc::Acquisition.Src
	"adjoint source functions"
	adjsrc::Acquisition.Src
	"acquisition geometry"
	acqgeom::Acquisition.Geom
	"acquisition geometry for adjoint propagation"
	adjacqgeom::Acquisition.Geom
	"Seismic model on modeling grid"
	modm::Models.Seismic
	"background Seismic model"
	modm0::Models.Seismic
	"Seismic on inversion grid"
	modi::Models.Seismic
	"initialize on this model"
	mod_initial::Models.Seismic
	"prior model on the inversion grid"
	priori::Models.Seismic
	"prior weights on the inversion grid"
	priorw::Models.Seismic
	"gradient Seismic model on modeling grid"
	gmodm::Models.Seismic
	"gradient Seismic on inversion grid"
	gmodi::Models.Seismic
	"parameterization"
	parameterization::Vector{Symbol}
	"model preconditioning"
	mprecon::AbstractArray{Float64,2}
	"compute data misfit"
	paTD::Tdatamisfit
	paminterp::Interpolation.Kernel_2D{Float64}
	verbose::Bool
	"synthetic or field inversion"
	attrib::Symbol
end

Base.print(pa::Param)=nothing
Base.show(pa::Param)=nothing
Base.display(pa::Param)=nothing

struct LS end # just cls inversion

struct LS_prior # cls inversion including prior 
	α::Vector{Float64} # weights
end

function LS_prior()
	return LS_prior([1.0, 0.5])
end

struct Migr end # returns first gradient

struct Migr_fd end # computes first gradient using FD

include("core.jl")
include("prior.jl")
include("xfwi.jl")
include("xwfwi.jl")

"""
Convert the data `TD` to `Src` after time reversal.
"""
function update_adjsrc!(adjsrc, δdat::Data.TD, adjacqgeom)
	(adjsrc.fields != δdat.fields) && error("dissimilar fields")
	nt=δdat.tgrid.nx
	for i in 1:adjacqgeom.nss
		for (j,field) in enumerate(δdat.fields)
			wav=adjsrc.wav[i,j] 
			dat=δdat.d[i,j]
			for ir in 1:δdat.acqgeom.nr[i]
				for it in 1:nt
					@inbounds wav[it,ir]=value_adjsrc(dat[nt-it+1,ir],eval(field)())
				end
			end
		end
	end
	return nothing
end

value_adjsrc(s, ::P) = s
value_adjsrc(s, ::Vx) = -1.0*s # see adj tests
value_adjsrc(s, ::Vz) = -1.0*s


function generate_adjsrc(fields, tgrid, adjacqgeom)
	adjsrc=Acquisition.Src_zeros(adjacqgeom, fields, tgrid)
	return adjsrc
end


"""
Constructor for `Param`

# Arguments

* `acqsrc::Acquisition.Src` : source time functions
* `acqgeom::Acquisition.Geom` : acquisition geometry
* `tgrid::Grid.M1D` : modelling time grid
* `attrib_mod::Union{ModFdtd, ModFdtdBorn, ModFdtdHBorn}` : modelling attribute
* `modm::Models.Seismic` : seismic model on modelling mesh 

# Optional Arguments
* `tgrid_obs::Grid.M1D` : time grid for observed data

* `igrid::Grid.M2D=modm.mgrid` : inversion grid if different from the modelling grid, i.e., `modm.mgrid`
* `igrid_interp_scheme` : interpolation scheme
  * `=:B1` linear 
  * `=:B2` second order 
* `mprecon_factor::Float64=1` : factor to control model preconditioner, always greater than 1
  * `=1` means the preconditioning is switched off
  * `>1` means the preconditioning is switched on, larger the mprecon_factor, stronger the applied preconditioner becomes
* `dobs::Data.TD` : observed data
* `dprecon::Data.TD=Data.TD_ones(1,dobs.tgrid,dobs.acqgeom)` : data preconditioning, defaults to one 
* `tlagssf_fracs=0.0` : maximum lagtime of unknown source filter
* `tlagrf_fracs=0.0` : maximum lagtime of unknown receiver filter
* `acqsrc_obs::Acquisition.Src=acqsrc` : source wavelets to generate *observed data*; can be different from `acqsrc`
* `modm_obs::Models.Seismic=modm` : actual seismic model to generate *observed data*
* `modm0::Models.Seismic=fill!(modm,0.0)` : background seismic model for Born modelling and inversion
* `parameterization::Vector{Symbol}` : subsurface parameterization
  * `=[:χvp, :χρ, :null]` for velocity and density inversion
  * `=[:null, :χρ, :null]` for density-only inversion
  * `=[:χvp, :null, :null]` for velocity-only inversion
* `verbose::Bool=false` : print level on STDOUT during inversion 
* `attrib::Symbol=:synthetic` : an attribute to control class
  * `=:synthetic` synthetic data inversion
  * `=:field` field data inversion
* nworker : number of workers (input nothing to use all available)
"""
function Param(
	       acqsrc::Acquisition.Src,
	       acqgeom::Acquisition.Geom,
	       tgrid::Grid.M1D,
	       attrib_mod::Union{ModFdtd, ModFdtdBorn, ModFdtdHBorn}, 
	       modm::Models.Seismic;
	       # other optional 
	       tgrid_obs::Grid.M1D=deepcopy(tgrid),
	       rfields=[:P],
	       igrid::Grid.M2D=nothing,
	       igrid_interp_scheme::Symbol=:B2,
	       mprecon_factor::Float64=1.0,
	       dobs::Data.TD=Data.TD_zeros(rfields,tgrid_obs,acqgeom),
	       dprecon=nothing,
	       tlagssf_fracs=0.0,
	       tlagrf_fracs=0.0,
	       acqsrc_obs::Acquisition.Src=deepcopy(acqsrc),
	       modm_obs::Models.Seismic=modm,
	       modm0=nothing,
	       parameterization::Vector{Symbol}=[:χvp, :χρ, :null],
	       verbose::Bool=false,
	       attrib::Symbol=:synthetic,
	       nworker=nothing,
	       )

	# make a copy of initial model
	modm=deepcopy(modm)

	if(igrid===nothing)
		igrid=deepcopy(modm.mgrid)
	end
	# igrid has to truncated because the gradient evaluation 
	# is inaccurate on boundaries
	mg=modm.mgrid
	igrid=Grid.M2D(max(mg.x[3],igrid.x[1]),
		min(igrid.x[end],mg.x[end-2]),
		max(igrid.z[1],mg.z[3]),
		min(igrid.z[end],mg.z[end-2]),
		 		igrid.nx,igrid.nz)


	# create modi according to igrid and interpolation of modm
	modi = Models.Seismic_zeros(igrid);
	# update bounds
	Models.adjust_bounds!(modi, modm)
	Models.interp_spray!(modm, modi, :interp)

	# save the initial model, as modi will be changed during the iterations
	mod_initial=deepcopy(modi)

	# allocate prior inputs
	priori=similar(modi)
	priorw=similar(modi)

	# use default background model modm0
	if(modm0 === nothing)
		modm0=deepcopy(modm) # some dummy
		fill!(modm0, 0.0)
	end
	isequal(modm0, modm_obs) && (@warn "modm0 == modm_obs")

	if(attrib == :synthetic) 
		isequal(modm, modm_obs) && error("initial model same as actual model, zero misfit?")
	end

	# acqgeom geometry for adjoint propagation
	adjacqgeom = AdjGeom(acqgeom)

	# generate adjoint sources
	adjsrc=generate_adjsrc(rfields, tgrid, adjacqgeom)

	# generating forward and adjoint modelling engines
	# to generate modelled data, border values, etc.
	# most of the parameters given to this are dummy
	paf=Fdtd.Param(npw=2, model=modm, born_flag=true,
		acqgeom=[acqgeom, adjacqgeom], acqsrc=[acqsrc, adjsrc], 
		rfields=rfields,
		sflags=[3, 2], rflags=[1, 1],
		backprop_flag=1, 
		tgridmod=tgrid, gmodel_flag=true, verbose=verbose, illum_flag=false, nworker=nworker)


	gmodm=similar(modm)
	gmodi=similar(modi)

	#println("added fake precon")
	#dprecon=deepcopy(dobs)
	#fill!(dprecon, 1.0)
	#Data.taper!(dprecon, 20.)

	# check dprecon
	if(!(dprecon===nothing))
		(!(isapprox(dprecon,dobs))) && error("invalid dprecon used")
	end
	
	# create Parameters for data misfit
	#coup=Coupling.TD_delta(dobs.tgrid, tlagssf_fracs, tlagrf_fracs, rfields, acqgeom)
	# choose the data misfit
	#	(iszero(tlagssf_fracs)) 
 	# tlagssf_fracs==[0.0] | tlagssf_fracs=[0.0])
	# paTD=Data.P_misfit(Data.TD_zeros(rfields,tgrid,acqgeom),dobs,w=dprecon,coup=coup, func_attrib=optims[1]);

	if(iszero(dobs))
		Random.randn!(dobs) # dummy dobs, update later
	end

	paTD=Data.P_misfit(Data.TD_zeros(rfields,tgrid,acqgeom),dobs,w=dprecon);

	paminterp=Interpolation.Kernel([modi.mgrid.x, modi.mgrid.z], [modm.mgrid.x, modm.mgrid.z], igrid_interp_scheme)

	mx=Inversion.X(modi.mgrid.nz*modi.mgrid.nx*count(parameterization.≠:null),2)
	mxm=Inversion.X(modm.mgrid.nz*modm.mgrid.nx*count(parameterization.≠:null),2)


	# put modm as a vector, according to parameterization in mxm.x
	Models.Seismic_get!(mxm.x,modm,parameterization)

	pa = Param(
	     mx, mxm,
	     paf,
	     deepcopy(acqsrc), 
	     Acquisition.Src_zeros(adjacqgeom, rfields, tgrid),
	     deepcopy(acqgeom), 
	     adjacqgeom,
	     deepcopy(modm), deepcopy(modm0), modi, mod_initial,
	     priori, priorw,
	     gmodm,gmodi,
	     parameterization,
	     zeros(2,2), # dummy, update mprecon later
	     paTD,
	     paminterp,
	     verbose, attrib)

	# generate observed data if attrib is synthetic
	if((attrib == :synthetic))
		update_observed_data!(pa, modm_obs, attrib_mod)
	end

	iszero(dobs) && ((attrib == :real) ? error("input observed data for real data inversion") : error("problem generating synthetic observed data"))

	# generate modelled data
	F!(pa, nothing, attrib_mod)

#	build_mprecon!(pa, Array(pa.paf.c.illum_stack), mprecon_factor)
	pa.paf.c.illum_flag=false # switch off illum flag for speed

	# default weights are 1.0
	fill!(pa.mx.w, 1.0)

	# update priori and priorw in pa to defaults
	update_prior!(pa)
	
	return pa
end

"""
update the *synthetic* observed data in `Param`
* allocated memory, don't use in inner loops
"""
function update_observed_data!(pa::Param, modm_obs, attrib_mod=ModFdtd())
	# save modm of pa to put it back later
	modm_copy=deepcopy(pa.modm)

	# change modm in pa to actual model
	copyto!(pa.modm, modm_obs)

	# update models in the forward engine and do modelling
	# ModFdtdBorn: modm0 will be used as a background model for Born modelling and modm will be perturbed model
	# ModFdtd: modm will be *the* model
	F!(pa, nothing, attrib_mod)

	# get observed data
	copyto!(pa.paTD.y, pa.paTD.x)

	# put back model and sources
	copyto!(pa.modm, modm_copy)

	Fdtd.initialize!(pa.paf.c)  # clear everything
end



"""
Return the number of inversion variables for FWI corresponding to `Param`.
This number of inversion variables depend on the size of inversion mesh.
"""
function xfwi_ninv(pa::Param)
	length(pa.mx.x)
end

"""
Return the number of inversion variables for source and receiver filter inversion 
corresponding to `Param`.
This number depends on the maximum lagtimes of the filters. 
"""
function wfwi_ninv(pa::Param)
	return  pa.paTD.coup.tgridssf.nx
end


"build a model precon"
function build_mprecon!(pa,illum::Array{Float64}, mprecon_factor=1.0)
	"illum cannot be zero when building mprecon"
	(any(illum .<= 0.0)) && error("illum cannot be negative or zero")
	(mprecon_factor < 1.0)  && error("invalid mprecon_factor")

	illumi = zeros(pa.modi.mgrid.nz, pa.modi.mgrid.nx)
	Interpolation.interp_spray!(illumi, illum, pa.paminterp, :spray)
	(any(illumi .<= 0.0)) && error("illumi cannot be negative or zero")

	# use log to fix scaling of illum 
	# (note that illum is not the exact diagonal of the Hessian matrix, so I am doing tricks here) 
	# larger the mprecon_factor, stronger the preconditioning
	# division by maximum is safe because illum is non zero
	minillumi = minimum(illumi)
	maxillumi = maximum(illumi)
	illumi -= minillumi; # remove mean 	
	(maxillumi ≠ 0.0) && (illumi ./= maxillumi) # normalize using maxillumi
	#illumi = 1.0 + log.(1. + (mprecon_factor-1) * illumi) # illumi now lies between 1.0 and some value depending on how large illumi is
	illumi = (1. + (mprecon_factor-1) * illumi) # illumi now lies between 1.0 and some value depending on how large illumi is

	illumi = vcat([vec(illumi) for i in 1:count(pa.parameterization.≠:null)]...)

	length(illumi) == xfwi_ninv(pa) ? nothing : error("something went wrong with creating mprecon")
	pa.mprecon = spdiagm(illumi,(0),xfwi_ninv(pa), xfwi_ninv(pa))
end


# convert x to modm by interpolation
# x is on the sparse grid, and is changed by parameterization
function x_to_modm!(pa, x)

	# put sparse x into modi (useless step, just for visualization)
	Models.Seismic_reparameterize!(pa.modi, x, pa.parameterization) 

	# get x on dense grid
	Interpolation.interp_spray!(x, pa.mxm.x, pa.paminterp, :interp, count(pa.parameterization.≠:null))

	# put dense x into modm, according to parameterization
	Models.Seismic_reparameterize!(pa.modm,pa.mxm.x,pa.parameterization) 
end


function δx_to_δmods!(pa, δx)
	# get x on dense grid
	Interpolation.interp_spray!(δx, pa.mxm.x, pa.paminterp, :interp, count(pa.parameterization.≠:null))

	# reparameterize accordingly to get [δKI, δρI]
	Models.pert_reparameterize!(pa.paf.c.δmod, pa.mxm.x, pa.modm0, pa.parameterization)

	# input [δKI, δρI] to the modeling engine
	Fdtd.update_δmods!(pa.paf.c, pa.paf.c.δmod)
end

"""
convert the gradient output from Fdtd to gx
"""
function spray_gradient!(gx,  pa, ::ModFdtd)

	# apply chain rule on gmodm to get gradient w.r.t. dense x
	Models.gradient_chainrule!(pa.mxm.gx, pa.paf.c.gradient, pa.modm, pa.parameterization)

	# spray gradient w.r.t. dense x to the space of sparse x
	Interpolation.interp_spray!(gx, pa.mxm.gx, pa.paminterp, :spray, count(pa.parameterization.≠:null))

	# visualize gmodi here?
	# visualize_gx!(pa.gmodm, pa.modm, pa.gmodi, pa.modi, gx, pa)
end
function spray_gradient!(gx, pa, ::ModFdtdBorn)

	# apply chain rule on gmodm to get gradient w.r.t. dense x (note that background model is used for Born modeling here)
	Models.gradient_chainrule!(pa.mxm.gx, pa.paf.c.gradient, pa.modm, pa.parameterization)

	# spray gradient w.r.t. dense x to the space of sparse x
	Interpolation.interp_spray!(gx, pa.mxm.gx, pa.paminterp, :spray, count(pa.parameterization.≠:null))

	# visualize gmodi here?
	# visualize_gx!(pa.gmodm, pa.modm, pa.gmodi, pa.modi, gx, pa)
end



"""
Return bound vectors for the `Seismic` model, 
depeding on paramaterization
"""
function Seismic_xbound!(lower_x, upper_x, pa)
	nznx = pa.modi.mgrid.nz*pa.modi.mgrid.nx;

	bound1 = similar(lower_x)
	# create a Seismic model with minimum possible values
	modbound = deepcopy(pa.modi);
	modbound.χvp = Models.χ(fill(modbound.vp0[1], size(modbound.χvp)), modbound.ref.vp)
	modbound.χvs = Models.χ(fill(modbound.vs0[1], size(modbound.χvs)), modbound.ref.vs)
	modbound.χρ = Models.χ(fill(modbound.ρ0[1], size(modbound.χρ)), modbound.ref.ρ)

	Models.Seismic_get!(bound1,modbound,pa.parameterization)


	bound2=similar(upper_x)
	# create a Seismic model with maximum possible values
	modbound.χvp = Models.χ(fill(modbound.vp0[2], size(modbound.χvp)), modbound.ref.vp)
	modbound.χvs = Models.χ(fill(modbound.vs0[2], size(modbound.χvs)), modbound.ref.vs)
	modbound.χρ = Models.χ(fill(modbound.ρ0[2], size(modbound.χρ)), modbound.ref.ρ)

	Models.Seismic_get!(bound2,modbound,pa.parameterization)

	# this sorting operation is important because 
	lower_x[:] = min.(bound1, bound2)
	upper_x[:] = max.(bound1, bound2)

	return modbound
end

# gmodi, gmodm <-- gx (just for visualizing gradient)
function visualize_gx!(gmodm, modm, gmodi, modi, gx, pa)
	# chain rule depending on re-parameterization
	Models.Seismic_chainrule!(gmodi, modi, gx, pa.parameterization, 1)

	# get gradient after interpolation (just for visualization, not exact)
	Models.interp_spray!(gmodi, gmodm, :interp, pa=pa.paminterp)
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

#	geomout.rx = geomin.sx; geomout.rz = geomin.sz;
#	geomout.nr = geomin.ns; 

	geomout.rx = geomin.rx; geomout.rz = geomin.rz;
	geomout.nr = geomin.nr; 

	return geomout
end


"""
Return functional and gradient of the LS objective 
"""
function func_grad_Coupling!(storage, x::Vector{Float64},pa::Param)

	pa.verbose ? println("computing gradient...") : nothing

	# x to w 
	Coupling_x!(x, pa, -1)

	if(storage === nothing)
		# compute misfit 
		f = Data.error!(pa.paTD)
		return f
	else
		f = Data.error!(pa.paTD, :dJssf)
		Coupling_gx!(storage, pa.paTD.dJssf)
		return storage
	end
end

"""
Convert coupling functions to x and vice versa
"""
function Coupling_x!(x::Vector{Float64}, 
		   pa::Param, 
		   flag::Int64)
	if(flag ==1) # convert w to x
		iss=1 # take only first source
		for ifield=1:length(pa.paTD.y.fields)
			w=pa.paTD.coup.ssf[iss,ifield]
			for i in eachindex(x)
				x[i]=w[i] # take only first source
			end
		end
	elseif(flag == -1) # convert x to w
		for iss=1:pa.paTD.y.acqgeom.nss, ifield=1:length(pa.paTD.y.fields)
			w=pa.paTD.coup.ssf[iss,ifield]
			for i in eachindex(x)
				w[i]=x[i]
			end
		end
	else
		error("invalid flag")
	end
end

"""
Convert dJssf to gradient vector
"""
function Coupling_gx!(gx, dJssf)
	gx[:]=0.0
	for ifield=1:size(dJssf,2), iss=1:size(dJssf,1)
		dwav=dJssf[iss,ifield]
		for i in eachindex(gx)
			gx[i]+=dwav[i]
		end
	end
end




end # module
