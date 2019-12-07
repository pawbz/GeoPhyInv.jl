"""
This module defines a type called `SeisInvExpt` that is to be constructed 
before performing 
and kind of inversion.
The following functionality is currently added in this module:
* a simple RTM
* full-waveform inversion
* source and receiver filter inversion
Once a `SeisInvExpt` object is constructed, the following routines update certain fields
of the `SeisInvExpt` after performing the inversion.
* `xfwi!` : updates the initial subsurface models in `SeisInvExpt` using its observed data
* `wfwi!` : updates the source and receiver filters in `SeisInvExpt` using its observed data
"""

include("X.jl")

# modeling types
struct ModFdtd end

struct ModFdtdBorn end

struct ModFdtdHBorn end


"""
FWI Parameters

# Fields

* `mgrid::Vector{StepRangeLen}` : modelling grid
* `igrid::Vector{StepRangeLen}` : inversion grid
* `acqsrc::SrcWav` : base source wavelet for modelling data
* `acqgeom::Geom` : acquisition geometry
* `tgrid::StepRangeLen` : 
* `attrib_mod`
* `model_obs` : model used for generating observed data
* `model0` : background velocity model (only used during Born modeling and inversion)
* `parameterization` : a vector of Symbols specifying parameterization of the inversion vector
* `verbose` : print level on STOUT
* `attrib` : synthetic or real

TODO: add an extra attribute for coupling functions inversion and modelling
"""
mutable struct SeisInvExpt{Tmod, Tdatamisfit}
	"model inversion variable"
	mx::X{Float64,1}
	mxm::X{Float64,1}
	"forward modelling parameters for"
	paf::Tmod
	"base source wavelet"
	acqsrc::SrcWav
	"adjoint source functions"
	adjsrc::SrcWav
	"acquisition geometry"
	acqgeom::Geom
	"acquisition geometry for adjoint propagation"
	adjacqgeom::Geom
	"Seismic model on modeling grid"
	modm::Medium
	"background Seismic model"
	modm0::Medium
	"Seismic on inversion grid"
	modi::Medium
	"initialize on this model"
	mod_initial::Medium
	"prior model on the inversion grid"
	priori::Medium
	"gradient Seismic model on modeling grid"
	gmodm::Medium
	"gradient Seismic on inversion grid"
	gmodi::Medium
	"parameterization"
	parameterization::Vector{Symbol}
	"model preconditioning"
	mprecon::AbstractArray{Float64,2}
	"compute data misfit"
	paTD::VNamedD_misfit
	paminterp::Interpolation.Kernel_2D{Float64}
	verbose::Bool
	"synthetic or field inversion"
	attrib::Symbol
end

Base.print(pa::SeisInvExpt)=nothing
Base.show(pa::SeisInvExpt)=nothing
Base.display(pa::SeisInvExpt)=nothing

# add data inverse covariance matrix here later
struct LS end # just cls inversion

struct LS_prior # cls inversion including prior 
	# add data inverse covariance matrix here later, instead of just Float64
	pdgls::Float64
	# pdgls
	pmgls::Misfits.P_gls{Float64}
end
LS_prior(α1::Float64, Q::LinearMap)=LS_prior(α1, Misfits.P_gls(Q))

"""
this method constructs prior term with Q=α*I
* `ninv` : number of inversion variables, use `xfwi_ninv` 
* `α` : scalar
"""
function LS_prior(ninv::Int, α=[1.0, 0.5])
	Q=LinearMap(
	     (y,x)->LinearAlgebra.mul!(y,x,α[2]), 
	     (y,x)->LinearAlgebra.mul!(y,x,α[2]), 
		      ninv, ninv, ismutating=true, issymmetric=true)
	return LS_prior(α[1],Misfits.P_gls(Q))
end


include("core.jl")
include("prior.jl")
include("xfwi.jl")
include("xwfwi.jl")

"""
Convert the data `Data` to `SrcWav` after time reversal.
"""
function update_adjsrc!(adjsrc, δdat::Data, adjacqgeom)
	(adjsrc.fields != δdat.fields) && error("dissimilar fields")
	nt=length(δdat.tgrid)
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
	adjsrc=SrcWav_zeros(adjacqgeom, fields, tgrid)
	return adjsrc
end


"""
Constructor for `SeisInvExpt`

# Arguments

* `acqsrc::SrcWav` : source time functions
* `acqgeom::Geom` : acquisition geometry
* `tgrid::StepRangeLen` : modelling time grid
* `attrib_mod::Union{ModFdtd, ModFdtdBorn, ModFdtdHBorn}` : modelling attribute
* `modm::Medium` : seismic model on modelling mesh 

# Optional Arguments
* `tgrid_obs::StepRangeLen` : time grid for observed data

* `igrid::Vector{StepRangeLen}=modm.mgrid` : inversion grid if different from the modelling grid, i.e., `modm.mgrid`
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
* `acqsrc_obs::SrcWav=acqsrc` : source wavelets to generate *observed data*; can be different from `acqsrc`
* `modm_obs::Medium=modm` : actual seismic model to generate *observed data*
* `modm0::Medium=fill!(modm,0.0)` : background seismic model for Born modelling and inversion
* `parameterization::Vector{Symbol}` : subsurface parameterization
  * `=[:χvp, :χrho, :null]` for velocity and density inversion
  * `=[:null, :χrho, :null]` for density-only inversion
  * `=[:χvp, :null, :null]` for velocity-only inversion
* `verbose::Bool=false` : print level on STDOUT during inversion 
* `attrib::Symbol=:synthetic` : an attribute to control class
  * `=:synthetic` synthetic data inversion
  * `=:field` field data inversion
* nworker : number of workers (input nothing to use all available)
"""
function SeisInvExpt(
	       acqsrc::SrcWav,
	       acqgeom::Geom,
	       tgrid::StepRangeLen,
	       attrib_mod::Union{ModFdtd, ModFdtdBorn, ModFdtdHBorn}, 
	       modm::Medium;
	       # other optional 
	       tgrid_obs::StepRangeLen=deepcopy(tgrid),
	       rfields=[:P],
	       igrid=nothing,
	       igrid_interp_scheme::Symbol=:B2,
	       mprecon_factor::Float64=1.0,
	       dobs::Data=Data(tgrid_obs, acqgeom, rfields),
	       dprecon=nothing,
	       tlagssf_fracs=0.0,
	       tlagrf_fracs=0.0,
	       acqsrc_obs::SrcWav=deepcopy(acqsrc),
	       modm_obs::Medium=modm,
	       modm0=nothing,
	       parameterization::Vector{Symbol}=[:χvp, :χrho, :null],
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
	amin=max(igrid[1][1],mg[1][3])
	amax=min(igrid[1][end],mg[1][end-2])
	if(length(igrid[1])==1)
		grid1=range(amin,step=min(step(igrid[1]), mg[1][end-2]-amin),length=1)  
	else
		grid1=range(amin,stop=amax,length=length(igrid[1]))
	end
	amin=max(mg[2][3],igrid[2][1])
	amax=min(igrid[2][end],mg[2][end-2])
	if(length(igrid[2])==1)
		grid2=range(amin,step=min(step(igrid[2]), mg[2][end-2]-amin),length=1)  
	else
		grid2=range(amin,stop=amax,length=length(igrid[2]))
	end
	igrid=[grid1, grid2]


	# create modi according to igrid and interpolation of modm
	modi = Medium(igrid, modm);
	# update bounds
	update!(modi,modm)
	interp_spray!(modm, modi, :interp)

	# save the initial model, as modi will be changed during the iterations
	mod_initial=deepcopy(modi)

	# allocate prior inputs
	priori=similar(modi)

	# use default background model modm0
	if(modm0 === nothing)
		modm0=deepcopy(modm) # some dummy
		fill!(modm0)
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
	paf=SeisForwExpt(npw=2, model=modm, born_flag=true,
		acqgeom=[acqgeom, adjacqgeom], acqsrc=[acqsrc, adjsrc], 
		rfields=rfields,
		sflags=[3, 2], rflags=[1, 1],
		backprop_flag=1, 
		tgridmod=tgrid, gmodel_flag=true, verbose=verbose, illum_flag=false, nworker=nworker)


	gmodm=Medium(modm.mgrid,[:χvp,:χrho])
	gmodi=Medium(modi.mgrid,[:χvp,:χrho])

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

	paTD=VNamedD_misfit(Data(tgrid,acqgeom,rfields),dobs,w=dprecon);

	paminterp=Interpolation.Kernel([modi.mgrid[2], modi.mgrid[1]], [modm.mgrid[2], modm.mgrid[1]], igrid_interp_scheme)

	mx=X(prod(length.(modi.mgrid))*count(parameterization.≠:null),2)
	mxm=X(prod(length.(modm.mgrid))*count(parameterization.≠:null),2)


	# put modm as a vector, according to parameterization in mxm.x
	copyto!(mxm.x,modm,parameterization)

	pa = SeisInvExpt(
	     mx, mxm,
	     paf,
	     deepcopy(acqsrc), 
	     SrcWav_zeros(adjacqgeom, rfields, tgrid),
	     deepcopy(acqgeom), 
	     adjacqgeom,
	     deepcopy(modm), deepcopy(modm0), modi, mod_initial,
	     priori,
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

	# update priori in pa to defaults
	update_prior!(pa)
	
	return pa
end

"""
update the *synthetic* observed data in `SeisInvExpt`
* allocated memory, don't use in inner loops
"""
function update_observed_data!(pa::SeisInvExpt, modm_obs, attrib_mod=ModFdtd())
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
Return the number of inversion variables for FWI corresponding to `SeisInvExpt`.
This number of inversion variables depend on the size of inversion mesh.
"""
function xfwi_ninv(pa::SeisInvExpt)
	length(pa.mx.x)
end

"""
Return the number of inversion variables for source and receiver filter inversion 
corresponding to `SeisInvExpt`.
This number depends on the maximum lagtimes of the filters. 
"""
function wfwi_ninv(pa::SeisInvExpt)
	return  pa.paTD.coup.tgridssf.nx
end


"build a model precon"
function build_mprecon!(pa,illum::Array{Float64}, mprecon_factor=1.0)
	"illum cannot be zero when building mprecon"
	(any(illum .<= 0.0)) && error("illum cannot be negative or zero")
	(mprecon_factor < 1.0)  && error("invalid mprecon_factor")

	illumi = zeros(length(pa.modi.mgrid[1]), length(pa.modi.mgrid[2]))
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
	copyto!(pa.modi, x, pa.parameterization) 

	# get x on dense grid
	Interpolation.interp_spray!(x, pa.mxm.x, pa.paminterp, :interp, count(pa.parameterization.≠:null))

	# put dense x into modm, according to parameterization
	copyto!(pa.modm,pa.mxm.x,pa.parameterization) 
end


function δx_to_δmods!(pa, δx)
	# get x on dense grid
	Interpolation.interp_spray!(δx, pa.mxm.x, pa.paminterp, :interp, count(pa.parameterization.≠:null))

	# reparameterize accordingly to get [δKI, δrhoI]
	copyto!(pa.paf.c.δmod, pa.mxm.x, pa.modm0, pa.parameterization)

	# input [δKI, δrhoI] to the modeling engine
	Fdtd.update_δmods!(pa.paf.c, pa.paf.c.δmod)
end

"""
convert the gradient output from Fdtd to gx
"""
function spray_gradient!(gx,  pa, ::ModFdtd)

	# apply chain rule on gmodm to get gradient w.r.t. dense x
	chainrule!(pa.mxm.gx, pa.paf.c.gradient, pa.modm, pa.parameterization)

	# spray gradient w.r.t. dense x to the space of sparse x
	Interpolation.interp_spray!(gx, pa.mxm.gx, pa.paminterp, :spray, count(pa.parameterization.≠:null))

	# visualize gmodi here?
	# visualize_gx!(pa.gmodm, pa.modm, pa.gmodi, pa.modi, gx, pa)
end
function spray_gradient!(gx, pa, ::ModFdtdBorn)

	# apply chain rule on gmodm to get gradient w.r.t. dense x (note that background model is used for Born modeling here)
	chainrule!(pa.mxm.gx, pa.paf.c.gradient, pa.modm, pa.parameterization)

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
	nznx = prod(length.(pa.modi.mgrid))

	bound1 = similar(lower_x)
	# create a Seismic model with minimum possible values
	modbound = deepcopy(pa.modi);
	copyto!(modbound[:vp], fill(modbound.bounds[:vp][1], size(modbound[:vp])))
	if (:vs ∈ names(modbound.m)[1]) 
	#	copyto!(modbound[:vs], fill(modbound.bounds[:vs][1], size(modbound[:vs])))
	end
	copyto!(modbound[:rho], fill(modbound.bounds[:rho][1], size(modbound[:rho])))

	copyto!(bound1,modbound,pa.parameterization)


	bound2=similar(upper_x)
	# create a Seismic model with maximum possible values
	copyto!(modbound[:vp], fill(modbound.bounds[:vp][2], size(modbound[:vp])))
	if (:vs ∈ names(modbound.m)[1]) 
	#	copyto!(modbound[:vs], fill(modbound.bounds[:vs][2], size(modbound[:vs])))
	end
	copyto!(modbound[:rho], fill(modbound.bounds[:rho][2], size(modbound[:rho])))

	copyto!(bound2,modbound,pa.parameterization)

	# this sorting operation is important because 
	lower_x[:] = min.(bound1, bound2)
	upper_x[:] = max.(bound1, bound2)

	return modbound
end

# gmodi, gmodm <-- gx (just for visualizing gradient)
function visualize_gx!(gmodm, modm, gmodi, modi, gx, pa)
	# chain rule depending on re-parameterization
	chainrule!(gmodi, modi, gx, pa.parameterization, 1)

	# get gradient after interpolation (just for visualization, not exact)
#	interp_spray!(gmodi, gmodm, :interp, :B2, [:χvp,:χrho], pa=pa.paminterp)
end

"""
Modify the input acquisition geometry 
such that the adjoint source time functions can 
be propagated from the receiver positions.
The number of supersources will remain the same.
All the recievers will be fired as simultaneous sources.
"""
function AdjGeom(geomin::Geom)
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
function func_grad_Coupling!(storage, x::Vector{Float64},pa::SeisInvExpt)

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
		   pa::SeisInvExpt, 
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




