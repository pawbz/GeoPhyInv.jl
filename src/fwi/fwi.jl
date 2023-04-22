"""
This module defines a type called `PFWI` that is to be constructed 
before performing 
and kind of inversion.
The following functionality is currently added in this module:
* a simple RTM
* full-waveform inversion
* source and receiver filter inversion
Once a `PFWI` object is constructed, the following routines update certain fields
of the `PFWI` after performing the inversion.
* `xfwi!` : updates the initial subsurface models in `PFWI` using its observed data
* `wfwi!` : updates the source and receiver filters in `PFWI` using its observed data
"""

include("X.jl")


# add data inverse covariance matrix here later
struct LS end # just cls inversion

struct LS_prior # cls inversion including prior 
	# add data inverse covariance matrix here later, instead of just Float64
	pdgls::Float64
	# pdgls
	# pmgls::P_gls{Float64}
end
# LS_prior(α1::Float64, Q::LinearMap)=LS_prior(α1, P_gls(Q))
struct Migr end
struct Migr_FD end



# create a timer object, used for inversion, see TimerOutputs.jl
global const fwi_to = TimerOutput();

mutable struct PFWI{Tmod, Tdatamisfit, Tinv}
	attrib_inv::Tinv
	"model inversion variable"
	mx::X{Float64,1}
	mxm::X{Float64,1}
	"forward modelling parameters for"
	paf::PFdtd{Tmod}
	"base source wavelet"
	srcwav::SrcWav
	"adjoint source functions"
	adjsrc::SrcWav
	"acquisition ageometry"
	ageom::AGeom
	"acquisition ageometry for adjoint propagation"
	adjageom::AGeom
	"Seismic model on modeling grid"
	mediumm::Medium
	"background Seismic model"
	mediumm0::Medium
	"Seismic on inversion grid"
	mediumi::Medium
	"initialize on this model"
	mod_initial::Medium
	"prior model on the inversion grid"
	priori::Medium
	"gradient Seismic model on modeling grid"
	gmediumm::Medium
	"gradient Seismic on inversion grid"
	gmediumi::Medium
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


Base.show(io::Base.IO, pa::PFWI)=nothing

"""
```julia
pa=SeisInvExpt(attrib_mod, attrib_inv; srcwav, ageom, tgrid, mediumm, mediumm_obs, dobs)
```

# Arguments

* `attrib_mod::Union{FdtdAcoustic,FdtdAcoustic{Born}}` : choose modeling attribute
* `attrib_inv::Union{LS,LS_prior,Migr,Migr_FD}` : choose inversion attribute
* `obj::Union{LS,LS_prior}` : which objective function?
  * `=LS()` : least-squares inversion 
  * `=LS_prior([1.0, 0.5])` regularized least-squares inversion
  * `=Migr()` : migration image
  * `=Migr_FD()` : finite-difference gradient (used for testing only)

# Keyword Arguments
* `srcwav::SrcWav` : source wavelets
* `ageom::AGeom` : acquisition 
* `tgrid::StepRangeLen` : modelling time grid
* `mediumm::Medium` : an instance of `Medium` (used to generate initial model on the inversion grid)
* `dobs::Records` : observed data (required only when `attrib==:field`)
* `mediumm_obs::Medium` : medium used for generating *observed data* (required only when `attrib==:synthetic`)

# Optional Keyword Arguments 

* `igrid` : inversion grid (defaults to modeling grid)
* `igrid_interp_scheme` : interpolation scheme
  * `=:B1` linear 
  * `=:B2` second order 
* `mediumm0` : background seismic model for Born modelling and inversion
* `parameterization::Vector{Symbol}` : subsurface parameterization
  * `=[:χvp, :χrho, :null]` for vp and density contrast inversion
  * `=[:null, :χrho, :null]` for density-only contrast inversion
  * `=[:χvp, :null, :null]` for vp-only contrast  inversion
* `tgrid_obs::StepRangeLen` : time grid for observed data
* `attrib::Symbol=:synthetic` : an attribute to control class
  * `=:synthetic` synthetic data inversion
  * `=:field` field data inversion
"""
SeisInvExpt(m1::Union{FdtdAcoustic, FdtdAcoustic{Born}},m2::Union{LS,LS_prior,Migr,Migr_FD};args1...)=PFWI(m1,m2;args1...) # change this when you have more inversion experiments 

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
	return LS_prior(α[1],P_gls(Q))
end


include("getprop.jl")
include("core.jl")
include("prior.jl")
include("xfwi.jl")
include("xwfwi.jl")
include("gallery.jl")

"""
Convert the data `Records` to `SrcWav` after time reversal.
"""
function update_adjsrc!(adjsrc, δdat::Records, adjageom)
	for i in 1:length(adjageom)
		for field in names(δdat[i].d)[1]
			wav=adjsrc[i].d[field] 
			dat=δdat[i].d[field]
			for ir in 1:δdat[i][:nr]
				nt=length(adjsrc[i].grid)  
				for it in 1:nt
					@inbounds wav[it,ir]=value_adjsrc(dat[nt-it+1,ir],eval(field)())
				end
			end
		end
	end
	return nothing
end

value_adjsrc(s, ::p) = s
value_adjsrc(s, ::vx) = -1.0*s # see adj tests
value_adjsrc(s, ::vz) = -1.0*s


function generate_adjsrc(fields, tgrid, adjageom)
	adjsrc=SrcWav(tgrid,adjageom,fields)
	return adjsrc
end


"""
Constructor for `PFWI`

# Arguments

* `attrib_mod::Union{FdtdAcoustic, FdtdAcoustic{Born}}` : modeling attribute
* `attrib_inv::Union{LS,LS_prior,Migr,Migr_FD} : inversion attribute


# Keyword Arguments
* `srcwav::SrcWav` : source time functions
* `ageom::AGeom` : acquisition ageometry
* `tgrid::StepRangeLen` : modeling time grid
* `mediumm::Medium` : seismic model on modeling mesh 
* `mediumm_obs::Medium=mediumm` : actual seismic model to generate *observed* data

# Optional Keyword Arguments 
* `tgrid_obs::StepRangeLen` : time grid for observed data
* `igrid::Vector{StepRangeLen}=mediumm.mgrid` : inversion grid if different from the modelling grid, i.e., `mediumm.mgrid`
* `igrid_interp_scheme` : interpolation scheme
  * `=:B1` linear 
  * `=:B2` second order 
* `mprecon_factor::Float64=1` : factor to control model preconditioner, always greater than 1
  * `=1` means the preconditioning is switched off
  * `>1` means the preconditioning is switched on, larger the mprecon_factor, stronger the applied preconditioner becomes
* `dobs::Records` : observed data
* `dprecon::Records=ones(Records)` : data preconditioning, defaults to one 
* `tlagssf_fracs=0.0` : maximum lagtime of unknown source filter
* `tlagrf_fracs=0.0` : maximum lagtime of unknown receiver filter
* `srcwav_obs::SrcWav=srcwav` : source wavelets to generate *observed data*; can be different from `srcwav`
* `mediumm0::Medium=fill!(mediumm,0.0)` : background seismic model for Born modelling and inversion
* `parameterization::Vector{Symbol}` : subsurface parameterization
  * `=[:χvp, :χrho, :null]` for velocity and density inversion
  * `=[:null, :χrho, :null]` for density-only inversion
  * `=[:χvp, :null, :null]` for velocity-only inversion
* `verbose::Bool=false` : print level on STDOUT during inversion 
* nworker : number of workers (input nothing to use all available)
"""
function PFWI(
	      attrib_mod::Union{FdtdAcoustic, FdtdAcoustic{Born}}, 
	      attrib_inv::Union{LS,LS_prior,Migr,Migr_FD};
	      srcwav::SrcWav=nothing,
	      ageom::AGeom=nothing,
	      tgrid::StepRangeLen=nothing,
	      mediumm::Medium=nothing,
	      # other optional 
	      tgrid_obs::StepRangeLen=deepcopy(tgrid),
	      rfields=[:p],
	      igrid=nothing,
	      igrid_interp_scheme::Symbol=:B2,
	      mprecon_factor::Float64=1.0,
	      dobs::Records=Records(tgrid_obs, ageom, rfields),
	      dprecon=nothing,
	      tlagssf_fracs=0.0,
	      tlagrf_fracs=0.0,
	      srcwav_obs::SrcWav=deepcopy(srcwav),
	      mediumm_obs::Medium=mediumm,
	      mediumm0=nothing,
	      parameterization::Vector{Symbol}=[:χvp, :χrho, :null],
	      verbose::Bool=false,
	      attrib::Symbol=:synthetic,
	      nworker=nothing,
	      )

	# make a copy of initial model
	mediumm=deepcopy(mediumm)

	if(igrid===nothing)
		igrid=deepcopy(mediumm.mgrid)
	end
	# igrid has to truncated because the gradient evaluation 
	# is inaccurate on boundaries
	mg=mediumm.mgrid
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


	# create mediumi according to igrid and interpolation of mediumm
	mediumi = Medium(igrid, mediumm);
	# update bounds
	update!(mediumi,mediumm)
	interp_spray!(mediumm, mediumi, :interp)

	# save the initial model, as mediumi will be changed during the iterations
	mod_initial=deepcopy(mediumi)

	# allocate prior inputs
	priori=similar(mediumi)

	# use default background model mediumm0
	if(mediumm0 === nothing)
		mediumm0=deepcopy(mediumm) # some dummy
		fill!(mediumm0)
	end
	isequal(mediumm0, mediumm_obs) && (@warn "mediumm0 == mediumm_obs")

	if(attrib == :synthetic) 
		isequal(mediumm, mediumm_obs) && error("initial model same as actual model, zero misfit?")
	end

	# ageom ageometry for adjoint propagation
	adjageom = AdjAGeom(ageom)

	# generate adjoint sources
	adjsrc=generate_adjsrc(rfields, tgrid, adjageom)


	# generating forward and adjoint modelling engines
	# to generate modelled data, border values, etc.
	# most of the parameters given to this are dummy
	paf=SeisForwExpt(attrib_mod,npw=2, medium=mediumm,
		ageom=[ageom, adjageom], srcwav=[srcwav, adjsrc], 
		rfields=rfields,
		sflags=[3, 2], rflags=[1, 1],
		backprop_flag=:save, 
		tgrid=tgrid, gmodel_flag=true, verbose=verbose, illum_flag=false, nworker=nworker)


	gmediumm=Medium(mediumm.mgrid,[:χvp,:χrho])
	gmediumi=Medium(mediumi.mgrid,[:χvp,:χrho])

	#println("added fake precon")
	#dprecon=deepcopy(dobs)
	#fill!(dprecon, 1.0)
	#Records.taper!(dprecon, 20.)

	# check dprecon
	if(!(dprecon===nothing))
		(!(isapprox(dprecon,dobs))) && error("invalid dprecon used")
	end
	
	# create Parameters for data misfit
	#coup=Coupling.TD_delta(dobs.tgrid, tlagssf_fracs, tlagrf_fracs, rfields, ageom)
	# choose the data misfit
	#	(iszero(tlagssf_fracs)) 
 	# tlagssf_fracs==[0.0] | tlagssf_fracs=[0.0])
	# paTD=Records.P_misfit(Records.TD_zeros(rfields,tgrid,ageom),dobs,w=dprecon,coup=coup, func_attrib=optims[1]);

	if((attrib == :synthetic))
		if(iszero(dobs))
			Random.randn!(dobs) # dummy dobs, update later
		end
	else
		if(iszero(dobs))
			error("input observed data is zero")
		end
	end

	paTD=VNamedD_misfit(Records(tgrid,ageom,rfields),dobs,w=dprecon);

	paminterp=Interpolation.Kernel([mediumi.mgrid[2], mediumi.mgrid[1]], [mediumm.mgrid[2], mediumm.mgrid[1]], igrid_interp_scheme)

	mx=X(prod(length.(mediumi.mgrid))*count(parameterization.≠:null),2)
	mxm=X(prod(length.(mediumm.mgrid))*count(parameterization.≠:null),2)


	# put mediumm as a vector, according to parameterization in mxm.x
	copyto!(mxm.x,mediumm,parameterization)

	pa = PFWI(
		  attrib_inv,
	     mx, mxm,
	     paf,
	     deepcopy(srcwav), 
	     SrcWav(tgrid,adjageom, rfields),
	     deepcopy(ageom), 
	     adjageom,
	     deepcopy(mediumm), deepcopy(mediumm0), mediumi, mod_initial,
	     priori,
	     gmediumm,gmediumi,
	     parameterization,
	     zeros(2,2), # dummy, update mprecon later
	     paTD,
	     paminterp,
	     verbose, attrib)

	# generate observed data if attrib is synthetic
	if((attrib == :synthetic))
		update_observed_data!(pa, mediumm_obs)
	end

	iszero(dobs) && ((attrib == :real) ? error("input observed data for real data inversion") : error("problem generating synthetic observed data"))

	# generate modelled data
	F!(pa, nothing)

#	build_mprecon!(pa, Array(pa.paf.c.illum_stack), mprecon_factor)
	pa.paf.c.illum_flag=false # switch off illum flag for speed

	# default weights are 1.0
	fill!(pa.mx.w, 1.0)

	# update priori in pa to defaults
	update_prior!(pa)
	
	return pa
end

"""
update the *synthetic* observed data in `PFWI`
* allocated memory, don't use in inner loops
"""
function update_observed_data!(pa::PFWI, mediumm_obs)
	# save mediumm of pa to put it back later
	mediumm_copy=deepcopy(pa.mediumm)

	# change mediumm in pa to actual model
	copyto!(pa.mediumm, mediumm_obs)

	# update models in the forward engine and do modelling
	# FdtdAcoustic{Born}: mediumm0 will be used as a background model for Born modelling and mediumm will be perturbed model
	# Fdtd: mediumm will be *the* model
	F!(pa, nothing)

	# get observed data
	copyto!(pa.paTD.y, pa.paTD.x)

	# put back model and sources
	copyto!(pa.mediumm, mediumm_copy)

	initialize!(pa.paf.c)  # clear everything
end



"""
Return the number of inversion variables for PFWI corresponding to `PFWI`.
This number of inversion variables depend on the size of inversion mesh.
"""
function xfwi_ninv(pa::PFWI)
	length(pa.mx.x)
end

"""
Return the number of inversion variables for source and receiver filter inversion 
corresponding to `PFWI`.
This number depends on the maximum lagtimes of the filters. 
"""
function wfwi_ninv(pa::PFWI)
	return  pa.paTD.coup.tgridssf.nx
end


"build a model precon"
function build_mprecon!(pa,illum::Array{Float64}, mprecon_factor=1.0)
	"illum cannot be zero when building mprecon"
	(any(illum .<= 0.0)) && error("illum cannot be negative or zero")
	(mprecon_factor < 1.0)  && error("invalid mprecon_factor")

	illumi = zeros(length(pa.mediumi.mgrid[1]), length(pa.mediumi.mgrid[2]))
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


# convert x to mediumm by interpolation
# x is on the sparse grid, and is changed by parameterization
function x_to_mediumm!(pa, x)

	# put sparse x into mediumi (useless step, just for visualization)
	copyto!(pa.mediumi, x, pa.parameterization) 

	# get x on dense grid
	Interpolation.interp_spray!(x, pa.mxm.x, pa.paminterp, :interp, count(pa.parameterization.≠:null))

	# put dense x into mediumm, according to parameterization
	copyto!(pa.mediumm,pa.mxm.x,pa.parameterization) 
end


function δx_to_δmods!(pa, δx)
	# get x on dense grid
	Interpolation.interp_spray!(δx, pa.mxm.x, pa.paminterp, :interp, count(pa.parameterization.≠:null))

	# reparameterize accordingly to get [δKI, δrhoI]
	copyto!(pa.paf.c.δmodall, pa.mxm.x, pa.mediumm0, pa.parameterization)

	# input [δKI, δrhoI] to the modeling engine
	update_δmods!(pa.paf.c, pa.paf.c.δmodall)
end

"""
convert the gradient output from Fdtd to gx
"""
function spray_gradient!(gx,  pa::PFWI{FdtdAcoustic,T1,T2}) where {T1,T2}

	# apply chain rule on gmediumm to get gradient w.r.t. dense x
	chainrule!(pa.mxm.gx, pa.paf.c.gradient, pa.mediumm, pa.parameterization)

	# spray gradient w.r.t. dense x to the space of sparse x
	Interpolation.interp_spray!(gx, pa.mxm.gx, pa.paminterp, :spray, count(pa.parameterization.≠:null))

	# visualize gmediumi here?
	# visualize_gx!(pa.gmediumm, pa.mediumm, pa.gmediumi, pa.mediumi, gx, pa)
end
function spray_gradient!(gx, pa::PFWI{FdtdAcoustic{Born},T1,T2}) where {T1,T2}   

	# apply chain rule on gmediumm to get gradient w.r.t. dense x (note that background model is used for Born modeling here)
	chainrule!(pa.mxm.gx, pa.paf.c.gradient, pa.mediumm, pa.parameterization)

	# spray gradient w.r.t. dense x to the space of sparse x
	Interpolation.interp_spray!(gx, pa.mxm.gx, pa.paminterp, :spray, count(pa.parameterization.≠:null))

	# visualize gmediumi here?
	# visualize_gx!(pa.gmediumm, pa.mediumm, pa.gmediumi, pa.mediumi, gx, pa)
end



"""
Return bound vectors for the `Medium` model, 
depeding on paramaterization
"""
function Seismic_xbound!(lower_x, upper_x, pa)
	nznx = prod(length.(pa.mediumi.mgrid))

	bound1 = similar(lower_x)
	# create a Seismic model with minimum possible values
	modbound = deepcopy(pa.mediumi);
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

# gmediumi, gmediumm <-- gx (just for visualizing gradient)
function visualize_gx!(gmediumm, mediumm, gmediumi, mediumi, gx, pa)
	# chain rule depending on re-parameterization
	chainrule!(gmediumi, mediumi, gx, pa.parameterization, 1)

	# get gradient after interpolation (just for visualization, not exact)
#	interp_spray!(gmediumi, gmediumm, :interp, :B2, [:χvp,:χrho], pa=pa.paminterp)
end

"""
Modify the input acquisition ageometry 
such that the adjoint source time functions can 
be propagated from the receiver positions.
The number of supersources will remain the same.
All the recievers will be fired as simultaneous sources.
"""
function AdjAGeom(ageomin::AGeom)
	ageomout = deepcopy(ageomin);
	for i in 1:length(ageomin)
		ageomout[i].s[:x] = ageomin[i].r[:x]; ageomout[i].s[:z] = ageomin[i].r[:z];
		ageomout[i].ns = ageomin[i].nr; 

	#	ageomout.r[:x] = ageomin.s[:x]; ageomout.r[:z] = ageomin.s[:z];
	#	ageomout.nr = ageomin.ns; 

		ageomout[i].r[:x] = ageomin[i].r[:x]; ageomout[i].r[:z] = ageomin[i].r[:z];
		ageomout[i].nr = ageomin[i].nr; 
	end

	return ageomout
end


"""
Return functional and gradient of the LS objective 
"""
function func_grad_Coupling!(storage, x::Vector{Float64},pa::PFWI)

	pa.verbose ? println("computing gradient...") : nothing

	# x to w 
	Coupling_x!(x, pa, -1)

	if(storage === nothing)
		# compute misfit 
		f = Records.error!(pa.paTD)
		return f
	else
		f = Records.error!(pa.paTD, :dJssf)
		Coupling_gx!(storage, pa.paTD.dJssf)
		return storage
	end
end

"""
Convert coupling functions to x and vice versa
"""
function Coupling_x!(x::Vector{Float64}, 
		   pa::PFWI, 
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
		for iss=1:pa.paTD.y.ageom.nss, ifield=1:length(pa.paTD.y.fields)
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




