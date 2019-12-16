# 
# # Credits 
# Author: Pawan Bharadwaj 
#         (bharadwaj.pawan@gmail.com)
# 
# * original code in FORTRAN90: March 2013
# * modified: 11 Sept 2013
# * major update: 25 July 2014
# * code optimization with help from Jan Thorbecke: Dec 2015
# * rewritten in Julia: June 2017
# * added parrallelization over supersources in Julia: July 2017
# * efficient parrallelization using distributed arrays: Sept 2017
# * optimized memory allocation: Oct 2017
# 









"""
Forward modeling using a finite-difference simulation of the acoustic wave-equation.
"""
struct Fdtd end


"""
Linearized forward modeling using a finite-difference simulation of the acoustic wave-equation.
"""
struct FdtdBorn end




global const npml = 50
global const nlayer_rand = 0

# 
#As forward modeling method, the 
#finite-difference method is employed. 
#It uses a discrete version of the two-dimensional isotropic acoustic wave equation.
#
#```math
#\pp[\tzero] - \pp[\tmo] = \dt \mB \left({\partial_x\vx}[\tmh]
# + \partial_z \vz[\tmh]  + \dt\sum_{0}^{\tmo}\sfo\right)
# ```
# ```math
#\pp[\tpo] - \pp[\tzero] = \dt \mB \left(\partial_x \vx[\tph]
# + {\partial_z \vz}[\tph]  + \dt\sum_{0}^{\tzero}\sfo\right)
# ```

"""
Modelling parameters common for all supersources
# Keyword Arguments that are modified by the method (some of them are returned as well)

* `gradient::Vector{Float64}=Medium(model.mgrid)` : gradient model modified only if `gmodel_flag`
* `TDout::Vector{Data.TD}=[Data.TD_zeros(rfields,tgridmod,ageom[ip]) for ip in 1:length(findn(rflags))]`
* `illum::Array{Float64,2}=zeros(length(model.mgrid[1]), length(model.mgrid[2]))` : source energy if `illum_flag`
* `boundary::Array{Array{Float64,4},1}` : stored boundary values for first propagating wavefield 
* `snaps::Array{Float64,4}=zeros(length(model.mgrid[1]),length(model.mgrid[2]),length(tsnaps),ageom[1].nss)` :snapshots saved at `tsnaps`

# Return (in order)

* modelled data for each propagating wavefield as `Vector{TD}`
* stored boundary values of the first propagating wavefield as `Array{Array{Float64,4},1}` (use for backpropagation)
* final conditions of the first propagating wavefield as `Array{Float64,4}` (use for back propagation)
* gradient model as `Seismic`
* stored snaps shots at tsnaps as Array{Float64,4} 
"""
mutable struct PFdtdc{T}
	jobname::Symbol
	attrib_mod::T
	npw::Int64 
	activepw::Vector{Int64}
	exmodel::Medium
	model::Medium
	ageom::Vector{AGeom}
	srcwav::Vector{SrcWav}
	abs_trbl::Vector{Symbol}
	sfields::Vector{Vector{Symbol}}
	isfields::Vector{Vector{Int64}}
	sflags::Vector{Int64} 
	irfields::Vector{Int64}
	rflags::Vector{Int64}
	δt::Float64
	δxI::Float64
	δzI::Float64
	nx::Int64
	nz::Int64
	nt::Int64
	δtI::Float64 
	δx24I::Float64
	δz24I::Float64
	a_x::Vector{Float64}
	b_x::Vector{Float64}
	k_xI::Vector{Float64}
	a_x_half::Vector{Float64}
	b_x_half::Vector{Float64}
	k_x_halfI::Vector{Float64} 
	a_z::Vector{Float64}
	b_z::Vector{Float64}
	k_zI::Vector{Float64}
	a_z_half::Vector{Float64}
	b_z_half::Vector{Float64}
	k_z_halfI::Vector{Float64}
	modtt::Matrix{Float64}
	modttI::Matrix{Float64} # just storing inv(modtt) for speed
	modrr::Matrix{Float64}
	modrrvx::Matrix{Float64}
	modrrvz::Matrix{Float64}
	δmodtt::Matrix{Float64}
	δmodrr::Matrix{Float64} 
	δmodrrvx::Matrix{Float64}
	δmodrrvz::Matrix{Float64}
	δmod::Vector{Float64} # perturbation vector (KI, rhoI)
	gradient::Vector{Float64}  # output gradient vector w.r.t (KI, rhoI)
	grad_modtt_stack::SharedArrays.SharedArray{Float64,2} # contains gmodtt
	grad_modrrvx_stack::SharedArrays.SharedArray{Float64,2}
	grad_modrrvz_stack::SharedArrays.SharedArray{Float64,2}
	grad_modrr_stack::Array{Float64,2}
	illum_flag::Bool
	illum_stack::SharedArrays.SharedArray{Float64,2}
	backprop_flag::Int64
	snaps_flag::Bool
	itsnaps::Vector{Int64}
	gmodel_flag::Bool
	ibx0::Int64
	ibz0::Int64
	ibx1::Int64
	ibz1::Int64
	isx0::Int64
	isz0::Int64
	datamat::SharedArrays.SharedArray{Float64,3}
	data::Vector{Data}
	verbose::Bool
end 

"""
Modelling parameters per every supersource for each worker
"""
mutable struct PFdtdss
	iss::Int64
	wavelets::Matrix{Matrix{Float64}}
	ssprayw::Vector{Matrix{Float64}}
	records::Vector{Array{Float64,3}}
	rinterpolatew::Vector{Matrix{Float64}}
	isx1::Vector{Vector{Int64}} 
	isx2::Vector{Vector{Int64}}
	isz1::Vector{Vector{Int64}} 
	isz2::Vector{Vector{Int64}}
	irx1::Vector{Vector{Int64}}
	irx2::Vector{Vector{Int64}}
	irz1::Vector{Vector{Int64}}
	irz2::Vector{Vector{Int64}}
        boundary::Vector{Array{Float64,3}}
	snaps::Array{Float64,3}
	illum::Matrix{Float64}
	grad_modtt::Matrix{Float64} 
	grad_modrrvx::Matrix{Float64}
	grad_modrrvz::Matrix{Float64}
end

"""
PFdtdeters per every worker, not necessarily for every supersource.
Note that a single worker can take care of multiple supersources.
"""
mutable struct PFdtdp
	ss::Vector{PFdtdss}
	p::Vector{Array{Float64,3}}
	pp::Vector{Array{Float64,3}}
	ppp::Vector{Array{Float64,3}}
	dpdx::Vector{Array{Float64,3}}
	dpdz::Vector{Array{Float64,3}}
	memory_dp_dx::Vector{Array{Float64,2}}
	memory_dp_dz::Vector{Array{Float64,2}}
	memory_dvx_dx::Vector{Array{Float64,2}}
	memory_dvz_dz::Vector{Array{Float64,2}}
	born_svalue_stack::Matrix{Float64} 
end


mutable struct PFdtd{T}
	sschunks::Vector{UnitRange{Int64}} # how supersources are distributed among workers
	p::DistributedArrays.DArray{PFdtdp,1,PFdtdp} # distributed parameters among workers
	c::PFdtdc{T} # common parameters
end


"""
```julia
pa = SeisForwExpt(attrib_mod; ageom, srcwav, medium, tgrid);
```
Method to create an instance of `SeisForwExpt`. 
The output of this method can be used as an input to the in-place method `update!`, to actually perform a
finite-difference modeling.

# Keyword Arguments

* `attrib_mod` : attribute to choose the type of modeling. Choose from 
  * `=Fdtd()` for full wavefield modeling  (finite-difference simulation of the acoustic wave-equation)
  * `=FdtdBorn()` for Born modeling 
* `model::Medium` : medium parameters 
* `tgrid` : modeling time grid, maximum time in `tgrid`should be greater than or equal to maximum source time, same sampling interval as in `srcwav`
* `ageom` :  acquisition 
* `srcwav` : source wavelets  

# Optional Keyword Arguments 

* `sflags=2` : source related flags 
  * `=0` inactive sources
  * `=1` sources with injection rate
  * `=2` volume injection sources
  * `=3` sources input after time reversal (use only during backpropagation)
* `rflags=1` : receiver related flags 
  * `=0` receivers do not record (or) inactive receivers
  * `=1` receivers are active only for the second propagating wavefield
* `rfields=[:P]` : multi-component receiver flag. Choose `Vector{Symbol}`
  * `=[:P]` record pressure 
  * `=[:Vx]` record horizontal component of particle velocity  
  * `=[:Vz]` record vertical component of particle velocity  
  * `=[:P, :Vx]` record both pressure and velocity 
* `tsnaps` : store snaps at these modeling times (defaults to half time)
  * `=[0.1,0.2,0.3]` record at these instances of tgrid
* `snaps_flag::Bool=false` : return snapshots or not
* `verbose::Bool=false` : verbose flag
"""
SeisForwExpt(attrib_mod::Union{Fdtd,FdtdBorn},args1...;args2...)=PFdtd(attrib_mod,args1...;args2...)


function initialize!(pap::PFdtdp)
	reset_per_ss!(pap)
	for issp in 1:length(pap.ss)
		pass=pap.ss[issp]
		for i in 1:length(pass.records)
			fill!(pass.records[i],0.0)
		end
		fill!(pass.snaps,0.0)
		fill!(pass.illum,0.0)
		fill!(pass.grad_modtt,0.0)
		fill!(pass.grad_modrrvx,0.0)
		fill!(pass.grad_modrrvz,0.0)
	end
end

function iszero_boundary(pa)
	result=false
	@sync begin
		for (ip, p) in enumerate(procs(pa.p))
			@async remotecall_wait(p) do 
				pap=localpart(pa.p)
				for issp in 1:length(pap.ss)
					pass=pap.ss[issp]
					for i in 1:length(pass.boundary)
						result = result | iszero(pass.boundary[i])
					end
				end
			end
		end
	end
	return result
end
		

function initialize_boundary!(pa)
	# zero out results stored per worker
	@sync begin
		for (ip, p) in enumerate(procs(pa.p))
			@async remotecall_wait(p) do 
				pap=localpart(pa.p)
				for issp in 1:length(pap.ss)
					pass=pap.ss[issp]
					for i in 1:length(pass.boundary)
						fill!(pass.boundary[i],0.0)
					end
				end
			end
		end
	end
end

function initialize!(pac::PFdtdc)
	fill!(pac.gradient,0.0)
	fill!(pac.grad_modtt_stack,0.0)
	fill!(pac.grad_modrrvx_stack,0.0)
	fill!(pac.grad_modrrvz_stack,0.0)
	fill!(pac.grad_modrr_stack,0.0)
	fill!(pac.illum_stack,0.0)
	for dat in pac.data
		fill!(dat, 0.0)
	end
	fill!(pac.datamat,0.0)
end

function reset_per_ss!(pap::PFdtdp)
	for ipw in 1:length(pap.p)
		fill!(pap.born_svalue_stack,0.0)
		fill!(pap.p[ipw],0.0)
		fill!(pap.pp[ipw],0.0)
		fill!(pap.ppp[ipw],0.0)
		fill!(pap.dpdx[ipw],0.0)
		fill!(pap.dpdz[ipw],0.0)
		fill!(pap.memory_dp_dz[ipw],0.0)
		fill!(pap.memory_dp_dx[ipw],0.0)
		fill!(pap.memory_dvx_dx[ipw],0.0)
		fill!(pap.memory_dvz_dz[ipw],0.0)
	end
end



"""
Primary method to generate Expt variable when Fdtd() and FdtdBorn() are used.

# Some internal arguments
* `jobname::String` : name
* `npw::Int64=1` : number of independently propagating wavefields in `model`
* `backprop_flag::Bool=Int64` : save final state variables and the boundary conditions for later use
  * `=1` save boundary and final values in `boundary` 
  * `=-1` use stored values in `boundary` for back propagation
* `gmodel_flag=false` : flag that is used to output gradient; there should be atleast two propagating wavefields in order to do so: 1) forward wavefield and 2) adjoint wavefield
* `illum_flag::Bool=false` : flag to output wavefield energy or source illumination; it can be used as preconditioner during inversion
* `abs_trbl::Vector{Symbol}=[:top, :bottom, :right, :left]` : use absorbing PML boundary conditions or not
  * `=[:top, :bottom]` apply PML conditions only at the top and bottom of the model 
  * `=[:bottom, :right, :left]` top is reflecting
"""
function PFdtd(attrib_mod;
	jobname::Symbol=:forward_propagation,
	npw::Int64=1, 
	medium::Medium=nothing,
	abs_trbl::Vector{Symbol}=[:top, :bottom, :right, :left],
	tgrid::StepRangeLen=nothing,
	ageom::Union{AGeom,Vector{AGeom}}=nothing,
	srcwav::Union{SrcWav,Vector{SrcWav}}=nothing,
	sflags::Union{Int,Vector{Int}}=fill(2,npw), 
	rflags::Union{Int,Vector{Int}}=fill(1,npw),
	rfields::Vector{Symbol}=[:P], 
	backprop_flag::Int64=0,  
	gmodel_flag::Bool=false,
	illum_flag::Bool=false,
	tsnaps::Vector{Float64}=fill(0.5*(tgrid[end]+tgrid[1]),1),
	snaps_flag::Bool=false,
	verbose::Bool=false,
	nworker=nothing)

	# convert to vectors 
	if(typeof(ageom)==AGeom); ageom=[ageom]; end
	if(typeof(srcwav)==SrcWav); srcwav=[srcwav]; end
	if(typeof(sflags)==Int); sflags=[sflags]; end
	if(typeof(rflags)==Int); rflags=[rflags]; end

	# alias
	model=medium
	tgridmod=tgrid

	#println("********PML Removed*************")
	#abs_trbl=[:null]

	# check sizes and errors based on input
	#(length(TDout) ≠ length(findn(rflags))) && error("TDout dimension")
	(length(ageom) ≠ npw) && error("ageom dimension")
	(length(srcwav) ≠ npw) && error("srcwav dimension")
	(length(sflags) ≠ npw) && error("sflags dimension")
	(length(rflags) ≠ npw) && error("rflags dimension")
	(maximum(tgridmod) < maximum(srcwav[1][1].grid)) && error("modeling time is less than source time")
	#(any([getfield(TDout[ip],:tgrid).δx < tgridmod.δx for ip=1:length(TDout)])) && error("output time grid sampling finer than modeling")
	#any([maximum(getfield(TDout[ip],:tgrid).x) > maximum(tgridmod) for ip=1:length(TDout)]) && error("output time > modeling time")

	#! no modeling if source wavelet is zero
	#if(maxval(abs(src_wavelets)) .lt. tiny(rzero_de)) 
	#        return
	#endif

	# all the propagating wavefields should have same supersources? check that?

	# check dimension of model

	# check if all sources are receivers are inside model
	any(.![(ageom[ip] ∈ model.mgrid) for ip=1:npw]) ? error("sources or receivers not inside model") : nothing


	length(ageom) != npw ? error("ageom size") : nothing
	length(srcwav) != npw ? error("srcwav size") : nothing
	all([issimilar(ageom[ip],srcwav[ip]) for ip=1:npw]) ? nothing : error("ageom and srcwav mismatch") 

	# necessary that nss and fields should be same for all nprop
	nss = length(ageom[1]);
	sfields=[names(srcwav[ipw][1].d)[1] for ipw in 1:npw]
	isfields = [Array{Int}(undef,length(sfields[ipw])) for ipw in 1:npw]
	for ipw in 1:npw
		for (i,field) in enumerate(sfields[ipw])
			isfields[ipw][i]=findfirst(isequal(field), [:P,:Vx,:Vz])
		end
	end
	fill(nss, npw) != [length(ageom[ip]) for ip=1:npw] ? error("different supersources") : nothing

	# create acquisition ageometry with each source shooting 
	# at every unique receiver position
	irfields = Array{Int}(undef,length(rfields))
	for (i,field) in enumerate(rfields)
		irfields[i]=findfirst(isequal(field), [:P,:Vx,:Vz])
	end


	#(verbose) &&	println(string("\t> number of super sources:\t",nss))	

	# find maximum and minimum frequencies in the source wavelets
	freqmin = Utils.findfreq(srcwav[1][1].d[1][:,1],srcwav[1][1].grid,attrib=:min) 
	freqmax = Utils.findfreq(srcwav[1][1].d[1][:,1],srcwav[1][1].grid,attrib=:max) 

	# minimum and maximum velocities
	vpmin = minimum(broadcast(minimum,[model.bounds[:vp]]))
	vpmax = maximum(broadcast(maximum,[model.bounds[:vp]]))
	#verbose && println("\t> minimum and maximum velocities:\t",vpmin,"\t",vpmax)


	check_fd_stability(vpmin, vpmax, step(model.mgrid[2]), step(model.mgrid[1]), freqmin, freqmax, step(tgridmod), verbose)


	# where to store the boundary values (careful, born scaterrers cannot be inside these!?)
	ibx0=npml-2; ibx1=length(model.mgrid[2])+npml+3
	ibz0=npml-2; ibz1=length(model.mgrid[1])+npml+3
#	ibx0=npml+1; ibx1=length(model.mgrid[2])+npml
#	ibz0=npml+1; ibz1=length(model.mgrid[1])+npml
#	println("**** Boundary Storage Changed **** ")

	# for snaps
	isx0, isz0=npml, npml

	# extend models in the PML layers
	exmodel = Medium_pml_pad_trun(model, nlayer_rand, npml);

	"density values on vx and vz stagerred grids"
	modrr=copy(exmodel[:rhoI])
	modrrvx = get_rhovxI(modrr)
	modrrvz = get_rhovzI(modrr)

	modttI = copy(exmodel[:K]) 
	modtt = copy(exmodel[:KI])

	if(typeof(attrib_mod)==FdtdBorn)
		(npw≠2) && error("born modeling needs npw=2")
	end

	#create some aliases
	nz, nx = length.(exmodel.mgrid)
	nzd, nxd = length.(model.mgrid)
	δx, δz = step(exmodel.mgrid[2]), step(exmodel.mgrid[1])
	δt = step(tgridmod)
	nt=length(tgridmod)

	δx24I, δz24I = inv(δx * 24.0), inv(δz * 24.0)
	δxI, δzI = inv(δx), inv(δz)
	δt = step(tgridmod)
	δtI = inv(δt)
	nt=length(tgridmod)
	mesh_x, mesh_z = exmodel.mgrid[2], exmodel.mgrid[1]


	# perturbation vectors are required on full space
	δmodtt = zeros(nz, nx)
	δmod = zeros(2*nzd*nxd)
	δmodrr = zeros(nz, nx)
	δmodrrvx = zeros(nz, nx)
	δmodrrvz = zeros(nz, nx)


	# pml_variables
	a_x, b_x, k_xI, a_x_half, b_x_half, k_x_halfI = pml_variables(nx, δt, δx, npml-3, vpmax, vpmin, freqmin, freqmax, 
							       [any(abs_trbl .== :left), any(abs_trbl .== :right)])
	a_z, b_z, k_zI, a_z_half, b_z_half, k_z_halfI = pml_variables(nz, δt, δz, npml-3, vpmax, vpmin, freqmin, freqmax,
							       [any(abs_trbl .== :top), any(abs_trbl .== :bottom)])

	gradient=zeros(2*nzd*nxd)
	grad_modtt_stack=SharedMatrix{Float64}(zeros(nz,nx))
	# these are full size, as rhoI -> rhovx and rhovz is also full size
	grad_modrrvx_stack=SharedMatrix{Float64}(zeros(nz,nx))
	grad_modrrvz_stack=SharedMatrix{Float64}(zeros(nz,nx))
	grad_modrr_stack=zeros(nz,nx)
	illum_stack=SharedMatrix{Float64}(zeros(nzd, nxd))

	itsnaps = [argmin(abs.(tgridmod .- tsnaps[i])) for i in 1:length(tsnaps)]

	nrmat=[ageom[ipw][iss].nr for ipw in 1:npw, iss in 1:nss]
	datamat=SharedArray{Float64}(nt,maximum(nrmat),nss)
	data=[Data(tgridmod,ageom[ip],rfields) for ip in 1:length(findall(!iszero, rflags))]

	# default is all prpagating wavefields are active
	activepw=[ipw for ipw in 1:npw]
	pac=PFdtdc(jobname,attrib_mod,npw,activepw,
	    exmodel,model,
	    ageom,srcwav,abs_trbl,sfields,isfields,sflags,
	    irfields,rflags,δt,δxI,δzI,
            nx,nz,nt,δtI,δx24I,δz24I,a_x,b_x,k_xI,a_x_half,b_x_half,k_x_halfI,a_z,b_z,k_zI,a_z_half,b_z_half,k_z_halfI,
	    modtt, modttI,modrr,modrrvx,modrrvz,
	    δmodtt,δmodrr,δmodrrvx,δmodrrvz,
	    δmod,
	    gradient,
	    grad_modtt_stack,
	    grad_modrrvx_stack,
	    grad_modrrvz_stack,grad_modrr_stack,
	    illum_flag,illum_stack,
	    backprop_flag,
	    snaps_flag,
	    itsnaps,
	    gmodel_flag,
	    ibx0,ibz0,ibx1,ibz1,
	    isx0,isz0,
	    datamat,
	    data,
	    verbose)	

	# dividing the supersources to workers
	if(nworker===nothing)
		nworker = min(nss, Distributed.nworkers())
	end
	work = Distributed.workers()[1:nworker]
	ssi=[round(Int, s) for s in range(0,stop=nss,length=nworker+1)]
	sschunks=Array{UnitRange{Int64}}(undef, nworker)
	for ib in 1:nworker       
		sschunks[ib]=ssi[ib]+1:ssi[ib+1]
	end

	# a distributed array of PFdtdp --- note that the parameters for each super source are efficiently distributed here
	papa=ddata(T=PFdtdp, init=I->PFdtdp(sschunks[I...][1],pac), pids=work);

	return PFdtd(sschunks, papa, pac)
end

"""
Create modeling parameters for each worker.
Each worker performs the modeling of supersources in `sschunks`.
The parameters common to all workers are stored in `pac`.
"""
function PFdtdp(sschunks::UnitRange{Int64},pac::PFdtdc)
	nx=pac.nx; nz=pac.nz; npw=pac.npw

	born_svalue_stack = zeros(nz, nx)

	p=[zeros(nz,nx,3) for ipw in 1:npw]; 
	dpdx=[zeros(nz,nx,3) for ipw in 1:npw]; 
	dpdz=[zeros(nz,nx,3) for ipw in 1:npw]; 
	pp=[zeros(nz,nx,3) for ipw in 1:npw]; 
	ppp=[zeros(nz,nx,3) for ipw in 1:npw]; 

	memory_dvx_dx=[zeros(nz,nx) for ipw in 1:npw]
	memory_dvx_dz=[zeros(nz,nx) for ipw in 1:npw]
	memory_dvz_dx=[zeros(nz,nx) for ipw in 1:npw]
	memory_dvz_dz=[zeros(nz,nx) for ipw in 1:npw]
	memory_dp_dx=[zeros(nz,nx) for ipw in 1:npw]
	memory_dp_dz=[zeros(nz,nx) for ipw in 1:npw]
	
	ss=[PFdtdss(iss, pac) for (issp,iss) in enumerate(sschunks)]

	pap=PFdtdp(ss,p,pp,ppp,dpdx,dpdz,memory_dp_dx,memory_dp_dz,memory_dvx_dx,memory_dvz_dz,
	    born_svalue_stack)

	return pap
end


"""
update the perturbation vector using the perturbed model
in this case, model will be treated as the background model 
* `δmod` is [δKI, δrhoI]
"""
function update_δmods!(pac::PFdtdc, δmod::Vector{Float64})
	nx=pac.nx; nz=pac.nz
	nznxd=prod(length.(pac.model.mgrid))
	copyto!(pac.δmod, δmod)
	fill!(pac.δmodtt,0.0)
	δmodtt=view(pac.δmodtt,npml+1:nz-npml,npml+1:nx-npml)
	for i in 1:nznxd
		# put perturbation due to KI as it is
		δmodtt[i] = δmod[i] 
	end
	fill!(pac.δmodrr,0.0)
	δmodrr=view(pac.δmodrr,npml+1:nz-npml,npml+1:nx-npml)
	for i in 1:nznxd
		# put perturbation due to rhoI here
		δmodrr[i] = δmod[nznxd+i]
	end
	# project δmodrr onto the vz and vx grids
	get_rhovxI!(pac.δmodrrvx, pac.δmodrr)
	get_rhovzI!(pac.δmodrrvz, pac.δmodrr)
	return nothing
end

"""
This method should be executed only after the updating the main model.
Update the `δmods` when a perturbed `model_pert` is input.
The model through which the waves are propagating 
is assumed to be the background model.
"""
function update_δmods!(pac::PFdtdc, model_pert::Medium)
	nznxd=prod(length.(pac.model.mgrid))
	fill!(pac.δmod,0.0)
	copyto!(pac.δmod, model_pert, [:KI, :rhoI])

	for i in 1:nznxd
		pac.δmod[i] -= pac.modtt[i] # subtracting the background model
		pac.δmod[nznxd+i] -= pac.modrr[i] # subtracting the background model
	end
	update_δmods!(pac, pac.δmod)
	return nothing
end

"""
```julia
update!(pa,medium_new)
```
Update `pa` with a new bundle of medium parameters `medium_new`, without additional memory allocation.
This routine is used during inversion, where medium parameters are iteratively updated. 
The ability to iteratively run the forward modeling task (with no additional memory allocation) on  
various subsurface models is necessary while implementing inversion 
algorithms.
"""
function update!(pa::PFdtd, model::Medium)
	return update!(pa.c, model)
end
function update!(pac::PFdtdc, model::Medium)

	copyto!(pac.model, model)

	Medium_pml_pad_trun!(pac.exmodel, pac.model, nlayer_rand)

	copyto!(pac.modttI, pac.exmodel, [:K]) 
	copyto!(pac.modtt, pac.exmodel, [:KI]) 
	copyto!(pac.modrr, pac.exmodel, [:rhoI])
	get_rhovxI!(pac.modrrvx, pac.modrr)
	get_rhovzI!(pac.modrrvz, pac.modrr)
	return nothing
end 

"""
```julia
update!(pa,srcwav_new,sflags)
```
Update `pa` with a new bundle of source wavelets `srcwav_new`, without additional memory allocation.
Optionally, `sflags` can be changed. 
"""
function update!(pa::PFdtd, srcwav::SrcWav, sflags::Any=nothing)
	update_srcwav!(pa,[srcwav],sflags)
end
function update_srcwav!(pa::PFdtd, srcwav::Vector{SrcWav}, sflags=nothing)
	# update srcwav in pa.c
	(length(srcwav) ≠ pa.c.npw) && error("cannot update")
	for i in 1:length(srcwav)
		copyto!(pa.c.srcwav[i], srcwav[i])
	end
	if(!(sflags===nothing))
		copyto!(pa.c.sflags, sflags)
	end

	# fill_wavelets for each supersource
	@sync begin
		for (ip, p) in enumerate(procs(pa.p))
			@async remotecall_wait(p) do 
				pap=localpart(pa.p)
				for is in 1:length(pap.ss)
					iss=pap.ss[is].iss
					wavelets=pap.ss[is].wavelets
					fill!.(wavelets, 0.0)
					fill_wavelets!(iss, wavelets, pa.c.srcwav, pa.c.sflags)
				end
			end
		end
	end
end

"""
Create modeling parameters for each supersource. 
Every worker models one or more supersources.
"""
function PFdtdss(iss::Int64, pac::PFdtdc)

	irfields=pac.irfields
	isfields=pac.isfields
	npw=pac.npw
	nt=pac.nt
	nx=pac.nx; nz=pac.nz
	nzd,nxd=length.(pac.model.mgrid)
	ageom=pac.ageom
	srcwav=pac.srcwav
	sflags=pac.sflags
	mesh_x, mesh_z = pac.exmodel.mgrid[2], pac.exmodel.mgrid[1]

	# records_output, distributed array among different procs
	records = [zeros(nt,pac.ageom[ipw][iss].nr,length(irfields)) for ipw in 1:npw]

	# gradient outputs
	grad_modtt = zeros(nz, nx)
	grad_modrrvx = zeros(nz, nx)
	grad_modrrvz = zeros(nz, nx)

	# saving illum
	illum =  (pac.illum_flag) ? zeros(nz, nx) : zeros(1,1)

	snaps = (pac.snaps_flag) ? zeros(nzd,nxd,length(pac.itsnaps)) : zeros(1,1,1)

	# source wavelets
	wavelets = [zeros(pac.ageom[ipw][iss].ns,length(isfields[ipw])) for ipw in 1:npw, it in 1:nt]
	fill_wavelets!(iss, wavelets, srcwav, sflags)

	

	# storing boundary values for back propagation
	nz1, nx1=length.(pac.model.mgrid)
	if(pac.backprop_flag ≠ 0)
		boundary=[zeros(3,nx1+6,nt),
		  zeros(nz1+6,3,nt),
		  zeros(3,nx1+6,nt),
		  zeros(nz1+6,3,nt),
		  zeros(nz1+2*npml,nx1+2*npml,3)]
	else
		boundary=[zeros(1,1,1) for ii in 1:5]
	end
	# source_spray_weights per supersource
	ssprayw = [zeros(4,ageom[ipw][iss].ns) for ipw in 1:npw]
	denomsI = [zeros(ageom[ipw][iss].ns) for ipw in 1:npw]
	isx1=[zeros(Int64,ageom[ipw][iss].ns) for ipw in 1:npw]
	isx2=[zeros(Int64,ageom[ipw][iss].ns) for ipw in 1:npw]
	isz1=[zeros(Int64,ageom[ipw][iss].ns) for ipw in 1:npw]
	isz2=[zeros(Int64,ageom[ipw][iss].ns) for ipw in 1:npw]
	for ipw=1:npw
		for is=1:ageom[ipw][iss].ns
			weights=ssprayw[ipw]
			denom=denomsI[ipw]
			Interpolation.get_spray_weights!(view(weights, :,is), view(denom,is), 
				    view(isx1[ipw],is), view(isx2[ipw],is),
				    view(isz1[ipw],is), view(isz2[ipw],is),
				    mesh_x, mesh_z, ageom[ipw][iss].sx[is], ageom[ipw][iss].sz[is])
		end
	end

	# receiver interpolation weights per sequential source
	rinterpolatew = [zeros(4,ageom[ipw][iss].nr) for ipw in 1:npw]
	denomrI = [zeros(ageom[ipw][iss].nr) for ipw in 1:npw]
	irx1=[zeros(Int64,ageom[ipw][iss].nr) for ipw in 1:npw]
	irx2=[zeros(Int64,ageom[ipw][iss].nr) for ipw in 1:npw]
	irz1=[zeros(Int64,ageom[ipw][iss].nr) for ipw in 1:npw]
	irz2=[zeros(Int64,ageom[ipw][iss].nr) for ipw in 1:npw]
	for ipw=1:npw
		for ir=1:ageom[ipw][iss].nr
			weights=rinterpolatew[ipw]
			denom=denomrI[ipw]
			Interpolation.get_interpolate_weights!(view(weights, :,ir), view(denom, ir), 
				  view(irx1[ipw],ir), view(irx2[ipw],ir),
				  view(irz1[ipw],ir), view(irz2[ipw],ir),
				  mesh_x, mesh_z, ageom[ipw][iss].rx[ir], ageom[ipw][iss].rz[ir])
		end
	end


	pass=PFdtdss(iss,wavelets,ssprayw,records,rinterpolatew,
	      isx1,isx2,isz1,isz2,irx1,irx2,irz1,irz2,boundary,snaps,illum, 
	      grad_modtt,grad_modrrvx,grad_modrrvz)


	return pass
end

include("source.jl")
include("receiver.jl")
include("core.jl")
include("rho_projection.jl")
include("gradient.jl")
include("born.jl")
include("boundary.jl")
include("gallery.jl")


"""
```julia
update!(pa)
```
In-place method to perform the experiment and update `pa` after wave propagation. After update, see
`pa[:data]` and `pa[:snaps]`.
"""
@fastmath function update!(pa::PFdtd)

	global to

	reset_timer!(to)

	@timeit to "initialize!" begin
		# zero out all the results stored in pa.c
		initialize!(pa.c) 
	
		# zero out results stored per worker
		@sync begin
			for (ip, p) in enumerate(procs(pa.p))
				@async remotecall_wait(p) do 
					initialize!(localpart(pa.p))
				end
			end
		end
	end


	@timeit to "mod_per_proc!" begin
		# all localparts of DArray are input to this method
		# parallelization over shots
		@sync begin
			for (ip, p) in enumerate(procs(pa.p))
				@async remotecall_wait(p) do 
					mod_per_proc!(pa.c, localpart(pa.p))
				end
			end
		end
	end

	@timeit to "stack_grads!" begin
		# stack gradients and illum over sources
		@sync begin
			for (ip, p) in enumerate(procs(pa.p))
				@sync remotecall_wait(p) do 
					(pa.c.gmodel_flag) && stack_grads!(pa.c, localpart(pa.p))
					(pa.c.illum_flag) && stack_illums!(pa.c, localpart(pa.p))
				end
			end
		end
	end

	@timeit to "update gradient" begin
		# update gradient model using grad_modtt_stack, grad_modrr_stack
		(pa.c.gmodel_flag) && update_gradient!(pa.c)
	end


	@timeit to "record data" begin
	for ipw in pa.c.activepw
		if(pa.c.rflags[ipw] ≠ 0) # record only if rflags is non-zero
			for ifield in 1:length(pa.c.irfields)

				fill!(pa.c.datamat, 0.0)
				@sync begin
					for (ip, p) in enumerate(procs(pa.p))
						@sync remotecall_wait(p) do 
							update_datamat!(ifield, ipw, pa.c, localpart(pa.p))
						end
					end
				end
				update_data!(ifield, ipw, pa.c)
			end
		end
	end
	end
	if(pa.c.verbose)
		show(to)	
		println("  ")
	end
	return nothing
end

function update_datamat!(ifield, ipw, pac::PFdtdc, pap::PFdtdp)
	datamat=pac.datamat
        pass=pap.ss
	for issp in 1:length(pass)
		iss=pass[issp].iss
		records=pass[issp].records[ipw]
		for ir in 1:pac.ageom[ipw][iss].nr
			for it in 1:pac.nt
				datamat[it,ir,iss]=records[it,ir,ifield]
			end
		end
        end
end

function update_data!(ifield, ipw, pac::PFdtdc)
	datamat=pac.datamat
	for iss in 1:length(pac.ageom[1])
		data=pac.data[ipw][iss].d[ifield]
		for ir in 1:pac.ageom[ipw][iss].nr
			for it in 1:pac.nt
				data[it,ir]=datamat[it,ir,iss]
			end
		end
	end
end


# update TDout after forming a vector and resampling
#	ipropout=0;
#	for iprop in 1:pac.npw
#		if(pac.rflags[iprop] ≠ 0)
#			ipropout += 1
##			Data.TD_resamp!(pac.data[ipropout], Data.TD_urpos((Array(records[:,:,iprop,:,:])), rfields, tgridmod, ageom[iprop],
##				ageom_urpos[1].nr[1],
##				(ageom_urpos[1].rz[1], ageom_urpos[1].rx[1])
##				)) 
#		end
#	end
	# return without resampling for testing
	#return [Data.TD(reshape(records[1+(iprop-1)*nd : iprop*nd],tgridmod.nx,recv_n,nss),
	#		       tgridmod, ageom[1]) for iprop in 1:npw]

function stack_illums!(pac::PFdtdc, pap::PFdtdp)
	nx, nz=pac.nx, pac.nz
	illums=pac.illum_stack
	pass=pap.ss
	for issp in 1:length(pass)
		gs=pass[issp].illum
		gss=view(gs,npml+1:nz-npml,npml+1:nx-npml)
		@. illums += gss
	end
end



# modelling for each processor
function mod_per_proc!(pac::PFdtdc, pap::PFdtdp) 
	# source_loop
	for issp in 1:length(pap.ss)
		reset_per_ss!(pap)

		iss=pap.ss[issp].iss

		if(pac.backprop_flag==-1)
			"initial conditions from boundary for first propagating field only"
			boundary_force_snap_p!(issp,pac,pap.ss,pap)
			boundary_force_snap_vxvz!(issp,pac,pap.ss,pap)
		end

		prog = Progress(pac.nt, dt=1, desc="\tmodeling supershot $iss/$(length(pac.ageom[1])) ", 
		  		color=:white) 
		# time_loop
		"""
		* don't use shared arrays inside this time loop, for speed when using multiple procs
		"""
		for it=1:pac.nt

			pac.verbose && next!(prog, :white)

			advance!(pac,pap)
		
			# force p[1] on boundaries
			(pac.backprop_flag==-1) && boundary_force!(it,issp,pac,pap.ss,pap)
	 
			add_source!(it, issp, iss, pac, pap.ss, pap, Source_B1())

			# no born flag for adjoint modelling
			if(!pac.gmodel_flag)
				(typeof(pac.attrib_mod)==FdtdBorn) && add_born_sources!(issp, pac, pap.ss, pap)
			end

			# record boundaries after time reversal already
			(pac.backprop_flag==1) && boundary_save!(pac.nt-it+1,issp,pac,pap.ss,pap)

			record!(it, issp, iss, pac, pap.ss, pap, Receiver_B1())

			(pac.gmodel_flag) && compute_gradient!(issp, pac, pap.ss, pap)

			(pac.illum_flag) && compute_illum!(issp, pap.ss, pap)

			if(pac.snaps_flag)
				iitsnaps=findall(x->==(x,it),pac.itsnaps)
				for itsnap in iitsnaps
					snaps_save!(itsnap,issp,pac,pap.ss,pap)
				end
			end

		end # time_loop
		"now pressure is at [nt], velocities are at [nt-1/2]"	

		"one more propagating step to save pressure at [nt+1] -- for time revarsal"
		advance!(pac,pap)

		"save last snap of pressure field"
		(pac.backprop_flag==1) && boundary_save_snap_p!(issp,pac,pap.ss,pap)

		"one more propagating step to save velocities at [nt+3/2] -- for time reversal"
		advance!(pac,pap)

		"save last snap of velocity fields with opposite sign for adjoint propagation"
		(pac.backprop_flag==1) && boundary_save_snap_vxvz!(issp,pac,pap.ss,pap)

		"scale gradients for each issp"
		(pac.gmodel_flag) && scale_gradient!(issp, pap.ss, step(pac.model.mgrid[2])*step(pac.model.mgrid[1]))
		

	end # source_loop
end # mod_per_shot




# Need illumination to estimate the approximate diagonal of Hessian
@inbounds @fastmath function compute_illum!(issp::Int64, pass::Vector{PFdtdss}, pap::PFdtdp)
	# saving illumination to be used as preconditioner 
	p=pap.p
	illum=pass[issp].illum
	pp=view(p[1],:,:,1)
	for i in eachindex(illum)
		illum[i] += abs2(pp[i])
	end
end


@fastmath @inbounds function snaps_save!(itsnap::Int64,issp::Int64,pac::PFdtdc,pass::Vector{PFdtdss},pap::PFdtdp)
	isx0=pac.isx0
	isz0=pac.isz0
	p=pap.p
	snaps=pass[issp].snaps
	for ix=1:size(snaps,2)
		@simd for iz=1:size(snaps,1)
			snaps[iz,ix,itsnap]=p[1][isz0+iz,isx0+ix,1]
		end
	end
end






function check_fd_stability(vpmin::Float64, vpmax::Float64, δx::Float64, δz::Float64, freqmin::Float64, freqmax::Float64, δt::Float64, verbose::Bool)


	# check spatial sampling
	δs_temp=round(vpmin/5.0/freqmax,digits=2);
	δs_max = maximum([δx, δz])
	str1=@sprintf("%0.2e",δs_max)
	str2=@sprintf("%0.2e",δs_temp)
	if(str1 ≠ str2)
		if(all(δs_max .> δs_temp)) 
			@warn "decrease spatial sampling ($str1) below $str2"
		else
			@info "spatial sampling ($str1) can be as high as $str2"
		end
	end

	# check time sampling
	δs_min = minimum([δx, δz])
	δt_temp=0.5*δs_min/vpmax
	str1=@sprintf("%0.2e", δt)
	str2=@sprintf("%0.2e", δt_temp)
	if(str1 ≠ str2)
		if(all(δt .> δt_temp))
			@warn "decrease time sampling ($str1) below $str2"
		else
			@info "time sampling ($str1) can be as high as $str2"
		end
	end

end

"""
Output vectors related to PML boundaries.
"""
function pml_variables(
		nx::Int64, 
		δt::Float64, 
		δx::Float64, mesh_na_pml::Int64, 
		velmax::Float64, velmin::Float64, 
		freqmin::Float64, freqmax::Float64,
	        flags::Vector{Bool}
		)

	NPOWER = 2.e0
	K_MAX_PML = 1.e0 # from Gedney page 8.11
	ALPHA_MAX_PML = 2.e0 * pi * ((freqmin + freqmax))/4.e0 # from Festa and Vilotte

	# thickness of the PML layer in meters
	thickness_PML_x = mesh_na_pml * δx

	"reflection coefficient (INRIA report section 6.1) http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf"
	Rcoef = 0.001e0

	#! check that NPOWER is okay
	#if(NPOWER < 1) stop "NPOWER must be greater than 1"

	# compute d0 from INRIA report section 6.1 http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
	d0_x = - (NPOWER + 1) * (velmax+velmin)/2.0 * log(Rcoef) / (2.e0 * thickness_PML_x)


	d_x=zeros(nx);	d_x_half=zeros(nx)
	k_x=ones(nx);	k_xI=ones(nx)
	k_x_half=ones(nx);	k_x_halfI=ones(nx)
	alpha_x=zeros(nx);	alpha_x_half=zeros(nx)
	a_x=zeros(nx);	a_x_half=zeros(nx)
	b_x=zeros(nx);	b_x_half=zeros(nx)


	"damping in the X direction"
	"origin of the PML layer (position of right edge minus thickness, in meters)"
	xoriginleft = thickness_PML_x
	xoriginright = (nx-1)*δx - thickness_PML_x

	for ix=1:nx
		# abscissa of current grid point along the damping profile
		xval = δx * real(ix-1)

		#---------- left edge
		if(flags[1]) 

			# define damping profile at the grid points
			abscissa_in_PML = xoriginleft - xval
			if(abscissa_in_PML >= 0.0) 
				abscissa_normalized = abscissa_in_PML / thickness_PML_x
				d_x[ix] = d0_x * abscissa_normalized.^NPOWER
				# this taken from Gedney page 8.2
				k_x[ix] = 1.e0 + (K_MAX_PML - 1.e0) * abscissa_normalized.^NPOWER
				alpha_x[ix] = ALPHA_MAX_PML * (1.e0 - abscissa_normalized) + 0.1e0 * ALPHA_MAX_PML
			end

			# define damping profile at half the grid points
			abscissa_in_PML = xoriginleft - (xval + δx/2.e0)
			if(abscissa_in_PML >= 0.0) 
				abscissa_normalized = abscissa_in_PML / thickness_PML_x
				d_x_half[ix] = d0_x * abscissa_normalized.^NPOWER
				# this taken from Gedney page 8.2
				k_x_half[ix] = 1.e0 + (K_MAX_PML - 1.e0) * abscissa_normalized.^NPOWER
				alpha_x_half[ix] = ALPHA_MAX_PML * (1.e0 - abscissa_normalized) + 0.1e0 * ALPHA_MAX_PML
			end
		end

		#---------- right edge
		if(flags[2]) 

			# define damping profile at the grid points
			abscissa_in_PML = xval - xoriginright
			if(abscissa_in_PML >= 0.0) 
				abscissa_normalized = abscissa_in_PML / thickness_PML_x
				d_x[ix] = d0_x * abscissa_normalized^NPOWER
				# this taken from Gedney page 8.2
				k_x[ix] = 1.e0 + (K_MAX_PML - 1.e0) * abscissa_normalized^NPOWER
				alpha_x[ix] = ALPHA_MAX_PML * (1.e0 - abscissa_normalized) + 0.1e0 * ALPHA_MAX_PML
			end

			# define damping profile at half the grid points
			abscissa_in_PML = xval + δx/2.e0 - xoriginright
			if(abscissa_in_PML >= 0.0) 
				abscissa_normalized = abscissa_in_PML / thickness_PML_x
				d_x_half[ix] = d0_x * abscissa_normalized^NPOWER
				# this taken from Gedney page 8.2
				k_x_half[ix] = 1.e0 + (K_MAX_PML - 1.e0) * abscissa_normalized^NPOWER
				alpha_x_half[ix] = ALPHA_MAX_PML * (1.e0 - abscissa_normalized) + 0.1e0 * ALPHA_MAX_PML
			end

		end
#
		# just in case, for -5 at the end
		(alpha_x[ix] < 0.0) ? alpha_x[ix] = 0.0 : nothing
		(alpha_x_half[ix] < 0.0) ? alpha_x_half[ix] = 0.0 : nothing

		b_x[ix] = exp(- (d_x[ix] / k_x[ix] + alpha_x[ix]) * δt)
		b_x_half[ix] = exp(- (d_x_half[ix] / k_x_half[ix] + alpha_x_half[ix]) * δt)

		# this to avoid division by zero outside the PML
		(abs(d_x[ix]) > 1.e-6) ? a_x[ix] = d_x[ix] * (b_x[ix] - 1.e0) / (k_x[ix] * (d_x[ix] + k_x[ix] * alpha_x[ix])) : nothing
		(abs(d_x_half[ix]) > 1.e-6) ? a_x_half[ix] = d_x_half[ix] * (b_x_half[ix] - 1.e0) / (k_x_half[ix] * (d_x_half[ix] + k_x_half[ix] * alpha_x_half[ix])) : nothing

	end
	k_xI = k_x.^(-1.e0)
	k_x_halfI = k_x_half.^(-1.e0)


	return a_x, b_x, k_xI, a_x_half, b_x_half, k_x_halfI
end


include("getprop.jl")

