module Fdtd

using Grid
using Interpolation
using Signals
import JuMIT.Models
import JuMIT.Acquisition
import JuMIT.Data
using ProgressMeter
using TimerOutputs
using Distributed
using DistributedArrays
using SharedArrays
using LinearAlgebra
using AxisArrays
using Printf

global const to = TimerOutput(); # create a timer object
global const npml = 50

#As forward modeling method, the 
#finite-difference method is employed. 
#It uses a discrete version of the two-dimensional isotropic acoustic wave equation.
#As shown in  
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

* `gradient::Vector{Float64}=Models.Seismic_zeros(model.mgrid)` : gradient model modified only if `gmodel_flag`
* `TDout::Vector{Data.TD}=[Data.TD_zeros(rfields,tgridmod,acqgeom[ip]) for ip in 1:length(findn(rflags))]`
* `illum::Array{Float64,2}=zeros(model.mgrid.nz, model.mgrid.nx)` : source energy if `illum_flag`
* `boundary::Array{Array{Float64,4},1}` : stored boundary values for first propagating wavefield 
* `snaps::Array{Float64,4}=zeros(model.mgrid.nz,model.mgrid.nx,length(tsnaps),acqgeom[1].nss)` :snapshots saved at `tsnaps`

# Return (in order)

* modelled data for each propagating wavefield as `Vector{TD}`
* stored boundary values of the first propagating wavefield as `Array{Array{Float64,4},1}` (use for backpropagation)
* final conditions of the first propagating wavefield as `Array{Float64,4}` (use for back propagation)
* gradient model as `Seismic`
* stored snaps shots at tsnaps as Array{Float64,4} 


"""
mutable struct Paramc
	jobname::Symbol
	npw::Int64 
	activepw::Vector{Int64}
	exmodel::Models.Seismic
	exmodel_pert::Models.Seismic
	model::Models.Seismic
	model_pert::Models.Seismic
	acqgeom::Vector{Acquisition.Geom}
	acqsrc::Vector{Acquisition.Src}
	abs_trbl::Vector{Symbol}
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
	modttI::Matrix{Float64}
	modrr::Matrix{Float64}
	modrrvx::Matrix{Float64}
	modrrvz::Matrix{Float64}
	δmodtt::Matrix{Float64}
	modrr_pert::Matrix{Float64}
	δmodrrvx::Matrix{Float64}
	δmodrrvz::Matrix{Float64}
	gradient::Vector{Float64}  # output gradient vector
	grad_modtt_stack::SharedArrays.SharedArray{Float64,2} # contains gmodtt
	grad_modrrvx_stack::SharedArrays.SharedArray{Float64,2}
	grad_modrrvz_stack::SharedArrays.SharedArray{Float64,2}
	grad_modrr_stack::SharedArrays.SharedArray{Float64,2}
	illum_flag::Bool
	illum_stack::SharedArrays.SharedArray{Float64,2}
	backprop_flag::Int64
	snaps_flag::Bool
	itsnaps::Vector{Int64}
	born_flag::Bool
	gmodel_flag::Bool
	ibx0::Int64
	ibz0::Int64
	ibx1::Int64
	ibz1::Int64
	isx0::Int64
	isz0::Int64
	datamat::SharedArrays.SharedArray{Float64,3}
	data::Vector{Data.TD}
	verbose::Bool
end 

"""
Modelling parameters per every supersource for each worker
"""
mutable struct Paramss
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
Parameters per every worker, not necessarily for every supersource.
Note that a single worker can take care of multiple supersources.
"""
mutable struct Paramp
	ss::Vector{Paramss}
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

mutable struct Param
	p::DistributedArrays.DArray{Paramp,1,Paramp} # distributed parameters among workers
	c::Paramc # common parameters
end

Base.print(pa::Param)=nothing
Base.show(pa::Param)=nothing

function initialize!(pap::Paramp)
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

function initialize!(pac::Paramc)
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

function reset_per_ss!(pap::Paramp)
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
Method to create `Fdtd` modeling parameters.
The output of this method can be used as an input to `mod!`, where the actual 
finite-difference modeling is performed.

# Keyword Arguments

* `npw::Int64=1` : number of independently propagating wavefields in `model`
* `model::Models.Seismic` : seismic medium parameters 
* `model_pert::Models.Seismic=model` : perturbed model, i.e., model + δmodel, used only for Born modeling 
* `tgridmod::Grid.M1D=` : modeling time grid, maximum time in tgridmod should be greater than or equal to maximum source time, same sampling interval as the wavelet
* `tgrid::Grid.M1D=tgridmod` : output records are resampled on this time grid
* `acqgeom::Vector{Acquisition.Geom}` :  acquisition geometry for each independently propagating wavefield
* `acqsrc::Vector{Acquisition.Src}` : source acquisition parameters for each independently propagating wavefield
* `sflags::Vector{Int64}=fill(2,npw)` : source related flags for each propagating wavefield
  * `=[0]` inactive sources
  * `=[1]` sources with injection rate
  * `=[2]` volume injection sources
  * `=[3]` sources input after time reversal (use only during backpropagation)
* `rflags::Vector{Int64}=fill(1,npw)` : receiver related flags for each propagating wavefield
  * `=[0]` receivers do not record (or) inactive receivers
  * `=[0,1]` receivers are active only for the second propagating wavefield
* `rfields::Vector{Symbol}=[:P]` : multi-component receiver flag; types fields the receivers record (to be changed later)
* `backprop_flag::Bool=Int64` : save final state variables and the boundary conditions for later use
  * `=1` save boundary and final values in `boundary` 
  * `=-1` use stored values in `boundary` for back propagation
* `abs_trbl::Vector{Symbol}=[:top, :bottom, :right, :left]` : use absorbing PML boundary conditions or not
  * `=[:top, :bottom]` apply PML conditions only at the top and bottom of the model 
  * `=[:bottom, :right, :left]` top is reflecting
* `born_flag::Bool=false` : do only Born modeling instead of full wavefield modelling (to be updated soon)
* `gmodel_flag=false` : flag that is used to output gradient; there should be atleast two propagating wavefields in order to do so: 1) forward wavefield and 2) adjoint wavefield
* `illum_flag::Bool=false` : flag to output wavefield energy or source illumination; it can be used as preconditioner during inversion
* `tsnaps::Vector{Float64}=fill(0.5*(tgridmod.x[end]+tgridmod.x[1]),1)` : store snaps at these modelling times
* `snaps_flag::Bool=false` : return snaps or not
* `verbose::Bool=false` : verbose flag

# Example

```julia
pa = JuMIT.Fdtd.Param(acqgeom=acqgeom, acqsrc=acqsrc, model=model, tgridmod=tgridmod);
JuMIT.Fdtd.mod!(pa);
```
# Credits 

Author: Pawan Bharadwaj 
        (bharadwaj.pawan@gmail.com)

* original code in FORTRAN90: March 2013
* modified: 11 Sept 2013
* major update: 25 July 2014
* code optimization with help from Jan Thorbecke: Dec 2015
* rewritten in Julia: June 2017
* added parrallelization over supersources in Julia: July 2017
* efficient parrallelization using distributed arrays: Sept 2017
* optimized memory allocation: Oct 2017
"""
function Param(;
	jobname::Symbol=:forward_propagation,
	npw::Int64=1, 
	model::Models.Seismic=nothing,
	abs_trbl::Vector{Symbol}=[:top, :bottom, :right, :left],
	born_flag::Bool=false,
	model_pert::Models.Seismic = model,
	tgridmod::Grid.M1D=nothing,
	acqgeom::Vector{Acquisition.Geom}=nothing,
	acqsrc::Array{Acquisition.Src}=nothing,
	sflags::Vector{Int64}=fill(2,npw), 
	rflags::Vector{Int64}=fill(1,npw),
	rfields::Vector{Symbol}=[:P], 
	backprop_flag::Int64=0,  
	gmodel_flag::Bool=false,
	illum_flag::Bool=false,
	tsnaps::Vector{Float64}=fill(0.5*(tgridmod.x[end]+tgridmod.x[1]),1),
	snaps_flag::Bool=false,
	verbose::Bool=false,
	nworker=nothing)

	#println("********PML Removed*************")
	#abs_trbl=[:null]
	

	# check sizes and errors based on input
	#(length(TDout) ≠ length(findn(rflags))) && error("TDout dimension")
	(length(acqgeom) ≠ npw) && error("acqgeom dimension")
	(length(acqsrc) ≠ npw) && error("acqsrc dimension")
	(length(sflags) ≠ npw) && error("sflags dimension")
	(length(rflags) ≠ npw) && error("rflags dimension")
	(maximum(tgridmod.x) < maximum(acqsrc[1].tgrid.x)) && error("modeling time is less than source time")
	#(any([getfield(TDout[ip],:tgrid).δx < tgridmod.δx for ip=1:length(TDout)])) && error("output time grid sampling finer than modeling")
	#any([maximum(getfield(TDout[ip],:tgrid).x) > maximum(tgridmod.x) for ip=1:length(TDout)]) && error("output time > modeling time")

	#! no modeling if source wavelet is zero
	#if(maxval(abs(src_wavelets)) .lt. tiny(rzero_de)) 
	#        return
	#endif

	# all the propagating wavefields should have same supersources? check that?

	# check dimension of model and model_pert

	# check if all sources are receivers are inside model
	any(.![Acquisition.Geom_check(acqgeom[ip], model.mgrid) for ip=1:npw]) ? error("sources or receivers not inside model") : nothing


	length(acqgeom) != npw ? error("acqgeom size") : nothing
	length(acqsrc) != npw ? error("acqsrc size") : nothing
	any([getfield(acqgeom[ip],:nss) != getfield(acqsrc[ip],:nss) for ip=1:npw])  ? error("different supersources") : nothing
	any([getfield(acqgeom[ip],:ns) != getfield(acqsrc[ip],:ns) for ip=1:npw])  ? error("different sources") : nothing

	# necessary that nss and fields should be same for all nprop
	nss = acqgeom[1].nss;
	sfields = [acqsrc[ipw].fields for ipw in 1:npw]
	isfields = [findall(in([:P, :Vx, :Vz]), acqsrc[ipw].fields) for ipw in 1:npw]
	fill(nss, npw) != [getfield(acqgeom[ip],:nss) for ip=1:npw] ? error("different supersources") : nothing

	# create acquisition geometry with each source shooting 
	# at every unique receiver position
	irfields = findall(in([:P, :Vx, :Vz]), rfields)


	#(verbose) &&	println(string("\t> number of super sources:\t",nss))	

	# find maximum and minimum frequencies in the source wavelets
	freqmin = Signals.DSP.findfreq(acqsrc[1].wav[1,1][:,1],acqsrc[1].tgrid,attrib=:min) 
	freqmax = Signals.DSP.findfreq(acqsrc[1].wav[1,1][:,1],acqsrc[1].tgrid,attrib=:max) 

	# minimum and maximum velocities
	vpmin = minimum(broadcast(minimum,[model.vp0, model_pert.vp0]))
	vpmax = maximum(broadcast(maximum,[model.vp0, model_pert.vp0]))
	#verbose && println("\t> minimum and maximum velocities:\t",vpmin,"\t",vpmax)


	check_fd_stability(vpmin, vpmax, model.mgrid.δx, model.mgrid.δz, freqmin, freqmax, tgridmod.δx, verbose)


	# where to store the boundary values (careful, born scaterrers cannot be inside these!?)
	ibx0=npml-2; ibx1=model.mgrid.nx+npml+3
	ibz0=npml-2; ibz1=model.mgrid.nz+npml+3
#	ibx0=npml+1; ibx1=model.mgrid.nx+npml
#	ibz0=npml+1; ibz1=model.mgrid.nz+npml
#	println("**** Boundary Storage Changed **** ")

	# for snaps
	isx0, isz0=npml, npml

	# extend models in the PML layers
	exmodel = Models.Seismic_pml_pad_trun(model);
	exmodel_pert = Models.Seismic_pml_pad_trun(model_pert);


	"density values on vx and vz stagerred grids"
	modrr=Models.Seismic_get(exmodel, :ρI)
	modrrvx = get_rhovxI(modrr)
	modrrvz = get_rhovzI(modrr)

	modttI = Models.Seismic_get(exmodel, :K) 
	modtt = Models.Seismic_get(exmodel, :KI)

	if(born_flag)
		(npw≠2) && error("born_flag needs npw=2")
	end

	#>>>>>>>>>>>>>### Dont need this, if born inversion is not performed ####
	isapprox(exmodel, exmodel_pert) || error("pertrubed model should be similar to background model")
	"inverse density contrasts for Born Modelling"
	modrr_pert = Models.Seismic_get(exmodel_pert, :ρI)
	δmodrrvx = get_rhovxI(modrr_pert)
	δmodrrvz = get_rhovzI(modrr_pert)
	for i in eachindex(δmodrrvx)
		δmodrrvx[i] -= modrrvx[i]
		δmodrrvz[i] -= modrrvz[i]
	end

	"inverse Bulk Modulus contrasts for Born Modelling"
	δmodtt = Models.Seismic_get(exmodel_pert, :KI)
	for i in eachindex(δmodtt)
		δmodtt[i] -= modtt[i]
	end
	#<<<<<<<<<<<<<<<<##############################################################



	#create some aliases
	nx, nz = exmodel.mgrid.nx, exmodel.mgrid.nz
	nxd, nzd = model.mgrid.nx, model.mgrid.nz
	δx, δz = exmodel.mgrid.δx, exmodel.mgrid.δz
	δt = tgridmod.δx 
	nt=tgridmod.nx

	δx24I, δz24I = (δx * 24.0)^(-1.0), (δz * 24.0)^(-1.0)
	δxI, δzI = (δx)^(-1.0), (δz)^(-1.0)
	δt = tgridmod.δx 
	δtI = (δt)^(-1.0)
	nt=tgridmod.nx
	mesh_x, mesh_z = exmodel.mgrid.x, exmodel.mgrid.z


	# pml_variables
	a_x, b_x, k_xI, a_x_half, b_x_half, k_x_halfI = pml_variables(nx, δt, δx, npml-5, vpmax, vpmin, freqmin, freqmax, 
							       [any(abs_trbl .== :left), any(abs_trbl .== :right)])
	a_z, b_z, k_zI, a_z_half, b_z_half, k_z_halfI = pml_variables(nz, δt, δz, npml-5, vpmax, vpmin, freqmin, freqmax,
							       [any(abs_trbl .== :top), any(abs_trbl .== :bottom)])

	gradient=zeros(2*nzd*nxd)
	grad_modtt_stack=SharedMatrix{Float64}(zeros(nzd,nxd))
	grad_modrrvx_stack=SharedMatrix{Float64}(zeros(nzd,nxd))
	grad_modrrvz_stack=SharedMatrix{Float64}(zeros(nzd,nxd))
	grad_modrr_stack=SharedMatrix{Float64}(zeros(nzd,nxd))
	illum_stack=SharedMatrix{Float64}(zeros(nzd, nxd))

	itsnaps = [argmin(abs.(tgridmod.x-tsnaps[i])) for i in 1:length(tsnaps)]

	nrmat=[acqgeom[ipw].nr[iss] for ipw in 1:npw, iss in 1:acqgeom[1].nss]
	datamat=SharedArray{Float64}(nt,maximum(nrmat),acqgeom[1].nss)
	data=[Data.TD_zeros(rfields,tgridmod,acqgeom[ip]) for ip in 1:length(findall(!iszero, rflags))]

	# default is all prpagating wavefields are active
	activepw=[ipw for ipw in 1:npw]
	pac=Paramc(jobname,npw,activepw,
	    exmodel,exmodel_pert,model,model_pert,
	    acqgeom,acqsrc,abs_trbl,isfields,sflags,
	    irfields,rflags,δt,δxI,δzI,
            nx,nz,nt,δtI,δx24I,δz24I,a_x,b_x,k_xI,a_x_half,b_x_half,k_x_halfI,a_z,b_z,k_zI,a_z_half,b_z_half,k_z_halfI,
	    modtt, modttI,modrr,modrrvx,modrrvz,
	    δmodtt,modrr_pert,δmodrrvx,δmodrrvz,
	    gradient,
	    grad_modtt_stack,
	    grad_modrrvx_stack,
	    grad_modrrvz_stack,grad_modrr_stack,
	    illum_flag,illum_stack,
	    backprop_flag,
	    snaps_flag,
	    itsnaps,born_flag,
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

	# a distributed array of Paramp --- note that the parameters for each super source are efficiently distributed here
	papa=ddata(T=Paramp, init=I->Paramp(sschunks[I...][1],pac), pids=work);

	return Param(papa, pac)
end


"""
Create modeling parameters for each worker.
Each worker performs the modeling of supersources in `sschunks`.
The parameters common to all workers are stored in `pac`.
"""
function Paramp(sschunks::UnitRange{Int64},pac::Paramc)
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
	
	ss=[Paramss(iss, pac) for (issp,iss) in enumerate(sschunks)]

	pap=Paramp(ss,p,pp,ppp,dpdx,dpdz,memory_dp_dx,memory_dp_dz,memory_dvx_dx,memory_dvz_dz,
	    born_svalue_stack)

	return pap
end

"""
Update the `Seismic` models in `Paramc` without additional memory allocation.
This routine is used during FWI, where medium parameters are itertively updated. 
"""
function update_model!(pac::Paramc, model::Models.Seismic, model_pert=nothing)

	copyto!(pac.model, model)

	Models.Seismic_pml_pad_trun!(pac.exmodel, pac.model)

	Models.Seismic_get!(pac.modttI, pac.exmodel, [:K]) 
	Models.Seismic_get!(pac.modtt, pac.exmodel, [:KI]) 

	Models.Seismic_get!(pac.modrr, pac.exmodel, [:ρI])
	get_rhovxI!(pac.modrrvx, pac.modrr)
	get_rhovzI!(pac.modrrvz, pac.modrr)


	if(!(model_pert === nothing))
		copyto!(pac.model_pert, model_pert)
		Models.Seismic_pml_pad_trun!(pac.exmodel_pert, pac.model_pert);
		Models.Seismic_get!(pac.modrr_pert, pac.exmodel_pert, [:ρI])
		get_rhovxI!(pac.δmodrrvx, pac.modrr_pert)
		get_rhovzI!(pac.δmodrrvz, pac.modrr_pert)
		# subtract background model to get the constrasts
		for i in eachindex(pac.δmodrrvx)
			pac.δmodrrvx[i] -= pac.modrrvx[i]
			pac.δmodrrvz[i] -= pac.modrrvz[i]
		end

		Models.Seismic_get!(pac.δmodtt, pac.exmodel_pert, [:KI])

		# subtract background model to get the constrasts
		for i in eachindex(pac.δmodtt)
			pac.δmodtt[i] -= pac.modtt[i]
		end
	end
end 

function update_acqsrc!(pa::Param, acqsrc::Vector{Acquisition.Src}, sflags=nothing)
	# update acqsrc in pa.c
	(length(acqsrc) ≠ pa.c.npw) && error("cannot update")
	for i in 1:length(acqsrc)
		copyto!(pa.c.acqsrc[i], acqsrc[i])
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
					fill_wavelets!(iss, wavelets, pa.c.acqsrc, pa.c.sflags)
				end
			end
		end
	end
end

"""
Create modeling parameters for each supersource. 
Every worker models one or more supersources.
"""
function Paramss(iss::Int64, pac::Paramc)

	irfields=pac.irfields
	isfields=pac.isfields
	npw=pac.npw
	nt=pac.nt
	nx=pac.nx; nz=pac.nz
	nxd=pac.model.mgrid.nx
	nzd=pac.model.mgrid.nz
	acqgeom=pac.acqgeom
	acqsrc=pac.acqsrc
	sflags=pac.sflags
	mesh_x, mesh_z = pac.exmodel.mgrid.x, pac.exmodel.mgrid.z

	# records_output, distributed array among different procs
	records = [zeros(nt,pac.acqgeom[ipw].nr[iss],length(irfields)) for ipw in 1:npw]

	# gradient outputs
	grad_modtt = zeros(nz, nx)
	grad_modrrvx = zeros(nz, nx)
	grad_modrrvz = zeros(nz, nx)

	# saving illum
	illum =  (pac.illum_flag) ? zeros(nz, nx) : zeros(1,1)

	snaps = (pac.snaps_flag) ? zeros(nzd,nxd,length(pac.itsnaps)) : zeros(1,1,1)

	# source wavelets
	wavelets = [zeros(pac.acqgeom[ipw].ns[iss],length(isfields[ipw])) for ipw in 1:npw, it in 1:nt]
	fill_wavelets!(iss, wavelets, acqsrc, sflags)

	

	# storing boundary values for back propagation
	nx1, nz1=pac.model.mgrid.nx, pac.model.mgrid.nz
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
	ssprayw = [zeros(4,acqgeom[ipw].ns[iss]) for ipw in 1:npw]
	denomsI = [zeros(acqgeom[ipw].ns[iss]) for ipw in 1:npw]
	isx1=[zeros(Int64,acqgeom[ipw].ns[iss]) for ipw in 1:npw]
	isx2=[zeros(Int64,acqgeom[ipw].ns[iss]) for ipw in 1:npw]
	isz1=[zeros(Int64,acqgeom[ipw].ns[iss]) for ipw in 1:npw]
	isz2=[zeros(Int64,acqgeom[ipw].ns[iss]) for ipw in 1:npw]
	for ipw=1:npw
		for is=1:acqgeom[ipw].ns[iss]
			weights=ssprayw[ipw]
			denom=denomsI[ipw]
			Interpolation.get_spray_weights!(view(weights, :,is), view(denom,is), 
				    view(isx1[ipw],is), view(isx2[ipw],is),
				    view(isz1[ipw],is), view(isz2[ipw],is),
			    mesh_x, mesh_z, acqgeom[ipw].sx[iss][is], acqgeom[ipw].sz[iss][is])
		end
	end

	# receiver interpolation weights per sequential source
	rinterpolatew = [zeros(4,acqgeom[ipw].nr[iss]) for ipw in 1:npw]
	denomrI = [zeros(acqgeom[ipw].nr[iss]) for ipw in 1:npw]
	irx1=[zeros(Int64,acqgeom[ipw].nr[iss]) for ipw in 1:npw]
	irx2=[zeros(Int64,acqgeom[ipw].nr[iss]) for ipw in 1:npw]
	irz1=[zeros(Int64,acqgeom[ipw].nr[iss]) for ipw in 1:npw]
	irz2=[zeros(Int64,acqgeom[ipw].nr[iss]) for ipw in 1:npw]
	for ipw=1:npw
		for ir=1:acqgeom[ipw].nr[iss]
			weights=rinterpolatew[ipw]
			denom=denomrI[ipw]
			Interpolation.get_interpolate_weights!(view(weights, :,ir), view(denom, ir), 
				  view(irx1[ipw],ir), view(irx2[ipw],ir),
				  view(irz1[ipw],ir), view(irz2[ipw],ir),
			    mesh_x, mesh_z, acqgeom[ipw].rx[iss][ir], acqgeom[ipw].rz[iss][ir])
		end
	end


	pass=Paramss(iss,wavelets,ssprayw,records,rinterpolatew,
	      isx1,isx2,isz1,isz2,irx1,irx2,irz1,irz2,boundary,snaps,illum, 
	      grad_modtt,grad_modrrvx,grad_modrrvz)


	return pass
end

include("source.jl")
include("receiver.jl")
include("core.jl")
include("gradient.jl")
include("born.jl")
include("boundary.jl")


"""
This method updated the input `Fdtd.Param` after the wave propagation.

# Arguments

* `pa::Param` : modelling parameters

# Useful fields in `pa` that are modified by the method

* `pa.c.TDout::Vector{Data.TD}` : seismic data at receivers after modeling, for each propagating wavefield
* `pa.c.snaps::Array{Float64,4}` : snaps with size `(nz,nx,length(tsnaps),nss)` saved at `tsnaps`
* `pa.c.gradient::Models.Seismic` : gradient model modified only if `gmodel_flag`
* `pa.c.illum_stack::Array{Float64,2}` source energy of size `(nz,nx)` if `illum_flag`

# Example

```julia
JuMIT.Fdtd.mod!(pa)
```
# Credits 

Author: Pawan Bharadwaj 
        (bharadwaj.pawan@gmail.com)

"""
@fastmath function mod!(pa::Param=Param())

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
	pa.c.verbose && show(to)	
	return nothing
end

function update_datamat!(ifield, ipw, pac::Paramc, pap::Paramp)
	datamat=pac.datamat
        pass=pap.ss
	for issp in 1:length(pass)
		iss=pass[issp].iss
		records=pass[issp].records[ipw]
		for ir in 1:pac.acqgeom[ipw].nr[iss]
			for it in 1:pac.nt
				datamat[it,ir,iss]=records[it,ir,ifield]
			end
		end
        end
end

function update_data!(ifield, ipw, pac::Paramc)
	datamat=pac.datamat
	for iss in 1:pac.acqgeom[1].nss
		data=pac.data[ipw].d[iss,ifield]
		for ir in 1:pac.acqgeom[ipw].nr[iss]
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
##			Data.TD_resamp!(pac.data[ipropout], Data.TD_urpos((Array(records[:,:,iprop,:,:])), rfields, tgridmod, acqgeom[iprop],
##				acqgeom_urpos[1].nr[1],
##				(acqgeom_urpos[1].rz[1], acqgeom_urpos[1].rx[1])
##				)) 
#		end
#	end
	# return without resampling for testing
	#return [Data.TD(reshape(records[1+(iprop-1)*nd : iprop*nd],tgridmod.nx,recv_n,nss),
	#		       tgridmod, acqgeom[1]) for iprop in 1:npw]

function stack_illums!(pac::Paramc, pap::Paramp)
	np=pac.model.mgrid.npml
	nx, nz=pac.nx, pac.nz
	illums=pac.illum_stack
	pass=pap.ss
	for issp in 1:length(pass)
		gs=pass[issp].illum
		gss=view(gs,np+1:nz-np,np+1:nx-np)
		@. illums += gss
	end
end



# modelling for each processor
function mod_per_proc!(pac::Paramc, pap::Paramp) 
	# source_loop
	for issp in 1:length(pap.ss)
		reset_per_ss!(pap)

		iss=pap.ss[issp].iss

		if(pac.backprop_flag==-1)
			"initial conditions from boundary for first propagating field only"
			boundary_force_snap_p!(issp,pac,pap.ss,pap)
			boundary_force_snap_vxvz!(issp,pac,pap.ss,pap)
		end

		prog = Progress(pac.nt, dt=1, desc="\tmodeling supershot $iss/$(pac.acqgeom[1].nss) ", 
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

			(pac.born_flag) && add_born_sources!(issp, pac, pap.ss, pap)

			# record boundaries after time reversal already
			(pac.backprop_flag==1) && boundary_save!(pac.nt-it+1,issp,pac,pap.ss,pap)

			record!(it, issp, iss, pac, pap.ss, pap, Receiver_B1())

			(pac.gmodel_flag) && compute_gradient!(issp, pac, pap.ss, pap)

			(pac.illum_flag) && compute_illum!(issp, pap.ss, pap)

			if(pac.snaps_flag)
				itsnap = findall(in(pac.itsnaps),it)
				(itsnap ≠ []) && (snaps_save!(itsnap[1],issp,pac,pap.ss,pap))
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
		(pac.gmodel_flag) && scale_gradient!(issp, pap.ss, pac.model.mgrid.δx*pac.model.mgrid.δz)
		

	end # source_loop
end # mod_per_shot




# Need illumination to estimate the approximate diagonal of Hessian
@inbounds @fastmath function compute_illum!(issp::Int64, pass::Vector{Paramss}, pap::Paramp)
	# saving illumination to be used as preconditioner 
	p=pap.p
	illum=pass[issp].illum
	pp=view(p[1],:,:,1)
	for i in eachindex(illum)
		illum[i] += abs2(pp[i])
	end
end


@fastmath @inbounds function snaps_save!(itsnap::Int64,issp::Int64,pac::Paramc,pass::Vector{Paramss},pap::Paramp)
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


"""
Project density on to a staggerred grid using simple interpolation
"""
function get_rhovxI(rhoI::Array{Float64})
	rhovxI = zero(rhoI);
	get_rhovxI!(rhovxI, rhoI)
	return rhovxI
end
function get_rhovxI!(rhovxI, rhoI::Array{Float64})
	for ix = 1:size(rhoI, 2)-1
		for iz = 1:size(rhoI,1)
			rhovxI[iz, ix] = 0.5e0 *(rhoI[iz,ix+1] + rhoI[iz,ix])
		end
	end
	return nothing
end # get_rhovxI

"""
Project density on to a staggerred grid using simple interpolation
"""
function get_rhovzI(rhoI::Array{Float64})
	rhovzI = zero(rhoI);
	get_rhovzI!(rhovzI, rhoI)
	return rhovzI
end
function get_rhovzI!(rhovzI, rhoI::Array{Float64})
	for ix = 1: size(rhoI, 2)
		for iz = 1: size(rhoI, 1)-1
			rhovzI[iz, ix] =  0.5e0 * (rhoI[iz+1,ix] + rhoI[iz,ix])
		end
	end
	return nothing
end # get_rhovzI





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


end # module
