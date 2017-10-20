__precompile__()

module Fdtd

import JuMIT.Grid
import JuMIT.Interpolation
import JuMIT.Models
import JuMIT.Acquisition
import JuMIT.Data
import JuMIT.Gallery
import JuMIT.DSP
using DistributedArrays

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

* `gmodel::Models.Seismic=Models.Seismic_zeros(model.mgrid)` : gradient model modified only if `gmodel_flag`
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
type Paramc
	jobname::Symbol
	npw::Int64 
	exmodel::Models.Seismic
	model::Models.Seismic
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
	modttI::Matrix{Float64}
	modrrvx::Matrix{Float64}
	modrrvz::Matrix{Float64}
	δmodtt::Matrix{Float64}
	δmodrrvx::Matrix{Float64}
	δmodrrvz::Matrix{Float64}
	grad_modtt_stack::SharedMatrix{Float64}
	grad_modrrvx_stack::SharedMatrix{Float64}
	grad_modrrvz_stack::SharedMatrix{Float64}
	grad_modrr_stack::Matrix{Float64}
	illum_flag::Bool
	illum_stack::SharedMatrix{Float64}
	backprop_flag::Int64
	snaps_flag::Bool
	itsnaps::Vector{Int64}
	born_flag::Bool
	gmodel::Models.Seismic
	gmodel_flag::Bool
	ibx0::Int64
	ibz0::Int64
	ibx1::Int64
	ibz1::Int64
	isx0::Int64
	isz0::Int64
	datamat::SharedArray{Float64,3}
	data::Vector{Data.TD}
	verbose::Bool
end 

"""
Modelling parameters per every supersource for each worker
"""
type Paramss
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
	born_svalue_stack::Matrix{Float64} 
	grad_modtt::Matrix{Float64} 
	grad_modrrvx::Matrix{Float64}
	grad_modrrvz::Matrix{Float64}
end

"""
Parameters per every worker, not necessarily for every supersource.
Note that a single worker can take care of multiple supersources.
"""
type Paramp
	ss::Vector{Paramss}
	p::Array{Float64,4}
	pp::Array{Float64,4}
	ppp::Array{Float64,4}
	dpdx::Array{Float64,4}
	dpdz::Array{Float64,4}
	memory_dp_dx::Array{Float64,3}
	memory_dp_dz::Array{Float64,3}
	memory_dvx_dx::Array{Float64,3}
	memory_dvz_dz::Array{Float64,3}
end

type Param
	p::DistributedArrays.DArray{Paramp,1,Paramp} # distributed parameters among workers
	c::Paramc # common parameters
end

function initialize!(pap::Paramp)
	reset_per_ss!(pap)
	for issp in 1:length(pap.ss)
		pass=pap.ss[issp]
		for i in 1:length(pass.boundary)
			pass.boundary[i][:]=0.0
		end
		for i in 1:length(pass.records)
			pass.records[i][:]=0.0
		end
		pass.snaps[:]=0.0
		pass.illum[:]=0.0
		pass.grad_modtt[:]=0.0
		pass.grad_modrrvx[:]=0.0
		pass.grad_modrrvz[:]=0.0
		pass.born_svalue_stack[:]=0.0
	end
end

function initialize!(pac::Paramc)
	pac.grad_modtt_stack[:]=0.0
	pac.grad_modrrvx_stack[:]=0.0
	pac.grad_modrrvz_stack[:]=0.0
	pac.grad_modrr_stack[:]=0.0
	pac.illum_stack[:]=0.0
	for ipw=1:pac.npw
		dat=pac.data[ipw]
		Data.TD_zeros!(dat)
	end
	pac.datamat[:]=0.0
	Models.Seismic_zeros!(pac.gmodel)
end

function reset_per_ss!(pap::Paramp)
	pap.p[:]=0.0
	pap.pp[:]=0.0
	pap.ppp[:]=0.0
	pap.dpdx[:]=0.0
	pap.dpdz[:]=0.0
	pap.memory_dp_dz[:]=0.0
	pap.memory_dp_dx[:]=0.0
	pap.memory_dvx_dx[:]=0.0
	pap.memory_dvz_dz[:]=0.0
end


"""
Method to create `Fdtd` modeling parameters.
The output of this method can be used as an input to `mod!`, where the actual 
finite-difference modeling is performed.

# Keyword Arguments

* `npw::Int64=1` : number of independently propagating wavefields in `model`
* `model::Models.Seismic=Gallery.Seismic(:acou_homo1)` : seismic medium parameters 
* `model_pert::Models.Seismic=model` : perturbed model, i.e., model + δmodel, used only for Born modeling 
* `tgridmod::Grid.M1D=Gallery.M1D(:acou_homo1)` : modeling time grid, maximum time in tgridmod should be greater than or equal to maximum source time, same sampling interval as the wavelet
* `tgrid::Grid.M1D=tgridmod` : output records are resampled on this time grid
* `acqgeom::Vector{Acquisition.Geom}=fill(Gallery.Geom(:acou_homo1),npw)` :  acquisition geometry for each independently propagating wavefield
* `acqsrc::Vector{Acquisition.Src}=fill(Gallery.Src(:acou_homo1),npw)` : source acquisition parameters for each independently propagating wavefield
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
	model::Models.Seismic=Gallery.Seismic(:acou_homo1),
	abs_trbl::Vector{Symbol}=[:top, :bottom, :right, :left],
	born_flag::Bool=false,
	model_pert::Models.Seismic = model,
	tgridmod::Grid.M1D = Gallery.M1D(:acou_homo1),
	acqgeom::Vector{Acquisition.Geom} = fill(Gallery.Geom(model.mgrid,:surf,nss=1),npw),
	acqsrc::Array{Acquisition.Src} = fill(Gallery.Src(:acou_homo1),npw),
	sflags::Vector{Int64}=fill(2,npw), 
	rflags::Vector{Int64}=fill(1,npw),
	rfields::Vector{Symbol}=[:P], 
	backprop_flag::Int64=0,  
	gmodel_flag::Bool=false,
	illum_flag::Bool=false,
	tsnaps::Vector{Float64}=fill(0.5*(tgridmod.x[end]+tgridmod.x[1]),1),
	snaps_flag::Bool=false,
	verbose::Bool=false)

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
	isfields = [findin([:P, :Vx, :Vz], acqsrc[ipw].fields) for ipw in 1:npw]
	fill(nss, npw) != [getfield(acqgeom[ip],:nss) for ip=1:npw] ? error("different supersources") : nothing

	# create acquisition geometry with each source shooting 
	# at every unique receiver position
	irfields = findin([:P, :Vx, :Vz], rfields)


	if(verbose)
		println(string("number of super sources:\t",nss))	
	end

	# find maximum and minimum frequencies in the source wavelets
	freqmin = DSP.findfreq(acqsrc[1].wav[1,1][:,1],acqsrc[1].tgrid,attrib=:min) 
	freqmax = DSP.findfreq(acqsrc[1].wav[1,1][:,1],acqsrc[1].tgrid,attrib=:max) 

	# minimum and maximum velocities
	vpmin = minimum(broadcast(minimum,[model.vp0, model_pert.vp0]))
	vpmax = maximum(broadcast(maximum,[model.vp0, model_pert.vp0]))
	verbose ? println("minimum and maximum velocities:\t",vpmin,"\t",vpmax) : nothing


	check_fd_stability(vpmin, vpmax, model.mgrid.δx, model.mgrid.δz, freqmin, freqmax, tgridmod.δx, verbose)


	# some constants for boundary
	ibx0=model.mgrid.npml-2; ibx1=model.mgrid.nx+model.mgrid.npml+3
	ibz0=model.mgrid.npml-2; ibz1=model.mgrid.nz+model.mgrid.npml+3

	# for snaps
	isx0, isz0=model.mgrid.npml, model.mgrid.npml

	# extend models in the PML layers
	exmodel = Models.Seismic_pml_pad_trun(model);
	exmodel_pert = Models.Seismic_pml_pad_trun(model_pert);


	"density values on vx and vz stagerred grids"
	modrrvx = get_rhovxI(Models.Seismic_get(exmodel, :ρI))
	modrrvz = get_rhovzI(Models.Seismic_get(exmodel, :ρI))

	modttI = Models.Seismic_get(exmodel, :K) 
	modtt = Models.Seismic_get(exmodel, :KI)

	if(born_flag)
		"inverse density contrasts for Born Modelling"
		δmodrrvx = get_rhovxI(Models.Seismic_get(exmodel_pert, :ρI)) - get_rhovxI(Models.Seismic_get(exmodel, :ρI))
		δmodrrvz = get_rhovzI(Models.Seismic_get(exmodel_pert, :ρI)) - get_rhovzI(Models.Seismic_get(exmodel, :ρI))
		"inverse Bulk Modulus contrasts for Born Modelling"
		δmodtt = Models.Seismic_get(exmodel_pert, :KI) - Models.Seismic_get(exmodel, :KI)
	else
		δmodtt, δmodrrvx, δmodrrvz = zeros(1,1), zeros(1,1), zeros(1,1)
	end



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
	a_x, b_x, k_xI, a_x_half, b_x_half, k_x_halfI = pml_variables(nx, δt, δx, model.mgrid.npml-5, vpmax, vpmin, freqmin, freqmax, 
							       [any(abs_trbl .== :left), any(abs_trbl .== :right)])
	a_z, b_z, k_zI, a_z_half, b_z_half, k_z_halfI = pml_variables(nz, δt, δz, model.mgrid.npml-5, vpmax, vpmin, freqmin, freqmax,
							       [any(abs_trbl .== :top), any(abs_trbl .== :bottom)])

	grad_modtt_stack=SharedMatrix{Float64}(zeros(nzd,nxd))
	grad_modrrvx_stack=SharedMatrix{Float64}(zeros(nzd,nxd))
	grad_modrrvz_stack=SharedMatrix{Float64}(zeros(nzd,nxd))
	grad_modrr_stack=zeros(nzd,nxd)
	gmodel=Models.Seismic_zeros(model.mgrid)
	illum_stack=SharedMatrix{Float64}(zeros(nzd, nxd))

	itsnaps = [indmin(abs.(tgridmod.x-tsnaps[i])) for i in 1:length(tsnaps)]

	nrmat=[acqgeom[ipw].nr[iss] for ipw in 1:npw, iss in 1:acqgeom[1].nss]
	datamat=SharedArray{Float64}(nt,maximum(nrmat),acqgeom[1].nss)
	data=[Data.TD_zeros(rfields,tgridmod,acqgeom[ip]) for ip in 1:length(findn(rflags))]

	pac=Paramc(jobname,npw,
	    exmodel,model,
	    acqgeom,acqsrc,abs_trbl,isfields,sflags,
	    irfields,rflags,δt,δxI,δzI,
            nx,nz,nt,δtI,δx24I,δz24I,a_x,b_x,k_xI,a_x_half,b_x_half,k_x_halfI,a_z,b_z,k_zI,a_z_half,b_z_half,k_z_halfI,
	    modttI,modrrvx,modrrvz,
	    δmodtt,δmodrrvx,δmodrrvz,
	    grad_modtt_stack,grad_modrrvx_stack,grad_modrrvz_stack,grad_modrr_stack,
	    illum_flag,illum_stack,
	    backprop_flag,
	    snaps_flag,
	    itsnaps,born_flag,
	    gmodel,
	    gmodel_flag,
	    ibx0,ibz0,ibx1,ibz1,
	    isx0,isz0,
	    datamat,
	    data,
	    verbose)	

	# dividing the supersources to workers
	nwork = min(nss, nworkers())
	work = workers()[1:nwork]
	ssi=[round(Int, s) for s in linspace(0,nss,nwork+1)]
	sschunks=Array{UnitRange{Int64}}(nwork)
	for ib in 1:nwork       
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

	p=zeros(nz,nx,3,npw); pp=zeros(p); ppp=zeros(p)
	dpdx=zeros(p); dpdz=zeros(p)

	memory_dvx_dx=zeros(nz,nx,npw)
	memory_dvx_dz=zeros(memory_dvx_dx)
	memory_dvz_dx=zeros(memory_dvx_dx)
	memory_dvz_dz=zeros(memory_dvx_dx)
	memory_dp_dx=zeros(memory_dvx_dx)
	memory_dp_dz=zeros(memory_dvx_dx)
	
	ss=[Paramss(iss, pac) for (issp,iss) in enumerate(sschunks)]

	pap=Paramp(ss,p,pp,ppp,dpdx,dpdz,memory_dp_dx,memory_dp_dz,memory_dvx_dx,memory_dvz_dz)

	return pap
end

"""
update the `Seismic` model in `Paramc`
"""
function update_model!(pac::Paramc, model::Models.Seismic)


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

	if(pac.born_flag)
		born_svalue_stack = zeros(nz, nx)
	else
		born_svalue_stack = zeros(1,1)
	end

	# storing boundary values for back propagation
	nx1, nz1=pac.model.mgrid.nx, pac.model.mgrid.nz
	npml=pac.model.mgrid.npml
	boundary=[zeros(3,nx1+6,nt),
	  zeros(nz1+6,3,nt),
	  zeros(3,nx1+6,nt),
	  zeros(nz1+6,3,nt),
	  zeros(nz1+2*npml,nx1+2*npml,3)
				    ]
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
	      isx1,isx2,isz1,isz2,irx1,irx2,irz1,irz2,boundary,snaps,illum,born_svalue_stack,
	      grad_modtt,grad_modrrvx,grad_modrrvz)


	return pass
end

"""
This method updated the input `Fdtd.Param` after the wave propagation.

# Arguments

* `pa::Param` : modelling parameters

# Useful fields in `pa` that are modified by the method

* `pa.c.TDout::Vector{Data.TD}` : seismic data at receivers after modeling, for each propagating wavefield
* `pa.c.snaps::Array{Float64,4}` : snaps with size `(nz,nx,length(tsnaps),nss)` saved at `tsnaps`
* `pa.c.gmodel::Models.Seismic` : gradient model modified only if `gmodel_flag`
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


	# all localparts of DArray are input to this method
	# parallelization over shots
	@sync begin
		for (ip, p) in enumerate(procs(pa.p))
			@async remotecall_wait(p) do 
				mod_per_proc!(pa.c, localpart(pa.p))
			end
		end
	end

	# stack gradients and illum over sources
	@sync begin
		for (ip, p) in enumerate(procs(pa.p))
			@sync remotecall_wait(p) do 
				(pa.c.gmodel_flag) && stack_grads!(pa.c, localpart(pa.p))
				(pa.c.illum_flag) && stack_illums!(pa.c, localpart(pa.p))
			end
		end
	end

	# update gradient model using grad_modtt_stack, grad_modrr_stack
	update_gmodel!(pa.c)

	for ipw in 1:pa.c.npw
		for ifield in 1:length(pa.c.irfields)

			pa.c.datamat[:]=0.0
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


function grad_modrr!(pac::Paramc)
	@simd for i in eachindex(pac.grad_modrr_stack)
		@inbounds pac.grad_modrr_stack[i] = (pac.grad_modrrvx_stack[i] + pac.grad_modrrvz_stack[i]) * (0.5)
	end
	grad_modrr_sprayrrvx!(pac.grad_modrr_stack,pac.grad_modrrvx_stack)
	grad_modrr_sprayrrvz!(pac.grad_modrr_stack,pac.grad_modrrvz_stack)
end
function grad_modrr_sprayrrvx!(grad_modrr_stack,grad_modrrvx_stack)
	for ix=2:size(grad_modrr_stack,2)-1
		for iz=2:size(grad_modrr_stack,1)-1
			@inbounds grad_modrr_stack[iz,ix+1] +=  0.5e0 * grad_modrrvx_stack[iz,ix]
		end
	end
end
function grad_modrr_sprayrrvz!(grad_modrr_stack,grad_modrrvz_stack)
	for ix=2:size(grad_modrr_stack,2)-1
		for iz=2:size(grad_modrr_stack,1)-1
			@inbounds grad_modrr_stack[iz+1,ix] +=  0.5e0 * grad_modrrvz_stack[iz,ix]
		end
	end
end

function stack_grads!(pac::Paramc, pap::Paramp)
	np=pac.model.mgrid.npml
	nx, nz=pac.nx, pac.nz

	# theses are SharedArrays
	gmodtt=pac.grad_modtt_stack
	gmodrrvx=pac.grad_modrrvx_stack
	gmodrrvz=pac.grad_modrrvz_stack
	pass=pap.ss
	for issp in 1:length(pass)
		gs=pass[issp].grad_modtt
		gss=view(gs,np+1:nz-np,np+1:nx-np)
		@. gmodtt += gss
		gs=pass[issp].grad_modrrvx
		gss=view(gs,np+1:nz-np,np+1:nx-np)
		@. gmodrrvx += gss
		gs=pass[issp].grad_modrrvz
		gss=view(gs,np+1:nz-np,np+1:nx-np)
		@. gmodrrvz += gss
	end
	# combine rrvx and rrvz
	grad_modrr!(pac::Paramc)
end

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


function update_gmodel!(pac::Paramc)
	KI0=Models.Seismic_get(pac.model,:KI0)
	ρI0=Models.Seismic_get(pac.model,:ρI0)
	gmodtt=view(pac.grad_modtt_stack,:)
	gmodrr=view(pac.grad_modrr_stack,:)
	
	Models.χg!(gmodtt,KI0,1)
	Models.χg!(gmodrr,ρI0,1)

	Models.Seismic_chainrule!(pac.gmodel, pac.model,gmodtt, gmodrr, [:χKI, :χρI], 1)

	Models.χg!(gmodtt,KI0,-1)
	Models.χg!(gmodrr,ρI0,-1)
end

# modelling for each processor
function mod_per_proc!(pac::Paramc, pap::Paramp) 
	# source_loop
	for issp in 1:length(pap.ss)
		reset_per_ss!(pap)

		iss=pap.ss[issp].iss

		if(pac.verbose)
			println("modelling supershot:\t", iss)
		end
		
		if(pac.backprop_flag==-1)
			"initial conditions from boundary for first propagating field only"
			boundary_force_snap_p!(issp,pac,pap.ss,pap)
			boundary_force_snap_vxvz!(issp,pac,pap.ss,pap)
		end

		# time_loop
		"""
		* don't use shared arrays inside this time loop, for speed when using multiple procs
		"""
		for it=1:pac.nt

			advance!(pac,pap)
		
			# force p[1] on boundaries
			(pac.backprop_flag==-1) && boundary_force!(it,issp,pac,pap.ss,pap)
	 
			add_source!(it, issp, iss, pac, pap.ss, pap)

			(pac.born_flag) && add_born_sources!(issp, pac, pap.ss, pap)

			# record boundaries after time reversal already
			(pac.backprop_flag==1) && boundary_save!(pac.nt-it+1,issp,pac,pap.ss,pap)

			record!(it, issp, iss, pac, pap.ss, pap)

			(pac.gmodel_flag) && compute_gradient!(issp, pac, pap.ss, pap)

			(pac.illum_flag) && compute_illum!(issp, pap.ss, pap)

			if(pac.snaps_flag)
				itsnap = findin(pac.itsnaps,it)
				(itsnap ≠ []) && (snaps_save!(itsnap[1],issp,pac,pap.ss,pap))
			end

		end # time_loop
		"now pressure is at [nt], velocities are at [nt-1/2]"	

		"one more propagating step to save pressure at [nt+1] -- for time revarsal"
		advance!(pac,pap)

		"save last snap of pressure field"
		boundary_save_snap_p!(issp,pac,pap.ss,pap)

		"one more propagating step to save velocities at [nt+3/2] -- for time reversal"
		advance!(pac,pap)

		"save last snap of velocity fields with opposite sign for adjoint propagation"
		boundary_save_snap_vxvz!(issp,pac,pap.ss,pap)

		"scale gradients for each issp"
		(pac.gmodel_flag) && scale_gradient!(issp, pap.ss, pac.model.mgrid.δx*pac.model.mgrid.δz)
		

	end # source_loop
end # mod_per_shot

# This routine ABSOLUTELY should not allocate any memory, called inside time loop.
@inbounds @fastmath function add_source!(it::Int64, issp::Int64, iss::Int64, pac::Paramc, pass::Vector{Paramss}, pap::Paramp)
	# aliases
	p=pap.p;
	wavelets=pass[issp].wavelets
	acqgeom=pac.acqgeom
	isx1=pass[issp].isx1
	isx2=pass[issp].isx2
	isz1=pass[issp].isz1
	isz2=pass[issp].isz2
	ssprayw=pass[issp].ssprayw
	modttI=pac.modttI
	"""
	adding source to pressure field at [it] 
	"""
	for ipw = 1:pac.npw
	for (ifields, ifield) in enumerate(pac.isfields[ipw])
	@simd for is = 1:acqgeom[ipw].ns[iss]
		"""
		use wavelets at [it], i.e., sum of source terms
		until [it-1]
		division of source term with δx and δz (see Jan's fdelmodc manual)
		"""
		source_term = wavelets[ipw,it][is, ifields] * pac.δt * pac.δxI * pac.δzI
		
		"""
		multiplication with modttI
		"""
		p[isz1[ipw][is], isx1[ipw][is],ifield, ipw] += 
			source_term * 
			ssprayw[ipw][1,is] * 
			modttI[isz1[ipw][is], isx1[ipw][is]]  
		p[isz1[ipw][is], isx2[ipw][is],ifield, ipw] += 
			source_term * 
			ssprayw[ipw][2,is] * 
			modttI[isz1[ipw][is], isx2[ipw][is]]
		p[isz2[ipw][is], isx1[ipw][is],ifield, ipw] += 
			source_term * 
			ssprayw[ipw][3,is] * 
			modttI[isz2[ipw][is], isx1[ipw][is]]
		p[isz2[ipw][is], isx2[ipw][is],ifield, ipw] += 
			source_term * 
			ssprayw[ipw][4,is] * 
			modttI[isz2[ipw][is], isx2[ipw][is]]
	end
	end
	end
end


# This routine ABSOLUTELY should not allocate any memory, called inside time loop.
@inbounds @fastmath function record!(it::Int64, issp::Int64, iss::Int64, pac::Paramc, pass::Vector{Paramss}, pap::Paramp)
	p=pap.p
	rinterpolatew=pass[issp].rinterpolatew
	irx1=pass[issp].irx1
	irx2=pass[issp].irx2
	irz1=pass[issp].irz1
	irz2=pass[issp].irz2

	for ipw = 1:pac.npw
		recs=pass[issp].records[ipw]
		for (ifieldr, ifield) in enumerate(pac.irfields)
			@simd for ir = 1:pac.acqgeom[ipw].nr[iss]
				recs[it,ir,ifieldr]= 
				(
				p[irz1[ipw][ir],irx1[ipw][ir],ifield,ipw]*
				rinterpolatew[ipw][1,ir]+
				p[irz1[ipw][ir],irx2[ipw][ir],ifield,ipw]*
				rinterpolatew[ipw][2,ir]+
				p[irz2[ipw][ir],irx1[ipw][ir],ifield,ipw]*
				rinterpolatew[ipw][3,ir]+
				p[irz2[ipw][ir],irx2[ipw][ir],ifield,ipw]*
				rinterpolatew[ipw][4,ir]
				)
		end
	end
	end
end

# This routine ABSOLUTELY should not allocate any memory, called inside time loop.
@inbounds @fastmath function compute_gradient!(issp::Int64, pac::Paramc, pass::Vector{Paramss}, pap::Paramp)
	# aliases
	p=pap.p
	pp=pap.pp
	ppp=pap.ppp
	dpdx=pap.dpdx
	dpdz=pap.dpdz
	δtI=pac.δtI
	grad_modtt=pass[issp].grad_modtt
	grad_modrrvx=pass[issp].grad_modrrvx
	grad_modrrvz=pass[issp].grad_modrrvz

	gmodtt!(grad_modtt,p,pp,ppp,pac.nx,pac.nz,δtI)
	gmodrrvx!(grad_modrrvx,dpdx,pac.nx,pac.nz)
	gmodrrvz!(grad_modrrvz,dpdz,pac.nx,pac.nz)
end
@inbounds @fastmath function gmodtt!(grad_modtt,p,pp,ppp,nx,nz,δtI,)
	for ix=1:nx
		@simd for iz=1:nz
			# p at [it], pp at [it-1]	# dpdx and dpdz at [it]	# gradients w.r.t inverse of modtt, i.e., 1/rho/c^2 
			@inbounds grad_modtt[iz,ix] += ((-1.0 * (ppp[iz, ix, 1,1] + p[iz, ix,  1,1] - 2.0 * pp[iz, ix,  1,1]) * δtI * δtI) *  pp[iz,ix,1,2])
		end
	end
end
@inbounds @fastmath function gmodrrvx!(grad_modrrvx,dpdx,nx,nz)
	for ix=1:nx
		@simd for iz=1:nz
			# gradient w.r.t. inverse of rho on vx and vz grids
			@inbounds grad_modrrvx[iz,ix] += (- dpdx[iz,ix,1,2]*dpdx[iz,ix,1,1])
		end
	end

end
@inbounds @fastmath function gmodrrvz!(grad_modrrvz,dpdz,nx,nz)
	for ix=1:nx
		@simd for iz=1:nz
			# gradient w.r.t. inverse of rho on vx and vz grids
			@inbounds grad_modrrvz[iz,ix] += (- dpdz[iz,ix,1,2]*dpdz[iz,ix,1,1])
		end
	end
end

@inbounds @fastmath function scale_gradient!(issp::Int64,pass::Vector{Paramss},δ)
	grad_modtt=pass[issp].grad_modtt
	grad_modrrvx=pass[issp].grad_modrrvx
	grad_modrrvz=pass[issp].grad_modrrvz
	"gradient is formed by intergration over time, hence multiply with δt, but why not?"
	"I don't completely understand where the factors δx and δz are coming from..."
	"probably the source term should not be multiplied by δxI and δzI during adjoint propagation"
	scale!(grad_modtt,δ)
	scale!(grad_modrrvx,δ)
	scale!(grad_modrrvz,δ)
end




@inbounds @fastmath function advance!(pac::Paramc, pap::Paramp)
	# aliases
	p=pap.p; pp=pap.pp; ppp=pap.ppp;
	dpdx=pap.dpdx; dpdz=pap.dpdz;
	memory_dp_dx=pap.memory_dp_dx; memory_dp_dz=pap.memory_dp_dz; memory_dvx_dx=pap.memory_dvx_dx; memory_dvz_dz=pap.memory_dvz_dz
	modttI=pac.modttI; modrrvx=pac.modrrvx; modrrvz=pac.modrrvz
	δx24I=pac.δx24I; δz24I=pac.δz24I; δt=pac.δt
	nx=pac.nx; nz=pac.nz
	a_x=pac.a_x; b_x=pac.b_x; k_xI=pac.k_xI; a_x_half=pac.a_x_half; b_x_half=pac.b_x_half; k_x_halfI=pac.k_x_halfI 
	a_z=pac.a_z; b_z=pac.b_z; k_zI=pac.k_zI; a_z_half=pac.a_z_half; b_z_half=pac.b_z_half; k_z_halfI=pac.k_z_halfI

	pppppp!(p,pp,ppp)

	"""
	compute dpdx and dpdz at [it-1] for all propagating fields
	"""
	update_dpdx!(p, dpdx, δx24I, memory_dp_dx, b_x_half, a_x_half, k_x_halfI, nx, nz)
	update_dpdz!(p, dpdz, δz24I, memory_dp_dz, b_z_half, a_z_half, k_z_halfI, nx, nz)

	"""
	update velocity at [it-1/2] using 
	velocity at [it-3/2] and dpdx and dpdz at [it-1] 
	"""
	update_vx!(p, dpdx, δt, modrrvx, nx, nz)
	update_vz!(p, dpdz, δt, modrrvz, nx, nz)


	dvdx!(dpdx,p,memory_dvx_dx,b_x,a_x,k_xI,nz,nx,δx24I)
	dvdz!(dpdz,p,memory_dvz_dz,b_z,a_z,k_zI,nz,nx,δz24I)

	"""
	compute pressure at [it] using p at [it-1] and dvxdx
	and dvzdz at [it-1/2]
	"""
	pvzvx!(p,dpdx,dpdz,modttI,nz,nx,δt)

end
function pppppp!(p,pp,ppp)
	@simd for i in eachindex(p)
		@inbounds ppp[i]=pp[i]
		@inbounds pp[i]=p[i]
	end
end

@inbounds @fastmath function dvdx!(dpdx,p,memory_dvx_dx,b_x,a_x,k_xI,nz,nx,δx24I)
	for ipw = 1:size(p,4)
		for ix=3:nx-2
		@simd for iz=3:nz-2
			@inbounds dpdx[iz,ix,2,ipw] = (27.e0*p[iz,ix,2,ipw]-27.e0*p[iz,ix-1,2,ipw]-p[iz,ix+1,2,ipw]+p[iz,ix-2,2,ipw]) * (δx24I)
			@inbounds memory_dvx_dx[iz,ix,ipw] = b_x[ix] * memory_dvx_dx[iz,ix,ipw] + a_x[ix] * dpdx[iz,ix,2,ipw] # pml 
			@inbounds dpdx[iz,ix,2,ipw] = dpdx[iz,ix,2,ipw] * k_xI[ix] + memory_dvx_dx[iz,ix,ipw] # pml
		end
		end
	end
end

@inbounds @fastmath function dvdz!(dpdz,p,memory_dvz_dz,b_z,a_z,k_zI,nz,nx,δz24I)
	for ipw = 1:size(p,4)
		for ix=3:nx-2
		@simd for iz=3:nz-2
			@inbounds dpdz[iz,ix,3,ipw] = (27.e0*p[iz,ix,3,ipw]-27.e0*p[iz-1,ix,3,ipw]-p[iz+1,ix,3,ipw]+p[iz-2,ix,3,ipw]) * (δz24I)
			@inbounds memory_dvz_dz[iz,ix,ipw] = b_z[iz] * memory_dvz_dz[iz,ix,ipw] + a_z[iz] * dpdz[iz,ix,3,ipw] # pml
			@inbounds dpdz[iz,ix,3,ipw] = dpdz[iz,ix,3,ipw] * k_zI[iz] + memory_dvz_dz[iz,ix,ipw] # pml
		end
		end
	end
end

@inbounds @fastmath function pvzvx!(p,dpdx,dpdz,modttI,nz,nx,δt)
	for ipw = 1:size(p,4)
		for ix=3:nx-2
		@simd for iz=3:nz-2
			@inbounds p[iz,ix,1,ipw] += (modttI[iz,ix] * (dpdx[iz,ix,2,ipw] + dpdz[iz,ix,3,ipw])) * δt #* boundary_p(iz,ix)
		end
		end
	end
end

@inbounds @fastmath function update_dpdx!(p, dpdx, δx24I, memory_dp_dx, b_x_half, a_x_half, k_x_halfI, nx, nz)
	for ipw = 1:size(p,4)
		for ix = 3:nx-2
		@simd for iz = 3:nz-2
			@inbounds dpdx[iz,ix,1,ipw] = (27.e0*p[iz,ix+1,1,ipw]-27.e0*p[iz,ix,1,ipw]-p[iz,ix+2,1,ipw]+p[iz,ix-1,1,ipw]) * (δx24I)
			@inbounds memory_dp_dx[iz,ix,ipw] = b_x_half[ix] * memory_dp_dx[iz,ix,ipw] + a_x_half[ix] * dpdx[iz,ix,1,ipw] # pml
			@inbounds dpdx[iz,ix,1,ipw] = dpdx[iz,ix,1,ipw] * k_x_halfI[ix] + memory_dp_dx[iz,ix,ipw] # pml
		end
		end
	end
end

@inbounds @fastmath function update_dpdz!(p, dpdz, δz24I, memory_dp_dz, b_z_half, a_z_half, k_z_halfI, nx, nz)
	for ipw = 1:size(p,4)
		for ix = 3:nx-2
		@simd for iz = 3:nz-2
			@inbounds dpdz[iz,ix,1,ipw] = (27.e0*p[iz+1,ix,1,ipw]-27.e0*p[iz,ix,1,ipw]-p[iz+2,ix,1,ipw]+p[iz-1,ix,1,ipw]) * (δz24I)
			@inbounds memory_dp_dz[iz,ix,ipw] = b_z_half[iz] * memory_dp_dz[iz,ix,ipw] + a_z_half[iz] * dpdz[iz,ix,1,ipw] # pml
			@inbounds dpdz[iz,ix,1,ipw] = dpdz[iz,ix,1,ipw] * k_z_halfI[iz] + memory_dp_dz[iz,ix,ipw] # pml
		end
		end
	end
end

@inbounds @fastmath function update_vx!(p, dpdx, δt, modrrvx,  nx, nz)
	for ipw = 1:size(p,4)
		for ix=3:nx-2
		@simd for iz=3:nz-2
			@inbounds p[iz,ix,2,ipw] += (dpdx[iz,ix,1,ipw]) * δt * modrrvx[iz,ix] #* boundary_vx(iz,ix)
		end
		end
	end
end

@inbounds @fastmath function update_vz!(p, dpdz, δt,  modrrvz, nx, nz)
	for ipw = 1:size(p,4)
		for ix=3:nx-2
		@simd for iz=3:nz-2
			@inbounds p[iz,ix,3,ipw] +=  (dpdz[iz,ix,1,ipw]) * δt * modrrvz[iz,ix] #* boundary_vz(iz,ix)
		end
		end
	end
end


# Need illumination to estimate the approximate diagonal of Hessian
@inbounds @fastmath function compute_illum!(issp::Int64, pass::Vector{Paramss}, pap::Paramp)
	# saving illumination to be used as preconditioner 
	p=pap.p
	illum=pass[issp].illum
	pp=view(p,:,:,1,1)
	@. illum += pp * pp
end


function add_born_sources!(issp::Int64, pac::Paramc, pass::Vector{Paramss}, pap::Paramp)

	born_svalue_stack=pass[issp].born_svalue_stack
	δx24I=pac.δx24I; δz24I=pac.δz24I; 
	δxI=pac.δxI; δzI=pac.δzI; 
	δt=pac.δt
	δtI=pac.δtI
	δmodtt=pac.δmodtt; modttI=pac.modttI;
	δmodrrvx=pac.δmodrrvx; δmodrrvz=pac.δmodrrvz
	nx=pac.nx; nz=pac.nz
	p=pap.p; pp=pap.pp; ppp=pap.ppp;
	dpdx=pap.dpdx; dpdz=pap.dpdz;

	# secondary sources for Born modeling
	# adding born sources from pressure(:,:,1) to pressure(:,:,2)
	# upto until [it-2]
	# lambdaI scatterrer source term at [it-1]
	# p is at [it], pp is at [it-1], ppp is at [it-2]
	# dpdx is at [it-1] and dpdz is at [it-1]
	# modrrvx scatterrer source term at [it-1]
	# modrrvz scatterrer source term at [it-1]

	born_stacktt!(born_svalue_stack,p,pp,ppp,δmodtt,nx,nz,δtI,δt)
	born_stackrrvx!(born_svalue_stack,dpdx,δmodrrvx,nx,nz,δx24I,δt)
	born_stackrrvz!(born_svalue_stack,dpdz,δmodrrvz,nx,nz,δz24I,δt)
	born_stack!(p,born_svalue_stack,modttI,nx,nz,δt)
end
@inbounds @fastmath function born_stacktt!(born_svalue_stack,p,pp,ppp,δmodtt,nx,nz,δtI,δt)
	for ix=3:nx-2
		@simd for iz=3:nz-2
			born_svalue_stack[iz,ix] += 
			δt * ((-1.0 * (ppp[iz, ix, 1,1] + p[iz, ix,  1,1] - 2.0 * pp[iz, ix,  1,1]) * δmodtt[iz, ix] * δtI * δtI)) 
		end
	end
end
@inbounds @fastmath function born_stackrrvx!(born_svalue_stack,dpdx,δmodrrvx,nx,nz,δx24I,δt)
	for ix=3:nx-2
		@simd for iz=3:nz-2
			born_svalue_stack[iz,ix] += 
				δt * ((27.e0*dpdx[iz,ix,1,1] * δmodrrvx[iz,ix] -27.e0*dpdx[iz,ix-1,1,1] * δmodrrvx[iz,ix-1] -dpdx[iz,ix+1,1,1] * δmodrrvx[iz,ix+1] +dpdx[iz,ix-2,1,1] * δmodrrvx[iz,ix-2] ) * (δx24I)) 
		end
	end

end
@inbounds @fastmath function born_stackrrvz!(born_svalue_stack,dpdz,δmodrrvz,nx,nz,δz24I,δt)
	for ix=3:nx-2
		@simd for iz=3:nz-2
			born_svalue_stack[iz,ix] += 
				δt * ((27.e0*dpdz[iz,ix,1,1] * δmodrrvz[iz,ix] -27.e0*dpdz[iz-1,ix,1,1] * δmodrrvz[iz-1,ix] -dpdz[iz+1,ix,1,1] * δmodrrvz[iz+1,ix] +dpdz[iz-2,ix,1,1] * δmodrrvz[iz-2,ix] ) * (δz24I))  

		end
	end
end
@inbounds @fastmath function born_stack!(p,born_svalue_stack,modttI,nx,nz,δt)
	for ix=3:nx-2
		@simd for iz=3:nz-2
			p[iz,ix,1,2] += born_svalue_stack[iz,ix] * δt * modttI[iz,ix] #* δxI * δzI 
		end
	end
end

@inbounds @fastmath function boundary_force_snap_p!(issp::Int64,pac::Paramc,pass::Vector{Paramss},pap::Paramp)
	p=pap.p
	boundary=pass[issp].boundary[5]
	ps=view(p,:,:,1,1)
	bs=view(boundary,:,:,1)
	copy!(ps,bs)
end
@inbounds @fastmath function boundary_force_snap_vxvz!(issp::Int64,pac::Paramc,pass::Vector{Paramss},pap::Paramp)
	# initial conditions from boundary for first propagating field only
	p=pap.p
	boundary=pass[issp].boundary[5]
	ps=view(p,:,:,2,1)
	bs=view(boundary,:,:,2)
	copy!(ps,bs)
	ps=view(p,:,:,3,1)
	bs=view(boundary,:,:,3)
	copy!(ps,bs)
end
@fastmath @inbounds function boundary_force!(it::Int64,issp::Int64,pac::Paramc,pass::Vector{Paramss},pap::Paramp)
	boundary=pass[issp].boundary
	p=pap.p
	ibx0=pac.ibx0; ibz0=pac.ibz0; ibx1=pac.ibx1; ibz1=pac.ibz1
	boundaryf_l!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	boundaryf_r!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	boundaryf_t!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	boundaryf_b!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
end
@fastmath @inbounds function boundaryf_l!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	for ix=1:3
		@simd for iz=1:ibz1-ibz0+1
			p[ibz0+iz-1,ibx0+ix-1,1,1] = boundary[4][iz,ix,it]
		end
	end
end
@fastmath @inbounds function boundaryf_r!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	for ix=1:3
		@simd for iz=1:ibz1-ibz0+1
			p[ibz0+iz-1,ibx1-ix+1,1,1] = boundary[2][iz,ix,it]
		end
	end
end
@fastmath @inbounds function boundaryf_t!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	@simd for ix=1:ibx1-ibx0+1
		for iz=1:3
			p[ibz0+iz-1,ibx0+ix-1,1,1] = boundary[1][iz,ix,it]
		end
	end
end
@fastmath @inbounds function boundaryf_b!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	@simd for ix=1:ibx1-ibx0+1
		for iz=1:3
			p[ibz1-iz+1,ibx0+ix-1,1,1] = boundary[3][iz,ix,it]
		end
	end
end


@inbounds @fastmath function boundary_save_snap_p!(issp::Int64,pac::Paramc,pass::Vector{Paramss},pap::Paramp)
	p=pap.p
	boundary=pass[issp].boundary[5]
	ps=view(p,:,:,1,1)
	bs=view(boundary,:,:,1)
	copy!(bs, ps)
end
@inbounds @fastmath function boundary_save_snap_vxvz!(issp::Int64,pac::Paramc,pass::Vector{Paramss},pap::Paramp)
	p=pap.p
	boundary=pass[issp].boundary[5]
	#vx
	bs=view(boundary,:,:,2)
	ps=view(p,:,:,2,1)
	copy!(bs, ps)
	scale!(bs,-1.)
	# vz
	bs=view(boundary,:,:,3)
	ps=view(p,:,:,3,1)
	copy!(bs, ps)
	scale!(bs,-1.)
end
@fastmath @inbounds function boundary_save!(it::Int64,issp::Int64,pac::Paramc,pass::Vector{Paramss},pap::Paramp)
	boundary=pass[issp].boundary
	p=pap.p
	ibx0=pac.ibx0; ibz0=pac.ibz0; ibx1=pac.ibx1; ibz1=pac.ibz1
	boundarys_l!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	boundarys_r!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	boundarys_t!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	boundarys_b!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
end
@fastmath @inbounds function boundarys_l!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	for ix=1:3
		@simd for iz=1:ibz1-ibz0+1
			boundary[4][iz,ix,it] = p[ibz0+iz-1,ibx0+ix-1,1,1] 
		end
	end
end
@fastmath @inbounds function boundarys_r!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	for ix=1:3
		@simd for iz=1:ibz1-ibz0+1
			boundary[2][iz,ix,it] = p[ibz0+iz-1,ibx1-ix+1,1,1]
		end
	end
end
@fastmath @inbounds function boundarys_t!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	@simd for ix=1:ibx1-ibx0+1
		for iz=1:3
			boundary[1][iz,ix,it] = p[ibz0+iz-1,ibx0+ix-1,1,1]
		end
	end
end
@fastmath @inbounds function boundarys_b!(it,ibx0,ibx1,ibz0,ibz1,p,boundary)
	@simd for ix=1:ibx1-ibx0+1
		for iz=1:3
			boundary[3][iz,ix,it] = p[ibz1-iz+1,ibx0+ix-1,1,1]
		end
	end
end


@fastmath @inbounds function snaps_save!(itsnap::Int64,issp::Int64,pac::Paramc,pass::Vector{Paramss},pap::Paramp)
	isx0=pac.isx0
	isz0=pac.isz0
	p=pap.p
	snaps=pass[issp].snaps
	for ix=1:size(snaps,2)
		@simd for iz=1:size(snaps,1)
			snaps[iz,ix,itsnap]=p[isz0+iz,isx0+ix,1,1]
		end
	end
end


"""
Project density on to a staggerred grid using simple interpolation
"""
function get_rhovxI(rhoI::Array{Float64})
	rhovxI = zeros(rhoI);
	for ix = 1:size(rhoI, 2)-1
		for iz = 1:size(rhoI,1)
			rhovxI[iz, ix] = 0.5e0 *(rhoI[iz,ix+1] + rhoI[iz,ix])
		end
	end
	return rhovxI
end # get_rhovxI

"""
Project density on to a staggerred grid using simple interpolation
"""
function get_rhovzI(rhoI::Array{Float64})
	rhovzI = zeros(rhoI);
	for ix = 1: size(rhoI, 2)
		for iz = 1: size(rhoI, 1)-1
				rhovzI[iz, ix] =  0.5e0 * (rhoI[iz+1,ix] + rhoI[iz,ix])
		end
	end
	return rhovzI
end # get_rhovzI


function fill_wavelets!(iss::Int64, wavelets::Array{Array{Float64,2},2}, acqsrc::Array{Acquisition.Src}, sflags::Vector{Int64})

	npw = size(wavelets,1)
	nt = size(wavelets,2)
	δt = acqsrc[1].tgrid.δx
	for ipw=1:npw
		ns, snfield = size(wavelets[ipw])
		for ifield=1:snfield, is=1:ns
			snt = acqsrc[ipw].tgrid.nx;
			if(sflags[ipw] == 0)
				nothing # just put zeros, no sources added
			elseif(sflags[ipw] == 1)
				"ϕ[t] = s[t]"
				for it=1:snt
					source_term = acqsrc[ipw].wav[iss,ifield][it,is]
					wavelets[ipw,it][is,ifield] = source_term
				end
			elseif(sflags[ipw] == 2)
				"ϕ[t] = ∑₁ᵗ⁻¹ s[t]"
				source_term_stack = 0.0;
				if(ifield == 1)
					for it=1:snt-1
						source_term_stack += (acqsrc[ipw].wav[iss,ifield][it,is] .* δt)
						wavelets[ipw,it+1][is,ifield] = source_term_stack
					end
				else
					for it=2:snt-1
						source_term_stack += (((acqsrc[ipw].wav[iss,ifield][it,is] .* δt) +
						   (acqsrc[ipw].wav[iss,ifield][it-1,is] .* δt)) * 0.5)
						wavelets[ipw,it+1][is,ifield] = source_term_stack
					end

				end
				if(nt > snt)
					wavelets[ipw,snt+1:end][is,ifield] = wavelets[ipw,snt][is,ifield]
				end
			elseif(sflags[ipw] == 3)
				"use this to add source sink: need during adjoint propagation from boundary"
				"multiplication with -1 for subtraction"
				"time reversal"
				"as the source wavelet has to be subtracted before the propagation step, I shift here by one sample"
				"ϕ[t] = "
				source_term_stack = 0.0;
				for it=1:snt-1
					source_term_stack += (acqsrc[ipw].wav[iss,ifield][it,is] .* δt)
					wavelets[ipw,nt-it+1][is,ifield] = -1.0 * source_term_stack
				end
				if(nt > snt)
					nt_diff = nt-snt
					wavelets[ipw,1:nt_diff+1][is,ifield] = wavelets[ipw,nt_diff+2][is,ifield]
				end
			end
		end
	end
end


function check_fd_stability(vpmin::Float64, vpmax::Float64, δx::Float64, δz::Float64, freqmin::Float64, freqmax::Float64, δt::Float64, verbose::Bool)

	# check spatial sampling
	δs_temp=vpmin/5.0/freqmax;
	δs_max = maximum([δx, δz])
	all(δs_max .> δs_temp) ? 
			warn(string("spatial sampling\t",δs_max,"\ndecrease spatial sampling below:\t",δs_temp), once=true) :
			verbose ? println("spatial sampling\t",δs_max,"\tcan be as high as:\t",δs_temp) : nothing 

	# check time sampling
	δs_min = minimum([δx, δz])
	δt_temp=0.5*δs_min/vpmax
	all(δt .> δt_temp) ? 
			warn(string("time sampling\t",δt,"\ndecrease time sampling below:\t",δt_temp), once=true) :
			verbose ? println("time sampling\t",δt,"\tcan be as high as:\t",δt_temp) : nothing

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
		(alpha_x_half[ix] < 0.0)? alpha_x_half[ix] = 0.0 : nothing

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
