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

function initialize_boundary(nx::Int64, nz::Int64, npml::Int64, nss::Int64, nt::Int64)
	nwork = min(nss, nworkers())
	work = workers()[1:nwork]
	return DArray((nss,), work,[nwork]) do I
	    			[[zeros(3,nx+6,nt),
				  zeros(nz+6,3,nt),
				  zeros(3,nx+6,nt),
				  zeros(nz+6,3,nt),
				  zeros(nz+2*npml,nx+2*npml,3)
							    ] for is in 1:nss]
				end
end

type Param
	f::Paramf # fields
	c::Paramc # common parameters
end

"""
Modelling parameters common for all supersources
"""
type Paramc
	jobname::AbstractString
	npropwav::Int64 
	abs_trbl::Vector{Symbol}

	npropwav::Int64
	src_ifields::Vector{Vector{Int64}}
	src_nsmul::Int64

	npropwav::Int64
	recv_ifields::Vector{Int64}
	recv_n::Int64

	δt::Float64
	δxI::Float64
	δzI::Float64
	nx::Int64
	nz::Int64
	δtI::Float64 
	δx24I::Float64
	δz24I::Float64
	δt::Float64

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


	modttI::Array{Float64}
	modrrvx::Array{Float64}
	modrrvz::Array{Float64}

	born_flag,
	model::Models.Seismic=Gallery.Seismic(:acou_homo1),
	model_pert::Models.Seismic = model,
	tgridmod::Grid.M1D = Gallery.M1D(:acou_homo1),
	acqgeom::Vector{Acquisition.Geom} = fill(Gallery.Geom(model.mgrid,:surf,nss=1),npropwav),
	src_flags::Vector{Int64}=fill(2,npropwav), 
	recv_flags::Vector{Int64}=fill(1,npropwav),
	recv_fields::Vector{Symbol}=[:P], 
	TDout::Vector{Data.TD}=[Data.TD_zeros(recv_fields,tgridmod,acqgeom[ip]) for ip in 1:length(findn(recv_flags)) ],
	backprop_flag::Int64=0,  
	boundary::DistributedArrays.DArray{Array{Array{Float64,3},1},1,Array{Array{Array{Float64,3},1},1}}=
			initialize_boundary(model.mgrid.nx, model.mgrid.nz, model.mgrid.npml, acqgeom[1].nss, tgridmod.nx),
	gmodel_flag::Bool=false,
	gmodel::Models.Seismic=Models.Seismic_zeros(model.mgrid),
	illum_flag::Bool=false,
	illum::Array{Float64,2}=zeros(model.mgrid.nz, model.mgrid.nx),
	tsnaps::Vector{Float64}=fill(0.5*(tgridmod.x[end]+tgridmod.x[1]),1),
	snaps::Array{Float64,4}=zeros(model.mgrid.nz,model.mgrid.nx,length(tsnaps),acqgeom[1].nss),
	snaps_flag::Bool=false,
	verbose::Bool=false

	grad_modtt_stack
	grad_modrrvx_stack 
	grad_modrrvz_stack 

end 

"""
Modelling parameters per supersource
"""
type Paramss
	boundary
	snaps_out
	illum_all
	born_svalue_stack::Array{Float64} 

	src_wav_mat::Array{Array{Float64,3},2}
	src_spray_weights::Array{Float64}

	recv_out::AbstractArray{Float64} 
	rec_interpolate_weights::Array{Float64}

	issmulx1::Array{Int64} 
	issmulx2::Array{Int64}
	issmulz1::Array{Int64} 
	issmulz2::Array{Int64}
	irecx1::Array{Int64}
	irecx2::Array{Int64}
	irecz1::Array{Int64}
	irecz2::Array{Int64}
end

"""
Parameters per every worker, not necessarily for every supersource.
Note that a single worker can take care of multiple supersources.
"""
type Paramp

	p::Array{Float64}
	pp::Array{Float64}
	ppp::Array{Float64}
	dpdx::Array{Float64}
	dpdz::Array{Float64}

	memory_dp_dx::Array{Float64}
	memory_dp_dz::Array{Float64}
	memory_dvx_dx::Array{Float64}
	memory_dvz_dz::Array{Float64}

	grad_modtt::Array{Float64} 
	grad_modrrvx::Array{Float64}
	grad_modrrvz::Array{Float64}
end




function Param(;
	jobname::AbstractString = "Hello",
	npropwav::Int64=1, 
	model::Models.Seismic=Gallery.Seismic(:acou_homo1),
	abs_trbl::Vector{Symbol}=[:top, :bottom, :right, :left],
	born_flag::Bool=false,
	model_pert::Models.Seismic = model,
	tgridmod::Grid.M1D = Gallery.M1D(:acou_homo1),
	acqgeom::Vector{Acquisition.Geom} = fill(Gallery.Geom(model.mgrid,:surf,nss=1),npropwav),
	acqsrc::Array{Acquisition.Src} = fill(Gallery.Src(:acou_homo1),npropwav),
	src_flags::Vector{Int64}=fill(2,npropwav), 
	recv_flags::Vector{Int64}=fill(1,npropwav),
	recv_fields::Vector{Symbol}=[:P], 
	TDout::Vector{Data.TD}=[Data.TD_zeros(recv_fields,tgridmod,acqgeom[ip]) for ip in 1:length(findn(recv_flags)) ],
	backprop_flag::Int64=0,  
	boundary::DistributedArrays.DArray{Array{Array{Float64,3},1},1,Array{Array{Array{Float64,3},1},1}}=
			initialize_boundary(model.mgrid.nx, model.mgrid.nz, model.mgrid.npml, acqgeom[1].nss, tgridmod.nx),
	gmodel_flag::Bool=false,
	gmodel::Models.Seismic=Models.Seismic_zeros(model.mgrid),
	illum_flag::Bool=false,
	illum::Array{Float64,2}=zeros(model.mgrid.nz, model.mgrid.nx),
	tsnaps::Vector{Float64}=fill(0.5*(tgridmod.x[end]+tgridmod.x[1]),1),
	snaps::Array{Float64,4}=zeros(model.mgrid.nz,model.mgrid.nx,length(tsnaps),acqgeom[1].nss),
	snaps_flag::Bool=false,
	verbose::Bool=false)

	# check sizes and errors based on input
	(length(TDout) ≠ length(findn(recv_flags))) && error("TDout dimension")
	(length(acqgeom) ≠ npropwav) && error("acqgeom dimension")
	(length(acqsrc) ≠ npropwav) && error("acqsrc dimension")
	(length(src_flags) ≠ npropwav) && error("src_flags dimension")
	(length(recv_flags) ≠ npropwav) && error("recv_flags dimension")
	(maximum(tgridmod.x) < maximum(acqsrc[1].tgrid.x)) && error("modeling time is less than source time")
	(any([getfield(TDout[ip],:tgrid).δx < tgridmod.δx for ip=1:length(TDout)])) && error("output time grid sampling finer than modeling")
	any([maximum(getfield(TDout[ip],:tgrid).x) > maximum(tgridmod.x) for ip=1:length(TDout)]) && error("output time > modeling time")

	#! no modeling if source wavelet is zero
	#if(maxval(abs(src_wavelets)) .lt. tiny(rzero_de)) 
	#        return
	#endif

	# all the propagating wavefields should have same supersources? check that?

	# check dimension of model and model_pert

	# check if all sources are receivers are inside model
	any(.![Acquisition.Geom_check(acqgeom[ip], model.mgrid) for ip=1:npropwav]) ? error("sources or receivers not inside model") : nothing



	length(acqgeom) != npropwav ? error("acqgeom size") : nothing
	length(acqsrc) != npropwav ? error("acqsrc size") : nothing
	any([getfield(acqgeom[ip],:nss) != getfield(acqsrc[ip],:nss) for ip=1:npropwav])  ? error("different supersources") : nothing
	any([getfield(acqgeom[ip],:ns) != getfield(acqsrc[ip],:ns) for ip=1:npropwav])  ? error("different sources") : nothing

	# necessary that nss and fields should be same for all nprop
	src_nseq = acqgeom[1].nss;
	src_fields = [acqsrc[ipropwav].fields for ipropwav in 1:npropwav]
	src_ifields = [findin([:P, :Vx, :Vz], acqsrc[ipropwav].fields) for ipropwav in 1:npropwav]
	fill(src_nseq, npropwav) != [getfield(acqgeom[ip],:nss) for ip=1:npropwav] ? error("different supersources") : nothing

	# create acquisition geometry with each source shooting 
	# at every unique receiver position
	acqgeom_urpos = Acquisition.Geom_get(acqgeom,:geomurpos);
	recv_n = acqgeom_urpos[1].nr[1] # same for all sources

	# same number of sources for all super sources
	acqgeom_uspos = Acquisition.Geom_get(acqgeom,:geomuspos);
	acqsrc_uspos = Acquisition.Src_uspos(acqsrc,acqgeom);
	src_nsmul = acqsrc_uspos[1].ns[1]; # same number of sources in all


	if(verbose)
		println(string("number of receivers:\t",recv_n))	
		println(string("number of sources:\t",src_nsmul))	
		println(string("number of super sources:\t",src_nseq))	
	end


	#create some aliases
	nx, nz = exmodel.mgrid.nx, exmodel.mgrid.nz
	nxd, nzd = model.mgrid.nx, model.mgrid.nz
	δx, δz = exmodel.mgrid.δx, exmodel.mgrid.δz
	δt = tgridmod.δx 
	tim_nt=tgridmod.nx


	
	
	pass = 
	pass DArray((nss,), work, [nwork]) do I
		[Paramss() for is in 1:nss]
	end
	 [zeros(3,nx+6,nt),
				  zeros(nz+6,3,nt),
				  zeros(3,nx+6,nt),
				  zeros(nz+6,3,nt),
				  zeros(nz+2*npml,nx+2*npml,3)
							    ] for is in 1:nss]


	# all localparts of DArrays are input to this method
	@sync begin
		for (ip, p) in enumerate(procs(boundary))
			@async remotecall_wait(p) do 
				Paramss_update_locals!(localindexes(boundary)[1], localpart(boundary),  
			  	               localpart(grad_modtt), localpart(grad_modrrvx), localpart(grad_modrrvz),
					       localpart(recv_out),  
					       localpart(snaps_out), localpart(illum_all),
					       snaps_flag,
					       illum_flag, backprop_flag, gmodel_flag, verbose, recv_n, recv_ifields, npropwav, src_nsmul, 
					       exmodel, model, tgridmod,
					       acqgeom_uspos, acqgeom_urpos,
					       modttI, modrrvx, modrrvz,
				               vpmax, vpmin, freqmin, freqmax, 
					       abs_trbl, tsnaps, src_wav_mat, src_ifields, 
					       born_flag,
					       localpart(born_svalue_stack), δmodtt,
					       δmodrrvx, δmodrrvz)
			end
		end
	end

end


function Paramc()
	# find maximum and minimum frequencies in the source wavelets
	freqmin = DSP.findfreq(acqsrc[1].wav[1,1][:,1],acqsrc[1].tgrid,attrib=:min) 
	freqmax = DSP.findfreq(acqsrc[1].wav[1,1][:,1],acqsrc[1].tgrid,attrib=:max) 

	# minimum and maximum velocities
	vpmin = minimum(broadcast(minimum,[model.vp0, model_pert.vp0]))
	vpmax = maximum(broadcast(maximum,[model.vp0, model_pert.vp0]))
	verbose ? println("minimum and maximum velocities:\t",vpmin,"\t",vpmax) : nothing


	check_fd_stability(vpmin, vpmax, model.mgrid.δx, model.mgrid.δz, freqmin, freqmax, tgridmod.δx, verbose)



	ibx0=model.mgrid.npml-2; ibx1=model.mgrid.nx+model.mgrid.npml+3
	ibz0=model.mgrid.npml-2; ibz1=model.mgrid.nz+model.mgrid.npml+3

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
		δmodtt, δmodrrvx, δmodrrvz = [0.0], [0.0], [0.0]
	end

	#create some aliases
	nx, nz = exmodel.mgrid.nx, exmodel.mgrid.nz
	nxd, nzd = model.mgrid.nx, model.mgrid.nz
	δx, δz = exmodel.mgrid.δx, exmodel.mgrid.δz
	δt = tgridmod.δx 
	tim_nt=tgridmod.nx

	# pml_variables
	a_x, b_x, k_xI, a_x_half, b_x_half, k_x_halfI = pml_variables(nx, δt, δx, model.mgrid.npml-5, vpmax, vpmin, freqmin, freqmax, 
							       [any(abs_trbl .== :left), any(abs_trbl .== :right)])
	a_z, b_z, k_zI, a_z_half, b_z_half, k_z_halfI = pml_variables(nz, δt, δz, model.mgrid.npml-5, vpmax, vpmin, freqmin, freqmax,
							       [any(abs_trbl .== :top), any(abs_trbl .== :bottom)])


	
	nwork = min(src_nseq, nworkers())
	work = workers()[1:nwork]

	return Paramc(ibx0, ibz0, exmodel, exmodel_pert, modrrvx, modrrvz, modttI, modtt, δmodtt, δmodrrvx, δmodrrvz,
	       freqmin, freqmax, vpmin, vpmax, a_x, b_x, k_xI, a_x_half, b_x_half, k_x_halfI,
	       a_z, b_z, k_zI, a_z_half, b_z_half, k_z_halfI
	       )
end

function Paramss()

	"""
	allocation of private variables
	"""
	p=zeros(nz,nx,3,npropwav); pp=zeros(p); ppp=zeros(p)
	dpdx=zeros(p); dpdz=zeros(p)

	memory_dvx_dx=zeros(nz,nx,npropwav)
	memory_dvx_dz=zeros(memory_dvx_dx)
	memory_dvz_dx=zeros(memory_dvx_dx)
	memory_dvz_dz=zeros(memory_dvx_dx)
	memory_dp_dx=zeros(memory_dvx_dx)
	memory_dp_dz=zeros(memory_dvx_dx)


function Paramss_update_locals!(sseqvec::UnitRange{Int64}, pass::Array{Paramss})

	# records_output, distributed array among different procs
	recv_ifields = findin([:P, :Vx, :Vz], recv_fields)
	recv_out = zeros(recv_n,length(recv_ifields),npropwav,tim_nt,src_nseq)


	# gradient outputs
	grad_modtt = zeros((nz, nx, src_nseq)
	grad_modrrvx = zeros((nz, nx, src_nseq)
	grad_modrrvz = zeros((nz, nx, src_nseq)


	"saving illum"
	illum_all =  (illum_flag) ? dzeros((nz, nx, src_nseq), work, [1,1,nwork]) : [0.0]

	snaps_out = (snaps_flag) ? dzeros(size(snaps), work, [1,1,1,nwork]) : [0.0]


	src_nseq=length(sseqvec)
	# source wavelets
	src_wav_mat = [zeros(src_nsmul,length(src_ifields[ipropwav]),tim_nt) for ipropwav in 1:npropwav, isseq in 1:src_nseq]
	fill_src_wav_mat!(src_wav_mat, acqsrc_uspos, src_flags)


	if(born_flag)
		born_svalue_stack = dzeros((exmodel.mgrid.nz, exmodel.mgrid.nx, src_nseq), work, [1,1,nwork])
	else
		born_svalue_stack = [0.0]
	end

	# source_loop
	for (iisseq,isseq) in enumerate(sseqvec)

		if(verbose)
			println("modelling supershot:\t", isseq)
		end
		
		"source_spray_weights per sequential source"
		src_spray_weights = zeros(4,src_nsmul,npropwav)
		denomsrcI = zeros(src_nsmul,npropwav)
		issmulx1 = fill(0,src_nsmul,npropwav); issmulx2=fill(0,src_nsmul,npropwav);
		issmulz1 = fill(0,src_nsmul,npropwav); issmulz2=fill(0,src_nsmul,npropwav);
		for ipropwav=1:npropwav
			for issmul=1:src_nsmul
				Interpolation.get_spray_weights!(view(src_spray_weights, :,issmul,ipropwav), view(denomsrcI,issmul,ipropwav), 
				    view(issmulx1,issmul,ipropwav), view(issmulx2,issmul,ipropwav),
				    view(issmulz1,issmul,ipropwav), view(issmulz2,issmul,ipropwav),
				    mesh_x, mesh_z, acqgeom_uspos[ipropwav].sx[isseq][issmul], acqgeom_uspos[ipropwav].sz[isseq][issmul])
			end
		end

		"receiver interpolation weights per sequential source"
		rec_interpolate_weights = zeros(4,recv_n,npropwav)
		denomrecI = zeros(recv_n,npropwav)
		irecx1 = fill(0,recv_n,npropwav); irecx2=fill(0,recv_n,npropwav);
		irecz1 = fill(0,recv_n,npropwav); irecz2=fill(0,recv_n,npropwav);
		for ipropwav=1:npropwav
			for irec=1:recv_n
				Interpolation.get_interpolate_weights!(view(rec_interpolate_weights, :,irec,ipropwav), view(denomrecI,irec,ipropwav), 
				    view(irecx1,irec,ipropwav), view(irecx2,irec,ipropwav),
				    view(irecz1,irec,ipropwav), view(irecz2,irec,ipropwav),
				    mesh_x, mesh_z, acqgeom_urpos[ipropwav].rx[isseq][irec], acqgeom_urpos[ipropwav].rz[isseq][irec])
			end
		end




end

"""
```math
\alpha=2
```
# Keyword Arguments

* `npropwav::Int64=1` : number of independently propagating wavefields in `model`
* `model::Models.Seismic=Gallery.Seismic(:acou_homo1)` : seismic medium parameters 
* `model_pert::Models.Seismic=model` : perturbed model, i.e., model + δmodel, used only for Born modeling 
* `tgridmod::Grid.M1D=Gallery.M1D(:acou_homo1)` : modeling time grid, maximum time in tgridmod should be greater than or equal to maximum source time, same sampling interval as the wavelet
* `tgrid::Grid.M1D=tgridmod` : output records are resampled on this time grid
* `acqgeom::Vector{Acquisition.Geom}=fill(Gallery.Geom(:acou_homo1),npropwav)` :  acquisition geometry for each independently propagating wavefield
* `acqsrc::Vector{Acquisition.Src}=fill(Gallery.Src(:acou_homo1),npropwav)` : source acquisition parameters for each independently propagating wavefield
* `src_flags::Vector{Int64}=fill(2,npropwav)` : source related flags for each propagating wavefield
  * `=[0]` inactive sources
  * `=[1]` sources with injection rate
  * `=[2]` volume injection sources
  * `=[3]` sources input after time reversal (use only during backpropagation) 
* `recv_flags::Vector{Int64}=fill(1,npropwav)` : receiver related flags for each propagating wavefield
  * `=[0]` receivers do not record (or) inactive receivers
  * `=[0,1]` receivers are active only for the second propagating wavefield
* `recv_fields::Vector{Symbol}=[:P]` : multi-component receiver flag; types fields the receivers record (to be changed later)
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

# Keyword Arguments that are modified by the method (some of them are returned as well)

* `gmodel::Models.Seismic=Models.Seismic_zeros(model.mgrid)` : gradient model modified only if `gmodel_flag`
* `TDout::Vector{Data.TD}=[Data.TD_zeros(recv_fields,tgridmod,acqgeom[ip]) for ip in 1:length(findn(recv_flags))]`
* `illum::Array{Float64,2}=zeros(model.mgrid.nz, model.mgrid.nx)` : source energy if `illum_flag`
* `boundary::Array{Array{Float64,4},1}` : stored boundary values for first propagating wavefield 
* `snaps::Array{Float64,4}=zeros(model.mgrid.nz,model.mgrid.nx,length(tsnaps),acqgeom[1].nss)` :snapshots saved at `tsnaps`

# Return (in order)

* modelled data for each propagating wavefield as `Vector{TD}`
* stored boundary values of the first propagating wavefield as `Array{Array{Float64,4},1}` (use for backpropagation)
* final conditions of the first propagating wavefield as `Array{Float64,4}` (use for back propagation)
* gradient model as `Seismic`
* stored snaps shots at tsnaps as Array{Float64,4} 

# Example

```julia
records, boundary, gmodel, snaps = mod(acqgeom=acqgeom, acqsrc=acqsrc, model=model, tgridmod=tgridmod);
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
"""
@fastmath function mod!(;
	jobname::AbstractString = "Hello",
	npropwav::Int64=1, 
	model::Models.Seismic=Gallery.Seismic(:acou_homo1),
	abs_trbl::Vector{Symbol}=[:top, :bottom, :right, :left],
	born_flag::Bool=false,
	model_pert::Models.Seismic = model,
	tgridmod::Grid.M1D = Gallery.M1D(:acou_homo1),
	acqgeom::Vector{Acquisition.Geom} = fill(Gallery.Geom(model.mgrid,:surf,nss=1),npropwav),
	acqsrc::Array{Acquisition.Src} = fill(Gallery.Src(:acou_homo1),npropwav),
	src_flags::Vector{Int64}=fill(2,npropwav), 
	recv_flags::Vector{Int64}=fill(1,npropwav),
	recv_fields::Vector{Symbol}=[:P], 
	TDout::Vector{Data.TD}=[Data.TD_zeros(recv_fields,tgridmod,acqgeom[ip]) for ip in 1:length(findn(recv_flags)) ],
	backprop_flag::Int64=0,  
	boundary::DistributedArrays.DArray{Array{Array{Float64,3},1},1,Array{Array{Array{Float64,3},1},1}}=
			initialize_boundary(model.mgrid.nx, model.mgrid.nz, model.mgrid.npml, acqgeom[1].nss, tgridmod.nx),
	gmodel_flag::Bool=false,
	gmodel::Models.Seismic=Models.Seismic_zeros(model.mgrid),
	illum_flag::Bool=false,
	illum::Array{Float64,2}=zeros(model.mgrid.nz, model.mgrid.nx),
	tsnaps::Vector{Float64}=fill(0.5*(tgridmod.x[end]+tgridmod.x[1]),1),
	snaps::Array{Float64,4}=zeros(model.mgrid.nz,model.mgrid.nx,length(tsnaps),acqgeom[1].nss),
	snaps_flag::Bool=false,
	verbose::Bool=false
	)
	

	# all localparts of DArrays are input to this method
	@sync begin
		for (ip, p) in enumerate(procs(boundary))
			@async remotecall_wait(p) do 
				mod_per_proc!(localindexes(boundary)[1], localpart(boundary),  
			  	               localpart(grad_modtt), localpart(grad_modrrvx), localpart(grad_modrrvz),
					       localpart(recv_out),  
					       localpart(snaps_out), localpart(illum_all),
					       snaps_flag,
					       illum_flag, backprop_flag, gmodel_flag, verbose, recv_n, recv_ifields, npropwav, src_nsmul, 
					       exmodel, model, tgridmod,
					       acqgeom_uspos, acqgeom_urpos,
					       modttI, modrrvx, modrrvz,
				               vpmax, vpmin, freqmin, freqmax, 
					       abs_trbl, tsnaps, src_wav_mat, src_ifields, 
					       born_flag,
					       localpart(born_svalue_stack), δmodtt,
					       δmodrrvx, δmodrrvz)
			end
		end
	end

	"summing gradient over all the sources after truncation in PML"
	grad_modtt_stack = Models.pml_pad_trun(squeeze(sum(Array(grad_modtt),3),3),model.mgrid.npml,-1);
	grad_modrrvx_stack = Models.pml_pad_trun(squeeze(sum(Array(grad_modrrvx),3),3),model.mgrid.npml,-1);
	grad_modrrvz_stack = Models.pml_pad_trun(squeeze(sum(Array(grad_modrrvz),3),3),model.mgrid.npml,-1);

	grad_modrr_stack = (grad_modrrvx_stack .+ grad_modrrvz_stack) .* (0.5)

	for ix=2:size(grad_modrr_stack,2)-1
		for iz=2:size(grad_modrr_stack,1)-1
			@inbounds grad_modrr_stack[iz,ix+1] +=  0.5e0 * grad_modrrvx_stack[iz,ix]
			@inbounds grad_modrr_stack[iz+1,ix] +=  0.5e0 * grad_modrrvz_stack[iz,ix]
		end
	end


	# update gradient model using grad_modtt_stack, grad_modrr_stack
	Models.Seismic_chainrule!(gmodel, model, 
				     vec(Models.χg(grad_modtt_stack, Models.Seismic_get(model,:KI0),1)),
				     vec(Models.χg(grad_modrr_stack, Models.Seismic_get(model,:ρI0),1)),
				     [:χKI, :χρI], 1
				     )
	# output illum
	(illum_flag) && (illum[:] = Models.pml_pad_trun(squeeze(sum(Array(illum_all),3),3),model.mgrid.npml,-1)) 

	# output snaps
	(snaps_flag) && (snaps[:] = Array(snaps_out))

	# update TDout after forming a vector and resampling
	ipropout=0;
	for iprop in 1:npropwav
		if(recv_flags[iprop] ≠ 0)
			ipropout += 1
			Data.TD_resamp!(TDout[ipropout], Data.TD_urpos((Array(recv_out[:,:,iprop,:,:])), recv_fields, tgridmod, acqgeom[iprop],
				acqgeom_urpos[1].nr[1],
				(acqgeom_urpos[1].rz[1], acqgeom_urpos[1].rx[1])
				)) 
		end
	end

	

	# return data depending on the receiver flags
	return TDout, boundary, gmodel, snaps

	# return without resampling for testing
	#return [Data.TD(reshape(recv_out[1+(iprop-1)*nd : iprop*nd],tgridmod.nx,recv_n,src_nseq),
	#		       tgridmod, acqgeom[1]) for iprop in 1:npropwav]
end

# modelling for each processor
function mod_per_proc!(sseqvec::UnitRange{Int64}, pass::Array{Paramss}) 
		       boundary::Array{Array{Array{Float64,3},1},1},
		       grad_modtt::Array{Float64}, grad_modrrvx::Array{Float64}, grad_modrrvz::Array{Float64},
		       recv_out::Array{Float64}, 
		       snaps_out::Array{Float64}, illum_all::Array{Float64}, snaps_flag::Bool,
		       illum_flag::Bool, backprop_flag::Int64, gmodel_flag::Bool, verbose::Bool, 
		       recv_n::Int64, recv_ifields::Vector{Int64}, npropwav::Int64, src_nsmul::Int64, 
		       exmodel::Models.Seismic, model::Models.Seismic, tgridmod::Grid.M1D, 
		       acqgeom_uspos::Array{Acquisition.Geom}, acqgeom_urpos::Array{Acquisition.Geom}, 
		       modttI::Array{Float64}, modrrvx::Array{Float64}, modrrvz::Array{Float64},
		       vpmax::Float64, vpmin::Float64, freqmin::Float64, freqmax::Float64, 
		       abs_trbl::Vector{Symbol}, tsnaps::Vector{Float64}, src_wav_mat::Array{Array{Float64,3},2}, 
		       src_ifields::Vector{Vector{Int64}},
		       born_flag::Bool,
		       born_svalue_stack::Array{Float64}, δmodtt::Array{Float64},
		       δmodrrvx::Array{Float64}, δmodrrvz::Array{Float64})

	#create some aliases
	nx, nz = exmodel.mgrid.nx, exmodel.mgrid.nz
	nxd, nzd = model.mgrid.nx, model.mgrid.nz
	δx, δz = exmodel.mgrid.δx, exmodel.mgrid.δz
	δx24I, δz24I = (δx * 24.0)^(-1.0), (δz * 24.0)^(-1.0)
	δxI, δzI = (δx)^(-1.0), (δz)^(-1.0)
	δt = tgridmod.δx 
	δtI = (δt)^(-1.0)
	tim_nt=tgridmod.nx
	ibx0, ibz0=model.mgrid.npml-2, model.mgrid.npml-2 
	ibx1, ibz1=model.mgrid.nx+model.mgrid.npml+3, model.mgrid.nz+model.mgrid.npml+3
	mesh_x, mesh_z = exmodel.mgrid.x, exmodel.mgrid.z

	# snap save
	itsnaps = [indmin(abs.(tgridmod.x-tsnaps[i])) for i in 1:length(tsnaps)]


	# for snap
	isx0, isz0 =model.mgrid.npml, model.mgrid.npml


	# source_loop
	for (iisseq,isseq) in enumerate(sseqvec)

		if(verbose)
			println("modelling supershot:\t", isseq)
		end
		
	


		# have to make a copy because using DArray inside time loop is very slow
		if(backprop_flag ≠ 0)
			btisseq = view(view(boundary, iisseq)[1],1)[1]
			brisseq = view(view(boundary, iisseq)[1],2)[1]
			bbisseq = view(view(boundary, iisseq)[1],3)[1]
			blisseq = view(view(boundary, iisseq)[1],4)[1]
		end

		if(backprop_flag==-1)
			"initial conditions from boundary for first propagating field only"
			p[:,:,:,1] = view(view(boundary, iisseq)[1],5)[1]
		end

		# time_loop
		"""
		* don't use shared arrays inside this time loop, for speed when using multiple procs
		"""
		for it=1:tim_nt

			advance!(p, pp, ppp, dpdx, dpdz, memory_dp_dx, memory_dp_dz, memory_dvx_dx, memory_dvz_dz,
				   modttI, modrrvx, modrrvz, 
				   δx24I, δz24I, δt, nx, nz,
				   a_x, b_x, k_xI, a_x_half, b_x_half, k_x_halfI, 
				   a_z, b_z, k_zI, a_z_half, b_z_half, k_z_halfI)
		
			# force p[1] on boundaries
			if(backprop_flag==-1) 
				boundary_force(it, view(p,:,:,1,1),
				   blisseq, brisseq,
				   btisseq, bbisseq,
		                        ibx0::Int64, ibz0::Int64, ibx1::Int64, ibz1::Int64)
			end
	 
			add_source!(it, isseq, p, src_wav_mat, src_spray_weights,
				modttI, issmulx1, issmulx2, issmulz1, issmulz2, npropwav, src_ifields, src_nsmul,
				δt::Float64, δxI::Float64, δzI::Float64
				)

			if(born_flag)
				add_born_sources!(iisseq, p, pp, ppp, dpdx, dpdz, δmodtt, modttI, 
				   δmodrrvx, δmodrrvz,
				      δtI, δt, born_svalue_stack, nx, nz, δxI, δzI, δx24I, δz24I)
			end

			# record boundaries
			if(backprop_flag==1)
				boundary_save(tim_nt-it+1, view(p,:,:,1,1),
				   blisseq, brisseq,
				   btisseq, bbisseq,
		                        ibx0::Int64, ibz0::Int64, ibx1::Int64, ibz1::Int64)
			end

			record!(it, iisseq, recv_out, p, rec_interpolate_weights, irecx1, irecx2, irecz1, irecz2, npropwav, recv_ifields, recv_n)

			if(gmodel_flag)
				compute_gradient!(iisseq, p, pp, ppp, dpdx, dpdz,
		      			grad_modtt, grad_modrrvx, grad_modrrvz, nx, nz, δtI)
			end

			if(illum_flag)
				compute_illum!(iisseq, p, illum_all)
			end

			if(snaps_flag)
				itsnapssave = findin(itsnaps,it)
				if(itsnapssave ≠ [])
					snap_save!(iisseq, itsnapssave[1],  p, snaps_out, isx0, isz0)
				end
			end

		end # time_loop
		"now pressure is at [tim_nt], velocities are at [tim_nt-1/2]"	

		"one more propagating step to save pressure at [tim_nt+1] -- for time revarsal"
		advance!(p, pp, ppp, dpdx, dpdz, memory_dp_dx, memory_dp_dz, memory_dvx_dx, memory_dvz_dz,
			   modttI, modrrvx, modrrvz, 
			   δx24I, δz24I, δt, nx, nz,
			   a_x, b_x, k_xI, a_x_half, b_x_half, k_x_halfI, 
			   a_z, b_z, k_zI, a_z_half, b_z_half, k_z_halfI)
		"save last snap of pressure field"
		boundary[iisseq][5][:,:,1] = p[:,:,1,1]

		"one more propagating step to save velocities at [tim_nt+3/2] -- for time reversal"
		advance!(p, pp, ppp, dpdx, dpdz, memory_dp_dx, memory_dp_dz, memory_dvx_dx, memory_dvz_dz,
			   modttI, modrrvx, modrrvz, 
			   δx24I, δz24I, δt, nx, nz,
			   a_x, b_x, k_xI, a_x_half, b_x_half, k_x_halfI, 
			   a_z, b_z, k_zI, a_z_half, b_z_half, k_z_halfI)
		"save last snap of velocity fields with opposite sign for adjoint propagation"
		boundary[iisseq][5][:,:,2:3] = (-1.0 .* p[:,:,2:3,1])

		"combine gradients for each individual isseq"
		if(gmodel_flag)
			"gradient is formed by intergration over time, hence multiply with δt, but why not?"
			"I don't completely understand where the factors δx and δz are coming from..."
			"probably the source term should not be multiplied by δxI and δzI during adjoint propagation"
			grad_modtt[:,:,iisseq] .*=  (δx * δz)
			grad_modrrvx[:,:,iisseq] .*=  (δx * δz)
			grad_modrrvz[:,:,iisseq] .*= (δx * δz)
		end

	end # source_loop
end # mod_per_shot

"""
This routine ABSOLUTELY should not allocate any memory, called inside time loop.
"""
@inline @inbounds @fastmath function add_source!(it::Int64, isseq::Int64, pass::Paramss, pac::Paramc)
	# aliases
	p=pass.p;
	src_wav_mat=pass.src_wav_mat
	issmulx1=pass.issmulx1
	issmulx2=pass.issmulx2
	issmulz1=pass.issmulz1
	issmulz2=pass.issmulz2
	src_spray_weights=pass.src_spray_weights
	modttI=pac.modttI
	"""
	adding source to pressure field at [it] 
	"""
	for ipropwav = 1:pac.npropwav
	for (ifieldsrc, ifield) in enumerate(pac.src_ifields[ipropwav])
	@simd for issmul = 1:pac.src_nsmul
		"""
		use src_wav_mat at [it], i.e., sum of source terms
		until [it-1]
		division of source term with δx and δz (see Jan's fdelmodc manual)
		"""
		source_term = src_wav_mat[ipropwav, isseq][issmul, ifieldsrc, it] * pac.δt * pac.δxI * pac.δzI
		
		"""
		multiplication with modttI
		"""
		p[issmulz1[issmul,ipropwav], issmulx1[issmul,ipropwav],ifield, ipropwav] += 
			source_term * 
			src_spray_weights[1,issmul,ipropwav] * 
			modttI[issmulz1[issmul,ipropwav], issmulx1[issmul,ipropwav]]  
		p[issmulz1[issmul,ipropwav], issmulx2[issmul,ipropwav],ifield, ipropwav] += 
			source_term * 
			src_spray_weights[2,issmul,ipropwav] * 
			modttI[issmulz1[issmul,ipropwav], issmulx2[issmul,ipropwav]]
		p[issmulz2[issmul,ipropwav], issmulx1[issmul,ipropwav],ifield, ipropwav] += 
			source_term * 
			src_spray_weights[3,issmul,ipropwav] * 
			modttI[issmulz2[issmul,ipropwav], issmulx1[issmul,ipropwav]]
		p[issmulz2[issmul,ipropwav], issmulx2[issmul,ipropwav],ifield, ipropwav] += 
			source_term * 
			src_spray_weights[4,issmul,ipropwav] * 
			modttI[issmulz2[issmul,ipropwav], issmulx2[issmul,ipropwav]]
	end
	end
	end
end


"""
This routine ABSOLUTELY should not allocate any memory, called inside time loop.
"""
@inline @inbounds @fastmath function record!(it::Int64, iisseq::Int64, pass::Paramss, pac::Paramc)
	recv_out=pass.recv_out
	p=pass.p
	rec_interpolate_weights=pass.recv_interpolate_weights

	for ipropwav = 1:pac.npropwav
	for (ifieldrecv, ifield) in enumerate(pac.recv_ifields)
	@simd for irec = 1:pac.recv_n
			recv_out[irec,ifieldrecv,ipropwav,it, iisseq] = 
			(
			p[irecz1[irec,ipropwav],irecx1[irec,ipropwav],ifield,ipropwav]*
			rec_interpolate_weights[1,irec,ipropwav]+
			p[irecz1[irec,ipropwav],irecx2[irec,ipropwav],ifield,ipropwav]*
			rec_interpolate_weights[2,irec,ipropwav]+
			p[irecz2[irec,ipropwav],irecx1[irec,ipropwav],ifield,ipropwav]*
			rec_interpolate_weights[3,irec,ipropwav]+
			p[irecz2[irec,ipropwav],irecx2[irec,ipropwav],ifield,ipropwav]*
			rec_interpolate_weights[4,irec,ipropwav]
			)
	end
	end
	end

end

"""
This routine ABSOLUTELY should not allocate any memory, called inside time loop.
"""
@inline @inbounds @fastmath function compute_gradient!(iisseq::Int64, pass::Paramss, pac::Paramc)
	# aliases
	p=pass.p
	pp=pass.pp
	ppp=pass.ppp
	δtI=pac.δtI
	grad_modtt=pass.grad_modtt
	grad_modrrvx=pass.grad_modrrvx
	grad_modrrvz=pass.grad_modrrvz
	dpdx=pass.dpdx
	dpdz=pass.dpdz

	for ix=1:pac.nx
		@simd for iz=1:pac.nz
			# p at [it], pp at [it-1]
			# dpdx and dpdz at [it]
			# gradients w.r.t inverse of modtt, i.e., 1/rho/c^2 
			@inbounds grad_modtt[iz,ix,iisseq] += ((-1.0 * (ppp[iz, ix, 1,1] + p[iz, ix,  1,1] - 2.0 * pp[iz, ix,  1,1]) * δtI * δtI) *  pp[iz,ix,1,2])
			       #(p(:,:,1,1) - pp(:,:,1,1)) * δtI * &
			       #        (p(:,:,1,2) - pp(:,:,1,2)) * δtI
			
			# gradient w.r.t. inverse of rho on vx and vz grids
			@inbounds grad_modrrvx[iz,ix,iisseq] += (- dpdx[iz,ix,1,2]*dpdx[iz,ix,1,1])
			@inbounds grad_modrrvz[iz,ix,iisseq] += (- dpdz[iz,ix,1,2]*dpdz[iz,ix,1,1])
		end
	end
end


@inbounds @fastmath function advance!(pass:Paramss, pac::Paramc)
	# aliases
	p=pass.p; pp=pass.pp; ppp=pass.ppp;
	dpdx=pass.dpdx; dpdz=pass.dpdz;
	memory_dp_dx=pass.memory_dp_dx; memory_dp_dz=pass.memory_dp_dz; memory_dvx_dx=pass.memory_dvx_dx; memory_dvz_dz=pass.memory_dvz_dz
	modttI=pac.modttI; modrrvx=pac.modrrvx; modrrvz=pac.modrrvz
	δx24I=pac.δx24I; δz24I=pac.δz24I; δt=pac.δt
	nx=pac.nx; nz=pac.nz
	a_x=pac.a_x; b_x=pac.b_x; k_xI=pac.k_xI; a_x_half=pac.a_x_half; b_x_half=pac.b_x_half; k_x_halfI=pac.k_x_halfI 
	a_z=pac.a_z; b_z=pac.b_z; k_zI=pac.k_zI; a_z_half=pac.a_z_half; b_z_half=pac.b_z_half; k_z_halfI=pac.k_z_halfI


	ppp .= pp
	pp .=  p

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

	for ipropwav = 1:size(p,4)

		"""
		compute p[it] using p[it-1] and velocity [it-1/2]
		"""
		for ix=3:nx-2
		@simd for iz=3:nz-2
			@inbounds dpdx[iz,ix,2,ipropwav] = (27.e0*p[iz,ix,2,ipropwav]-27.e0*p[iz,ix-1,2,ipropwav]-p[iz,ix+1,2,ipropwav]+p[iz,ix-2,2,ipropwav]) * (δx24I)
			@inbounds dpdz[iz,ix,3,ipropwav] = (27.e0*p[iz,ix,3,ipropwav]-27.e0*p[iz-1,ix,3,ipropwav]-p[iz+1,ix,3,ipropwav]+p[iz-2,ix,3,ipropwav]) * (δz24I)
			@inbounds memory_dvx_dx[iz,ix,ipropwav] = b_x[ix] * memory_dvx_dx[iz,ix,ipropwav] + a_x[ix] * dpdx[iz,ix,2,ipropwav] # pml 
			@inbounds memory_dvz_dz[iz,ix,ipropwav] = b_z[iz] * memory_dvz_dz[iz,ix,ipropwav] + a_z[iz] * dpdz[iz,ix,3,ipropwav] # pml
			@inbounds dpdx[iz,ix,2,ipropwav] = dpdx[iz,ix,2,ipropwav] * k_xI[ix] + memory_dvx_dx[iz,ix,ipropwav] # pml
			@inbounds dpdz[iz,ix,3,ipropwav] = dpdz[iz,ix,3,ipropwav] * k_zI[iz] + memory_dvz_dz[iz,ix,ipropwav] # pml

			"""
			compute pressure at [it] using p at [it-1] and dvxdx
			and dvzdz at [it-1/2]
			"""
			@inbounds p[iz,ix,1,ipropwav] += (modttI[iz,ix] * (dpdx[iz,ix,2,ipropwav] + dpdz[iz,ix,3,ipropwav])) * δt #* boundary_p(iz,ix)
		end
		end
	end
end

@inbounds @fastmath function update_dpdx!(p, dpdx, δx24I, memory_dp_dx, b_x_half, a_x_half, k_x_halfI, nx, nz)
	for ipropwav = 1:size(p,4)
		for ix = 3:nx-2
		@simd for iz = 3:nz-2
			@inbounds dpdx[iz,ix,1,ipropwav] = (27.e0*p[iz,ix+1,1,ipropwav]-27.e0*p[iz,ix,1,ipropwav]-p[iz,ix+2,1,ipropwav]+p[iz,ix-1,1,ipropwav]) * (δx24I)
			@inbounds memory_dp_dx[iz,ix,ipropwav] = b_x_half[ix] * memory_dp_dx[iz,ix,ipropwav] + a_x_half[ix] * dpdx[iz,ix,1,ipropwav] # pml
			@inbounds dpdx[iz,ix,1,ipropwav] = dpdx[iz,ix,1,ipropwav] * k_x_halfI[ix] + memory_dp_dx[iz,ix,ipropwav] # pml
		end
		end
	end
end

@inbounds @fastmath function update_dpdz!(p, dpdz, δz24I, memory_dp_dz, b_z_half, a_z_half, k_z_halfI, nx, nz)
	for ipropwav = 1:size(p,4)
		for ix = 3:nx-2
		@simd for iz = 3:nz-2
			@inbounds dpdz[iz,ix,1,ipropwav] = (27.e0*p[iz+1,ix,1,ipropwav]-27.e0*p[iz,ix,1,ipropwav]-p[iz+2,ix,1,ipropwav]+p[iz-1,ix,1,ipropwav]) * (δz24I)
			@inbounds memory_dp_dz[iz,ix,ipropwav] = b_z_half[iz] * memory_dp_dz[iz,ix,ipropwav] + a_z_half[iz] * dpdz[iz,ix,1,ipropwav] # pml
			@inbounds dpdz[iz,ix,1,ipropwav] = dpdz[iz,ix,1,ipropwav] * k_z_halfI[iz] + memory_dp_dz[iz,ix,ipropwav] # pml
		end
		end
	end
end

@inbounds @fastmath function update_vx!(p, dpdx, δt, modrrvx,  nx, nz)
	for ipropwav = 1:size(p,4)
		for ix=3:nx-2
		@simd for iz=3:nz-2
			@inbounds p[iz,ix,2,ipropwav] += (dpdx[iz,ix,1,ipropwav]) * δt * modrrvx[iz,ix] #* boundary_vx(iz,ix)
		end
		end
	end
end

@inbounds @fastmath function update_vz!(p, dpdz, δt,  modrrvz, nx, nz)
	for ipropwav = 1:size(p,4)
		for ix=3:nx-2
		@simd for iz=3:nz-2
			@inbounds p[iz,ix,3,ipropwav] +=  (dpdz[iz,ix,1,ipropwav]) * δt * modrrvz[iz,ix] #* boundary_vz(iz,ix)
		end
		end
	end
end


"""
Need illumination to estimate the approximate diagonal of Hessian
"""
@inline function compute_illum!(iisseq::Int64, p::Array{Float64}, illum_all::Array{Float64})
	for ix=1:size(illum_all,2)
		@simd for iz=1:size(illum_all,1)
			# saving illumination to be used as preconditioner 
			illum_all[iz,ix,iisseq] += (p[iz,ix,1,1] * p[iz,ix,1,1])
		end
	end
end

function add_born_sources!(iisseq::Int64, pass::Paramss, pac::Paramc)

 	born_svalue_stack=pass.born_svalue_stack
	δx24I=pac.δx24I; δz24I=pac.δz24I; 
	δxI=pac.δxI; δzI=pac.δzI; 
	δt=pac.δt
	δtI=pac.δtI
	δmodtt=pac.δmodtt; modttI=pac.modttI;
	δmodrrvx=pac.δmodrrvx; δmodrrvz=pac.δmodrrvz
	nx=pac.nx; nz=pac.nz
	p=pass.p; pp=pass.pp; ppp=pass.ppp;
	dpdx=pass.dpdx; dpdz=pass.dpdz;

	# secondary sources for Born modeling
	# adding born sources from pressure(:,:,1) to pressure(:,:,2)
	# upto until [it-2]
	# lambdaI scatterrer source term at [it-1]
	# p is at [it], pp is at [it-1], ppp is at [it-2]
	# dpdx is at [it-1] and dpdz is at [it-1]
	# modrrvx scatterrer source term at [it-1]
	# modrrvz scatterrer source term at [it-1]

	for ix=3:nx-2
		@simd for iz=3:nz-2

			@inbounds born_svalue_stack[iz,ix, iisseq] += 
				δt * ((-1.0 * (ppp[iz, ix, 1,1] + p[iz, ix,  1,1] - 2.0 * pp[iz, ix,  1,1]) * δmodtt[iz, ix] * δtI * δtI) + (
			  (27.e0*dpdx[iz,ix,1,1] * δmodrrvx[iz,ix] -27.e0*dpdx[iz,ix-1,1,1] * δmodrrvx[iz,ix-1] -dpdx[iz,ix+1,1,1] * δmodrrvx[iz,ix+1] +dpdx[iz,ix-2,1,1] * δmodrrvx[iz,ix-2] ) * (δx24I)) + (
			  (27.e0*dpdz[iz,ix,1,1] * δmodrrvz[iz,ix] -27.e0*dpdz[iz-1,ix,1,1] * δmodrrvz[iz-1,ix] -dpdz[iz+1,ix,1,1] * δmodrrvz[iz+1,ix] +dpdz[iz-2,ix,1,1] * δmodrrvz[iz-2,ix] ) * (δz24I)) ) 

			@inbounds p[iz,ix,1,2] += born_svalue_stack[iz,ix, iisseq] * δt * modttI[iz,ix] #* δxI * δzI 
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


@inbounds function boundary_force(it::Int64, p::AbstractArray{Float64}, 
		     bl::AbstractArray{Float64}, br::AbstractArray{Float64},
		     bt::AbstractArray{Float64}, bb::AbstractArray{Float64},
		     ibx0::Int64, ibz0::Int64, ibx1::Int64, ibz1::Int64)

	nlayer=3
	for ix=1:nlayer
		@simd for iz=1:ibz1-ibz0+1
			p[ibz0+iz-1, ibx0+ix-1] = bl[iz,ix, it]
			p[ibz0+iz-1, ibx1-ix+1] = br[iz,ix, it]
		end
	end

	@simd for ix=1:ibx1-ibx0+1
		for iz=1:nlayer
			p[ibz0+iz-1,ibx0+ix-1] = bt[iz,ix, it]
			p[ibz1-iz+1,ibx0+ix-1] = bb[iz,ix, it]
		end
	end

end

@inbounds function boundary_save(it::Int64, p::AbstractArray{Float64}, 
		     bl::AbstractArray{Float64}, br::AbstractArray{Float64},
		     bt::AbstractArray{Float64}, bb::AbstractArray{Float64},
		     ibx0::Int64, ibz0::Int64, ibx1::Int64, ibz1::Int64)

	nlayer=3
	for ix=1:nlayer
		@simd for iz=1:ibz1-ibz0+1
			bl[iz,ix,it] = p[ibz0+iz-1, ibx0+ix-1] 
			br[iz,ix,it] = p[ibz0+iz-1, ibx1-ix+1]
		end
	end

	@simd for ix=1:ibx1-ibx0+1
		for iz=1:nlayer
			bt[iz,ix,it] = p[ibz0+iz-1,ibx0+ix-1]
			bb[iz,ix,it] = p[ibz1-iz+1,ibx0+ix-1]
		end
	end

end

function snap_save!(iisseq::Int64, itsnapssave::Int64, 
		    p::AbstractArray{Float64}, snap::AbstractArray{Float64}, isx0::Int64, isz0::Int64)
	for ix=1:size(snap,2)
		for iz=1:size(snap,1)
			snap[iz,ix, itsnapssave, iisseq] = p[isz0+iz, isx0+ix,1,1]
		end
	end
end

function fill_src_wav_mat!(src_wav_mat::Array{Array{Float64, 3}, 2}, acqsrc_uspos::Array{Acquisition.Src},src_flags::Vector{Int64})

	npropwav, src_nseq = size(src_wav_mat)
	δt = acqsrc_uspos[1].tgrid.δx
	for isseq=1:src_nseq
	for ipropwav=1:npropwav
		src_nsmul, src_nfield, tim_nt = size(src_wav_mat[ipropwav,isseq])
		for ifield=1:src_nfield, issmul=1:src_nsmul
			src_tim_nt = acqsrc_uspos[ipropwav].tgrid.nx;
			if(src_flags[ipropwav] == 0)
				nothing # just put zeros, no sources added
			elseif(src_flags[ipropwav] == 1)
				"ϕ[t] = s[t]"
				for it=1:src_tim_nt
					source_term = acqsrc_uspos[ipropwav].wav[isseq,ifield][it,issmul]
					src_wav_mat[ipropwav, isseq][issmul,ifield,it] = source_term
				end
			elseif(src_flags[ipropwav] == 2)
				"ϕ[t] = ∑₁ᵗ⁻¹ s[t]"
				source_term_stack = 0.0;
				if(ifield == 1)
					for it=1:src_tim_nt-1
						source_term_stack += (acqsrc_uspos[ipropwav].wav[isseq,ifield][it,issmul] .* δt)
						src_wav_mat[ipropwav, isseq][issmul,ifield,it+1] = source_term_stack
					end
				else
					for it=2:src_tim_nt-1
						source_term_stack += (((acqsrc_uspos[ipropwav].wav[isseq,ifield][it,issmul] .* δt) +
						   (acqsrc_uspos[ipropwav].wav[isseq,ifield][it-1,issmul] .* δt)) * 0.5)
						src_wav_mat[ipropwav, isseq][issmul,ifield,it+1] = source_term_stack
					end

				end
				if(tim_nt > src_tim_nt)
					src_wav_mat[ipropwav, isseq][issmul,ifield,src_tim_nt+1:end] = src_wav_mat[ipropwav, isseq][issmul,ifield,src_tim_nt]
				end
			elseif(src_flags[ipropwav] == 3)
				"use this to add source sink: need during adjoint propagation from boundary"
				"multiplication with -1 for subtraction"
				"time reversal"
				"as the source wavelet has to be subtracted before the propagation step, I shift here by one sample"
				"ϕ[t] = "
				source_term_stack = 0.0;
				for it=1:src_tim_nt-1
					source_term_stack += (acqsrc_uspos[ipropwav].wav[isseq,ifield][it,issmul] .* δt)
					src_wav_mat[ipropwav, isseq][issmul,ifield,tim_nt-it+1] = -1.0 * source_term_stack
				end
				if(tim_nt > src_tim_nt)
					tim_nt_diff = tim_nt-src_tim_nt
					src_wav_mat[ipropwav, isseq][issmul,ifield,1:tim_nt_diff+1] = src_wav_mat[ipropwav, isseq][issmul,ifield,tim_nt_diff+2]
				end
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

end # module
