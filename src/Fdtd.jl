__precompile__()

module Fdtd

import SIT.Grid
import SIT.Interpolation
import SIT.Models
import SIT.F90libs
import SIT.Acquisition
import SIT.Data
import SIT.Gallery
import SIT.DSP

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
* `recv_nfield::Int64=1` : multi-component receiver flag; types fields the receivers record (to be changed later)
* `backprop_flag::Bool=Int64` : save final state variables and the boundary conditions for later use
  * `=1` save boundary values in `boundary` and final values in `initp`
  * `=-1` force boundary values that are in `boundary` and use initial values in `initp` for back propagation
* `abs_trbl::Vector{Symbol}=[:top, :bottom, :right, :left]` : use absorbing PML boundary conditions or not
  * `=[:top, :bottom]` apply PML conditions only at the top and bottom of the model 
  * `=[:bottom, :right, :left]` top is reflecting
* `born_flag::Bool=false` : do only Born modeling instead of full wavefield modelling (to be updated soon)
* `gmodel_flag=false` : flag that is used to output gradient; there should be atleast two propagating wavefields in order to do so: 1) forward wavefield and 2) adjoint wavefield
* `illum_flag::Bool=false` : flag to output wavefield energy or source illumination; it can be used as preconditioner during inversion
* `tsnaps::Vector{Float64}=fill(0.5*(tgridmod.x[end]+tgridmod.x[1]),1)` : store snaps at these modelling times
* `verbose::Bool=false` : verbose flag

# Keyword Arguments that are modified by the method (some of them are returned as well)

* `gmodel::Models.Seismic=Models.Seismic_zeros(model.mgrid)` : gradient model modified only if `gmodel_flag`
* `TDout::Vector{Data.TD}=[Data.TD_zeros(recv_nfield,tgridmod,acqgeom[ip]) for ip in 1:length(findn(recv_flags))]`
* `illum_out::Array{Float64,2}=zeros(model.mgrid.nz, model.mgrid.nx)` : source energy if `illum_flag`
* `boundary::Array{Array{Float64,4},1}` : stored boundary values for first propagating wavefield 
* `initp::Array{Float64,4}` : stored final values for the first propagating wavefield
* `snaps::Array{Float64,4}=zeros(model.mgrid.nz,model.mgrid.nx,length(tsnaps),acqgeom[1].nss)` :snapshots saved at `tsnaps`

# Return (in order)

* modelled data for each propagating wavefield as `Vector{TD}`
* stored boundary values of the first propagating wavefield as `Array{Array{Float64,4},1}` (use for backpropagation)
* final conditions of the first propagating wavefield as `Array{Float64,4}` (use for back propagation)
* gradient model as `Seismic`
* stored snaps shots at tsnaps as Array{Float64,4} 

# Example

```julia
records, boundary, initp, gmodel, snaps = mod(acqgeom=acqgeom, acqsrc=acqsrc, model=model, tgridmod=tgridmod);
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
	acqgeom::Vector{Acquisition.Geom} = fill(Gallery.Geom(:acou_homo1),npropwav),
	acqsrc::Array{Acquisition.Src} = fill(Gallery.Src(:acou_homo1),npropwav),
	src_flags::Vector{Int64}=fill(2,npropwav), 
	recv_flags::Vector{Int64}=fill(1,npropwav),
	recv_nfield::Int64=1, 
	TDout::Vector{Data.TD}=[Data.TD_zeros(recv_nfield,tgridmod,acqgeom[ip]) for ip in 1:length(findn(recv_flags)) ],
	backprop_flag::Int64=0,  
	boundary::Array{Array{Float64,4},1}=[
						zeros(3,model.mgrid.nx+6,tgridmod.nx,acqgeom[1].nss),
						zeros(model.mgrid.nz+6,3,tgridmod.nx,acqgeom[1].nss),
						zeros(3,model.mgrid.nx+6,tgridmod.nx,acqgeom[1].nss),
						zeros(model.mgrid.nz+6,3,tgridmod.nx,acqgeom[1].nss)],
	initp::Array{Float64,4}=
	        zeros(model.mgrid.nz+2*model.mgrid.npml,model.mgrid.nx+2*model.mgrid.npml,3,acqgeom[1].nss),
	gmodel_flag=false,
	gmodel::Models.Seismic=Models.Seismic_zeros(model.mgrid),
	illum_flag::Bool=false,
	illum_out::Array{Float64,2}=zeros(model.mgrid.nz, model.mgrid.nx),
	tsnaps::Vector{Float64}=fill(0.5*(tgridmod.x[end]+tgridmod.x[1]),1),
	snaps::Array{Float64,4}=zeros(model.mgrid.nz,model.mgrid.nx,length(tsnaps),acqgeom[1].nss),
	verbose::Bool=false
	)
	

	# check size TDout
	length(TDout) == length(findn(recv_flags)) ? nothing : error("TDout dimension")

	maximum(tgridmod.x) < maximum(acqsrc[1].tgrid.x) ? error("modeling time is less than source time") : nothing

	any([getfield(TDout[ip],:tgrid).δx < tgridmod.δx for ip=1:length(TDout)]) ? 
				error("output time grid sampling finer than modeling") : nothing
	any([maximum(getfield(TDout[ip],:tgrid).x) > maximum(tgridmod.x) for ip=1:length(TDout)]) ?
	                        error("output time > modeling time") : nothing

	# find maximum and minimum frequencies in the source wavelets
	freqmin = DSP.findfreq(acqsrc[1].wav[1,1][:,1],acqsrc[1].tgrid,attrib=:min) 
	freqmax = DSP.findfreq(acqsrc[1].wav[1,1][:,1],acqsrc[1].tgrid,attrib=:max) 

	# minimum and maximum velocities
	vpmin = minimum(broadcast(minimum,[model.vp0, model_pert.vp0]))
	vpmax = maximum(broadcast(maximum,[model.vp0, model_pert.vp0]))
	verbose ? println("minimum and maximum velocities:\t",vpmin,"\t",vpmax) : nothing


	check_fd_stability(vpmin, vpmax, model.mgrid.δx, model.mgrid.δz, freqmin, freqmax, tgridmod.δx, verbose)

	#! no modeling if source wavelet is zero
	#if(maxval(abs(src_wavelets)) .lt. tiny(rzero_de)) 
	#        return
	#endif

	# all the propagating wavefield should have same sources, receivers, ? check that?

	# check dimension of model and model_pert

	# check if all sources are receivers are inside model

	length(acqgeom) != npropwav ? error("acqgeom size") : nothing
	length(acqsrc) != npropwav ? error("acqsrc size") : nothing
	any([getfield(acqgeom[ip],:nss) != getfield(acqsrc[ip],:nss) for ip=1:npropwav])  ? error("different supersources") : nothing
	any([getfield(acqgeom[ip],:ns) != getfield(acqsrc[ip],:ns) for ip=1:npropwav])  ? error("different sources") : nothing

	# necessary that nss and nfield should be same for all nprop
	src_nseq = acqgeom[1].nss;
	src_nfield = acqsrc[1].nfield
	fill(src_nseq, npropwav) != [getfield(acqgeom[ip],:nss) for ip=1:npropwav] ? error("different supersources") : nothing
	fill(src_nfield, npropwav) != [getfield(acqsrc[ip],:nfield) for ip=1:npropwav] ? error("different fields") : nothing

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


	# extend models in the PML layers
	exmodel = Models.Seismic_pad_trun(model);
	exmodel_pert = Models.Seismic_pad_trun(model_pert);

	"""
	create some aliases
	"""
	nx = exmodel.mgrid.nx
	nz = exmodel.mgrid.nz
	nxd = model.mgrid.nx
	nzd = model.mgrid.nz
	δx = exmodel.mgrid.δx
	δz = exmodel.mgrid.δz
	δx24I = (exmodel.mgrid.δx * 24.0)^(-1.0)
	δxI = (exmodel.mgrid.δx)^(-1.0)
	δz24I = (exmodel.mgrid.δz * 24.0)^(-1.0)
	δzI = (exmodel.mgrid.δz)^(-1.0)
	δt = tgridmod.δx 
	δtI = (tgridmod.δx)^(-1.0)
	tim_nt=tgridmod.nx
	mesh_x=exmodel.mgrid.x
	mesh_z=exmodel.mgrid.z

	# records_output, shared array among different procs
	recv_out = SharedArray(Float64, (tim_nt,recv_n,recv_nfield,npropwav,src_nseq))

	# snap save
	itsnaps = [indmin(abs(tgridmod.x-tsnaps[i])) for i in 1:length(tsnaps)]

	# gradient outputs
	grad_modtt = SharedArray(Float64, (nz, nx, src_nseq))
	grad_modtt .= 0.0
	grad_modrrvx = SharedArray(Float64, (nz, nx, src_nseq))
	grad_modrrvx .= 0.0
	grad_modrrvz = SharedArray(Float64, (nz, nx, src_nseq))
	grad_modrrvz .= 0.0


	# pml_variables
	a_x, b_x, k_xI, a_x_half, b_x_half, k_x_halfI = pml_variables(nx, δt, δx, model.mgrid.npml-5, vpmax, vpmin, freqmin, freqmax, 
							       [any(abs_trbl .== :left), any(abs_trbl .== :right)])
	a_z, b_z, k_zI, a_z_half, b_z_half, k_z_halfI = pml_variables(nz, δt, δz, model.mgrid.npml-5, vpmax, vpmin, freqmin, freqmax,
							       [any(abs_trbl .== :top), any(abs_trbl .== :bottom)])

	# boundary coordinates
	ibx0=model.mgrid.npml-2; ibx1=model.mgrid.nx+model.mgrid.npml+3
	ibz0=model.mgrid.npml-2; ibz1=model.mgrid.nz+model.mgrid.npml+3

	# for snap
	isx0=model.mgrid.npml; isz0=model.mgrid.npml


	"saving boundary values for all super sources"
	if(backprop_flag == 1)
		btout = SharedArray(Float64,3,ibx1-ibx0+1,tim_nt,src_nseq)
		brout = SharedArray(Float64,ibz1-ibz0+1,3,tim_nt,src_nseq)
		bbout = SharedArray(Float64,3,ibx1-ibx0+1,tim_nt,src_nseq)
		blout = SharedArray(Float64,ibz1-ibz0+1,3,tim_nt,src_nseq)
	end

	"p_out will save final snaps that can be used for time reversal"
	p_out = SharedArray(Float64, size(initp))


	if(illum_flag)
		"saving illum"
		illum_all = SharedArray(Float64, (nz, nx, src_nseq))
		illum_all .= 0.0
	end

	snapsout = SharedArray(Float64, size(snaps))

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
		born_svalue_stack = zeros(δmodtt)
	end


	# source wavelets
	src_wav_mat = zeros(src_nsmul,src_nfield,npropwav,tim_nt,src_nseq)
	fill_src_wav_mat!(src_wav_mat, acqsrc_uspos, src_flags)

	# source_loop
	@sync @parallel for isseq = 1:src_nseq

		if(verbose)
			println("modelling supershot:\t", isseq)
		end
		
		"allocate receiver matrix per sequential source"
		rec_mat=zeros(recv_n, recv_nfield, npropwav, tim_nt)


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

		"calling first time will be slow -- dummy advance with all zeros"
		advance!(p, pp, ppp, dpdx, dpdz, memory_dp_dx, memory_dp_dz, memory_dvx_dx, memory_dvz_dz,
				   modttI, modrrvx, modrrvz, 
				   δx24I, δz24I, δt, nx, nz,
				   a_x, b_x, k_xI, a_x_half, b_x_half, k_x_halfI, 
				   a_z, b_z, k_zI, a_z_half, b_z_half, k_z_halfI)

		gradisseq_modtt=zeros(nz,nx)
		gradisseq_modrrvx=zeros(gradisseq_modtt)
		gradisseq_modrrvz=zeros(gradisseq_modtt)

		if(illum_flag)
			illumisseq_out=zeros(nz,nx)
		end


		if(backprop_flag == 1)
			"saving boundary values per issseq"
			brisseq = zeros(ibz1-ibz0+1,3,tim_nt)
			btisseq = zeros(3,ibx1-ibx0+1,tim_nt)
			bbisseq = zeros(3,ibx1-ibx0+1,tim_nt)
			blisseq = zeros(ibz1-ibz0+1,3,tim_nt)
		elseif(backprop_flag==-1)
			"initial conditions from initp for first propagating field only"
			p[:,:,:,1].=initp[:,:,:,isseq]
			btisseq = view(view(boundary,1)[1],:,:,:,isseq)
			brisseq = view(view(boundary,2)[1],:,:,:,isseq)
			bbisseq = view(view(boundary,3)[1],:,:,:,isseq)
			blisseq = view(view(boundary,4)[1],:,:,:,isseq)
		end

		# private variable for snapshot
		snapsisseq = zeros(nzd, nxd, length(itsnaps))

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
				boundary_force(view(p,:,:,1,1),
				   view(blisseq,:,:,it), view(brisseq,:,:,it),
				   view(btisseq,:,:,it), view(bbisseq,:,:,it),
		                        ibx0::Int64, ibz0::Int64, ibx1::Int64, ibz1::Int64)
			end
	 
			add_source!(p, view(src_wav_mat,:,:,:,it,isseq), src_spray_weights,
				modttI, issmulx1, issmulx2, issmulz1, issmulz2, npropwav, src_nfield, src_nsmul,
				δt::Float64, δxI::Float64, δzI::Float64
				)


			# record boundaries
			if(backprop_flag==1)
				boundary_save(view(p,:,:,1,1),
				   view(blisseq,:,:,it), view(brisseq,:,:,it),
				   view(btisseq,:,:,it), view(bbisseq,:,:,it),
		                        ibx0::Int64, ibz0::Int64, ibx1::Int64, ibz1::Int64)
			end

			record!(view(rec_mat, :, :, :, it), p, rec_interpolate_weights, irecx1, irecx2, irecz1, irecz2, npropwav, recv_nfield, recv_n)

			if(gmodel_flag)
				compute_gradient!(p, pp, ppp, dpdx, dpdz,
		      			gradisseq_modtt, gradisseq_modrrvx, gradisseq_modrrvz,nx,nz,δtI)
			end

			if(illum_flag)
				compute_illum!(p, illumisseq_out)
			end

			if(findin(itsnaps,it)!=[])
				snap_save!(view(p,:,:,1,1), view(snapsisseq,:,:,findin(itsnaps,it)[1]), isx0, isz0)
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
		p_out[:,:,1,isseq] .= p[:,:,1,1]

		"one more propagating step to save velocities at [tim_nt+3/2] -- for time reversal"
		advance!(p, pp, ppp, dpdx, dpdz, memory_dp_dx, memory_dp_dz, memory_dvx_dx, memory_dvz_dz,
			   modttI, modrrvx, modrrvz, 
			   δx24I, δz24I, δt, nx, nz,
			   a_x, b_x, k_xI, a_x_half, b_x_half, k_x_halfI, 
			   a_z, b_z, k_zI, a_z_half, b_z_half, k_z_halfI)
		"save last snap of velocity fields with opposite sign for adjoint propagation"
		p_out[:,:,2:3,isseq] .= -1.0 * p[:,:,2:3,1]

		"rec_mat to recv_out"
		for ipropwav=1:npropwav
			for ifield=1:recv_nfield
				for irec=1:recv_n
					for it=1:tim_nt
						@inbounds recv_out[it,irec,ifield,ipropwav,isseq]=rec_mat[irec,ifield,ipropwav,it]
					end
				end
			end
		end

		"combine gradients for each individual isseq"
		if(gmodel_flag)
			"gradient is formed by intergration over time, hence multiply with δt, but why not?"
			"I don't completely understand where the factors δx and δz are coming from..."
			"probably the source term should not be multiplied by δxI and δzI during adjoint propagation"
			grad_modtt[:,:,isseq] .=  gradisseq_modtt .* (δx * δz)
			grad_modrrvx[:,:,isseq] .=  gradisseq_modrrvx .* (δx * δz)
			grad_modrrvz[:,:,isseq] .= gradisseq_modrrvz .* (δx * δz)
		end
		# update boundary_out after time reversal
		if(backprop_flag ==1)
			btout[:,:,:,isseq] .= flipdim(btisseq,3)
			brout[:,:,:,isseq] .= flipdim(brisseq,3)
			bbout[:,:,:,isseq] .= flipdim(bbisseq,3)
			blout[:,:,:,isseq] .= flipdim(blisseq,3)
		end

		# update snapshot
		snapsout[:,:,:,isseq] .= snapsisseq[:,:,:]

		if(illum_flag)
			# output illum_out
			illum_all[:,:,isseq] .= illumisseq_out
		end

	end # source_loop

	"check if recv_out is zeros"
	isapprox(maximum(abs(recv_out)),0.0) && warn("recv_out are zeros")

	"summing gradient over all the sources after truncation in PML"
	grad_modtt_stack = Models.pad_trun(squeeze(sum(Array(grad_modtt),3),3),model.mgrid.npml,-1);
	grad_modrrvx_stack = Models.pad_trun(squeeze(sum(Array(grad_modrrvx),3),3),model.mgrid.npml,-1);
	grad_modrrvz_stack = Models.pad_trun(squeeze(sum(Array(grad_modrrvz),3),3),model.mgrid.npml,-1);

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
	if(illum_flag)
		# output illum_out
		illum_out[:,:] = Models.pad_trun(squeeze(sum(Array(illum_all),3),3),model.mgrid.npml,-1)
	end

	snaps[:,:,:] = Array(snapsout)

	# update TDout after forming a vector and resampling
	ipropout=0;
	for iprop in 1:npropwav
		if(recv_flags[iprop] ≠ 0)
			ipropout += 1
			Data.TD_resamp!(
			       Data.TD_urpos(
			       recv_out[:,:,:,iprop,:],
				recv_nfield,
				tgridmod, acqgeom[iprop],
				acqgeom_urpos[1].nr[1],
				(acqgeom_urpos[1].rz[1], acqgeom_urpos[1].rx[1])
				), TDout[ipropout]) 
		end
	end

	
	# update initp and boundary outputs
	if(backprop_flag==1)
		initp .= Array(p_out) 
		boundary[:]=[Array(btout), Array(brout), Array(bbout), Array(blout)]
	end


	# return data depending on the receiver flags
	return TDout, boundary, initp, gmodel, snaps

	# return without resampling for testing
	#return [Data.TD(reshape(recv_out[1+(iprop-1)*nd : iprop*nd],tgridmod.nx,recv_n,src_nseq),
	#		       tgridmod, acqgeom[1]) for iprop in 1:npropwav]
end



@inline @inbounds @fastmath function add_source!( 
				p::Array{Float64}, src_wav_mat::AbstractArray{Float64},
				src_spray_weights::Array{Float64},
				modttI::Array{Float64},
				issmulx1::Array{Int64}, issmulx2::Array{Int64},
				issmulz1::Array{Int64}, issmulz2::Array{Int64},
				npropwav::Int64, src_nfield::Int64, src_nsmul::Int64,
				δt::Float64, δxI::Float64, δzI::Float64
				)
	"""
	adding source to pressure field at [it] 
	"""
	for ipropwav = 1:npropwav
	for ifield = 1:src_nfield
	@simd for issmul = 1:src_nsmul
		"""
		use src_wav_mat at [it], i.e., sum of source terms
		until [it-1]
		division of source term with δx and δz (see Jan's fdelmodc manual)
		"""
		source_term = src_wav_mat[issmul, ifield, ipropwav] * δt * δxI * δzI
		
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


@inline @inbounds @fastmath function record!(
		 rec_mat::AbstractArray{Float64}, 
		 p::Array{Float64}, 
		 rec_interpolate_weights::Array{Float64},
		 irecx1::Array{Int64}, irecx2::Array{Int64},
		 irecz1::Array{Int64}, irecz2::Array{Int64},
		 npropwav::Int64, recv_nfield::Int64, 
		 recv_n::Int64)

	for ipropwav = 1:npropwav
	for ifield = 1:recv_nfield
	@simd for irec = 1:recv_n
			rec_mat[irec,ifield,ipropwav] = 
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

@inline @inbounds @fastmath function compute_gradient!(p::Array{Float64}, pp::Array{Float64}, ppp::Array{Float64},
			   dpdx::Array{Float64}, dpdz::Array{Float64},
			   gradisseq_modtt::Array{Float64}, gradisseq_modrrvx::Array{Float64},
			   gradisseq_modrrvz::Array{Float64},
			   nx::Int64, nz::Int64, δtI::Float64 
			  )

	for ix=1:nx
		@simd for iz=1:nz
			# p at [it], pp at [it-1]
			# dpdx and dpdz at [it]
			# gradients w.r.t inverse of modtt, i.e., 1/rho/c^2 
			@inbounds gradisseq_modtt[iz,ix] += ((-1.0 * (ppp[iz, ix, 1,1] + p[iz, ix,  1,1] - 2.0 * pp[iz, ix,  1,1]) * δtI * δtI) *  pp[iz,ix,1,2])
			       #(p(:,:,1,1) - pp(:,:,1,1)) * δtI * &
			       #        (p(:,:,1,2) - pp(:,:,1,2)) * δtI
			
			# gradient w.r.t. inverse of rho on vx and vz grids
			@inbounds gradisseq_modrrvx[iz,ix] += (- dpdx[iz,ix,1,2]*dpdx[iz,ix,1,1])
			@inbounds gradisseq_modrrvz[iz,ix] += (- dpdz[iz,ix,1,2]*dpdz[iz,ix,1,1])
		end
	end
end


@inline @inbounds @fastmath function advance!(
		  p::Array{Float64}, pp::Array{Float64}, ppp::Array{Float64},  dpdx::Array{Float64}, dpdz::Array{Float64},  
		  memory_dp_dx::Array{Float64}, memory_dp_dz::Array{Float64}, memory_dvx_dx::Array{Float64}, memory_dvz_dz::Array{Float64},
		  modttI::Array{Float64}, modrrvx::Array{Float64}, modrrvz::Array{Float64},
		  δx24I::Float64, δz24I::Float64, δt::Float64,
		  nx::Int64, nz::Int64,
		  a_x::Vector{Float64}, b_x::Vector{Float64}, k_xI::Vector{Float64}, a_x_half::Vector{Float64}, b_x_half::Vector{Float64}, k_x_halfI::Vector{Float64}, 
		  a_z::Vector{Float64}, b_z::Vector{Float64}, k_zI::Vector{Float64}, a_z_half::Vector{Float64}, b_z_half::Vector{Float64}, k_z_halfI::Vector{Float64}
		 )



	ppp .= pp
	pp .=  p

	"""
	compute dpdx and dpdz at [it-1] for all propagating fields
	"""
	for ipropwav = 1:size(p,4)
		for ix = 3:nx-2
		@simd for iz = 3:nz-2
			@inbounds dpdx[iz,ix,1,ipropwav] = (27.e0*p[iz,ix+1,1,ipropwav]-27.e0*p[iz,ix,1,ipropwav]-p[iz,ix+2,1,ipropwav]+p[iz,ix-1,1,ipropwav]) * (δx24I)
			@inbounds memory_dp_dx[iz,ix,ipropwav] = b_x_half[ix] * memory_dp_dx[iz,ix,ipropwav] + a_x_half[ix] * dpdx[iz,ix,1,ipropwav] # pml
			@inbounds dpdx[iz,ix,1,ipropwav] = dpdx[iz,ix,1,ipropwav] * k_x_halfI[ix] + memory_dp_dx[iz,ix,ipropwav] # pml
		end
		end

		for ix = 3:nx-2
		@simd for iz = 3:nz-2
			@inbounds dpdz[iz,ix,1,ipropwav] = (27.e0*p[iz+1,ix,1,ipropwav]-27.e0*p[iz,ix,1,ipropwav]-p[iz+2,ix,1,ipropwav]+p[iz-1,ix,1,ipropwav]) * (δz24I)
			@inbounds memory_dp_dz[iz,ix,ipropwav] = b_z_half[iz] * memory_dp_dz[iz,ix,ipropwav] + a_z_half[iz] * dpdz[iz,ix,1,ipropwav] # pml
			@inbounds dpdz[iz,ix,1,ipropwav] = dpdz[iz,ix,1,ipropwav] * k_z_halfI[iz] + memory_dp_dz[iz,ix,ipropwav] # pml
		end
		end


		"""
		update velocity at [it-1/2] using 
		velocity at [it-3/2] and dpdx and dpdz at [it-1] 
		"""
		for ix=3:nx-2
		@simd for iz=3:nz-2
			@inbounds p[iz,ix,2,ipropwav] += (dpdx[iz,ix,1,ipropwav]) * δt * modrrvx[iz,ix] #* boundary_vx(iz,ix)
		end
		end

		for ix=3:nx-2
		@simd for iz=3:nz-2
			@inbounds p[iz,ix,3,ipropwav] +=  (dpdz[iz,ix,1,ipropwav]) * δt * modrrvz[iz,ix] #* boundary_vz(iz,ix)
		end
		end

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


"""
Need illumination to estimate the approximate diagonal of Hessian
"""
@inline function compute_illum!(p::Array{Float64}, illumisseq_out::Array{Float64})
	for ix=1:size(illumisseq_out,2)
		@simd for iz=1:size(illumisseq_out,1)
			# saving illumination to be used as preconditioner 
			illumisseq_out[iz,ix] += (p[iz,ix,1,1] * p[iz,ix,1,1])
		end
	end
end

function born_sources()
	# secondary sources for Born modeling
	if(born_flag)
	# adding born sources from pressure(:,:,1) to pressure(:,:,2)
	for ix=3:nx-2
	for iz=3:nz-2

		# upto until [it-2]
		# lambdaI scatterrer source term at [it-1]
		# p is at [it], pp is at [it-1], ppp is at [it-2]
		# dpdx is at [it-1] and dpdz is at [it-1]
		# modrrvx scatterrer source term at [it-1]
		# modrrvz scatterrer source term at [it-1]
		born_svalue_stack[iz,ix] += δt * ((-1.0 * (ppp[iz, ix, 1,1] + p[iz, ix,  1,1] - 2.0 * pp[iz, ix,  1,1]) * δmodtt[iz, ix] * δtI * δtI) )#+ (
#							  (27.e0*dpdx[iz,ix,1,1] * δmodrrvx[iz,ix] -27.e0*dpdx[iz,ix-1,1,1] * δmodrrvx[iz,ix-1] -dpdx[iz,ix+1,1,1] * δmodrrvx[iz,ix+1] +dpdx[iz,ix-2,1,1] * δmodrrvx[iz,ix-2] ) * (δx24I)) + (
#							  (27.e0*dpdz[iz,ix,1,1] * δmodrrvz[iz,ix] -27.e0*dpdz[iz-1,ix,1,1] * δmodrrvz[iz-1,ix] -dpdz[iz+1,ix,1,1] * δmodrrvz[iz+1,ix] +dpdz[iz-2,ix,1,1] * δmodrrvz[iz-2,ix] ) * (δz24I)) ) 

		p[iz,ix,1,2] += born_svalue_stack[iz,ix] * δt * δxI * δzI * modttI[iz,ix]  
	end
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


@inbounds function boundary_force(p::AbstractArray{Float64}, 
		     bl::AbstractArray{Float64}, br::AbstractArray{Float64},
		     bt::AbstractArray{Float64}, bb::AbstractArray{Float64},
		     ibx0::Int64, ibz0::Int64, ibx1::Int64, ibz1::Int64)

	nlayer=3
	for ix=1:nlayer
		@simd for iz=1:ibz1-ibz0+1
			p[ibz0+iz-1, ibx0+ix-1] = bl[iz,ix]
			p[ibz0+iz-1, ibx1-ix+1] = br[iz,ix]
		end
	end

	@simd for ix=1:ibx1-ibx0+1
		for iz=1:nlayer
			p[ibz0+iz-1,ibx0+ix-1] = bt[iz,ix]
			p[ibz1-iz+1,ibx0+ix-1] = bb[iz,ix]
		end
	end

end

@inbounds function boundary_save(p::AbstractArray{Float64}, 
		     bl::AbstractArray{Float64}, br::AbstractArray{Float64},
		     bt::AbstractArray{Float64}, bb::AbstractArray{Float64},
		     ibx0::Int64, ibz0::Int64, ibx1::Int64, ibz1::Int64)

	nlayer=3
	for ix=1:nlayer
		@simd for iz=1:ibz1-ibz0+1
			bl[iz,ix] = p[ibz0+iz-1, ibx0+ix-1] 
			br[iz,ix] = p[ibz0+iz-1, ibx1-ix+1]
		end
	end

	@simd for ix=1:ibx1-ibx0+1
		for iz=1:nlayer
			bt[iz,ix] = p[ibz0+iz-1,ibx0+ix-1]
			bb[iz,ix] = p[ibz1-iz+1,ibx0+ix-1]
		end
	end

end

function snap_save!(p::AbstractArray{Float64}, snap::AbstractArray{Float64}, isx0::Int64, isz0::Int64)
	nz, nx=size(snap)
	for ix=1:nx
		for iz=1:nz
			snap[iz,ix] = p[isz0+iz, isx0+ix]
		end
	end
end

function fill_src_wav_mat!(src_wav_mat::Array{Float64}, acqsrc_uspos::Array{Acquisition.Src},src_flags::Vector{Int64})

	src_nsmul, src_nfield, npropwav, tim_nt, src_nseq = size(src_wav_mat)
	δt = acqsrc_uspos[1].tgrid.δx
	for ipropwav=1:npropwav, ifield=1:src_nfield, isseq=1:src_nseq, issmul=1:src_nsmul
		src_tim_nt = acqsrc_uspos[ipropwav].tgrid.nx;
		if(src_flags[ipropwav] == 0)
			nothing # just put zeros, no sources added
		elseif(src_flags[ipropwav] == 1)
			"ϕ[t] = s[t]"
			for it=1:src_tim_nt
				source_term = acqsrc_uspos[ipropwav].wav[isseq,ifield][it,issmul]
				src_wav_mat[issmul,ifield,ipropwav,it,isseq] = source_term
			end
		elseif(src_flags[ipropwav] == 2)
			"ϕ[t] = ∑₁ᵗ⁻¹ s[t]"
			source_term_stack = 0.0;
			for it=1:src_tim_nt-1
				source_term_stack += (acqsrc_uspos[ipropwav].wav[isseq,ifield][it,issmul] .* δt)
				src_wav_mat[issmul,ifield,ipropwav,it+1,isseq] = source_term_stack
			end
			if(tim_nt > src_tim_nt)
				src_wav_mat[issmul,ifield,ipropwav,src_tim_nt+1:end,isseq] = src_wav_mat[issmul,ifield,ipropwav,src_tim_nt,isseq]
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
				src_wav_mat[issmul,ifield,ipropwav,tim_nt-it+1,isseq] = -1.0 * source_term_stack
			end
			if(tim_nt > src_tim_nt)
				tim_nt_diff = tim_nt-src_tim_nt
				src_wav_mat[issmul,ifield,ipropwav,1:tim_nt_diff+1,isseq] = src_wav_mat[issmul,ifield,ipropwav,tim_nt_diff+2,isseq]
			end
		end
	end
end


function check_fd_stability(vpmin::Float64, vpmax::Float64, δx::Float64, δz::Float64, freqmin::Float64, freqmax::Float64, δt::Float64, verbose::Bool)

	# check spatial sampling
	δs_temp=vpmin/5.0/freqmax;
	δs_max = maximum([δx, δz])
	all(δs_max .> δs_temp) ? 
			warn(string("spatial sampling\t",δs_max,"\ndecrease spatial sampling below:\t",δs_temp)) :
			verbose ? println("spatial sampling\t",δs_max,"\tcan be as high as:\t",δs_temp) : nothing 

	# check time sampling
	δs_min = minimum([δx, δz])
	δt_temp=0.5*δs_min/vpmax
	all(δt .> δt_temp) ? 
			warn(string("time sampling\t",δt,"\ndecrease time sampling below:\t",δt_temp)) :
			verbose ? println("time sampling\t",δt,"\tcan be as high as:\t",δt_temp) : nothing

end

end # module