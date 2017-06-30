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

#function mod_alloc(;
#	nrecwav::Int64=1, 
#	model::Models.Seismic=Gallery.Seismic(:acou_homo1),
#	tgridmod::Grid.M1D=Gallery.M1D(:acou_homo1),
#	tgrid::Grid.M1D=tgridmod,
#	acqgeom::Array{Acquisition.Geom}=fill(Gallery.Geom(:acou_homo1),npropwav),
#	recv_nfield::Int64=1, 
#	)
#
#	initp=zeros(model.mgrid.nz+2*model.mgrid.npml,model.mgrid.nx+2*model.mgrid.npml,3,acqgeom[1].nss)
#	TDout=[Data.TD_zeros(recv_nfield,tgrid,acqgeom[ip]) for ip in 1:nrecwav]
#	boundary=zeros(Grid.M2D_boundary(model.mgrid, 3, :outer, onlycount=true),tgridmod.nx,acqgeom[1].nss)
#	gmodel=Models.Seismic_zeros(model.mgrid)
#
#	return TDout, boundary, initp, gmodel
#end

"""
As forward modeling method, the 
finite-difference method is employed. 
It uses a discrete version of the two-dimensional isotropic acoustic wave equation.
As shown in  
Figure~ef{fdshmesh_acou}, a staggered-grid mesh is used 
# Keyword Arguments
* `npropwav::Int64=1` : number of independently propagating wavefields in `model`
* `model::Models.Seismic` : seismic medium parameters 
* `model_pert::Models.Seismic` : perturbed model, i.e., model + δmodel, used for Born modeling 
* `tgridmod::Grid.M1D` : modelling time vector
* `tgrid::Grid.M1D` : output records are interpolated on this time vector
* `recv_nfield::Int64=1` : number of fields that receivers record 
* `boundary_save_flag::Bool=Int64` : save final state variables and the boundary conditions for later use
* `boundary_in::Any=[0]` : input final state variables and boundary
* `abs_trbl::Vector{Symbol}=[:top, :bottom, :right, :left]` : PML boundary conditions 
* `born_flag::Bool=false` : Born modeling flag
* `gmodel_flag=false` : output gradient or not
* `verbose::Bool=false` : verbose flag

# Example
```julia
julia> records, boundary_save  = mod();
```
Credits: Pawan Bharadwaj, 2017
"""
@fastmath function mod!(;
	jobname::AbstractString = "Hello",
	npropwav::Int64=1, 
	model::Models.Seismic = Gallery.Seismic(:acou_homo1),
	abs_trbl::Vector{Symbol}=[:top, :bottom, :right, :left],
	born_flag::Bool=false,
	model_pert::Models.Seismic = model,
	tgridmod::Grid.M1D = Gallery.M1D(:acou_homo1),
	acqgeom::Array{Acquisition.Geom} = fill(Gallery.Geom(:acou_homo1),npropwav),
	acqsrc::Array{Acquisition.Src} = fill(Gallery.Src(:acou_homo1),npropwav),
	src_flags::Vector{Int64}=fill(2,npropwav), # bilinear for all
	recv_flags::Vector{Int64}=fill(1,npropwav), # bilinear for all
	recv_nfield::Int64=1, 
	TDout::Vector{Data.TD}=[Data.TD_zeros(recv_nfield,tgridmod,acqgeom[ip]) for ip in 1:length(findn(recv_flags)) ],
	backprop_flag::Int64=0,  
	boundary::Array{Float64}=
	        zeros(Grid.M2D_boundary(model.mgrid, 3, :outer, onlycount=true),tgridmod.nx,acqgeom[1].nss),
	initp::Array{Float64}=
	        zeros(model.mgrid.nz+2*model.mgrid.npml,model.mgrid.nx+2*model.mgrid.npml,3,acqgeom[1].nss),
	gmodel_flag=false,
	gmodel::Models.Seismic=Models.Seismic_zeros(model.mgrid),
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


	# check spatial sampling
	ds_temp=Models.χ([minimum(model.χvp)],model.vp0,-1)[1]/5.0/freqmax;
	ds_max = maximum([model.mgrid.δx, model.mgrid.δz])
	all(ds_max .> ds_temp) ? 
			warn(string("spatial sampling\t",ds_max,"\ndecrease spatial sampling below:\t",ds_temp)) :
			verbose ? println("spatial sampling\t",ds_max,"\tcan be as high as:\t",ds_temp) : nothing 

	# check time sampling
	ds_min = minimum([model.mgrid.δx, model.mgrid.δz])
	dt_temp=0.5*ds_min/Models.χ([(maximum(model.χvp))],model.vp0,-1)[1]
	all(tgridmod.δx .> dt_temp) ? 
			warn(string("time sampling\t",tgridmod.δx,"\ndecrease time sampling below:\t",dt_temp)) :
			verbose ? println("time sampling\t",tgridmod.δx,"\tcan be as high as:\t",dt_temp) : nothing



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

	# gradient outputs
	grad_modtt = SharedArray(Float64, (nz, nx, src_nseq))
	grad_modtt[:,:,:] = 0.0
	grad_modrr = SharedArray(Float64, (nz, nx, src_nseq))
	grad_modrr[:,:,:] = 0.0


	# pml_variables
	a_x, b_x, k_xI, a_x_half, b_x_half, k_x_halfI = pml_variables(nx, δt, δx, model.mgrid.npml-5, vpmax, vpmin, freqmin, freqmax, 
							       [any(abs_trbl .== :left), any(abs_trbl .== :right)])
	a_z, b_z, k_zI, a_z_half, b_z_half, k_z_halfI = pml_variables(nz, δt, δz, model.mgrid.npml-5, vpmax, vpmin, freqmin, freqmax,
							       [any(abs_trbl .== :top), any(abs_trbl .== :bottom)])

	# boundary coordinates
	boundary_z, boundary_x, boundary_n = Grid.M2D_boundary(model.mgrid, 3, :outer)
	iboundary_x = broadcast(x->indmin(abs(mesh_x-x)), boundary_x)
	iboundary_z = broadcast(x->indmin(abs(mesh_z-x)), boundary_z)

	"saving boundary values"
	if(backprop_flag == 1)
		boundary_out = SharedArray(Float64,size(boundary))
	end

	"p_out will save final snapshots that can be used for time reversal"
	p_out = SharedArray(Float64, size(initp))


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
	for ipropwav=1:npropwav
		for ifield=1:src_nfield
			for isseq=1:src_nseq
				for issmul=1:src_nsmul
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
		end
	end

	# source_loop
	@sync @parallel for isseq = 1:src_nseq

		if(verbose)
			println("modelling supershot:\t", isseq)
		end
		
		"allocate receiver matrix per sequential source"
		rec_mat=zeros(recv_n, recv_nfield, npropwav, tim_nt)


		# source_spray_weights
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
		allocations
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
		advance!(p, dpdx, dpdz, memory_dp_dx, memory_dp_dz, memory_dvx_dx, memory_dvz_dz,
				   modttI, modrrvx, modrrvz, 
				   δx24I, δz24I, δt, nx, nz,
				   a_x, b_x, k_xI, a_x_half, b_x_half, k_x_halfI, 
				   a_z, b_z, k_zI, a_z_half, b_z_half, k_z_halfI)

		gradis_modtt=zeros(nz,nx)
		gradis_modrrvx=zeros(gradis_modtt)
		gradis_modrrvz=zeros(gradis_modtt)

		"initial conditions from initp for first propagating field only"
		if(backprop_flag==-1)
			p[:,:,:,1].=initp[:,:,:,isseq]
		end

		# time_loop
		for it=1:tim_nt
			ppp .= pp
			pp .=  p

			advance!(p, dpdx, dpdz, memory_dp_dx, memory_dp_dz, memory_dvx_dx, memory_dvz_dz,
				   modttI, modrrvx, modrrvz, 
				   δx24I, δz24I, δt, nx, nz,
				   a_x, b_x, k_xI, a_x_half, b_x_half, k_x_halfI, 
				   a_z, b_z, k_zI, a_z_half, b_z_half, k_z_halfI)
		
			# force p[1] on boundarys
			if(backprop_flag==-1) 
				for iboundary = 1:boundary_n
					@inbounds p[iboundary_z[iboundary], iboundary_x[iboundary],1,1] = boundary[iboundary,it,isseq] 
				end
			end
	 
			add_source!(it, isseq, p, src_wav_mat, src_spray_weights,
				modttI, issmulx1, issmulx2, issmulz1, issmulz2, npropwav, src_nfield, src_nsmul,
				δt::Float64, δxI::Float64, δzI::Float64
				)


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

			# record boundarys
			if(backprop_flag==1)
				for iboundary = 1:boundary_n
				       @inbounds boundary_out[iboundary,it,isseq] = p[iboundary_z[iboundary],iboundary_x[iboundary],1,1]
				end
			end

			record!(it,rec_mat, p, rec_interpolate_weights, irecx1, irecx2, irecz1, irecz2, npropwav, recv_nfield, recv_n)

			if(gmodel_flag)
				compute_gradient!(p, pp, ppp, dpdx, dpdz,
		      			gradis_modtt, gradis_modrrvx, gradis_modrrvz,nx,nz,δtI)
			end



		end # time_loop
		"now pressure is at [tim_nt], velocities are at [tim_nt-1/2]"	

		"one more propagating step to save pressure at [tim_nt+1] -- for time revarsal"
		advance!(p, dpdx, dpdz, memory_dp_dx, memory_dp_dz, memory_dvx_dx, memory_dvz_dz,
			   modttI, modrrvx, modrrvz, 
			   δx24I, δz24I, δt, nx, nz,
			   a_x, b_x, k_xI, a_x_half, b_x_half, k_x_halfI, 
			   a_z, b_z, k_zI, a_z_half, b_z_half, k_z_halfI)
		"save last snapshot of pressure field"
		p_out[:,:,1,isseq] = p[:,:,1,1]

		"one more propagating step to save velocities at [tim_nt+3/2] -- for time reversal"
		advance!(p, dpdx, dpdz, memory_dp_dx, memory_dp_dz, memory_dvx_dx, memory_dvz_dz,
			   modttI, modrrvx, modrrvz, 
			   δx24I, δz24I, δt, nx, nz,
			   a_x, b_x, k_xI, a_x_half, b_x_half, k_x_halfI, 
			   a_z, b_z, k_zI, a_z_half, b_z_half, k_z_halfI)
		"save last snapshot of velocity fields with opposite sign for adjoint propagation"
		p_out[:,:,2:3,isseq] = -1.0 * p[:,:,2:3,1]

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
			gradis_modtt .*=  (δx * δz)
			gradis_modrrvx .*= (δx * δz)
			gradis_modrrvz .*= (δx * δz)


			grad_modtt[:,:,isseq] = gradis_modtt

			for ix=2:nx-2
				for iz=2:nz-2
					@inbounds grad_modrr[iz,ix, isseq] += 0.5e0 * gradis_modrrvx[iz,ix] 
					@inbounds grad_modrr[iz,ix, isseq] +=  0.5e0 * gradis_modrrvz[iz,ix]
				end
			end
			for ix=2:nx-2
				for iz=2:nz-2
					@inbounds grad_modrr[iz,ix+1, isseq] +=  0.5e0 * gradis_modrrvx[iz,ix]
				end
			end
			for ix=2:nx-2
				for iz=2:nz-2
					@inbounds grad_modrr[iz+1,ix, isseq] +=  0.5e0 * gradis_modrrvz[iz,ix]
				end
			end


		end

	end # source_loop

	"check if recv_out is zeros"
	isapprox(maximum(abs(recv_out)),0.0) && warn("recv_out are zeros")

	"summing gradient over all the sources after truncation in PML"
	grad_modtt_stack = Models.pad_trun(squeeze(sum(Array(grad_modtt),3),3),model.mgrid.npml,-1);
	grad_modrr_stack = Models.pad_trun(squeeze(sum(Array(grad_modrr),3),3),model.mgrid.npml,-1);


	# update gradient model using grad_modtt_stack, grad_modrr_stack
	Models.Seismic_chainrule!(gmodel, model, 
				     vec(Models.χg(grad_modtt_stack, Models.Seismic_get(model,:KI0),1)),
				     vec(Models.χg(grad_modrr_stack, Models.Seismic_get(model,:ρI0),1)),
				     [:χKI, :χρI], 1
				     )

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

	# update boundary after time reversal
	if(backprop_flag ==1)
		boundary .= Array(flipdim(boundary_out,2))
	end

	# update initp
	if(backprop_flag==1)
		initp .= Array(p_out) 
	end

	# return data depending on the receiver flags
	return TDout, boundary, initp, gmodel

	# return without resampling for testing
	#return [Data.TD(reshape(recv_out[1+(iprop-1)*nd : iprop*nd],tgridmod.nx,recv_n,src_nseq),
	#		       tgridmod, acqgeom[1]) for iprop in 1:npropwav]
end

@inbounds @fastmath function add_source!(it::Int64, isseq::Int64, 
				p::Array{Float64}, src_wav_mat::Array{Float64},
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
		source_term = src_wav_mat[issmul, ifield, ipropwav, it, isseq] * δt * δxI * δzI
		
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


@inbounds @fastmath function record!(it::Int64,
		 rec_mat::Array{Float64}, 
		 p::Array{Float64}, 
		 rec_interpolate_weights::Array{Float64},
		 irecx1::Array{Int64}, irecx2::Array{Int64},
		 irecz1::Array{Int64}, irecz2::Array{Int64},
		 npropwav::Int64, recv_nfield::Int64, 
		 recv_n::Int64)

	for ipropwav = 1:npropwav
	for ifield = 1:recv_nfield
	@simd for irec = 1:recv_n
			rec_mat[irec,ifield,ipropwav, it] = 
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

@inbounds @fastmath function compute_gradient!(p::Array{Float64}, pp::Array{Float64}, ppp::Array{Float64},
			   dpdx::Array{Float64}, dpdz::Array{Float64},
			   gradis_modtt::Array{Float64}, gradis_modrrvx::Array{Float64},
			   gradis_modrrvz::Array{Float64},
			   nx::Int64, nz::Int64, δtI::Float64 
			  )

	for ix=1:nx
		@simd for iz=1:nz
			# p at [it], pp at [it-1]
			# dpdx and dpdz at [it]
			# gradients w.r.t inverse of modtt, i.e., 1/rho/c^2 
			@inbounds gradis_modtt[iz,ix] += ((-1.0 * (ppp[iz, ix, 1,1] + p[iz, ix,  1,1] - 2.0 * pp[iz, ix,  1,1]) * δtI * δtI) *  pp[iz,ix,1,2])
			       #(p(:,:,1,1) - pp(:,:,1,1)) * δtI * &
			       #        (p(:,:,1,2) - pp(:,:,1,2)) * δtI
			
			# gradient w.r.t. inverse of rho on vx and vz grids
			@inbounds gradis_modrrvx[iz,ix] += (- dpdx[iz,ix,1,2]*dpdx[iz,ix,1,1])
			@inbounds gradis_modrrvz[iz,ix] += (- dpdz[iz,ix,1,2]*dpdz[iz,ix,1,1])
		end
	end
end


@inbounds @fastmath function advance!(
		  p::Array{Float64},  dpdx::Array{Float64}, dpdz::Array{Float64},  
		  memory_dp_dx::Array{Float64}, memory_dp_dz::Array{Float64}, memory_dvx_dx::Array{Float64}, memory_dvz_dz::Array{Float64},
		  modttI::Array{Float64}, modrrvx::Array{Float64}, modrrvz::Array{Float64},
		  δx24I::Float64, δz24I::Float64, δt::Float64,
		  nx::Int64, nz::Int64,
		  a_x::Vector{Float64}, b_x::Vector{Float64}, k_xI::Vector{Float64}, a_x_half::Vector{Float64}, b_x_half::Vector{Float64}, k_x_halfI::Vector{Float64}, 
		  a_z::Vector{Float64}, b_z::Vector{Float64}, k_zI::Vector{Float64}, a_z_half::Vector{Float64}, b_z_half::Vector{Float64}, k_z_halfI::Vector{Float64}
		 )


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
		! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
		! update velocity at [it-1/2] using 
		! velocity at [it-3/2] and dpdx and dpdz at [it-1] 
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
		! compute p[it] using p[it-1] and velocity [it-1/2]
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
			! compute pressure at [it] using p at [it-1] and dvxdx
			! and dvzdz at [it-1/2]
			"""
			@inbounds p[iz,ix,1,ipropwav] += (modttI[iz,ix] * (dpdx[iz,ix,2,ipropwav] + dpdz[iz,ix,3,ipropwav])) * δt #* boundary_p(iz,ix)
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

end # module
