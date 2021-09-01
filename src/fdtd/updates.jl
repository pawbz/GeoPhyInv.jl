# Methods to update FDTD structs in several ways


"""
update the perturbation vector using the perturbed medium
in this case, medium will be treated as the background medium 
* `δmod` is [δKI, δrhoI]
"""
function update_δmods!(pac::P_common, δmod::Vector{Float64})
	nx=pac.ic[:nx]; nz=pac.ic[:nz]
	nznxd=prod(length.(pac.medium.mgrid))
	copyto!(pac.δmodall, δmod)
	fill!(pac.δmod[:tt],0.0)
	δmodtt=view(pac.δmod[:tt],npml+1:nz-npml,npml+1:nx-npml)
	for i in 1:nznxd
		# put perturbation due to KI as it is
		δmodtt[i] = δmod[i] 
	end
	fill!(pac.δmod[:rr],0.0)
	δmodrr=view(pac.δmod[:rr],npml+1:nz-npml,npml+1:nx-npml)
	for i in 1:nznxd
		# put perturbation due to rhoI here
		δmodrr[i] = δmod[nznxd+i]
	end
	# project δmodrr onto the vz and vx grids
	get_rhovxI!(pac.δmod[:rrvx], pac.δmod[:rr])
	get_rhovzI!(pac.δmod[:rrvz], pac.δmod[:rr])
	return nothing
end

"""
This method should be executed only after the updating the main medium.
Update the `δmods` when a perturbed `medium_pert` is input.
The medium through which the waves are propagating 
is assumed to be the background medium.
"""
function update_δmods!(pac::P_common, medium_pert::Medium)
	nznxd=prod(length.(pac.medium.mgrid))
	fill!(pac.δmodall,0.0)
	copyto!(pac.δmodall, medium_pert, [:KI, :rhoI])

	for i in 1:nznxd
		pac.δmodall[i] -= pac.mod[:tt][i] # subtracting the background medium
		pac.δmodall[nznxd+i] -= pac.mod[:rr][i] # subtracting the background medium
	end
	update_δmods!(pac, pac.δmodall)
	return nothing
end

"""
```julia
update!(pa,medium_new)
```
Update `pa` with a new bundle of medium parameters `medium_new`, without additional memory allocation.
This routine is used during inversion, where medium parameters are iteratively updated. 
The ability to iteratively run the forward mediuming task (with no additional memory allocation) on  
various subsurface mediums is necessary while implementing inversion 
algorithms.
"""
function update!(pa::PFdtd, medium::Medium)
	return update!(pa.c, medium)
end
function update!(pac::P_common, medium::Medium)

	copyto!(pac.medium, medium)

	Medium_pml_pad_trun!(pac.exmedium, pac.medium, nlayer_rand)

	copyto!(pac.mod[:ttI], pac.exmedium, [:K]) 
	copyto!(pac.mod[:tt], pac.exmedium, [:KI]) 
	copyto!(pac.mod[:rr], pac.exmedium, [:rhoI])
	get_rhovxI!(pac.mod[:rrvx], pac.mod[:rr])
	get_rhovzI!(pac.mod[:rrvz], pac.mod[:rr])
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
	(length(srcwav) ≠ pa.c.ic[:npw]) && error("cannot update")
	for i in 1:length(srcwav)
		copyto!(pa.c.srcwav[i], srcwav[i])
	end
	if(!(sflags===nothing))
		copyto!(pa.c.sflags, sflags)
	end

	for ipw in 1:pa.c.ic[:npw]
		# fill_wavelets for each supersource
		@sync begin
			for (ip, p) in enumerate(procs(pa.p))
				@async remotecall_wait(p) do 
					pap=localpart(pa.p)[ipw]
					for is in 1:length(pap.ss)
						iss=pap.ss[is].iss
						wavelets=pap.ss[is].wavelets
						broadcast(x->fill!.(x,0.0), wavelets)
						fill_wavelets!(ipw, iss, wavelets, pa.c.srcwav, pa.c.sflags)
					end
				end
			end
		end
	end
end

function update!(pass::P_x_worker_x_pw_x_ss, ipw, iss, ageomss::AGeomss, pac)
	@assert ageomss.ns == pac.ageom[ipw][iss].ns
	@assert ageomss.nr == pac.ageom[ipw][iss].nr

	mesh_x, mesh_z = pac.exmedium.mgrid[2], pac.exmedium.mgrid[1]
	ssprayw=pass.ssprayw
	rinterpolatew=pass.rinterpolatew
	fill!(ssprayw,0.0)
	fill!(rinterpolatew,0.0)
	sindices=pass.sindices
	rindices=pass.rindices

	for is=1:ageomss.ns
		weights=ssprayw
		Interpolation.get_spray_weights!(view(weights, :,is),  
			    view(sindices[:x1],is), view(sindices[:x2],is),
			    view(sindices[:z1],is), view(sindices[:z2],is),
			    mesh_x, mesh_z, ageomss.s[:x][is], ageomss.s[:z][is])
	end
	for ir=1:ageomss.nr
		weights=rinterpolatew
		Interpolation.get_interpolate_weights!(view(weights, :,ir),
			  view(rindices[:x1],ir), view(rindices[:x2],ir),
			  view(rindices[:z1],ir), view(rindices[:z2],ir),
			  mesh_x, mesh_z, ageomss.r[:x][ir], ageomss.r[:z][ir])
	end

end

# if just one propagating field
update!(pa::PFdtd, ageom::AGeom)=update!(pa,[ageom])

function update!(pa::PFdtd, ageom::Vector{AGeom})
	for ipw in 1:pa.c.ic[:npw]
		copyto!(pa.c.ageom[ipw], ageom[ipw])
		@sync begin
			for (ip, p) in enumerate(procs(pa.p))
				@async remotecall_wait(p) do 
					pap=localpart(pa.p)[ipw]
					for is in 1:length(pap.ss)
						iss=pap.ss[is].iss
						update!(localpart(pap).ss[is],ipw,iss,ageom[ipw][iss],pa.c)
					end
				end
			end
		end
	end
end


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


	@timeit to "mod_x_proc!" begin
		# all localparts of DArray are input to this method
		# parallelization over shots
		@sync begin
			for (ip, p) in enumerate(procs(pa.p))
				@async remotecall_wait(p) do 
					mod_x_proc!(pa.c, localpart(pa.p))
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
		# update gradient medium using grad_modtt_stack, grad_modrr_stack
		(pa.c.gmodel_flag) && update_gradient!(pa.c)
	end


	@timeit to "record data" begin
	for ipw in pa.c.activepw
		if(pa.c.rflags[ipw] ≠ 0) # record only if rflags is non-zero
			for rfield in pa.c.rfields

				fill!(pa.c.datamat, 0.0)
				@sync begin
					for (ip, p) in enumerate(procs(pa.p))
						@sync remotecall_wait(p) do 
							update_datamat!(rfield, ipw, pa.c, localpart(pa.p))
						end
					end
				end
				update_data!(rfield, ipw, pa.c)
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

function update_datamat!(rfield, ipw, pac::P_common, pap::P_x_worker)
	datamat=pac.datamat
	pass=pap[ipw].ss
	for issp in 1:length(pass)
		iss=pass[issp].iss
		records=pass[issp].records
		for ir in 1:pac.ageom[ipw][iss].nr
			for it in 1:pac.ic[:nt]
				datamat[it,ir,iss]=records[rfield][it,ir]
			end
		end
        end
end

function update_data!(rfield, ipw, pac::P_common)
	datamat=pac.datamat
	for iss in 1:length(pac.ageom[1])
		data=pac.data[ipw][iss].d[rfield]
		for ir in 1:pac.ageom[ipw][iss].nr
			for it in 1:pac.ic[:nt]
				data[it,ir]=datamat[it,ir,iss]
			end
		end
	end
end

