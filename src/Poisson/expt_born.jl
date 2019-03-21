
# born forward modeling to generate δψ
# update psi using LP, has a loop over snapshots
# most expensive routine
function bornmod!(pa::ParamExpt, δσ, σ0)
	updateA!(pa.paσ, δσ)
	updateA!(pa.paσ0, σ0)
	qrA=factorize(pa.paσ.A)
	qrA0=factorize(pa.paσ0.A)
	fill!(pa.g, 0.0)
	@showprogress "time loop ψ\t" for it in 1:length(pa.tgrid)
		snap_in=view(pa.LP,:,:,it)
		# generate ψ0
		forwψ!(pa.ψ0,snap_in,pa.paσ0,qrA)

		# compute A(δσ)*ψ0, store this in adjψ
		applyA!(pa.adjψ, pa.ψ0, pa.paσ; A=qrA)

		# compute A(σ0)⁻¹ * adjψ
		applyinvA!(pa.ψ, pa.adjψ, pa.paσ0; A=qrA0)

		# record
		dat_slice=view(pa.data,it,:)
		mul!(dat_slice,pa.ACQ,pa.ψ) # record

	end
	return pa
end



function operator_Born(pa)
	fw=(y,x)->Fborn_map!(y, x, pa)
	bk=(y,x)->Fadj_map!(y, x, pa)

	return LinearMap(fw, bk, 
		  length(pa.paTD.dJx),  # length of output
		  xfwi_ninv(pa), # length of input
		  ismutating=true)
end

function Fborn_map!(δd, δσ, pa, σ0)
	bornmod!(pa, δσ, σ0)
	copyto!(δd, pa.data)
end

function Fadj_map!(δy, δx, pa)
	updateA!(pa.paσ, δσ)
	updateA!(pa.paσ0, σ0)
	qrA=factorize(pa.paσ.A)
	qrA0=factorize(pa.paσ0.A)
	fill!(pa.g, 0.0)
	@showprogress "time loop ψ\t" for it in 1:length(pa.tgrid)
		snap_in=view(pa.LP,:,:,it)
		# generate ψ0
		forwψ!(pa.ψ0,snap_in,pa.paσ0,qrA)

		# compute A(δσ)*ψ0, store this in adjψ
		applyA!(pa.adjψ, pa.ψ0, pa.paσ; A=qrA)

		# compute A(σ0)⁻¹ * adjψ
		applyinvA!(pa.ψ, pa.adjψ, pa.paσ0; A=qrA0)

		# record
		dat_slice=view(pa.data,it,:)
		mul!(dat_slice,pa.ACQ,pa.ψ) # record


		adjψ!(pa.adjψ, pa.adjψ2, pa.ψ0, 
		      pa.adjsrc, pa.g, pa.gtemp,
		      pa.ACQ, pa.data_misfit, pa.paσ, qrAt)

	end
	return pa


	adjψ_core!(adjψ, adjψ2, ψ, adjsrc, g, gtemp, ACQ, data_misfit, pa, qrAt)
updateA!(pa.paσ, σ)
	qrA=factorize(pa.paσ.A)
	qrAt=factorize(pa.paσ.At)
	fill!(pa.g, 0.0)
	@showprogress "time loop ψ\t" for it in 1:length(pa.tgrid)
		snap_in=view(pa.LP,:,:,it)
		forwψ!(pa.ψ,snap_in,pa.paσ,qrA)

		# record
		dat_slice=view(pa.data,it,:)
		mul!(dat_slice,pa.ACQ,pa.ψ) # record

		# wanna do adjoint modeling? depends on mode
		dat_slice_obs=view(pa.data_obs,it,:)
		adjψ!(pa.adjψ, pa.adjψ2, pa.ψ, 
		      pa.adjsrc, pa.g, pa.gtemp,
		      pa.ACQ, pa.data_misfit, dat_slice_obs, pa.paσ, qrAt, mode)
	end
	return pa

end


