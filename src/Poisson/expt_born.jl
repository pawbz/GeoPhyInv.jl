using LinearMaps
"""
A(σ₀)δψ=-A(δσ)ψ₀
Therefore,
perturbed data: δψ=-A⁻¹(σ₀)A(δσ)ψ₀
"""
function bornmod!(pa::ParamExpt, δσ, σ0)
	updateA!(pa.paσ, δσ)
	updateA!(pa.paσ0, σ0)
	# pa.paσ.A[end,:] .= 0.0 # is it necessary? A=A(σ₀)+A(δσ), except for the last row
	qrA0=factorize(pa.paσ0.A)
	@showprogress 5 "time loop ψ\t" for it in 1:length(pa.tgrid)
		snap_in=view(pa.LP,:,:,it)
		# generate ψ0
		forwψ!(pa.ψ0,snap_in,pa.paσ0,qrA0)

		# compute A(δσ)*ψ0, store this in adjψ
		applyA!(pa.adjψ, pa.ψ0, pa.paσ; A=pa.paσ.A)

		# compute A(σ0)⁻¹ * adjψ
		# do we need to mute boundary source?! --- otherwise it seems to bug!
		applyinvA!(pa.ψ, pa.adjψ, pa.paσ0; A=qrA0, mute_boundary_source=false)

		# record
		dat_slice=view(pa.data,it,:)
		mul!(dat_slice,pa.ACQ,pa.ψ) # record
		rmul!(dat_slice, -1.0)
	end
	return nothing
end


# just want to extract the data out into δd
function Fborn_map!(δd, δσ, pa, σ0)
	bornmod!(pa, δσ, σ0)
	copyto!(δd, pa.data)
end

function Fadj_map!(δσ, δd, pa, σ0)
	updateA!(pa.paσ0, σ0)
	qrAt0=factorize(pa.paσ0.At)
	qrA0=factorize(pa.paσ0.A)
	fill!(pa.g, 0.0)
	@showprogress 10 "time loop ψ\t" for it in 1:length(pa.tgrid)
		snap_in=view(pa.LP,:,:,it)
		# generate ψ0
		forwψ!(pa.ψ0,snap_in,pa.paσ0,qrA0)

		# prepare the adjoint source
		for ir in 1:size(pa.ACQ,1)
			pa.data_misfit[ir]=δd[it+(ir-1)*length(pa.tgrid)]
		end
		# migrate the data!; using this routine will also test its consistency
		adjψ_core!(pa.adjψ, pa.adjψ2, pa.ψ0, 
		      pa.adjsrc, pa.g, pa.gtemp,
		      pa.ACQ, pa.data_misfit, pa.paσ0, qrAt0)

	end
	copyto!(δσ,pa.g)
	return nothing
end




"""
Return a `LinearMap` object to perform linearized modeling and its transpose.

# Arguments
* `pa::PoissonExpt` 
* `σ0::Array` 
"""
function operator_Born(pa::ParamExpt, σ0)
	@assert length(σ0)==prod(length.(pa.mgrid))
	fw=(y,x)->Fborn_map!(y, x, pa, σ0)
	bk=(y,x)->Fadj_map!(y, x, pa, σ0)

	return LinearMap(fw, bk, 
		  length(pa.data),  # length of output
		  length(pa.σ), # length of input
		  ismutating=true)
end

