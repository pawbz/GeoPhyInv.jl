"""
Compute coeffs that are required to solve the memory variable eq.37 from Robertsson et. al.
"""
function get_memcoeff(mod::Medium)
	tau_sigma = mod[:tau_sigma]
	tau_epsilon = mod[:tau_epsilon]

	@assert size(tau_epsilon) == size(tau_sigma)

	nsls, nz, nx = size(tau_sigma)
	# compute factor to convert from relaxed (ω -> 0) to unrelaxed (ω -> ∞)

	K_relaxation_factor=relaxed_conversion_factor(tau_sigma, tau_epsilon)

	memcoeff2=inv.(tau_sigma) .* (1.0 .- tau_epsilon .* inv.(tau_sigma))
	for ix in 1:nx
		for iz in 1:nz
			for i in 1:nsls
				memcoeff2[i,iz,ix] *= inv(K_relaxation_factor[iz,ix])
			end
		end
	end

	memcoeff1=inv.(tau_sigma)

	return memcoeff1, memcoeff2
end
