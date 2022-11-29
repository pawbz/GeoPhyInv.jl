"""
complex bulk modulus K(ω)
then Q is given by Re(K(ω))/Im(K(ω))
* `nsls` : number of standard linear solids
* 
"""
function get_relaxation_times(Q, nsls, fmin, fmax)
	tau_epsilon=zeros(nsls, size(Q)...)
	tau_sigma=zeros(nsls, size(Q)...)

	nom=100
	om=range(fmin, stop=fmax, length=nom);

	# unifrmly distribute tau_sigma over the pass band (eq. 27)
	τ_σ0=inv.(exp.(range(log(fmin), stop=log(fmax), length=nsls)))

	for ix in 1:size(Q,2)
		for iz in 1:size(Q,1)
			Q1=Q[iz,ix]

			for i in 1:nsls
				tau_sigma[i, iz, ix] = τ_σ0[i]
			end


			# see eq. 24
			A=hcat([om .* inv.(1. .+ om.*om .* tau_sigma[i,iz,ix].*tau_sigma[i,iz,ix]) for i in 1:nsls]...)

			# our objective is to optimize Q such that it is close to Q_0, for all ω	
			τ_ϵ0 = A\fill(inv(Q1),nom)

			for i in 1:nsls
				tau_epsilon[i,iz,ix] = tau_sigma[i,iz,ix] .+ τ_ϵ0[i]
			end
		end
    end
    return tau_sigma, tau_epsilon
end


"""
Compute that factor that converts unrelaxed bulk modulus to relaxed bulk modulus.
"""
function relaxed_conversion_factor(tau_sigma, tau_epsilon)
	@assert size(tau_epsilon) == size(tau_sigma)
	nsls, nz, nx = size(tau_sigma)
	K_relaxation_factor=zeros(nz,nx)
	for ix in 1:nx
		for iz in 1:nz
			stau=0.0
			for i in 1:nsls
				stau  += (tau_epsilon[i,iz,ix]*inv(tau_sigma[i,iz,ix]) - 1)
			end
			K_relaxation_factor[iz,ix] = 1.0 - stau 
		end
	end
	return K_relaxation_factor
end


# return complex bulk modulus, for a frequency domain solver (currently implemented for homogeneous media)
function complexK(Ku, freqs, tau_sigma, tau_epsilon)
    fac=relaxed_conversion_factor(tau_sigma, tau_epsilon)[1]
    nsls=size(tau_sigma,1)

    Kc=complex.(zeros(length(freqs)), zeros(length(freqs)))
    for (ifreq,freq) in enumerate(freqs)
        s=complex(0.0,0.0)
        for isls in 1:nsls
            s+=(1.0 + im*freq*tau_epsilon[isls,1,1])*inv(1.0 + im*freq*tau_sigma[isls,1,1])
        end
		Kc[ifreq]=(1.0 + inv(nsls)*(-nsls + s))*Ku#*inv(fac)
		s=0.0
		for isls in 1:nsls
			s+=tau_epsilon[isls,1,1]*inv(tau_sigma[isls,1,1])
		end
		Kc[ifreq] *= inv(1 + inv(nsls)*(-nsls + s))
	end
	return Kc
	
end