
"""
complex bulk modulus K(ω)
then Q is given by Re(K(ω))/Im(K(ω))
* `n_sls` : number of standard linear solids
* 
"""
function calculate_relaxation_times(Q, n_sls, fmin, fmax)
	τ_ϵ=zeros(n_sls, size(Q)...)
	τ_σ=zeros(n_sls, size(Q)...)

	nom=100
	om=range(fmin, stop=fmax, length=nom);

	# unifrmly distribute τ_σ over the pass band (eq. 27)
	τ_σ0=inv.(exp.(range(log(fmin), stop=log(fmax), length=n_sls)))

	for ix in 1:size(Q,2)
		for iz in 1:size(Q,1)
			Q1=Q[iz,ix]

			for i in 1:n_sls
				τ_σ[i, iz, ix] = τ_σ0[i]
			end


			# see eq. 24
			A=hcat([om .* inv.(1. .+ om.*om .* τ_σ[i,iz,ix].*τ_σ[i,iz,ix]) for i in 1:n_sls]...)

			# our objective is to optimize Q such that it is close to Q_0, for all ω	
			τ_ϵ0 = A\fill(inv(Q1),nom)


			for i in 1:n_sls
				τ_ϵ[i,iz,ix] = τ_σ[i,iz,ix] .+ τ_ϵ0[i]
			end
		end
	end

	# compute factor to convert from relaxed (ω -> 0) to unrelaxed (ω -> ∞)

	K_relaxation_factor=zeros(size(Q))
	for ix in 1:size(Q,2)
		for iz in 1:size(Q,1)
			stau=0.0
			for i in 1:n_sls
				stau  += (1. - τ_ϵ[i,iz,ix]*inv(τ_σ[i,iz,ix]))
			end
			K_relaxation_factor[iz,ix] = 1.0 - stau 
		end
	end

	Rmemory_coeff2=inv.(τ_σ) .* (1.0 .- τ_ϵ .* inv.(τ_σ))
	for ix in 1:size(Q,2)
		for iz in 1:size(Q,1)
			for i in 1:n_sls
				Rmemory_coeff2[i,iz,ix] *= inv(K_relaxation_factor[iz,ix])
			end
		end
	end

	Rmemory_coeff1=inv.(τ_σ)


	return Rmemory_coeff1, Rmemory_coeff2
end








#=
  call compute_attenuation_coeffs(N_SLS,QKappa,f0,f_min_attenuation,f_max_attenuation,tau_epsilon_kappa,tau_sigma_kappa)


subroutine compute_attenuation_coeffs(N,Qref,f0,f_min,f_max,tau_epsilon,tau_sigma)

  implicit none

! pi
  double precision, parameter :: PI = 3.141592653589793d0
  double precision, parameter :: TWO_PI = 2.d0 * PI

  integer, intent(in) :: N
  double precision, intent(in) :: Qref,f_min,f_max,f0
  double precision, dimension(1:N), intent(out) :: tau_epsilon,tau_sigma

  integer i
  double precision, dimension(1:N) :: point,weight

! nonlinear optimization with constraints
  call nonlinear_optimization(N,Qref,f0,point,weight,f_min,f_max)

  do i = 1,N
    tau_sigma(i) = 1.d0 / point(i)
    tau_epsilon(i) = tau_sigma(i) * (1.d0 + N * weight(i))
  enddo

! print *,'points = '
! do i = 1,N
!   print *,point(i)
! enddo
! print *

! print *,'weights = '
! do i = 1,N
!   print *,weight(i)
! enddo
! print *

  print *,'tau_epsilon computed by SolvOpt() = '
  do i = 1,N
    print *,tau_epsilon(i)
  enddo
  print *

  print *,'tau_sigma computed by SolvOpt() = '
  do i = 1,N
    print *,tau_sigma(i)
  enddo
  print *

end subroutine compute_attenuation_coeffs

=#
