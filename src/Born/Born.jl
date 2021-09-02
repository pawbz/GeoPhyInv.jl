module Born

import GeoPhyInv: AGeom, SrcWav
import GeoPhyInv: Medium, complexK
import GeoPhyInv.Records
using SpecialFunctions
using Statistics
using FFTW
using DSP

mutable struct Param
end

"""
k = wavenumber = 2pif/v0
"""
function G0_homo_acou(x::T, z::T, k, rho0::T) where {T<:AbstractFloat}
	sqrt(x*x + z*z) == 0.0 ? error("distance cannot be zero") : nothing
	k == 0.0 ? error("wavenumber cannot be zero") : nothing
	G0 = -0.25 * rho0 * im * hankelh2(0,k*sqrt(x*x + z*z))
	return G0
end


function mod(medium::Medium;
		 Q=nothing,
	     medium_pert::Medium=nothing,
             born_flag::Bool=false,
	     tgridmod::StepRangeLen=nothing,
	     tgrid::StepRangeLen = tgridmod,
	     ageom::AGeom=nothing,
	     srcwav::SrcWav=nothing,
	     src_flag::Int64=2,
		  )
	(length(ageom) != length(srcwav))  ? error("different supersources") : nothing
	@info "fix this"
#	(length(ageom) != getfield(srcwav,:ns))  ? error("different sources") : nothing


	nt = (length(tgridmod) == length(srcwav[1].grid)) ? length(tgrid) : error("srcwav tgrid")
	np2 = nextpow(2, 2*nt);	
	fnpow2grid = FFTW.fftfreq(np2,inv(step(tgridmod)));


	if(:Q ∈ names(medium.m)[1])
		vp0=medium.ref[:vp]
		rho0=medium.ref[:rho]
		Ku=rho0*vp0*vp0
		tau_epsilon=medium[:tau_epsilon]
		tau_sigma=medium[:tau_sigma]

		Kc=complexK(Ku, abs.(2. * pi .* fnpow2grid), tau_sigma, tau_epsilon)

		rho0=fill(rho0, length(fnpow2grid))
		vp0=sqrt.(Kc .* inv.(rho0))
	else
		vp0=fill(medium.ref[:vp], length(fnpow2grid))
		rho0=fill(medium.ref[:rho], length(fnpow2grid))
	end

		

	if(born_flag)
		mesh_x=medium_pert.mgrid[2]
		mesh_z=medium_pert.mgrid[1]
		nz=length(mesh_z)
		nx=length(mesh_x)
		δx = step(mesh_x)
		δz = step(mesh_z)

		δmodtt = medium_pert[:KI] - (vp0 * vp0 * rho0)^(-1)
		δmodrr = medium_pert[:rhoI] - (rho0)^(-1)
	end

	

	nss = length(ageom)
	data=Records(tgridmod, ageom, [:p])

	for iss in 1:nss, ifield in 1:length(data[iss].d)
		sx = ageom[iss].s[:x]
		rx = ageom[iss].r[:x]
		sz = ageom[iss].s[:z]
		rz = ageom[iss].r[:z]

		for ir = 1:ageom[iss].nr

		dtemp=zeros(nt)
		dpow2all=complex.(zeros(np2), zeros(np2));
		wpow2=complex.(zeros(np2), zeros(np2)); 
		
		for is=1:ageom[iss].ns
			# zero pad wavelet
			for it in 1:nt
				wpow2[it]=complex(srcwav[iss].d[ifield][it,is])
			end
			FFTW.fft!(wpow2) # source wavelet in the frequency domain

			x = sx[is] - rx[ir]
			z = sz[is] - rz[ir]
			# analytical expression for every frequency, except zero
			dpow2 = complex.(zeros(np2), zeros(np2)); 
			dpow2[1] = complex.(0.0, 0.0)
			for iω in 2:np2
				ω = 2. * pi * abs(fnpow2grid[iω])
				k = ω * inv(vp0[iω])

				if(born_flag)
					term = complex(0., 0.)
					for ix=1:nx
						@simd for iz=1:nz
							if(δmodtt[iz,ix] ≠ 0.0)
								term += (G0_homo_acou(sx[is]-mesh_x[ix], sz[is]-mesh_z[iz], k, rho0)[1] 
									 * G0_homo_acou(rx[ir]-mesh_x[ix], rz[ir]-mesh_z[iz], k, rho0)[1] .* ω .* ω .* δmodtt[iz,ix]) * δx * δz
								# factor due to integration
							end
						end
					end

				else
					if(src_flag == 2)
						term = G0_homo_acou(x, z, k, rho0[iω]);
					elseif(src_flag == 1)
						dpow2[iω] = G0_homo_acou(x, z, k, rho0[iω]) * im * abs(fnpow2grid[iω])
					else
						error("invalid src_flag")
					end
				end
				if(fnpow2grid[iω] > 0) 
					dpow2[iω] = term;
				elseif(fnpow2grid[iω] < 0)
					dpow2[iω] = conj(term);
				end

			end
			# convolution and stack over simultaneous sources
			dpow2all += dpow2 .* wpow2;
		end

		# back to time domain
		FFTW.ifft!(dpow2all)

		# truncate
		for it in 1:nt
			data[iss].d[ifield][it,ir] = real(dpow2all[it])
		end
		end
	end
	return data
end


end # module 
