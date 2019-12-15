module Born

import GeoPhyInv: AGeom, SrcWav
import GeoPhyInv: Medium
import GeoPhyInv.Data
using SpecialFunctions
using FFTW
using DSP

mutable struct Param
end

"""
k = wavenumber = 2pif/v0
"""
function G0_homo_acou(x::T, z::T, k::T, rho0::T) where {T<:AbstractFloat}
	sqrt(x*x + z*z) == 0.0 ? error("distance cannot be zero") : nothing
	k == 0.0 ? error("wavenumber cannot be zero") : nothing
	G0 = -0.25 * rho0 * im * complex.(besselj0(k*sqrt(x*x + z*z)), -bessely0(k*sqrt(x*x + z*z)))
	return G0
end


function mod(;
	     vp0::Float64=2000.0,
	     rho0::Float64=2000.0,
	     model_pert::Medium=nothing,
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
		

	if(born_flag)
		mesh_x=model_pert.mgrid[2]
		mesh_z=model_pert.mgrid[1]
		nz=length(mesh_z)
		nx=length(mesh_x)
		δx = step(mesh_x)
		δz = step(mesh_z)

		δmodtt = model_pert[:KI] - (vp0 * vp0 * rho0)^(-1)
		δmodrr = model_pert[:rhoI] - (rho0)^(-1)
	end

	
	fnpow2grid = FFTW.fftfreq(np2,inv(step(tgridmod)));

	nss = length(ageom)
	data=Data(tgridmod, ageom, [:P])

	for iss in 1:nss, ifield in 1:length(data[iss].d)
		sx = ageom[iss].sx
		rx = ageom[iss].rx
		sz = ageom[iss].sz
		rz = ageom[iss].rz

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
				k = ω / vp0

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
						term = G0_homo_acou(x, z, k, rho0);
					elseif(src_flag == 1)
						dpow2[iω] = G0_homo_acou(x, z, k, rho0) * im * abs(fnpow2grid[iω])
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
