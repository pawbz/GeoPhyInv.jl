module Born

using Grid
import JuMIT.Acquisition
import JuMIT.Models
import JuMIT.Data
using SpecialFunctions
using FFTW

mutable struct Param
end

"""
k = wavenumber = 2pif/v0
"""
function G0_homo_acou(x::T, z::T, k::T, ρ0::T) where {T<:AbstractFloat}
	sqrt(x*x + z*z) == 0.0 ? error("distance cannot be zero") : nothing
	k == 0.0 ? error("wavenumber cannot be zero") : nothing
	G0 = -0.25 * ρ0 * im * complex.(besselj0(k*sqrt(x*x + z*z)), -bessely0(k*sqrt(x*x + z*z)))
	return G0
end


function mod(;
	     vp0::Float64=2000.0,
	     ρ0::Float64=2000.0,
	     model_pert::Models.Seismic=nothing,
             born_flag::Bool=false,
	     tgridmod::Grid.M1D=nothing,
	     tgrid::Grid.M1D = tgridmod,
	     acqgeom::Acquisition.Geom=nothing,
	     acqsrc::Acquisition.Src=nothing,
	     src_flag::Int64=2,
		  )
	(getfield(acqgeom,:nss) != getfield(acqsrc,:nss))  ? error("different supersources") : nothing
	(getfield(acqgeom,:ns) != getfield(acqsrc,:ns))  ? error("different sources") : nothing


	nt = (length(tgridmod.x) == length(acqsrc.tgrid.x)) ? tgrid.nx : error("acqsrc tgrid")
	np2 = nextpow(2, 2*nt);	
		

	if(born_flag)
		mesh_x=model_pert.mgrid.x
		mesh_z=model_pert.mgrid.z
		nz=model_pert.mgrid.nz
		nx=model_pert.mgrid.nx
		δx = model_pert.mgrid.δx
		δz = model_pert.mgrid.δz

		δmodtt = Models.Seismic_get(model_pert, :KI) - (vp0 * vp0 * ρ0)^(-1)
		δmodrr = Models.Seismic_get(model_pert, :ρI) - (ρ0)^(-1)
	end

	# npow2 grid for time
	tnpow2grid = Grid.M1D_fft(np2, tgridmod.δx);
	# corresponding npow2 frequency grid 
	fnpow2grid = Grid.M1D_fft(tnpow2grid);

	nss = acqgeom.nss
	nr = acqgeom.nr
	ns = acqgeom.ns
	data = Data.TD(
	      [zeros(nt,nr[iss]) for iss=1:nss, ifield=1:1],
	      [:P],
	      tgridmod,acqgeom)

	for ifield = 1:length(data.fields), iss = 1:nss
		sx = acqgeom.sx[iss][:]
		rx = acqgeom.rx[iss][:]
		sz = acqgeom.sz[iss][:]
		rz = acqgeom.rz[iss][:]

		for ir = 1:nr[iss]

		dtemp=zeros(nt)
		dpow2all=complex.(zeros(np2), zeros(np2));
		wpow2=complex.(zeros(np2), zeros(np2)); 
		
		for is=1:acqgeom.ns[iss]
			# zero pad wavelet
			for it in 1:nt
				wpow2[it]=complex(acqsrc.wav[iss,ifield][it,is])
			end
			FFTW.fft!(wpow2) # source wavelet in the frequency domain

			x = sx[is] - rx[ir]
			z = sz[is] - rz[ir]
			# analytical expression for every frequency, except zero
			dpow2 = complex.(zeros(np2), zeros(np2)); 
			dpow2[1] = complex.(0.0, 0.0)
			for iω in 2:np2
				ω = 2. * pi * abs(fnpow2grid.x[iω])
				k = ω / vp0

				if(born_flag)
					term = complex(0., 0.)
					for ix=1:nx
						@simd for iz=1:nz
							if(δmodtt[iz,ix] ≠ 0.0)
								term += (G0_homo_acou(sx[is]-mesh_x[ix], sz[is]-mesh_z[iz], k, ρ0)[1] 
									 * G0_homo_acou(rx[ir]-mesh_x[ix], rz[ir]-mesh_z[iz], k, ρ0)[1] .* ω .* ω .* δmodtt[iz,ix]) * δx * δz
								# factor due to integration
							end
						end
					end

				else
					if(src_flag == 2)
						term = G0_homo_acou(x, z, k, ρ0);
					elseif(src_flag == 1)
						dpow2[iω] = G0_homo_acou(x, z, k, ρ0) * im * abs(fnpow2grid.x[iω])
					else
						error("invalid src_flag")
					end
				end
				if(fnpow2grid.x[iω] > 0) 
					dpow2[iω] = term;
				elseif(fnpow2grid.x[iω] < 0)
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
			data.d[iss,ifield][it,ir] = real(dpow2all[it])
		end
		end
	end
	return data
end


end # module 
