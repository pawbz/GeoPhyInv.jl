module Analytic

import SIT.Grid
import SIT.Acquisition
import SIT.Data
import SIT.DSP
import SIT.Gallery


"""
k = wavenumber = 2pif/v0
"""
function G0_homo_acou{T<:AbstractFloat}(x::T, z::T, k::T)
	sqrt(x*x + z*z) == 0.0 ? error("distance cannot be zero") : nothing
	k == 0.0 ? error("wavenumber cannot be zero") : nothing
	G0 = complex(T(0.0),T(0.25)) * hankelh1(0,k*sqrt(x*x + z*z))
	G0z = (complex(T(0.0),T(-0.25))*complex(k*z,T(0.0))*hankelh1(1,(k*sqrt(x*x + z*z)))
			/complex(sqrt(x*x + z*z),T(0.0)))
	G0x = (complex(T(0.0),T(-0.25))*complex(k*x,T(0.0))*hankelh1(1, (k*sqrt(x*x + z*z)))
			/complex(sqrt(x*x + z*z),T(0.0)))
	return G0
end


function mod(;
	     vp0::Float64=2000.0,
	     ρ0::Float64=1000.0,
	     tgridmod::Grid.M1D = Gallery.M1D(:acou_homo1),
	     tgrid::Grid.M1D = tgridmod,
	     acqgeom::Acquisition.Geom=Gallery.Geom(:acou_homo1),
	     acqsrc::Acquisition.Src=Gallery.Src(:acou_homo1),
	     src_flags::Float64=2.0,
		  )
	(getfield(acqgeom,:nss) != getfield(acqsrc,:nss))  ? error("different supersources") : nothing
	(getfield(acqgeom,:ns) != getfield(acqsrc,:ns))  ? error("different sources") : nothing


	nt = (length(tgridmod.x) == length(acqsrc.tgrid.x)) ? tgrid.nx : error("acqsrc tgrid")
	np2 = nextpow2(2*nt);	
		
	wpow2=complex(zeros(np2), zeros(np2)); 
	dpow2all=complex(zeros(np2), zeros(np2));
	dtemp=zeros(nt)

	# npow2 grid for time
	tnpow2grid = Grid.M1D_npow2(np2, tgridmod.δx);
	# corresponding npow2 frequency grid 
	fnpow2grid = Grid.M1D_npow2_tf(tnpow2grid);

	nss = acqgeom.nss
	nr = acqgeom.nr
	data = Data.TD(
	      [zeros(nt,nr[iss]) for iss=1:nss, ifield=1:1],
	      1,tgridmod,acqgeom)

	for ifield = 1:data.nfield, iss = 1:nss, ir = 1:nr[iss]
		
		dpow2 = complex(zeros(np2), zeros(np2)); 
		for is=1:acqgeom.ns[iss]
			# zero pad wavelet
			DSP.nlag_npow2_pad_truncate!(acqsrc.wav[iss,ifield][:,is], wpow2, nt-1, 0, np2, 1)
			fft!(wpow2) # source wavelet in the frequency domain

			x = acqgeom.sx[iss][is] - acqgeom.rx[iss][ir]
			z = acqgeom.sz[iss][is] - acqgeom.rz[iss][ir]
			for iω in 1:np2
				k = 2. * pi * fnpow2grid.x[iω] / vp0
				G0 = G0_homo_acou(x, z, k)
				dpow2[iω] = G0[1]
			end
			# convolution and stack over simultaneous sources
			dpow2all += dpow2 .* wpow2;
		end

		# back to time domain
		ifft!(dpow2all)

		# truncate
		DSP.nlag_npow2_pad_truncate!(dtemp, dpow2all, nt-1, 0, np2, -1)
		data.d[iss,ifield][:,ir] = dtemp
	end
	return Data.TD_resamp(data, tgrid) 
end

	#

end # module 
