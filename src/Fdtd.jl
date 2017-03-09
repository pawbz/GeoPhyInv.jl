module Fdtd

import SIT.Grid
import SIT.Models
import SIT.F90libs
import SIT.Acquisition


"""
test function
"""
function __init__()

end


"""
modeling
"""
function fdtd_mod(;
		  jobname::AbstractString = "Hello",
		  npropwav::Int64 = 1, 
		  mgrid::Grid.M2D = Grid.M2D(:samp1),
		  model::Models.Seismic = Models.Seismic("test_homo_acoustic"),
		  tgrid::Grid.M1D = Grid.M1D(:timesamp1),
		  acqgeom::Acquisition.Geom = Acquisition.Geom(),
		  acqsrc::Acquisition.Src = Acquisition.Src(),
		  src_nsmul::Int64 = 1,
		  src_flags::AbstractString = "[BILINEAR]",
		  recv_flags = "[BILINEAR]"
		 )

	#if(tgrid.nt .lt. src_nt) call abort_msg("fd2_mod: tgrid.nt .lt. src_nt")

	#! no modeling if source wavelet is zero
	#if(maxval(abs(src_wavelets)) .lt. tiny(rzero_de)) then
	#        return
	#endif


	recv_out = zeros(tgrid.nt*acqgeom.nr*npropwav*acqgeom.ns)
	snaps_in = zeros(tgrid.nt)
	snaps_out = ones(tgrid.nt)
	grad_modtt = zeros(mgrid.nz, mgrid.nx)
	grad_modrr = zeros(mgrid.nz, mgrid.nx)


	ccall( (:fdtd_mod, F90libs.fdtd), Void,
	      (Ptr{UInt8}, Ref{Int64},
	       Ptr{Float64}, Ptr{Float64},
	       Ref{Int64}, Ref{Int64},
	       Ref{Float64}, Ref{Float64},
	       Ref{Float64}, Ref{Float64},
	       Ref{Int64}, Ptr{UInt8},
	       Ref{Int64}, Ref{Float64},
	       Ref{Int64}, Ptr{Float64},
	       Ref{Int64}, Ref{Int64},
	       Ptr{Float64}, Ptr{Float64},
	       Ptr{UInt8},
	       Ref{Int64}, Ptr{Float64}, Ptr{Float64},
	       Ptr{Float64},
	       Ptr{UInt8},
	       Ptr{Float64},
	       Ptr{Float64},
	       Ptr{Float64},
	       Ptr{Float64}
	       ),
	      jobname, npropwav,
	      Models.χ(model.χKI, model.K0I,-1),
	      Models.χ(model.χρI, model.ρ0I,-1),
	      mgrid.nx, mgrid.nz,
	      mgrid.δx, mgrid.δz,
	      mgrid.x, mgrid.z,
	      mgrid.npml, "mgrid_abs_trbl",
	      tgrid.nt, tgrid.δt,
	      acqsrc.tgrid.nt, acqsrc.wav,
	      acqgeom.ns, src_nsmul,
	      acqgeom.sx, acqgeom.sz,
	      src_flags,
	      acqgeom.nr,
	      acqgeom.rx, acqgeom.rz,
	      recv_out,
	      recv_flags,
	      snaps_in,
	      snaps_out,
	      grad_modtt,
	      grad_modrr,
	     )


	isapprox(maximum(abs(recv_out)),0.0) && warn("receiver records are zeros")

	return reshape(recv_out,(tgrid.nt,acqgeom.nr,npropwav,acqgeom.ns))
end

end # module
