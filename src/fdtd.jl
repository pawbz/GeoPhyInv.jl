module fdtd

using SeismicInversion.mesh
using SeismicInversion.models
using SeismicInversion.f90libs


function __init__()

end

"""
modeling
"""
function fdtd_mod(
		  jobname="Hello",
		  npropwav=1, 
		  mesh::Mesh2D = Mesh2D(), 
		  model::SeismicModel = SeismicModel("test_homo_acoustic", mesh),
		  tim_nt = 2000, tim_del = 0.001,
		  src_nt = 2000,
		  src_wavelets = randn(2000),
		  src_nseq = 1, src_nsmul = 1,
		  src_flags = "Hello",
		  src_x = zeros(1), src_z = zeros(1),
		  recv_n = 1,
		  recv_x = zeros(1),
		  recv_z = zeros(1),
		  recv_flags = "Hello"
		 )

#if(tim_nt .lt. src_nt) call abort_msg("fd2_mod: tim_nt .lt. src_nt")

#! no modeling if source wavelet is zero
#if(maxval(abs(src_wavelets)) .lt. tiny(rzero_de)) then
#        return
#endif


	recv_out = zeros(tim_nt)
	snaps_in = zeros(tim_nt)
	snaps_out = ones(tim_nt)
	grad_modtt = zeros(mesh.nz, mesh.nx)
	grad_modrr = zeros(mesh.nz, mesh.nx)


ccall( (:fdtd_mod, f90libs.fdtd), Void, 
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
      χ(model.χKI, model.K0I,-1), 
      χ(model.χρI, model.ρ0I,-1), 
      mesh.nx, mesh.nz,
      mesh.δx, mesh.δz, 
      mesh.x, mesh.z,
      mesh.npml, "mesh_abs_trbl",
      tim_nt, tim_del, 
      src_nt, src_wavelets,
      src_nseq, src_nsmul,
      src_x, src_z,
      src_flags,
      recv_n, recv_x, recv_z,
      recv_out,
      recv_flags,
      snaps_in,
      snaps_out,
      grad_modtt,
      grad_modrr,
     )

end

end # module
