module fdtd




function __init__()

end

"""
modeling
"""

function fdtd_mod(
		  jobname="Hello",
		  npropwav=1, 
		  mesh_nx = 200,
		  mesh_nz = 200,
		  mesh_X = linspace(0,2000,mesh_nx),
		  mesh_Z = linspace(0,2000,mesh_nz),
		  modTT = 2000.0 * ones(mesh_nz,mesh_nx),
		  modRR = 2000.0 * ones(mesh_nz,mesh_nx),
		 )





ccall( (:fdtd_mod, SeismicInversion.f90libs.fdtd), Void, 
      (Ptr{UInt8}, Ptr{Float64}, Ptr{Float64}, Ref{Int64}, 
							       Ref{Int64}), 
      jobname, B, n1, n2)

end

end # module
