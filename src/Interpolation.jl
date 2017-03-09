module Interpolation







function interp_model(;modelin::Seismic=Seismic("test_homo_acou"), modelout=Seis
	 meshgridin::Grid.Mod2D = Grid.Mod2D(), 
	 meshgridout::Grid.Mod2D = Grid.Mod2D(), 


	ccall( (:resample2D, F90libs.fdtd), Void, 
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

	      meshgridin.x, 
	      meshgridin.z,  
	      meshgridin.δx,  
	      meshgridin.δz,  
	      meshgridin.nx,  
	      meshgridin.nz,  
	      meshgridout.x, 
	      meshgridout.z,  
	      meshgridout.δx,  
	      meshgridout.δz,  
	      meshgridout.nx,  
	      meshgridout.nz,  
	      y_in, 
	      y_out, 
	      zero_flag, 
	      flag 
	     )

      jobname, npropwav,
      Models.χ(model.χKI, model.K0I,-1), 
      Models.χ(model.χρI, model.ρ0I,-1), 
      mgrid.nx, mgrid.nz,
      mgrid.δx, mgrid.δz, 
 


end # module
