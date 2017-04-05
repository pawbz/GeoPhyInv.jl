module Models_Interpolation



"""
this subroutine interpolates or sprays by nearest neighbour
interpolation returns y using y11, y12, y22, y21
spraying returns y11, y12, y21, y22 using y

"""
function interpolate_spray_B0_2D( &
        x1, x2, &
        z1, z2, &
        z,  x,  &
        y11, y12, &
        y21, y22, &
        y, &
        flag &
        )
implicit none
real, intent(inout)                     :: y11, y12, y21, y22, y
real, intent(in)                        :: x, z, x1, z1, x2, z2
real                                    :: d(4), ytemp(4)
integer, intent(in)                     :: flag

d(1) = sqrt((x-x1)**2 + (z-z1)**2)
d(2) = sqrt((x-x2)**2 + (z-z1)**2)
d(3) = sqrt((x-x2)**2 + (z-z2)**2)
d(4) = sqrt((x-x1)**2 + (z-z2)**2)
ytemp(4) = rzero_de;

if(flag .eq. 1) then
        ytemp = rzero_de
        ytemp(minloc(d)) = y;
        y11 = ytemp(1); y12 = ytemp(2); y22 = ytemp(3); y21 = ytemp(4)
        return
elseif(flag .eq. -1) then
        ytemp(1) = y11; ytemp(2) = y12; ytemp(3) = y22; ytemp(4) = y21;
        y = sum(ytemp(minloc(d)))
        return
endif

return 

end subroutine interpolate_spray_B0_2D








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
