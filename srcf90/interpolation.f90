module interpolation
use error_messages! USEINDEPFILE
use precision_mod! USEINDEPFILE

implicit none

contains

subroutine resample1D( &
        ! this subroutine interpolates discrete y = f(x), the input and the output
        ! meshes need not have same minimum and maximum values. They have different number of 
        ! samples depending on the sampling interval
        ! only the cases allowed: these cases will result in truncation
                !minval(mesh_out) .gt. minval(mesh_in)
                !maxval(mesh_out) .lt. maxval(mesh_in)
        mesh_in_X, & ! input mesh
        mesh_in_dx, &
        mesh_in_nx, &
        mesh_out_X, & ! output mesh
        mesh_out_dx, &
        mesh_out_nx, &
        y_in, & ! y - values correponding to the input [nx_in]
        y_out, & ! y - values corresponding to the output [mesh_out_nx]
        flag, & ! flag is 1 for bilinear interpolation
                ! flag is 2 for smoothing cubic spline interpolation
        zero_flag& ! input y-values which are zero are not used for interpolation 
                ! For example if the velocity model has solid-fluid interface.. 
        ) bind(c, name="resample1D")
implicit none
integer, intent(in)                     :: mesh_in_nx, mesh_out_nx
real, intent(in)                        :: mesh_in_X(mesh_in_nx), mesh_out_X(mesh_out_nx)
real, intent(in)                        :: mesh_in_dx, mesh_out_dx
real, intent(in)                        :: y_in(mesh_in_nx)
real, intent(out)                       :: y_out(mesh_out_nx)
logical, intent(in), optional           :: zero_flag
integer, intent(in)                     :: flag

logical                                 :: zflag
real                                    :: y1, y2, x1, x2, x, x3, x4, y3, y4
integer                                 :: ixmin, ixmax, it

zflag = .false.
if(present(zero_flag)) then
        if(zero_flag) zflag = .true.
endif

! output to zero
y_out(:) = rzero_de;

if(flag .eq. 0) then
        do it = 1, mesh_out_nx
                call index_mesh1D(x0=mesh_out_X(it),mesh_X=mesh_in_X,mesh_dx=mesh_in_dx,&
                        mesh_nx=mesh_in_nx,ixmin=ixmin, ixmax=ixmax) 
                y1 = y_in(ixmin);
                y2 = y_in(ixmax);
                x1 = mesh_in_X(ixmin); x2 = mesh_in_X(ixmax); x = mesh_out_X(it)
                call interpolate_spray_B0_1D(x=x, y=y_out(it), x1=x1, y1=y1, x2=x2, y2=y2, flag=-1)
        enddo
elseif(flag .eq. 1) then
        do it = 1, mesh_out_nx
                call index_mesh1D(x0=mesh_out_X(it),mesh_X=mesh_in_X,mesh_nx=mesh_in_nx, &
                        mesh_dx=mesh_in_dx,ixmin=ixmin, ixmax=ixmax) 
                ! bilinear interpolation
                ! http://www.ajdesigner.com/phpinterpolation/bilinear_interpolation_equation.php
                y1 = y_in(ixmin);
                y2 = y_in(ixmax);
                x1 = mesh_in_X(ixmin); x2 = mesh_in_X(ixmax); x = mesh_out_X(it)
                if( ((y1 * y2) .eq. rzero_de) .and. (zflag) ) then
                        y_out(it) = rzero_de
                else
                        call interpolate_spray_B1_1D(x=x, y=y_out(it), x1=x1, y1=y1, x2=x2, y2=y2, flag=-1)
                endif
        enddo


! smoothing cubic spline interpolation
elseif(flag .eq. 2) then
        do it = 1, mesh_out_nx
                call index_mesh1D(x0=mesh_out_X(it),mesh_X=mesh_in_X,mesh_nx=mesh_in_nx, &
                        mesh_dx=mesh_in_dx,ixmin=ixmin, ixmax=ixmax) 

                y2 = y_in(ixmin);
                y3 = y_in(ixmax);
                x2 = mesh_in_X(ixmin); x = mesh_out_X(it)
                x3 = mesh_in_X(ixmax);
                x1 = x2 - (mesh_in_dx) 
                x4 = x3 + (mesh_in_dx); 
                if(ixmin .eq. 1) then
                        y4 = y_in(ixmax + 1);
                        y1 = rzero_de
                        
                        call interpolate_spray_B3_1D(x=x, y=y_out(it), x1=x1, y1=y1, x2=x2, y2=y2,&
                                x3=x3, x4=x4, y3=y3, y4=y4, flag=-2)

                elseif(ixmax .eq. mesh_in_nx) then
                        y4 = rzero_de
                        y1 = y_in(ixmin - 1);
                        
                        call interpolate_spray_B3_1D(x=x, y=y_out(it), x1=x1, y1=y1, x2=x2, y2=y2,&
                                x3=x3, x4=x4, y3=y3, y4=y4, flag=-3)
                else
                        y1 = y_in(ixmin - 1);
                        y4 = y_in(ixmax + 1);
                        
                        call interpolate_spray_B3_1D(x=x, y=y_out(it), x1=x1, y1=y1, x2=x2, y2=y2,&
                                x3=x3, x4=x4, y3=y3, y4=y4, flag=-1)
                endif

        enddo
endif

end subroutine resample1D

subroutine spray1D( &
        ! this subroutine sprays y_in on y_out
        ! x-vector have same minimum and maximum values, but have different number of 
        ! samples depending on the sampling interval
        mesh_in_X, & ! input mesh
        mesh_in_dx, &
        mesh_in_nx, &
        mesh_out_X, & ! output mesh
        mesh_out_dx, &
        mesh_out_nx, &
        y_in, & ! y - values correponding to the input [nx_in]
        y_out, & ! y - values corresponding to the output [mesh_out_nx]
        flag & ! flag is 1 for bilinear spray
               ! flag is 2 for smoothing cubic spline spray
        ) bind(c, name="spray1D")
implicit none
integer, intent(in)                     :: mesh_in_nx, mesh_out_nx
real, intent(in)                        :: mesh_in_X(mesh_in_nx), mesh_out_X(mesh_out_nx)
real, intent(in)                        :: mesh_in_dx, mesh_out_dx
real, intent(in)                        :: y_in(mesh_in_nx)
real, intent(out)                       :: y_out(mesh_out_nx)
integer, intent(in)                     :: flag

real                                    :: y1, y2, x1, x2, x, y, x3, x4, y3, y4
integer                                 :: ixmin, ixmax, it
real, allocatable                       :: y_out_extend(:)

! output to rzero_de
y_out(:) = rzero_de;

if(flag .eq. 0) then
        do it = 1, mesh_in_nx
                call index_mesh1D(x0=mesh_in_X(it),mesh_X=mesh_out_X, &
                        mesh_dx=mesh_out_dx,mesh_nx=mesh_out_nx,ixmin=ixmin, ixmax=ixmax) 

                ! find y1 and y2
                x1 = mesh_out_X(ixmin); x2 = mesh_out_X(ixmax); x = mesh_in_X(it)
                y = y_in(it);
                call interpolate_spray_B0_1D(x=x, y=y, x1=x1, y1=y1, x2=x2, y2=y2, flag=1)

                ! spray y1 and y2
                y_out(ixmin) = y_out(ixmin) + y1
                y_out(ixmax) = y_out(ixmax) + y2
        enddo
elseif(flag .eq. 1) then
        do it = 1, mesh_in_nx
                call index_mesh1D(x0=mesh_in_X(it),mesh_X=mesh_out_X, &
                        mesh_dx=mesh_out_dx,mesh_nx=mesh_out_nx,ixmin=ixmin, ixmax=ixmax) 

                ! find y1 and y2
                x1 = mesh_out_X(ixmin); x2 = mesh_out_X(ixmax); x = mesh_in_X(it)
                y = y_in(it);
                call interpolate_spray_B1_1D(x=x, y=y, x1=x1, y1=y1, x2=x2, y2=y2, flag=1)

                ! spray y1 and y2
                y_out(ixmin) = y_out(ixmin) + y1
                y_out(ixmax) = y_out(ixmax) + y2
        enddo
elseif(flag .eq. 2) then
        do it = 1, mesh_in_nx
                call index_mesh1D(x0=mesh_in_X(it),mesh_X=mesh_out_X, &
                        mesh_dx=mesh_out_dx,mesh_nx=mesh_out_nx,ixmin=ixmin, ixmax=ixmax) 

                x2 = mesh_out_X(ixmin); x = mesh_in_X(it); x3 = mesh_out_X(ixmax);
                x1 = x2 - (mesh_out_dx) 
                x4 = x3 + (mesh_out_dx); 
                y = y_in(it);

                if(ixmin .eq. 1) then
                        ! find y1, y2, y3 and y4
                        call interpolate_spray_B3_1D(x=x, y=y, x1=x1, y1=y1, x2=x2, y2=y2, flag=2, &
                                x3=x3, x4=x4, y3=y3, y4=y4)

                        ! spray y1, y2, y3  and y4
                        y_out_extend(ixmax + 1) = y_out_extend(ixmax + 1) + y4
                        y_out_extend(ixmin) = y_out_extend(ixmin) + y2
                        y_out_extend(ixmax) = y_out_extend(ixmax) + y3
                elseif(ixmax .eq. mesh_out_nx) then
                        ! find y1, y2, y3 and y4
                        call interpolate_spray_B3_1D(x=x, y=y, x1=x1, y1=y1, x2=x2, y2=y2, flag=3, &
                                x3=x3, x4=x4, y3=y3, y4=y4)

                        ! spray y1, y2, y3  and y4
                        y_out_extend(ixmin - 1) = y_out_extend(ixmin - 1) + y1
                        y_out_extend(ixmin) = y_out_extend(ixmin) + y2
                        y_out_extend(ixmax) = y_out_extend(ixmax) + y3
                else
                        ! find y1, y2, y3 and y4
                        call interpolate_spray_B3_1D(x=x, y=y, x1=x1, y1=y1, x2=x2, y2=y2, flag=1, &
                                x3=x3, x4=x4, y3=y3, y4=y4)

                        ! spray y1, y2, y3  and y4
                        y_out_extend(ixmin - 1) = y_out_extend(ixmin - 1) + y1
                        y_out_extend(ixmax + 1) = y_out_extend(ixmax + 1) + y4
                        y_out_extend(ixmin) = y_out_extend(ixmin) + y2
                        y_out_extend(ixmax) = y_out_extend(ixmax) + y3
                endif
        enddo
endif

end subroutine spray1D
!
subroutine interpolate_spray_B0_1D( &
! this subroutine interpolates or sprays to the nearest neighbour
! interpolation returns y using y1, y2
! spraying returns y1, y2 using y
        x1, x2, &
        x,  &
        y1, y2, &
        y, &
        flag &
        )
implicit none
real, intent(inout)                     :: x1, x2, y1, y2
real, intent(inout)                     :: x, y
integer, intent(in)                     :: flag
real                                    :: d1, d2

d1 = abs(x1 - x)
d2 = abs(x2 - x)

if(flag .eq. 1) then
        if(d1 .le. d2) then
                y1 = y; y2 = rzero_de
        else
                y2 = y; y1 = rzero_de
        endif
        return
elseif(flag .eq. -1) then
        if(d1 .le. d2) then
                y = y1;
        else
                y = y2;
        endif
        return
endif

end subroutine interpolate_spray_B0_1D
!
!
subroutine interpolate_spray_B1_1D( &
! this subroutine interpolates or sprays bilinearly [one is adjoint of another]
! interpolation returns y using y1, y2
! spraying returns y1, y2 using y
! 
!                        +                      
!                        |                      
!    y1= f(x1)           |      y2= f(x2)   
!      +-----------------x--------+             
!                 y=f(x) |                      
!                        +                      
        x1, x2, &
        x,  &
        y1, y2, &
        y, &
        flag &
        )
implicit none
real, intent(inout)                     :: x1, x2, y1, y2
real, intent(inout)                     :: x, y
integer, intent(in)                     :: flag
real                                    :: denom

! bilinear interpolation
! http://www.ajdesigner.com/phpinterpolation/bilinear_interpolation_equation.php
denom = (x2 - x1)
if(flag .eq. 1) then
        if((x2-x) .eq. rzero_de) then
        else
                y1 = y * (x2-x)/denom
        endif
        if((x-x1) .eq. rzero_de) then
                y2 = rzero_de;
        else
                y2 = y * (x-x1)/denom
        endif
        return
elseif(flag .eq. -1) then
        y = (y1*(x2-x) + y2*(x-x1)) / denom
        return
endif

end subroutine interpolate_spray_B1_1D
!
!
!
subroutine interpolate_spray_B3_1D( &
! this subroutine interpolates or sprays 
! using cubic bspline
! interpolation returns y using y1, y2, y3, y4
! spraying returns y1, y2, y3, y4 using y
! 
!                        +                      
!                        |                      
!    y1      y2          |        y3       y4
!      +-------+---------x--------+--------+             
!                 y=f(x) |                      
!                        +                      
        x1, x2, &
        x3, x4, &
        x,  &
        y1, y2, &
        y3, y4, &
        y, &
        flag &
        )
implicit none
real, intent(in)                        :: x1, x2, x3, x4
real, intent(inout)                     :: y1, y2, y3, y4
real, intent(in)                        :: x
real, intent(inout)                     :: y
integer, intent(in)                     :: flag
real                                    :: h, c1, c2, c3, c4
real                                    :: hcube, hcubeI


h = (x2 - x1);
hcube = h*h*h;
hcubeI = hcube**(-1.0)
if(h .eq. rzero_de) then
         if(flag .eq. 1) then
                y1 = y 
                y2 = y 
                y3 = y 
                y4 = y 
                return
        elseif(flag .eq. -1) then
                y = 0.25 * (y1 + y2 + y3 + y4)

                return
        endif
else
        c1 = (1.0 / 6.0 * (2.0 * h - (x - x1))**3) * hcubeI
        c4 = (1.0 / 6.0 * (2.0 * h - (x4 - x))**3) * hcubeI
        c2 = ((2.0 / 3.0 * hcube - (0.5 * (x2 - x)**2.0 * (2.0 * h - (x - x2)) ) )) * hcubeI
        c3 = ((2.0 / 3.0 * hcube - (0.5 * (x3 - x)**2.0 * (2.0 * h + (x - x3)) ) )) * hcubeI
        ! center spray
        if(flag .eq. 1) then
                y1 = y * c1 
                y2 = y * c2 
                y3 = y * c3 
                y4 = y * c4 
                return
        ! center interpolate
        elseif(flag .eq. -1) then
                y = (y1 * c1 + y2 * c2 + y3 * c3 + y4 * c4)
                return
        ! left edge spray
        elseif(flag .eq. 2) then
                y2 = y * (c2 + 2.0*c1) 
                y3 = y * (c3 - c1)
                y4 = y * c4 
                return
        ! left edge interpolate
        elseif(flag .eq. -2) then
                y = (y2 * (c2 + 2.0*c1) + y3 * (c3 - c1) + y4 * c4)
                return
        ! right edge spray        
        elseif(flag .eq. 3) then
                y1 = y * c1 
                y2 = y * (c2 - c4) 
                y3 = y * (c3+ 2.0 * c4) 
                return
        ! right edge interpolate
        elseif(flag .eq. -3) then
                y = (y1 * c1 + y2 * (c2 - c4) + y3 * (c3 + 2.0*c4))
                return
        endif
endif

end subroutine interpolate_spray_B3_1D
!
subroutine resample2D( &
        ! y_out = interpolation(y_in)
        ! this subroutine interpolates a discrete y = f(z, x), the input and the output
        ! model matrices have same minimum and maximum values, but have different number of 
        ! samples depending on the sampling interval in both directions
        mesh_in_X, & ! input mesh
        mesh_in_Z, &
        mesh_in_dx, &
        mesh_in_dz, &
        mesh_in_nx, &
        mesh_in_nz, &
        mesh_out_X, & ! output mesh
        mesh_out_Z, &
        mesh_out_dx, & 
        mesh_out_dz, &
        mesh_out_nx, &
        mesh_out_nz, &
        y_in, & ! y - values correponding to the input [nz_in, nx_in]
        y_out, & ! y - values corresponding to the output [nz_out, nx_out]
        zero_flag, & ! input y-values which are zero are not used for interpolation 
                ! For example if the velocity model has solid-fluid interface.. 
        flag & ! bilinear interpolation or cubic spline interpolation
        ) bind(c, name="resample2D")
implicit none
integer, intent(in)                     :: mesh_in_nx, mesh_in_nz
real, intent(in)                        :: mesh_in_dx, mesh_in_dz
real, intent(in)                        :: mesh_in_X(mesh_in_nx), mesh_in_Z(mesh_in_nz)
integer, intent(in)                     :: mesh_out_nx, mesh_out_nz
real, intent(in)                        :: mesh_out_dx, mesh_out_dz
real, intent(in)                        :: mesh_out_X(mesh_out_nx), mesh_out_Z(mesh_out_nz)
logical, intent(in), optional           :: zero_flag
logical                                 :: zflag
real, intent(in)                        :: y_in(mesh_in_nz,mesh_in_nx)
real, intent(out)                       :: y_out(mesh_out_nz,mesh_out_nx)
real                                    :: y11, y12, y22, y21, x1, x2, z1, z2, z, x
real                                    :: y1, y2, y3, y4
real                                    :: x3, x4, z3, z4
integer                                 :: ixmin, ixmax, izmax, izmin, ix, iz, iix, iiz
real, allocatable                       :: ymat(:,:), y_in_x(:,:)
integer, intent(in)                     :: flag

zflag = .false.
if(present(zero_flag)) then
        if(zero_flag) zflag = .true.
endif


! output to zero
y_out(:, :) = rzero_de;

if(flag .eq. 0) then
        do ix = 1, mesh_out_nx
                do iz = 1, mesh_out_nz
                        call index_mesh2D(mesh_X=mesh_in_X,mesh_Z=mesh_in_Z,&
                                mesh_nx=mesh_in_nx,mesh_nz=mesh_in_nz,&
                                x0=mesh_out_X(ix),z0=mesh_out_Z(iz),&
                                ixmin=ixmin,ixmax=ixmax,izmin=izmin,izmax=izmax)

                        y11 = y_in(izmin, ixmin);
                        y12 = y_in(izmin, ixmax);
                        y21 = y_in(izmax, ixmin);
                        y22 = y_in(izmax, ixmax);

                        x1 = mesh_in_X(ixmin); x2 = mesh_in_X(ixmax); x = mesh_out_X(ix)
                        z1 = mesh_in_Z(izmin); z2 = mesh_in_Z(izmax); z = mesh_out_Z(iz)

                        call interpolate_spray_B0_2D(y11=y11,y12=y12,y21=y21,y22=y22,y=y_out(iz,ix),&
                                x=x,z=z,x1=x1,x2=x2,z1=z1,z2=z2,flag=-1)

               enddo
        enddo
elseif(flag .eq. 1) then
        do ix = 1, mesh_out_nx
                do iz = 1, mesh_out_nz
                        call index_mesh2D(mesh_X=mesh_in_X,mesh_Z=mesh_in_Z,&
                                mesh_nx=mesh_in_nx,mesh_nz=mesh_in_nz,&
                                x0=mesh_out_X(ix),z0=mesh_out_Z(iz),&
                                ixmin=ixmin,ixmax=ixmax,izmin=izmin,izmax=izmax)

                        y11 = y_in(izmin, ixmin);
                        y12 = y_in(izmin, ixmax);
                        y21 = y_in(izmax, ixmin);
                        y22 = y_in(izmax, ixmax);

                        x1 = mesh_in_X(ixmin); x2 = mesh_in_X(ixmax); x = mesh_out_X(ix)
                        z1 = mesh_in_Z(izmin); z2 = mesh_in_Z(izmax); z = mesh_out_Z(iz)

                        if( ((y11 * y12 * y21 *y22) .eq. rzero_de) .and. (zflag) ) then
                                y_out(iz,ix) = rzero_de;
                        else
                                ! bilinear interpolation
                                call interpolate_spray_B1_2D(y11=y11,y12=y12,y21=y21,y22=y22,y=y_out(iz,ix),&
                                        x=x,z=z,x1=x1,x2=x2,z1=z1,z2=z2,flag=-1)
                        endif

               enddo
        enddo

elseif(flag .eq. 2) then
        ! y interpolation along x-direction
        allocate(y_in_x(mesh_in_nz, mesh_out_nx)); y_in_x = rzero_de
        y1 = rzero_de; y2 = rzero_de; y3 = rzero_de; y4 = rzero_de;
        if(mesh_in_nx .eq. 1 .or. mesh_in_nz .eq. 1) call abort_msg("resample2D: cubic spline interpolation not valid in this case")
        do iz = 1, mesh_in_nz
                do ix = 1, mesh_out_nx
                        call index_mesh2D(mesh_X=mesh_in_X,mesh_Z=mesh_in_Z,&
                                mesh_nx=mesh_in_nx,mesh_nz=mesh_in_nz,&
                                x0=mesh_out_X(ix),z0=mesh_in_Z(iz),&
                                ixmin=ixmin,ixmax=ixmax,izmin=izmin,izmax=izmax)
                        ! smoothing cubic spline interpolation
                        y2 = y_in(iz,ixmin);
                        y3 = y_in(iz,ixmax);

                        x2 = mesh_in_X(ixmin); x = mesh_out_X(ix)
                        x3 = mesh_in_X(ixmax);
                        x1 = x2 - (mesh_in_dx) 
                        x4 = x3 + (mesh_in_dx); 

                        if(ixmin .eq. 1) then
                                y1 = rzero_de;
                                y4 = y_in(iz,ixmax + 1);
                                call interpolate_spray_B3_1D(x=x, y=y_in_x(iz, ix), x1=x1, y1=y1, x2=x2, y2=y2,&
                                        x3=x3, x4=x4, y3=y3, y4=y4, flag=-2)
                        elseif(ixmax .eq. mesh_in_nx) then
                                y1 = y_in(iz,ixmin - 1);
                                y4 = rzero_de;
                                call interpolate_spray_B3_1D(x=x, y=y_in_x(iz, ix), x1=x1, y1=y1, x2=x2, y2=y2,&
                                        x3=x3, x4=x4, y3=y3, y4=y4, flag=-3)
                        else
                                y1 = y_in(iz,ixmin - 1);
                                y4 = y_in(iz,ixmax + 1);
                                call interpolate_spray_B3_1D(x=x, y=y_in_x(iz, ix), x1=x1, y1=y1, x2=x2, y2=y2,&
                                        x3=x3, x4=x4, y3=y3, y4=y4, flag=-1)
                        endif
                enddo
        enddo

        y1 = rzero_de; y2 = rzero_de; y3 = rzero_de; y4 = rzero_de;
        ! interpolation along z direction
        do iz = 1, mesh_out_nz
                do ix = 1, mesh_out_nx

                        call index_mesh2D(mesh_X=mesh_in_X,mesh_Z=mesh_in_Z,&
                                mesh_nx=mesh_in_nx,mesh_nz=mesh_in_nz,&
                                x0=mesh_out_X(ix),z0=mesh_out_Z(iz),&
                                ixmin=ixmin,ixmax=ixmax,izmin=izmin,izmax=izmax)

                        z2 = mesh_in_Z(izmin); x = mesh_out_Z(iz)
                        z3 = mesh_in_Z(izmax);
                        z1 = z2 - (mesh_in_dz) 
                        z4 = z3 + (mesh_in_dz); 
                      
                        y2 = y_in_x(izmin,ix);
                        y3 = y_in_x(izmax,ix);

                        if(izmin .eq. 1) then
                                y1 = rzero_de;
                                y4 = y_in_x(izmax + 1,ix);

                                call interpolate_spray_B3_1D(x=x, y=y_out(iz, ix), x1=z1, y1=y1, x2=z2, y2=y2,&
                                        x3=z3, x4=z4, y3=y3, y4=y4, flag=-2)
                        elseif(izmax .eq. mesh_in_nz) then
                                y1 = y_in_x(izmin - 1,ix);
                                y4 = rzero_de;

                                call interpolate_spray_B3_1D(x=x, y=y_out(iz, ix), x1=z1, y1=y1, x2=z2, y2=y2,&
                                        x3=z3, x4=z4, y3=y3, y4=y4, flag=-3)
                        else
                                y1 = y_in_x(izmin - 1,ix);
                                y4 = y_in_x(izmax + 1,ix);

                                call interpolate_spray_B3_1D(x=x, y=y_out(iz, ix), x1=z1, y1=y1, x2=z2, y2=y2,&
                                        x3=z3, x4=z4, y3=y3, y4=y4, flag=-1)
                        endif
                enddo
        enddo
        deallocate(y_in_x)
endif


end subroutine resample2D

subroutine spray2D( &
        ! y_out = spray(y_in)
        ! this subroutine sprays input onto an output mesh 
        ! model matrices have same minimum and maximum values, but have different number of 
        ! samples depending on the sampling interval in both directions
        mesh_in_X, & ! input mesh
        mesh_in_Z, &
        mesh_in_dx, &
        mesh_in_dz, &
        mesh_in_nx, &
        mesh_in_nz, &
        mesh_out_X, & ! output mesh
        mesh_out_Z, &
        mesh_out_dx, & 
        mesh_out_dz, &
        mesh_out_nx, &
        mesh_out_nz, &
        y_in, & ! y - values correponding to the input [nz_in, nx_in]
        y_out, & ! y - values corresponding to the output [nz_out, nx_out]
        flag & 
        ) bind(c, name="spray2D")
implicit none
integer, intent(in)                     :: mesh_in_nx, mesh_in_nz
real, intent(in)                        :: mesh_in_dx, mesh_in_dz
real, intent(in)                        :: mesh_in_X(mesh_in_nx), mesh_in_Z(mesh_in_nz)
integer, intent(in)                     :: mesh_out_nx, mesh_out_nz
real, intent(in)                        :: mesh_out_dx, mesh_out_dz
real, intent(in)                        :: mesh_out_X(mesh_out_nx), mesh_out_Z(mesh_out_nz)
real, intent(in)                        :: y_in(:,:)
real, intent(out)                       :: y_out(:,:)
real                                    :: y11, y12, y22, y21, x1, x2, z1, z2, z, x, y
real                                    :: x3, x4, z3, z4
real                                    :: y1, y2, y3, y4
real, allocatable                       :: ymat(:,:), y_out_z(:,:)
integer                                 :: ixmin, ixmax, izmax, izmin, ix, iz
integer, intent(in)                     :: flag

!if(minval(mesh_in_X) .ne. minval(mesh_out_X)) call abort_msg('spray2D: meshes have different Xmin values')
!if(maxval(mesh_in_X) .ne. maxval(mesh_out_X)) call abort_msg('spray2D: meshes have different Xmax values')
!if(minval(mesh_in_Z) .ne. minval(mesh_out_Z)) call abort_msg('spray2D: meshes have different Zmin values')
!if(maxval(mesh_in_Z) .ne. maxval(mesh_out_Z)) call abort_msg('spray2D: meshes have different Zmax values')

call check_dimension('spray2D: y_in size',y_in, mesh_in_nz,mesh_in_nx)
call check_dimension('spray2D: y_out size',y_out, mesh_out_nz,mesh_out_nx)

! output to zero
y_out(:, :) = rzero_de;

if(flag .eq. 0) then
        do ix = 1, mesh_in_nx
                do iz = 1, mesh_in_nz
                        call index_mesh2D(mesh_X=mesh_out_X,mesh_Z=mesh_out_Z,&
                                mesh_nx=mesh_out_nx,mesh_nz=mesh_out_nz,&
                                x0=mesh_in_X(ix),z0=mesh_in_Z(iz),&
                                ixmin=ixmin,ixmax=ixmax,izmin=izmin,izmax=izmax)

                        x1 = mesh_out_X(ixmin); x2 = mesh_out_X(ixmax); x = mesh_in_X(ix)
                        z1 = mesh_out_Z(izmin); z2 = mesh_out_Z(izmax); z = mesh_in_Z(iz)

                        ! find y11, y12, y21, y22
                        y = y_in(iz,ix);
                        call interpolate_spray_B0_2D(y11=y11,y12=y12,y21=y21,y22=y22,y=y,&
                                x=x,z=z,x1=x1,x2=x2,z1=z1,z2=z2,flag=1)

                        ! spray them
                        y_out(izmin, ixmin) = y_out(izmin, ixmin) + y11;
                        y_out(izmin, ixmax) = y_out(izmin, ixmax) + y12;
                        y_out(izmax, ixmin) = y_out(izmax, ixmin) + y21;
                        y_out(izmax, ixmax) = y_out(izmax, ixmax) + y22;
               enddo
        enddo
elseif(flag .eq. 1) then
        do ix = 1, mesh_in_nx
                do iz = 1, mesh_in_nz
                        call index_mesh2D(mesh_X=mesh_out_X,mesh_Z=mesh_out_Z,&
                                mesh_nx=mesh_out_nx,mesh_nz=mesh_out_nz,&
                                x0=mesh_in_X(ix),z0=mesh_in_Z(iz),&
                                ixmin=ixmin,ixmax=ixmax,izmin=izmin,izmax=izmax)

                        x1 = mesh_out_X(ixmin); x2 = mesh_out_X(ixmax); x = mesh_in_X(ix)
                        z1 = mesh_out_Z(izmin); z2 = mesh_out_Z(izmax); z = mesh_in_Z(iz)

                        ! find y11, y12, y21, y22
                        y = y_in(iz,ix);
                        call interpolate_spray_B1_2D(y11=y11,y12=y12,y21=y21,y22=y22,y=y,&
                                x=x,z=z,x1=x1,x2=x2,z1=z1,z2=z2,flag=1)

                        ! spray them
                        y_out(izmin, ixmin) = y_out(izmin, ixmin) + y11;
                        y_out(izmin, ixmax) = y_out(izmin, ixmax) + y12;
                        y_out(izmax, ixmin) = y_out(izmax, ixmin) + y21;
                        y_out(izmax, ixmax) = y_out(izmax, ixmax) + y22;
               enddo
        enddo

elseif(flag .eq. 2) then

        if(mesh_in_nx .eq. 1 .or. mesh_in_nz .eq. 1) call abort_msg("spray2D: cubic spline interpolation not valid in this case")
        ! spray values along z direction
        allocate(y_out_z(mesh_out_nz,mesh_in_nx)); y_out_z = rzero_de
        do ix = 1, mesh_in_nx
                do iz = 1, mesh_in_nz

                        call index_mesh2D(mesh_X=mesh_out_X,mesh_Z=mesh_out_Z,&
                                mesh_nx=mesh_out_nx,mesh_nz=mesh_out_nz,&
                                x0=mesh_in_X(ix),z0=mesh_in_Z(iz),&
                                ixmin=ixmin,ixmax=ixmax,izmin=izmin,izmax=izmax)

                        z2 = mesh_out_Z(izmin); z = mesh_in_Z(iz); z3 = mesh_out_Z(izmax);
                        z1 = z2 - (mesh_out_dz) 
                        z4 = z3 + (mesh_out_dz); 
                        y = y_in(iz,ix);


                        if(izmin .eq. 1) then
                                ! find y1, y2, y3 and y4
                                call interpolate_spray_B3_1D(x=z, y=y, x1=z1, y1=y1, x2=z2, y2=y2, flag=2, &
                                        x3=z3, x4=z4, y3=y3, y4=y4)
                                ! spray y1, y2, y3  and y4
                                y_out_z(izmax + 1, ix) = y_out_z(izmax + 1, ix) + y4
                                y_out_z(izmin, ix) = y_out_z(izmin, ix) + y2
                                y_out_z(izmax, ix) = y_out_z(izmax, ix) + y3

                        elseif(izmax .eq. mesh_out_nz) then
                                ! find y1, y2, y3 and y4
                                call interpolate_spray_B3_1D(x=z, y=y, x1=z1, y1=y1, x2=z2, y2=y2, flag=3, &
                                        x3=z3, x4=z4, y3=y3, y4=y4)
                                ! spray y1, y2, y3  and y4
                                y_out_z(izmin - 1, ix) = y_out_z(izmin - 1, ix) + y1
                                y_out_z(izmin, ix) = y_out_z(izmin, ix) + y2
                                y_out_z(izmax, ix) = y_out_z(izmax, ix) + y3

                        else
                                ! find y1, y2, y3 and y4
                                call interpolate_spray_B3_1D(x=z, y=y, x1=z1, y1=y1, x2=z2, y2=y2, flag=1, &
                                        x3=z3, x4=z4, y3=y3, y4=y4)
                                ! spray y1, y2, y3  and y4
                                y_out_z(izmin - 1, ix) = y_out_z(izmin - 1, ix) + y1
                                y_out_z(izmax + 1, ix) = y_out_z(izmax + 1, ix) + y4
                                y_out_z(izmin, ix) = y_out_z(izmin, ix) + y2
                                y_out_z(izmax, ix) = y_out_z(izmax, ix) + y3
                        endif
                enddo
        enddo

        ! spray values along x direction
        do iz = 1, mesh_out_nz
                do ix = 1, mesh_in_nx

                        call index_mesh2D(mesh_X=mesh_out_X,mesh_Z=mesh_out_Z,&
                                mesh_nx=mesh_out_nx,mesh_nz=mesh_out_nz,&
                                x0=mesh_in_X(ix),z0=mesh_out_Z(iz),&
                                ixmin=ixmin,ixmax=ixmax,izmin=izmin,izmax=izmax)

                        x2 = mesh_out_X(ixmin); x = mesh_in_X(ix); x3 = mesh_out_X(ixmax);
                        x1 = x2 - (mesh_out_dx) 
                        x4 = x3 + (mesh_out_dx); 
                        y = y_out_z(iz,ix);

                        if(ixmin .eq. 1) then
                                ! find y1, y2, y3 and y4
                                call interpolate_spray_B3_1D(x=x, y=y, x1=x1, y1=y1, x2=x2, y2=y2, flag=2, &
                                        x3=x3, x4=x4, y3=y3, y4=y4)

                                ! spray y1, y2, y3  and y4
                                y_out(iz, ixmax + 1) = y_out(iz, ixmax + 1) + y4
                                y_out(iz, ixmin) = y_out(iz, ixmin) + y2
                                y_out(iz, ixmax) = y_out(iz, ixmax) + y3
                        elseif(ixmax .eq. mesh_out_nx) then
                                ! find y1, y2, y3 and y4
                                call interpolate_spray_B3_1D(x=x, y=y, x1=x1, y1=y1, x2=x2, y2=y2, flag=3, &
                                        x3=x3, x4=x4, y3=y3, y4=y4)

                                ! spray y1, y2, y3  and y4
                                y_out(iz, ixmin - 1) = y_out(iz, ixmin - 1) + y1
                                y_out(iz, ixmin) = y_out(iz, ixmin) + y2
                                y_out(iz, ixmax) = y_out(iz, ixmax) + y3
                        else
                                ! find y1, y2, y3 and y4
                                call interpolate_spray_B3_1D(x=x, y=y, x1=x1, y1=y1, x2=x2, y2=y2, flag=1, &
                                        x3=x3, x4=x4, y3=y3, y4=y4)

                                ! spray y1, y2, y3  and y4
                                y_out(iz, ixmin - 1) = y_out(iz, ixmin - 1) + y1
                                y_out(iz, ixmax + 1) = y_out(iz, ixmax + 1) + y4
                                y_out(iz, ixmin) = y_out(iz, ixmin) + y2
                                y_out(iz, ixmax) = y_out(iz, ixmax) + y3
                        endif
                enddo
        enddo
        deallocate(y_out_z)

endif

end subroutine spray2D
!

subroutine interpolate_spray_B1_2D( &
! this subroutine interpolates or sprays bilinearly [one is adjoint of another]
! interpolation returns y using y11, y12, y22, y21
! spraying returns y11, y12, y21, y22 using y
!                        +                      
!                        |                      
!    y11= f(z1,x1)       |      y12= f(z1,x2)   
!      +--------------------------+             
!      |                 |        |             
!      |                 |        |             
!      |                 |        |             
!      |                 |        |             
!      |                 |        |             
!      |                 |        |             
!      |                 |        |             
!      |                 |        |             
!      |        y= f(z,x)|        |             
!+-----------------------*---------------------+
!      |                 |        |             
!      |                 |        |             
!      +--------------------------+             
!  y21= f(z2,x1)         |       y22= f(z2,x2)  
!                        |                      
!                        |                      
!                        +                      
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
real                                    :: denom
integer, intent(in)                     :: flag

if((x2 - x1)*(z2 - z1) .ne. rzero_de) then
        ! bilinear interpolation
        ! http://www.ajdesigner.com/phpinterpolation/bilinear_interpolation_equation.php
        denom = (x2 - x1)*(z2 - z1)
        if(flag .eq. 1) then
                y11 = y * (x2-x)*(z2-z)/denom
                y12 = y * (x-x1)*(z2-z)/denom
                y21 = y * (x2-x)*(z-z1)/denom
                y22 = y * (x-x1)*(z-z1)/denom
                return
        elseif(flag .eq. -1) then
                y = (y11*(x2-x)*(z2-z) + y12*(x-x1)*(z2-z) + &
                    y21*(x2-x)*(z-z1) + y22*(x-x1)*(z-z1)) / denom
                return
        endif
elseif((x2 - x1) .eq. rzero_de) then
        denom = (z2 - z1)
        if(flag .eq. 1) then
                y11 = y * (z2-z)/denom
                y12 = y * (z2-z)/denom
                y21 = y * (z-z1)/denom
                y22 = y * (z-z1)/denom
                return
        elseif(flag .eq. -1) then
                y = (y11*(z2-z) + &
                    y21*(z-z1)) / denom
                return
        endif
elseif((z2 - z1) .eq. rzero_de) then
        denom = (x2 - x1)
        if(flag .eq. 1) then
                y11 = y * (x2-x)/denom
                y12 = y * (x-x1)/denom
                y21 = y * (x2-x)/denom
                y22 = y * (x-x1)/denom
                return
        elseif(flag .eq. -1) then
                y = (y11*(x2-x) + y12*(x-x1)) / denom
                return
        endif
endif

end subroutine interpolate_spray_B1_2D

subroutine interpolate_spray_B0_2D( &
! this subroutine interpolates or sprays by nearest neighbour
! interpolation returns y using y11, y12, y22, y21
! spraying returns y11, y12, y21, y22 using y
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

end subroutine interpolate_spray_B0_2D

subroutine index_mesh1D( &
         ! subroutine to find indices of a points around x0 on mesh1D
                 mesh_X, &
                 mesh_dx, &
                 mesh_nx, &
                 x0, & ! coordinates of the point 
                         ixmin, ixmax & ! 
                       )
implicit none
integer, intent(in)                             :: mesh_nx
real, intent(in)                                :: mesh_X(mesh_nx), mesh_dx
integer, intent(out)                            :: ixmin, ixmax
real, intent(in)                                :: x0
logical, allocatable                            :: mask1D(:)
integer                                         :: i1, i2
!                        +                      
!                        |                      
!         (tmin)         |           (tmax)   
!      +-----------------x--------+             
!                   (x0) |                      
!                        +                      

allocate(mask1D(mesh_nx)); mask1D = .true.
i1 = minloc(array = abs(mesh_X - x0), dim = 1, mask = mask1D)
mask1D(i1) = .false.
i2 = minloc(array = abs(mesh_X - x0), dim = 1, mask = mask1D)
ixmin = min(i1, i2); ixmax = max(i1, i2)
deallocate(mask1D)
end subroutine index_mesh1D


subroutine index_mesh2D( &
                        mesh_X, mesh_Z, & ! mesh2D
                        mesh_nx, mesh_nz, &
                        x0, z0, & ! coordinates of the point 
                        ixmin, ixmax, izmin, izmax & ! 
                       )
implicit none
integer, intent(in)                     :: mesh_nx, mesh_nz
real, intent(in)                        :: mesh_X(mesh_nx), mesh_Z(mesh_nz)
integer, intent(out)                    :: ixmin, ixmax, izmin, izmax
real, intent(in)                        :: x0, z0
logical, allocatable                    :: mask1D(:)
integer                                 :: i1, i2
!                        +                      
!                        |                      
!          (z1,x1)       |            (z1,x2)   
!      +--------------------------+             
!      |                 |        |             
!      |                 |        |             
!      |                 |        |             
!      |                 |        |             
!      |                 |        |             
!      |                 |        |             
!      |                 |        |             
!      |                 |        |             
!      |          (z0,x0)|        |             
!+-----------------------*---------------------+
!      |                 |        |             
!      |                 |        |             
!      +--------------------------+             
!        (z2,x1)         |             (z2,x2)  
!                        |                      
!                        |                      
!                        +                      
if(mesh_nx .eq. 1) then
    239         ixmin = 1; ixmax = 1;
else
        allocate(mask1D(mesh_nx)); mask1D = .true.
        i1 = minloc(array = abs(mesh_X(:) - x0), dim = 1, mask = mask1D)
        mask1D(i1) = .false.
        i2 = minloc(array = abs(mesh_X(:) - x0), dim = 1, mask = mask1D)
        ixmin = min(i1, i2); ixmax = max(i1, i2)
        deallocate(mask1D)
endif

if(mesh_nz .eq. 1) then
    250         izmin = 1; izmax = 1;
else
        allocate(mask1D(mesh_nz)); mask1D = .true.
        i1 = minloc(array = abs(mesh_Z(:) - z0), dim = 1, mask = mask1D)
        mask1D(i1) = .false.
        i2 = minloc(array = abs(mesh_Z(:) - z0), dim = 1, mask = mask1D)
        izmin = min(i1, i2); izmax = max(i1, i2)
        deallocate(mask1D)
endif
end subroutine index_mesh2D


end module interpolation

