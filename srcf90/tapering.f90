module tapering
use precision_mod! USEINDEPFILE
use error_messages! USEINDEPFILE
implicit none

contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! subroutines for tapering and windowing
! go_window
function window2D(nz, nx, nx0, nz0, Nxlen, Nzlen, alphax1, alphax2, alphaz1, alphaz2)

! this function returns a matrix with a 2D window in it.
! the cosine tapered window is centered at [nz0, nx0] and has total width of Nlen
! the output has the matrix has dimensions [nz, nx]
! alpha controls the steepness of the edges
! if Nlen is large such that the window exceeds the window, then window is cut
implicit none
integer                                 :: ix, iz, ixx, izz
integer, intent(in)                     :: nz, nx, nz0, nx0, Nxlen, Nzlen
real, intent(in)                        :: alphax1, alphax2, alphaz1, alphaz2
real                                    :: alpx1, alpx2, alpz1, alpz2
real, dimension(nz,nx)                  :: window2D
real                                    :: xvalue, zvalue
real, parameter                                                         :: pi=3.141592654
optional                                :: alphax2, alphaz1, alphaz2

if(nz0 .gt. nz) then
        write(*,*) 'nz, nz0:', nz, nz0
        stop 'window2D: center of the window not inside matrix'
endif
if(nx0 .gt. nx) then
        write(*,*) 'nx, nx0:', nx, nx0
        stop 'window2D: center of the window not inside matrix'
endif


if(Nxlen .ge. 3* nx) then
        stop 'window2D: Nxlen too large'
endif
if(Nzlen .ge. 3*nz) then
        stop 'window2D: Nzlen too large'
endif

alpx1 = alphax1
if(present(alphaz1)) then
        alpz1 = alphaz1
else
        alpz1 = alpx1
endif
if(present(alphaz2)) then
        alpz2 = alphaz2
else
        alpz2 = alpz1
endif
if(present(alphax2)) then
        alpx2 = alphax2
else
        alpx2 = alpx1
endif

        window2D(:,:) = rzero_de;
        do ix = 1, nx
                do iz = 1, nz
                        ! x
                        ixx = ix - (nx0 - int(Nxlen * 0.5e0)); xvalue = rzero_de;
                        if((ixx .ge. 0) .and. (ixx .le. alpx1 * (Nxlen - 1) / 2)) then
                                if(alpx1 .eq. 0) then
                                        xvalue = 0.0;
                                else
                                        xvalue = (0.5e0 * (1.0 + cos(pi * ((2 * ixx / (alpx1 *(Nxlen - 1))) - 1.0))))
                                endif
                        elseif((ixx .ge. alpx1 * (Nxlen - 1) / 2) .and. (ixx .lt. (Nxlen -1) * &
                                                         (1.0 - alpx2 * 0.5e0))) then
                                xvalue = 1.0 ;
                        elseif((ixx .ge. (Nxlen - 1) *(1.0 - alpx2 * 0.5e0)) .and. (ixx .le.(Nxlen - 1))) then
                                if(alpx2 .eq. 0) then
                                        xvalue = 0.0;
                                else
                                        xvalue = (0.5e0 * (1.0 + cos(pi * ((2 * ixx / (alpx2 *(Nxlen - 1))) - &
                                                        (2.0e0 / alpx2) + 1.0))))       
                                endif
                        endif
                        ! z
                        izz = iz - (nz0 - int(Nzlen * 0.5e0)); zvalue = rzero_de;
                        if((izz .ge. 0) .and. (izz .le. alpz1 * (Nzlen - 1) / 2)) then
                                if(alpz1 .eq. 0) then
                                        zvalue = 0.0;
                                else
                                        zvalue = (0.5e0 * (1.0 + cos(pi * ((2 * izz / (alpz1 *(Nzlen - 1))) - 1.0))))
                                endif
                        elseif((izz .ge. alpz1 * (Nzlen - 1) / 2) .and. (izz .lt. (Nzlen -1) * &
                                                         (1.0 - alpz2 * 0.5e0))) then
                                       zvalue = 1.0 ;
                        elseif((izz .ge. (Nzlen -1) * (1.0 - alpz2 * 0.5e0)) .and. (izz .le.(Nzlen - 1))) then
                                if(alpz2 .eq. 0) then
                                        zvalue = 0.0;
                                else
                                        zvalue = (0.5e0 * (1.0 + cos(pi * ((2 * izz / (alpz2 *(Nzlen - 1))) - &
                                                (2.0e0 / alpz2) + 1.0))))
                                endif
                        endif
                        window2D(iz,ix) = window2D(iz,ix) + xvalue * zvalue
                enddo
        enddo
        window2D = (1.0e0 - window2D)

end function window2D

function window1D(nz, nz0, Nzlen, alphaz)

! this function returns a vector with a 1D window in it.
! the cosine tapered window is centered at [nz0] and has total width of Nzlen
! the output vector has dimensions [nz]
! alpha controls the steepness of the edges
! if Nlen is large such that the window exceeds the window, then window is cut
implicit none
integer                                 :: iz, izz
integer, intent(in)                     :: nz, nz0, Nzlen
real, intent(in)                        :: alphaz
real, dimension(nz)                     :: window1D
real                                    :: zvalue
real, parameter                                                         :: pi=3.141592654

if(nz0 .gt. nz) then
        write(*,*) 'nz, nz0:', nz, nz0
        stop 'window1D: center of the window not inside vector'
endif

if(Nzlen .ge. 3*nz) then
        stop 'window1D: Nzlen too large'
endif


        window1D(:) = rzero_de;
                do iz = 1, nz
                        ! z
                        izz = iz - (nz0 - int(Nzlen * 0.5e0)); zvalue = rzero_de;
                        if((izz .ge. 0) .and. (izz .le. alphaz * (Nzlen - 1) / 2)) then
                                if(alphaz .eq. 0) then
                                        zvalue = 0.0;
                                else
                                        zvalue = (0.5e0 * (1.0 + cos(pi * ((2 * izz / (alphaz *(Nzlen - 1))) - 1.0))))
                                endif
                        elseif((izz .ge. alphaz * (Nzlen - 1) / 2) .and. (izz .lt. (Nzlen -1) * &
                                                         (1.0 - alphaz * 0.5e0))) then
                                       zvalue = 1.0 ;
                        elseif((izz .ge. (Nzlen -1) * (1.0 - alphaz * 0.5e0)) .and. (izz .le.(Nzlen - 1))) then
                                if(alphaz .eq. 0) then
                                        zvalue= 0.0;
                                else
                                        zvalue = (0.5e0 * (1.0 + cos(pi * ((2 * izz / (alphaz *(Nzlen - 1))) - &
                                                (2.0e0 / alphaz) + 1.0))))
                                endif
                        endif
                        window1D(iz) = window1D(iz) + zvalue
                enddo
        window1D = (1.0e0 - window1D)

end function window1D



function insert_smooth_2Dbox(mat, x1, x2, z1, z2, alphax1, alphax2, alphaz1, alphaz2, value)
implicit none
integer, intent(in)                     :: x1, x2, z1, z2
real, intent(in), dimension(:,:)        :: mat
real, intent(in)                        :: alphax1, alphax2, alphaz1, alphaz2, value
real                                    :: alpx1, alpx2, alpz1, alpz2
real, dimension(size(mat, 1), size(mat, 2))                     :: insert_smooth_2Dbox
optional                                :: alphax2, alphaz1, alphaz2


! this function inserts a smooth rectangle in a 2D model
! uses the function window2D which can generate tapered rectangle 
! inside the box, mat will be equal to value

if(x1 .le. 0) stop 'insert_smooth_box: x1 out of range'
if(x2 .gt. size(mat,2)) stop 'insert_smooth_box: x2 out of range'
if(z2 .gt. size(mat,1)) stop 'insert_smooth_box: z2 out of range'
if(z1 .le. 0) stop 'insert_smooth_box: z1 out of range'
!

alpx1 = alphax1
if(present(alphaz1)) then
        alpz1 = alphaz1
else
        alpz1 = alpx1
endif
if(present(alphaz2)) then
        alpz2 = alphaz2
else
        alpz2 = alpz1
endif
if(present(alphax2)) then
        alpx2 = alphax2
else
        alpx2 = alpx1
endif


insert_smooth_2Dbox = rzero_de;
                insert_smooth_2Dbox = mat + (value - mat) * &
                        (1.0-window2D(nz = size(mat,1), nx= size(mat,2), &
                                nx0= int(0.5e0*(x1 + x2)) +1, &
                               nz0=  int(0.5e0*(z1 + z2)) + 1, &
                               Nzlen = z2 - z1 + 5 , alphaz1=alpz1, alphaz2=alpz2, &
                               Nxlen = x2 - x1 + 5 , alphax1=alpx1, alphax2=alpx2))

end function insert_smooth_2Dbox

subroutine taper_trc( &
!  Apply cosine tapers to first frctnfront*nz points of beginning
!  and last frctnback*nz points of time series array z.
                frctnfront, & ! fraction of points to be tapered in front 
                        ! input zero do not apply taper
                frctnback, &  ! fraction of points to be tapered at end
                        ! input zero do not apply taper
                z, & !the vector
                nz) ! size of vector
implicit none
real, intent(in)                   :: frctnfront,frctnback
real                               :: z(:), pi,  f, sf
integer                            :: nz, i, lz

call check_dimension('taper_trc: z',z, nz)
pi = 4.0*atan(1.0)
lz = int(nz*frctnfront,i_de)
if (lz >= 1) then
        sf = pi/lz
        do i=1,lz
                f = 0.5*(1.0-cos(sf*(i-1)))
                z(i) = z(i)*f
        end do
end if
lz = int(nz*frctnback,i_de)
if (lz >= 1) then
        sf = pi/lz
        do i=nz,nz-lz+1,-1
                f = 0.5*(1.0-cos(sf*(nz-i)))
                z(i) = z(i)*f
        end do
end if
return
end subroutine taper_trc


subroutine windows(window,n,iwindow) 
!   Routine WINDOW: To Obtain Window Function. 
!   Input parameters: 
!    n      : the length of window data. 
!    iwindow: window type desired. 
!    if     : iwindow=1: rectangular window ,  =2: triangular window , 
!                    =3: cosin window ,        =4: Hanning window , 
!                    =5: Hamming window ,      =6: Blackman window , 
!                    =7: Papoulis window . 
!   Output parameters: 
!    window     : N dimensioned real array.the result is in w(0) to w(n-1). 
!    ierror:IF IERROR=0: no error,    =1: Iwindow out of range. 
!                                       in chapter 8 
! ---------------------------------------------------------------------- 
IMPLICIT NONE
        INTEGER, INTENT(IN) :: n,iwindow
INTEGER :: i
        REAL, DIMENSION(0:n-1) :: window
        REAL :: pn,a,b,c, pi
        pi=3.1416
        pn=2.*pi/float(n)
SELECT CASE(iwindow)
CASE(1)
        
     DO i=0,n-1 
         window(i)=1. 
     END DO 
CASE (2)
     DO i=0,n-1 
         window(i)=1.-abs(1.-2.*i/float(n)) 
     END DO
CASE(3) 
    DO i=0,n-1 
        window(i)=sin(pn*i/2.) 
    END DO
     write(*,*)'test' 
CASE(4)
      DO i=0,n-1 
      window(i)=0.5*(1.0-cos(pn*i)) 
      END DO 
CASE(5)
      DO i=0,n-1 
         window(i)=0.54-0.46*cos(pn*i) 
      END DO 
CASE(6)
      DO i=0,n-1 
        window(i)=0.42-0.5*cos(pn*i)+0.08*cos(2.*pn*i) 
       END DO
CASE(7) 
      DO i=0,n-1 
           a=abs(sin(pn*i))/pi 
           b=1.-2.*(abs(i-n/2.))/float(n) 
           c=cos(pn*i) 
          window(i)=a-b*c 
           
      END DO

CASE DEFAULT
WRITE(*,*)'INVALID window type :: window' 
STOP
END SELECT
end subroutine windows

subroutine taper(nz,gather,gather1d,sample,&
        percent1, & ! percentage taper in the beginning (default=10)
        percent2, & ! percentage taper towards the end (default=10)
        dimflag)
! this function returns windows and also tapers the ends if the input
! Author : Pawan Bharadwaj
!          p.b.pisupati@tudelft.nl
! Input: 
!       nz              -- number of samples along 1st dimension
!       gather          -- 1st dimension is tapered of this gather if dimflag='1' (default)
!                       -- 2nd dimension is tapered if dimflag='2'
!       sample          -- returns a sample window with both ends tapered of size nz

!       flag            -- 1  only tapered in the begginning
!                       -- 2  only tapered at end
!                       -- dont use flag to taper both sides

       implicit none
        real,intent(inout),dimension(:,:), optional             :: gather
        real,intent(inout),dimension(:), optional               :: gather1d
        real,intent(in),optional                                :: percent1, percent2
        real                                                    :: perc1, perc2
        real,intent(out),dimension(nz),optional                 :: sample
        integer                                                 :: n,inn, dflag
        integer,intent(in)                                      :: nz, dimflag
        optional                                                :: dimflag
if(present(percent1)) then
        perc1=percent1
else 
        perc1=10.0
end if
if(present(percent2)) then
        perc2=percent2
else 
        perc2=10.0
end if

if(present(dimflag)) then
        dflag = dimflag
else
        dflag = 1 ! default
endif
if((dflag .ne. 1) .and. (dflag .ne. 2)) stop 'taper: invalid dflag'
if (present(gather)) then
        if(size(gather, 1) .ne. nz) stop 'taper: gather first dimension size .ne. nz'
        if(dflag .eq. 1) then
                n=size(gather,2);
                do inn=1,n
                        call taper_trc(perc1/1e2,perc2/1e2,gather(:,inn),nz)
                end do
        elseif(dflag .eq. 2) then
                n=size(gather,1);
                do inn=1,n
                        call taper_trc(perc1/1e2,perc2/1e2,gather(inn,:),size(gather, 2))
                end do
        endif
end if
if (present(gather1d)) then
        if(size(gather1d, 1) .ne. nz) stop 'taper: gather1d first dimension size .ne. nz'
        call taper_trc(perc1/1e2,perc2/1e2,gather1d(:),nz)
end if



if (present(sample)) then
        sample(:)=1.0
        call taper_trc(perc1/1e2,perc2/1e2,sample,nz)
end if
END subroutine taper

subroutine lr_window(window,flag,factor,nz)
! this subroutine returns either a right or left window depending on factor
! Input
!               flag            -- 1 (from 0 to 1)
!               factor          -- controling steepness
!               nz              -- number of samples in window
! Output
!               window          -- output window.. either left (0 to 1) ot right (1 to 0)
implicit none
                real, intent(out),dimension(nz)         :: window
                real, intent(in)                        :: factor
                integer,intent(in)                      :: nz,flag
                real,allocatable,dimension(:)           :: sample
                integer                                 :: nzt

window(1:nz)=0.0
nzt=int(nz/factor);
allocate(sample(nzt))
sample(:)=1.0;
call taper(nzt,sample=sample,percent1=factor*1e2,percent2=factor*1e2)
if(flag.eq.1) then
        window=sample(1:nz);
endif
if(flag.eq.2) then
        window=sample(nzt-nz+1:nzt)
endif

end subroutine lr_window


end module tapering
