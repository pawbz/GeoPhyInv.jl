module fft
use error_messages! USEINDEPFILE
use interpolation! USEINDEPFILE
use precision_mod! USEINDEPFILE
use, intrinsic :: iso_c_binding
use omp_lib

implicit none

contains

subroutine lenpow2(n,    & ! any integer 
                   ni    & ! next nearest power of 2 to 2 * n
                  ) bind(c, name="lenpow2")
implicit none
integer, intent(in)                             :: n
integer, intent(out)                            :: ni
integer                                         :: i
do i = 1, 50
        if(2**i .ge. 2*n) then
                ni = 2**i
                exit
        end if
end do
end subroutine lenpow2

function findfreq(x,& ! time series
	dt, & ! time sampling interval
	str, &  ! str  = PEAK returns peak frequency
		! str  = MIN returns minumum frequency
		! str  = MAX returns maximum frequency
	threshold & ! threshold to select maximum and minumum frequency
		    ! compared to the normalized power spectrum 
)
implicit none
real, dimension(:), intent(in)          :: x
real, intent(in)                        :: dt, threshold
character(len=*), intent(in)		:: str
optional				:: threshold

real					:: thresholdin
integer					:: npow2, nlag, ind,n
real                                    :: findfreq
real					:: df
complex, dimension(:),allocatable       :: cx
real, dimension(:), allocatable         :: ax, fvec
logical, dimension(:), allocatable      :: mask

if(present(threshold)) then; thresholdin = threshold; else; thresholdin = 1e-3; endif

n=size(x); 
call lenpow2(2*n,npow2)
nlag = (npow2/2);
df = 1.0 / npow2 / dt 
allocate(cx(npow2), ax(npow2), fvec(npow2))
cx = cmplx(rzero_de, rzero_de);
ax = rzero_de;

cx(1:n) = cmplx(x, rzero_de)

call cfft(cx, npow2, 1)
ax = real(cx * conjg(cx))
fvec = npow2_vector(d=df, npow2=npow2)

! normalizing
if(maxval(ax) .eq. rzero_de) then
        findfreq = 0.0
        return
else
	ax = ax / maxval(ax)
endif

if(str .eq. 'PEAK') then
	ind = maxloc(ax(1:nlag+1),dim=1)
	findfreq = fvec(ind) 
elseif(str .eq. 'MIN') then
	allocate(mask(1:nlag+1))
	mask = ax(1:nlag+1) .gt. thresholdin 
	findfreq = minval(fvec(1:nlag+1), dim=1, mask=mask)
	deallocate(mask)
elseif(str .eq. 'MAX') then
	allocate(mask(1:nlag+1))
	mask = ax(1:nlag+1) .gt. thresholdin 
	findfreq = maxval(fvec(1:nlag+1), dim=1, mask=mask)
	deallocate(mask)
else
        findfreq = 0.0
endif
end function findfreq


pure subroutine cfft( dat, n, iforw )
! df = 1/ n / dt
! complex fft
implicit none
integer                 :: i1, i2a, i2b, i3, i3rev, ip1, ip2,  isign
integer, intent(in)     :: iforw, n
real                    :: theta, sinth
complex, intent(inout)  :: dat(:)
complex                 :: temp, w, wstp
isign = -iforw
i3rev = 1
do i3 = 1, n
        if ( i3 < i3rev ) then   ! switch values
                temp          = dat( i3 )
                dat( i3 )    = dat( i3rev )
                dat( i3rev ) = temp
        endif
        ! following loop is just to compute i3rev
        ip1 = n / 2
        do while (  i3rev > ip1 )
                if ( ip1 <= 1 ) exit
                i3rev = i3rev - ip1
                ip1   = ip1 / 2
        end do
        i3rev = i3rev + ip1
end do
ip1 = 1
do while  ( ip1 < n )
        ip2   = ip1 * 2
        theta = 6.283185307 / float( isign * ip2 )
        sinth = sin( theta / 2.)
        wstp  = cmplx( -2. * sinth * sinth, sin( theta ) )
        w     = cone_de
        do i1 = 1, ip1
                do i3 = i1, n, ip2
                        i2a = i3
                        i2b = i2a + ip1
                        temp      = w * dat( i2b )
                        dat( i2b ) = dat( i2a ) - temp
                        dat( i2a ) = dat( i2a ) + temp
                end do
                w = w * wstp + w
        end do
        ip1 = ip2
end do
return
end subroutine cfft

subroutine cfft2( dat, n1, n2, iforw )
! apply a 2 dimensional forward and inverse 
! fourier transform using cfft
! df = 1/ n / dt
! complex fft
implicit none

integer, intent(in)             :: n1, n2, iforw
complex, intent(inout)          :: dat(:,:)
integer                         :: i1, i2

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i1, i2) NUM_THREADS(OMP_GET_MAX_THREADS())&
!$OMP SHARED(n2, n1, dat, iforw)
!$OMP DO
do i1 = 1, n1
        call cfft(dat(i1, :), n2, iforw)
enddo 
!$OMP END DO

!$OMP DO
do i2 = 1, n2
        call cfft(dat(:, i2), n1, iforw)
enddo 
!$OMP END DO
!$OMP END PARALLEL

end subroutine cfft2

pure subroutine nlag_npow2_pad_truncate( &
                nplags, & ! number of positive lags
                nnlags, & ! number of negative lags
                npow2, &! number of samples in xpow2
                x, & ! real signal with dimension nplag + nnlag + 1
                     ! first has decreasing negative lags
                     ! then has zero lag at nnlags + 1
                     ! then has increasing positive nplags lags 
                     ! signal contains only positive lags and zero lag if nnlag
                     ! =0 and vice versa
                xpow2, &! output lenpow2 real vector with dimension npow2 
                xpow2_complex, &! output lenpow2 vector with dimension npow2
                flag &
                        ! flag = 1 means xpow2 is returned using x (zero pad or truncation accordingly)
                        ! flag = -1 means x is returned using xpow2 (zero pad or truncation accordingly)
                )
implicit none
real, intent(inout)                     :: x(nplags + nnlags + 1)
real, intent(inout)                     :: xpow2(npow2)
complex, intent(inout)                  :: xpow2_complex(npow2)
integer, intent(in)                     :: nplags, nnlags, npow2
integer, intent(in)                     :: flag
optional                                :: xpow2, xpow2_complex
integer                                 :: i

!if(iand(npow2,npow2-1) .ne. 0) call abort_msg('n_npow2_pad_truncate: npow2 must be power of 2')
!if(n .gt. npow2/2 + 1) call abort_msg('n_npow2_pad_truncate: n .gt. npow2/2 + 1')

if(flag .eq. 1) then
        if(present(xpow2)) then
                xpow2 = rzero_de;
                ! zero lag        
                xpow2(1) = x(nnlags + 1);
                if(nplags .gt. 0) then
                        ! positive lags
                        xpow2(2 : nplags + 1) = x(nnlags + 2 : nnlags +1 + nplags);
                endif
                if(nnlags .gt. 0) then
                        ! negative lags
                        forall(i=1:nnlags)
                                xpow2(npow2 - i + 1) = x(nnlags + 1 - i)
                        endforall
                endif
        endif
        if(present(xpow2_complex)) then
                xpow2_complex = czero_de
                ! zero lag        
                xpow2_complex(1) = cmplx(x(nnlags + 1), rzero_de);
                if(nplags .gt. 0) then
                        ! positive lags
                        xpow2_complex(2 : nplags + 1) = cmplx(x(nnlags + 2 : nnlags +1 + nplags), rzero_de);
                endif
                if(nnlags .gt. 0) then
                        ! negative lags
                        forall(i=1:nnlags)
                                xpow2_complex(npow2 - i + 1) = cmplx(x(nnlags + 1 - i), rzero_de);
                        endforall
                endif
        endif
elseif(flag .eq. -1) then
        x = rzero_de;
        if(present(xpow2)) then
                x(nnlags + 1) = xpow2(1);
                if(nplags .gt. 0) then
                        forall(i=1:nplags)
                                x(nnlags + 1 + i) = xpow2(1 + i)
                        endforall
                endif
                if(nnlags .gt. 0) then
                        forall(i=1:nnlags)
                                x(nnlags + 1 - i) = xpow2(npow2 - i + 1)
                        endforall
                endif
        elseif(present(xpow2_complex)) then
                x(nnlags + 1) = real(xpow2_complex(1));
                if(nplags .gt. 0) then
                        forall(i=1:nplags)
                                x(nnlags + 1 + i) = real(xpow2_complex(1 + i))
                        endforall
                endif
                if(nnlags .gt. 0) then
                        forall(i=1:nnlags)
                                x(nnlags + 1 - i) = real(xpow2_complex(npow2 - i + 1))
                        endforall
                endif
        endif
endif
end subroutine nlag_npow2_pad_truncate

function npow2_vector( &
        ! this subroutine returns a time or frequency vector for npow2 
                d, &     ! change 
                npow2 &  ! dimension
                ) result(vec) !vector of size npow2
implicit none
real                                    :: vec(npow2)
real, intent(in)                        :: d
integer, intent(in)                     :: npow2
integer                                 :: i

vec = rzero_de;

! zero lag
vec(1) = rzero_de;

! positive lags
do i = 1, npow2/2
        vec (1 + i) = d * i
enddo

! negative lags -- one less than positive lags
do i = 1, npow2/2-1 
        vec(npow2 -i +1) = -d * i
enddo

end function npow2_vector


!subroutine zerophase_bp_mapping(&
!! band pass mapping of the spectrum of a zerophase signal
!! the spectrum is real
!!  ^                                                           
!!  |                                                           
!!  |                                                           
!!  |       Spectrum before bandpass mapping  - xb
!!  |                                                           
!!  |        +-----------------+                                
!!  |        |                 |                                
!!  |        |                 |                                
!!  |        |                 |                                
!!  +--------+-----------------+--------------------------->    
!!          f1                 f2                               
!!  ^                                                           
!!  |                                                           
!!  |                                                           
!!  |       Spectrum after bandpass mapping - xf                    
!!  |                                                           
!!  +-----------------------------------------+                 
!!  |                                         |                 
!!  |                                         |                 
!!  |                                         |                 
!!  +-----------------------------------------+------------->   
!!  f1                                        f2                
!!
!                npow2, & ! number of samples in the spectrum 
!                df, & ! frequency sampling
!                f1, & ! first frequency
!                f2, & ! second frequency
!                xspec, & ! amplitude spectrum [npow2] [real]
!                flag & ! flag .eq. 1 means output xf using xb
!                       ! flag .eq. -1 means output xb using xf
!                )
!
!implicit none
!integer, intent(in)                             :: npow2
!real, intent(in)                                :: df, f1, f2
!real, intent(inout)                             :: xspec(:)
!integer, intent(in)                             :: flag
!
!type(mesh1D)                                    :: fb, ff
!integer                                         :: nlag, ifreq, ifreq1, ifreq2
!real, allocatable                               :: xbtemp(:), xftemp(:)
!real                                            :: fn, x, x1, x2, y1, y2
!
!        
!call check_dimension('zerophase_bp_mapping: size of xspec', xspec, npow2)
!! check if npow2 is actually a power of two
!if(iand(npow2,npow2-1) .ne. 0) call abort_msg('zerophase_bp_mapping: npow2 must be power of 2')
!
!! nlag
!nlag = (npow2/2);
!
!fn = 0.5 *  npow2 * df; ! nyquist
!
!! temp vectors for only positive lags
!allocate(xbtemp(1:nlag+1) )
!allocate(xftemp(1:nlag+1) )
!xbtemp = rzero_de; xftemp = rzero_de; 
!
!call create_mesh1D(mesh=fb, Tmin=0.0, Tmax=fn, nt=nlag+1)
!call create_mesh1D(mesh=ff, Tmin=f1, Tmax=f2, nt=nlag+1)
!
!call is_nan('zerophase_bp_mapping: input xspec', xspec)
!
!if(flag .eq. -1) then
!        ! input, output to zeros
!        xftemp = xspec(1:nlag + 1)
!        xspec = rzero_de;
!        xbtemp  = rzero_de;
!        call spray(y_in=xftemp, y_out=xbtemp, mesh_in=ff, mesh_out=fb, flag = 1)
!        xspec(1 : nlag + 1) = xbtemp;
!        xspec(npow2: nlag +2 : -1) = xbtemp(2:nlag);
!                
!elseif(flag .eq. 1) then
!        ! input, output to zeros
!        xbtemp = xspec(1:nlag + 1)
!        xspec = rzero_de;
!        call resample(y_in=xbtemp, y_out=xftemp, mesh_in=fb, mesh_out=ff, flag = 1)
!        xspec(1 : nlag + 1) = xftemp;
!        xspec(npow2: nlag +2 : -1) = xftemp(2:nlag);
!endif
!
!call is_nan('zerophase_bp_mapping: output xspec', xspec)
!
!end subroutine zerophase_bp_mapping

end module fft
