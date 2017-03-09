module filters
use string_routines! USEINDEPFILE
use io! USEINDEPFILE
use fft! USEINDEPFILE
use precision_mod! USEINDEPFILE
use tapering! USEINDEPFILE

! Author : Pawan Bharadwaj
!          p.b.pisupati@tudelft.nl
!          March 2015

implicit none
contains

subroutine ZPspec_npow2( &
! ZP stands for zerophase
! this subroutine returns a butterworth zerophase filter (in frequency domain)
! of dimension npow2 
                        ZPspec,& ! filter in freq domain (complex but imag part is zero)
                        order, & ! order of the butterworth filter
                        npow2, & ! dimension of the filter in the frequency domain
                        df, &    ! sampling of frequency 
                        fmin,&   ! lowcut
                        fmax &   ! highcut 
                ! if fmin .eq. 0 then lowpass filter at fmax
                ! if fmax .eq. 0 then highpass filter at fmin
                ! otherwise; bandpass filter fmin--fmax
                                )
! Author : Pawan Bharadwaj
!          p.b.pisupati@tudelft.nl

implicit none
real, intent(in)                        :: fmin, fmax, df
integer, intent(in)                     :: order
complex, intent(out)                    :: ZPspec(:)
integer, intent(in)                     :: npow2
real, allocatable                       :: fvec(:), x1(:)

if(fmin .lt. 0.0) call abort_msg('ZPspec_npow2: fmin cannot be .lt. zero')
if(fmax .lt. 0.0) call abort_msg('ZPspec_npow2: fmax cannot be .lt. zero')
if(fmax .ne. 0.0) then
        if(fmax .le. fmin) call abort_msg('ZPspec_npow2: fmax .le. fmin')
endif
call check_dimension('ZPspec_npow2: size of x',ZPspec, npow2)

! create a freq vector
allocate(fvec(npow2)); fvec = rzero_de;
fvec =  npow2_vector(d= df, npow2 = npow2)

if(fmin .eq. 0.0) then
        ! lowpass filter
        ZPspec  = cmplx((1.0 + (fvec/fmax)**(2.0*real(order,r_de)))**(-1.0), 0.0)
elseif(fmax .eq. 0.0) then
        ! highpass filter
        ZPspec  = cmplx((1.0 + (fmin/fvec)**(2.0*real(order,r_de)))**(-1.0), 0.0)
elseif((fmin .ne. 0.0) .and. (fmax .ne. 0.0)) then
        ! bandpass filter
        allocate(x1(npow2)); x1 = 0.0;
        x1 = sqrt(fmin*fmax) / (fmax-fmin) * (fvec/sqrt(fmin*fmax) + sqrt(fmax*fmin)/fvec);
        ZPspec  = cmplx((1.0 + (x1)**(2.0*real(order,r_de)))**(-1.0), 0.0)
        deallocate(x1)
endif
! DC component is always forced zero
ZPspec(1) = cmplx(0.0, 0.0);
! deallocate frequency vector
deallocate(fvec)

end subroutine ZPspec_npow2

subroutine MPspec_npow2( &
! first a zerophase butterworth filter is created
! then it is converted to minimum phase filter in the cepstral domain
! more info: http://www.katjaas.nl/minimumphase/minimumphase.html
! positive quefrencies correspond to minimum phase component of the signal
! negetive quefrencies correspond to maximum phase component of the signal
! signal = minimum phase [convolved with] maximum phase
                        MPspec,& ! minimum phase filter in freq domain (note that it is complex)
                        order, & ! order of the butterworth filter
                        npow2, & ! dimension of the filter in the frequency domain
                        df, &    ! sampling of frequency 
                        fmin,&   ! lowcut
                        fmax &   ! highcut 
                        )
                ! if fmin .eq. 0 then lowpass filter at fmax
                ! if fmax .eq. 0 then highpass filter at fmin
                ! otherwise; bandpass filter fmin--fmax
! Author : Pawan Bharadwaj
!          p.b.pisupati@tudelft.nl

implicit none
real, intent(in)                        :: fmin, fmax, df
integer, intent(in)                     :: order
complex, intent(out)                    :: MPspec(:)
complex, allocatable                    :: X(:)
integer, intent(in)                     :: npow2

allocate(X(npow2)); X = czero_de;
call ZPspec_npow2(ZPspec = X, fmin=fmin, fmax=fmax, df=df, npow2=npow2, order=order)
call ZPspec_to_MPspec_npow2(ZPspec = X, MPspec = MPspec)
deallocate(X)

end subroutine MPspec_npow2

subroutine ZPspec_to_MPspec_npow2( &
! more info: http://www.katjaas.nl/minimumphase/minimumphase.html
                ZPspec,&! frequency spectrum of a zerophase signal 
                        ! (complex, but only abs is used discarding phase)
                MPspec &! Minimumphase spectrum with both phase and amplitude
                        ! (complex)
                        )
! Author : Pawan Bharadwaj
!          p.b.pisupati@tudelft.nl

complex, intent(in)                     :: ZPspec(:)
complex, intent(out)                    :: MPspec(:)
complex, allocatable                    :: X(:)
integer                                 :: npow2
real                                    :: damp 

! dimension check
call check_dimension('ZPspec_to_MPspec', MPspec, size(ZPspec,1))
npow2 = size(ZPspec, 1);
if(iand(npow2,npow2-1) .ne. 0) call abort_msg('ZPspec_to_MPspec: dimension must be power of 2')

allocate(X(npow2)); X = cmplx(0.0, 0.0);

! prevent log(0)
damp = 1e-20 * maxval(abs(ZPspec));
! logarithm 
X = log(cmplx(abs(ZPspec) + damp, 0.0)) 
! to cepstral domain - IFFT
call cfft(X, npow2, -1)
! only real part output
X = cmplx(real(X), 0.0);
! scaling
X = X / cmplx(real(npow2), 0.0)

! positive cepstrum x 2
X(2 : npow2/2 + 1) = cmplx(2.0, 0.0) * X(2 : npow2/2 + 1)
! remove negative quefrencies
X(npow2/2 + 2 : npow2) = cmplx(0.0, 0.0) 

! FFT
call cfft(X, npow2, 1)

! exponential
MPspec = exp(X)

MPspec = MPspec/cmplx(maxval(abs(MPspec)), 0.0)
deallocate(X)

end subroutine ZPspec_to_MPspec_npow2


function apply_butterfilter(gatherin,& ! input gather; first dimension is time
                         dt, & ! time sampling
                         fmin, & ! minimum frequency
                         fmax, & ! maximum frequency
                         order, & ! order of butterworth filter
                         n1, &  ! time dimension
                         n2, &  ! second dimension
                         n3, &  ! third dimension
                         filter_type & ! filter_type
                                ! 0 means zerophase filter
                                ! 1 means minimum phase filter
                                ) result(gatherout)
                                        ! result after appling a ZP butterworth
                                        ! filter

implicit none
integer, intent(in)                     :: n1, n2, n3, order, filter_type
real, intent(in)                        :: gatherin(:), dt, fmin, fmax
real                                    :: gatherout(n1*n2*n3)
complex, allocatable                    :: X(:)

integer                                 :: npow2

call check_dimension('apply_butterfilter: gatherin', gatherin, n1*n2*n3)
gatherout(:) = 0.0;

! filter length
if(filter_type .eq. 0) then
        call lenpow2(n=2*n1 + 1, ni=npow2)
elseif(filter_type .eq. 1) then
        ! large length for minimum phase filter
        ! to prevent cepstral aliasing
        call lenpow2(n=10*n1, ni=npow2)
else
        call abort_msg('apply_butterfilter: unknown filter_type')
endif
       
! design a filter - X
allocate(X(npow2));
if(filter_type .eq. 0) then
        call ZPspec_npow2(ZPspec = X, fmin=fmin, fmax=fmax, df=1.0 / npow2 / dt, npow2=npow2, order=order)
elseif(filter_type .eq. 1) then
        call MPspec_npow2(MPspec = X, fmin=fmin, fmax=fmax, df=1.0 / npow2 / dt, npow2=npow2, order=order)
endif
! apply filter
gatherout = apply_filter(gatherin = gatherin, X = X, npow2=npow2, n1=n1, n2=n2, n3=n3)
! deallocate X
deallocate(X)
end function apply_butterfilter

function apply_filter(gatherin,& ! input gather; first dimension is time
                         X, & ! filter in freq domain [npow2]
                         npow2, & ! dimension after zero padding
                         n1, &  ! time dimension before zero padding
                         n2, &  ! second dimension
                         n3  &  ! third dimension
                                ) result(gatherout)
                                        ! result after appling X
implicit none
integer, intent(in)                     :: n1, n2, n3, npow2
real, intent(in)                        :: gatherin(:)
complex, intent(in)                     :: X(:)
real                                    :: gatherout(n1*n2*n3)

integer                                 :: i3, i2, i1, a1, a2
complex, allocatable                    :: GIN(:)
real, allocatable                       :: temp1d(:)

call check_dimension('apply_filter: gatherin', gatherin, n1*n2*n3)
call check_dimension('apply_filter: X', X, npow2)
gatherout(:) = 0.0;

do i3 = 1, n3
        ! allocations
        allocate(GIN(npow2)); 
        allocate(temp1d(n1)); 

        do i2 = 1, n2
                a1 = 1 + (i2-1) * n1 + (i3-1) * n1 * n2;   a2 =   (i2) * n1 + (i3-1) * n1 * n2 ;

                GIN = cmplx(rzero_de, rzero_de); 
                ! input time series is only causal, assume 
                temp1d = (gatherin(a1:a2));
                ! zero pad
                call nlag_npow2_pad_truncate(x = temp1d, &
                                 xpow2_complex = GIN, nplags = n1-1, nnlags=0, npow2 = npow2, flag = 1)
                ! FFT
                call cfft(GIN, npow2, 1)
                ! apply the filter
                GIN = GIN * X
                ! IFFT
                call cfft(GIN, npow2, -1)
                ! scaling
                GIN = GIN / cmplx(real(npow2), 0.0)
                ! truncate
                call nlag_npow2_pad_truncate(x = gatherout(a1:a2), &
                                 xpow2_complex = GIN, nplags = n1-1,nnlags=0, npow2 = npow2, flag = -1)
        enddo
        ! deallocate rest
        deallocate(GIN, temp1d)
enddo

end function apply_filter

subroutine bandfilter2d(gather_out,gather_in,f1,f2,dt,nt,nx)
! (not used -- use apply_butterfilter )

! this subroutine causal hicut filter of traces in gather_in
! Input:
!       gather_in                       :: input time series
!       dt                              :: time sampling interval
!       fh                              :: frequencies above fh are removed
! Output:
!       gather_out                      :: filtered trace

        implicit none
        integer,intent(in)   :: nt,nx
        integer              :: ix,ntpad,ntnew
        real, dimension(:,:),allocatable &
                             :: gin,gout
        real, intent(in)     :: gather_in(nt,nx),f1,f2,dt
        real, intent(out)    :: gather_out(nt,nx)
        real,allocatable     :: envel(:)
if(size(gather_out,1).ne.nt .or. size(gather_out,2).ne.nx) then
        write(*,*)'do_bw_filter: output gather size different from input gather'
        STOP
end if

ntpad=0
ntnew=nt+ntpad*2
allocate(gin(1:ntnew,nx),gout(1:ntnew,nx))
allocate(envel(ntnew))
gin(:,:)=0.0;
gout(:,:)=0.0;
gin(ntpad+1:ntpad+nt,:)=gather_in;
gather_out(:,:)=0.0

do ix=1,nx
        gout(:,ix)=gin(:,ix)
        call band(gout(:,ix),nd=ntnew,f1=f1,f2=f2,delt=dt,nroll=2,icaus=1)
enddo
gather_out(:,:)=gout(ntpad+1:ntpad+nt,:);
end subroutine bandfilter2d

subroutine bandfilter1d(trc_out,trc_in,f1,f2,dt,nt)
! (not used -- use apply_butterfilter )

! this subroutine applies a CAUSAL hicut filter to  trc_in.

! Input:
!       trc_in                  :: input time series
!       dt                      :: time sampling interval
!       fh                      :: frequencies above fh are removed
! Output:
!       trc_out                 :: filtered trace

        implicit none
        integer,intent(in)   :: nt
        real, intent(in)     :: trc_in(nt),f1,f2,dt
        real, intent(inout)  :: trc_out(:)
if(size(trc_out,1).ne.nt .or.size(trc_in,1).ne.nt) then
        write(*,*)'do_bw_filter: output gather size different from input gather'
        STOP
end if
trc_out=trc_in;
call band(trc_out,nd=nt,f1=f1,f2=f2,delt=dt,nroll=2,icaus=1)
end subroutine bandfilter1d


subroutine band(s,nd,f1,f2,delt,nroll,icaus)
! (not used -- use apply_butterfilter )
!  Butterworth bandpass filter order 2*nroll (nroll<=8) (see Kanasewich, 
!    Time Sequence Analysis in Geophysics, Third Edition, 
!    University of Alberta Press, 1981)
!  written by W.B. Joyner 01/07/97
!  causal if icaus.eq.1 - zero phase shift otherwise
!  s(j) input = the time series to be filtered - output = the 
!    filtered series
!  dimension of s(j) must be as least as large as the larger of 
!    the following:
!    nd+3.0*float(nroll)/(f1*delt)
!    nd+6.0*float(nroll)/((f2-f1)*delt)
!  nd = the number of points in the time series
!  f1, f2 = the cutoff frequencies
!  delt = the timestep
! Dates: xx/xx/xx - Written by Bill Joyner
!        09/12/00 - Changed "n" to "nroll" to eliminate confusion with
!                   Kanesewich, who uses "n" as the order (=2*nroll), and
!                   cleaned up appearance of code by using spaces, indents, etc.
!        09/12/00 - double precision statements added so that the output
!                   series has sufficient precision for double integration.
!        11/08/00 - Increased dimension of s from 50000 to 100000
!        02/14/01 - Modified/corrected real numbers function calls to 
!                  double precision - cds
!        02/22/01 - Removed dimension of s (it is up to the user to specify
!                   it properly)
        implicit none
        integer, intent(in)     :: nd
        real, intent(inout)     :: s(nd)
        real, intent(in)        :: f1,f2, delt
        integer, intent(in)     :: icaus, nroll
        integer                 :: k,i,j,np1,np2,npad
        
      real(r_dp) fact(16),b1(16),b2(16)
      real(r_dp) pi,w1,xp,yp,x1,x2,y1,y2, w2
      real(r_dp) pre, pim, argre, argim, rho, theta, sjre, sjim
      real(r_dp) bj, cj, con
      if(f1.eq.0..or.f1.eq.f2) return
      pi=4.0d0*datan(1.0d0)
      w1=2.0d0*pi*f1
      w1=2.0d0*dtan(w1*delt/2.0d0)/real(delt,r_dp)
      w2=2.0d0*pi*f2
      w2=2.0d0*dtan(w2*delt/2.0d0)/real(delt,r_dp)
      do k=1,nroll
        pre=-dsin(pi*dfloat(2*k-1)/dfloat(4*nroll))
        pim=dcos(pi*dfloat(2*k-1)/dfloat(4*nroll))
        argre=(pre**2-pim**2)*(w2-w1)**2/4.0d0-w1*w2
        argim=2.0d0*pre*pim*(w2-w1)**2/4.0d0
        rho=(argre**2+argim**2)**(1.0d0/4.0d0)
        theta=pi+datan2(argim,argre)/2.0d0
        do i=1,2
          sjre=pre*(w2-w1)/2.0d0+(-1)**i*rho*(-dsin(theta-pi/2.0d0))
          sjim=pim*(w2-w1)/2.0d0+(-1)**i*rho*dcos(theta-pi/2.0d0)
          bj=-2.0d0*sjre
          cj=sjre**2+sjim**2
          con=1.0d0/(2.0d0/delt+bj+cj*delt/2.0d0)
          fact(2*k+i-2)=(w2-w1)*con
          b1(2*k+i-2)=(cj*delt-4.0d0/delt)*con
          b2(2*k+i-2)=(2.0d0/delt-bj+cj*delt/2.0d0)*con
        end do
      end do
      np2=nd
      if(icaus.ne.1) then
        npad=3*int(float(nroll)/(f1*delt),i_de)
        if( npad .lt. int(6.0*float(nroll)/((f2-f1)*delt),i_de) ) then
          npad=int(6.0*float(nroll)/((f2-f1)*delt),i_de)
        end if
        np1=nd+1
        np2=nd+npad
        do j=np1,np2
          s(j)=0.0
        end do
      end if
      do k=1,2*nroll
        x1=0.0d0
        x2=0.0d0
        y1=0.0d0
        y2=0.0d0
        do j=1,np2
          xp=real(s(j),r_dp)
          yp=fact(k)*(xp-x2)-b1(k)*y1-b2(k)*y2
          s(j)=real(yp,r_de)
          y2=y1
          y1=yp
          x2=x1
          x1=xp
        end do
      end do
      if(icaus.ne.1) then
        do k=1,2*nroll
          x1=0.0d0
          x2=0.0d0
          y1=0.0d0
          y2=0.0d0
          do j=1,np2
            xp=real(s(np2-j+1),r_dp)
            yp=fact(k)*(xp-x2)-b1(k)*y1-b2(k)*y2
            s(np2-j+1)=real(yp,r_de)
            y2=y1
            y1=yp
            x2=x1
            x1=xp
          end do
        end do
      end if
      return
end subroutine band

end module filters
