module source_wavelets
use filters! USEINDEPFILE
use error_messages! USEINDEPFILE

implicit none
private
public                                                          :: ricker, ormsby, MPfiltered_spike, fricker
real,parameter,private                                          :: zero = 0.d0

contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! wavelet related fucntions
! go_wavelet
subroutine fricker(wfreqs,fpeak,freqs)
! returns a zero phase ricker wavelet in frequency domain for frequencies freqs
! 
implicit none 
real, intent(in)                :: freqs(:), fpeak
complex, intent(out)            :: wfreqs(:)
integer                         :: ifreq, nfreq
real                            :: fmax, pi
complex                         :: a,b


pi=3.14159265358979e0

if(size(freqs) .ne. size(wfreqs)) stop 'fricker: size mismatch freqs, wfreqs'
fmax = 6.45 * fpeak; nfreq = size(freqs);

do ifreq =1, nfreq
        a = (1.0,0.0) * freqs(ifreq)  + (0.0,1.0) * 0.0e0   
        b = (a/fpeak)**2;
        wfreqs(ifreq) = (b * exp(-b)) * cmplx((2.0e0 / (fpeak * pi**0.5e0)),0.0)
enddo
        
end subroutine fricker

subroutine ricker( &
! this subroutine returns a ricker wavelet with a peak at fqdom 
                  wavelet, &! [ntw] a ricker wavelet towards beginning of arrray
                  wavelet0, &! [ntw0] wavelet such that max is a center of array
                  ntw, &! number of samples in wavelet
                  ntw0, &! number of samples in wavelet0 (must be odd) 
                  fqdom, &! dominant frequency
                  dt, &! time sampling interval
                  diff_flag, & ! time derivative of the ricker is diff_flag ==1
                  norm_flag & ! if eq 1, then normalizes the wavelet to max 1
                        )
implicit none
real,intent(in)                         :: fqdom, dt
integer,intent(in)                      :: ntw, ntw0, diff_flag
logical, intent(in)                     :: norm_flag
real                                    :: pf
real                                    :: t, tsquare
real, parameter                         :: pi = 3.1416
integer                                 :: it, loc, nt, nthalf
real, intent(out)                       :: wavelet(:), wavelet0(:)
real, allocatable                       :: wav(:)
optional                                :: wavelet, ntw, wavelet0, ntw0, diff_flag, norm_flag

nt = 0
if(dt .eq. 0) call abort_msg('ricker: dt is zero')
if(present(wavelet)) then
        if(present(ntw)) then
                call check_dimension('ricker: wavelet', wavelet, ntw)
                nt = maxval((/ntw, nt/))
        else
                call abort_msg('ricker: need ntw for wavelet')
        endif
endif

if(present(wavelet0)) then
        if(present(ntw0)) then
                if(mod(ntw0, 2) .eq. 0) call abort_msg('ricker: time samples of wavelet0 must be odd')
                call check_dimension('ricker: wavelet0', wavelet0, ntw0)
                nt = maxval((/ntw0, nt/))
        else
                call abort_msg('ricker: need ntw0 for wavelet0')
        endif
endif

! some constants
pf = (pi**2)*(fqdom**2)

! a vector is odd number of samples (nt + 1 corresponds to time zero)
allocate(wav(2*nt+1)); wav = zero;
! k = (1 - 2* pf * t^2) * Exp[-pf *t^2]
! Simplify[D[k,t]] 
! FortranForm[Simplify[D[k,t]]]
if(present(diff_flag)) then
        if(diff_flag .eq. 1) then
                ! derivative of the ricker
                do it = 1, 2*nt+ 1
                        tsquare = ((nt + 1 - it) * dt)**2
                        t       = -1.0 * ((nt + 1 - it) * dt)
                        wav(it) = (2.0 * pf * t * (-3.0 + 2.0 * pf * tsquare)) * exp(-1.0 * pf * tsquare)
                enddo
        endif
else
        ! ricker wavelet
        do it = 1, 2*nt + 1
                tsquare = ((nt + 1 - it) * dt)**2
                wav(it) = (1.0 - 2.0 * pf * tsquare) * exp(-1.0e0 * pf * tsquare)
        enddo
endif
call is_nan('ricker: wav', wav)
        
! normalize
if(present(norm_flag)) then
        if(maxval(abs(wav)) .eq. zero) call abort_msg('ricker: maxval(abs(wav)) == 0 ?')
        if(norm_flag) then
                wav = wav / maxval(abs(wav));
        endif
endif
call is_nan('ricker: wav', wav)

if(present(wavelet)) then
        do loc = 1, 2*nt +1 
                if(abs(wav(loc)) .ge. 1e-6) then
                        wavelet = wav(loc:loc+ntw-1)
                        exit
                end if
        enddo
        if(maxval(abs(wavelet)) .eq. zero) call abort_msg('ricker: maxval(abs(wavlet)) == 0 ?')
        call is_nan('ricker: wavelet', wavelet)

endif

if(present(wavelet0)) then
        nthalf = (ntw0 - 1)/2
        wavelet0 = wav(nt + 1 - nthalf : nt + 1 + nthalf)
        if(maxval(abs(wavelet0)) .eq. zero) call abort_msg('ricker: maxval(abs(wavlet0)) == 0 ?')
        call is_nan('ricker: wavelet0', wavelet0)
endif
deallocate(wav)

end subroutine ricker


subroutine ormsby( &
! this subroutine returns a ormsby wavelet f1-f2 ------ f3-f4 
                  wavelet, &! [ntw] an ormsby wavelet towards beginning of arrray
                  wavelet0, &! [ntw0] wavelet such that max is at the center of array (only odd samples)
                  ntw, &! number of samples in wavelet
                  ntw0, &! number of samples in wavelet0 (must be odd) 
                  f1, &! lower bound min freq 
                  f2, &! lower bound max freq
                  f3, &! higher bound min freq
                  f4, &! higher bound max freq
                  dt, &! time sampling interval
                  norm_flag & ! if eq 1, then normalizes the output wavelet to max 1
                        )
implicit none
real,intent(in)                         :: f1, f2, f3, f4, dt
integer,intent(in)                      :: ntw, ntw0
logical, intent(in)                     :: norm_flag
real, intent(out)                       :: wavelet(:), wavelet0(:)
real, allocatable                       :: wav(:)

real                                    :: t, S1, S2, S3, S4, A43, A34, A21, A12
real, parameter                         :: pi = 3.1416

integer                                 :: it, loc, nt, nthalf
optional                                :: wavelet, ntw, wavelet0, ntw0, norm_flag

nt = 0
if(dt .eq. 0) call abort_msg('ormbsy: dt is zero')
if(present(wavelet)) then
        if(present(ntw)) then
                call check_dimension('ormbsy: wavelet', wavelet, ntw)
                nt = maxval((/ntw, nt/))
        else
                call abort_msg('ormbsy: need ntw for wavelet')
        endif
endif

if(present(wavelet0)) then
        if(present(ntw0)) then
                if(mod(ntw0, 2) .eq. 0) call abort_msg('ormbsy: time samples of wavelet0 must be odd')
                call check_dimension('ormbsy: wavelet0', wavelet0, ntw0)
                nt = maxval((/ntw0, nt/))
        else
                call abort_msg('ormbsy: need ntw0 for wavelet0')
        endif
endif

! some constants
A43 = (pi*f4)**2 / (pi*f4 - pi*f3);
A34 = (pi*f3)**2 / (pi*f4 - pi*f3);
A21 = (pi*f2)**2 / (pi*f2 - pi*f1);
A12 = (pi*f1)**2 / (pi*f2 - pi*f1);

! a vector is odd number of samples (nt + 1 corresponds to time zero)
allocate(wav(2*nt+1)); wav = zero;
! ormsby wavelet
do it = 1, 2*nt + 1
        t = ((nt + 1 - it) * dt)
        S4 = (sinc(pi*f4*t))**2
        S3 = (sinc(pi*f3*t))**2
        S2 = (sinc(pi*f2*t))**2
        S1 = (sinc(pi*f1*t))**2

        wav(it) =  (A43 * S4 - A34 * S3) - (A21 * S2 - A12 * S1)
enddo
call is_nan('ormsby: wav', wav)

! normalize
if(present(norm_flag)) then
        if(maxval(abs(wav)) .eq. zero) call abort_msg('ormsby: maxval(abs(wav)) == 0 ?')
        if(norm_flag) then
                wav = wav / maxval(abs(wav));
        endif
endif
call is_nan('ormsby: wav', wav)
! finding the first "non-zero" location and returning the wavelet
if(present(wavelet)) then
        do loc = 1, 2*nt +1
                if(abs(wav(loc)) .ge. 5e-3) then
                        wavelet = wav(loc:loc+ntw-1)
                        exit
                end if
        enddo
        if(maxval(abs(wavelet)) .eq. zero) call abort_msg('ormsby: maxval(abs(wavlet)) == 0 ?')
        call is_nan('ormsby: wavelet', wavelet)
        wavelet(1) = zero; ! first sample is zero
endif

! wavelet0 has peak amplitude at the center
if(present(wavelet0)) then
        nthalf = (ntw0 - 1)/2
        wavelet0 = wav(nt + 1 - nthalf : nt + 1 + nthalf)
        if(maxval(abs(wavelet0)) .eq. zero) call abort_msg('ormsby: maxval(abs(wavlet0)) == 0 ?')
        call is_nan('ormsby: wavelet0', wavelet0)
endif
deallocate(wav)

end subroutine ormsby

function sinc(x)
implicit none
real, intent(in)                        :: x
real                                    :: sinc

! sinc(x) = sin(x)/x
if(x .eq. zero) then
        sinc =  1.0;
else
        sinc = sin(x) / x
endif

end function sinc

subroutine MPfiltered_spike( &
        ! a butterworth minimumphase filter is 
        ! applied to a spike
                  wavelet, &! [ntw] a wavelet minimum phase
                  ntw, &! number of samples in wavelet
                  fmin, &! min. freq
                  fmax, &! max. freq
                  order, & ! order
                  dt, &! time sampling interval
                  norm_flag & ! if eq 1, then normalizes the wavelet to max 1
                        )
implicit none
real,intent(in)                         :: fmin, fmax, dt
integer,intent(in)                      :: ntw, order
logical, intent(in)                     :: norm_flag
real, intent(out)                       :: wavelet(:)
real, allocatable                       :: wav(:)

if(dt .eq. 0) call abort_msg('MPfiltered_spike: dt is zero')
call check_dimension('MPfiltered_spike: wavelet', wavelet, ntw)

allocate(wav(ntw));
! delta function
wav = 0.0; wav(1) = 1.0;
! apply minimum phase butterworth hicut filter
wavelet = apply_butterfilter(gatherin=wav, dt=dt, &
                fmin=fmin, fmax=fmax, order=order,&
                n1=ntw, n2=1, n3=1, filter_type=1)
deallocate(wav)

! normalize
if(maxval(abs(wavelet)) .eq. zero) &
        call abort_msg('MPfiltered_spike: maxval(abs(wavelet)) == 0 ?')
if(norm_flag) then
        wavelet = wavelet / maxval(abs(wavelet));
endif
call is_nan('MPfiltered_spike: wavelet', wavelet)


end subroutine MPfiltered_spike


end module source_wavelets
