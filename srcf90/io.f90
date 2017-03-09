module io
use precision_mod! USEINDEPFILE
use error_messages! USEINDEPFILE
use omp_lib
use unique! USEINDEPFILE
use sorting! USEINDEPFILE
use string_routines! USEINDEPFILE

implicit none
private
integer, parameter, private                                  :: i4=4 ! a 4 byte number = how many record length units?
 ! On GNU,   recl unit = 1 byte, therefore I4 = 4
 ! On Intel, recl unit = 4 bytes, so I4 = 1

! segy headers byte counts
integer, parameter, private                                  :: nbytes_header=       240
integer, parameter, private                                  :: byte0_tracl  =       1
integer, parameter, private                                  :: byte0_tracr  =       5
integer, parameter, private                                  :: byte0_fldr   =       9
integer, parameter, private                                  :: byte0_tracf  =       13      
integer, parameter, private                                  :: byte0_ep     =       17      
integer, parameter, private                                  :: byte0_cdp    =       21      
integer, parameter, private                                  :: byte0_cdpt   =       25      
integer, parameter, private                                  :: byte0_trid   =       29      
integer, parameter, private                                  :: byte0_nvs    =       31      
integer, parameter, private                                  :: byte0_nhs    =       33      
integer, parameter, private                                  :: byte0_duse   =       35      
integer, parameter, private                                  :: byte0_offset =       37      
integer, parameter, private                                  :: byte0_gelev  =       41      
integer, parameter, private                                  :: byte0_selev  =       45      
integer, parameter, private                                  :: byte0_sdepth =       49      
integer, parameter, private                                  :: byte0_gdel   =       53      
integer, parameter, private                                  :: byte0_sdel   =       57      
integer, parameter, private                                  :: byte0_swdep  =       61      
integer, parameter, private                                  :: byte0_gwdep  =       65      
integer, parameter, private                                  :: byte0_scalel =       69      
integer, parameter, private                                  :: byte0_scalco =       71      
integer, parameter, private                                  :: byte0_sx     =       73      
integer, parameter, private                                  :: byte0_sy     =       77      
integer, parameter, private                                  :: byte0_gx     =       81      
integer, parameter, private                                  :: byte0_gy     =       85      
integer, parameter, private                                  :: byte0_counit =       89      
integer, parameter, private                                  :: byte0_wevel  =       91      
integer, parameter, private                                  :: byte0_swevel =       93      
integer, parameter, private                                  :: byte0_sut    =       95      
integer, parameter, private                                  :: byte0_gut    =       97      
integer, parameter, private                                  :: byte0_sstat  =       99
integer, parameter, private                                  :: byte0_gstat  =       101     
integer, parameter, private                                  :: byte0_tstat  =       103     
integer, parameter, private                                  :: byte0_laga   =       105     
integer, parameter, private                                  :: byte0_lagb   =       107     
integer, parameter, private                                  :: byte0_delrt  =       109     
integer, parameter, private                                  :: byte0_muts   =       111     
integer, parameter, private                                  :: byte0_mute   =       113     
integer, parameter, private                                  :: byte0_ns     =       115     
integer, parameter, private                                  :: byte0_dt     =       117     
integer, parameter, private                                  :: byte0_gain   =       119     
integer, parameter, private                                  :: byte0_igc    =       121     
integer, parameter, private                                  :: byte0_igi    =       123     
integer, parameter, private                                  :: byte0_corr   =       125     
integer, parameter, private                                  :: byte0_sfs    =       127     
integer, parameter, private                                  :: byte0_sfe    =       129     
integer, parameter, private                                  :: byte0_slen   =       131     
integer, parameter, private                                  :: byte0_styp   =       133     
integer, parameter, private                                  :: byte0_stas   =       135     
integer, parameter, private                                  :: byte0_stae   =       137     
integer, parameter, private                                  :: byte0_tatyp  =       139     
integer, parameter, private                                  :: byte0_afilf  =       141     
integer, parameter, private                                  :: byte0_afils  =       143     
integer, parameter, private                                  :: byte0_nofilf =       145     
integer, parameter, private                                  :: byte0_nofils =       147     
integer, parameter, private                                  :: byte0_lcf    =       149     
integer, parameter, private                                  :: byte0_hcf    =       151     
integer, parameter, private                                  :: byte0_lcs    =       153     
integer, parameter, private                                  :: byte0_hcs    =       155     
integer, parameter, private                                  :: byte0_year   =       157     
integer, parameter, private                                  :: byte0_day    =       159     
integer, parameter, private                                  :: byte0_hour   =       161     
integer, parameter, private                                  :: byte0_minute =       163     
integer, parameter, private                                  :: byte0_sec    =       165     
integer, parameter, private                                  :: byte0_timbas =       167     
integer, parameter, private                                  :: byte0_trwf   =       169     
integer, parameter, private                                  :: byte0_grnors =       171     
integer, parameter, private                                  :: byte0_grnofr =       173     
integer, parameter, private                                  :: byte0_grnlof =       175     
integer, parameter, private                                  :: byte0_gaps   =       177     
integer, parameter, private                                  :: byte0_ofrav  =       179     
integer, parameter, private                                  :: byte0_d1     =       181     
integer, parameter, private                                  :: byte0_f1     =       185     
integer, parameter, private                                  :: byte0_d2     =       189     
integer, parameter, private                                  :: byte0_f2     =       193     
integer, parameter, private                                  :: byte0_ungpow =       197 
integer, parameter, private                                  :: byte0_unscale =      201


public                                                       :: makesu, readbin, writebin, &
                                                                writedat, &
                                                                readsu_common, readsu_trace, readsu_data_fast
interface writebin
        module procedure         writebinreal
        module procedure         writebincomplex
end interface writebin
interface writedat
        module procedure         writedatreal
end interface writedat
interface readbin
        module procedure         readbinreal
        module procedure         readbincomplex
end interface readbin


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! INPUT - OUTPUT fucntions - binary and su format
contains

subroutine ioerr_msg(ioerr, fname, msg)
implicit none
integer, intent(in)             :: ioerr
character(len=*), intent(in)    :: fname, msg
logical                         :: stop_flag

stop_flag = .false.

if(ioerr.gt.0.0) then
        write(*,*) 'error in accessing a file'
        stop_flag = .true.
endif
if(ioerr.lt.0.0) then
        write(*,*)'end of a file reached'
        stop_flag = .true.
endif

if(stop_flag) then
        write(*,*) 'filename: ', fname
        call abort_msg(msg)
endif
end subroutine ioerr_msg

subroutine readbinreal(fname, & ! name of the binary file
                 dat, & ! data to be read
                 n1, & ! dimension of the data (maximum records to be read from file)
                 n1max, & ! maximum number of records present in the file
                 fileid) ! optional fileid
! subroutine reads a binary file from the disk

character(len=*), intent(in)                    :: fname
integer, optional, intent(in)                   :: n1
integer, optional, intent(out)                  :: n1max
real, optional, intent(out)                     :: dat(:)
integer, intent(in), optional                   :: fileid

real                                            :: tempdat
integer                                         :: i1, ioerr, fileidin, ioerr2, fsize 

if(present(fileid)) then; fileidin = fileid; else; fileidin = 10; endif

! inquire file size in bytes
open(fileidin, file=trim(fname), access='direct', recl=1*i4, status='old', iostat=ioerr)
inquire(file=trim(fname), size = fsize)
if(present(n1max)) n1max = fsize/4

! output to zero
if(present(dat)) then
        if(size(dat,1) .gt. fsize/4) call abort_msg('readbinreal: &
                & not enough 4 byte records in the file, '//trim(fname))
        if(present(n1)) then
                call check_dimension('readbinreal', dat, n1)
        endif
        dat = rzero_de;
        open(fileidin, file=trim(fname), access='direct', recl=1*i4, status='old', iostat=ioerr)
        call ioerr_msg(ioerr= ioerr, fname = fname, msg = 'readbinreal')
        do i1 = 1, size(dat,1)
                read(fileidin,rec=i1, iostat = ioerr2) tempdat
                if(ioerr2 .eq. 0) then
                        dat(i1) = tempdat
                else
                        call abort_msg('readbinreal: error in reading record')
                endif
        enddo
endif
close(fileidin)

end subroutine readbinreal

subroutine readbincomplex(fname, dat, n1, fileid)
implicit  none
character(len=*), intent(in)                    :: fname
integer,          intent(in)                    :: n1
complex,             intent(out)                :: dat(:)
integer, intent(in), optional                   :: fileid

integer                                         :: i1,  fileidin

real, allocatable, dimension(:)                 :: datR, datI

! output to zero
call check_dimension('readbinreal', real(dat), n1)
dat = cmplx(0.0, 0.0);

if(present(fileid)) then; fileidin = fileid; else; fileidin = 10; endif



allocate(datR(n1))
allocate(datI(n1))
        
call readbinreal(fname=fname//'R',dat=datR,n1=n1, fileid= fileidin)
call readbinreal(fname=fname//'I',dat=datI,n1=n1, fileid= fileidin)

do i1=1, n1
        dat(i1) = cmplx(datR(i1), datI(i1))
enddo
 
deallocate(datR, datI)

end subroutine readbincomplex


subroutine writebinreal(fname, & ! name of the binaryfile
         dat, & ! data to be written 
         n1, & ! dimension of the data, optional
         fileid & ! optional fileid
                )
! subroutine writes binary data to disk

character(len=*), intent(in)                    :: fname
integer, intent(in), optional                   :: n1
real, intent(in)                                :: dat(:)
integer, intent(in), optional                   :: fileid

integer                                         :: n1in, ioerr, fileidin

if(present(n1)) then; n1in = n1; else; n1in = size(dat, 1); endif
if(present(fileid)) then; fileidin = fileid; else; fileidin = 10; endif
call check_dimension('writebinreal', dat, n1in)

open(fileidin, file=fname, access='direct', recl=n1in*i4, status='replace', iostat=ioerr)
call ioerr_msg(ioerr= ioerr, fname = fname, msg = 'writebinreal')
write(fileidin,rec=1) real(dat(1:n1in),r_sp)
close(fileidin)

end subroutine writebinreal


subroutine writedatreal(fname, & ! name of the binaryfile
         dat, & ! data to be written 
         n1, & ! dimension of the data, optional
         fileid & ! optional fileid
                )
! subroutine writes a vector to disk

character(len=*), intent(in)                    :: fname
integer, intent(in), optional                   :: n1
real, intent(in)                                :: dat(:)
integer, intent(in), optional                   :: fileid

integer                                         :: n1in,  fileidin, ioerr

if(present(n1)) then; n1in = n1; else; n1in = size(dat, 1); endif
if(present(fileid)) then; fileidin = fileid; else; fileidin = 10; endif
call check_dimension('writedatreal', dat, n1in)

open(fileidin, file=fname, action='write', status='replace', iostat=ioerr)
call ioerr_msg(ioerr= ioerr, fname = fname, msg = 'writedatreal')
write(fileidin,*) real(dat(1:n1in),r_sp)
close(fileidin)

end subroutine writedatreal


subroutine writebincomplex(fname, dat, n1, fileid)
implicit none
character(len=*), intent(in)                    :: fname
integer, intent(in), optional                   :: n1
complex, intent(in)                             :: dat(:)
integer, intent(in), optional                   :: fileid

integer                                         :: n1in,  fileidin

real, dimension(:), allocatable                 :: datR, datI


if(present(n1)) then; n1in = n1; else; n1in = size(dat, 1); endif
if(present(fileid)) then; fileidin = fileid; else; fileidin = 10; endif
call check_dimension('writebincomplex', real(dat), n1in)

allocate(datR(n1in))
allocate(datI(n1in))
datR = real(dat); datI = aimag(dat) 

call writebinreal(fname=fname//'R',dat=datR,n1=n1in,fileid=fileidin)
call writebinreal(fname=fname//'I',dat=datI,n1=n1in,fileid=fileidin)

deallocate(datR, datI)

end subroutine writebincomplex


subroutine readsu_common( &
        ! this subroutines reads headers from a SU file common to all traces
                        fname,  & ! file name of the SU file
                        deltat, & ! time sampling
                        nt, & !   number of samples in each trace
                        ntrc, & ! number of traces/ records present in the SU file
                        nsource, & ! total number of sources
                        nreceiver,  & ! total number of receivers
                        fileid, & ! fileid of the file to be opened
                        success_flag, & ! true if the all headers present are read 
			verbose_flag & ! info print on screen
                        )
implicit none
integer, intent(out), optional                  :: nt, nsource, nreceiver, ntrc
real, intent(out), optional                     :: deltat
integer, intent(in), optional                   :: fileid
character(len=*), intent(in)                    :: fname
logical, intent(out)                            :: success_flag
logical, intent(in),optional                    :: verbose_flag
                
integer                                         :: ioerr, fileidin, ntot4, ioerr2, ntot_records, itrc
integer                                         :: nsrc, nrec

integer, allocatable                            :: sxall(:), gxall(:), gyall(:), syall(:)
integer, allocatable                            :: gelevall(:), selevall(:)
real, allocatable                               :: tempx(:), tempz(:)
integer, allocatable                            :: traclall(:), tracfall(:), fldrall(:), ntrcall(:)

integer*2, allocatable                          :: scalcoall(:), scalelall(:)
real, allocatable                               :: scalcoall_factor(:), scalelall_factor(:)
logical, allocatable                            :: temp_logical(:)
logical						:: vflag

! segy headers
integer(i_dp)                                   :: tracl,tracr,fldr,tracf,ep,cdp,cdpt
integer(i_sp)                                   :: trid,nvs,nhs,duse 
integer(i_dp)                                   :: offset,gelev,selev,sdepth,gdel,sdel,swdep,gwdep
integer(i_sp)                                   :: scalel,scalco
integer(i_dp)                                   :: sx,sy,gx,gy
integer(i_sp)                                   :: counit,wevel,swevel,sut,gut,sstat,gstat, &
                                                   tstat,laga,lagb,delrt,muts,mute,ns,dt 
integer(i_sp)                                   :: gain,igc,igi,corr,sfs,sfe,slen,styp,stas, &
                                                   stae,tatyp,afilf,afils,nofilf,nofils,lcf,hcf, &
                                                   lcs,hcs,year,day,hour,minute,sec,timbas, &
                                                   trwf,grnors,grnofr,grnlof,gaps,otrav 
real(r_sp)                                      :: d1,f1,d2,f2,ungpow,unscale
integer(i_sp)                                   :: mark,unass(17)


if(present(fileid)) then; fileidin = fileid; else; fileidin = OMP_get_thread_num(); endif
if(present(verbose_flag)) then; vflag = verbose_flag; else; vflag = .false.; endif

! inquire total 4-byte records present in the SU file
call readbinreal(fname = trim(fname), n1max = ntot4)

! opening fil with recl = 2 bytes for ns and dt
open(fileidin, file=trim(fname), access='direct', recl=1*i4/2, status='old', iostat=ioerr)
call ioerr_msg(ioerr= ioerr, fname = fname, msg = 'readsu_common: opening file')
read(fileidin, rec=(byte0_ns+1)/2, iostat = ioerr2) ns
call ioerr_msg(ioerr= ioerr2, fname = fname, msg = 'readsu_common: cannot read ns')
read(fileidin, rec=(byte0_dt+1)/2, iostat = ioerr2) dt
call ioerr_msg(ioerr= ioerr2, fname = fname, msg = 'readsu_common: cannot read dt')
close(fileidin)

! by default success_flag is true
! if any of the headers not present.. flag will be false
success_flag = .true.

if(vflag) write(*,*)'readsu_common: headers in SU filename ', trim(fname)

if(present(nt)) then
        nt = int(ns,i_dp)
        if(nt .eq. 0) then
                if(vflag) write(*,*) '     header ns not set '
                success_flag = .false.
	else
		if(vflag) write(*,*) '     ns: ', nt
        endif
endif

if(present(deltat)) then
        deltat = dt * 1e-6
        if(deltat .eq. 0) then
                if(vflag) write(*,*) '     header dt not set in SU file ',trim(fname)
                success_flag = .false.
	else
		if(vflag) write(*,*) '     dt: ', deltat
        endif
endif

ntot_records = ntot4 / (nbytes_header/4 + ns)
if(present(ntrc)) then
        ntrc = ntot_records
        if(ntrc .eq. 0) then
                if(vflag) write(*,*) '     cannot determine ntrc from SU file ',trim(fname)
                success_flag = .false.
	else
		if(vflag) write(*,*) '     ntrc: ', ntrc
        endif
endif

if(present(nsource) .or. present(nreceiver)) then

        allocate( &
                traclall(ntot_records), fldrall(ntot_records),&
                tracfall(ntot_records), &
                scalcoall(ntot_records), scalelall(ntot_records),&
                scalcoall_factor(ntot_records), scalelall_factor(ntot_records),&
                ntrcall(ntot_records))

        ! open SU file to read headers of 4 byte recl
        open(fileidin, file=trim(fname), access='direct', recl=1*i4, status='old', iostat=ioerr)
        call ioerr_msg(ioerr= ioerr, fname = fname, msg = 'readsu_common: opening file')
        do itrc = 1, ntot_records
                read(fileidin, rec=((byte0_tracl + (itrc-1)*(nbytes_header+ns*4))+3)/4, &
                                                iostat = ioerr2) traclall(itrc)
                read(fileidin, rec=((byte0_tracf + (itrc-1)*(nbytes_header+ns*4))+3)/4, &
                                                iostat = ioerr2) tracfall(itrc)
                read(fileidin, rec=((byte0_fldr + (itrc-1)*(nbytes_header+ns*4))+3)/4, &
                                                iostat = ioerr2) fldrall(itrc)
        enddo
        close(fileidin)

        ! open SU file to read headers of 2 byte recl
        open(fileidin, file=trim(fname), access='direct', recl=1*i4/2, status='old', iostat=ioerr)
        call ioerr_msg(ioerr= ioerr, fname = fname, msg = 'readsu_common: opening file')
        do itrc = 1, ntot_records
                read(fileidin, rec=((byte0_scalco + (itrc-1)*(nbytes_header+ns*4))+1)/2, &
                                                iostat = ioerr2) scalcoall(itrc)
                read(fileidin, rec=((byte0_scalel + (itrc-1)*(nbytes_header+ns*4))+1)/2, &
                                                iostat = ioerr2) scalelall(itrc)
        enddo
        close(fileidin)

        ! scaling factors present or not?
        if(any(abs(scalcoall) .eq. izero_sp)) then
                scalcoall_factor = 1.e0
        else    
                scalcoall_factor = abs(real(scalcoall,r_de))**(scalcoall/abs(scalcoall))
        endif
        if(any(abs(scalelall) .eq. izero_sp)) then
                scalelall_factor = 1.e0
        else
                scalelall_factor = abs(real(scalelall,r_de))**(scalelall/abs(scalelall))
        endif

endif

if(present(nsource)) then
        allocate( &
                sxall(ntot_records), syall(ntot_records), &
                selevall(ntot_records))

        ! open SU file to read headers of 4 byte recl
        open(fileidin, file=trim(fname), access='direct', recl=1*i4, status='old', iostat=ioerr)
        call ioerr_msg(ioerr= ioerr, fname = fname, msg = 'readsu_common: opening file')
        do itrc = 1, ntot_records
                read(fileidin, rec=((byte0_sx + (itrc-1)*(nbytes_header+ns*4))+3)/4, &
                                                iostat = ioerr2) sxall(itrc)
                read(fileidin, rec=((byte0_sy + (itrc-1)*(nbytes_header+ns*4))+3)/4, &
                                                iostat = ioerr2) syall(itrc)
                read(fileidin, rec=((byte0_selev + (itrc-1)*(nbytes_header+ns*4))+3)/4, &
                                                iostat = ioerr2) selevall(itrc)
        enddo
        close(fileidin)

        allocate(temp_logical(ntot_records))
        allocate(tempx(ntot_records), tempz(ntot_records))
        ! applying scales
        tempx = real(sxall,r_de) * scalcoall_factor 
        tempz = real(selevall,r_de) * scalelall_factor 
        call unique_mask(array1=tempx, array2=tempz, unique_mask=temp_logical)
        nsrc = count(temp_logical)

        if(nsrc .ne. izero_de) then
                nsource = nsrc
		if(vflag) write(*,*) '     number of source positions: ', nsrc
        else
                success_flag = .false.
                if(vflag) write(*,*) '     nsource cannot be estimated '
                nsource = 0
        endif
        deallocate(tempx, tempz, temp_logical)
endif

if(present(nreceiver)) then
        allocate( &
                gxall(ntot_records), gyall(ntot_records), &
                gelevall(ntot_records))

        ! open SU file to read headers of 4 byte recl
        open(fileidin, file=trim(fname), access='direct', recl=1*i4, status='old', iostat=ioerr)
        call ioerr_msg(ioerr= ioerr, fname = fname, msg = 'readsu_common: opening file')
        do itrc = 1, ntot_records
                read(fileidin, rec=((byte0_gx + (itrc-1)*(nbytes_header+ns*4))+3)/4, &
                                                iostat = ioerr2) gxall(itrc)
                read(fileidin, rec=((byte0_gy + (itrc-1)*(nbytes_header+ns*4))+3)/4, &
                                                iostat = ioerr2) gyall(itrc)
                read(fileidin, rec=((byte0_gelev + (itrc-1)*(nbytes_header+ns*4))+3)/4, &
                                                iostat = ioerr2) gelevall(itrc)
        enddo
        close(fileidin)

        allocate(temp_logical(ntot_records))
        allocate(tempx(ntot_records), tempz(ntot_records))
        ! applying scales
        tempx = real(gxall,r_de) * scalcoall_factor 
        tempz = real(gelevall,r_de) * scalelall_factor 
        call unique_mask(array1=tempx, array2=tempz, unique_mask=temp_logical)
        nrec = count(temp_logical)

        if(nrec .ne. izero_de) then
                nreceiver = nrec
		if(vflag) write(*,*) '     number of receiver positions: ', nrec
        else
                success_flag = .false.
                if(vflag) write(*,*) '     nsource cannot be estimated '
                nreceiver = 0
        endif
        deallocate(tempx, tempz, temp_logical)
endif
if(allocated(traclall)) deallocate(traclall)
if(allocated(fldrall)) deallocate(fldrall)
if(allocated(tracfall)) deallocate(tracfall)
if(allocated(scalcoall)) deallocate(scalcoall)
if(allocated(scalelall_factor)) deallocate(scalelall_factor)
if(allocated(scalcoall_factor)) deallocate(scalcoall_factor)
if(allocated(scalelall)) deallocate(scalelall)
if(allocated(ntrcall)) deallocate(ntrcall)
if(allocated(sxall)) deallocate(sxall)
if(allocated(syall)) deallocate(syall)
if(allocated(gxall)) deallocate(gxall)
if(allocated(gyall)) deallocate(gyall)
if(allocated(selevall)) deallocate(selevall)
if(allocated(gelevall)) deallocate(gelevall)

end subroutine readsu_common


subroutine readsu_trace( &
        ! this subroutines reads trace headers and data from a su file
                        fname,  & ! file name of the su file
                        ntmax, & ! maximum number of time samples in each trace to be read
                                 ! it should be less than ns of the su file, ofc time sampling is same
                        nsourcemax, & ! max source positions (< total source positions in su file)
                                      ! the source positions are sorted and first nsourcemax sources will be selected 
                        nreceivermax,  & ! max receiver positions (< total receiver positions in su file)
                        sr_pos, & ! source- receiver acquisition geometry [nrec, nsrc, 5]
                        csgs, & ! common shot gathers [ntmax *  nreceivermax * nsourcemax]
                        fileid & ! fileid of the file to be opened
                        )
implicit none
integer, intent(in), optional                   :: ntmax
integer, intent(in)                             :: nsourcemax, nreceivermax
integer, intent(in), optional                   :: fileid
character(len=*), intent(in)                    :: fname
real, intent(out), optional                     :: sr_pos(:,:,:)
real, intent(out), optional                     :: csgs(:)
                
integer                                         :: ioerr, fileidin, ntot4, ioerr2, ntot_records, itrc, it
integer                                         :: isrc, irec, k

integer, allocatable                            :: sxall(:), gxall(:), gyall(:), syall(:)
integer, allocatable                            :: gelevall(:), selevall(:)
real, allocatable                               :: tempx(:), tempz(:), temp_src(:), temp_rec(:), sr_posWtemp(:,:)
real, allocatable                               :: sourcesx(:), sourcesz(:), receiversx(:), receiversz(:)
real(r_sp), allocatable				:: datatrace(:,:)
integer, allocatable                            :: traclall(:), tracfall(:), fldrall(:)
integer*2, allocatable                          :: scalcoall(:), scalelall(:), tridall(:)
real, allocatable                               :: scalcoall_factor(:), scalelall_factor(:)
logical, allocatable                            :: temp_logical(:)

! segy headers
integer(i_dp)                                   :: tracl,tracr,fldr,tracf,ep,cdp,cdpt
integer(i_sp)                                   :: trid,nvs,nhs,duse 
integer(i_dp)                                   :: offset,gelev,selev,sdepth,gdel,sdel,swdep,gwdep
integer(i_sp)                                   :: scalel,scalco
integer(i_dp)                                   :: sx,sy,gx,gy
integer(i_sp)                                   :: counit,wevel,swevel,sut,gut,sstat,gstat, &
                                                   tstat,laga,lagb,delrt,muts,mute,ns,dt 
integer(i_sp)                                   :: gain,igc,igi,corr,sfs,sfe,slen,styp,stas, &
                                                   stae,tatyp,afilf,afils,nofilf,nofils,lcf,hcf, &
                                                   lcs,hcs,year,day,hour,minute,sec,timbas, &
                                                   trwf,grnors,grnofr,grnlof,gaps,otrav 
real(r_sp)                                      :: d1,f1,d2,f2,ungpow,unscale
integer(i_sp)                                   :: mark,unass(17)


if(present(fileid)) then; fileidin = fileid; else; fileidin = omp_get_thread_num(); endif
if(present(sr_pos)) then
        call check_dimension('readsu_trace: size of sr_pos', sr_pos, nreceivermax, nsourcemax, 5)
        sr_pos = 0.0;
endif
if(present(csgs)) then
        if(.not.(present(ntmax))) call abort_msg('readsu_trace: ntmax needed for csgs output')
        call check_dimension('readsu_trace: size of csgs', csgs, ntmax * nreceivermax * nsourcemax )      
        csgs = rzero_de
endif

! inquire total 4-byte records present in the su file
call readbinreal(fname = trim(fname), n1max = ntot4)


! opening fil with recl = 2 bytes for ns and dt
open(fileidin, file=trim(fname), access='direct', recl=1*i4/2, status='old', iostat=ioerr)
call ioerr_msg(ioerr= ioerr, fname = fname, msg = 'readsu_trace: opening file')
read(fileidin, rec=(byte0_ns+1)/2, iostat = ioerr2) ns
call ioerr_msg(ioerr= ioerr2, fname = fname, msg = 'readsu_trace: cannot read ns')
read(fileidin, rec=(byte0_dt+1)/2, iostat = ioerr2) dt
call ioerr_msg(ioerr= ioerr2, fname = fname, msg = 'readsu_trace: cannot read dt')
close(fileidin)
ntot_records = ntot4 / (nbytes_header/4 + ns)
if(present(ntmax)) then
        if(ntmax .gt. ns) call abort_msg('readsu_trace: ns of the su file '//trim(fname)//' < ntmax')
endif


if(present(csgs) .or. present(sr_pos)) then
        allocate( &
                sxall(ntot_records), syall(ntot_records), &
                gxall(ntot_records), gyall(ntot_records), &
                selevall(ntot_records), gelevall(ntot_records), &
                traclall(ntot_records), fldrall(ntot_records),&
                tracfall(ntot_records), &
                scalcoall(ntot_records), scalelall(ntot_records), &
                scalcoall_factor(ntot_records), scalelall_factor(ntot_records), &
                tridall(ntot_records), &
		datatrace(ns,ntot_records)&
                )

	open(fileidin, file=trim(adjustl(fname)), access='direct', recl=(ns + 60)*i4, &
			status='old', iostat = ioerr)
	do itrc = 1, ntot_records
		! zero to all headers
                tracl=0; tracr=0; fldr=0; tracf=0; ep=0; cdp=0; cdpt=0; trid=0;
                nvs=0; nhs=0; duse=0; offset=0; gelev=0; selev=0;
                sdepth=0; gdel=0; sdel=0; swdep=0; gwdep=0; scalel=0; scalco=0; sx=0;
                sy=0; gx=0; gy=0; counit=0; wevel=0; swevel=0; sut=0; 
                gut=0; sstat=0; gstat=0; tstat=0; laga=0
                lagb=0; delrt=0; muts=0; mute=0; ns=0; dt=0; gain=0; igc=0; igi=0; corr=0; sfs=0; sfe=0
                slen=0; styp=0; stas=0; stae=0; tatyp=0; afilf=0; afils=0; nofilf=0; 
                nofils=0; lcf=0; hcf=0; lcs=0
                hcs=0; year=0; day=0; hour=0; minute=0; sec=0; timbas=0; trwf=0; grnors=0; grnofr=0
                grnlof=0; gaps=0; otrav=0; d1=0.0; f1=0.0; d2=0.0; f2=0.0; ungpow=0.0; unscale=0.0
                mark=0; unass(:)=0;

		read(fileidin, rec=itrc, iostat = ioerr2) tracl,tracr,fldr,tracf,ep,cdp,cdpt, &
			trid,nvs,nhs,duse,offset,gelev,selev,sdepth, &
			gdel,sdel,swdep,gwdep,scalel,scalco, &
			sx,sy,gx,gy,counit,wevel,swevel,sut,gut,sstat,gstat, &
			tstat,laga,lagb,delrt,muts,mute,ns,dt, &
			gain,igc,igi,corr,sfs,sfe,slen,styp,stas, &
			stae,tatyp,afilf,afils,nofilf,nofils,lcf,hcf, &
			lcs,hcs,year,day,hour,minute,sec,timbas, &
			trwf,grnors,grnofr,grnlof,gaps,otrav, &
			d1,f1,d2,f2,ungpow,unscale,mark,(unass(k),k=1,17), &
			(datatrace(it, itrc),it=1,int(ns,i_dp))
		
                        tridall(itrc) = trid
			sxall(itrc) = sx
			syall(itrc) = sy
			selevall(itrc) = selev
			gxall(itrc) = gx
			gyall(itrc) = gy
			gelevall(itrc) = gelev
			traclall(itrc) = tracl
			tracfall(itrc) = tracf
			fldrall(itrc) = fldr
			scalcoall(itrc) = scalco
			scalelall(itrc) = scalel
	enddo
	close(fileidin)

        ! scaling factors present or not?
        if(any(abs(scalcoall) .eq. izero_sp)) then
                scalcoall_factor = 1.e0
        else    
                scalcoall_factor = abs(real(scalcoall,r_de))**(scalcoall/abs(scalcoall))
        endif
        if(any(abs(scalelall) .eq. izero_sp)) then
                scalelall_factor = 1.e0
        else
                scalelall_factor = abs(real(scalelall,r_de))**(scalelall/abs(scalelall))
        endif

        ! sourcesx sourcesz receiersx receiversz
        allocate(sourcesx(nsourcemax), sourcesz(nsourcemax))
        allocate(receiversx(nreceivermax), receiversz(nreceivermax))
        allocate(temp_logical(ntot_records))
        allocate(tempx(ntot_records), tempz(ntot_records))
        tempx = real(sxall,r_de) * scalcoall_factor
        tempz = real(selevall,r_de) * scalelall_factor 
        call unique_mask(array1=tempx, array2=tempz, unique_mask=temp_logical)
        if(nsourcemax .gt. count(temp_logical)) call abort_msg('readsu_trace: &
                        &sources in su file '//trim(fname)//' < nsourcemax'//i2str(count(temp_logical)))
	isrc = 0; sourcesx = rzero_de; sourcesz = rzero_de
        do itrc = 1, ntot_records
                if(temp_logical(itrc)) then
                        isrc = isrc + 1
                        if(isrc .gt. nsourcemax) exit
                        sourcesx(isrc) = tempx(itrc)
                        sourcesz(isrc) = tempz(itrc)
                endif
        enddo

        tempx = real(gxall,r_de) * scalcoall_factor 
        tempz = real(gelevall,r_de) * scalelall_factor 
        call unique_mask(array1=tempx, array2=tempz, unique_mask=temp_logical)
        if(nreceivermax .gt. count(temp_logical)) call abort_msg('readsu_trace: &
                        &receivers in su file '//trim(fname)//' < nreceivermax'//i2str(count(temp_logical)))
	irec = 0; receiversx = rzero_de; receiversz = rzero_de;
        do itrc = 1, ntot_records
                if(temp_logical(itrc)) then
                        irec = irec + 1
                        if(irec .gt. nreceivermax) exit
                        receiversx(irec) = tempx(itrc)
                        receiversz(irec) = tempz(itrc)
                endif
        enddo
        deallocate(tempx, tempz, temp_logical)

        ! sr_pos(:,:,1:4)
        if(present(sr_pos)) then
                do isrc=1,nsourcemax 
                        sr_pos(:,isrc,1) = sourcesx(isrc)
                        sr_pos(:,isrc,2) = sourcesz(isrc)
                enddo
                do irec=1,nreceivermax
                        sr_pos(irec,:,3) = receiversx(irec)
                        sr_pos(irec,:,4) = receiversz(irec)
                enddo
        endif


        ! csgs and/or sr_pos(:,:,5)
        ! open su file to read headers of 4 byte recl
        allocate(temp_src(nsourcemax))
        allocate(temp_rec(nreceivermax))
        allocate(sr_posWtemp(nreceivermax,nsourcemax))
        ! default is all the weights are zero
        sr_posWtemp  = 0.0;
        do itrc = 1, ntot_records
                if(tridall(itrc) .ne. itwo_sp) then
                        temp_src = (sourcesx - real(sxall(itrc),r_de) * scalcoall_factor(itrc) )**2 &
                                + (sourcesz - real(selevall(itrc),r_de) * scalelall_factor(itrc))**2
                        temp_rec = (receiversx - real(gxall(itrc),r_de) * scalcoall_factor(itrc) )**2 &
                                + (receiversz - real(gelevall(itrc),r_de) * scalelall_factor(itrc) )**2
                        if((minval(temp_src) .le. 1e-6) .and. minval(temp_rec) .le. 1e-6) then
                                isrc = minloc(temp_src, dim=1)
                                irec = minloc(temp_rec, dim=1)
                                if(sr_posWtemp(irec, isrc) .eq. 1.e0) then
                                        call abort_msg(&
                                'readsu_trace: multiple traces for same source & receiver in su file '//trim(fname))
                                else
                                        sr_posWtemp(irec, isrc) = 1.e0
                                endif
                                if(present(csgs)) then
                                        csgs(1 + (irec-1)*ntmax + (isrc-1)*ntmax*nreceivermax : &
                                        ntmax + (irec-1)*ntmax + (isrc-1)*ntmax*nreceivermax) = real(datatrace(1:ntmax,itrc),r_de)
                                endif
                        endif
                endif
        enddo
        if(present(sr_pos)) then
                sr_pos(:, :, 5) = sr_posWtemp;
        endif
        deallocate(sr_posWtemp)
        deallocate(temp_src, temp_rec)
endif
deallocate(receiversx, receiversz, sourcesx, sourcesz)
if(allocated(traclall)) deallocate(traclall)
if(allocated(fldrall)) deallocate(fldrall)
if(allocated(tracfall)) deallocate(tracfall)
if(allocated(scalcoall)) deallocate(scalcoall)
if(allocated(scalelall)) deallocate(scalelall)
if(allocated(scalcoall_factor)) deallocate(scalcoall_factor)
if(allocated(scalelall_factor)) deallocate(scalelall_factor)
if(allocated(sxall)) deallocate(sxall)
if(allocated(syall)) deallocate(syall)
if(allocated(gxall)) deallocate(gxall)
if(allocated(gyall)) deallocate(gyall)
if(allocated(selevall)) deallocate(selevall)
if(allocated(gelevall)) deallocate(gelevall)
if(allocated(tridall)) deallocate(tridall)
if(allocated(datatrace)) deallocate(datatrace)

end subroutine readsu_trace

subroutine readsu_nt_nrecords( &
        ! this subroutines reads data from a su file quickly.. headers ignored
                        fname,  & ! file name of the su file
                        nt, & ! output ns from the headers
                        ntrc  &! output number of traces
                        ) bind(c, name="readsu_nt_nrecords") 
implicit none
integer(i_dp), intent(out)                      :: nt, ntrc
character(len=1, kind=C_char), intent(in)       :: fname(*)
integer                                         :: ioerr, ioerr2, fileidin, ntot4, ntot_records
real(r_sp), allocatable				:: datatrace(:,:)
integer(i_sp)                                   :: ns


fileidin = OMP_get_thread_num();

! inquire total 4-byte records present in the su file
call readbinreal(fname = trim(ctofstr(fname)), n1max = ntot4)

! opening fil with recl = 2 bytes for ns and dt
open(fileidin, file=trim(ctofstr(fname)), access='direct', recl=1*i4/2, status='old', iostat=ioerr)
call ioerr_msg(ioerr= ioerr, fname = ctofstr(fname), msg = "readsu_nt_nrecords: opening file")
read(fileidin, rec=(byte0_ns+1)/2, iostat = ioerr2) ns
call ioerr_msg(ioerr= ioerr2, fname = ctofstr(fname), msg = "readsu_data_nt_nrecords: cannot read ns")
close(fileidin)
ntot_records = ntot4 / (nbytes_header/4 + ns)

! output
ntrc = int(ntot_records); nt = int(ns)

end subroutine readsu_nt_nrecords




subroutine readsu_data_fast( &
        ! this subroutines reads data from a su file quickly.. headers ignored
                        fname,  & ! file name of the su file
                        nt, & ! ns from the headers
                        ntrc, &! number of traces
                        dat  & ! data  [nt * ntrc]
                        ) bind(c, name="readsu_data_fast") 
implicit none
integer, intent(in)                             :: nt, ntrc
character(len=1, kind=C_char), intent(in)       :: fname(*)
real(r_de), intent(out)                         :: dat(nt * ntrc)
integer                                         :: ioerr, fileidin, ntot4, ioerr2, ntot_records, itrc, it, k
real(r_sp), allocatable				:: datatrace(:,:)

! segy headers
integer(i_dp)                                   :: tracl,tracr,fldr,tracf,ep,cdp,cdpt
integer(i_sp)                                   :: trid,nvs,nhs,duse 
integer(i_dp)                                   :: offset,gelev,selev,sdepth,gdel,sdel,swdep,gwdep
integer(i_sp)                                   :: scalel,scalco
integer(i_dp)                                   :: sx,sy,gx,gy
integer(i_sp)                                   :: counit,wevel,swevel,sut,gut,sstat,gstat, &
                                                   tstat,laga,lagb,delrt,muts,mute,ns,dt 
integer(i_sp)                                   :: gain,igc,igi,corr,sfs,sfe,slen,styp,stas, &
                                                   stae,tatyp,afilf,afils,nofilf,nofils,lcf,hcf, &
                                                   lcs,hcs,year,day,hour,minute,sec,timbas, &
                                                   trwf,grnors,grnofr,grnlof,gaps,otrav 
real(r_sp)                                      :: d1,f1,d2,f2,ungpow,unscale
integer(i_sp)                                   :: mark,unass(17)


fileidin = OMP_get_thread_num();

! inquire total 4-byte records present in the su file
call readbinreal(fname = trim(ctofstr(fname)), n1max = ntot4)


! opening fil with recl = 2 bytes for ns and dt
open(fileidin, file=trim(ctofstr(fname)), access='direct', recl=1*i4/2, status='old', iostat=ioerr)
call ioerr_msg(ioerr= ioerr, fname = ctofstr(fname), msg = 'readsu_data_fast: opening file')
read(fileidin, rec=(byte0_ns+1)/2, iostat = ioerr2) ns
call ioerr_msg(ioerr= ioerr2, fname = ctofstr(fname), msg = 'readsu_data_fast: cannot read ns')
close(fileidin)
ntot_records = ntot4 / (nbytes_header/4 + ns)
if(nt .ne. ns) call abort_msg("readsu_data_fast: ns of the su file '//trim(ctofstr(fname))//' .ne. nt")
if(ntot_records .ne. ntrc) call abort_msg("readsu_data_fast: ntrc .ne. ntot_records")


allocate(datatrace(ns,ntot_records))
open(fileidin, file=trim(adjustl(ctofstr(fname))), access='direct', recl=(ns + 60)*i4, &
                status='old', iostat = ioerr)
do itrc = 1, ntot_records
        ! zero to all headers
        tracl=0; tracr=0; fldr=0; tracf=0; ep=0; cdp=0; cdpt=0; trid=0;
        nvs=0; nhs=0; duse=0; offset=0; gelev=0; selev=0;
        sdepth=0; gdel=0; sdel=0; swdep=0; gwdep=0; scalel=0; scalco=0; sx=0;
        sy=0; gx=0; gy=0; counit=0; wevel=0; swevel=0; sut=0; 
        gut=0; sstat=0; gstat=0; tstat=0; laga=0
        lagb=0; delrt=0; muts=0; mute=0; ns=0; dt=0; gain=0; igc=0; igi=0; corr=0; sfs=0; sfe=0
        slen=0; styp=0; stas=0; stae=0; tatyp=0; afilf=0; afils=0; nofilf=0; 
        nofils=0; lcf=0; hcf=0; lcs=0
        hcs=0; year=0; day=0; hour=0; minute=0; sec=0; timbas=0; trwf=0; grnors=0; grnofr=0
        grnlof=0; gaps=0; otrav=0; d1=0.0; f1=0.0; d2=0.0; f2=0.0; ungpow=0.0; unscale=0.0
        mark=0; unass(:)=0;

        read(fileidin, rec=itrc, iostat = ioerr2) tracl,tracr,fldr,tracf,ep,cdp,cdpt, &
                trid,nvs,nhs,duse,offset,gelev,selev,sdepth, &
                gdel,sdel,swdep,gwdep,scalel,scalco, &
                sx,sy,gx,gy,counit,wevel,swevel,sut,gut,sstat,gstat, &
                tstat,laga,lagb,delrt,muts,mute,ns,dt, &
                gain,igc,igi,corr,sfs,sfe,slen,styp,stas, &
                stae,tatyp,afilf,afils,nofilf,nofils,lcf,hcf, &
                lcs,hcs,year,day,hour,minute,sec,timbas, &
                trwf,grnors,grnofr,grnlof,gaps,otrav, &
                d1,f1,d2,f2,ungpow,unscale,mark,(unass(k),k=1,17), &
                (datatrace(it, itrc),it=1,int(ns,i_dp))
        
enddo
close(fileidin)
! return with default real data type
dat = real(pack(datatrace,.true.), r_de)
if(allocated(datatrace)) deallocate(datatrace)

end subroutine readsu_data_fast



subroutine makesu(dat,& ! data vector (regularly sampled) 
                        ! the data vector is first sorted in time, then by receivers and finally by sources
                 nsrc, & ! number of sources (optional, default is nsrc = 1)
                 nrec, & ! number of receivers (optional, default is nrec =1)
                 SR_pos, & ! a 3D matrix of positions of sources and receivers 
                           ! the ordering of the sources and receivers should match to that of the 
                           ! data vector [nrec, nsrc, 5] -- optional
                 nt, & ! number of time samples (optional, default is length of the data vector)
                 deltat, & ! time sampling interval
                 ! some non-seismic headers
                 dz, & ! sampling in first dimension
                 dx, & ! sampling in second dimension
                 fz, & ! first sample of z direction
                 fx, & ! first sample of x direction
                 fname, & ! name of the file to be saved
                 fileid & !
                         )
implicit none
real, intent(in)                :: dat(:)
character(len=*), intent(in)    :: fname
integer,intent(in)              :: fileid
integer, intent(in)             :: nt, nsrc, nrec
real, intent(in)                :: deltat, SR_pos(:,:,:) 
real, intent(in)                :: fx, fz, dz, dx

optional                        :: fx, fz ! default = 0
optional                        :: dx, dz ! default = 1
real(r_sp)                      :: fxin, fzin, dzin, dxin

optional                        :: nsrc, nrec, nt, deltat, SR_pos
optional                        :: fileid ! default = 10 [any number]

integer                         :: fileidin
integer                         :: ntin, nrecin, nsrcin
integer(i_sp)                   :: deltatin
real(r_sp), allocatable         :: SR_posin(:,:,:)
integer                         :: isrc, irec, itrc, it, ioerr, ioerr2, k

! segy headers
integer(i_dp)                                   :: tracl,tracr,fldr,tracf,ep,cdp,cdpt
integer(i_sp)                                   :: trid,nvs,nhs,duse 
integer(i_dp)                                   :: offset,gelev,selev,sdepth,gdel,sdel,swdep,gwdep
integer(i_sp)                                   :: scalel,scalco
integer(i_dp)                                   :: sx,sy,gx,gy
integer(i_sp)                                   :: counit,wevel,swevel,sut,gut,sstat,gstat, &
                                                   tstat,laga,lagb,delrt,muts,mute,ns,dt 
integer(i_sp)                                   :: gain,igc,igi,corr,sfs,sfe,slen,styp,stas, &
                                                   stae,tatyp,afilf,afils,nofilf,nofils,lcf,hcf, &
                                                   lcs,hcs,year,day,hour,minute,sec,timbas, &
                                                   trwf,grnors,grnofr,grnlof,gaps,otrav 
real(r_sp)                                      :: d1,f1,d2,f2,ungpow,unscale
integer(i_sp)                                   :: mark,unass(17)

if(present(fileid)) then ; fileidin = fileid; else; fileidin = OMP_get_thread_num(); endif
if(present(nsrc)) then ; nsrcin = nsrc; else; nsrcin = 1; endif
if(present(nrec)) then ; nrecin = nrec; else; nrecin = 1; endif
if(present(nt)) then ; ntin = nt; else; ntin = (size(dat,1) / nrecin / nsrcin); endif
if(ntin * nrecin * nsrcin .ne. size(dat,1)) call abort_msg('makesu: nt * nrec * nsrcin .ne. size of dat')
if(present(deltat)) then ; deltatin = int(deltat * 1e6, i_sp); else; deltatin = 0; endif
if(present(fx)) then ; fxin = real(fx,r_sp); else; fxin = 0.0; endif
if(present(fz)) then ; fzin = real(fz,r_sp); else; fzin = 0.0; endif
if(present(dx)) then ; dxin = real(dx,r_sp); else; dxin = 0.0; endif
if(present(dz)) then ; dzin = real(dz,r_sp); else; dzin = 0.0; endif

allocate(SR_posin(nrecin, nsrcin,5));
if(present(SR_pos)) then
        call check_dimension('makesu: size of SR_pos', SR_pos,  nrecin, nsrcin, 5)
        SR_posin = real(SR_pos, r_sp);
else
	SR_posin = 0.0
endif

! check dimensions data
call check_dimension('makesu: size of dat', dat, ntin * nrecin * nsrcin)

! default values of all the segy headers
tracl=0; tracr=0; fldr=0; tracf=0; ep=0; cdp=0; cdpt=0; trid=0;
nvs=0; nhs=0; duse=0; offset=0; gelev=0; selev=0;
sdepth=0; gdel=0; sdel=0; swdep=0; gwdep=0; scalel=0; scalco=0; sx=0;
sy=0; gx=0; gy=0; counit=0; wevel=0; swevel=0; sut=0; 
gut=0; sstat=0; gstat=0; tstat=0; laga=0
lagb=0; delrt=0; muts=0; mute=0; ns=0; dt=0; gain=0; igc=0; igi=0; corr=0; sfs=0; sfe=0
slen=0; styp=0; stas=0; stae=0; tatyp=0; afilf=0; afils=0; nofilf=0; 
nofils=0; lcf=0; hcf=0; lcs=0
hcs=0; year=0; day=0; hour=0; minute=0; sec=0; timbas=0; trwf=0; grnors=0; grnofr=0
grnlof=0; gaps=0; otrav=0; d1=0.0; f1=0.0; d2=0.0; f2=0.0; ungpow=0.0; unscale=0.0
mark=0; unass(:)=0;

! headers common for all traces
dt = deltatin ! dt should be in microseconds
ns = int(ntin, i_sp)
scalco = -100 ! positions stored in cm
scalel = -100 ! positions stored in cm
counit = 1
! non seismic headers if present
d1 = dzin;
d2 = dxin;
f1 = fzin;
f2 = fxin;

! the segy headers which are actually added
if(deltatin .eq. 0) then
        trid = 3 ! dummy meaning non-seismic data 
else
        trid = 1
endif

open(fileidin, file=trim(adjustl(ctofstr(fname))), access='direct', recl=(ntin+60)*i4, &
                status='replace', iostat = ioerr)
do isrc = 1, nsrcin
        ! source specific headers
        fldr = isrc
        do irec = 1, nrecin
                itrc = irec + (isrc - 1) * nrecin
                ! receiver specific headers
                tracf = irec
                tracl = itrc            
                tracr = itrc            
                sx = int(SR_posin(irec, isrc, 1) * abs(real(scalco,r_sp))**(-ione_sp*scalco/abs(scalco)),i_dp) ;
                selev = int(SR_posin(irec, isrc, 2) * abs(real(scalel,r_sp))**(-ione_sp*scalel/abs(scalel)),i_dp) ;
                gx = int(SR_posin(irec, isrc, 3) * abs(real(scalco,r_sp))**(-ione_sp*scalco/abs(scalco)),i_dp) ;
                gelev = int(SR_posin(irec, isrc, 4) * abs(real(scalel,r_sp))**(-ione_sp*scalel/abs(scalel)),i_dp) ;

                write(fileidin,rec=itrc, iostat = ioerr2) tracl,tracr,fldr,tracf,ep,cdp,cdpt, &
                        trid,nvs,nhs,duse,offset,gelev,selev,sdepth, &
                        gdel,sdel,swdep,gwdep,scalel,scalco, &
                        sx,sy,gx,gy,counit,wevel,swevel,sut,gut,sstat,gstat, &
                        tstat,laga,lagb,delrt,muts,mute,ns,dt, &
                        gain,igc,igi,corr,sfs,sfe,slen,styp,stas, &
                        stae,tatyp,afilf,afils,nofilf,nofils,lcf,hcf, &
                        lcs,hcs,year,day,hour,minute,sec,timbas, &
                        trwf,grnors,grnofr,grnlof,gaps,otrav, &
                        d1,f1,d2,f2,ungpow,unscale,mark,(unass(k),k=1,17), &
                        ! only single precision data can be stored
                        (real(dat(it + (itrc-1) * ntin),r_sp),it=1,ntin)
        enddo
enddo
close(fileidin)


end subroutine makesu

end module io
