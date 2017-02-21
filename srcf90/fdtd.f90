!\begin{comment}
! TODO 
! diagonal hessian for velocity and density gradients.

! USE_PML_LAYERS  (or) ABS_TRBL

!
!\end{comment}
!!
!\begin{tikzpicture}
!  \matrix (m) [matrix of math nodes,row sep=3em,column sep=4em,minimum width=2em]
!  {
!     P & P \\
!     P & P \\};
!  \path[-stealth]
!    (m-1-1) edge node [left] {${V}_z$} (m-2-1)
!            edge [] node [below] {${V}_x$} (m-1-2)
!    (m-2-1.east|-m-2-2) edge node [below] {${V}_x$}
!            node [above] {} (m-2-2)
!    (m-1-2) edge node [right] {${V}_z$} (m-2-2)
!            edge [] (m-2-1);
!\end{tikzpicture}
!
!\begin{comment}
module fdtd
use precision_mod! USEINDEPFILE
use staggered_grid! USEINDEPFILE
use error_messages! USEINDEPFILE
use fft! USEINDEPFILE
use omp_lib

type model
              real, dimension(:,:), allocatable               :: vel
               real, dimension(:,:), allocatable               :: rho
        end type model

!\end{comment}
!\subsection{Description}
! As forward modeling method, the 
! finite-difference method is employed. 
! It uses a discrete version of the two-dimensional isotropic acoustic wave equation.
! As shown in  
! Figure~\ref{fdshmesh_acou}, a staggered-grid mesh is used 
! with particle velocity, $\vy$, density, $\mR$ and wave-speed, $\mV$,
! on the main grid. The stress components, $\sigmayx$ and $\sigmayz$,  
! are on the secondary grids.
! As part of the input, material parameters need to be given. 
! The spatial grid size is 
! important to make sure the data 
! are computed with no numerical dispersion.
! The time sampling interval is important to prevent the code 
! from being unstable.
! We used Convolutional Perfectly Matched
! Layer (C-PML) absorbing conditions \citep{KoMa07,RoGe00}.
! Fourth-order operators are used in space and second-order operators are 
! used in time. 
! The absorbing conditions are adapted from SEISMIC\_CPML, 
! a set of eleven open-source FORTRAN90 programs\footnote
! {\url{https://geodynamics.org/cig/software/seismic_cpml/}}.
! The partial differential equations governing the SH wave propagation 
! are:

! $$ \mM = \mR \mV^2$$
! %\mL = \mR * (cp**2 - 2 * cs**2)

! $$\mR \frac{\partial\vy}{\partial t} = \frac{\partial\sigmayx}{\partial x} 
! + \frac{\partial \sigmayz}{\partial z}  + \sfo$$
! $$\frac{\partial \sigmayx}{\partial t} = \mM \frac{\partial \vy}{\partial x}$$
! $$\frac{\partial \sigmayz}{\partial t} = \mM \frac{\partial \vy}{\partial z}$$

!\begin{figure}[hb]
!\begin{center}
!\begin{tikzpicture}[auto, thick]
!\definecolor{lavander}{cmyk}{0,0.48,0,0}
!\definecolor{violet}{cmyk}{0.79,0.88,0,0}
!\definecolor{burntorange}{cmyk}{0,0.52,1,0}
!\definecolor{amber}{rgb}{1.0, 0.75, 0.0}
!\definecolor{asparagus}{rgb}{0.53, 0.66, 0.42}
!\definecolor{babyblue}{rgb}{0.54, 0.81, 0.94}
!
!
!\tikzstyle{vx}=[draw,circle,black,bottom color=asparagus,
!                  top color= white, text=black,minimum width=10pt]
!\tikzstyle{vz}=[draw,circle,black,bottom color=white,
!                  top color= babyblue, text=black,minimum width=10pt]
!\tikzstyle{pressure}=[draw,circle,black, left color=amber,
!                       text=black,minimum width=14pt]
!\tikzstyle{legendsp}=[rectangle, draw, black, rounded corners,
!                     thin,left color=amber, 
!                     text=black, minimum width=2.5cm, minimum
!height=0.5cm]
!\tikzstyle{legendp}=[rectangle, draw, black, rounded corners, thin,
!                     bottom color=white, top color= babyblue,
!                     text= black, minimum width= 2.5cm, minimum height=0.5cm]
!\tikzstyle{legend_general}=[rectangle, rounded corners, thin,
!                           black, top color= white, bottom color=asparagus, draw,
!text=black,
!                           minimum width=2.5cm, minimum height=0.5cm]
!
!  % Place super vx and connect them
!  \foreach \place/\name in {{(0,1)/b}, {(0,3)/c}, {(0,5)/d}}
!    \node[pressure] (\name) at \place {};
!   %
!   % Place normal vx
!  \foreach \pos/\i in {right of/1}
!    \node[vx, \pos = b] (b\i) {};
!  \foreach \pos/\i in {right of/1}
!    \node[vx, \pos = c] (c\i) {};
!  \foreach \pos/\i in {right of/1}
! \node[vx, \pos = d] (d\i) {};
!  \foreach \pos/\i in {right of/1};
!
!  \foreach \pos/\i in {below of/1}
!    \node[vz, \pos = d] (d\i) {};
! \foreach \pos/\i in {below of/1}
!    \node[vz, \pos = c] (c\i) {};
! \foreach \pos/\i in {below of/1}
!    \node[vz, \pos = b] (b\i) {};
!
!
!  \foreach \place/\name in {{(2,1)/f}, {(2,3)/g}, {(2,5)/h}}
!    \node[pressure] (\name) at \place {};
!
!   %
!   % Place normal vx
!  \foreach \pos/\i in {right of/1}
!    \node[vx, \pos = f] (f\i) {};
!  \foreach \pos/\i in {right of/1}
!    \node[vx, \pos = g] (g\i) {};
!  \foreach \pos/\i in {right of/1}
! \node[vx, \pos = h] (h\i) {};
!  \foreach \pos/\i in {right of/1};
!
!  \foreach \pos/\i in {below of/1}
!    \node[vz, \pos = h] (h\i) {};
! \foreach \pos/\i in {below of/1}
!    \node[vz, \pos = g] (g\i) {};
! \foreach \pos/\i in {below of/1}
!    \node[vz, \pos = f] (f\i) {};
!
! \foreach \place/\name in {{(4,1)/f}, {(4,3)/g}, {(4,5)/h}}
!    \node[pressure] (\name) at \place {};
!  \foreach \pos/\i in {right of/1}
!    \node[vx, \pos = f] (f\i) {};
!  \foreach \pos/\i in {right of/1}
!    \node[vx, \pos = g] (g\i) {};
!  \foreach \pos/\i in {right of/1}
! \node[vx, \pos = h] (h\i) {};
!  \foreach \pos/\i in {right of/1};
!  \foreach \pos/\i in {below of/1}
!    \node[vz, \pos = h] (h\i) {};
! \foreach \pos/\i in {below of/1}
!    \node[vz, \pos = g] (g\i) {};
! \foreach \pos/\i in {below of/1}
!    \node[vz, \pos = f] (f\i) {};
!   %\foreach \speer/\peer in {e/e1,e/e2,e/e3}
!   % \path (\speer) edge (\peer);
!   %
!   %%%%%%%%
!   % Legends
!   \node[legend_general] at (8,2) {\small{\sigmayx}};
!   \node[legendp] at (8,3) {\small{\sigmayz}};
!   \node[legendsp] at (8,4) {\small{\vy, \mR, \mV}};
!\end{tikzpicture}
!\caption{Staggered-grid finite-difference modeling mesh}
!\label{fdshmesh}
!\end{center}
!\end{figure}


!\end{comment}
!\begin{comment}
! to solve the two-dimensional acoutic wave equation
! using a finite-difference method with Convolutional Perfectly Matched
! Layer (C-PML) conditions.

! 2D acoustic finite-difference code in velocity and stress formulation
! with Convolutional-PML (C-PML) absorbing conditions for an isotropic medium

! C-PML adapted from : Dimitri Komatitsch, University of Pau, France, April 2007.
! Fourth-order implementation by Dimitri Komatitsch and Roland Martin, University of Pau, France, August 2007.


! Author           : Pawan Bharadwaj
!                    bharadwaj.pawan@gmail.com

! this module solves three different cases 
! Acoustic 2D  ...

! rho = density
! mu, lambda = lame"
! cp, cs = P and S wave velocities

! mu = rho * cs * cs = zero
! lambda = rho * (cp**2 - 2 * cs**2)

! their description given below now ;)

! rho * dvx_dt = dp_dx 
! rho * dvz_dt = dp_dz 
! dp_dt = (lambda ) * dvx_dx + lambda * dvz_dz + lambda * integral(pressure_source))

! pressure source is added after an integraion (dividing with iomega in freq domain) in order not to 
! input injection rate as we are solving first order systems ..

! main grid  --- pressure is the main field that is recorded
! and all the parameters (lambda and density) are on the pressure grid

! The staggered-grid formulation of Madariaga (1976) and Virieux (1986) is used:
!
!            ^ z
!            |
!            |
!
!            +-------------------+
!            |                   |
!            |                   |
!            |                   |
!            |                   |
!            |        v_z        |
!            +---------+         |
!            |         |         |
!            |         |         |
!            |         |         |
!            |         |         |
!            |         |         |
!            +---------+---------+  ---> x
!           v_x        p        v_x
!                    lambda
!                     rho  
!
! but a fourth-order spatial operator is used
!
! The C-PML implementation is based in part on formulas given in Roden and Gedney (2000)
! total number of grid points in each direction of the grid

!\begin{comment}

! thickness of the PML layer in grid points (typically 10 or 20 points are sufficient)

contains

subroutine fdtd_mod(&
        jobname, & ! optional argument used as prefix while saving ()
        nprop, & ! number of independent wavefields propagating independently in the same medium
! >>>>>>> model parameters
        modTT, & ! model parameter matrix before second order time derivative [mesh_nz][mesh_nx]
        modRR, & ! model parameter vector before  [mesh_nz][mesh_nx]
        mesh_X, & ! modeling mesh X vector
        mesh_Z, & ! modeling mesh Z vector
        mesh_dx, & ! 
        mesh_dz, & ! 
        mesh_nx, & !
        mesh_nz, & !
        na_pml, & ! number of PML layers to be added
        abs_trbl, & ! absorbing boundary conditions Top, Right, Bottom, Left
! >>>>>>> time parameters
        ntmod,& ! total time samples
        dt, & ! time sampling
! >>>>>>> source parameters
        ntsrc, & ! number time samples in source wavelet (ntsrc <= ntmod)
        src_wavelets, & ! source wavelets [ntsrc] * [src_nsmul] * [src_nseq]
        src_nseq, & ! number of sequential sources (OMP loop)
        src_nsmul, & ! number of simultaneous sources
        src_x, & ! x coordinates of all the sources [src_nsmul] * [src_nseq]
        src_z, & ! z coordinated of all the sources [src_nsmul] * [src_nseq]
        source_flag, & ! parameters determining source 
                        ! contains [IRATE], then injection rate is ON
                        ! contains [TIME_REVERSE], then inject time reverse of the src_wavelets
                        ! contains [BILINEAR], then bilinear interpolation
                        ! contains [NEAREST_NEIGH], then nearest neighbour interpolation
! >>>>>>> receiver parameters
        recv_n, & ! number of receivers for each source 
        recv_x, & ! x coord for receivers [recv_n] * [src_nseq]
        recv_z, & ! z  "                  [recv_n] * [src_nseq]
        records, & ! recorded pressure at receivers [ntmod] * [recv_n] * [src_nseq]
        records_flag, & ! '[P]' records contain only pressure
                        ! '[P][Vz]' records contain pressure and Vz (need to be implemented)  
                        ! '[DLDVEL_P]' 
                        ! 
! >>>>>>> snapsshots and illumination and gradients
        snaps_in, & ! input snaps shots in the order [mesh_nz,mesh_nx,ntmod,src_nseq] 
        snaps_out, & ! output snaps shots in the order [mesh_nz,mesh_nx,ntmod,src_nseq]
        grad,  &  ! output gradient models w.r.t velocity and density  [src_nseq]      
        varpert, & ! background model, optional for testing only
        snap_save_at& ! a snap shot at this time is saved 
        ) !bind(c, name="fdtd_mod")
! Author: Pawan Bharadwaj
!         p.b.pisupati@tudelft.nl
!         modified: 11 Sept 2013
!         updated on 25 July 2014
                ! parallelized w.r.t src_nseq
                ! added ntsrc so that source wavelets can have different
                ! dimensions 

        implicit none
character(len=5), intent(in)                            :: jobname
real, intent(in)                                        :: dt, snap_save_at
integer, intent(in)                                     :: nprop
integer, intent(in)                                     :: ntmod, ntsrc, src_nsmul, src_nseq, recv_n
type(model), intent(in)                                 :: var, varpert
optional                                                :: varpert
type(model), intent(inout)                              :: grad(:)
real, dimension(:), intent(in)                          :: src_wavelets, src_x, src_z, recv_x, recv_z

integer, intent(in)                                     :: mesh_nx, mesh_nz
real, intent(in)                                        :: mesh_X(mesh_nx), mesh_Z(mesh_nz), mesh_dx, mesh_dz
real, dimension(:), intent(out)                         :: records
character(len=*), intent(in)                            :: records_flag
integer, intent(in)                                     :: na_pml  !number of PML layers
character(len=4), intent(in)                            :: abs_trbl

real, dimension(:,:,:,:), intent(in)                    :: snaps_in
real, dimension(:,:,:,:), intent(out)                   :: snaps_out
character(len=*), intent(in)                            :: source_flag

! which are optional
optional                                                :: jobname, snap_save_at,  snaps_in, snaps_out
optional                                                :: records, records_flag, recv_n, recv_x, recv_z
optional                                                :: grad

! other variables for this only this subroutine
! total number of grid points in each direction of the grid
integer                                                 :: nx
integer                                                 :: nz

! size of a grid cell
real                                                    :: deltax, deltax24I 
real                                                    :: deltaz, deltaz24I
real                                                    :: deltat, deltatI 

integer                                                 :: ix,iz,it,ixx,izz, irec, issmul, isseq
real, dimension(:,:,:), allocatable                     :: src_wav_mat
real                                                    :: source_term

! source position variables                            
integer, allocatable, dimension(:)                      :: issmulx1, issmulx2, issmulz1, issmulz2
integer                                                 :: issmulxtemp1, issmulxtemp2, issmulztemp1, issmulztemp2
real, allocatable, dimension(:)                         :: denomsrcI
real, allocatable, dimension(:,:)                       :: src_spray_weights
real                                                    :: svalue
! receiver position variables           
integer, allocatable, dimension(:)                      :: irecx1, irecx2, irecz1, irecz2
integer                                                 :: irecxtemp1, irecxtemp2, irecztemp1, irecztemp2
real, dimension(:), allocatable                         :: denomrecI
real, dimension(:,:), allocatable                       :: rec_interpolate_weights

! recording data
real, dimension(:,:), allocatable                       :: rec_mat(:,:)

! to find indices
logical, allocatable                                    :: maskX(:), maskZ(:)

! extended models 
real, dimension(:,:), allocatable                       :: lambda, rho, rhovxI, rhovzI, rhoI

! models for background velocity                        
real, dimension(:,:), allocatable                       :: rhoIpert, rhovIpert

! for boundary conditions on pressure grid .. arrays either one or zero
real, dimension(:,:), allocatable                       :: boundary_p, boundary_vx, boundary_vz

! minimum and maximum frequency                                
real                                                    :: freqmin, freqmax

! minimum and maximum velocities                                
real                                                    :: velmin, velmax

! for stability
real                                                    :: Courant_number

! power to compute d0 profile
real, parameter                                         :: NPOWER = 2.e0

real, parameter                                         :: K_MAX_PML = 1.e0 ! from Gedney page 8.11
real                                                    :: ALPHA_MAX_PML  ! from Festa and Vilotte
! 1D arrays for the damping profiles 
real, dimension(:),allocatable                          :: d_x,K_x,K_xI,alpha_x,a_x,b_x,d_x_half,&
                                                           K_x_half,K_x_halfI,alpha_x_half,a_x_half,b_x_half
real, dimension(:),allocatable                          :: d_z,K_z,K_zI,alpha_z,a_z,b_z,d_z_half,&
                                                           K_z_half,K_z_halfI,alpha_z_half,a_z_half,b_z_half

real                                                    :: thickness_PML_x,thickness_PML_z,&
                                                           xoriginleft,xoriginright,&
                                                           zoriginbottom,zorigintop
real                                                    :: Rcoef,d0_x,d0_z,&
                                                             xval,zval,abscissa_in_PML,&
                                                             abscissa_normalized


! main arrays while wave propagation
real, dimension(:,:), allocatable                       :: vx,vz,p
real, dimension(:,:), allocatable                       :: pp ! for time differential of pressure
real, dimension(:,:), allocatable                       :: dpFdx,dpFdz,dpdx,dpdz ! spatial derivatives of pressure in vx and vz grids
real, dimension(:,:), allocatable                       :: dvxdx,dvzdz

real, dimension(:,:), allocatable                       :: grad_modTT
real, dimension(:,:), allocatable                       :: grad_rhovxI, grad_rhovzI, grad_modRR

real, allocatable                                       :: saved_snap(:,:)

real(r_dp)                                              :: cpu_time1, cpu_time2

! arrays for the memory variables
! could declare these arrays in PML only to save a lot of memory, but proof of concept only here
real, dimension(:,:), allocatable                       :: memory_dvx_dx, &
                                                                   memory_dvx_dz, &
                                                                   memory_dvz_dx, &
                                                                   memory_dvz_dz, &
                                                                   memory_dp_dx, &
                                                                   memory_dp_dz
! flags to add PML layers to the edges of the grid
logical                                                 :: USE_PML_XMIN
logical                                                 :: USE_PML_XMAX
logical                                                 :: USE_PML_ZMIN
logical                                                 :: USE_PML_ZMAX

real, parameter                                         :: pi=3.141592654


cpu_time1  = omp_get_wtime();

! dimension checks
call check_dimension("fd2_mod:",var%vel,mesh_nz,mesh_nx)
call check_dimension("fd2_mod:",var%rho,mesh_nz,mesh_nx)
if(ntmod .lt. ntsrc) call abort_msg("fd2_mod: ntmod .lt. ntsrc")
call check_dimension("fd2_mod: size src_wavelets", src_wavelets, ntsrc*src_nsmul*src_nseq)
call check_dimension("fd2_mod: size src_x", src_x, src_nsmul*src_nseq)
call check_dimension("fd2_mod: size src_z", src_z, src_nsmul*src_nseq)
if(.not.allocated(var%rho)) then
        call abort_msg("fd2_mod: input model doesnt have rho")
endif
if(.not.allocated(var%vel)) then
        call abort_msg("fd2_mod: input model type model doesnt have vel")
endif
if(present(grad)) then
        if(size(grad,1) .ne. src_nseq) call abort_msg("size(grad) .ne. src_nseq")
endif
if(present(recv_n)) then
        if(.not.(present(recv_x))) call abort_msg("fd2_mod: need recv_x")
        if(.not.(present(recv_z))) call abort_msg("fd2_mod: need recv_z")
endif
if(present(recv_x)) then
        call check_dimension("fd2_mod: size recv_x", recv_x, recv_n*src_nseq)
endif
if(present(recv_z)) then
        call check_dimension("fd2_mod: size recv_z", recv_z, recv_n*src_nseq)
endif
if(present(records)) then
        if(.not.(present(recv_n))) call abort_msg("fd2_mod: need recv_n")
        if(.not.(present(records_flag))) call abort_msg("fd2_mod: need records_flag")
        call check_dimension("fd2_mod: size records", records, ntmod*recv_n*src_nseq)
        records = rzero_de
endif
if(present(snaps_out)) then
        call check_dimension("fd2_mod: size snaps_out", snaps_out,  mesh_nz,mesh_nx,ntmod, src_nseq)
        snaps_out = rzero_de
endif
if(present(snaps_in)) then
        call check_dimension("fd2_mod: size snaps_in", snaps_in, mesh_nz, mesh_nx,ntmod, src_nseq)
endif

! no modeling if source wavelet is zero
if(maxval(abs(src_wavelets)) .lt. tiny(rzero_de)) then
        return
endif
! convert src_wavelets to src_wav_mat
allocate(src_wav_mat(src_nsmul, 0:ntmod+1, src_nseq))
src_wav_mat = rzero_de;

do isseq = 1, src_nseq
        do issmul = 1, src_nsmul
                source_term = rzero_de
                do it = 1, ntsrc-1
                        if(index(source_flag,'[TIME_REVERSE]') .ne. 0) then
                                source_term = &
                                src_wavelets(ntsrc - it  + (issmul-1)*ntsrc + (isseq-1)*src_nsmul * ntsrc)
                        else
                                source_term = &
                                src_wavelets(it + (issmul-1)*ntsrc + (isseq-1)*src_nsmul * ntsrc)
                        endif
                        if(index(source_flag,'[IRATE]') .eq. 0) then
                                src_wav_mat(issmul, it + 1, isseq) = 2.0 * source_term + &
                                        src_wav_mat(issmul, it - 1, isseq)
                        else
                                src_wav_mat(issmul, it, isseq) = source_term
                        endif
                enddo
        enddo
enddo
if(maxval(src_wav_mat) .ge. huge(1.0e0)) call abort_msg("fd2_mod: huge src_wav_mat?")

! allocate masks
allocate(maskX(mesh_nx)); maskX = .true.
allocate(maskZ(mesh_nz)); maskZ = .true.

if(na_pml .le. 0) call abort_msg("fd2_mod: invalid na_pml")

! extended dimensions of the modeling grid
nx    = mesh_nx+2*na_pml
nz    = mesh_nz+2*na_pml

! assign variables
deltax=mesh_dx
deltax24I = (deltax * 24.e0)**(-1.e0)
deltaz=mesh_dz
deltaz24I = (deltaz * 24.e0)**(-1.e0)
deltat=dt
deltatI= (deltat)**(-1.e0)

! minimum and maximum velocities
velmin = minval(var%vel); velmax = maxval(var%vel)

! minimum and maximum frequencies
freqmin = huge(1.e0); freqmax = 0.0;
do isseq = 1, src_nseq
        do issmul = 1, src_nsmul
                if(maxval(abs(src_wav_mat(issmul,:,isseq))) .gt. tiny(rzero_de)) then
                        freqmin = &
                        min(freqmin,findfreq(x=pack(src_wav_mat(issmul,:,isseq),.true.),dt=dt,str="MIN",threshold=2e-2))
                        freqmax = &
                        max(freqmax,findfreq(x=pack(src_wav_mat(issmul,:,isseq),.true.),dt=dt,str="MAX",threshold=2e-2))
                endif
        enddo
enddo
ALPHA_MAX_PML = 2.e0*PI*((freqmin + freqmax))/4.e0 ! from Festa and Vilotte

! (type model var) to (rho and lambda), which are global variables
allocate(lambda(0:nz+1,0:nx+1),rho(0:nz+1,0:nx+1))
allocate(rhovxI(0:nz+1,0:nx+1),rhovzI(0:nz+1,0:nx+1),rhoI(0:nz+1,0:nx+1))
lambda(:,:)=rzero_de; rho(:,:)= rzero_de; rhoI = rzero_de; rhovxI(:,:) = rzero_de; rhovzI(:,:) = rzero_de
lambda = extend((var%vel * var%vel * var%rho), na_pml, mesh_nx, mesh_nz)
rho = extend((var%rho), na_pml, mesh_nx, mesh_nz)

! rhoI on p frid
call get_rhoI(rho, rhoI)

! rho on vx grid
call get_rhovxI(rhoI, rhovxI)

! rho on vz grid
call get_rhovzI(rhoI, rhovzI)


call is_nan("fd2_mod: lambda", lambda)
call is_nan("fd2_mod: rho", rho)
call is_nan("fd2_mod: rhoI", rhoI)
call is_nan("fd2_mod: rhovzI", rhovzI)
call is_nan("fd2_mod: rhovxI", rhovxI)

! boundary conditions on p
allocate(boundary_p(0:nz+1,0:nx+1)); allocate(boundary_vx(0:nz+1,0:nx+1));
allocate(boundary_vz(0:nz+1,0:nx+1))
! default
boundary_p = rone_de; boundary_vx = rone_de; boundary_vz = rone_de;
USE_PML_ZMIN = .true.; USE_PML_XMAX = .true.; USE_PML_ZMAX = .true.; USE_PML_XMIN = .true.

! free surface boundary conditions on sides of the model.
if(abs_trbl(1:1) .eq. "0") then
        boundary_p(1:na_pml+1,:) = rzero_de
        !boundary_vz(1:na_pml,:) = rzero_de
!        boundary_vx(1:na_pml,:) = rzero_de
!        USE_PML_ZMIN = .false.
endif
if(abs_trbl(2:2) .eq. "0") then
        boundary_p(1:nz,na_pml+mesh_nx:nx) = rzero_de
        USE_PML_XMAX = .false.
endif
if(abs_trbl(3:3) .eq. "0") then 
        boundary_p(na_pml+mesh_nz,1:nx) = rzero_de
        USE_PML_ZMAX = .false.
endif
if(abs_trbl(4:4) .eq. "0") then
        boundary_p(1:nz,na_pml+1) = rzero_de
        USE_PML_XMIN = .false.
endif


! Dirichlet conditions (rigid boundaries) on the edges or at the bottom of the PML layers
boundary_vx(1,:) = rzero_de
boundary_vx(nz,:) = rzero_de
boundary_vx(:,1) = rzero_de
boundary_vx(:,nx) = rzero_de
boundary_vz(1,:) = rzero_de
boundary_vz(nz,:) = rzero_de
boundary_vz(:,1) = rzero_de
boundary_vz(:,nx) = rzero_de


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! for damping profiles
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allocate         (d_x(nx),K_x(nx),K_xI(nx),alpha_x(nx),&
                          a_x(nx),b_x(nx),d_x_half(nx),&
                          K_x_half(nx),K_x_halfI(nx),alpha_x_half(nx),a_x_half(nx),b_x_half(nx))

allocate         (d_z(nz),K_z(nz),K_zI(nz),alpha_z(nz),&
                          a_z(nz),b_z(nz),d_z_half(nz),&
                          K_z_half(nz),K_z_halfI(nz),alpha_z_half(nz),a_z_half(nz),b_z_half(nz))
! thickness of the PML layer in meters
thickness_PML_x = na_pml * deltax
thickness_PML_z = na_pml * deltaz

! reflection coefficient (INRIA report section 6.1) http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
Rcoef = 0.001e0

! check that NPOWER is okay
if(NPOWER < 1) stop "NPOWER must be greater than 1"

! compute d0 from INRIA report section 6.1 http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
d0_x = - (NPOWER + 1) * (velmax+velmin)/2.0 * log(Rcoef) / (2.e0 * thickness_PML_x)
d0_z = - (NPOWER + 1) * (velmax+velmin)/2.0 * log(Rcoef) / (2.e0 * thickness_PML_z)

d_x(:) = rzero_de
d_x_half(:) = rzero_de
K_x(:) = 1.e0
K_xI(:) = 1.e0
K_x_half(:) = 1.e0
K_x_halfI(:) = 1.e0
alpha_x(:) = rzero_de
alpha_x_half(:) = rzero_de
a_x(:) = rzero_de
a_x_half(:) = rzero_de

d_z(:) = rzero_de
d_z_half(:) = rzero_de
K_z(:) = 1.e0
K_zI(:) = 1.e0
K_z_half(:) = 1.e0
K_z_halfI(:) = 1.e0
alpha_z(:) = rzero_de
alpha_z_half(:) = rzero_de
a_z(:) = rzero_de
a_z_half(:) = rzero_de

! damping in the X direction
! origin of the PML layer (position of right edge minus thickness, in meters)
xoriginleft = thickness_PML_x
xoriginright = (nx-1)*deltax - thickness_PML_x

do ix = 1,nx
        ! abscissa of current grid point along the damping profile
        xval = deltax * real(ix-1)

        !---------- left edge
        if(USE_PML_XMIN) then

                ! define damping profile at the grid points
                abscissa_in_PML = xoriginleft - xval
                if(abscissa_in_PML >= rzero_de) then
                abscissa_normalized = abscissa_in_PML / thickness_PML_x
                d_x(ix) = d0_x * abscissa_normalized**NPOWER
                ! this taken from Gedney page 8.2
                K_x(ix) = 1.e0 + (K_MAX_PML - 1.e0) * abscissa_normalized**NPOWER
                alpha_x(ix) = ALPHA_MAX_PML * (1.e0 - abscissa_normalized) + 0.1e0 * ALPHA_MAX_PML
                endif

                ! define damping profile at half the grid points
                abscissa_in_PML = xoriginleft - (xval + deltax/2.e0)
                if(abscissa_in_PML >= rzero_de) then
                abscissa_normalized = abscissa_in_PML / thickness_PML_x
                d_x_half(ix) = d0_x * abscissa_normalized**NPOWER
                ! this taken from Gedney page 8.2
                K_x_half(ix) = 1.e0 + (K_MAX_PML - 1.e0) * abscissa_normalized**NPOWER
                alpha_x_half(ix) = ALPHA_MAX_PML * (1.e0 - abscissa_normalized) + 0.1e0 * ALPHA_MAX_PML
                endif

        endif

        !---------- right edge
        if(USE_PML_XMAX) then

                ! define damping profile at the grid points
                abscissa_in_PML = xval - xoriginright
                if(abscissa_in_PML >= rzero_de) then
                abscissa_normalized = abscissa_in_PML / thickness_PML_x
                d_x(ix) = d0_x * abscissa_normalized**NPOWER
                ! this taken from Gedney page 8.2
                K_x(ix) = 1.e0 + (K_MAX_PML - 1.e0) * abscissa_normalized**NPOWER
                alpha_x(ix) = ALPHA_MAX_PML * (1.e0 - abscissa_normalized) + 0.1e0 * ALPHA_MAX_PML
                endif

                ! define damping profile at half the grid points
                abscissa_in_PML = xval + deltax/2.e0 - xoriginright
                if(abscissa_in_PML >= rzero_de) then
                abscissa_normalized = abscissa_in_PML / thickness_PML_x
                d_x_half(ix) = d0_x * abscissa_normalized**NPOWER
                ! this taken from Gedney page 8.2
                K_x_half(ix) = 1.e0 + (K_MAX_PML - 1.e0) * abscissa_normalized**NPOWER
                alpha_x_half(ix) = ALPHA_MAX_PML * (1.e0 - abscissa_normalized) + 0.1e0 * ALPHA_MAX_PML
                endif

        endif

        ! just in case, for -5 at the end
        if(alpha_x(ix) < rzero_de) alpha_x(ix) = rzero_de
        if(alpha_x_half(ix) < rzero_de) alpha_x_half(ix) = rzero_de

        b_x(ix) = exp(- (d_x(ix) / K_x(ix) + alpha_x(ix)) * deltat)
        b_x_half(ix) = exp(- (d_x_half(ix) / K_x_half(ix) + alpha_x_half(ix)) * deltat)

        ! this to avoid division by zero outside the PML
        if(abs(d_x(ix)) > 1.e-6) a_x(ix) = d_x(ix) * (b_x(ix) - 1.e0) / (K_x(ix) * (d_x(ix) + K_x(ix) * alpha_x(ix)))
        if(abs(d_x_half(ix)) > 1.e-6) a_x_half(ix) = d_x_half(ix) * &
        (b_x_half(ix) - 1.e0) / (K_x_half(ix) * (d_x_half(ix) + K_x_half(ix) * alpha_x_half(ix)))

enddo
K_xI = K_x**(-1.e0)
K_x_halfI = K_x_half**(-1.e0)

! damping in the Y direction
! origin of the PML layer (position of right edge minus thickness, in meters)
zoriginbottom = thickness_PML_z
zorigintop = nz*deltaz - thickness_PML_z

do iz = 1,nz

        ! abscissa of current grid point along the damping profile
        zval = deltaz * real(iz-1)

        !---------- bottom edge
        if(USE_PML_ZMIN) then

                ! define damping profile at the grid points
                abscissa_in_PML = zoriginbottom - zval
                if(abscissa_in_PML >= rzero_de) then
                abscissa_normalized = abscissa_in_PML / thickness_PML_z
                d_z(iz) = d0_z * abscissa_normalized**NPOWER
                ! this taken from Gedney page 8.2
                K_z(iz) = 1.e0 + (K_MAX_PML - 1.e0) * abscissa_normalized**NPOWER
                alpha_z(iz) = ALPHA_MAX_PML * (1.e0 - abscissa_normalized) + 0.1e0 * ALPHA_MAX_PML
                endif

                ! define damping profile at half the grid points
                abscissa_in_PML = zoriginbottom - (zval + deltaz/2.e0)
                if(abscissa_in_PML >= rzero_de) then
                abscissa_normalized = abscissa_in_PML / thickness_PML_z
                d_z_half(iz) = d0_z * abscissa_normalized**NPOWER
                ! this taken from Gedney page 8.2
                K_z_half(iz) = 1.e0 + (K_MAX_PML - 1.e0) * abscissa_normalized**NPOWER
                alpha_z_half(iz) = ALPHA_MAX_PML * (1.e0 - abscissa_normalized) + 0.1e0 * ALPHA_MAX_PML
                endif

        endif

        !---------- top edge
        if(USE_PML_ZMAX) then

                ! define damping profile at the grid points
                abscissa_in_PML = zval - zorigintop
                if(abscissa_in_PML >= rzero_de) then
                abscissa_normalized = abscissa_in_PML / thickness_PML_z
                d_z(iz) = d0_z * abscissa_normalized**NPOWER
                ! this taken from Gedney page 8.2
                K_z(iz) = 1.e0 + (K_MAX_PML - 1.e0) * abscissa_normalized**NPOWER
                alpha_z(iz) = ALPHA_MAX_PML * (1.e0 - abscissa_normalized) + 0.1e0 * ALPHA_MAX_PML
                endif

                ! define damping profile at half the grid points
                abscissa_in_PML = zval + deltaz/2.e0 - zorigintop
                if(abscissa_in_PML >= rzero_de) then
                abscissa_normalized = abscissa_in_PML / thickness_PML_z
                d_z_half(iz) = d0_z * abscissa_normalized**NPOWER
                ! this taken from Gedney page 8.2
                K_z_half(iz) = 1.e0 + (K_MAX_PML - 1.e0) * abscissa_normalized**NPOWER
                alpha_z_half(iz) = ALPHA_MAX_PML * (1.e0 - abscissa_normalized) + 0.1e0 * ALPHA_MAX_PML
                endif

        endif

        b_z(iz) = exp(- (d_z(iz) / K_z(iz) + alpha_z(iz)) * deltat)
        b_z_half(iz) = exp(- (d_z_half(iz) / K_z_half(iz) + alpha_z_half(iz)) * deltat)

        ! this to avoid division by zero outside the PML
        if(abs(d_z(iz)) > 1.e-6) a_z(iz) = d_z(iz) * (b_z(iz) - 1.e0) / (K_z(iz) * (d_z(iz) + K_z(iz) * alpha_z(iz)))
        if(abs(d_z_half(iz)) > 1.e-6) a_z_half(iz) = d_z_half(iz) * &
        (b_z_half(iz) - 1.e0) / (K_z_half(iz) * (d_z_half(iz) + K_z_half(iz) * alpha_z_half(iz)))

enddo
K_zI = K_z**(-1.e0)
K_z_halfI = K_z_halfI**(-1.e0)

! check the Courant stability condition for the explicit time scheme
! R. Courant et K. O. Friedrichs et H. Lewy (1928)
courant_number = velmin * deltat * sqrt(1.e0/deltax**2 + 1.e0/deltaz**2)
if(courant_number > 1.e0) call abort_msg("fd2_mod: time step is too large, simulation will be unstable")


!$OMP PARALLEL NUM_THREADS(OMP_GET_MAX_THREADS()) DEFAULT(NONE) &
!$OMP PRIVATE(isseq, it, issmul, ixx, izz, maskX, maskZ, denomsrcI, denomrecI, svalue, &
!$OMP issmulxtemp1,issmulxtemp2,issmulztemp1,issmulztemp2,issmulx1,issmulz1,issmulx2,issmulz2,src_spray_weights, &
!$OMP irecxtemp1,irecxtemp2,irecztemp1,irecztemp2,irecx1,irecz1,irecx2,irecz2,rec_interpolate_weights, &
!$OMP pp,dpFdx,dpFdz,dpdx,dpdz,dvxdx,dvzdz,grad_modTT,grad_rhovxI,grad_rhovzI,grad_modRR, &
!$OMP vx,vz,p, &
!$OMP memory_dvx_dx, &
!$OMP memory_dvx_dz, &
!$OMP memory_dvz_dx, &
!$OMP memory_dvz_dz, &
!$OMP memory_dp_dx, &
!$OMP memory_dp_dz, &
!$OMP saved_snap, rec_mat) &
!$OMP SHARED(src_nseq,src_wav_mat,source_flag, var, grad, mesh_nz, mesh_nx, nz, nx, snaps_in, &
!$OMP mesh_X, mesh_Z, &
!$OMP snaps_out, ntmod, ntsrc, src_nsmul, na_pml, deltax, deltaz, deltax24I, deltaz24I, deltat, deltatI, & 
!$OMP src_x, src_z, recv_n, recv_z, recv_x, snap_save_at, records, records_flag, &
!$OMP jobname, lambda, rho,rhoI,rhovxI,rhovzI, &
!$OMP b_x_half, a_x_half, b_z_half, a_z_half, b_x, b_z, a_x, a_z, &
!$OMP K_x, K_z, K_xI, K_zI, varpert, rhoIpert, rhovIpert, &
!$OMP K_x_half, K_z_half, K_x_halfI, K_z_halfI, &
!$OMP boundary_p, boundary_vx, boundary_vz) 

!$OMP DO 
src_par_loop: do isseq = 1, src_nseq

        allocate         (p(0:nz+1,0:nx+1), pp(0:nz+1,0:nx+1), &
                          vx(0:nz+1,0:nx+1),vz(0:nz+1,0:nx+1))
        allocate         (dpdx(0:nz+1,0:nx+1), dpdz(0:nz+1,0:nx+1))
        allocate         (dvxdx(0:nz+1,0:nx+1), dvzdz(0:nz+1,0:nx+1))

        allocate         (memory_dvx_dx(0:nz+1,0:nx+1), &
                                 memory_dvx_dz(0:nz+1,0:nx+1), &
                                 memory_dvz_dx(0:nz+1,0:nx+1), &
                                 memory_dvz_dz(0:nz+1,0:nx+1), &
                                 memory_dp_dx(0:nz+1,0:nx+1), &
                                 memory_dp_dz(0:nz+1,0:nx+1))

        ! initialize main arrays
        p(:,:) = rzero_de; pp(:,:) = rzero_de;  vx(:,:) = rzero_de;  vz(:,:) = rzero_de
        dpdx = rzero_de; dpdz = rzero_de;
        dvxdx = rzero_de; dvzdz = rzero_de;

        ! temp array per thread to store records
        if(present(recv_n)) then
                allocate   (rec_mat(recv_n,ntmod))
                rec_mat = rzero_de
        endif

        ! initialize PML arrays
        memory_dp_dx(:,:) = rzero_de
        memory_dp_dz(:,:) = rzero_de
        memory_dvx_dx(:,:) = rzero_de
        memory_dvx_dz(:,:) = rzero_de
        memory_dvz_dx(:,:) = rzero_de
        memory_dvz_dz(:,:) = rzero_de

        if(present(grad)) then
                allocate(grad_rhovxI(mesh_nz, mesh_nx), grad_rhovzI(mesh_nz, mesh_nx), grad_modTT(mesh_nz,mesh_nx))
                grad_rhovxI =  rzero_de; grad_rhovzI = rzero_de; grad_modTT = rzero_de;
                allocate(grad_modRR(mesh_nz, mesh_nx)); grad_modRR = rzero_de;
                allocate(dpFdx(mesh_nz, mesh_nx),dpFdz(mesh_nz, mesh_nx)) 
                dpFdx(:,:) = rzero_de ; dpFdz(:,:) = rzero_de; ! spatial derivatives of forward propagated field

                if(.not.(present(snaps_in))) stop "fd2_mod: need snaps_in for gradient"
        endif

        if(present(snap_save_at)) allocate(saved_snap(mesh_nz, mesh_nx))

        ! source indices 
        allocate(issmulx1(src_nsmul),issmulx2(src_nsmul),issmulz1(src_nsmul),issmulz2(src_nsmul))
        issmulx1 = izero_de; issmulz1 = izero_de; issmulx2 = izero_de; issmulz2 = izero_de 
        allocate(denomsrcI(src_nsmul)); denomsrcI = rzero_de;
        allocate(src_spray_weights(src_nsmul,4)); src_spray_weights = rzero_de;
        do issmul = 1, src_nsmul
                if(index(source_flag,'[BILINEAR]') .ne. 0) then
                        maskX = .true.
                        issmulxtemp1 = &
                        minloc(array = abs(mesh_X - src_x(issmul + (isseq-1)*src_nsmul)), dim = 1, mask = maskX)
                        maskX(issmulxtemp1) = .false.
                        issmulxtemp2 = &
                        minloc(array = abs(mesh_X - src_x(issmul + (isseq-1)*src_nsmul)), dim = 1, mask = maskX)
                        issmulx1(issmul) = min(issmulxtemp1, issmulxtemp2); issmulx2(issmul) = max(issmulxtemp1, issmulxtemp2)

                        maskZ = .true.
                        issmulztemp1 = &
                        minloc(array = abs(mesh_Z - src_z(issmul + (isseq-1)*src_nsmul)), dim = 1, mask = maskZ)
                        maskZ(issmulztemp1) = .false.
                        issmulztemp2 = &
                        minloc(array = abs(mesh_Z - src_z(issmul + (isseq-1)*src_nsmul)), dim = 1, mask = maskZ)
                        issmulz1(issmul) = min(issmulztemp1, issmulztemp2); issmulz2(issmul) = max(issmulztemp1, issmulztemp2)

                        denomsrcI(issmul) = ((mesh_X(issmulx2(issmul)) - mesh_X(issmulx1(issmul)))*(mesh_Z(issmulz2(issmul))&
                                         - mesh_Z(issmulz1(issmul))))**(-1.e0)
                        ! for issmulz1, issmulx1
                        src_spray_weights(issmul,1) = (mesh_X(issmulx2(issmul))-src_x(issmul + (isseq-1)*src_nsmul))*&
                        (mesh_Z(issmulz2(issmul))-src_z(issmul + (isseq-1)*src_nsmul))*denomsrcI(issmul) 
                        ! for issmulz1, issmulx2
                        src_spray_weights(issmul,2) = (src_x(issmul + (isseq-1)*src_nsmul)-mesh_X(issmulx1(issmul)))*&
                        (mesh_Z(issmulz2(issmul))-src_z(issmul + (isseq-1)*src_nsmul))*denomsrcI(issmul)
                        ! for issmulz2, issmulx1
                        src_spray_weights(issmul,3) = (mesh_X(issmulx2(issmul))-src_x(issmul + (isseq-1)*src_nsmul))*&
                        (src_z(issmul + (isseq-1)*src_nsmul)-mesh_Z(issmulz1(issmul)))*denomsrcI(issmul)
                        ! for issmulz2, issmulx2 
                        src_spray_weights(issmul,4) = (src_x(issmul + (isseq-1)*src_nsmul)-mesh_X(issmulx1(issmul)))*&
                        (src_z(issmul + (isseq-1)*src_nsmul)-mesh_Z(issmulz1(issmul)))*denomsrcI(issmul)

                elseif(index(source_flag,'[NEAREST_NEIGH]') .ne. 0) then
                        issmulx1(issmul) = &
                        minloc(array = abs(mesh_X - src_x(issmul + (isseq-1)*src_nsmul)), dim = 1)

                        issmulz1(issmul) = &
                        minloc(array = abs(mesh_Z - src_z(issmul + (isseq-1)*src_nsmul)), dim = 1)
                        
                        ! dummy issmulz2 and issmulx2
                        issmulz2(issmul) = int(mesh_nz*0.5)
                        issmulx2(issmul) = int(mesh_nx*0.5)
                        
                        ! only add source at nearest grid point
                        src_spray_weights(issmul,1) = rone_de;
                endif
        enddo

        ! receiver indices
        if(present(recv_n)) then
                allocate(irecx1(recv_n),irecx2(recv_n),irecz1(recv_n),irecz2(recv_n))
                irecx1 = izero_de; irecz1 = izero_de; irecx2 = izero_de; irecz2 = izero_de 
                allocate(denomrecI(recv_n)); denomrecI = rzero_de;
                allocate(rec_interpolate_weights(recv_n,4)); rec_interpolate_weights = rzero_de;
                do irec = 1, recv_n
                        if(index(records_flag, '[BILINEAR]') .ne. 0) then
                                maskX = .true.
                                irecxtemp1 = &
                                minloc(array = &
                                abs(mesh_X - recv_x(irec + (isseq-1)*recv_n)), dim = 1, mask = maskX)
                                maskX(irecxtemp1) = .false.
                                irecxtemp2 = &
                                minloc(array = &
                                abs(mesh_X - recv_x(irec + (isseq-1)*recv_n)), dim = 1, mask = maskX)
                                irecx1(irec) = min(irecxtemp1, irecxtemp2); irecx2(irec) = max(irecxtemp1, irecxtemp2)

                                maskZ = .true.
                                irecztemp1 = &
                                minloc(array = &
                                abs(mesh_Z - recv_z(irec + (isseq-1)*recv_n)), dim = 1, mask = maskZ)
                                maskZ(irecztemp1) = .false.
                                irecztemp2 = &
                                minloc(array = &
                                abs(mesh_Z - recv_z(irec + (isseq-1)*recv_n)), dim = 1, mask = maskZ)
                                irecz1(irec) = min(irecztemp1, irecztemp2); irecz2(irec) = max(irecztemp1, irecztemp2)

                                denomrecI(irec) = ((mesh_X(irecx2(irec)) - mesh_X(irecx1(irec)))*(mesh_Z(irecz2(irec)) &
                                                - mesh_Z(irecz1(irec))))**(-1.e0)

                                ! irecz1, irecx1
                                rec_interpolate_weights(irec,1) = &
                                        ((mesh_X(irecx2(irec))-recv_x(irec + (isseq-1) * recv_n))*&
                                        (mesh_Z(irecz2(irec))-recv_z(irec + (isseq-1) * recv_n)))*denomrecI(irec)
                                ! irecz1, irecx2
                                rec_interpolate_weights(irec,2) = &
                                        ((recv_x(irec + (isseq-1) * recv_n)-mesh_X(irecx1(irec)))*&
                                        (mesh_Z(irecz2(irec))-recv_z(irec + (isseq-1) * recv_n)))*denomrecI(irec)
                                ! irecz2, irecx1
                                rec_interpolate_weights(irec,3) = &
                                        ((mesh_X(irecx2(irec))-recv_x(irec + (isseq-1) * recv_n))*&
                                        (recv_z(irec + (isseq-1) * recv_n)-mesh_Z(irecz1(irec))))*denomrecI(irec)
                                ! irecz2, irecx2
                                rec_interpolate_weights(irec,4) = &
                                        ((recv_x(irec + (isseq-1) * recv_n)-mesh_X(irecx1(irec)))*&
                                        (recv_z(irec + (isseq-1) * recv_n)-mesh_Z(irecz1(irec))))*denomrecI(irec)
                        elseif(index(records_flag,'[NEAREST_NEIGH]') .ne. 0) then
                                irecx1(irec) = &
                                minloc(array = &
                                abs(mesh_X - recv_x(irec + (isseq-1)*recv_n)), dim = 1)

                                irecz1(irec) = &
                                minloc(array = &
                                abs(mesh_Z - recv_z(irec + (isseq-1)*recv_n)), dim = 1)

                                ! dummy issmulz2 and issmulx2
                                irecz2(issmul) = int(mesh_nz*0.5)
                                irecx2(issmul) = int(mesh_nx*0.5)

                                ! record at nearest grid point
                                rec_interpolate_weights(irec,1) = rone_de;
                        endif
                enddo
        endif


        time_loop: do it = 1, ntmod
!                if((wtd .eq. "[MODELING]") .and. (mod(it, nint(ntmod/10.0)) .eq. 0)) then
!                        write(*,*) "MODELING: ",nint((it*100.0)/ntmod), &
!                        " % complete"
!                endif

                ! save p at [it-1] to pp
                pp = p

                !$OMKP PARALLEL NUM_THREADS(OMP_GET_MAX_THREADS()) DEFAULT(NONE) &
                !$OMKP PRIVATE(ix,iz) &
                !$OMKP SHARED(dvxdx, dvzdz, p,deltax24I,deltaz24I,a_x,a_z,memory_dvz_dz,lambda,&
                !$OMKP memory_dvx_dx,K_xI,K_zI,vx,vz,deltat,rhovxI,boundary_p,nx,nz,b_x,b_z)
                !$OMKP DO
                !forall(ix=2:nx,iz=2:nz)
                do ix=2,nx
                do iz=2,nz
                        dvxdx(iz,ix) = (27.e0*vx(iz,ix)-27.e0*vx(iz,ix-1)-vx(iz,ix+1)+vx(iz,ix-2)) * (deltax24I)
                        dvzdz(iz,ix) = (27.e0*vz(iz,ix)-27.e0*vz(iz-1,ix)-vz(iz+1,ix)+vz(iz-2,ix)) * (deltaz24I)
                        memory_dvx_dx(iz,ix) = b_x(ix) * memory_dvx_dx(iz,ix) + a_x(ix) * dvxdx(iz,ix) ! pml 
                        memory_dvz_dz(iz,ix) = b_z(iz) * memory_dvz_dz(iz,ix) + a_z(iz) * dvzdz(iz,ix) ! pml
                        dvxdx(iz,ix) = dvxdx(iz,ix) * k_xI(ix) + memory_dvx_dx(iz,ix) !pml
                        dvzdz(iz,ix) = dvzdz(iz,ix) * k_zI(iz) + memory_dvz_dz(iz,ix) !pml

                        ! compute pressure at [it] using p at [it-1] and dvxdx
                        ! and dvzdz at [it-1/2]
                        p(iz,ix) = p(iz,ix) + &
                                 (lambda(iz,ix) * (dvxdx(iz,ix) + dvzdz(iz,ix))) * deltat * boundary_p(iz,ix)
                enddo
                enddo
                !endforall
                !$OMKP ENDDO
                !$OMKP END PARALLEL

		! adding source to pressure at [it] 
		do issmul = 1, src_nsmul
                        ! compute source wavelet at [it-1/2] by simple
                        ! interpolation
                        ! division of source term with deltax and deltaz (see Jan's fdelmodc manual)
                        svalue  = 0.5e0 * (src_wav_mat(issmul, it, isseq) + src_wav_mat(issmul, it-1, isseq)) * deltat &
                                        * deltax24I * deltaz24I
                        p(na_pml+issmulz1(issmul), na_pml+issmulx1(issmul)) = p(na_pml+issmulz1(issmul), na_pml+issmulx1(issmul)) &
                                + svalue * src_spray_weights(issmul,1) * lambda(na_pml+issmulz1(issmul), na_pml+issmulx1(issmul))  
                        p(na_pml+issmulz1(issmul), na_pml+issmulx2(issmul)) = p(na_pml+issmulz1(issmul), na_pml+issmulx2(issmul)) &
                                + svalue * src_spray_weights(issmul,2) * lambda(na_pml+issmulz1(issmul), na_pml+issmulx2(issmul))
                        p(na_pml+issmulz2(issmul), na_pml+issmulx1(issmul)) = p(na_pml+issmulz2(issmul), na_pml+issmulx1(issmul)) &
                                + svalue * src_spray_weights(issmul,3) * lambda(na_pml+issmulz2(issmul), na_pml+issmulx1(issmul))
                        p(na_pml+issmulz2(issmul), na_pml+issmulx2(issmul)) = p(na_pml+issmulz2(issmul), na_pml+issmulx2(issmul)) &
                                + svalue * src_spray_weights(issmul,4) * lambda(na_pml+issmulz2(issmul), na_pml+issmulx2(issmul))
		enddo

                ! compute dpdx and dpdz at [it]
                !forall(ix=1:nx-1,iz=1:nz-1)
                do ix=1,nx-1
                do iz=1,nz-1
                        dpdx(iz,ix) = (27.e0*p(iz,ix+1)-27.e0*p(iz,ix)-p(iz,ix+2)+p(iz,ix-1)) * (deltax24I)
                        memory_dp_dx(iz,ix) = b_x_half(ix) * memory_dp_dx(iz,ix) + a_x_half(ix) * dpdx(iz,ix) ! pml
                        dpdx(iz,ix) = dpdx(iz,ix) * K_x_halfI(ix) + memory_dp_dx(iz,ix) ! pml
                enddo
                enddo
                !endforall

                !forall(ix=1:nx-1,iz=1:nz-1)
                do ix=1,nx-1
                do iz=1,nz-1
                        dpdz(iz,ix) = (27.e0*p(iz+1,ix)-27.e0*p(iz,ix)-p(iz+2,ix)+p(iz-1,ix)) * (deltaz24I)
                        memory_dp_dz(iz,ix) = b_z_half(iz) * memory_dp_dz(iz,ix) + a_z_half(iz) * dpdz(iz,ix) !pml
                        dpdz(iz,ix) = dpdz(iz,ix) * K_z_halfI(iz) + memory_dp_dz(iz,ix) !pml
                enddo
                enddo
                !endforall


                ! gradient calculation
                if(present(grad) .and. (it .ne. ntmod) .and. (it .ne. 1)) then

                        ! p at [it], pp at [it-1]
                        ! dpdx and dpdz at [it]
                        ! gradients w.r.t inverse of lambda, i.e., 1/rho/c^2 
                        forall(ix=1:mesh_nx,iz=1:mesh_nz)
                                grad_modTT(iz,ix) = grad_modTT(iz,ix) + &
                                        (snaps_in(iz,ix,ntmod-it,isseq)-snaps_in(iz,ix,ntmod-it+1,isseq)) * deltatI * &
                                        (p(iz+na_pml,ix+na_pml)-pp(iz+na_pml,ix+na_pml)) * deltatI
                        endforall
                        

                        if(allocated(grad(isseq)%rho)) then
                                forall(ix=2:mesh_nx-2,iz=2:mesh_nz-2)
                                        ! calculating dpFdx are on vx grid
                                        dpFdx(iz,ix) = &
                                                ((27.e0*snaps_in(iz,ix+1,ntmod-it,isseq)-27.e0*snaps_in(iz,ix,ntmod-it,isseq)-&
                                                snaps_in(iz,ix+2,ntmod-it,isseq)+snaps_in(iz,ix-1,ntmod-it,isseq)) &
                                                * (deltax24I))     
                                        ! calculating dpFdz are on vz grid
                                        dpFdz(iz,ix) = &
                                                ((27.e0*snaps_in(iz+1,ix,ntmod-it,isseq)-27.e0*snaps_in(iz,ix,ntmod-it,isseq)-&
                                                snaps_in(iz+2,ix,ntmod-it,isseq)+snaps_in(iz-1,ix,ntmod-it,isseq)) &
                                                * (deltaz24I))     
                                endforall
                                ! gradient w.r.t. rho_inv
                                forall(ix=2:mesh_nx,iz=2:mesh_nz)
!                                        grad_rho_inv(iz,ix) = grad_rho_inv(iz,ix) - &
!                                                        ((rhalf * (dpFdx(iz,ix) + dpFdx(iz,ix-1)) * &
!                                        rhalf * (dpdx(iz+na_pml,ix+na_pml) + dpdx(iz+na_pml,ix-1+na_pml))) + &
!                                                        (rhalf * (dpFdz(iz,ix) + dpFdz(iz-1,ix)) *&
!                                        rhalf * (dpdz(iz+na_pml,ix+na_pml) + dpdz(iz-1+na_pml,ix+na_pml))))
                                        grad_rhovxI(iz,ix) = grad_rhovxI(iz,ix) - dpdx(iz+na_pml,ix+na_pml)*dpFdx(iz,ix) 
                                        grad_rhovzI(iz,ix) = grad_rhovzI(iz,ix) - dpdz(iz+na_pml,ix+na_pml)*dpFdz(iz,ix) 
                                endforall
                        endif
                                
                endif

                ! saving a field snapshot (commented for speed)
                !if(allocated(saved_snap) ) then
                !        if(it .eq. int(ntmod * snap_save_at)) then
                !                saved_snap(:,:) = p(na_pml+1:mesh_nz+na_pml,1+na_pml:mesh_nx+na_pml) ;
                !        endif
                !endif

                ! store pressure seismograms to "records" 
                if(present(recv_n)) then
                        do irec = 1, recv_n
!                                        if(index(records_flag, '[P]') .ne. 0) then
                                rec_mat(irec,it) = &
                                        (&
                                        p(na_pml+irecz1(irec),na_pml+irecx1(irec))*rec_interpolate_weights(irec,1)+&
                                        p(na_pml+irecz1(irec),na_pml+irecx2(irec))*rec_interpolate_weights(irec,2)+&
                                        p(na_pml+irecz2(irec),na_pml+irecx1(irec))*rec_interpolate_weights(irec,3)+&
                                        p(na_pml+irecz2(irec),na_pml+irecx2(irec))*rec_interpolate_weights(irec,4)&
                                        )
!                                       endif
                        enddo
                endif

                if(present(snaps_out)) then
                        forall(ix=1:mesh_nx,iz=1:mesh_nz) 
                                snaps_out(iz,ix,it,isseq) = p(na_pml+iz,na_pml+ix)
                        endforall
                end if

                ! update velocity at [it+1/2] using 
                ! velcity at [it-1/2] and dpdx and dpdz at [it] 
                !forall(ix=1:nx-1,iz=1:nz-1)
                do ix=1,nx-1
                do iz=1,nz-1
                        vx(iz,ix) = vx(iz,ix) + (dpdx(iz,ix)) * deltat * rhovxI(iz,ix) * boundary_vx(iz,ix)
                enddo
                enddo
                !endforall

                !forall(ix=1:nx-1,iz=1:nz-1)
                do ix=1,nx-1
                do iz=1,nz-1
                        vz(iz,ix) = vz(iz,ix) + (dpdz(iz,ix)) * deltat * rhovzI(iz,ix) * boundary_vz(iz,ix)
                enddo
                enddo
                !endforall

        enddo time_loop

        ! output records
        if(present(records)) then
                records(1 + (isseq-1) * recv_n * ntmod : &
                        recv_n * ntmod + (isseq-1) * recv_n *ntmod) = pack(transpose(rec_mat),.true.)
        endif

        ! output gradient
        ! edges of the gradient zero
        if(present(grad)) then
                ! grad of vel
                grad(isseq)%vel = rzero_de;
                forall(ix=2:mesh_nx-2,iz=2:mesh_nz-2)
                        grad(isseq)%vel(iz,ix) = -2.0e0 * &
                                (var%vel(iz,ix))**(-3e0) * var%rho(iz,ix)**(-1e0) * grad_modTT(iz,ix)
                endforall
        

                grad_modRR = rzero_de;
                forall(ix=2:mesh_nx-2,iz=2:mesh_nz-2)
                        grad_modRR(iz,ix) = grad_modRR(iz,ix) + 0.5e0 * grad_rhovxI(iz,ix)
                        grad_modRR(iz,ix+1) = grad_modRR(iz,ix+1) + 0.5e0 * grad_rhovxI(iz,ix)
                        grad_modRR(iz,ix) = grad_modRR(iz,ix) + 0.5e0 * grad_rhovzI(iz,ix)
                        grad_modRR(iz+1,ix) = grad_modRR(iz+1,ix) + 0.5e0 * grad_rhovzI(iz,ix)
                endforall

                if(present(varpert)) then

                        allocate(rhoIpert(mesh_nz,mesh_nx), rhovIpert(mesh_nz,mesh_nx))
                        call get_rhoI(rho=varpert%rho, rhoI=rhoIpert)

                        call get_rhovxI(rhoI=rhoIpert, rhovxI=rhovIpert)
                        write(*,*) "adjoint validation: <deltarhovxI, grad_rhovxI>", &
                                dot_product(pack(rhovIpert-rhovxI,.true.), pack(grad_rhovxI,.true.))

                        call get_rhovzI(rhoI=rhoIpert, rhovzI=rhovIpert)
                        write(*,*) "adjoint validation: <deltarhovzI, grad_rhovzI>", &
                                dot_product(pack(rhovIpert-rhovzI,.true.), pack(grad_rhovzI,.true.))

                        write(*,*) "adjoint validation: <deltarhoI, grad_modRR>", &
                                dot_product(pack(rhoI-rhoIpert,.true.), pack(grad_modRR,.true.))

                        deallocate(rhoIpert, rhovIpert)


                        write(*,*) "adjoint validation: <deltalambdaI, grad_modTT>", &
                                dot_product(pack( &
                                (varpert%vel**2.0 * varpert%rho)**(-1.0) - (var%vel**2.0 * var%rho)**(-1.0), &
                                .true.), pack(grad_modTT,.true.))


                endif

 
                grad(isseq)%rho = rzero_de;
                forall(ix=2:mesh_nx-2,iz=2:mesh_nz-2)
                        grad(isseq)%rho(iz,ix) = -1.e0 * (var%rho(iz,ix))**(-2e0) * grad_modRR(iz,ix) &
                                - (var%vel(iz,ix) * var%rho(iz,ix))**(-2e0) * grad_modTT(iz,ix)
                endforall
        endif

        !if(allocated(saved_snap) .and. isseq .eq. 1) then
        !        if(present(jobname))    then
        !                junk = trim(out_dir)//"/snapshot"//trim(jobname)
        !        else
        !                junk = trim(out_dir)//"/snapshot"
        !        endif
        !        call makesu(dat=pack(saved_snap,.true.),fname=junk,&
        !                dx=deltax,dz=deltaz,nrec=mesh_nx,nt=mesh_nz, fileid=OMP_get_thread_num())
        !        deallocate(saved_snap)
        !endif


        deallocate                      (p, pp,&
                                         vx,vz)
        deallocate			(dpdx, dpdz)
        deallocate			(dvxdx, dvzdz)

        deallocate                      (memory_dp_dx,&
                                         memory_dp_dz,&
                                         memory_dvx_dx,&
                                         memory_dvx_dz,&
                                         memory_dvz_dx,&
                                         memory_dvz_dz)
        deallocate(issmulx1, issmulx2, issmulz1, issmulz2)
        deallocate(denomsrcI, src_spray_weights)
        if(present(recv_n)) then
                deallocate(irecx1, irecx2, irecz1, irecz2)
                deallocate(denomrecI, rec_interpolate_weights)
        endif

        if(allocated(rec_mat)) deallocate  (rec_mat)
        if(allocated(dpFdz)) deallocate(dpFdz)
        if(allocated(dpFdx)) deallocate(dpFdx)
        if(allocated(saved_snap)) deallocate(saved_snap)
        if(allocated(grad_modRR)) deallocate(grad_modRR)
        if(allocated(grad_rhovxI)) deallocate(grad_rhovxI)
        if(allocated(grad_rhovzI)) deallocate(grad_rhovzI)
        if(allocated(grad_modTT)) deallocate(grad_modTT)


enddo src_par_loop
!$OMP ENDDO 
!$OMP END PARALLEL

deallocate                      (d_x,K_x,K_xI,alpha_x,&
                                 a_x,b_x,d_x_half,&
                                 K_x_half,K_x_halfI,alpha_x_half,a_x_half,b_x_half)

deallocate                      (d_z,K_z,K_zI,alpha_z,&
                                 a_z,b_z,d_z_half,&
                                 K_z_half,K_z_halfI,alpha_z_half,a_z_half,b_z_half)
deallocate                      (lambda, rho, rhoI, rhovxI, rhovzI)

deallocate(maskX); deallocate(maskZ)
if(allocated(src_wav_mat)) deallocate(src_wav_mat)

! testing the speed -- uncomment if u want to print out cpu time
!cpu_time2  = omp_get_wtime(); write(*,*) "fd2_mod: ",(cpu_time2-cpu_time1),"seconds"


end subroutine fdtd_mod

end module fdtd
!\end{comment}
