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
use string_routines! USEINDEPFILE
use omp_lib

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
! >>>>>>> general inputs
        jobname, & ! name of the subroutine call
        npropwav, & ! number of wavefields propagating independently in the same medium
                    ! inside the time loop
! >>>>>>> model parameters
        modtt, & ! model parameter matrix before second order time derivative [mesh_nz][mesh_nx]
        modrr, & ! model parameter vector before  [mesh_nz][mesh_nx]
! >>>>>>> background model parameters
        modtt0, & ! model parameter matrix before second order time derivative [mesh_nz][mesh_nx]
        modrr0, & ! model parameter vector before  [mesh_nz][mesh_nx]
! >>>>>>> modeling mesh parameters
        mesh_nx, & !
        mesh_nz, & !
        mesh_dx, & ! 
        mesh_dz, & ! 
        mesh_x, & ! modeling mesh X vector
        mesh_z, & ! modeling mesh Z vector
        mesh_na_pml, & ! number of PML layers to be added
        mesh_abs_trbl, & ! absorbing boundary conditions Top, Right, Bottom, Left
! >>>>>>> time parameters
        tim_nt,& ! total time samples
        tim_del, & ! time sampling
! >>>>>>> source parameters
        src_nt, & ! number time samples in source wavelet (src_nt <= tim_nt)
        src_wavelets, & ! source wavelets [src_nt] * [src_nsmul] * [src_nseq] * [src_nfield] * [npropwav]
        src_nseq, & ! number of sequential sources; OMP loop
        src_nsmul, & ! number of simultaneous sources
        src_nfield, & ! number of source fields
        src_x, & ! x coordinates of all the sources [src_nsmul] * [src_nseq] * [npropwav]
        src_z, & ! z coordinated of all the sources [src_nsmul] * [src_nseq] * [npropwav]
        src_flags, & ! parameters determining source 
                        ! contains [IRATE], then injection rate is ON
                        ! contains [TIME_REVERSE], then inject time reverse of the src_wavelets
                        ! contains [BILINEAR], then bilinear interpolation
                        ! contains [NEAREST_NEIGH], then nearest neighbour interpolation
! >>>>>>> receiver parameters
        recv_n, & ! number of receivers for each source 
        recv_nfield, & ! number of field recorded        
        recv_x, & ! x coord for receivers [recv_n] * [src_nseq] * [npropwav]
        recv_z, & ! z  "                  [recv_n] * [src_nseq] * [npropwav]
        recv_out, & ! recorded pressure at receivers [tim_nt] * [recv_n] * [src_nseq] * [recv_nfield] * [npropwav]
        recv_flags, & ! 0 means no records; 1 means nearest neighbour; 2 means bilinear [npropwav]
! >>>>>>> flags between propagating wavefields
        prop_flags, &
! >>>>>>> indices of borders
        border_n, & ! total number of grid points on the boundary
        border_x, & ! x indices of boundary 
        border_z, & ! z indices of boundary
! >>>>>>> boundary storage -- input; only for the first propagating field
        border_in_flag, & ! border_in used when true
        border_in, & !  
! >>>>>>> boundary storage -- output; only for the first propagating field
        border_out_flag, & ! border_out used when true
        border_out, & ! 
! >>>>>>> initial values; only for the first propagating field
        p_in, &
        p_out, &
! >>>>>>> gradients; assumes that the first propagating wavefield is 
                        !forward state and the second is adjoint state
        grad_out_flag, & 
        grad_modtt,  &  ! output gradient models w.r.t velocity and density  [src_nseq]      
        grad_modrr, & ! output gradient 
! >>>>>>>
        born_flag &
        ) bind(c, name="fdtd_mod")
! Author: Pawan Bharadwaj
!         p.b.pisupati@tudelft.nl
!         modified: 11 Sept 2013
!         updated on 25 July 2014
                ! parallelized w.r.t src_nseq
                ! added src_nt so that source wavelets can have different
                ! dimensions 

        implicit none
character(len=1, kind=C_char), intent(in)               :: jobname(*), mesh_abs_trbl(*)
real, intent(in)                                        :: recv_flags(npropwav), src_flags(npropwav)
character(len=1, kind=C_char), intent(in)               :: prop_flags(*)
integer, intent(in)                                     :: npropwav
real, dimension(mesh_nz,mesh_nx), intent(in)            :: modtt, modrr, modtt0, modrr0


integer, intent(in)                                     :: mesh_nx, mesh_nz
real, intent(in)                                        :: mesh_x(mesh_nx), mesh_z(mesh_nz), mesh_dx, mesh_dz
integer, intent(in)                                     :: mesh_na_pml  !number of PML layers


real, intent(in)                                        :: tim_del
integer, intent(in)                                     :: tim_nt
integer, intent(in)                                     :: src_nfield
integer, intent(in)                                     :: recv_nfield


integer, intent(in)                                     :: src_nt, src_nsmul, src_nseq
real, dimension(src_nt*src_nsmul*src_nseq*src_nfield*npropwav), intent(in) &
                                                        :: src_wavelets
real, dimension(src_nsmul*src_nseq*npropwav), intent(in):: src_x, src_z

integer, intent(in)                                     :: recv_n
real, dimension(tim_nt*recv_n*src_nseq*recv_nfield*npropwav), intent(out) &
                                                        :: recv_out
real, dimension(recv_n*src_nseq*npropwav), intent(in) &
                                                        :: recv_x, recv_z


integer, intent(in)                                     :: border_n
real, dimension(border_n,tim_nt,src_nseq),   intent(in) :: border_in
real, dimension(border_n,tim_nt,src_nseq),   intent(out):: border_out
logical, intent(in)                                     :: border_in_flag, border_out_flag

real, dimension(border_n), intent(in) &
                                                        :: border_x, border_z


real, intent(in)                                        :: p_in(0:mesh_nz+1,0:mesh_nx+1,3,src_nseq,2)
real, intent(out)                                       :: p_out(0:mesh_nz+1,0:mesh_nx+1,3,src_nseq,2)                                                


real, dimension(1:mesh_nz,1:mesh_nx,src_nseq), intent(inout) &
                                                        :: grad_modtt, grad_modrr
logical, intent(in)                                     :: born_flag, grad_out_flag 

! other variables for this only this subroutine
! total number of grid points in each direction of the grid
integer                                                 :: nx
integer                                                 :: nz

! size of a grid cell
real                                                    :: deltax, deltax24I 
real                                                    :: deltaz, deltaz24I
real                                                    :: tim_delI 

integer                                                 :: ix, iz, it, ixx, izz, irec, issmul, isseq, ipropwav
integer                                                 :: ifield, iborder
real, dimension(:,:,:,:,:), allocatable                 :: src_wav_mat
real                                                    :: source_term

! source position variables                            
integer, allocatable, dimension(:,:)                    :: issmulx1, issmulx2, issmulz1, issmulz2
integer                                                 :: issmulxtemp1, issmulxtemp2, issmulztemp1, issmulztemp2, &
                                                                issmultemp
real, allocatable, dimension(:,:)                       :: denomsrcI
real, allocatable, dimension(:,:,:)                     :: src_spray_weights
real                                                    :: svalue, src_xtemp, src_ztemp
! receiver position variables           
integer, allocatable, dimension(:,:)                    :: irecx1, irecx2, irecz1, irecz2
integer                                                 :: irecxtemp1, irecxtemp2, irecztemp1, irecztemp2, irectemp
real, dimension(:,:), allocatable                       :: denomrecI
real, dimension(:,:,:), allocatable                     :: rec_interpolate_weights
real                                                    :: recv_xtemp, recv_ztemp


integer, allocatable, dimension(:)                      :: iborderx, iborderz
! recording data
real, dimension(:,:,:,:), allocatable                   :: rec_mat

! to find indices
logical, allocatable                                    :: maskX(:), maskZ(:)

! extended models 
real, dimension(:,:), allocatable                       :: modttI, modrrvx, modrrvz
real, dimension(:,:), allocatable                       :: deltamodtt, deltamodrrvx, deltamodrrvz

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
real, dimension(:,:,:,:), allocatable                   :: p
real, dimension(:,:,:,:), allocatable                   :: pp, ppp ! for time differential of pressure
real, dimension(:,:,:), allocatable                     :: dpdx,dpdz ! spatial derivatives of pressure in vx and vz grids
real, dimension(:,:,:), allocatable                     :: dvxdx,dvzdz

real, dimension(:,:), allocatable                       :: gradis_modrrvx, gradis_modrrvz, gradis_modtt

real, allocatable                                       :: saved_snap(:,:)

real(r_dp)                                              :: cpu_time1, cpu_time2

! arrays for the memory variables
! could declare these arrays in PML only to save a lot of memory, but proof of concept only here
real, dimension(:,:,:), allocatable                     :: memory_dvx_dx, &
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

! output to zeros
recv_out = rzero_de
border_out = rzero_de
p_out = rzero_de
grad_modtt =  rzero_de;
grad_modrr = rzero_de;

cpu_time1  = omp_get_wtime();

! convert src_wavelets to src_wav_mat
allocate(src_wav_mat(src_nsmul, src_nfield, npropwav, 0:tim_nt+1, src_nseq))
src_wav_mat = rzero_de;
do ipropwav = 1, npropwav
        do ifield = 1, src_nfield
        do isseq = 1, src_nseq
                do issmul = 1, src_nsmul
                        source_term = rzero_de
                        do it = 1, src_nt-1
                              !  if(index(ctofstr(src_flags),'[TIME_REVERSE]') .ne. 0) then
                              !          source_term = &
                              !          src_wavelets(src_nt - it  &
                              !                  + (issmul-1)*src_nt &
                              !                  + (isseq-1)*src_nsmul*src_nt &
                              !                  + (ifield-1)*src_nseq*src_nsmul*src_nt &
                              !                  + (ipropwav-1)*src_nfield*src_nseq*src_nsmul*src_nt &
                              !                          ) 
                              !  else
                                        source_term = &
                                        src_wavelets(it &
                                                + (issmul-1)*src_nt &
                                                + (isseq-1)*src_nsmul*src_nt &
                                                + (ifield-1)*src_nseq*src_nsmul*src_nt &
                                                + (ipropwav-1)*src_nseq*src_nfield*src_nsmul*src_nt &
                                                )
                              !  endif
                                if(src_flags(ipropwav) .lt. 0) then
                                        src_wav_mat(issmul, ifield, ipropwav, it + 1, isseq) = 2.0 * source_term + &
                                                src_wav_mat(issmul, ifield, ipropwav, it - 1, isseq)
                                else
                                        src_wav_mat(issmul, ifield, ipropwav, it, isseq) = source_term
                                endif
                        enddo
                enddo
        enddo
        enddo
enddo
if(maxval(src_wav_mat) .ge. huge(1.0e0)) call abort_msg("fd2_mod: huge src_wav_mat?")

if(born_flag .or. grad_out_flag) then
        if(npropwav .lt. 2) stop "fdtd_mod: need atleast 2 propagating wavefields"
endif

! allocate masks
allocate(maskX(mesh_nx)); maskX = .true.
allocate(maskZ(mesh_nz)); maskZ = .true.

if(mesh_na_pml .le. 0) call abort_msg("fd2_mod: invalid mesh_na_pml")

! extended dimensions of the modeling grid
nx    = mesh_nx
nz    = mesh_nz

! assign variables
deltax=mesh_dx
deltax24I = (deltax * 24.e0)**(-1.e0)
deltaz=mesh_dz
deltaz24I = (deltaz * 24.e0)**(-1.e0)
tim_delI= (tim_del)**(-1.e0)

! minimum and maximum velocities
velmin = minval(sqrt(modrr*modtt**(-rone_de))); 
velmax = maxval(sqrt(modrr*modtt**(-rone_de))); 
write(*,*) "minimum and maximum velocities", velmin, velmax

! minimum and maximum frequencies
freqmin = huge(1.e0); freqmax = 0.0;
do ipropwav = 1, npropwav
        do ifield=1, src_nfield
        do isseq = 1, src_nseq
                do issmul = 1, src_nsmul
                if(maxval(abs(src_wav_mat(issmul,ifield,ipropwav,:,isseq))) .gt. tiny(rzero_de)) then
                                freqmin = &
                                        min(freqmin,findfreq(x=pack(src_wav_mat(issmul,ifield,ipropwav,:,isseq),.true.),&
                                dt=tim_del,str="MIN",threshold=2e-2))
                                freqmax = &
                                        max(freqmax,findfreq(x=pack(src_wav_mat(issmul,ifield,ipropwav,:,isseq),.true.),&
                                dt=tim_del, str="MAX",threshold=2e-2))
                        endif
                enddo
        enddo
        enddo
enddo
ALPHA_MAX_PML = 2.e0*PI*((freqmin + freqmax))/4.e0 ! from Festa and Vilotte

! model parameters - modttI
allocate(modttI(0:nz+1,0:nx+1), deltamodtt(0:nz+1,0:nx+1)); modttI=rzero_de; deltamodtt = rzero_de;
modttI = modtt**(-rone_de) 
deltamodtt = modtt-modtt0

! model parameters - modrrvx and modrrvz
allocate(modrrvx(0:nz+1,0:nx+1),modrrvz(0:nz+1,0:nx+1)); 
allocate(deltamodrrvx(0:nz+1,0:nx+1),deltamodrrvz(0:nz+1,0:nx+1));
modrrvx(:,:) = rzero_de; modrrvz(:,:) = rzero_de
deltamodrrvx = rzero_de; deltamodrrvz = rzero_de;

modrrvx = get_rhovxI(modrr)
modrrvz = get_rhovzI(modrr)

deltamodrrvx = modrrvx - get_rhovxI(modrr0)
deltamodrrvz = modrrvz - get_rhovzI(modrr0)

call is_nan("fd2_mod: modttI", modttI)
call is_nan("fd2_mod: modrrvz", modrrvz)
call is_nan("fd2_mod: modrrvx", modrrvx)

! boundary conditions on p
allocate(boundary_p(0:nz+1,0:nx+1)); allocate(boundary_vx(0:nz+1,0:nx+1));
allocate(boundary_vz(0:nz+1,0:nx+1))
! default
boundary_p = rone_de; boundary_vx = rone_de; boundary_vz = rone_de;
USE_PML_ZMIN = .true.; USE_PML_XMAX = .true.; USE_PML_ZMAX = .true.; USE_PML_XMIN = .true.

if(index(ctofstr(mesh_abs_trbl),"[T]") .eq. 0) USE_PML_ZMIN = .false.
if(index(ctofstr(mesh_abs_trbl),"[R]") .eq. 0) USE_PML_XMAX = .false.
if(index(ctofstr(mesh_abs_trbl),"[L]") .eq. 0) USE_PML_XMIN = .false.
if(index(ctofstr(mesh_abs_trbl),"[B]") .eq. 0) USE_PML_ZMAX = .false.


! free surface boundary conditions on sides of the model.
if(ctofstr(mesh_abs_trbl) .eq. "0") then
        boundary_p(1:mesh_na_pml+1,:) = rzero_de
        !boundary_vz(1:mesh_na_pml,:) = rzero_de
!        boundary_vx(1:mesh_na_pml,:) = rzero_de
!        USE_PML_ZMIN = .false.
endif
if(ctofstr(mesh_abs_trbl) .eq. "0") then
        boundary_p(1:nz,mesh_na_pml+mesh_nx:nx) = rzero_de
        USE_PML_XMAX = .false.
endif
if(ctofstr(mesh_abs_trbl(3:3)) .eq. "0") then 
        boundary_p(mesh_na_pml+mesh_nz,1:nx) = rzero_de
        USE_PML_ZMAX = .false.
endif
if(ctofstr(mesh_abs_trbl(4:4)) .eq. "0") then
        boundary_p(1:nz,mesh_na_pml+1) = rzero_de
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
thickness_PML_x = mesh_na_pml * deltax
thickness_PML_z = mesh_na_pml * deltaz

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

        b_x(ix) = exp(- (d_x(ix) / K_x(ix) + alpha_x(ix)) * tim_del)
        b_x_half(ix) = exp(- (d_x_half(ix) / K_x_half(ix) + alpha_x_half(ix)) * tim_del)

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

        b_z(iz) = exp(- (d_z(iz) / K_z(iz) + alpha_z(iz)) * tim_del)
        b_z_half(iz) = exp(- (d_z_half(iz) / K_z_half(iz) + alpha_z_half(iz)) * tim_del)

        ! this to avoid division by zero outside the PML
        if(abs(d_z(iz)) > 1.e-6) a_z(iz) = d_z(iz) * (b_z(iz) - 1.e0) / (K_z(iz) * (d_z(iz) + K_z(iz) * alpha_z(iz)))
        if(abs(d_z_half(iz)) > 1.e-6) a_z_half(iz) = d_z_half(iz) * &
        (b_z_half(iz) - 1.e0) / (K_z_half(iz) * (d_z_half(iz) + K_z_half(iz) * alpha_z_half(iz)))

enddo
K_zI = K_z**(-1.e0)
K_z_halfI = K_z_halfI**(-1.e0)

! check the Courant stability condition for the explicit time scheme
! R. Courant et K. O. Friedrichs et H. Lewy (1928)
courant_number = velmin * tim_del * sqrt(1.e0/deltax**2 + 1.e0/deltaz**2)
if(courant_number > 1.e0) call abort_msg("fd2_mod: time step is too large, simulation will be unstable")


allocate(iborderx(border_n), iborderz(border_n))
! border indices
do iborder = 1, border_n
        iborderx(iborder) = &
                minloc(array = &
                abs(mesh_x - border_x(iborder)), dim = 1)

        iborderz(iborder) = &
                minloc(array = &
                abs(mesh_z - border_z(iborder)), dim = 1)
enddo


!$OMP PARALLEL NUM_THREADS(OMP_GET_MAX_THREADS()) DEFAULT(NONE) &
!$OMP PRIVATE(ipropwav, isseq, it, issmul, ixx, izz, maskX, maskZ, denomsrcI, denomrecI, svalue, &
!$OMP issmulxtemp1,issmulxtemp2,issmulztemp1,issmulztemp2,issmulx1,issmulz1,issmulx2,issmulz2,src_spray_weights, &
!$OMP issmultemp, src_ztemp, src_xtemp, &
!$OMP irecxtemp1,irecxtemp2,irecztemp1,irecztemp2,irecx1,irecz1,irecx2,irecz2,rec_interpolate_weights, &
!$OMP irectemp, recv_xtemp, recv_ztemp, iborder, &
!$OMP pp,ppp,dpdx,dpdz,dvxdx,dvzdz,gradis_modrrvx,gradis_modrrvz, gradis_modtt,&
!$OMP p, &
!$OMP memory_dvx_dx, &
!$OMP memory_dvx_dz, &
!$OMP memory_dvz_dx, &
!$OMP memory_dvz_dz, &
!$OMP memory_dp_dx, &
!$OMP memory_dp_dz, &
!$OMP saved_snap, rec_mat) &
!$OMP SHARED(npropwav,src_nseq,src_nfield,src_wav_mat,src_flags, mesh_nz, mesh_nx, nz, nx, border_in, border_in_flag,&
!$OMP border_out_flag, mesh_x, mesh_z, p_in, p_out, grad_modtt, grad_modrr, &
!$OMP border_out, grad_out_flag, tim_nt, src_nt, src_nsmul, deltax, deltaz, deltax24I, deltaz24I, tim_del, tim_delI, & 
!$OMP src_x, src_z, recv_n, recv_z, recv_x, recv_out, recv_flags, recv_nfield, &
!$OMP iborderx, iborderz, border_n, born_flag,&
!$OMP jobname, modttI, modrrvx,modrrvz, &
!$OMP deltamodtt, deltamodrrvx, deltamodrrvz, &
!$OMP b_x_half, a_x_half, b_z_half, a_z_half, b_x, b_z, a_x, a_z, &
!$OMP K_x, K_z, K_xI, K_zI,  &
!$OMP K_x_half, K_z_half, K_x_halfI, K_z_halfI, &
!$OMP boundary_p, boundary_vx, boundary_vz) 

!$OMP DO 
src_par_loop: do isseq = 1, src_nseq

        allocate (p(0:nz+1,0:nx+1,3,npropwav), pp(0:nz+1,0:nx+1,3,npropwav), &
                ppp(0:nz+1,0:nx+1,3,npropwav))
        allocate (dpdx(0:nz+1,0:nx+1,npropwav), dpdz(0:nz+1,0:nx+1,npropwav))
        allocate (dvxdx(0:nz+1,0:nx+1,npropwav), dvzdz(0:nz+1,0:nx+1,npropwav))

        allocate (memory_dvx_dx(0:nz+1,0:nx+1,npropwav), &
                memory_dvx_dz(0:nz+1,0:nx+1,npropwav), &
                memory_dvz_dx(0:nz+1,0:nx+1,npropwav), &
                memory_dvz_dz(0:nz+1,0:nx+1,npropwav), &
                memory_dp_dx(0:nz+1,0:nx+1,npropwav), &
                memory_dp_dz(0:nz+1,0:nx+1,npropwav))

        ! initialize main arrays
        p = rzero_de; pp = rzero_de; ppp =rzero_de
        ! initial conditions for propagating field 1
        p(:,:,:,1) = p_in(:,:,:,isseq,1); pp(:,:,:,1) = p_in(:,:,:,isseq,2); ppp=rzero_de;
        dpdx = rzero_de; dpdz = rzero_de;
        dvxdx = rzero_de; dvzdz = rzero_de;

        ! temp array per thread to store recv_out
        allocate   (rec_mat(recv_n, recv_nfield, npropwav, tim_nt)); rec_mat = rzero_de

        ! initialize PML arrays
        memory_dp_dx  = rzero_de
        memory_dp_dz  = rzero_de
        memory_dvx_dx = rzero_de
        memory_dvx_dz = rzero_de
        memory_dvz_dx = rzero_de
        memory_dvz_dz = rzero_de

        if(grad_out_flag) then
                allocate(gradis_modrrvx(0:nz+1, 0:nx+1), gradis_modrrvz(0:nz+1, 0:nx+1))
                allocate(gradis_modtt(0:nz+1, 0:nx+1))
                gradis_modrrvx =  rzero_de; gradis_modrrvz = rzero_de;
                gradis_modtt = rzero_de;
        endif


        ! source indices 
        allocate(issmulx1(src_nsmul,npropwav),issmulx2(src_nsmul,npropwav), &
                issmulz1(src_nsmul,npropwav),issmulz2(src_nsmul,npropwav))
        issmulx1 = izero_de; issmulz1 = izero_de; issmulx2 = izero_de; issmulz2 = izero_de 
        allocate(denomsrcI(src_nsmul,npropwav)); denomsrcI = rzero_de;
        allocate(src_spray_weights(4,src_nsmul,npropwav)); src_spray_weights = rzero_de;
        do ipropwav = 1, npropwav
                do issmul = 1, src_nsmul
                        issmultemp = issmul+(isseq-1)*src_nsmul+(ipropwav-1)*src_nsmul*src_nseq;
                        src_xtemp = src_x(issmultemp);
                        src_ztemp = src_z(issmultemp);
                        if(abs(src_flags(ipropwav)) .eq. 2) then ! [BILINEAR]
                                maskX = .true.
                                issmulxtemp1 = &
                                minloc(array = abs(mesh_x - &
                                src_xtemp), &
                                dim = 1, mask = maskX)
                                maskX(issmulxtemp1) = .false.
                                issmulxtemp2 = &
                                minloc(array = abs(mesh_x - &
                                src_xtemp), &
                                dim = 1, mask = maskX)
                                issmulx1(issmul,ipropwav) = min(issmulxtemp1, issmulxtemp2); 
                                issmulx2(issmul,ipropwav) = max(issmulxtemp1, issmulxtemp2)

                                maskZ = .true.
                                issmulztemp1 = &
                                minloc(array = abs(mesh_z - &
                                src_ztemp), &
                                dim = 1, mask = maskZ)
                                maskZ(issmulztemp1) = .false.
                                issmulztemp2 = &
                                minloc(array = abs(mesh_z - &
                                src_ztemp), &
                                dim = 1, mask = maskZ)
                                issmulz1(issmul,ipropwav) = min(issmulztemp1, issmulztemp2); 
                                issmulz2(issmul,ipropwav) = max(issmulztemp1, issmulztemp2)

                                denomsrcI(issmul,ipropwav) = ((mesh_x(issmulx2(issmul,ipropwav)) - &
                                                    mesh_x(issmulx1(issmul,ipropwav)))*&
                                                    (mesh_z(issmulz2(issmul,ipropwav))&
                                                    - mesh_z(issmulz1(issmul,ipropwav))))**(-1.e0)
                                ! for issmulz1, issmulx1
                                src_spray_weights(1,issmul,ipropwav) = &
                                (mesh_x(issmulx2(issmul,ipropwav))-&
                                src_xtemp)*&
                                (mesh_z(issmulz2(issmul,ipropwav))-&
                                src_ztemp)*&
                                denomsrcI(issmul,ipropwav) 
                                ! for issmulz1, issmulx2
                                src_spray_weights(2,issmul,ipropwav) = &
                                (src_xtemp-&
                                mesh_x(issmulx1(issmul,ipropwav)))*&
                                (mesh_z(issmulz2(issmul,ipropwav))-&
                                src_ztemp)*&
                                denomsrcI(issmul,ipropwav)
                                ! for issmulz2, issmulx1
                                src_spray_weights(3,issmul,ipropwav) = &
                                (mesh_x(issmulx2(issmul,ipropwav))-&
                                src_xtemp)*&
                                (src_ztemp-&
                                mesh_z(issmulz1(issmul,ipropwav)))*&
                                denomsrcI(issmul,ipropwav)
                                ! for issmulz2, issmulx2 
                                src_spray_weights(4,issmul,ipropwav) = &
                                (src_xtemp-&
                                mesh_x(issmulx1(issmul,ipropwav)))*&
                                (src_ztemp-&
                                mesh_z(issmulz1(issmul,ipropwav)))*&
                                denomsrcI(issmul,ipropwav)

                        elseif(abs(src_flags(ipropwav)) .eq. 1) then ! [NEAREST_NEIGH]
                                issmulx1(issmul,ipropwav) = &
                                minloc(array = abs(mesh_x - &
                                src_xtemp), dim = 1)

                                issmulz1(issmul,ipropwav) = &
                                minloc(array = abs(mesh_z - &
                                src_ztemp), dim = 1)
                                
                                ! dummy issmulz2 and issmulx2
                                issmulz2(issmul,ipropwav) = int(mesh_nz*0.5)
                                issmulx2(issmul,ipropwav) = int(mesh_nx*0.5)
                                
                                ! only add source at nearest grid point
                                src_spray_weights(1,issmul,ipropwav) = rone_de;
                        elseif(abs(src_flags(ipropwav)) .eq. 0) then ! [OFF]
                                ! dummy 
                                issmulz1(issmul,ipropwav) = int(mesh_nz*0.5)
                                issmulx1(issmul,ipropwav) = int(mesh_nx*0.5)
                                issmulz2(issmul,ipropwav) = int(mesh_nz*0.5)
                                issmulx2(issmul,ipropwav) = int(mesh_nx*0.5)
                                src_spray_weights(1,issmul,ipropwav) = rzero_de;

                        endif
                enddo
        enddo

        ! receiver indices
        allocate(irecx1(recv_n,npropwav),irecx2(recv_n,npropwav),&
                irecz1(recv_n,npropwav),irecz2(recv_n,npropwav))
        irecx1 = izero_de; irecz1 = izero_de; irecx2 = izero_de; irecz2 = izero_de 
        allocate(denomrecI(recv_n,npropwav)); denomrecI = rzero_de;
        allocate(rec_interpolate_weights(4,recv_n,npropwav)); rec_interpolate_weights = rzero_de;
        do ipropwav = 1, npropwav
                do irec = 1, recv_n
                        irectemp = irec+(isseq-1)*recv_n+(ipropwav-1)*recv_n*src_nseq
                        recv_xtemp = recv_x(irectemp);
                        recv_ztemp = recv_z(irectemp);
                        if(recv_flags(ipropwav) .eq. 2) then
                                maskX = .true.
                                irecxtemp1 = &
                                minloc(array = &
                                abs(mesh_x - recv_xtemp), dim = 1, mask = maskX)
                                maskX(irecxtemp1) = .false.
                                irecxtemp2 = &
                                minloc(array = &
                                abs(mesh_x - recv_xtemp), dim = 1, mask = maskX)
                                irecx1(irec,ipropwav) = min(irecxtemp1, irecxtemp2); 
                                irecx2(irec,ipropwav) = max(irecxtemp1, irecxtemp2)

                                maskZ = .true.
                                irecztemp1 = &
                                minloc(array = &
                                abs(mesh_z - recv_ztemp), &
                                dim = 1, mask = maskZ)
                                maskZ(irecztemp1) = .false.
                                irecztemp2 = &
                                minloc(array = &
                                abs(mesh_z - recv_ztemp), &
                                dim = 1, mask = maskZ)
                                irecz1(irec,ipropwav) = min(irecztemp1, irecztemp2); 
                                irecz2(irec,ipropwav) = max(irecztemp1, irecztemp2)

                                denomrecI(irec,ipropwav) = ((mesh_x(irecx2(irec,ipropwav)) - &
                                        mesh_x(irecx1(irec,ipropwav)))*(mesh_z(irecz2(irec,ipropwav)) &
                                                - mesh_z(irecz1(irec,ipropwav))))**(-1.e0)

                                ! irecz1, irecx1
                                rec_interpolate_weights(1,irec,ipropwav) = &
                                        ((mesh_x(irecx2(irec,ipropwav))-&
                                        recv_xtemp)*&
                                        (mesh_z(irecz2(irec,ipropwav))-&
                                        recv_ztemp))*&
                                        denomrecI(irec,ipropwav)
                                ! irecz1, irecx2
                                rec_interpolate_weights(2,irec,ipropwav) = &
                                        ((recv_xtemp-&
                                        mesh_x(irecx1(irec,ipropwav)))*&
                                        (mesh_z(irecz2(irec,ipropwav))-&
                                        recv_ztemp))*&
                                        denomrecI(irec,ipropwav)
                                ! irecz2, irecx1
                                rec_interpolate_weights(3,irec,ipropwav) = &
                                        ((mesh_x(irecx2(irec,ipropwav))-&
                                        recv_xtemp)*&
                                        (recv_ztemp-&
                                        mesh_z(irecz1(irec,ipropwav))))*&
                                        denomrecI(irec,ipropwav)
                                ! irecz2, irecx2
                                rec_interpolate_weights(4,irec,ipropwav) = &
                                        ((recv_xtemp-&
                                        mesh_x(irecx1(irec,ipropwav)))*&
                                        (recv_ztemp-&
                                        mesh_z(irecz1(irec,ipropwav))))*&
                                        denomrecI(irec,ipropwav)
                        elseif(recv_flags(ipropwav) .eq. 1) then
                                irecx1(irec,ipropwav) = &
                                minloc(array = &
                                abs(mesh_x - recv_xtemp), dim = 1)

                                irecz1(irec,ipropwav) = &
                                minloc(array = &
                                abs(mesh_z - recv_ztemp), dim = 1)

                                ! dummies
                                irecz2(irec,ipropwav) = int(mesh_nz*0.5)
                                irecx2(irec,ipropwav) = int(mesh_nx*0.5)

                                ! record at nearest grid point
                                rec_interpolate_weights(1,irec,ipropwav) = rone_de;
                        elseif(recv_flags(ipropwav) .eq. 0) then
                                ! dummies
                                irecz1(irec,ipropwav) = int(mesh_nz*0.5)
                                irecx1(irec,ipropwav) = int(mesh_nx*0.5)
                                irecz2(irec,ipropwav) = int(mesh_nz*0.5)
                                irecx2(irec,ipropwav) = int(mesh_nx*0.5)
                                ! record at nearest grid point
                                rec_interpolate_weights = rzero_de;
                        endif
                enddo
        enddo


        time_loop: do it = 1, tim_nt
!                if((wtd .eq. "[MODELING]") .and. (mod(it, nint(tim_nt/10.0)) .eq. 0)) then
!                        write(*,*) "MODELING: ",nint((it*100.0)/tim_nt), &
!                        " % complete"
!                endif

                ! save p at [it-1] to pp and [it-2] to ppp
                ppp = pp
                pp = p

                !$OMK PARALLEL NUM_THREADS(OMP_GET_MAX_THREADS()) DEFAULT(NONE) &
                !$OMK PRIVATE(ix,iz) &
                !$OMK SHARED(dvxdx, dvzdz, p,deltax24I,deltaz24I,a_x,a_z,memory_dvz_dz,modttI,&
                !$OMK memory_dvx_dx,K_xI,K_zI,vx,vz,tim_del,modrrvx,boundary_p,nx,nz,b_x,b_z)
                !$OMK DO
                !forall(ix=2:nx,iz=2:nz)
                do ipropwav = 1,npropwav
                do ix=2,nx
                do iz=2,nz
                        dvxdx(iz,ix,ipropwav) = &
                        (27.e0*p(iz,ix,2,ipropwav)-&
                        27.e0*p(iz,ix-1,2,ipropwav)-&
                            p(iz,ix+1,2,ipropwav)+&
                            p(iz,ix-2,2,ipropwav)) * (deltax24I)
                        dvzdz(iz,ix,ipropwav) = &
                        (27.e0*p(iz,ix,3,ipropwav)-&
                        27.e0*p(iz-1,ix,3,ipropwav)-&
                            p(iz+1,ix,3,ipropwav)+&
                            p(iz-2,ix,3,ipropwav)) * (deltaz24I)
                        memory_dvx_dx(iz,ix,ipropwav) = b_x(ix) * memory_dvx_dx(iz,ix,ipropwav) + &
                                                        a_x(ix) * dvxdx(iz,ix,ipropwav) ! pml 
                        memory_dvz_dz(iz,ix,ipropwav) = b_z(iz) * memory_dvz_dz(iz,ix,ipropwav) + &
                                                        a_z(iz) * dvzdz(iz,ix,ipropwav) ! pml
                        dvxdx(iz,ix,ipropwav) = dvxdx(iz,ix,ipropwav) * k_xI(ix) + &
                                                memory_dvx_dx(iz,ix,ipropwav) !pml
                        dvzdz(iz,ix,ipropwav) = dvzdz(iz,ix,ipropwav) * k_zI(iz) + &
                                                memory_dvz_dz(iz,ix,ipropwav) !pml

                        ! compute pressure at [it] using p at [it-1] and dvxdx
                        ! and dvzdz at [it-1/2]
                        p(iz,ix,1,ipropwav) = p(iz,ix,1,ipropwav) + &
                                (modttI(iz,ix) * (dvxdx(iz,ix,ipropwav) + &
                                                dvzdz(iz,ix,ipropwav))) * &
                                                        tim_del * boundary_p(iz,ix)
                enddo
                enddo
                enddo
                !endforall
                !$OMK ENDDO
                !$OMK END PARALLEL
                !$OMK END PARALLEL

                ! adding source to pressure at [it] 
                do ipropwav = 1,npropwav
                do ifield = 1, src_nfield
                do issmul = 1, src_nsmul
                        ! compute source wavelet at [it-1/2] by simple
                        ! interpolation
                        ! division of source term with deltax and deltaz (see Jan's fdelmodc manual)
                        svalue  = 0.5e0 * (src_wav_mat(issmul,ifield,ipropwav,it, isseq) + &
                                src_wav_mat(issmul,ifield,ipropwav,it-1,isseq)) * tim_del &
                                        * deltax24I * deltaz24I


                                p(issmulz1(issmul,ipropwav), issmulx1(issmul,ipropwav),ifield, ipropwav) = &
                        p(issmulz1(issmul,ipropwav),issmulx1(issmul,ipropwav),ifield, ipropwav) &
                                + svalue * src_spray_weights(1,issmul,ipropwav) * &
                                modttI(issmulz1(issmul,ipropwav), issmulx1(issmul,ipropwav))  
                        p(issmulz1(issmul,ipropwav), issmulx2(issmul,ipropwav),ifield, ipropwav) = &
                                p(issmulz1(issmul,ipropwav),issmulx2(issmul,ipropwav),ifield, ipropwav) &
                                + svalue * src_spray_weights(2,issmul,ipropwav) * &
                                modttI(issmulz1(issmul,ipropwav), issmulx2(issmul,ipropwav))
                        p(issmulz2(issmul,ipropwav), issmulx1(issmul,ipropwav),ifield, ipropwav) = &
                                p(issmulz2(issmul,ipropwav),issmulx1(issmul,ipropwav),ifield, ipropwav) &
                                + svalue * src_spray_weights(3,issmul,ipropwav) * &
                                modttI(issmulz2(issmul,ipropwav), issmulx1(issmul,ipropwav))
                        p(issmulz2(issmul,ipropwav), issmulx2(issmul,ipropwav),ifield, ipropwav) = &
                                p(issmulz2(issmul,ipropwav),issmulx2(issmul,ipropwav),ifield, ipropwav) &
                                + svalue * src_spray_weights(4,issmul,ipropwav) * &
                                modttI(issmulz2(issmul,ipropwav), issmulx2(issmul,ipropwav))
                enddo
                enddo
                enddo

                ! force p[1] on borders
                if(border_in_flag) then
                do iborder = 1, border_n
                        p(iborderz(iborder), iborderx(iborder),1,1) = &
                        border_in(iborder,it,isseq) 
                enddo
                endif

                ! compute dpdx and dpdz at [it] only for the first propagating field, to use for born modeling
                !forall(ix=1:nx-1,iz=1:nz-1)
                ipropwav = 1
                do ix=1,nx-1
                do iz=1,nz-1
                        dpdx(iz,ix,ipropwav) = &
                                (27.e0*p(iz,ix+1,1,ipropwav)-&
                                27.e0*p(iz,ix,1,ipropwav)-&
                                    p(iz,ix+2,1,ipropwav)+&
                                    p(iz,ix-1,1,ipropwav)) * (deltax24I)
                        memory_dp_dx(iz,ix,ipropwav) = b_x_half(ix) * memory_dp_dx(iz,ix,ipropwav) + &
                                                    a_x_half(ix) * dpdx(iz,ix,ipropwav) ! pml
                        dpdx(iz,ix,ipropwav) = dpdx(iz,ix,ipropwav) * K_x_halfI(ix) + &
                                            memory_dp_dx(iz,ix,ipropwav) ! pml
                enddo
                enddo
                !endforall

                !forall(ix=1:nx-1,iz=1:nz-1)
                do ix=1,nx-1
                do iz=1,nz-1
                        dpdz(iz,ix,ipropwav) = &
                                (27.e0*p(iz+1,ix,1,ipropwav)-&
                                27.e0*p(iz,ix,1,ipropwav)-&
                                    p(iz+2,ix,1,ipropwav)+&
                                    p(iz-1,ix,1,ipropwav)) * (deltaz24I)
                        memory_dp_dz(iz,ix,ipropwav) = b_z_half(iz) * memory_dp_dz(iz,ix,ipropwav) + &
                                                    a_z_half(iz) * dpdz(iz,ix,ipropwav) !pml
                        dpdz(iz,ix,ipropwav) = dpdz(iz,ix,ipropwav) * K_z_halfI(iz) + &
                                            memory_dp_dz(iz,ix,ipropwav) !pml
                enddo
                enddo


                if(born_flag) then
                ! adding born sources from pressure(:,:,1) to pressure(:,:,2)
                do ix=1,nx-1
                do iz=1,nz-1
                        !p(iz,ix,2) = p(iz,ix,1);
                        p(iz,ix,1,2) = p(iz,ix,1,2) + & 
                                ! lambdaI scatterrer
                                (-1.0 * &
                                (ppp(iz, ix, 1,1) + p(iz, ix,  1,1) - 2.0 * pp(iz, ix,  1,1)) &
                                                 * deltamodtt(iz, ix) * tim_delI * tim_delI) + &
                                ! modrrvx scatterrer
                                (27.e0*dpdx(iz,ix,1) * deltamodrrvx(iz,ix) &
                                        -27.e0*dpdx(iz,ix-1,1) * deltamodrrvx(iz,ix-1) &
                                        -dpdx(iz,ix+1,1) * deltamodrrvx(iz,ix+1) &
                                        +dpdx(iz,ix-2,1) * deltamodrrvx(iz,ix-2) &
                                        ) * (deltax24I) + &
                                ! modrrvz scatterrer
                                (27.e0*dpdz(iz,ix,1) * deltamodrrvz(iz,ix) &
                                        -27.e0*dpdz(iz-1,ix,1) * deltamodrrvz(iz-1,ix) &
                                        -dpdz(iz+1,ix,1) * deltamodrrvz(iz+1,ix) &
                                        +dpdz(iz-2,ix,1) * deltamodrrvz(iz-2,ix) &
                                        ) * (deltaz24I)
        
                enddo
                enddo
                endif


                ! compute dpdx and dpdz at [it] for all other propagating fields
                !forall(ix=1:nx-1,iz=1:nz-1)
                do ipropwav = 2, npropwav
                do ix=1,nx-1
                do iz=1,nz-1
                        dpdx(iz,ix,ipropwav) = &
                                (27.e0*p(iz,ix+1,1,ipropwav)-&
                                27.e0*p(iz,ix,1,ipropwav)-&
                                    p(iz,ix+2,1,ipropwav)+&
                                    p(iz,ix-1,1,ipropwav)) * (deltax24I)
                        memory_dp_dx(iz,ix,ipropwav) = b_x_half(ix) * memory_dp_dx(iz,ix,ipropwav) + &
                                                    a_x_half(ix) * dpdx(iz,ix,ipropwav) ! pml
                        dpdx(iz,ix,ipropwav) = dpdx(iz,ix,ipropwav) * K_x_halfI(ix) + &
                                            memory_dp_dx(iz,ix,ipropwav) ! pml
                enddo
                enddo
                !endforall

                !forall(ix=1:nx-1,iz=1:nz-1)
                do ix=1,nx-1
                do iz=1,nz-1
                        dpdz(iz,ix,ipropwav) = &
                                (27.e0*p(iz+1,ix,1,ipropwav)-&
                                27.e0*p(iz,ix,1,ipropwav)-&
                                    p(iz+2,ix,1,ipropwav)+&
                                    p(iz-1,ix,1,ipropwav)) * (deltaz24I)
                        memory_dp_dz(iz,ix,ipropwav) = b_z_half(iz) * memory_dp_dz(iz,ix,ipropwav) + &
                                                    a_z_half(iz) * dpdz(iz,ix,ipropwav) !pml
                        dpdz(iz,ix,ipropwav) = dpdz(iz,ix,ipropwav) * K_z_halfI(iz) + &
                                            memory_dp_dz(iz,ix,ipropwav) !pml
                enddo
                enddo
                enddo
                !endforall



                ! gradient calculation
                if(grad_out_flag) then

                        ! p at [it], pp at [it-1]
                        ! dpdx and dpdz at [it]
                        ! gradients w.r.t inverse of modttI, i.e., 1/rho/c^2 
                        gradis_modtt = gradis_modtt + (p(:,:,1,1) - pp(:,:,1,1)) * tim_delI * &
                                        (p(:,:,1,2) - pp(:,:,1,2)) * tim_delI
                        

                        ! gradient w.r.t. inverse of rho on vx and vz grids
                        gradis_modrrvx = gradis_modrrvx - dpdx(:,:,2)*dpdx(:,:,1) 
                        gradis_modrrvz = gradis_modrrvz - dpdz(:,:,2)*dpdz(:,:,1) 

                endif

                ! saving a field snapshot (commented for speed)
                !if(allocated(saved_snap) ) then
                !        if(it .eq. int(tim_nt * snap_save_at)) then
                !                saved_snap(:,:) = p(1:mesh_nz+mesh_na_pml,1+mesh_na_pml:mesh_nx+mesh_na_pml) ;
                !        endif
                !endif

                ! update velocity at [it+1/2] using 
                ! velcity at [it-1/2] and dpdx and dpdz at [it] 
                !forall(ix=1:nx-1,iz=1:nz-1)
                do ipropwav = 1, npropwav
                do ix=1,nx-1
                do iz=1,nz-1
                        p(iz,ix,2,ipropwav) = p(iz,ix,2,ipropwav) + (dpdx(iz,ix,ipropwav)) * tim_del * &
                                modrrvx(iz,ix) * boundary_vx(iz,ix)
                enddo
                enddo
                !endforall

                !forall(ix=1:nx-1,iz=1:nz-1)
                do ix=1,nx-1
                do iz=1,nz-1
                        p(iz,ix,3,ipropwav) = p(iz,ix,3,ipropwav) + (dpdz(iz,ix,ipropwav)) * tim_del * &
                                modrrvz(iz,ix) * boundary_vz(iz,ix)
                enddo
                enddo
                enddo
                !endforall

                ! store pressure seismograms to rec_mat 
                do ipropwav = 1, npropwav
                do ifield = 1, recv_nfield
                do irec = 1, recv_n
!                                if(index(recv_flags, '[P]') .ne. 0) then
                                rec_mat(irec,ifield,ipropwav, it) = &
                                (&
                                p(irecz1(irec,ipropwav),irecx1(irec,ipropwav),ifield,ipropwav)*&
                                rec_interpolate_weights(1,irec,ipropwav)+&
                                p(irecz1(irec,ipropwav),irecx2(irec,ipropwav),ifield,ipropwav)*&
                                rec_interpolate_weights(2,irec,ipropwav)+&
                                p(irecz2(irec,ipropwav),irecx1(irec,ipropwav),ifield,ipropwav)*&
                                rec_interpolate_weights(3,irec,ipropwav)+&
                                p(irecz2(irec,ipropwav),irecx2(irec,ipropwav),ifield,ipropwav)*&
                                rec_interpolate_weights(4,irec,ipropwav)&
                                )
!                               endif
                enddo
                enddo
                enddo


                ! record borders
                if(border_out_flag) then
                !do ipropwav = 1, npropwav
                !do ifield = 1, 3
                do iborder = 1, border_n
                        border_out(iborder,it,isseq) = &
                        p(iborderz(iborder),iborderx(iborder),1,1)
                enddo
                !enddo
                !enddo
                endif

        enddo time_loop

        do ipropwav = 1, npropwav
        do ifield=1, recv_nfield
        do irec = 1, recv_n
        do it = 1, tim_nt 
                ! output recv_out
                recv_out(it + & 
                         (irec-1)*tim_nt + &
                         (isseq-1)*recv_n*tim_nt + &
                         (ifield-1)*recv_n*tim_nt*src_nseq + &
                         (ipropwav-1)*recv_nfield*recv_n*tim_nt*src_nseq &
                        ) = &
                        rec_mat(irec,ifield,ipropwav,it)
        enddo
        enddo
        enddo
        enddo

        ! output p_out of the first propagating field
        p_out(:,:,:,isseq,1) = p(:,:,:,1) 
        p_out(:,:,:,isseq,2) = pp(:,:,:,1)


        ! edges of the gradient zero
        !if(present(grad)) then
        !        ! grad of vel
        !        grad(isseq)%vel = rzero_de;
        !        forall(ix=2:mesh_nx-2,iz=2:mesh_nz-2)
        !                grad(isseq)%vel(iz,ix) = -2.0e0 * &
        !                        (var%vel(iz,ix))**(-3e0) * var%rho(iz,ix)**(-1e0) * grad_modtt(iz,ix)
        !        endforall
        !

        ! output gradient
        if(grad_out_flag) then
                grad_modtt(:,:,isseq) = gradis_modtt(1:nz,1:nx);

                forall(ix=2:nx-2,iz=2:nz-2)
                        grad_modrr(iz,ix, isseq) = grad_modrr(iz,ix,isseq) + 0.5e0 * gradis_modrrvx(iz,ix)
                        grad_modrr(iz,ix+1, isseq) = grad_modrr(iz,ix+1,isseq) + 0.5e0 * gradis_modrrvx(iz,ix)
                        grad_modrr(iz,ix, isseq) = grad_modrr(iz,ix,isseq) + 0.5e0 * gradis_modrrvz(iz,ix)
                        grad_modrr(iz+1,ix, isseq) = grad_modrr(iz+1,ix,isseq) + 0.5e0 * gradis_modrrvz(iz,ix)
                endforall
        endif


        !        grad(isseq)%rho = rzero_de;
        !        forall(ix=2:mesh_nx-2,iz=2:mesh_nz-2)
        !                grad(isseq)%rho(iz,ix) = -1.e0 * (var%rho(iz,ix))**(-2e0) * grad_modrr(iz,ix) &
        !                        - (var%vel(iz,ix) * var%rho(iz,ix))**(-2e0) * grad_modtt(iz,ix)
        !        endforall
        !endif


        deallocate                      (p, pp, ppp&
                                        )
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
        deallocate(irecx1, irecx2, irecz1, irecz2)
        deallocate(denomrecI, rec_interpolate_weights)

        if(allocated(rec_mat)) deallocate  (rec_mat)
        if(allocated(gradis_modrrvx)) deallocate(gradis_modrrvx)
        if(allocated(gradis_modrrvz)) deallocate(gradis_modrrvz)
        if(allocated(gradis_modtt)) deallocate(gradis_modtt)


enddo src_par_loop
!$OMP ENDDO 
!$OMP END PARALLEL

deallocate                      (iborderx, iborderz)
deallocate                      (d_x,K_x,K_xI,alpha_x,&
                                a_x,b_x,d_x_half,&
                                K_x_half,K_x_halfI,alpha_x_half,a_x_half,b_x_half)

deallocate                      (d_z,K_z,K_zI,alpha_z,&
                                a_z,b_z,d_z_half,&
                                K_z_half,K_z_halfI,alpha_z_half,a_z_half,b_z_half)
deallocate                      (modttI,  modrrvx, modrrvz, deltamodtt, deltamodrrvx, deltamodrrvz)

deallocate(maskX); deallocate(maskZ)
if(allocated(src_wav_mat)) deallocate(src_wav_mat)

! testing the speed -- uncomment if u want to print out cpu time
cpu_time2  = omp_get_wtime(); write(*,*) "fd2_mod: ",(cpu_time2-cpu_time1),"seconds"

end subroutine fdtd_mod

end module fdtd
!\end{comment}
