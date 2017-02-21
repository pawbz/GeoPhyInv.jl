
module staggered_grid

use precision_mod! USEINDEPFILE
use error_messages! USEINDEPFILE
public

contains


subroutine get_musigmayx(mu, musigmayx)
implicit none
real, intent(in)                                :: mu(:,:)
real, intent(out)                               :: musigmayx(:,:)
integer                                         :: iz, ix 
call check_dimension("fd1sh: get_musigmayx", mu, musigmayx)
musigmayx = rzero_de;
do ix = 1, size(mu, 2)-1
        do iz = 1, size(mu, 1)
                if(mu(iz,ix+1) * mu(iz,ix) .ne. rzero_de) then
                        musigmayx(iz, ix) = 0.5e0 * (mu(iz,ix+1) + mu(iz,ix))
                else
                        musigmayx(iz, ix) = rzero_de
                endif
        enddo
enddo
end subroutine get_musigmayx

subroutine get_musigmayz(mu, musigmayz)
implicit none
real, intent(in)                                :: mu(:,:)
real, intent(out)                               :: musigmayz(:,:)
integer                                         :: iz, ix 
call check_dimension("fd1sh: get_musigmayz", mu, musigmayz)
musigmayz = rzero_de;

do ix = 1, size(mu, 2)
        do iz = 1, size(mu, 1)-1
                if(mu(iz+1,ix) * mu(iz,ix) .ne. rzero_de) then
                        musigmayz(iz, ix) =  0.5e0 * (mu(iz+1,ix) + mu(iz,ix))
                else
                        musigmayz(iz, ix) = rzero_de
                endif
        enddo
enddo
end subroutine get_musigmayz


function extend(input, np, nxprev, nzprev)
! this subroutine extends the input model on all four sides by np points
!       input                           : input (nxprev , nzprev) 
!       np                              : number of pml layers to be added to extend models
!       nxprev                          : size in x direction before extending
!       nzprev                          : size in z direction before extending
!       
! output 
!       extend                          : extended model of diemension (0:(nxprev + 2* np) +1 , 0:(nzprev + 2*np)+1)

implicit none
integer, intent(in)                     :: np, nxprev, nzprev
real, intent(in), dimension(:,:)        :: input
real                                    :: extend(0: (nzprev + 2*np) + 1, 0: (nxprev + 2*np) + 1)
integer                                 :: ix,iz

if(size(input,1) .ne. nzprev) stop "extend: size mismatch, nz" 
if(size(input,2) .ne. nxprev) stop "extend: size mismatch, nx" 
! initialize
extend(:,:)=rzero_de;

! padding 
extend(np+1:nzprev+np,np+1:nxprev+np)=input;
do ix=0,nxprev + 2*np +1
        do iz=0,np
                extend(iz,ix)=extend(np+1,ix)
                extend(nzprev +2*np-iz+1,ix)=extend(nzprev+np,ix)
        end do
end do
do iz=0,nzprev+2*np+1
        do ix=0,np
                extend(iz,ix)=extend(iz,np+1)
                extend(iz,nxprev+2*np+1-ix)=extend(iz,nxprev+np)
        end do
end do
end function extend

subroutine get_rhoI(rho, rhoI)
implicit none
real, intent(in)                                :: rho(:,:)
real, intent(out)                               :: rhoI(:,:)
integer                                         :: iz, ix 
call check_dimension("fd1acou: get_rhoI", rho, rhoI)
rhoI = rzero_de;
do ix = 1, size(rho, 2)-1
        do iz = 1, size(rho, 1)
                if(rho(iz,ix) .ne. rzero_de) then
                        rhoI(iz,ix) = (rho(iz,ix))**(-1.e0)
                else
                        call abort_msg("get_rhoI: cannot be zero")
                endif
        enddo
enddo
end subroutine get_rhoI


subroutine get_rhovxI(rhoI, rhovxI)
implicit none
real, intent(in)                                :: rhoI(:,:)
real, intent(out)                               :: rhovxI(:,:)
integer                                         :: iz, ix 
call check_dimension("fd1acou: get_rhovxI", rhoI, rhovxI)
rhovxI = rzero_de;
do ix = 1, size(rhoI, 2)-1
        do iz = 1, size(rhoI, 1)
                rhovxI(iz, ix) = 0.5e0 *(rhoI(iz,ix+1) + rhoI(iz,ix))
        enddo
enddo
end subroutine get_rhovxI

subroutine get_rhovzI(rhoI, rhovzI)
implicit none
real, intent(in)                                :: rhoI(:,:)
real, intent(out)                               :: rhovzI(:,:)
integer                                         :: iz, ix 
call check_dimension("fd1acou: get_rhovzI", rhoI, rhovzI)
rhovzI = rzero_de;

do ix = 1, size(rhoI, 2)
        do iz = 1, size(rhoI, 1)-1
                        rhovzI(iz, ix) =  0.5e0 * (rhoI(iz+1,ix) + rhoI(iz,ix))
        enddo
enddo
end subroutine get_rhovzI

end module staggered_grid
