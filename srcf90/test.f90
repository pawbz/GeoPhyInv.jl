module test
use precision_mod! USEINDEPFILE
use string_routines! USEINDEPFILE
use iso_c_binding
implicit none

contains

! subroutine simple print 
subroutine sprint(sayhello) bind(c, name="sprint")
        implicit none
        character(len=1, kind=C_char), intent(in)            :: sayhello(*)

        write(*,*) repeat(" ",10),"sprint success in test.f90 "
        write(*,*) repeat(" ",10),"Hello from ", ctofstr(sayhello)
end subroutine sprint

! subroutine simple print using many cores
subroutine sprint_core() bind(c, name="sprint_core")
        use omp_lib

        implicit none

        !$OMP PARALLEL
        write(*,*) repeat(" ",10),"sprint success core", OMP_GET_THREAD_NUM(), "in test.f90"
        !$OMP END PARALLEL

end subroutine sprint_core

subroutine matrix_in_out(A, B, n1, n2) bind(c, name="matrix_in_out")
        implicit  none
        real(r_dp), intent(in)              :: A(n1,n2)
        real(r_dp), intent(inout)           :: B(n1)
        integer(i_sp), intent(in)           :: n1, n2

        write(*,*) repeat(" ",10),"in fortran n1 n2", n1, n2
        write(*,*) repeat(" ",10),"in fortran A", size(A)
        write(*,*) repeat(" ",10),"in fortran A", A
        write(*,*) repeat(" ",10),"in fortran B", size(B)
        write(*,*) repeat(" ",10),"in fortran B", B
        B(1) = 0.0

end subroutine matrix_in_out

!subroutine derived_data_type(model_in, model_out) bind(c, name="derived_data_type")
!        implicit none
!        type(model), intent(in)         :: model_in
!        type(model), intent(inout)        :: model_out
!        write(*,*) repeat(" ",10),"in fortran model%A", model_in%A
!        write(*,*) repeat(" ",10),"in fortran model%B", model_in%B
!end subroutine derived_data_type

end module test
