module test
use precision_mod! USEINDEPFILE

        type, bind(c) :: model
                real            :: m1, m2, m3
        endtype
contains

! subroutine simple print 
subroutine sprint() bind(c, name="sprint")
        implicit none

        write(*,*) "sprint success in test.f90"

end subroutine sprint

! subroutine simple print using many cores
subroutine sprint_core() bind(c, name="sprint_core")
        use omp_lib

        implicit none

        !$OMP PARALLEL
        write(*,*) "sprint success core", OMP_GET_THREAD_NUM(), "in test.f90"
        !$OMP END PARALLEL

end subroutine sprint_core

subroutine matrix_in_out(A, B) bind(c, name="matrix_in_out")
        implicit  none
        real(r_dp), intent(in)           :: A(:)
        real(r_dp), intent(inout)           :: B(:)

        write(*,*) "in fortran ", size(A)

end subroutine matrix_in_out

subroutine derived_data_type(model_in, model_out) bind(c, name="derived_data_type")
        implicit none
        type(model), intent(in)         :: model_in
        type(model), intent(out)        :: model_out
end subroutine derived_data_type

end module test
