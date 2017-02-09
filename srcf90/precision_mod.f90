module precision_mod
implicit none
public

integer, parameter                             :: r_sp = selected_real_kind(p=6,r=37) ! single precision real
integer, parameter                             :: r_dp = selected_real_kind(p=12) ! double precision real
integer, parameter                             :: i_sp = selected_int_kind(3) ! single precision int
integer, parameter                             :: i_dp = selected_int_kind(9) ! double precision int
real, parameter                                :: rzero_de = real(0.0)
real, parameter                                :: rhalf_de = real(0.5)
real, parameter                                :: rtwo_de  = real(2.0)
real, parameter                                :: rone_de  = real(1.0)


integer, parameter                             :: izero_de = int(0.0)
integer, parameter                             :: itwo_de  = int(2.0)
integer, parameter                             :: ione_de  = int(1.0)



integer, parameter                             :: r_de = kind(rzero_de) ! default real kind
integer, parameter                             :: i_de = kind(izero_de) ! default integer kind



real(r_dp), parameter                          :: rzero_dp = real(0,r_dp)
real(r_dp), parameter                          :: rhalf_dp = real(0.5,r_dp)
real(r_dp), parameter                          :: rtwo_dp  = real(2.0,r_dp)
real(r_dp), parameter                          :: rone_dp  = real(1.0,r_dp)
real(r_sp), parameter                          :: rzero_sp = real(0,r_sp)
real(r_sp), parameter                          :: rhalf_sp = real(0.5,r_sp)
real(r_sp), parameter                          :: rtwo_sp  = real(2.0,r_sp)
real(r_sp), parameter                          :: rone_sp  = real(1.0,r_sp)


integer(i_dp), parameter                       :: izero_dp = int(0,i_dp)
integer(i_dp), parameter                       :: itwo_dp  = int(2.0,i_dp)
integer(i_dp), parameter                       :: ione_dp  = int(1.0,i_dp)
integer(i_sp), parameter                       :: izero_sp = int(0,i_sp)
integer(i_sp), parameter                       :: itwo_sp  = int(2.0,i_sp)
integer(i_sp), parameter                       :: ione_sp  = int(1.0,i_sp)



integer, parameter                             :: c_sp = kind((rone_sp,rone_sp))
integer, parameter                             :: c_dp = kind((rone_dp,rone_dp))
complex, parameter                             :: czero_de = cmplx(rzero_de,rzero_de)
complex, parameter                             :: cone_de = cmplx(rone_de,rzero_de)
complex(c_dp), parameter                       :: czero_dp = cmplx(rzero_dp,rzero_dp,kind=c_dp)
complex(c_sp), parameter                       :: czero_sp = cmplx(rzero_sp,rzero_sp,kind=c_sp)
complex(c_dp), parameter                       :: cone_dp = cmplx(rone_dp,rzero_dp,kind=c_dp)
complex(c_sp), parameter                       :: cone_sp = cmplx(rone_sp,rzero_sp,kind=c_sp)

integer, parameter                             :: c_de = kind(czero_de) ! default complex kind

! maximum lengths of string variables used throughout the program
! length of string
integer, parameter                                                      :: max_str_len=400
! max length of first dimension
integer, parameter                                                      :: max_str_dim=200

end module precision_mod


