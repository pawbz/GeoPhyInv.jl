module sorting
implicit none
private

interface swap
        module procedure swapreal, swapint
end interface swap
interface findminimum
        module procedure findminimumreal, findminimumint
end interface findminimum
interface sort
        module procedure sortreal, sortint
end interface sort



public                                                  :: sort

contains

! --------------------------------------------------------------------
! integer function  findminimum():
!    this function returns the location of the minimum in the section
! between start and end.
! --------------------------------------------------------------------

integer function  findminimumint(x, start, end)
implicit  none
integer, dimension(1:), intent(in) :: x
integer, intent(in)                :: start, end
integer                            :: minimum
integer                            :: location
integer                            :: i

minimum  = x(start)             ! assume the first is the min
location = start                        ! record its position
do i = start+1, end             ! start with next elements
 if (x(i) < minimum) then       !   if x(i) less than the min?
    minimum  = x(i)             !      yes, a new minimum found
    location = i                !      record its position
 end if
end do
findminimumint = location          ! return the position
end function  findminimumint

integer function  findminimumreal(x, start, end)
implicit  none
real, dimension(1:), intent(in)    :: x
integer, intent(in)                :: start, end
real                               :: minimum
integer                            :: location
integer                            :: i

minimum  = x(start)             ! assume the first is the min
location = start                        ! record its position
do i = start+1, end             ! start with next elements
 if (x(i) < minimum) then       !   if x(i) less than the min?
    minimum  = x(i)             !      yes, a new minimum found
    location = i                !      record its position
 end if
end do
findminimumreal = location          ! return the position
end function  findminimumreal



! --------------------------------------------------------------------
! subroutine  swap():
!    this subroutine swaps the values of its two formal arguments.
! --------------------------------------------------------------------

subroutine  swapint(a, b)
implicit  none
integer, intent(inout) :: a, b
integer                :: temp
temp = a
a    = b
b    = temp
end subroutine  swapint
subroutine  swapreal(a, b)
implicit  none
real, intent(inout)     :: a, b
real                    :: temp
temp = a
a    = b
b    = temp
end subroutine  swapreal


! --------------------------------------------------------------------
! subroutine  sort():
!    this subroutine receives an array x() and sorts it into ascending
! order.
! --------------------------------------------------------------------

subroutine  sortint(x, size)
implicit  none
integer, dimension(1:), intent(inout) :: x
integer, intent(in)                   :: size
integer                               :: i
integer                               :: location

do i = 1, size-1                        ! except for the last
        location = findminimum(x, i, size)     ! find min from this to last
        call  swap(x(i), x(location))  ! swap this and the minimum
end do
end subroutine  sortint

subroutine  sortreal(x, size)
implicit  none
real, dimension(1:), intent(inout)    :: x
integer, intent(in)                   :: size
integer                               :: i
integer                               :: location

do i = 1, size-1                        ! except for the last
        location = findminimum(x, i, size)     ! find min from this to last
        call  swap(x(i), x(location))  ! swap this and the minimum
end do
end subroutine  sortreal


end module sorting

