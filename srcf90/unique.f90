module unique
use error_messages! USEINDEPFILE
implicit none
private

public                                          :: unique_mask 

interface unique_mask
        module procedure unique_maskreal1, unique_maskreal2, unique_maskinteger1, unique_maskinteger2
end interface unique_mask

contains

subroutine unique_maskreal1(array1, unique_mask)
implicit none
real, intent(in)                        :: array1(:)
logical, intent(out)                    :: unique_mask(:)
integer                                 :: nz, iz
nz = size(array1)
call check_dimension('unique_maskreal1: output size', unique_mask, nz)
unique_mask = .true.
do iz = nz,2,-1
   unique_mask(iz) = .not.(any(array1(:iz-1) .eq. array1(iz)))
end do
end subroutine unique_maskreal1

subroutine unique_maskreal2(array1, array2, unique_mask)
implicit none
real, intent(in)                        :: array1(:), array2(:)
logical, intent(out)                    :: unique_mask(:)
integer                                 :: nz, iz
nz = size(array1)
call check_dimension('unique_maskreal2: size if array2', array2, nz)
call check_dimension('unique_maskreal2: output size', unique_mask, nz)
unique_mask = .true.
do iz = nz,2,-1
        unique_mask(iz) = .not.(any(array1(:iz-1)==array1(iz).and.array2(:iz-1)==array2(iz)))
end do
end subroutine unique_maskreal2

subroutine unique_maskinteger1(array1, unique_mask)
implicit none
integer, intent(in)                     :: array1(:)
logical, intent(out)                    :: unique_mask(:)
integer                                 :: nz, iz
nz = size(array1)
call check_dimension('unique_maskinteger1: output size', unique_mask, nz)
unique_mask = .true.
do iz = nz,2,-1
   unique_mask(iz) = .not.(any(array1(:iz-1) .eq. array1(iz)))
end do
end subroutine unique_maskinteger1

subroutine unique_maskinteger2(array1, array2, unique_mask)
implicit none
integer, intent(in)                     :: array1(:), array2(:)
logical, intent(out)                    :: unique_mask(:)
integer                                 :: nz, iz
nz = size(array1)
call check_dimension('unique_maskinteger2: size if array2', array2, nz)
call check_dimension('unique_maskinteger2: output size', unique_mask, nz)
unique_mask = .true.
do iz = nz,2,-1
        unique_mask(iz) = .not.(any(array1(:iz-1)==array1(iz).and.array2(:iz-1)==array2(iz)))
end do
end subroutine unique_maskinteger2



end module unique
