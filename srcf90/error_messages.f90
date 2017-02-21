module error_messages
use precision_mod! USEINDEPFILE
use string_routines! USEINDEPFILE
implicit none

interface check_dimension
        module procedure         check_sizereal1D, check_sizereal2D, check_sizereal3D, check_sizereal4D
        module procedure         check_sizeint1D, check_sizeint2D, check_sizeint3D
        module procedure         check_sizelogical1D
        module procedure         check_sizecomplex1D
        module procedure         check_size_equalreal1D, check_size_equalreal2D, check_size_equalreal3D
        module procedure         check_size_equalrealint1D
end interface check_dimension

interface is_nan
        module procedure         is_nanint, is_nanint1D, is_nanint2D, is_nanint3D
        module procedure         is_nanreal, is_nanreal1D, is_nanreal2D, is_nanreal3D
        module procedure         is_nancomplex2D
end interface is_nan 

contains

subroutine abort_msg(message)
character(len=*), intent(in)                    :: message
write(*,*) trim(message)
write(*,*) 'program aborted with internal error'
stop
end subroutine abort_msg


subroutine keyword_input_error(message,keyword1)
implicit none
character(len=*), intent(in)                    :: keyword1, message
write(*,*) 'error input keyword: ',trim(keyword1) 
write(*,*) '     ',trim(message)
stop 'program aborted'
end subroutine keyword_input_error

subroutine keyword_input_alert(message,keyword1)
implicit none
character(len=*), intent(in)                    :: keyword1, message
write(*,*) 'alert input keyword: ',trim(keyword1) 
write(*,*) '     ',trim(message)
end subroutine keyword_input_alert



subroutine keyword_input_required(keyword1, keyword2, keyword3, keyword4, keyword5, keyword6, keyword7, keyword8)
implicit none
character(len=*), intent(in)                    :: keyword1
character(len=*), optional, intent(in)          :: keyword2, keyword3, keyword4, keyword5, keyword6, keyword7, &
                                                        keyword8
integer                                         :: nn

nn = 0
write(*,*) 'one of the following input keywords required:'
        
if(trim(keyword1) .ne. "") then
        nn = nn + 1
        write(*,*) '     '//trim(i2str(nn))//') ', trim(keyword1)      
endif
if(present(keyword2)) then
        if(trim(keyword2) .ne. "") then
                nn = nn + 1
                write(*,*) '     '//trim(i2str(nn))//') ', trim(keyword2) 
        endif
endif
if(present(keyword3)) then
        if(trim(keyword3) .ne. "") then
                nn = nn + 1
                write(*,*) '     '//trim(i2str(nn))//') ', trim(keyword3)    
        endif
endif
if(present(keyword4)) then
        if(trim(keyword4) .ne. "") then
                nn = nn + 1
                write(*,*) '     '//trim(i2str(nn))//') ', trim(keyword4)    
        endif
endif
if(present(keyword5)) then
        if(trim(keyword5) .ne. "") then
                nn = nn + 1
                write(*,*) '     '//trim(i2str(nn))//') ', trim(keyword5)    
        endif
endif
if(present(keyword6)) then
        if(trim(keyword6) .ne. "") then
                nn = nn + 1
                write(*,*) '     '//trim(i2str(nn))//') ', trim(keyword6)    
        endif
endif
if(present(keyword7)) then
        if(trim(keyword7) .ne. "") then
                nn = nn + 1
                write(*,*) '     '//trim(i2str(nn))//') ', trim(keyword7)    
        endif
endif
if(present(keyword8)) then
        if(trim(keyword8) .ne. "") then
                nn = nn + 1
                write(*,*) '     '//trim(i2str(nn))//') ', trim(keyword8)    
        endif
endif
stop 'program aborted'
end subroutine keyword_input_required

subroutine check_sizelogical1D(message, x, n1)
implicit none
logical, intent(in)                     :: x(:)
integer, intent(in)                     :: n1
character(len=*)                        :: message

if(size(x,1) .ne. n1) then
        write(*,*) message
        write(*,*) 'size mismatch'
        write(*,*) 'shape of x: ', shape(x)
        write(*,*) 'n1: ',n1 
        stop 'program aborted with internal error'
endif
end subroutine check_sizelogical1D

subroutine check_sizecomplex1D(message, x, n1)
implicit none
complex, intent(in)                     :: x(:)
integer, intent(in)                     :: n1
character(len=*)                        :: message

if(size(x,1) .ne. n1) then
        write(*,*) message
        write(*,*) 'size mismatch'
        write(*,*) 'shape of x: ', shape(x)
        write(*,*) 'n1: ',n1 
        stop 'program aborted with internal error'
        
endif
end subroutine check_sizecomplex1D

subroutine check_sizereal1D(message, x, n1)
implicit none
real, intent(in)                        :: x(:)
integer, intent(in)                     :: n1
character(len=*)                        :: message

if(size(x,1) .ne. n1) then
        write(*,*) message
        write(*,*) 'size mismatch'
        write(*,*) 'shape of x: ', shape(x)
        write(*,*) 'n1: ',n1 
        stop 'program aborted with internal error'
        
endif
end subroutine check_sizereal1D

subroutine check_sizereal2D(message, x, n1, n2)
implicit none
real, intent(in)                        :: x(:,:)
integer, intent(in)                     :: n1, n2
character(len=*)                        :: message

if(size(x,1) .ne. n1 .or. size(x,2) .ne. n2 ) then
        write(*,*) message
        write(*,*) 'size mismatch'
        write(*,*) 'shape of x: ', shape(x)
        write(*,*) 'n1,n2: ',n1, n2
        stop 'program aborted with internal error'
        
endif
end subroutine check_sizereal2D

subroutine check_sizereal3D(message, x, n1, n2, n3)
implicit none
real, intent(in)                        :: x(:,:,:)
integer, intent(in)                     :: n1, n2, n3
character(len=*)                        :: message

if(size(x,1) .ne. n1 .or. size(x,2) .ne. n2 .or. size(x,3) .ne. n3) then
        write(*,*) message
        write(*,*) 'size mismatch'
        write(*,*) 'shape of x: ', shape(x)
        write(*,*) 'n1,n2,n3: ',n1, n2, n3
        stop 'program aborted with internal error'
        
endif
end subroutine check_sizereal3D

subroutine check_sizereal4D(message, x, n1, n2, n3,n4)
implicit none
real, intent(in)                        :: x(:,:,:,:)
integer, intent(in)                     :: n1, n2, n3, n4
character(len=*)                        :: message

if(size(x,1) .ne. n1 .or. size(x,2) .ne. n2 .or. size(x,3) .ne. n3 .or. size(x,4) .ne. n4) then
        write(*,*) message
        write(*,*) 'size mismatch'
        write(*,*) 'shape of x: ', shape(x)
        write(*,*) 'n1,n2,n3, n4: ',n1, n2, n3, n4
        stop 'program aborted with internal error'
        
endif
end subroutine check_sizereal4D



subroutine check_sizeint1D(message, x, n1)
implicit none
integer, intent(in)                        :: x(:)
integer, intent(in)                     :: n1
character(len=*)                        :: message

if(size(x,1) .ne. n1) then
        write(*,*) message
        write(*,*) 'size mismatch'
        write(*,*) 'shape of x: ', shape(x)
        write(*,*) 'n1: ',n1 
        stop 'program aborted with internal error'
        
endif
end subroutine check_sizeint1D

subroutine check_sizeint2D(message, x, n1, n2)
implicit none
integer, intent(in)                        :: x(:,:)
integer, intent(in)                     :: n1, n2
character(len=*)                        :: message

if(size(x,1) .ne. n1 .or. size(x,2) .ne. n2 ) then
        write(*,*) message
        write(*,*) 'size mismatch'
        write(*,*) 'shape of x: ', shape(x)
        write(*,*) 'n1,n2: ',n1, n2
        stop 'program aborted with internal error'
        
endif
end subroutine check_sizeint2D

subroutine check_sizeint3D(message, x, n1, n2, n3)
implicit none
integer, intent(in)                        :: x(:,:,:)
integer, intent(in)                     :: n1, n2, n3
character(len=*)                        :: message

if(size(x,1) .ne. n1 .or. size(x,2) .ne. n2 .or. size(x,3) .ne. n3) then
        write(*,*) message
        write(*,*) 'size mismatch'
        write(*,*) 'shape of x: ', shape(x)
        write(*,*) 'n1,n2,n3: ',n1, n2, n3
        stop 'program aborted with internal error'
        
endif
end subroutine check_sizeint3D


subroutine check_size_equalreal1D(message, x,y,z)
implicit none
real, intent(in)                        :: x(:),y(:)
real, intent(in), optional              :: z(:)
character(len=*)                        :: message

if(any(shape(x) .ne. shape(y))) then
        write(*,*) message
        write(*,*) 'size mismatch'
        write(*,*) 'shape of x:', shape(x)
        write(*,*) 'shape of y:', shape(y)
        stop 'program aborted with internal error'
endif
if(present(z)) then
        if(any(shape(x) .ne. shape(z))) then
                write(*,*) message
                write(*,*) 'size mismatch'
                write(*,*) 'shape of x:', shape(x)
                write(*,*) 'shape of z:', shape(y)
                stop 'program aborted with internal error'
        endif
endif
end subroutine check_size_equalreal1D

subroutine check_size_equalrealint1D(message, x,y,z)
implicit none
real, intent(in)                        :: x(:)
integer, intent(in)                     :: y(:)
real, intent(in), optional              :: z(:)
character(len=*)                        :: message

if(any(shape(x) .ne. shape(y))) then
        write(*,*) message
        write(*,*) 'size mismatch'
        write(*,*) 'shape of x:', shape(x)
        write(*,*) 'shape of y:', shape(y)
        stop 'program aborted with internal error'
endif
if(present(z)) then
        if(any(shape(x) .ne. shape(z))) then
                write(*,*) message
                write(*,*) 'size mismatch'
                write(*,*) 'shape of x:', shape(x)
                write(*,*) 'shape of z:', shape(y)
                stop 'program aborted with internal error'
        endif
endif
end subroutine check_size_equalrealint1D


subroutine check_size_equalreal2D(message, x,y,z)
implicit none
real, intent(in)                        :: x(:,:),y(:,:)
real, intent(in), optional              :: z(:,:)
character(len=*)                        :: message

if(any(shape(x) .ne. shape(y))) then
        write(*,*) message
        write(*,*) 'size mismatch'
        write(*,*) 'shape of x:', shape(x)
        write(*,*) 'shape of y:', shape(y)
        stop 'program aborted with internal error'
endif
if(present(z)) then
        if(any(shape(x) .ne. shape(z))) then
                write(*,*) message
                write(*,*) 'size mismatch'
                write(*,*) 'shape of x:', shape(x)
                write(*,*) 'shape of z:', shape(y)
                stop 'program aborted with internal error'
        endif
endif
end subroutine check_size_equalreal2D

subroutine check_size_equalreal3D(message, x,y,z)
implicit none
real, intent(in)                        :: x(:,:,:),y(:,:,:)
real, intent(in), optional              :: z(:,:,:)
character(len=*)                        :: message

if(any(shape(x) .ne. shape(y))) then
        write(*,*) message
        write(*,*) 'size mismatch'
        write(*,*) 'shape of x:', shape(x)
        write(*,*) 'shape of y:', shape(y)
        stop 'program aborted with internal error'
endif
if(present(z)) then
        if(any(shape(x) .ne. shape(z))) then
                write(*,*) message
                write(*,*) 'size mismatch'
                write(*,*) 'shape of x:', shape(x)
                write(*,*) 'shape of y:', shape(y)
                write(*,*) 'shape of z:', shape(y)
                stop 'program aborted with internal error'
        endif
endif
end subroutine check_size_equalreal3D

! check nans
subroutine is_nanint1D(message, x)
! prints message and abort it x is nan
implicit none
integer, intent(in)                     :: x(:)
character(len=*)                        :: message
if(any(.not.(x .eq. x))) then
        write(*,*) 'NaN values'
        write(*,*) message
endif
end subroutine is_nanint1D
subroutine is_nanint(message, x)
implicit none
integer, intent(in)                     :: x
character(len=*)                        :: message
if((.not.(x .eq. x))) then
        write(*,*) 'NaN values'
        write(*,*) message
endif
end subroutine is_nanint


subroutine is_nanint2D(message, x)
implicit none
integer, intent(in)                     :: x(:,:)
character(len=*)                        :: message
if(any(.not.(x .eq. x))) then
        write(*,*) 'NaN values'
        write(*,*) message
endif
end subroutine is_nanint2D
subroutine is_nanint3D(message, x)
implicit none
integer, intent(in)                     :: x(:,:,:)
character(len=*)                        :: message
if(any(.not.(x .eq. x))) then
        write(*,*) 'NaN values'
        write(*,*) message
endif
end subroutine is_nanint3D
subroutine is_nanreal3D(message, x)
implicit none
real, intent(in)                        :: x(:,:,:)
character(len=*)                        :: message
if(any(.not.(x .eq. x))) then
        write(*,*) 'NaN values'
        write(*,*) message
endif
end subroutine is_nanreal3D
subroutine is_nanreal2D(message, x)
implicit none
real, intent(in)                        :: x(:,:)
character(len=*)                        :: message
if(any(.not.(x .eq. x))) then
        write(*,*) 'NaN values'
        write(*,*) message
endif
end subroutine is_nanreal2D
subroutine is_nanreal1D(message, x)
implicit none
real, intent(in)                        :: x(:)
character(len=*)                        :: message
if(any(.not.(x .eq. x))) then
        write(*,*) 'NaN values'
        write(*,*) message
endif
end subroutine is_nanreal1D
subroutine is_nanreal(message, x)
implicit none
real, intent(in)                        :: x
character(len=*)                        :: message
if((.not.(x .eq. x))) then
        write(*,*) 'NaN values'
        write(*,*) message
endif
end subroutine is_nanreal

subroutine is_nancomplex2D(message, x)
implicit none
complex, intent(in)                     :: x(:,:)
character(len=*)                        :: message
if(any(.not.(x .eq. x))) then
        write(*,*) 'NaN values'
        write(*,*) message
endif
end subroutine is_nancomplex2D


end module error_messages

