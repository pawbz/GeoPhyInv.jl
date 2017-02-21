module string_routines
use precision_mod! USEINDEPFILE
use iso_c_binding
implicit none

public
contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! subroutines for searching a keyword, making up a filename and others related with characters.
! go_string
! 
!  filename(fname,str1,num,str2)
!               #output a string fname which is consisted of
!               #str1, num and str2
!  get_num(x,nn,str,line,nl)
!               #output nn numbers to array x
!               #from a line containning the key word str and '='
!  get_str(str_new,len_new,str,line,nl) 
!               #output a string str_new and its length len_new
!               #from a line containning the key word str and '='
!  strlen(str,lenn)
!               #output the length of string str without any blank in the
!               #middle of it 
!  num2str(string,num,len,len_new,itb,ch)
!               #output a string string which is converted from the number 
!               #num with a maximum length len 


function ctofstr(s) result(str)
implicit none
character(kind=c_char,len=1), intent(in) :: s(*)
character(len=:), allocatable            :: str
integer                                  :: i, nchars
i = 1
do
        if (s(i) == c_null_char) exit
        i = i + 1
end do
nchars = i - 1  ! Exclude null character from Fortran string
allocate(character(len=nchars) :: str)
str = transfer(s(1:nchars), str)
end function ctofstr

! 
function n2str(num, formatspec)
! this fucntion converts real to string
implicit none
real, intent(in)                        :: num
character(len=max_str_len)              :: n2str
character(len=*), intent(in), optional  :: formatspec

if(.not.(present(formatspec))) then
        write(n2str,*)num
        n2str=trim(adjustl(n2str))
else
        write(n2str,formatspec)num
        n2str=trim(adjustl(n2str))
endif
end function n2str

function i2str(num)
! this fucntion converts integer to string
implicit none
integer, intent(in)                     :: num
character(len=max_str_len)              :: i2str
write(i2str,*)num
i2str=trim(adjustl(i2str))
end function i2str


!pure subroutine itoa(number, str, str_len)
!!     ITOA converts an interger 'number' to a string 'str'.
!implicit none
!integer number, str_len
!character(*) str
!integer i, indx, n
!if (number .eq. 0) then
!str(1:1) = '0'
!str_len = 1
!return
!endif
!str_len = int(log10(dble(number))+1.0d0)
!if (str_len.gt.0) then
!do i=1,str_len
!  indx = str_len-i+1
!  n = (number/10)*10
!  str(indx:indx) = char(48+(number-n))
!  number = n/10
!enddo
!endif
!write(*,*) str(1:str_len)
!end subroutine itoa
!
pure subroutine get_real_new(x,nn,str,line,is_suc)
!this subroutine searches for 'str' in line(1:nl) and reads nn numbers after '=' and 
! assigns the numbers to x 

implicit none
integer                                 :: lenn, i, ii, j, ind, indeq, k, nl, ioerr
integer, intent(in)                     :: nn
integer, intent(out)                    :: is_suc
real, intent(out)                       :: x(:)
character(len=*), intent(in)            :: line(:)
character(len=len(line))                :: line_tmp
character(len=max_str_len)              :: streq
character(len=*), intent(in)            :: str

nl = size(line)
lenn=len(str)
streq = ''
is_suc = 0
! loop over the dimension of line
do i=1,nl
        line_tmp=line(i)
        ! deleting everything after hash
        ii=index(line_tmp,'#')
        if(ii.ne.0)then
                do j=ii,len(line)
                        line_tmp(j:j)=' '
                enddo
        endif
        ind=index(line_tmp,str(1:lenn))
        if(ind .ne. 0)then
                indeq=index(line_tmp,'=')
                streq = trim(adjustl(streq))//trim(adjustl(line_tmp(indeq+1:len(line))))
                is_suc=1
        endif
enddo
! read
read(streq,*,iostat=ioerr) (x(k),k=1,nn)
if(ioerr .ne. 0) then
        is_suc = 0
        return
endif
end subroutine get_real_new

pure subroutine get_real_new_separator(separator,x,nn,str,line,is_suc)
! the string line is converted to a vector line_vec using parse
! this subroutine searches for 'str' in line_vec(1:) and reads 1 numbers after '=' and 
! assigns the numbers to x 
implicit none
real,intent(out)                :: x(:)
integer, intent(in)             :: nn
character(len=*), intent(in)    :: line, separator
character(len=max_str_len)      :: line_vec(max_str_dim)
character(len=*), intent(in)    :: str
integer, intent(out)            :: is_suc
integer                         :: nargs
call parse(line,trim(adjustl(separator)),line_vec,nargs)
call get_real_new(x,nn,trim(str),line_vec,is_suc)
end subroutine get_real_new_separator

pure subroutine get_str_new(str_new,str,line,len_new)
!this subroutine searches for 'str' in every 'line' and assigns 'str_new'
!the string after '=' sign 
! nl is the total number of lines in line(1:nl)
! if nothing found then len_new.eq.0 or else it will some value .eq. len(str_new)

implicit none
integer                                   :: lenn, nl, i, ii, j, is_found, indeq, ioerr
integer, intent(out)                      :: len_new
character(len=*), intent(in)              :: line(:)
character(len=len(line))                  :: line_tmp
character(len=max_str_len)                :: streq
character(len=*), intent(out)             :: str_new
character(len=max_str_len)                :: str_tmp
character(len=*), intent(in)              :: str

nl = size(line)
str_new =' '

lenn=len_trim(str)
len_new=0
streq = ''
! loop over line dimension
do i=1,nl
        line_tmp=(line(i))
        if(len_trim(line_tmp) .ne. 0) then
                ! deleting everything after hash
                ii=index(line_tmp,'#')
                if(ii.ne.0)then
                        do j=ii,max_str_len
                        line_tmp(j:j)=' '
                        enddo
                endif
                is_found=index(line_tmp,str(1:lenn))
                if(is_found.ne.0) then
                        indeq = index(line_tmp,'=')
                        streq = trim(adjustl(streq))//trim(adjustl(line_tmp(indeq+1:max_str_len)))
                endif
        endif
enddo
read(streq,'(a)',iostat=ioerr) str_tmp
if(ioerr .ne. 0) then
        len_new=0
        str_new = ''
else
        len_new = len_trim(str_tmp)
        str_new = (str_tmp)
endif
end  subroutine get_str_new

pure subroutine get_str_new_separator(separator,str_new,str,line,is_suc)
! the string line is converted to a vector line_vec using parse
! a separator is required to do so
! this subroutine searches for 'str' in line_vec(1:) and reads 1 numbers after '=' and 
! assigns the numbers to x 
implicit none
character(len=*), intent(in)              :: line, separator
character(len=max_str_len)                :: line_vec(max_str_dim)
character(len=*), intent(out)             :: str_new
character(len=*), intent(in)              :: str
integer, intent(out)                      :: is_suc
integer                                   :: nargs
call parse(line,trim(adjustl(separator)),line_vec,nargs)
call get_str_new(str_new,str,line_vec,is_suc)
end subroutine get_str_new_separator


pure subroutine removebksl(str)
! Removes backslash (\) characters. Double backslashes (\\) are replaced
! by a single backslash.
implicit none
character(len=*), intent(inout)         :: str
character(len=1):: ch
character(len=len_trim(str))::outstr
integer                         :: lenstr, ibsl, k, i

str=adjustl(str)
lenstr=len_trim(str)
outstr=' '
k=0
ibsl=0                        ! backslash initially inactive
do i=1,lenstr
	ch=str(i:i)
	if(ibsl == 1) then          ! backslash active
		k=k+1
		outstr(k:k)=ch
		ibsl=0
		cycle
	end if
	if(ch == '\\') then          ! backslash with backslash inactive
		ibsl=1
		cycle
	end if
	k=k+1
	outstr(k:k)=ch              ! non-backslash with backslash inactive
end do
str=adjustl(outstr)
end subroutine removebksl

pure subroutine replacechar(str, char1, char2)
! replaces all occurences char1 with char2 in the string
! Removes backslash (\) characters. Double backslashes (\\) are replaced
! by a single backslash.
implicit none
character(len=*), intent(inout)		:: str
character(len=1), intent(in)		:: char1, char2
character(len=1)			:: ch
character(len=len(str))			:: strin
integer                         	:: lenstr, i

strin = str;
lenstr=len(strin);
str = ' ';
do i=1, lenstr
	ch=strin(i:i)
	if(ch == char1) then
		str(i:i) = char2
	else
		str(i:i) = ch
	endif
enddo
end subroutine replacechar



!pure subroutine strlen(fname,lenn)
!
!implicit none
!character*(*) fname
!integer nll, lenn, i
!
!nll=len(fname)
!lenn=0
!do i=1,nll
!if(fname(i:i).ne.' ')then
!  lenn=lenn+1
!else if(lenn.ne.0)then
!  return
!endif
!enddo
!end subroutine strlen

!pure subroutine num2str(string,num,lenstr,len_new,itb,ch)
!
!implicit none
!character ch
!character(len=max_str_len) string
!integer len_new, nnn, num, itb, iii, ii, iic, i, lenstr, i0
!
!if(lenstr.eq.0.or.num.lt.0)then
!len_new=0
!return
!endif
!if(num.eq.0)then
!nnn=1
!else
!nnn=int((log10(dble(num)))+0.000000001d0+1.0d0)
!endif
!iii=num
!if(lenstr.eq.999)then
!len_new=nnn
!else
!len_new=lenstr
!endif
!i0=ichar('0')
!do i=1,len_new
!string(i:i)=ch
!enddo
!do i=1,nnn
!ii=iii/10
!iic=iii-ii*10
!iii=ii
!if(itb.eq.1)then
!  string((nnn-i+1):(nnn-i+1))=char(i0+iic)
!else
!  if(len_new.ge.i)then
!    string((len_new-i+1):(len_new-i+1))=char(i0+iic)
!  endif
!endif
!enddo
!
!end subroutine num2str

pure subroutine split(str,delims,before,sep)

! Routine finds the first instance of a character from 'delims' in the
! the string 'str'. The characters before the found delimiter are
! output in 'before'. The characters after the found delimiter are
! output in 'str'. The optional output character 'sep' contains the 
! found delimiter. A delimiter in 'str' is treated like an ordinary 
! character if it is preceded by a backslash (\). If the backslash 
! character is desired in 'str', then precede it with another backslash.

implicit none
character(len=*), intent(out)            :: str,before
character(len=*), intent(in)             :: delims
character,optional, intent(out)          :: sep
logical :: pres
character :: ch,cha
integer                 :: ibsl, ipos, k, i, lenstr, iposa

pres=present(sep)
str=adjustl(str)
call compact(str)
lenstr=len_trim(str)
if(lenstr == 0) return        ! string str is empty
k=0
ibsl=0                        ! backslash initially inactive
before=' '
do i=1,lenstr
   ch=str(i:i)
   if(ibsl == 1) then          ! backslash active
      k=k+1
      before(k:k)=ch
      ibsl=0
      cycle
   end if
   if(ch == '\\') then          ! backslash with backslash inactive
      k=k+1
      before(k:k)=ch
      ibsl=1
      cycle
   end if

   ipos=index(delims,ch)         
   if(ipos == 0) then          ! character is not a delimiter
      k=k+1
      before(k:k)=ch
      cycle
   end if
   if(ch /= ' ') then          ! character is a delimiter that is not a space
      str=str(i+1:)
      if(pres) sep=ch
      exit
   end if
   cha=str(i+1:i+1)            ! character is a space delimiter
   iposa=index(delims,cha)

   if(iposa > 0) then          ! next character is a delimiter
      str=str(i+2:)
      if(pres) sep=cha
      exit
   else
      str=str(i+1:)
      if(pres) sep=ch
      exit
   end if
end do

if(i >= lenstr) str=''
str=adjustl(str)              ! remove initial spaces
return
end subroutine split

pure subroutine compact(str)
! Converts multiple spaces and tabs to single spaces; deletes control characters;
! removes initial spaces.
implicit none
character(len=*), intent(inout) :: str
character(len=1):: ch
character(len=len_trim(str)):: outstr
integer                                 :: lenstr, k, isp, i, ich 
str=adjustl(str)
lenstr=len_trim(str)
outstr=' '
isp=0
k=0
do i=1,lenstr
  ch=str(i:i)
  ich=iachar(ch)
  select case(ich)
    case(9,32)     ! space or tab character
      if(isp==0) then
        k=k+1
        outstr(k:k)=' '
      end if
      isp=1
    case(33:)      ! not a space, quote, or control character
      k=k+1
      outstr(k:k)=ch
      isp=0
  end select
end do
str=adjustl(outstr)
end subroutine compact

pure subroutine parse(str,delims,args,nargs)
! Parses the string 'str' into arguments args(1), ..., args(nargs) based on
! the delimiters contained in the string 'delims'. Preceding a delimiter in
! 'str' by a backslash (\) makes this particular instance not a delimiter.
! The integer output variable nargs contains the number of arguments found.
implicit none
character(len=*), intent(in)                        :: str
character(len=*), intent(in)                        :: delims
character(len=len_trim(str))                        :: strsav
character(len=*),dimension(:), intent(out)          :: args
integer, intent(out)                                :: nargs
integer                                             :: lenstr,   na, i, k

strsav=str
call compact(strsav)
na=size(args)
do i=1,na
  args(i)=' '
end do  
nargs=0
lenstr=len_trim(strsav)
if(lenstr==0) return
k=0
do
   if(len_trim(strsav) == 0) exit
   nargs=nargs+1
   call split(strsav,delims,args(nargs))
   call removebksl(args(nargs))
end do   
end subroutine parse

end module string_routines

