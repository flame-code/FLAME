subroutine parsearray_real(substring,nsub,string,nstr,array,narr,found)
!This routine will fill the array with narr values if string contains the keyword substring
use String_Utility
implicit none
integer:: narr,nstr,nsub,j,k,k1,k2,nstr1,nstr2
character(nstr):: string
character(nsub):: substring,substring_int
real(8):: array(narr)
logical:: found
found=.false.
string=adjustl(string)
!Find Comment line
nstr1 = len_trim(string)
nstr2 = SCAN(string(1:nstr1), "!#")
if(nstr2.ne.0) nstr1=nstr2-1
if(nstr1.lt.nsub) return



!First make substring to upcase
!write(*,*) substring
 substring_int=StrLowCase ( substring )
   k1 = index(string(1:nsub),trim(substring_int))
 substring_int=StrUpCase ( substring )
   k2 = index(string(1:nsub),trim(substring_int))
!write(*,*) "RUNNING PARSER REAL ARRAY"
    if(k1.ne.0.or.k2.ne.0) then
      found=.true.
      k=max(k1,k2) + LEN (substring)
      nstr2 = SCAN(string(1:nstr1), "=")
      if(nstr2.ne.0) k=nstr2+1
!The string contains narr values
      read(string(k:nstr1),*,err=10,end=20) (array(j),j=1,narr) 
!If successful, then return
      return
      
      10 continue
      write(*,'(a,a)') "ERROR occured in reading variable ",trim(substring) 
      stop
      
      20 continue
      write(*,'(a,i5,a,a)') "Expecting ", narr, " real numbers but not enough data provided for the variable ",trim(substring) 
      stop
    endif
end subroutine

subroutine parsearray_int(substring,nsub,string,nstr,array,narr,found)
!This routine will fill the array with narr values if string contains the keyword substring
use String_Utility
implicit none
integer:: narr,nstr,nsub,j,k,k1,k2,nstr1,nstr2
character(nstr):: string
character(nsub):: substring,substring_int
integer:: array(narr)
logical:: found
found=.false.
string=adjustl(string)
!Find Comment line
nstr1 = len_trim(string)
nstr2 = SCAN(string(1:nstr1), "!#")
if(nstr2.ne.0) nstr1=nstr2-1
if(nstr1.lt.nsub) return



!First make substring to upcase
!write(*,*) substring
 substring_int=StrLowCase ( substring )
   k1 = index(string(1:nsub),trim(substring_int))
 substring_int=StrUpCase ( substring )
   k2 = index(string(1:nsub),trim(substring_int))
!write(*,*) "RUNNING PARSER INTEGER ARRAY"
    if(k1.ne.0.or.k2.ne.0) then
      found=.true.
      k=max(k1,k2) + LEN (substring)
      nstr2 = SCAN(string(1:nstr1), "=")
      if(nstr2.ne.0) k=nstr2+1
!The string contains narr values
      read(string(k:nstr1),*,err=10,end=20) (array(j),j=1,narr) 
!If successful, then return
      return
      
      10 continue
      write(*,'(a,a)') "ERROR occured in reading variable ",trim(substring) 
      stop
      
      20 continue
      write(*,'(a,i5,a,a)') "Expecting ", narr, " real numbers but not enough data provided for the variable ",trim(substring) 
      stop
    endif
end subroutine

subroutine parsescalar_real(substring,nsub,string,nstr,realval,found)
!This routine will fill the array with narr values if string contains the keyword substring
use String_Utility
implicit none
integer:: narr,nstr,nsub,j,k,k1,k2,nstr1,nstr2
character(nstr):: string
character(nsub):: substring,substring_int
real(8):: realval
logical:: found
found=.false.
string=adjustl(string)
!Find Comment line
nstr1 = len_trim(string)
nstr2 = SCAN(string(1:nstr1), "!#")
if(nstr2.ne.0) nstr1=nstr2-1
if(nstr1.lt.nsub) return



!First make substring to upcase
!write(*,*) substring
 substring_int=StrLowCase ( substring )
   k1 = index(string(1:nsub),trim(substring_int))
 substring_int=StrUpCase ( substring )
   k2 = index(string(1:nsub),trim(substring_int))
!write(*,*) "RUNNING PARSER REAL SCALAR"
    if(k1.ne.0.or.k2.ne.0) then
      found=.true.
      k=max(k1,k2) + LEN (substring)
      nstr2 = SCAN(string(1:nstr1), "=")
      if(nstr2.ne.0) k=nstr2+1
!The string contains narr values
      read(string(k:nstr1),*,err=10,end=20) realval
!If successful, then return
      return
      
      10 continue
      write(*,'(a,a)') "ERROR occured in reading variable ",trim(substring) 
      stop
      
      20 continue
      write(*,'(a,i5,a,a)') "Expecting ", 1 , " real numbers but not enough data provided for the variable ",trim(substring) 
      stop
    endif
end subroutine

subroutine parsescalar_int(substring,nsub,string,nstr,intval,found)
!This routine will fill the array with narr values if string contains the keyword substring
use String_Utility
implicit none
integer:: narr,nstr,nsub,j,k,k1,k2,nstr1,nstr2
character(nstr):: string
character(nsub):: substring,substring_int
integer:: intval
logical:: found
found=.false.
string=adjustl(string)
!Find Comment line
nstr1 = len_trim(string)
nstr2 = SCAN(string(1:nstr1), "!#")
if(nstr2.ne.0) nstr1=nstr2-1
if(nstr1.lt.nsub) return



!First make substring to upcase
!write(*,*) substring
 substring_int=StrLowCase ( substring )
   k1 = index(string(1:nsub),trim(substring_int))
 substring_int=StrUpCase ( substring )
   k2 = index(string(1:nsub),trim(substring_int))
!write(*,*) "RUNNING PARSER REAL SCALAR"
    if(k1.ne.0.or.k2.ne.0) then
      found=.true.
      k=max(k1,k2) + LEN (substring)
      nstr2 = SCAN(string(1:nstr1), "=")
      if(nstr2.ne.0) k=nstr2+1
!The string contains narr values
      read(string(k:nstr1),*,err=10,end=20) intval
!If successful, then return
      return
      
      10 continue
      write(*,'(a,a)') "ERROR occured in reading variable ",trim(substring) 
      stop
      
      20 continue
      write(*,'(a,i5,a,a)') "Expecting ", 1 , " integer number but not enough data provided for the variable ",trim(substring) 
      stop
    endif
end subroutine

subroutine parse_logical(substring,nsub,string,nstr,logval,found)
!This routine will fill the array with narr values if string contains the keyword substring
use String_Utility
implicit none
integer:: narr,nstr,nsub,j,k,k1,k2,nstr1,nstr2
character(nstr):: string
character(nsub):: substring,substring_int
logical:: logval
logical:: found
found=.false.
string=adjustl(string)
!Find Comment line
nstr1 = len_trim(string)
nstr2 = SCAN(string(1:nstr1), "!#")
if(nstr2.ne.0) nstr1=nstr2-1
if(nstr1.lt.nsub) return


!First make substring to upcase
!write(*,*) substring
 substring_int=StrLowCase ( substring )
   k1 = index(string(1:nsub),trim(substring_int))
 substring_int=StrUpCase ( substring )
   k2 = index(string(1:nsub),trim(substring_int))
!write(*,*) "RUNNING PARSER REAL SCALAR"
    if(k1.ne.0.or.k2.ne.0) then
      found=.true.
      k=max(k1,k2) + LEN (substring)
      nstr2 = SCAN(string(1:nstr1), "=")
      if(nstr2.ne.0) k=nstr2+1
!The string contains narr values
      read(string(k:nstr1),*,err=10,end=20) logval
!If successful, then return
      return
      
      10 continue
      write(*,'(a,a)') "ERROR occured in reading variable ",trim(substring) 
      stop
      
      20 continue
      write(*,'(a,i5,a,a)') "Expecting ", 1 , " integer number but not enough data provided for the variable ",trim(substring) 
      stop
    endif
end subroutine

subroutine parsescalar_string(substring,nsub,string,nstr,stringval,nstringval,found)
!This routine will fill the array with narr values if string contains the keyword substring
use String_Utility
implicit none
integer:: narr,nstr,nsub,j,k,k1,k2,nstr1,nstr2,nstringval
character(nstr):: string
character(nsub):: substring,substring_int
character(nstringval):: stringval1,stringval
logical:: found
found=.false.
string=adjustl(string)
!Find Comment line
nstr1 = len_trim(string)
nstr2 = SCAN(string(1:nstr1), "!#")
if(nstr2.ne.0) nstr1=nstr2-1
if(nstr1.lt.nsub) return



!First make substring to upcase
!write(*,*) substring
 substring_int=StrLowCase ( substring )
   k1 = index(string(1:nsub),trim(substring_int))
 substring_int=StrUpCase ( substring )
   k2 = index(string(1:nsub),trim(substring_int))
!write(*,*) "RUNNING PARSER STRING"
    if(k1.ne.0.or.k2.ne.0) then
      found=.true.
      k=max(k1,k2) + LEN (substring)
      nstr2 = SCAN(string(1:nstr1), "=")
      if(nstr2.ne.0) k=nstr2+1
!The string contains narr values
      read(string(k:nstr1),*,err=10,end=20) stringval1
      stringval = stringval1
!      stringval =StrUpCase(stringval1)
!If successful, then return
      return
      
      10 continue
      write(*,'(a,a)') "ERROR occured in reading variable ",trim(substring) 
      stop
      
      20 continue
      write(*,'(a,i5,a,a)') "Expecting ", 1 , " string value but not enough data provided for the variable ",trim(substring) 
      stop
    endif
end subroutine


subroutine parsearray_string(substring,nsub,string,nstr,string_array,nstringval,narr,found)
!This routine will fill the array with narr values if string contains the keyword substring
use String_Utility
implicit none
integer:: narr,nstr,nsub,j,k,k1,k2,nstr1,nstr2,nstringval
character(nstr):: string
character(nsub):: substring,substring_int
character(nstringval):: stringval1(narr),string_array(narr)
logical:: found
found=.false.
string=adjustl(string)
!Find Comment line
nstr1 = len_trim(string)
nstr2 = SCAN(string(1:nstr1), "!#")
if(nstr2.ne.0) nstr1=nstr2-1
if(nstr1.lt.nsub) return



!First make substring to upcase
!write(*,*) substring
 substring_int=StrLowCase ( substring )
   k1 = index(string(1:nsub),trim(substring_int))
 substring_int=StrUpCase ( substring )
   k2 = index(string(1:nsub),trim(substring_int))
!write(*,*) "RUNNING PARSER STRING"
    if(k1.ne.0.or.k2.ne.0) then
      found=.true.
      k=max(k1,k2) + LEN (substring)
      nstr2 = SCAN(string(1:nstr1), "=")
      if(nstr2.ne.0) k=nstr2+1
!The string contains narr values
      read(string(k:nstr1),*,err=10,end=20) (stringval1(j),j=1,narr)
      do j=1,narr
         string_array(j) = StrUpCase(stringval1(j))
      enddo
!If successful, then return
      return
      
      10 continue
      write(*,'(a,a)') "ERROR occured in reading variable ",trim(substring) 
      stop
      
      20 continue
      write(*,'(a,i5,a,a)') "Expecting ", narr , " string values but not enough data provided for the variable ",trim(substring) 
      stop
    endif
end subroutine

subroutine exist_string(substring,nsub,string,nstr,found)
!This routine will fill the array with narr values if string contains the keyword substring
use String_Utility
implicit none
integer:: nstr,nsub,j,k,k1,k2,nstr1,nstr2,nstringval
character(nstr):: string
character(nsub):: substring,substring_int
logical:: found
found=.false.
string=adjustl(string)
!Find Comment line
nstr1 = len_trim(string)
nstr2 = SCAN(string(1:nstr1), "!#")
if(nstr2.ne.0) nstr1=nstr2-1
if(nstr1.lt.nsub) return



!First make substring to upcase
!write(*,*) substring
 substring_int=StrLowCase ( substring )
   k1 = index(string(1:nsub),trim(substring_int))
 substring_int=StrUpCase ( substring )
   k2 = index(string(1:nsub),trim(substring_int))
!write(*,*) "RUNNING PARSER STRING"
    if(k1.ne.0.or.k2.ne.0) then
      found=.true.
    endif
end subroutine
