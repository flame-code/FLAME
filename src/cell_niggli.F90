!program diver_niggli
!implicit none
!integer:: option,i,imax,nat
!real(8):: pi,latvec_in(3,3),latvec_out(3,3),dist_ang(6),epslen,epsvol,nigmat(6),eps,transmat(3,3),a(3,3),vol,energy,enthalpy,dist_ang_init(6)
!character(40)::filename,file_pre
!character(5):: fn5
!logical:: file_exists
!pi=acos(-1.d0)
!
!epsvol=1.d-6
!imax=1000
!file_pre="poslow"
!!file_pre="poslocm_"
!!file_pre="pos_problem"
!do i=1,imax
!write(fn5,'(i5.5)') i
!  filename=trim(file_pre)//fn5//".ascii"
!  file_exists=.false.
!  INQUIRE(FILE=trim(filename), EXIST=file_exists)
!  if(file_exists) then
!
!    write(888,*) trim(filename)
!    open(unit=2,file=trim(filename))
!    read(2,*) nat,enthalpy,energy
!    latvec_in=0.d0
!    read(2,*) latvec_in(1,1),latvec_in(1,2),latvec_in(2,2) 
!    read(2,*) latvec_in(1,3),latvec_in(2,3),latvec_in(3,3) 
!    close(2)
!    call latvec2dist_ang(dist_ang_init,latvec_in,pi)
!
!    a=latvec_in
!    vol= a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)+&
!         a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)
!!    write(*,'(8(es18.10))') enthalpy,dist_ang_init,vol
!    eps=epsvol*vol**(1.d0/3.d0)
!    call niggli(latvec_in,latvec_out,transmat,eps)
!    a=latvec_out
!    vol= a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)+&
!         a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)
!    call latvec2dist_ang(dist_ang,latvec_out,pi)
!    write(*,'(8(es18.10))') enthalpy,dist_ang,vol 
!    endif
!enddo
!end program


subroutine fixcell_niggli(nat,latvec,xred)
implicit none
!This routine will apply the niggli reduction to the cell and to the reduced coordinates and transform them back into the cell
integer:: nat,iat
real(8):: latvec(3,3),xred(3,nat),latvec_out(3,3),eps,transmat(3,3),imat(3,3),epsvol,a(3,3),vol

epsvol=1.d-6
    a=latvec
    vol= a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)+&
         a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)
    eps=epsvol*vol**(1.d0/3.d0)
    call niggli(latvec,latvec_out,transmat,eps)
    call invertmat(transmat,imat,3)
    do iat=1,nat
       xred(:,iat)=matmul(imat,xred(:,iat))
       xred(1,iat)=modulo(xred(1,iat),1.d0)
       xred(2,iat)=modulo(xred(2,iat),1.d0)
       xred(3,iat)=modulo(xred(3,iat),1.d0)
    enddo
latvec=latvec_out
end subroutine


subroutine niggli(latvec_in,latvec_out,transmat,eps)
!This soubroutine will produce the niggli reduced cell based on the improved algorithm of R. W. Grosse-Kunstleve,* N. K. Sauter and P. D. Adams
!On exit, latvec_out will contain the reduced cell, and transmat will provide the transformation operator in latvec_in
implicit none
real(8):: latvec_in(3,3),latvec_out(3,3),transmat(3,3),eps,tmpmat(3,3)
real(8):: pi,nigmat(6),dist_ang(6),nigmat_check(6)
integer:: iout,l,m,n
!logical:: def_gt_0,feq,flt,fgt,is_niggli_cell,
logical:: debug
logical:: def_gt_0, fgt, feq, flt, is_niggli_cell

!debug=.true.
debug=.false.

pi=acos(-1.d0)
transmat=0.d0
transmat(1,1)=1.d0
transmat(2,2)=1.d0
transmat(3,3)=1.d0
!0 step
!Initiallize the niggli matrix
call init_nigmat(latvec_in,nigmat,l,m,n,eps,pi)
if(debug) write(888,'(a,6(1x,es15.7))') "step 0",nigmat
!1 step
1000 continue
    ! A1
    if (fgt(nigmat(1),nigmat(2),eps) .or. (feq(nigmat(1), nigmat(2), eps) .and. fgt(abs(nigmat(4)), abs(nigmat(5)),eps))) then
      call a1_action(nigmat,tmpmat,eps)
      transmat=matmul(transmat,tmpmat)
if(debug) then
      write(888,'(a,6(1x,es15.7))') "step 1",nigmat
!Check the transformation matrix
            latvec_out=matmul(latvec_in,transmat)
            call init_nigmat(latvec_out,nigmat_check,l,m,n,eps,pi)
            if(.not.all(abs(nigmat_check-nigmat).lt.1.d-7)) then
            write(*,*) "Transformation matrix wrong"
            write(888,'(a,6(1x,es15.7))') "step 1",nigmat_check
            stop
            endif
endif
    endif
    ! A2
    if (fgt(nigmat(2), nigmat(3),eps) .or. (feq(nigmat(2),nigmat(3),eps) .and. fgt(abs(nigmat(5)), abs(nigmat(6)),eps))) then
      call a2_action(nigmat,tmpmat,eps)
      transmat=matmul(transmat,tmpmat)
if(debug) then
      write(888,'(a,6(1x,es15.7))') "step 2",nigmat
!Check the transformation matrix
            latvec_out=matmul(latvec_in,transmat)
            call init_nigmat(latvec_out,nigmat_check,l,m,n,eps,pi)
            if(.not.all(abs(nigmat_check-nigmat).lt.1.d-7)) then
            write(*,*) "Transformation matrix wrong"
            write(888,'(a,6(1x,es15.7))') "step 2",nigmat_check
            stop
            endif
endif
      goto 1000
    endif
    ! A3
    if (def_gt_0(nigmat,eps)) then
      call a3_action(nigmat,tmpmat,eps)
      transmat=matmul(transmat,tmpmat)
if(debug) then
      write(888,'(a,6(1x,es15.7))') "step 3",nigmat
!Check the transformation matrix
            latvec_out=matmul(latvec_in,transmat)
            call init_nigmat(latvec_out,nigmat_check,l,m,n,eps,pi)
            if(.not.all(abs(nigmat_check-nigmat).lt.1.d-7)) then
            write(*,*) "Transformation matrix wrong"
            write(888,'(a,6(1x,es15.7))') "step 3",nigmat_check
            stop
            endif
endif
    else
    ! A4
      call a4_action(nigmat,tmpmat,eps)
      transmat=matmul(transmat,tmpmat)
if(debug) then
      write(888,'(a,6(1x,es15.7))') "step 4",nigmat
!Check the transformation matrix
            latvec_out=matmul(latvec_in,transmat)
            call init_nigmat(latvec_out,nigmat_check,l,m,n,eps,pi)
            if(.not.all(abs(nigmat_check-nigmat).lt.1.d-7)) then
            write(*,*) "Transformation matrix wrong"
            write(888,'(a,6(1x,es15.7))') "step 4",nigmat_check
            stop
            endif
endif
    endif 
    ! A5
    if (fgt(abs(nigmat(4)), nigmat(2),eps) &
        &.or. (feq(nigmat(4), nigmat(2),eps) .and. flt(nigmat(5)+nigmat(5), nigmat(6),eps))&
        &.or. (feq(nigmat(4),-nigmat(2),eps) .and. flt(nigmat(6), 0.d0,eps))) then
      call a5_action(nigmat,tmpmat,eps)
      transmat=matmul(transmat,tmpmat)
if(debug) then
      write(888,'(a,6(1x,es15.7))') "step 5",nigmat
!Check the transformation matrix
            latvec_out=matmul(latvec_in,transmat)
            call init_nigmat(latvec_out,nigmat_check,l,m,n,eps,pi)
            if(.not.all(abs(nigmat_check-nigmat).lt.1.d-7)) then
            write(*,*) "Transformation matrix wrong"
            write(888,'(a,6(1x,es15.7))') "step 5",nigmat_check
            stop
            endif
endif
      goto 1000
    endif
    ! A6
    if (fgt(abs(nigmat(5)), nigmat(1),eps)&
        &.or. (feq(nigmat(5), nigmat(1),eps) .and. flt(nigmat(4)+nigmat(4), nigmat(6),eps))&
        &.or. (feq(nigmat(5),-nigmat(1),eps) .and. flt(nigmat(6), 0.d0,eps))) then
      call a6_action(nigmat,tmpmat,eps)
      transmat=matmul(transmat,tmpmat)
if(debug) then
      write(888,'(a,6(1x,es15.7))') "step 6",nigmat
!Check the transformation matrix
            latvec_out=matmul(latvec_in,transmat)
            call init_nigmat(latvec_out,nigmat_check,l,m,n,eps,pi)
            if(.not.all(abs(nigmat_check-nigmat).lt.1.d-7)) then
            write(*,*) "Transformation matrix wrong"
            write(888,'(a,6(1x,es15.7))') "step 6",nigmat_check
            stop
            endif
endif
      goto 1000
    endif
    ! A7
    if (fgt(abs(nigmat(6)),nigmat(1),eps)&
        &.or. (feq(nigmat(6), nigmat(1),eps) .and. flt(nigmat(4)+nigmat(4),nigmat(5),eps))&
        &.or. (feq(nigmat(6),-nigmat(1),eps) .and. flt(nigmat(5), 0.d0,eps))) then
      call a7_action(nigmat,tmpmat,eps)
      transmat=matmul(transmat,tmpmat)
if(debug) then
      write(888,'(a,6(1x,es15.7))') "step 7",nigmat
!Check the transformation matrix
            latvec_out=matmul(latvec_in,transmat)
            call init_nigmat(latvec_out,nigmat_check,l,m,n,eps,pi)
            if(.not.all(abs(nigmat_check-nigmat).lt.1.d-7)) then
            write(*,*) "Transformation matrix wrong"
            write(888,'(a,6(1x,es15.7))') "step 7",nigmat_check
            stop
endif
            endif
      goto 1000
    endif
    ! A8
    if (flt(nigmat(4)+nigmat(5)+nigmat(6)+nigmat(1)+nigmat(2),0.d0,eps)&
        &.or. (feq(nigmat(4)+nigmat(5)+nigmat(6)+nigmat(1)+nigmat(2),0.d0,eps)) .and.&
        & fgt(nigmat(1)+nigmat(1)+nigmat(5)+nigmat(5)+nigmat(6),0.d0,eps))then
      call a8_action(nigmat,tmpmat,eps)
      transmat=matmul(transmat,tmpmat)
if(debug) then
      write(888,'(a,6(1x,es15.7))') "step 8",nigmat
!Check the transformation matrix
            latvec_out=matmul(latvec_in,transmat)
            call init_nigmat(latvec_out,nigmat_check,l,m,n,eps,pi)
            if(.not.all(abs(nigmat_check-nigmat).lt.1.d-7)) then
            write(*,*) "Transformation matrix wrong"
            write(888,'(a,6(1x,es15.7))') "step 8",nigmat_check
            stop
            endif
endif
      goto 1000
    endif
latvec_out=matmul(latvec_in,transmat)
if(debug) then
!Check the transformation matrix
    if(.not.is_niggli_cell(nigmat,eps)) stop "We did not generate a niggli cell"
    call init_nigmat(latvec_out,nigmat_check,l,m,n,eps,pi)
    if(.not.is_niggli_cell(nigmat_check,eps))then
    write(*,*) "No niggli cell in check!"
    endif
    if(.not.all(abs(nigmat_check-nigmat).lt.1.d-7)) then
    write(*,*) "Transformation matrix wrong"
    nigmat=nigmat_check
    stop
    endif
endif
end subroutine

subroutine  a1_action(nigmat,tmpmat,eps)
implicit none
real(8):: nigmat(6),tmpmat(3,3),eps
call n1_action(nigmat,tmpmat,eps)
end subroutine
subroutine  a2_action(nigmat,tmpmat,eps)
implicit none
real(8):: nigmat(6),tmpmat(3,3),eps
call n2_action(nigmat,tmpmat,eps)
end subroutine
subroutine  a3_action(nigmat,tmpmat,eps)
implicit none
real(8):: nigmat(6),tmpmat(3,3),eps
call n3_true_action(nigmat,tmpmat,eps)
end subroutine
subroutine  a4_action(nigmat,tmpmat,eps)
implicit none
real(8):: nigmat(6),tmpmat(3,3),eps
call n3_false_action(nigmat,tmpmat,eps)
end subroutine
subroutine  a5_action(nigmat,tmpmat,eps)
implicit none
real(8):: nigmat(6),tmpmat(3,3),eps
    if (nigmat(4) > 0.d0) then
      tmpmat(1,:)=(/1.d0,0.d0, 0.d0/)
      tmpmat(2,:)=(/0.d0,1.d0,-1.d0/)
      tmpmat(3,:)=(/0.d0,0.d0, 1.d0/)
      nigmat(3)=nigmat(3)+ nigmat(2)- nigmat(4)
      nigmat(4)=nigmat(4)-(nigmat(2)+nigmat(2))
      nigmat(5)=nigmat(5)-(nigmat(6))
    else
      tmpmat(1,:)=(/1.d0,0.d0, 0.d0/)
      tmpmat(2,:)=(/0.d0,1.d0, 1.d0/)
      tmpmat(3,:)=(/0.d0,0.d0, 1.d0/)
      nigmat(3)=nigmat(3)+ nigmat(2)+ nigmat(4)
      nigmat(4)=nigmat(4)+(nigmat(2)+nigmat(2))
      nigmat(5)=nigmat(5)+(nigmat(6))
    endif
    if(nigmat(3).lt.0.d0) stop "Error in action a5"
end subroutine

subroutine  a6_action(nigmat,tmpmat,eps)
implicit none
real(8):: nigmat(6),tmpmat(3,3),eps
logical:: flt,fgt,feq
    if (nigmat(5) > 0.d0) then
      tmpmat(1,:)=(/1.d0,0.d0,-1.d0/)
      tmpmat(2,:)=(/0.d0,1.d0, 0.d0/)
      tmpmat(3,:)=(/0.d0,0.d0, 1.d0/)
      nigmat(3)=nigmat(3)+ nigmat(1)- nigmat(5)
      nigmat(4)=nigmat(4)-(nigmat(6))
      nigmat(5)=nigmat(5)-(nigmat(1)+nigmat(1))
    else
      tmpmat(1,:)=(/1.d0,0.d0, 1.d0/)
      tmpmat(2,:)=(/0.d0,1.d0, 0.d0/)
      tmpmat(3,:)=(/0.d0,0.d0, 1.d0/)
      nigmat(3)=nigmat(3)+ nigmat(1)+ nigmat(5)
      nigmat(4)=nigmat(4)+(nigmat(6))
      nigmat(5)=nigmat(5)+(nigmat(1)+nigmat(1))
    endif
    if(nigmat(3).lt.0.d0) stop "Error in action a6" 
end subroutine

subroutine  a7_action(nigmat,tmpmat,eps)
implicit none
real(8):: nigmat(6),tmpmat(3,3),eps
    if (nigmat(6) > 0.d0) then
      tmpmat(1,:)=(/1.d0,-1.d0,0.d0/)
      tmpmat(2,:)=(/0.d0,1.d0, 0.d0/)
      tmpmat(3,:)=(/0.d0,0.d0, 1.d0/)
      nigmat(2)=nigmat(2)+ nigmat(1)- nigmat(6)
      nigmat(4)=nigmat(4)-(nigmat(5))
      nigmat(6)=nigmat(6)-(nigmat(1)+nigmat(1))
    else
      tmpmat(1,:)=(/1.d0, 1.d0,0.d0/)
      tmpmat(2,:)=(/0.d0,1.d0, 0.d0/)
      tmpmat(3,:)=(/0.d0,0.d0, 1.d0/)
      nigmat(2)=nigmat(2)+ nigmat(1)+ nigmat(6)
      nigmat(4)=nigmat(4)+(nigmat(5))
      nigmat(6)=nigmat(6)+(nigmat(1)+nigmat(1))
    endif
    if(nigmat(2).lt.0.d0) stop "Error in action a7"
end subroutine


subroutine  a8_action(nigmat,tmpmat,eps)
implicit none
real(8):: nigmat(6),tmpmat(3,3),eps
      tmpmat(1,:)=(/1.d0,0.d0, 1.d0/)
      tmpmat(2,:)=(/0.d0,1.d0, 1.d0/)
      tmpmat(3,:)=(/0.d0,0.d0, 1.d0/)
      nigmat(3)=nigmat(3)+ nigmat(1)+ nigmat(2)+ nigmat(4)+ nigmat(5)+ nigmat(6)
      nigmat(4)=nigmat(4)+ nigmat(2)+ nigmat(2)+ nigmat(6)
      nigmat(5)=nigmat(5)+ nigmat(1)+ nigmat(1)+ nigmat(6)
    if(nigmat(3).lt.0.d0) stop "Error in action a7"
end subroutine

subroutine  n1_action(nigmat,tmpmat,eps)
implicit none
real(8):: nigmat(6),tmpmat(3,3),eps,tmp
      tmpmat(1,:)=(/ 0.d0,-1.d0, 0.d0/)
      tmpmat(2,:)=(/-1.d0, 0.d0, 0.d0/)
      tmpmat(3,:)=(/ 0.d0, 0.d0,-1.d0/)
    tmp=nigmat(1) 
    nigmat(1)=nigmat(2)
    nigmat(2)=tmp
    tmp=nigmat(4)
    nigmat(4)=nigmat(5)
    nigmat(5)=tmp
end subroutine

subroutine  n2_action(nigmat,tmpmat,eps)
implicit none
real(8):: nigmat(6),tmpmat(3,3),eps,tmp
      tmpmat(1,:)=(/-1.d0, 0.d0, 0.d0/)
      tmpmat(2,:)=(/ 0.d0, 0.d0,-1.d0/)
      tmpmat(3,:)=(/ 0.d0,-1.d0, 0.d0/)
    tmp=nigmat(2) 
    nigmat(2)=nigmat(3)
    nigmat(3)=tmp
    tmp=nigmat(5)
    nigmat(5)=nigmat(6)
    nigmat(6)=tmp
end subroutine

subroutine  n3_true_action(nigmat,tmpmat,eps)
implicit none
real(8):: nigmat(6),tmpmat(3,3),eps,tmp
logical:: flt
real(8):: i,j,k
i=1.d0
j=1.d0
k=1.d0
    if (flt(nigmat(4), 0.d0,eps)) i = -1.d0
    if (flt(nigmat(5), 0.d0,eps)) j = -1.d0
    if (flt(nigmat(6), 0.d0,eps)) k = -1.d0
      tmpmat(1,:)=(/ i,   0.d0, 0.d0/)
      tmpmat(2,:)=(/ 0.d0, j  , 0.d0/)
      tmpmat(3,:)=(/ 0.d0,0.d0,    k/)
    nigmat(4)=abs(nigmat(4))
    nigmat(5)=abs(nigmat(5))
    nigmat(6)=abs(nigmat(6))
end subroutine

subroutine  n3_false_action(nigmat,tmpmat,eps)
implicit none
real(8):: nigmat(6),tmpmat(3,3),eps,tmp
logical:: flt,fgt
real(8):: i,j,k,f(3)
integer:: z
    f = 1.d0
    z = -1
    if (fgt(nigmat(4), 0.d0,eps)) then 
        f(1) = -1.d0
    elseif (.not. flt(nigmat(4), 0.d0,eps))then
        z = 1
    endif
    if (fgt(nigmat(5), 0.d0,eps)) then 
        f(2) = -1.d0
    elseif (.not. flt(nigmat(5), 0.d0,eps))then 
        z = 2
    endif
    if (fgt(nigmat(6), 0.d0,eps)) then 
        f(3) = -1.d0
    elseif (.not. flt(nigmat(6), 0.d0,eps))then 
        z = 3
    endif
    if (f(1)*f(2)*f(3) < 0.d0) then
       if(z==-1) stop "Wrong z in n4"
       f(z) = -1.d0
    endif
    tmpmat(1,:)=(/ f(1), 0.d0, 0.d0/)
    tmpmat(2,:)=(/ 0.d0, f(2), 0.d0/)
    tmpmat(3,:)=(/ 0.d0,0.d0,  f(3)/)
    nigmat(4)=-abs(nigmat(4))
    nigmat(5)=-abs(nigmat(5))
    nigmat(6)=-abs(nigmat(6))
end subroutine

function def_gt_0(nigmat,eps)
implicit none
real(8):: nigmat(6),eps
integer:: nzero,npositive,i
logical:: flt,def_gt_0
def_gt_0=.false.
nzero=0
npositive=0
do i=4,6
!if(0.d0.lt.nigmat(i)) then
if(flt(0.d0,nigmat(i),eps)) then
  npositive=npositive+1
  elseif (.not.flt(nigmat(i),0.d0,eps)) then
  nzero=nzero+1
endif
enddo
if(npositive==3.or.(nzero==0.and.npositive==1)) def_gt_0=.true.
end function

function fpositive(nigmat,eps)
implicit none
real(8):: nigmat(6),eps
integer:: nzero,npositive,i
logical:: fpositive,flt
fpositive=.false.
nzero=0
npositive=0
do i=4,6
!if(0.d0.lt.nigmat(i)) then
if(flt(0.d0,nigmat(i),eps)) then
  npositive=npositive+1
  elseif (.not.flt(nigmat(i),0.d0,eps)) then
  nzero=nzero+1
endif
enddo
if(npositive==3.or.(nzero==0.and.npositive==1)) fpositive=.true.
end function

function flt(a,b,eps) 
implicit none
real(8)::a,b,eps
logical:: flt
if(a.lt.b-eps) then
  flt=.true.
else
  flt=.false.
endif
end function

function fgt(a,b,eps) 
implicit none
real(8)::a,b,eps
logical:: fgt,flt
fgt=flt(b,a,eps)
end function

function fle(a,b,eps) 
implicit none
real(8)::a,b,eps
logical:: fle
if(.not.(b.lt.a-eps)) then
  fle=.true.
else
  fle=.false.
endif
end function

function fge(a,b,eps) 
implicit none
real(8)::a,b,eps
logical:: fge
if(.not.(a.lt.b-eps)) then
  fge=.true.
else
  fge=.false.
endif
end function

function feq(a,b,eps) 
implicit none
real(8)::a,b,eps
logical:: feq,flt
feq=(.not. (flt(a,b,eps).or.flt(b, a,eps)))
end function

!subroutine init_nigmat(dist_ang,nigmat,l,m,n,eps,pi)
subroutine init_nigmat(latvec,nigmat,l,m,n,eps,pi)
!Initiallization step (step 0I of the Krivy Gruber algo 
implicit none
real(8):: dist_ang(6),nigmat(6),pi,convang,eps,latvec(3,3)
integer:: l,m,n
convang=pi/180.d0
nigmat(1)=dot_product(latvec(:,1),latvec(:,1))!dist_ang(1)**2
nigmat(2)=dot_product(latvec(:,2),latvec(:,2))!dist_ang(2)**2
nigmat(3)=dot_product(latvec(:,3),latvec(:,3))!dist_ang(3)**2
nigmat(4)=2.d0*dot_product(latvec(:,2),latvec(:,3))!2.d0*dist_ang(2)*dist_ang(3)*cos(convang*dist_ang(4))
nigmat(5)=2.d0*dot_product(latvec(:,3),latvec(:,1))!2.d0*dist_ang(3)*dist_ang(1)*cos(convang*dist_ang(5))
nigmat(6)=2.d0*dot_product(latvec(:,1),latvec(:,2))!2.d0*dist_ang(1)*dist_ang(2)*cos(convang*dist_ang(6))
l=0;m=0;n=0
end subroutine

subroutine init_lmn(nigmat,l,m,n,eps)
implicit none
real(8):: dist_ang(6),nigmat(6),eps
integer:: l,m,n
l=0
m=0
n=0
if(nigmat(4).lt.-eps) l=-1
if(nigmat(4).gt. eps) l= 1
if(nigmat(5).lt.-eps) m=-1
if(nigmat(5).gt. eps) m= 1
if(nigmat(6).lt.-eps) n=-1
if(nigmat(6).gt. eps) n= 1
end subroutine


!subroutine dist_ang2latvec(dist_ang,latvec,pi)
!!This subroutine will generate the lattice vector representation of the cell
!!from the length/angle representation
!implicit none
!real(8):: dist_ang(6),latvec(3,3),convang,pi
!convang=pi/180.d0
!
!latvec(1,1)=dist_ang(1)
!latvec(2,1)=0.d0
!latvec(3,1)=0.d0
!latvec(1,2)=dist_ang(2)*cos(convang*dist_ang(6))
!latvec(2,2)=dist_ang(2)*sin(convang*dist_ang(6))
!latvec(3,2)=0.d0
!latvec(1,3)=dist_ang(3)*cos(convang*dist_ang(5))
!latvec(2,3)=(dist_ang(2)*dist_ang(3)*cos(convang*dist_ang(4))-latvec(1,2)*latvec(1,3))/latvec(2,2)
!latvec(3,3)=(dist_ang(3)**2-latvec(1,3)**2-latvec(2,3)**2)
!if(latvec(3,3).lt.0.d0) stop "Bad angles in dist_ang2latvec"
!latvec(3,3)=sqrt(latvec(3,3))
!end subroutine

subroutine latvec2dist_ang(dist_ang,latvec,pi)
!This subroutine will generate the distance and angles of a cell starting from latvec
implicit none
real(8):: dist_ang(6),latvec(3,3),pi,convang
convang=180.d0/pi
dist_ang(1)=sqrt(dot_product(latvec(:,1),latvec(:,1)))
dist_ang(2)=sqrt(dot_product(latvec(:,2),latvec(:,2)))
dist_ang(3)=sqrt(dot_product(latvec(:,3),latvec(:,3)))
dist_ang(4)=convang*acos(dot_product(latvec(:,2),latvec(:,3))/(dist_ang(2)*dist_ang(3)))
dist_ang(5)=convang*acos(dot_product(latvec(:,3),latvec(:,1))/(dist_ang(3)*dist_ang(1)))
dist_ang(6)=convang*acos(dot_product(latvec(:,1),latvec(:,2))/(dist_ang(1)*dist_ang(2)))
end subroutine

function is_niggli_cell(nigmat,eps) 
implicit none
real(8):: nigmat(6),eps
logical::is_niggli_cell,is_buerger_cell,feq,fgt
is_niggli_cell=.true.
    if (.not.is_buerger_cell(nigmat,eps)) then
        is_niggli_cell=.false.
        return
    endif
    if (feq(nigmat(4),nigmat(2),eps)) then
      if (fgt(nigmat(6),nigmat(5)+nigmat(5),eps)) then 
        is_niggli_cell=.false.
        return
      endif
    endif
    if (feq(nigmat(5),nigmat(1),eps)) then
      if (fgt(nigmat(6),nigmat(4)+nigmat(4),eps)) then 
        is_niggli_cell=.false.
        return
      endif
    endif
    if (feq(nigmat(6),nigmat(1),eps)) then
      if (fgt(nigmat(5),nigmat(4)+nigmat(4),eps)) then 
        is_niggli_cell=.false.
        return
      endif
    endif
    if (feq(nigmat(4),-nigmat(2),eps)) then
      if (.not.feq(nigmat(6),0.d0,eps)) then 
        is_niggli_cell=.false.
        return
      endif
    endif
    if (feq(nigmat(5),-nigmat(1),eps)) then
      if (.not.feq(nigmat(6),0.d0,eps)) then 
        is_niggli_cell=.false.
        return
      endif
    endif
    if (feq(nigmat(6),-nigmat(1),eps)) then
      if (.not.feq(nigmat(5),0.d0,eps)) then 
        is_niggli_cell=.false.
        return
      endif
    endif
    if (feq(nigmat(4)+nigmat(5)+nigmat(6)+nigmat(1)+nigmat(2),0.d0,eps)) then
      if (fgt(nigmat(1)+nigmat(1)+nigmat(5)+nigmat(5)+nigmat(6),0.d0,eps)) then 
        is_niggli_cell=.false.
        return
      endif
    endif
end function

function  meets_primary_conditions(nigmat,eps) 
implicit none
real(8):: nigmat(6),eps
logical::meets_primary_conditions,fgt
meets_primary_conditions=.true.
    if (fgt(nigmat(1), nigmat(2),eps)) then
        meets_primary_conditions=.false.
        return
    endif
    if (fgt(nigmat(2), nigmat(3),eps)) then
        meets_primary_conditions=.false.
        return
    endif
    if (fgt(abs(nigmat(4)), nigmat(2),eps)) then
        meets_primary_conditions=.false.
        return
    endif
    if (fgt(abs(nigmat(5)), nigmat(1),eps)) then
        meets_primary_conditions=.false.
        return
    endif
    if (fgt(abs(nigmat(6)), nigmat(1),eps)) then
        meets_primary_conditions=.false.
        return
    endif
end function

function  meets_main_conditions(nigmat,eps) 
implicit none
integer:: typer
real(8):: nigmat(6),eps
logical::meets_main_conditions,meets_primary_conditions,flt
meets_main_conditions=.true.
    if (.not. meets_primary_conditions(nigmat,eps)) then
      meets_main_conditions=.false.
      return 
    endif
    
    if (typer(nigmat,eps) == 0) then 
      meets_main_conditions=.false.
      return
    endif
    if (typer(nigmat,eps) == 2) then
      if (flt(nigmat(4)+nigmat(5)+nigmat(6)+nigmat(1)+nigmat(2), 0.d0,eps)) then
          meets_main_conditions=.false.
          return
      endif
    endif
end function

function is_buerger_cell(nigmat,eps) 
implicit none
real(8):: nigmat(6),eps
logical::is_buerger_cell,meets_main_conditions,feq,fgt

    if (.not. meets_main_conditions(nigmat,eps)) then
        is_buerger_cell=.false.
        return
    endif
    if (feq(nigmat(1),nigmat(2),eps)) then
      if (fgt(abs(nigmat(4)),abs(nigmat(5)),eps)) then
          is_buerger_cell=.false.
          return
      endif
    endif
    if (feq(nigmat(2), nigmat(3),eps)) then
      if (fgt(abs(nigmat(5)), abs(nigmat(6)),eps))then 
          is_buerger_cell=.false.
          return
      endif
   endif
  is_buerger_cell=.true. 

end function

function typer(nigmat,eps) 
implicit none
real(8):: nigmat(6),eps
integer:: nzero,npositive,i,typer
logical:: flt
typer=0
nzero=0
npositive=0
do i=4,6
!if(0.d0.lt.nigmat(i)) then
if(flt(0.d0,nigmat(i),eps)) then
  npositive=npositive+1
  elseif (.not.flt(nigmat(i),0.d0,eps)) then
  nzero=nzero+1
endif
enddo
    typer=0
    if (npositive == 3) typer=1
    if (npositive == 0) typer=2
end function

