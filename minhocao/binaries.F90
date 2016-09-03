program binaries
!This program computes all possible binary compositions up to a maximum of NMAX
implicit none
integer:: NMAX,i,j,k,l,n
real(8):: con_comp(3,10000),con
real(8),allocatable:: printarray(:,:)
character(30):: formatting

con_comp=1.d10
write(*,*) "Maximum number of atoms per cell"
read(*,*) NMAX
allocate(printarray(3,NMAX))
n=0
do i=0,NMAX
   do j=NMAX-i,0,-1
     if(.not.(i==0.and.j==0)) then
     k=0
     con=real(i,8)/(real(i,8)+real(j,8))
     if (n.ge.1) then
       call hunt(con_comp(1,1:n),n,con,k)
       do l=n,k+1,-1
         con_comp(:,l+1)=con_comp(:,l)
       enddo
     endif
     n=n+1
     con_comp(1,k+1)=con
     con_comp(2,k+1)=real(i,8)
     con_comp(3,k+1)=real(j,8)
     endif
enddo
enddo

k=0
l=1
printarray(:,1)=con_comp(:,1)
do i=1,n
  if(abs(con_comp(1,i)-con_comp(1,i+1)).gt.1.d-15) then
     write(formatting,'(a,i5,a)') "(e12.4,",l,"(5x,2(i3)))"
     write(*,formatting) printarray(1,1),int(printarray(2:3,1:l))
     l=1
     printarray(:,l)=con_comp(:,i+1)
   else
     l=l+1
     printarray(:,l)=con_comp(:,i+1)
   endif
enddo

end program


!************************************************************************************

subroutine hunt(xx,n,x,jlo)
  implicit none
  !C x is in interval [xx(jlo),xx(jlow+1)[ ; xx(0)=-Infinity ; xx(n+1) = Infinity
  !Arguments
  integer :: jlo,n
  real(kind=8) :: x,xx(n)
  !Local variables
  integer :: inc,jhi,jm
  logical :: ascnd
  if (n.le.0) stop 'hunt'
  if (n == 1) then
     if (x.ge.xx(1)) then
        jlo=1
     else
        jlo=0
     endif
     return
  endif
  ascnd=xx(n).ge.xx(1)
  if(jlo.le.0.or.jlo.gt.n)then
     jlo=0
     jhi=n+1
     goto 3
  endif
  inc=1
  if(x.ge.xx(jlo).eqv.ascnd)then
1    continue
     jhi=jlo+inc
     if(jhi.gt.n)then
        jhi=n+1
     else if(x.ge.xx(jhi).eqv.ascnd)then
        jlo=jhi
        inc=inc+inc
        goto 1
     endif
  else
     jhi=jlo
2    continue
     jlo=jhi-inc
     if(jlo.lt.1)then
        jlo=0
     else if(x.lt.xx(jlo).eqv.ascnd)then
        jhi=jlo
        inc=inc+inc
        goto 2
     endif
  endif
3 continue
  if(jhi-jlo == 1)then
     if(x == xx(n))jlo=n
     if(x == xx(1))jlo=1
     return
  endif
  jm=(jhi+jlo)/2
  if(x.ge.xx(jm).eqv.ascnd)then
     jlo=jm
  else
     jhi=jm
  endif
  goto 3
END SUBROUTINE hunt

