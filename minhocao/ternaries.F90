program binaries
!This program computes all possible binary compositions up to a maximum of NMAX
implicit none
integer:: NMAX,i,j,k,l,n,h
real(8):: con12_comp(3,10000),con12,pcon1
real(8):: con23_comp(3,10000),con23,pcon2
real(8):: con31_comp(3,10000),con31,pcon3
real(8):: pcon_lst(8,400,1000)
real(8),allocatable:: printarray(:,:)
character(30):: formatting
real(8):: vec12(2),vec13(2),convec23(2)
real(8):: rotmat(3,3),ep(2,3),p(2),pdist
logical:: found,inside
integer:: corner(3,3)
real(8):: pcorner(2,3)

!Composition 1 is at the left bottom corner
!Composition 2 is at the right bottom corner
!Composition 3 is at the top corner

con12_comp=1.d10
con23_comp=1.d10
con31_comp=1.d10
pcon_lst=0.d0

write(*,*) "Maximum number of atoms per cell"
read(*,*) NMAX
write(*,*) "Define segment with 3 corners of a triangle (0 0 0 for elementary corners)"
read(*,*) corner(:,1)
read(*,*) corner(:,2)
read(*,*) corner(:,3)
call getpoint(1,1,1,p,ep)
do i=1,3
if(all(corner(:,i).le.0)) then 
  pcorner(:,i)=ep(:,i)
else
  call getpoint(corner(1,i),corner(2,i),corner(3,i),p,ep)
  pcorner(:,i)=p
endif 
enddo
!write(*,*) pcorner(:,1)
!write(*,*) pcorner(:,2)
!write(*,*) pcorner(:,3)

allocate(printarray(3,NMAX))
n=0
do i=0,NMAX
   do j=NMAX-i,0,-1
      do h=NMAX-j-i,0,-1
        if(.not.((i==0.and.j==0.and.h==0))) then
        call getpoint(i,j,h,p,ep)
        call getconc(i,j,h,con12,con23,con31) 
        call getpconc(i,j,h,pcon1,pcon2,pcon3,ep)
        call checkpoint(p,pcorner(:,1),pcorner(:,2),pcorner(:,3),inside)
!Check if this concentration has already been found
!pcon_list: 1-3=pcon,4-6=integer concentration,7=number of compositions found of given pcon
        found=.false.
        if(inside) then
        do k=1,n
           pdist=0.d0
           pdist=pdist+(pcon_lst(1,k,1)-pcon1)**2
           pdist=pdist+(pcon_lst(2,k,1)-pcon2)**2
           pdist=pdist+(pcon_lst(3,k,1)-pcon3)**2
           if((pdist.lt.1.d-15)) then
              pcon_lst(7,  k,1)= pcon_lst(7,k,1)+1.d0 
              pcon_lst(1:3,k,int(pcon_lst(7,k,1)))=(/pcon1,pcon2,pcon3/)
              pcon_lst(4:6,k,int(pcon_lst(7,k,1)))=(/i,j,h/)
              found=.true.
           endif
        enddo
        if(.not.found) then
           n=n+1
           pcon_lst(7  ,n,1)= pcon_lst(7,n,1)+1.d0 
           pcon_lst(1:3,n,int(pcon_lst(7,n,1)))=(/pcon1,pcon2,pcon3/)
           pcon_lst(4:6,n,int(pcon_lst(7,n,1)))=(/i,j,h/)
        endif
!        write(*,'(3i5,2(es15.7),6(es10.3))') i,j,h,p,con12,con23,con31,pcon1/(ep(2,3)-ep(2,1)),pcon2/(ep(2,3)-ep(2,1)),pcon3/(ep(2,3)-ep(2,1))
        endif
        endif
enddo
enddo
enddo


open(unit=44,file='ternary.plt')
k=0
l=1
printarray(:,1)=con12_comp(:,1)
do i=1,n
     write(formatting,'(a,i5,a)') "(3(e12.4),",int(pcon_lst(7,i,1)),"(5x,3(i3)))"
     write(*,formatting) pcon_lst(1:3,i,1),int(pcon_lst(4:6,i,1:int(pcon_lst(7,i,1))))
     call getpoint(int(pcon_lst(4,i,1)),int(pcon_lst(5,i,1)),int(pcon_lst(6,i,1)),p,ep)
     write(44,'(3(e12.4))') p
enddo
close(44)

end program

subroutine getconc(i,j,h,con12,con23,con31) 
implicit none
real(8):: con12,con23,con31
integer:: i,j,h
con12=0.d0
con23=0.d0
con31=0.d0
   if((real(i,8)+real(j,8).gt.1.d-10))     con12=real(i,8)/(real(i,8)+real(j,8))
   if((real(j,8)+real(h,8).gt.1.d-10))     con23=real(j,8)/(real(j,8)+real(h,8))
   if((real(i,8)+real(h,8).gt.1.d-10))     con31=real(h,8)/(real(i,8)+real(h,8))
end subroutine

subroutine getpoint(i,j,h,point,ep)
implicit none
integer:: i,j,h
real(8):: k,rotmat(3,3),p(3),point(2),ep(2,3)
k=real(i+j+h,8)
p=(/real(i,8),real(j,8),real(h,8)/)/k
!Rotate plane normal to z maxis
call getrotmat(rotmat,ep)
p=matmul(rotmat,p)
point=p(1:2)
end subroutine

subroutine getpconc(i,j,h,pcon1,pcon2,pcon3,ep)
implicit none
integer:: i,j,h
real(8):: k,rotmat(3,3),p(3),point(2),ep(2,3)
real(8):: pcon1,pcon2,pcon3,fac
k=real(i+j+h,8)
p=(/real(i,8),real(j,8),real(h,8)/)/k
fac=sqrt(3.d0*0.5d0)
pcon1=p(1)*fac
pcon2=p(2)*fac
pcon3=p(3)*fac
pcon1=pcon1/(ep(2,3)-ep(2,1))
pcon2=pcon2/(ep(2,3)-ep(2,1))
pcon3=pcon3/(ep(2,3)-ep(2,1))
end subroutine


subroutine checkpoint(point,point1,point2,point3,inside)
!Check if point is in the trianle defined by point1,2,3
implicit none
integer:: i,j,h
real(8),dimension(2):: point,point1,point2,point3
real(8),dimension(2):: t1,t2,t3
real(8),dimension(3):: p1,p2,p3,c1,c2
real(8):: angle,p(2,3)
real(8):: pi
logical:: inside
inside=.true.
pi=acos(-1.d0)
!Compute and sum up the angles
angle=0.d0
p(:,1)=point1
p(:,2)=point2
p(:,3)=point3
do i=1,3
   j=mod(i,3)+1
   h=mod(i+1,3)+1
   t1=point-p(:,i)
   t2=p(:,j)-p(:,i)
   t3=p(:,h)-p(:,i)
   p1=(/t1,0.d0/)
   p2=(/t2,0.d0/)
   p3=(/t3,0.d0/)
   call cross_product(p2,p1,c1)
   call cross_product(p2,p3,c2)
   if(dot_product(c1,c2).lt.-1.d-10) inside=.false.
enddo
end subroutine


subroutine getrotmat(rotmat,ep)
implicit none
real(8):: k,angle,ep(2,3)
real(8),dimension(3):: a,b,c,normal,axe,tmp
real(8),dimension(3,3):: rotmat1,rotmat2,rotmat
a=(/1.d0,0.d0,0.d0/)
b=(/0.d0,1.d0,0.d0/)
c=(/0.d0,0.d0,1.d0/)
!Rotate plane normal to z maxis
normal=(/1.d0,1.d0,1.d0/)
normal=normal/sqrt(dot_product(normal,normal))
call cross_product(normal,c,axe)
axe=axe/sqrt(dot_product(axe,axe))
angle=acos(dot_product(normal,c))
call rotation(rotmat1,angle,axe)
axe=-c
tmp=matmul(rotmat1,b-a)
tmp=tmp/sqrt(dot_product(tmp,tmp))
angle=acos(dot_product(tmp,a))
call rotation(rotmat2,angle,axe)
rotmat=matmul(rotmat2,rotmat1)
a=matmul(rotmat,a)
b=matmul(rotmat,b)
c=matmul(rotmat,c)
ep(:,1)=a(1:2)
ep(:,2)=b(1:2)
ep(:,3)=c(1:2)
end subroutine


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

!************************************************************************************

 subroutine rotation(rotmat,angle,axe)
 !This subroutine will calculate the rotational matrix rotmat for a
 !3-dim vector around an axis 'axe' by the angle 'angle'.
 implicit none
 real(8),INTENT(IN) :: angle
 real(8),INTENT(IN) :: axe(3)
 real(8):: rotator(3,3)
 real(8):: rotmat(3,3)

 !Define Rotation Matrix
 rotator(1,1)=dcos(angle)+(axe(1)**2)*(1.d0-dcos(angle))
 rotator(1,2)=axe(1)*axe(2)*(1.d0-dcos(angle))-axe(3)*dsin(angle)
 rotator(1,3)=axe(1)*axe(3)*(1.d0-dcos(angle))+axe(2)*dsin(angle)

 rotator(2,1)=axe(2)*axe(1)*(1.d0-dcos(angle))+axe(3)*dsin(angle)
 rotator(2,2)=dcos(angle)+(axe(2)**2)*(1.d0-dcos(angle))
 rotator(2,3)=axe(2)*axe(3)*(1.d0-dcos(angle))-axe(1)*dsin(angle)

 rotator(3,1)=axe(3)*axe(1)*(1.d0-dcos(angle))-axe(2)*dsin(angle)
 rotator(3,2)=axe(3)*axe(2)*(1.d0-dcos(angle))+axe(1)*dsin(angle)
 rotator(3,3)=dcos(angle)+(axe(3)**2)*(1.d0-dcos(angle))
 rotmat(:,:)=rotator(:,:)

 !do i=1,3
 !   vector2(i)=rotator(i,1)*vector(1)+rotator(i,2)*vector(2)+rotator(i,3)*vector(3)
 !enddo
 !vector(:)=vector2(:)
 end subroutine rotation

!************************************************************************************

 subroutine cross_product(a,b,crossp)
 !a very simple implementation of the cross product
 implicit none
 real(8)::a(3),b(3)
 real(8)::crossp(3)
 crossp(1)=a(2)*b(3)-a(3)*b(2)
 crossp(2)=a(3)*b(1)-a(1)*b(3)
 crossp(3)=a(1)*b(2)-a(2)*b(1)
 return
 end subroutine

 subroutine dot_p(a,b,dotp)
 !a very simple implementation of the dot product
 implicit none
 real(8)::a(3),b(3)
 real(8)::dotp(3)
 integer::i
 dotp=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
 end subroutine

!************************************************************************************
