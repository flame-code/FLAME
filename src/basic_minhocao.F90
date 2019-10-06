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

!**********************************************************************************************

 subroutine dot_p(a,b,dotp)
 !a very simple implementation of the dot product
 implicit none
 real(8)::a(3),b(3)
 real(8)::dotp(3)
 integer::i
 dotp=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
 end subroutine

!**********************************************************************************************

subroutine bin_write(filename,array,n)
implicit none
integer:: n
real(8):: array(n)
character(40):: filename
open(unit = 120, status='replace',file=trim(filename),form='unformatted')  ! create a new file, or overwrite an existing on
write(120) array
close(120) ! close the file
!!open(unit = 120, status='replace',file=trim(filename),form='formatted')  ! create a new file, or overwrite an existing on
!!write(120,*) array
!!close(120) ! close the file
end subroutine

!**********************************************************************************************

subroutine bin_read(filename,array,n)
implicit none
integer:: n
real(8):: array(n)
character(40):: filename
open(unit = 120, status='old',file=trim(filename),form='unformatted')  ! open an existing file
read(120) array ! read the data into array x, of the appropriate data type
close(120) ! close the file
end subroutine

function round(enerd,accur)
  implicit none
  real*8 enerd,accur,round
  integer*8 ii
  ii=enerd/accur
  round=ii*accur
  !           write(*,'(a,1pe24.17,1x,i17,1x,1pe24.17)') 'enerd,ii,round',enerd,ii,round
  return
end function round
 subroutine rotation(rotmat,angle,axe)
 !This subroutine will calculate the rotational matrix rotmat for a
 !3-dim vector around an axis 'axe' by the angle 'angle'.
 implicit none
 real(8),INTENT(IN) :: angle
 real(8),INTENT(IN) :: axe(3)
 real(8):: rotator(3,3)
 real(8):: rotmat(3,3),cosang,sinang
 cosang=dcos(angle)
 sinang=dsin(angle)


 !Define Rotation Matrix
 rotator(1,1)=cosang+(axe(1)**2)*(1.d0-cosang)
 rotator(1,2)=axe(1)*axe(2)*(1.d0-cosang)-axe(3)*sinang
 rotator(1,3)=axe(1)*axe(3)*(1.d0-cosang)+axe(2)*sinang
                                                
 rotator(2,1)=axe(2)*axe(1)*(1.d0-cosang)+axe(3)*sinang
 rotator(2,2)=cosang+(axe(2)**2)*(1.d0-cosang) 
 rotator(2,3)=axe(2)*axe(3)*(1.d0-cosang)-axe(1)*sinang
                                              
 rotator(3,1)=axe(3)*axe(1)*(1.d0-cosang)-axe(2)*sinang
 rotator(3,2)=axe(3)*axe(2)*(1.d0-cosang)+axe(1)*sinang
 rotator(3,3)=cosang+(axe(3)**2)*(1.d0-cosang)
 rotmat(:,:)=rotator(:,:)

 !do i=1,3
 !   vector2(i)=rotator(i,1)*vector(1)+rotator(i,2)*vector(2)+rotator(i,3)*vector(3)
 !enddo
 !vector(:)=vector2(:)
 end subroutine rotation

!**********************************************************************************************

 subroutine invertmat(mat,matinv,n)
 implicit none
 real(8),intent(in) :: mat(n,n)
 integer               :: n
 real(8),allocatable   :: WORK(:)
 real(8)               :: matinv(n,n),det(3),a(n,n),div
 integer               :: IPIV(n), INFO
 integer               :: LDWORK
 !Here only for a 3*3 matrix
 if (n==3) then
 a=mat
 div=(a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)+a(1,2)*a(2,3)*a(3,1)+&
 &a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)) 
 div=1.d0/div
      matinv(1,1) = (a(2,2)*a(3,3)-a(2,3)*a(3,2))*div
      matinv(1,2) =-(a(1,2)*a(3,3)-a(1,3)*a(3,2))*div
      matinv(1,3) = (a(1,2)*a(2,3)-a(1,3)*a(2,2))*div
      matinv(2,1) =-(a(2,1)*a(3,3)-a(2,3)*a(3,1))*div
      matinv(2,2) = (a(1,1)*a(3,3)-a(1,3)*a(3,1))*div
      matinv(2,3) =-(a(1,1)*a(2,3)-a(1,3)*a(2,1))*div
      matinv(3,1) = (a(2,1)*a(3,2)-a(2,2)*a(3,1))*div
      matinv(3,2) =-(a(1,1)*a(3,2)-a(1,2)*a(3,1))*div
      matinv(3,3) = (a(1,1)*a(2,2)-a(1,2)*a(2,1))*div
 else
 !General n*n matrix 
 matinv=mat
 allocate(WORK(n))
 call  DGETRF( n, n, matinv, n, IPIV, INFO )
 if (info.ne.0) stop "Error in DGETRF"
 LDWORK=-1
 call  DGETRI( n, matinv, n, IPIV, WORK,LDWORK , INFO )
 LDWORK=WORK(1)
 deallocate(WORK)
 allocate(WORK(LDWORK))
 call  DGETRI( n, matinv, n, IPIV, WORK,LDWORK , INFO )
 if (info.ne.0) stop "Error in DGETRI"
 endif
 end subroutine

!**********************************************************************************************

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
