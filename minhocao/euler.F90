!!!program test
!!!implicit none
!!!integer:: i
!!!real(8):: A(3,3),B(3,3),C(3,3),alpha,beta,gamma,vec(3),pi,phi,tmp1(3),vec2(3),rot1(3),rot2(3),quat(4),quatb(4),quatout(4)
!!!pi=acos(-1.d0)
!!!
!!!call random_number(alpha)
!!!call random_number(beta)
!!!call random_number(gamma)
!!!
!!!!!!alpha=pi/4.d0
!!!!!!beta=0.d0!pi/4.d0
!!!!!!gamma=0.d0
!!!
!!!!write(*,*) alpha,beta,gamma
!!!call euler2rotmat(A, alpha,beta,gamma)
!!!write(*,*) A(1,:)
!!!write(*,*) A(2,:)
!!!write(*,*) A(3,:)
!!!write(*,'(a,3(3(es25.15),a))') "A=[",A(1,:),";",A(2,:),";",A(3,:),"]"
!!!
!!!call random_number(alpha)
!!!call random_number(beta)
!!!call random_number(gamma)
!!!call euler2rotmat(B, alpha,beta,gamma)
!!!write(*,*) B(1,:)
!!!write(*,*) B(2,:)
!!!write(*,*) B(3,:)
!!!write(*,'(a,3(3(es25.15),a))') "B=[",B(1,:),";",B(2,:),";",B(3,:),"]"
!!!
!!!c=matmul(A,B)
!!!
!!!call rotmat2quat(A,quat)
!!!write(*,*) "QUAT"
!!!write(*,'(4(es25.15))')quat
!!!call quat2rotmat(A,quat)
!!!write(*,*) "ROTMAT"
!!!write(*,'(3(es25.15))') A(1,:)
!!!write(*,'(3(es25.15))') A(2,:)
!!!write(*,'(3(es25.15))') A(3,:)
!!!
!!!call rotmat2quat(B,quatb)
!!!write(*,*) "QUAT_B"
!!!write(*,'(4(es25.15))')quatb
!!!
!!!call quatxquat(quat,quatb,quatout)
!!!write(*,*) "QUATxQUAT_B"
!!!write(*,'(4(es25.15))')quatout
!!!write(*,'(a,3(3(es25.15),a))') "C    =[",B(1,:),";",B(2,:),";",B(3,:),"]"
!!!call quat2rotmat(C,quatout)
!!!write(*,'(a,3(3(es25.15),a))') "Cquat=[",B(1,:),";",B(2,:),";",B(3,:),"]"
!!!
!!!
!!!
!!!call rotmat2euler(A,alpha,beta,gamma)
!!!write(*,*) alpha,beta,gamma
!!!
!!!vec=(/1,0,0/)
!!!do i=1,1
!!!!  call uniform_euler(A,alpha,beta,gamma)
!!!  call rand_sphere(vec,pi)
!!!!  vec(3)=abs(vec(3))
!!!  call random_number(alpha)
!!!  alpha=(alpha-0.5d0)*4.d0*pi
!!!!  alpha=pi 
!!!  call rotmat2vecang(A,alpha,vec)
!!!  rot1=matmul(A,vec)
!!!  write(*,'(a,4(es25.15))') "Vec,alpha 1 : ",vec,alpha
!!!  write(*,'(a,3(3(es20.8),a))') "A=[",A(1,:),";",A(2,:),";",A(3,:),"]"
!!!  write(14,'(a,3(3(es20.8),a))') "A=[",A(1,:),";",A(2,:),";",A(3,:),"];"
!!!  write(14,*) "vrrotmat2vec(A)"
!!!  call rotmat2vecang(A,vec2,phi)
!!!  write(13,*) vec(:)
!!!  write(*,'(a,4(es25.15))') "Vec,alpha 2 : ",vec2,phi
!!!  beta= dot_product(vec2-vec,vec2-vec)
!!!  gamma=dot_product(vec2+vec,vec2+vec)
!!!  call vecang2rotmat(A,phi,vec2)
!!!  rot2=matmul(A,vec)
!!!  if(dot_product(rot1-rot2,rot1-rot2).gt.1.d-10) stop "Its the end of it!"
!!!!!  if((abs(alpha-phi).gt.1.d-10.and.abs(alpha+phi).gt.1.d-10).or.(beta.gt.1.d-10.and.gamma.gt.1.d-10)) then
!!!!!!  if(beta.gt.1.d-12.and.gamma.gt.1.d-12) then
!!!!!    write(*,'(a,3(es15.7))') "Wrong angle ",abs(alpha-phi),beta,gamma
!!!!!    stop
!!!!!  endif
!!!enddo
!!!!
!!!!call rotmat2vecang(A,vec,phi)
!!!!write(*,*) vec
!!!!write(*,*) phi
!!!end program

subroutine rotmat2vecang(A,axis,phi)
implicit none
integer:: i
real(8):: A(3,3),axis(3),phi,Atrc,eps,A_upper(3),signs(3),flip(3),shifted(3),pi,den
pi=acos(-1.d0)
!VRROTMAT2VEC Convert rotation from matrix to axis-angle representation.
!   R = VRROTMAT2VEC(M) returns an axis-angle representation of  
!   rotation defined by the rotation matrix M.
!
!   R = VRROTMAT2VEC(M, OPTIONS) converts the rotation with the default 
!   algorithm parameters replaced by values defined in the structure
!   OPTIONS.
!
!   The OPTIONS structure contains the following parameters:
!
!     'epsilon'
!        Minimum value to treat a number as zero. 
!        Default value of 'epsilon' is 1e-12.
!
!   The result R is a 4-element axis-angle rotation row vector.
!   First three elements specify the rotation axis, the last element
!   defines the angle of rotation.
!
!   See also VRROTVEC2MAT, VRROTVEC.

!   Copyright 1998-2008 HUMUSOFT s.r.o. and The MathWorks, Inc.
!   $Revision: 1.1.6.1 $ $Date: 2008/10/31 07:09:44 $ $Author: batserve $

! Rotation axis elements flipping in the singular case phi == pi 
! is discussed here:
! www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToAngle

! test input arguments

  eps = 1.d-12
!! make the conversion
Atrc = A(1,1)+A(2,2)+A(3,3)
if (abs(Atrc - 3.d0) .le. eps) then
  axis=(/0.d0,1.d0,0.d0/)
  phi=0.d0
elseif  (abs(Atrc + 1) .le. eps) then
!   phi == pi
!   This singularity requires elaborate sign ambiguity resolution
!  
!   compute axis of rotation, make sure all elements are >= 0
!   real signs are obtained by the flipping algorithm below
  axis(1)=sqrt(max(0.5d0*(A(1,1)+1.d0),0.d0))
  axis(2)=sqrt(max(0.5d0*(A(2,2)+1.d0),0.d0))
  axis(3)=sqrt(max(0.5d0*(A(3,3)+1.d0),0.d0))
!
!  axis elements that are <= epsilon set to zero
!  
   do i=1,3
      if(axis(i).le.eps) axis(i)=0.d0
   enddo
!   Flipping
!   
!   The algorithm uses the elements above diagonal to determine the signs 
!   of rotation axis coordinate in the singular case Phi = pi. 
!   All valid combinations of 0, positive and negative values lead to 
!   3 different cases:
!   If (Sum(signs)) >= 0 ... leave all coordinates positive
!   If (Sum(signs)) == -1 and all values are non-zero 
!     ... flip the coordinate that is missing in the term that has + sign, 
!         e.g. if 2AyAz is positive, flip x
!   If (Sum(signs)) == -1 and 2 values are zero 
!     ... flip the coord next to the one with non-zero value 
!     ... ambiguous, we have chosen shift right
!  
!   construct vector [M23 M13 M12] ~ [2AyAz 2AxAz 2AxAy]
!   (in the order to facilitate flipping):    ^
!                                    [no_x  no_y  no_z ]
!  
   A_upper = (/A(2,3),A(1,3),A(1,2)/)
!   elements with || smaller than epsilon are considered to be 0
   signs=0.d0
   do i=1,3
    if(abs(A_upper(i)) > eps) signs(i) = sign(1.d0,A_upper(i)) 
   enddo
!  
  if (sum(signs) >= 0.d0) then
!   none of the signs is negative
!    
!     don't flip any axis element
    flip = (/1.d0, 1.d0, 1.d0/)
!    
  elseif (.not.any(signs.eq.0.d0)) then
!   none of the signs is zero, 2 negative 1 positive
!    
!     flip the coordinate that is missing in the term that has + sign
    flip = -signs;
!    
  else
write(*,*) "In if 1c"
stop
!   2 signs are 0, 1 negative
!      flip the coord to the right of the one with non-zero value 
!      [-1 0 0]->[1 -1 1], [0 -1 0]->[1 1 -1], [0 0 -1]->[-1 1 1]
     shifted = (/signs(3), signs(1), signs(2)/)
     flip=shifted
     do i=1,3
      if(shifted(i)==0.d0) flip(i) = shifted(i) + 1.d0
     enddo
  endif
!   flip the axis
  axis = axis * flip
  axis = axis
  phi= pi
else
!   general case
  phi = acos((Atrc - 1.d0)*0.5d0);
  den = 2.d0*sin(phi)
  axis = (/A(3,2)-A(2,3), A(1,3)-A(3,1), A(2,1)-A(1,2)/) / den 
  axis = axis
  phi = phi
endif
end subroutine


subroutine uniform_euler(A,alpha,beta,gamma,vec,theta)
!Using the algorithm to generate "Fast random rotation matrices"
!http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.53.1357
implicit none
real(8):: alpha,beta,gamma,R(3,3),pi,c,s,A(3,3),rand(3),sq,v(3),vec(3),theta
pi=acos(-1.d0)
A=0.d0
R=0.d0
call random_number(rand)
c=cos(2.d0*pi*rand(1))
s=sin(2.d0*pi*rand(1))
R(3,3)=1.d0
R(1,1)=c
R(2,2)=c
R(1,2)=-s
R(2,1)=s
c=cos(2.d0*pi*rand(2))
s=sin(2.d0*pi*rand(2))
sq=sqrt(rand(3))
v(1)=c*sq
v(2)=s*sq
v(3)=sqrt(1.d0-rand(3))
A(1,1)=v(1)**2;A(1,2)=v(1)*v(2);A(1,3)=v(1)*v(3)
A(2,1)=A(1,2) ;A(2,2)=v(2)**2  ;A(2,3)=v(2)*v(3)
A(3,1)=A(1,3) ;A(3,2)=A(2,3)   ;A(3,3)=v(3)**2
A=-A*2.d0
A(1,1)=A(1,1)+1.d0
A(2,2)=A(2,2)+1.d0
A(3,3)=A(3,3)+1.d0
A=-matmul(A,R)
call rotmat2euler(A, alpha,beta, gamma)

end subroutine



subroutine euler2rotmat(A, alpha,beta,gamma)
!This routine will compute the matrix A to perform a euler rotation according to YZY'
!http://en.citizendium.org/wiki/Euler_angles
implicit none
integer:: i
real(8):: A(3,3),c1,c2,c3,s1,s2,s3,alpha,beta,gamma
c1=cos(alpha)
c2=cos(beta)
c3=cos(gamma)
s1=sin(alpha)
s2=sin(beta)
s3=sin(gamma)
!zyz convention
A(1,1)=c1*c2*c3-s1*s3
A(2,1)=s1*c2*c3+c1*s3
A(3,1)=-s2*c3
A(1,2)=-c1*c2*s3-s1*c3
A(2,2)=-s1*c2*s3+c1*c3
A(3,2)=s2*s3
A(1,3)=c1*s2
A(2,3)=s1*s2
A(3,3)=c2
!!zyz convention
!A(1,1)=-s1*s3+c1*c2*c3
!A(2,1)=-s1*c3-c1*c2*s3
!A(3,1)=c1*s2
!A(1,2)=c1*s3+s1*c2*c3
!A(2,2)=c1*c3-s1*c2*s3
!A(3,2)=s1*s2
!A(1,3)=-s2*c3
!A(2,3)=s2*s3
!A(3,3)=c2
end subroutine

subroutine rotmat2euler(A, alpha,beta, gamma)
!     Determine Euler angles of the proper rotation matrix A. YZY'
!     All conventions according to Biedenharn & Louck
!     Author P.E.S. Wormer 1985
      implicit none
      real*8:: alpha,beta,gamma
      real*8:: A(3,3),AA(3,3)

      AA=transpose(A)
      AA=A

!--Check if matrix is orthogonal with unit determinant.
      call Check(AA)

!--Get polar angles of third column
      call get_polar(AA(:,3), alpha, beta)

!--Compute gamma
      call Comgamma(AA, alpha, beta, gamma)
end subroutine rotmat2euler

      subroutine Check(A)

!     Check if matrix is orthogonal with unit determinant.

      implicit none
!      implicit real*8(a-h,o-z)
      integer:: i,j,k
      real*8,parameter:: thresh = 1.d-12
      real*8:: A(3,3),t,det
      do i=1,3
         do j=1,i-1
            t = 0.d0
            do k=1,3
               t = t + A(i,k)*A(j,k)
            enddo
            if (abs(t) .gt. thresh) then
               write(*,'(1x,a,i1,a,i1,a,d12.5 )') ' Row ',i, ' and ', j, 'non-orthogonal',abs(t)
               stop 'non-orthogonal'
            endif
         enddo
         t = 0.d0
         do k=1,3
            t = t + A(i,k)*A(i,k)
         enddo
         if (abs(t-1.d0) .gt. thresh) then   
            write(*,'(1x,a,i1,a,d25.15 )')       ' Row ',i, ' non-normalized:' , t
            stop 'non-normalized'
         endif
      enddo

      t = det(A)
      if (abs(t-1.d0) .gt. thresh) then
         if (abs(t+1.d0) .lt. thresh) then
            write(*,'(//1x,a,1x,d14.6,a)')      'Non-unit determinant:',t, ' will interchange column 1 and 2'
            do i=1,3
               T = A(i,2)
               A(i,2) = A(i,1)
               A(i,1) = T
            enddo
         else
         write(*,'(//1x,a,d20.10,a)')      'Non-unit determinant:',t
            stop ' Determinant'
         endif
      endif

      end subroutine
      real*8 function det(A)
      implicit none
      integer:: i
!     determinant of A
      real*8:: A(3,3), B(3),pi

      B(1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
      B(2) = A(3,2)*A(1,3) - A(1,2)*A(3,3)
      B(3) = A(1,2)*A(2,3) - A(2,2)*A(1,3)
      det = 0.d0
      do i=1,3
         det = det + A(i,1)*B(i)
      enddo
      end function
      Subroutine get_polar(A, alpha, beta)

!     Get polar angles of vector A.
      implicit none
      integer:: i
      real*8, parameter:: thresh = 1.d-50
      real*8:: A(3),alpha,beta,R,cb,pi
      R = sqrt( A(1)**2 + A(2)**2 + A(3)**2 )
      if ( abs(R) .lt. thresh) then
         write(*,'(a)') ' zero vector'
         alpha = 0.d0
         beta  = 0.d0
         return
      endif

      beta = acos(A(3)/R);write(*,*) "beta",A(3)/R
      cb   = abs (A(3)/R)
      if ( abs(cb-1.d0) .lt. thresh ) then
         alpha = 0.d0
         return
      endif
      alpha = acos( A(1) / (sqrt( A(1)**2 + A(2)**2 ) ) )
      if ( A(2) .lt. 0.d0 ) then
         pi = acos(-1.d0)
         alpha = 2.d0*pi - alpha
      endif

      end subroutine

      Subroutine Comgamma(A, alpha,beta, gamma)

!     Compute third Euler angle gamma.
      implicit none
      integer:: i
      real*8, parameter:: thresh = 1.d-50
      real*8:: A(3,3), B1(3), B2(3),alpha,beta,gamma,cb,sb,ca,sa,sg,pi,cg

      cb = cos(beta)
      sb = sin(beta)
      ca = cos(alpha)
      sa = sin(alpha)
      B1(1) = ca*cb
      B1(2) = sa*cb
      B1(3) = -sb
      B2(1) = -sa
      B2(2) =  ca
      B2(3) =  0.d0
      cg = 0.d0
      sg = 0.d0
      do i=1,3
        cg = cg + B1(i)*A(i,1)
        sg = sg + B2(i)*A(i,1)
      enddo
      gamma = acos(cg)
      if ( sg  .lt. 0.d0 ) then
         pi = acos(-1.d0)
         gamma = 2.d0*pi - gamma
      endif

      end subroutine  

!************************************************************************************

subroutine rand_sphere(radxyz,pi)
!random spherical distribution
!This routine will compute a random uniformli distributed point on a sphere of radius 1
!Uniform distribution of points on a sphere:
!http://mathworld.wolfram.com/SpherePointPicking.html
implicit none
real(8), intent(in):: pi
real(8):: u,v,theta,phi,radxyz(3),sinphi,costheta,cosphi,sintheta


  call RANDOM_NUMBER(v)
  call RANDOM_NUMBER(u)

  theta=u*2.d0*pi
  phi=acos(2.d0*v-1.d0)

  sinphi=sin(phi)
  costheta=cos(theta)
  cosphi=cos(phi)
  sintheta=sin(theta)

  radxyz(1)=costheta*sinphi
  radxyz(2)=sintheta*sinphi
  radxyz(3)=cosphi
end subroutine


!************************************************************************************

 subroutine vecang2rotmat(rotmat,angle,axe)
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
 end subroutine vecang2rotmat

!************************************************************************************

subroutine vecang2quat(axe,angle,quat)
!This subroutine will construct a quaternion representing the rotation
!of angle around an axis axe and store it in quat
implicit none
real(8):: axe(3),angle,quat(4)
quat(2:4)=axe/sqrt(dot_product(axe,axe))*sin(0.5d0*angle)
quat(1)=cos(0.5d0*angle)
end subroutine

!************************************************************************************

subroutine quat2rotmat(rotmat,quat)
!This subroutine will generate a rotation matrix from a quaternion reperesentation
implicit none
real(8):: quat(0:3),rotmat(3,3)

!call q2m ( quat,rotmat)
!!!     local variables
real(8)::     q01
real(8)::     q02
real(8)::     q03
real(8)::     q12
real(8)::     q13
real(8)::     q23
real(8)::     q1s
real(8)::     q2s
real(8)::     q3s
real(8)::     l2
real(8)::     sharpn
 
 
 
      q01  =  quat(0) * quat(1)
      q02  =  quat(0) * quat(2)
      q03  =  quat(0) * quat(3)
 
      q12  =  quat(1) * quat(2)
      q13  =  quat(1) * quat(3)
 
      q23  =  quat(2) * quat(3)
 
      q1s  =  quat(1) * quat(1)
      q2s  =  quat(2) * quat(2)
      q3s  =  quat(3) * quat(3)
 
     l2   =  quat(0) * quat(0) + q1s + q2s + q3s

     if ( l2 .ne. 1.0d0 .and. l2 .ne. 0.0d0 ) then

        sharpn = 1.0d0 / l2

        q01    = q01 * sharpn
        q02    = q02 * sharpn
        q03    = q03 * sharpn

        q12    = q12 * sharpn
        q13    = q13 * sharpn

        q23    = q23 * sharpn

        q1s    = q1s * sharpn
        q2s    = q2s * sharpn
        q3s    = q3s * sharpn

     end if



     rotmat(1,1) =  1.d0  -  2.d0 * ( q2s + q3s )
     rotmat(2,1) =           2.d0 * ( q12 + q03 )
     rotmat(3,1) =           2.d0 * ( q13 - q02 )

     rotmat(1,2) =           2.d0 * ( q12 - q03 )
     rotmat(2,2) =  1.d0  -  2.d0 * ( q1s + q3s )
     rotmat(3,2) =           2.d0 * ( q23 + q01 )

     rotmat(1,3) =           2.d0 * ( q13 + q02 )
     rotmat(2,3) =           2.d0 * ( q23 - q01 )
     rotmat(3,3) =  1.d0  -  2.d0 * ( q1s + q2s )

     return
     end subroutine



subroutine rotmat2quat(rotmat,quat)
!This subroutine will compute the quaternion from a rotation matrix
implicit none
real(8):: rotmat(3,3),quat(0:3)
!!!!call m2q(rotmat,quat)
!!!!end subroutine
!     local variables
      real(8)     c
      real(8)     cc4
      real(8)     factor
      real(8)     l2
      real(8)     mtrace
      real(8)     polish
      real(8)     s(3)
      real(8)     s114
      real(8)     s224
      real(8)     s334
      real(8)     trace
 
 
 
      trace  = rotmat(1,1) + rotmat(2,2) + rotmat(3,3)
      mtrace = 1.0d0 - trace
 
      cc4    = 1.0d0  + trace
      s114   = mtrace + 2.0d0*rotmat(1,1)
      s224   = mtrace + 2.0d0*rotmat(2,2)
      s334   = mtrace + 2.0d0*rotmat(3,3)
 
!     note that if you simply add cc4 + s114 + s224 + s334
!     you get four. thus at least one of the 4 terms is greater than 1.
!
      if ( 1.0d0 .le. cc4 ) then
         c      =  dsqrt  ( cc4 * 0.25d0 )
         factor =  1.0d0 /( c   * 4.0d0  )
 
         s(1)   = ( rotmat(3,2) - rotmat(2,3) )*factor
         s(2)   = ( rotmat(1,3) - rotmat(3,1) )*factor
         s(3)   = ( rotmat(2,1) - rotmat(1,2) )*factor
 
      else if ( 1.0d0 .le. s114 ) then
 
         s(1)   = dsqrt  ( s114 * 0.25d0 )
         factor = 1.0d0 /( s(1) * 4.0d0  )
 
         c      = ( rotmat(3,2) - rotmat(2,3) ) * factor
         s(2)   = ( rotmat(1,2) + rotmat(2,1) ) * factor
         s(3)   = ( rotmat(1,3) + rotmat(3,1) ) * factor
 
 
      else if ( 1.0d0 .le. s224 ) then
 
         s(2)   = dsqrt  ( s224 * 0.25d0 )
         factor = 1.0d0 /( s(2) * 4.0d0  )
 
         c      = ( rotmat(1,3) - rotmat(3,1) ) * factor
         s(1)   = ( rotmat(1,2) + rotmat(2,1) ) * factor
         s(3)   = ( rotmat(2,3) + rotmat(3,2) ) * factor
 
      else
 
         s(3)   = dsqrt  ( s334 * 0.25d0 )
         factor = 1.0d0 /( s(3) * 4.0d0  )
 
         c      = ( rotmat(2,1) - rotmat(1,2) ) * factor
         s(1)   = ( rotmat(1,3) + rotmat(3,1) ) * factor
         s(2)   = ( rotmat(2,3) + rotmat(3,2) ) * factor
 
      end if
!
!     if the magnitude of this quaternion is not one, we polish it
!     up a bit.
!
      l2 = c*c + s(1)*s(1) + s(2)*s(2) + s(3)*s(3)
 
      if ( l2 .ne. 1.0d0 ) then
         polish = 1.0d0/dsqrt(l2)
         c      =    c*polish
         s(1)   = s(1)*polish
         s(2)   = s(2)*polish
         s(3)   = s(3)*polish
      end if
 
      if ( c .gt. 0.0d0 ) then
         quat(0) = c
         quat(1) = s(1)
         quat(2) = s(2)
         quat(3) = s(3)
      else
         quat(0) = -c
         quat(1) = -s(1)
         quat(2) = -s(2)
         quat(3) = -s(3)
      end if
 
      return
      end subroutine


subroutine quatxquat(quat1,quat2,quatout)
!This routine will compute the product of two quaternions
real(8):: quat1(0:3),quat2(0:3),quatout(0:3),crossp(3)
quatout(0)=quat1(0)*quat2(0)-dot_product(quat1(1:3),quat2(1:3))
call cross_product(quat1(1:3),quat2(1:3),crossp)
quatout(1:3)=quat1(0)*quat2(1:3)+quat2(0)*quat1(1:3)+crossp
end subroutine

subroutine qinv(quat_in,quat_out)
!This routine will revert the rotation from quat_in to quat_out
implicit none
real(8):: quat_in(4),quat_out(4)
quat_out(1)=quat_in(1)
quat_out(2:4)=-quat_in(2:4)
end subroutine

