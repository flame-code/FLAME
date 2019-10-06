!Various methods to initialize the velocities for the MD part of Minima Hopping
!GAUSSIAN DISTRIBUTION**********************************************************
      subroutine gausdist(nat,vxyz,amass)
!generates 3*nat random numbers distributed according to  exp(-.5*vxyz**2)
      implicit none!real*8 (a-h,o-z)
      real:: s1,s2
      real(8):: t1,t2,tt,amass(nat)
! On Intel the random_number can take on the values 0. and 1.. To prevent overflow introduce eps
      real(8),parameter:: eps=1.d-8
      real(8),dimension(3*nat)::  vxyz
      integer:: nat,i
      do i=1,3*nat-1,2
        call random_number(s1)
        t1=eps+(1.d0-2.d0*eps)*dble(s1)
        call random_number(s2)
        t2=dble(s2)
        tt=sqrt(-2.d0*log(t1))
        vxyz(i)=tt*cos(6.28318530717958648d0*t2)
        vxyz(i+1)=tt*sin(6.28318530717958648d0*t2)
      enddo
        call random_number(s1)
        t1=eps+(1.d0-2.d0*eps)*dble(s1)
        call random_number(s2)
        t2=dble(s2)
        tt=sqrt(-2.d0*log(t1))
        vxyz(3*nat)=tt*cos(6.28318530717958648d0*t2)

        call elim_moment(nat,vxyz,amass)
      return
      end subroutine

!************************************************************************************

!GAUSSIAN DISTRIBUTION FOR THE CELL VECTORS***************************************
      subroutine gausdist_cell(latvec,vlat)
! generates 3*3 random numbers distributed according to  exp(-.5*vxyz**2) for the cell vectors
      implicit none
      integer:: i
      real:: s1,s2
      real(8) :: t1,t2,tt
! On Intel the random_number can take on the values 0. and 1.. To prevent overflow introduce eps
      real(8),parameter:: eps=1.d-8
      real(8)::  vlat(9),latvec(9)

      do i=1,3*3-1,2
        call random_number(s1)
        t1=eps+(1.d0-2.d0*eps)*dble(s1)
        call random_number(s2)
        t2=dble(s2)
        tt=sqrt(-2.d0*log(t1))
        vlat(i)=tt*cos(6.28318530717958648d0*t2)
        vlat(i+1)=tt*sin(6.28318530717958648d0*t2)
      enddo
        call random_number(s1)
        t1=eps+(1.d0-2.d0*eps)*dble(s1)
        call random_number(s2)
        t2=dble(s2)
        tt=sqrt(-2.d0*log(t1))
        vlat(3*3)=tt*cos(6.28318530717958648d0*t2)

        call elim_torque_cell(latvec,vlat)
      return
      end subroutine








!Various methods to initialize the velocities for the MD part of Minima Hopping
!GAUSSIAN DISTRIBUTION**********************************************************
      subroutine gausdist_simple(nat,vxyz,amass)
!generates 3*nat random numbers distributed according to  exp(-.5*vxyz**2)
      use mod_utils
      implicit none!real*8 (a-h,o-z)
      real(8):: s1,s2
      real(8):: t1,t2,tt,amass(nat)
! On Intel the random_number can take on the values 0. and 1.. To prevent overflow introduce eps
      real(8),parameter:: eps=1.d-8
      real(8),dimension(3*nat)::  vxyz
      integer:: nat,i
      do i=1,3*nat-1,2
        call random_number_generator_simple(s1)
        t1=eps+(1.d0-2.d0*eps)*dble(s1)
        call random_number_generator_simple(s2)
        t2=dble(s2)
        tt=sqrt(-2.d0*log(t1))
        vxyz(i)=tt*cos(6.28318530717958648d0*t2)
        vxyz(i+1)=tt*sin(6.28318530717958648d0*t2)
      enddo
        call random_number_generator_simple(s1)
        t1=eps+(1.d0-2.d0*eps)*dble(s1)
        call random_number_generator_simple(s2)
        t2=dble(s2)
        tt=sqrt(-2.d0*log(t1))
        vxyz(3*nat)=tt*cos(6.28318530717958648d0*t2)

        call elim_moment(nat,vxyz,amass)
      return
      end subroutine

!************************************************************************************

!GAUSSIAN DISTRIBUTION FOR THE CELL VECTORS***************************************
      subroutine gausdist_cell_simple(latvec,vlat)
! generates 3*3 random numbers distributed according to  exp(-.5*vxyz**2) for the cell vectors
      use mod_utils
      implicit none
      integer:: i
      real(8):: s1,s2
      real(8) :: t1,t2,tt
! On Intel the random_number can take on the values 0. and 1.. To prevent overflow introduce eps
      real(8),parameter:: eps=1.d-8
      real(8)::  vlat(9),latvec(9)

      do i=1,3*3-1,2
        t1=eps+(1.d0-2.d0*eps)*dble(s1)
        call random_number_generator_simple(s2)
        t2=dble(s2)
        tt=sqrt(-2.d0*log(t1))
        vlat(i)=tt*cos(6.28318530717958648d0*t2)
        vlat(i+1)=tt*sin(6.28318530717958648d0*t2)
      enddo
        call random_number_generator_simple(s1)
        t1=eps+(1.d0-2.d0*eps)*dble(s1)
        call random_number_generator_simple(s2)
        t2=dble(s2)
        tt=sqrt(-2.d0*log(t1))
        vlat(3*3)=tt*cos(6.28318530717958648d0*t2)

        call elim_torque_cell(latvec,vlat)
      return
      end subroutine










