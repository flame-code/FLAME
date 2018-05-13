!!  module cell_utils
!!    implicit none
!!  
!!    private
!!    public ::          &
!!      cell_volume,     &
!!      backtocell_cart, &
!!      find_kpt,        &
!!      k_expansion
!!  
!!  contains
!!  
!!  
!!    !************************************************************************************
!!    real(8) function cell_volume(v)
!!      real(8), intent(in) :: v(:,:)
!!      
!!      cell_volume = &
!!        +v(1,1)*v(2,2)*v(3,3) &
!!        -v(1,1)*v(2,3)*v(3,2) &
!!        -v(1,2)*v(2,1)*v(3,3) &
!!        +v(1,2)*v(2,3)*v(3,1) &
!!        +v(1,3)*v(2,1)*v(3,2) &
!!        -v(1,3)*v(2,2)*v(3,1)
!!  
!!    end function cell_volume
!!  
!!  
!!    !************************************************************************************
!!    ! This subroutine will transform back all atoms into the periodic cell
!!    ! defined by the 3 lattice vectors in latvec=[v1.v2.v3]
!!    !************************************************************************************
!!    subroutine backtocell_cart(nat, latvec, rxyz)
!!      integer, intent(in) :: nat
!!      integer :: i,iat,j
!!      real(8) :: latvec(3,3), rxyz(3,nat), crossp(3),a(3),b(3), nvec(3,3), dist(6),eps,count,rxyz_red(3,nat)
!!      logical:: neccesary
!!  
!!      !First check if the volume is positive
!!      if(cell_volume(latvec).le.0.d0) stop "Negative volume during backtocell"
!!      call rxyz_cart2int(latvec,rxyz_red,rxyz,nat)
!!      do iat=1,nat
!!          rxyz_red(1,iat)=modulo(modulo(rxyz_red(1,iat),1.d0),1.d0)
!!          rxyz_red(2,iat)=modulo(modulo(rxyz_red(2,iat),1.d0),1.d0)
!!          rxyz_red(3,iat)=modulo(modulo(rxyz_red(3,iat),1.d0),1.d0)
!!      enddo
!!      call rxyz_int2cart(latvec,rxyz_red,rxyz,nat)
!!  ! v=latvec
!!  ! vol= v(1,1)*v(2,2)*v(3,3)-v(1,1)*v(2,3)*v(3,2)-v(1,2)*v(2,1)*v(3,3)+&
!!  !      v(1,2)*v(2,3)*v(3,1)+v(1,3)*v(2,1)*v(3,2)-v(1,3)*v(2,2)*v(3,1)
!!  ! if(vol.le.0.d0) stop "Negative volume during backtocell"
!!  ! call rxyz_int2cart(latvec,rxyz_red,rxyz,nat)
!!   
!!  
!!  !    ! To really be on the safe side, the translation vector can be shortened by a factor eps in order
!!  !    ! to get the atom into the cell. 
!!  !    eps = 1.d-15
!!  !    count = 0.d0
!!  !    neccesary = .true.
!!  !    do while(neccesary)
!!  !      neccesary = .false.
!!  !      count = count + 1.d0
!!  !
!!  !      ! generate 3 normal vectors of the 3 planes
!!  !      call nveclatvec(latvec,nvec)
!!  !      do iat = 1, nat
!!  !        ! 3 planes through origin (xy,yz,zx)
!!  !        do i = 1, 3
!!  !          dist(i) = DOT_PRODUCT(rxyz(:,iat),nvec(:,i))
!!  !          if(dist(i).lt.-abs(dist(i))*eps) then
!!  !            rxyz(:,iat) = rxyz(:,iat) + latvec(:,mod(i + 1, 3) + 1) !*eps
!!  !            neccesary = .true.
!!  !          endif
!!  !
!!  !          ! 3 planes on top/side/back (xy,yz,zx)
!!  !          dist(i + 3) = DOT_PRODUCT(rxyz(:,iat) - latvec(:,mod(i + 1,3) + 1),nvec(:,i))
!!  !          if(dist(i + 3).gt.abs(dist(i + 3))*eps) then
!!  !            rxyz(:,iat) = rxyz(:,iat) - latvec(:,mod(i + 1,3) + 1) !*eps
!!  !            neccesary = .true.
!!  !          endif
!!  !        enddo
!!  !      enddo
!!  !
!!  !      if(count.gt.1.d6) stop "Too many iterations in back-to-cell"
!!  !    enddo
!!   end subroutine
!!  
!!  
!!   ! This code will define the KPT mesh based on the desired grid density
!!   subroutine find_kpt(k1, k2, k3, lat, gridden)
!!     integer, intent(out) :: k1,k2,k3
!!     real(8), intent(in)  :: lat(3,3), gridden
!!     
!!     integer :: i, j
!!     real(8) :: lat1, lat2, lat3
!!     real(8) :: angles(3), cos_arr(3)
!!     real(8) :: glat(3,3), crossp(3), vol, glen(3), a(3,3)
!!     real(8) :: pi
!!     
!!     pi = acos(-1.d0)
!!     
!!     lat1 = sqrt(lat(1,1)**2 + lat(2,1)**2 + lat(3,1)**2)
!!     lat2 = sqrt(lat(1,2)**2 + lat(2,2)**2 + lat(3,2)**2)
!!     lat3 = sqrt(lat(1,3)**2 + lat(2,3)**2 + lat(3,3)**2)
!!     
!!     cos_arr(1) = (lat(1,2)*lat(1,3) + lat(2,2)*lat(2,3) + lat(3,2)*lat(3,3))/(lat2*lat3)
!!     cos_arr(2) = (lat(1,1)*lat(1,3) + lat(2,1)*lat(2,3) + lat(3,1)*lat(3,3))/(lat1*lat3)
!!     cos_arr(3) = (lat(1,1)*lat(1,2) + lat(2,1)*lat(2,2) + lat(3,1)*lat(3,2))/(lat1*lat2)
!!     
!!     angles(1) = (acos(cos_arr(1))/pi)*180.d0
!!     angles(2) = (acos(cos_arr(2))/pi)*180.d0
!!     angles(3) = (acos(cos_arr(3))/pi)*180.d0
!!     
!!     vol = cell_volume(lat)
!!     
!!     call cross_product(lat(:,2), lat(:,3), crossp(:))
!!     glat(:,1) = 2.d0*pi*crossp(:)/vol
!!     call cross_product(lat(:,3), lat(:,1), crossp(:))
!!     glat(:,2) = 2.d0*pi*crossp(:)/vol
!!     call cross_product(lat(:,1), lat(:,2), crossp(:))
!!     glat(:,3) = 2.d0*pi*crossp(:)/vol
!!     
!!     !Compute the correct kpts
!!     glen(1) = sqrt(glat(1,1)**2 + glat(2,1)**2 + glat(3,1)**2)
!!     glen(2) = sqrt(glat(1,2)**2 + glat(2,2)**2 + glat(3,2)**2)
!!     glen(3) = sqrt(glat(1,3)**2 + glat(2,3)**2 + glat(3,3)**2)
!!     
!!     call track_kpt(gridden, glen(1), k1)
!!     call track_kpt(gridden, glen(2), k2)
!!     call track_kpt(gridden, glen(3), k3)
!!     
!!   contains
!!     
!!     subroutine track_kpt(gridden, glen, kpt)
!!       real(8), intent(in) :: gridden, glen
!!       
!!       integer :: kpt,j
!!       real(8) :: d_test
!!       
!!       kpt = int(glen/(gridden*2.d0*pi))
!!       if (kpt == 0) kpt = 1
!!       d_test=glen/(kpt*2.d0*pi)
!!       if (d_test.ge.gridden) then
!!         do j = 1, 25
!!           kpt = kpt + j
!!           d_test = glen/(kpt*2.d0*pi)
!!           if (d_test.le.gridden) exit
!!         enddo
!!       endif
!!     end subroutine track_kpt
!!     
!!   end subroutine find_kpt
!!  
!!  subroutine k_expansion(latvec,xred,ka,kb,kc,k_latvec,k_xcart)
!!  !This routine expands the cell defined by latevec to a supercell
!!  !of dimension ka,kb,kc. The atomic positions in real space
!!  !will be returned in k_xcart. k_nat will then be the number of
!!  !all atoms in the supercell and is ka*kb*kc*nat
!!  use global, only: nat
!!  implicit none
!!  real(8):: latvec(3,3),k_latvec(3,3),k_xcart(3,nat,ka,kb,kc),xred(3,nat) 
!!  integer:: iat,k,l,m,ka,kb,kc
!!  do k=1,ka
!!  do l=1,kb
!!  do m=1,kc
!!  do iat=1,nat
!!     k_xcart(:,iat,k,l,m)=matmul(latvec,xred(:,iat))+&
!!     &real(k-1,8)*latvec(:,1)+real(l-1,8)*latvec(:,2)+real(m-1,8)*latvec(:,3)
!!  enddo
!!  enddo
!!  enddo
!!  enddo
!!  k_latvec(:,1)=real(ka,8)*latvec(:,1)
!!  k_latvec(:,2)=real(kb,8)*latvec(:,2)
!!  k_latvec(:,3)=real(kc,8)*latvec(:,3)
!!  end subroutine
!!  
!!  end module cell_utils
