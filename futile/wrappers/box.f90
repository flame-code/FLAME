!> @file
!!    Modulefile for handling fundamental data structed and methods of the simulation box
!!
!! @author
!!    G. Fisicaro, L. Genovese (September 2015)
!!    Copyright (C) 2016-2016 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module box

  use f_precisions, gp=>f_double

  private

  !>parameter for the definition of the bc
  integer, parameter :: FREE=0
  integer, parameter :: PERIODIC=1

  type, public :: cell
     logical :: orthorhombic !<true if the cell is orthorhombic
     integer, dimension(3) :: bc
     integer, dimension(3) :: ndims
     real(gp), dimension(3) :: hgrids
     real(gp), dimension(3) :: angrad !<angles between the dimensions in radiant (alpha_bc,beta_ac,gamma_bc)
     !derived data
     integer(f_long) :: ndim !< product of the dimension, long integer to avoid overflow
     real(gp) :: volume_element
     real(gp), dimension(3,3) :: habc !<primitive volume elements in the translation vectors direction
     real(gp), dimension(3,3) :: gd !<covariant metric needed for non-orthorhombic operations
     real(gp), dimension(3,3) :: gu !<controvariant metric needed for non-orthorhombic operations
     real(gp) :: detgd !<determinant of the covariant matrix
  end type cell

  type, public :: box_iterator
     integer :: i3s !<starting point in the dimension z
     integer :: i3e !<ending point in the dimension z
     integer :: i23 !<collapsed index in 23 dimension
     integer :: ind !<one-dimensional index for arrays
     !< 3D indices in absolute coordinates in the given box specified by boxat
     integer, dimension(3)  :: ibox  
     !> actual index inside the box
     integer :: i,j,k !better as scalars
     !> Sub-box to iterate over the points (ex. around atoms)
     !! start and end points for each direction
     integer, dimension(2,3) :: nbox 
     real(gp), dimension(3) :: oxyz !<origin of the coordinate system
     real(gp), dimension(3) :: rxyz !<coordinates of the grid point
     real(gp), dimension(3) :: tmp !< size 3 array buffer to avoid the creation of temporary arrays
     logical :: whole !<to assess if we run over the entire box or not (no check over the internal point)
     !>reference mesh from which it starts
     type(cell), pointer :: mesh
  end type box_iterator

!!$  interface box_iter
!!$     module procedure box_iter_c,box_iter_base
!!$  end interface box_iter

  interface dotp
     module procedure dotp,dotp_add1,dotp_add2
  end interface dotp

  interface square
     module procedure square,square_add
  end interface square

  public :: cell_r,cell_periodic_dims,distance,closest_r,square,cell_new,box_iter,box_next_point
  public :: cell_geocode,box_next_x,box_next_y,box_next_z,dotp,cell_null

contains

  !> Nullify the cell type
  pure function cell_null() result(me)
   implicit none
   type(cell) :: me
   me%orthorhombic=.true.
   me%bc=0
   me%ndims=0
   me%hgrids=0.0_gp
   me%angrad=0.0_gp
   !derived data
   me%ndim=0
   me%volume_element=0.0_gp
   me%habc=0.0_gp
   me%gd=0.0_gp
   me%gu=0.0_gp
   me%detgd=0.0_gp
  end function cell_null

  !> Nullify the iterator dpbox type
  pure subroutine nullify_box_iterator(boxit)
    implicit none
    type(box_iterator), intent(out) :: boxit
    boxit%i3s =-1 
    boxit%i3e =-1
    boxit%i23 =-1
    boxit%ind =-1
    boxit%i=-1
    boxit%j=-1
    boxit%k=-1
    boxit%ibox=0
    boxit%ibox(1)=-1
    boxit%nbox=-1 
    boxit%oxyz=-1.0_gp
    boxit%rxyz=-1.0_gp
    nullify(boxit%mesh)
    boxit%whole=.false.
  end subroutine nullify_box_iterator

!!$  function box_iter_c(mesh,origin) result(boxit)
!!$    type(cell), intent(in), target :: mesh
!!$    !> starting point of the box in the z direction
!!$    integer, intent(in), optional :: i3s
!!$    !> number of planes of the box to be considered
!!$    integer, intent(in), optional :: n3p
!!$    !> Box of start and end points which have to be considered
!!$    integer, dimension(2,3), intent(in), optional :: nbox
!!$    !> real coordinates of the origin in the reference frame of the 
!!$    !box (the first point has the 000 coordinate)
!!$    real(gp), dimension(3), intent(in), optional :: origin
!!$    type(box_iterator) :: boxit
!!$
!!$  end function box_iter_c

  !>define an iterator over the cell points
  function box_iter(mesh,nbox,origin,i3s,n3p,centered,cutoff) result(boxit)
    use f_utils, only: f_zero
    implicit none
    type(cell), intent(in), target :: mesh
    !>when true the origin is placed at the center of the box, origin is ignored
    logical, intent(in), optional :: centered
    !> starting point of the box in the z direction
    integer, intent(in), optional :: i3s
    !> number of planes of the box to be considered
    integer, intent(in), optional :: n3p
    !> Box of start and end points which have to be considered
    integer, dimension(2,3), intent(in), optional :: nbox
    real(gp), intent(in), optional :: cutoff !< determine the box around the origin
    !> real coordinates of the origin in the reference frame of the 
    !! box (the first point has the 000 coordinate)
    real(gp), dimension(3), intent(in), optional :: origin

    type(box_iterator) :: boxit
    
    call nullify_box_iterator(boxit)
    
    !associate the mesh
    boxit%mesh => mesh

    call f_zero(boxit%oxyz)
    if (present(origin)) boxit%oxyz=origin
    if (present(centered)) then
       if (centered) boxit%oxyz=0.5_gp*real(boxit%mesh%ndims)*boxit%mesh%hgrids
    end if

    if(present(nbox)) then
       boxit%nbox=nbox
       boxit%whole=.false.
    else if (present(cutoff)) then
       boxit%nbox(1,:)=floor((boxit%oxyz-cutoff)/boxit%mesh%hgrids)
       boxit%nbox(2,:)=ceiling((boxit%oxyz+cutoff)/boxit%mesh%hgrids)
       boxit%whole=.false.
    else
       boxit%whole=.true.
       boxit%nbox(1,:)=1
       boxit%nbox(2,:)=mesh%ndims
    end if

    if (present(i3s)) then
       boxit%i3s=i3s
    else
       boxit%i3s=1
    end if

    if (present(n3p)) then
       boxit%i3e=boxit%i3s+n3p-1
    else
       boxit%i3e=boxit%i3s+mesh%ndims(3)-1
    end if
    if (boxit%whole) boxit%whole=boxit%i3s == 1 .and. boxit%i3e==mesh%ndims(3)
    call set_starting_point(boxit)

    call probe_iterator(boxit)

  end function box_iter

  !>verify if the iterator can be used as expected
  subroutine probe_iterator(bit)
    use f_precisions
    use yaml_strings
    use dictionaries
    use dynamic_memory
    implicit none
    type(box_iterator), intent(inout) :: bit
    !local variables
    integer :: iz,iy,ix
    integer(f_long) :: icnt,itgt
    logical(f_byte), dimension(:), allocatable :: lxyz
    integer, dimension(:), allocatable :: lx,ly,lz

    !first, count if the iterator covers all the points
    itgt=int(bit%mesh%ndims(1),f_long)*int(bit%mesh%ndims(2),f_long)*&
         int(bit%i3e-bit%i3s+1,f_long)

    !allocate array of values corresponding to the expected grid
    lx=f_malloc(bit%mesh%ndims(1),id='lx')
    ly=f_malloc(bit%mesh%ndims(2),id='ly')
    lz=f_malloc(bit%i3e-bit%i3s+1,id='lz')

    lxyz=f_malloc0(itgt,id='lxyz')

    do iz=bit%i3s,bit%i3e
       lz(iz-bit%i3s+1)=iz
    end do
    do iy=1,bit%mesh%ndims(2)
       ly(iy)=iy
    end do
    do ix=1,bit%mesh%ndims(1)
       lx(ix)=ix
    end do

    !separable mode
    iz=0
    icnt=0
    do while(box_next_z(bit))
       iz=iz+1
       iy=0
       do while(box_next_y(bit))
          iy=iy+1
          ix=0
          do while(box_next_x(bit))
             ix=ix+1
             icnt=icnt+1
             if (lx(bit%i) /= bit%i) &
                  call f_err_throw('Error value, ix='+bit%i+&
                  ', expected='+lx(bit%i))
             !convert the value of the logical array
             if (lxyz(bit%ind)) call f_err_throw('Error point ind='+bit%ind+&
               ', i,j,k='+yaml_toa([bit%i,bit%j,bit%k]))
             lxyz(bit%ind)=f_T
          end do
          if (ix /= bit%mesh%ndims(1)) &
               call f_err_throw('Error boxit, ix='+ix+&
               ', itgtx='+bit%mesh%ndims(1))
          if (ly(bit%j) /= bit%j) &
               call f_err_throw('Error value, iy='+bit%j+&
               ', expected='+ly(bit%j))
       end do
       if (iy /= bit%mesh%ndims(2)) &
            call f_err_throw('Error boxit, iy='+iy+&
            ', itgty='+bit%mesh%ndims(2))
       if (lz(bit%k-bit%i3s+1) /= bit%k) &
            call f_err_throw('Error value, iz='+bit%k+&
            ', expected='+lz(bit%k-bit%i3s+1))
    end do
    if (iz /= bit%i3e-bit%i3s+1) call f_err_throw('Error boxit, iz='+iz+&
         ', itgtz='+(bit%i3e-bit%i3s+1))
    if (icnt /= itgt) call f_err_throw('Error sep boxit, icnt='+icnt+&
         ', itgt='+itgt)

    !complete mode
    icnt=int(0,f_long)
    do while(box_next_point(bit))
       icnt=icnt+1
       !here we might see if there are points from which 
       !we passed twice
       !print *,bit%i,bit%j,bit%k
       if (.not. lxyz(bit%ind)) &
            call f_err_throw('Error point (2) ind='+bit%ind+&
            ', i,j,k='+yaml_toa([bit%i,bit%j,bit%k]))
       lxyz(bit%ind)=f_F
    end do
    if (icnt /= itgt) call f_err_throw('Error boxit, icnt='+icnt+&
         ', itgt='+itgt)

    if (any(lxyz)) call f_err_throw('Error boxit, points not covered')

    call f_free(lxyz)
    call f_free(lx,ly,lz)
    
  end subroutine probe_iterator

  pure subroutine set_starting_point(bit)
    implicit none
    type(box_iterator), intent(inout) :: bit
    !local variables
    logical :: test

    bit%ibox=bit%nbox(1,:)
!!$    if (bit%whole) then
!!$       bit%i=1
!!$       bit%j=1
!!$       bit%k=1
!!$    else
!!$       bit%i=-1
!!$       bit%j=-1
!!$       bit%k=-1
!!$    end if
    bit%k=bit%nbox(1,3)-1
    bit%ind=0
    bit%i23=0
  end subroutine set_starting_point

  !find the first z value which is available from the starting point
  function box_next_z(bit) result(ok)
    implicit none
    type(box_iterator), intent(inout) :: bit
    logical :: ok

    ok = bit%i3e >= bit%i3s
    if (.not. ok) return !there is nothing to explore
    ok= bit%ibox(3) <= bit%nbox(2,3)
    do while(ok)
       if (bit%whole) then
          bit%k=bit%ibox(3)
       else 
          call internal_point(bit%mesh%bc(3),bit%ibox(3),bit%mesh%ndims(3),&
               bit%k,bit%i3s,bit%i3e,ok)
          if (.not. ok) bit%ibox(3)=bit%ibox(3)+1
       end if
       if (ok) then
          bit%ibox(3)=bit%ibox(3)+1
          exit
       end if
       ok = bit%ibox(3) <= bit%nbox(2,3)
    end do
    !reset x and y
    if (ok) then
       call update_boxit_z(bit)
       bit%ibox(2)=bit%nbox(1,2)
       bit%ibox(1)=bit%nbox(1,1)
    end if

    !in the case the z_direction is over, make the iterator ready for new use
    if (.not. ok) call set_starting_point(bit)
    
  end function box_next_z

  !find the first z value which is available from the starting point
  function box_next_y(bit) result(ok)
    implicit none
    type(box_iterator), intent(inout) :: bit
    logical :: ok

    ok= bit%ibox(2) <= bit%nbox(2,2)
    do while(ok)
       if (bit%whole) then
          bit%j=bit%ibox(2)
       else 
          call internal_point(bit%mesh%bc(2),bit%ibox(2),bit%mesh%ndims(2),&
               bit%j,1,bit%mesh%ndims(2),ok)
          if (.not. ok) bit%ibox(2)=bit%ibox(2)+1
       end if
       if (ok) then
          bit%ibox(2)=bit%ibox(2)+1
          exit
       end if
       ok = bit%ibox(2) <= bit%nbox(2,2)
    end do
    !reset x
    if (ok) then
       call update_boxit_y(bit)
       bit%ibox(1)=bit%nbox(1,1)
    end if


  end function box_next_y

  !find the first z value which is available from the starting point
  function box_next_x(bit) result(ok)
    implicit none
    type(box_iterator), intent(inout) :: bit
    logical :: ok

    ok= bit%ibox(1) <= bit%nbox(2,1)
    do while(ok)
       if (bit%whole) then
          bit%i=bit%ibox(1)
       else 
          call internal_point(bit%mesh%bc(1),bit%ibox(1),bit%mesh%ndims(1),&
               bit%i,1,bit%mesh%ndims(1),ok)
          if (.not. ok) bit%ibox(1)=bit%ibox(1)+1
       end if
       if (ok) then
          !increment for after
          bit%ibox(1)=bit%ibox(1)+1
          exit
       end if
       ok = bit%ibox(1) <= bit%nbox(2,1)
    end do
    if (ok) call update_boxit_x(bit)
  end function box_next_x

  pure subroutine update_boxit_x(boxit)
    implicit none
    type(box_iterator), intent(inout) :: boxit 

    !one dimensional index
    boxit%ind = boxit%i+boxit%mesh%ndims(1)*boxit%i23

    !the position associated to the coordinates
    boxit%rxyz(1)=cell_r(boxit%mesh,boxit%i,1)-boxit%oxyz(1)

  end subroutine update_boxit_x

  pure subroutine update_boxit_y(boxit)
    implicit none
    type(box_iterator), intent(inout) :: boxit 
    !here we have the indices      boxit%ibox as well as boxit%ixyz
    !we might then calculate the related quantities
    !two dimensional index, last two elements
    boxit%i23=(boxit%j-1)+&
         boxit%mesh%ndims(2)*(boxit%k-boxit%i3s)

    !the position associated to the coordinates
    boxit%rxyz(2)=cell_r(boxit%mesh,boxit%j,2)-boxit%oxyz(2)

  end subroutine update_boxit_y

  pure subroutine update_boxit_z(boxit)
    implicit none
    type(box_iterator), intent(inout) :: boxit 

    !the position associated to the coordinates
    boxit%rxyz(3)=cell_r(boxit%mesh,boxit%k,3)-boxit%oxyz(3)

  end subroutine update_boxit_z

  !this routine should not use ibox as it is now prepared for the next step
  pure subroutine update_boxit(boxit)
    implicit none
    type(box_iterator), intent(inout) :: boxit 

    !here we have the indices      boxit%ibox as well as boxit%ixyz
    !we might then calculate the related quantities
    !two dimensional index, last two elements
    boxit%i23=(boxit%j-1)+&
         boxit%mesh%ndims(2)*(boxit%k-boxit%i3s)
    !one dimensional index
    boxit%ind = boxit%i+boxit%mesh%ndims(1)*boxit%i23

    !the position associated to the coordinates
    boxit%rxyz(1)=cell_r(boxit%mesh,boxit%i,1)-boxit%oxyz(1)
    boxit%rxyz(2)=cell_r(boxit%mesh,boxit%j,2)-boxit%oxyz(2)
    boxit%rxyz(3)=cell_r(boxit%mesh,boxit%k,3)-boxit%oxyz(3)
  end subroutine update_boxit

  function box_next_point(boxit)
    implicit none
    type(box_iterator), intent(inout) :: boxit
    logical :: box_next_point
    !local variables
    logical :: go

    box_next_point=associated(boxit%mesh)
    if (.not. box_next_point) return

    !this put the starting point
    if (boxit%k==boxit%nbox(1,3)-1) then
       go=box_next_z(boxit)
       if (go) go=box_next_y(boxit)
       !if this one fails then there are no slices available
       if (.not. go) then
          box_next_point=.false.
       end if
    end if
    !simulate loop
    flattened_loop: do 
       if (box_next_x(boxit)) exit flattened_loop
       if (box_next_y(boxit)) cycle flattened_loop !and then redo the check for x
       box_next_point =box_next_z(boxit)
       if (box_next_point) box_next_point =box_next_y(boxit)
       if (.not. box_next_point) exit flattened_loop
    end do flattened_loop
    
  end function box_next_point

  pure subroutine internal_point(bc,ipoint,npoint,jpoint,ilow,ihigh,go)
    implicit none
    integer, intent(in) :: bc
    integer, intent(in) :: npoint,ilow,ihigh,ipoint
    logical, intent(out) :: go
    integer, intent(out) :: jpoint

    if (bc == PERIODIC) then
       jpoint=modulo(ipoint-1,npoint)+1
    else
       jpoint=ipoint
    end if
    go=jpoint >= ilow
    if (go) go= jpoint <= ihigh

  end subroutine internal_point

  function cell_new(geocode,ndims,hgrids,alpha_bc,beta_ac,gamma_ab,abc) result(mesh)
    use numerics, only: onehalf,pi
    use wrapper_linalg, only: det_3x3
    use f_utils, only: f_assert
    use dictionaries, only: f_err_throw
    implicit none
    character(len=1), intent(in) :: geocode
    integer, dimension(3), intent(in) :: ndims
    real(gp), dimension(3), intent(in) :: hgrids
    !real(gp), dimension(3), intent(in), optional :: angrad
    real(gp), intent(in), optional :: alpha_bc,beta_ac,gamma_ab
    !> arrays of the unit vectors of the cell. Normalized, in fortran order a_i=abc(i,1), b_i=abc(i,2)
    real(gp), dimension(3,3), intent(in), optional :: abc
    type(cell) :: mesh
    !local variables
    real(gp) :: aa,cc,a2,cosang
    integer :: i,j


    select case(geocode)
    case('P')
       mesh%bc=PERIODIC
    case('S')
       mesh%bc=PERIODIC
       mesh%bc(2)=FREE
    case('F')
       mesh%bc=FREE
    case('W')
       mesh%bc=FREE
       mesh%bc(3)=PERIODIC
    case default
       call f_err_throw('Invalid specification of the variable "geocode"')
    end select
    mesh%ndims=ndims
    mesh%hgrids=hgrids
    mesh%ndim=product(int(ndims,f_long))

    !default orthorhombic
    mesh%angrad=onehalf*pi

    if (present(alpha_bc)) mesh%angrad(1)=alpha_bc
    if (present(beta_ac)) mesh%angrad(2)=beta_ac
    if (present(gamma_ab)) mesh%angrad(3)=gamma_ab

    call f_assert(all(mesh%angrad > 0.0_gp),'Error, Cell new, some of the angles are not positive')

    if (geocode == 'S') then
       call f_assert(mesh%angrad(1)-onehalf*pi,id='Alpha angle invalid')
       call f_assert(mesh%angrad(3)-onehalf*pi,id='Gamma angle invalid')
    end if

    mesh%orthorhombic=all(mesh%angrad==onehalf*pi)

    if ((geocode == 'F' .or. geocode== 'W') .and. (.not. mesh%orthorhombic)) &
         call f_err_throw('For geocode="F","W" the cell must be orthorhombic')

    if (.not. mesh%orthorhombic) then
       !some consistency check on the angles should be performed
       !1) sum(angrad) < twopi
       if (all(mesh%angrad==mesh%angrad(1))) then
          !Treat the case of equal angles (except all right angles) :
          !generates trigonal symmetry wrt third axis
          cosang=cos(mesh%angrad(1))
          a2=2.0_gp/3.0_gp*(1.0_gp-cosang)
          aa=sqrt(a2)
          cc=sqrt(1.0_gp-a2)
          mesh%habc(1,1)=aa; mesh%habc(2,1)=0.0_gp; mesh%habc(3,1)=cc
          mesh%habc(1,2)=-0.5_gp*aa ; mesh%habc(2,2)=sqrt(3.0_gp)*0.5_gp*aa ; mesh%habc(3,2)=cc
          mesh%habc(1,3)=-0.5_gp*aa ; mesh%habc(2,3)=-sqrt(3.0_gp)*0.5_gp*aa ; mesh%habc(3,3)=cc
          !Set the covariant metric
          mesh%gd(1,1) = 1.0_gp
          mesh%gd(1,2) = cos(mesh%angrad(3)) !gamma_ab
          mesh%gd(1,3) = cos(mesh%angrad(2)) !beta_ac
          mesh%gd(2,2) = 1.0_gp
          mesh%gd(2,3) = cos(mesh%angrad(1)) !alpha_bc
          mesh%gd(3,3) = 1.0_gp 
          !Set the determinant of the covariant metric
          mesh%detgd = 1.0_gp - cos(mesh%angrad(1))**2 - cos(mesh%angrad(2))**2 - cos(mesh%angrad(3))**2 +&
               2.0_gp*cos(mesh%angrad(1))*cos(mesh%angrad(2))*cos(mesh%angrad(3))
          !Set the contravariant metric
          mesh%gu(1,1) = (sin(mesh%angrad(1))**2)/mesh%detgd
          mesh%gu(1,2) = (cos(mesh%angrad(2))*cos(mesh%angrad(1))-cos(mesh%angrad(3)))/mesh%detgd
          mesh%gu(1,3) = (cos(mesh%angrad(3))*cos(mesh%angrad(1))-cos(mesh%angrad(2)))/mesh%detgd
          mesh%gu(2,2) = (sin(mesh%angrad(2))**2)/mesh%detgd
          mesh%gu(2,3) = (cos(mesh%angrad(3))*cos(mesh%angrad(2))-cos(mesh%angrad(1)))/mesh%detgd
          mesh%gu(3,3) = (sin(mesh%angrad(3))**2)/mesh%detgd
       else if (geocode == 'P') then
          mesh%habc=0.0_gp
          mesh%habc(1,1)=1.0_gp
          mesh%habc(1,2)=cos(mesh%angrad(3))
          mesh%habc(2,2)=sin(mesh%angrad(3))
          mesh%habc(1,3)=cos(mesh%angrad(2))
          mesh%habc(2,3)=(cos(mesh%angrad(1))-mesh%habc(1,2)*mesh%habc(1,3))/mesh%habc(2,2)
          mesh%habc(3,3)=sqrt(1.0_gp-mesh%habc(1,3)**2-mesh%habc(2,3)**2)
          !Set the covariant metric
          mesh%gd(1,1) = 1.0_gp
          mesh%gd(1,2) = cos(mesh%angrad(3)) !gamma_ab
          mesh%gd(1,3) = cos(mesh%angrad(2)) !beta_ac
          mesh%gd(2,2) = 1.0_gp
          mesh%gd(2,3) = cos(mesh%angrad(1)) !alpha_bc
          mesh%gd(3,3) = 1.0_gp 
          !Set the determinant of the covariant metric
          mesh%detgd = 1.0_gp - cos(mesh%angrad(1))**2 - cos(mesh%angrad(2))**2 - cos(mesh%angrad(3))**2 +&
               2.0_gp*cos(mesh%angrad(1))*cos(mesh%angrad(2))*cos(mesh%angrad(3))
          !Set the contravariant metric
          mesh%gu(1,1) = (sin(mesh%angrad(1))**2)/mesh%detgd
          mesh%gu(1,2) = (cos(mesh%angrad(2))*cos(mesh%angrad(1))-cos(mesh%angrad(3)))/mesh%detgd
          mesh%gu(1,3) = (cos(mesh%angrad(3))*cos(mesh%angrad(1))-cos(mesh%angrad(2)))/mesh%detgd
          mesh%gu(2,2) = (sin(mesh%angrad(2))**2)/mesh%detgd
          mesh%gu(2,3) = (cos(mesh%angrad(3))*cos(mesh%angrad(2))-cos(mesh%angrad(1)))/mesh%detgd
          mesh%gu(3,3) = (sin(mesh%angrad(3))**2)/mesh%detgd
       else !only Surfaces is possible here
          mesh%habc=0.0_gp
          mesh%habc(1,1)=1.0_gp
          mesh%habc(2,2)=1.0_gp
          mesh%habc(1,3)=cos(mesh%angrad(2))
          mesh%habc(3,3)=sin(mesh%angrad(2))
          !Set the covariant metric
          mesh%gd=0.0_gp
          mesh%gd(1,1) = 1.0_gp
          mesh%gd(1,3) = cos(mesh%angrad(2)) !beta_ac
          mesh%gd(2,2) = 1.0_gp
          mesh%gd(3,3) = 1.0_gp 
          !Set the determinant of the covariant metric
          mesh%detgd = sin(mesh%angrad(2))**2 
          !Set the contravariant metric
          mesh%gu=0.0_gp
          mesh%gu(1,1) = 1.0_gp/mesh%detgd
          mesh%gu(1,3) = -cos(mesh%angrad(2))/mesh%detgd
          mesh%gu(2,2) = 1.0_gp!/mesh%detgd
          mesh%gu(3,3) = 1.0_gp/mesh%detgd
       end if
       !Rescale habc using hgrid
       mesh%habc(:,1)=hgrids*mesh%habc(:,1)
       mesh%habc(:,2)=hgrids*mesh%habc(:,2)
       mesh%habc(:,3)=hgrids*mesh%habc(:,3)
       !the volume element
       !Compute unit cell volume
       mesh%volume_element=det_3x3(mesh%habc)
    else
       mesh%habc=0.0_gp
       do i=1,3
          mesh%habc(i,i)=hgrids(i)
       end do
       mesh%angrad=onehalf*pi
       mesh%volume_element=product(mesh%hgrids)
       mesh%gd(1,1) = 1.0_gp
       mesh%gd(1,2) = 0.0_gp
       mesh%gd(1,3) = 0.0_gp
       mesh%gd(2,2) = 1.0_gp
       mesh%gd(2,3) = 0.0_gp
       mesh%gd(3,3) = 1.0_gp 
       mesh%detgd = 1.0_gp 
       !Set the contravariant metric
       mesh%gu(1,1) = 1.0_gp
       mesh%gu(1,2) = 0.0_gp
       mesh%gu(1,3) = 0.0_gp
       mesh%gu(2,2) = 1.0_gp
       mesh%gu(2,3) = 0.0_gp 
       mesh%gu(3,3) = 1.0_gp
    end if
    mesh%gd(2,1) = mesh%gd(1,2)
    mesh%gd(3,1) = mesh%gd(1,3)
    mesh%gd(3,2) = mesh%gd(2,3)

    mesh%gu(2,1) = mesh%gu(1,2)
    mesh%gu(3,1) = mesh%gu(1,3)
    mesh%gu(3,2) = mesh%gu(2,3)
    do i=1,3
       do j=1,3
          if (abs(mesh%gd(i,j)).lt.1.0d-15) mesh%gd(i,j)=0.0_gp
          if (abs(mesh%gu(i,j)).lt.1.0d-15) mesh%gu(i,j)=0.0_gp
       end do
    end do

    !here we should verify that the the inverse metric times the metric is the identity



  end function cell_new

  !> returns a logical array of size 3 which is .true. for all the periodic dimensions
  pure function cell_periodic_dims(mesh) result(peri)
    implicit none
    type(cell), intent(in) :: mesh
    logical, dimension(3) :: peri
    !local variables

    peri= mesh%bc == PERIODIC

  end function cell_periodic_dims

  !>give the associated geocode, 'X' for unknown
  pure function cell_geocode(mesh)
    implicit none
    type(cell), intent(in) :: mesh
    character(len=1) :: cell_geocode
    !local variables
    logical, dimension(3) :: peri

    peri=cell_periodic_dims(mesh)
    if (all(peri)) then
       cell_geocode='P'
    else if (.not. any(peri)) then
       cell_geocode='F'
    else if (peri(1) .and. .not. peri(2) .and. peri(3)) then
       cell_geocode='S'
    else if (.not. peri(1) .and. .not. peri(2) .and. peri(3)) then
       cell_geocode='W'
    else
       cell_geocode='X'
    end if

  end function cell_geocode


  !>gives the value of the coordinate from the grid point
  elemental pure function cell_r(mesh,i,dim) result(t)
    implicit none
    integer, intent(in) :: i
    type(cell), intent(in) :: mesh
    integer, intent(in) :: dim
    real(gp) :: t

    t=mesh%hgrids(dim)*(i-1)
  end function cell_r

  function distance(mesh,v1,v2) result(d)
    use dictionaries, only: f_err_throw
    implicit none
    real(gp), dimension(3), intent(in) :: v1,v2
    type(cell), intent(in) :: mesh
    real(gp) :: d
    !local variables
    integer :: i
    real(gp) :: d2

    if (mesh%orthorhombic) then
       d2=0.0_gp
       do i=1,3
          d2=d2+r_wrap(mesh%bc(i),mesh%hgrids(i)*mesh%ndims(i),&
               v1(i),v2(i))**2
       end do
       d=sqrt(d2)
    else
       call f_err_throw('Distance not yet implemented for nonorthorhombic cells')
    end if

  end function distance

!!$  pure function min_dist(bc,alat,r,r_old)
!!$    implicit none
!!$    integer, intent(in) :: bc
!!$    real(gp), intent(in) :: r,r_old,alat
!!$    real(gp) :: min_dist
!!$
!!$    !for periodic BC calculate mindist only if the center of mass can be defined without the modulo
!!$    min_dist=abs(r-r_old)
!!$    if (bc==PERIODIC) then
!!$       if (min_dist > 0.5_gp*alat) then
!!$          if (r < 0.5_gp*alat) then
!!$             min_dist=abs(r+alat-r_old)
!!$          else
!!$             min_dist=abs(r-alat-r_old)
!!$          end if
!!$       end if
!!$    end if
!!$
!!$  end function min_dist

  !> Calculates the minimum difference between two coordinates
  pure function r_wrap(bc,alat,r,c)
    implicit none
    integer, intent(in) :: bc
    real(gp), intent(in) :: r,c,alat
    real(gp) :: r_wrap

    !for periodic BC calculate mindist only if the center of mass can be defined without the modulo
    r_wrap=r-c
    if (bc==PERIODIC) then
       if (abs(r_wrap) > 0.5_gp*alat) then
          if (r < 0.5_gp*alat) then
             r_wrap=r+alat-c
          else
             r_wrap=r-alat-c
          end if
       end if
    end if

  end function r_wrap

  !>find the closest center according to the periodiciy of the
  !! box and provide the vector
  pure function closest_r(mesh,v,center) result(r)
    implicit none
    real(gp), dimension(3), intent(in) :: v,center
    type(cell), intent(in) :: mesh
    real(gp), dimension(3) :: r
    !local variables
    integer :: i

    if (mesh%orthorhombic) then
       do i=1,3
          r(i)=r_wrap(mesh%bc(i),mesh%hgrids(i)*mesh%ndims(i),&
               v(i),center(i))
       end do
    end if

  end function closest_r

  !> Calculates the square of the vector r in the cell defined by mesh
  !! Takes into account the non-orthorhombicity of the box
  pure function square(mesh,v)
    implicit none
    !> array of coordinate in the mesh reference frame
    real(gp), dimension(3), intent(in) :: v
    type(cell), intent(in) :: mesh !<definition of the cell
    real(gp) :: square
    integer :: i,j

    if (mesh%orthorhombic) then
       square=v(1)**2+v(2)**2+v(3)**2
    else
       square=dotp(mesh,v,v)
    end if

  end function square

  function square_add(mesh,v_add) result(square)
    implicit none
    !> array of coordinate in the mesh reference frame
    real(gp) :: v_add
    type(cell), intent(in) :: mesh !<definition of the cell
    real(gp) :: square

    if (mesh%orthorhombic) then
       call dotp_external_ortho(v_add,v_add,square)
    end if

  end function square_add


  pure function dotp(mesh,v1,v2)
    implicit none
    real(gp), dimension(3), intent(in) :: v1,v2
    type(cell), intent(in) :: mesh !<definition of the cell
    real(gp) :: dotp
    !local variables
    integer :: i,j

    if (mesh%orthorhombic) then
       dotp=v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)
    else
       dotp=0.0_gp
       do i=1,3
          do j=1,3
             dotp=dotp+mesh%gu(i,j)*v1(i)*v2(j)
          end do
       end do
    end if

  end function dotp

  function dotp_add2(mesh,v1,v2_add) result(dotp)
    implicit none
    real(gp), dimension(3), intent(in) :: v1
    real(gp) :: v2_add !<intent in, cannot be declared as such
    type(cell), intent(in) :: mesh !<definition of the cell
    real(gp) :: dotp

    if (mesh%orthorhombic) then
       call dotp_external_ortho(v1,v2_add,dotp)
    end if

  end function dotp_add2

  function dotp_add1(mesh,v1_add,v2) result(dotp)
    implicit none
    real(gp), dimension(3), intent(in) :: v2
    real(gp) :: v1_add !<intent in, cannot be declared as such
    type(cell), intent(in) :: mesh !<definition of the cell
    real(gp) :: dotp

    if (mesh%orthorhombic) then
       call dotp_external_ortho(v1_add,v2,dotp)
    end if

  end function dotp_add1

end module box

subroutine dotp_external_ortho(v1,v2,dotp)  
  use f_precisions, only: gp=>f_double
  implicit none
  real(gp), dimension(3), intent(in) :: v1,v2
  real(gp), intent(out) :: dotp

  dotp=v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)
end subroutine dotp_external_ortho
