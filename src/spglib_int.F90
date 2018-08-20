!program convert
!implicit none
!real(8), allocatable:: rxyz(:,:),xred(:,:)
!real(8):: dproj(6),acell(3),rprim(3,3),latvec(3,3),angdeg(6)
!integer:: nat, iat,count,jat,ntype,spg
!real(8):: angbohr,pressure_gpa,pressconv,tol
!parameter(angbohr=1.889725989d0)
!parameter(pressconv=29421.033d0)
!character(3), allocatable:: typat(:)
!character(3), allocatable:: typatt(:)
!integer, allocatable:: kinds(:),nkindsat(:)
!character(40):: filename
!character(40):: type_format,formatt
!logical::new
!write(*,*) "Input file:"
!read(*,*) filename
!write(*,*) "Tolerance"
!read(*,*) tol
!open(unit=12,file=trim(filename))
!read(12,*) nat
!allocate(kinds(nat),rxyz(3,nat),xred(3,nat),typat(nat))
!read(12,*) dproj(1:3)
!read(12,*) dproj(4:6)
!
!do iat=1,nat
!read(12,*) rxyz(:,iat),typat(iat)
!enddo
!close(12)
!
!call dproj2latvec(dproj,latvec)
!!write(*,*) "Latvec"
!!write(*,'(3(1x,f15.10))') latvec(:,1)
!!write(*,'(3(1x,f15.10))') latvec(:,2)
!!write(*,'(3(1x,f15.10))') latvec(:,3)
!
!call backtocell(nat,latvec,rxyz)
!
!
!call latvec2dist_ang(latvec,angdeg)
!!write(*,*) "a,b,c,alpha,beta,gamma"
!!write(*,'(6(1x,es25.16))') angdeg
!!
!!write(type_format,'(a,i0,a)') "(",nat,"(1x,i0))"
!
!
!kinds(1)=1
!ntype=1
!do iat=2,nat
!new=.true.
! do jat=1,iat-1
! if(trim(typat(iat))==trim(typat(jat))) then
!     new=.false.
!     kinds(iat)=kinds(jat)
! endif
! enddo
! if(new) then
! ntype=ntype+1
! kinds(iat)=ntype
! endif
!enddo
!allocate(nkindsat(ntype),typatt(ntype))
!nkindsat=0
!
!do iat=1,ntype
!new=.false.
!do jat=1,nat
!  if(kinds(jat)==iat) then
!      nkindsat(iat)=nkindsat(iat)+1
!      typatt(iat)=trim(typat(jat))
!      new=.true.
!  endif
!  if(new.and.kinds(jat).lt.iat) stop "Atoms must be ordered!!!"
!enddo
!enddo
!!write(*,*) nkindsat
!
!
!
!!write(*,*) "Atom kinds"
!!write(*,type_format) kinds(:)
!!
!!
!!write(*,*) "Cartesian Coordinates"
!!do iat=1,nat
!!write(*,'(3(1x,es25.16),a4)') rxyz(:,iat),typat(iat)
!!enddo 
!
!call rxyz_cart2int(latvec,xred,rxyz,nat)
!
!call get_spg(nat,xred,latvec,kinds,tol,spg)
!!write(*,*) "Reduced Coordinates"
!!do iat=1,nat
!!write(*,'(3(1x,f15.10))') xred(:,iat)
!!enddo 
!end program convert
!
subroutine get_spg(num_atom,positions,lattice,atom_types,symprec,spg)
  use spglib_f08
  use yaml_output
implicit none
integer:: nat, typat(num_atom), spg
  ! Arguments ------------------------------------
  ! scalars
  integer, intent(in) :: num_atom!, max_num_sym, is_time_reversal
  real(8), intent(in) :: symprec
  ! arrays
  integer, intent(in), dimension(num_atom) :: atom_types
!  integer, intent(in), dimension(3) :: mesh, is_shift
  real(8), intent(in), dimension(3, 3) :: lattice
  real(8), dimension(3,3):: lattrans
  real(8), intent(in), dimension(3, num_atom) :: positions
  ! Local variables-------------------------------
  ! scalars
  integer :: i, j, counter, weight, space_group, num_sym, indent
  integer :: num_ir_grid
  character(len=21) :: international
  character(len=10) :: schoenflies
  character(len=30) :: space
  character(len=128) :: fmt
  ! arrays
!  integer, dimension(3, 3, max_num_sym) :: rotations
!  integer, dimension(3, mesh(1)*mesh(2)*mesh(3)) :: grid_point
!  integer, dimension(mesh(1)*mesh(2)*mesh(3)) :: map
!  real(8), dimension(3, max_num_sym) :: translations

  type(SpglibDataset) :: dset
!   The allocatable components of dset get deallocated on scope exit
  lattrans=transpose(lattice)
  dset = spg_get_dataset(lattrans, positions, atom_types, num_atom, symprec)
  
  num_sym = dset % n_operations
  
  indent = 1
  if (dset % spacegroup_number /= 0) then
     schoenflies = ' '
     space_group = spg_get_schoenflies( schoenflies, & 
          & lattrans, positions, atom_types, num_atom, symprec );

     !write(*,'(a,i5)') " # SPGLIB: space_group ", dset % spacegroup_number
     call yaml_map('space_group',dset%spacegroup_number,fmt='(i5)')
  !   write(*,'(a, a)') " # SPGLIB: international ", trim(dset % international_symbol)
  !   write(*,'(a, a)') " # SPGLIB: schoenflies ", trim(schoenflies)
  else
     write(*,'(a)') " # SPGLIB: Space group could not be found"
  end if
  spg=dset % spacegroup_number  

end subroutine


!************************************************************************************

subroutine spg_cell_refine(nat_in,nat_out,nat_max,positions,lattice,atom_types,symprec,spg)
  use spglib_f08
implicit none
integer:: nat, spg
  ! Arguments ------------------------------------
  ! scalars
  integer:: nat_in,nat_max!, max_num_sym, is_time_reversal
  integer:: nat_out
  real(8), intent(in) :: symprec
  ! arrays
  integer, intent(inout), dimension(nat_max) :: atom_types
!  integer, intent(in), dimension(3) :: mesh, is_shift
  real(8), intent(inout), dimension(3, 3) :: lattice
  real(8), dimension(3,3):: lattrans
  real(8), intent(inout), dimension(3, nat_max) :: positions
  ! Local variables-------------------------------
  ! scalars
  integer :: i, j, counter, weight, space_group, num_sym, indent
  integer :: num_ir_grid
  character(len=21) :: international
  character(len=10) :: schoenflies
  character(len=30) :: space
  character(len=128) :: fmt
  type(SpglibDataset) :: dset
!   The allocatable components of dset get deallocated on scope exit
  lattrans=transpose(lattice)
  dset = spg_get_dataset(lattrans, positions, atom_types, nat_in, symprec)
  
  num_sym = dset % n_operations
  
  indent = 1
  if (dset % spacegroup_number /= 0) then
     schoenflies = ' '
     space_group = spg_get_schoenflies( schoenflies, & 
          & lattrans, positions, atom_types, nat_in, symprec );

     write(*,'(a,i5)') " # SPGLIB: space_group ", dset % spacegroup_number
  !   write(*,'(a, a)') " # SPGLIB: international ", trim(dset % international_symbol)
  !   write(*,'(a, a)') " # SPGLIB: schoenflies ", trim(schoenflies)
  else
     write(*,'(a)') " # SPGLIB: Space group could not be found"
  end if
  spg=dset % spacegroup_number  

  nat_out=spg_refine_cell( lattrans, positions, atom_types, nat_in, symprec)
  lattice=transpose(lattrans)
end subroutine

!************************************************************************************

subroutine spg_cell_primitive(nat_in,nat_out,nat_max,positions,lattice,atom_types,symprec,spg)
  use spglib_f08
implicit none
integer:: nat, spg
  ! Arguments ------------------------------------
  ! scalars
  integer:: nat_in,nat_max!, max_num_sym, is_time_reversal
  integer:: nat_out
  real(8), intent(in) :: symprec
  ! arrays
  integer, intent(inout), dimension(nat_max) :: atom_types
!  integer, intent(in), dimension(3) :: mesh, is_shift
  real(8), intent(inout), dimension(3, 3) :: lattice
  real(8), dimension(3,3):: lattrans
  real(8), intent(inout), dimension(3, nat_max) :: positions
  ! Local variables-------------------------------
  ! scalars
  integer :: i, j, counter, weight, space_group, num_sym, indent
  integer :: num_ir_grid
  character(len=21) :: international
  character(len=10) :: schoenflies
  character(len=30) :: space
  character(len=128) :: fmt
  type(SpglibDataset) :: dset
!   The allocatable components of dset get deallocated on scope exit
  lattrans=transpose(lattice)
  dset = spg_get_dataset(lattrans, positions, atom_types, nat_in, symprec)
  
  num_sym = dset % n_operations
  
  indent = 1
  if (dset % spacegroup_number /= 0) then
     schoenflies = ' '
     space_group = spg_get_schoenflies( schoenflies, & 
          & lattrans, positions, atom_types, nat_in, symprec );

     write(*,'(a,i5)') " # SPGLIB: space_group ", dset % spacegroup_number
  !   write(*,'(a, a)') " # SPGLIB: international ", trim(dset % international_symbol)
  !   write(*,'(a, a)') " # SPGLIB: schoenflies ", trim(schoenflies)
  else
     write(*,'(a)') " # SPGLIB: Space group could not be found"
  end if
  spg=dset % spacegroup_number  

  nat_out=spg_find_primitive( lattrans, positions, atom_types, nat_in, symprec)
  lattice=transpose(lattrans)
end subroutine


!!************************************************************************************

!!************************************************************************************
!
! subroutine rxyz_cart2int(latvec,rxyzint,rxyzcart,nat)
! !This subrouine will convert the internal coordinates into cartesian coordinates
! implicit none
! real(8):: rxyzint(3,nat), rxyzcart(3,nat),latvec(3,3),latvecinv(3,3)
! integer:: nat,iat
! call invertmat(latvec,latvecinv,3)
! do iat=1,nat
!  rxyzint(:,iat)=matmul(latvecinv,rxyzcart(:,iat))
! enddo
! end subroutine rxyz_cart2int
!
!
!!************************************************************************************
!
! subroutine dproj2latvec(dproj,latvec)
! !This subroutine will convert the distance and projective representation of 
! !a periodic cell (dxx,dyx,dyy,dzx,dzy,dzz) into a 
! !lattice vektor format (vec1(:,1),vec2(:,2),vec3(:,3)) with dxx oriented into x direction
! implicit none
! real*8:: dproj(6),latvec(3,3)
!
! latvec(:,:)=0.d0
! latvec(1,1)=dproj(1)
! latvec(1,2)=dproj(2)
! latvec(2,2)=dproj(3)
! latvec(1,3)=dproj(4)
! latvec(2,3)=dproj(5)
! latvec(3,3)=dproj(6)
! return
! end subroutine
!
!!************************************************************************************
!!********************************************************************************
!
!subroutine latvec2dist_ang(latvec,dist_ang)
!implicit none
!real(8):: latvec(3,3),dist_ang(6),pi,dist_ang_tmp(6)
!integer:: i
!pi=acos(-1.d0)
!dist_ang(1)=sqrt(latvec(1,1)**2+latvec(2,1)**2+latvec(3,1)**2)
!dist_ang(2)=sqrt(latvec(1,2)**2+latvec(2,2)**2+latvec(3,2)**2)
!dist_ang(3)=sqrt(latvec(1,3)**2+latvec(2,3)**2+latvec(3,3)**2)
!dist_ang_tmp(4)=dot_product(latvec(:,2),latvec(:,3))/dist_ang(2)/dist_ang(3)
!dist_ang_tmp(5)=dot_product(latvec(:,3),latvec(:,1))/dist_ang(3)/dist_ang(1)
!dist_ang_tmp(6)=dot_product(latvec(:,1),latvec(:,2))/dist_ang(1)/dist_ang(2)
!dist_ang(4)=180.d0/pi*acos(max(min(dist_ang_tmp(4),1.d0),-1.d0))
!dist_ang(5)=180.d0/pi*acos(max(min(dist_ang_tmp(5),1.d0),-1.d0))
!dist_ang(6)=180.d0/pi*acos(max(min(dist_ang_tmp(6),1.d0),-1.d0))
!
!
!do i=4,6
!if(isnan(dist_ang(i))) then
! write(*,*) i,dist_ang_tmp(i),dist_ang(i),acos(dist_ang(i))
!! stop
!endif
!enddo
!
!end subroutine
!
! subroutine invertmat(mat,matinv,n)
! implicit none
! real(8),intent(in) :: mat(n,n)
! integer               :: n
! real(8),allocatable   :: WORK(:)
! real(8)               :: matinv(n,n),det(3),a(n,n),div
! integer               :: IPIV(n), INFO 
! integer               :: LDWORK
! !Here only for a 3*3 matrix
! a=mat
! div=(a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)+a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)) 
! div=1.d0/div
!      matinv(1,1) = (a(2,2)*a(3,3)-a(2,3)*a(3,2))*div
!      matinv(1,2) =-(a(1,2)*a(3,3)-a(1,3)*a(3,2))*div
!      matinv(1,3) = (a(1,2)*a(2,3)-a(1,3)*a(2,2))*div
!      matinv(2,1) =-(a(2,1)*a(3,3)-a(2,3)*a(3,1))*div
!      matinv(2,2) = (a(1,1)*a(3,3)-a(1,3)*a(3,1))*div
!      matinv(2,3) =-(a(1,1)*a(2,3)-a(1,3)*a(2,1))*div
!      matinv(3,1) = (a(2,1)*a(3,2)-a(2,2)*a(3,1))*div
!      matinv(3,2) =-(a(1,1)*a(3,2)-a(1,2)*a(3,1))*div
!      matinv(3,3) = (a(1,1)*a(2,2)-a(1,2)*a(2,1))*div
! end subroutine 
!
! subroutine backtocell(nat,latvec,rxyz)
! !This subroutine will transform back all atoms into the periodic cell
! !defined by the 3 lattice vectors in latvec=[v1.v2.v3]
! implicit none
! integer:: nat,i,iat,j
! real(8) :: latvec(3,3), rxyz(3,nat), crossp(3),a(3),b(3), nvec(3,3), dist(6),eps,count
! logical:: neccesary
! !To really be on the safe side, the translation vector can be shortened by  a factor eps in order
! !to get the atom into the cell. 
!!  eps=1.d-10
!  eps=1.d-15
! ! eps=1.d0-eps
! count=0.d0
! neccesary=.true.
! do while(neccesary)
! neccesary=.false.
! count=count+1.d0
! !generate 3 normal vectors of the 3 planes
! call nveclatvec(latvec,nvec)
! do iat=1,nat
! !3 planes through origin (xy,yz,zx)
! do i=1,3
! dist(i)=DOT_PRODUCT(rxyz(:,iat),nvec(:,i))
!! if(dist(i).lt.0.d0) then
! if(dist(i).lt.-abs(dist(i))*eps) then
!! write(*,*) "unten 1",i,iat
! rxyz(:,iat)=rxyz(:,iat)+latvec(:,mod(i+1,3)+1)!*eps
! neccesary=.true.
! endif
!
! !3 planes on top/side/back (xy,yz,zx)
! dist(i+3)=DOT_PRODUCT(rxyz(:,iat)-latvec(:,mod(i+1,3)+1),nvec(:,i))
!! if(dist(i+3).gt.0.d0) then
! if(dist(i+3).gt.abs(dist(i+3))*eps) then
!! write(*,*) "unten 1",i,iat
!! if(dist(i+3).gt.eps) then
!! write(*,*) "oben 2",i,iat
! rxyz(:,iat)=rxyz(:,iat)-latvec(:,mod(i+1,3)+1)!*eps
! neccesary=.true.
! endif
! enddo
! enddo
! if(count.gt.1.d6) stop "Too many iterations in back-to-cell"
! enddo
! end subroutine
!
!
!!************************************************************************************
!!************************************************************************************
!
! subroutine nveclatvec(latvec,nvec)
! !Will calculate the normalized normal vector to the 3 planes of the cell
! implicit none
! real*8, intent(in) :: latvec(3,3)
! real*8, intent(out):: nvec(3,3)
! real*8             :: a(3),b(3),crossp(3),norm
! integer:: i
! do i=1,3
! a=latvec(:,i)
! b=latvec(:,mod(i,3)+1)
! call cross_product(a,b,crossp)
! norm=dsqrt(crossp(1)*crossp(1)+crossp(2)*crossp(2)+crossp(3)*crossp(3))
! nvec(:,i)=crossp(:)/norm
! enddo
! end
!!************************************************************************************
!
! subroutine cross_product(a,b,crossp)
! !a very simple implementation of the cross product
! implicit none
! real(8)::a(3),b(3)
! real(8)::crossp(3)
! crossp(1)=a(2)*b(3)-a(3)*b(2)
! crossp(2)=a(3)*b(1)-a(1)*b(3)
! crossp(3)=a(1)*b(2)-a(2)*b(1)
! return
! end subroutine
!
!!************************************************************************************
