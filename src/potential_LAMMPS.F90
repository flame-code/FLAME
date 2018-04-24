module interface_lammps
  use global
  use defs_basis
#if defined(HAVE_LAMMPS)
  use mpi
  use LAMMPS
  use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int
#endif

  implicit none
  private
  public :: &
   init_lammps,&
   call_lammps,&
   count_lammps
   ! The following line is unnecessary, as I have included these three entities
   ! with the LAMMPS module, but I leave them in anyway to remind people where
   ! they came from


   ! Notes:
   !  * If LAMMPS returns a scalar that is allocated by the library interface
   !     (see library.cpp), then that memory is deallocated automatically and
   !     the argument to lammps_extract_fix must be a SCALAR.
   !  * If LAMMPS returns a pointer to an array, consisting of internal LAMMPS
   !     data, then the argument must be an interoperable Fortran pointer.
   !     Interoperable means it is of type INTEGER (C_INT) or of type
   !     REAL (C_DOUBLE) in this context.
   !  * Pointers should NEVER be deallocated, as that would deallocate internal
   !     LAMMPS data!
   !  * Note that just because you can read the values of, say, a compute at
   !     any time does not mean those values represent the "correct" values.
   !     LAMMPS will abort you if you try to grab a pointer to a non-current
   !     entity, but once it's bound, it's your responsibility to check that
   !     it's current before evaluating.
   !  * IMPORTANT:  Two-dimensional arrays (such as 'x' from extract_atom)
   !     will be transposed from what they might look like in C++.  This is
   !     because of different bookkeeping conventions between Fortran and C
   !     that date back to about 1970 or so (when C was written).
   !  * Arrays start from 1, EXCEPT for mass from extract_atom, which
   !     starts from 0.  This is because the C array actually has a blank
   !     first element (and thus mass[1] corresponds to the mass of type 1)
   integer:: count_lammps = 0
   integer:: nat_lammps
   real(8):: cut_lammps,cut_lammps2   !Maximal cutoff used in lammps, since only first image convention is used!
   type (C_ptr) :: lmp
   real (C_double), pointer :: compute => NULL()
   real (C_double) :: fix, fix2
   integer (C_int), pointer :: iint
   real (C_double),pointer::  boxxlo
   real (C_double),pointer::  boxxhi
   real (C_double),pointer::  boxylo
   real (C_double),pointer::  boxyhi
   real (C_double),pointer::  boxzlo
   real (C_double),pointer::  boxzhi
   real (C_double),pointer:: xy
   real (C_double),pointer:: xz
   real (C_double),pointer:: yz
   real (C_double), dimension(:), pointer :: compute_v => NULL()
   real (C_int), dimension(:), pointer :: compute_vi => NULL()
   real (C_double), dimension(:,:), pointer :: x => NULL()
   real (C_int), dimension(:), pointer :: d => NULL()
   real (C_double), dimension(:), pointer :: mass => NULL()
   integer, dimension(:), allocatable :: types
   double precision, dimension(:), allocatable :: r
   double precision, allocatable :: bxxlo
   integer :: error, narg, me, nprocs,  iat
   character (len=1024) :: command_line
   character(40):: filename
   logical:: file_exists
   logical:: use_backtocell
 

contains

subroutine init_lammps(parini,nnat)
use mod_parini, only: typ_parini
implicit none
type(typ_parini), intent(in):: parini
integer:: iat,nnat
write(*,'(a)') " # Initiallizing lammps"
! And here's how to to it with a string constant of your choice
   call lammps_open_no_mpi ('lmp -log log.lammps', lmp)
   use_backtocell=.true.
!If there is a file called "in.lammps_standalon", then this is the only file
!that will ever be read by MHM to set up the calculation. The user is
!responsible that the atomic number and types correspond to what we have in our
!poscar file!!!
filename="in.lammps_standalone"
INQUIRE(FILE=trim(filename), EXIST=file_exists)
if(file_exists) then
  write(*,*) "Lammps is reading from standalonte file! MHM and LAMMPS files must be consistent!"
  write(*,*) "The following LAMMPS commands will be executed after reading in.lammps_standalone:"
  write(*,*) "LMP:   compute thermo_virial all pressure thermo_temp virial"
  write(*,*) "LMP:   fix fix_nph all nve"
  write(*,*) "LMP:   thermo 0"
  cut_lammps=0.d0
  cut_lammps2=0.d0
  call lammps_file (lmp, 'in.lammps_standalone')
  nat_lammps=lammps_get_natoms(lmp)
  if (nnat.ne.nat_lammps) stop "Number of atoms in LAMMPS not consistend with MHM!"
  use_backtocell=.false.
  goto 2000
endif  


if(any(parini%znucl(:)==201).or.any(parini%znucl(:)==202)) then
!we are dealing with a LJ system
   call lammps_command (lmp, 'units lj')
else
!simple angstroem units
   call lammps_command (lmp, 'units metal')
endif
!!!Set atomic type
call lammps_file (lmp, 'in.lammps_init')
!   call lammps_command (lmp, 'atom_style charge')
if(parini%bc==1) then
   call lammps_command (lmp, 'boundary p p p')
elseif(parini%bc==2) then
   call lammps_command (lmp, 'boundary p p f')
else
   call lammps_command (lmp, 'boundary f f f')
endif
call lammps_command (lmp, 'atom_modify map array')
!Set up cell
call lammps_command (lmp, 'lattice none 1.')
call lammps_command (lmp, 'region simbox prism 0 1 0 1 0 1 0 0 0 side in units box')
!Set up atom types
write(command_line,'(a,i5,a)') "create_box ",parini%ntypat," simbox"
call lammps_command (lmp, trim(command_line))
do iat=1,parini%ntypat
write(command_line,'(a,i5,es25.15)') "mass ",iat,amu_emass*parini%amu(iat)
call lammps_command (lmp, trim(command_line))
enddo
!Setup all atoms
nat_lammps=nnat
do iat=1,nnat
write(command_line,'(a,i5,a,3(es25.15),a)') "create_atoms ",parini%typat_global(modulo(iat-1,parini%nat)+1), " single " , 0.d0,0.d0,0.d0," units box"
write(*,*) trim(command_line)
call lammps_command (lmp, trim(command_line))
enddo
!Get set commands
filename="in.lammps_set"
INQUIRE(FILE=trim(filename), EXIST=file_exists)
if(file_exists) then
  call lammps_file (lmp, 'in.lammps_set')
endif
   ! Extracts pointer to atom types
!   call lammps_gather_atoms (lmp, 'type', 1, types)
!   print *, 'types is ', types(1:)
!Read in coeficcients from file for lammps
filename="in.lammps"
call get_cut(filename,cut_lammps)
cut_lammps2=2.d0*cut_lammps
call lammps_file (lmp, trim(filename))
write(*,*) trim(filename)//" read"
2000 continue
!Setup fixes
call lammps_command (lmp, "compute thermo_virial all pressure thermo_temp virial")
call lammps_command (lmp, "fix fix_nph all nve")
!call lammps_command (lmp, "neigh_modify every 2 delay 10 check yes page 100000")
!call lammps_command (lmp, "neighbor 2 bin")
!!!!!call lammps_command (lmp, "thermo_style custom step temp pe etotal press vol")
!!command_line="dump myDump all custom 1 dump.lammps  id type x y z fx fy fz"
!!call lammps_command (lmp, trim(command_line))
!!command_line='dump_modify myDump format "%d %d %25.20g %25.20g %25.20g %25.20g %25.20g %25.20g"'
!!call lammps_command (lmp, trim(command_line))
call lammps_command (lmp, "thermo 0")
!!command_line='thermo_style custom step pe press fmax fnorm       xlo xhi ylo yhi zlo zhi xy xz yz'
!!call lammps_command (lmp, trim(command_line))
!!command_line='thermo_modify format 2 %25.20g'
!!call lammps_command (lmp, trim(command_line))


end subroutine

!!!********************************************************************************
subroutine call_lammps(parini,latvec,xred0,fcart,energy,strten)
!!!All H-Atoms have to be at the end of the rxyz-array
use mod_parini, only: typ_parini
implicit none
type(typ_parini), intent(in):: parini
integer:: i,n_in(3),offset_in(3),do_kpt_in,cellsearch,tb_or_meam,n_silicon,iprec,iat,nat_lammps_new
integer:: nec1,nec2,nec3,nec1t,nec2t,nec3t,jat,k,l,m
real(8):: xred(3,parini%nat),xred0(3,parini%nat),fcart(3,parini%nat),strten(6),stress(3,3),energy,count
real(8):: vol,alpha,beta,gamma,a,b,c,rotmat(3,3),fcartblj(3,parini%nat),strtenblj(6),energyblj,strtenvir(6)
real(8):: tilts(6),tiltsm(6),ftot(3),randmov(3,parini%nat)
real(8):: latvec(3,3),latvec_ang(3,3),latvec_tilt(3,3),latvec_tilt_inv(3,3),k_latvec(3,3)
real(8),allocatable:: k_xcart(:,:,:,:,:),k_xred(:,:,:,:,:)
character(5)::dnat1,dnat2
!!call random_number(randmov)
!!randmov=randmov*0.01d0
xred=xred0!+randmov
!Setup the cell, atoms etc and feed it into lammps
latvec_ang=latvec*Bohr_Ang
!Rotate cell according to the tilt conditions
a=sqrt(dot_product(latvec_ang(:,1),latvec_ang(:,1)))
b=sqrt(dot_product(latvec_ang(:,2),latvec_ang(:,2)))
c=sqrt(dot_product(latvec_ang(:,3),latvec_ang(:,3)))
alpha=dot_product(latvec_ang(:,2),latvec_ang(:,3))/b/c
beta=dot_product(latvec_ang(:,3),latvec_ang(:,1))/c/a
gamma=dot_product(latvec_ang(:,1),latvec_ang(:,2))/a/b
tilts(1)=a
tilts(4)=b*gamma
tilts(5)=c*beta
tilts(2)=sqrt(b**2-tilts(4)**2)
tilts(6)=(b*c*alpha-tilts(4)*tilts(5))/tilts(2)
tilts(3)=sqrt(c**2-tilts(5)**2-tilts(6)**2)
latvec_tilt=0.d0
latvec_tilt(1,1)=tilts(1)
latvec_tilt(2,2)=tilts(2)
latvec_tilt(3,3)=tilts(3)
latvec_tilt(1,2)=tilts(4)
latvec_tilt(1,3)=tilts(5)
latvec_tilt(2,3)=tilts(6)
if(use_backtocell) call backtocell(parini%nat,latvec_tilt,xred)

!Extended cell
call n_rep_dim(latvec_tilt,cut_lammps2,nec1,nec2,nec3)
if(parini%verb.gt.0) write(*,'(a,i5,i5,i5)') " #Expanding cell with periodic images to: ",nec1,nec2,nec3
nat_lammps_new=nec1*nec2*nec3*parini%nat
allocate(k_xcart(3,parini%nat,nec1,nec2,nec3),k_xred(3,parini%nat,nec1,nec2,nec3))
call k_expansion(parini,latvec_tilt,xred,nec1,nec2,nec3,k_latvec,k_xcart)

tilts(1)=k_latvec(1,1)
tilts(2)=k_latvec(2,2)
tilts(3)=k_latvec(3,3)
tilts(4)=k_latvec(1,2)
tilts(5)=k_latvec(1,3)
tilts(6)=k_latvec(2,3)

!!open(unit=222,file="poslammps_init.ascii")
!!write(222,*) "Comment"
!!write(222,*) tilts(1),tilts(4),tilts(2)
!!write(222,*) tilts(5),tilts(6),tilts(3)
!!  do k=1,nec1
!!  do l=1,nec2
!!  do m=1,nec3
!!  do iat=1,nat
!!   write(222,*) k_xcart(:,iat,k,l,m),"Si"
!!  enddo
!!  enddo
!!  enddo
!!  enddo
!!close(222)

!Correct the tilt
tiltsm=tilts
if(tiltsm(4).ge.0.5d0*tilts(1).or.tiltsm(4).lt.-0.5d0*tilts(1)) then
tiltsm(4)=modulo(tilts(4)-0.5d0*tilts(1),tilts(1))-0.5d0*tilts(1)
endif
if(tiltsm(5).ge.0.5d0*tilts(1).or.tiltsm(5).lt.-0.5d0*tilts(1)) then
tiltsm(5)=modulo(tilts(5)-0.5d0*tilts(1),tilts(1))-0.5d0*tilts(1)
endif
if(tiltsm(6).ge.0.5d0*tilts(2).or.tiltsm(6).lt.-0.5d0*tilts(2)) then
tiltsm(6)=modulo(tilts(6)-0.5d0*tilts(2),tilts(2))-0.5d0*tilts(2)
tiltsm(5)=tiltsm(5)+(tiltsm(6)-tilts(6))/tiltsm(2)*tiltsm(4)
  if(tiltsm(5).ge.0.5d0*tilts(1).or.tiltsm(5).lt.-0.5d0*tilts(1)) then
    tiltsm(5)=modulo(tiltsm(5)-0.5d0*tiltsm(1),tiltsm(1))-0.5d0*tiltsm(1)
  endif
endif

!Old way to get the tilts
!tiltsm(4)=mod(tilts(4),0.5d0*tilts(1))
!tiltsm(5)=mod(tilts(5),0.5d0*tilts(1))
!tiltsm(6)=mod(tilts(6),0.5d0*tilts(2))

!Check tilts on latvec
if(any(tilts(:).ne.tiltsm(:))) then
  k_latvec=0.d0
  k_latvec(1,1)=tiltsm(1)
  k_latvec(2,2)=tiltsm(2)
  k_latvec(3,3)=tiltsm(3)
  k_latvec(1,2)=tiltsm(4)
  k_latvec(1,3)=tiltsm(5)
  k_latvec(2,3)=tiltsm(6)
  if(use_backtocell)  call rxyz_cart2int(k_latvec,k_xred,k_xcart,nat_lammps_new)
  if(use_backtocell)  call backtocell(nat_lammps_new,k_latvec,k_xred)
  if(use_backtocell)  call rxyz_int2cart(k_latvec,k_xred,k_xcart,nat_lammps_new)
!Check if necs are still fine
  call n_rep_dim(k_latvec,cut_lammps,nec1t,nec2t,nec3t)
  if(nec1t.gt.1.or.nec2t.gt.1.or.nec3t.gt.1) then
    write(*,*) "Not enough nec in lammps"
    write(*,*) nec1,nec2,nec3
    write(*,*) nec1t,nec2t,nec3t
    stop
  endif
endif
tilts=tiltsm

latvec_tilt(:,1)=k_latvec(:,1)/real(nec1,8)
latvec_tilt(:,2)=k_latvec(:,2)/real(nec2,8)
latvec_tilt(:,3)=k_latvec(:,3)/real(nec3,8)


!Delete/create atoms depending on wether we need more or less atoms to get the correct expansions
if(nat_lammps_new-nat_lammps.lt.0) then
!Close and reinitiallize lammps
!!  call lammps_close (lmp)
!!  call init_lammps(nat_lammps_new)
!group atoms to delete
 write(dnat1,'(i5)') nat_lammps_new+1
 write(dnat2,'(i5)') nat_lammps;dnat2=adjustl(dnat2)
 write(command_line,'(a)') "group delat id    "//dnat1//":"//dnat2
! write(*,*) "now "//trim(command_line)
 call lammps_command (lmp, trim(command_line))
 call lammps_command (lmp, "delete_atoms group delat compress no")
! write(*,*) trim("now delete_atoms group delat compress no")
elseif(nat_lammps_new-nat_lammps.gt.0) then
 do iat=nat_lammps+1,nat_lammps_new
 write(command_line,'(a,i5,a,3(es25.15),a)') "create_atoms ",parini%typat_global(modulo(iat-1,parini%nat)+1), " single " , 0.d0,0.d0,0.d0," units box"
 call lammps_command (lmp, trim(command_line))
 enddo
!Get set commands
filename="in.lammps_set"
INQUIRE(FILE=trim(filename), EXIST=file_exists)
if(file_exists) then
  call lammps_file (lmp, 'in.lammps_set')
endif
endif
nat_lammps=nat_lammps_new 


!open(unit=222,file="poslammps.ascii")
!write(222,*) "Comment"
!write(222,*) tilts(1),tilts(4),tilts(2)
!write(222,*) tilts(5),tilts(6),tilts(3)
!  do k=1,nec1
!  do l=1,nec2
!  do m=1,nec3
!  do iat=1,nat
!   write(222,*) k_xcart(:,iat,k,l,m),"Si"
!  enddo
!  enddo
!  enddo
!  enddo
!close(222)

!Setup the cell
!call lammps_command (lmp, "neigh_modify every 2 delay 10 check yes page 100000")
!call lammps_command (lmp, "neighbor 2 bin")
!Perform it in multiple steps
!First expand xx,yy,zz
write(command_line,'(a,3(a,es8.1,es25.15),a)')&!,3(a,es25.15),a)') &
     &'change_box all ', ' x final ',0.d0,5.d0*tilts(1),' y final ',0.d0,5.d0*tilts(2),' z final ',0.d0,5.d0*tilts(3),&
     &' units box'
call lammps_command (lmp, trim(command_line))
!Now update xy,xz,yz
write(command_line,'(a,3(a,es25.15),a)') &
     &'change_box all ', ' xy final ',tilts(4),' xz final ',tilts(5),' yz final ',tilts(6),&
     &' units box'
call lammps_command (lmp, trim(command_line))
!Now update back the correct xx,yy,zz
write(command_line,'(a,3(a,es8.1,es25.15),a)')&!,3(a,es25.15),a)') &
     &'change_box all ', ' x final ',0.d0,tilts(1),' y final ',0.d0,tilts(2),' z final ',0.d0,tilts(3),&
     &' units box'
call lammps_command (lmp, trim(command_line))

!   call lammps_extract_global(iint, lmp, 'natoms')
!   print *, 'we have',iint,'atoms'
!Setup the atomic positions
! Allocates an array and assigns all positions to it
  call lammps_gather_atoms (lmp, 'x', 3, r)
  jat=1
  do k=1,nec1
  do l=1,nec2
  do m=1,nec3
  do iat=1,parini%nat
!   write(*,'(a,3es25.15,a)') "create_atoms 1 single   ",k_xcart(:,iat,k,l,m),"     units box"
   r(jat)=k_xcart(1,iat,k,l,m);jat=jat+1
   r(jat)=k_xcart(2,iat,k,l,m);jat=jat+1
   r(jat)=k_xcart(3,iat,k,l,m);jat=jat+1
  enddo
  enddo
  enddo
  enddo

!Puts those position data back
   call lammps_scatter_atoms (lmp, 'x', r)
!   write(*,*) "SIZE R",size(r)
!Puts data back for velocity
   r=0.d0
   call lammps_scatter_atoms (lmp, 'v', r)
!Run the calculation
   call lammps_command (lmp, 'run 0')
!This extracts the vector stress tensor of compute thermo_temp
   call lammps_extract_compute (compute_v, lmp, 'thermo_virial', 0, 1)
!   print *, 'Stress is ', compute_v(1:6)
   strten(1)=-compute_v(1)
   strten(2)=-compute_v(2)
   strten(3)=-compute_v(3)
   strten(4)=-compute_v(6)
   strten(5)=-compute_v(5)
   strten(6)=-compute_v(4)
!write(*,*) "strten", strten
!This is not very clear yet...
if((any(parini%znucl(:)==201).or.any(parini%znucl(:)==202))) then
   strten=strten/Ha_eV*Bohr_Ang**3  !From ev/ang^3 to Ha/Bohr^3
!   strtenvir=strtenvir/Ha_eV*Bohr_Ang**3
else
   strten=strten*1.d-4/HaBohr3_GPa  !From bar to gigapascal
!   strtenvir=strtenvir*1.d9/HaBohr3_GPa
endif

   
! Extracts a pointer to the arrays of positions for all atoms
   call lammps_extract_atom (x, lmp, 'x')
   if ( .not. associated (x) ) print *, 'x is not associated'
!   print *,"coordinates"
!   do iat=1,nat
!   write(*,'(a,3(es15.7))') 'rxyz is ', x(1:3,iat)  ! Prints x, y, z for atom 1
!   enddo
! Extracts a pointer to the arrays of forces for all atoms
!   call lammps_extract_atom (x, lmp, 'f')
!   if ( .not. associated (x) ) print *, 'x is not associated'
!   call lammps_gather_atoms (lmp, 'x', 3, r)
!do iat=1,nat
!  if(abs(k_xcart(1,iat,1,1,1)-r((iat-1)*3+1)).gt.1.d-10.or.&
!     abs(k_xcart(2,iat,1,1,1)-r((iat-1)*3+2)).gt.1.d-10.or.&
!     abs(k_xcart(3,iat,1,1,1)-r((iat-1)*3+3)).gt.1.d-10) then
!  do jat=1,nat  
!     write(*,'(6es25.15)') k_xcart(:,jat,1,1,1),r((jat-1)*3+1:(jat-1)*3+3)
!  enddo
!  stop
!  endif
!enddo
   call lammps_gather_atoms (lmp, 'f', 3, r)
!Rotate forces back to original cell
   call rotmat_fcart_stress(latvec_ang,latvec_tilt,rotmat)
!!   call invertmat(latvec_tilt,latvec_tilt_inv,3)
!!   rotmat=matmul(latvec_ang,latvec_tilt_inv)
   ftot=0.d0
   do iat=1,parini%nat
!Transform forces back to original cell
     fcart(:,iat)=matmul(rotmat,r((iat-1)*3+1:(iat-1)*3+3))/Ha_eV*Bohr_Ang
     ftot=ftot+fcart(:,iat)
   enddo
   call rotate_stresstensor(strten,rotmat)
!   write(*,'(a,3es25.15)') "FTOT",ftot
   if(abs(sum(ftot)).gt.1.d-10) stop "Total forces not zero!!!"
! This extracts the vector compute of compute thermo_temp
   call lammps_extract_compute (compute, lmp, 'thermo_pe', 0, 0)
!   print *, 'Energy is ', compute
   energy=compute/Ha_eV/real(nec1*nec2*nec3,8)
  deallocate(k_xcart,k_xred)


end subroutine
!********************************************************************************
subroutine n_rep_dim(latvec,cut,nec1,nec2,nec3)
!This subroutine will return how many periodic expansions for each lattice vector direction are necessary for the periodic boundary conditions
!with for the given cut. nec1,nec2,nec3 for latvec(:,1),latvec(:,2),latvec(:,3)
implicit none
real*8 :: latvec(3,3),cut,nvec(3,3),point(3),point0(3),dist(3),eps,dd
integer:: i
integer:: nec1,nec2,nec3
! eps=1.d-6
nec1=0
nec2=0
nec3=0
call nveclatvec(latvec,nvec)
point0=(/0.d0,0.d0,0.d0/)
do i=1,3
call dist2plane(latvec(:,mod(i+1,3)+1),nvec(:,i),point0,dist(i))
! write(*,*) "cut",i,cut, dist
enddo
dist=abs(dist)
nec1=int(cut/dist(2))+1
nec2=int(cut/dist(3))+1
nec3=int(cut/dist(1))+1
end subroutine

!*********************************************************************************

subroutine get_cut(filename,cut_lammps)
implicit none
integer:: n,k,io
character(40)::filename
real(8):: cut_lammps
logical:: found
character(150)::all_line
found=.false.
cut_lammps=0.d0
open(unit=778,file=trim(filename))
do while(.true.)
   read(778,'(a150)',end=99)all_line
   n = len_trim(all_line)
   k = index(all_line(1:n),"cutoff")
    if(k.ne.0) then
      k = scan(all_line(1:n),"f",.true.)
      found=.true.
      read(all_line(k+1:n),*,iostat=io) cut_lammps
      if(io.lt.0) stop "Error while reading cutoff in lammps"
      exit
    endif
enddo
99 continue
close(778)
all_line= "The keyword 'cutoff' must be provided in "//trim(filename)//" followed by the maximal atomic cutoff of all potentials used in lammps"
if(.not.found) then
write(*,*)  trim(all_line)
stop
endif
if(cut_lammps.le.0.d0) stop "cut_lammps must be larger than 0!!!"
end subroutine

!!!********************************************************************************
end module
