!*****************************************************************************************
subroutine single_point_task(parini)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr, typ_file_info
    use mod_potential, only: fcalls, perfstatus, potential
    use mod_processors, only: iproc
    use mod_const, only: ev2ha, ang2bohr
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    type(typ_atoms_arr):: atoms_arr
    type(typ_file_info):: file_info
    real(8):: tt1, tt2, fxyz(3)
    integer:: iconf, iat
    logical:: acf_exists
    inquire(file='posinp.acf',exist=acf_exists)
    if(acf_exists) then
        call acf_read_new(parini,'posinp.acf',10000,atoms_arr)
    else
        atoms_arr%nconf=1
        allocate(atoms_arr%atoms(atoms_arr%nconf))
        call read_poscar_for_single_point(parini,atoms_arr%atoms(1))
    endif
    do iconf=1,atoms_arr%nconf
        call set_ndof(atoms_arr%atoms(iconf))
    enddo
    potential=trim(parini%potential_potential)
    file_info%filename_positions='posout.acf'
    file_info%print_force=parini%print_force_single_point
    file_info%file_position='new'
    if(trim(parini%frmt_single_point)/='unknown') then
        file_info%frmt=trim(parini%frmt_single_point)
    endif
    do iconf=1,atoms_arr%nconf
        if(trim(potential)/='netsock' .or. iconf==1) then 
            call init_potential_forces(parini,atoms_arr%atoms(iconf))
        endif
        call cal_potential_forces(parini,atoms_arr%atoms(iconf))
        !call atoms_all%fatall(1:3,1:atoms_all%atoms%nat,iconf)=atoms_all%atoms%fat(1:3,1:atoms_all%atoms%nat)
        if(iconf==1) then
            tt1=0.d0
            tt2=0.d0
        else
            tt1=atoms_arr%atoms(iconf)%epot-atoms_arr%atoms(1)%epot
            tt2=atoms_arr%atoms(iconf)%epot-atoms_arr%atoms(iconf-1)%epot
        endif
        write(*,'(a,i7,e24.15,2f10.5)') 'EPOT',iconf,atoms_arr%atoms(iconf)%epot,tt1,tt2
        !if(parini%print_force_single_point) then
        !    do iat=1,atoms_all%atoms%nat
        !        fxyz(1)=atoms_all%atoms%fat(1,iat)
        !        fxyz(2)=atoms_all%atoms%fat(2,iat)
        !        fxyz(3)=atoms_all%atoms%fat(3,iat)
        !        write(*,'(3f12.8)') fxyz(1),fxyz(2),fxyz(3)
        !    enddo
        !endif
        if(trim(potential)/='netsock' .or. iconf==atoms_arr%nconf) then 
            call final_potential_forces(parini,atoms_arr%atoms(iconf))
        endif
        if (iconf==2)  file_info%file_position='append'
        call acf_write(file_info,atoms=atoms_arr%atoms(iconf),strkey='posout') !to be fixed later by atoms_arr
    enddo
    !call atom_all_deallocate(atoms_all,ratall=.true.,fatall=.true.,epotall=.true.)
end subroutine single_point_task
!*****************************************************************************************
subroutine read_poscar_for_single_point(parini,atoms)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use global, only: nat, ntypat, typat, znucl, char_type, units, amu, rcov
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms):: atoms
    !local variables
    real(8), allocatable:: xred(:,:)
    real(8), allocatable:: fcart(:,:)
    logical, allocatable:: fixat(:)
    integer, allocatable:: fragarr(:)
    real(8):: strten(6), printval1, printval2
    logical:: fixlat(7), readfix, readfrag
    integer:: iat
    call poscar_getsystem_alborz('POSCAR')
    allocate(xred(3,nat),source=0.d0)
    allocate(fcart(3,nat),source=0.d0)
    if(.not.allocated(fixat)) allocate(fixat(nat),source=.false.)
    if(.not.allocated(fragarr)) allocate(fragarr(nat),source=0)
    atoms%nat=nat
    atoms%boundcond='bulk'
    call atom_allocate_old(atoms,nat,0,0)
    call read_atomic_file_poscar_alborz('POSCAR',atoms%nat,units,xred,atoms%cellvec,fcart,strten, &
        fixat,fixlat,readfix,fragarr,readfrag,printval1,printval2)
    call rxyz_int2cart_alborz(atoms%nat,atoms%cellvec,xred,atoms%rat)
    do iat=1,nat
        atoms%sat(iat)=trim(char_type(typat(iat)))
    enddo
    deallocate(fixat)
    deallocate(xred)
    deallocate(fcart)
    deallocate(fragarr)
    if(allocated(znucl)) deallocate(znucl)
    if(allocated(char_type)) deallocate(char_type)
    if(allocated(amu)) deallocate(amu)
    if(allocated(rcov)) deallocate(rcov)
    if(allocated(typat)) deallocate(typat)
end subroutine read_poscar_for_single_point
!*****************************************************************************************
subroutine poscar_getsystem_alborz(filename)
!This subroutine simply returns the 
!Number of atoms
!Types of atoms
!Elemental symbol or Z-values
!and allocates the necessary arrays
!Allocations are done on:
!znucl,char_type,amu,rcov,typat,
use global, only: nat,ntypat,typat,znucl,char_type,units,amu,rcov
implicit none
integer:: i,j,k,n,iat
integer, allocatable:: nitype(:)
real(8):: scaling_tmp,latvec_tmp(3,3)
logical:: vasp_5,red
character(250):: firstline,line_vaspversion,all_line
character(*) :: filename
character(1)  :: CHARAC
logical       :: nat_found,ntypat_found,znucl_found
!Vasp only uses angstroem
units="angstroem"
!Get the two filenames
vasp_5=.true. 
!Get nat and allocate some arrays
open(unit=46,file=trim(filename))
read(46,'(a250)') firstline  !this line is either comment or contains the atomic elements, if the format is vasp4. So lets keep it in memory
!Read the first structure to initialize stuff
!We only read temporary stuff here, which is of no significance. Essentially,
!we are only interested in the number of atoms, types, etc.
read(46,*) scaling_tmp  
read(46,*) latvec_tmp(:,1)
read(46,*) latvec_tmp(:,2)
read(46,*) latvec_tmp(:,3)
!Here we determine if we are dealing with VASP5 format or not
! 6th line, contains either the number of ions or their type
read(46,'(a250)') line_vaspversion
!Here we get how many types we have. We do not care if the elements in this
!line are characters or numbers, so if vasp version 4 or 5
all_line=" "//trim(line_vaspversion)
!write(*,*) all_line
ntypat=0
do i=1,250-1!we cannot have more than 250 characters per line, sorry...
if(all_line(i:i)==" ".and.all_line(i+1:i+1).ne." ") then
ntypat=ntypat+1
endif
enddo
if(ntypat.gt.0) ntypat_found=.true.
!Here we allocate the arrays of types, and the character stuff
if(.not.allocated(znucl)) allocate(znucl(ntypat))
if(.not.allocated(char_type)) allocate(char_type(ntypat))
if(.not.allocated(nitype)) allocate(nitype(ntypat))          
if(.not.allocated(amu)) allocate(amu(ntypat))
if(.not.allocated(rcov)) allocate(rcov(ntypat))
!Check if the character is a number or a real character
      READ(line_vaspversion,*,ERR=99,END=99) CHARAC
      IF (.NOT.(CHARAC>='0' .AND. CHARAC<='9')) THEN
             !We are deling with vasp5 format
             vasp_5=.true.
             read(line_vaspversion,*) char_type(:)
             read(46,*) nitype(:)
      else
             !We are dealing with vasp4 format
             vasp_5=.false.
             read(line_vaspversion,*) nitype(:)
             !Get the chemical character from the first line
             read(firstline,*) char_type(:)
     endif
!Compute total number of atoms and types and shit
nat=0
do i=1,ntypat
      nat=nat+nitype(i)
enddo
if(nat.gt.0) nat_found=.true.
if(.not.allocated(typat)) allocate(typat(nat))          
iat=0
do i=1,ntypat
    do j=1,nitype(i)
      iat=iat+1
      typat(iat)=i
    enddo
    call symbol2znucl(amu(i),rcov(i),char_type(i),znucl(i))
enddo
znucl_found=.true.

!Be careful: in VASP, you can have the same element multiple times not in sequence,
!but since the POTCAR the also must contain multiple instances of the same element
!the elements will be treated to be "different" also in MHM
99 continue
close(46)
end subroutine

!************************************************************************************
subroutine read_atomic_file_poscar_alborz(filename,nat,units,xred,latvec,fcart,strten,&
           &fixat,fixlat,readfix,fragarr,readfrag,printval1,printval2)
use defs_basis, only: Ha_eV
use mod_const, only: bohr2ang
use global, only: reduced
!This file will read the atomic files latvec and xred, but NOT the types etc.
implicit none
integer:: i,ntypat_tmp,nat,natin,iat,ierror,io,n,k,fragarr(nat),fragarr_tmp,lhead(nat),llist(nat),nmol
integer, allocatable:: nitype(:)
logical:: fixat(nat),fixlat(7),readfix,reduced_tmp,readfrag
character(1) :: isfixed,ch_tmp,tmp_ch
character(2) :: element
character(*):: filename,units
real(8):: pos(3,nat),xred(3,nat),latvec(3,3),dproj(6),strten(6),fcart(3,nat)
real(8):: angbohr,evhartree,enthalpy_at,printval1,printval2,scaling
character(200):: line
character(150):: all_line,all_line_tmp,line_vaspversion
logical:: vasp_5,red
character(1)  :: CHARAC

!No way of defining fixlat 
fixlat=.false.

!No fragarr
fragarr = -1

!Vasp only uses angstroem
units="angstroem"
if(trim(units)=="angstroem") then
  angbohr=1.d0/bohr2ang
  evhartree=1.d0/Ha_eV
elseif(trim(units)=="bohr") then
  angbohr=1.d0
  evhartree=1.d0
else
  stop "Wrong file unit format"
endif
open(unit=46,file=trim(filename))
read(46,'(a150)',end=99) all_line_tmp  !this line is either comment or contains the atomic elements, if the format is vasp4. So lets keep it in memory
!Read the first structure to initialize stuff
!We only read temporary stuff here, which is of no significance. Essentially,
!we are only interested in the number of atoms, types, etc.
read(46,*,end=99) scaling  
read(46,*,end=99) latvec(:,1)
read(46,*,end=99) latvec(:,2)
read(46,*,end=99) latvec(:,3)
!Here we determine if we are dealing with VASP5 format or not
! 6th line, contains either the number of ions or their type
read(46,'(a150)',end=99) line_vaspversion
!Here we get how many types we have. We do not care if the elements in this
!line are characters or numbers, so if vasp version 4 or 5
all_line=" "//trim(line_vaspversion)
!write(*,*) all_line
ntypat_tmp=0
do i=1,150-1!we cannot have more than 250 characters per line, sorry...
if(all_line(i:i)==" ".and.all_line(i+1:i+1).ne." ") then
ntypat_tmp=ntypat_tmp+1
endif
enddo
allocate(nitype(ntypat_tmp))
!Check if the character is a number or a real character
      READ(line_vaspversion,*,ERR=99,END=99) CHARAC
      IF (.NOT.(CHARAC>='0' .AND. CHARAC<='9')) THEN
             !We are deling with vasp5 format
             vasp_5=.true.
!             read(line_vaspversion,*) char_type(:)
             read(46,*) nitype(:)
      else
             !We are dealing with vasp4 format
             vasp_5=.false.
             read(line_vaspversion,*) nitype(:)
             !Get the chemical character from the first line
!             read(firstline,*) char_type(:)
     endif
!Compute total number of atoms and types and shit
natin=0
do i=1,ntypat_tmp
      natin=natin+nitype(i)
enddo
!Check if the total number of atoms is the same as in input
if(nat.ne.natin) stop "POSCAR does not contain the same number of atoms as declared!"
!Read positions
read(46,'(a150)',end=99)all_line
!write(*,*) all_line
n = len_trim(all_line)
k = index(all_line(1:n),"S") !In selective dynamics
  if(k.ne.0) then  !ignore this line
  read(46,'(a150)',end=99)all_line
  endif 
all_line=adjustl(all_line)
all_line=trim(all_line)
all_line=all_line(1:1)
if(trim(all_line).eq."D".or.trim(all_line).eq."d") then
  red=.true.
elseif(trim(all_line).eq."C".or.trim(all_line).eq."c") then
  red=.false.
else
  stop "Coordinates must be either Direct or Cartesian"
endif
!!!
!!!if(trim(adjustl(all_line)).eq."Direct".or.trim(adjustl(all_line)).eq."direct") then
!!!  red=.true.
!!!elseif(trim(adjustl(all_line)).eq."Cartesian".or.trim(adjustl(all_line)).eq."cartesian") then
!!!  red=.false.
!!!else
!!!  stop "Coordinates must be either Direct or Cartesian"
!!!endif

!Read positions
do iat=1,nat
   read(46,*,end=99) xred(:,iat)
enddo
  latvec=latvec*scaling
  if(.not.red) then
   pos=xred
   call rxyz_cart2int(latvec,xred,pos,nat)
  endif
!Convert everything to atomic units
  latvec=latvec*angbohr  
99 continue
close(46)
end subroutine

!************************************************************************************
