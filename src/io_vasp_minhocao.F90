subroutine poscar_getsystem(parini,filename)
!This subroutine simply returns the 
!Number of atoms
!Types of atoms
!Elemental symbol or Z-values
!and allocates the necessary arrays
!Allocations are done on:
!znucl,char_type,amu,rcov,typat,
use mod_parini, only: typ_parini
use global, only: units
implicit none
type(typ_parini), intent(inout):: parini
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
parini%ntypat_global=0
do i=1,250-1!we cannot have more than 250 characters per line, sorry...
if(all_line(i:i)==" ".and.all_line(i+1:i+1).ne." ") then
parini%ntypat_global=parini%ntypat_global+1
endif
enddo
if(parini%ntypat_global.gt.0) ntypat_found=.true.
!Here we allocate the arrays of types, and the character stuff
if(.not.allocated(parini%znucl)) allocate(parini%znucl(parini%ntypat_global))
if(.not.allocated(parini%char_type)) allocate(parini%char_type(parini%ntypat_global))
if(.not.allocated(nitype)) allocate(nitype(parini%ntypat_global))          
if(.not.allocated(parini%amu)) allocate(parini%amu(parini%ntypat_global))
if(.not.allocated(parini%rcov)) allocate(parini%rcov(parini%ntypat_global))
!Check if the character is a number or a real character
      READ(line_vaspversion,*,ERR=99,END=99) CHARAC
      IF (.NOT.(CHARAC>='0' .AND. CHARAC<='9')) THEN
             !We are deling with vasp5 format
             vasp_5=.true.
             read(line_vaspversion,*) parini%char_type(:)
             read(46,*) nitype(:)
      else
             !We are dealing with vasp4 format
             vasp_5=.false.
             read(line_vaspversion,*) nitype(:)
             !Get the chemical character from the first line
             read(firstline,*) parini%char_type(:)
     endif
!Compute total number of atoms and types and all
parini%nat=0
do i=1,parini%ntypat_global
      parini%nat=parini%nat+nitype(i)
enddo
if(parini%nat.gt.0) nat_found=.true.
if(.not.allocated(parini%typat_global)) allocate(parini%typat_global(parini%nat))          
iat=0
do i=1,parini%ntypat_global
    do j=1,nitype(i)
      iat=iat+1
      parini%typat_global(iat)=i
    enddo
    call symbol2znucl(parini%amu(i),parini%rcov(i),parini%char_type(i),parini%znucl(i))
enddo
znucl_found=.true.

!Be careful: in VASP, you can have the same element multiple times not in sequence,
!but since the POTCAR the also must contain multiple instances of the same element
!the elements will be treated to be "different" also in MHM
99 continue
close(46)
end subroutine

!************************************************************************************

subroutine write_atomic_file_poscar(parini,filename,nat,units,xred,latvec0,fcart,strten,char_type,&
           &ntypat,typat,fixat,fixlat,energy,pressure,printval1,printval2)
!This routine will write the file "filename" in vasp file format
!The input unit will always be in atomic units (bohr, hartree), but the output can be specified by the vaule in "units"
!So if units==angstroem, the file will be converted to angstroem
!   if units==bohr, the positions will not be changed
use defs_basis, only: Ha_eV,Bohr_Ang,HaBohr3_GPa
use global, only: reduced
use mod_parini, only: typ_parini
implicit none
type(typ_parini), intent(in):: parini
integer:: nat,natin,iat,ntypat,typat(nat)
character(*):: filename,units
real(8):: pos(3,nat),xred(3,nat),latvec(3,3),latvec0(3,3),dproj(6),rotmat(3,3),v(3,3),ucvol,fcart(3,nat),strten(6)
real(8):: angbohr,hartree2ev,in_GPA,in_ang3,int_press
real(8):: energy, etotal, enthalpy, enthalpy_at,pressure,printval1,printval2
character(2):: char_type(ntypat)
logical:: fixat(nat),fixlat(7)
integer:: count,jat,ntype
real(8):: pressure_gpa,pressconv
character(3), allocatable:: typatt(:)
integer, allocatable:: nkindsat(:)
character(40):: type_format,formatt
logical::new

units="angstroem"

if(trim(units)=="angstroem") then
  angbohr=Bohr_Ang
  hartree2ev=Ha_eV
  in_GPA=HaBohr3_GPa
  in_ang3=Bohr_Ang**3
  int_press=HaBohr3_GPa/Ha_eV*Bohr_Ang**3
elseif(trim(units)=="bohr") then
  angbohr=1.d0
  hartree2ev=1.d0
  in_GPA=1.d0
  in_ang3=1.d0
  int_press=1.d0
else
  stop "Wrong file unit format"
endif

latvec=latvec0

ntype=ntypat
allocate(nkindsat(ntype),typatt(ntype))
nkindsat=0

do iat=1,ntype
new=.false.
do jat=1,nat
  if(typat(jat)==iat) then
      nkindsat(iat)=nkindsat(iat)+1
      typatt(iat)=trim(char_type(typat(jat)))
      new=.true.
  endif
  if(new.and.typat(jat).lt.iat) stop "Atoms must be ordered!!!"
enddo
enddo

!Must convert latvec to anstroem
latvec=latvec*angbohr

open(unit=46,file=trim(filename))
write(46,*) trim(filename)
write(46,*) 1
write(46,'(3(1x,f15.10))') latvec(:,1)
write(46,'(3(1x,f15.10))') latvec(:,2)
write(46,'(3(1x,f15.10))') latvec(:,3)
write(formatt,'(a,i5,a)') "(",ntype,"(a5))"
write(46,trim(formatt)) typatt(:)
write(formatt,'(a,i5,a)') "(",ntype,"(i5))"
!write(2,trim(formatt)) (nkindsat(iat),  iat=1,ntype)
write(46,trim(formatt)) nkindsat(:)
write(46,*) "Direct"
do iat=1,nat
write(46,*)xred(:,iat)
enddo
write(46,*) " "
if(parini%verb.ge.3) then
do iat=1,nat
write(46,*)fcart(:,iat)*hartree2ev/angbohr
enddo
endif
close(46)

end subroutine

!************************************************************************************

subroutine read_atomic_file_poscar(filename,nat,units,xred,latvec,fcart,strten,&
           &fixat,fixlat,readfix,fragarr,readfrag,printval1,printval2)
use defs_basis, only: Bohr_Ang,Ha_eV
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
  angbohr=1.d0/Bohr_Ang
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

