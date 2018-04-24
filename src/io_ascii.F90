subroutine ascii_getsystem(parini,filename)
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
integer:: i,j,k,n,iat,jat
integer, allocatable:: nitype(:)
real(8):: scaling_tmp,dproj_tmp(6),pos(3)
logical:: new
character(250):: firstline,all_line
character(40) :: filename
character(1)  :: CHARAC
character(2),allocatable:: char_type_tmp(:)
!uses angstroem, but here it actually does not matter
units="angstroem"
!Get nat and allocate some arrays
open(unit=46,file=trim(filename))
read(46,*) parini%nat !The first line contains the number of atoms
read(46,*) dproj_tmp(1:3)
read(46,*) dproj_tmp(4:6)
!Here we allocate a temporary array for the atomic character 
if(.not.allocated(char_type_tmp)) allocate(char_type_tmp(parini%nat))
if(.not.allocated(parini%typat_global)) allocate(parini%typat_global(parini%nat))
!There are some comment lines possible, reduced and/or fixlat
k=0
do iat=1,parini%nat
1010 continue 
    read(46,'(a250)') all_line
    n = len_trim(all_line)
    k = index(all_line(1:n),"reduced")
    if(k.ne.0) then 
         goto 1010
    endif
    k = index(all_line(1:n),"fixlat",back=.true.)
    if(k.ne.0) then
         goto 1010
    endif
    read(all_line,*) pos(:),char_type_tmp(iat)
enddo

!Count how many different atom kinds there are
if(.not.allocated(parini%typat_global)) allocate(parini%typat_global(parini%nat))          
parini%typat_global(1)=1
parini%ntypat_global=1
do iat=2,parini%nat
new=.true.
 do jat=1,iat-1
 if(trim(char_type_tmp(iat))==trim(char_type_tmp(jat))) then
     new=.false.
     parini%typat_global(iat)=parini%typat_global(jat)
 endif
 enddo
 if(new) then
 parini%ntypat_global=parini%ntypat_global+1
 parini%typat_global(iat)=parini%ntypat_global
 endif
enddo

!Here we allocate the arrays of types, and the character stuff
if(.not.allocated(parini%znucl)) allocate(parini%znucl(parini%ntypat_global))
if(.not.allocated(parini%char_type)) allocate(parini%char_type(parini%ntypat_global))
if(.not.allocated(parini%amu)) allocate(parini%amu(parini%ntypat_global))
if(.not.allocated(parini%rcov)) allocate(parini%rcov(parini%ntypat_global))

!Now get the sole atomic characters
do iat=1,parini%nat
  parini%char_type(parini%typat_global(iat))=char_type_tmp(iat)
enddo

!Assign Znucl here
do i=1,parini%ntypat_global
      call symbol2znucl(parini%amu(i),parini%rcov(i),parini%char_type(i),parini%znucl(i))
enddo

99 continue
close(46)
end subroutine

!************************************************************************************


subroutine read_atomic_file_ascii(filename,nat,units,xred,latvec,fcart,strten,&
           &fixat,fixlat,readfix,fragarr,readfrag,printval1,printval2)
!subroutine read_atomic_file_ascii(filename,nat,units,xred,latvec,acell,rprim,printval1,printval2)
!This routine will read the file "filename" ins ascii file format
!The output unit will always be in atomic units (bohr)
!So if units==angstroem, the file will be converted
!   if units==bohr, the positions will not be changed
use defs_basis, only: Bohr_Ang,Ha_eV
use global, only: reduced
implicit none
integer:: nat,natin,iat,ierror,io,n,k,fragarr(nat),fragarr_tmp,lhead(nat),llist(nat),nmol,m,l
logical:: fixat(nat),fixlat(7),readfix,reduced_tmp,readfrag
character(1) :: isfixed,ch_tmp
character(2) :: element
character(40):: filename,units
real(8):: pos(3,nat),xred(3,nat),latvec(3,3),dproj(6),strten(6),fcart(3,nat)
real(8):: angbohr,evhartree,enthalpy_at,printval1,printval2
character(200):: line
fragarr_tmp=-1
reduced_tmp=.false.
if(readfrag) fragarr=-1
if(trim(units)=="angstroem") then
!  angbohr=1.889725989d0
  angbohr=1.d0/Bohr_Ang
  evhartree=1.d0/Ha_eV
elseif(trim(units)=="bohr") then
  angbohr=1.d0
  evhartree=1.d0
else
  stop "Wrong file unit format"
endif

open(unit=46,file=trim(filename),iostat=ierror)
 if (ierror /= 0) then
    write(*,*) ' COULD not read file ',filename
    stop
 end if

!If the input file is called "poscur.ascii", then dont read the input enthalpy and energy, since they might not be provided
if(trim(filename)=="poscur.ascii") then
read(46,*)natin
else
read(46,*)natin, printval1,printval2
endif
if(natin.ne.nat) stop "Number of atoms not consistent with params.in input"
read(46,*) dproj(1:3)
read(46,*) dproj(4:6)
do iat=1,nat
1010  read(46,'(a200)') line

  if(iat==1) then  
    n = len_trim(line)
    k = index(line(1:n),"reduced")
    if(k.ne.0) then 
         reduced_tmp=.true.
         goto 1010
    endif
  endif

  if(iat==1) then  
    n = len_trim(line)
    k = index(line(1:n),"fixlat",back=.true.)
    if(k.ne.0) then
      k = scan(line(1:n),"t",.true.)
      if(readfix) then
         read(line(k+1:n),*,iostat=io) fixlat(:)
         if(io.lt.0) stop "Could not read lattice constraints"
      endif 
      goto 1010
    endif
  endif

  if(readfix) then
     read(line,*,iostat=io) pos(:,iat),element,isfixed
     if(io.lt.0)  isfixed="a"
     if(isfixed=="f" .or. isfixed=="F")  then
         fixat(iat)=.true.
     else
         fixat(iat)=.false.
     endif
     read(line,*,iostat=io) pos(:,iat),element,fragarr_tmp,isfixed
     if(io.lt.0)  isfixed="a"
     if(isfixed=="f" .or. isfixed=="F")  then
         fixat(iat)=.true.
     else
         fixat(iat)=.false.
     endif
  endif
  if(readfrag) then
     read(line,*,iostat=io) pos(:,iat),element,fragarr_tmp
     if(io.lt.0) then 
         fragarr_tmp=-1 
     else
         fragarr(iat)=fragarr_tmp
     endif
     read(line,*,iostat=io) pos(:,iat),element,isfixed,fragarr_tmp
     if(io.lt.0) then 
         fragarr_tmp=-1 
     else
         fragarr(iat)=fragarr_tmp
     endif
  endif
  read(line,*,iostat=io) pos(:,iat)
      if(io.lt.0) then 
         write(*,*) "File read error in",trim(filename)
         stop
      endif 
enddo
if(readfrag) then
     if(any(fixat).and.any(fragarr(:).ge.1)) then
         write(*,*) fragarr
         stop "Fragmented structure not supported with fixed atoms"
     endif
endif
strten=0.d0
fcart=0.d0
if(trim(filename).ne."poscur.ascii") then
  read(46,*,end=99) ch_tmp,strten(:)
  do iat=1,nat
  read(46,'(a200)',end=99) line
             n=len_trim(line)
     if(iat==1) then 
             k=SCAN(line(1:n),"[") 
     else
             k=SCAN(line(1:n),"#") 
     endif
             do l=1,3 
                 m=SCAN(line(k:n),";") 
                 read(line(k+1:k+m-2),*) fcart(l,iat)
                 k=k+m
             enddo
  enddo
  fcart=fcart*evhartree/angbohr
endif
99 close(46)


dproj=dproj*angbohr
call dproj2latvec(dproj,latvec)
if(reduced_tmp) then
xred=pos
else
pos=pos*angbohr
call rxyz_cart2int(latvec,xred,pos,nat)
endif

if(trim(filename)=="poscur.ascii") reduced=reduced_tmp
!!if(readfrag) then
!!call refragment(fragarr,nat)
!!call make_linked_list(fragarr,lhead,llist,nat,nmol)
!!endif
!stop
end subroutine

!************************************************************************************

subroutine write_atomic_file_ascii(parini,filename,nat,units,xred,latvec0,fcart,strten,char_type,&
           &ntypat,typat,fixat,fixlat,energy,pressure,printval1,printval2)
!This routine will write the file "filename" in ascii file format
!The input unit will always be in atomic units (bohr, hartree), but the output can be specified by the vaule in "units"
!So if units==angstroem, the file will be converted to angstroem
!   if units==bohr, the positions will not be changed
use defs_basis, only: Ha_eV,Bohr_Ang,HaBohr3_GPa
use global, only: reduced
use mod_parini, only: typ_parini

implicit none
type(typ_parini), intent(in):: parini
integer:: nat,natin,iat,ntypat,typat(nat),j
character(40):: filename,units
real(8):: pos(3,nat),xred(3,nat),latvec(3,3),latvec0(3,3),dproj(6),rotmat(3,3),v(3,3),ucvol,fcart(3,nat),strten(6)
real(8):: angbohr,hartree2ev,in_GPA,in_ang3,int_press
real(8):: energy, etotal, enthalpy, enthalpy_at,pressure,printval1,printval2,tmp(3)
character(2):: char_type(ntypat)
character(1):: endline
logical:: fixat(nat),fixlat(7)
if(trim(units)=="angstroem") then
!  angbohr=1.d0/1.889725989d0
!  hartree2ev=27.211396132d0
!  in_GPA=29421.033d0
!  in_ang3=0.148184743d0
!  int_press=160.217646200d0
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

do iat=1,nat
   pos(:,iat)=matmul(latvec,xred(:,iat))
enddo

call latvec2dproj(dproj,latvec,rotmat,pos,nat)

!Compute cell volume
 v=latvec
 ucvol= v(1,1)*v(2,2)*v(3,3)-v(1,1)*v(2,3)*v(3,2)-v(1,2)*v(2,1)*v(3,3)+&
      v(1,2)*v(2,3)*v(3,1)+v(1,3)*v(2,1)*v(3,2)-v(1,3)*v(2,2)*v(3,1)
!Convert units
!!!The second entry in the first line is the enthalpy in atomic units
!!enthalpy_at=energy+pressure*ucvol
dproj=dproj*angbohr
pos=pos*angbohr
if(reduced) pos=xred
etotal=energy*hartree2ev
enthalpy=energy*hartree2ev+pressure*in_GPA/int_press*ucvol*in_ang3
ucvol=ucvol*in_ang3

open(unit=46,file=trim(filename))
if(trim(units)=="angstroem") then
  write(46,'(i4,1x,es25.15,1x,es25.15,1x,a20,es15.7,a14,es15.7,a13,es15.7)')&
  & nat,printval1,printval2,"         energy(eV)=",etotal," enthalpy(eV)=",enthalpy,&
  &" ucvol(ang3)=",ucvol
elseif(trim(units)=="bohr") then
  write(46,'(i4,1x,es25.15,1x,es25.15,1x,a20,es15.7,a14,es15.7,a13,es15.7)')&
  & nat,printval1,printval2,"         energy(Ha)=",etotal," enthalpy(Ha)=",enthalpy,&
  &" ucvol(bhr3)=",ucvol
endif
write(46,*) dproj(1:3)
write(46,*) dproj(4:6)
if(reduced) write(46,'(a)') "#keyword:reduced"
if(any(fixlat)) write(46,'(a,7(1x,L))') "#keyword:fixlat ",fixlat(:)
do iat=1,nat
  if(any(fixat(:))) then
      if(fixat(iat))then
         write(46,'(3(1x,es25.15),2x,a2,2x,a1)') pos(:,iat),trim(char_type(typat(iat))),"F"
      else
         write(46,'(3(1x,es25.15),2x,a2)') pos(:,iat),trim(char_type(typat(iat)))
      endif 
  elseif(all(parini%fragarr.gt.0)) then 
      write(46,'(3(1x,es25.15),2x,a2,2x,i5)')       pos(:,iat),trim(char_type(typat(iat))),parini%fragarr(iat)
  else 
      write(46,'(3(1x,es25.15),2x,a2)')       pos(:,iat),trim(char_type(typat(iat)))
  endif
enddo
if(parini%verb.ge.3) then
write(46,'(a,6(1x,es25.15),a)')  "# ",strten(:)," Stress tensor in GPa"

!write the first position
  iat=1
  if (nat==iat) then
     endline=']'
  else
     endline=char(92)
  end if
tmp(:)=matmul(rotmat,fcart(:,1))
if(trim(units)=="angstroem") then
  write(46, "(A,3(1pe25.15,A),a)") "#metaData: forces (eV/Angstroem) =[",(tmp(j)*hartree2ev/angbohr, ";",j=1,3),' '//endline
elseif(trim(units)=="bohr") then
  write(46, "(A,3(1pe25.15,A),a)") "#metaData: forces (Ha/Bohr) =[",(tmp(j)*hartree2ev/angbohr, ";",j=1,3),' '//endline
endif
!then the rest until the second-last
  do iat=2,nat
     if (nat==iat) then
        endline=']'
     else
        endline=char(92)
     end if
     tmp(:)=matmul(rotmat,fcart(:,iat))
     write(46, "(A,3(1pe25.15,A),a)") "# ",(tmp(j)*hartree2ev/angbohr, ";",j=1,3),' '//endline
  end do
!!write(46,'(a,3(1x,es25.15),a)')  "# ",matmul(rotmat,fcart(:,1))*hartree2ev/angbohr," Forces in file units"
!!do iat=2,nat
!!   write(46,'(a,3(1x,es25.15),a)')  "# ",matmul(rotmat,fcart(:,iat))*hartree2ev/angbohr
!!enddo
endif
close(46)
end subroutine

!************************************************************************************
