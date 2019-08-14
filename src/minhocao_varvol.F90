subroutine varvol(parini,parres,latvec,xred,tolmin,tolmax,ntol,findsym)
 use global, only: units
 use defs_basis
 use interface_code
 use mod_parini, only: typ_parini
! Main program to test potential subroutines
       implicit none
       type(typ_parini), intent(in):: parini
       type(typ_parini), intent(inout):: parres
       integer::  ka,kb,kc,iprec
! nat: number of atoms (e.g. vacancy)
!       parameter(nint=2*32*1+1)
       real(8):: latvec(3,3),xred(3,parini%nat),latvec0(3,3),xred0(3,parini%nat)
       real(8):: vel_in(3,parini%nat),vel_lat_in(3,3),fcart(3,parini%nat),strten(6)
       real(8):: counter,count_geopt,enthalpy,energy,vol,ext_press
       character(2)::atom
       character(40):: filename,folder
       character(400):: command
       character(12):: fn
       logical:: getwfk,reuse_str
       real(8):: tolmin,tolmax,spgtol
       integer:: ntol,spgint,itime
       logical:: findsym,is_percentage
       real(8):: vfmin,vfmax,vfsteps,vmin,vmax,vsteps,vcur,vol0,vinit,vfinal,convert,alpha
       open(unit=88,file="varvol.in")
       read(88,*) vfmin,vfmax,vfsteps,is_percentage 
!Convert units: usually the unis are given as in the poscur. Assume everything is in volume/atom
if(trim(units)=="angstroem") then
  convert=Bohr_Ang
else
  convert =1.d0
endif

!Define current directory
folder=""
latvec0=latvec
xred0=xred
call getvol(latvec0,vol0)
  if(.not.is_percentage) then
     vfmin=parini%nat*vfmin/vol0/convert**3
     vfmax=parini%nat*vfmax/vol0/convert**3
     vfsteps=parini%nat*vfsteps/vol0/convert**3
  endif
   

open(unit=2,file="VolEnergies")
write(2,'(a)')" #   TrgtPressure(GPa)        Enthalpy(eV)              Energy(eV)&
                &                 Energy(eV/atom)           ExtPressure(GPa)          Volume(ang^3)             Volume(ang^3/atom)        Volume(%)                 SPG"

close(2)
!Now run over all volumes
vinit=vfmin*vol0
vfinal=vfmax*vol0
vsteps=vfsteps*vol0
vcur=vinit
itime=1
do while (vcur.le.vfinal+1.d-10)
   write(*,'(a,es15.7)') " # Now running at volume per atom ",vcur*convert**3/parini%nat
!Translate vcur to a lattice vector
   alpha=vcur/vol0
   alpha=alpha**(1.d0/3.d0)
   latvec=latvec0*alpha
!First remove all posgeopt files
   call system('rm -f posgeopt.*.ascii')
       if(itime==1)  then
           getwfk=.false.
       else
           getwfk=.true.
       endif
   call get_energyandforces_single(parini,parres,latvec,xred0,fcart,strten,energy,iprec,getwfk)
!Check for symmetry
   if(findsym) then
      call find_symmetry(parini,parini%nat,xred,latvec,parini%typat_global,tolmin,tolmax,ntol,spgtol,spgint)
   else
      spgint=0
      spgtol=0.d0
   endif
   call get_enthalpy(latvec,energy,parres%target_pressure_habohr,enthalpy)
   call getvol(latvec,vol)
   ext_press=strten(1)+strten(2)+strten(3)
   ext_press=-ext_press/3.d0*HaBohr3_GPa
   if(is_percentage) then
      write(fn,'(f12.3)') vcur/vol0
   else
      write(fn,'(f12.3)') vcur*convert**3/parini%nat
   endif
   fn = repeat( '0', 12-len_trim(adjustl(fn))) // adjustl(fn)
   filename="posvol."//fn//".ascii"
   call write_atomic_file_ascii(parini,filename,parini%nat,units,xred,latvec,fcart,strten,&
         &parini%char_type,parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,energy,parres%target_pressure_habohr,ext_press,enthalpy)
   open(unit=2,file="VolEnergies",access="append")
   write(2,'(7(1x,es25.15),1x,es25.8,3x,i5)') parres%target_pressure_gpa,enthalpy*Ha_ev,energy*Ha_ev,energy*Ha_ev/parini%nat,ext_press,vol*Bohr_Ang**3,vol*Bohr_Ang**3/parini%nat,vol/vol0,spgint
   close(2)
   vcur=vcur+vsteps
   itime=itime+1
enddo
write(*,'(a)') " # Finished relaxing structures. Output in GPa folders and in the file 'VolEnergies'"
write(*,'(a)') " # Have a nice day!"


end subroutine
