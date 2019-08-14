subroutine enthalpyrelax(parini,parres,latvec,xred,tolmin,tolmax,ntol,findsym)
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
       real(8):: latvec(3,3),xred(3,parini%nat),pinit,pfinal,psteps,pcur,latvec0(3,3),xred0(3,parini%nat),vel_vol_in
       real(8):: vel_in(3,parini%nat),vel_lat_in(3,3),fcart(3,parini%nat),strten(6)
       real(8):: counter,count_geopt,enthalpy,energy,vol,ext_press
       character(2)::atom
       character(40):: filename,folder
       character(400):: command
       character(7):: fn
       logical:: getwfk,reuse_str
       real(8):: tolmin,tolmax,spgtol
       integer:: ntol,spgint
       logical:: findsym
       open(unit=88,file="enthalpyparams.in")
       read(88,*) pinit,pfinal,psteps,reuse_str !Read which degrees of freedom should be tested, usually only one of them
!Define current directory
  folder=""


latvec0=latvec
xred0=xred
count_geopt=0.d0
open(unit=2,file="Enthalpies")
write(2,'(a)')" #   TrgtPressure(GPa)        Enthalpy(eV)              Energy(eV)&
                &                 ExtPressure(GPa)          Volume(ang^3)             SPG"

close(2)
!Now run over all pressures
pcur=pinit
do while (pcur.le.pfinal)
!pcur=pinit,pfinal,psteps

write(*,'(a,es15.7)') " # Now running at pressure GPa ",pcur
if(.not.reuse_str) latvec=latvec0
if(.not.reuse_str) xred=xred0
!First remove all posgeopt files
call system('rm -f posgeopt.*.ascii')
 parres%target_pressure_gpa=pcur                                   !Target pressure in GPA
 parres%target_pressure_habohr=parres%target_pressure_gpa/HaBohr3_GPA
!Call geometry optimizer 
  vel_in=0.d0
  vel_lat_in=0.d0
   iprec=1
   if(parini%geopt_ext) then
     call geopt_external(parini,latvec,xred,fcart,strten,energy,iprec,ka,kb,kc,counter)
   else
     if(parini%paropt_geopt%approach=="RBFGS")  call GEOPT_RBFGS_MHM(parini,parres,latvec,xred,fcart,strten,energy,iprec,counter,folder)
     if(parini%paropt_geopt%approach=="MBFGS")  call GEOPT_MBFGS_MHM(parini,parres,latvec,xred,fcart,strten,energy,iprec,counter,folder)
     if(parini%paropt_geopt%approach=="FIRE")   call GEOPT_FIRE_MHM(parini,parres,latvec,xred,fcart,strten,vel_in,vel_lat_in,vel_vol_in,energy,iprec,counter,folder)
     if(parini%paropt_geopt%approach=="SQNM")   call     GEOPT_SQNM(parini,parres,latvec,xred,fcart,strten,energy,iprec,counter,folder)
     if(parini%paropt_geopt%approach=="QBFGS")  call    GEOPT_QBFGS(parini,parres,latvec,xred,fcart,strten,energy,iprec,counter,folder)
     if(parini%paropt_geopt%approach=="SD")     call     GEOPT_SD  (parini,parres,latvec,xred,fcart,strten,energy,iprec,counter,folder)
   endif
!Check for symmetry
   if(findsym) then
      call find_symmetry(parini,parini%nat,xred,latvec,parini%typat_global,tolmin,tolmax,ntol,spgtol,spgint)
   else
      spgint=0
      spgtol=0.d0
   endif
!Update GEOPT counter
    count_geopt=count_geopt+counter     
    write(*,'(a,i7)') " # Counter of GEOPT updated: ", int(count_geopt)
    call get_enthalpy(latvec,energy,parres%target_pressure_habohr,enthalpy)
    call getvol(latvec,vol)
    ext_press=strten(1)+strten(2)+strten(3)
    ext_press=-ext_press/3.d0*HaBohr3_GPa
    write(fn,'(f7.2)') pcur
    fn = repeat( '0', 7-len_trim(adjustl(fn))) // adjustl(fn)
    write(*,*) fn
    filename="posenth."//fn//".ascii"
    call write_atomic_file_ascii(parini,filename,parini%nat,units,xred,latvec,fcart,strten,&
          &parini%char_type,parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,energy,parres%target_pressure_habohr,ext_press,enthalpy)
    open(unit=2,file="Enthalpies",access="append")
    write(2,'(5(1x,es25.15),3x,i5)') parres%target_pressure_gpa,enthalpy*Ha_ev,energy*Ha_ev,ext_press,vol*Bohr_Ang**3,spgint
    close(2)
    command="mkdir GPa"//fn
    call system(trim(command)) 
    command="cp posgeopt.*.ascii GPa"//fn
    call system(trim(command)) 
pcur=pcur+psteps
enddo
write(*,'(a)') " # Finished relaxing structures. Output in GPa folders and in the file 'Enthalpies'"
write(*,'(a)') " # Have a nice day!"


end subroutine
