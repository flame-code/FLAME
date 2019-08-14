subroutine poslowrelax(parini,parres,latvec,xred,tolmin,tolmax,ntol)
 use global, only: units
 use defs_basis
 use interface_code
 use mod_parini, only: typ_parini
! Main program to test potential subroutines
       implicit none
       type(typ_parini), intent(in):: parini
       type(typ_parini), intent(inout):: parres
!       use parameters 
       integer::  iprec,nstruct,i,ntol,spgint
! nat: number of atoms (e.g. vacancy)
!       parameter(nint=2*32*1+1)
       real(8):: latvec(3,3),xred(3,parini%nat),pinit,pfinal,psteps,latvec0(3,3),xred0(3,parini%nat),vel_vol_in
       real(8):: vel_in(3,parini%nat),vel_lat_in(3,3),fcart(3,parini%nat),strten(6)
       real(8):: counter,count_geopt,enthalpy,energy,vol,ext_press,tolmin,tolmax,spgtol_pos,spg_pos
       character(2)::atom
       character(40):: filename,filenameout,folder
       character(400):: command
       character(7):: fn
       character(5):: fn5
       logical:: getwfk,reuse_str,file_exists,readfix,readfrag
!Define current directory
  folder=""

open(unit=4,file="poslowrelax.in")
read(4,*) nstruct
close(4)

count_geopt=0.d0
open(unit=2,file="Enthalpies_poslow")
write(2,'(a)')" # poslow.index     TrgtPressure(GPa)        Enthalpy(eV)              Energy(eV)&
              &                 ExtPressure(GPa)          Volume(ang^3)    spg"
close(2)
!Now run over all poslows found in the folder
do i=1,nstruct
        write(fn5,'(i5.5)') i
        filename = 'poslow'//fn5//'.ascii'
        INQUIRE(FILE=trim(filename), EXIST=file_exists)  
        readfix=.false.
        readfrag=.false.
if(file_exists) then
   call read_atomic_file_ascii(filename,parini%nat,units,xred,latvec,fcart,strten,parini%fixat,parini%fixlat,readfix,parini%fragarr,readfrag,enthalpy,energy)
   write(*,'(a,es15.7)') " # Now running "//trim(filename)//" at pressure GPa ",parres%target_pressure_gpa
   !First remove all posgeopt files
   call system('rm -f posgeopt.*.ascii')
   !Call geometry optimizer 
     vel_in=0.d0
     vel_vol_in=0.d0
     vel_lat_in=0.d0
      iprec=1
      if(parini%geopt_ext) then
        call geopt_external(parini,latvec,xred,fcart,strten,energy,iprec,parres%ka,parres%kb,parres%kc,counter)
      else
        if(parini%paropt_geopt%approach=="RBFGS")  call GEOPT_RBFGS_MHM(parini,parres,latvec,xred,fcart,strten,energy,iprec,counter,folder)
        if(parini%paropt_geopt%approach=="MBFGS")  call GEOPT_MBFGS_MHM(parini,parres,latvec,xred,fcart,strten,energy,iprec,counter,folder)
        if(parini%paropt_geopt%approach=="FIRE")   call GEOPT_FIRE_MHM(parini,parres,latvec,xred,fcart,strten,vel_in,vel_lat_in,vel_vol_in,energy,iprec,counter,folder)
        if(parini%paropt_geopt%approach=="SQNM")   call     GEOPT_SQNM(parini,parres,latvec,xred,fcart,strten,energy,iprec,counter,folder)
        if(parini%paropt_geopt%approach=="QBFGS")  call    GEOPT_QBFGS(parini,parres,latvec,xred,fcart,strten,energy,iprec,counter,folder)
        if(parini%paropt_geopt%approach=="SD")     call     GEOPT_SD  (parini,parres,latvec,xred,fcart,strten,energy,iprec,counter,folder)
      endif
!Check for symmetry
   if(parini%findsym) then
      call find_symmetry(parini,parini%nat,xred,latvec,parini%typat_global,tolmin,tolmax,ntol,spgtol_pos,spgint)
      spg_pos=real(spgint,8)
   else
      spg_pos=0.d0
      spgtol_pos=0.d0
   endif
!Update GEOPT counter
       count_geopt=count_geopt+counter     
       write(*,'(a,i7)') " # Counter of GEOPT updated: ", int(count_geopt)
       call get_enthalpy(latvec,energy,parini%target_pressure_habohr,enthalpy)
       call getvol(latvec,vol)
       ext_press=strten(1)+strten(2)+strten(3)
       ext_press=-ext_press/3.d0*HaBohr3_GPa
       filenameout=trim(filename)//"relaxed.ascii"
       call write_atomic_file_ascii(parini,filenameout,parini%nat,units,xred,latvec,fcart,strten,&
             &parini%char_type,parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,energy,parini%target_pressure_habohr,ext_press,enthalpy)
       open(unit=2,file="Enthalpies_poslow",access="append")
       write(2,'(a,5(1x,es25.15),i5)')trim(filename),parres%target_pressure_gpa,&
             &enthalpy*Ha_ev,energy*Ha_ev,ext_press,vol*Bohr_Ang**3,spgint
       close(2)
       command="mkdir "//trim(filename)//"_DIR"
       call system(trim(command)) 
       command="cp posgeopt.*.ascii "//trim(filename)//"_DIR"
       call system(trim(command)) 
endif
enddo
write(*,'(a)') " # Finished relaxing structures. Output in _DIR folders and in the file 'Enthalpies'"
write(*,'(a)') " # Have a nice day!"


end subroutine
