subroutine rotate_like_crazy(parini,parres,latvec,xred,tolmin,tolmax,ntol)
 use global, only: units
 use defs_basis
 use interface_code
! Main program to test potential subroutines
!       use parameters 
 use mod_parini, only: typ_parini
       implicit none
       type(typ_parini), intent(in):: parini
       type(typ_parini), intent(inout):: parres
       integer::  iprec,nstruct,i,ntol,spgint
! nat: number of atoms (e.g. vacancy)
!       parameter(nint=2*32*1+1)
       real(8):: latvec(3,3),xred(3,parini%nat),pinit,pfinal,psteps,latvec0(3,3),xred0(3,parini%nat),vel_vol_in
       real(8):: vel_in(3,parini%nat),vel_lat_in(3,3),fcart(3,parini%nat),strten(6),axis(3),rotmat(3,3),angle
       real(8):: counter,count_geopt,enthalpy,energy,vol,ext_press,tolmin,tolmax,spgtol_pos,spg_pos
       character(2)::atom
       character(40):: filename,filenameout,folder
       character(400):: command
       character(7):: fn
       character(5):: fn5
       logical:: getwfk,reuse_str,file_exists,readfix,readfrag
!Define current directory
  folder=""
!Read from the rotation file
  open(unit=45,file="rotate.in")
  read(45,*) nstruct
  close(45)

!Now run over all poslows found in the folder
do i=1,nstruct

call random_number(axis(:))
call random_number(angle)
        axis=axis-0.5d0
        angle=(angle-0.5d0)*10.d0
        axis=axis/sqrt(sum(axis(:)**2))
        call rotation(rotmat,angle,axis)
        latvec=matmul(rotmat,latvec)
!Compute energy
        call get_energyandforces_single(parini,parres,latvec,xred,fcart,strten,energy,iprec,getwfk)
        call get_enthalpy(latvec,energy,parini%target_pressure_habohr,enthalpy)
        call getvol(latvec,vol)
        ext_press=strten(1)+strten(2)+strten(3)
        ext_press=-ext_press/3.d0*HaBohr3_GPa
        write(fn5,'(i5.5)') i
        filename = 'poslow'//fn5
        filenameout=trim(filename)//".rotated.ascii"
        call write_atomic_file_ascii(parini,filenameout,parini%nat,units,xred,latvec,fcart,strten,&
             &parini%char_type,parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,energy,parini%target_pressure_habohr,ext_press,enthalpy)
        open(unit=2,file="Enthalpies_rotated",access="append")
        write(2,'(a,5(1x,es25.15),i5)')trim(filename),parres%target_pressure_gpa,&
             &enthalpy*Ha_ev,energy*Ha_ev,ext_press,vol*Bohr_Ang**3,spgint
        close(2)
enddo
write(*,'(a)') " # Finished rotating structures. Output in _DIR folders and in the file 'Enthalpies'"
write(*,'(a)') " # Have a nice day!"

end subroutine
