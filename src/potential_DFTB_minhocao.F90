module interface_dftb
  use global
  use defs_basis

  implicit none

  private
  public :: &
 make_input_dftb,                     &
 get_output_dftb,                     &
 get_dos_dftb,                        &
 dftb_geopt,                          &
 make_input_dftb_geopt,               &
 get_output_dftb_geopt

contains

  subroutine make_input_dftb(parini,latvec,xred0,iprec,ka,kb,kc,getwfk,dos)
  !This routine will append some informations to a file already containing some informations about the abininit runs
  !The informations appended are:
  !-The atomic informations
  !A file KPOINTS is generated. There are two options available:
  !Automated generation if vasp_kpt_mode==1
  !Monkhorst pack if vasp_kpt_mode==2 
  !ATTENTION:
  !The meaning of dkpt1 and dkpt2 is different depending on vasp_kpt_mode:
  !accuracy is given by the integer length of dkpt for vasp_kpt_mode==1 (10 for insulators, 100 for metals)
  !accuracy is 2pi/bohr*dkpt for vasp_kpt_mode==2 
  use mod_parini, only: typ_parini
  use defs_basis, only: Bohr_Ang
  implicit none
  type(typ_parini), intent(in):: parini
  logical, intent(in), optional :: dos
  real(8):: xred(3,parini%nat),xcart(3,parini%nat),xred0(3,parini%nat)
  real(8):: dproj(6),acell(3),rprim(3,3),latvec(3,3),dkpt,angbohr
  integer:: iat,iprec,ka,kb,kc,itype
  integer:: nat_type(parini%ntypat_global)
  logical:: getwfk
  character(1):: fn
  character(150):: command,all_atoms
  angbohr=1.d0/Bohr_Ang
  xred=xred0
  call backtocell(parini%nat,latvec,xred)
  !getwfk=.false.
  
  if(iprec==1) then
  dkpt=parini%dkpt1
  else
  dkpt=parini%dkpt2
  endif
  
  write(fn,'(i1.1)') iprec
  command = "cp -f dftb_in."//fn//".hsd dftb_in.hsd"
  call system(command)
  
  open(unit=87,file="input_restart.gen")
  if(getwfk) then
  !write(87,'(a)') "# Use previous charge file "
  !write(87,'(a)') " ReadInitialCharges = Yes"
  write(87,'(a)') "Yes"
  else
  !write(87,'(a)') "# Dont use previous charge file "
  !write(87,'(a)') " ReadInitialCharges = No"
  write(87,'(a)') "No"
  endif
  close(87)
  
  !Kpoint mesh
  open(unit=87,file="input_kpoints.gen")
  write(87,'(3(1x,i5),a)') ka,0,0," # a-direction"
  write(87,'(3(1x,i5),a)') 0,kb,0," # b-direction"
  write(87,'(3(1x,i5),a)') 0,0,kc," # c-direction"
  write(87,'(3(1x,i5),a)') 0,0,0,"  # Shifts"
  close(87)
  
  !The Geometry section
  !write(87,'(a)') "#Geometry section: lattice and atoms"
  open(unit=87,file="input_geometry.gen")
    if(parini%bc==2) then
      write(87,'(i5,a)') parini%nat, " C"
    else
      write(87,'(i5,a)') parini%nat, " F"
    endif
      write(87,*) (parini%char_type(iat)(1:2)//" ", iat=1,parini%ntypat_global)
    if(parini%bc==2) then
      call rxyz_int2cart(latvec,xred,xcart,parini%nat)
      do iat = 1, parini%nat
        write(87,'(i5,1x,i5,3(1x,es25.15))') iat, parini%typat_global(iat), xcart(:, iat)*Bohr_Ang
      end do
    else
      do iat = 1, parini%nat
        write(87,'(i5,1x,i5,3(1x,es25.15))') iat, parini%typat_global(iat), xred(:, iat)
      end do
      write(87,'(3(1x,es25.15))') 0.d0,0.d0,0.d0
      write(87,'(3(1x,es25.15))') latvec(:, 1)*Bohr_Ang
      write(87,'(3(1x,es25.15))') latvec(:, 2)*Bohr_Ang
      write(87,'(3(1x,es25.15))') latvec(:, 3)*Bohr_Ang
    endif
  close(87)
  end  subroutine
  
  subroutine get_output_dftb(parini,fcart,energy,strten)
  use mod_parini, only: typ_parini
  use defs_basis
  !Since its a single call, we only have forces and stresses from one configuration!
  implicit none
  type(typ_parini), intent(in):: parini
  integer:: io,i,iat,n,k,l,m,int_tmp,nat_cell
  real(8):: fcart(3,parini%nat),energy,strten(6),value,latvec(3,3),xred(3,parini%nat),str_matrix(3,3),vol,a(3,3),scaling
  character(11):: ch_tmp
  character(250)::all_line
  
  energy=1.d10
  fcart=1.d10
  strten=1.d10
  nat_cell=0
  
  ch_tmp="old"
  open(unit=32,file="detailed.out")
  do while(.true.)
   read(32,'(a250)',end=99)all_line
  !!write(*,*) all_line
   n = len_trim(all_line)
  
   k = index(all_line(1:n),"Total Mermin free energy")
    if(k.ne.0) then
      m = len_trim(all_line)
      l = scan(all_line(1:m),":",.true.)
      read(all_line(l+1:m),*) energy
      cycle
    endif
   k = index(all_line(1:n),"Total Forces")
    if(k.ne.0) then
      do iat=1,parini%nat
      read(32,*) fcart(:,iat) 
      enddo
      cycle
    endif
   k = index(all_line(1:n),"Total stress tensor")
    if(k.ne.0) then
      do iat=1,3
      read(32,*) str_matrix(:,iat)
      enddo
       strten(1)=-str_matrix(1,1)
       strten(2)=-str_matrix(2,2)
       strten(3)=-str_matrix(3,3)
       strten(6)=-str_matrix(1,2)
       strten(5)=-str_matrix(1,3)
       strten(4)=-str_matrix(2,3)
      cycle
    endif
  enddo
  
  99 continue 
  close(32)
  if(parini%bc==2) strten=0.d0
  if(energy==1.d10.or.strten(1)==1.d10.or.fcart(1,1)==1.d10) stop "Could not find all requested variables"
  !Transform all to bohr
  energy=energy!/real(nat_cell,8)*real(nat,8)/Ha_eV
  strten=strten!/HaBohr3_GPa
  !strten=0.d0
  fcart=fcart!/Ha_eV*Bohr_Ang
  end subroutine
  
  ! Only implemented for closed shell!!!
  subroutine get_dos_dftb(fdos,efermi)
    character(250):: all_line
    integer:: n, k, l, m, i
    real(8):: fdos, density, efermi, energy, diff, Emax, Emin, number_of_steps, t
    efermi=1.d10
    fdos=1.d10
!Get the efermi from the normal output of dftb 
    open(unit=32,file="detailed.out")
    do while(.true.)
     read(32,'(a250)',end=98)all_line
    !!write(*,*) all_line
     n = len_trim(all_line)
     k = index(all_line(1:n),"Fermi energy:")
      if(k.ne.0) then
        m = len_trim(all_line)
        l = scan(all_line(1:m),":",.true.)
        read(all_line(l+1:m),*) efermi
        exit
      endif
    enddo

98 continue 
    if(efermi==1.d10) stop "Fermi energy not found in detailed.out"  
    efermi=efermi*Ha_eV
    close(32)
!Run the python script to convert 
    call system("dp_dos band.out dos_total.dat")
    diff = 1.d10
    open(unit=32,file="dos_total.dat")
    do
      read(32,'(a250)', end=99) all_line
      read(all_line,*) energy, density
      if(abs(efermi-energy).lt.diff) then
        fdos = density 
        diff = abs(efermi-energy)
      endif
    enddo

99  continue
    close(32)
    call system("mv dos_total.dat dos_total.dat.bak")
    write(*,*) efermi, fdos
    if(fdos==1.d10) stop "Fermi level not found in dos_total.dat"  
  end subroutine get_dos_dftb

  
  subroutine dftb_geopt(parini,latvec,xred,fcart,strten,energy,iprec,ka,kb,kc,counter)
  !This routine will setup the input file for a vasp geometry optimization
  !It will also call the run script and harvest the output
  use mod_parini, only: typ_parini
  implicit none
  type(typ_parini), intent(in):: parini
  real(8):: xred(3,parini%nat),fcart(3,parini%nat),strten(6),energy,counter,tmp
  real(8):: dproj(6),acell(3),rprim(3,3),latvec(3,3),fmax
  integer:: iat,iprec,ka,kb,kc,itype
  logical:: getwfk
  character(4):: tmp_char
  getwfk=.false.
  !Set up the input file to perform geometry optimization
  !First run with low precision
  iprec=2
  call make_input_dftb_geopt(parini,latvec,xred,iprec,ka,kb,kc,getwfk)
  !Run the job NOW!
  write(*,'(a,i1)')' # GEOPT START NATIVE DFTB OPTIMIZER RUN 1, IPREC=',iprec
  call system("./runjob.sh")
  !Collect data
  call get_output_dftb_geopt(parini,latvec,xred,fcart,energy,strten,fmax)
  !Write intermediate data
  call system("grep 'Geometry step:' dftb.log|wc -l>tmp_count")
  open(unit=32,file="tmp_count")
  read(32,*) tmp
  counter=counter+tmp
  close(32)
  write(*,'(a,i5,a,es15.7)')' # GEOPT INTERMEDIATE 1 FINISHED: ITER=',int(tmp),", FMAX=",fmax
  !Now run with high precision
  iprec=1
  call make_input_dftb_geopt(parini,latvec,xred,iprec,ka,kb,kc,getwfk)
  !Run the job NOW!
  write(*,'(a,i1)')' # GEOPT START NATIVE DFTB OPTIMIZER RUN 2, IPREC=',iprec
  call system("./runjob.sh")
  !Collect data
  call get_output_dftb_geopt(parini,latvec,xred,fcart,energy,strten,fmax)
  !Write intermediate data
  call system("grep 'Geometry step:' dftb.log|wc -l>tmp_count")
  open(unit=32,file="tmp_count")
  read(32,*) tmp
  counter=counter+tmp
  close(32)
  write(*,'(a,i5,a,es15.7)')' # GEOPT INTERMEDIATE 2 FINISHED: ITER=',int(tmp),", FMAX=",fmax
  call make_input_dftb_geopt(parini,latvec,xred,iprec,ka,kb,kc,getwfk)
  !Run the job NOW!
  write(*,'(a,i1)')' # GEOPT START NATIVE DFTB OPTIMIZER FINAL, IPREC=',iprec
  call system("./runjob.sh")
  !Collect data
  call get_output_dftb_geopt(parini,latvec,xred,fcart,energy,strten,fmax)
  !Write intermediate data
  call system("grep 'Geometry step:' dftb.log|wc -l>tmp_count")
  open(unit=32,file="tmp_count")
  read(32,*) tmp
  counter=counter+tmp
  close(32)
  write(*,'(a,i5,a,es15.7)')' # GEOPT CONVERGED: ITER=',int(tmp),", FMAX=",fmax
  call system("rm -f tmp_count")
  end subroutine
  
  
  subroutine make_input_dftb_geopt(parini,latvec,xred,iprec,ka,kb,kc,getwfk)
  use mod_parini, only: typ_parini
  use defs_basis,only: Bohr_Ang
  implicit none
  type(typ_parini), intent(in):: parini
  real(8):: xred(3,parini%nat),xcart(3,parini%nat)
  real(8):: dproj(6),acell(3),rprim(3,3),latvec(3,3),dkpt,angbohr
  integer:: iat,iprec,ka,kb,kc,itype,nat_type(parini%ntypat_global)
  logical:: getwfk
  character(1):: fn
  character(150):: command,all_atoms
  angbohr=1.d0/Bohr_Ang
  
  getwfk=.false.
  
  if(iprec==1) then
  dkpt=parini%dkpt1
  else
  dkpt=parini%dkpt2
  endif
  
  write(fn,'(i1.1)') iprec
  command = "cp -f dftb_ingeo."//fn//".hsd dftb_in.hsd"
  call system(command)
  
  open(unit=87,file="input_restart.gen")
  if(getwfk) then
  !write(87,'(a)') "# Use previous charge file "
  !write(87,'(a)') " ReadInitialCharges = Yes"
  write(87,'(a)') "Yes"
  else
  !write(87,'(a)') "# Dont use previous charge file "
  !write(87,'(a)') " ReadInitialCharges = No"
  write(87,'(a)') "No"
  endif
  close(87)
  
  !Kpoint mesh
  open(unit=87,file="input_kpoints.gen")
  if(dkpt==0.d0) then
  write(87,'(3(1x,i5),a)') ka,0,0," # a-direction"
  write(87,'(3(1x,i5),a)') 0,kb,0," # b-direction"
  write(87,'(3(1x,i5),a)') 0,0,kc," # c-direction"
  write(87,'(3(1x,i5),a)') 0,0,0,"  # Shifts"
  else
  call find_kpt(ka,kb,kc,latvec,dkpt)
  write(87,'(3(1x,i5),a)') ka,0,0," # a-direction"
  write(87,'(3(1x,i5),a)') 0,kb,0," # b-direction"
  write(87,'(3(1x,i5),a)') 0,0,kc," # c-direction"
  write(87,'(3(1x,i5),a)') 0,0,0, " # Shifts"
  endif
  close(87)
  write(*,'(a,3(1x,i5))') " # KPT mesh set up as follows: ",ka,kb,kc
  
  !The Geometry section
  open(unit=87,file="input_geometry.gen")
    if(parini%bc==2) then
      write(87,'(i5,a)') parini%nat, " C"
    else
      write(87,'(i5,a)') parini%nat, " F"
    endif
      write(87,*) (parini%char_type(iat)(1:2)//" ", iat=1,parini%ntypat_global)
    if(parini%bc==2) then
      call rxyz_int2cart(latvec,xred,xcart,parini%nat)
      do iat = 1, parini%nat
        write(87,'(i5,1x,i5,3(1x,es25.15))') iat, parini%typat_global(iat), xcart(:, iat)*Bohr_Ang
      end do
    else
      do iat = 1, parini%nat
        write(87,'(i5,1x,i5,3(1x,es25.15))') iat, parini%typat_global(iat), xred(:, iat)
      end do
      write(87,'(3(1x,es25.15))') 0.d0,0.d0,0.d0
      write(87,'(3(1x,es25.15))') latvec(:, 1)*Bohr_Ang
      write(87,'(3(1x,es25.15))') latvec(:, 2)*Bohr_Ang
      write(87,'(3(1x,es25.15))') latvec(:, 3)*Bohr_Ang
    endif
  close(87)

  !The Geometry optimizer section
  open(unit=87,file="input_driver.gen")
  if(iprec==2) then
  write(87,'(a,es25.15)') " MaxForceComponent = ", parini%paropt_geopt%fmaxtol*10.d0
  else
  write(87,'(a,es25.15)') " MaxForceComponent = ", parini%paropt_geopt%fmaxtol
  endif
  write(87,'(a,i5)') " MaxSteps = ",parini%paropt_geopt%nit 
  if(((all(parini%fixlat(1:6))).and.(.not.parini%fixlat(7))).or.parini%bc==2) then
    write(87,'(a)') " LatticeOpt = No"
  else
    write(87,'(a)') " LatticeOpt = Yes"
  endif
  write(87,'(a,es25.15)') " Pressure = ", parini%target_pressure_habohr
  close(87)

  open(unit=87,file="dftb_in.hsd",access="append")
  write(87,'(a)') " Driver = ConjugateGradient{"
  if(iprec==2) then
  write(87,'(a,es25.15)') " MaxForceComponent = ", parini%paropt_geopt%fmaxtol*10.d0
  else
  write(87,'(a,es25.15)') " MaxForceComponent = ", parini%paropt_geopt%fmaxtol
  endif
  write(87,'(a,i5)') " MaxSteps = ",parini%paropt_geopt%nit 
  if(((all(parini%fixlat(1:6))).and.(.not.parini%fixlat(7))).or.parini%bc==2) then
    write(87,'(a)') " LatticeOpt = No"
  else
    write(87,'(a)') " LatticeOpt = Yes"
  endif
  if(parini%bc==1) write(87,'(a,es25.15)') " Pressure = ", parini%target_pressure_habohr
  if(any(parini%fixat(:))) then
     write(87,'(a)') "Constraints = {"
     do iat=1,parini%nat
       if(parini%fixat(iat)) then
          write(87,'(i5,1x,a)') iat, " 1.0, 0.0, 0.0 "
          write(87,'(i5,1x,a)') iat, " 0.0, 1.0, 0.0 "
          write(87,'(i5,1x,a)') iat, " 0.0, 0.0, 1.0 "
       endif
     enddo
     write(87,'(a)') "}"
  endif
  write(87,'(a)') " } "
  close(87)
 
  end  subroutine
  
  
  
  subroutine get_output_dftb_geopt(parini,latvec,xred,fcart,energy,strten,fmax)
  use mod_parini, only: typ_parini
  use defs_basis
  !Since its a single call, we only have forces and stresses from one configuration!
  implicit none
  type(typ_parini), intent(in):: parini
  integer:: io,i,iat,n,k,l,m,int_tmp,istr
  real(8):: fcart(3,parini%nat),energy,strten(6),value,latvec(3,3),xred(3,parini%nat),str_matrix(3,3),vol,a(3,3),scaling,r_tmp
  real(8):: fmax,strtarget(6),dstr(6),xtmp(3,parini%nat)
  character(11):: ch_tmp
  character(150)::all_line
  logical:: cartesian
  call get_output_dftb(parini,fcart,energy,strten)
  open(unit=87,file="geo_end.gen")
  read(87,*) int_tmp,ch_tmp
  if(trim(ch_tmp)=="C") then
    cartesian=.true.
  else
    cartesian=.false.
  endif
  read(87,*) ch_tmp
  do iat=1,parini%nat
    read(87,*) int_tmp,int_tmp,xtmp(:,iat)
  enddo
  if(cartesian) then
     xtmp=xtmp/Bohr_Ang
     call rxyz_cart2int(latvec,xred,xtmp,parini%nat)
  else
     xred=xtmp     
     read(87,*) r_tmp
     do iat=1,3
       read(87,*) latvec(:,iat)
     enddo
     latvec=latvec/Bohr_Ang
  endif

  close(87)
  !!Compute maximal component of forces, EXCLUDING any fixed components
   fmax=0.0d0
   do iat=1,parini%nat
     do i=1,3
  !     if (dtsets(1)%iatfix(i,iat) /= 1) then
         if( abs(fcart(i,iat)) >= fmax ) fmax=abs(fcart(i,iat))
  !     end if
     end do
   end do
   strtarget=0.d0
   strtarget(1:3)=-parini%target_pressure_habohr
   dstr(:)=strten(:)-strtarget(:)
  !Eventually take into account the stress
   do istr=1,6
       if(abs(dstr(istr))*parini%paropt_geopt%strfact >= fmax ) fmax=abs(dstr(istr))*parini%paropt_geopt%strfact
   end do
  end subroutine
  
  end module interface_dftb
  
