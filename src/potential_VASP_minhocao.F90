module interface_vasp
  use global
  use defs_basis
  !use cell_utils

  implicit none

  private
  public :: &
    make_input_vasp,        &
    get_output_vasp,        &
    get_dos_vasp,           & 
    vasp_geopt,             &
    make_input_vasp_geopt,  &
    get_output_vasp_geopt   
 

contains

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
  subroutine make_input_vasp(parini,latvec, xred, iprec, ka, kb, kc, getwfk, dos)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
    real(8), intent(in) :: latvec(3,3)
    real(8), intent(in) :: xred(3,parini%nat)
    integer, intent(inout) :: ka, kb, kc
    logical, intent(in) :: getwfk
    logical, intent(in), optional :: dos
 

    real(8):: dproj(6), acell(3), rprim(3,3), angbohr, dkpt
    integer:: iat, iprec, itype
    integer:: nat_type(parini%ntypat_global)
    character(1):: fn

    if(iprec == 1) then
      dkpt = parini%dkpt1
    else
      dkpt = parini%dkpt2
    endif

    call system("cp  vasprun.xml vasprun.xml.bak")
    call system("rm -f INCAR POSCAR KPOINTS vasprun.xml")
    write(fn,'(i1.1)') iprec
    call system("cp -f INCAR."//fn//" INCAR")

    open(unit=87, file="INCAR", ACCESS="APPEND")
    if(.not.getwfk) then
      call system("./rm_wavecar.sh")
      write(87,'(a)') "LWAVE  = .FALSE."
      write(87,'(a)') "LCHARG = .FALSE."
    endif
    !Setup for only one force call
!    write(87,'(a,es25.15)') "PSTRESS = ",target_pressure_gpa*10.d0
    write(87,'(a)') "NSW    = 0"
    write(87,'(a)') "IBRION = 2"
    if(((all(parini%fixlat(1:6))).and.(.not.parini%fixlat(7))).or.parini%bc==2) then
        write(87,'(a)') "ISIF   = 0"
    else
        write(87,'(a)') "ISIF   = 2"
    endif
    close(87)
    
    if(present(dos)) then
      write(87,'(a)') "NEDOS = 2000"
    end if

    !Kpoint mesh
    open(unit=87,file="KPOINTS")
    write(87,'(a,i5)') "# Definition of the k-point mesh ",parini%vasp_kpt_mode
    write(87,'(i5)') 0
!    if(dkpt==0.d0.or.vasp_kpt_mode==2) then
    if(parini%vasp_kpt_mode==2) then
      write(87,'(a)') "Gamma"!"Monkhorst Pack"
      write(87,'(3(1x,i5),a)') ka, kb, kc,"  # Number of gridpoints in each dimension"
      write(87,'(3(1x,i5),a)') 0,0,0,"  # Shifts"
    elseif(parini%vasp_kpt_mode==1) then
      write(87,'(a)') "Auto"
      write(87,'(i5,a)') dkpt," # K-mesh length"
    else
      stop "Wrong kpt option"
    endif
!    write(*,'(a,3(1x,i5))') " # KPT mesh set up as follows: ", ka, kb, kc

    close(87)

    ! The atomic positions and the unit cell is defined here
    ! The data are read from a file, called posinp.ascii
    ! The units are expected to be in angstroem
    ! This particular part is made for Silicon and Hydrogen cells
    ! In VASP, the atomic types need to follow in groups, and each number of
    ! a specific type must be given
    nat_type=0
    do itype=1,parini%ntypat_global
      do iat=1,parini%nat
        if(parini%typat_global(iat)==itype) nat_type(itype) = nat_type(itype) + 1
      enddo
    enddo

    open(unit=87,file="POSCAR")
    write(87,'(a,i5)') "Accuracy level ",iprec 
    write(87,'(es15.7)') 1.d0
    write(87,*) latvec(:,1)*Bohr_Ang
    write(87,*) latvec(:,2)*Bohr_Ang
    write(87,*) latvec(:,3)*Bohr_Ang
    write(87,*) nat_type(:)
    write(87,'(a)') "Direct"
    do iat=1,parini%nat
      write(87,'(3(1x,es25.15),i5)') xred(:,iat)
    enddo
    close(87)

  end subroutine make_input_vasp


  ! Only implemented for closed shell!!!
  subroutine get_dos_vasp(fdos,efermi)
    character(150):: all_line
    integer:: n, k, l, m, i
    real(8):: fdos, density, efermi, energy, diff, Emax, Emin, number_of_steps, t

    diff = 1.d10
    open(unit=32, file="DOSCAR")
    ! first line is the total number of atoms and 3 numbers i don't know what they mean
    read(32,*) i, i, i, i
    ! line 2 to 5 are useless
    do i = 2, 5
      read (32,*)
    end do
    read(32,*) Emax, Emin, number_of_steps, efermi, t

    do
      read(32,'(a150)', end=99) all_line
      read(all_line,*) energy, density
      if(abs(efermi-energy).lt.diff) then
        fdos = density 
        diff = abs(efermi-energy)
      endif
    enddo

99  continue
    close(32)

  end subroutine get_dos_vasp


  subroutine get_output_vasp(parini,fcart,energy,strten)
  !use global, only: nat,target_pressure_gpa
  !use defs_basis
  !Since its a single call, we only have forces and stresses from one configuration!
  use mod_parini, only: typ_parini
  implicit none
  type(typ_parini), intent(in):: parini
  integer:: io,i,iat,n,k,l,m,i_tmp
  real(8):: fcart(3,parini%nat),energy,strten(6),value,latvec(3,3),xred(3,parini%nat),str_matrix(3,3),vol,a(3,3),scaling
  character(11):: ch_tmp
  character(150)::all_line
  logical:: vasp_5
character(40):: filename,units_tmp
real(8):: printval1,printval2
logical:: readfix,readfrag
logical:: fixat_tmp(parini%nat),fixlat_tmp(7)
  !Set to true uf you are using vasp version 5.x
  vasp_5=.true.
  
  energy=1.d10
  fcart=1.d10
  strten=1.d10
  
  ch_tmp="old"
  open(unit=32,file="vasprun.xml")
  do while(.true.)
   read(32,'(a150)',end=99)all_line
  !!write(*,*) all_line
   n = len_trim(all_line)
  
   k = index(all_line(1:n),"rec_basis")
    if(k.ne.0) then
      read(32,*)ch_tmp
      read(32,*)ch_tmp
      read(32,*)ch_tmp
      cycle
    endif
   k = index(all_line(1:n),"volumeweight")
    if(k.ne.0) then
      cycle
    endif
   
   k = index(all_line(1:n),"basis")
    if(k.ne.0) then
    !write(*,*) "Latvec found"
      read(32,'(a150)',end=99)all_line
      m = len_trim(all_line)
      l = scan(all_line(1:m),"/",.true.)
      read(all_line(1:l-2),*) ch_tmp,latvec(:,1)
      read(32,'(a150)',end=99)all_line
      m = len_trim(all_line)
      l = scan(all_line(1:m),"/",.true.)
      read(all_line(1:l-2),*) ch_tmp,latvec(:,2)
      read(32,'(a150)',end=99)all_line
      m = len_trim(all_line)
      l = scan(all_line(1:m),"/",.true.)
      read(all_line(1:l-2),*) ch_tmp,latvec(:,3)
    !write(*,*) latvec
    cycle
    endif
  
   k = index(all_line(1:n),"volume")
    if(k.ne.0) then
    l = VERIFY(all_line(1:n)," ",.false.)
    read(all_line(l:n),'(a17,F16.8,a)')ch_tmp,vol,ch_tmp
    endif
  
!   k = index(all_line(1:n),"e_fr_energy")
   k = index(all_line(1:n),"e_wo_entrp")
    if(k.ne.0) then
    l = VERIFY(all_line(1:n)," ",.false.)
    !write(*,*) "ETOT found"
      read(all_line(l:n),'(a22,F16.8,a)')ch_tmp,energy,ch_tmp
    !write(*,*) energy
    endif
  
   k = index(all_line(1:n),"forces")
    if(k.ne.0) then
    !write(*,*) "Forces found"
      do iat=1,parini%nat
      read(32,'(a150)',end=99)all_line
      m = len_trim(all_line)
      l = scan(all_line(1:m),"/",.true.)
      read(all_line(1:l-2),*) ch_tmp,fcart(:,iat)
      enddo
    !write(*,*) xred
    endif
   k = index(all_line(1:n),"stress")
    if(k.ne.0) then
    !write(*,*) "stress found"
      read(32,'(a150)',end=99)all_line
      m = len_trim(all_line)
      l = scan(all_line(1:m),"/",.true.)
      read(all_line(1:l-2),*) ch_tmp,str_matrix(:,1)
      read(32,'(a150)',end=99)all_line
      m = len_trim(all_line)
      l = scan(all_line(1:m),"/",.true.)
      read(all_line(1:l-2),*) ch_tmp,str_matrix(:,2)
      read(32,'(a150)',end=99)all_line
      m = len_trim(all_line)
      l = scan(all_line(1:m),"/",.true.)
      read(all_line(1:l-2),*) ch_tmp,str_matrix(:,3)
       strten(1)=-str_matrix(1,1)
       strten(2)=-str_matrix(2,2)
       strten(3)=-str_matrix(3,3)
       strten(6)=-str_matrix(1,2)
       strten(5)=-str_matrix(1,3)
       strten(4)=-str_matrix(2,3)
    !write(*,*) strten
    endif
  enddo
  
  99 continue 
  close(32)
  if(((all(parini%fixlat(1:6))).and.(.not.parini%fixlat(7))).or.parini%bc==2) strten=0.d0
  if(energy==1.d10.or.strten(1)==1.d10.or.fcart(1,1)==1.d10) stop "Could not find all requested variables"
  
  if(parini%target_pressure_gpa.ne.0.d0) then
  !In the vasprun.xml file at the end you will have the enthalpy instead of the total energy in the file, so 
  !we need to transform it back, remember pressures are in kilobar in vasp
  !        energy=energy-target_pressure_gpa*10.d0/1.60217733d-19/1.d22*vol
!          energy=energy-target_pressure_gpa*10.d0/(1.60217733d3)*vol
    !write(*,*) "ETOT found"
    !write(*,*) energy
  endif
  
  !Since in CONTCAR the cell and atomic positions are written with higher accuracy, get it from there:
  filename="CONTCAR"
  units_tmp="angstrom"
  readfix=.false.
  readfrag=.false.  
  call read_atomic_file_poscar(filename,parini%nat,units_tmp,xred,latvec,fcart,strten,&
           &fixat_tmp,fixlat_tmp,readfix,parini%fragarr,readfrag,printval1,printval2)
  latvec=latvec*Bohr_Ang !Internally already converted
!  open(unit=32,file="CONTCAR")
!  read(32,*)ch_tmp
!  read(32,*)scaling
!  read(32,*) latvec(:,1)
!  read(32,*) latvec(:,2)
!  read(32,*) latvec(:,3)
!  latvec=latvec*scaling
!  if(vasp_5) read(32,*)ch_tmp
!  read(32,*)i_tmp
!  read(32,*)ch_tmp
!      do iat=1,nat
!        read(32,*)xred(:,iat)
!      enddo
!  close(32)
  
  !Transform all to bohr
  latvec=latvec/Bohr_Ang
  energy=energy/Ha_eV
  strten=strten*0.1d0/HaBohr3_GPa
  fcart=fcart/Ha_eV*Bohr_Ang
  end subroutine
  
  subroutine vasp_geopt(parini,latvec,xred,fcart,strten,energy,iprec,ka,kb,kc,counter)
  !This routine will setup the input file for a vasp geometry optimization
  !It will also call the run script and harvest the output
  !use global, only: nat
  use mod_parini, only: typ_parini
  implicit none
  type(typ_parini), intent(in):: parini
  real(8):: xred(3,parini%nat),fcart(3,parini%nat),strten(6),energy,counter,tmp
  real(8):: dproj(6),acell(3),rprim(3,3),latvec(3,3)
  integer:: iat,iprec,ka,kb,kc,itype
  logical:: getwfk
  character(4):: tmp_char
  getwfk=.false.
  !Set up the input file to perform geometry optimization
   call make_input_vasp_geopt(parini,latvec,xred,iprec,ka,kb,kc,getwfk)
  ! call system("sleep 1")
  !Run the job NOW!
   call system("./runjob_geovasp.sh")
  ! call system("sleep 1")
  !Now harvest the structure, energy, forces, etc
   call get_output_vasp_geopt(parini,latvec,xred,fcart,energy,strten)
  !Check how many iterations have been needed
!   call system("grep Conjugate OUTCAR_geo_a |wc -l>tmp_count") 
!   call system("grep Conjugate OUTCAR_geo_b |wc -l>>tmp_count") 
!   call system("grep Conjugate OUTCAR_geo_c |wc -l>>tmp_count") 
   call system("rm -f tmp_count")
   call system("grep ' 1)  ---------------------------------------' OUTCAR_geo_a |wc -l>tmp_count") 
   call system("grep ' 1)  ---------------------------------------' OUTCAR_geo_b |wc -l>>tmp_count") 
   call system("grep ' 1)  ---------------------------------------' OUTCAR_geo_c |wc -l>>tmp_count") 
  ! call system("sleep 1")
   open(unit=32,file="tmp_count")
   read(32,*) tmp
   counter=counter+tmp
   read(32,*) tmp
   counter=counter+tmp
   read(32,*) tmp
   counter=counter+tmp
   close(32)
  !We create a backup of the geometry optimization file
   call system("cp vasprun.xml vasprun.xml.bak") 
  end subroutine
  
  subroutine make_input_vasp_geopt(parini,latvec,xred,iprec,ka,kb,kc,getwfk)
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
  !use global, only: nat,ntypat,znucl,typat,dkpt1,dkpt2,char_type,ntime_geopt,tolmxf,target_pressure_gpa,vasp_kpt_mode
  !use defs_basis,only: Bohr_Ang
  use mod_parini, only: typ_parini
  implicit none
  type(typ_parini), intent(in):: parini
  real(8):: xred(3,parini%nat)
  real(8):: dproj(6),acell(3),rprim(3,3),latvec(3,3),dkpt,angbohr
  real(8):: HaBohr_eVAng
  integer:: iat,iprec,ka,kb,kc,itype,nat_type(parini%ntypat_global)
  logical:: getwfk
  character(1):: fn
  character(150):: command
  angbohr=1.d0/Bohr_Ang
  
  getwfk=.false.
  
  if(iprec==1) then
  dkpt=parini%dkpt1
  else
  dkpt=parini%dkpt2
  endif
  
  !Rescale tolmxf to units of eV/Ang from Ha/Bohr
  HaBohr_eVAng=Ha_eV/Bohr_Ang

  
  
  command= "cp  vasprun.xml vasprun.xml.bak"
  call system(command)
  command= "rm -f INCAR POSCAR KPOINTS vasprun.xml OUTCAR_geo_a* OUTCAR_geo_b* OUTCAR_geo_c*"
  call system(command)
  write(fn,'(i1.1)') iprec
  command = "cp -f INCAR."//fn//".a INCAR_geo_a"
  call system(command)
  command = "cp -f INCAR."//fn//".b INCAR_geo_b"
  call system(command)
  command = "cp -f INCAR."//fn//" INCAR_geo_c"
  call system(command)
  
  open(unit=87,file="INCAR_geo_a",ACCESS="APPEND")
  !Setup for only one force call
  write(87,'(a)') ""
  !Setup for only a sequence of geopt
  write(87,'(a,i5)') "NSW = ",int(parini%paropt_geopt%nit*0.75d0)
  write(87,'(a,es25.15)') "PSTRESS = ",parini%target_pressure_gpa*10.d0
  write(87,'(a,es25.15)') "EDIFFG = ",-parini%paropt_geopt%fmaxtol*8.d0*HaBohr_eVAng
  !write(87,'(a)') "IBRION = 2"
  if(((all(parini%fixlat(1:6))).and.(.not.parini%fixlat(7))).or.parini%bc==2) then
     write(87,'(a)') "ISIF   = 0"
  else
     write(87,'(a)') "ISIF   = 3"
  endif 
  close(87)
  
  close(87)
  
  open(unit=87,file="INCAR_geo_b",ACCESS="APPEND")
  !Setup for only one force call
  write(87,'(a)') ""
  !Setup for only a sequence of geopt
  write(87,'(a,i5)') "NSW = ",int(parini%paropt_geopt%nit*0.25d0)
  write(87,'(a,es25.15)') "PSTRESS = ",parini%target_pressure_gpa*10.d0
  write(87,'(a,es25.15)') "EDIFFG = ",-parini%paropt_geopt%fmaxtol*HaBohr_eVAng
  !write(87,'(a)') "IBRION = 2"
  if(((all(parini%fixlat(1:6))).and.(.not.parini%fixlat(7))).or.parini%bc==2) then
     write(87,'(a)') "ISIF   = 0"
  else
     write(87,'(a)') "ISIF   = 3"
  endif 
  close(87)
  
  open(unit=87,file="INCAR_geo_c",ACCESS="APPEND")
  !Setup for only one force call
  write(87,'(a)') ""
  !Setup for only a sequence of geopt
  write(87,'(a,i5)') "NSW = ",0
!  write(87,'(a,es25.15)') "PSTRESS = ",target_pressure_gpa*10.d0
  write(87,'(a)') "IBRION = 2"
  if(((all(parini%fixlat(1:6))).and.(.not.parini%fixlat(7))).or.parini%bc==2) then
     write(87,'(a)') "ISIF   = 0"
  else
     write(87,'(a)') "ISIF   = 2"
  endif 
  close(87)
  
  !Kpoint mesh
  open(unit=87,file="KPOINTS")
  write(87,'(a,i5)') "# Definition of the k-point mesh ",parini%vasp_kpt_mode
  write(87,'(i5)') 0
  if(dkpt==0.d0) then
  write(87,'(a)') "Gamma"!"Monkhorst Pack"
  write(87,'(3(1x,i5),a)') ka,kb,kc,"  # Number of gridpoints in each dimension"
  write(87,'(3(1x,i5),a)') 0,0,0,"  # Shifts"
  elseif(parini%vasp_kpt_mode==2) then
  call find_kpt(ka,kb,kc,latvec,dkpt)
  write(87,'(a)') "Gamma"!"Monkhorst Pack"
  write(87,'(3(1x,i5),a)') ka,kb,kc,"  # Number of gridpoints in each dimension"
  write(87,'(3(1x,i5),a)') 0,0,0,"  # Shifts"
  elseif(parini%vasp_kpt_mode==1) then
  write(87,'(a)') "Auto"
  write(87,'(i5,a)') dkpt," # K-mesh length"
  else
  stop "Wrong kpt option"
  endif
  write(*,'(a,3(1x,i5))') " # KPT mesh set up as follows: ",ka,kb,kc
  close(87)
  
  
  !The atomic positions and the unit cell is defined here
  !The data are read from a file, called posinp.ascii
  !The units are expected to be in angstroem
  !This particular part is made for Silicon and Hydrogen cells
  !In VASP, the atomic types need to follow in groups, and each number of
  !a specific type must be given
  nat_type=0
  do itype=1,parini%ntypat_global
  do iat=1,parini%nat
  if(parini%typat_global(iat)==itype) nat_type(itype)=nat_type(itype)+1
  enddo
  enddo
  
  
  open(unit=87,file="POSCAR")
  write(87,'(a)') "For geopt initially"
  write(87,'(es15.7)') 1.d0
  write(87,*)latvec(:,1)/angbohr 
  write(87,*)latvec(:,2)/angbohr 
  write(87,*)latvec(:,3)/angbohr 
  write(87,*) nat_type(:) 
  if(any(parini%fixat(:))) write(87,'(a)') "Selective dynamics"
  write(87,'(a)') "Direct"
  do iat=1,parini%nat
  if(any(parini%fixat(:))) then
    if(parini%fixat(iat)) then
         write(87,'(3(1x,es25.15),a)') xred(:,iat), " F F F "  
    else
         write(87,'(3(1x,es25.15),a)') xred(:,iat), " T T T "  
    endif 
  else
    write(87,'(3(1x,es25.15),i5)') xred(:,iat)
  endif
  enddo 
  close(87)
  end  subroutine
  
  subroutine get_output_vasp_geopt(parini,latvec,xred,fcart,energy,strten)
  !use global, only: nat,target_pressure_gpa
  !use defs_basis
  !Since its a single call, we only have forces and stresses from one configuration!
  use mod_parini, only: typ_parini
  implicit none
  type(typ_parini), intent(in):: parini
  integer:: io,i,iat,n,k,l,m
  real(8):: fcart(3,parini%nat),energy,strten(6),value,latvec(3,3),xred(3,parini%nat),str_matrix(3,3),vol,a(3,3),scaling
  character(11):: ch_tmp
  character(150)::all_line
  logical:: vasp_5
character(40):: filename,units_tmp
real(8):: printval1,printval2
logical:: readfix,readfrag
logical:: fixat_tmp(parini%nat),fixlat_tmp(7)
  !if vasp is version 5.x, use vasp_5=.true.
  
  vasp_5=.true.
  
  energy=1.d10
  fcart=1.d10
  strten=1.d10
  
  ch_tmp="old"
  open(unit=32,file="vasprun.xml")
  do while(.true.)
   read(32,'(a150)',end=99)all_line
  !!write(*,*) all_line
   n = len_trim(all_line)
  
   k = index(all_line(1:n),"rec_basis")
    if(k.ne.0) then
      read(32,*)ch_tmp
      read(32,*)ch_tmp
      read(32,*)ch_tmp
      cycle
    endif
   k = index(all_line(1:n),"volumeweight")
    if(k.ne.0) then
      cycle
    endif
   
   k = index(all_line(1:n),"basis")
    if(k.ne.0) then
    !write(*,*) "Latvec found"
    !write(*,*) all_line(1:n)
      read(32,'(a150)',end=99)all_line
      m = len_trim(all_line)
      l = scan(all_line(1:m),"/",.true.)
      read(all_line(1:l-2),*) ch_tmp,latvec(:,1)
      read(32,'(a150)',end=99)all_line
      m = len_trim(all_line)
      l = scan(all_line(1:m),"/",.true.)
      read(all_line(1:l-2),*) ch_tmp,latvec(:,2)
      read(32,'(a150)',end=99)all_line
      m = len_trim(all_line)
      l = scan(all_line(1:m),"/",.true.)
      read(all_line(1:l-2),*) ch_tmp,latvec(:,3)
    !write(*,*) latvec
    cycle
    endif
  
  
   k = index(all_line(1:n),"volume")
    if(k.ne.0) then
    l = VERIFY(all_line(1:n)," ",.false.)
    !write(*,*) "Volume found" 
    !write(*,*) all_line(l:n)
    read(all_line(l:n),'(a17,F16.8,a)')ch_tmp,vol,ch_tmp
    !write(*,*) vol
    endif
  
!   k = index(all_line(1:n),"e_fr_energy")
   k = index(all_line(1:n),"e_wo_entrp")
    if(k.ne.0) then
    l = VERIFY(all_line(1:n)," ",.false.)
    !write(*,*) "ETOT found"
      read(all_line(l:n),'(a22,F16.8,a)')ch_tmp,energy,ch_tmp
    !write(*,*) energy
    endif
  
   k = index(all_line(1:n),"forces")
    if(k.ne.0) then
    !write(*,*) "Forces found"
      do iat=1,parini%nat
      read(32,'(a150)',end=99)all_line
      m = len_trim(all_line)
      l = scan(all_line(1:m),"/",.true.)
      read(all_line(1:l-2),*) ch_tmp,fcart(:,iat)
      enddo
    !write(*,*) xred
    endif
   k = index(all_line(1:n),"stress")
    if(k.ne.0) then
    !write(*,*) "stress found"
      read(32,'(a150)',end=99)all_line
      m = len_trim(all_line)
      l = scan(all_line(1:m),"/",.true.)
      read(all_line(1:l-2),*) ch_tmp,str_matrix(:,1)
      read(32,'(a150)',end=99)all_line
      m = len_trim(all_line)
      l = scan(all_line(1:m),"/",.true.)
      read(all_line(1:l-2),*) ch_tmp,str_matrix(:,2)
      read(32,'(a150)',end=99)all_line
      m = len_trim(all_line)
      l = scan(all_line(1:m),"/",.true.)
      read(all_line(1:l-2),*) ch_tmp,str_matrix(:,3)
       strten(1)=-str_matrix(1,1)
       strten(2)=-str_matrix(2,2)
       strten(3)=-str_matrix(3,3)
       strten(6)=-str_matrix(1,2)
       strten(5)=-str_matrix(1,3)
       strten(4)=-str_matrix(2,3)
    !write(*,*) strten
    endif
  enddo
  
  99 continue 
  close(32)
  if(((all(parini%fixlat(1:6))).and.(.not.parini%fixlat(7))).or.parini%bc==2) strten=0.d0
  if(energy==1.d10.or.strten(1)==1.d10.or.fcart(1,1)==1.d10) stop "Could not find all requested variables"
  
  if(parini%target_pressure_gpa.ne.0.d0) then
  !In the vasprun.xml file at the end you will have the enthalpy instead of the total energy in the file, so 
  !we need to transform it back, remember pressures are in kilobar in vasp
  !        energy=energy-target_pressure_gpa*10.d0/1.60217733d-19/1.d22*vol
  !        energy=energy-target_pressure_gpa*10.d0/1.60217733d3*vol
    !write(*,*) "ETOT found"
    !write(*,*) energy
  endif
  
  !Since in CONTCAR the cell and atomic positions are written with higher accuracy, get it from there:
  filename="CONTCAR"
  units_tmp="angstrom"
  readfix=.false.
  readfrag=.false.  
  call read_atomic_file_poscar(filename,parini%nat,units_tmp,xred,latvec,fcart,strten,&
           &fixat_tmp,fixlat_tmp,readfix,parini%fragarr,readfrag,printval1,printval2)
  latvec=latvec*Bohr_Ang !Internally already converted
!  open(unit=32,file="CONTCAR")
!  read(32,*)ch_tmp
!  read(32,*)scaling
!  read(32,*) latvec(:,1)
!  read(32,*) latvec(:,2)
!  read(32,*) latvec(:,3)
!  latvec=latvec*scaling
!  read(32,*)ch_tmp
!  if(vasp_5) read(32,*)ch_tmp
!  read(32,*)ch_tmp
!      do iat=1,nat
!        read(32,*)xred(:,iat)
!      enddo
!  close(32)
  
  !Transform all to bohr
  latvec=latvec/Bohr_Ang
  energy=energy/Ha_eV
  strten=strten*0.1d0/HaBohr3_GPa
  fcart=fcart/Ha_eV*Bohr_Ang
  end subroutine
  end module interface_vasp
