module interface_siesta
  use defs_basis
  !use cell_utils

  implicit none

  private
  public :: &
    make_input_siesta,      &
    get_output_siesta,      &
    siesta_geopt,           &
    make_input_siesta_geopt,&
    get_output_siesta_geopt

contains

  
!This routine will append some informations to a file already containing some informations about the abininit runs
!The informations appended are:
! - Read/dont read wavefunction from file is ignored for siesta
! - The kpoint mesh
! - The atomic informations
subroutine make_input_siesta(parini,latvec, xred, iprec, ka, kb, kc, getwfk, dos)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
    real(8), intent(in) :: latvec(3,3)
    real(8), intent(in) :: xred(3,parini%nat)
    integer, intent(inout) :: ka, kb, kc
    logical, intent(in) :: getwfk
    logical, intent(in), optional :: dos

    real(8) :: dproj(6), acell(3), rprim(3,3), dkpt
    integer :: iat, iprec, itype
    character(1):: fn
    if(iprec == 1) then
      dkpt = parini%dkpt1
    else
      dkpt = parini%dkpt2
    endif
    !For siesta, there are two modes for setting up the k-point mesh:
    !mode 1: setting up of the mesh done by seista, using kgrid_cutoff
    !        then dkptx is directly the cutoff length in Bohr
    !mode 2: setting up the monkhorst pack mesh by setting up the mesh according to 
    !        dkptx as the density of the mesh in the B-Zone

    call system("rm -f siestarun.fdf")
    call system("rm -f FORCE_STRESS.bak")
    call system("cp  FORCE_STRESS FORCE_STRESS.bak")
    write(fn,'(i1.1)') iprec
    call system("cp -f siestarun."//fn//".fdf siestarun.fdf")

    open(unit=87,file="siestarun.fdf",ACCESS="APPEND")

    !Setup for only one force call
    write(87,'(a)') "MD.TypeOfRun         cg"
    write(87,'(a)') "MD.NumCGsteps        0"

    !Kpoint mesh
    write(87,'(a)') "# Definition of the k-point mesh"

    if(parini%siesta_kpt_mode==1.and.dkpt.ne.0.d0) then
      write(87,'(a,es25.15,a)') "kgrid_cutoff    ", dkpt,"  Bohr"
    else
      write(*,'(a,3(1x,i5))') " # KPT mesh set up as follows: ",ka,kb,kc
      write(87,'(a)') "%block kgrid_Monkhorst_Pack"
      write(87,'(3(1x,i5),1x,f7.3)') ka, 0, 0, 0.5
      write(87,'(3(1x,i5),1x,f7.3)') 0, kb, 0, 0.5
      write(87,'(3(1x,i5),1x,f7.3)') 0, 0, kc, 0.5
      write(87,'(a)') "%endblock kgrid_Monkhorst_Pack"
    endif

    write(87,'(a)') " " 
    write(87,'(a)') " " 

    !The atomic positions and the unit cell is defined here
    !The data are read from a file, called posinp.ascii
    !The units are expected to be in angstroem
    !This particular part is made for Silicon and Hydrogen cells
    write(87,'(a)') "# Definition of the unit cell"
    write(87,*) "NumberOfSpecies", parini%ntypat_global !Number of atom types
    write(87,*) "NumberOfAtoms",      parini%nat !Number of atoms

    write(87,'(a)') "%block ChemicalSpeciesLabel"
    do itype = 1, parini%ntypat_global
      write(87,'(i5,1x,i5,1x,a2)') itype, int(parini%znucl(itype)), trim(parini%char_type(itype)) 
    enddo
    write(87,'(a)') "%endblock ChemicalSpeciesLabel"
    write(87,'(a,f7.3,a)') "LatticeConstant  ",1.d0,"  Bohr"
    write(87,'(a)') "%block LatticeVectors"
    write(87,'(3(1x,es25.15))') latvec(:,1) 
    write(87,'(3(1x,es25.15))') latvec(:,2) 
    write(87,'(3(1x,es25.15))') latvec(:,3) 
    write(87,'(a)') "%endblock LatticeVectors"
    write(87,*) " "
    write(87,*) " "
    write(87,*) "AtomicCoordinatesFormat Fractional" 
    write(87,*) "%block AtomicCoordinatesAndAtomicSpecies"
    do iat=1,parini%nat
      write(87,'(3(1x,es25.15),i5)') xred(:,iat),parini%typat_global(iat)
    enddo
    write(87,*) "%endblock AtomicCoordinatesAndAtomicSpecies"
    close(87)
  end subroutine make_input_siesta


  subroutine get_output_siesta(parini,fcart,energy,strten)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
    real(8), intent(out) ::  fcart(3, parini%nat), energy, strten(6)

    integer:: io,i,iat
    real(8):: value, str_matrix(3,3)
    character(11):: ch_tmp

    open(unit=32,file="FORCE_STRESS")
    !First value is energy in Ry
    read(32,*) value   
    energy = value*0.5d0  !Convert to Hartree
    !Now get the stress tensor in Ry/Bohr**3
    read(32,*) str_matrix(:,1) 
    read(32,*) str_matrix(:,2) 
    read(32,*) str_matrix(:,3) 

    strten(1) = str_matrix(1,1)
    strten(2) = str_matrix(2,2)
    strten(3) = str_matrix(3,3)
    strten(6) = str_matrix(1,2)
    strten(5) = str_matrix(1,3)
    strten(4) = str_matrix(2,3)
    strten = strten*0.5d0  !Convert to Hartree
    read(32,*) value !Temporary

    !Now get atomic forces
    do iat=1,parini%nat
      read(32,*) value, value, fcart(:,iat)
    enddo
    !Convert back to Ha/Bohr instead of Ry/Bohr
    fcart=fcart*0.5d0
    close(32)

  end subroutine get_output_siesta

  subroutine siesta_geopt(parini,latvec,xred,fcart,strten,energy,iprec,ka,kb,kc,counter)
  use mod_parini, only: typ_parini
  !This routine will setup the input file for a siesta geometry optimization
  !It will also call the run script and harvest the output
  implicit none
  type(typ_parini), intent(in):: parini
  real(8):: xred(3,parini%nat),fcart(3,parini%nat),strten(6),energy,counter
  real(8):: dproj(6),acell(3),rprim(3,3),latvec(3,3)
  integer:: iat,iprec,ka,kb,kc,itype
  logical:: getwfk
  character(4):: tmp_char
  getwfk=.false.
  !Set up the input file to perform geometry optimization
   call make_input_siesta_geopt(parini,latvec,xred,iprec,ka,kb,kc,getwfk)
  ! call system("sleep 1")
  !Run the job NOW!
   call system("./runjob.sh")
  ! call system("sleep 1")
  !Now harvest the structure, energy, forces, etc
   call get_output_siesta_geopt(parini,latvec,xred,fcart,energy,strten)
  !Check how many iterations have been needed
   call system("grep Begin siestarun.out|tail -1>tmp_count")
  ! call system("sleep 1")
   open(unit=32,file="tmp_count")
   read(32,*) tmp_char,tmp_char,tmp_char,tmp_char,counter
   close(32)
   call system("rm -f tmp_count")
  !We create a backup of the geometry optimization file
   call system("cp siestarun.out siestarun_geoptbak.out") 
  end subroutine
  
  
  subroutine make_input_siesta_geopt(parini,latvec,xred,iprec,ka,kb,kc,getwfk)
  use mod_parini, only: typ_parini
  !This routine will append some informations to a file already containing some informations about the siesta runs
  !The informations appended are:
  !-Geometry optimization parameters
  !-The kpoint mesh
  !-The atomic informations
  !use global, only: nat,ntypat,znucl,typat,dkpt1,dkpt2,char_type,ntime_geopt,tolmxf,target_pressure_gpa,siesta_kpt_mode
  implicit none
  type(typ_parini), intent(in):: parini
  real(8):: xred(3,parini%nat)
  real(8):: dproj(6),acell(3),rprim(3,3),latvec(3,3),dkpt
  integer:: iat,iprec,ka,kb,kc,itype
  logical:: getwfk
  character(1):: fn
  character(150):: command
  !For siesta, there are two modes for setting up the k-point mesh:
  !mode 1: setting up of the mesh done by seista, using kgrid_cutoff
  !        then dkptx is directly the cutoff length in Bohr
  !mode 2: setting up the monkhorst pack mesh by setting up the mesh according to 
  !        dkptx as the density of the mesh in the B-Zone
!  siesta_kpt_mode=2
  
  getwfk=.false.
  
  if(iprec==1) then
  dkpt=parini%dkpt1
  else
  dkpt=parini%dkpt2
  endif
  
  
  command= "rm -f siestarun.fdf"
  call system(command)
  command= "rm -f FORCE_STRESS.bak"
  call system(command)
  command= "mv  FORCE_STRESS FORCE_STRESS.bak"
  call system(command)
  command= "rm -f siestarun.1.STRUCT_OUT.bak"
  call system(command)
  command= "mv -f siestarun.1.STRUCT_OUT siestarun.1.STRUCT_OUT.bak"
  call system(command)
  write(fn,'(i1.1)') iprec
  command = "cp -f siestarun."//fn//".fdf siestarun.fdf"
  call system(command)
  
  open(unit=87,file="siestarun.fdf",ACCESS="APPEND")
  !Geometry optimization
  !Setup for only one force call
  write(87,'(a)')            "MD.TypeOfRun         CG"
  !write(87,'(a)')           "MD.TypeOfRun         FIRE"
  !write(87,'(a)')           "MD.TypeOfRun         Broyden"
  write(87,'(a,i5)')          "MD.NumCGsteps  ",parini%paropt_geopt%nit
  write(87,'(a)')            "MD.VariableCell .true."
  write(87,'(a,es25.15,a)')        "MD.MaxForceTol  ",  2.d0*parini%paropt_geopt%fmaxtol ,"  Ry/Bohr"
  write(87,'(a)')            "MD.MaxStressTol 0.1 GPa"
  write(87,'(a,es25.15,a)')  "MD.TargetPressure  ",parini%target_pressure_gpa, " GPa"
  
  !Kpoint mesh
  write(87,'(a)') "# Definition of the k-point mesh"
  
  if(dkpt==0.d0) then
  write(87,'(a,3(1x,i5),a)') "#",ka,kb,kc,"  # Number of gridpoints in each dimension"
  elseif(parini%siesta_kpt_mode==2) then
  call find_kpt(ka,kb,kc,latvec,dkpt)
  write(87,'(a,3(1x,i5),a)') "#",ka,kb,kc,"  # Number of gridpoints in each dimension"
  endif
  
  if(parini%siesta_kpt_mode==1.and.dkpt.ne.0.d0) then
  write(87,'(a,es25.15,a)') "kgrid_cutoff    ", dkpt,"  Bohr"
  else
  write(*,'(a,3(1x,i5))') " # KPT mesh set up as follows: ",ka,kb,kc
  write(87,'(a)') "%block kgrid_Monkhorst_Pack"
  write(87,'(3(1x,i5),1x,f7.3)') ka,0,0,0.5
  write(87,'(3(1x,i5),1x,f7.3)') 0,kb,0,0.5
  write(87,'(3(1x,i5),1x,f7.3)') 0,0,kc,0.5
  write(87,'(a)') "%endblock kgrid_Monkhorst_Pack"
  endif
  
  !**************************************************************************
  if(getwfk) write(87,'(a)') " getwfk -1           # This is to speed up the calculation, by restarting"
  if(getwfk) write(87,'(a)') "                     # from previous wavefunctions, transferred from the old "
  if(getwfk) write(87,'(a)') "                     # to the new k-points."
  write(87,'(a)') " " 
  write(87,'(a)') " " 
  
  
  !The atomic positions and the unit cell is defined here
  !The data are read from a file, called posinp.ascii
  !The units are expected to be in angstroem
  !This particular part is made for Silicon and Hydrogen cells
  write(87,'(a)') "# Definition of the unit cell"
  write(87,*) "NumberOfSpecies", parini%ntypat_global !Number of atom types
  write(87,*) "NumberOfAtoms",      parini%nat !Number of atoms
  
  write(87,'(a)') "%block ChemicalSpeciesLabel"
  do itype=1,parini%ntypat_global
     write(87,'(i5,1x,i5,1x,a2)') itype,int(parini%znucl(itype)),trim(parini%char_type(itype)) 
  enddo
  write(87,'(a)') "%endblock ChemicalSpeciesLabel"
  write(87,'(a,f7.3,a)') "LatticeConstant  ",1.d0,"  Bohr"
  write(87,'(a)') "%block LatticeVectors"
  write(87,'(3(1x,es25.15))') latvec(:,1) 
  write(87,'(3(1x,es25.15))') latvec(:,2) 
  write(87,'(3(1x,es25.15))') latvec(:,3) 
  write(87,'(a)') "%endblock LatticeVectors"
  write(87,*) " "
  write(87,*) " "
  write(87,*) "AtomicCoordinatesFormat Fractional" 
  write(87,*) "%block AtomicCoordinatesAndAtomicSpecies"
  do iat=1,parini%nat
  write(87,'(3(1x,es25.15),i5)') xred(:,iat),parini%typat_global(iat)
  enddo 
  write(87,*) "%endblock AtomicCoordinatesAndAtomicSpecies"
  close(87)
  end  subroutine
  
  subroutine get_output_siesta_geopt(parini,latvec,xred,fcart,energy,strten)
  !use global, only: nat
  !use defs_basis
  use mod_parini, only: typ_parini
  implicit none
  type(typ_parini), intent(in):: parini
  integer:: io,i,iat,itmp1,itmp2
  real(8):: fcart(3,parini%nat),energy,strten(6),value,str_matrix(3,3),latvec(3,3),xred(3,parini%nat)
  character(11):: ch_tmp
  open(unit=32,file="FORCE_STRESS")
  !First value is energy in Ry
  read(32,*) value   
  energy=value*0.5d0  !Convert to Hartree
  !Now get the stress tensor in Ry/Bohr**3
  read(32,*) str_matrix(:,1) 
  read(32,*) str_matrix(:,2) 
  read(32,*) str_matrix(:,3) 
       strten(1)=str_matrix(1,1)
       strten(2)=str_matrix(2,2)
       strten(3)=str_matrix(3,3)
       strten(6)=str_matrix(1,2)
       strten(5)=str_matrix(1,3)
       strten(4)=str_matrix(2,3)
  strten=strten*5.d0  !Convert to Hartree
  read(32,*) value !Temporary
  !Now get atomic forces
  do iat=1,parini%nat
  read(32,*) value,value,fcart(:,iat)
  enddo
  !Convert back to Ha/Bohr instead of Ry/Bohr
  fcart=fcart*0.5d0
  close(32)
  !write(*,*) "Energy", energy
  !write(*,*) "Stress", strten
  !write(*,*) "Forces"
  !write(*,*) fcart 
  open(unit=32,file="siestarun.1.STRUCT_OUT")
  read(32,*) latvec(:,1)
  read(32,*) latvec(:,2)
  read(32,*) latvec(:,3)
  latvec=latvec/Bohr_Ang
  read(32,*) iat
  if(iat.ne.parini%nat) stop "Something wrong with atomic number!"
  do iat=1,parini%nat
  read(32,*) itmp1,itmp2,xred(:,iat)
  enddo
  close(32)
  end subroutine
  end module interface_siesta
