module interface_espresso
  use global
  use defs_basis
  !use cell_utils

  implicit none

  private
  public :: &
    make_input_espresso,        &
    get_output_espresso,        &
    espresso_geopt,             &
    make_input_espresso_geopt,  &
    get_output_espresso_geopt   
 

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
  subroutine make_input_espresso(parini,latvec, xred, iprec, ka, kb, kc, getwfk, dos)
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

    write(fn,'(i1.1)') iprec
    call system("cp -f espresso."//fn//".CONTROL   espresso.CONTROL")
    call system("cp -f espresso."//fn//".SYSTEM    espresso.SYSTEM")
    call system("cp -f espresso."//fn//".ELECTRONS espresso.ELECTRONS")

    open(unit=87, file="espresso.CONTROL", ACCESS="APPEND")
        write(87,'(a)') 'calculation =   "scf"    ,'
        write(87,'(a)') 'tstress     =   .true.   ,'
        write(87,'(a)') 'tprnfor     =   .true.   ,'
!Close block    
        write(87,'(a)') "/"
    close(87)
    open(unit=87, file="espresso.SYSTEM", ACCESS="APPEND")
        write(87,'(a,i5,a)') 'nat =  ',parini%nat,' ,'
        write(87,'(a,i5,a)') 'ntyp =  ',parini%ntypat_global,' ,'
        write(87,'(a)')      'nosym = .true. ,'

!Close block    
        write(87,'(a)') "/"
    close(87)
    open(unit=87, file="espresso.ELECTRONS", ACCESS="APPEND")
!Close block    
        write(87,'(a)') "/"
    close(87)
!Create Block with atoms, structure, etc.
    open(unit=87,file="espresso.STRUCT")
        write(87,'(a)') "ATOMIC_SPECIES"
        do itype=1,parini%ntypat_global
           write(87,'(a,2x,f10.5,2x,a)') trim(parini%char_type(itype)),parini%amu(itype),trim(parini%char_type(itype))//".PSP"
        enddo
        write(87,'(a)') "ATOMIC_POSITIONS crystal"
        do iat=1,parini%nat
           write(87,'(a,2x,3(es25.15))') trim(parini%char_type(parini%typat_global(iat))),xred(:,iat)
        enddo
        write(87,'(a)') "CELL_PARAMETERS bohr"
           write(87,'(3(es25.15))') latvec(:,1)
           write(87,'(3(es25.15))') latvec(:,2)
           write(87,'(3(es25.15))') latvec(:,3)
    close(87)
    open(unit=87,file="espresso.KPOINTS")
        write(87,'(a)') "K_POINTS automatic"
           write(87,'(i5,i5,i5,5x,i5,i5,i5)') ka,kb,kc,0,0,0
    close(87)
!    if(.not.getwfk) then
!      call system("./rm_wavecar.sh")
!      write(87,'(a)') "LWAVE  = .FALSE."
!      write(87,'(a)') "LCHARG = .FALSE."
!    endif
    call system("cat espresso.CONTROL espresso.SYSTEM espresso.ELECTRONS espresso.STRUCT espresso.KPOINTS>espresso.in")
  end subroutine make_input_espresso




  subroutine get_output_espresso(parini,fcart,energy,strten)
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
  
  energy=1.d10
  fcart=1.d10
  strten=1.d10
  
  ch_tmp="old"
  open(unit=32,file="espresso.out")
  do while(.true.)
   read(32,'(a150)',end=99)all_line
  !!write(*,*) all_line
   n = len_trim(all_line)
  
   k = index(all_line(1:n),"Forces acting on atoms")
    if(k.ne.0) then
      write(*,*) "FORCES FOUND"
      read(32,*)
      do iat=1,parini%nat
      read(32,*)ch_tmp,i_tmp,ch_tmp,i_tmp,ch_tmp,ch_tmp,fcart(:,iat)
      write(*,*) fcart(:,iat)
      enddo
      cycle
    endif
   k = index(all_line(1:n),"total   stress ")
    if(k.ne.0) then
      write(*,*) "STRESSES FOUND"
      read(32,*)str_matrix(:,1)
      read(32,*)str_matrix(:,2)
      read(32,*)str_matrix(:,3)
      write(*,*)str_matrix(:,1)
      write(*,*)str_matrix(:,2)
      write(*,*)str_matrix(:,3)
       strten(1)=-str_matrix(1,1)
       strten(2)=-str_matrix(2,2)
       strten(3)=-str_matrix(3,3)
       strten(6)=-str_matrix(1,2)
       strten(5)=-str_matrix(1,3)
       strten(4)=-str_matrix(2,3)
      cycle
    endif
   
   k = index(all_line(1:n),"!    total energy              =")
    if(k.ne.0) then
            write(*,*) "ENERGY FOUND"
            read(all_line(1:n),*) ch_tmp,ch_tmp,ch_tmp,ch_tmp,energy
            write(*,*) energy
    cycle
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
  
  !Transform all to bohr
  energy=energy*0.5d0
  strten=strten*0.5d0
  fcart=fcart*0.5d0
  end subroutine
  
  subroutine espresso_geopt(parini,latvec,xred,fcart,strten,energy,iprec,ka,kb,kc,counter)
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
   call make_input_espresso_geopt(parini,latvec,xred,iprec,ka,kb,kc,getwfk)
  ! call system("sleep 1")
  !Run the job NOW!
   call system("./runjob_geoespresso.sh")
  ! call system("sleep 1")
  !Now harvest the structure, energy, forces, etc
   call get_output_espresso_geopt(parini,latvec,xred,fcart,energy,strten)
  !Check how many iterations have been needed
!   call system("grep Conjugate OUTCAR_geo_a |wc -l>tmp_count") 
!   call system("grep Conjugate OUTCAR_geo_b |wc -l>>tmp_count") 
!   call system("grep Conjugate OUTCAR_geo_c |wc -l>>tmp_count") 
   call system("rm -f tmp_count")
   call system("grep '!    total energy              =     ' espresso_geo_a.out   |wc -l>tmp_count") 
   call system("grep '!    total energy              =     ' espresso_geo_b.out   |wc -l>>tmp_count") 
   call system("grep '!    total energy              =     ' espresso_geo_c.out   |wc -l>>tmp_count") 
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
  
  subroutine make_input_espresso_geopt(parini,latvec,xred,iprec,ka,kb,kc,getwfk)
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
  
  !Rescale tolmxf to units of eV/Ang from Ha/Bohr
  HaBohr_eVAng=Ha_eV/Bohr_Ang
  
    if(iprec == 1) then
      dkpt = parini%dkpt1
    else
      dkpt = parini%dkpt2
    endif

    write(fn,'(i1.1)') iprec
    call system("cp -f espresso."//fn//".a.CONTROL   espresso_geo_a.CONTROL")
    call system("cp -f espresso."//fn//".a.SYSTEM    espresso_geo_a.SYSTEM")
    call system("cp -f espresso."//fn//".a.ELECTRONS espresso_geo_a.ELECTRONS")
    call system("cp -f espresso."//fn//".b.CONTROL   espresso_geo_b.CONTROL")
    call system("cp -f espresso."//fn//".b.SYSTEM    espresso_geo_b.SYSTEM")
    call system("cp -f espresso."//fn//".b.ELECTRONS espresso_geo_b.ELECTRONS")
    call system("cp -f espresso."//fn//".CONTROL    espresso_geo_c.CONTROL")
    call system("cp -f espresso."//fn//".SYSTEM     espresso_geo_c.SYSTEM")
    call system("cp -f espresso."//fn//".ELECTRONS  espresso_geo_c.ELECTRONS")

!BLOCK A_GEOPT----------------------------------------------
    open(unit=87, file="espresso_geo_a.CONTROL", ACCESS="APPEND")
        if(((all(parini%fixlat(1:6))).and.(.not.parini%fixlat(7))).or.parini%bc==2) then
            write(87,'(a)')          'calculation =   "relax"    ,'
        else
            write(87,'(a)')          'calculation =   "vc-relax"  ,'
        endif
        write(87,'(a)')           'tstress     =   .true.   ,'
        write(87,'(a)')           'tprnfor     =   .true.   ,'
        write(87,'(a,i5,a)')      'nstep       = ',int(parini%paropt_geopt%nit*0.5d0),' ,'
        write(87,'(a,es15.7,a)')  'forc_conv_thr = ',parini%paropt_geopt%fmaxtol*8.d0*2.d0,' ,'
!Close block    
        write(87,'(a)') "/"
    close(87)
    open(unit=87, file="espresso_geo_a.SYSTEM", ACCESS="APPEND")
        write(87,'(a,i5,a)') 'nat =  ',parini%nat,' ,'
        write(87,'(a,i5,a)') 'ntyp = ',parini%ntypat_global,' ,'
        write(87,'(a)')      'nosym = .true. ,'
!Close block    
        write(87,'(a)') "/"
    close(87)
    open(unit=87, file="espresso_geo_a.ELECTRONS", ACCESS="APPEND")
!Close block    
        write(87,'(a)') "/"
    close(87)
    open(unit=87, file="espresso_geo_a.IONS")
        write(87,'(a)') "&IONS"
        write(87,'(a)') "ion_dynamics = 'bfgs' ,"
!Close block    
        write(87,'(a)') "/"
    close(87)
    open(unit=87, file="espresso_geo_a.CELL")
        write(87,'(a)') "&CELL"
        if(((all(parini%fixlat(1:6))).and.(.not.parini%fixlat(7))).or.parini%bc==2) then
             continue
        else
        write(87,'(a)')          "cell_dynamics   = 'bfgs' ,"
        write(87,'(a,es15.7,a)') "press           = ",parini%target_pressure_gpa*10.d0
        write(87,'(a,es15.7,a)') "press_conv_thr  = ",parini%paropt_geopt%fmaxtol/parini%paropt_geopt%strfact*10.d0*HaBohr3_GPa*8.d0
        write(87,'(a)')          "cell_factor     = 4.d0"
        endif
!Close block    
        write(87,'(a)') "/"
    close(87)


!BLOCK B_GEOPT----------------------------------------------
    open(unit=87, file="espresso_geo_b.CONTROL", ACCESS="APPEND")
        if(((all(parini%fixlat(1:6))).and.(.not.parini%fixlat(7))).or.parini%bc==2) then
            write(87,'(a)')          'calculation =   "relax"    ,'
        else
            write(87,'(a)')          'calculation =   "vc-relax"  ,'
        endif
        write(87,'(a)')           'tstress     =   .true.   ,'
        write(87,'(a)')           'tprnfor     =   .true.   ,'
        write(87,'(a,i5,a)')      'nstep       = ',int(parini%paropt_geopt%nit*0.5d0),' ,'
        write(87,'(a,es15.7,a)')  'forc_conv_thr = ',parini%paropt_geopt%fmaxtol*2.d0,' ,'
!Close block    
        write(87,'(a)') "/"
    close(87)
    open(unit=87, file="espresso_geo_b.SYSTEM", ACCESS="APPEND")
        write(87,'(a,i5,a)') 'nat =  ',parini%nat,' ,'
        write(87,'(a,i5,a)') 'ntyp = ',parini%ntypat_global,' ,'
        write(87,'(a)')      'nosym = .true. ,'
!Close block    
        write(87,'(a)') "/"
    close(87)
    open(unit=87, file="espresso_geo_b.ELECTRONS", ACCESS="APPEND")
!Close block    
        write(87,'(a)') "/"
    close(87)
    open(unit=87, file="espresso_geo_b.IONS")
        write(87,'(a)') "&IONS"
        write(87,'(a)') "ion_dynamics = 'bfgs' ,"
!Close block    
        write(87,'(a)') "/"
    close(87)
    open(unit=87, file="espresso_geo_b.CELL")
        write(87,'(a)') "&CELL"
        if(((all(parini%fixlat(1:6))).and.(.not.parini%fixlat(7))).or.parini%bc==2) then
           continue
        else
        write(87,'(a)')          "cell_dynamics   = 'bfgs' ,"
        write(87,'(a,es15.7,a)') "press           = ",parini%target_pressure_gpa*10.d0
        write(87,'(a,es15.7,a)') "press_conv_thr  = ",parini%paropt_geopt%fmaxtol/parini%paropt_geopt%strfact*10.d0*HaBohr3_GPa
        write(87,'(a)')          "cell_factor     = 4.d0"
        endif
!Close block    
        write(87,'(a)') "/"
    close(87)


!BLOCK C_GEOPT----------------------------------------------
    open(unit=87, file="espresso_geo_c.CONTROL", ACCESS="APPEND")
        write(87,'(a)')           'calculation =   "scf"    ,'
        write(87,'(a)')           'tstress     =   .true.   ,'
        write(87,'(a)')           'tprnfor     =   .true.   ,'
!Close block    
        write(87,'(a)') "/"
    close(87)
    open(unit=87, file="espresso_geo_c.SYSTEM", ACCESS="APPEND")
        write(87,'(a,i5,a)') 'nat =  ',parini%nat,' ,'
        write(87,'(a,i5,a)') 'ntyp = ',parini%ntypat_global,' ,'
        write(87,'(a)')      'nosym = .true. ,'
!Close block    
        write(87,'(a)') "/"
    close(87)
    open(unit=87, file="espresso_geo_c.ELECTRONS", ACCESS="APPEND")
!Close block    
        write(87,'(a)') "/"
    close(87)
    open(unit=87, file="espresso_geo_c.IONS")
        write(87,'(a)') "&IONS"
!Close block    
        write(87,'(a)') "/"
    close(87)
    open(unit=87, file="espresso_geo_c.CELL")
        write(87,'(a)') "&CELL"
!Close block    
        write(87,'(a)') "/"
    close(87)

!Create Block with atoms, structure, etc.
    open(unit=87,file="espresso.STRUCT")
        write(87,'(a)') "ATOMIC_SPECIES"
        do itype=1,parini%ntypat_global
           write(87,'(a,2x,f10.5,2x,a)') trim(parini%char_type(itype)),parini%amu(itype),trim(parini%char_type(itype))//".PSP"
        enddo
        write(87,'(a)') "ATOMIC_POSITIONS crystal"
        do iat=1,parini%nat
        if(parini%fixat(iat)) then
              write(87,'(a,2x,3(es25.15),a)') trim(parini%char_type(parini%typat_global(iat))),xred(:,iat),' 0 0 0 '
        else
              write(87,'(a,2x,3(es25.15),a)') trim(parini%char_type(parini%typat_global(iat))),xred(:,iat),' 1 1 1 '
        endif
        enddo
        write(87,'(a)') "CELL_PARAMETERS bohr"
           write(87,'(3(es25.15))') latvec(:,1)
           write(87,'(3(es25.15))') latvec(:,2)
           write(87,'(3(es25.15))') latvec(:,3)
    close(87)
!Create Block with atoms, structure, etc.
    open(unit=87,file="espresso.KPOINTS")
        write(87,'(a)') "K_POINTS automatic"
           write(87,'(i5,i5,i5,5x,i5,i5,i5)') ka,kb,kc,0,0,0
    close(87)
!    if(.not.getwfk) then
!      call system("./rm_wavecar.sh")
!      write(87,'(a)') "LWAVE  = .FALSE."
!      write(87,'(a)') "LCHARG = .FALSE."
!    endif
!    call system("cat espresso.CONTROL espresso.SYSTEM espresso.ELECTRONS espresso.STRUCT>espresso.in")

  end  subroutine
  
  subroutine get_output_espresso_geopt(parini,latvec,xred,fcart,energy,strten)
  !use global, only: nat,target_pressure_gpa
  !use defs_basis
  !Since its a single call, we only have forces and stresses from one configuration!
  use mod_parini, only: typ_parini
  implicit none
  type(typ_parini), intent(in):: parini
  integer:: io,i,iat,n,k,l,m,i_tmp
  real(8):: fcart(3,parini%nat),energy,strten(6),value,latvec(3,3),xred(3,parini%nat),str_matrix(3,3),vol,a(3,3),scaling,xcart(3,parini%nat),alat
  character(11):: ch_tmp
  character(150)::all_line
  logical:: vasp_5
  !if vasp is version 5.x, use vasp_5=.true.
  
  vasp_5=.true.
  
  energy=1.d10
  fcart=1.d10
  strten=1.d10
  
  ch_tmp="old"
  open(unit=32,file="espresso_geo_c.out")
  do while(.true.)
   read(32,'(a150)',end=99)all_line
  !!write(*,*) all_line
   n = len_trim(all_line)
  
   k = index(all_line(1:n),"Forces acting on atoms")
    if(k.ne.0) then
      write(*,*) "FORCES FOUND"
      read(32,*)
      do iat=1,parini%nat
      read(32,*)ch_tmp,i_tmp,ch_tmp,i_tmp,ch_tmp,ch_tmp,fcart(:,iat)
      write(*,*) fcart(:,iat)
      enddo
      cycle
    endif
   k = index(all_line(1:n),"total   stress ")
    if(k.ne.0) then
      write(*,*) "STRESSES FOUND"
      read(32,*)str_matrix(:,1)
      read(32,*)str_matrix(:,2)
      read(32,*)str_matrix(:,3)
      write(*,*)str_matrix(:,1)
      write(*,*)str_matrix(:,2)
      write(*,*)str_matrix(:,3)
       strten(1)=-str_matrix(1,1)
       strten(2)=-str_matrix(2,2)
       strten(3)=-str_matrix(3,3)
       strten(6)=-str_matrix(1,2)
       strten(5)=-str_matrix(1,3)
       strten(4)=-str_matrix(2,3)
      cycle
    endif
   
   k = index(all_line(1:n),"!    total energy              =")
    if(k.ne.0) then
            write(*,*) "ENERGY FOUND"
            read(all_line(1:n),*) ch_tmp,ch_tmp,ch_tmp,ch_tmp,energy
            write(*,*) energy
    cycle
    endif

   k = index(all_line(1:n),"lattice parameter (alat)")
    if(k.ne.0) then
    write(*,*) "alat found"
    write(*,*) all_line(1:n)
      l = scan(all_line(1:n),"=",.true.)
      read(all_line(l+1:n),*) alat
    write(*,*) alat
    cycle
    endif
   
   k = index(all_line(1:n),"crystal axes: (cart. coord. in units of alat)")
    if(k.ne.0) then
      read(32,'(a150)',end=99)all_line
      m = len_trim(all_line)
      l = scan(all_line(1:m),"(",.true.)
      read(all_line(l+1:m),*) latvec(:,1)
      read(32,'(a150)',end=99)all_line
      read(all_line(l+1:m),*) latvec(:,2)
      read(32,'(a150)',end=99)all_line
      read(all_line(l+1:m),*) latvec(:,3)
    write(*,*) latvec
    latvec=latvec*alat
    cycle
    endif

   k = index(all_line(1:n),"site n.     atom                  positions (alat units)")
    if(k.ne.0) then
      write(*,*) "xcart found"
      do iat=1,parini%nat
      read(32,'(a150)',end=99)all_line
      m = len_trim(all_line)
      l = scan(all_line(1:m),"(",.true.)
      read(all_line(l+1:m),*) xcart(:,iat)
      enddo
      write(*,*) xcart
      xcart=xcart*alat
      call rxyz_cart2int(latvec,xred,xcart,parini%nat)
    cycle
    endif

    k = index(all_line(1:n),"CELL_PARAMETERS (bohr)")
    if(k.ne.0) then
      read(32,*,end=99) latvec(:,1)
      read(32,*,end=99) latvec(:,2)
      read(32,*,end=99) latvec(:,3)
    write(*,*) latvec
    cycle
    endif

   k = index(all_line(1:n),"ATOMIC_POSITIONS (crystal)")
    if(k.ne.0) then
      write(*,*) "xcart found"
      do iat=1,parini%nat
      read(32,*,end=99) ch_tmp, xred(:,iat)
      enddo
      write(*,*) xred
    cycle
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
  
  
  !Transform all to bohr
  energy=energy*0.5d0
  strten=strten*0.5d0
  fcart=fcart*0.5d0
  end subroutine
  end module interface_espresso

!!    subroutine rewrite_struct_espresso(filename_in,filename_out_struct,filename_out_kpt)
!!    !use global, only: nat,target_pressure_gpa
!!    !use defs_basis
!!    !Since its a single call, we only have forces and stresses from one configuration!
!!    implicit none
!!    integer:: io,i,iat,n,k,l,m,i_tmp,nat,ntypat,itype
!!    real(8):: energy,strten(6),value,latvec(3,3),str_matrix(3,3),vol,a(3,3),scaling,alat,r_tmp
!!    real(8),allocatable:: fcart(:,:),xred(:,:),xcart(:,:),amu(:)
!!    integer,allocatable:: fix(:,:)
!!    character(2),allocatable:: char_type(:),typat_char(:)
!!    character(11):: ch_tmp
!!    character(150)::all_line
!!    character(40) ::filename_in,filename_out_struct,filename_out_kpt
!!    real(8):: dkpt1,dkpt2,dkpt_12(2)
!!    integer:: ka,kb,kc,kpt_abc(3)
!!    logical:: tmpl,auto_kpt,found
!!    !if vasp is version 5.x, use vasp_5=.true.
!!    
!!    energy=1.d10
!!    fcart=1.d10
!!    strten=1.d10
!!    
!!    ch_tmp="old"
!!    open(unit=32,file=trim(filename_in))
!!    !First loop to get the number of atoms, number of types, etc  
!!    do while(.true.)
!!     read(32,'(a150)',end=99)all_line
!!    !!write(*,*) all_line
!!     n = len_trim(all_line)
!!     k = index(all_line(1:n),"number of atoms/cell      =")
!!      if(k.ne.0) then
!!        write(*,*) "NAT"
!!        l = scan(all_line(1:n),"=",.true.)
!!        read(all_line(l+1:n),*) nat
!!        cycle
!!      endif
!!      k = index(all_line(1:n),"number of atomic types    =")
!!      if(k.ne.0) then
!!        write(*,*) "NTYPAT"
!!        l = scan(all_line(1:n),"=",.true.)
!!        read(all_line(l+1:n),*) ntypat
!!        cycle
!!      endif
!!     enddo
!!     close(99)
!!     allocate(fix(3,nat),fcart(3,nat),xred(3,nat),xcart(3,nat),char_type(ntypat),amu(ntypat),typat_char(nat))
!!  
!!  
!!  !Do second loop to get the other stuff
!!    open(unit=32,file=trim(filename_in))
!!    !First loop to get the number of atoms, number of types, etc  
!!    do while(.true.)
!!     read(32,'(a150)',end=99)all_line
!!    !!write(*,*) all_line
!!     n = len_trim(all_line)
!!     k = index(all_line(1:n),"total   stress ")
!!      if(k.ne.0) then
!!        write(*,*) "STRESSES FOUND"
!!        read(32,*)str_matrix(:,1)
!!        read(32,*)str_matrix(:,2)
!!        read(32,*)str_matrix(:,3)
!!        write(*,*)str_matrix(:,1)
!!        write(*,*)str_matrix(:,2)
!!        write(*,*)str_matrix(:,3)
!!         strten(1)=-str_matrix(1,1)
!!         strten(2)=-str_matrix(2,2)
!!         strten(3)=-str_matrix(3,3)
!!         strten(6)=-str_matrix(1,2)
!!         strten(5)=-str_matrix(1,3)
!!         strten(4)=-str_matrix(2,3)
!!        cycle
!!      endif
!!     
!!     k = index(all_line(1:n),"!    total energy              =")
!!      if(k.ne.0) then
!!              write(*,*) "ENERGY FOUND"
!!              read(all_line(1:n),*) ch_tmp,ch_tmp,ch_tmp,ch_tmp,energy
!!              write(*,*) energy
!!      cycle
!!      endif
!!  
!!     k = index(all_line(1:n),"lattice parameter (alat)")
!!      if(k.ne.0) then
!!      write(*,*) "alat found"
!!      write(*,*) all_line(1:n)
!!        l = scan(all_line(1:n),"=",.true.)
!!        read(all_line(l+1:n),*) alat
!!      write(*,*) alat
!!      cycle
!!      endif
!!     
!!     k = index(all_line(1:n),"crystal axes: (cart. coord. in units of alat)")
!!      if(k.ne.0) then
!!        read(32,'(a150)',end=99)all_line
!!        m = len_trim(all_line)
!!        l = scan(all_line(1:m),"(",.true.)
!!        read(all_line(l+1:m),*) latvec(:,1)
!!        read(32,'(a150)',end=99)all_line
!!        read(all_line(l+1:m),*) latvec(:,2)
!!        read(32,'(a150)',end=99)all_line
!!        read(all_line(l+1:m),*) latvec(:,3)
!!      write(*,*) latvec
!!      latvec=latvec*alat
!!      cycle
!!      endif
!!  
!!     k = index(all_line(1:n),"atomic species   valence    mass     pseudopotential")
!!      if(k.ne.0) then
!!              do itype=1,ntypat
!!              read(32,*,end=99) char_type(itype),r_tmp,amu(itype)
!!              enddo 
!!      cycle
!!      endif
!!  
!!  
!!  
!!     k = index(all_line(1:n),"site n.     atom                  positions (alat units)")
!!      if(k.ne.0) then
!!        write(*,*) "xcart found"
!!        do iat=1,nat
!!        read(32,'(a150)',end=99)all_line
!!        m = len_trim(all_line)
!!        l = scan(all_line(1:m),"(",.true.)
!!        read(all_line(1:m),*)i_tmp,typat_char(iat)
!!        read(all_line(l+1:m),*) xcart(:,iat)
!!        enddo
!!        write(*,*) xcart
!!        xcart=xcart*alat
!!        call rxyz_cart2int(latvec,xred,xcart,nat)
!!      cycle
!!      endif
!!  
!!      k = index(all_line(1:n),"CELL_PARAMETERS (bohr)")
!!      if(k.ne.0) then
!!        read(32,*,end=99) latvec(:,1)
!!        read(32,*,end=99) latvec(:,2)
!!        read(32,*,end=99) latvec(:,3)
!!      write(*,*) latvec
!!      cycle
!!      endif
!!  
!!     k = index(all_line(1:n),"ATOMIC_POSITIONS (crystal)")
!!      if(k.ne.0) then
!!        fix=1
!!        write(*,*) "xcart found"
!!        do iat=1,nat
!!        read(32,'(a150)',end=99)all_line
!!        read(all_line,*,iostat=io) ch_tmp, xred(:,iat),fix(:,iat)
!!        if(io.lt.0)  fix(:,iat)=1
!!        enddo
!!        write(*,*) xred
!!        write(*,*) fix
!!      cycle
!!      endif
!!    enddo
!!    
!!    99 continue 
!!    close(32)
!!    if(energy==1.d10.or.strten(1)==1.d10.or.fcart(1,1)==1.d10) stop "Could not find all requested variables"
!!    
!!    
!!    
!!    !Transform all to bohr
!!    energy=energy*0.5d0
!!    strten=strten*0.5d0
!!    fcart=fcart*0.5d0
!!  
!!  
!!  !Get infor from params_new.in
!!  open(unit=12,file="params_new.in")
!!    do while(.true.)
!!     read(12,'(a150)',end=97)all_line
!!     n = len_trim(all_line)
!!  !Block KPT****************
!!  !AUTO_KPT
!!     call parse_logical("AUTO_KPT",8,all_line(1:n),n,auto_kpt,found)
!!     if(found) cycle
!!  !KPTMESH
!!     call parsearray_int("KPTMESH",7,all_line(1:n),n,kpt_abc(1:3),3,found)
!!     if(found) then
!!       ka=kpt_abc(1)
!!       kb=kpt_abc(2)
!!       kc=kpt_abc(3)
!!     endif
!!     if(found) cycle
!!  !KPTDEN
!!     call parsearray_real("KPTDEN",6,all_line(1:n),n,dkpt_12(1:2),2,found)
!!     if(found) then
!!       dkpt1=dkpt_12(1)
!!       dkpt2=dkpt_12(2)
!!     endif
!!     if(found) cycle
!!  !Block KPT****************
!!  enddo
!!   97 continue
!!    if(AUTO_KPT) then
!!      ka=0;kb=0;kc=0
!!    else
!!      dkpt1=0.d0
!!      dkpt2=0.d0
!!    endif
!!  
!!   close(12)
!!  if(dkpt1.ne.0.d0) then
!!  call find_kpt(ka,kb,kc,latvec,dkpt1)
!!  endif
!!  
!!  !Create Block with atoms, structure, etc.
!!      open(unit=87,file=trim(filename_out_struct))
!!          write(87,'(a)') "ATOMIC_SPECIES"
!!          do itype=1,ntypat
!!             write(87,'(a,2x,f10.5,2x,a)') trim(char_type(itype)),amu(itype),trim(char_type(itype))//".PSP"
!!          enddo
!!          write(87,'(a)') "ATOMIC_POSITIONS crystal"
!!          do iat=1,nat
!!                write(87,'(a,2x,3(es25.15),2x,3(i5))') trim(typat_char(iat)),xred(:,iat),fix(:,iat)
!!          enddo
!!          write(87,'(a)') "CELL_PARAMETERS bohr"
!!             write(87,'(3(es25.15))') latvec(:,1)
!!             write(87,'(3(es25.15))') latvec(:,2)
!!             write(87,'(3(es25.15))') latvec(:,3)
!!      close(87)
!!      open(unit=87,file=trim(filename_out_kpt))
!!          write(87,'(a)') "K_POINTS automatic"
!!             write(87,'(i5,i5,i5,5x,i5,i5,i5)') ka,kb,kc,0,0,0
!!      close(87)
!!    end subroutine
