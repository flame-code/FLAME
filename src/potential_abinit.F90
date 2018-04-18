module interface_abinit
  use global
  use defs_basis
!  !use cell_utils

  implicit none

  private
  public :: &
    make_input_abinit, &
    get_output_abinit, &
    get_dos_abinit

contains

  !This routine will append some informations to a file already containing some informations about the abininit runs
  !The informations appended are:
  ! - Read/dont read wavefunction from file
  ! - The kpoint mesh
  ! - The atomic informations
  subroutine make_input_abinit(parini,latvec, xred, iprec, ka, kb, kc, getwfk, dos)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
    real(8), intent(in) :: latvec(3,3)
    real(8), intent(in) :: xred(3,parini%nat)
    integer, intent(inout) :: ka, kb, kc
    logical, intent(in) :: getwfk
    logical, intent(in), optional :: dos

    real(8):: dproj(6), acell(3), rprim(3,3), dkpt
    integer:: iat, iprec
    character(1):: fn
    character(150):: command

    if(iprec == 1) then
      dkpt = parini%dkpt1
    else
      dkpt = parini%dkpt2
    endif

    call system("rm -f abirun.in")
    call system("rm -f abirun_out.bak")
    call system("cp  abirun.out abirun_out.bak")
    write(fn,'(i1.1)') iprec
    call system("cp -f abirun."//fn//".in abirun.in")

    open(unit=87,file="abirun.in",ACCESS="APPEND")

    if(present(dos)) then
      write(87,'(a)') "prtdos 2         # Print the density of states"
    end if

    ! Kpoint mesh
    write(87,'(a)') "# Definition of the k-point mesh"
    write(87,'(a)') "kptopt   1   # Option for the automatic generation of k points, taking into account the symmetry"
    if(dkpt==0.d0.or.parini%abinit_kpt_mode==2) then
       write(87,'(a,3(1x,i5),a)') "ngkpt  ", ka, kb, kc, "  # Number of gridpoints in each dimension"
       write(87,'(a)') "nshiftk  1"
       write(87,'(a)') "shiftk  0.5 0.5 0.5  # These shifts will be the same for all grids"
    elseif(parini%abinit_kpt_mode==1) then
      write(87,'(a)') "# Automatic generation of the kpt mesh"
      write(87,'(a,es25.15,a)') "kptrlen ",dkpt," # K-mesh length"
    else
      stop "Wrong kpt option"
    endif

    if(getwfk) then
     write(87,'(a)') " irdwfk 1            # This is to speed up the calculation, by restarting"
     write(87,'(a)') "                     # from previous wavefunctions, transferred from the old "
     write(87,'(a)') "                     # to the new k-points."
     write(87,'(a)') "                     # Careful: be sure to move the previous WFK file to the"
     write(87,'(a)') "                     # input file!!!"
    else
     call system("./rm_WFK.sh")
    endif
    write(87,'(a)') " " 
    write(87,'(a)') " " 

    !The atomic positions and the unit cell is defined here
    !The data are read from a file, called posinp.ascii
    !The units are expected to be in angstroem
    !This particular part is made for Silicon and Hydrogen cells
    write(87,'(a)') "# Definition of the unit cell"
    call latvec2acell_rprim(latvec, acell, rprim)
    write(87,*) "natom"
    write(87,*)  parini%nat
    write(87,*) "ntypat"
    write(87,*)  parini%ntypat_global !Number of atom types        
    write(87,*) "znucl"
    write(87,*)  int(parini%znucl)  !Atom type
    write(87,*) "typat" 
    write(87,*)  parini%typat_global
    write(87,*) " "
    write(87,*) " "
    write(87,*) "acell"
    write(87,*) acell
    write(87,*) "rprim"
    write(87,'(3(1x,es25.15))')  rprim(:,1)
    write(87,'(3(1x,es25.15))')  rprim(:,2)
    write(87,'(3(1x,es25.15))')  rprim(:,3)
    write(87,*) "xred"
    do iat=1,parini%nat
      write(87,'(3(1x,es25.15))') xred(:,iat)
    enddo
    close(87)
  end subroutine make_input_abinit


  subroutine get_output_abinit(parini,fcart, energy, strten)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
    real(8), intent(out) ::  fcart(3, parini%nat), energy, strten(6)
    
    integer :: io, iat
    real(8) :: value
    character(11)  :: ch_tmp
    character(150) ::all_line
    
    ch_tmp="old"
    open(unit=32,file="abirun.out")
    do
      read(32,'(a150)',iostat=io) all_line
      
      ! write(*,*) all_line
      if(io.lt.0) stop "File ended unexpectedly"
      read(all_line, *, iostat=io) ch_tmp, value
      
      if(trim(ch_tmp).eq."etotal") then
        write(*,*) "etotal found", value
        energy=value
        
      elseif (trim(ch_tmp).eq."fcart") then
        write(*,*) "fcart found"
        read(all_line,*,iostat=io)ch_tmp, fcart(:,1)
        if(parini%nat.gt.1) then
          do iat = 2, parini%nat
            read(32,*) fcart(:,iat)
          enddo
        endif
        
      elseif(trim(ch_tmp).eq."strten") then
        write(*,*) "strten found"
        read(all_line,*,iostat=io) ch_tmp, strten(1:3)
        read(32,*) strten(4:6)
        exit ! nothing else to read
      endif
    enddo
    
    close(32)
  end subroutine get_output_abinit


  subroutine get_dos_abinit(fdos, efermi)
    real(8), intent(out) :: fdos, efermi

    character(150) :: all_line
    integer :: n, k, l, m
    real(8) :: density, energy, diff

    diff = 1.d10
    open(unit=32,file="abirun__xo_DOS")
    do
      read(32,'(a150)',end=99) all_line

      n = len_trim(all_line)
      k = index(all_line(1:n), "Fermi")
      if(k.ne.0) then
        m = len_trim(all_line)
        l = scan(all_line(1:m), ":", .true.)
        read(all_line(l+1:m),*) efermi
      endif

      n = len_trim(all_line)
      k = index(all_line(1:n), "energy(Ha)     DOS  integrated DOS")
      if(k.ne.0) then
        exit
      endif
    enddo

    do 
      read(32,'(a150)',end=99) all_line
      read(all_line,*) energy,density
      if(abs(efermi-energy).lt.diff) then
        fdos = density 
        diff = abs(efermi-energy)
      endif
    enddo

99  continue
    close(32)
  end subroutine get_dos_abinit


end module interface_abinit
