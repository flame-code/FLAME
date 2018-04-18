module interface_cp2k
  use global
  use defs_basis
  !use cell_utils

  implicit none

  private
  public :: &
    make_input_cp2k, &
    get_output_cp2k


contains


  subroutine make_input_cp2k(parini,latvec,xred,iprec,ka,kb,kc,getwfk)
  !This routine will append some informations to a file already containing some informations about the cp2k runs
  !The informations appended are:
  !-Read/dont read wavefunction from file
  !-The kpoint mesh
  !-The atomic informations
  !use global, only: nat,ntypat,znucl,typat,dkpt1,dkpt2,ka1,kb1,kc1,max_kpt,reuse_kpt,char_type
  !use defs_basis, only: Bohr_Ang
  use mod_parini, only: typ_parini
  implicit none
  type(typ_parini), intent(in):: parini
  real(8):: xred(3,parini%nat)
  real(8):: dproj(6),acell(3),rprim(3,3),latvec(3,3),dkpt,k_latvec(3,3)
  real(8), allocatable:: k_xcart(:,:,:,:,:)
  integer:: iat,iprec,ka,kb,kc,k,l,m
  logical:: getwfk
  character(1):: fn
  character(150):: command
  
  
  if(iprec==1) then
  dkpt=parini%dkpt1
  else
  dkpt=parini%dkpt2
  endif
  
  
  command= "rm -f cp2k.in"
  call system(command)
  command= "rm -f cp2k_out.bak"
  call system(command)
  command= "cp  cp2k.out cp2k_out.bak"
  call system(command)
  write(fn,'(i1.1)') iprec
  command = "cp -f cp2k."//fn//".in cp2k.in"
  call system(command)
  
  open(unit=87,file="cp2k_subsys.in")
  !Kpoint mesh
  write(87,'(a,3(1x,i5),a)') "# Expansion with respect to k-mesh ",ka,kb,kc,"  # Number of gridpoints in each dimension"
  write(87,'(a)') " " 
  
  
  !The atomic positions and the unit cell is defined here
  allocate(k_xcart(3,parini%nat,ka,kb,kc))
  call k_expansion(parini,latvec,xred,ka,kb,kc,k_latvec,k_xcart)
  write(87,'(a)') "# Definition of the unit cell"
  write(87,'(a)') "  &CELL"
  write(87,'(a,3(1x,es25.15))') "      A ",k_latvec(:,1)*Bohr_Ang
  write(87,'(a,3(1x,es25.15))') "      B ",k_latvec(:,2)*Bohr_Ang
  write(87,'(a,3(1x,es25.15))') "      C ",k_latvec(:,3)*Bohr_Ang
  write(87,'(a)') "  &END CELL"
  write(87,'(a)') "  &COORD"
  do k=1,ka
  do l=1,kb
  do m=1,kc
  do iat=1,parini%nat
  write(87,'(a2,3(1x,es25.15))') trim(parini%char_type(parini%typat_global(iat))),k_xcart(:,iat,k,l,m)*Bohr_Ang
  enddo
  enddo
  enddo
  enddo
  write(87,'(a)') "  &END COORD"
  close(87)
  deallocate(k_xcart)
  end subroutine
  
  subroutine get_output_cp2k(parini,fcart,energy,strten)
  use defs_basis
  !Since its a single call, we only have forces and stresses from one configuration!
  use mod_parini, only: typ_parini
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
  open(unit=32,file="cp2k.out")
  do while(.true.)
   read(32,'(a250)',end=99)all_line
  !!write(*,*) all_line
   n = len_trim(all_line)
  
   k = index(all_line(1:n),"Total energy [eV]:")
    if(k.ne.0) then
      m = len_trim(all_line)
      l = scan(all_line(1:m),":",.true.)
      read(all_line(l+1:m),*) energy
      cycle
    endif
   k = index(all_line(1:n),"Atoms:")
    if(k.ne.0) then
      m = len_trim(all_line)
      l = scan(all_line(1:m),":",.true.)
      read(all_line(l+1:m),*) nat_cell
      cycle
    endif
   k = index(all_line(1:n),"Atom   Kind   Element")
    if(k.ne.0) then
      do iat=1,parini%nat
      read(32,*) int_tmp,int_tmp,ch_tmp,fcart(:,iat) 
      enddo
      cycle
    endif
   k = index(all_line(1:n),"STRESS TENSOR [GPa]")
    if(k.ne.0) then
    read(32,'(a250)',end=99)all_line
    read(32,'(a250)',end=99)all_line
      do iat=1,3
      read(32,*) ch_tmp,str_matrix(:,iat)
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
  if(energy==1.d10.or.strten(1)==1.d10.or.fcart(1,1)==1.d10.or.nat_cell==0) stop "Could not find all requested variables"
  !Transform all to bohr
  energy=energy/real(nat_cell,8)*real(parini%nat,8)/Ha_eV
  strten=strten/HaBohr3_GPa
  !strten=0.d0
  fcart=fcart!/Ha_eV*Bohr_Ang
  end subroutine
  end module interface_cp2k
  
