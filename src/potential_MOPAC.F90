module interface_mopac
  use global
  use defs_basis
  !use cell_utils

  implicit none

  private
  public :: &
    make_input_mopac, &
    get_output_mopac


contains


  subroutine make_input_mopac(parini,latvec,xred,iprec,ka,kb,kc,getwfk)
  !This routine will append some informations to a file already containing some informations about the mopac runs
  !The informations appended are:
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
  
  
  command= "rm -f mopac.in"
  call system(command)
  command= "rm -f mopac_out.bak"
  call system(command)
  command= "cp  mopac.out mopac_out.bak"
  call system(command)
  command= "mv  mopac.aux mopac_aux.bak"
  call system(command)
  write(fn,'(i1.1)') iprec
  command = "cp -f mopac."//fn//".in mopac.mop"
  call system(command)
  
  open(unit=87,file="mopac.mop",ACCESS="APPEND")

  !Kpoint mesh
  write(87,'(a,3(1x,i5),a)') "# Expansion with respect to k-mesh ",ka,kb,kc,"  # Number of gridpoints in each dimension"
  write(87,'(a)') " " 
  
  
  !The atomic positions and the unit cell is defined here
  allocate(k_xcart(3,parini%nat,ka,kb,kc))
  call k_expansion(parini,latvec,xred,ka,kb,kc,k_latvec,k_xcart)
  do k=1,ka
  do l=1,kb
  do m=1,kc
  do iat=1,parini%nat
  if(.not.parini%fixat(iat)) then
      write(87,'(a2,3(1x,es25.15,1x,i5))') trim(parini%char_type(parini%typat_global(iat))),k_xcart(1,iat,k,l,m)*&
      &Bohr_Ang,1,k_xcart(2,iat,k,l,m)*Bohr_Ang,1,k_xcart(3,iat,k,l,m)*Bohr_Ang,1
  else
      write(87,'(a2,3(1x,es25.15,1x,i5))') trim(parini%char_type(parini%typat_global(iat))),k_xcart(1,iat,k,l,m)*&
      &Bohr_Ang,0,k_xcart(2,iat,k,l,m)*Bohr_Ang,0,k_xcart(3,iat,k,l,m)*Bohr_Ang,0
  endif
  enddo
  enddo
  enddo
  enddo
  write(87,'(a,3(1x,es25.15,1x,i5))') "Tv",k_latvec(1,1)*Bohr_Ang,1,k_latvec(2,1)*Bohr_Ang,1,k_latvec(3,1)*Bohr_Ang,1
  write(87,'(a,3(1x,es25.15,1x,i5))') "Tv",k_latvec(1,2)*Bohr_Ang,1,k_latvec(2,2)*Bohr_Ang,1,k_latvec(3,2)*Bohr_Ang,1
  write(87,'(a,3(1x,es25.15,1x,i5))') "Tv",k_latvec(1,3)*Bohr_Ang,1,k_latvec(2,3)*Bohr_Ang,1,k_latvec(3,3)*Bohr_Ang,1
  close(87)
  deallocate(k_xcart)
  end subroutine
  
  subroutine get_output_mopac(parini,fcart,energy,strten)
  use defs_basis
  !Since its a single call, we only have forces and stresses from one configuration!
  use mod_parini, only: typ_parini
  implicit none
  type(typ_parini), intent(in):: parini
  integer:: io,i,iat,n,k,l,m,int_tmp,nat_cell
  real(8):: fcart(3,parini%nat),energy,strten(6),value,latvec(3,3),xred(3,parini%nat),str_matrix(3,3),vol,a(3,3),scaling,transmat(3,3)
  real(8),allocatable:: fcart_tmp(:,:),tmp_grad(:)
  character(11):: ch_tmp
  character(250)::all_line,formatting
  
  energy=1.d10
  fcart=1.d10
  strten=1.d10
  nat_cell=0
  
  ch_tmp="old"
  open(unit=32,file="mopac.aux")
  do while(.true.)
   read(32,'(a250)',end=99)all_line
  !!write(*,*) all_line
   n = len_trim(all_line)
  
   k = index(all_line(1:n),"HEAT_OF_FORMATION:KCAL/MOL=")
    if(k.ne.0) then
      m = len_trim(all_line)
      l = scan(all_line(1:m),"=",.true.)
      read(all_line(l+1:m),*) energy
      write(*,*) "ENERGY", energy
      cycle
    endif
   k = index(all_line(1:n),"EMPIRICAL_FORMULA")
    if(k.ne.0) then
      m = len_trim(all_line)
      l = scan(all_line(1:m),"=",.true.)
      read(all_line(l+1:m),*) nat_cell
      write(*,*) "NATCELL", nat_cell
      cycle
    endif
   k = index(all_line(1:n),"TRANS_VECTS")
    if(k.ne.0) then
      read(32,*)latvec(:,1)
      read(32,*)latvec(:,2)
      read(32,*)latvec(:,3)
     call getvol(latvec,vol)
      write(*,*) "VOL", VOL
      cycle
    endif
   k = index(all_line(1:n),"GRADIENTS:KCAL")
    if(k.ne.0) then
       if(.not.allocated(fcart_tmp)) allocate(fcart_tmp(3,nat_cell))
       allocate(tmp_grad(10*ceiling(real((3*nat_cell+9),8)/10.d0)))
       do i=1,ceiling(real((3*nat_cell+9),8)/10.d0)
         read(32,'(10(f18.13))') tmp_grad((i-1)*10+1:i*10)
       enddo
       !transform fcart_tmp into Hartree/Bohr
       tmp_grad=tmp_grad/Ha_kcalmol*Bohr_Ang
       vol=vol/Bohr_ang**3
       fcart_tmp(:,:)=-reshape(tmp_grad(1:3*nat_cell),(/3,nat_cell/))
       str_matrix(:,:)=-reshape(tmp_grad(3*nat_cell+1:3*nat_cell+9),(/3,3/))
       write(*,*) fcart_tmp(:,:),str_matrix(:,1),str_matrix(:,2),str_matrix(:,3)

       fcart=fcart_tmp(:,1:parini%nat)
       deallocate(fcart_tmp)
       str_matrix=-matmul(str_matrix,transpose(latvec)/Bohr_ang)/vol
       write(*,*) "Trans Vects" 
       write(*,*) latvec(:,1)
       write(*,*) latvec(:,2)
       write(*,*) latvec(:,3)
       write(*,*) "Stress Matrix" 
       write(*,*) str_matrix(:,1)
       write(*,*) str_matrix(:,2)
       write(*,*) str_matrix(:,3)
       !Symmetrize the matrix, since noisy
       str_matrix(1,2)=0.5d0*(str_matrix(1,2)+str_matrix(2,1))
       str_matrix(1,3)=0.5d0*(str_matrix(1,3)+str_matrix(3,1))
       str_matrix(2,3)=0.5d0*(str_matrix(2,3)+str_matrix(3,2))
       strten(1)=str_matrix(1,1)
       strten(2)=str_matrix(2,2)
       strten(3)=str_matrix(3,3)
       strten(6)=str_matrix(1,2)
       strten(5)=str_matrix(1,3)
       strten(4)=str_matrix(2,3)
       write(*,*) "Pressure", HaBohr3_GPa*(strten(1)+strten(2)+strten(3))/3.d0
      cycle
    endif
  enddo
  
  99 continue 
  close(32)
  if(energy==1.d10.or.strten(1)==1.d10.or.fcart(1,1)==1.d10.or.nat_cell==0) stop&
   & "Could not find all requested variables"
  !Transform all to bohr
  energy=energy/real(nat_cell,8)*real(parini%nat,8)/Ha_kcalmol
  end subroutine
  end module interface_mopac
  
