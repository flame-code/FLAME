module interface_msock
  use global
  use defs_basis
  USE F90SOCKETS, ONLY : create_socket, open_socket, writebuffer, readbuffer
  use modsocket
  implicit none
  CHARACTER(LEN=MSGLEN) :: header
  LOGICAL :: hasdata=.false.

  private
  public :: &
    init_msock, &
    evaluate_msock
 
contains

  subroutine evaluate_msock(parini,latvec, xred, fcart, strten, energy, ka, kb, kc, iprec)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
    real(8), intent(in) :: latvec(3,3)
    real(8), intent(in) :: xred(3,parini%nat)
    integer, intent(inout) :: ka, kb, kc
    integer, intent(in) :: iprec
    real(8), intent(out):: fcart(3,parini%nat),energy,strten(6)
    real(8):: dist_ang(6), latvec_ang(3,3), latvec_rot(3,3), rxyz(3,parini%nat), rotmat(3,3), vol
    integer:: i,iat,itype
    integer:: nat_type(parini%ntypat_global)
    character(1):: fn
!Socket buffers
    integer:: nmsg
    character*1024:: msg
    CHARACTER(LEN=20)   :: KPT_STRING
    CHARACTER(LEN=3)    :: KKA,KKB,KKC
    CHARACTER(LEN=6)    :: ECUTWFSCL
    integer:: cbuf, last_reset=1000
    integer, save:: k_msock(3)=0,k_old_msock(3)=0,iprec_msock
!Additional data can and will be passed down to the slave, as long as it requests a reinitiallization
!when obtaining the atomic position for the next force evaluation
!Check if kpt has been updated
if((k_msock(1)/=ka.or.k_msock(2)/=kb.or.k_msock(3)/=kc).and.last_reset.gt.10) then
    k_msock(1)=max(ka,k_msock(1))
    k_msock(2)=max(kb,k_msock(2))
    k_msock(3)=max(kc,k_msock(3))
    if(any(k_msock.ne.k_old_msock)) then
      write(KKA,'(i3)') k_msock(1)
      write(KKB,'(i3)') k_msock(2)
      write(KKC,'(i3)') k_msock(3)
      KPT_STRING="UPKPT"//KKA//KKB//KKC
      sock_extra_string=trim(sock_extra_string)//trim(KPT_STRING)//" "
      last_reset=0
      k_old_msock=k_msock
    endif
endif

!Check if precision has changed
if(iprec/=iprec_msock) then
      sock_extra_string=trim(sock_extra_string)//" CRESET "
      write(ECUTWFSCL,'(f6.4)') parini%sock_ecutwf(iprec)
      sock_extra_string=trim(sock_extra_string)//" ECUTWF "//ECUTWFSCL
      iprec_msock=iprec
      last_reset=1000
      k_msock=0
      k_old_msock=0
endif
    msg=trim(sock_extra_string)
    nmsg=60
    call send_data(xred,latvec,parini%nat,parini%nat,msg,nmsg,latvec_rot)
    call get_data(energy,fcart,strten,latvec,latvec_rot,parini%nat)

  end subroutine

  subroutine init_msock(parini)
  use mod_parini, only: typ_parini
  implicit none
  type(typ_parini), intent(in):: parini
  character*1024:: host
      write(*,'(a,a,a,i6)') " # SOCKET MASTER: Trying to open the socket ", trim(parini%sock_host), " on port ",parini%sock_port
      host = TRIM(parini%sock_host)//achar(0)
      call create_socket( sock_socket, parini%sock_inet, parini%sock_port, host )
  end subroutine
  
  subroutine send_data(pos,latvec,nat,repid,msg,nmsg,latvec_rot)
  implicit none
  integer,intent(in):: nat, nmsg, repid
  real(8),intent(in):: latvec(3,3),pos(3,nat)
  real(8),intent(out)::latvec_rot(3,3)
  real(8):: latvec_inv(3,3),pos_back(3,nat),pos_cart(3,nat),dist_ang(6)
  real(8),parameter :: pi = 3.141592653589793239d0
  real(8),allocatable:: send_array(:)
  character*1024:: msg
!Get the status of slave and send the data
  do
!Ask for the status
      call writebuffer(sock_socket,"STATUS      ",MSGLEN)
      call readbuffer (sock_socket, header, MSGLEN)
      write(*,'(a,a)') " # SOCKET MASTER: header returned ", trim(header)
      if (trim(header)=='READY       ') then
              exit
      elseif (trim(header)=='NEEDINIT') then
!Needs initiallization, send the string
              call writebuffer(sock_socket,'INIT        ',MSGLEN)
              call writebuffer(sock_socket,repid)    !Replica ID, not really used here
              call writebuffer(sock_socket,nmsg)     !Message size
              call writebuffer(sock_socket,msg,nmsg) !Message contents
              exit
      endif
  enddo
!Tell the slave to prepare for receiving the data
            call writebuffer(sock_socket,"POSDATA     ",MSGLEN)
            write(*,*) latvec
!Now the slave is ready to receive the positions
!Arrays need to be transposed, according to ipi, but in fact it would not be
!necessary since MHM is acting as a server here
!Prepare data for atomic positions
            pos_back=pos
            call backtocell(nat,latvec,pos_back)
!Rotate the cell to upper triangle
            call latvec2dist_ang(dist_ang,latvec,pi)
            call dist_ang2latvec(dist_ang,latvec_rot,pi)
            call rxyz_int2cart(latvec_rot,pos_back,pos_cart,nat)
!Data must be re-formatted (and units converted) in the units and shapes used in the wrapper
            allocate(send_array(9))
            send_array=reshape(transpose(latvec_rot),(/9/))
            call writebuffer(sock_socket,send_array,9)               ! Writing the lattice vector
            call invertmat (latvec_rot, latvec_inv, 3)
            send_array=reshape(transpose(latvec_inv),(/9/))
            call writebuffer(sock_socket,send_array,9)               ! Writing the inverse lattice vectors
            deallocate(send_array)  
            allocate(send_array(3*nat))
            send_array=reshape(pos_cart,(/3*nat/))
            call writebuffer(sock_socket,nat)                        ! Writing the number of atoms
            call writebuffer(sock_socket,send_array,3*nat)           ! Writing the positions
            deallocate(send_array)
            write(*,'(a)') " # SOCKET MASTER: positions and lattice sent"
  end subroutine

  subroutine get_data(etot,fcart,strten,latvec,latvec_rot,nat)
  implicit none
  integer,intent(in) :: nat
  real(8),intent(in) :: latvec(3,3),latvec_rot(3,3)
  real(8),intent(out):: fcart(3,nat),etot,strten(6)
  real(8):: rotmat(3,3),strmat(3,3),vol
  integer:: iat,nat_get
  real(8),allocatable:: get_array(:)
  integer:: nmsg
  character*1024:: msg
  write(*,'(a)') " # SOCKET MASTER: requesting energy, forces and stress"
!Get the status of slave and get the data
  do
!Ask for the status
      call writebuffer(sock_socket,"STATUS      ",MSGLEN)      
      call readbuffer (sock_socket, header,       MSGLEN)
      write(*,'(a,a)') " # SOCKET MASTER: header returned ", trim(header)
      if (trim(header)=='HAVEDATA') then
!The forces are ready, ask for them now
              call writebuffer(sock_socket,"GETFORCE    ",MSGLEN)    
              call readbuffer (sock_socket, header,       MSGLEN)
              if (trim(header)/='FORCEREADY') then
                      stop "Socket did not reply with FORCEREADY when sending &
                      &GETFORCE"
              endif
              exit
      endif
    enddo
    !Read potential energy
    call readbuffer(sock_socket, etot)
    !Read the array size, in terms of nat
    call readbuffer(sock_socket, nat_get)
    if (nat.ne.nat_get) stop "Array sizes not consistent"
    allocate(get_array(3*nat_get))
    call readbuffer(sock_socket, get_array, 3*nat_get)  !Getting the forces
    fcart=reshape(get_array,(/3,nat_get/))
    !Get the stress tensor in matrix form
    deallocate(get_array)
    allocate(get_array(9))
    call readbuffer(sock_socket, get_array, 9)  ! Getting stress tensor
    strmat = RESHAPE(get_array, (/3,3/))
    !Get additional comment, if any
    call readbuffer(sock_socket, nmsg)
    if(nmsg.gt.0) then
        call readbuffer(sock_socket, msg, nmsg)
    endif
    deallocate(get_array)
!Convert the forces and stresses with rotation matrix
!Rotate forces back to original cell
    call rotmat_fcart_stress(latvec,latvec_rot,rotmat)
    do iat=1,nat
       fcart(:,iat)=matmul(rotmat,fcart(:,iat))
    enddo
!Rotate stress back to original cell
    strten(1) = -strmat(1,1)
    strten(2) = -strmat(2,2)
    strten(3) = -strmat(3,3)
    strten(6) = -strmat(2,1)
    strten(5) = -strmat(3,1)
    strten(4) = -strmat(3,2)
    call getvol(latvec_rot,vol)
    strten=strten/vol
    call rotate_stresstensor(strten,rotmat)
    write(*,'(a)') " # SOCKET MASTER: force, stress and energy received"
  end subroutine



end module
  subroutine socket_stop()
  USE F90SOCKETS, ONLY : create_socket, open_socket, writebuffer, readbuffer
  use modsocket
  implicit none
  character*1024:: msg
!Send kill signal
      write(*,'(a)') " # SOCKET MASTER: Sending STOP command" 
      call writebuffer(sock_socket,"STOP        ",MSGLEN)
  end subroutine
