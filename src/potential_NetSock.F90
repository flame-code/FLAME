!module interface_msock
!  use global
!  use defs_basis
!  USE F90SOCKETS, ONLY : create_socket, open_socket, writebuffer, readbuffer
!  use modsocket
!  implicit none
!  CHARACTER(LEN=MSGLEN) :: header
!  LOGICAL :: hasdata=.false.
!
!  private
!  public :: &
!    init_msock, &
!    evaluate_msock
! 
!contains

!  subroutine evaluate_msock(iproc,nat,latvec, xred, fcart, strten, energy, ka, kb, kc, iprec)
  subroutine cal_potential_forces_netsock(atoms)
!  use defs_basis
  USE F90SOCKETS, ONLY : create_socket, open_socket, writebuffer, readbuffer
  use mod_potential, only: sock_socket, sock_inet, sock_port,sock_host,MSGLEN,sock_extra_string,reset
  use mod_atoms, only: typ_atoms, update_ratp
  implicit none
  type(typ_atoms), intent(inout):: atoms
  integer:: iproc
  integer:: i,iat
  CHARACTER(LEN=MSGLEN) :: header
  LOGICAL :: hasdata=.false.
  real(8):: latvec(3,3)
  real(8):: xred(3,atoms%nat)
  real(8):: fcart(3,atoms%nat),energy,strten(6)
  real(8):: dist_ang(6), latvec_ang(3,3), latvec_rot(3,3), rotmat(3,3), vol
  character(1):: fn
!Socket buffers
    integer:: nmsg
    character*1024:: msg
    CHARACTER(LEN=20)   :: KPT_STRING
    CHARACTER(LEN=3)    :: KKA,KKB,KKC
    CHARACTER(LEN=6)    :: ECUTWFSCL
    integer:: cbuf, last_reset=1000
!Additional data can and will be passed down to the slave, as long as it requests a reinitiallization
!when obtaining the atomic position for the next force evaluation
!Check if kpt has been updated

!Check if precision has changed or any other challenges
if(reset) then
      sock_extra_string=trim(sock_extra_string)//" CRESET "
      last_reset=1000
endif
!Send message
    msg=trim(sock_extra_string)
    nmsg=60
!Copy the lattice vectors
    latvec=atoms%cellvec
!Transform to reduced coordinates
!    call rxyz_cart2int(latvec,xred,atoms%rat,atoms%nat)
    call update_ratp(atoms)
    call rxyz_cart2int_alborz(atoms%nat,latvec,atoms%ratp,xred)
    call send_data(xred,latvec,atoms%nat,atoms%nat,msg,nmsg,latvec_rot)
    call get_data(energy,fcart,strten,latvec,latvec_rot,atoms%nat)
!We have all info here, but we will only pass etot and fcart
    atoms%fat=fcart
    atoms%epot=energy
!Echo the forces and energy
    write(*,*) "Energy: ", energy
    write(*,*) "Forces: "
    do iat=1,atoms%nat
      write(*,'(3(es20.12))') fcart(:,iat)
    enddo
  end subroutine

  subroutine init_netsock(parini)
  USE F90SOCKETS, ONLY : create_socket, open_socket, writebuffer, readbuffer
  use mod_parini, only: typ_parini
  use mod_potential, only: sock_socket, sock_inet, sock_port,sock_host,MSGLEN,sock_extra_string
  implicit none
  type(typ_parini), intent(in):: parini
  character*1024:: host
      sock_inet=parini%inisock_inet
      sock_port=parini%inisock_port
      sock_host=parini%inisock_host
      write(*,'(a,a,a,i6)') " # SOCKET MASTER: Trying to open the socket ", trim(sock_host), " on port ",sock_port
      host = TRIM(sock_host)//achar(0)
      call create_socket(sock_socket, sock_inet, sock_port, host )
  end subroutine
  
  subroutine send_data(pos,latvec,nat,repid,msg,nmsg,latvec_rot)
  USE F90SOCKETS, ONLY : create_socket, open_socket, writebuffer, readbuffer
  use mod_potential, only: sock_socket, sock_inet, sock_port,sock_host,MSGLEN,sock_extra_string
  implicit none
  integer,intent(in):: nat, nmsg, repid
  real(8),intent(in):: latvec(3,3),pos(3,nat)
  real(8),intent(out)::latvec_rot(3,3)
  CHARACTER(LEN=MSGLEN) :: header
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
!Now the slave is ready to receive the positions
!Arrays need to be transposed, according to ipi, but in fact it would not be
!necessary since MHM is acting as a server here
!Prepare data for atomic positions
            pos_back=pos
            call backtocell_alborz(nat,latvec,pos_back)
!Rotate the cell to upper triangle
            latvec_rot=latvec
            call latvec2dist_ang(dist_ang,latvec,pi)
            call dist_ang2latvec(dist_ang,latvec_rot,pi)
!!!            call rxyz_int2cart(latvec_rot,pos_back,pos_cart,nat)
            call rxyz_int2cart_alborz(nat,latvec_rot,pos_back,pos_cart)
!Data must be re-formatted (and units converted) in the units and shapes used in the wrapper
            allocate(send_array(9))
            send_array=reshape(transpose(latvec_rot),(/9/))
            call writebuffer(sock_socket,send_array,9)               ! Writing the lattice vector
!            call invertmat_alborz (latvec_rot, latvec_inv, 3)
            call invertmat_alborz (latvec_rot, latvec_inv)
            send_array=reshape(transpose(latvec_inv),(/9/))
            call writebuffer(sock_socket,send_array,9)               ! Writing the inverse lattice vectors
            deallocate(send_array)  
            allocate(send_array(3*nat))
            send_array=reshape(pos_cart,(/3*nat/))
            call writebuffer(sock_socket,nat)                        ! Writing the number of atoms
            call writebuffer(sock_socket,send_array,3*nat)           ! Writing the positions
            deallocate(send_array)
            write(*,'(a)') " # SOCKET MASTER: positions and lattice sent"
!!  contains
!!       subroutine latvec2dist_ang(dist_ang,latvec,pi)
!!       !This subroutine will generate the distance and angles of a cell starting from latvec
!!       implicit none
!!       real(8):: dist_ang(6),latvec(3,3),pi,convang
!!       convang=180.d0/pi
!!       dist_ang(1)=sqrt(dot_product(latvec(:,1),latvec(:,1)))
!!       dist_ang(2)=sqrt(dot_product(latvec(:,2),latvec(:,2)))
!!       dist_ang(3)=sqrt(dot_product(latvec(:,3),latvec(:,3)))
!!       dist_ang(4)=convang*acos(dot_product(latvec(:,2),latvec(:,3))/(dist_ang(2)*dist_ang(3)))
!!       dist_ang(5)=convang*acos(dot_product(latvec(:,3),latvec(:,1))/(dist_ang(3)*dist_ang(1)))
!!       dist_ang(6)=convang*acos(dot_product(latvec(:,1),latvec(:,2))/(dist_ang(1)*dist_ang(2)))
!!       end subroutine
!**********************************************************************************************

  end subroutine

  subroutine get_data(etot,fcart,strten,latvec,latvec_rot,nat)
  USE F90SOCKETS, ONLY : create_socket, open_socket, writebuffer, readbuffer
  use mod_potential, only: sock_socket, sock_inet, sock_port,sock_host,MSGLEN,sock_extra_string
  implicit none
  integer,intent(in) :: nat
  real(8),intent(in) :: latvec(3,3),latvec_rot(3,3)
  real(8),intent(out):: fcart(3,nat),etot,strten(6)
  real(8):: rotmat(3,3),strmat(3,3),vol
  integer:: iat,nat_get
  real(8),allocatable:: get_array(:)
  integer:: nmsg
  character*1024:: msg
  CHARACTER(LEN=MSGLEN) :: header
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
    call rotmat_fcart_stress_other(latvec,latvec_rot,rotmat)
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
    call getvol_alborz(latvec_rot,vol)
    strten=strten/vol
!If cell is ortho    call rotate_stresstensor(strten,rotmat)
    write(*,'(a)') " # SOCKET MASTER: force, stress and energy received"
    !contains
  end subroutine
       !************************************************************************************
       subroutine rotmat_fcart_stress_other(latvec_init,latvec_trans,rotmat)
       !This subroutine will compute a rotation matrix, which transforms
       !fcart_trans into the original orientation forces fcart by fcart=matmul(rotmat,fcart_trans)
       !stress_trans into the original orientation stress by stress=rotmat*stress_trans*rotnat^T
       implicit none
       real(8):: latvec_init(3,3),latvec_trans(3,3),latvec_trans_inv(3,3),rotmat(3,3)
!       call invertmat_alborz(latvec_trans,latvec_trans_inv,3)
       call invertmat_alborz(latvec_trans,latvec_trans_inv)
       rotmat=matmul(latvec_init,latvec_trans_inv)
       end subroutine rotmat_fcart_stress_other
       !**********************************************************************************************
       subroutine rotate_stresstensor_other(strten,rotmat)
       !This subroutine will rotate the stress tensor by rotmat according to rotmat*stress*rotmat^T
       implicit none
       real(8):: strten(6),rotmat(3,3),stress(3,3)
               stress(1,1) =  strten(1) 
               stress(2,2) =  strten(2) 
               stress(3,3) =  strten(3) 
               stress(2,1) =  strten(6) 
               stress(3,1) =  strten(5) 
               stress(3,2) =  strten(4) 
               stress(1,2) =  stress(2,1)
               stress(1,3) =  stress(3,1)
               stress(2,3) =  stress(3,2)
                  stress=matmul(rotmat,matmul(stress,transpose(rotmat)))
               strten(1) =  stress(1,1)
               strten(2) =  stress(2,2)
               strten(3) =  stress(3,3)
               strten(6) =  stress(2,1)
               strten(5) =  stress(3,1)
               strten(4) =  stress(3,2)
       end subroutine rotate_stresstensor_other

  subroutine final_netsock()
  USE F90SOCKETS, ONLY : create_socket, open_socket, writebuffer, readbuffer
  use mod_potential, only: sock_socket, sock_inet, sock_port,sock_host,MSGLEN,sock_extra_string
  implicit none
  character*1024:: host
       call writebuffer(sock_socket,"STOP        ",MSGLEN)
  end subroutine
