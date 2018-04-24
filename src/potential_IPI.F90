module interface_ipi
  use global
  use defs_basis
  USE F90SOCKETS, ONLY : open_socket, writebuffer, readbuffer
  use modsocket

  implicit none
      ! SOCKET VARIABLES
!      INTEGER, PARAMETER :: MSGLEN=12   ! length of the headers of the driver/wrapper communication protocol
!      INTEGER socket, inet, port        ! socket ID & address of the server
!      CHARACTER(LEN=1024) :: host

      ! SOCKET COMMUNICATION BUFFERS
      CHARACTER(LEN=12) :: header
      LOGICAL :: isinit=.false., hasdata=.false.
      INTEGER cbuf


  private
  public :: &
    init_ipi, &
    evaluate_ipi
 

contains

  subroutine init_ipi(parini,nat)
  use mod_parini, only: typ_parini
  implicit none
  type(typ_parini), intent(in):: parini
  integer, intent(in):: nat
  real(8):: vir_tmp(3*3),vir_tmp_inv(3*3),msgbuffer(3*nat)
  CHARACTER(LEN=2048) :: initbuffer      ! it's unlikely a string this large will ever be passed...
  integer:: cbuf,i
  CHARACTER(LEN=1024) :: sock_host_open
      ! Calls the interface to the POSIX sockets library to open a communication channel
!      ipi_inet=0 !0 for unix socket, 1 for tcp
!      ipi_port=3141
!      ipi_host="mh-driver"//achar(0)
      write(*,'(a,a,a,i6)') " # IPI: Trying to open the socket ", trim(parini%sock_host), " on port ",parini%sock_port
      sock_host_open=trim(adjustl(parini%sock_host))//achar(0)
      CALL open_socket(sock_socket, parini%sock_inet, parini%sock_port, sock_host_open)
!Dummy call to check mpi
!Get the status of ipi to receive data
    do
      call readbuffer(sock_socket, header, MSGLEN)
      if (trim(header)/='STATUS') exit
      call writebuffer(sock_socket,"READY       ",MSGLEN)      
    enddo
    if (trim(header) == "POSDATA") then 
    !receives the positions & the cell data
    ! first reads cell and the number of atoms
    call readbuffer(sock_socket, vir_tmp, 9)
    call readbuffer(sock_socket, vir_tmp_inv, 9)
    call readbuffer(sock_socket, cbuf)
    !Check if the number of atoms coincide!
    if(nat.ne.cbuf) stop "Number of atoms not the same as in i-PI"
    CALL readbuffer(sock_socket, msgbuffer, nat*3)
    else
    stop "Did not get data"
    endif
  end subroutine



  subroutine evaluate_ipi(parini,nat,latvec, xred0, fcart, strten, energy, ka, kb, kc, iprec)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
    integer, intent(in):: nat
    real(8), intent(in) :: latvec(3,3)
    real(8), intent(in) :: xred0(3,nat)
    integer, intent(inout) :: ka, kb, kc
    integer, intent(in) :: iprec
    real(8), intent(out):: fcart(3,nat),energy,strten(6)
    real(8):: dist_ang(6), latvec_ang(3,3), latvec_rot(3,3), rxyz(3,nat), rotmat(3,3), vol, xred(3,nat)
    integer:: i,iat,itype
    integer:: nat_type(parini%ntypat_global)
    character(1):: fn
    logical:: isinit=.true.
!Ipi stuff
    CHARACTER(LEN=2048) :: initbuffer      ! it's unlikely a string this large will ever be passed...
    CHARACTER(LEN=20)   :: KPT_STRING
    CHARACTER(LEN=3)    :: KKA,KKB,KKC
    CHARACTER(LEN=6)    :: ECUTWFSCL
    real(8):: msgbuffer(3*nat), mtxbuf(9), virial(3,3),inverse(3,3)
    integer:: cbuf, last_reset=1000
    integer, save:: k_ipi(3)=0,k_old_ipi(3)=0,iprec_ipi

!Check if kpt has been updated
if((k_ipi(1)/=ka.or.k_ipi(2)/=kb.or.k_ipi(3)/=kc).and.last_reset.gt.10) then
   k_ipi(1)=max(ka,k_ipi(1))
   k_ipi(2)=max(kb,k_ipi(2))
   k_ipi(3)=max(kc,k_ipi(3))
   if(any(k_ipi.ne.k_old_ipi)) then
   write(KKA,'(i3)') k_ipi(1)
   write(KKB,'(i3)') k_ipi(2)
   write(KKC,'(i3)') k_ipi(3)
      KPT_STRING="UPKPT"//KKA//KKB//KKC
      sock_extra_string=trim(sock_extra_string)//trim(KPT_STRING)//" "
   last_reset=0
   k_old_ipi=k_ipi
   endif
endif

!Check if precision has changed
if(iprec/=iprec_ipi) then
      sock_extra_string=trim(sock_extra_string)//" CRESET "
      write(ECUTWFSCL,'(f6.4)') parini%sock_ecutwf(iprec)
      sock_extra_string=trim(sock_extra_string)//" ECUTWF "//ECUTWFSCL
      iprec_ipi=iprec
      last_reset=1000
      k_ipi=0
      k_old_ipi=0
endif

xred=xred0
call backtocell(nat,latvec,xred)
!Convert everything from "internal atomic" units to the real "angstrom-ev" units
if(any(parini%znucl(:).gt.200)) then   !When calling the driver software of Ceriotti
latvec_ang=latvec*Bohr_Ang !Only for lennard jones
else
latvec_ang=latvec
endif

!Rotate the cell to upper triangle
call latvec2dist_ang(dist_ang,latvec_ang,pi)
call dist_ang2latvec(dist_ang,latvec_rot,pi)
call rxyz_int2cart(latvec_rot,xred,rxyz,nat)
! Data must be re-formatted (and units converted) in the units and shapes used in the wrapper
            DO i = 1, nat
               msgbuffer(3*(i-1)+1:3*i) = rxyz(:,i)
            ENDDO
            virial = transpose(latvec_rot)

!Get the status of ipi and send the data
    do
      call readbuffer(sock_socket, header, MSGLEN)
      if (trim(header)/='STATUS') exit
      call writebuffer(sock_socket,"HAVEDATA    ",MSGLEN)
    enddo
    if (trim(header)/='GETFORCE') stop 'Error in socket communication!'
!Now that ipi is ready to get the data, send them
!The positions are ready, write into the buffer to send it to the ipi
            CALL writebuffer(sock_socket,"FORCEREADY  ",MSGLEN)
            CALL writebuffer(sock_socket,0.d0)  ! Writing the potential
            CALL writebuffer(sock_socket,nat)  ! Writing the number of atoms
            CALL writebuffer(sock_socket,msgbuffer,3*nat) ! Writing the forces
            CALL writebuffer(sock_socket,reshape(virial,(/9/)),9)  ! Writing the virial tensor, NOT divided by the volume
            cbuf = 60 ! Size of the "extras" string
            CALL writebuffer(sock_socket,cbuf) ! This would write out the "extras" string, but in this case we only use a dummy string.
  !          sock_extra_string="No info from MHM"
            CALL writebuffer(sock_socket,sock_extra_string,cbuf) !Sends and extra string, such as to reset the cell etc

!Get the status of ipi to receive data
    do
      call readbuffer(sock_socket, header, MSGLEN)
      if (trim(header)/='STATUS') exit
      call writebuffer(sock_socket,"READY       ",MSGLEN)      
    enddo
    if (trim(header) == "POSDATA") then 
    !receives the positions & the cell data
    ! first reads cell and the number of atoms
    CALL readbuffer(sock_socket, mtxbuf, 9)  ! Cell matrix
    virial = RESHAPE(mtxbuf, (/3,3/))
    CALL readbuffer(sock_socket, mtxbuf, 9)  ! Cell matrix
    inverse = RESHAPE(mtxbuf, (/3,3/))
    call readbuffer(sock_socket, cbuf)
    CALL readbuffer(sock_socket, msgbuffer, nat*3)
    else
    stop "Did not get data"
    endif
    DO i = 1, nat
        fcart(:,i) = msgbuffer(3*(i-1)+1:3*i)
    ENDDO
!Rotate forces back to original cell
    call rotmat_fcart_stress(latvec_ang,latvec_rot,rotmat)
    do iat=1,nat
       fcart(:,iat)=matmul(rotmat,fcart(:,iat))
    enddo
write(*,*) virial
write(*,*) fcart
if(any(parini%znucl(:).gt.200)) fcart=fcart/Ha_eV*Bohr_Ang
energy=virial(1,3)
write(*,*) energy
if(any(parini%znucl(:).gt.200)) energy=energy/Ha_eV
strten(1) = -virial(1,1)
strten(2) = -virial(2,2)
strten(3) = -virial(3,3)
strten(6) = -virial(2,1)
strten(5) = -virial(3,1)
strten(4) = -virial(3,2)
call getvol(latvec_rot,vol)
strten=strten/vol
if(any(parini%znucl(:).gt.200)) strten=strten/Ha_eV*Bohr_Ang**3
call rotate_stresstensor(strten,rotmat)
last_reset=last_reset+1
  end subroutine evaluate_ipi

  end module interface_ipi
