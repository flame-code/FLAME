program f2fslave

! Test fortran-to-fortran socket communication
! JMS, Mar.2015

  use f90sockets, only: open_socket, writebuffer, readbuffer
  implicit none
  integer,parameter    :: MSGLEN=12
  ! this should be the IP address or fully-qualified host name on which the slave should connect
  ! in order to use UNIX sockets, this should be the name of the socket.
  ! these (and port in case of TCP sockets) should match the address and port on which the master is running
  character(len=1032)  :: host='127.0.0.1'
  integer              :: inet = 1   ! this has to be one to use a TCP/IP socket, zero to have a UNIX socket
  integer              :: port=21211 ! this is the port number (only used for TCP/IP sockets) 
  integer              :: socket
  character(len=MSGLEN):: header
  integer              :: nat,repid
  !Controlling variables
  CHARACTER*1024 :: msg
  INTEGER :: nmsg
  INTEGER:: nat_get
  LOGICAL :: isinit=.false., hasdata=.false., firststep=.true.
  !Lattice vectors and others  
  REAL*8  :: latvec(3,3), latvec_inv(3,3), strmat(3,3),etot
  REAL*8, ALLOCATABLE :: fcart(:,:),pos(:,:),get_array(:),send_array(:)


  read(*,*) port
  host = TRIM(host)//achar(0)
  call open_socket( socket, inet, port, host )
  
  do    ! receive-send iteration
    call readbuffer(socket, header, MSGLEN)
    write(*,'(a,a)') ' # SOCKET SLAVE: header received ',trim(header)
      if (trim(header) == "STATUS") then
         call send_status(header, MSGLEN, isinit)
      else if (trim(header) == "INIT") then
         call get_init   (header, MSGLEN, repid, isinit)
!Although the number of atoms and the arrays assiciated with that side should be
!allocate already, we allocate fcart and pos here         
         nat=repid
         if (.not.allocated(fcart)) then
           allocate(fcart(3,nat))
         endif
         if (.not.allocated(pos)) then
           allocate(pos(3,nat))
         endif
      else if (trim(header) == "POSDATA") then
         call get_data(pos,latvec,nat,nat_get)
!!!Compute the forces and stress here!!!
         strmat=0.d0
         fcart=0.d0
         etot=0.d0
         hasdata=.true.
      else if (trim(header)=="GETFORCE") then
         nmsg = 0
         call send_data(fcart,strmat,nat,nmsg,msg)
         isinit = .false. ! resets init so that we will get replica index again at next step!
         hasdata= .false.
    elseif (trim(header)=="STOP") then
      exit
    elseif (trim(header)=="WAIT") then
      cycle
    endif
  enddo

contains
  subroutine send_status(header, MSGLEN, isinit)
  !Report the status to the master
  implicit none
  integer:: MSGLEN
  character(MSGLEN):: header
  logical:: isinit
            if (hasdata) then
               header="HAVEDATA    "
            else if (isinit) then
               header="READY       "
            else 
               header="NEEDINIT    "
            endif
            call writebuffer(socket,header,MSGLEN)               
            write(*,'(a,a)')   " # SOCKET SLAVE: header sent ", trim(header)
  end subroutine
  subroutine get_init(header, MSGLEN, repid, isinit)
  !Get initiallization string plus a repid
  implicit none
  integer:: MSGLEN,repid
  character(MSGLEN):: header
  logical:: isinit
         write(*,'(a)')   " # SOCKET SLAVE: initiallizing... "
         call readbuffer(socket, repid)        ! actually this reads the replica id         
         call readbuffer(socket, nmsg)         ! length of parameter string -- ignored at present!
         call readbuffer(socket, msg, nmsg)    ! the actual message
         write(*,'(a,a)') " # SOCKET SLAVE: initiallization string ", msg(1:nmsg)
         isinit=.true.
  end subroutine
  subroutine get_data(pos,latvec,nat,nat_get)
  ! Receives the positions & the cell data
  implicit none
  integer:: nat,nat_get,i
  real(8):: pos(3,nat),pos_cart(3,nat),latvec(3,3)
  real(8),allocatable:: get_array(:)
! first reads cell and the number of atoms
         write(*,'(a)')   " # SOCKET SLAVE: waiting for positions "
         allocate(get_array(9))
         call readbuffer(socket, get_array , 9)
         latvec = transpose(RESHAPE(get_array, (/3,3/)))       !cell vector      
         call readbuffer(socket, get_array, 9)
         latvec_inv = transpose(RESHAPE(get_array, (/3,3/)))   !inverse cell vector
         deallocate(get_array)
         call readbuffer(socket, nat_get)                      !number of atoms
         if (nat.ne.nat_get) stop "Received NAT not the same as the &
         local NAT"
         allocate(get_array(3*nat_get))
         call readbuffer(socket, get_array, nat_get*3)
         pos_cart = RESHAPE(get_array, (/ 3 , nat /) ) 
         call rxyz_cart2int(latvec,pos,pos_cart,nat)
         deallocate(get_array)
         write(21,*) 'lat'
         write(21,*)  latvec(:,1)  
         write(21,*)  latvec(:,2)  
         write(21,*)  latvec(:,3)  
         write(21,*) 'pos'
         do i=1,nat
            write(21,*) pos(:,i)
         enddo

         write(*,'(a)')   " # SOCKET SLAVE: received positions "
  end subroutine
  subroutine send_data(fcart,strmat,nat,nmsg,msg)
  ! Sends the energy, forces and stresses
  implicit none
  integer::nat,nmsg
  real(8):: fcart(3,nat),strmat(3,3)
  real(8),allocatable:: get_array(:)
  CHARACTER*1024 :: msg
 ! communicates energy info back to master
         write(*,'(a)')   " # SOCKET SLAVE: sending energy, forces and stress "
         call writebuffer(socket,"FORCEREADY  ",MSGLEN)         
         call writebuffer(socket,etot)
         allocate(send_array(3*nat))
         call writebuffer(socket,nat)            
         call writebuffer(socket,send_array,3*nat)
         deallocate(send_array)
         allocate(send_array(9))
         send_array=reshape(strmat,(/9/))
         call writebuffer(socket,send_array,9)
         deallocate(send_array)
         ! i-pi can also receive an arbitrary string, that will be printed out to the "extra" 
         ! trajectory file. this is useful if you want to return additional information, e.g.
         ! atomic charges, wannier centres, etc. one must return the number of characters, then
         ! the string. here we just send back zero characters.
         call writebuffer(socket,nmsg)
         if(nmsg.gt.0) then
                 call writebuffer(socket,msg,nmsg)
         endif
  end subroutine
end program f2fslave

!************************************************************************************

 subroutine rxyz_cart2int(latvec,rxyzint,rxyzcart,nat)
 !This subrouine will convert the internal coordinates into cartesian coordinates
 implicit none
 real(8):: rxyzint(3,nat), rxyzcart(3,nat),latvec(3,3),latvecinv(3,3)
 integer:: nat,iat
 call invertmat(latvec,latvecinv,3)
 do iat=1,nat
  rxyzint(:,iat)=matmul(latvecinv,rxyzcart(:,iat))
 enddo
 end subroutine rxyz_cart2int

!************************************************************************************

 subroutine invertmat(mat,matinv,n)
 implicit none
 real(8),intent(in) :: mat(n,n)
 integer               :: n
 real(8),allocatable   :: WORK(:)
 real(8)               :: matinv(n,n),det(3),a(n,n),div
 integer               :: IPIV(n), INFO
 integer               :: LDWORK
 !Here only for a 3*3 matrix
 if (n==3) then
 a=mat
 div=(a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)+a(1,2)*a(2,3)*a(3,1)+&
 &a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)) 
 div=1.d0/div
      matinv(1,1) = (a(2,2)*a(3,3)-a(2,3)*a(3,2))*div
      matinv(1,2) =-(a(1,2)*a(3,3)-a(1,3)*a(3,2))*div
      matinv(1,3) = (a(1,2)*a(2,3)-a(1,3)*a(2,2))*div
      matinv(2,1) =-(a(2,1)*a(3,3)-a(2,3)*a(3,1))*div
      matinv(2,2) = (a(1,1)*a(3,3)-a(1,3)*a(3,1))*div
      matinv(2,3) =-(a(1,1)*a(2,3)-a(1,3)*a(2,1))*div
      matinv(3,1) = (a(2,1)*a(3,2)-a(2,2)*a(3,1))*div
      matinv(3,2) =-(a(1,1)*a(3,2)-a(1,2)*a(3,1))*div
      matinv(3,3) = (a(1,1)*a(2,2)-a(1,2)*a(2,1))*div
 else
 !General n*n matrix 
 matinv=mat
 allocate(WORK(n))
 call  DGETRF( n, n, matinv, n, IPIV, INFO )
 if (info.ne.0) stop "Error in DGETRF"
 LDWORK=-1
 call  DGETRI( n, matinv, n, IPIV, WORK,LDWORK , INFO )
 LDWORK=WORK(1)
 deallocate(WORK)
 allocate(WORK(LDWORK))
 call  DGETRI( n, matinv, n, IPIV, WORK,LDWORK , INFO )
 if (info.ne.0) stop "Error in DGETRI"
 endif
 end subroutine

!************************************************************************************
