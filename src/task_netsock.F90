!*****************************************************************************************
subroutine netsock_task(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr, typ_file_info, set_ndof, update_rat
    use mod_potential, only: fcalls, perfstatus, potential
    use mod_processors, only: iproc
    use mod_const, only: ev2ha, ang2bohr, bohr2ang
    USE F90SOCKETS, ONLY : create_socket, open_socket, writebuffer, readbuffer
    use mod_potential, only: sock_socket, sock_inet, sock_port,sock_host,MSGLEN,sock_extra_string,reset
    use mod_acf, only: acf_read_new
    use mod_yaml_conf, only: read_yaml_conf
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variables
    type(typ_atoms_arr):: atoms_arr
    type(typ_file_info):: file_info
    real(8):: tt1, tt2, fxyz(3)
    integer:: iconf, iat
    logical:: hasdata,isinit
    logical:: yaml_exists,acf_exists,vasp_exists
    character(len=MSGLEN):: header
    character(len=60)::     msg
    integer :: nmsg,repid,str_index,PsiId,nat_get
    real(8) :: latvec(3,3),strmat(3,3),etot,vol
    real(8),allocatable:: pos(:,:),fcart(:,:)
    isinit = .true. 
    hasdata= .false.
!Open sockets if we want to communicate with server
    if(parini%usesocket) then
       sock_inet=parini%inisock_inet
       sock_port=parini%inisock_port
       sock_host=parini%inisock_host
       sock_host = TRIM(sock_host)//achar(0)
       call open_socket(sock_socket, sock_inet, sock_port, sock_host )
    endif
!----------------------------
    inquire(file='posinp.yaml',exist=yaml_exists)
    inquire(file='posinp.acf',exist=acf_exists)
    inquire(file='POSCUR',exist=vasp_exists)
    if(yaml_exists) then
        call read_yaml_conf(parini,'posinp.yaml',10000,atoms_arr)
    elseif(acf_exists) then
        call acf_read_new(parini,'posinp.acf',10000,atoms_arr)
    elseif(vasp_exists) then
        atoms_arr%nconf=1
        allocate(atoms_arr%atoms(atoms_arr%nconf))
        call read_poscar_for_single_point(parini,atoms_arr%atoms(1))
    else
        stop "A structure file of the atomic configuration, &
            &including types,&
            &must be provided when FLAME is run as a socket client"
    endif
!----------------------------
!If sockets are use, only run if a single configuration is present
    if(parini%usesocket.and.atoms_arr%nconf.gt.1) stop "Only single configurations allowed with socket communication"
!If sockets are used and netsock is the potential the code will bomb
    if(parini%usesocket.and.trim(potential)=='netsock') stop "Potential cannot be netsock if socket communication requested"

    do iconf=1,atoms_arr%nconf
        call set_ndof(atoms_arr%atoms(iconf))
    enddo
    potential=trim(parini%potential_potential)
    file_info%filename_positions='posout.acf'
    file_info%print_force=parini%print_force_single_point
    file_info%file_position='new'
    if(trim(parini%frmt_single_point)/='unknown') then
        file_info%frmt=trim(parini%frmt_single_point)
    endif
    do iconf=1,atoms_arr%nconf
        if(trim(potential)/='netsock' .or. iconf==1) then 
            call init_potential_forces(parini,atoms_arr%atoms(iconf))
        endif
        allocate(pos(3,atoms_arr%atoms(iconf)%nat),fcart(3,atoms_arr%atoms(iconf)%nat))
9001    continue !Socket loop
!receive-send iteration
        if(parini%usesocket) then 
            call readbuffer(sock_socket, header, MSGLEN)
            write(*,'(i6,a33,a100)') iproc, ' # SOCKET SLAVE: header received ',trim(header)
            if (trim(header) == "STATUS") then 
                call send_status(header, MSGLEN, isinit)
            else if (trim(header) == "INIT") then 
                call get_init   (header, MSGLEN, repid, isinit)
                isinit=.true.
!Although the number of atoms and the arrays assiciated with that side should be
!allocated already, we allocate fcart and pos here         
                if(atoms_arr%atoms(iconf)%nat.ne.repid) stop "The number of atoms on the master and slave side should be the same"
!Check if cell should be reset, in this case forget previous WF
                str_index = index(msg(1:msglen),"CRESET",.true.)
                if(str_index.gt.0) then
                    write(*,*) " MESSAGE: CRESET"
                    PsiId=0
                else
                    PsiId=1
                endif
            else if (trim(header) == "POSDATA") then 
                call get_data_slave(pos,latvec,atoms_arr%atoms(iconf)%nat,nat_get)
                atoms_arr%atoms(iconf)%cellvec=latvec
                if(potential=='vcblj')atoms_arr%atoms(iconf)%cellvec=latvec*bohr2ang
!Convert the units of positions and lattice vectors in a format that alborz will understand
                call rxyz_int2cart_alborz(atoms_arr%atoms(iconf)%nat,atoms_arr%atoms(iconf)%cellvec,pos,atoms_arr%atoms(iconf)%ratp)
                call update_rat(atoms_arr%atoms(iconf),upall=.true.)
!Compute the forces and stress here!!!
                call cal_potential_forces(parini,atoms_arr%atoms(iconf))
!Get Force
                fcart=atoms_arr%atoms(iconf)%fat
                if(potential=='vcblj')fcart=atoms_arr%atoms(iconf)%fat*ev2ha/ang2bohr
!Get Stress
                call getvol_alborz(atoms_arr%atoms(iconf)%cellvec,vol)
                strmat=atoms_arr%atoms(iconf)%stress!*vol
                if(potential=='vcblj')strmat=-atoms_arr%atoms(iconf)%stress*ev2ha*vol
!Get Energy
                etot=atoms_arr%atoms(iconf)%epot
                if(potential=='vcblj')etot=atoms_arr%atoms(iconf)%epot*ev2ha
!Now I have the data
                hasdata=.true.
            else if (trim(header)=="GETFORCE") then 
                nmsg = 0
                call send_data_slave(etot,fcart,strmat,atoms_arr%atoms(iconf)%nat,nmsg,msg)
                isinit = .false. ! resets init so that we will get replica index again at next step!
                hasdata= .false.
            elseif (trim(header)=="STOP") then 
                goto 9002
            elseif (trim(header)=="WAIT") then 
                goto 9001
            endif
            goto 9001
        else
            stop "Task netsock needs to be set"
        endif 
        if(parini%usesocket) goto 9001
        deallocate(pos,fcart)
    enddo
9002 continue
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
            call writebuffer(sock_socket,header,MSGLEN)               
            write(*,'(a,a)')   " # SOCKET SLAVE: header sent ", trim(header)
  end subroutine
  subroutine get_init(header, MSGLEN, repid, isinit)
  !Get initiallization string plus a repid
  implicit none
  integer:: MSGLEN,repid
  character(MSGLEN):: header
  logical:: isinit
         write(*,'(a)')   " # SOCKET SLAVE: initiallizing... "
         call readbuffer(sock_socket, repid)        ! actually this reads the replica id         
         call readbuffer(sock_socket, nmsg)  ; write(*,*) nmsg       ! length of parameter string -- ignored at present!
         call readbuffer(sock_socket, msg, nmsg)    ! the actual message
         write(*,'(a,a)') " # SOCKET SLAVE: initiallization string ", msg(1:nmsg)
         isinit=.true.
  end subroutine
  subroutine get_data_slave(pos,latvec,nat,nat_get)
  ! Receives the positions & the cell data
  implicit none
  integer:: nat,nat_get,i
  real(8):: pos(3,nat),pos_cart(3,nat),latvec(3,3),latvec_inv(3,3)
  real(8),allocatable:: get_array(:)
! first reads cell and the number of atoms
         write(*,'(a)')   " # SOCKET SLAVE: waiting for positions "
         allocate(get_array(9))
         call readbuffer(sock_socket, get_array , 9)
         latvec = transpose(RESHAPE(get_array, (/3,3/)))       !cell vector      
         call readbuffer(sock_socket, get_array, 9)
         latvec_inv = transpose(RESHAPE(get_array, (/3,3/)))   !inverse cell vector
         deallocate(get_array)
         call readbuffer(sock_socket, nat_get)                      !number of atoms
         if (nat.ne.nat_get) stop "Received NAT not the same as the &
         local NAT"
         allocate(get_array(3*nat_get))
         call readbuffer(sock_socket, get_array, nat_get*3)
         pos_cart = RESHAPE(get_array, (/ 3 , nat /) ) 
!         call rxyz_cart2int(latvec,pos,pos_cart,nat)
         call rxyz_cart2int_alborz(nat,latvec,pos_cart,pos)
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
  subroutine send_data_slave(etot,fcart,strmat,nat,nmsg,msg)
  ! Sends the energy, forces and stresses
  implicit none
  integer::nat,nmsg
  real(8):: fcart(3,nat),strmat(3,3),etot
  real(8),allocatable :: get_array(:)
  REAL(8),allocatable :: send_array(:)
  character(len=MSGLEN):: msg
 ! communicates energy info back to master
         write(*,'(a)')   " # SOCKET SLAVE: sending energy, forces and stress "
         call writebuffer(sock_socket,"FORCEREADY  ",MSGLEN)         
         call writebuffer(sock_socket,etot)
         allocate(send_array(3*nat))
         send_array=reshape(fcart,(/3*nat/))
         call writebuffer(sock_socket,nat)            
         call writebuffer(sock_socket,send_array,3*nat)
         deallocate(send_array)
         allocate(send_array(9))
         send_array=reshape(strmat,(/9/))
         call writebuffer(sock_socket,send_array,9)
         deallocate(send_array)
         ! i-pi can also receive an arbitrary string, that will be printed out to the "extra" 
         ! trajectory file. this is useful if you want to return additional information, e.g.
         ! atomic charges, wannier centres, etc. one must return the number of characters, then
         ! the string. here we just send back zero characters.
         call writebuffer(sock_socket,nmsg)
         if(nmsg.gt.0) then
                 call writebuffer(sock_socket,msg,nmsg)
         endif
  end subroutine
  subroutine str2arr(string,strarr,n)
  implicit none
  integer:: n,i
  character(n)::string
  character(1)::strarr(n)
  do i=1,n
    strarr(i)=string(i:i)
  enddo
  end subroutine
  subroutine arr2str(string,strarr,n)
  implicit none
  integer:: n,i
  character(n)::string
  character(1)::strarr(n)
  do i=1,n
    string(i:i)=strarr(i)
  enddo
  end subroutine
!*****************************************************************************************
subroutine read_poscar_for_single_point(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, atom_allocate_old, update_rat
    use global, only: units
    implicit none
    type(typ_parini), intent(inout):: parini !poscar_getsystem must be called from parser
    !so parini can be intent(in) in future.
    type(typ_atoms):: atoms
    !local variables
    real(8), allocatable:: xred(:,:)
    real(8), allocatable:: fcart(:,:)
    logical, allocatable:: fixat(:)
    integer, allocatable:: fragarr(:)
    real(8):: strten(6), printval1, printval2
    logical:: fixlat(7), readfix, readfrag
    integer:: iat
    character(40):: filename
    logical:: file_exists
    filename='posinp.vasp'
    inquire(file=trim(filename),exist=file_exists)
    if (.not. file_exists) then
        filename='POSCAR'
        inquire(file=trim(filename),exist=file_exists)
        if (.not. file_exists) stop "VASP file not found"
    endif
    write (*,*) "Reading structure from ",trim(filename)

    call poscar_getsystem(parini,trim(filename))
    allocate(xred(3,parini%nat),source=0.d0)
    allocate(fcart(3,parini%nat),source=0.d0)
    if(.not.allocated(fixat)) allocate(fixat(parini%nat),source=.false.)
    if(.not.allocated(fragarr)) allocate(fragarr(parini%nat),source=0)
    atoms%nat=parini%nat
    atoms%boundcond='bulk'
    call atom_allocate_old(atoms,parini%nat,0,0)
    call read_atomic_file_poscar(filename,atoms%nat,units,xred,atoms%cellvec,fcart,strten, &
        fixat,fixlat,readfix,fragarr,readfrag,printval1,printval2)
    call rxyz_int2cart_alborz(atoms%nat,atoms%cellvec,xred,atoms%ratp)
    call update_rat(atoms,upall=.true.)
    do iat=1,parini%nat
        atoms%sat(iat)=trim(parini%char_type(parini%typat_global(iat)))
    enddo
    deallocate(fixat)
    deallocate(xred)
    deallocate(fcart)
    deallocate(fragarr)
    if(allocated(parini%znucl)) deallocate(parini%znucl)
    if(allocated(parini%char_type)) deallocate(parini%char_type)
    if(allocated(parini%amu)) deallocate(parini%amu)
    if(allocated(parini%rcov)) deallocate(parini%rcov)
    if(allocated(parini%typat_global)) deallocate(parini%typat_global)
end subroutine read_poscar_for_single_point
!*****************************************************************************************
end subroutine netsock_task
!************************************************************************************
