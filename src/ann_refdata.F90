!*****************************************************************************************
module mod_refdata
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr
    use mod_symfunc, only: typ_symfunc_arr
    use mod_ann, only: typ_ann_arr
    implicit none
    private
    public:: typ_refdata
    type typ_refdata
        logical:: initialized=.false.
        type(typ_atoms_arr):: atoms_train
        type(typ_atoms_arr):: atoms_valid
        type(typ_symfunc_arr):: symfunc_train
        type(typ_symfunc_arr):: symfunc_valid
        contains
        procedure, public, pass(self):: init_refdata
        procedure, public, pass(self):: fini_refdata
        procedure, public, pass(self):: read_refdata
        procedure, public, pass(self):: prepare_refdata
    end type typ_refdata
contains
!*****************************************************************************************
subroutine init_refdata(self)
    implicit none
    class(typ_refdata), intent(inout):: self
    !type(typ_atoms_arr), intent(in), target:: atoms_train
    !type(typ_atoms_arr), intent(in), target:: atoms_valid
    !type(typ_symfunc_arr), intent(in), target:: symfunc_train
    !type(typ_symfunc_arr), intent(in), target:: symfunc_valid
    !local variables
end subroutine init_refdata
!*****************************************************************************************
subroutine fini_refdata(self)
    implicit none
    class(typ_refdata), intent(inout):: self
    !type(typ_atoms_arr), intent(in), target:: atoms_train
    !type(typ_atoms_arr), intent(in), target:: atoms_valid
    !type(typ_symfunc_arr), intent(in), target:: symfunc_train
    !type(typ_symfunc_arr), intent(in), target:: symfunc_valid
    !local variables
end subroutine fini_refdata
!*****************************************************************************************
subroutine read_refdata(self,parini,ann_arr,nann)
    use mod_ann_io_yaml, only: read_data_yaml
    use mod_ann_io_yaml, only: read_input_ann_yaml, read_ann_yaml
    use mod_processors, only: iproc
    use yaml_output
    implicit none
    class(typ_refdata), intent(inout):: self
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    integer, intent(in):: nann
    !local variables
    logical:: file_exists
    character (50)::fname
    !Reading configurations and their energies and forces
    inquire(file="list_posinp_train.yaml",exist=file_exists)
    if(file_exists) then
        call read_data_yaml(parini,'list_posinp_train.yaml',self%atoms_train,ann_arr=ann_arr)
    else
        call read_data_old(parini,'list_posinp_train',self%atoms_train)
    endif
    inquire(file="list_posinp_valid.yaml",exist=file_exists)
    if(file_exists) then
        call read_data_yaml(parini,'list_posinp_valid.yaml',self%atoms_valid,ann_arr=ann_arr)
    else
        call read_data_old(parini,'list_posinp_valid',self%atoms_valid)
    endif
    call ann_arr%set_number_of_ann(nann)
    if(ann_arr%nann==0) stop 'ERROR: number of type of atoms zero in ann_train'
    call yaml_map('number of ann',ann_arr%nann)
    !write(*,*) 'Here', ann_arr%nann
    allocate(ann_arr%ann(ann_arr%nann))
    ann_arr%approach=trim(parini%approach_ann)
    fname=trim(parini%stypat(1))//'.ann.input.yaml'
    inquire(file=trim(fname),exist=ann_arr%exists_yaml_file)
    if(ann_arr%exists_yaml_file) then
        if(parini%restart_param) then
            call read_ann_yaml(parini,ann_arr)
        else
            call read_input_ann_yaml(parini,iproc,ann_arr)
        endif
    else
        call read_input_ann(parini,iproc,ann_arr)
    endif
end subroutine read_refdata
!*****************************************************************************************
subroutine prepare_refdata(self,parini,ann_arr)
    use mod_processors, only: iproc
    use yaml_output
    implicit none
    class(typ_refdata), intent(inout):: self
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    !local variables
    real(8):: time1, time2, time3
    !-------------------------------------------------------
    if(iproc==0) then
        call yaml_map('number of ANN wights',ann_arr%nwtot)
        call yaml_map('number of training data points',self%atoms_train%nconf)
        call yaml_map('number of validating data points',self%atoms_valid%nconf)
    endif
    call set_conf_inc_random(parini,self%atoms_train)
    call set_conf_inc_random(parini,self%atoms_valid)
    call prepare_atoms_arr(parini,self%atoms_train)
    call prepare_atoms_arr(parini,self%atoms_valid)
    !-------------------------------------------------------
    call cpu_time(time1)
    call set_gbounds(parini,ann_arr,self%atoms_valid,'bounds_valid',self%symfunc_valid)
    call cpu_time(time2)
    call set_gbounds(parini,ann_arr,self%atoms_train,'bounds_train',self%symfunc_train)
    call cpu_time(time3)
    !write(*,'(a,2f10.1)') 'TIMING: evaluation symmetry functions: ',time2-time1,time3-time2
    call yaml_map('TIMING evaluation symmetry functions for valid',time2-time1)
    call yaml_map('TIMING evaluation symmetry functions for train',time3-time2)
    !-------------------------------------------------------------------------------------
    !IMPORTANT: The following must be done after set_gbounds is called for training set.
    !if(trim(parini%symfunc)/='do_not_save') then
    call apply_gbounds_atom(parini,ann_arr,self%atoms_valid,self%symfunc_valid)
    call apply_gbounds_atom(parini,ann_arr,self%atoms_train,self%symfunc_train)
    !endif
end subroutine prepare_refdata
!*****************************************************************************************
subroutine prepare_atoms_arr(parini,atoms_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr, update_ratp, update_rat
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms_arr), intent(inout):: atoms_arr
    !local variables
    integer:: iconf, iat, i
    real(8), allocatable:: ratred(:,:)
    do iconf=1,atoms_arr%nconf
        if(trim(atoms_arr%atoms(iconf)%boundcond)=='bulk') then
        allocate(ratred(1:3,1:atoms_arr%atoms(iconf)%nat))
        call update_ratp(atoms_arr%atoms(iconf))
        call rxyz_cart2int_alborz(atoms_arr%atoms(iconf)%nat,atoms_arr%atoms(iconf)%cellvec,atoms_arr%atoms(iconf)%ratp,ratred)
        call backtocell_alborz(atoms_arr%atoms(iconf)%nat,atoms_arr%atoms(iconf)%cellvec,ratred)
        call rxyz_int2cart_alborz(atoms_arr%atoms(iconf)%nat,atoms_arr%atoms(iconf)%cellvec,ratred,atoms_arr%atoms(iconf)%ratp)
        call update_rat(atoms_arr%atoms(iconf),upall=.true.)
        deallocate(ratred)
        endif
        do iat=1,atoms_arr%atoms(iconf)%nat
            do i=1,parini%ntypat
                if(trim(atoms_arr%atoms(iconf)%sat(iat))==trim(parini%stypat(i))) then
                    atoms_arr%atoms(iconf)%itypat(iat)=parini%ltypat(i)
                    exit
                endif
            enddo
        enddo
    enddo
end subroutine prepare_atoms_arr
!*****************************************************************************************
subroutine set_conf_inc_random(parini,atoms_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr
    use mod_utils
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms_arr), intent(inout):: atoms_arr
    !local variables
    integer:: iconf, irand, ierr
    integer, allocatable:: ind_list(:)
    real(8):: tt
#if defined(MPI)
    include 'mpif.h'
#endif
    if(parini%nconf_rmse==0) then
        write(*,*) 'ERROR: parini%nconf_rmse=0'
        stop
    endif
    if(parini%nconf_rmse>=atoms_arr%nconf) then
        allocate(atoms_arr%conf_inc(atoms_arr%nconf),source=.true.)
        atoms_arr%nconf_inc=atoms_arr%nconf
        return
    endif
    allocate(atoms_arr%conf_inc(atoms_arr%nconf),source=.false.)
    atoms_arr%nconf_inc=parini%nconf_rmse
    allocate(ind_list(atoms_arr%nconf_inc))
    if(parini%mpi_env%iproc==0) then
        irand=0
        do
            if(irand==atoms_arr%nconf_inc) exit
            if(trim(parini%rng_type)=='only_for_tests') then
                call random_number_generator_simple(tt)
            else
                call random_number(tt)
            endif
            tt=tt*real(atoms_arr%nconf)
            iconf=int(tt)+1
            if(atoms_arr%conf_inc(iconf)) cycle
            atoms_arr%conf_inc(iconf)=.true.
            irand=irand+1
            ind_list(irand)=iconf
        enddo
    endif
    call MPI_BARRIER(parini%mpi_env%mpi_comm,ierr)
    call MPI_BCAST(ind_list,atoms_arr%nconf_inc,MPI_INTEGER,0,parini%mpi_env%mpi_comm,ierr)
    if(parini%mpi_env%iproc>0) then
        do irand=1,atoms_arr%nconf_inc
            atoms_arr%conf_inc(ind_list(irand))=.true.
        enddo
    endif
    deallocate(ind_list)
    call MPI_BARRIER(parini%mpi_env%mpi_comm,ierr)
end subroutine set_conf_inc_random
!*****************************************************************************************
subroutine set_gbounds(parini,ann_arr,atoms_arr,strmess,symfunc_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc_arr
    use mod_atoms, only: typ_atoms_arr
    use mod_linked_lists, only: typ_pia_arr
    use mod_processors, only: iproc, nproc
    use wrapper_MPI, only: fmpi_allreduce, FMPI_SUM
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms_arr), intent(inout):: atoms_arr
    character(*), intent(in):: strmess
    type(typ_symfunc_arr), intent(inout):: symfunc_arr
    !local variables
    integer:: iconf
    if(.not. allocated(symfunc_arr%symfunc)) then
        symfunc_arr%nconf=atoms_arr%nconf
        allocate(symfunc_arr%symfunc(symfunc_arr%nconf))
    endif
    !write(*,'(a,i3,i6)') 'iproc,nconf ',iproc,atoms_arr%nconf
    do iconf=1,atoms_arr%nconf
        call symfunc_arr%symfunc(iconf)%init_symfunc(parini%mpi_env,parini%iverbose,parini%bondbased_ann,parini%symfunc_type_ann)
        symfunc_arr%symfunc(iconf)%ng=ann_arr%ann(1)%nn(0) !HERE
        symfunc_arr%symfunc(iconf)%nat=atoms_arr%atoms(iconf)%nat
    enddo
    !configuration: do iconf=1+iproc,atoms_arr%nconf,nproc
    configuration: do iconf=1,atoms_arr%nconf
        if(trim(parini%symfunc)/='read') then
            call symfunc_arr%symfunc(iconf)%get_symfunc(ann_arr,atoms_arr%atoms(iconf),.false.)
            if(parini%mpi_env%nproc>1) then
                call fmpi_allreduce(symfunc_arr%symfunc(iconf)%y(1,1), &
                    ann_arr%ann(1)%nn(0)*atoms_arr%atoms(iconf)%nat,op=FMPI_SUM,comm=parini%mpi_env%mpi_comm)
            endif
            if(.not. parini%save_symfunc_force_ann) then
                call f_free(symfunc_arr%symfunc(iconf)%y0d)
            endif
            call f_free(symfunc_arr%symfunc(iconf)%y0dr)
            deallocate(symfunc_arr%symfunc(iconf)%linked_lists%prime_bound)
            deallocate(symfunc_arr%symfunc(iconf)%linked_lists%bound_rad)
            deallocate(symfunc_arr%symfunc(iconf)%linked_lists%bound_ang)
        !elseif(trim(parini%symfunc)/='read') then
        !    stop 'ERROR: arini%symfunc contains none of the three acceptable possibilies'
        endif
        if(trim(parini%symfunc)=='write') then
            call write_symfunc(parini,iconf,atoms_arr,strmess,symfunc_arr)
        elseif(trim(parini%symfunc)=='read') then
            call read_symfunc(parini,iconf,ann_arr,atoms_arr,strmess,symfunc_arr)
            deallocate(symfunc_arr%symfunc(iconf)%linked_lists%prime_bound)
            deallocate(symfunc_arr%symfunc(iconf)%linked_lists%bound_rad)
            deallocate(symfunc_arr%symfunc(iconf)%linked_lists%bound_ang)
        !elseif(trim(parini%symfunc)=='do_not_save') then
        !    call f_free(symfunc_arr%symfunc(iconf)%y0d)
        !    call f_free(symfunc_arr%symfunc(iconf)%y0dr)
        !    deallocate(symfunc_arr%symfunc(iconf)%linked_lists%prime_bound)
        !    deallocate(symfunc_arr%symfunc(iconf)%linked_lists%bound_ang)
        endif
    enddo configuration
    call save_gbounds(parini,ann_arr,atoms_arr,strmess,symfunc_arr)
end subroutine set_gbounds
!*****************************************************************************************
subroutine write_symfunc(parini,iconf,atoms_arr,strmess,symfunc_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc_arr
    use mod_atoms, only: typ_atoms_arr, update_ratp
    use mod_linked_lists, only: typ_pia_arr
    use mod_processors, only: iproc, nproc
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iconf
    type(typ_atoms_arr), intent(inout):: atoms_arr
    character(*), intent(in):: strmess
    type(typ_symfunc_arr), intent(inout):: symfunc_arr
    !local variables
    real(8), allocatable:: wa(:)
    integer:: nwa
    character(30):: filename
    integer:: i, ig, iat, ios, n, ib
    !Symmetry functions are written into files to be used for
    !subsequent training runs.
    if(trim(strmess)=='bounds_train') then
        write(filename,'(a25,i5.5)') '../symfunc/train.symfunc.',iconf
    elseif(trim(strmess)=='bounds_valid') then
        write(filename,'(a25,i5.5)') '../symfunc/valid.symfunc.',iconf
    else
        stop 'ERROR: invalid content in strmess in gset_bounds '
    endif
    open(unit=311,file=trim(filename),status='replace',form='unformatted', &
        access='stream',iostat=ios)
    if(ios/=0) then
        write(*,'(2a)') 'ERROR: failure openning file ',trim(filename)
        stop
    endif
    associate(nb=>symfunc_arr%symfunc(iconf)%linked_lists%maxbound_rad)
    associate(ng=>symfunc_arr%symfunc(iconf)%ng)
    associate(nat=>atoms_arr%atoms(iconf)%nat)
    if(parini%save_symfunc_force_ann) then
        nwa=3+nat*(3+ng)+ng*3*nb
    else
        nwa=3+nat*(3+ng)
    endif
    allocate(wa(1:nwa))
    wa(1)=real(nat,8)
    wa(2)=real(ng,8)
    wa(3)=real(nb,8)
    call update_ratp(atoms_arr%atoms(iconf))
    do iat=1,atoms_arr%atoms(iconf)%nat
        wa(3+iat*3-2)=atoms_arr%atoms(iconf)%ratp(1,iat)
        wa(3+iat*3-1)=atoms_arr%atoms(iconf)%ratp(2,iat)
        wa(3+iat*3-0)=atoms_arr%atoms(iconf)%ratp(3,iat)
    enddo
    n=3+3*atoms_arr%atoms(iconf)%nat
    do iat=1,atoms_arr%atoms(iconf)%nat
        do ig=1,symfunc_arr%symfunc(iconf)%ng
            n=n+1
            wa(n)=symfunc_arr%symfunc(iconf)%y(ig,iat)
        enddo
    enddo
    if(parini%save_symfunc_force_ann) then
        do ib=1,nb
            do i=1,3
                do ig=1,symfunc_arr%symfunc(iconf)%ng
                    n=n+1
                    wa(n)=symfunc_arr%symfunc(iconf)%y0d(ig,i,ib)
                enddo
            enddo
        enddo
    !else
    !    call f_free(symfunc_arr%symfunc(iconf)%y0d)
    endif
    end associate
    end associate
    end associate
    write(311,iostat=ios) wa
    if(ios/=0) then
        write(*,'(2a)') 'ERROR: failure writing to file ',trim(filename)
        stop
    endif
    !write(311,'(2i6)') atoms_arr%atoms(iconf)%nat,symfunc_arr%symfunc(iconf)%ng
    !do iat=1,atoms_arr%atoms(iconf)%nat
    !    do ig=1,symfunc_arr%symfunc(iconf)%ng
    !        write(311,'(es25.16)') symfunc_arr%symfunc(iconf)%y(ig,iat)
    !    enddo
    !enddo
    close(311)
    deallocate(wa)
end subroutine write_symfunc
!*****************************************************************************************
subroutine read_symfunc(parini,iconf,ann_arr,atoms_arr,strmess,symfunc_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc_arr
    use mod_atoms, only: typ_atoms_arr, update_ratp
    use mod_linked_lists, only: typ_pia_arr
    use mod_processors, only: iproc, nproc
    use mod_linkedlists, only: typ_linkedlists
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iconf
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms_arr), intent(inout):: atoms_arr
    character(*), intent(in):: strmess
    type(typ_symfunc_arr), intent(inout):: symfunc_arr
    !local variables
    real(8), allocatable:: wa(:)
    integer:: nwa
    character(30):: filename
    integer:: i, ig, iat, ios, nat_t, ng_t, nb_t, n, ib
    type(typ_pia_arr):: pia_arr_tmp
    real(8):: ttx, tty, ttz
    real(8):: eps=epsilon(1.d0)
    character(100):: smsg
    type(typ_linkedlists):: linkedlists
    !Symmetry functions which are previously calculated and written by
    !some other run is going to be read from files
    if(trim(strmess)=='bounds_train') then
        write(filename,'(a25,i5.5)') '../symfunc/train.symfunc.',iconf
    elseif(trim(strmess)=='bounds_valid') then
        write(filename,'(a25,i5.5)') '../symfunc/valid.symfunc.',iconf
    else
        stop 'ERROR: invalid content in strmess in gset_bounds '
    endif
    !open(unit=311,file=trim(filename),status='old',iostat=ios)
    open(unit=311,file=trim(filename),status='old',form='unformatted', &
        access='stream',iostat=ios)
    if(ios/=0) then
        write(*,'(2a)') 'ERROR: failure openning file ',trim(filename)
        stop
    endif
    !read(311,*) nat_t,ng_t
    associate(nb=>symfunc_arr%symfunc(iconf)%linked_lists%maxbound_rad)
    associate(ng=>symfunc_arr%symfunc(iconf)%ng)
    associate(nat=>atoms_arr%atoms(iconf)%nat)
    symfunc_arr%symfunc(iconf)%linked_lists%rcut=ann_arr%rcut
    symfunc_arr%symfunc(iconf)%linked_lists%triplex=.true.
    call linkedlists%calc_linkedlists(atoms_arr%atoms(iconf),.true.,symfunc_arr%symfunc(iconf)%linked_lists,pia_arr_tmp,&
        parini%mpi_env,parini%iverbose,parini%bondbased_ann)
    deallocate(pia_arr_tmp%pia)
    symfunc_arr%symfunc(iconf)%y=f_malloc0((/1.to.ng,1.to.nat/),id='symfunc%y')
    if(parini%save_symfunc_force_ann) then
        symfunc_arr%symfunc(iconf)%y0d=f_malloc0((/1.to.ng,1.to.3,1.to.nb/),id='symfunc%y0d')
        nwa=3+nat*(3+ng)+ng*3*nb
    else
        nwa=3+nat*(3+ng)
    endif
    wa=f_malloc([1.to.nwa],id='wa')
    read(311,iostat=ios) wa
    if(ios/=0) then
        write(*,'(2a)') 'ERROR: failure reading from file ',trim(filename)
        stop
    endif
    nat_t=nint(wa(1))
    ng_t=nint(wa(2))
    nb_t=nint(wa(3))
    if(parini%save_symfunc_behnam) then
    eps=eps*1.d4
    if(nat_t/=nat .or. ng_t/=ng) then
        write(*,'(a,7i6)') 'ERROR: inconsistent nat or ng ',iconf,nat_t,nat,ng_t,ng,nb_t,nb
        stop
    endif
    else
    if(nat_t/=nat .or. ng_t/=ng .or. nb_t/=nb) then
        write(*,'(a,7i6)') 'ERROR: inconsistent nat or ng ',iconf,nat_t,nat,ng_t,ng,nb_t,nb
        stop
    endif
    endif
    call update_ratp(atoms_arr%atoms(iconf))
    do iat=1,nat
        ttx=abs(wa(3+iat*3-2)-atoms_arr%atoms(iconf)%ratp(1,iat))
        tty=abs(wa(3+iat*3-1)-atoms_arr%atoms(iconf)%ratp(2,iat))
        ttz=abs(wa(3+iat*3-0)-atoms_arr%atoms(iconf)%ratp(3,iat))
        if(ttx>eps .or. tty>eps .or. ttz>eps) then
            smsg='ERROR: inconsistency of configuration in symmetry functions file. '
            !write(*,*) wa(3+iat*3-2),wa(3+iat*3-1),wa(3+iat*3-0)
            !write(*,*) atoms_arr%atoms(iconf)%rat(1,iat),atoms_arr%atoms(iconf)%rat(2,iat),atoms_arr%atoms(iconf)%rat(3,iat)
            !write(*,*) 'IAT',iat
            write(*,'(a,3es14.5,i7,a,i5)') trim(smsg),ttx,tty,ttz,iconf,trim(atoms_arr%fn(iconf)),atoms_arr%lconf(iconf)
            stop
        endif
    enddo
    n=3+3*atoms_arr%atoms(iconf)%nat
    do iat=1,nat
        do ig=1,ng
            n=n+1
            symfunc_arr%symfunc(iconf)%y(ig,iat)=wa(n)
        enddo
    enddo
    if(parini%save_symfunc_force_ann) then
        do ib=1,nb
            do i=1,3
                do ig=1,ng
                    n=n+1
                    symfunc_arr%symfunc(iconf)%y0d(ig,i,ib)=wa(n)
                enddo
            enddo
        enddo
    endif
    end associate
    end associate
    end associate
    !do iat=1,atoms_arr%atoms(iconf)%nat
    !    do ig=1,symfunc_arr%symfunc(iconf)%ng
    !        read(311,*) symfunc_arr%symfunc(iconf)%y(ig,iat)
    !    enddo
    !enddo
    close(311)
    call f_free(wa)
    write(*,*) "Reading symmetry functions done."
end subroutine read_symfunc
!*****************************************************************************************
subroutine save_gbounds(parini,ann_arr,atoms_arr,strmess,symfunc_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc_arr
    use mod_atoms, only: typ_atoms_arr
    use mod_linked_lists, only: typ_pia_arr
    use mod_processors, only: iproc, nproc
    use dynamic_memory
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms_arr), intent(inout):: atoms_arr
    character(*), intent(in):: strmess
    type(typ_symfunc_arr), intent(inout):: symfunc_arr
    !local variables
    integer:: iconf, ib, ig, i, iat, jat, i0
    real(8), allocatable:: gminarr(:,:), gmaxarr(:,:) !, poll_period
    integer, allocatable:: iatmin(:,:), iatmax(:,:), iconfmin(:,:), iconfmax(:,:)
    integer, allocatable:: ibmin(:), ibmax(:)
    integer:: ngmax
    ngmax=0
    do i=1,ann_arr%nann
        ngmax=max(ngmax,ann_arr%ann(i)%nn(0))
    enddo
    allocate(ibmin(ngmax),ibmax(ngmax))
    allocate(gminarr(1:ngmax,1:parini%ntypat))
    gminarr=huge(1.d20)
    allocate(gmaxarr(1:ngmax,1:parini%ntypat))
    gmaxarr=-huge(1.d20)
    allocate(iatmin(1:ngmax,1:parini%ntypat))
    iatmin=0.d0
    allocate(iatmax(1:ngmax,1:parini%ntypat))
    iatmax=0.d0
    allocate(iconfmin(1:ngmax,1:parini%ntypat))
    iconfmin=0.d0
    allocate(iconfmax(1:ngmax,1:parini%ntypat))
    iconfmax=0.d0
    ibmin(1:ngmax)=0 ; ibmax(1:ngmax)=0
    do iconf=1,atoms_arr%nconf
        !if(mod(iconf-1,nproc)==iproc) cycle
        !write(41,'(i6,i3)',advance='no') mod(iconf-1,nproc),iproc
        do iat=1,atoms_arr%atoms(iconf)%nat
            i=atoms_arr%atoms(iconf)%itypat(iat)
            do ig=1,symfunc_arr%symfunc(iconf)%ng
                if(symfunc_arr%symfunc(iconf)%y(ig,iat)<gminarr(ig,i)) then
                    gminarr(ig,i)=symfunc_arr%symfunc(iconf)%y(ig,iat)
                    iatmin(ig,i)=iat
                    iconfmin(ig,i)=iconf
                endif
                if(symfunc_arr%symfunc(iconf)%y(ig,iat)>gmaxarr(ig,i)) then
                    gmaxarr(ig,i)=symfunc_arr%symfunc(iconf)%y(ig,iat)
                    iatmax(ig,i)=iat
                    iconfmax(ig,i)=iconf
                endif
            enddo
        enddo
    enddo
    call yaml_mapping_open('symfunc bounds')
    if(iproc==0) then
    do i=1,ann_arr%nann
        do ig=1,ann_arr%ann(1)%nn(0) !HERE
            call yaml_sequence(advance='no')
            call yaml_mapping_open(trim(strmess),flow=.true.)
            !write(*,'(2(i7,2i4,es20.10),1x,a)') &
            !    iconfmin(ig,i),atoms_arr%atoms(iconfmin(ig,i))%nat,ibmin(ig),gminarr(ig,i), &
            !    iconfmax(ig,i),atoms_arr%atoms(iconfmax(ig,i))%nat,ibmax(ig),gmaxarr(ig,i),trim(strmess)
            call yaml_map('iconfmin',iconfmin(ig,i),fmt='(i7)')
            call yaml_map('nat',atoms_arr%atoms(iconfmin(ig,i))%nat,fmt='(i4)')
            call yaml_map('ibmin',ibmin(ig),fmt='(i4)')
            call yaml_map('gminarr',gminarr(ig,i),fmt='(es20.10)')
            call yaml_map('iconfmax',iconfmax(ig,i),fmt='(i7)')
            call yaml_map('nat',atoms_arr%atoms(iconfmax(ig,i))%nat,fmt='(i4)')
            call yaml_map('ibmax',ibmax(ig),fmt='(i4)')
            call yaml_map('gmaxarr',gmaxarr(ig,i),fmt='(es20.10)')
            call yaml_map('strmess',trim(strmess))
            !write(*,'(2(a50,i6,1x))') trim(atoms_arr%fn(iconfmin(ig,i))),atoms_arr%lconf(iconfmin(ig,i)), &
            !    trim(atoms_arr%fn(iconfmax(ig,i))),atoms_arr%lconf(iconfmax(ig,i))
            call yaml_map('fn_min',trim(atoms_arr%fn(iconfmin(ig,i))))
            call yaml_map('lconf_min',atoms_arr%lconf(iconfmin(ig,i)),fmt='(i6)')
            call yaml_map('fn_max',trim(atoms_arr%fn(iconfmax(ig,i))))
            call yaml_map('lconf_max',atoms_arr%lconf(iconfmax(ig,i)),fmt='(i6)')
            call yaml_mapping_close()
        enddo
    enddo
    endif
    call yaml_mapping_close()
    do i=1,ann_arr%nann
    do i0=1,ann_arr%ann(i)%nn(0)
        !if(abs(gminarr(i0))<epsilon(1.d0)
        if(gminarr(i0,i)==0.d0) then
            ann_arr%ann(i)%gbounds(1,i0)=-epsilon(1.d0)
        else
            ann_arr%ann(i)%gbounds(1,i0)=gminarr(i0,i)
        endif
        if(gmaxarr(i0,i)==0.d0) then
            ann_arr%ann(i)%gbounds(2,i0)=epsilon(1.d0)
        else
            ann_arr%ann(i)%gbounds(2,i0)=gmaxarr(i0,i)
        endif
        ann_arr%ann(i)%two_over_gdiff(i0)=2.d0/(ann_arr%ann(i)%gbounds(2,i0)-ann_arr%ann(i)%gbounds(1,i0))
    enddo
    enddo
    !if(trim(parini%symfunc)=='do_not_save') then
    !    do iconf=1,atoms_arr%nconf
    !        call f_free(symfunc_arr%symfunc(iconf)%y)
    !        deallocate(symfunc_arr%symfunc(iconf)%linked_lists%bound_rad)
    !    enddo
    !endif
    deallocate(gminarr)
    deallocate(gmaxarr)
    deallocate(iatmin)
    deallocate(iatmax)
    deallocate(iconfmin)
    deallocate(iconfmax)
    deallocate(ibmin,ibmax)
end subroutine save_gbounds
!*****************************************************************************************
subroutine apply_gbounds_atom(parini,ann_arr,atoms_arr,symfunc_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc_arr
    use mod_atoms, only: typ_atoms_arr
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(in):: ann_arr
    type(typ_atoms_arr), intent(inout):: atoms_arr
    type(typ_symfunc_arr), intent(inout):: symfunc_arr
    !local variables
    integer:: iconf, iat, ig, ib, i, i0, isat
    real(8):: tt
    do iconf=1,atoms_arr%nconf
        do iat=1,atoms_arr%atoms(iconf)%nat
            i=atoms_arr%atoms(iconf)%itypat(iat)
            do ig=1,symfunc_arr%symfunc(iconf)%ng
                tt=symfunc_arr%symfunc(iconf)%y(ig,iat)
                tt=(tt-ann_arr%ann(i)%gbounds(1,ig))*ann_arr%ann(i)%two_over_gdiff(ig)-1.d0
                symfunc_arr%symfunc(iconf)%y(ig,iat)=tt
            enddo
        enddo
        if(parini%save_symfunc_force_ann) then
            do ib=1,symfunc_arr%symfunc(iconf)%linked_lists%maxbound_rad
                iat=symfunc_arr%symfunc(iconf)%linked_lists%bound_rad(1,ib)
                isat=atoms_arr%atoms(iconf)%itypat(iat)
                do i0=1,ann_arr%ann(isat)%nn(0)
                    !normalization of y0d
                    tt=ann_arr%ann(isat)%two_over_gdiff(i0)
                    symfunc_arr%symfunc(iconf)%y0d(i0,1,ib)=symfunc_arr%symfunc(iconf)%y0d(i0,1,ib)*tt
                    symfunc_arr%symfunc(iconf)%y0d(i0,2,ib)=symfunc_arr%symfunc(iconf)%y0d(i0,2,ib)*tt
                    symfunc_arr%symfunc(iconf)%y0d(i0,3,ib)=symfunc_arr%symfunc(iconf)%y0d(i0,3,ib)*tt
                    !normalization of y0dr
                    !symfunc%y0dr(i0,1:9,ib)=symfunc%y0dr(i0,1:9,ib)*ann_arr%ann(isat)%two_over_gdiff(i0)
                enddo
            enddo
        endif
    enddo
end subroutine apply_gbounds_atom
!*****************************************************************************************
end module mod_refdata
!*****************************************************************************************
