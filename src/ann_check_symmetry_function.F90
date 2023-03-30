!*****************************************************************************************
subroutine ann_check_symmetry_function(parini,path)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc
    use mod_ann_io_yaml, only: read_input_ann_yaml, read_data_yaml
    use mod_yaml_conf, only: write_yaml_conf
    use mod_atoms, only: typ_atoms_arr, atom_deallocate_old, atom_copy_old
    use mod_atoms, only: typ_file_info
    use mod_flm_futile
    implicit none
    type(typ_parini), intent(in):: parini
    character(len=*), intent(in):: path
    !local variables
    type(typ_ann_arr) :: ann_arr
    type(typ_atoms_arr):: atoms_check
    type(typ_atoms_arr):: atoms_sel
    type(typ_symfunc):: symfunc
    type(mpi_environment):: mpi_env
    type(typ_file_info):: file_info
    logical:: file_exists, conf_new
    character(30):: fnout, fnout1
    character(50):: fname
    integer:: iat, jat, i, ig, iconf, jconf, i0, nconf_sel
    integer:: ios, ipair, npair, ipairs, ipaire, mpair
    integer:: ierr, icolor, ikey, jproc, mproc, mconf, iconfs, iconfe
    real(8):: tt, distance, distance2, de
    real(8):: time1, time2, time3
    real(8), allocatable:: gminarr(:)
    real(8), allocatable:: gmaxarr(:)
    real(8), allocatable:: diff(:), c(:,:)
    real(8), allocatable:: yall(:,:,:)
    real(8), allocatable:: distance_all(:)
    real(8), allocatable:: dist(:,:)
    integer, allocatable:: iconf_sel(:)
    integer, allocatable:: F(:)
    integer, allocatable:: ncounts(:)
    integer, allocatable:: idispls(:)
    integer, allocatable:: ind_pairs(:,:)
#if defined(MPI)
    include 'mpif.h'
#endif
    call f_routine(id='ann_check_symmetry_function')
    !----------------------------------------------------------
    !write(*,*) trim(parini%stypat_ann)
    !call count_words(parini%stypat_ann,ann_arr%nann)
    !read(parini%stypat_ann,*) ann_arr%stypat(1:ann_arr%nann)
    if(parini%mpi_env%nproc>1) then
        icolor=parini%mpi_env%iproc+1
        ikey=1
        call MPI_COMM_SPLIT(MPI_COMM_WORLD,icolor,ikey,mpi_env%mpi_comm,ierr)
    endif
    mpi_env%nproc=1
    mpi_env%iproc=0
    call symfunc%init_symfunc(mpi_env,parini%iverbose,parini%bondbased_ann,parini%symfunc_type_ann)
    ann_arr%nann=parini%ntypat
    if(ann_arr%nann==0) stop 'ERROR: number of type of atoms zero in check_symmetry_function'
    allocate(ann_arr%ann(ann_arr%nann))
    ann_arr%approach=trim(parini%approach_ann)
    fname=trim(path)//'/'//trim(parini%stypat(1))//'.ann.input.yaml'
    inquire(file=trim(fname),exist=ann_arr%exists_yaml_file)
    if(ann_arr%exists_yaml_file) then
        call read_input_ann_yaml(parini,parini%mpi_env%iproc,ann_arr,path)
    else
        call read_input_ann(parini,parini%mpi_env%iproc,ann_arr)
    endif
    inquire(file="list_posinp_check.yaml",exist=file_exists)
    if(file_exists) then
        call read_data_yaml(parini,'list_posinp_check.yaml',atoms_check)
    else
        call read_data_old(parini,'list_posinp_check',atoms_check)
    endif
    !----------------------------------------------------------
    if(parini%mpi_env%iproc==0 .and. parini%iverbose>1) then
        write(*,'(a34,i8)') 'number of checking data points:   ',atoms_check%nconf
    endif
    do iconf=1,atoms_check%nconf
        do iat=1,atoms_check%atoms(iconf)%nat
            do i=1,ann_arr%nann
                if(trim(atoms_check%atoms(iconf)%sat(iat))==trim(parini%stypat(i))) then
                    atoms_check%atoms(iconf)%itypat(iat)=parini%ltypat(i)
                    exit
                endif
            enddo
        enddo
    enddo
    !----------------------------------------------------------
    do i=2,ann_arr%nann
        if(ann_arr%ann(i)%nn(0)/=ann_arr%ann(1)%nn(0)) then
            write(*,*) 'ERROR: current implementation for identical length of descriptor for all elements'
            stop
        endif
    enddo
    do iconf=2,atoms_check%nconf
        if(atoms_check%atoms(iconf)%nat/=atoms_check%atoms(1)%nat) then
            write(*,*) 'ERROR: this subtask can be used for systems having the same size'
            stop
        endif
    enddo
    associate(ng=>ann_arr%ann(1)%nn(0))
    associate(nat=>atoms_check%atoms(1)%nat)
    allocate(diff(ng),c(nat,nat),F(nat))
    allocate(yall(ng,nat,atoms_check%nconf),source=0.d0)
    allocate(gminarr(ng),gmaxarr(ng)) 
    gminarr(1:ng)=huge(1.d0) ; gmaxarr(1:ng)=-huge(1.d0)
    !----------------- Compute symmetry functions ----------
    call cpu_time(time1)
    if(parini%mpi_env%nproc>1) then
        allocate(ncounts(0:parini%mpi_env%nproc-1))
        allocate(idispls(0:parini%mpi_env%nproc-1))
        do jproc=0,parini%mpi_env%nproc-1
            mconf=atoms_check%nconf/parini%mpi_env%nproc
            iconfs=jproc*mconf+1
            mproc=mod(atoms_check%nconf,parini%mpi_env%nproc)
            iconfs=iconfs+max(0,jproc-parini%mpi_env%nproc+mproc)
            if(jproc>parini%mpi_env%nproc-mproc-1) mconf=mconf+1
            iconfe=iconfs+mconf-1
            ncounts(jproc)=mconf*ng*nat
            if(jproc==0) then
                idispls(0)=0
            else
                idispls(jproc)=idispls(jproc-1)+ncounts(jproc-1)
            endif
        enddo
        mconf=atoms_check%nconf/parini%mpi_env%nproc
        iconfs=parini%mpi_env%iproc*mconf+1
        mproc=mod(atoms_check%nconf,parini%mpi_env%nproc)
        iconfs=iconfs+max(0,parini%mpi_env%iproc-parini%mpi_env%nproc+mproc)
        if(parini%mpi_env%iproc>parini%mpi_env%nproc-mproc-1) mconf=mconf+1
        iconfe=iconfs+mconf-1
    else
        iconfs=1
        iconfe=atoms_check%nconf
    endif
    do iconf=iconfs,iconfe
        call symfunc%get_symfunc(ann_arr,atoms_check%atoms(iconf),.false.)
        do iat=1,atoms_check%atoms(iconf)%nat
            do ig=1,ng
                yall(ig,iat,iconf)=symfunc%y(ig,iat)
            enddo
        enddo
        if(parini%symfunc_type_ann=='behler') then
            deallocate(symfunc%linked_lists%prime_bound)
            deallocate(symfunc%linked_lists%bound_rad)
            deallocate(symfunc%linked_lists%bound_ang)
        endif
        call f_free(symfunc%y)
        call f_free(symfunc%y0d)
        call f_free(symfunc%y0dr)
    enddo
    if(parini%mpi_env%nproc>1) then
    call MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,yall,ncounts,idispls,MPI_DOUBLE_PRECISION,parini%mpi_env%mpi_comm,ierr)
    endif
    !----------------- Normalize if required ---------------
    if(parini%normalization_ann) then
        do iconf=1,atoms_check%nconf
            do iat=1,atoms_check%atoms(iconf)%nat
                do ig=1,ng
                    if(yall(ig,iat,iconf)<gminarr(ig)) then
                        gminarr(ig)=yall(ig,iat,iconf)
                    endif
                    if(yall(ig,iat,iconf)>gmaxarr(ig)) then
                        gmaxarr(ig)=yall(ig,iat,iconf)
                    endif
                enddo
            enddo
        enddo
        do i=1,ann_arr%nann
            do i0=1,ann_arr%ann(i)%nn(0)
                ann_arr%ann(i)%gbounds(1,i0)=gminarr(i0)
                ann_arr%ann(i)%gbounds(2,i0)=gmaxarr(i0)
                ann_arr%ann(i)%two_over_gdiff(i0)=2.d0/(ann_arr%ann(i)%gbounds(2,i0)-ann_arr%ann(i)%gbounds(1,i0))
            enddo
        enddo
        do iconf=1,atoms_check%nconf
            do iat=1,atoms_check%atoms(iconf)%nat
                i=atoms_check%atoms(iconf)%itypat(iat)
                do ig=1,ng
                    tt=yall(ig,iat,iconf)
                    tt=(tt-ann_arr%ann(i)%gbounds(1,ig))*ann_arr%ann(i)%two_over_gdiff(ig)-1.d0
                    yall(ig,iat,iconf)=tt
                    !write(88,'(2i4,es24.15)') ig,iat,symfunc_check%symfunc(1)%y(ig,iat)
                enddo
            enddo
        enddo
    endif
    call cpu_time(time2)
    !----------------------------------------------------------
    npair=atoms_check%nconf*(atoms_check%nconf-1)/2
    allocate(ind_pairs(2,npair))
    allocate(distance_all(npair))
    ipair=0
    do iconf=1,atoms_check%nconf
        do jconf=iconf+1,atoms_check%nconf
            ipair=ipair+1
            ind_pairs(1,ipair)=iconf
            ind_pairs(2,ipair)=jconf
        enddo !over jconf
    enddo !over iconf
    if(ipair/=npair) stop 'ERROR: ipair/=npair'
    if(parini%mpi_env%nproc>1) then
        do jproc=0,parini%mpi_env%nproc-1
            mpair=npair/parini%mpi_env%nproc
            ipairs=jproc*mpair+1
            mproc=mod(npair,parini%mpi_env%nproc)
            ipairs=ipairs+max(0,jproc-parini%mpi_env%nproc+mproc)
            if(jproc>parini%mpi_env%nproc-mproc-1) mpair=mpair+1
            ipaire=ipairs+mpair-1
            ncounts(jproc)=mpair
            if(jproc==0) then
                idispls(0)=0
            else
                idispls(jproc)=idispls(jproc-1)+ncounts(jproc-1)
            endif
        enddo
        mpair=npair/parini%mpi_env%nproc
        ipairs=parini%mpi_env%iproc*mpair+1
        mproc=mod(npair,parini%mpi_env%nproc)
        ipairs=ipairs+max(0,parini%mpi_env%iproc-parini%mpi_env%nproc+mproc)
        if(parini%mpi_env%iproc>parini%mpi_env%nproc-mproc-1) mpair=mpair+1
        ipaire=ipairs+mpair-1
    else
        ipairs=1
        ipaire=npair
    endif
    do ipair=ipairs,ipaire
        iconf=ind_pairs(1,ipair)
        jconf=ind_pairs(2,ipair)
        do iat=1,atoms_check%atoms(iconf)%nat
            do jat=1,atoms_check%atoms(jconf)%nat
                do ig=1,ng
                    diff(ig)=yall(ig,iat,iconf)-yall(ig,jat,jconf)
                enddo
                c(iat,jat)=dot_product(diff,diff)
            enddo !over jat
        enddo !over iat
        call hung(atoms_check%atoms(iconf)%nat,c,F,distance2)
        distance_all(ipair)=sqrt(distance2)
    enddo
    if(parini%mpi_env%nproc>1) then
    call MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,distance_all,ncounts,idispls,MPI_DOUBLE_PRECISION,parini%mpi_env%mpi_comm,ierr)
    deallocate(ncounts)
    deallocate(idispls)
    endif
    if(.not. parini%pickdiffconfs) then
    if(parini%mpi_env%iproc==0) then
    fnout='incompatible'
    fnout1='distall'
    open(unit=111,file=fnout1,status='replace',iostat=ios)
    open(unit=11,file=fnout,status='replace',iostat=ios)
    if(ios/=0) then
        write(*,'(a)') 'ERROR: failure openning output file'
        stop
    endif
    do ipair=1,npair
        iconf=ind_pairs(1,ipair)
        jconf=ind_pairs(2,ipair)
        distance=distance_all(ipair)
        de=abs(atoms_check%atoms(iconf)%epot-atoms_check%atoms(jconf)%epot)
        write(111,'(2(a25,i6,1x),2es20.10)') adjustl(trim(atoms_check%fn(iconf))),atoms_check%lconf(iconf), &
            adjustl(trim(atoms_check%fn(jconf))),atoms_check%lconf(jconf),distance,de
        if(distance.le.parini%dtol_ann .and. de.ge.parini%etol_ann) then
            write(11,'(2(a25,i6,1x),2es20.10)') adjustl(trim(atoms_check%fn(iconf))),atoms_check%lconf(iconf), &
                adjustl(trim(atoms_check%fn(jconf))),atoms_check%lconf(jconf),distance,de
            cycle
        endif
    enddo
    close(11)
    close(111)
    deallocate(ind_pairs)
    deallocate(distance_all)
    endif !end of if mpi_env%iproc==0
    endif !end of if .not. parini%pickdiffconfs

    if(parini%pickdiffconfs) then
    allocate(dist(atoms_check%nconf,atoms_check%nconf),source=0.d0)
    do ipair=1,npair
        iconf=ind_pairs(1,ipair)
        jconf=ind_pairs(2,ipair)
        dist(iconf,jconf)=distance_all(ipair)
        dist(jconf,iconf)=distance_all(ipair)
    enddo

    allocate(iconf_sel(atoms_check%nconf),source=0)
    nconf_sel=0
    do iconf=1,atoms_check%nconf
        if(iconf==1) then
            nconf_sel=nconf_sel+1
            iconf_sel(nconf_sel)=iconf
            continue
        endif
        conf_new=.true.
        do jconf=1,nconf_sel
            if(dist(iconf,iconf_sel(jconf))<parini%dtol_pickdiffconfs) then
                conf_new=.false.
                exit
            endif
        enddo
        if(conf_new) then
            nconf_sel=nconf_sel+1
            iconf_sel(nconf_sel)=iconf
        endif
    enddo
    atoms_sel%nconf=nconf_sel
    allocate(atoms_sel%atoms(atoms_sel%nconf))
    do iconf=1,nconf_sel
        call atom_copy_old(atoms_check%atoms(iconf_sel(iconf)),atoms_sel%atoms(iconf),'atoms_check->atoms_sel')
    enddo
    file_info%filename_positions='posout.yaml'
    do iconf=1,atoms_sel%nconf
        if(iconf==1) then
            file_info%file_position='new'
        else
            file_info%file_position='append'
        endif
        call write_yaml_conf(file_info,atoms_sel%atoms(iconf),strkey='posout')
    enddo
    do iconf=1,nconf_sel
        call atom_deallocate_old(atoms_sel%atoms(iconf))
    enddo
    deallocate(atoms_sel%atoms)
    deallocate(dist)
    deallocate(iconf_sel)
    endif !end of if parini%pickdiffconfs

    call cpu_time(time3)
    if(parini%iverbose>2) then
        write(*,*) 'timing ',time2-time1,time3-time2
    endif
    do iconf=1,atoms_check%nconf
        call atom_deallocate_old(atoms_check%atoms(iconf))
    enddo
    deallocate(atoms_check%atoms)
    deallocate(atoms_check%fn)
    deallocate(atoms_check%lconf)
    deallocate(diff,c,F)
    deallocate(yall)
    deallocate(gminarr,gmaxarr)
    call symfunc%fini_symfunc()
    end associate
    end associate
    if(parini%mpi_env%nproc>1) then
        call MPI_COMM_FREE(mpi_env%mpi_comm,ierr)
    endif
    call f_release_routine()
end subroutine ann_check_symmetry_function
!*****************************************************************************************
