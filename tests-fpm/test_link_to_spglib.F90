!*****************************************************************************************
subroutine test_link_to_spglib()
    use iso_fortran_env, only: error_unit, output_unit
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_allocate_old, atom_deallocate_old
    use mod_atoms, only: atom_copy_old, set_rat, typ_file_info
    use mod_yaml_conf, only: write_yaml_conf
    use mod_colors, only: green_passed, red_failed
    use mod_const, only: bohr2ang
    use mod_gensymcrys, only: typ_gensymcrys
    use mod_flm_futile
    implicit none
    !local variables
    type(typ_parini):: parini
    integer:: nsym_tot, NCELLS, nconf, ntry, NAT_CELL_MAX, iconf, NKINDS, i, iat, ispg_err
    integer:: i_spacegroup(10)=[186,216,136,225,194,1,221,160,36,44], i_spacegroup_t(10)
    real(8):: target_vol_per_atom, tol
    !real(8):: time1, time2
    logical:: succeeded
    type(typ_gensymcrys):: gensymcrys
    integer, allocatable:: cryssys_all(:), brav_all(:), nsymp_all(:), nsym_all(:), ind_rsym_all(:)
    integer, allocatable:: NAT_CELL_ALL(:)
    integer, allocatable:: NAT_KINDS(:)
    integer, allocatable:: itypat(:)
    real(8), allocatable:: KINDSDIST_MIN(:,:)
    real(8), allocatable:: rsym_all(:,:,:)
    real(8), allocatable:: rxyz_all(:,:,:)
    real(8), allocatable:: cv_all(:,:,:)
    real(8), allocatable:: ratred(:,:)
    character(5), allocatable:: stypat(:)
    character(5), allocatable:: sat_all(:,:)
    type(typ_atoms):: atoms
    type(typ_file_info):: file_info
    !call cpu_time(time1)
    parini%iverbose=0
    parini%task='genconf'
    parini%task='gensymcrys'
    parini%types_main='Zn O'
    parini%datafilesdir='./datafiles'
    call set_atomc_types_info(parini)
    parini%mpi_env%nproc=1
    parini%mpi_env%iproc=0
    call gensymcrys%init_gensymcrys(parini%iverbose)
    nconf=10
    ntry=20
    NCELLS=1
    NAT_CELL_MAX=16
    NKINDS=parini%ntypat
    target_vol_per_atom=12.d0
    allocate(stypat(NKINDS))
    do i=1,NKINDS
        stypat(i)=parini%stypat(i)
    enddo
    allocate(NAT_KINDS(NKINDS))
    do i=1,NKINDS
        NAT_KINDS(i)=8 !eight Zn atoms and 8 O atoms
    enddo
    allocate(KINDSDIST_MIN(NKINDS,NKINDS))
    KINDSDIST_MIN(1,1)=1.22d0+1.22d0 ; KINDSDIST_MIN(1,2)=1.22d0+0.66d0
    KINDSDIST_MIN(2,1)=1.22d0+0.66d0 ; KINDSDIST_MIN(2,2)=0.66d0+0.66d0
    KINDSDIST_MIN=KINDSDIST_MIN*0.9d0
    call set_nsym_tot(nsym_tot)
    allocate(cryssys_all(230),brav_all(230),nsymp_all(230),nsym_all(230),ind_rsym_all(230))
    allocate(rsym_all(4,4,nsym_tot))
    allocate(NAT_CELL_ALL(nconf),cv_all(3,3,nconf))
    do iconf=1,nconf
    NAT_CELL_ALL(iconf)=sum(NAT_KINDS(1:NKINDS))
    enddo
    allocate(rxyz_all(3,NAT_CELL_MAX,nconf),sat_all(NAT_CELL_MAX,nconf))
    if(parini%mpi_env%iproc==0) then
    call read_coeffs_gensymcrys(parini,nsym_tot,rsym_all,cryssys_all,brav_all,nsymp_all,nsym_all,ind_rsym_all)
    endif
    do iconf=1,nconf
    call gensymcrys_single(gensymcrys,nsym_tot,rsym_all,cryssys_all,brav_all,nsymp_all, &
        nsym_all,ind_rsym_all,target_vol_per_atom,stypat,NCELLS,NAT_CELL_ALL(iconf),ntry, &
        i_spacegroup(iconf),NKINDS,NAT_KINDS,KINDSDIST_MIN,cv_all(1,1,iconf),rxyz_all(1,1,iconf),sat_all(1,iconf),succeeded)
    enddo
    !call cpu_time(time2)
    !write(*,*) 'time= ',time2-time1
    tol=1.d-3
    allocate(itypat(NAT_CELL_ALL(1)))
    itypat(1:8)=1
    itypat(9:16)=2
    allocate(ratred(3,NAT_CELL_ALL(1)))
    ispg_err=0
    do iconf=1,nconf
    call rxyz_cart2int_alborz(NAT_CELL_ALL(iconf),cv_all(1,1,iconf),rxyz_all(1,1,iconf),ratred)
    call get_spg(NAT_CELL_ALL(iconf),ratred,cv_all(1,1,iconf),itypat,tol,i_spacegroup_t(iconf))
    if(i_spacegroup_t(iconf)/=i_spacegroup(iconf)) then
        ispg_err=ispg_err+1
    endif
    enddo
    deallocate(ratred)

    !if(parini%mpi_env%iproc==0) then
    !file_info%filename_positions='posout.yaml'
    !do iconf=1,nconf !nconf_allproc
    !    call atom_allocate_old(atoms,NAT_CELL_ALL(iconf),0,0)
    !    rxyz_all(1:3,1:atoms%nat,iconf)=rxyz_all(1:3,1:atoms%nat,iconf)/bohr2ang
    !    atoms%cellvec(1:3,1:3)=cv_all(1:3,1:3,iconf)/bohr2ang
    !    call set_rat(atoms,rxyz_all(1,1,iconf),.true.)
    !    atoms%boundcond='bulk'
    !    atoms%units_length_io='angstrom'
    !    do iat=1,atoms%nat
    !        atoms%sat(iat)=sat_all(iat,iconf)
    !    enddo
    !    if(iconf==1) then
    !        file_info%file_position='new'
    !    else
    !        file_info%file_position='append'
    !    endif
    !    call write_yaml_conf(file_info,atoms=atoms,strkey='posout')
    !    call atom_deallocate_old(atoms)
    !enddo
    !endif

    if(ispg_err==0) then
        write(output_unit,'(2a)') green_passed,' in test_link_to_spglib: space group of all 10 configurations are OK!'
    else
        write(error_unit,'(2a,i6)') red_failed,' in test_link_to_spglib: ispg_err=  ',ispg_err
        call exit(1)
    end if
end subroutine test_link_to_spglib
!*****************************************************************************************
