!*****************************************************************************************
subroutine test_symfunc_atom_behler()
    use iso_fortran_env, only: error_unit, output_unit
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_deallocate
    use mod_symfunc, only: typ_symfunc
    use mod_yaml_conf, only: read_yaml_conf
    use mod_ann_io_yaml, only: read_input_ann_yaml
    use mod_colors, only: green_passed, red_failed
    use mod_flm_futile
    implicit none
    !local variables
    real(8):: errmax
    type(typ_parini):: parini
    type(typ_symfunc):: symfunc
    type(typ_ann_arr):: ann_arr
    type(typ_atoms_arr):: atoms_arr
    type(mpi_environment):: mpi_env
    integer:: iconf, iat, iann, ng, nat, nb, ig, ii, jj, ios, ib
    real(8), allocatable:: y(:,:)
    real(8), allocatable:: y0d(:,:,:)
    real(8), allocatable:: y0dr(:,:,:)
    parini%iverbose=3
    parini%symfunc_type_ann='behler'
    parini%types_main='Ca F'
    call set_atomc_types_info(parini)
    mpi_env%nproc=1
    mpi_env%iproc=0
    call read_yaml_conf(parini,'tests-fpm/data/posinp.yaml',10000,atoms_arr)
    !write(*,*) atoms_arr%atoms(1)%nat
    !write(*,*) trim(atoms_arr%atoms(1)%boundcond)
    associate(atoms=>atoms_arr%atoms(1))
    ann_arr%approach=trim(parini%approach_ann)
    call ann_arr%set_number_of_ann(parini%ntypat)
    if(ann_arr%nann==0) stop 'ERROR: number of type of atoms zero in test_symfunc_atom_behler'
    allocate(ann_arr%ann(ann_arr%nann))
    do iat=1,atoms%nat
        do iann=1,ann_arr%nann
            if(trim(adjustl(atoms%sat(iat)))==trim(parini%stypat(iann))) then
                atoms%itypat(iat)=parini%ltypat(iann)
                exit
            endif
        enddo
    enddo
    call read_input_ann_yaml(parini,mpi_env%iproc,ann_arr,'tests-fpm/data')
    call symfunc%init_symfunc(mpi_env,parini%iverbose,.false.,parini%symfunc_type_ann)

    !ann_arr%ann(1)%g2eta(1)=ann_arr%ann(1)%g2eta(1)+1.d-3
    !ann_arr%ann(1)%g2eta(2)=ann_arr%ann(1)%g2eta(2)+1.d-3
    !ann_arr%ann(2)%g2eta(1)=ann_arr%ann(2)%g2eta(1)+1.d-3
    !ann_arr%ann(2)%g2eta(2)=ann_arr%ann(2)%g2eta(2)+1.d-3

    call symfunc%get_symfunc(ann_arr,atoms,.false.)
    open(unit=761,file='tests-fpm/data/symfunc_atom_behler_free.ref',status='old',iostat=ios)
    read(761,*) ng,nat,nb
    !write(*,*) ng,nat,nb
    !write(*,*) ann_arr%ann(1)%nn(0),atoms%nat,symfunc%linked_lists%maxbound_rad
    if(ng/=ann_arr%ann(1)%nn(0) .or. nat/=atoms%nat .or. nb/=symfunc%linked_lists%maxbound_rad) then
        write(error_unit,'(2a)') red_failed,' in test_symfunc_atom_behler: '
        write(error_unit,'(a,3i6)') 'ng,nat,nb from ref. file: ',ng,nat,nb
        write(error_unit,'(a,3i6)') '          current values: ', &
            ann_arr%ann(1)%nn(0),atoms%nat,symfunc%linked_lists%maxbound_rad
        call exit(1)
    endif
    allocate(y(ng,nat))
    allocate(y0d(ng,3,nb))
    allocate(y0dr(ng,9,nb))
    do iat=1,atoms%nat
        do ig=1,ng
            read(761,*) ii,jj,y(ig,iat)
        enddo
    enddo
    do ib=1,nb
        do ig=1,ng
            read(761,*) ii,jj,y0d(ig,1,ib),y0d(ig,2,ib),y0d(ig,3,ib)
        enddo
    enddo
    do ib=1,nb
        do ig=1,ng
            read(761,*) ii,jj, &
                y0dr(ig,1,ib),y0dr(ig,2,ib),y0dr(ig,3,ib), &
                y0dr(ig,4,ib),y0dr(ig,5,ib),y0dr(ig,6,ib), &
                y0dr(ig,7,ib),y0dr(ig,8,ib),y0dr(ig,9,ib)
        enddo
    enddo
    close(761)
    errmax=0.d0
    do iat=1,atoms%nat
        do ig=1,ng
            errmax=max(errmax,abs(symfunc%y(ig,iat)-y(ig,iat)))
        enddo
    enddo
    do ib=1,nb
        do ig=1,ng
            errmax=max(errmax,abs(symfunc%y0d(ig,1,ib)-y0d(ig,1,ib)))
            errmax=max(errmax,abs(symfunc%y0d(ig,2,ib)-y0d(ig,2,ib)))
            errmax=max(errmax,abs(symfunc%y0d(ig,3,ib)-y0d(ig,3,ib)))
            errmax=max(errmax,abs(symfunc%y0dr(ig,1,ib)-y0dr(ig,1,ib)))
            errmax=max(errmax,abs(symfunc%y0dr(ig,2,ib)-y0dr(ig,2,ib)))
            errmax=max(errmax,abs(symfunc%y0dr(ig,3,ib)-y0dr(ig,3,ib)))
            errmax=max(errmax,abs(symfunc%y0dr(ig,4,ib)-y0dr(ig,4,ib)))
            errmax=max(errmax,abs(symfunc%y0dr(ig,5,ib)-y0dr(ig,5,ib)))
            errmax=max(errmax,abs(symfunc%y0dr(ig,6,ib)-y0dr(ig,6,ib)))
            errmax=max(errmax,abs(symfunc%y0dr(ig,7,ib)-y0dr(ig,7,ib)))
            errmax=max(errmax,abs(symfunc%y0dr(ig,8,ib)-y0dr(ig,8,ib)))
            errmax=max(errmax,abs(symfunc%y0dr(ig,9,ib)-y0dr(ig,9,ib)))
        enddo
    enddo
    if(errmax<1.d-13) then
        write(output_unit,'(2a)') green_passed,' in test_symfunc_atom_behler: errmax'
    else
        write(error_unit,'(2a,es14.5)') red_failed,' in test_symfunc_atom_behler: errmax=  ',errmax
        call exit(1)
    end if
    call symfunc%fini_symfunc()

    deallocate(y,y0d,y0dr)
    end associate
    
    do iconf=1,atoms_arr%nconf
        call atom_deallocate(atoms_arr%atoms(iconf))
    enddo
    deallocate(atoms_arr%atoms)
end subroutine test_symfunc_atom_behler
!*****************************************************************************************
