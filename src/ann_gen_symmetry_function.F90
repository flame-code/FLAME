!*****************************************************************************************
subroutine ann_gen_symmetry_function(parini)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, ann_arr_deallocate
    use mod_symfunc, only: typ_symfunc, typ_symfunc_arr
    use mod_atoms, only: typ_atoms_arr
    use mod_processors, only: iproc, mpi_comm_abz
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    type(typ_ann_arr) :: ann_arr
    type(typ_atoms_arr):: atoms_gen
    type(typ_symfunc_arr):: symfunc_gen
    type(typ_symfunc):: symfunc
    character(400):: fnarr(100000), fn_tmp, filename
    integer:: iat,jat,i,j,ig,jg,iconf,jconf, ng, nat,a,l,ll
    real(8) :: distance
    real(8), allocatable:: diff(:),tt1(:),tt2(:)
    integer:: nconftot, ios, k
    character (50)::fname
    logical:: file_exists
    
    !write(*,*) trim(parini%stypat_ann)
    !call count_words(parini%stypat_ann,ann_arr%nann)
    !read(parini%stypat_ann,*) ann_arr%stypat(1:ann_arr%nann)
    ann_arr%nann=parini%ntypat
    
    !do i=1,ann_arr%nann
    !    ann_arr%ltypat(i)=i
    !    write(*,*) i,ann_arr%stypat(i)
    !enddo
    
    if(ann_arr%nann==0) stop 'ERROR: number of type of atoms zero in gen_symmetry_function'
    allocate(ann_arr%ann(ann_arr%nann))
    ann_arr%approach=trim(parini%approach_ann)
    
    fname = trim(parini%stypat(1))//'.ann.input.yaml'
    inquire(file=trim(fname),exist=ann_arr%exists_yaml_file)
    if( ann_arr%exists_yaml_file) then
        call read_input_ann_yaml(parini,iproc,ann_arr)
    else
        call read_input_ann(parini,iproc,ann_arr)
    endif
    
    inquire(file="list_posinp_gen.yaml",exist=file_exists)
    if(file_exists) then
        call read_data_yaml(parini,'list_posinp_gen.yaml',atoms_gen)
    else
        call read_data_old(parini,'list_posinp_gen',atoms_gen)
    endif
    !---------------------------------------------------------- 
    open(unit=1,file='list_posinp_gen',status='old',iostat=ios)
    if(ios/=0) then;write(*,'(a)') 'ERROR: failure openning list_posinp_gen__';stop;endif
    nconftot=0
    do
        read(1,'(a)',iostat=k) fn_tmp
        if(k<0) exit
        !fn_tmp=adjustl(trim(filename))
        if(fn_tmp(1:1)=='#') cycle
        nconftot=nconftot+1
        fnarr(nconftot)=fn_tmp(9:)
        !write(*,'(a)') trim(fn_tmp)
        !write(*,'(a)') trim(fnarr(nconftot))
        !stop
    enddo
    close(1)
    !----------------------------------------------------------   
    if(iproc==0) then
        write(*,'(a36,i8)') 'number of generating data points:   ',atoms_gen%nconf
    endif
    do iconf=1,atoms_gen%nconf
        do iat=1,atoms_gen%atoms(iconf)%nat
            do i=1,ann_arr%nann 
                if(trim(atoms_gen%atoms(iconf)%sat(iat))==trim(parini%stypat(i))) then
                atoms_gen%atoms(iconf)%itypat(iat)=parini%ltypat(i)
                exit
                endif
            enddo
        enddo
    enddo
    !----------------------------------------------------------   
    do iconf=1,atoms_gen%nconf
        filename=trim(fnarr(iconf))//'.fp'
        write(*,'(a)') trim(filename)
        open(unit=1,file=trim(filename),status='replace',iostat=ios)
        if(ios/=0) then;write(*,'(2a)') 'ERROR: failure openning ',trim(filename);stop;endif
        call symmetry_functions(parini,ann_arr,atoms_gen%atoms(iconf),symfunc,.false.)
        do iat=1,atoms_gen%atoms(iconf)%nat
            do ig=1,ann_arr%ann(1)%nn(0)
                write(1,'(2i4,es24.15)') ig,iat,symfunc%y(ig,iat)
                !write(999,'(2i4,es24.15,i4)') ig,iat, ann_arr%yall(ig,iat),atoms_gen%nconf
            enddo
        enddo
        close(1)
    enddo
    call ann_arr_deallocate(ann_arr)
end subroutine ann_gen_symmetry_function 
!*****************************************************************************************
