!*****************************************************************************************
subroutine ann_check_symmetry_function(parini)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc, typ_symfunc_arr
    use mod_atoms, only: typ_atoms_arr
    use mod_processors, only: iproc, mpi_comm_abz
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    type(typ_ann_arr) :: ann_arr
    type(typ_atoms_arr):: atoms_check
    type(typ_symfunc_arr):: symfunc_check
    type(typ_symfunc):: symfunc
    character(400):: fnarr(100000), fn_tmp, filename
    character(30):: fnout,fnout1
    integer:: iat,jat,i,j,ig,jg,iconf,jconf, ng, nat,a,i0
    integer:: iatmin(140), iatmax(140), iconfmin(140), iconfmax(140)
    real(8) :: tt,distance,distance2,de
    real(8), allocatable:: diff(:),c(:,:)
    real(8):: gminarr(140), gmaxarr(140)
    integer,allocatable:: F(:)
    integer:: nconftot, ios, k
    character (50)::fname
    logical:: file_exists
    call f_routine(id='ann_check_symmetry_function')
    associate(etol=>parini%etol_ann,dtol=>parini%dtol_ann)
    !---------------------------------------------------------- 
    !write(*,*) trim(parini%stypat_ann)
    !call count_words(parini%stypat_ann,ann_arr%nann)
    !read(parini%stypat_ann,*) ann_arr%stypat(1:ann_arr%nann)
    ann_arr%nann=parini%ntypat
    !do i=1,ann_arr%nann
    !    ann_arr%ltypat(i)=i
    !    write(*,*) i,ann_arr%stypat(i)
    !enddo
    if(ann_arr%nann==0) stop 'ERROR: number of type of atoms zero in check_symmetry_function'
    allocate(ann_arr%ann(ann_arr%nann))
    ann_arr%approach=trim(parini%approach_ann)
    fname = trim(parini%stypat(1))//'.ann.input.yaml'
    inquire(file=trim(fname),exist=ann_arr%exists_yaml_file)
    if( ann_arr%exists_yaml_file) then
        call read_input_ann_yaml(parini,iproc,ann_arr)
    else 
        call read_input_ann(parini,iproc,ann_arr) 
    endif
    inquire(file="list_posinp_check.yaml",exist=file_exists)
    if(file_exists) then
        call read_data_yaml(parini,'list_posinp_check.yaml',atoms_check)
    else
        call read_data_old(parini,'list_posinp_check',atoms_check)
    endif
    !---------------------------------------------------------- 
    !open(unit=1,file='list_posinp_check',status='old',iostat=ios)
    !if(ios/=0) then;write(*,'(a)') 'ERROR: failure openning list_posinp_check__';stop;endif
    !nconftot=0
    !do
    !    read(1,'(a)',iostat=k) fn_tmp
    !    if(k<0) exit
    !    fn_tmp=adjustl(trim(filename))
    !    if(fn_tmp(1:1)=='#') cycle
    !    nconftot=nconftot+1
    !enddo
    !close(1)
    !----------------------------------------------------------   
    if(iproc==0) then
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
    gminarr(1:140)=1.d20 ; gmaxarr(1:140)=-1.d20
    iatmin(1:140)=0 ; iatmax(1:140)=0
    iconfmin(1:140)=0 ; iconfmax(1:140)=0
    if(.not. allocated(symfunc_check%symfunc)) then
        symfunc_check%nconf=atoms_check%nconf
        allocate(symfunc_check%symfunc(symfunc_check%nconf))
    endif
    do iconf=1,atoms_check%nconf
        symfunc_check%symfunc(iconf)%ng=ann_arr%ann(1)%nn(0) 
        symfunc_check%symfunc(iconf)%nat=atoms_check%atoms(iconf)%nat
        associate(ng=>symfunc_check%symfunc(iconf)%ng)
        associate(nat=>symfunc_check%symfunc(iconf)%nat)
        allocate(symfunc_check%symfunc(iconf)%y(ng,nat))
        end associate
        end associate
    enddo
!-----------------Compute symmetry functions with/without normalization-------------------------
    if(parini%normalization_ann) then
        configurations: do iconf=1,atoms_check%nconf
            call symmetry_functions(parini,ann_arr,atoms_check%atoms(iconf),symfunc,.false.)
            if(parini%symfunc_type_ann=='behler') then
                deallocate(symfunc%linked_lists%prime_bound)
                deallocate(symfunc%linked_lists%bound_rad)
                deallocate(symfunc%linked_lists%bound_ang)
            endif
            do iat=1,atoms_check%atoms(iconf)%nat
                do ig=1,symfunc_check%symfunc(iconf)%ng
                    symfunc_check%symfunc(iconf)%y(ig,iat)=symfunc%y(ig,iat)
                    if(symfunc_check%symfunc(iconf)%y(ig,iat)<gminarr(ig)) then
                        gminarr(ig)=symfunc_check%symfunc(iconf)%y(ig,iat)
                        iatmin(ig)=iat
                        iconfmin(ig)=iconf
                    endif
                    if(symfunc_check%symfunc(iconf)%y(ig,iat)>gmaxarr(ig)) then
                        gmaxarr(ig)=symfunc_check%symfunc(iconf)%y(ig,iat)
                        iatmax(ig)=iat
                        iconfmax(ig)=iconf
                    endif
                enddo
            enddo
            !write(*,*) allocated(symfunc%y),allocated(symfunc%y0d),allocated(symfunc%y0dr)
            call f_free(symfunc%y)
            call f_free(symfunc%y0d)
            call f_free(symfunc%y0dr)
        enddo configurations
        do i=1,ann_arr%nann
            do i0=1,ann_arr%ann(i)%nn(0)
                ann_arr%ann(i)%gbounds(1,i0)=gminarr(i0)
                ann_arr%ann(i)%gbounds(2,i0)=gmaxarr(i0)
                ann_arr%ann(i)%two_over_gdiff(i0)=2.d0/(ann_arr%ann(i)%gbounds(2,i0)-ann_arr%ann(i)%gbounds(1,i0))
            enddo
        enddo
        do iconf=1,atoms_check%nconf
            do iat=1,atoms_check%atoms(iconf)%nat
                ng=symfunc_check%symfunc(iconf)%ng
                i=atoms_check%atoms(iconf)%itypat(iat)
                do ig=1,ng
                    tt=symfunc_check%symfunc(iconf)%y(ig,iat)
                    tt=(tt-ann_arr%ann(i)%gbounds(1,ig))*ann_arr%ann(i)%two_over_gdiff(ig)-1.d0
                    symfunc_check%symfunc(iconf)%y(ig,iat)=tt
                    write(88,'(2i4,es24.15)') ig,iat,symfunc_check%symfunc(1)%y(ig,iat)
                enddo
            enddo
        enddo
    else
        do iconf=1,atoms_check%nconf
            call symmetry_functions(parini,ann_arr,atoms_check%atoms(iconf),symfunc,.false.)
            do iat=1,atoms_check%atoms(iconf)%nat
                do ig=1,symfunc_check%symfunc(iconf)%ng
                    symfunc_check%symfunc(iconf)%y(ig,iat)=symfunc%y(ig,iat)
                enddo
            enddo
        enddo
    endif
    !----------------------------------------------------------   
        fnout='incompatible'
        fnout1='distall'
        open(unit=111,file=fnout1,status='replace',iostat=ios)
        open(unit=11,file=fnout,status='replace',iostat=ios)
        if(ios/=0) then
            write(*,'(a)') 'ERROR: failure openning output file'
            stop
        endif
    !---------------------------------------------------------- 
    do iconf=1,atoms_check%nconf
    allocate(diff(symfunc_check%symfunc(iconf)%ng),c(atoms_check%atoms(iconf)%nat,atoms_check%atoms(iconf)%nat),F(atoms_check%atoms(iconf)%nat)) 
        do jconf=iconf+1,atoms_check%nconf
                do iat=1,atoms_check%atoms(iconf)%nat
                    do jat=1,atoms_check%atoms(jconf)%nat
                        do ig=1,symfunc_check%symfunc(iconf)%ng
                            diff(ig)=symfunc_check%symfunc(iconf)%y(ig,iat)-symfunc_check%symfunc(jconf)%y(ig,jat)
                        enddo
                        c(iat,jat)=dot_product(diff,diff)
                    enddo !over jat
                enddo !over iat
                call hung(atoms_check%atoms(iconf)%nat,c,F,distance2)
                distance=sqrt(distance2)
                de=abs(atoms_check%atoms(iconf)%epot-atoms_check%atoms(jconf)%epot)
                write(111,'(2(a25,i6,1x),2es20.10)') adjustl(trim(atoms_check%fn(iconf))),atoms_check%lconf(iconf), &
                    adjustl(trim(atoms_check%fn(jconf))),atoms_check%lconf(jconf),distance,de
                if(distance.le.dtol .and. de.ge.etol) then
                    write(11,'(2(a25,i6,1x),2es20.10)') adjustl(trim(atoms_check%fn(iconf))),atoms_check%lconf(iconf), &
                        adjustl(trim(atoms_check%fn(jconf))),atoms_check%lconf(jconf),distance,de
                    cycle
                else
                    write(12,'(2(a25,i6,1x),2es20.10)') trim(atoms_check%fn(iconf)),atoms_check%lconf(iconf), &
                        trim(atoms_check%fn(jconf)),atoms_check%lconf(jconf),distance,de
                endif
        enddo !over jconf
    deallocate(diff,c,F) 
    enddo !over iconf
    end associate
    close(11)
    close(111)
    call f_release_routine()
end subroutine ann_check_symmetry_function 
!*****************************************************************************************
