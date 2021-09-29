!*****************************************************************************************
subroutine read_input_ann(parini,iproc,ann_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iproc
    type(typ_ann_arr), intent(inout):: ann_arr
    !local variables
    integer:: ios, iann, i, j
    character(256):: fn_fullpath
    character(5):: stypat
    real(8):: rcut
    do iann=1,ann_arr%nann
        if(parini%bondbased_ann) then
            stypat=parini%stypat(1)
        else
            stypat=parini%stypat(iann)
        endif
        fn_fullpath=trim(parini%cwd)//'/'//'input.ann.'//trim(stypat)
        open(unit=1,file=trim(fn_fullpath),status='old',iostat=ios)
        if(ios/=0) then
            write(*,'(2a)') 'ERROR: failure openning ','input.ann.'//trim(stypat)
            stop
        endif
        call read_symmetry_functions(parini,iproc,1,ann_arr%ann(iann),rcut)
        if(iann==1) then
            ann_arr%rcut=rcut
        else
            if(abs(rcut-ann_arr%rcut)>epsilon(1.d0)) then
                write(*,'(a,2f6.3)') 'ERROR: rcut is different in input.ann.* files: ',ann_arr%rcut,rcut
                stop
            endif
        endif
        close(1)
    enddo
    !if(.not. (parini%bondbased_ann .and. trim(ann_arr%approach)=='tb')) then
    !    do i=1,ann_arr%nann
    !        do j=i,ann_arr%nann
    !            ann_arr%reprcut(i,j)=ann_arr%ann(i)%rionic+ann_arr%ann(j)%rionic
    !            ann_arr%reprcut(j,i)=ann_arr%ann(i)%rionic+ann_arr%ann(j)%rionic
    !        enddo
    !    enddo
    !endif
end subroutine read_input_ann
!*****************************************************************************************
subroutine read_symmetry_functions(parini,iproc,ifile,ann,rcut)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iproc, ifile
    type(typ_ann), intent(inout):: ann
    real(8), intent(out):: rcut
    !local variables
    integer:: ig, i, icol, ios, i0, iline
    character(256):: strline
    character(50):: str1, str2, str3
    character(5):: sat1, sat2
    !logical:: all_read
    !read(ifile,*) ann%nn(1),ann%nn(2)
    do iline=1,10
        read(ifile,'(a)') ann%hlines(iline)
    enddo
    strline=trim(ann%hlines(1))
    ann%nl=0
    !all_read=.false.
    do icol=1,10-1
        call elim_white_space(strline)
        if(strline(1:1)=='#') exit
        read(strline,*,iostat=ios) ann%nn(icol)
        write(*,'(i8,5x,a)') ios,trim(strline)
        if(ios<0) then
            exit
        else
            ann%nl=ann%nl+1
            !write(*,*) ann%nl,ann%nn(icol)
        endif
        if(ann%nn(icol)==0) then 
            write(*,'(a,i5,a)') 'ERROR: number of nodes at',icol,'th hidden layer is zero'
            stop
        endif
        do i=1,256
            if(strline(i:i)==' ' .or. strline(i:i)=='' .or. strline(i:i)=='#') then
                exit
            else
                strline(i:i)=' '
            endif
        enddo
        !if(all_read) exit
    enddo
    ann%nl=ann%nl+1 !adding the output layer to total number of layers
    ann%nn(ann%nl)=1 !setting the output layer
    read(ann%hlines(2),*) str1,rcut
    if(trim(parini%approach_ann)/='atombased' .and. trim(parini%approach_ann)/='tb' ) then
        read(ann%hlines(3),*) str1,ann%ampl_chi,str2,ann%prefactor_chi
        read(ann%hlines(4),*) str1,ann%zion,str2,ann%gausswidth_ion,str3,ann%ener_ref
        read(ann%hlines(5),*) str1,ann%gausswidth,str2,ann%hardness,str3,ann%chi0
        read(ann%hlines(6),*) str1,ann%spring_const,str2,ann%qinit
    elseif(trim(parini%approach_ann)=='tb') then
        read(ann%hlines(4),*) str3,ann%ener_ref
    endif
    !read(ann%hlines(6),*) str1,ann%rionic
    !---------------------------------------------
    i0=0

    read(ifile,'(a)') strline
    read(strline,*) ann%ng1
    if(.not. parini%bondbased_ann .and. ann%ng1>0) then
        stop 'ERROR: symmetry function of type G3 not implemented yet.'
    endif
    do ig=1,ann%ng1
        if(.not. parini%bondbased_ann) then
            stop 'ERROR: g1 is not ready.'
        endif
        i0=i0+1
        read(ifile,'(a)') strline
        read(strline,*,iostat=ios) ann%g1eta(ig),ann%g1rs(ig),ann%gbounds(1,i0),ann%gbounds(2,i0)
        if(ios<0) then
            write(*,'(a)') 'ERROR: 5 columns are required for each of G1 symmetry functions,'
            write(*,'(a)') '       including 2 values for bounds (neglected for ANN training).'
            stop
        endif
    enddo
    read(ifile,'(a)') strline
    read(strline,*) ann%ng2
    do ig=1,ann%ng2
        i0=i0+1
        read(ifile,'(a)') strline
        read(strline,*,iostat=ios) sat1,ann%g2eta(ig),ann%g2rs(ig),ann%gbounds(1,i0),ann%gbounds(2,i0)
        call set_radial_atomtype(parini,sat1,ann%g2i(ig))
        if(ios<0) then
            write(*,'(a)') 'ERROR: 5 columns are required for each of G2 symmetry functions,'
            write(*,'(a)') '       including 2 values for bounds (neglected for ANN training).'
            stop
        endif
    enddo
    read(ifile,'(a)') strline
    read(strline,*) ann%ng3
    if(ann%ng3>0) stop 'ERROR: symmetry function of type G3 not implemented yet.'
    read(ifile,'(a)') strline
    read(strline,*) ann%ng4
    do ig=1,ann%ng4
        stop 'ERROR: g4 is not ready.'
        i0=i0+1
        read(ifile,'(a)') strline
        read(strline,*,iostat=ios) ann%g4eta(ig),ann%g4zeta(ig), &
            ann%g4lambda(ig),ann%gbounds(1,i0),ann%gbounds(2,i0)
        if(ios<0) then
            write(*,'(a)') 'ERROR: 6 columns are required for each of G4 symmetry functions,'
            write(*,'(a)') '       including 2 values for bounds (neglected for ANN training).'
            stop
        endif
    enddo
    read(ifile,'(a)') strline
    read(strline,*) ann%ng5
    do ig=1,ann%ng5
        i0=i0+1
        read(ifile,'(a)') strline
        read(strline,*,iostat=ios) sat1,sat2,ann%g5eta(ig),ann%g5zeta(ig), &
            ann%g5lambda(ig),ann%gbounds(1,i0),ann%gbounds(2,i0)
        call set_angular_atomtype(parini,sat1,sat2,ann%g5i(1,ig))
        if(ios<0) then
            write(*,'(a)') 'ERROR: 6 columns are required for each of G5 symmetry functions,'
            write(*,'(a)') '       including 2 values for bounds (neglected for ANN training).'
            stop
        endif
    enddo
    !read(ifile,'(a)') strline
    !read(strline,*) ann%ng6
    !if(mod(ann%ng6,3)/=0) stop 'ERROR: ng6 must be multiple of three.'
    !do ig=1,ann%ng6/3
    !    stop 'ERROR: g6 is not ready.'
    !    i0=i0+1
    !    read(ifile,'(a)') strline
    !    read(strline,*,iostat=ios) ann%g6eta(ig),ann%gbounds(1,i0),ann%gbounds(2,i0)
    !    if(ios<0) then
    !        write(*,'(a)') 'ERROR: 4 columns are required for each of G6 symmetry functions,'
    !        write(*,'(a)') '       including 2 values for bounds (neglected for ANN training).'
    !        stop
    !    endif
    !    i0=i0+1
    !    read(ifile,'(a)') strline
    !    read(strline,*,iostat=ios) ann%gbounds(1,i0),ann%gbounds(2,i0)
    !    if(ios<0) then
    !        write(*,'(a)') 'ERROR: 2 columns are required for each of G6 symmetry functions,'
    !        write(*,'(a)') '       including 2 values for bounds (neglected for ANN training).'
    !        stop
    !    endif
    !    i0=i0+1
    !    read(ifile,'(a)') strline
    !    read(strline,*,iostat=ios) ann%gbounds(1,i0),ann%gbounds(2,i0)
    !    if(ios<0) then
    !        write(*,'(a)') 'ERROR: 2 columns are required for each of G6 symmetry functions,'
    !        write(*,'(a)') '       including 2 values for bounds (neglected for ANN training).'
    !        stop
    !    endif
    !enddo
    ann%ng6=0
    ann%nn(0)=ann%ng1+ann%ng2+ann%ng3+ann%ng4+ann%ng5+ann%ng6
    if(iproc==0) then
        do i=0,ann%nl
            write(*,'(a,i1,a,i4)') 'ann%(',i,')=',ann%nn(i)
        enddo
        !write(*,'(a,3i4)') 'n0,n1,n2 ',ann%nn(0),ann%nn(1),ann%nn(2)
        write(*,'(a,6i4)') 'ng1,ng2,ng3,ng4,ng5,ng6 ',ann%ng1,ann%ng2,ann%ng3,ann%ng4,ann%ng5,ann%ng6
    endif
end subroutine read_symmetry_functions
!*****************************************************************************************
subroutine set_radial_atomtype(parini,sat1,ityp)
    use mod_parini, only: typ_parini
    implicit none
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: sat1
    integer, intent(out):: ityp(1)
    !local variables
    integer:: i
    character(200):: str
    str='ERROR: set_radial_atomtype: types in input.ini and input.ann.* are inconsistent:'
    do i=1,parini%ntypat
        !write(*,*) trim(parini%stypat(i)),trim(sat1)
        if(trim(parini%stypat(i))==trim(sat1)) then
            ityp(1)=parini%ltypat(i)
            exit
        endif
        if(i==parini%ntypat) then
            write(*,'(a,2(1x,a))') trim(str),trim(parini%stypat(i)),trim(sat1)
            stop
        endif
    enddo
end subroutine set_radial_atomtype
!*****************************************************************************************
subroutine set_angular_atomtype(parini,sat1,sat2,ityp)
    use mod_parini, only: typ_parini
    implicit none
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: sat1, sat2
    integer, intent(out):: ityp(2)
    !local variables
    integer:: i, j
    character(200):: str
    str='ERROR: set_angular_atomtype: types in input.ini and input.ann.* are inconsistent.'
    firstloop: do i=1,parini%ntypat
        do j=1,parini%ntypat
            if(trim(parini%stypat(i))==trim(sat1) .and. trim(parini%stypat(j))==trim(sat2)) then
                ityp(1)=parini%ltypat(i)
                ityp(2)=parini%ltypat(j)
                exit firstloop
            endif
        enddo
        if(i==parini%ntypat .and. j==parini%ntypat) then
            write(*,'(a)') trim(str)
            stop
        endif
    enddo firstloop
end subroutine set_angular_atomtype
!*****************************************************************************************
subroutine write_ann_all(parini,ann_arr,iter)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(in):: ann_arr
    integer, intent(in):: iter
    !local variables
    character(16):: fn
    character(1):: fn_tt
    character(50):: filename
    integer:: i
    if(iter==-1) then
        write(fn,'(a10)') '.ann.param'
    else
        write(fn,'(a11,i5.5)') '.ann.param.',iter
    endif
    if(parini%bondbased_ann .and. trim(ann_arr%approach)=='tb') then
        if(parini%ntypat>1) then
            stop 'ERROR: writing ANN parameters for tb available only ntypat=1'
        endif
        do i=1,ann_arr%nann
            write(fn_tt,'(i1)') i
            filename=trim(parini%stypat(1))//fn_tt//trim(fn)
            write(*,'(a)') trim(filename)
            call write_ann(parini,filename,ann_arr%ann(i))
        enddo
    elseif(trim(ann_arr%approach)=='eem1' .or. trim(ann_arr%approach)=='cent1' &
        .or. trim(ann_arr%approach)=='cent2' .or. trim(ann_arr%approach)=='centt') then
        do i=1,ann_arr%nann
            filename=trim(parini%stypat(i))//trim(fn)
            write(*,'(a)') trim(filename)
            call write_ann(parini,filename,ann_arr%ann(i))
        enddo
    else
        stop 'ERROR: writing ANN parameters is only for cent1,centt,tb'
    endif
end subroutine write_ann_all
!*****************************************************************************************
subroutine write_ann(parini,filename,ann)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann
    implicit none
    type(typ_parini), intent(in):: parini
    character(*):: filename
    type(typ_ann), intent(in):: ann
    !local variables
    !integer:: 
    !real(8):: 
    integer:: i, j, k, l, ios, ialpha, i0, iline
    character(5):: sat1, sat2
    open(unit=1,file=trim(filename),status='replace',iostat=ios)
    if(ios/=0) then;write(*,'(2a)') 'ERROR: failure openning ',trim(filename);stop;endif
    !write(1,'(2i6,2x,a)') ann%nn(1),ann%nn(2),'#number of nodes in 1st and 2nd hidden layers.'
    do iline=1,10
        write(1,'(a)') trim(ann%hlines(iline))
    enddo
    i0=0
    !-------------------------------------------------------
    write(1,'(i6,2x,a)') ann%ng1,'#ng1'
    do i=1,ann%ng1
        i0=i0+1
        write(1,'(2f8.4,2es24.15)') ann%g1eta(i),ann%g1rs(i),ann%gbounds(1,i0),ann%gbounds(2,i0)
    enddo
    !-------------------------------------------------------
    write(1,'(i6,2x,a)') ann%ng2,'#ng2'
    do i=1,ann%ng2
        i0=i0+1
        sat1=parini%stypat(ann%g2i(i))
        write(1,'(1a5,2f8.4,2es24.15)') trim(sat1),ann%g2eta(i),ann%g2rs(i),ann%gbounds(1,i0),ann%gbounds(2,i0)
    enddo
    !-------------------------------------------------------
    write(1,'(i6,2x,a)') ann%ng3,'#ng3'
    do i=1,ann%ng3
        i0=i0+1
        write(1,'(1f8.4)') ann%g3kappa(i)
    enddo
    !-------------------------------------------------------
    write(1,'(i6,2x,a)') ann%ng4,'#ng4'
    do i=1,ann%ng4
        i0=i0+1
        write(1,'(3f8.4,2es24.15)') ann%g4eta(i),ann%g4zeta(i),ann%g4lambda(i),ann%gbounds(1,i0),ann%gbounds(2,i0)
    enddo
    !-------------------------------------------------------
    write(1,'(i6,2x,a)') ann%ng5,'#ng5'
    do i=1,ann%ng5
        i0=i0+1
        sat1=parini%stypat(ann%g5i(1,i))
        sat2=parini%stypat(ann%g5i(2,i))
        write(1,'(2a5,3f8.4,2es24.15)') trim(sat1),trim(sat2),ann%g5eta(i),ann%g5zeta(i),ann%g5lambda(i),ann%gbounds(1,i0),ann%gbounds(2,i0)
    enddo
    !-------------------------------------------------------
    !write(1,'(i6,2x,a)') ann%ng6,'#ng6'
    !do i=1,ann%ng6/3
    !    i0=i0+1
    !    write(1,'(1f8.4,2es24.15)') ann%g6eta(i),ann%gbounds(1,i0),ann%gbounds(2,i0)
    !    i0=i0+1
    !    write(1,'(16x,2es24.15)') ann%gbounds(1,i0),ann%gbounds(2,i0)
    !    i0=i0+1
    !    write(1,'(16x,2es24.15)') ann%gbounds(1,i0),ann%gbounds(2,i0)
    !enddo
    !-------------------------------------------------------
    do ialpha=1,ann%nl
        write(1,'(2(a,i1))') '#main nodes weights connecting layers ',ialpha,' and ',ialpha-1
        do j=1,ann%nn(ialpha)
            do i=1,ann%nn(ialpha-1)
                write(1,'(es24.15)') ann%a(i,j,ialpha)
            enddo
        enddo
        write(1,'(a,i1)') '#bias nodes weights for layer ',ialpha
        do i=1,ann%nn(ialpha)
            write(1,'(es24.15)') ann%b(i,ialpha)
        enddo
    enddo
    !-------------------------------------------------------
    close(1)
end subroutine write_ann
!*****************************************************************************************
subroutine read_ann(parini,ann_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_processors, only: iproc
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    !local variables
    !integer:: 
    integer:: i, j, k, l, ios, i0, ifile, ialpha, iann
    !character(100):: ttstr
    real(8):: bound_l, bound_u, rcut
    character(16):: fn
    character(1):: fn_tt
    character(50):: filename
    do iann=1,ann_arr%nann
        write(fn,'(a10)') '.ann.param'
        if(parini%bondbased_ann .and. trim(ann_arr%approach)=='tb') then
            if(parini%ntypat>1) then
                stop 'ERROR: writing ANN parameters for tb available only ntypat=1'
            endif
            write(fn_tt,'(i1)') iann
            filename=trim(parini%stypat(1))//fn_tt//trim(fn)
            write(*,'(a)') trim(filename)
        elseif(trim(ann_arr%approach)=='eem1' .or. trim(ann_arr%approach)=='cent1' &
            .or. trim(ann_arr%approach)=='cent2'.or. trim(ann_arr%approach)=='centt') then
            filename=trim(parini%stypat(iann))//trim(fn)
        else
            stop 'ERROR: reading ANN parameters is only for cent1,centt,tb'
        endif
        open(unit=1,file=trim(filename),status='old',iostat=ios)
        if(ios/=0) then
            write(*,'(2a)') 'ERROR: failure openning ',trim(filename)
            stop
        endif
        !-------------------------------------------------------
        !call read_symmetry_functions(parini,iproc,1,ann)
        ifile=1
        call read_symmetry_functions(parini,iproc,1,ann_arr%ann(iann),rcut)
        if(iann==1) then
            ann_arr%rcut=rcut
        else
            if(abs(rcut-ann_arr%rcut)>epsilon(1.d0)) then
                write(*,'(a,2f6.3)') 'ERROR: rcut is different in input.ann.* files: ',ann_arr%rcut,rcut
                stop
            endif
        endif
        do i0=1,ann_arr%ann(iann)%nn(0)
            bound_l=ann_arr%ann(iann)%gbounds(1,i0)
            bound_u=ann_arr%ann(iann)%gbounds(2,i0)
            ann_arr%ann(iann)%two_over_gdiff(i0)=2.d0/(bound_u-bound_l)
        enddo
        !-------------------------------------------------------
        do ialpha=1,ann_arr%ann(iann)%nl
            read(1,*)
            do j=1,ann_arr%ann(iann)%nn(ialpha)
                do i=1,ann_arr%ann(iann)%nn(ialpha-1)
                    read(1,*) ann_arr%ann(iann)%a(i,j,ialpha)
                enddo
            enddo
            read(1,*)
            do i=1,ann_arr%ann(iann)%nn(ialpha)
                read(1,*) ann_arr%ann(iann)%b(i,ialpha)
            enddo
        enddo
        !-------------------------------------------------------
        close(1)
    enddo
    !if(.not. (parini%bondbased_ann .and. trim(ann_arr%approach)=='tb')) then
    !    do i=1,ann_arr%nann
    !        do j=i,ann_arr%nann
    !            ann_arr%reprcut(i,j)=ann_arr%ann(i)%rionic+ann_arr%ann(j)%rionic
    !            ann_arr%reprcut(j,i)=ann_arr%ann(i)%rionic+ann_arr%ann(j)%rionic
    !        enddo
    !    enddo
    !endif
end subroutine read_ann
!*****************************************************************************************
subroutine read_data_old(parini,filename_list,atoms_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr, atom_copy_old, set_rat_atoms
    use mod_atoms, only: atom_allocate_old, atom_deallocate, atom_deallocate_old
    use mod_acf, only: acf_read_new
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: filename_list
    type(typ_atoms_arr), intent(inout):: atoms_arr
    !local variables
    integer:: i, iat, ios, k, iconf
    character(256):: filename, fn_tmp, filename_force
    type(typ_atoms_arr):: atoms_arr_of !configuration of one file
    type(typ_atoms_arr):: atoms_arr_t
    real(8):: fx, fy, fz
    integer:: nconfmax, ind, len_filename
    call f_routine(id='read_data_old')
    nconfmax=1*10**5
    allocate(atoms_arr_t%atoms(nconfmax))
    allocate(atoms_arr_t%fn(nconfmax))
    allocate(atoms_arr_t%lconf(nconfmax))
    open(unit=1,file=filename_list,status='old',iostat=ios)
    if(ios/=0) then;write(*,'(a)') 'ERROR: failure openning list_posinp';stop;endif
    atoms_arr_t%nconf=0
    over_files: do
        read(1,'(a)',iostat=k) filename
        if(k<0) exit
        fn_tmp=adjustl(trim(filename))
        if(fn_tmp(1:1)=='#') cycle
        call acf_read_new(parini,filename,10000,atoms_arr_of)
        ind=index(filename,'/',back=.true.)
        len_filename=len(filename)
        filename_force=filename(1:ind)//'force_'//filename(ind+1:len_filename)
        if(parini%read_forces_ann) then
            open(unit=2,file=trim(filename_force),status='old',iostat=ios)
        endif
        if(ios/=0) then;write(*,'(a)') 'ERROR: failure openning force of list_posinp';stop;endif
        over_iconf: do iconf=1,atoms_arr_of%nconf
            atoms_arr_t%nconf=atoms_arr_t%nconf+1
            if(atoms_arr_t%nconf>nconfmax) then
                stop 'ERROR: too many configurations, change parameter nconfmax.'
            endif
            call atom_allocate_old(atoms_arr_t%atoms(atoms_arr_t%nconf),atoms_arr_of%atoms(iconf)%nat,0,0)
            atoms_arr_t%atoms(atoms_arr_t%nconf)%epot=atoms_arr_of%atoms(iconf)%epot
            atoms_arr_t%atoms(atoms_arr_t%nconf)%qtot=atoms_arr_of%atoms(iconf)%qtot
            atoms_arr_t%atoms(atoms_arr_t%nconf)%boundcond=trim(atoms_arr_of%atoms(iconf)%boundcond)
            atoms_arr_t%atoms(atoms_arr_t%nconf)%cellvec(1:3,1:3)=atoms_arr_of%atoms(iconf)%cellvec(1:3,1:3)
            if(parini%read_forces_ann) read(2,*)
            call set_rat_atoms(atoms_arr_t%atoms(atoms_arr_t%nconf),atoms_arr_of%atoms(iconf),setall=.true.)
            do iat=1,atoms_arr_of%atoms(iconf)%nat
                atoms_arr_t%atoms(atoms_arr_t%nconf)%sat(iat)=atoms_arr_of%atoms(iconf)%sat(iat)
                if(parini%read_forces_ann) then
                    read(2,*) fx,fy,fz
                    atoms_arr_t%atoms(atoms_arr_t%nconf)%fat(1,iat)=fx
                    atoms_arr_t%atoms(atoms_arr_t%nconf)%fat(2,iat)=fy
                    atoms_arr_t%atoms(atoms_arr_t%nconf)%fat(3,iat)=fz
                else
                    atoms_arr_t%atoms(atoms_arr_t%nconf)%fat(1,iat)=0.d0
                    atoms_arr_t%atoms(atoms_arr_t%nconf)%fat(2,iat)=0.d0
                    atoms_arr_t%atoms(atoms_arr_t%nconf)%fat(3,iat)=0.d0
                endif
            enddo
            atoms_arr_t%fn(atoms_arr_t%nconf)=trim(filename)
            atoms_arr_t%lconf(atoms_arr_t%nconf)=iconf
            call atom_deallocate(atoms_arr_of%atoms(iconf))
        enddo over_iconf
        deallocate(atoms_arr_of%atoms)
        if(parini%read_forces_ann) close(2)
    enddo over_files
    close(1)
    atoms_arr%nconf=atoms_arr_t%nconf
    allocate(atoms_arr%atoms(atoms_arr%nconf))
    allocate(atoms_arr%fn(atoms_arr%nconf))
    allocate(atoms_arr%lconf(atoms_arr%nconf))
    do iconf=1,atoms_arr%nconf
        call atom_copy_old(atoms_arr_t%atoms(iconf),atoms_arr%atoms(iconf), &
            'atoms_arr_t%atoms(iconf)->atoms_arr%atoms(iconf)')
        atoms_arr%fn(iconf)=trim(atoms_arr_t%fn(iconf))
        atoms_arr%lconf(iconf)=atoms_arr_t%lconf(iconf)
    enddo
    do iconf=1,atoms_arr_t%nconf
        call atom_deallocate_old(atoms_arr_t%atoms(iconf))
    enddo
    deallocate(atoms_arr_t%atoms)
    deallocate(atoms_arr_t%fn)
    deallocate(atoms_arr_t%lconf)
    call f_release_routine()
end subroutine read_data_old
!*****************************************************************************************
