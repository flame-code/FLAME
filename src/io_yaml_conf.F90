!*****************************************************************************************
module mod_yaml_conf
    implicit none
    private
    public:: read_yaml_conf, write_yaml_conf
    public:: read_yaml_conf_getdict, read_yaml_conf_getatoms
    public:: atoms2dict
contains
!*****************************************************************************************
subroutine read_yaml_conf(parini,filename,nconfmax,atoms_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr
    use dictionaries
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: filename
    integer, intent(in):: nconfmax
    type(typ_atoms_arr), intent(inout):: atoms_arr
    !local variables
    type(dictionary), pointer :: confs_list=>null()
    if(nconfmax<1) then
        write(*,'(a)') 'ERROR: why do you call acf_read_new with nconfmax<1 ?'
        stop
    endif
    call read_yaml_conf_getdict(parini,filename,confs_list)
    call read_yaml_conf_getatoms(confs_list,nconfmax,atoms_arr)
    call dict_free(confs_list)
    nullify(confs_list)
    call yaml_mapping_open('Number of configurations read',flow=.true.)
    call yaml_map('filename',trim(filename))
    call yaml_map('nconf',atoms_arr%nconf)
    call yaml_mapping_close()
end subroutine read_yaml_conf
!*****************************************************************************************
subroutine read_yaml_conf_getdict(parini,filename,confs_list)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr
    use dictionaries
    use yaml_parse
    use dynamic_memory
    !use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: filename
    type(dictionary), pointer, intent(out) :: confs_list
    !local variables
    character(256):: fn_tmp
    character(256):: fn_fullpath
    character, dimension(:), allocatable :: fbuf
    integer(kind = 8) :: cbuf, cbuf_len
    fn_tmp=adjustl(trim(filename))
    if(fn_tmp(1:1)=='/') then
        fn_fullpath=trim(filename)
    else
        fn_fullpath=trim(parini%cwd)//'/'//trim(filename)
    endif
    call getFileContent(cbuf,cbuf_len,filename,len_trim(filename))
    fbuf=f_malloc0_str(1,int(cbuf_len),id='fbuf')
    call copyCBuffer(fbuf,cbuf,cbuf_len)
    call freeCBuffer(cbuf)
    !then parse the user's input file
    call yaml_parse_from_char_array(confs_list,fbuf)
    call f_free_str(1,fbuf)
end subroutine read_yaml_conf_getdict
!*****************************************************************************************
subroutine read_yaml_conf_getatoms(confs_list,nconfmax,atoms_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr, atom_allocate, update_rat
    use mod_const, only: bohr2ang
    use dictionaries
    use yaml_parse
    use dynamic_memory
    use mod_acf, only: str_motion2bemoved
    use yaml_output
    implicit none
    type(dictionary), pointer, intent(in):: confs_list
    integer, intent(in):: nconfmax
    type(typ_atoms_arr), intent(inout):: atoms_arr
    !local variables
    integer:: iconf, nconf, iat, nat, ii, iiconf
    character(3):: str_motion
    character(10):: str_units_length
    real(8):: x, y, z, cf_length, fx, fy, fz 
    real(8):: cvax, cvay, cvaz
    real(8):: cvbx, cvby, cvbz
    real(8):: cvcx, cvcy, cvcz
    type(dictionary), pointer :: dict1=>null()
    type(dictionary), pointer :: dict2=>null()
    nconf=dict_len(confs_list)
    atoms_arr%nconf=nconf
    if(nconf<1) stop 'ERROR: nconf<1 in read_yaml_conf'
    if(nconf>nconfmax) stop 'ERROR: nconf>nconfmax in read_yaml_conf'
    allocate(atoms_arr%atoms(nconf))

    do iconf=1,nconf
        iiconf=iconf-1
        dict1=>confs_list//iiconf//'conf'
        nat=dict1//'nat'
        call atom_allocate(atoms_arr%atoms(iconf),nat,0,0)
        atoms_arr%atoms(iconf)%boundcond=dict1//'bc'
        cf_length=1.d0
        if(has_key(dict1,"units_length")) then
            str_units_length=dict1//'units_length'
            atoms_arr%atoms(iconf)%units_length_io=trim(str_units_length)
            if(trim(str_units_length)=='angstrom') then
                cf_length=1.d0/bohr2ang
            elseif(trim(str_units_length)=='atomic') then
                cf_length=1.d0
            else
                write(*,*) 'ERROR: only atomic and angstrom are allowed in structural yaml input.'
                write(*,*) '       configuration= ',iconf
                stop
            endif
        endif
        dict2=>dict1//'cell'
        cvax=dict2//0//0 ; cvay=dict2//0//1 ; cvaz=dict2//0//2
        cvbx=dict2//1//0 ; cvby=dict2//1//1 ; cvbz=dict2//1//2
        cvcx=dict2//2//0 ; cvcy=dict2//2//1 ; cvcz=dict2//2//2
        nullify(dict2)
        atoms_arr%atoms(iconf)%cellvec(1,1)=cvax*cf_length
        atoms_arr%atoms(iconf)%cellvec(2,1)=cvay*cf_length
        atoms_arr%atoms(iconf)%cellvec(3,1)=cvaz*cf_length
        atoms_arr%atoms(iconf)%cellvec(1,2)=cvbx*cf_length
        atoms_arr%atoms(iconf)%cellvec(2,2)=cvby*cf_length
        atoms_arr%atoms(iconf)%cellvec(3,2)=cvbz*cf_length
        atoms_arr%atoms(iconf)%cellvec(1,3)=cvcx*cf_length
        atoms_arr%atoms(iconf)%cellvec(2,3)=cvcy*cf_length
        atoms_arr%atoms(iconf)%cellvec(3,3)=cvcz*cf_length
        if(has_key(dict1,"qtot")) then
            atoms_arr%atoms(iconf)%qtot=dict1//'qtot'
        endif
        if(has_key(dict1,"dpm")) then
            atoms_arr%atoms(iconf)%dpm=dict1//'dpm'
        endif
        if(has_key(dict1,"elecfield")) then
            atoms_arr%atoms(iconf)%elecfield=dict1//'elecfield'
        endif
        if(has_key(dict1,"epot")) then
            atoms_arr%atoms(iconf)%epot=dict1//'epot'
        endif
        if(dict_len(dict1//'coord'//0)==5) then
            do iat=1,nat
                ii=iat-1
                dict2=>dict1//'coord'//ii
                x=dict2//0
                y=dict2//1
                z=dict2//2
                atoms_arr%atoms(iconf)%ratp(1,iat)=x*cf_length
                atoms_arr%atoms(iconf)%ratp(2,iat)=y*cf_length
                atoms_arr%atoms(iconf)%ratp(3,iat)=z*cf_length
                atoms_arr%atoms(iconf)%sat(iat)=dict1//'coord'//ii//3
                str_motion=dict1//'coord'//ii//4
                call str_motion2bemoved(str_motion,atoms_arr%atoms(iconf)%bemoved(1,iat))
                nullify(dict2)
            enddo
            call update_rat(atoms_arr%atoms(iconf),upall=.true.)
        else
            do iat=1,nat
                ii=iat-1
                dict2=>dict1//'coord'//ii
                x=dict2//0
                y=dict2//1
                z=dict2//2
                atoms_arr%atoms(iconf)%ratp(1,iat)=x*cf_length
                atoms_arr%atoms(iconf)%ratp(2,iat)=y*cf_length
                atoms_arr%atoms(iconf)%ratp(3,iat)=z*cf_length
                atoms_arr%atoms(iconf)%sat(iat)=dict1//'coord'//ii//3
                nullify(dict2)
                atoms_arr%atoms(iconf)%bemoved(1,iat)=.true.
                atoms_arr%atoms(iconf)%bemoved(2,iat)=.true.
                atoms_arr%atoms(iconf)%bemoved(3,iat)=.true.
            enddo
            call update_rat(atoms_arr%atoms(iconf),upall=.true.)
        endif
        if(has_key(dict1,"force")) then
            do iat=1,nat
                ii=iat-1
                dict2=>dict1//'force'//ii
                fx=dict2//0
                fy=dict2//1
                fz=dict2//2
                atoms_arr%atoms(iconf)%fat(1,iat)=fx
                atoms_arr%atoms(iconf)%fat(2,iat)=fy
                atoms_arr%atoms(iconf)%fat(3,iat)=fz
                nullify(dict2)
            enddo
        endif
        nullify(dict1)
    enddo
end subroutine read_yaml_conf_getatoms
!*****************************************************************************************
subroutine write_yaml_conf(file_info,atoms,strkey)
    use mod_atoms, only: typ_file_info, typ_atoms
    use dictionaries
    use futile
    use yaml_parse
    use dynamic_memory
    use yaml_output
    implicit none
    type(typ_file_info), intent(inout):: file_info
    type(typ_atoms), intent(in):: atoms
    character(*), optional, intent(in):: strkey
    !local variables
    integer:: iunit, ierr
    character(256):: str_msg
    type(dictionary), pointer :: dict1=>null()

    dict1=>dict_new()

    call set(dict1//'conf',dict_new())

    call atoms2dict(file_info,atoms,dict1)

    iunit=f_get_free_unit(10**5)

    if(trim(file_info%file_position)=='new') then
        call yaml_set_stream(unit=iunit,filename=trim(file_info%filename_positions),&
             record_length=92,istat=ierr,setdefault=.false.,tabbing=0,position='rewind')
    elseif(trim(file_info%file_position)=='append') then
        call yaml_set_stream(unit=iunit,filename=trim(file_info%filename_positions),&
             record_length=92,istat=ierr,setdefault=.false.,tabbing=0,position='append')
    endif
    if(ierr/=0) then
        str_msg='Failed to create'//trim(file_info%filename_positions)
        str_msg=trim(str_msg)//'error code='//trim(yaml_toa(ierr))
       call yaml_warning(trim(str_msg))
    end if
    call yaml_release_document(unit=iunit)
    call yaml_new_document(unit=iunit)

    call yaml_dict_dump(dict1,unit=iunit)
    call dict_free(dict1)
    nullify(dict1)
    call yaml_close_stream(unit=iunit)
end subroutine write_yaml_conf
!*****************************************************************************************
subroutine atoms2dict(file_info,atoms,dict1)
    use mod_atoms, only: typ_file_info, typ_atoms, get_rat
    use mod_const, only: bohr2ang
    use dictionaries
    use futile
    use yaml_parse
    use dynamic_memory
    use yaml_output
    implicit none
    type(typ_file_info), intent(inout):: file_info
    type(typ_atoms), intent(in):: atoms
    type(dictionary), pointer, intent(inout):: dict1
    !local variables
    integer:: iat, ii
    character(3):: str_motion
    real(8):: x, y, z, cf_length, fx, fy, fz
    real(8):: cvax, cvay, cvaz
    real(8):: cvbx, cvby, cvbz
    real(8):: cvcx, cvcy, cvcz
    type(dictionary), pointer :: single_atom_list=>null()
    type(dictionary), pointer :: coord_list=>null()
    type(dictionary), pointer :: conf_dict=>null()
    type(dictionary), pointer :: force_list=>null()
    type(dictionary), pointer :: stress_list=>null()
    real(8), allocatable:: rat(:,:)
    allocate(rat(3,atoms%nat))

    conf_dict=>dict1//'conf'

    call set(conf_dict//'nat',atoms%nat)
    call set(conf_dict//'bc',trim(atoms%boundcond))
    cf_length=1.d0
    if(trim(atoms%units_length_io)=='angstrom') then
        call set(conf_dict//'units_length',trim(atoms%units_length_io))
        if(trim(atoms%units_length_io)=='angstrom') then
            cf_length=bohr2ang
        endif
    endif
    cvax=atoms%cellvec(1,1)*cf_length
    cvay=atoms%cellvec(2,1)*cf_length
    cvaz=atoms%cellvec(3,1)*cf_length
    cvbx=atoms%cellvec(1,2)*cf_length
    cvby=atoms%cellvec(2,2)*cf_length
    cvbz=atoms%cellvec(3,2)*cf_length
    cvcx=atoms%cellvec(1,3)*cf_length
    cvcy=atoms%cellvec(2,3)*cf_length
    cvcz=atoms%cellvec(3,3)*cf_length
    call set(conf_dict//'cell'//0,(/cvax,cvay,cvaz/))
    call set(conf_dict//'cell'//1,(/cvbx,cvby,cvbz/))
    call set(conf_dict//'cell'//2,(/cvcx,cvcy,cvcz/))

    if(abs(atoms%qtot)>1.d-16) then
        call set(conf_dict//'qtot',atoms%qtot)
    endif
    call set(conf_dict//'epot',atoms%epot)
    call set(conf_dict//'dpm',(/atoms%dpm(1),atoms%dpm(2),atoms%dpm(3)/))
    call set(conf_dict//'coord',list_new(.item. "it will be overwritten"))
    coord_list=>conf_dict//'coord'
    call get_rat(atoms,rat)
    do iat=1,atoms%nat
        ii=iat-1
        x=rat(1,iat)*cf_length
        y=rat(2,iat)*cf_length
        z=rat(3,iat)*cf_length
        call set(coord_list//ii,list_new(.item. "it will be overwritten"))
        single_atom_list=>coord_list//ii
        call set(single_atom_list,(/x,y,z/))
        call set(single_atom_list//3,trim(atoms%sat(iat)))
        write(str_motion,'(3l1)') atoms%bemoved(1,iat),atoms%bemoved(2,iat),atoms%bemoved(3,iat)
        call set(single_atom_list//4,str_motion)
        nullify(single_atom_list)
    enddo
    nullify(coord_list)

    if(file_info%print_force) then !CORRECT_IT
        force_list=>conf_dict//'force'
        do iat=1,atoms%nat
            ii=iat-1
            fx=atoms%fat(1,iat)
            fy=atoms%fat(2,iat)
            fz=atoms%fat(3,iat)
            call set(force_list//ii,list_new(.item. "it will be overwritten"))
            single_atom_list=>force_list//ii
            call set(single_atom_list,(/fx,fy,fz/))
            nullify(single_atom_list)
        enddo
        nullify(force_list)
        if(trim(atoms%boundcond)=='bulk') then
            stress_list=>conf_dict//'stress'
            call set(stress_list,(/atoms%stress(1,1), &
                                   atoms%stress(2,1), &
                                   atoms%stress(3,1), &
                                   atoms%stress(1,2), &
                                   atoms%stress(2,2), &
                                   atoms%stress(3,2), &
                                   atoms%stress(1,3), &
                                   atoms%stress(2,3), &
                                   atoms%stress(3,3)/))
            nullify(stress_list)
        endif
    endif
    nullify(conf_dict)
    deallocate(rat)
end subroutine atoms2dict
!*****************************************************************************************
end module mod_yaml_conf
!*****************************************************************************************
