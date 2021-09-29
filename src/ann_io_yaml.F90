!*****************************************************************************************
subroutine read_input_ann_yaml(parini,iproc,ann_arr)
    use futile
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iproc
    type(typ_ann_arr), intent(inout):: ann_arr
    !local variables
    integer:: ios, iann, i, j
    character(50):: fname,str1
    character(5):: stypat
    real(8)::rcut
    !call f_lib_initialize()
    !call yaml_new_document()
    call yaml_sequence_open('ann input files') !,flow=.true.)
    do iann=1,ann_arr%nann
        if(parini%bondbased_ann) then
            stypat=parini%stypat(1)
        else
            stypat=parini%stypat(iann)
        endif
        !-------------------------------------------------------
        fname = trim(stypat)//'.ann.input.yaml'
        call yaml_sequence(advance='no')
        call yaml_comment(trim(fname))
        !write(*,*)trim(fname)
        call get_symfunc_parameters_yaml(parini,iproc,fname,ann_arr%ann(iann),rcut)
        ann_arr%rcut = rcut
        !-------------------------------------------------------
        call dict_free(ann_arr%ann(iann)%dict)
        nullify(ann_arr%ann(iann)%dict)
    enddo
    call yaml_sequence_close()
    !call f_lib_finalize()
end subroutine read_input_ann_yaml
!*****************************************************************************************
subroutine get_symfunc_parameters_yaml(parini,iproc,fname,ann,rcut)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann
    use dictionaries
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann), intent(inout):: ann
    integer, intent(in):: iproc
    real(8), intent(out):: rcut
    !local variables
    integer:: ig, ios, i0, i, il 
    character(50):: fname, sat1, sat2
    character(250):: tt, str1
    integer :: count1, count2, count3, count4, count5, count6
    character(5):: stypat
    type(dictionary), pointer :: subdict_ann=>null()
    type(dictionary), pointer :: dict_tmp=>null()
    character(50):: str_out_ann
    character(5):: str_out_ann_tt
    call set_dict_ann(ann,fname,stypat)
    !call yaml_comment('USER INPUT FILE',hfill='~')
    if(parini%iverbose>=2) call yaml_dict_dump(ann%dict)
    !call yaml_comment('',hfill='~')
    !logical:: all_read

    subdict_ann => ann%dict//"main"

    ann%nl=dict_len(subdict_ann//"nodes")
    do il=0,ann%nl-1
        ann%nn(il+1)=subdict_ann//"nodes"//il
    enddo
    ann%nl=ann%nl+1 !adding the output layer to total number of layers
    ann%nn(ann%nl)=1 !setting the output layer
    !if(trim(parini%subtask_ann)=='check_symmetry_function') then
    !    rcut=subdict_ann//"rcut"
    !endif
    !if(trim(parini%approach_ann)=='atombased') then
    !    rcut=subdict_ann//"rcut"
    !endif
    rcut               =  subdict_ann//"rcut"
    ann%method         =  subdict_ann//"method"
    ann%ener_ref       =  subdict_ann//"ener_ref" 
    if(trim(parini%approach_ann)=='eem1' .or. trim(parini%approach_ann)=='cent1' .or. &
         trim(parini%approach_ann)=='centt' .or. trim(parini%approach_ann)=='cent2' &
        .or. trim(parini%approach_ann)=='cent3') then
        ann%ampl_chi       =  subdict_ann//"ampl_chi" 
        ann%prefactor_chi  =  subdict_ann//"prefactor_chi" 
        ann%gausswidth     =  subdict_ann//"gausswidth" 
        ann%hardness       =  subdict_ann//"hardness" 
        ann%chi0           =  subdict_ann//"chi0" 
        ann%qinit          =  subdict_ann//"qinit"
    endif
    if(trim(parini%approach_ann)=='centt' .or. trim(parini%approach_ann)=='cent2' .or. trim(parini%approach_ann)=='cent3') then
        ann%zion           =  subdict_ann//"zion" 
        ann%gausswidth_ion =  subdict_ann//"gausswidth_ion" 
        ann%spring_const   =  subdict_ann//"spring_const"
    endif
    !if(trim(parini%approach_ann)=='tb') then
    !    ann%ener_ref       =  subdict_ann//"ener_ref" 
    !endif
    !ann%rionic    = subdict_ann//"rionic"
    nullify(subdict_ann)
    !---------------------------------------------
    subdict_ann=>ann%dict//"symfunc"
    dict_tmp=>dict_iter(subdict_ann)
    ann%ng1 = 0 
    ann%ng2 = 0
    ann%ng3 = 0
    ann%ng4 = 0
    ann%ng5 = 0
    ann%ng6 = 0

    i0 = 0
    do while(associated(dict_tmp))
        tt=trim(dict_key(dict_tmp))
        if (tt(1:3)=="g01") ann%ng1 = ann%ng1+1
        if (tt(1:3)=="g02") ann%ng2 = ann%ng2+1
        if (tt(1:3)=="g03") ann%ng3 = ann%ng3+1
        if (tt(1:3)=="g04") ann%ng4 = ann%ng4+1
        if (tt(1:3)=="g05") ann%ng5 = ann%ng5+1
        if (tt(1:3)=="g06") ann%ng6 = ann%ng6+1
        dict_tmp=>dict_next(dict_tmp)
    end do

    ann%nn(0)=ann%ng1+ann%ng2+ann%ng3+ann%ng4+ann%ng5+ann%ng6
    call ann%ann_allocate(maxval(ann%nn(0:ann%nl)),ann%nl)
    if(iproc==0) then
        str_out_ann=''
        do i=0,ann%nl-1
            write(str_out_ann_tt,'(i5)') ann%nn(i)
            str_out_ann=trim(str_out_ann)//trim(adjustl(str_out_ann_tt))//'-'
        enddo
        write(str_out_ann_tt,'(i5)') ann%nn(i)
        str_out_ann=trim(str_out_ann)//trim(adjustl(str_out_ann_tt))
        call yaml_map('ANN architecture',trim(str_out_ann))
        !write(*,*) 'ANN architecture',trim(str_out_ann)
        !do i=0,ann%nl
        !    write(*,'(a,i1,a,i4)') 'ann%(',i,')=',ann%nn(i)
        !enddo
        !write(*,'(a,3i4)') 'n0,n1,n2 ',ann%nn(0),ann%nn(1),ann%nn(2)
        call yaml_mapping_open('symfunc numbers from input',flow=.true.)
        call yaml_map('ng1',ann%ng1)
        call yaml_map('ng2',ann%ng2)
        call yaml_map('ng3',ann%ng3)
        call yaml_map('ng4',ann%ng4)
        call yaml_map('ng5',ann%ng5)
        call yaml_map('ng6',ann%ng6)
        call yaml_mapping_close()
        !write(*,'(a,6i4)') 'ng1,ng2,ng3,ng4,ng5,ng6 ',ann%ng1,ann%ng2,ann%ng3,ann%ng4,ann%ng5,ann%ng6
    endif
    if(.not. parini%bondbased_ann .and. ann%ng1>0) then
        stop 'ERROR: symmetry function of type G3 not implemented yet.'
    endif
    if(ann%ng3>0) stop 'ERROR: symmetry function of type G3 not implemented yet.'
    if(mod(ann%ng6,3)/=0) stop 'ERROR: ng6 must be multiple of three.'
    
    count1=0; count2=0; count3=0; count4=0; count5=0; count6=0

    dict_tmp=>dict_iter(subdict_ann)
    do while(associated(dict_tmp))
        tt=trim(dict_key(dict_tmp))

        if (tt(1:3)=="g01") then
            if(.not. parini%bondbased_ann) then
                stop 'ERROR: g1 is not ready.'
            endif
            count1=count1+1
            i0=i0+1
            str1=subdict_ann//tt
            read(str1,*,iostat=ios)ann%g1eta(count1),ann%g1rs(count1),ann%gbounds(1,i0),ann%gbounds(2,i0)
            if(ios<0) then
                write(*,'(a)') 'ERROR: 5 columns are required for each of G1 symmetry functions,'
                write(*,'(a)') '       including 2 values for bounds (neglected for ANN training).'
                stop
            endif
        endif
        !method =  ann%dict//"main"//"method" 
        if (tt(1:3)=="g02") then
            count2=count2+1
            i0=i0+1
            str1=subdict_ann//tt
            if (trim(ann%method) == "behler") then
                read(str1,*,iostat=ios)ann%g2eta(count2),ann%g2rs(count2),ann%gbounds(1,i0),ann%gbounds(2,i0),sat1
                call set_radial_atomtype(parini,sat1,ann%g2i(count2))
            else
                read(str1,*,iostat=ios)ann%g2eta(count2),ann%g2rs(count2),ann%gbounds(1,i0),ann%gbounds(2,i0)
            endif
            if(ios<0) then
                write(*,*)ios
                write(*,'(a)') 'ERROR: 4 columns are required for each of G2 symmetry functions,'
                write(*,'(a)') '       including 2 values for bounds (neglected for ANN training).'
                stop
            endif
        endif

        if (tt(1:3)=="g04") then
            stop 'ERROR: g4 is not ready.'
            count4=count4+1
            i0=i0+1
            str1=subdict_ann//tt
            read(str1,*,iostat=ios) ann%g4eta(count4),ann%g4zeta(count4), &
                ann%g4lambda(count4),ann%gbounds(1,i0),ann%gbounds(2,i0)
            call set_radial_atomtype(parini,sat1,ann%g2i(count4))
            if(ios<0) then
                write(*,'(a)') 'ERROR: 6 columns are required for each of G4 symmetry functions,'
                write(*,'(a)') '       including 2 values for bounds (neglected for ANN training).'
                stop
            endif
        endif

        if (tt(1:3)=="g05") then
            count5=count5+1
            i0=i0+1
            str1=subdict_ann//tt
            if (trim(ann%method) == "behler") then
                read(str1,*,iostat=ios) ann%g5eta(count5),ann%g5zeta(count5), &
                    ann%g5lambda(count5),ann%gbounds(1,i0),ann%gbounds(2,i0),sat1,sat2
                call set_angular_atomtype(parini,sat1,sat2,ann%g5i(1,count5))
            else
                read(str1,*,iostat=ios) ann%g5eta(count5),ann%g5zeta(count5), &
                    ann%g5lambda(count5),ann%gbounds(1,i0),ann%gbounds(2,i0)
            endif
                if(ios<0) then
                    write(*,'(a)') 'ERROR: 5 columns are required for each of G5 symmetry functions,'
                write(*,'(a)') '       including 2 values for bounds (neglected for ANN training).'
                stop
            endif
        endif

        dict_tmp=>dict_next(dict_tmp)
    end do
    nullify(dict_tmp)
    nullify(subdict_ann)



!    do ig=1,ann%ng6/3
!        stop 'ERROR: g6 is not ready.'
!        i0=i0+1
!        read(ifile,'(a)') strline
!        read(strline,*,iostat=ios) ann%g6eta(ig),ann%gbounds(1,i0),ann%gbounds(2,i0)
!        if(ios<0) then
!            write(*,'(a)') 'ERROR: 4 columns are required for each of G6 symmetry functions,'
!            write(*,'(a)') '       including 2 values for bounds (neglected for ANN training).'
!            stop
!        endif
!        i0=i0+1
!        read(ifile,'(a)') strline
!        read(strline,*,iostat=ios) ann%gbounds(1,i0),ann%gbounds(2,i0)
!        if(ios<0) then
!            write(*,'(a)') 'ERROR: 2 columns are required for each of G6 symmetry functions,'
!            write(*,'(a)') '       including 2 values for bounds (neglected for ANN training).'
!            stop
!        endif
!        i0=i0+1
!        read(ifile,'(a)') strline
!        read(strline,*,iostat=ios) ann%gbounds(1,i0),ann%gbounds(2,i0)
!        if(ios<0) then
!            write(*,'(a)') 'ERROR: 2 columns are required for each of G6 symmetry functions,'
!            write(*,'(a)') '       including 2 values for bounds (neglected for ANN training).'
!            stop
!        endif
!    enddo

  !  call dict_free(ann%dict)
  !  nullify(ann%dict)
end subroutine get_symfunc_parameters_yaml
!*****************************************************************************************
subroutine write_ann_all_yaml(parini,ann_arr,iter)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(in):: ann_arr
    integer, intent(in):: iter
    !local variables
    character(21):: fn
    character(1):: fn_tt
    character(50):: filename
    integer:: i
    if(iter==-1) then
        write(fn,'(a15)') '.ann.param.yaml'
    else
        write(fn,'(a16,i5.5)') '.ann.param.yaml.',iter
    endif
    if(parini%bondbased_ann .and. trim(ann_arr%approach)=='tb') then
        if(parini%ntypat>1) then
            stop 'ERROR: writing ANN parameters for tb available only ntypat=1'
        endif
        do i=1,ann_arr%nann
            write(fn_tt,'(i1)') i
            filename=trim(parini%stypat(1))//fn_tt//trim(fn)
            call yaml_comment(trim(filename))
            !write(*,'(a)') trim(filename)
            call write_ann_yaml(parini,filename,ann_arr%ann(i),ann_arr%rcut)
        enddo
    elseif(trim(ann_arr%approach)=='atombased' .or. trim(ann_arr%approach)=='eem1' .or. &
        trim(ann_arr%approach)=='cent1' .or. trim(ann_arr%approach)=='centt' &
        .or. trim(ann_arr%approach)=='cent2' .or. trim(ann_arr%approach)=='cent3') then
        do i=1,ann_arr%nann
            filename=trim(parini%stypat(i))//trim(fn)
            !write(*,'(a)') trim(filename)
            call yaml_comment(trim(filename))
            call write_ann_yaml(parini,filename,ann_arr%ann(i),ann_arr%rcut)
        enddo
    else
        stop 'ERROR: writing ANN parameters is only for cent1,centt,cent3,tb'
    endif
end subroutine write_ann_all_yaml
!*****************************************************************************************
subroutine write_ann_yaml(parini,filename,ann,rcut)
    use futile
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann
    use dictionaries
    implicit none
    type(typ_parini), intent(in):: parini
    character(*):: filename
    type(typ_ann), intent(in):: ann
    real(8), intent(in):: rcut
    !local variables
    !integer:: 
    !real(8):: 
    integer:: i, j, k, l, ios, ialpha, i0, iline, ierr, iunit, il
    character(5):: sat1, sat2
    character(8):: key1
    character(250):: str1
    type(dictionary), pointer :: dict_ann=>null()
    type(dictionary), pointer :: subdict_ann=>null()
    dict_ann=>dict_new()
    !-------------------------------------------------------
    call set(dict_ann//"main","nodes")
    subdict_ann=>dict_ann//"main"
    !ann%nl-1 is output layer so following loop goes up to ann%nl-2
    do il=0,ann%nl-2
        call set(subdict_ann//"nodes"//il,ann%nn(il+1))
    enddo
    call set(subdict_ann//"rcut",rcut)
    call set(subdict_ann//"method",ann%method)
    call set(subdict_ann//"ener_ref",ann%ener_ref)
    if(trim(parini%approach_ann)=='eem1' .or. trim(parini%approach_ann)=='cent1' .or. &
       trim(parini%approach_ann)=='cent2' .or. &
        trim(parini%approach_ann)=='centt' .or. trim(parini%approach_ann)=='cent3') then
        call set(subdict_ann//"ampl_chi",ann%ampl_chi)
        call set(subdict_ann//"prefactor_chi",ann%prefactor_chi)
        call set(subdict_ann//"gausswidth",ann%gausswidth)
        call set(subdict_ann//"hardness",ann%hardness)
        call set(subdict_ann//"chi0",ann%chi0)
        call set(subdict_ann//"qinit",ann%qinit)
    endif
    if(trim(parini%approach_ann)=='centt' .or. trim(parini%approach_ann)=='cent2' .or. &
        trim(parini%approach_ann)=='cent3') then
        call set(subdict_ann//"zion",ann%zion)
        call set(subdict_ann//"gausswidth_ion",ann%gausswidth_ion)
        call set(subdict_ann//"spring_const",ann%spring_const)
    endif
    !if(trim(parini%approach_ann)=='tb') then
    !    !ann%ener_ref       =  subdict_ann//"ener_ref" 
    !    call set(subdict_ann//"ener_ref",ann%ener_ref)
    !endif
    nullify(subdict_ann)
    !-------------------------------------------------------
    subdict_ann=>dict_ann//"symfunc"
    i0=0
    do i=1,ann%ng1
        write(key1,'(a,i3.3)')"g01_",i 
        i0=i0+1
        write(str1,'(2f8.4,2es24.15)') ann%g1eta(i),ann%g1rs(i),ann%gbounds(1,i0),ann%gbounds(2,i0)
        call set(subdict_ann//key1,str1)
    enddo
    !-------------------------------------------------------
    do i=1,ann%ng2
        write(key1,'(a,i3.3)')"g02_",i 
        i0=i0+1
        if (trim(ann%method) == "behler") then
            sat1=parini%stypat(ann%g2i(i))
            write(str1,'(2f8.4,2es24.15,1a5)') ann%g2eta(i),ann%g2rs(i),ann%gbounds(1,i0),ann%gbounds(2,i0),trim(sat1)
        else
            write(str1,'(2f10.6,2es24.15)') ann%g2eta(i),ann%g2rs(i),ann%gbounds(1,i0),ann%gbounds(2,i0)
        endif
        call set(subdict_ann//key1,str1)
    enddo
    !-------------------------------------------------------
    do i=1,ann%ng3
        write(key1,'(a,i3.3)')"g03_",i 
        i0=i0+1
        write(str1,'(1f8.4)') ann%g3kappa(i)
        call set(subdict_ann//key1,str1)
    enddo
    do i=1,ann%ng4
        write(key1,'(a,i3.3)')"g04_",i 
        i0=i0+1
        write(str1,'(3f8.4,2es24.15)') ann%g4eta(i),ann%g4zeta(i),ann%g4lambda(i),ann%gbounds(1,i0),ann%gbounds(2,i0)
        call set(subdict_ann//key1,str1)
    enddo
    !-------------------------------------------------------
    do i=1,ann%ng5
        write(key1,'(a,i3.3)')"g05_",i 
        i0=i0+1
        if (trim(ann%method) == "behler") then
            sat1=parini%stypat(ann%g5i(1,i))
            sat2=parini%stypat(ann%g5i(2,i))
            write(str1,'(3f8.4,2es24.15,2a5)') ann%g5eta(i),ann%g5zeta(i),ann%g5lambda(i),ann%gbounds(1,i0), &
                                               ann%gbounds(2,i0),trim(sat1),trim(sat2)
        else
            write(str1,'(3f10.6,2es24.15)') ann%g5eta(i),ann%g5zeta(i),ann%g5lambda(i),ann%gbounds(1,i0),ann%gbounds(2,i0)
        endif
        call set(subdict_ann//key1,str1)
    enddo
    nullify(subdict_ann)
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
    subdict_ann=>dict_ann//"weights"
    i0=0
    do ialpha=1,ann%nl
        write(str1,'(2(a,i1))') '#main nodes weights connecting layers ',ialpha,' and ',ialpha-1
        call set(subdict_ann//i0//"comment",str1)

        do j=1,ann%nn(ialpha)
            do i=1,ann%nn(ialpha-1)
                write(str1,'(es24.15)') ann%a(i,j,ialpha)
                call set(subdict_ann//i0,str1)
                i0=i0+1
            enddo
        enddo
        write(str1,'(a,i1)') '#bias nodes weights for layer ',ialpha
        call set(subdict_ann//i0//"comment",str1)
        do i=1,ann%nn(ialpha)
            write(str1,'(es24.15)') ann%b(i,ialpha)
            call set(subdict_ann//i0,str1)
            i0=i0+1
        enddo
    enddo
    nullify(subdict_ann)
    !-------------------------------------------------------
    iunit=f_get_free_unit(10**5)
    call yaml_set_stream(unit=iunit,filename=trim(filename),&
         record_length=92,istat=ierr,setdefault=.false.,tabbing=0,position='rewind')
    if (ierr==0) then
        call yaml_dict_dump(dict_ann,unit=iunit)
       call yaml_close_stream(unit=iunit)
    else
       call yaml_warning('Failed to create'//trim(filename)//', error code='//trim(yaml_toa(ierr)))
    end if
    !-------------------------------------------------------
    call dict_free(dict_ann)
    nullify(dict_ann)
end subroutine write_ann_yaml
!*****************************************************************************************
subroutine read_ann_yaml(parini,ann_arr)
    use futile
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_processors, only: iproc
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    !local variables
    !integer:: 
    integer:: i, j, k, l, ios, i0, i1, ifile, ialpha, iann
    !character(100):: ttstr
    real(8):: bound_l, bound_u, rcut
    character(16):: fn
    character(1):: fn_tt
    character(50):: filename
    do iann=1,ann_arr%nann
        if (parini%restart_param) then
            write(fn,'(a15)') '.ann.input.yaml'
        else
            write(fn,'(a15)') '.ann.param.yaml'
        endif
        if(parini%bondbased_ann .and. trim(ann_arr%approach)=='tb') then
            if(parini%ntypat>1) then
                stop 'ERROR: writing ANN parameters for tb available only ntypat=1'
            endif
            write(fn_tt,'(i1)') iann
            filename=trim(parini%stypat(1))//fn_tt//trim(fn)
            write(*,'(a)') trim(filename)
        elseif(trim(ann_arr%approach)=='eem1' .or. trim(ann_arr%approach)=='cent1' .or. &
            trim(ann_arr%approach)=='centt' .or. trim(ann_arr%approach)=='cent3' .or. &
            trim(ann_arr%approach)=='cent2' .or. trim(ann_arr%approach)=='atombased') then
            filename=trim(parini%stypat(iann))//trim(fn)
        else
            stop 'ERROR: reading ANN parameters is only for cent1,centt,cent3,tb'
        endif
        !-------------------------------------------------------
        call get_symfunc_parameters_yaml(parini,iproc,filename,ann_arr%ann(iann),rcut)
        ann_arr%rcut = rcut
        do i0=1,ann_arr%ann(iann)%nn(0)
            bound_l=ann_arr%ann(iann)%gbounds(1,i0)
            bound_u=ann_arr%ann(iann)%gbounds(2,i0)
            ann_arr%ann(iann)%two_over_gdiff(i0)=2.d0/(bound_u-bound_l)
        enddo
        !-------------------------------------------------------
        i1=0
        do ialpha=1,ann_arr%ann(iann)%nl
            do j=1,ann_arr%ann(iann)%nn(ialpha)
                do i=1,ann_arr%ann(iann)%nn(ialpha-1)
                    ann_arr%ann(iann)%a(i,j,ialpha)=ann_arr%ann(iann)%dict//"weights"//i1
                    i1=i1+1
                enddo
            enddo
            do i=1,ann_arr%ann(iann)%nn(ialpha)
                ann_arr%ann(iann)%b(i,ialpha)=ann_arr%ann(iann)%dict//"weights"//i1
                i1=i1+1
            enddo
        enddo
        !-------------------------------------------------------
        call dict_free(ann_arr%ann(iann)%dict)
        nullify(ann_arr%ann(iann)%dict)
    enddo
end subroutine read_ann_yaml
!*****************************************************************************************
subroutine set_dict_ann(ann,fname,stypat)
    use dictionaries
    use yaml_parse
    use dynamic_memory
    use mod_ann, only: typ_ann
    implicit none
    !local variales
    type(dictionary), pointer :: dict=>null()
    type(typ_ann), intent(inout):: ann
    character, dimension(:), allocatable :: fbuf
    character(5):: stypat
    character(len=*):: fname 
    integer(kind = 8) :: cbuf, cbuf_len

    call getFileContent(cbuf,cbuf_len,fname,len_trim(fname))
    fbuf=f_malloc0_str(1,int(cbuf_len),id='fbuf')
    call copyCBuffer(fbuf,cbuf,cbuf_len)
    call freeCBuffer(cbuf)
    !then parse the user's input file
    call yaml_parse_from_char_array(dict,fbuf)
    call f_free_str(1,fbuf)
    call dict_copy(ann%dict,dict//0)

    call dict_free(dict)
    nullify(dict)
end subroutine set_dict_ann
!*****************************************************************************************
subroutine read_data_yaml(parini,filename_list,atoms_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr, atom_allocate_old, atom_deallocate, atom_copy_old
    use mod_atoms, only: atom_deallocate_old, set_rat_atoms
    use mod_yaml_conf, only: read_yaml_conf
    use mod_bin, only: read_bin_conf
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: filename_list
    type(typ_atoms_arr), intent(inout):: atoms_arr
    !local variables
    integer:: i, iat, ios, k, iconf
    character(256):: filename, fn_tmp, fn_ext
    type(typ_atoms_arr):: atoms_arr_of !configuration of one file
    type(typ_atoms_arr):: atoms_arr_t
    real(8):: fx, fy, fz
    integer:: nconfmax, ind, len_filename, nfiles, nfiles_max, ifile
    character(256), allocatable:: fn_list(:)
    call f_routine(id='read_data_yaml')
    nconfmax=1*10**5
    nfiles_max=5*10**4
    allocate(atoms_arr_t%atoms(nconfmax))
    allocate(atoms_arr_t%fn(nconfmax))
    allocate(atoms_arr_t%lconf(nconfmax))
    atoms_arr_t%nconf=0
    allocate(fn_list(nfiles_max))
    call read_list_files_yaml(filename_list,nfiles_max,fn_list,nfiles)
    over_files: do ifile=1,nfiles
        filename=trim(fn_list(ifile))
        fn_tmp=adjustl(trim(filename))
        if(fn_tmp(1:1)=='#') cycle
        !call acf_read_new(parini,filename,10000,atoms_arr_of)
        ind=index(filename,'.',back=.true.)
        len_filename=len(filename)
        fn_ext=filename(ind+1:len_filename)
        if(trim(fn_ext)=='bin') then
            call read_bin_conf(parini,filename,atoms_arr_of)
        elseif(trim(fn_ext)=='yaml') then
            call read_yaml_conf(parini,filename,10000,atoms_arr_of)
        else
            write(*,*) 'ERROR: only binary and yaml files can be read in read_data_yaml'
            stop
        endif
        !if(parini%read_forces_ann) then
        !    open(unit=2,file=trim(filename_force),status='old',iostat=ios)
        !endif
        !if(ios/=0) then;write(*,'(a)') 'ERROR: failure openning force of list_posinp';stop;endif
        over_iconf: do iconf=1,atoms_arr_of%nconf
            atoms_arr_t%nconf=atoms_arr_t%nconf+1
            if(atoms_arr_t%nconf>nconfmax) then
                stop 'ERROR: too many configurations, change parameter nconfmax.'
            endif
            if(atoms_arr_of%atoms(iconf)%ntrial>0) then
                call atom_allocate_old(atoms_arr_t%atoms(atoms_arr_t%nconf),atoms_arr_of%atoms(iconf)%nat,0,0, &
                    ntrial=atoms_arr_of%atoms(iconf)%ntrial)
            else
                call atom_allocate_old(atoms_arr_t%atoms(atoms_arr_t%nconf),atoms_arr_of%atoms(iconf)%nat,0,0)
            endif
            atoms_arr_t%atoms(atoms_arr_t%nconf)%epot=atoms_arr_of%atoms(iconf)%epot
            atoms_arr_t%atoms(atoms_arr_t%nconf)%qtot=atoms_arr_of%atoms(iconf)%qtot
            atoms_arr_t%atoms(atoms_arr_t%nconf)%dpm(1:3)=atoms_arr_of%atoms(iconf)%dpm(1:3)
            atoms_arr_t%atoms(atoms_arr_t%nconf)%elecfield(1:3)=atoms_arr_of%atoms(iconf)%elecfield(1:3)
            atoms_arr_t%atoms(atoms_arr_t%nconf)%boundcond=trim(atoms_arr_of%atoms(iconf)%boundcond)
            atoms_arr_t%atoms(atoms_arr_t%nconf)%cellvec(1:3,1:3)=atoms_arr_of%atoms(iconf)%cellvec(1:3,1:3)
            !if(parini%read_forces_ann) read(2,*)
            do iat=1,atoms_arr_of%atoms(iconf)%ntrial
                atoms_arr_t%atoms(atoms_arr_t%nconf)%trial_ref_energy(iat)=atoms_arr_of%atoms(iconf)%trial_ref_energy(iat)
                atoms_arr_t%atoms(atoms_arr_t%nconf)%trial_ref_nat(iat)=atoms_arr_of%atoms(iconf)%trial_ref_nat(iat)
                atoms_arr_t%atoms(atoms_arr_t%nconf)%trial_ref_disp(1,iat)=atoms_arr_of%atoms(iconf)%trial_ref_disp(1,iat)
                atoms_arr_t%atoms(atoms_arr_t%nconf)%trial_ref_disp(2,iat)=atoms_arr_of%atoms(iconf)%trial_ref_disp(2,iat)
                atoms_arr_t%atoms(atoms_arr_t%nconf)%trial_ref_disp(3,iat)=atoms_arr_of%atoms(iconf)%trial_ref_disp(3,iat)
            end do
            call set_rat_atoms(atoms_arr_t%atoms(atoms_arr_t%nconf),atoms_arr_of%atoms(iconf),setall=.true.)
            do iat=1,atoms_arr_of%atoms(iconf)%nat
                atoms_arr_t%atoms(atoms_arr_t%nconf)%sat(iat)=atoms_arr_of%atoms(iconf)%sat(iat)
                if(parini%read_forces_ann) then
                    !read(2,*) fx,fy,fz
                    atoms_arr_t%atoms(atoms_arr_t%nconf)%fat(1,iat)=atoms_arr_of%atoms(iconf)%fat(1,iat)
                    atoms_arr_t%atoms(atoms_arr_t%nconf)%fat(2,iat)=atoms_arr_of%atoms(iconf)%fat(2,iat)
                    atoms_arr_t%atoms(atoms_arr_t%nconf)%fat(3,iat)=atoms_arr_of%atoms(iconf)%fat(3,iat)
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
        !if(parini%read_forces_ann) close(2)
    enddo over_files
    !close(1)
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
    deallocate(fn_list)
    call f_release_routine()
end subroutine read_data_yaml
!*****************************************************************************************
