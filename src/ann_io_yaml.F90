!*****************************************************************************************
subroutine read_input_ann_yaml(parini,iproc,ann_arr)
    use mod_interface
    use futile
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
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
    do iann=1,ann_arr%n
        if(parini%bondbased_ann) then
            stypat=parini%stypat(1)
        else
            stypat=parini%stypat(iann)
        endif
        !-------------------------------------------------------
        fname = trim(stypat)//'.ann.input.yaml'
        write(*,*)trim(fname)
        call get_symfunc_parameters_yaml(parini,iproc,fname,ann_arr%ann(iann),rcut)
        ann_arr%rcut = rcut
        !-------------------------------------------------------

    enddo
    !call f_lib_finalize()
end subroutine read_input_ann_yaml
!*****************************************************************************************
subroutine get_symfunc_parameters_yaml(parini,iproc,fname,ann,rcut)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann
    use dictionaries
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    type(dictionary), pointer :: subdict_ann
    type(dictionary), pointer :: dict_tmp
    type(typ_ann), intent(inout):: ann
    integer, intent(in):: iproc
    !local variables
    integer:: ig, ios, i0, i, il 
    character(50):: fname, method, sat1, sat2
    character(250):: tt, str1
    integer :: count1, count2, count3, count4, count5, count6
    character(5):: stypat
    real(8)::rcut
    call set_dict_ann(ann,fname,stypat)
    call yaml_comment('USER INPUT FILE',hfill='~')
    call yaml_dict_dump(ann%dict_ann)
    call yaml_comment('',hfill='~')
    !logical:: all_read

    ann%nl=dict_len(ann%dict_ann//"main"//"nodes")
    do il=0,ann%nl-1
        ann%nn(il+1)=ann%dict_ann//"main"//"nodes"//il
    enddo
    ann%nl=ann%nl+1 !adding the output layer to total number of layers
    ann%nn(ann%nl)=1 !setting the output layer
    if(trim(parini%subtask_ann)=='check_symmetry_function') then
        rcut=ann%dict_ann//"main"//"rcut"
    endif
    if(trim(parini%approach_ann)/='atombased' .and. trim(parini%approach_ann)/='tb' ) then
        ann%ampl_chi       =  ann%dict_ann//"main"//"ampl_chi" 
        ann%prefactor_chi  =  ann%dict_ann//"main"//"prefactor_chi" 
        ann%zion           =  ann%dict_ann//"main"//"zion" 
        ann%gausswidth_ion =  ann%dict_ann//"main"//"gausswidth_ion" 
        ann%ener_ref       =  ann%dict_ann//"main"//"ener_ref" 
        ann%gausswidth     =  ann%dict_ann//"main"//"gausswidth" 
        ann%hardness       =  ann%dict_ann//"main"//"hardness" 
        ann%chi0           =  ann%dict_ann//"main"//"chi0" 
        ann%spring_const   =  ann%dict_ann//"main"//"spring_const"
        ann%qinit          =  ann%dict_ann//"main"//"qinit"
        rcut               =  ann%dict_ann//"main"//"rcut"
    elseif(trim(parini%approach_ann)=='tb') then
        ann%ener_ref       =  ann%dict_ann//"main"//"ener_ref" 
    endif
    !ann%rionic    = ann%dict_ann//"main"//"rionic"
    !---------------------------------------------
    dict_tmp=>dict_iter(ann%dict_ann//"symfunc")
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
    if(iproc==0) then
        do i=0,ann%nl
            write(*,'(a,i1,a,i4)') 'ann%(',i,')=',ann%nn(i)
        enddo
        !write(*,'(a,3i4)') 'n0,n1,n2 ',ann%nn(0),ann%nn(1),ann%nn(2)
        write(*,'(a,6i4)') 'ng1,ng2,ng3,ng4,ng5,ng6 ',ann%ng1,ann%ng2,ann%ng3,ann%ng4,ann%ng5,ann%ng6
    endif
    if(.not. parini%bondbased_ann .and. ann%ng1>0) then
        stop 'ERROR: symmetry function of type G3 not implemented yet.'
    endif
    if(ann%ng3>0) stop 'ERROR: symmetry function of type G3 not implemented yet.'
    if(mod(ann%ng6,3)/=0) stop 'ERROR: ng6 must be multiple of three.'
    
    count1=0; count2=0; count3=0; count4=0; count5=0; count6=0

    dict_tmp=>dict_iter(ann%dict_ann//"symfunc")
    do while(associated(dict_tmp))
        tt=trim(dict_key(dict_tmp))

        if (tt(1:3)=="g01") then
            if(.not. parini%bondbased_ann) then
                stop 'ERROR: g1 is not ready.'
            endif
            count1=count1+1
            i0=i0+1
            str1=ann%dict_ann//"symfunc"//tt
            read(str1,*,iostat=ios)ann%g1eta(count1),ann%g1rs(count1),ann%gbounds(1,i0),ann%gbounds(2,i0)
            if(ios<0) then
                write(*,'(a)') 'ERROR: 5 columns are required for each of G1 symmetry functions,'
                write(*,'(a)') '       including 2 values for bounds (neglected for ANN training).'
                stop
            endif
        endif
        method =  ann%dict_ann//"main"//"method" 
        if (tt(1:3)=="g02") then
            count2=count2+1
            i0=i0+1
            str1=ann%dict_ann//"symfunc"//tt
            if (trim(method) == "behler") then
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
            str1=ann%dict_ann//"symfunc"//tt
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
            str1=ann%dict_ann//"symfunc"//tt
            if (trim(method) == "behler") then
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

  !  call dict_free(ann%dict_ann)
  !  nullify(ann%dict_ann)
end subroutine get_symfunc_parameters_yaml
!*****************************************************************************************
subroutine write_ann_all_yaml(parini,ann_arr,iter)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
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
        do i=1,ann_arr%n
            write(fn_tt,'(i1)') i
            filename=trim(parini%stypat(1))//fn_tt//trim(fn)
            write(*,'(a)') trim(filename)
            call write_ann_yaml(parini,filename,ann_arr%ann(i))
        enddo
    elseif(trim(ann_arr%approach)=='eem1' .or. trim(ann_arr%approach)=='cent2') then
        do i=1,ann_arr%n
            filename=trim(parini%stypat(i))//trim(fn)
            write(*,'(a)') trim(filename)
            call write_ann_yaml(parini,filename,ann_arr%ann(i))
        enddo
    else
        stop 'ERROR: writing ANN parameters is only for eem1,cent2,tb'
    endif
end subroutine write_ann_all_yaml
!*****************************************************************************************
subroutine write_ann_yaml(parini,filename,ann)
    use mod_interface
    use futile
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann
    implicit none
    type(typ_parini), intent(in):: parini
    character(*):: filename
    type(typ_ann), intent(in):: ann
    !local variables
    !integer:: 
    !real(8):: 
    integer:: i, j, k, l, ios, ialpha, i0, iline, ierr, iunit
    character(5):: sat1, sat2
    character(8):: key1
    character(250):: str1
    character(50)::  method
    method =  ann%dict_ann//"main"//"method"
    i0=0
    do i=1,ann%ng1
        write(key1,'(a,i3.3)')"g01_",i 
        i0=i0+1
        write(str1,'(2f8.4,2es24.15)') ann%g1eta(i),ann%g1rs(i),ann%gbounds(1,i0),ann%gbounds(2,i0)
        call set(ann%dict_ann//"symfunc"//key1,str1)
    enddo
    !-------------------------------------------------------
    do i=1,ann%ng2
        write(key1,'(a,i3.3)')"g02_",i 
        i0=i0+1
        if (trim(method) == "behler") then
            sat1=parini%stypat(ann%g2i(i))
            write(str1,'(2f8.4,2es24.15,1a5)') ann%g2eta(i),ann%g2rs(i),ann%gbounds(1,i0),ann%gbounds(2,i0),trim(sat1)
        else
            write(str1,'(2f8.4,2es24.15)') ann%g2eta(i),ann%g2rs(i),ann%gbounds(1,i0),ann%gbounds(2,i0)
        endif
        call set(ann%dict_ann//"symfunc"//key1,str1)
    enddo
    !-------------------------------------------------------
    do i=1,ann%ng3
        write(key1,'(a,i3.3)')"g03_",i 
        i0=i0+1
        write(str1,'(1f8.4)') ann%g3kappa(i)
        call set(ann%dict_ann//"symfunc"//key1,str1)
    enddo
    do i=1,ann%ng4
        write(key1,'(a,i3.3)')"g04_",i 
        i0=i0+1
        write(str1,'(3f8.4,2es24.15)') ann%g4eta(i),ann%g4zeta(i),ann%g4lambda(i),ann%gbounds(1,i0),ann%gbounds(2,i0)
        call set(ann%dict_ann//"symfunc"//key1,str1)
    enddo
    !-------------------------------------------------------
    do i=1,ann%ng5
        write(key1,'(a,i3.3)')"g05_",i 
        i0=i0+1
        if (trim(method) == "behler") then
            sat1=parini%stypat(ann%g5i(1,i))
            sat2=parini%stypat(ann%g5i(2,i))
            write(str1,'(3f8.4,2es24.15,2a5)') ann%g5eta(i),ann%g5zeta(i),ann%g5lambda(i),ann%gbounds(1,i0), &
                                               ann%gbounds(2,i0),trim(sat1),trim(sat2)
        else
            write(str1,'(3f8.4,2es24.15)') ann%g5eta(i),ann%g5zeta(i),ann%g5lambda(i),ann%gbounds(1,i0),ann%gbounds(2,i0)
        endif
        call set(ann%dict_ann//"symfunc"//key1,str1)
    enddo
!    !-------------------------------------------------------
!    write(1,'(i6,2x,a)') ann%ng6,'#ng6'
!    do i=1,ann%ng6/3
!        i0=i0+1
!        write(1,'(1f8.4,2es24.15)') ann%g6eta(i),ann%gbounds(1,i0),ann%gbounds(2,i0)
!        i0=i0+1
!        write(1,'(16x,2es24.15)') ann%gbounds(1,i0),ann%gbounds(2,i0)
!        i0=i0+1
!        write(1,'(16x,2es24.15)') ann%gbounds(1,i0),ann%gbounds(2,i0)
!    enddo
    !-------------------------------------------------------
    i0=0
    do ialpha=1,ann%nl
        write(str1,'(2(a,i1))') '#main nodes weights connecting layers ',ialpha,' and ',ialpha-1
        call set(ann%dict_ann//"weights"//i0//"comment",str1)

        do j=1,ann%nn(ialpha)
            do i=1,ann%nn(ialpha-1)
                write(str1,'(es24.15)') ann%a(i,j,ialpha)
                call set(ann%dict_ann//"weights"//i0,str1)
                i0=i0+1
            enddo
        enddo
        write(str1,'(a,i1)') '#bias nodes weights for layer ',ialpha
        call set(ann%dict_ann//"weights"//i0//"comment",str1)
        do i=1,ann%nn(ialpha)
            write(str1,'(es24.15)') ann%b(i,ialpha)
            call set(ann%dict_ann//"weights"//i0,str1)
            i0=i0+1
        enddo
    enddo
    iunit=f_get_free_unit(10**5)
    call yaml_set_stream(unit=iunit,filename=trim(filename),&
         record_length=92,istat=ierr,setdefault=.false.,tabbing=0,position='rewind')
    if (ierr==0) then
        call yaml_dict_dump(ann%dict_ann,unit=iunit)
       call yaml_close_stream(unit=iunit)
    else
       call yaml_warning('Failed to create'//trim(filename)//', error code='//trim(yaml_toa(ierr)))
    end if
    !-------------------------------------------------------
end subroutine write_ann_yaml
!*****************************************************************************************
subroutine read_ann_yaml(parini,ann_arr)
    use mod_interface
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
    do iann=1,ann_arr%n
        write(fn,'(a15)') '.ann.param.yaml'
        if(parini%bondbased_ann .and. trim(ann_arr%approach)=='tb') then
            if(parini%ntypat>1) then
                stop 'ERROR: writing ANN parameters for tb available only ntypat=1'
            endif
            write(fn_tt,'(i1)') iann
            filename=trim(parini%stypat(1))//fn_tt//trim(fn)
            write(*,'(a)') trim(filename)
        elseif(trim(ann_arr%approach)=='eem1' .or. trim(ann_arr%approach)=='cent2') then
            filename=trim(parini%stypat(iann))//trim(fn)
        else
            stop 'ERROR: reading ANN parameters is only for eem1,cent2,tb'
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
                    ann_arr%ann(iann)%a(i,j,ialpha)=ann_arr%ann(iann)%dict_ann//"weights"//i1
                    i1=i1+1
                enddo
            enddo
            do i=1,ann_arr%ann(iann)%nn(ialpha)
                ann_arr%ann(iann)%b(i,ialpha)=ann_arr%ann(iann)%dict_ann//"weights"//i1
                i1=i1+1
            enddo
        enddo
        !-------------------------------------------------------
    enddo
end subroutine read_ann_yaml
!*****************************************************************************************
subroutine set_dict_ann(ann,fname,stypat)
    use mod_interface
    use dictionaries
    use yaml_parse
    use dynamic_memory
    use mod_ann, only: typ_ann
    implicit none
    !local variales
    type(dictionary), pointer :: dict
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
    call dict_copy(ann%dict_ann,dict//0)

    call dict_free(dict)
    nullify(dict)
end subroutine set_dict_ann
!*****************************************************************************************
