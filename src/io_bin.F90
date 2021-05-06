!*****************************************************************************************
module mod_bin
    implicit none
    private
    public:: read_bin_conf, write_bin_conf
contains
!*****************************************************************************************
subroutine read_bin_conf(parini,filename,atoms_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr, atom_allocate
    use dynamic_memory
    use futile
    implicit none
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: filename
    type(typ_atoms_arr), intent(inout):: atoms_arr
    !local variables
    integer:: ios, iconf, i, nconf, iat, k, nat, ii, ndigit, iiconf
    real(8):: ver
    integer:: iver, iunit
    iunit=f_get_free_unit(10**5)
    open(unit=iunit,file=trim(filename),status='old',form='unformatted', &
        access='stream',iostat=ios)
    if(ios/=0) then
        write(*,'(2a)') 'ERROR: failure openning file ',trim(filename)
        stop
    endif
    read(iunit,iostat=ios) ver
    if(ios/=0) then
        write(*,'(2a)') 'ERROR: failure reading to file ',trim(filename)
        stop
    endif
    iver=int(ver)
    if(iver==1) then
        call read_bin_conf_v1(parini,filename,iunit,atoms_arr)
    else
        write(*,*) 'ERROR: unknown version of binary file, ver=',iver
        stop
    endif
    close(iunit)
    call yaml_mapping_open('Number of configurations read',flow=.true.)
    call yaml_map('filename',trim(filename))
    call yaml_map('nconf',atoms_arr%nconf)
    call yaml_mapping_close()
end subroutine read_bin_conf
!*****************************************************************************************
subroutine read_bin_conf_v1(parini,filename,iunit,atoms_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr, iatom_to_sat, atom_allocate, update_rat
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: filename
    integer, intent(in):: iunit
    type(typ_atoms_arr), intent(inout):: atoms_arr
    !local variables
    integer:: ios, iconf, i, nconf, iat, k, ii, ndigit, iiconf
    integer:: nwa, nread, isat, nat, irr, iwa, ipos, ibc
    real(8), allocatable:: wa(:)
    real(8):: r_nconf, rr
    logical:: cell_present, epot_present, fat_present
    logical:: qtot_present, bemoved_present, vat_present
    read(iunit,iostat=ios) r_nconf
    nconf=int(r_nconf)
    if(ios/=0) then
        write(*,'(2a)') 'ERROR: failure reading to file ',trim(filename)
        stop
    endif
    nwa=10**8
    wa=f_malloc([1.to.nwa],id='wa_read_bin_conf_v1')
    read(iunit,iostat=ios) wa
    !if(ios/=0) then
    !    inquire(unit=iunit,pos=ipos) !can be used in future
    !    nwa=ipos-ipos_ref
    !    exit
    !endif
    !ipos_ref=ipos_ref+nwa
    if(ios==0) then
        write(*,*) 'ERROR: read_bin_conf_v1 cannot read such a large file.',trim(filename)
        stop
    endif
    inquire(unit=iunit,pos=ipos)
    if(mod(ipos-1,8)/=0) then
        write(*,*) 'ERROR: read_bin_conf_v1: something wrong in file= ',trim(filename)
        stop
    endif
    nread=(ipos-1)/8-2
    atoms_arr%nconf=nconf
    allocate(atoms_arr%atoms(nconf))
    iwa=0
    do iconf=1,nconf
        cell_present=.false.
        epot_present=.false.
        fat_present=.false.
        qtot_present=.false.
        bemoved_present=.false.
        vat_present=.false.
        iwa=iwa+1 ; if(wa(iwa)==1) cell_present=.true.
        iwa=iwa+1 ; if(wa(iwa)==1) epot_present=.true.
        iwa=iwa+1 ; if(wa(iwa)==1) fat_present=.true.
        iwa=iwa+1 ; if(wa(iwa)==1) qtot_present=.true.
        iwa=iwa+1 ; if(wa(iwa)==1) bemoved_present=.true.
        iwa=iwa+1 ; if(wa(iwa)==1) vat_present=.true.
        iwa=iwa+1 ; ibc=int(wa(iwa))
        if(ibc==0) then
            atoms_arr%atoms(iconf)%boundcond='free'
        elseif(ibc==1) then
            atoms_arr%atoms(iconf)%boundcond='wire'
        elseif(ibc==2) then
            atoms_arr%atoms(iconf)%boundcond='slab'
        elseif(ibc==3) then
            atoms_arr%atoms(iconf)%boundcond='bulk'
        else
            write(*,*) "ERROR: invalid coding for boundcond in read_bin_conf_v1, ibc=",ibc
            stop
        endif
        iwa=iwa+1 ; nat=int(wa(iwa))
        call atom_allocate(atoms_arr%atoms(iconf),nat,0,0)
        do iat=1,atoms_arr%atoms(iconf)%nat
            iwa=iwa+1 ; atoms_arr%atoms(iconf)%ratp(1,iat)=wa(iwa)
            iwa=iwa+1 ; atoms_arr%atoms(iconf)%ratp(2,iat)=wa(iwa)
            iwa=iwa+1 ; atoms_arr%atoms(iconf)%ratp(3,iat)=wa(iwa)
        enddo
        call update_rat(atoms_arr%atoms(iconf),upall=.true.)
        do iat=1,atoms_arr%atoms(iconf)%nat
            iwa=iwa+1 ; isat=int(wa(iwa))
            call iatom_to_sat(isat,atoms_arr%atoms(iconf)%sat(iat))
        enddo
        if(cell_present) then
            iwa=iwa+1 ; atoms_arr%atoms(iconf)%cellvec(1,1)=wa(iwa)
            iwa=iwa+1 ; atoms_arr%atoms(iconf)%cellvec(2,1)=wa(iwa)
            iwa=iwa+1 ; atoms_arr%atoms(iconf)%cellvec(3,1)=wa(iwa)
            iwa=iwa+1 ; atoms_arr%atoms(iconf)%cellvec(1,2)=wa(iwa)
            iwa=iwa+1 ; atoms_arr%atoms(iconf)%cellvec(2,2)=wa(iwa)
            iwa=iwa+1 ; atoms_arr%atoms(iconf)%cellvec(3,2)=wa(iwa)
            iwa=iwa+1 ; atoms_arr%atoms(iconf)%cellvec(1,3)=wa(iwa)
            iwa=iwa+1 ; atoms_arr%atoms(iconf)%cellvec(2,3)=wa(iwa)
            iwa=iwa+1 ; atoms_arr%atoms(iconf)%cellvec(3,3)=wa(iwa)
        endif
        if(epot_present) then
            iwa=iwa+1 ; atoms_arr%atoms(iconf)%epot=wa(iwa)
        else
            atoms_arr%atoms(iconf)%epot=1.d100
        endif
        if(fat_present) then
            do iat=1,atoms_arr%atoms(iconf)%nat
                iwa=iwa+1 ; atoms_arr%atoms(iconf)%fat(1,iat)=wa(iwa)
                iwa=iwa+1 ; atoms_arr%atoms(iconf)%fat(2,iat)=wa(iwa)
                iwa=iwa+1 ; atoms_arr%atoms(iconf)%fat(3,iat)=wa(iwa)
            enddo
        endif
        if(qtot_present) then
            iwa=iwa+1 ; atoms_arr%atoms(iconf)%qtot=wa(iwa)
        else
            atoms_arr%atoms(iconf)%qtot=0.d0
        endif
        atoms_arr%atoms(iconf)%bemoved=.true.
        if(bemoved_present) then
            do iat=1,atoms_arr%atoms(iconf)%nat
                iwa=iwa+1 ; rr=wa(iwa)
                irr=int(rr)
                select case(irr)
                case(0)
                    atoms_arr%atoms(iconf)%bemoved(1,iat)=.false.
                    atoms_arr%atoms(iconf)%bemoved(2,iat)=.false.
                    atoms_arr%atoms(iconf)%bemoved(3,iat)=.false.
                case(1)
                    atoms_arr%atoms(iconf)%bemoved(2,iat)=.false.
                    atoms_arr%atoms(iconf)%bemoved(3,iat)=.false.
                case(2)
                    atoms_arr%atoms(iconf)%bemoved(1,iat)=.false.
                    atoms_arr%atoms(iconf)%bemoved(3,iat)=.false.
                case(3)
                    atoms_arr%atoms(iconf)%bemoved(3,iat)=.false.
                case(4)
                    atoms_arr%atoms(iconf)%bemoved(1,iat)=.false.
                    atoms_arr%atoms(iconf)%bemoved(2,iat)=.false.
                case(5)
                    atoms_arr%atoms(iconf)%bemoved(2,iat)=.false.
                case(6)
                    atoms_arr%atoms(iconf)%bemoved(1,iat)=.false.
                case(7)
                    !already done
                case default
                    stop 'ERROR: wrong value for bemoved in read_bin_conf_v1'
                end select
            enddo
        endif
        if(vat_present) then
            do iat=1,atoms_arr%atoms(iconf)%nat
                iwa=iwa+1 ; atoms_arr%atoms(iconf)%vat(1,iat)=wa(iwa)
                iwa=iwa+1 ; atoms_arr%atoms(iconf)%vat(2,iat)=wa(iwa)
                iwa=iwa+1 ; atoms_arr%atoms(iconf)%vat(3,iat)=wa(iwa)
            enddo
        endif
    enddo
    if(iwa/=nread) then
        write(*,*) 'ERROR: inconsistent values in file ',trim(filename)
    endif
    !write(*,*) 'nread,iwa',nread,iwa
    call f_free(wa)
end subroutine read_bin_conf_v1
!*****************************************************************************************
subroutine write_bin_conf(file_info,atoms,strkey)
    use mod_atoms, only: typ_atoms, typ_file_info
    use dynamic_memory
    use futile
    implicit none
    type(typ_file_info), intent(inout):: file_info
    type(typ_atoms), intent(in):: atoms
    character(*), optional, intent(in):: strkey
    !local variables
    integer:: ios, iconf, i, nconf, iat, k, nat, ii, ndigit, iiconf
    real(8):: ver
    integer:: iver, iunit
    iunit=f_get_free_unit(10**5)
    if(trim(file_info%file_position)=='new') then
        open(unit=iunit,file=trim(file_info%filename_positions),status='replace',form='unformatted', &
            access='stream',action='write',iostat=ios)
    elseif(trim(file_info%file_position)=='append') then
        open(unit=iunit,file=trim(file_info%filename_positions),status='old',position='append',form='unformatted', &
            access='stream',action='readwrite',iostat=ios)
    endif
    if(ios/=0) then
        write(*,'(2a)') 'ERROR: failure openning file ',trim(file_info%filename_positions)
        stop
    endif
    call write_bin_conf_v1(file_info%filename_positions,file_info%file_position,iunit,atoms)
    close(iunit)
end subroutine write_bin_conf
!*****************************************************************************************
subroutine write_bin_conf_v1(filename,file_position,iunit,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, sat_to_iatom, get_rat
    use dynamic_memory
    implicit none
    character(*), intent(in):: filename
    character(*), intent(in):: file_position
    integer, intent(in):: iunit
    type(typ_atoms), intent(in):: atoms
    !local variables
    integer:: ios, iconf, i, nconf, iat, k, ii, ndigit, iiconf
    integer:: nwa, nread, isat, irr, iwa, ipos, ibc, ibm, i1, i2, i3
    real(8), allocatable:: wa(:)
    real(8):: r_nconf, rr, ver
    logical:: cell_present, epot_present, fat_present
    logical:: qtot_present, bemoved_present, vat_present
    real(8), allocatable:: rat(:,:)
    allocate(rat(3,atoms%nat))
    nwa=10**8
    wa=f_malloc([1.to.nwa],id='wa_write_bin_conf_v1')
    iwa=0
    if(file_position=='new') then
        ver=1.d0
        iwa=iwa+1 ; wa(iwa)=ver
        ipos=1
        nconf=1
        iwa=iwa+1 ; wa(iwa)=real(nconf,8)
    elseif(file_position=='append') then
        inquire(unit=iunit,pos=ipos)
        read(iunit,iostat=ios,pos=9) r_nconf
        if(ios/=0) then
            write(*,*) 'ERROR: cannot read the number of already written confs in the file= ',trim(filename)
            stop
        endif
        nconf=int(r_nconf)
        nconf=nconf+1
        write(iunit,iostat=ios,pos=9) real(nconf,8)
    else
        write(*,*) 'ERROR: unknown file_position in write_bin_conf_v1, file_position',trim(file_position)
    endif

    cell_present=.true.
    epot_present=.true.
    fat_present=.true.
    qtot_present=.true.
    bemoved_present=.true.
    vat_present=.false.
    if(cell_present) then
        iwa=iwa+1 ; wa(iwa)=1.d0
    else
        iwa=iwa+1 ; wa(iwa)=0.d0
    endif
    if(epot_present) then
        iwa=iwa+1 ; wa(iwa)=1.d0
    else
        iwa=iwa+1 ; wa(iwa)=0.d0
    endif
    if(fat_present) then
        iwa=iwa+1 ; wa(iwa)=1.d0
    else
        iwa=iwa+1 ; wa(iwa)=0.d0
    endif
    if(qtot_present) then
        iwa=iwa+1 ; wa(iwa)=1.d0
    else
        iwa=iwa+1 ; wa(iwa)=0.d0
    endif
    if(bemoved_present) then
        iwa=iwa+1 ; wa(iwa)=1.d0
    else
        iwa=iwa+1 ; wa(iwa)=0.d0
    endif
    if(vat_present) then
        iwa=iwa+1 ; wa(iwa)=1.d0
    else
        iwa=iwa+1 ; wa(iwa)=0.d0
    endif
    if(trim(atoms%boundcond)=='free') then
        iwa=iwa+1 ; wa(iwa)=0.d0
    elseif(trim(atoms%boundcond)=='wire') then
        iwa=iwa+1 ; wa(iwa)=1.d0
    elseif(trim(atoms%boundcond)=='slab') then
        iwa=iwa+1 ; wa(iwa)=2.d0
    elseif(trim(atoms%boundcond)=='bulk') then
        iwa=iwa+1 ; wa(iwa)=3.d0
    else
        write(*,*) "ERROR: unknown boundary conditions in write_bin_conf_v1, ",trim(atoms%boundcond)
        stop
    endif
    iwa=iwa+1 ; wa(iwa)=real(atoms%nat,8)
    call get_rat(atoms,rat)
    do iat=1,atoms%nat
        iwa=iwa+1 ; wa(iwa)=rat(1,iat)
        iwa=iwa+1 ; wa(iwa)=rat(2,iat)
        iwa=iwa+1 ; wa(iwa)=rat(3,iat)
    enddo
    do iat=1,atoms%nat
        call sat_to_iatom(atoms%sat(iat),isat)
        iwa=iwa+1 ; wa(iwa)=real(isat,8)
    enddo
    if(cell_present) then
        iwa=iwa+1 ; wa(iwa)=atoms%cellvec(1,1)
        iwa=iwa+1 ; wa(iwa)=atoms%cellvec(2,1)
        iwa=iwa+1 ; wa(iwa)=atoms%cellvec(3,1)
        iwa=iwa+1 ; wa(iwa)=atoms%cellvec(1,2)
        iwa=iwa+1 ; wa(iwa)=atoms%cellvec(2,2)
        iwa=iwa+1 ; wa(iwa)=atoms%cellvec(3,2)
        iwa=iwa+1 ; wa(iwa)=atoms%cellvec(1,3)
        iwa=iwa+1 ; wa(iwa)=atoms%cellvec(2,3)
        iwa=iwa+1 ; wa(iwa)=atoms%cellvec(3,3)
    endif
    if(epot_present) then
        iwa=iwa+1 ; wa(iwa)=atoms%epot
    endif
    if(fat_present) then
        do iat=1,atoms%nat
            iwa=iwa+1 ; wa(iwa)=atoms%fat(1,iat)
            iwa=iwa+1 ; wa(iwa)=atoms%fat(2,iat)
            iwa=iwa+1 ; wa(iwa)=atoms%fat(3,iat)
        enddo
    endif
    if(qtot_present) then
        iwa=iwa+1 ; wa(iwa)=atoms%qtot
    endif
    if(bemoved_present) then
        do iat=1,atoms%nat
            if(atoms%bemoved(1,iat)) then
                i1=1
            else
                i1=0
            endif
            if(atoms%bemoved(2,iat)) then
                i2=1
            else
                i2=0
            endif
            if(atoms%bemoved(3,iat)) then
                i3=1
            else
                i3=0
            endif
            ibm=i1+i2*2+i3*4
            iwa=iwa+1 ; wa(iwa)=real(ibm,8)
        enddo
    endif
    if(vat_present) then
        do iat=1,atoms%nat
            iwa=iwa+1 ; wa(iwa)=atoms%vat(1,iat)
            iwa=iwa+1 ; wa(iwa)=atoms%vat(2,iat)
            iwa=iwa+1 ; wa(iwa)=atoms%vat(3,iat)
        enddo
    endif
    write(iunit,iostat=ios,pos=ipos) wa(1:iwa)
    if(ios/=0) then
        write(*,*) 'ERROR: cannot write wa in write_bin_conf_v1, ios= ',ios
        stop
    endif
    call f_free(wa)
    deallocate(rat)
end subroutine write_bin_conf_v1
!*****************************************************************************************
end module mod_bin
!*****************************************************************************************
