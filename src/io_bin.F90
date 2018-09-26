!*****************************************************************************************
subroutine read_bin_conf(parini,filename,atoms_arr)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr
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
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr
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
            iwa=iwa+1 ; atoms_arr%atoms(iconf)%rat(1,iat)=wa(iwa)
            iwa=iwa+1 ; atoms_arr%atoms(iconf)%rat(2,iat)=wa(iwa)
            iwa=iwa+1 ; atoms_arr%atoms(iconf)%rat(3,iat)=wa(iwa)
        enddo
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
