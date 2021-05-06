!*****************************************************************************************
!This file contains routines to read/write ACF files.
!ACF: Alborz Configuration Format
!*****************************************************************************************
module mod_acf
    implicit none
    private
    public:: acf_write, acf_write_new, rotate4acf
    public:: acf_force_write, acf_read, acf_read_new
    public:: str_motion2bemoved
contains
!*****************************************************************************************
subroutine acf_write(file_info,atoms,atoms_all,strkey)
    use mod_atoms, only: typ_file_info, typ_atoms, typ_atoms_all, atom_copy_old
    use mod_atoms, only: atom_deallocate_old, update_ratp, set_rat_atoms, set_rat
    use mod_atoms, only: update_rat
    use mod_const, only: bohr2ang, ha2ev
    implicit none
    type(typ_file_info), intent(inout):: file_info
    type(typ_atoms), optional, intent(in):: atoms
    type(typ_atoms_all), optional, intent(in):: atoms_all
    character(*), optional, intent(in):: strkey
    !local variables
    integer:: iat, ios, nconf, iconf
    real(8):: x, y, z, epot,epot2, xt, yt, zt, hinv(3,3)
    real(8):: anrm, bnrm, tt1, tt2, tt3, dproj(6), cv(3,3), rotmat(3,3)
    logical:: l1, l2, l3
    type(typ_atoms):: atoms_t
    character(5):: tch1
    character(40):: tch2
    character(10):: units
    character(10):: coord
    character(256):: str_line6
    if(present(atoms) .and. present(atoms_all)) then
        write(*,'(a)') 'ERROR: in acf_write, either atoms or atoms_all can be present.'
        stop
    else if(.not. present(atoms) .and. .not. present(atoms_all)) then
        write(*,'(a)') 'ERROR: in acf_write, atoms or atoms_all must be present.'
        stop
    endif
    if(present(strkey)) then
        if(len_trim(strkey)>35) then
            stop 'ERROR: strkey too long in acf_write'
        endif
    endif
    if(trim(file_info%filename_positions)=='unknown') then
        stop 'ERROR: in acf_force_write filename is unknown'
    endif
    if(trim(file_info%file_position)=='new') then
        open(unit=1358,file=trim(file_info%filename_positions),status='replace',iostat=ios)
        file_info%nconf=0
    elseif(trim(file_info%file_position)=='append') then
        open(unit=1358,file=(file_info%filename_positions),status='old',position='append',iostat=ios)
    else
        write(*,'(2a)') 'ERROR: in acf_write file_position is unknown, ',trim(file_info%file_position)
        stop
    endif
    if(ios/=0) then
        write(*,'(2a)') 'ERROR: failure openning ',trim(file_info%filename_positions)
        stop
    endif
    if(present(atoms)) then
        nconf=1
        call atom_copy_old(atoms,atoms_t,'atoms->atoms_t')
    else if(present(atoms_all)) then 
        nconf=atoms_all%nconf
        call atom_copy_old(atoms_all%atoms,atoms_t,'atoms_all%atoms->atoms_t')
    else
        write(*,'(a)') 'ERROR: in acf_write, what is going on?'
        stop
    endif
    if(present(atoms)) then
        coord=trim(atoms%coordinates_type)
    endif
    if(present(atoms_all)) then
        coord=trim(atoms_all%atoms%coordinates_type)
    endif
    cv(1:3,1:3)=atoms_t%cellvec(1:3,1:3)
    if(trim(atoms_t%boundcond)=='bulk') then
        !write(*,'(a,3es9.0)') 'BEFORE ',cv(1,1),cv(1,2),cv(1,3)
        !write(*,'(a,3es9.0)') 'BEFORE ',cv(2,1),cv(2,2),cv(2,3)
        !write(*,'(a,3es9.0)') 'BEFORE ',cv(3,1),cv(3,2),cv(3,3)
        call update_ratp(atoms_t)
        call latvec2dproj_alborz(dproj,cv,rotmat,atoms_t%ratp,atoms_t%nat)
        call update_rat(atoms_t,upall=.true.)
        !write(*,'(a,3es9.0)') 'AFTER  ',cv(1,1),cv(1,2),cv(1,3)
        !write(*,'(a,3es9.0)') 'AFTER  ',cv(2,1),cv(2,2),cv(2,3)
        !write(*,'(a,3es9.0)') 'AFTER  ',cv(3,1),cv(3,2),cv(3,3)
        !write(*,'(a)') '--------------------------------------------------------------'
        anrm=sqrt(cv(1,1)**2+cv(2,1)**2+cv(3,1)**2)
        bnrm=sqrt(cv(1,2)**2+cv(2,2)**2+cv(3,2)**2)
        tt1=abs(cv(2,1))/anrm
        tt2=abs(cv(3,1))/anrm
        tt3=abs(cv(3,2))/bnrm
        if(tt1>1.d-8 .or. tt2>1.d-8 .or. tt3>1.d-8) then
            write(*,'(a,3es9.0)') 'ERROR: cell vectors not properly rotated for acf format',tt1,tt2,tt3
            stop
        endif
    endif
    if(trim(coord)=='reduced') call invertmat_alborz(cv,hinv)
    if(present(atoms)) then
        units=trim(atoms%units)
    endif
    if(present(atoms_all)) then
        units=trim(atoms_all%atoms%units)
    endif
    if(trim(units)=='angstrom') then
        !atoms_t%cellvec(1:3,1:3)=atoms_t%cellvec(1:3,1:3)*bohr2ang
        cv(1:3,1:3)=cv(1:3,1:3)*bohr2ang
    endif
    do iconf=1,nconf
        write(1358,'(a)') '#1st line comment'
        write(1358,'(a)') '#2nd line comment'
        write(1358,'(a)') '#3rd line comment'
        write(1358,*)
        write(1358,*)
        write(str_line6,'(a)') 'c5=bemoved'
        if(trim(units)=='atomic') then
            write(str_line6,'(a,1x,a,1x,a)') trim(str_line6),'units=',trim(units)
            !write(91,*) trim(str_line6)
        endif
        write(1358,'(a)') trim(str_line6)
        if(present(atoms)) then
            epot=atoms%epot
            epot2=epot-atoms%ebattery
            tch2=trim(strkey)
        else if(present(atoms_all)) then 
            epot=atoms_all%epotall(iconf)
            write(tch1,'(i5.5)') iconf
            tch2=trim(strkey)//tch1
        else
            write(*,'(a)') 'ERROR: in acf_write, what is going on?'
            stop
        endif
        if(present(strkey)) then
            write(1358,'(2(a,1x))',advance='no') 'label=',trim(tch2)
        endif
        if(trim(units)=='angstrom') then
            epot=epot*ha2ev
            epot2=epot2*ha2ev
        endif
            !write(*,*)"externalwork2", atoms%ebattery
        if(atoms%ebattery > 0.d0 .or. atoms%ebattery < 0.d0 ) then
            write(1358,'(a,e24.15,a10,e24.15)') 'epot=',epot ,'epot2=',epot2
        else
            write(1358,'(a,e24.15)') 'epot=',epot
        endif
        file_info%nconf=file_info%nconf+1
        write(1358,'(i7,2x,a,2x,i7)') atoms_t%nat,trim(atoms_t%boundcond),file_info%nconf
        if(present(atoms)) then
            call set_rat_atoms(atoms_t,atoms,setall=.true.)
        else if(present(atoms_all)) then 
            call set_rat(atoms_t,atoms_all%ratall(1,1,iconf),setall=.true.)
        else
            write(*,'(a)') 'ERROR: in acf_write, what is going on?'
            stop
        endif
        if(trim(atoms_t%boundcond)=='bulk') then
            call update_ratp(atoms_t)
            call rotate4acf(atoms_t%nat,atoms_t%ratp(1,1),atoms_t%cellvec,cv)
            call update_rat(atoms_t,upall=.true.)
        endif
        write(1358,'(3es24.15)') cv(1,1),cv(1,2),cv(2,2)
        write(1358,'(3es24.15)') cv(1,3),cv(2,3),cv(3,3)
        call update_ratp(atoms_t)
        do iat=1,atoms_t%nat
            xt=atoms_t%ratp(1,iat)
            yt=atoms_t%ratp(2,iat)
            zt=atoms_t%ratp(3,iat)
            if(trim(coord)=='cartesian') then
                if(trim(units)=='angstrom' .and. trim(atoms_t%boundcond)/='bulk') then
                    x=xt*bohr2ang
                    y=yt*bohr2ang
                    z=zt*bohr2ang
                else
                    x=xt
                    y=yt
                    z=zt
                endif
            else
                !It must be reduced coordinates so applying r=hs
                x=hinv(1,1)*xt+hinv(1,2)*yt+hinv(1,3)*zt
                y=hinv(2,1)*xt+hinv(2,2)*yt+hinv(2,3)*zt
                z=hinv(3,1)*xt+hinv(3,2)*yt+hinv(3,3)*zt
            endif
            l1=atoms_t%bemoved(1,iat);l2=atoms_t%bemoved(2,iat);l3=atoms_t%bemoved(3,iat)
            !write(1358,'(a5,2x,3e22.13,2x,3l1)') trim(atoms_t%sat(iat)),x,y,z,l1,l2,l3
            write(1358,'(a5,2x,3es24.15,2x,3l1)') trim(atoms_t%sat(iat)),x,y,z,l1,l2,l3
        enddo
    enddo
    close(1358)
    call atom_deallocate_old(atoms_t)
    if(file_info%print_force) then
        if(present(atoms)) then
            if(present(strkey)) then
                call acf_force_write(file_info,atoms=atoms,strkey=strkey)
            else
                call acf_force_write(file_info,atoms=atoms)
            endif
        else if(present(atoms_all)) then 
            if(present(strkey)) then
                call acf_force_write(file_info,atoms_all=atoms_all,strkey=strkey)
            else
                call acf_force_write(file_info,atoms_all=atoms_all)
            endif
        else
            write(*,'(a)') 'ERROR: in acf_write, what is going on?'
            stop
        endif
    endif
end subroutine acf_write
!*****************************************************************************************
subroutine acf_write_new(file_info,atoms_arr,strkey)
    use mod_atoms, only: typ_file_info, typ_atoms, typ_atoms_arr, atom_copy_old
    use mod_atoms, only: atom_deallocate_old, update_ratp, update_rat
    use mod_const, only: bohr2ang, ha2ev
    implicit none
    type(typ_file_info), intent(inout):: file_info
    type(typ_atoms_arr),intent(in):: atoms_arr
    character(*), optional, intent(in):: strkey
    !local variables
    integer:: iat, ios, nconf, iconf
    real(8):: x, y, z, epot, xt, yt, zt, hinv(3,3)
    real(8):: anrm, bnrm, tt1, tt2, tt3, dproj(6), cv(3,3), rotmat(3,3)
    logical:: l1, l2, l3
    type(typ_atoms):: atoms_t
    character(5):: tch1
    character(40):: tch2
    character(10):: units
    character(10):: coord
    character(256):: str_line6
    !if(present(atoms) .and. present(atoms_all)) then
    !    write(*,'(a)') 'ERROR: in acf_write, either atoms or atoms_all can be present.'
    !    stop
    !else if(.not. present(atoms) .and. .not. present(atoms_all)) then
    !    write(*,'(a)') 'ERROR: in acf_write, atoms or atoms_all must be present.'
    !    stop
    !endif
    if(present(strkey)) then
        if(len_trim(strkey)>35) then
            stop 'ERROR: strkey too long in acf_write'
        endif
    endif
    if(trim(file_info%filename_positions)=='unknown') then
        stop 'ERROR: in acf_force_write filename is unknown'
    endif
    if(trim(file_info%file_position)=='new') then
        open(unit=1358,file=trim(file_info%filename_positions),status='replace',iostat=ios)
        file_info%nconf=0
    elseif(trim(file_info%file_position)=='append') then
        open(unit=1358,file=(file_info%filename_positions),status='old',position='append',iostat=ios)
    else
        write(*,'(2a)') 'ERROR: in acf_write file_position is unknown, ',trim(file_info%file_position)
        stop
    endif
    if(ios/=0) then
        write(*,'(2a)') 'ERROR: failure openning ',trim(file_info%filename_positions)
        stop
    endif
    configurations: do iconf=1,atoms_arr%nconf    
        call atom_copy_old(atoms_arr%atoms(iconf),atoms_t,'atoms_arr%atoms(iconf)->atoms_t')
        coord=trim(atoms_arr%atoms(iconf)%coordinates_type)
        cv(1:3,1:3)=atoms_t%cellvec(1:3,1:3)
    if(trim(atoms_t%boundcond)=='bulk') then
        !write(*,'(a,3es9.0)') 'BEFORE ',cv(1,1),cv(1,2),cv(1,3)
        !write(*,'(a,3es9.0)') 'BEFORE ',cv(2,1),cv(2,2),cv(2,3)
        !write(*,'(a,3es9.0)') 'BEFORE ',cv(3,1),cv(3,2),cv(3,3)
        call update_ratp(atoms_t)
        call latvec2dproj_alborz(dproj,cv,rotmat,atoms_t%ratp,atoms_t%nat)
        call update_rat(atoms_t,upall=.true.)
        !write(*,'(a,3es9.0)') 'AFTER  ',cv(1,1),cv(1,2),cv(1,3)
        !write(*,'(a,3es9.0)') 'AFTER  ',cv(2,1),cv(2,2),cv(2,3)
        !write(*,'(a,3es9.0)') 'AFTER  ',cv(3,1),cv(3,2),cv(3,3)
        !write(*,'(a)') '--------------------------------------------------------------'
        anrm=sqrt(cv(1,1)**2+cv(2,1)**2+cv(3,1)**2)
        bnrm=sqrt(cv(1,2)**2+cv(2,2)**2+cv(3,2)**2)
        tt1=abs(cv(2,1))/anrm
        tt2=abs(cv(3,1))/anrm
        tt3=abs(cv(3,2))/bnrm
        if(tt1>1.d-8 .or. tt2>1.d-8 .or. tt3>1.d-8) then
            write(*,'(a,3es9.0)') 'ERROR: cell vectors not properly rotated for acf format',tt1,tt2,tt3
            stop
        endif
    endif
    if(trim(coord)=='reduced') call invertmat_alborz(cv,hinv)
    units=trim(atoms_arr%atoms(iconf)%units)
    if(trim(units)=='angstrom') then
        !atoms_t%cellvec(1:3,1:3)=atoms_t%cellvec(1:3,1:3)*bohr2ang
        cv(1:3,1:3)=cv(1:3,1:3)*bohr2ang
    endif
    !do iconf=1,atoms_arr%nconf
        write(1358,'(a)') '#1st line comment'
        write(1358,'(a)') '#2nd line comment'
        write(1358,'(a)') '#3rd line comment'
        write(1358,*)
        write(1358,*)
        write(str_line6,'(a)') 'c5=bemoved'
        if(trim(units)=='atomic') then
            write(str_line6,'(a,1x,a,1x,a)') trim(str_line6),'units=',trim(units)
            !write(91,*) trim(str_line6)
        endif
        write(1358,'(a)') trim(str_line6)
        !epot=atoms_arr%epotall(iconf)
            epot=atoms_t%epot
            write(tch1,'(i5.5)') iconf
            tch2=trim(strkey)//tch1
        if(present(strkey)) then
            write(1358,'(2(a,1x))',advance='no') 'label=',trim(tch2)
        endif
        if(trim(units)=='angstrom') then
            epot=epot*ha2ev
        endif
        write(1358,'(a,e24.15)') 'epot=',epot
        file_info%nconf=file_info%nconf+1
        write(1358,'(i7,2x,a,2x,i7)') atoms_t%nat,trim(atoms_t%boundcond),file_info%nconf
        !if(present(atoms)) then
        !    do iat=1,atoms_t%nat
        !        atoms_t%rat(1,iat)=atoms%rat(1,iat)
        !        atoms_t%rat(2,iat)=atoms%rat(2,iat)
        !        atoms_t%rat(3,iat)=atoms%rat(3,iat)
        !    enddo
        !else if(present(atoms_all)) then 
        !    do iat=1,atoms_t%nat
        !        atoms_t%rat(1,iat)=atoms_all%ratall(1,iat,iconf)
        !        atoms_t%rat(2,iat)=atoms_all%ratall(2,iat,iconf)
        !        atoms_t%rat(3,iat)=atoms_all%ratall(3,iat,iconf)
        !    enddo
        !else
        !    write(*,'(a)') 'ERROR: in acf_write, what is going on?'
        !    stop
        !endif
        if(trim(atoms_t%boundcond)=='bulk') then
            call update_ratp(atoms_t)
            call rotate4acf(atoms_t%nat,atoms_t%ratp(1,1),atoms_t%cellvec,cv)
            call update_rat(atoms_t,upall=.true.)
        endif
        write(1358,'(3es24.15)') cv(1,1),cv(1,2),cv(2,2)
        write(1358,'(3es24.15)') cv(1,3),cv(2,3),cv(3,3)
        atoms: do iat=1,atoms_t%nat
            xt=atoms_t%ratp(1,iat)
            yt=atoms_t%ratp(2,iat)
            zt=atoms_t%ratp(3,iat)
            if(trim(coord)=='cartesian') then
                if(trim(units)=='angstrom' .and. trim(atoms_t%boundcond)/='bulk') then
                    x=xt*bohr2ang
                    y=yt*bohr2ang
                    z=zt*bohr2ang
                else
                    x=xt
                    y=yt
                    z=zt
                endif
            else
                !It must be reduced coordinates so applying r=hs
                x=hinv(1,1)*xt+hinv(1,2)*yt+hinv(1,3)*zt
                y=hinv(2,1)*xt+hinv(2,2)*yt+hinv(2,3)*zt
                z=hinv(3,1)*xt+hinv(3,2)*yt+hinv(3,3)*zt
            endif
            l1=atoms_t%bemoved(1,iat);l2=atoms_t%bemoved(2,iat);l3=atoms_t%bemoved(3,iat)
            !write(1358,'(a5,2x,3e22.13,2x,3l1)') trim(atoms_t%sat(iat)),x,y,z,l1,l2,l3
            write(1358,'(a5,2x,3es24.15,2x,3l1)') trim(atoms_t%sat(iat)),x,y,z,l1,l2,l3
        enddo atoms
    !if(file_info%print_force) then
    !    if(present(atoms)) then
    !        if(present(strkey)) then
    !            call acf_force_write(file_info,atoms=atoms,strkey=strkey)
    !        else
    !            call acf_force_write(file_info,atoms=atoms)
    !        endif
    !    else if(present(atoms_all)) then 
    !        if(present(strkey)) then
    !            call acf_force_write(file_info,atoms_all=atoms_all,strkey=strkey)
    !        else
    !            call acf_force_write(file_info,atoms_all=atoms_all)
    !        endif
    !    else
    !        write(*,'(a)') 'ERROR: in acf_write, what is going on?'
    !        stop
    !    endif
    !endif
    enddo configurations
    close(1358)
    call atom_deallocate_old(atoms_t)
end subroutine acf_write_new
!*****************************************************************************************
subroutine rotate4acf(nat,rat,cv,cvrot)
    implicit none
    integer, intent(in):: nat
    real(8), intent(inout):: rat(3,nat)
    real(8), intent(in):: cv(3,3), cvrot(3,3)
    !local variables
    real(8):: cvinv(3,3), x, y, z
    integer:: iat
    !real(8), allocatable:: ratred(:,:)
    call invertmat_alborz(cv,cvinv)
    do iat=1,nat
        x=cvinv(1,1)*rat(1,iat)+cvinv(1,2)*rat(2,iat)+cvinv(1,3)*rat(3,iat)
        y=cvinv(2,1)*rat(1,iat)+cvinv(2,2)*rat(2,iat)+cvinv(2,3)*rat(3,iat)
        z=cvinv(3,1)*rat(1,iat)+cvinv(3,2)*rat(2,iat)+cvinv(3,3)*rat(3,iat)
        x=modulo(modulo(x,1.d0),1.d0)
        y=modulo(modulo(y,1.d0),1.d0)
        z=modulo(modulo(z,1.d0),1.d0)
        rat(1,iat)=cvrot(1,1)*x+cvrot(1,2)*y+cvrot(1,3)*z
        rat(2,iat)=cvrot(2,1)*x+cvrot(2,2)*y+cvrot(2,3)*z
        rat(3,iat)=cvrot(3,1)*x+cvrot(3,2)*y+cvrot(3,3)*z
    enddo
end subroutine rotate4acf
!*****************************************************************************************
subroutine acf_force_write(file_info,atoms,atoms_all,strkey)
    use mod_atoms, only: typ_file_info, typ_atoms, typ_atoms_all, atom_copy_old
    use mod_atoms, only: atom_deallocate_old
    implicit none
    type(typ_file_info), intent(inout):: file_info
    type(typ_atoms), optional, intent(in):: atoms
    type(typ_atoms_all), optional, intent(in):: atoms_all
    character(*), optional, intent(in):: strkey
    !local variables
    integer:: iat, ios, nconf, iconf
    real(8):: x, y, z, epot
    logical:: l1, l2, l3
    type(typ_atoms):: atoms_t
    character(5):: tch1
    character(40):: tch2
    if(present(atoms) .and. present(atoms_all)) then
        write(*,'(a)') 'ERROR: in acf_force_write, either atoms or atoms_all can be present.'
        stop
    else if(.not. present(atoms) .and. .not. present(atoms_all)) then
        write(*,'(a)') 'ERROR: in acf_force_write, atoms or atoms_all must be present.'
        stop
    endif
    if(present(strkey)) then
        if(len_trim(strkey)>29) then
            stop 'ERROR: strkey too long in acf_force_write'
        endif
    endif
    if(trim(file_info%filename_forces)=='unknown') then
        if(trim(file_info%filename_positions)=='unknown') then
            stop 'ERROR: in acf_force_write filename is unknown'
        endif
        file_info%filename_forces='force_'//trim(file_info%filename_positions)
    endif
    if(trim(file_info%file_position)=='new') then
        open(unit=1358,file=trim(file_info%filename_forces),status='replace',iostat=ios)
        file_info%nconf=0
    elseif(trim(file_info%file_position)=='append') then
        open(unit=1358,file=trim(file_info%filename_forces),status='old',position='append',iostat=ios)
    else
        write(*,'(2a)') 'ERROR: in acf_force_write file_position is unknown, ',trim(file_info%file_position)
        stop
    endif
    if(ios/=0) then
        write(*,'(2a)') 'ERROR: failure openning ',trim(file_info%filename_forces)
        stop
    endif
    if(present(atoms)) then
        nconf=1
        call atom_copy_old(atoms,atoms_t,'atoms->atoms_t')
    else if(present(atoms_all)) then 
        nconf=atoms_all%nconf
        call atom_copy_old(atoms_all%atoms,atoms_t,'atoms_all%atoms->atoms_t')
    else
        write(*,'(a)') 'ERROR: in acf_force_write, what is going on?'
        stop
    endif
    do iconf=1,nconf
        write(1358,*)
        write(1358,*)
        write(1358,*)
        write(1358,*)
        write(1358,*)
        write(1358,*)
        if(present(atoms)) then
            epot=atoms%epot
            tch2=trim(strkey)
        else if(present(atoms_all)) then 
            epot=atoms_all%epotall(iconf)
            write(tch1,'(i5.5)') iconf
            tch2=trim(strkey)//tch1
        else
            write(*,'(a)') 'ERROR: in acf_force_write, what is going on?'
            stop
        endif
        if(present(strkey)) then
            write(1358,'(2(a,1x))',advance='no') 'label=',trim(tch2)
        endif
        write(1358,'(a,e24.15)') 'epot=',epot
        file_info%nconf=file_info%nconf+1
        write(1358,'(i7,2x,a,2x,i7)') atoms_t%nat,trim(atoms_t%boundcond),file_info%nconf
        write(1358,*)
        write(1358,*)
        do iat=1,atoms_t%nat
            if(present(atoms)) then
                x=atoms%fat(1,iat)
                y=atoms%fat(2,iat)
                z=atoms%fat(3,iat)
            else if(present(atoms_all)) then 
                x=atoms_all%fatall(1,iat,iconf)
                y=atoms_all%fatall(2,iat,iconf)
                z=atoms_all%fatall(3,iat,iconf)
            else
                write(*,'(a)') 'ERROR: in acf_force_write, what is going on?'
                stop
            endif
            l1=atoms_t%bemoved(1,iat);l2=atoms_t%bemoved(2,iat);l3=atoms_t%bemoved(3,iat)
            !write(1358,'(a5,2x,3e22.13,2x,3l1)') trim(atoms_t%sat(iat)),x,y,z,l1,l2,l3
            !write(1358,'(a5,2x,3es24.15,2x,3l1)') trim(atoms_t%sat(iat)),x,y,z,l1,l2,l3
            !write(1358,'(3f12.8)') x,y,z
            write(1358,trim(file_info%frmt)) trim(atoms_t%sat(iat)),x,y,z,l1,l2,l3
        enddo
    enddo
    close(1358)
    call atom_deallocate_old(atoms_t)
end subroutine acf_force_write
!*****************************************************************************************
subroutine acf_read(parini,filename,nconfmax,atoms,atoms_all)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_atoms_all, atom_all_allocate, atom_copy_old
    use mod_atoms, only: atom_allocate_old, update_ratp, get_rat, update_rat
    use mod_const, only: bohr2ang, ha2ev
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: filename
    integer, intent(in):: nconfmax
    type(typ_atoms), optional, intent(inout):: atoms
    type(typ_atoms_all), optional, intent(inout):: atoms_all
    !local variables
    type(typ_atoms):: atoms_t
    integer:: ios, iconf, i, nconf, iat, k, nat
    character(256):: str
    character(256):: fn_fullpath
    character(3):: str_motion
    character(5):: s
    character(50):: c5
    real(8):: t1, t2, t3, x, y, z
    if(present(atoms) .and. present(atoms_all)) then
        write(*,'(a)') 'ERROR: in acf_read, either atoms or atoms_all can be present.'
        stop
    else if(.not. present(atoms) .and. .not. present(atoms_all)) then
        write(*,'(a)') 'ERROR: in acf_read, atoms or atoms_all must be present.'
        stop
    endif
    if(nconfmax<1) then
        write(*,'(a)') 'ERROR: why do you call acf_read with nconfmax<1 ?'
        stop
    endif
    if(present(atoms) .and. nconfmax/=1) then
        write(*,'(a)') 'ERROR: in acf_read: present(atoms) and nconfmax/=1'
    endif
    !call atom_copy(atoms,atoms_t,'atoms->atoms_t')
    !call atom_copy(atoms_all%atoms,atoms_t,'atoms_all%atoms->atoms_t')
    !atom_allocate_old(atoms,sat,rat,vat,amass,fat,bemoved,qat,rcov,fp)
    fn_fullpath=trim(parini%cwd)//'/'//trim(filename)
    open(unit=1358, file=trim(fn_fullpath),status='old',iostat=ios)
    if(ios/=0) then;write(*,'(2a)') 'ERROR: failure openning ',trim(fn_fullpath);stop;endif
    iconf=0
    do
        read(1358,*,iostat=k)
        if(k<0) exit
        read(1358,*)
        read(1358,*)

        c5='nothing'
        do i=4,7
            read(1358,'(a256)') str
            call str_parse(str,i,atoms_t,c5)
        enddo
        read(1358,*) nat,atoms_t%boundcond,nconf
        if(nconf/=iconf+1) then
            write(*,'(3a,2i6)') 'ERROR: inconsistent configuration index: ', &
                'filename,nconf,iconf ',trim(filename),nconf,iconf
        endif
        !if(iconf==0) nconf_p=nconf
        !if(iconf>0 .and. nconf/=nconf_p-1) then
        !    write(*,'(3a,2i6)') 'ERROR: inconsistent configuration index: ', &
        !        'filename,nconf,nconf_p ',trim(filename),nconf,nconf_p
        !endif
        read(1358,*) atoms_t%cellvec(1,1),atoms_t%cellvec(1,2),atoms_t%cellvec(2,2)
        read(1358,*) atoms_t%cellvec(1,3),atoms_t%cellvec(2,3),atoms_t%cellvec(3,3)
        atoms_t%cellvec(2,1)=0.d0 ; atoms_t%cellvec(3,1)=0.d0 ; atoms_t%cellvec(3,2)=0.d0
        if(trim(atoms_t%units)=='angstrom') then
            atoms_t%cellvec(1:3,1:3)=atoms_t%cellvec(1:3,1:3)/bohr2ang
            atoms_t%epot=atoms_t%epot/ha2ev
        endif
        !write(*,'(i7,2x,a,2x,i7)') nat,trim(atoms_t%boundcond),nconf
        !write(*,'(3es24.15)') atoms_t%cellvec(1,1),atoms_t%cellvec(1,2),atoms_t%cellvec(2,2)
        !write(*,'(3es24.15)') atoms_t%cellvec(1,3),atoms_t%cellvec(2,3),atoms_t%cellvec(3,3)
        if(iconf==0) then
            call atom_allocate_old(atoms_t,nat,0,0)
        endif
        do iat=1,atoms_t%nat
            read(1358,'(a256)') str
            !write(*,*) trim(str)
            !write(*,*) allocated(atoms_t%rat)
            read(str,*) atoms_t%sat(iat),x,y,z
            if(trim(atoms_t%coordinates_type)=='cartesian') then
                if(trim(atoms_t%units)=='angstrom') then
                    atoms_t%ratp(1,iat)=x/bohr2ang
                    atoms_t%ratp(2,iat)=y/bohr2ang
                    atoms_t%ratp(3,iat)=z/bohr2ang
                else
                    atoms_t%ratp(1,iat)=x
                    atoms_t%ratp(2,iat)=y
                    atoms_t%ratp(3,iat)=z
                endif
            else
                !It must be reduced coordinates so applying r=hs
                atoms_t%ratp(1,iat)=atoms_t%cellvec(1,1)*x+atoms_t%cellvec(1,2)*y+atoms_t%cellvec(1,3)*z
                atoms_t%ratp(2,iat)=atoms_t%cellvec(2,1)*x+atoms_t%cellvec(2,2)*y+atoms_t%cellvec(2,3)*z
                atoms_t%ratp(3,iat)=atoms_t%cellvec(3,1)*x+atoms_t%cellvec(3,2)*y+atoms_t%cellvec(3,3)*z
            endif
            if(trim(c5)=='bemoved') then
                read(str,*) s,t1,t2,t3,str_motion
                !write(*,*) trim(str_motion)
                call str_motion2bemoved(str_motion,atoms_t%bemoved(1,iat))
            endif
            !write(*,'(a5,2x,3es24.15,2x,3l1)') trim(atoms_t%sat(iat)), &
            !    atoms_t%rat(1,iat),atoms_t%rat(2,iat),atoms_t%rat(3,iat), &
            !    atoms_t%bemoved(1,iat),atoms_t%bemoved(2,iat),atoms_t%bemoved(3,iat)
        enddo
        call update_rat(atoms_t,upall=.true.)

        if(present(atoms)) then
            call atom_copy_old(atoms_t,atoms,'atoms_t->atoms')
        endif
        if(iconf==0 .and. present(atoms_all)) then
            !atoms_all%atoms%nat is set in next line. It is needed by atom_all_allocate
            call atom_copy_old(atoms_t,atoms_all%atoms,'atoms_t->atoms_all%atoms')
            atoms_all%nconfmax=nconfmax
            call atom_all_allocate(atoms_all,ratall=.true.,fatall=.true.,epotall=.true.,qtotall=.true.)
        endif
        if(present(atoms_all)) then
            !atoms_all%ratall(1:3,1:atoms_t%nat,iconf+1)=atoms_t%rat(1:3,1:atoms_t%nat)
            call get_rat(atoms_t,atoms_all%ratall(1,1,iconf+1))
            atoms_all%epotall(iconf+1)=atoms_t%epot
            atoms_all%qtotall(iconf+1)=atoms_t%qtot
        endif
        !call atom_allocate_old(atoms_t,sat,rat,vat,amass,fat,bemoved,qat,rcov,fp)
        iconf=iconf+1
        if(iconf>=nconfmax) exit
        !nconf_p=nconf
    enddo !end of loop over iconf
    if(present(atoms_all)) then
        atoms_all%nconf=iconf
    endif
    call yaml_mapping_open('reading structures')
    call yaml_map('filename',trim(filename))
    call yaml_map('nconf',iconf)
    call yaml_mapping_close()
    !write(*,'(3(a,1x),i4)') 'Number of configurations read from',trim(filename),'is',iconf
    close(1358)
end subroutine acf_read
!*****************************************************************************************
subroutine acf_read_new(parini,filename,nconfmax,atoms_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_allocate, atom_copy
    use mod_atoms, only: atom_deallocate, update_rat
    use mod_const, only: bohr2ang, ha2ev
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: filename
    integer, intent(in):: nconfmax
    type(typ_atoms_arr), intent(inout):: atoms_arr
    !local variables
    type(typ_atoms):: atoms_t
    integer:: ios, iconf, i, nconf, iat, k, nat
    character(256):: str, fn_tmp
    character(256):: fn_fullpath
    character(3):: str_motion
    character(5):: s
    character(50):: c5
    real(8):: t1, t2, t3, x, y, z
    if(nconfmax<1) then
        write(*,'(a)') 'ERROR: why do you call acf_read_new with nconfmax<1 ?'
        stop
    endif
    fn_tmp=adjustl(trim(filename))
    if(fn_tmp(1:1)=='/') then
        fn_fullpath=trim(filename)
    else
        fn_fullpath=trim(parini%cwd)//'/'//trim(filename)
    endif
    open(unit=1358, file=trim(fn_fullpath),status='old',iostat=ios)
    if(ios/=0) then;write(*,'(2a)') 'ERROR: failure openning ',trim(fn_fullpath);stop;endif
    allocate(atoms_arr%atoms(nconfmax))
    iconf=0
    do
        read(1358,*,iostat=k)
        if(k<0) exit
        read(1358,*)
        read(1358,*)

        c5='nothing'
        do i=4,7
            read(1358,'(a256)') str
            call str_parse(str,i,atoms_t,c5)
        enddo
        read(1358,*) nat,atoms_t%boundcond,nconf
        if(nconf/=iconf+1) then
            write(*,'(3a,2i6)') 'ERROR: inconsistent configuration index: ', &
                'filename,nconf,iconf ',trim(filename),nconf,iconf
        endif
        !if(iconf==0) nconf_p=nconf
        !if(iconf>0 .and. nconf/=nconf_p-1) then
        !    write(*,'(3a,2i6)') 'ERROR: inconsistent configuration index: ', &
        !        'filename,nconf,nconf_p ',trim(filename),nconf,nconf_p
        !endif
        read(1358,*) atoms_t%cellvec(1,1),atoms_t%cellvec(1,2),atoms_t%cellvec(2,2)
        read(1358,*) atoms_t%cellvec(1,3),atoms_t%cellvec(2,3),atoms_t%cellvec(3,3)
        atoms_t%cellvec(2,1)=0.d0 ; atoms_t%cellvec(3,1)=0.d0 ; atoms_t%cellvec(3,2)=0.d0
        if(trim(atoms_t%units)=='angstrom') then
            atoms_t%cellvec(1:3,1:3)=atoms_t%cellvec(1:3,1:3)/bohr2ang
            atoms_t%epot=atoms_t%epot/ha2ev
        endif
        !write(*,'(i7,2x,a,2x,i7)') nat,trim(atoms_t%boundcond),nconf
        !write(*,'(3es24.15)') atoms_t%cellvec(1,1),atoms_t%cellvec(1,2),atoms_t%cellvec(2,2)
        !write(*,'(3es24.15)') atoms_t%cellvec(1,3),atoms_t%cellvec(2,3),atoms_t%cellvec(3,3)
        if(iconf==0) then
            call atom_allocate(atoms_t,nat,0,0)
        endif
        do iat=1,atoms_t%nat
            read(1358,'(a256)') str
            !write(*,*) trim(str)
            !write(*,*) allocated(atoms_t%rat)
            read(str,*) atoms_t%sat(iat),x,y,z
            if(trim(atoms_t%coordinates_type)=='cartesian') then
                if(trim(atoms_t%units)=='angstrom') then
                    atoms_t%ratp(1,iat)=x/bohr2ang
                    atoms_t%ratp(2,iat)=y/bohr2ang
                    atoms_t%ratp(3,iat)=z/bohr2ang
                else
                    atoms_t%ratp(1,iat)=x
                    atoms_t%ratp(2,iat)=y
                    atoms_t%ratp(3,iat)=z
                endif
            else
                !It must be reduced coordinates so applying r=hs
                atoms_t%ratp(1,iat)=atoms_t%cellvec(1,1)*x+atoms_t%cellvec(1,2)*y+atoms_t%cellvec(1,3)*z
                atoms_t%ratp(2,iat)=atoms_t%cellvec(2,1)*x+atoms_t%cellvec(2,2)*y+atoms_t%cellvec(2,3)*z
                atoms_t%ratp(3,iat)=atoms_t%cellvec(3,1)*x+atoms_t%cellvec(3,2)*y+atoms_t%cellvec(3,3)*z
            endif
            if(trim(c5)=='bemoved') then
                read(str,*) s,t1,t2,t3,str_motion
                !write(*,*) trim(str_motion)
                call str_motion2bemoved(str_motion,atoms_t%bemoved(1,iat))
            endif
            !write(*,'(a5,2x,3es24.15,2x,3l1)') trim(atoms_t%sat(iat)), &
            !    atoms_t%rat(1,iat),atoms_t%rat(2,iat),atoms_t%rat(3,iat), &
            !    atoms_t%bemoved(1,iat),atoms_t%bemoved(2,iat),atoms_t%bemoved(3,iat)
        enddo
        call update_rat(atoms_t,upall=.true.)
        iconf=iconf+1
        call atom_copy(atoms_t,atoms_arr%atoms(iconf),'atoms_t->atoms_arr%atoms')
        if(iconf>=nconfmax) exit
        !nconf_p=nconf
    enddo !end of loop over iconf
    atoms_arr%nconf=iconf
    call atom_deallocate(atoms_t)
    call yaml_mapping_open('Number of configurations read',flow=.true.)
    call yaml_map('filename',trim(filename))
    call yaml_map('nconf',iconf)
    call yaml_mapping_close()
    !write(*,'(3(a,1x),i4)') 'Number of configurations read from',trim(filename),'is',iconf
    close(1358)
end subroutine acf_read_new
!*****************************************************************************************
subroutine str_parse(str,line,atoms_t,c5)
    use mod_atoms, only: typ_atoms, typ_atoms_all
    implicit none
    character(*), intent(in):: str
    integer, intent(in):: line
    type(typ_atoms), intent(inout):: atoms_t
    character(*), intent(inout):: c5
    !character(*), intent(inout):: coordinates
    !local variables
    integer:: idx
    character(256):: str_t
    character(10):: units
    character(20):: coord
    !-------------------------------------------------------
    !Is column five specified?
    idx=index(str,'c5=')
    !write(*,*) idx,trim(str)
    if(idx>0) then
        if(line==6) then
            str_t=str(idx+3:256)
            !write(*,*) trim(str_t)
            read(str_t,*) c5
        else
            write(*,'(a,i2)') 'ERROR: c5 can only be in 6th line for each conf. line=',line
            stop
        endif
    endif
    !-------------------------------------------------------
    !Is the input units specified?
    idx=index(str,'units=')
    if(idx>0) then
        if(line==6) then
            str_t=str(idx+6:256)
            read(str_t,*) units
            if(trim(units)=='angstrom' .or. trim(units)=='atomic') then
                atoms_t%units=trim(units)
            else
                stop 'ERROR: only angstrom and atomic are accepted in acf files.'
            endif
        else
            write(*,'(a,i2)') 'ERROR: units can only be in 6th line for each conf. line=',line
            stop
        endif
    endif
    !-------------------------------------------------------
    !Is the input coordinates specified?
    idx=index(str,'coord=')
    if(idx>0) then
        if(line==6) then
            str_t=str(idx+6:256)
            read(str_t,*) coord
            if(trim(coord)=='cartesian' .or. trim(coord)=='reduced') then
                atoms_t%coordinates_type=trim(coord)
            else
                stop 'ERROR: only angstrom and atomic are accepted in acf files.'
            endif
        else
            write(*,'(a,i2)') 'ERROR: coordinates can only be in 6th line for each conf. line=',line
            stop
        endif
    endif
    !-------------------------------------------------------
    !Is potential energy specified?
    idx=index(str,'epot=')
    !write(*,*) idx,trim(str)
    if(idx>0) then
        if(line==7) then
            str_t=str(idx+5:256)
            !write(*,*) trim(str_t)
            read(str_t,*) atoms_t%epot
        else
            write(*,'(a,i2)') 'ERROR: epot can only be in 7th line for each conf. line=',line
            stop
        endif
    endif
    !-------------------------------------------------------
    !Is total charge of the system specified?
    idx=index(str,'qtot=')
    !write(*,*) idx,trim(str)
    if(idx>0) then
        if(line==7) then
            str_t=str(idx+5:256)
            !write(*,*) trim(str_t)
            read(str_t,*) atoms_t%qtot
        else
            write(*,'(a,i2)') 'ERROR: qtot can only be in 7th line for each conf. line=',line
            stop
        endif
    endif
end subroutine str_parse
!*****************************************************************************************
subroutine str_motion2bemoved(str_motion,bemoved)
    implicit none
    character(*), intent(in):: str_motion
    logical, intent(inout):: bemoved(3)
    !local variables
    if(str_motion(1:1)=='T') then
        bemoved(1)=.true.
    elseif(str_motion(1:1)=='F') then
        bemoved(1)=.false.
    else
        stop 'ERROR: incorrect symbol for atomic motion in input configuration'
    endif
    if(str_motion(2:2)=='T') then
        bemoved(2)=.true.
    elseif(str_motion(2:2)=='F') then
        bemoved(2)=.false.
    else
        stop 'ERROR: incorrect symbol for atomic motion in input configuration'
    endif
    if(str_motion(3:3)=='T') then
        bemoved(3)=.true.
    elseif(str_motion(3:3)=='F') then
        bemoved(3)=.false.
    else
        stop 'ERROR: incorrect symbol for atomic motion in input configuration'
    endif
end subroutine str_motion2bemoved
!*****************************************************************************************
end module mod_acf
!*****************************************************************************************
