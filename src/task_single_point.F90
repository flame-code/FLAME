!*****************************************************************************************
subroutine single_point_task(parini)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr, typ_file_info
    use mod_potential, only: fcalls, perfstatus, potential
    use mod_processors, only: iproc
    use mod_const, only: ev2ha, ang2bohr
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    type(typ_atoms_arr):: atoms_arr
    type(typ_file_info):: file_info
    real(8):: tt1, tt2, fxyz(3)
    integer:: iconf, iat
    call acf_read_new(parini,'posinp.acf',10000,atoms_arr)
    do iconf=1,atoms_arr%nconf
        call set_ndof(atoms_arr%atoms(iconf))
    enddo
    potential=trim(parini%potential_potential)
    file_info%filename_positions='posout.acf'
    file_info%print_force=parini%print_force_single_point
    file_info%file_position='new'
    if(trim(parini%frmt_single_point)/='unknown') then
        file_info%frmt=trim(parini%frmt_single_point)
    endif
    do iconf=1,atoms_arr%nconf
        if(trim(potential)/='netsock' .or. iconf==1) then 
            call init_potential_forces(parini,atoms_arr%atoms(iconf))
        endif
        call cal_potential_forces(parini,atoms_arr%atoms(iconf))
        !call atoms_all%fatall(1:3,1:atoms_all%atoms%nat,iconf)=atoms_all%atoms%fat(1:3,1:atoms_all%atoms%nat)
        if(iconf==1) then
            tt1=0.d0
            tt2=0.d0
        else
            tt1=atoms_arr%atoms(iconf)%epot-atoms_arr%atoms(1)%epot
            tt2=atoms_arr%atoms(iconf)%epot-atoms_arr%atoms(iconf-1)%epot
        endif
        write(*,'(a,i7,e24.15,2f10.5)') 'EPOT',iconf,atoms_arr%atoms(iconf)%epot,tt1,tt2
        !if(parini%print_force_single_point) then
        !    do iat=1,atoms_all%atoms%nat
        !        fxyz(1)=atoms_all%atoms%fat(1,iat)
        !        fxyz(2)=atoms_all%atoms%fat(2,iat)
        !        fxyz(3)=atoms_all%atoms%fat(3,iat)
        !        write(*,'(3f12.8)') fxyz(1),fxyz(2),fxyz(3)
        !    enddo
        !endif
        if(trim(potential)/='netsock' .or. iconf==atoms_arr%nconf) then 
            call final_potential_forces(parini,atoms_arr%atoms(iconf))
        endif
        if (iconf==2)  file_info%file_position='append'
        call acf_write(file_info,atoms=atoms_arr%atoms(iconf),strkey='posout') !to be fixed later by atoms_arr
    enddo
    !call atom_all_deallocate(atoms_all,ratall=.true.,fatall=.true.,epotall=.true.)
end subroutine single_point_task
!*****************************************************************************************
