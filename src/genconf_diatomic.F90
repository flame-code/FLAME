!*****************************************************************************************
subroutine genconf_diatomic(parini,genconf)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_all, typ_file_info, atom_all_allocate, atom_all_deallocate
    use mod_atoms, only: atom_allocate_old, atom_deallocate_old, set_rat
    use mod_genconf, only: typ_genconf
    use mod_processors, only: iproc
    use mod_potential, only: potential
    use mod_const, only: bohr2ang
    use mod_acf, only: acf_write
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_genconf), intent(in):: genconf
    !local variables
    type(typ_atoms_all):: atoms_all
    type(typ_file_info):: file_info
    integer:: iconf, iat
    real(8):: dx, r, cv(3,3)
    call atom_allocate_old(atoms_all%atoms,2,0,0)
    atoms_all%nconfmax=parini%npoint_genconf
    call atom_all_allocate(atoms_all,ratall=.true.,epotall=.true.)
    atoms_all%atoms%sat(1)='Na'
    atoms_all%atoms%sat(2)='Cl'
    !call set_rcov(atoms_all%atoms)
    atoms_all%atoms%boundcond='free'
    !r0=atoms_all%atoms%rcov(1)+atoms_all%atoms%rcov(2)
    !setting the cell dimensions
    cv(1:3,1:3)=0.d0
    cv(1,1)=11.d0+parini%dmax_genconf
    cv(2,2)=11.d0
    cv(3,3)=11.d0
    atoms_all%atoms%cellvec(1:3,1:3)=cv(1:3,1:3)/bohr2ang
    if(trim(genconf%cal_pot)=='yes') then
        potential=trim(parini%potential_potential)
        call init_potential_forces(parini,atoms_all%atoms)
    endif

    !write(*,'(a,i6)') 'nat ',atoms%nat
    !dx=2.5d0*r0/real(atoms_all%nconfmax-1,8)
    dx=(parini%dmax_genconf-parini%dmin_genconf)/real(atoms_all%nconfmax-1,8)
    do iconf=1,atoms_all%nconfmax
        atoms_all%ratall(1:3,1,iconf)=4.d0/bohr2ang
        !atoms_all%ratall(1,2,iconf)=4.d0+r0+i*dx
        atoms_all%ratall(1,2,iconf)=(4.d0+parini%dmin_genconf+(iconf-1)*dx)/bohr2ang
        atoms_all%ratall(2:3,2,iconf)=4.d0/bohr2ang
        if(trim(genconf%cal_pot)=='yes') then
            call set_rat(atoms_all%atoms,atoms_all%ratall(1,1,iconf),setall=.true.)
            call cal_potential_forces(parini,atoms_all%atoms)
            atoms_all%epotall(iconf)=atoms_all%atoms%epot
            r=atoms_all%ratall(1,2,iconf)-atoms_all%ratall(1,1,iconf)
            write(21,'(i5,f10.4,es24.15)') iconf,r*bohr2ang,atoms_all%epotall(iconf)
        endif
    enddo
    atoms_all%nconf=atoms_all%nconfmax
    !writing all confs. to posout.acf
    file_info%filename_positions='posout.acf'
    file_info%file_position='new'
    call acf_write(file_info,atoms_all=atoms_all,strkey='posout')
    !finalizing
    if(trim(genconf%cal_pot)=='yes') then
        call final_potential_forces(parini,atoms_all%atoms)
    endif
    call atom_deallocate_old(atoms_all%atoms)
    call atom_all_deallocate(atoms_all,ratall=.true.,epotall=.true.)
end subroutine genconf_diatomic
!*****************************************************************************************
