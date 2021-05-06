!*****************************************************************************************
subroutine genconf_trimer(parini,genconf)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_all, typ_file_info, atom_all_allocate, atom_all_deallocate
    use mod_atoms, only: atom_allocate_old, atom_deallocate_old, set_rcov
    use mod_genconf, only: typ_genconf
    use mod_processors, only: iproc
    use mod_potential, only: potential
    use mod_acf, only: acf_write
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_genconf), intent(in):: genconf
    !local variables
    type(typ_atoms_all):: atoms_all
    type(typ_file_info):: file_info
    integer:: nrho, nphi, irho2, irho3, iphi, iconf, irho2_inc, irho3_inc
    real(8):: r0, r12, r13, r23, drho, dphi, pi, dx, dy
    real(8):: de
    logical:: conf_repeat
    pi=4.d0*atan(1.d0)
    call atom_allocate_old(atoms_all%atoms,3,0,0)
    nrho=20
    nphi=20
    atoms_all%nconfmax=1000
    call atom_all_allocate(atoms_all,ratall=.true.,epotall=.true.)
    atoms_all%atoms%sat(1)='Au'
    atoms_all%atoms%sat(2)='Au'
    atoms_all%atoms%sat(3)='Au'
    call set_rcov(atoms_all%atoms)
    !atoms_all%atoms%rcov(1:2)=1.d0
    r0=atoms_all%atoms%rcov(1)+atoms_all%atoms%rcov(2)
    !setting the cell dimensions
    atoms_all%atoms%cellvec(1:3,1:3)=0.d0
    atoms_all%atoms%cellvec(1,1)=20.d0
    atoms_all%atoms%cellvec(2,2)=20.d0
    atoms_all%atoms%cellvec(3,3)=10.d0
    !initialize potential if it must be calculated
    if(trim(genconf%cal_pot)=='yes') then
        potential=trim(parini%potential_potential)
        call init_potential_forces(parini,atoms_all%atoms)
    endif
    !write(*,'(a,i6)') 'nat ',atoms%nat
    drho=3.d0*r0/real(nrho,8)
    dphi=pi/real(nphi,8)
    atoms_all%nconf=0
    irho2_inc=1
    irho3_inc=1
    irho2=0
    irho3=0
    do iphi=1,nphi-1
        do
            do
                !write(*,*) iphi,irho2,irho3,atoms_all%nconf
                !first atom
                atoms_all%atoms%ratp(1,1)=6.d0
                atoms_all%atoms%ratp(2,1)=6.d0
                atoms_all%atoms%ratp(3,1)=6.d0
                !second atoms
                atoms_all%atoms%ratp(1,2)=6.d0+irho2*drho
                atoms_all%atoms%ratp(2,2)=6.d0
                atoms_all%atoms%ratp(3,2)=6.d0
                !third atom
                atoms_all%atoms%ratp(1,3)=6.d0+irho3*drho*cos(iphi*dphi)
                atoms_all%atoms%ratp(2,3)=6.d0+irho3*drho*sin(iphi*dphi)
                atoms_all%atoms%ratp(3,3)=6.d0
                !calculating distances
                r12=irho2*drho
                r13=irho3*drho
                dx=atoms_all%atoms%ratp(1,3)-atoms_all%atoms%ratp(1,2)
                dy=atoms_all%atoms%ratp(2,3)-atoms_all%atoms%ratp(2,2)
                r23=sqrt(dx**2+dy**2)
                if(r12<0.7d0*r0 .or. r13<0.7d0*r0 .or. r23<0.7d0*r0) goto 1000
                if(r12>2.0d0*r0 .or. r13>2.0d0*r0 .or. r23>2.0d0*r0) goto 1000
                atoms_all%nconf=atoms_all%nconf+1
                if(atoms_all%nconf>atoms_all%nconfmax) then
                    stop 'ERROR: too many configuration in genconf_trimer'
                endif
                if(trim(genconf%cal_pot)=='yes') then
                    call cal_potential_forces(parini,atoms_all%atoms)
                    conf_repeat=.false.
                    do iconf=1,atoms_all%nconf-1
                        de=abs(atoms_all%atoms%epot-atoms_all%epotall(iconf))
                        if(de<1.d-14) then
                            write(*,'(a,es15.5,3f10.4)') 'WARNING: two confs. having tiny energy difference: ',de,r12,r13,r23
                            conf_repeat=.true.
                            atoms_all%nconf=atoms_all%nconf-1
                            exit
                        endif
                    enddo
                    if(conf_repeat) goto 1000
                    atoms_all%epotall(atoms_all%nconf)=atoms_all%atoms%epot
                    write(21,'(i5,3f10.4,es24.15)') atoms_all%nconf,r12,r13,r23,atoms_all%atoms%epot
                endif
                atoms_all%ratall(1:3,1:3,atoms_all%nconf)=atoms_all%atoms%ratp(1:3,1:3)
1000            continue
                !if(irho2==0) exit
                if((irho3==irho2 .and. irho3_inc==1) .or. (irho3==0 .and. irho3_inc==-1)) then
                    irho3_inc=-1*irho3_inc
                    exit
                else
                    irho3=irho3+irho3_inc
                endif
            enddo !end of loop over irho3
            if((irho2==nrho .and. irho2_inc==1) .or. (irho2==0 .and. irho2_inc==-1)) then
                irho2_inc=-1*irho2_inc
                exit
            else
                irho2=irho2+irho2_inc
            endif
        enddo !end of loop over irho2
    enddo
    !writing all confs. to posout.acf
    atoms_all%atoms%boundcond='free'
    file_info%filename_positions='posout.acf'
    file_info%file_position='new'
    call acf_write(file_info,atoms_all=atoms_all,strkey='posout')
    !finalizing
    if(trim(genconf%cal_pot)=='yes') then
        call final_potential_forces(parini,atoms_all%atoms)
    endif
    call atom_deallocate_old(atoms_all%atoms)
    call atom_all_deallocate(atoms_all,ratall=.true.,epotall=.true.)
end subroutine genconf_trimer
!*****************************************************************************************
