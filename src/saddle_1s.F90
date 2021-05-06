!*****************************************************************************************
subroutine dimer_method(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_atoms_arr, typ_file_info, atom_deallocate_old
    use mod_atoms, only: atom_copy_old, atom_deallocate, set_ndof, atom_calmaxforcecomponent
    use mod_atoms, only: update_ratp, update_rat
    use mod_potential, only: potential, fcalls
    use mod_saddle, only: dmconverged, str_moving_atoms_rand, dimsep, ampl
    use mod_opt, only: typ_paropt
    use mod_yaml_conf, only: write_yaml_conf, read_yaml_conf
    use mod_processors, only: iproc
    use mod_const, only: ang2bohr
    use dynamic_memory
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    type(typ_paropt):: paropt, paropt_m, paropt_m_prec
    type(typ_atoms):: atoms_s
    type(typ_atoms_arr):: atoms_arr
    type(typ_file_info):: file_info
    integer:: iat, ios
    real(8):: epot_t
    real(8):: epot_m0, fmax_m0 !, fnrmtol
    real(8):: fnrm, curv, fmax_s
    real(8), allocatable:: uvn(:,:)
    !logical, allocatable:: atoms_move_random(:)
    character(10):: filename
    character(100):: comment
    !character(256):: comment1, comment2
    call f_routine(id='dimer_method')
    !----------------------------------------------------------------
    potential=trim(parini%potential_potential)
    paropt_m=parini%paropt_geopt
    if(parini%two_level_geopt) then
        paropt_m_prec=parini%paropt_geopt_prec
    endif
    str_moving_atoms_rand=parini%str_moving_atoms_rand_saddle
    dimsep=parini%dimsep_saddle*ang2bohr
    ampl=parini%ampl_saddle*ang2bohr

    paropt=parini%paropt_saddle_opt
    !call acf_read(parini,'posinp.acf',1,atoms=atoms_s)
    call read_yaml_conf(parini,'posinp.yaml',10000,atoms_arr)
    if(atoms_arr%nconf/=1) stop 'ERROR: atoms_arr%nconf/=1 in dimer_method'
    call atom_copy_old(atoms_arr%atoms(1),atoms_s,'atoms_arr%atoms(iconf)->atoms_s')
    call atom_deallocate(atoms_arr%atoms(1))
    deallocate(atoms_arr%atoms)
    call set_ndof(atoms_s)
    !read(comment2,*) atoms_s%boundcond,atoms_s%cellvec(1,1),atoms_s%cellvec(2,2),atoms_s%cellvec(3,3),atoms_s%cellvec(1,2),atoms_s%cellvec(1,3),atoms_s%cellvec(2,3)
    !atoms_s%cellvec(2,1)=0.d0 ; atoms_s%cellvec(3,1)=0.d0 ; atoms_s%cellvec(3,2)=0.d0
    !write(*,*) 'nat ',atoms_s%nat
    uvn=f_malloc([1.to.3,1.to.atoms_s%nat],id='uvn')
    !----------------------------------------------------------------
    !call dm_init(dimsep,degreefreedom,alphax,fnrmtol,maxitec,maxitsd,maxitcg)
    !call dm_init(1.d-3,0,1.d-3,1.d-3,100,2000,500)
    call read_input(atoms_s) !,paropt)
    !call init_potential_forces(cell)
    paropt%approach='SDCG'
    paropt%feedback=2
    paropt%nsatur=1
    paropt%fnrmtolsatur=10.d-1
    paropt%alpha0=2.d0*paropt%alphax
    !call initminimize(paropt)
    call pot_initialize(parini,atoms_s,paropt,paropt_m)
    !---------------------------------------------------------------------------
    call cal_potential_forces(parini,atoms_s)
    epot_m0=atoms_s%epot
    call atom_calmaxforcecomponent(atoms_s%nat,atoms_s%bemoved,atoms_s%fat,fmax_m0)
    !write(*,'(a,es24.15,es13.3)') 'ENERGY: epot_m0,fmax_m0 ',epot_m0,fmax_m0
    call yaml_mapping_open('ENERGY_m0',flow=.true.)
    call yaml_map('epot',epot_m0,fmt='(es24.15)')
    call yaml_map('fmax',fmax_m0,fmt='(es12.3)')
    call yaml_mapping_close()
    file_info%filename_positions='posout.yaml'
    file_info%file_position='new'
    !call acf_write(file_info,atoms=atoms_s,strkey='posout')
    call write_yaml_conf(file_info,atoms_s,strkey='posout')
    !---------------------------------------------------------------------------
    call update_ratp(atoms_s)
    call random_move_atoms(parini,atoms_s%nat,atoms_s%bemoved,atoms_s%cellvec,atoms_s%ratp)
    call update_rat(atoms_s)
    !call cal_potential_forces(parini,iproc,3*atoms_s%nat,atoms_s%rat,atoms_s%fat,epot_m0)
    !call atom_calmaxforcecomponent(atoms_s%nat,atoms_s%bemoved,atoms_s%fat,fmax_m0)
    !write(*,'(a,es24.15,es13.3)') 'ENERGY: epot_m0,fmax_m0 ',epot_m0,fmax_m0
    file_info%filename_positions='posout.yaml'
    file_info%file_position='append'
    !call acf_write(file_info,atoms=atoms_s,strkey='posout')
    call write_yaml_conf(file_info,atoms_s,strkey='posout')
    !stop
    !---------------------------------------------------------------------------
    call update_ratp(atoms_s)
    call dimmethimproved(parini,iproc,atoms_s,atoms_s%nat,atoms_s%ndof,atoms_s%ratp,atoms_s%epot,atoms_s%fat,curv,uvn,paropt)
    call update_rat(atoms_s)
    file_info%filename_positions='posout.yaml'
    file_info%file_position='append'
    !call acf_write(file_info,atoms=atoms_s,strkey='posout')
    call write_yaml_conf(file_info,atoms_s,strkey='posout')

    write(filename,'(a10)') 'modout.xyz'
    write(comment,'(a11,es24.15,a7,i5)') ' angstroem ',curv,' iter= ',0
    call writexyz(filename,'new',atoms_s%nat,uvn,atoms_s%bemoved,atoms_s%sat,atoms_s%cellvec,atoms_s%boundcond,comment)
    !call random_number(uvn)
    !call normalizevector(3*atoms_s%nat,uvn)
    if(dmconverged) then
        !write(*,'(a,f10.3)') 'energy difference between saddle point and intial minimum ', &
        !    atoms_s%epot-epot_m0
        call atom_calmaxforcecomponent(atoms_s%nat,atoms_s%bemoved,atoms_s%fat,fmax_s)
        !write(*,'(a,es24.15,es13.3)') 'ENERGY: epot_s,fmax_s   ',atoms_s%epot,fmax_s
        call yaml_mapping_open('ENERGY_s',flow=.true.)
        call yaml_map('epot',atoms_s%epot,fmt='(es24.15)')
        call yaml_map('fmax',fmax_s,fmt='(es12.3)')
        call yaml_mapping_close()
        !paropt_m%alphax=paropt%alphax
        !paropt_m%fmaxtol=paropt%fmaxtol
        !paropt_m%nit=1000
        !paropt_m%dt_start=5.d-3
        !paropt_m%dtmax=10.d0*paropt_m%dt_start
        !paropt_m%approach='FIRE'
        !paropt_m%lprint=.true.
        call find_minima(parini,iproc,atoms_s,paropt_m,paropt_m_prec,uvn,curv,epot_m0)
    endif
    call yaml_map('force_call',int(fcalls))
    !write(*,'(a,i7)') 'force_call',int(fcalls)
    call f_free(uvn) !,atoms_move_random)
    call atom_deallocate_old(atoms_s,sat=.true.,rat=.true.,fat=.true.,bemoved=.true.)
    !call deallocateatomsarrays
    call f_release_routine()
end subroutine dimer_method
!*****************************************************************************************
subroutine read_input(atoms_s) !,paropt)
    use mod_saddle, only: dimsep, ampl, do_elim_trans, do_elim_rot, &
        maxitsd, maxitcg, sdconverged, cgconverged, &
        nmatr, moving_atoms_rand, str_moving_atoms_rand
    use mod_atoms, only: typ_atoms
    use yaml_output
    !use mod_opt, only: typ_paropt
    implicit none
    type(typ_atoms), intent(inout):: atoms_s
    !type(typ_paropt), intent(inout):: paropt
    !integer, intent(out)::
    !local variables
    integer:: ios, iat
    !character(256):: strline
    !open(unit=1,file='input.dimer',status='old',iostat=ios)
    !read(1,*) ampl
    !read(1,*) dimsep
    !!read(1,*) fnrmtol
    !read(1,*) paropt%fmaxtol
    !!read(1,*) alphax
    !read(1,*) paropt%alphax
    !!read(1,*) maxitsd
    !read(1,*) paropt%nit
    !!read(1,*) maxitcg
    !!read(1,*) paropt%nitcg
    !read(1,'(a)') strline
    !close(1)
    read(str_moving_atoms_rand,*) nmatr
    if(nmatr==0 .or. nmatr>10) stop 'ERROR: nmatr=0 or nmatr>10'
    !write(*,'(a)') trim(str_moving_atoms_rand)
    read(str_moving_atoms_rand,*) nmatr,moving_atoms_rand(1:nmatr)
    call yaml_map('moving_atoms_rand',moving_atoms_rand(1:nmatr))
    !write(*,'(a)',advance='no') 'moving_atoms_rand '
    !do iat=1,nmatr
    !    write(*,'(i5)',advance='no') moving_atoms_rand(iat)
    !enddo
    !write(*,*)
    sdconverged=.false.
    cgconverged=.false.
    if(atoms_s%ndof/=3*atoms_s%nat) then
        do_elim_trans=.false.
        do_elim_rot=.false.
    elseif(trim(atoms_s%boundcond)=='free') then
        do_elim_trans=.true.
        do_elim_rot=.true.
    elseif(trim(atoms_s%boundcond)=='bulk' .and. trim(atoms_s%boundcond)=='slab') then
        do_elim_trans=.true.
        do_elim_rot=.false.
    else
        stop 'ERROR: does not know what to do with translation and rotation.'
    endif
end subroutine read_input
!*****************************************************************************************
subroutine random_move_atoms(parini,nat,atom_motion,cellvec,rat)
    use mod_parini, only: typ_parini
    use mod_saddle, only: ampl, moving_atoms_rand
    use mod_const, only: ang2bohr
    use mod_utils
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: nat
    logical, intent(in):: atom_motion(3,nat)
    real(8), intent(in):: cellvec(3,3)
    real(8), intent(inout):: rat(3,nat)
    !local variables
    integer:: iat, ix, iy, iz
    real(8):: x, y, z, dx, dy, dz, r, rt, trand(3), rc, weight
    character(10):: filename
    character(100):: comment
    rc=4.0d0*ang2bohr
    do iat=1,nat
        r=1.d20
        do iz=-1,1
        do iy=-1,1
        do ix=-1,1
        x=rat(1,iat)+ix*cellvec(1,1)+iy*cellvec(1,2)+iz*cellvec(1,3)
        y=rat(2,iat)+ix*cellvec(2,1)+iy*cellvec(2,2)+iz*cellvec(2,3)
        z=rat(3,iat)+ix*cellvec(3,1)+iy*cellvec(3,2)+iz*cellvec(3,3)
        dx=x-rat(1,moving_atoms_rand(1))
        dy=y-rat(2,moving_atoms_rand(1))
        dz=z-rat(3,moving_atoms_rand(1))
        rt=sqrt(dx**2+dy**2+dz**2)
        if(rt<r) r=rt
        enddo
        enddo
        enddo
        !if(r>2.0d0) cycle 
        if(r<rc) then
            weight=ampl*(1.d0-(r/rc)**2)**4
            if(trim(parini%rng_type)=='only_for_tests') then
                call random_number_generator_simple(3,trand)
            else
                call random_number(trand)
            endif
            trand(1:3)=(trand(1:3)-0.5d0)*2.d0
        else
            !weight=0.d0
            cycle
        endif
        !write(*,'(i5,f10.5)') iat,weight
        if(atom_motion(1,iat)) rat(1,iat)=rat(1,iat)+weight*trand(1)
        if(atom_motion(2,iat)) rat(2,iat)=rat(2,iat)+weight*trand(2)
        if(atom_motion(3,iat)) rat(3,iat)=rat(3,iat)+weight*trand(3)
    enddo
end subroutine random_move_atoms
!*****************************************************************************************
subroutine find_minima(parini,iproc,atoms_s,paropt_m,paropt_m_prec,uvn,curv,epot_m0)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_file_info, atom_deallocate_old, atom_allocate_old
    use mod_atoms, only: atom_calmaxforcecomponent
    use mod_atoms, only: atom_copy_old, get_rat, update_rat
    use mod_yaml_conf, only: write_yaml_conf
    use mod_opt, only: typ_paropt
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iproc !, nat, ndof
    type(typ_atoms), intent(in):: atoms_s
    type(typ_paropt), intent(inout):: paropt_m, paropt_m_prec
    real(8), intent(in):: uvn(3,atoms_s%nat), curv, epot_m0
    !local variables
    type(typ_atoms):: atoms_m1, atoms_m2
    type(typ_file_info):: file_info
    integer:: iat
    real(8):: zeta, fmax1, fmax2
    !real(8), allocatable:: rat1(:,:), rat2(:,:)
    !real(8), allocatable:: fat1(:,:), fat2(:,:)
    character(100):: comment
    !allocate(rat1(3,atoms_s%nat),rat2(3,atoms_s%nat))
    !allocate(fat1(3,atoms_s%nat),fat2(3,atoms_s%nat))
    call atom_allocate_old(atoms_m1,atoms_s%nat,0,0,sat=.true.,fat=.true.,bemoved=.true.)
    call atom_allocate_old(atoms_m2,atoms_s%nat,0,0,sat=.true.,fat=.true.,bemoved=.true.)
    call atom_copy_old(atoms_s,atoms_m1,'atoms_s->atoms_m1')
    call atom_copy_old(atoms_s,atoms_m2,'atoms_s->atoms_m2')
    !atoms_m1%rat(1:3,1:atoms_s%nat)=atoms_s%rat(1:3,1:atoms_s%nat)
    !atoms_m2%rat(1:3,1:atoms_s%nat)=atoms_s%rat(1:3,1:atoms_s%nat)
    zeta=min(0.2d0,0.1d0/abs(curv))
    call get_rat(atoms_s,atoms_m1%ratp)
    do iat=1,atoms_s%nat
        if(atoms_s%bemoved(1,iat)) atoms_m1%ratp(1,iat)=atoms_m1%ratp(1,iat)+zeta*uvn(1,iat)
        if(atoms_s%bemoved(2,iat)) atoms_m1%ratp(2,iat)=atoms_m1%ratp(2,iat)+zeta*uvn(2,iat)
        if(atoms_s%bemoved(3,iat)) atoms_m1%ratp(3,iat)=atoms_m1%ratp(3,iat)+zeta*uvn(3,iat)
    enddo
    call update_rat(atoms_m1)
    if(parini%two_level_geopt) then
        call minimize(parini,iproc,atoms_m1,paropt_m_prec)
    endif
    call minimize(parini,iproc,atoms_m1,paropt_m)
    call atom_calmaxforcecomponent(atoms_s%nat,atoms_s%bemoved,atoms_m1%fat,fmax1)
    !if(paropt_m%converged) write(*,'(a,es24.15,es13.3)') 'ENERGY: epot1,fmax1     ',atoms_m1%epot,fmax1
    if(paropt_m%converged) then
        call yaml_mapping_open('ENERGY_m1',flow=.true.)
        call yaml_map('epot',atoms_m1%epot,fmt='(es24.15)')
        call yaml_map('fmax',fmax1,fmt='(es12.3)')
        call yaml_mapping_close()
    endif
    file_info%filename_positions='posmin1.yaml'
    file_info%file_position='new'
    !call acf_write(file_info,atoms=atoms_m1,strkey='minimum1')
    call write_yaml_conf(file_info,atoms_m1,strkey='minimum1')
    call get_rat(atoms_s,atoms_m2%ratp)
    do iat=1,atoms_s%nat
        if(atoms_s%bemoved(1,iat)) atoms_m2%ratp(1,iat)=atoms_m2%ratp(1,iat)-zeta*uvn(1,iat)
        if(atoms_s%bemoved(2,iat)) atoms_m2%ratp(2,iat)=atoms_m2%ratp(2,iat)-zeta*uvn(2,iat)
        if(atoms_s%bemoved(3,iat)) atoms_m2%ratp(3,iat)=atoms_m2%ratp(3,iat)-zeta*uvn(3,iat)
    enddo
    call update_rat(atoms_m2)
    call atom_calmaxforcecomponent(atoms_s%nat,atoms_s%bemoved,atoms_m2%fat,fmax2)
    if(parini%two_level_geopt) then
        call minimize(parini,iproc,atoms_m2,paropt_m_prec)
    endif
    call minimize(parini,iproc,atoms_m2,paropt_m)
    !if(paropt_m%converged) write(*,'(a,es24.15,es13.3)') 'ENERGY: epot2,fmax2     ',atoms_m2%epot,fmax2
    if(paropt_m%converged) then
        call yaml_mapping_open('ENERGY_m2',flow=.true.)
        call yaml_map('epot',atoms_m2%epot,fmt='(es24.15)')
        call yaml_map('fmax',fmax2,fmt='(es12.3)')
        call yaml_mapping_close()
    endif
    file_info%filename_positions='posmin2.yaml'
    file_info%file_position='new'
    !call acf_write(file_info,atoms=atoms_m2,strkey='minimum2')
    call write_yaml_conf(file_info,atoms_m2,strkey='minimum2')
    !write(*,'(a,4f13.3)') 'bh0,curv,bh1,bh2 ', &
    !    atoms_s%epot-epot_m0,curv,atoms_s%epot-atoms_m1%epot,atoms_s%epot-atoms_m2%epot
    call yaml_mapping_open('barriers',flow=.true.)
    call yaml_map('bh0',atoms_s%epot-epot_m0,fmt='(f13.3)')
    call yaml_map('curv',curv,fmt='(f13.3)')
    call yaml_map('bh1',atoms_s%epot-atoms_m1%epot,fmt='(f13.3)')
    call yaml_map('bh2',atoms_s%epot-atoms_m2%epot,fmt='(f13.3)')
    call yaml_mapping_close()

    call atom_deallocate_old(atoms_m1,sat=.true.,rat=.true.,fat=.true.,bemoved=.true.)
    call atom_deallocate_old(atoms_m2,sat=.true.,rat=.true.,fat=.true.,bemoved=.true.)
    !deallocate(rat1,rat2)
    !deallocate(fat1,fat2)
end subroutine find_minima
!*****************************************************************************************
subroutine alongnegcurvature(iproc,atoms,uvn,c)
    use mod_atoms, only: typ_atoms, atom_deallocate_old
    use mod_saddle, only: dimsep
    use mod_atoms, only: atom_copy_old, update_rat, update_ratp
    implicit none
    integer, intent(in):: iproc
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(in):: uvn(3,atoms%nat) !unit vector along the dimer, \hat{n}.
    real(8), intent(inout):: c !curvature
    !local variables
    type(typ_atoms):: atoms_t
    integer:: istat, i, ncount
    real(8):: pi, fnrm, c0, t1, t2, t3, DDOT, trand(10)
    stop 'ERROR: check this routine before use'
    !if(iproc==0) write(*,'(a,f15.5)') 'dimsep ',dimsep
    !ncount=icount
    pi=4.d0*atan(1.d0)
    call atom_copy_old(atoms,atoms_t,'atoms->atoms_t')
    !call cal_potential_forces(parini,iproc,n,x,f,epotprime)
    !-------------------------------
    !do i=-5,5
    !    alpha=i*1.d-3
    !    rat1(1:3,natr+1:nat)=rat(1:3,natr+1:nat)
    !    rat1(1:3,1:natr)=rat(1:3,1:natr)+alpha*uvn(1:3,1:natr)
    !    call cal_potential_forces(parini,iproc,3*nat,rat1,fat,epot1)
    !    write(91,*) alpha,epot1
    !enddo
    !call random_number(trand)
    !call random_number(uvn)
    !write(*,*) 'do_elim_trans,do_elim_rot ',do_elim_trans,do_elim_rot
    !if(do_elim_trans) call projtransout(nr,uvn)
    !if(do_elim_rot) call projrotout(n,nr,x,uvn)
    !call normalizevector(nr,uvn)
    !call lowestcurvature(iproc,n/3,nr,x,uvn,f,2.d0,10,c0,c,1)
    call update_ratp(atoms_t)
    do i=-5,5
        atoms_t%ratp(1:3,1:atoms%nat)=atoms%ratp(1:3,1:atoms%nat)+i*4.d-2*uvn(1:3,1:atoms%nat)
        call update_rat(atoms_t,upall=.true.)
        stop 'ERROR: line was commented due to missing parini in this routine'
        !call cal_potential_forces(parini,atoms_t)
        call calnorm(3*atoms_t%nat,atoms_t%fat,fnrm)
        write(31,*) i*4.d-2,fnrm,atoms_t%epot
    enddo
    call atom_deallocate_old(atoms_t)
    !-------------------------------
end subroutine alongnegcurvature
!*****************************************************************************************
