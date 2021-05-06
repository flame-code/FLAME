!*****************************************************************************************
subroutine cal_hessian_4p(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_atoms_arr, typ_file_info, atom_copy_old
    use mod_atoms, only: atom_copy, atom_deallocate, set_atomic_mass
    use mod_atoms, only: update_ratp, update_rat, set_rat, set_rat_iat, get_rat
    use mod_processors, only: iproc
    use mod_potential, only: potential
    use mod_yaml_conf, only: write_yaml_conf, read_yaml_conf
    use futile
    use dynamic_memory
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    type(typ_atoms):: atoms
    type(typ_atoms_arr):: atoms_arr
    type(typ_file_info):: file_info
    real(8), allocatable:: hess(:,:), freq(:)
    real(8), allocatable:: rat_center(:), eval(:), work(:)
    real(8):: h, rlarge, twelfth, twothird, shift, dm, tt, s, alpha, ttmass, ttnorm
    real(8):: xyz(3)
    integer:: istat, lwork, info, jj !, nat_yes
    integer:: i, iat, ixyz, j, jat, jxyz, imode, ii, iconf, iunit
    character(5):: fn
    character(10):: strkey
    !call acf_read(parini,'posinp.acf',1,atoms=atoms)
    call read_yaml_conf(parini,'posinp.yaml',10000,atoms_arr)
    if(atoms_arr%nconf/=1) stop 'ERROR: atoms_arr%nconf/=1 in cal_hessian_4p'
    call atom_copy(atoms_arr%atoms(1),atoms,'atoms_arr%atoms(1)->atoms')
    call atom_deallocate(atoms_arr%atoms(1))
    deallocate(atoms_arr%atoms)
    potential=trim(parini%potential_potential)
    !logical, allocatable:: yes(:)
    allocate(hess(3*atoms%nat,3*atoms%nat),stat=istat)
    allocate(rat_center(3*atoms%nat),stat=istat)
    allocate(eval(3*atoms%nat),stat=istat)
    lwork=100*atoms%nat
    allocate(work(lwork),stat=istat)
    !h=1.d-1
    !h=7.5d-2
    !h=5.d-2
    h=4.d-2
    !h=2.d-2
    rlarge=1.d0*1.d6
    twelfth=-1.d0/(12.d0*h)
    twothird=-2.d0/(3.d0*h)
    !-------------------------------------------------------
    call init_potential_forces(parini,atoms)
    call get_rat(atoms,rat_center)
    ii=1
    do i=ii,3*atoms%nat
        iat=(i-1)/3+1
        ixyz=mod(i-1,3)+1
        !if(.not. yes(iat)) cycle
        call set_rat(atoms,rat_center,setall=.true.)
        !-----------------------------------------
        call update_ratp(atoms)
        xyz=atoms%ratp(1:3,iat)
        xyz(ixyz)=atoms%ratp(ixyz,iat)-2*h
        call set_rat_iat(atoms,iat,xyz)
        call cal_potential_forces(parini,atoms)
        do j=1,3*atoms%nat
            jat=(j-1)/3+1
            jxyz=mod(j-1,3)+1
            hess(j,i)=twelfth*atoms%fat(jxyz,jat)
        enddo
        !-----------------------------------------
        call update_ratp(atoms)
        xyz=atoms%ratp(1:3,iat)
        xyz(ixyz)=atoms%ratp(ixyz,iat)+h
        call set_rat_iat(atoms,iat,xyz)
        call cal_potential_forces(parini,atoms)
        do j=1,3*atoms%nat
            jat=(j-1)/3+1
            jxyz=mod(j-1,3)+1
            hess(j,i)=hess(j,i)-twothird*atoms%fat(jxyz,jat)
        enddo
        !-----------------------------------------
        call update_ratp(atoms)
        xyz=atoms%ratp(1:3,iat)
        xyz(ixyz)=atoms%ratp(ixyz,iat)+2*h
        call set_rat_iat(atoms,iat,xyz)
        call cal_potential_forces(parini,atoms)
        do j=1,3*atoms%nat
            jat=(j-1)/3+1
            jxyz=mod(j-1,3)+1
            hess(j,i)=hess(j,i)+twothird*atoms%fat(jxyz,jat)
        enddo
        !-----------------------------------------
        call update_ratp(atoms)
        xyz=atoms%ratp(1:3,iat)
        xyz(ixyz)=atoms%ratp(ixyz,iat)+h
        call set_rat_iat(atoms,iat,xyz)
        call cal_potential_forces(parini,atoms)
        do j=1,3*atoms%nat
            jat=(j-1)/3+1
            jxyz=mod(j-1,3)+1
            hess(j,i)=hess(j,i)-twelfth*atoms%fat(jxyz,jat)
        enddo
        !-----------------------------------------
        !if(iproc==0) then
        !    open(unit=1388,file='restart.dat',status='replace')
        !    write(1388,*) i
        !    do j=1,3*atoms%nat
        !    do k=1,3*atoms%nat
        !        write(1388,'(1es25.16)') hess(k,j)
        !    enddo
        !    enddo
        !    close(1388)
        !endif
        !-----------------------------------------
    enddo
    call set_rat(atoms,rat_center,setall=.true.)
    call final_potential_forces(parini,atoms)
    !-------------------------------------------------------
    !deallocate(yes)
    !check symmetry
    dm=0.d0
    do i=1,3*atoms%nat
        do j=1,i-1
            s=.5d0*(hess(i,j)+hess(j,i))
            tt=abs(hess(i,j)-hess(j,i))/(1.d0+abs(s))
            dm=max(dm,tt)
            hess(i,j)=s
            hess(j,i)=s
        enddo
    enddo
    if(dm>1.d-1) write(*,*) 'max dev from sym',dm
    !do j=1,3*atoms%nat
    !    do i=1,3*atoms%nat
    !        write(201,*) hess(i,j)
    !    enddo
    !enddo
    !-------------------------------------------------------
    !project out rotations
    if(trim(atoms%boundcond)=='free') then
        call projectout_rotation(atoms,hess,rlarge)
    endif
    !-------------------------------------------------------
    !project out translations
    shift=rlarge/atoms%nat
    do j=1,3*atoms%nat-2,3
        do i=1,3*atoms%nat-2,3
            hess(i+0,j+0)=hess(i+0,j+0)+shift
            hess(i+1,j+1)=hess(i+1,j+1)+shift
            hess(i+2,j+2)=hess(i+2,j+2)+shift
        enddo
    enddo
    !-------------------------------------------------------
    !do j=1,3*atoms%nat
    !    do i=1,3*atoms%nat
    !        write(202,*) hess(i,j)
    !    enddo
    !enddo
    !-------------------------------------------------------
    !check
    call set_atomic_mass(atoms)
    do j=1,3*atoms%nat
        jat=(j-1)/3+1
        ttmass=1.d0/sqrt(atoms%amass(jat))
        do i=1,3*atoms%nat
            hess(i,j)=hess(i,j)*ttmass
        enddo
        do i=1,3*atoms%nat
            hess(j,i)=hess(j,i)*ttmass
        enddo
    enddo
    call DSYEV('V','L',3*atoms%nat,hess,3*atoms%nat,eval,work,lwork,info)
    do j=1,3*atoms%nat
        if(hess(1,j)<0.d0) then
            do i=1,3*atoms%nat
                hess(i,j)=-hess(i,j)
            enddo
        endif
    enddo
    do j=1,3*atoms%nat
        jat=(j-1)/3+1
        jxyz=mod(j-1,3)+1
        ttmass=1.d0/sqrt(atoms%amass(jat))
        ttnorm=0.d0
        do i=1,3*atoms%nat
            hess(i,j)=hess(i,j)*ttmass
            ttnorm=ttnorm+hess(i,j)**2
        enddo
        ttnorm=sqrt(ttnorm)
        do i=1,3*atoms%nat
            hess(i,j)=hess(i,j)/ttnorm
        enddo
    enddo
    if(info/=0) stop 'DSYEV'
    if(iproc==0) then
        freq=f_malloc0([1.to.3*atoms%nat],id='freq')
        !iunit=f_get_free_unit(10**5)
        !open(unit=iunit,file='phonons.dat',status='replace')
        do i=1,3*atoms%nat
            if(eval(i)<0.d0) then
                freq(i)=-sqrt(-eval(i))*219474.63068d0
            else
                freq(i)=sqrt(eval(i))*219474.63068d0
            endif
        enddo
        call yaml_map('vibrational frequencies',freq,fmt='(e23.15)')
        !call yaml_map('vibrational frequencies',sqrt(eval)*219474.63068d0,fmt='(e23.15)')
        !write(iunit,*) '---  TB eigenvalues in a.u. -------------'
        atoms_arr%nconf=2*10+1
        allocate(atoms_arr%atoms(atoms_arr%nconf))
        call yaml_sequence_open('vibrational eigenmodes')
        call update_ratp(atoms)
        do imode=1,3*atoms%nat
            call yaml_sequence(advance='no')
            write(fn,'(i5.5)') imode
            strkey='mode_'//fn
            call yaml_map(strkey,hess(1:3*atoms%nat,imode),fmt='(e23.15)')
            file_info%filename_positions='mode_'//fn//'.yaml'
            file_info%file_position='new'
            file_info%print_force=.false.
            !tt=eval(imode)/(.529d0**2/27.2114d0)
            !write(iunit,'(a,i6,2e15.5)') 'eval (eV/A^2), a.u. ',imode,tt,eval(imode)
            !eval(imode)=tt
            !write(iunit,'(i6,e15.5)') imode,eval(imode)
            if(trim(atoms%boundcond)=='free' .and. imode>3*atoms%nat-6) cycle
            if(trim(atoms%boundcond)/='free' .and. imode>3*atoms%nat-3) cycle
            iconf=0
            do j=-10,10
                iconf=iconf+1
                call atom_copy_old(atoms,atoms_arr%atoms(iconf),'atoms->atoms_arr%atoms(iconf)')
                alpha=j*5.d-2 !1.d-4/eval(imode)
                do iat=1,atoms%nat
                    atoms_arr%atoms(iconf)%ratp(1,iat)=atoms%ratp(1,iat)+alpha*hess(3*iat-2,imode)
                    atoms_arr%atoms(iconf)%ratp(2,iat)=atoms%ratp(2,iat)+alpha*hess(3*iat-1,imode)
                    atoms_arr%atoms(iconf)%ratp(3,iat)=atoms%ratp(3,iat)+alpha*hess(3*iat-0,imode)
                enddo
                call update_rat(atoms_arr%atoms(iconf),upall=.true.)
            enddo
            !call acf_write_new(file_info,atoms_arr=atoms_arr,strkey='mode'//fn)
            do iconf=1,atoms_arr%nconf
                if(iconf==1) then
                    file_info%file_position='new'
                else
                    file_info%file_position='append'
                endif
                call write_yaml_conf(file_info,atoms_arr%atoms(iconf),strkey='mode'//fn)
            enddo
        enddo
        call yaml_sequence_close()
        deallocate(atoms_arr%atoms)
        call f_free(freq)
    endif
    call atom_deallocate(atoms)
end subroutine cal_hessian_4p
!*****************************************************************************************
subroutine projectout_rotation(atoms,hess,rlarge)
    use mod_atoms, only: typ_atoms, update_ratp
    use dynamic_memory
    implicit none
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(inout):: hess(3*atoms%nat,3*atoms%nat)
    real(8), intent(in):: rlarge
    !local variables
    real(8):: cmx, cmy, cmz
    real(8):: t1, t2, t3
    real(8), allocatable:: vrot(:,:)
    integer:: ixyz, iat, i, j
    vrot=f_malloc0([1.to.3*atoms%nat,1.to.3],id='vrot')
    call update_ratp(atoms)
    call calc_rotation_eigenvectors(atoms%nat,atoms%ratp,vrot)
    do j=1,3*atoms%nat
        do i=1,3*atoms%nat
            hess(i,j)=hess(i,j)+(rlarge)*vrot(i,1)*vrot(j,1)
            hess(i,j)=hess(i,j)+(rlarge)*vrot(i,2)*vrot(j,2)
            hess(i,j)=hess(i,j)+(rlarge)*vrot(i,3)*vrot(j,3)
        enddo
    enddo
    call f_free(vrot)
end subroutine projectout_rotation
!*****************************************************************************************
