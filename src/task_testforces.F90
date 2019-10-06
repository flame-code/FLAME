!*****************************************************************************************
subroutine task_testforces(parini)
    use mod_parini, only: typ_parini
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    if(trim(parini%testforces_method)=='finite_difference') then
        call testforces_fd(parini)
    else if(trim(parini%testforces_method)=='stress_fd') then
        call teststress_fd(parini)
    else if(trim(parini%testforces_method)=='stress_fd_cellvec') then
        call teststress_fd_cellvec(parini)
    else
        write(*,'(2a)') 'ERROR: unknown parini%testforces_method: ',trim(parini%testforces_method)
        stop
    endif
end subroutine task_testforces
!*****************************************************************************************
subroutine testforces_fd(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_copy_old, atom_deallocate
    use mod_atoms, only: get_rat_iat, set_rat_iat
    use mod_potential, only: potential
    use mod_processors, only: iproc
    use mod_const, only: bohr2ang
    use mod_acf, only: acf_read
    use mod_yaml_conf, only: read_yaml_conf
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    type(typ_atoms):: atoms , atoms_center
    type(typ_atoms_arr):: atoms_arr
    integer:: iat
    real(8):: h, fnrm, fd, epot_r, epot_l, fsum(3), xyz(3)
    h=1.d-5/bohr2ang
    !call acf_read(parini,'posinp.acf',1,atoms=atoms)
    call read_yaml_conf(parini,'posinp.yaml',1,atoms_arr)
    if(atoms_arr%nconf/=1) stop 'ERROR: atoms_arr%nconf/=1 in testforces_fd'
    call atom_copy_old(atoms_arr%atoms(1),atoms,'atoms_arr%atoms(iconf)->atoms_s')
    call atom_deallocate(atoms_arr%atoms(1))
    deallocate(atoms_arr%atoms)
    potential=trim(parini%potential_potential)
    call init_potential_forces(parini,atoms)
    call cal_potential_forces(parini,atoms)
    call calnorm(3*atoms%nat,atoms%fat,fnrm)
    call yaml_map('finite difference step size',h)
    call yaml_mapping_open('input structure info',label='id001')
    call yaml_map('epot',atoms%epot)
    call yaml_map('fnrm',fnrm)
    !write(*,*) 'epot=',atoms%epot
    !write(*,*) 'fnrm=',fnrm
    call atom_copy_old(atoms,atoms_center,'atoms->atoms_center')
    fsum(1)=0.d0 ; fsum(2)=0.d0 ; fsum(3)=0.d0
    do iat=1,atoms%nat
        fsum(1)=fsum(1)+atoms_center%fat(1,iat)
        fsum(2)=fsum(2)+atoms_center%fat(2,iat)
        fsum(3)=fsum(3)+atoms_center%fat(3,iat)
    enddo
    !write(*,'(a,3es14.5)') 'fsum=',fsum(1),fsum(2),fsum(3)
    call yaml_map('fsum_x',fsum(1))
    call yaml_map('fsum_y',fsum(2))
    call yaml_map('fsum_z',fsum(3))
    !write(*,'(a,es14.5)') 'finite difference step size: ',h
    call yaml_mapping_close()
    call yaml_sequence_open('force_error')
    do iat=1,atoms%nat
        call yaml_sequence(advance='no')
        call get_rat_iat(atoms,iat,xyz)
        !testing x-component of the force
        xyz(1)=xyz(1)-h
        call set_rat_iat(atoms,iat,xyz)
        call cal_potential_forces(parini,atoms)
        epot_l=atoms%epot
        xyz(1)=xyz(1)+2.d0*h
        call set_rat_iat(atoms,iat,xyz)
        call cal_potential_forces(parini,atoms)
        epot_r=atoms%epot
        fd=-(epot_r-epot_l)/(2.d0*h)
        call yaml_map('epot_l_x',epot_l)
        call yaml_map('epot_r_x',epot_r)
        call yaml_map('epot_diff_x',epot_r-epot_l)
        !write(*,'(a,2es24.15,es14.5)') 'EPOTs: ',epot_l,epot_r,epot_r-epot_l
        call yaml_map('iat',iat)
        call yaml_map('F_x',atoms_center%fat(1,iat))
        call yaml_map('err_x',fd-atoms_center%fat(1,iat),fmt='(f20.12)')
        !write(*,'(a,i,a,es11.2,2x,es19.10)') 'F_x error of atom ',iat,' is', &
        !    fd-atoms_center%fat(1,iat),atoms_center%fat(1,iat)
        xyz(1)=xyz(1)-h
        call set_rat_iat(atoms,iat,xyz)

        !testing y-component of the force
        xyz(2)=xyz(2)-h
        call set_rat_iat(atoms,iat,xyz)
        call cal_potential_forces(parini,atoms)
        epot_l=atoms%epot
        xyz(2)=xyz(2)+2.d0*h
        call set_rat_iat(atoms,iat,xyz)
        call cal_potential_forces(parini,atoms)
        epot_r=atoms%epot
        fd=-(epot_r-epot_l)/(2.d0*h)
        call yaml_map('epot_l_y',epot_l)
        call yaml_map('epot_r_y',epot_r)
        call yaml_map('epot_diff_y',epot_r-epot_l)
        !write(*,'(a,2es24.15,es14.5)') 'EPOTs: ',epot_l,epot_r,epot_r-epot_l
        call yaml_map('F_y',atoms_center%fat(2,iat))
        call yaml_map('err_y',fd-atoms_center%fat(2,iat),fmt='(f20.12)')
        !write(*,'(a,i,a,es11.2,2x,es19.10)') 'F_y error of atom ',iat,' is', &
        !    fd-atoms_center%fat(2,iat),atoms_center%fat(2,iat)
        xyz(2)=xyz(2)-h
        call set_rat_iat(atoms,iat,xyz)

        !testing z-component of the force
        xyz(3)=xyz(3)-h
        call set_rat_iat(atoms,iat,xyz)
        call cal_potential_forces(parini,atoms)
        epot_l=atoms%epot
        xyz(3)=xyz(3)+2.d0*h
        call set_rat_iat(atoms,iat,xyz)
        call cal_potential_forces(parini,atoms)
        epot_r=atoms%epot
        fd=-(epot_r-epot_l)/(2.d0*h)
        call yaml_map('epot_l_z',epot_l)
        call yaml_map('epot_r_z',epot_r)
        call yaml_map('epot_diff_z',epot_r-epot_l)
        !write(*,'(a,2es24.15,es14.5)') 'EPOTs: ',epot_l,epot_r,epot_r-epot_l
        call yaml_map('F_z',atoms_center%fat(3,iat))
        call yaml_map('err_z',fd-atoms_center%fat(3,iat),fmt='(f20.12)')
        !write(*,'(a,i,a,es11.2,2x,es19.10)') 'F_z error of atom ',iat,' is', &
        !    fd-atoms_center%fat(3,iat),atoms_center%fat(3,iat)
        xyz(3)=xyz(3)-h
        call set_rat_iat(atoms,iat,xyz)
    enddo
    call yaml_sequence_close()
    call final_potential_forces(parini,atoms)
end subroutine testforces_fd
!*****************************************************************************************
subroutine teststress_fd(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_potential, only: potential
    use mod_processors, only: iproc
    use mod_acf, only: acf_read
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    integer :: i,j,k,iat
    type(typ_atoms):: atoms 
    real(8):: hij,vol, tt
    real(8), allocatable::rat_int(:,:)
    real(8):: stress(3,3),stress_center(3,3),Ed(3,3)
    integer, parameter:: m=3
    real(8), parameter:: h=5.d-5
    real(8):: ener(-m:m), c(-m:m)=(/-1.d0,9.d0,-45.d0,0.d0,45.d0,-9.d0,1.d0/)
    potential=trim(parini%potential_potential)
   
    call acf_read(parini,'posinp.acf',1,atoms=atoms)
    allocate(rat_int(3,atoms%nat))
    call update_ratp(atoms)
    call rxyz_cart2int_alborz(atoms%nat,atoms%cellvec,atoms%ratp,rat_int)
    call vc_init_potential_forces(atoms)
    call cal_potential_forces_vc(iproc,atoms%nat,rat_int,atoms%cellvec,atoms%pressure,atoms%fat,atoms%celldv,atoms%stress,atoms%epot,atoms%enth)
    write(*,*) "epot=",atoms%epot
    stress_center(1:3,1:3)=atoms%stress(1:3,1:3)
    write (*,*) "Old_Stress:"
    write (*,'(3es14.5)') stress_center(1:3,1:3)
    call cell_vol(atoms%nat,atoms%cellvec,vol)
    vol=vol*real(atoms%nat,8)
!Calculating the stress tensor by using celldv
    do i=1,3
        do j=1,3
            tt=0.d0
            do k=1,3
                tt=tt+(atoms%celldv(i,k)*atoms%cellvec(j,k))
            enddo
            stress(i,j)=-tt/vol
!            write (*,'(2i2,es14.5,es10.1)') i,j,stress(i,j),stress(i,j)+stress_center(i,j)
        enddo
    enddo
    write(*,*) "The new stress calculated by using celldv:"
    write (*,'(3es14.5,8x,3es10.1)') stress(1:3,1:3),stress(1:3,1:3)+stress_center(1:3,1:3)

!Calculating the derivation of energy (7 points)    
     do i=1,3
         do j=1,3
            tt=0.d0
            hij=atoms%cellvec(i,j)
                do k=-m,m
                    atoms%cellvec(i,j)=hij+k*h
                    call cal_potential_forces_vc(iproc,atoms%nat,rat_int,atoms%cellvec,atoms%pressure,atoms%fat,atoms%celldv,atoms%stress,atoms%epot,atoms%enth)
                    
                    ener(k)=atoms%epot
                    tt=tt+(ener(k)*c(k))
                enddo
            Ed(i,j)=tt/(60.d0*h)
            atoms%cellvec(i,j)=hij
        enddo
    enddo
!Calculating the stress tensor by using the derivaton of energy respect to h.
    do i=1,3
        do j=1,3
            tt=0.d0
            do k=1,3
                tt=tt+(Ed(i,k)*atoms%cellvec(j,k))
            enddo
            stress(i,j)=-tt/vol
!            write (*,'(2i2,es14.5,es10.1)') i,j,stress(i,j),stress(i,j)+stress_center(i,j)
        enddo
    enddo
    write(*,*) "fat"
        do iat=1,atoms%nat
            write(*,'(3es12.4)') atoms%fat(1,iat),atoms%fat(2,iat),atoms%fat(3,iat)
        enddo
    write(*,*) "This is the volume in test_task:  ",vol
    !write(*,*) "The new stress calculated by finite difference:"
    !write (*,'(3es14.5,8x,3es10.1)') stress(1:3,1:3),stress(1:3,1:3)+stress_center(1:3,1:3)
write(*,*)" ------------------------------------------------------"
write(*,'(a11,8x,a11,5x,a11)') " FD-stress:",'VC_BLJ:','Diff=VC_BLJ-FD:'
do i=1,3
    do j=1,3
     write(*,"(es12.3,5x,es12.3,5x,es12.3)") stress(i,j),stress_center(i,j),stress_center(i,j)+stress(i,j)
    enddo
enddo
write(*,*) "------------------------------------------------------"
end subroutine teststress_fd
!*****************************************************************************************
subroutine teststress_fd_cellvec(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp, update_rat
    use mod_potential, only: potential
    use mod_processors, only: iproc
    use mod_acf, only: acf_read
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    integer :: i,j,k
    type(typ_atoms):: atoms 
    real(8):: hij,vol, tt
    real(8), allocatable::rat_int(:,:)
    real(8):: stress(3,3),stress_center(3,3),Ed(3,3)
    integer, parameter:: m=3
    real(8), parameter:: h=5.d-5
    real(8):: ener(-m:m), c(-m:m)=(/-1.d0,9.d0,-45.d0,0.d0,45.d0,-9.d0,1.d0/)
    potential=trim(parini%potential_potential)
   
    call acf_read(parini,'posinp.acf',1,atoms=atoms) !Cartezy
    allocate(rat_int(3,atoms%nat))
    call update_ratp(atoms)
    call rxyz_cart2int_alborz(atoms%nat,atoms%cellvec,atoms%ratp,rat_int)
    call init_potential_forces(parini,atoms)
    call cell_vol(atoms%nat,atoms%cellvec,vol)   !to calculate the volume per atom
    vol=vol*real(atoms%nat,8)                    !the volume of the cell
    call cal_potential_forces(parini,atoms)
    stress_center(1:3,1:3)=atoms%stress(1:3,1:3)

!Calculating the derivation of energy (7 points)    
     do i=1,3
         do j=1,3
            tt=0.d0
            hij=atoms%cellvec(i,j)
                do k=-m,m
                    atoms%cellvec(i,j)=hij+k*h
                     call rxyz_int2cart_alborz(atoms%nat,atoms%cellvec,rat_int,atoms%ratp)  !to convert the reduced coordinates to cartezian
                     call update_rat(atoms)
                     call cal_potential_forces(parini,atoms)
                    ener(k)=atoms%epot
                    tt=tt+(ener(k)*c(k))
                enddo
            Ed(i,j)=tt/(60.d0*h)
            atoms%cellvec(i,j)=hij
        enddo
    enddo
!Calculating the stress tensor by using the derivaton of energy respect to h.
    do i=1,3
        do j=1,3
            tt=0.d0
            do k=1,3
                tt=tt+(Ed(i,k)*atoms%cellvec(j,k))
            enddo
            !test stress fails due to inconsistency of definition: per volume
            stress(i,j)=-tt !/vol
!            write (*,'(2i2,es14.5,es10.1)') i,j,stress(i,j),stress(i,j)+stress_center(i,j)
        enddo
    enddo
write(*,*)" ------------------------------------------------------"
write(*,'(a11,8x,a11,5x,a11)') " FD-stress:",'NN-stress:','Diff=NN-FD:'
do i=1,3
    do j=1,3
     write(*,"(es12.3,5x,es12.3,5x,es12.3)") stress(i,j),stress_center(i,j),stress_center(i,j)-stress(i,j)
    enddo
enddo
write(*,*) "------------------------------------------------------"
end subroutine teststress_fd_cellvec
!*****************************************************************************************
