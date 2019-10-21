!*****************************************************************************************
subroutine dynamics(parini)
    use mod_parini, only: typ_parini
    use mod_acf, only: acf_read, acf_write
    use mod_potential, only: potential, perfstatus
    use mod_atoms, only: typ_atoms, typ_file_info, set_ndof, atom_deallocate_old
    use mod_dynamics, only: dt, nmd,nfreq,md_method
    use mod_processors, only: iproc
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variables
    type(typ_atoms):: atoms
    character(56):: comment
    
    potential=trim(parini%potential_potential)
    md_method=trim(parini%md_method_dynamics)
    dt=parini%dt_dynamics
    dt=dt*41.341373336493 ! convert fs to atomic unit
    nmd=parini%nmd_dynamics
    nfreq=parini%nfreq_dynamics
    
        write(*,'(a,es15.5)') 'ERROR: dt must be set in input file and dt>0, dt= ',dt 
    if(.not. dt>0.d0) then
        write(*,'(a,es15.5)') 'ERROR: dt must be set in input file and dt>0, dt= ',dt 
        stop
    endif
    !dt=4.d-3 !good for LJ
    call acf_read(parini,'posinp.acf',1,atoms=atoms)
    call set_ndof(atoms)
    
    if(trim(md_method)=='nve') then
        call md_nve(parini,atoms)
    else if(trim(md_method)=='nvt_langev') then
        call md_nvt_langevin(parini,atoms)
    else if(trim(md_method)=='nvt_nosecp') then
        call md_nvt_nose_hoover_cp(parini,atoms)
    else if(trim(md_method)=='nvt_nose') then
        call md_nvt_nose_hoover_chain(parini,atoms)
        !call md_nvt_nose_hoover(parini,atoms)
    else if(trim(md_method)=='nph') then
        call md_nph(parini,atoms)
    else 
        write(*,'(2a)') 'ERROR: unknown parini%md_method:',trim(md_method)
        stop
    endif
    call atom_deallocate_old(atoms)
end subroutine dynamics
!*****************************************************************************************
subroutine md_nve(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_potential, only: potential, perfstatus
    use mod_atoms, only: typ_atoms, typ_file_info, atom_copy_old, set_atomic_mass
    use mod_atoms, only: atom_deallocate_old, update_ratp, update_rat
    use mod_acf, only: acf_write
    use mod_dynamics, only: dt, nmd,nfreq
    use mod_velocity, only: set_velocities
    use mod_processors, only: iproc
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    type(typ_atoms):: atoms, atoms_old
    type(typ_file_info):: file_info
    integer:: iat, i
    integer:: imd
    real(8):: etot, epotold, etotold, ekin_target
    real(8):: DNRM2, fnrm, t1, aboltzmann, totmass, temp, etotavg,mass_conv
    real(8):: rcm(3), vcm(3)
    character(56):: comment
    !aboltzmann=8.6173324d-5 !eV/K
    aboltzmann=8.6173324d-5/27.211385d0 !a.u./K
    
    call atom_copy_old(atoms,atoms_old,'atoms->atoms_old')
    call init_potential_forces(parini,atoms)
    
    call set_atomic_mass(atoms)
    totmass=0.d0
    do iat=1,atoms%nat
        totmass=totmass+atoms%amass(iat)
    enddo
    ! temp in K and ekin_target in a.u
    ekin_target=1.5d0*atoms%nat*aboltzmann*parini%init_temp_dynamics 
    call set_velocities(atoms,ekin_target)
    epotold=atoms%epot
    call cal_potential_forces(parini,atoms)
    write(*,'(a,2e20.10)') 'epotold,epot',epotold,atoms%epot
    file_info%filename_positions='posout.acf'
    file_info%file_position='new'
    file_info%print_force=parini%print_force_dynamics
    call acf_write(file_info,atoms=atoms,strkey='posout')
    atoms%ekin=0.d0
    vcm(1:3)=0.d0
    rcm(1:3)=0.d0
    call update_ratp(atoms)
    do iat=1,atoms%nat
        t1=atoms%amass(iat)
        atoms%ekin=atoms%ekin+t1*(atoms%vat(1,iat)**2+atoms%vat(2,iat)**2+atoms%vat(3,iat)**2)
        vcm(1:3)=vcm(1:3)+t1*atoms%vat(1:3,iat)
        rcm(1:3)=rcm(1:3)+t1*atoms%ratp(1:3,iat)
    enddo
    vcm(1:3)=vcm(1:3)/totmass
    rcm(1:3)=rcm(1:3)/totmass
    atoms%ekin=0.5d0*atoms%ekin
    etot=atoms%epot+atoms%ekin
    temp=atoms%ekin/(1.5d0*atoms%nat*aboltzmann)
    etotold=etot
    write(21,*) "#imd, epot, ekin, etot, etot-etotold, temp"
    write(21,'(i9,5es25.15)') 0,atoms%epot,atoms%ekin,etot,etot-etotold,temp
    write(22,'(i9,6es20.10)') 0,rcm(1:3),vcm(1:3)
    do imd=1,nmd
        call update_ratp(atoms)
        do iat=1,atoms%nat
            t1=dt**2/(2.d0*atoms%amass(iat))
            atoms%ratp(1:3,iat)=atoms%ratp(1:3,iat)+dt*atoms%vat(1:3,iat)+t1*atoms%fat(1:3,iat)
        enddo
        call update_rat(atoms,upall=.true.)
        call cal_potential_forces(parini,atoms)
        do iat=1,atoms%nat
            t1=dt/(2.d0*atoms%amass(iat))
            atoms%vat(1:3,iat)=atoms%vat(1:3,iat)+t1*(atoms_old%fat(1:3,iat)+atoms%fat(1:3,iat))
        enddo
        atoms_old%fat(1:3,1:atoms%nat)=atoms%fat(1:3,1:atoms%nat)
        atoms%ekin=0.d0
        vcm(1:3)=0.d0
        rcm(1:3)=0.d0
        call update_ratp(atoms)
        do iat=1,atoms%nat
            t1=atoms%amass(iat)
            atoms%ekin=atoms%ekin+t1*(atoms%vat(1,iat)**2+atoms%vat(2,iat)**2+atoms%vat(3,iat)**2)
            vcm(1:3)=vcm(1:3)+t1*atoms%vat(1:3,iat)
            rcm(1:3)=rcm(1:3)+t1*atoms%ratp(1:3,iat)
        enddo
        vcm(1:3)=vcm(1:3)/totmass
        rcm(1:3)=rcm(1:3)/totmass
        atoms%ekin=0.5d0*atoms%ekin
        etot=atoms%epot+atoms%ekin
        !if(imd>100) etotavg=etotavg+etot
        temp=atoms%ekin/(1.5d0*atoms%nat*aboltzmann)
        if(mod(imd,1)==0) then
            write(21,'(i9,5es25.15)') imd,atoms%epot,atoms%ekin,etot,etot-etotold,temp
            write(22,'(i9,6es20.10)') imd,rcm(1:3),vcm(1:3)
        endif
        !if(mod(imd,10)==0) then
        if(mod(imd,nfreq)==0) then
            file_info%file_position='append'
            write(23,'(i9,5es25.15)') imd,atoms%epot,atoms%ekin,etot,etot-etotold,temp
            call acf_write(file_info,atoms=atoms,strkey='posout')
        endif
        etotold=etot
    enddo !end of loop over imd
    call final_potential_forces(parini,atoms)
    !call atom_deallocate_old(atoms)
    call atom_deallocate_old(atoms_old)
end subroutine md_nve
!*****************************************************************************************
subroutine md_nph(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_potential, only: potential, perfstatus
    use mod_atoms, only: typ_atoms, typ_file_info, set_rat, update_ratp, update_rat
    use mod_acf, only: acf_write
    use mod_velocity, only: set_velocities
    use mod_dynamics, only: dt, nmd
    use mod_processors, only: iproc
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    type(typ_atoms):: atoms
    type(typ_file_info):: file_info
    integer:: iat,i,j,k,imd
    real(8):: etot, epotold, etotold,vol
    !real(8):: enthalpy, epotold, enthalpyold,vol
    real(8):: enthalpy, enthalpyold,trc_hdthd,trc,enth,enthold
    real(8):: t1,aboltzmann,totmass,temp,w,cell_ekin
    real(8):: rcm(3),vcm(3)
    real(8):: cellvec_old(3,3),cellvec_new(3,3),cellvec_inv(3,3),hd_new(3,3),hd(3,3),hdd(3,3),eta(3,3),hdthd(1:3,1:3)
    real(8):: stress_tot(3,3),vvt(3,3),g(3,3),ginv(3,3),gd(3,3),ginv_gd(3,3),F(3,3)
    real(8), allocatable::rat_new(:,:),rxyz_red(:,:),s(:,:),s_new(:,:),s_old(:,:),vat_new(:,:)
    real(8), allocatable::sd(:,:),sdd(:,:),hsd(:,:)
    real(8), allocatable::fat_int(:,:),tt1(:,:),tt2(:,:)
    real(8)::alpha,beta,gama,a,b,c
    character(56):: comment
    aboltzmann=8.6173324d-5
    etot=0.d0;enth=0.d0
    allocate(s(3,atoms%nat),s_new(3,atoms%nat),s_old(3,atoms%nat),sd(3,atoms%nat),sdd(3,atoms%nat), &
            tt1(3,atoms%nat),tt2(3,atoms%nat),fat_int(3,atoms%nat),rat_new(3,atoms%nat) & 
            ,rxyz_red(3,atoms%nat),hsd(3,atoms%nat),vat_new(3,atoms%nat))
    call init_potential_forces(parini,atoms)
    atoms%amass(1:atoms%nat)=1.d0
    totmass=0.d0
    do iat=1,atoms%nat
        totmass=totmass+atoms%amass(iat)
    enddo
    w=totmass/atoms%nat
    call set_velocities(atoms)
    epotold=atoms%epot
    call update_ratp(atoms)
    call rxyz_cart2int_alborz(atoms%nat,atoms%cellvec,atoms%ratp,s)
    call backtocell_alborz(atoms%nat,atoms%cellvec,s)
    do iat=1,atoms%nat
        atoms%ratp(1:3,iat)=matmul(atoms%cellvec,s(1:3,iat))
    enddo
    call update_rat(atoms,upall=.true.)
    call cal_potential_forces(parini,atoms)
    write(*,'(a,2e20.10)') 'epotold,epot',epotold,atoms%epot
    file_info%filename_positions='posout.acf'
    file_info%file_position='new'
    file_info%print_force=parini%print_force_dynamics
    call acf_write(file_info,atoms=atoms,strkey='posout')
    atoms%ekin=0.d0
    vcm(1:3)=0.d0
    rcm(1:3)=0.d0
    call update_ratp(atoms)
    do iat=1,atoms%nat
        t1=atoms%amass(iat)
        atoms%ekin=atoms%ekin+t1*(atoms%vat(1,iat)**2+atoms%vat(2,iat)**2+atoms%vat(3,iat)**2)
        vcm(1:3)=vcm(1:3)+t1*atoms%vat(1:3,iat)
        rcm(1:3)=rcm(1:3)+t1*atoms%ratp(1:3,iat)
    enddo
    vcm(1:3)=vcm(1:3)/totmass
    rcm(1:3)=rcm(1:3)/totmass
    atoms%ekin=0.5d0*atoms%ekin
    etot=atoms%epot+atoms%ekin 
    temp=atoms%ekin/(1.5d0*atoms%nat*aboltzmann)
    call cell_vol(atoms%nat,atoms%cellvec,vol)
    vol=vol*real(atoms%nat,8)
    etot=etot+vol*atoms%pressure
    enth=atoms%epot+vol*atoms%pressure
    etotold=etot
    enthold=enth
    write(21,*) '#    imd,  atoms_epot,  atoms_ekin, etot  ,enth, enth-enthold,  '
    write(21,'(i9,5es25.15)') 0,atoms%epot,atoms%ekin,etot,enth,enth-enth
    
    !-----------------------Start MD-----------------------------------------
    call invertmat_alborz(atoms%cellvec,cellvec_inv)
    do iat=1,atoms%nat
        sd(1:3,iat)=matmul(cellvec_inv,atoms%vat(1:3,iat))
    enddo
    gd=0.d0
    hd=0.d0
    do imd=1,nmd
        g(1:3,1:3)=matmul(transpose(atoms%cellvec),atoms%cellvec)
        call invertmat_alborz(g,ginv)
        gd(1:3,1:3)=matmul(transpose(hd),atoms%cellvec)+matmul(transpose(atoms%cellvec),hd)
        ginv_gd(1:3,1:3)=matmul(ginv,gd)
        tt2=0.d0
        do iat=1,atoms%nat
            t1=atoms%amass(iat)
            do i=1,3
                do j=1,3
                    tt2(i,iat)=tt2(i,iat)+ginv_gd(i,j)*sd(j,iat)
                enddo
            enddo
        enddo
        call invertmat_alborz(atoms%cellvec,cellvec_inv)
        fat_int=0.d0
        do iat=1,atoms%nat
            do i=1,3
                do j=1,3
                    fat_int(i,iat)=fat_int(i,iat)+cellvec_inv(i,j)*atoms%fat(j,iat)
                enddo
            enddo
        enddo
        do iat=1,atoms%nat
            t1=atoms%amass(iat)
            tt1(1:3,iat)=fat_int(1:3,iat)/t1
            s_new(1:3,iat)=s(1:3,iat)+dt*sd(1:3,iat)+tt1(1:3,iat)*dt**2-tt2(1:3,iat)*dt**2
        enddo
        call cell_vol(atoms%nat,atoms%cellvec,vol)
        vol=vol*real(atoms%nat,8)
        !eta=vol*transpose(cellvec_inv) 
        call cross_product_alborz(atoms%cellvec(1,2),atoms%cellvec(1,3),eta(1,1))
        call cross_product_alborz(atoms%cellvec(1,3),atoms%cellvec(1,1),eta(1,2))
        call cross_product_alborz(atoms%cellvec(1,1),atoms%cellvec(1,2),eta(1,3))
        atoms%vat=0.d0
        do iat=1,atoms%nat
            do i=1,3
                do j=1,3
                    atoms%vat(i,iat)=atoms%vat(i,iat)+atoms%cellvec(i,j)*sd(j,iat)
                enddo
            enddo
        enddo
        vvt(1:3,1:3)=0.d0
        do iat=1,atoms%nat
            t1=atoms%amass(iat)
            do i=1,3
                do j=1,3
                    vvt(i,j)=vvt(i,j)+t1*(atoms%vat(i,iat)*atoms%vat(j,iat))
                enddo
            enddo
        enddo
        stress_tot(1:3,1:3)=(vvt(1:3,1:3)/vol)-atoms%stress(1:3,1:3)
        F(1:3,1:3)=matmul(stress_tot,eta)-atoms%pressure*eta(1:3,1:3) 
        hdd(1:3,1:3)=F(1:3,1:3)/w
        cellvec_new(1:3,1:3)=atoms%cellvec(1:3,1:3)+dt*hd(1:3,1:3)+hdd(1:3,1:3)*dt**2
        do iat=1,atoms%nat
            do i=1,3
                rxyz_red(i,iat)=modulo(modulo(s_new(i,iat),1.d0),1.d0)
            enddo
        enddo
        rat_new=0.d0
        do iat=1,atoms%nat
            do i=1,3
                do j=1,3
                    rat_new(i,iat)=rat_new(i,iat)+cellvec_new(i,j)*rxyz_red(j,iat)
                enddo
            enddo
        enddo
        s_old(1:3,1:atoms%nat)=s(1:3,1:atoms%nat)
        s(1:3,1:atoms%nat)=s_new(1:3,1:atoms%nat)
        call set_rat(atoms,rat_new,setall=.true.)
        cellvec_old(1:3,1:3)=atoms%cellvec(1:3,1:3)
        atoms%cellvec(1:3,1:3)=cellvec_new(1:3,1:3)
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        call cal_potential_forces(parini,atoms)
        write(*,*) 'FNRM ',sqrt(sum(atoms%fat(1:3,1:atoms%nat)**2))
        atoms%vat=0.d0
        do iat=1,atoms%nat
           sd(1:3,iat)=(s(1:3,iat)-s_old(1:3,iat))/dt
           do i=1,3
                do j=1,3
                    atoms%vat(i,iat)=atoms%vat(i,iat)+atoms%cellvec(i,j)*sd(j,iat)
                enddo
            enddo
        enddo
        atoms%ekin=0.d0
        vcm(1:3)=0.d0
        rcm(1:3)=0.d0
        call update_ratp(atoms)
        do iat=1,atoms%nat
            t1=atoms%amass(iat)
            atoms%ekin=atoms%ekin+t1*(atoms%vat(1,iat)**2+atoms%vat(2,iat)**2+atoms%vat(3,iat)**2)
            vcm(1:3)=vcm(1:3)+t1*atoms%vat(1:3,iat)
            rcm(1:3)=rcm(1:3)+t1*atoms%ratp(1:3,iat)
        enddo
        vcm(1:3)=vcm(1:3)/totmass
        rcm(1:3)=rcm(1:3)/totmass
        atoms%ekin=0.5d0*atoms%ekin
        temp=atoms%ekin/(1.5d0*atoms%nat*aboltzmann)
        hd(1:3,1:3)=(atoms%cellvec(1:3,1:3)-cellvec_old(1:3,1:3))/dt
        hdthd(1:3,1:3)=matmul(transpose(hd),hd)
        trc=hdthd(1,1)+hdthd(2,2)+hdthd(3,3)
        cell_ekin=0.5d0*w*trc
        etot=atoms%epot+atoms%ekin+cell_ekin+vol*atoms%pressure
        enth=atoms%epot+atoms%pressure*vol
        if(mod(imd,1)==0) then
            write(21,'(i9,5es25.15)') imd,atoms%epot,atoms%ekin,etot, enth,enth-enthold
            write(22,'(i9,6es20.10)') imd,rcm(1:3),vcm(1:3)
        endif
        if(mod(imd,10)==0) then
            file_info%file_position='append'
            call acf_write(file_info,atoms=atoms,strkey='posout')
        endif
        etotold=etot
        enthold=enth
        !enthalpyold=atoms%enth
    enddo !end of loop over imd
    call final_potential_forces(parini,atoms)
    deallocate(s,s_new,s_old,sd,tt1,tt2,fat_int,rat_new,rxyz_red)
end subroutine md_nph
!*****************************************************************************************
subroutine ekin_temprature(atoms,temp,vcm,rcm,totmass) 
    use mod_atoms, only: typ_atoms, update_ratp
    implicit none
    type(typ_atoms):: atoms
    integer :: iat
    real(8) :: vcm(3), rcm(3), t1, tmp
    real(8) :: aboltzmann, totmass, temp 

    aboltzmann= 3.1668113916289087d-6
    atoms%ekin=0.d0
    vcm(1:3)=0.d0
    rcm(1:3)=0.d0
    totmass= 0.d0
    call update_ratp(atoms)
    do iat=1,atoms%nat
        t1=atoms%amass(iat)
        totmass=totmass + t1
        atoms%ekin=atoms%ekin+t1*(atoms%vat(1,iat)**2+atoms%vat(2,iat)**2+atoms%vat(3,iat)**2)
        vcm(1:3)=vcm(1:3)+t1*atoms%vat(1:3,iat)
        rcm(1:3)=rcm(1:3)+t1*atoms%ratp(1:3,iat)
    enddo
    vcm(1:3)=vcm(1:3)/totmass
    rcm(1:3)=rcm(1:3)/totmass
    atoms%ekin=0.5d0*atoms%ekin
    temp=atoms%ekin/(1.5d0*atoms%nat*aboltzmann)

end subroutine ekin_temprature
!********************************************************************************************
