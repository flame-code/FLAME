!*****************************************************************************************
subroutine md_nvt_langevin(parini,atoms)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_potential, only: potential, perfstatus
    use mod_atoms, only: typ_atoms, typ_file_info
    use mod_dynamics, only: dt, nmd
    use mod_processors, only: iproc
    use mod_potential, only: bias 
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variables
    type(typ_atoms):: atoms
    type(typ_file_info):: file_info
    integer:: iat, ierr, nat_t, i
    integer:: imd, ff
    real(8):: etot, epotold, etotold
    real(8):: DNRM2, fnrm, t1, aboltzmann, totmass, temp, etotavg
    real(8):: t2,t3,t4 
    real(8):: rl, ru 
    real(8):: rcm(3), vcm(3)
    real(8):: scale_vat, temp_trget, ekin_target
    real(8):: eta(3,atoms%nat)
    real(8):: gama, factor, factor_iat
    real(8) :: sum1, sum2, sum3
    real(8) :: kt, temp_prev, tol, tolerance 
    character(56):: comment
    real(8):: langev(atoms%nat), forces_langevin(3,atoms%nat)
    real(8):: rat_next(3,atoms%nat), vat_old(3,atoms%nat)
    real(8):: rat_init(3,atoms%nat)
    real(8):: r, dx(3) , rsq, msd1, msd2, msd3

    call random_seed() 
    rat_init=atoms%rat
    !  ___________parameters_______________________________________
    gama=1.d-3
    aboltzmann= 3.1668113916289087d-6
    temp_trget = parini%temp_dynamics
    kt = aboltzmann*temp_trget
    tolerance = 1.d-2

    open(unit=1000,file="velocity",status='replace')
    open(unit=1111,file="displace.txt",status='replace')
    call init_potential_forces(parini,atoms)
    call get_atomic_mass(atoms,totmass)
    langev(:)=sqrt(2*gama*atoms%amass(:)*kt/dt)
    ekin_target=1.5d0*atoms%nat*aboltzmann*parini%init_temp_dynamics
    !_______________________initial velocity __________________________
    if (parini%restart_dynamics )then
        open(unit=1001,file="velocity_r",status='old')
        read(1001,*)
        read(1001,*) 
        read(1001,*)
        do iat=1,atoms%nat
            read(1001,*) atoms%vat(1,iat),atoms%vat(2,iat),atoms%vat(3,iat)
        enddo
    else
        if ( parini%init_temp_dynamics==0.d0) then
            atoms%vat(:,:)=0.d0
        else
            call set_velocities(atoms, ekin_target)
        endif
    endif
    !____________________________________________________________________
    epotold=atoms%epot
    call cal_potential_forces(parini,atoms)
    if (trim(bias)=='yes') then
        call plane_repulsion(atoms)
    endif

    write(*,'(a,2e20.10)') 'epotold,epot',epotold,atoms%epot
    file_info%filename_positions='posout.acf'
    file_info%file_position='new'
    file_info%print_force=parini%print_force_dynamics
    call acf_write(file_info,atoms=atoms,strkey='posout')
    !____________________________________________________________________
    call ekin_temprature(atoms,temp,vcm,rcm,totmass) 
    etot=atoms%epot+atoms%ekin
    etotold=etot
    write(21,'(i9,4es25.15)') 0,etot,atoms%epot,atoms%ekin,temp
    write(22,'(i9,6es20.10)') 0,rcm(1:3),vcm(1:3)
   !_________________________ The first md step _________________________
    call set_langevin_randforce(eta,atoms%nat)
    if (parini%restart_dynamics )then
        t1=dt*dt
    else
        t1=0.5*dt*dt
    endif
    do iat=1,atoms%nat
        forces_langevin(1:3,iat)=atoms%fat(1:3,iat)+langev(iat)*eta(1:3,iat)-gama*atoms%amass(iat)*atoms%vat(1:3,iat)
        rat_next(1:3,iat)=atoms%rat(1:3,iat) + t1*forces_langevin(1:3,iat)/atoms%amass(iat) + dt*atoms%vat(1:3,iat)
        atoms%vat(1:3,iat)=(rat_next(1:3,iat)-atoms%rat(1:3,iat))/dt
        atoms%rat(1:3,iat)=rat_next(1:3,iat)
    enddo
    !call back_to_cell(atoms)

    !_____________________________________________________________________
    do imd=2,nmd
        parini%time_dynamics = (imd-1)*dt

        vat_old =atoms%vat
        parini%cal_charge = .false.
        !if(mod(imd,10)==0) then
        !    parini%cal_charge = .true.
        !endif
        epotold=atoms%epot
        call cal_potential_forces(parini,atoms)
        if (trim(bias)=='yes') then
            call plane_repulsion(atoms)
        endif
        call ekin_temprature(atoms,temp,vcm,rcm,totmass) 
        etot=atoms%epot+atoms%ekin
        etotold=etot
        write(21,'(i9,4es25.15)') imd-1,etot,atoms%epot,atoms%ekin,temp
        write(22,'(i9,6es20.10)') imd-1,rcm(1:3),vcm(1:3)

        call set_langevin_randforce(eta,atoms%nat)
        do iat=1,atoms%nat
            forces_langevin(1:3,iat)=atoms%fat(1:3,iat)+langev(iat)*eta(1:3,iat)-gama*atoms%amass(iat)*atoms%vat(1:3,iat)
            rat_next(1:3,iat)= atoms%rat(1:3,iat)+atoms%vat(1:3,iat)*dt &
            &       + dt*dt*forces_langevin(1:3,iat)/atoms%amass(iat)
        enddo
        tol=1.d0
        call ekin_temprature(atoms,temp_prev,vcm,rcm,totmass) 
        do while (tol>tolerance)
            atoms%vat = 0.5*((rat_next-atoms%rat)/dt+vat_old)
            call ekin_temprature(atoms,temp,vcm,rcm,totmass) 
            do iat=1,atoms%nat
                forces_langevin(1:3,iat)=atoms%fat(1:3,iat)+langev(iat)*eta(1:3,iat)-gama*atoms%amass(iat)*atoms%vat(1:3,iat)
                rat_next(1:3,iat)= atoms%rat(1:3,iat)+ vat_old(1:3,iat)*dt &
                &       + dt*dt*forces_langevin(1:3,iat)/atoms%amass(iat)
            enddo
            tol=dabs(temp-temp_prev)/dabs(temp_prev)
            temp_prev=temp
        enddo

        atoms%vat = (rat_next - atoms%rat )/dt
        atoms%rat = rat_next
        !call back_to_cell(atoms)
        if(mod(imd,100)==0) then
            file_info%file_position='append'
            call acf_write(file_info,atoms=atoms,strkey='posout')
            write(1111,*) '#'
            write(1111,*) '#    imd = ',imd, parini%time_dynamics
            write(1111,*) '#'
            msd1= 0.d0
            msd2= 0.d0
            msd3= 0.d0
            do iat=1,atoms%nat

                dx(1:3)=atoms%rat(:,iat)-rat_init(:,iat)
                rsq=(dx(1)**2+dx(2)**2+dx(3)**2)
                r=sqrt(dx(1)**2+dx(2)**2+dx(3)**2)
                msd1 = msd1 + rsq                !all directions
                msd2 = msd2 + dx(1)**2+dx(2)**2  !x,y directions
                msd3 = msd3 + dx(3)**2           !z   direction
                write(1111,'(i5,2a5,4es25.17)')iat," ",atoms%sat(iat), dx, r
            enddo
            write(1111,*) '  MSD    = ', msd1/atoms%nat 
            write(1111,*) '  MSD_xy = ', msd2/atoms%nat 
            write(1111,*) '  MSD_z  = ', msd3/atoms%nat 
            write(1111,*) '# -----------------------------------------------'
        
        endif
        etotold=etot
    enddo !end of loop over imd
    close(1000)
    call final_potential_forces(parini,atoms)
end subroutine md_nvt_langevin
!*****************************************************************************************
subroutine md_nvt_nose_hoover_cp(parini,atoms)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_potential, only: potential, perfstatus
    use mod_atoms, only: typ_atoms, typ_file_info
    use mod_dynamics, only: dt, nmd
    use mod_processors, only: iproc
    use mod_potential, only: bias 
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variables
    type(typ_atoms):: atoms
    type(typ_file_info):: file_info
    integer:: iat, ierr, nat_t, i
    integer:: imd, ff,ntherm,tt
    real(8):: etot, epotold, etotold
    real(8):: DNRM2, fnrm, t1, aboltzmann, totmass, temp, etotavg
    real(8):: t2,t3,t4 
    real(8):: rl, ru 
    real(8):: rcm(3), vcm(3)
    real(8):: scale_vat, temp_trget, ekin_target
    real(8) :: sum1, sum2, sum3
    real(8) :: kt, temp_prev, tol, tolerance 
    character(56):: comment
    real(8):: forces_nosehoover(3,atoms%nat)
    real(8):: rat_next(3,atoms%nat), vat_old(3,atoms%nat)
    real(8), allocatable ::zeta_next(:,:,:) ,zeta(:,:,:) ,zeta_prev(:,:,:) ,dzeta(:,:,:),mass_q(:),dzeta_old(:,:,:)
    real(8):: rat_init(3,atoms%nat)
    real(8):: r, dx(3) , rsq, msd1, msd2, msd3

    call random_seed() 
    rat_init=atoms%rat

    call init_potential_forces(parini,atoms)
    open(unit=1000,file="velocity",status='replace')
    open(unit=1111,file="displace.txt",status='replace')
    file_info%filename_positions='posout.acf'
    file_info%file_position='new'
    file_info%print_force=parini%print_force_dynamics
    call acf_write(file_info,atoms=atoms,strkey='posout')

    !  ___________parameters_______________________________________
    ntherm=3
    aboltzmann= 3.1668113916289087d-6
    temp_trget = parini%temp_dynamics
    kt = aboltzmann*temp_trget
    tolerance = 1.d-4

    allocate(zeta_next(3,atoms%nat,ntherm), zeta(3,atoms%nat,ntherm),&
             zeta_prev(3,atoms%nat,ntherm), dzeta(3,atoms%nat,ntherm),&
             mass_q(ntherm),dzeta_old(3,atoms%nat,ntherm))
    zeta      = 0.d0
    zeta_next = 0.d0
    zeta_prev = 0.d0
    dzeta     = 0.d0
    mass_q    =kt*((50*dt)**2)
    ekin_target=1.5d0*atoms%nat*aboltzmann*parini%init_temp_dynamics

    call get_atomic_mass(atoms,totmass)
    !_______________________initial velocity __________________________
    if (parini%restart_dynamics )then
        open(unit=1001,file="velocity_r",status='old')
        read(1001,*)
        read(1001,*) 
        read(1001,*)
        do iat=1,atoms%nat
            read(1001,*) atoms%vat(1,iat),atoms%vat(2,iat),atoms%vat(3,iat)
        enddo
    else
        if ( parini%init_temp_dynamics==0.d0) then
            atoms%vat(:,:)=0.d0
        else
            call set_velocities(atoms, ekin_target)
        endif
    endif
    !____________________________________________________________________
    call ekin_temprature(atoms,temp,vcm,rcm,totmass) 

    epotold=atoms%epot
    call cal_potential_forces(parini,atoms)
    if (trim(bias)=='yes')  call plane_repulsion(atoms)
    etot=atoms%epot+atoms%ekin
    etotold=etot

    write(*,'(a,2e20.10)') 'epotold,epot',epotold,atoms%epot
    write(21,'(i9,4es25.15)') 0,etot,atoms%epot,atoms%ekin,temp
    write(22,'(i9,6es20.10)') 0,rcm(1:3),vcm(1:3)
   !_________________________ The first md step _________________________
    imd=1
    t1=0.5*dt*dt
    do iat=1,atoms%nat
        forces_nosehoover(1:3,iat)=atoms%fat(1:3,iat)-atoms%amass(iat)*atoms%vat(1:3,iat)*dzeta(1:3,iat,1)
    enddo
    do iat=1,atoms%nat
        rat_next(1:3,iat)=atoms%rat(1:3,iat) + t1*forces_nosehoover(1:3,iat)/atoms%amass(iat) + dt*atoms%vat(1:3,iat)
    enddo
    call thermostat_evolution(atoms,zeta_next,zeta,zeta_prev,dzeta,mass_q,kt,ntherm,imd)
    do iat=1,atoms%nat
        atoms%vat(1:3,iat)=(rat_next(1:3,iat)-atoms%rat(1:3,iat))/dt
        atoms%rat(1:3,iat)=rat_next(1:3,iat)
    enddo
    call back_to_cell(atoms)
    dzeta=(zeta_next-zeta)/dt
    zeta_prev=zeta
    zeta=zeta_next
    !_____________________________________________________________________
    do imd=2,nmd
        parini%time_dynamics = (imd-1)*dt
        vat_old =atoms%vat
        dzeta_old=dzeta
        epotold=atoms%epot
        call cal_potential_forces(parini,atoms)
        if (trim(bias)=='yes')  call plane_repulsion(atoms)

        call ekin_temprature(atoms,temp,vcm,rcm,totmass) 
        etot=atoms%epot+atoms%ekin
        etotold=etot
        write(21,'(i9,4es25.15)') imd-1,etot,atoms%epot,atoms%ekin,temp
        write(22,'(i9,6es20.10)') imd-1,rcm(1:3),vcm(1:3)

        do iat=1,atoms%nat
            forces_nosehoover(1:3,iat)=atoms%fat(1:3,iat)-atoms%amass(iat)*atoms%vat(1:3,iat)*dzeta(1:3,iat,1)
            rat_next(1:3,iat)= atoms%rat(1:3,iat)+atoms%vat(1:3,iat)*dt &
            &       + dt*dt*forces_nosehoover(1:3,iat)/atoms%amass(iat)
        enddo
        call thermostat_evolution(atoms,zeta_next,zeta,zeta_prev,dzeta,mass_q,kt,ntherm,imd)
        tol=1.d0
        call ekin_temprature(atoms,temp_prev,vcm,rcm,totmass) 
        do while (tol>tolerance)
            atoms%vat = 0.5*((rat_next-atoms%rat)/dt+vat_old)
            dzeta=0.5*((zeta_next - zeta )/dt+dzeta_old)
            call ekin_temprature(atoms,temp,vcm,rcm,totmass) 
            do iat=1,atoms%nat
                forces_nosehoover(1:3,iat)=atoms%fat(1:3,iat)-atoms%amass(iat)*atoms%vat(1:3,iat)*dzeta(1:3,iat,1)
                rat_next(1:3,iat)= atoms%rat(1:3,iat)+atoms%vat(1:3,iat)*dt &
                &       + dt*dt*forces_nosehoover(1:3,iat)/atoms%amass(iat)
            enddo
            call thermostat_evolution(atoms,zeta_next,zeta,zeta_prev,dzeta,mass_q,kt,ntherm,imd)
            tol=dabs(temp-temp_prev)/dabs(temp_prev)
            temp_prev=temp
        enddo

        atoms%vat = (rat_next - atoms%rat )/dt
        atoms%rat = rat_next
        dzeta = (zeta_next - zeta )/dt
        zeta_prev=zeta
        zeta=zeta_next
        call back_to_cell(atoms)
        if(mod(imd,100)==0) then
            file_info%file_position='append'
            call acf_write(file_info,atoms=atoms,strkey='posout')
            write(1111,*) '#'
            write(1111,*) '#    imd = ',imd, parini%time_dynamics
            write(1111,*) '#'
            msd1= 0.d0
            msd2= 0.d0
            msd3= 0.d0
            do iat=1,atoms%nat

                dx(1:3)=atoms%rat(:,iat)-rat_init(:,iat)
                rsq=(dx(1)**2+dx(2)**2+dx(3)**2)
                r=sqrt(dx(1)**2+dx(2)**2+dx(3)**2)
                msd1 = msd1 + rsq                !all directions
                msd2 = msd2 + dx(1)**2+dx(2)**2  !x,y directions
                msd3 = msd3 + dx(3)**2           !z   direction
                write(1111,'(i5,2a5,4es25.17)')iat," ",atoms%sat(iat), dx, r
            enddo
            write(1111,*) '  MSD    = ', msd1/atoms%nat 
            write(1111,*) '  MSD_xy = ', msd2/atoms%nat 
            write(1111,*) '  MSD_z  = ', msd3/atoms%nat 
            write(1111,*) '# -----------------------------------------------'
         !   write(1000,*) '#'
         !   write(1000,*) '#    imd = ', imd
         !   write(1000,*) '#'
         !   do iat=1,atoms%nat
         !       write(1000,'(3es25.17)') atoms%vat(1,iat),atoms%vat(2,iat),atoms%vat(3,iat)
         !   enddo
        endif
        write(221,'(i9,4es25.15)') imd,etot,atoms%epot,atoms%ekin,temp
        etotold=etot
        
    enddo !end of loop over imd
    close(1000)
    call final_potential_forces(parini,atoms)
end subroutine md_nvt_nose_hoover_cp
!*****************************************************************************************
subroutine md_nvt_nose_hoover_chain(parini,atoms)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_potential, only: potential, perfstatus
    use mod_atoms, only: typ_atoms, typ_file_info
    use mod_dynamics, only: dt, nmd
    use mod_processors, only: iproc
    use mod_potential, only: bias 
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variables
    type(typ_atoms):: atoms
    type(typ_file_info):: file_info
    integer:: iat, ierr, nat_t, i
    integer:: imd, ff,ntherm
    real(8):: etot, epotold, etotold
    real(8):: DNRM2, fnrm, t1, aboltzmann, totmass, temp, etotavg
    real(8):: t2,t3,t4 
    real(8):: rl, ru 
    real(8):: rcm(3), vcm(3)
    real(8):: scale_vat, temp_trget, ekin_target
    real(8) :: sum1, sum2, sum3
    real(8) :: kt, temp_prev, tol, tolerance 
    character(56):: comment
    real(8):: forces_nosehoover(3,atoms%nat)
    real(8):: rat_next(3,atoms%nat), vat_old(3,atoms%nat)
    real(8):: omega
    real(8), allocatable ::zeta_next(:) ,zeta(:) ,zeta_prev(:) ,dzeta(:),mass_q(:),dzeta_old(:)
    real(8):: rat_init(3,atoms%nat)
    real(8):: r, dx(3) , rsq, msd1, msd2, msd3

    call random_seed() 
    rat_init=atoms%rat

    call init_potential_forces(parini,atoms)

    open(unit=1000,file="velocity",status='replace')
    open(unit=1111,file="displace.txt",status='replace')
    file_info%filename_positions='posout.acf'
    file_info%file_position='new'
    file_info%print_force=parini%print_force_dynamics
    call acf_write(file_info,atoms=atoms,strkey='posout')

    !  ___________parameters_______________________________________
    ntherm=3
    aboltzmann= 3.1668113916289087d-6
    temp_trget = parini%temp_dynamics
    kt = aboltzmann*temp_trget
    tolerance = 1.d-4
    ekin_target=1.5d0*atoms%nat*aboltzmann*parini%init_temp_dynamics

    allocate(zeta_next(ntherm), zeta(ntherm),&
             zeta_prev(ntherm), dzeta(ntherm),&
             mass_q(ntherm), dzeta_old(ntherm))
    zeta      = 0.d0
    zeta_next = 0.d0
    zeta_prev = 0.d0
    dzeta     = 0.d0
    !omega = 1.d0/40
    mass_q   = kt*dt**2*400
    mass_q(1)= 3.d0*atoms%nat*kt*dt**2*400

    call get_atomic_mass(atoms,totmass)
    !_______________________initial velocity __________________________
    if (parini%restart_dynamics )then
        open(unit=1001,file="velocity_r",status='old')
        read(1001,*)
        read(1001,*) 
        read(1001,*)
        do iat=1,atoms%nat
            read(1001,*) atoms%vat(1,iat),atoms%vat(2,iat),atoms%vat(3,iat)
        enddo
    else
        if ( parini%init_temp_dynamics==0.d0) then
            atoms%vat(:,:)=0.d0
        else
            call set_velocities(atoms, ekin_target)
        endif
    endif
    !____________________________________________________________________
    epotold=atoms%epot
    call cal_potential_forces(parini,atoms)
    if (trim(bias)=='yes')   call plane_repulsion(atoms)

    call ekin_temprature(atoms,temp,vcm,rcm,totmass) 
    etot=atoms%epot+atoms%ekin
    etotold=etot

    write(*,'(a,2e20.10)') 'epotold,epot',epotold,atoms%epot
    write(21,'(i9,4es25.15)') 0,etot,atoms%epot,atoms%ekin,temp
    write(22,'(i9,6es20.10)') 0,rcm(1:3),vcm(1:3)
   !_________________________ The first md step _________________________
    imd=1
    t1=0.5*dt*dt
    do iat=1,atoms%nat
        forces_nosehoover(1:3,iat)=atoms%fat(1:3,iat)-atoms%amass(iat)*atoms%vat(1:3,iat)*dzeta(1)
    enddo
    call thermostat_evolution_2(atoms,zeta_next,zeta,zeta_prev,dzeta,mass_q,kt,ntherm,imd)
    do iat=1,atoms%nat
        rat_next(1:3,iat)=atoms%rat(1:3,iat) + dt*atoms%vat(1:3,iat) + t1*forces_nosehoover(1:3,iat)/atoms%amass(iat) 
        atoms%vat(1:3,iat)=(rat_next(1:3,iat)-atoms%rat(1:3,iat))/dt
        atoms%rat(1:3,iat)=rat_next(1:3,iat)
    enddo
    call back_to_cell(atoms)
    dzeta=(zeta_next-zeta)/dt
    zeta_prev=zeta
    zeta=zeta_next
    !_____________________________________________________________________
    do imd=2,nmd
        parini%time_dynamics = (imd-1)*dt
    write(110,*) imd,dzeta
        vat_old =atoms%vat
        dzeta_old=dzeta

        epotold=atoms%epot
        call cal_potential_forces(parini,atoms)
        if (trim(bias)=='yes') call plane_repulsion(atoms)
        call ekin_temprature(atoms,temp,vcm,rcm,totmass) 
        etot=atoms%epot+atoms%ekin
        etotold=etot

        write(21,'(i9,4es25.15)') imd-1,etot,atoms%epot,atoms%ekin,temp
        write(22,'(i9,6es20.10)') imd-1,rcm(1:3),vcm(1:3)

        do iat=1,atoms%nat
            forces_nosehoover(1:3,iat)=atoms%fat(1:3,iat)-atoms%amass(iat)*atoms%vat(1:3,iat)*dzeta(1)
            rat_next(1:3,iat)= atoms%rat(1:3,iat)+atoms%vat(1:3,iat)*dt &
            &       + dt*dt*forces_nosehoover(1:3,iat)/atoms%amass(iat)
        enddo
        call thermostat_evolution_2(atoms,zeta_next,zeta,zeta_prev,dzeta,mass_q,kt,ntherm,imd)
        call ekin_temprature(atoms,temp_prev,vcm,rcm,totmass) 
        tol=1.d0
       ! do while (tol>tolerance)
       !     atoms%vat = 0.5*((rat_next-atoms%rat)/dt+vat_old)
       !     dzeta=0.5*((zeta_next - zeta )/dt+dzeta_old)
       !     call ekin_temprature(atoms,temp,vcm,rcm,totmass) 
       !     do iat=1,atoms%nat
       !         forces_nosehoover(1:3,iat)=atoms%fat(1:3,iat)-atoms%amass(iat)*atoms%vat(1:3,iat)*dzeta(1)
       !         rat_next(1:3,iat)= atoms%rat(1:3,iat)+atoms%vat(1:3,iat)*dt &
       !         &       + dt*dt*forces_nosehoover(1:3,iat)/atoms%amass(iat)
       !     enddo
       !     call thermostat_evolution_2(atoms,zeta_next,zeta,zeta_prev,dzeta,mass_q,kt,ntherm,imd)
       !     tol=dabs(temp-temp_prev)/dabs(temp_prev)
       !     temp_prev=temp
       ! enddo

       ! atoms%vat =1.5d0*(rat_next - atoms%rat )/dt -0.5d0* vat_old
        atoms%vat = (rat_next - atoms%rat )/dt
        atoms%rat = rat_next
       ! dzeta = 0.5d0*(3.d0*zeta_next-4.d0*zeta+zeta_prev)/dt 
        dzeta =(zeta_next-zeta)/dt 
        zeta_prev=zeta
        zeta=zeta_next
        etotold=etot

        call back_to_cell(atoms)
        if(mod(imd,1)==0) then
            file_info%file_position='append'
            call acf_write(file_info,atoms=atoms,strkey='posout')
            write(1111,*) '#'
            write(1111,*) '#    imd = ',imd, parini%time_dynamics
            write(1111,*) '#'
            msd1= 0.d0
            msd2= 0.d0
            msd3= 0.d0
            do iat=1,atoms%nat

                dx(1:3)=atoms%rat(:,iat)-rat_init(:,iat)
                rsq=(dx(1)**2+dx(2)**2+dx(3)**2)
                r=sqrt(dx(1)**2+dx(2)**2+dx(3)**2)
                msd1 = msd1 + rsq                !all directions
                msd2 = msd2 + dx(1)**2+dx(2)**2  !x,y directions
                msd3 = msd3 + dx(3)**2           !z   direction
                write(1111,'(i5,2a5,4es25.17)')iat," ",atoms%sat(iat), dx, r
            enddo
            write(1111,*) '  MSD    = ', msd1/atoms%nat 
            write(1111,*) '  MSD_xy = ', msd2/atoms%nat 
            write(1111,*) '  MSD_z  = ', msd3/atoms%nat 
            write(1111,*) '# -----------------------------------------------'
         !   write(1000,*) '#'
         !   write(1000,*) '#    imd = ', imd
         !   write(1000,*) '#'
         !   do iat=1,atoms%nat
         !       write(1000,'(3es25.17)') atoms%vat(1,iat),atoms%vat(2,iat),atoms%vat(3,iat)
         !   enddo
        
        endif
        etotold=etot
    enddo !end of loop over imd
    close(1000)
    call final_potential_forces(parini,atoms)
end subroutine md_nvt_nose_hoover_chain
!*****************************************************************************************
subroutine set_langevin_randforce(eta,nat)
    implicit none
    integer :: nat
    real(8) ::eta(3,nat), sum1, sum2, sum3

    call gausdist_alborz(nat,eta(1,:))
    call gausdist_alborz(nat,eta(2,:))
    call gausdist_alborz(nat,eta(3,:))
    sum1=sum(eta(1,:))
    sum2=sum(eta(2,:))
    sum3=sum(eta(2,:))
    eta(1,:)= eta(1,:)- sum1/nat
    eta(2,:)= eta(2,:)- sum2/nat
    eta(3,:)= eta(3,:)- sum3/nat

end subroutine set_langevin_randforce
!*****************************************************************************************
subroutine back_to_cell(atoms)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms):: atoms
    integer :: iat

        if (trim(atoms%boundcond)=='free') then
            continue
        else if (trim(atoms%boundcond)=='bulk') then
            do iat=1,atoms%nat
                atoms%rat(1,iat)=modulo(modulo(atoms%rat(1,iat),atoms%cellvec(1,1)),atoms%cellvec(1,1))
                atoms%rat(2,iat)=modulo(modulo(atoms%rat(2,iat),atoms%cellvec(2,2)),atoms%cellvec(2,2))
                atoms%rat(3,iat)=modulo(modulo(atoms%rat(2,iat),atoms%cellvec(3,3)),atoms%cellvec(3,3))
            enddo
        else if (trim(atoms%boundcond)=='slab') then
            do iat=1,atoms%nat
                atoms%rat(1,iat)=modulo(modulo(atoms%rat(1,iat),atoms%cellvec(1,1)),atoms%cellvec(1,1))
                atoms%rat(2,iat)=modulo(modulo(atoms%rat(2,iat),atoms%cellvec(2,2)),atoms%cellvec(2,2))
            enddo
        else
            write(*,*)"periodic BC is just modified for slab "
            stop
        endif

end subroutine back_to_cell
!***********************************************************************************************
subroutine plane_repulsion(atoms)
    use mod_atoms, only: typ_atoms
    implicit none
    integer :: iat
    real(8) rl, ru,t1,t2
    type(typ_atoms):: atoms
    do iat=1,atoms%nat
        rl=atoms%rat(3,iat)
        if (rl<7.d0) then
            t1=100.d0*exp(-2.d0*rl)
            atoms%fat(3,iat)=atoms%fat(3,iat)+t1
            atoms%epot = atoms%epot+t1/3.d0
        endif
        ru=atoms%cellvec(3,3)-atoms%rat(3,iat)
        if (ru<7.d0) then
            t2=100.d0*exp(-2.d0*ru)
            atoms%fat(3,iat)=atoms%fat(3,iat)-t2
            atoms%epot = atoms%epot+t2/3.d0
        endif
    enddo
end subroutine plane_repulsion
!*************************************************************************************************************
subroutine thermostat_evolution(atoms,zeta_next,zeta,zeta_prev,dzeta,mass_q,kt,ntherm,imd)
    use mod_interface
    use mod_atoms, only: typ_atoms, typ_file_info
    use mod_dynamics, only: dt, nmd
    implicit none
    type(typ_atoms):: atoms
    integer :: ntherm, imd 
    integer :: i, iat, ith 
    real(8) :: kt, t1
    real(8) :: zeta_next(3,atoms%nat,ntherm), zeta(3,atoms%nat,ntherm),zeta_prev(3,atoms%nat,ntherm)
    real(8) :: dzeta(3,atoms%nat,ntherm), mass_q(ntherm)
    real(8) :: force_therm(3,atoms%nat,ntherm)
    
    do ith=1,ntherm
    if (ith==1) then
        do iat=1,atoms%nat
        do i=1,3
            force_therm(i,iat,1)=atoms%amass(iat)*atoms%vat(i,iat)**2-kt
        enddo
        enddo
    else if (ith==ntherm-1) then
        do iat=1,atoms%nat
        do i=1,3
            force_therm(i,iat,ith)=(mass_q(ntherm-2)*(dzeta(i,iat,ntherm-2))**2-kt)+&
                                   (mass_q(ntherm  )*(dzeta(i,iat,ntherm  ))**2-kt)
        enddo
        enddo
    else if (ith==ntherm) then
        do iat=1,atoms%nat
        do i=1,3
            force_therm(i,iat,ntherm)=mass_q(ntherm-1)*(dzeta(i,iat,ntherm-1)**2)-kt
        enddo
        enddo
    else    
        do iat=1,atoms%nat
        do i=1,3
            force_therm(i,iat,ith)=mass_q(ith-1)*(dzeta(i,iat,ith-1)**2)-kt
        enddo        
        enddo        
    endif
    enddo        

    t1=dt**2
    if (imd==1) then  !taylor
        do ith=1,ntherm
        if (ith==1) then
            do iat=1,atoms%nat
            do i=1,3
                zeta_next(i,iat,1)=zeta(i,iat,1)+dzeta(i,iat,1)*dt+&
                                (force_therm(i,iat,1)-mass_q(1)*dzeta(i,iat,1)*dzeta(i,iat,2))*t1/(2.d0*mass_q(1))
            enddo
            enddo
        else if (ith==ntherm-1)then
            do iat=1,atoms%nat
            do i=1,3
                zeta_next(i,iat,ith)=zeta(i,iat,ith)+dzeta(i,iat,ith)*dt+&
                                (force_therm(i,iat,ith)-mass_q(ith)*dzeta(i,iat,ith)*dzeta(i,iat,ntherm))*t1/(2.d0*mass_q(ith))
            enddo
            enddo
        else if (ith==ntherm) then
            do iat=1,atoms%nat
            do i=1,3
                zeta_next(i,iat,ntherm)=zeta(i,iat,ntherm)+dzeta(i,iat,ntherm)*dt+&
                                (force_therm(i,iat,ntherm)-mass_q(ntherm)*dzeta(i,iat,ntherm)*dzeta(i,iat,ntherm-1))*t1/(2.d0*mass_q(ntherm))
            enddo
            enddo
        else
            do iat=1,atoms%nat
            do i=1,3
                zeta_next(i,iat,ith)=zeta(i,iat,ith)+dzeta(i,iat,ith)*dt+&
                                (force_therm(i,iat,ith)-mass_q(ith)*dzeta(i,iat,ith)*dzeta(i,iat,ith+1))*t1/(2.d0*mass_q(ith))
            enddo
            enddo
        endif
        enddo
    else   ! verlet
        do ith=1,ntherm
        if (ith==1) then
            do iat=1,atoms%nat
            do i=1,3
                zeta_next(i,iat,1)=zeta(i,iat,1)+dzeta(i,iat,1)*dt+&
                                (force_therm(i,iat,1)-mass_q(1)*dzeta(i,iat,1)*dzeta(i,iat,2))*t1/(mass_q(1))
            enddo
            enddo
        else if (ith==ntherm-1)then
            do iat=1,atoms%nat
            do i=1,3
                zeta_next(i,iat,ith)=zeta(i,iat,ith)+dzeta(i,iat,ith)*dt+&
                                (force_therm(i,iat,ith)-mass_q(ith)*dzeta(i,iat,ith)*dzeta(i,iat,ntherm))*t1/(mass_q(ith))
            enddo
            enddo
        else if (ith==ntherm) then
            do iat=1,atoms%nat
            do i=1,3
                zeta_next(i,iat,ntherm)=zeta(i,iat,ntherm)+dzeta(i,iat,ntherm)*dt+&
                                (force_therm(i,iat,ntherm)-mass_q(ntherm)*dzeta(i,iat,ntherm)*dzeta(i,iat,ntherm-1))*t1/(mass_q(ntherm))
            enddo
            enddo
        else
            do iat=1,atoms%nat
            do i=1,3
                zeta_next(i,iat,ith)=zeta(i,iat,ith)+dzeta(i,iat,ith)+&
                                (force_therm(i,iat,ith)-mass_q(ith)*dzeta(i,iat,ith)*dzeta(i,iat,ith+1))*t1/(mass_q(ith))
            enddo
            enddo
        endif
        enddo
    endif
        !write(110,*)"************************************************************"
        !write(110,*)imd,"force_therm"
        !write(110,*)force_therm(:,19,1)
        !write(110,*)imd,"zeta"
        !write(110,*)zeta(:,19,1)
        !write(110,*)zeta_next(:,19,1)

end subroutine thermostat_evolution
!*************************************************************************************************************
subroutine thermostat_evolution_2(atoms,zeta_next,zeta,zeta_prev,dzeta,mass_q,kt,ntherm,imd)
    use mod_interface
    use mod_atoms, only: typ_atoms, typ_file_info
    use mod_dynamics, only: dt, nmd
    implicit none
    type(typ_atoms):: atoms
    integer :: ntherm, imd 
    integer :: i, iat, ith 
    real(8) :: kt, t1, tt
    real(8) :: zeta_next(ntherm), zeta(ntherm),zeta_prev(ntherm)
    real(8) :: dzeta(ntherm), mass_q(ntherm)
    real(8) :: force_therm(ntherm)
    
    tt=0.d0
    do iat=1,atoms%nat
        t1=atoms%amass(iat)
        tt=tt+t1*(atoms%vat(1,iat)**2+atoms%vat(2,iat)**2+atoms%vat(3,iat)**2)
    enddo

    if (ntherm==1) then
        force_therm(1)=tt-3*atoms%nat*kt
    else
        do ith=1,ntherm
            if (ith==1) then
                force_therm(1)=tt-3*atoms%nat*kt- mass_q(1)*dzeta(1)*dzeta(2)
            elseif (ith==ntherm) then
                force_therm(ntherm)=mass_q(ntherm-1)*(dzeta(ntherm-1)**2) - kt
            else
                force_therm(ith)=mass_q(ith-1)*(dzeta(ith-1)**2) - kt- mass_q(ith)*dzeta(ith)*dzeta(ith+1)
            endif
        enddo        
    endif

    t1=dt**2
    if (imd==1) then
        do ith=1,ntherm
            zeta_next(ith)=zeta(ith)+dzeta(ith)*dt+&
                        (force_therm(ith))*t1/(2.d0*mass_q(ith))
        enddo
    else
        do ith=1,ntherm
            zeta_next(ith)=zeta(ith)+dzeta(ith)*dt+&
                        (force_therm(ith))*t1/(mass_q(ith))
        enddo
    endif

end subroutine thermostat_evolution_2
!************************************************************************
subroutine md_nvt_nose_hoover(parini,atoms)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_potential, only: potential, perfstatus
    use mod_atoms, only: typ_atoms, typ_file_info
    use mod_dynamics, only: dt, nmd
    use mod_processors, only: iproc
    use mod_potential, only: bias 
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variables
    type(typ_atoms):: atoms
    type(typ_file_info):: file_info
    integer:: iat, ierr, nat_t, i
    integer:: imd, ff
    real(8):: etot, epotold, etotold
    real(8):: DNRM2, fnrm, t1, aboltzmann, totmass, temp, etotavg
    real(8):: t2,t3,t4 
    real(8):: rl, ru 
    real(8):: rcm(3), vcm(3)
    real(8):: scale_vat, temp_trget, ekin_target
    real(8):: eta(3,atoms%nat)
    real(8):: gama, factor, factor_iat
    real(8) :: sum1, sum2, sum3
    real(8) :: kt, temp_prev, tol, tolerance 
    character(56):: comment
    real(8)::  forces_nose(3,atoms%nat)
    real(8):: rat_next(3,atoms%nat), vat_old(3,atoms%nat)
    real(8):: zeta, dzeta, mass_q, force_q, zeta_next, dzeta_old 
    real(8):: rat_init(3,atoms%nat)
    real(8):: r, dx(3) , rsq, msd1, msd2, msd3

    call random_seed() 
    rat_init=atoms%rat

    file_info%filename_positions='posout.acf'
    file_info%file_position='new'
    file_info%print_force=parini%print_force_dynamics
    call acf_write(file_info,atoms=atoms,strkey='posout')
    open(unit=1000,file="velocity",status='replace')
    open(unit=1111,file="displace.txt",status='replace')
    call init_potential_forces(parini,atoms)

    !  ___________parameters_______________________________________
    aboltzmann= 3.1668113916289087d-6
    temp_trget = parini%temp_dynamics
    kt = aboltzmann*temp_trget
    tolerance = 1.d-4
    mass_q= (3*atoms%nat)*kt*(80*dt)**2
    ekin_target=1.5d0*atoms%nat*aboltzmann*parini%init_temp_dynamics
    dzeta = 0.d0   ;    zeta  = 0.d0

    call get_atomic_mass(atoms,totmass)
    !_______________________initial velocity __________________________
    if (parini%restart_dynamics )then
        open(unit=1001,file="velocity_r",status='old')
        read(1001,*)
        read(1001,*) 
        read(1001,*)
        do iat=1,atoms%nat
            read(1001,*) atoms%vat(1,iat),atoms%vat(2,iat),atoms%vat(3,iat)
        enddo
    else
        if ( parini%init_temp_dynamics==0.d0) then
            atoms%vat(:,:)=0.d0
        else
            call set_velocities(atoms, ekin_target)
        endif
    endif
    !____________________________________________________________________

    epotold=atoms%epot
    call cal_potential_forces(parini,atoms)
    if (trim(bias)=='yes')  call plane_repulsion(atoms)

    call ekin_temprature(atoms,temp,vcm,rcm,totmass) 
    etot=atoms%epot+atoms%ekin
    etotold=etot
    write(*,'(a,2e20.10)') 'epotold,epot',epotold,atoms%epot
    write(21,'(i9,4es25.15)') 0,etot,atoms%epot,atoms%ekin,temp
    write(22,'(i9,6es20.10)') 0,rcm(1:3),vcm(1:3)
   !_________________________ The first md step _________________________
    if (parini%restart_dynamics )then
        t1=dt*dt
    else
        t1=0.5*dt*dt
    endif
    force_q = 0.d0
    do iat=1,atoms%nat
        forces_nose(1:3,iat)=atoms%fat(1:3,iat)-dzeta*atoms%vat(1:3,iat)*atoms%amass(iat)
        force_q = force_q + atoms%amass(iat)*(atoms%vat(1,iat)**2+atoms%vat(2,iat)**2+atoms%vat(3,iat)**2)
    enddo
        force_q = (force_q -(3*atoms%nat)*kt)!/mass_q
        zeta_next = zeta + t1* force_q/mass_q +dt*dzeta
    do iat=1,atoms%nat
        rat_next(1:3,iat)=atoms%rat(1:3,iat) + t1*forces_nose(1:3,iat)/atoms%amass(iat) + dt*atoms%vat(1:3,iat)
        atoms%vat(1:3,iat)=(rat_next(1:3,iat)-atoms%rat(1:3,iat))/dt
        atoms%rat(1:3,iat)=rat_next(1:3,iat)
    enddo
    dzeta = (zeta_next-zeta)/dt 
    call back_to_cell(atoms)
    !_____________________________________________________________________
    zeta=zeta_next

    do imd=2,nmd
        parini%time_dynamics = (imd-1)*dt

        vat_old =atoms%vat
        dzeta_old = dzeta
        parini%cal_charge = .false.
        epotold=atoms%epot
        call cal_potential_forces(parini,atoms)
        if (trim(bias)=='yes')  call plane_repulsion(atoms)
        call ekin_temprature(atoms,temp,vcm,rcm,totmass) 
        etot=atoms%epot+atoms%ekin
        etotold=etot
        write(21,'(i9,4es25.15)') imd-1,etot,atoms%epot,atoms%ekin,temp
        write(22,'(i9,6es20.10)') imd-1,rcm(1:3),vcm(1:3)

        force_q = 0.d0
        do iat=1,atoms%nat
            forces_nose(1:3,iat)=atoms%fat(1:3,iat)-dzeta*atoms%vat(1:3,iat)*atoms%amass(iat)
            rat_next(1:3,iat)= atoms%rat(1:3,iat)+atoms%vat(1:3,iat)*dt &
            &       + dt*dt*forces_nose(1:3,iat)/atoms%amass(iat)
            force_q = force_q + atoms%amass(iat)*(atoms%vat(1,iat)**2+atoms%vat(2,iat)**2+atoms%vat(3,iat)**2)
        enddo
        force_q = (force_q -(3*atoms%nat)*kt)
        zeta_next = zeta + dt*dt* force_q/mass_q +dt*dzeta
        tol=1.d0
        call ekin_temprature(atoms,temp_prev,vcm,rcm,totmass) 
        do while (tol>tolerance)
            atoms%vat = 0.5*((rat_next-atoms%rat)/dt+vat_old)
            dzeta     = 0.5*((zeta_next-zeta)/dt+dzeta_old)
            call ekin_temprature(atoms,temp,vcm,rcm,totmass) 
            force_q = 0.d0
            do iat=1,atoms%nat
                forces_nose(1:3,iat)=atoms%fat(1:3,iat)-dzeta*atoms%vat(1:3,iat)*atoms%amass(iat)
                rat_next(1:3,iat)= atoms%rat(1:3,iat)+ atoms%vat(1:3,iat)*dt &
                &       + dt*dt*forces_nose(1:3,iat)/atoms%amass(iat)
                force_q = force_q + atoms%amass(iat)*(atoms%vat(1,iat)**2+atoms%vat(2,iat)**2+atoms%vat(3,iat)**2)
            enddo
            force_q = (force_q -(3*atoms%nat)*kt)!/mass_q
            zeta_next = zeta + dt*dt* force_q/mass_q +dt*dzeta
            tol=dabs(temp-temp_prev)/dabs(temp_prev)
            temp_prev=temp
        enddo

        atoms%vat = (rat_next - atoms%rat )/dt
        atoms%rat =  rat_next
        dzeta = (zeta_next - zeta )/dt
        zeta = zeta_next
        etotold=etot
        call back_to_cell(atoms)
        if(mod(imd,100)==0) then
            file_info%file_position='append'
            call acf_write(file_info,atoms=atoms,strkey='posout')
            write(1111,*) '#'
            write(1111,*) '#    imd = ',imd, parini%time_dynamics
            write(1111,*) '#'
            msd1= 0.d0
            msd2= 0.d0
            msd3= 0.d0
            do iat=1,atoms%nat

                dx(1:3)=atoms%rat(:,iat)-rat_init(:,iat)
                rsq=(dx(1)**2+dx(2)**2+dx(3)**2)
                r=sqrt(dx(1)**2+dx(2)**2+dx(3)**2)
                msd1 = msd1 + rsq                !all directions
                msd2 = msd2 + dx(1)**2+dx(2)**2  !x,y directions
                msd3 = msd3 + dx(3)**2           !z   direction
                write(1111,'(i5,2a5,4es25.17)')iat," ",atoms%sat(iat), dx, r
            enddo
            write(1111,*) '  MSD    = ', msd1/atoms%nat 
            write(1111,*) '  MSD_xy = ', msd2/atoms%nat 
            write(1111,*) '  MSD_z  = ', msd3/atoms%nat 
            write(1111,*) '# -----------------------------------------------'
        !    write(1000,*) '#'
        !    write(1000,*) '#    imd = ', imd
        !    write(1000,*) '#'
        !    do iat=1,atoms%nat
        !        write(1000,'(3es25.17)') atoms%vat(1,iat),atoms%vat(2,iat),atoms%vat(3,iat)
        !    enddo
        
        endif
    enddo !end of loop over imd
    close(1000)
    call final_potential_forces(parini,atoms)
end subroutine md_nvt_nose_hoover
!*****************************************************************************************
subroutine get_atomic_mass(atoms,totmass)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms):: atoms
    real(8):: totmass
    integer:: iat
    totmass=0.d0
    do iat=1,atoms%nat
        if(trim(atoms%sat(iat))=='Na') then 
            atoms%amass(iat)= 22.98976928*1822.888485540949556d0
        else if(trim(atoms%sat(iat))=='Cl') then
            atoms%amass(iat)= 35.45*1822.888485540949556d0
        else if(trim(atoms%sat(iat))=='Si') then
            atoms%amass(iat)= 28.085*1822.888485540949556d0
        else if(trim(atoms%sat(iat))=='Zr') then
            atoms%amass(iat)= 91.224*1822.888485540949556d0
        else if(trim(atoms%sat(iat))=='Y') then
            atoms%amass(iat)= 88.90585*1822.888485540949556d0
        else if(trim(atoms%sat(iat))=='O') then
            atoms%amass(iat)= 15.9994*1822.888485540949556d0
        else if(trim(atoms%sat(iat))=='K') then 
            atoms%amass(iat)= 22.98976928*1822.888485540949556d0
        else if(trim(atoms%sat(iat))=='Br') then
            atoms%amass(iat)= 35.45*1822.888485540949556d0
        else
            write(*,*)"unknown atomic type"
            stop
        endif
        totmass=totmass+atoms%amass(iat)
    enddo
end subroutine get_atomic_mass
!*****************************************************************************************
