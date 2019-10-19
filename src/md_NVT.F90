!*****************************************************************************************
subroutine md_nvt_langevin(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_potential, only: potential, perfstatus
    use mod_atoms, only: typ_atoms, typ_file_info
    use mod_atoms, only: get_rat, update_ratp, update_rat, set_rat
    use mod_acf, only: acf_write
    use mod_velocity, only: set_velocities
    use mod_dynamics, only: dt, nmd ,nfreq
    use mod_processors, only: iproc
    !use mod_potential, only: bias 
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variables
    type(typ_atoms):: atoms
    type(typ_file_info):: file_info
    integer:: iat, ierr, nat_t, i
    integer:: imd, ff , rmd
    real(8):: etot, epotold, etotold
    real(8):: DNRM2, fnrm, t1, aboltzmann, totmass, temp, etotavg
    real(8):: t2,t3,t4 ,tt
    real(8):: rl, ru 
    real(8):: rcm(3), vcm(3)
    real(8):: scale_vat, temp_trget, ekin_target
    real(8):: eta(3,atoms%nat)
    real(8):: gama, factor, factor_iat
    real(8) :: sum1, sum2, sum3
    real(8) :: kt, temp_prev, tol, tolerance 
    character(56):: comment , nn
    real(8):: langev(atoms%nat), forces_langevin(3,atoms%nat)
    real(8):: rat_next(3,atoms%nat), vat_old(3,atoms%nat)
    real(8):: rat_init(3,atoms%nat)
    real(8):: r, dx(3) , rsq, msd1, msd2, msd3,pi,dipole
    pi=4.d0*atan(1.d0)

    call random_seed() 
    call get_rat(atoms,rat_init)
    !  ___________parameters_______________________________________
    gama=1.d-3
    aboltzmann= 3.1668139952584056d-06
    temp_trget = parini%temp_dynamics
    kt = aboltzmann*temp_trget
    tolerance = 1.d-2

    if (parini%restart_dynamics )then
        open(unit=21,file="md_out.dat",status='old',Access = 'append')
        open(unit=1111,file="displacement.dat",status='old',Access = 'append')
    else
        open(unit=21,file="md_out.dat",status='replace')
        open(unit=1111,file="displacement.dat",status='replace')
        write(21,'(a9,4a25)') "imd","E_tot",'E_pot','E_kin','Temp'
    endif
    call init_potential_forces(parini,atoms)
    call get_atomic_mass(atoms,totmass)
    langev(:)=sqrt(2*gama*atoms%amass(:)*kt/dt)
    ekin_target=1.5d0*atoms%nat*aboltzmann*parini%init_temp_dynamics
    !_______________________initial velocity __________________________
    if (parini%restart_dynamics )then
        open(unit=1003,file="restart.dat",status='old')
        read(1003,*) nn  , rmd
        read(1003,*)  
        do iat=1,atoms%nat
            read(1003,'(3es25.17)') atoms%ratp(1,iat),atoms%ratp(2,iat),atoms%ratp(3,iat)
        enddo
        read(1003,*) 
        do iat=1,atoms%nat
            read(1003,'(3es25.17)') atoms%vat(1,iat),atoms%vat(2,iat),atoms%vat(3,iat)
        enddo
        close(1003)
        call update_rat(atoms,upall=.true.)
        call update_ratp(atoms)
        write(21,"(a,i8,a)") "#   *********************** restart from imd:",rmd,"  **************************"
    else
        rmd = 1
        if ( parini%init_temp_dynamics==0.d0) then
            atoms%vat(:,:)=0.d0
        else
            call set_velocities(atoms, ekin_target)
        endif
    endif
    !____________________________________________________________________
    epotold=atoms%epot
    call update_ratp(atoms)
    call cal_potential_forces(parini,atoms)

    write(*,'(a,2e20.10)') 'epotold,epot',epotold,atoms%epot
    file_info%filename_positions='trajectory.acf'
    if (parini%restart_dynamics )then
        file_info%file_position='append'
    else
        file_info%file_position='new'
        file_info%print_force=parini%print_force_dynamics
        call acf_write(file_info,atoms=atoms,strkey='trajectory')
    endif
    !____________________________________________________________________
    call set_langevin_randforce(eta,atoms%nat)
    if (parini%restart_dynamics )then
        t1=dt*dt
    else
        t1=0.5*dt*dt
    endif
    !____________________________________________________________________
    call ekin_temprature(atoms,temp,vcm,rcm,totmass) 
    etot=atoms%epot+atoms%ekin
    etotold=etot
   !_________________________ The first md step _________________________
    if (rmd==1) then
        write(21,'(i9,4es25.15)') 0,etot,atoms%epot,atoms%ekin,temp
       ! write(22,'(i9,6es20.10)') 0,rcm(1:3),vcm(1:3)
        call update_ratp(atoms)
        do iat=1,atoms%nat
            forces_langevin(1:3,iat)=atoms%fat(1:3,iat)+langev(iat)*eta(1:3,iat)-gama*atoms%amass(iat)*atoms%vat(1:3,iat)
            rat_next(1:3,iat)=atoms%ratp(1:3,iat) + t1*forces_langevin(1:3,iat)/atoms%amass(iat) + dt*atoms%vat(1:3,iat)
            atoms%vat(1:3,iat)=(rat_next(1:3,iat)-atoms%ratp(1:3,iat))/dt
            atoms%ratp(1:3,iat)=rat_next(1:3,iat)
        enddo
        call update_rat(atoms,upall=.true.)
    endif
    !_____________________________________________________________________
    do imd=1+rmd,nmd+rmd
        parini%time_dynamics = (imd-1)*dt

        dipole=0.d0
        call update_ratp(atoms)
        do iat=1,atoms%nat
            dipole=dipole+atoms%qat(iat)*atoms%ratp(3,iat)
        enddo
        write(31,*) "dtime",parini%time_dynamics,parini%vu_ac_ewald*sin(2*pi*parini%frequency_ewald*parini%time_dynamics),dipole

        vat_old =atoms%vat
        parini%cal_charge = .false.
        epotold=atoms%epot
        call cal_potential_forces(parini,atoms)
        call ekin_temprature(atoms,temp,vcm,rcm,totmass) 
        etot=atoms%epot+atoms%ekin
        etotold=etot
        write(21,'(i9,4es25.15)') imd-1,etot,atoms%epot,atoms%ekin,temp
       ! write(22,'(i9,6es20.10)') imd-1,rcm(1:3),vcm(1:3)

        call set_langevin_randforce(eta,atoms%nat)
        call update_ratp(atoms)
        do iat=1,atoms%nat
            forces_langevin(1:3,iat)=atoms%fat(1:3,iat)+langev(iat)*eta(1:3,iat)-gama*atoms%amass(iat)*atoms%vat(1:3,iat)
            rat_next(1:3,iat)= atoms%ratp(1:3,iat)+atoms%vat(1:3,iat)*dt &
            &       + dt*dt*forces_langevin(1:3,iat)/atoms%amass(iat)
        enddo
        tol=1.d0
        call ekin_temprature(atoms,temp_prev,vcm,rcm,totmass) 
        do while (tol>tolerance)
            call update_ratp(atoms)
            atoms%vat = 0.5*((rat_next-atoms%ratp)/dt+vat_old)
            call ekin_temprature(atoms,temp,vcm,rcm,totmass) 
            do iat=1,atoms%nat
                forces_langevin(1:3,iat)=atoms%fat(1:3,iat)+langev(iat)*eta(1:3,iat)-gama*atoms%amass(iat)*atoms%vat(1:3,iat)
                call update_ratp(atoms)
                rat_next(1:3,iat)= atoms%ratp(1:3,iat)+ vat_old(1:3,iat)*dt &
                &       + dt*dt*forces_langevin(1:3,iat)/atoms%amass(iat)
            enddo
            tol=dabs(temp-temp_prev)/dabs(temp_prev)
            temp_prev=temp
        enddo

        call update_ratp(atoms)
        atoms%vat = (rat_next - atoms%ratp )/dt
        call set_rat(atoms,rat_next,setall=.true.)
        if(mod(imd,nfreq)==0) then
            file_info%file_position='append'
            call acf_write(file_info,atoms=atoms,strkey='trajectory')
            call cal_potential_forces(parini,atoms)

            open(unit=1003,file="restart.dat",status='replace')
                write(1003,*) "MD_STEP:" , imd
                write(1003,*) "ATOM:" 
                do iat=1,atoms%nat
                    write(1003,'(3es25.17)') atoms%ratp(1,iat),atoms%ratp(2,iat),atoms%ratp(3,iat)
                enddo
                write(1003,*) "VELOCITY:" 
                do iat=1,atoms%nat
                    write(1003,'(3es25.17)') atoms%vat(1,iat),atoms%vat(2,iat),atoms%vat(3,iat)
                enddo
            write(1003,*) "ENER:" ,  atoms%epot
            close(1003)



            call ekin_temprature(atoms,temp,vcm,rcm,totmass) 
            tt=atoms%epot+atoms%ekin


            write(1111,*) '#'
            write(1111,*) '#    imd = ',imd, parini%time_dynamics
            write(1111,*) '#'
            msd1= 0.d0
            msd2= 0.d0
            msd3= 0.d0
            call update_ratp(atoms)
            do iat=1,atoms%nat
                dx(1:3)=atoms%ratp(:,iat)-rat_init(:,iat)
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
    close(21)
    call final_potential_forces(parini,atoms)
end subroutine md_nvt_langevin
!*****************************************************************************************
subroutine md_nvt_nose_hoover_cp(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_potential, only: potential, perfstatus
    use mod_atoms, only: typ_atoms, typ_file_info, set_rat, get_rat, update_ratp
    use mod_velocity, only: set_velocities
    use mod_acf, only: acf_write
    use mod_dynamics, only: dt, nmd
    use mod_processors, only: iproc
    !use mod_potential, only: bias 
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
    real(8):: rat_next(3,atoms%nat), rat_prev(3,atoms%nat),vat_old(3,atoms%nat) 
    real(8), allocatable ::zeta_next(:,:,:) ,zeta(:,:,:) ,zeta_prev(:,:,:) ,dzeta(:,:,:),mass_q(:),dzeta_old(:,:,:)
    real(8):: rat_init(3,atoms%nat)
    real(8):: r, dx(3) , rsq, msd1, msd2, msd3

    !call random_seed() 
    call get_rat(atoms,rat_init)

    call init_potential_forces(parini,atoms)
    open(unit=1000,file="velocity",status='replace')
    open(unit=1111,file="displacement.dat",status='replace')
    file_info%filename_positions='trajectory.acf'
    file_info%file_position='new'
    file_info%print_force=parini%print_force_dynamics
    call acf_write(file_info,atoms=atoms,strkey='trajectory')

    !  ___________parameters_______________________________________
    ntherm=3
    aboltzmann= 3.1668139952584056d-06
    temp_trget = parini%temp_dynamics
    kt = aboltzmann*temp_trget
    tolerance = 1.d-9
!  HERE

    allocate(zeta_next(3,atoms%nat,ntherm), zeta(3,atoms%nat,ntherm),&
             zeta_prev(3,atoms%nat,ntherm), dzeta(3,atoms%nat,ntherm),&
             mass_q(ntherm),dzeta_old(3,atoms%nat,ntherm))
    zeta      = 0.d0
    zeta_next = 0.d0
    zeta_prev = 0.d0
    dzeta     = 0.d0
    mass_q    = 1.d0/kt
    !mass_q    =kt*((4000)**2)
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
    !if (trim(bias)=='yes')  call plane_repulsion(atoms)
    etot=atoms%epot+atoms%ekin
    etotold=etot
    do iat=1,atoms%nat
        forces_nosehoover(1:3,iat)=atoms%fat(1:3,iat)-atoms%amass(iat)*atoms%vat(1:3,iat)*dzeta(1:3,iat,1)
    enddo

    write(*,'(a,2e20.10)') 'epotold,epot',epotold,atoms%epot
    write(21,'(i9,4es25.15)') 0,etot,atoms%epot,atoms%ekin,temp
   ! write(22,'(i9,6es20.10)') 0,rcm(1:3),vcm(1:3)
   !_________________________ The first md step _________________________
    imd=1
    t1=0.5*dt*dt
    call update_ratp(atoms)
    do iat=1,atoms%nat
        rat_next(1:3,iat)=atoms%ratp(1:3,iat) + t1*forces_nosehoover(1:3,iat)/atoms%amass(iat) + dt*atoms%vat(1:3,iat)
    enddo
    do iat=1,atoms%nat
        atoms%vat(1:3,iat)=(rat_next(1:3,iat)-atoms%ratp(1:3,iat))/dt
    enddo
    call ekin_temprature(atoms,temp,vcm,rcm,totmass) 
    call thermostat_evolution(atoms,zeta_next,zeta,zeta_prev,dzeta,mass_q,kt,ntherm,imd)
    dzeta=(zeta_next-zeta)/dt
    zeta_prev = zeta
    zeta      = zeta_next
    call get_rat(atoms,rat_prev)
    call set_rat(atoms,rat_next,setall=.true.)
    !_____________________________________________________________________
    do imd=2,nmd
        parini%time_dynamics = (imd-1)*dt
        epotold=atoms%epot
        call cal_potential_forces(parini,atoms)
        !if (trim(bias)=='yes')  call plane_repulsion(atoms)

        call ekin_temprature(atoms,temp,vcm,rcm,totmass) 
        etot=atoms%epot+atoms%ekin
        etotold=etot
        write(21,'(i9,4es25.15)') imd-1,etot,atoms%epot,atoms%ekin,temp
      !  write(22,'(i9,6es20.10)') imd-1,rcm(1:3),vcm(1:3)
        call update_ratp(atoms)
        do iat=1,atoms%nat
            forces_nosehoover(1:3,iat)=atoms%fat(1:3,iat)-atoms%amass(iat)*atoms%vat(1:3,iat)*dzeta(1:3,iat,1)
            rat_next(1:3,iat)= 2.d0*atoms%ratp(1:3,iat)-rat_prev(1:3,iat) &
            &       + dt*dt*forces_nosehoover(1:3,iat)/atoms%amass(iat)
        enddo
        call thermostat_evolution(atoms,zeta_next,zeta,zeta_prev,dzeta,mass_q,kt,ntherm,imd)
        tol=1.d0
        call ekin_temprature(atoms,temp,vcm,rcm,totmass) 
        temp_prev = temp 
        do while (tol>tolerance)
            atoms%vat = 0.5*((rat_next  - rat_prev)/dt)
            dzeta     = 0.5*((zeta_next - zeta_prev)/dt)
            call ekin_temprature(atoms,temp,vcm,rcm,totmass) 
            call update_ratp(atoms)
            do iat=1,atoms%nat
                forces_nosehoover(1:3,iat)=atoms%fat(1:3,iat)-atoms%amass(iat)*atoms%vat(1:3,iat)*dzeta(1:3,iat,1)
                rat_next(1:3,iat)= 2.d0*atoms%ratp(1:3,iat)-rat_prev(1:3,iat) &
                &       + dt*dt*forces_nosehoover(1:3,iat)/atoms%amass(iat)
            enddo
            call thermostat_evolution(atoms,zeta_next,zeta,zeta_prev,dzeta,mass_q,kt,ntherm,imd)
            tol=dabs(temp-temp_prev)/dabs(temp_prev)
            temp_prev=temp
        enddo
        call update_ratp(atoms)
        atoms%vat = 0.5*(3.d0*rat_next-4.d0*atoms%ratp+rat_prev)/dt
        !atoms%vat = (rat_next-atoms%rat)/dt
        call get_rat(atoms,rat_prev)
        call set_rat(atoms,rat_next,setall=.true.)
        dzeta=0.5d0*(3.d0*zeta_next-4.d0*zeta+zeta_prev)/dt
        zeta_prev=zeta
        zeta=zeta_next
       ! call back_to_cell(atoms)
        if(mod(imd,100)==0) then
            file_info%file_position='append'
            call acf_write(file_info,atoms=atoms,strkey='trajectory')
            write(1111,*) '#'
            write(1111,*) '#    imd = ',imd, parini%time_dynamics
            write(1111,*) '#'
            msd1= 0.d0
            msd2= 0.d0
            msd3= 0.d0
            call update_ratp(atoms)
            do iat=1,atoms%nat
                dx(1:3)=atoms%ratp(:,iat)-rat_init(:,iat)
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
        write(221,'(i9,4es25.15)') imd+1,etot,atoms%epot,atoms%ekin,temp
        etotold=etot
        
    enddo !end of loop over imd
    close(1000)
    call final_potential_forces(parini,atoms)
end subroutine md_nvt_nose_hoover_cp
!*****************************************************************************************
subroutine md_nvt_nose_hoover_chain(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_potential, only: potential, perfstatus
    use mod_atoms, only: typ_atoms, typ_file_info, get_rat, update_ratp, update_rat
    use mod_velocity, only: set_velocities
    use mod_acf, only: acf_write
    use mod_dynamics, only: dt, nmd, nfreq
    use mod_processors, only: iproc
    !use mod_potential, only: bias 
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variables
    type(typ_atoms):: atoms
    type(typ_file_info):: file_info
    integer:: iat, ierr, nat_t, i, j
    integer:: imd, ff,ntherm, ith, rmd
    real(8):: etot, epotold, etotold
    real(8):: DNRM2, fnrm, t1, aboltzmann, totmass, temp, etotavg
    real(8):: t2,t3,t4 
    real(8):: rl, ru 
    real(8):: rcm(3), vcm(3)
    real(8):: drcm(3), rcm_init(3)
    real(8):: scale_vat, temp_trget, ekin_target
    real(8) :: sum1, sum2, sum3
    real(8) :: kt, temp_prev, tol, tolerance 
    character(56):: comment , nn
    real(8):: omega, tt
    real(8), allocatable :: zeta(:), dzeta(:), mass_q(:), azeta(:)
    real(8):: rat_init(3,atoms%nat)
    real(8):: r, dx(3) , rsq, msd1, msd2, msd3
    real(8):: nof, enhc 
    real(8):: dt2, dt4, dt8 
    integer:: vfile
    real(8):: temp1, temp2
    real(8):: sumf1, sumf2, sumf3
    real(8):: tau, time_unit
    logical:: lfist= .true.
    call random_seed() 
    call get_rat(atoms,rat_init)
    time_unit=41.341373336
!   dt=parini%dt_dynamics

    call init_potential_forces(parini,atoms)

    if (parini%restart_dynamics )then
        open(unit=21,file="md_out.dat",status='old',Access = 'append')
        open(unit=1111,file="displacement.dat",status='old',Access = 'append')
        open(unit=1112,file="MSD.dat",status='old',Access = 'append')
    else
        open(unit=21,file="md_out.dat",status='replace')
        write(21,'(a9,4a25)') "imd","E_tot",'E_pot','E_kin','Temp'
        open(unit=1111,file="displacement.dat",status='replace')
        open(unit=1112,file="MSD.dat",status='replace')
        write(1112,'(a15,3a25)') "imd " , " MSD       "  , " MSD_xy       " , " MSD_z       " 
    endif

    file_info%filename_positions='trajectory.acf'
    if (parini%restart_dynamics )then
        file_info%file_position='append'
    else
        file_info%file_position='new'
        file_info%print_force=parini%print_force_dynamics
        call acf_write(file_info,atoms=atoms,strkey='trajectory')
    endif

    !  ___________parameters_______________________________________
    ntherm = parini%ntherm
    if (ntherm < 2) then
        stop "ERRROR: nose-hoover chain works with more than 2 thermostat"
    endif
    aboltzmann= 3.1668139952584056d-06
    temp_trget = parini%temp_dynamics
    omega = parini%highest_frequency ! THz
    !omega = 20.d0 !THz

    dt2 = 0.5d0*dt
    dt4 = 0.5d0*dt2
    dt8 = 0.5d0*dt4

    kt = aboltzmann*temp_trget
    nof=0
    do iat=1,atoms%nat
        if(atoms%bemoved(1,iat)) nof=nof+1
        if(atoms%bemoved(2,iat)) nof=nof+1
        if(atoms%bemoved(3,iat)) nof=nof+1
    enddo
    !nof = (3.d0*atoms%nat)
    !nof = (3.d0*atoms%nat+ntherm)
    !ekin_target=0.5d0*nof*aboltzmann*parini%init_temp_dynamics
    ekin_target=1.5d0*atoms%nat*aboltzmann*parini%init_temp_dynamics

    allocate(zeta(ntherm), dzeta(ntherm),&
             mass_q(ntherm), azeta(ntherm))
    zeta  = 0.d0
    dzeta = 0.d0
    azeta = 0.d0
    tau = 1.d3/ omega *time_unit !Time in atomic unit

    !mass_q   = kt*tau**2
    mass_q   = kt*(dt/time_unit)**2*tau**2
    !mass_q(1)= 3.d0*atoms%nat*kt*tau**2
    !mass_q(1)= nof*kt*tau**2
    mass_q(1)= nof*kt*(dt/time_unit)**2*tau**2

    call get_atomic_mass(atoms,totmass)
    !_______________________initial velocity __________________________
    if (parini%restart_dynamics )then
        open(unit=1003,file="restart.dat",status='old')
        read(1003,*) nn  , rmd
        read(1003,*)  
        do iat=1,atoms%nat
            read(1003,'(3es25.17)') atoms%ratp(1,iat),atoms%ratp(2,iat),atoms%ratp(3,iat)
        enddo
        read(1003,*) 
        do iat=1,atoms%nat
            read(1003,'(3es25.17)') atoms%vat(1,iat),atoms%vat(2,iat),atoms%vat(3,iat)
        enddo
        read(1003,*) 
        do ith=1,ntherm
            read(1003,*) zeta(ith),dzeta(ith)
        enddo
        close(1003)
        call update_rat(atoms,upall=.true.)
        call update_ratp(atoms)
        write(21,"(a,i8,a)") "#   *********************** restart from imd:",rmd,"  **************************"
    else
        rmd = 0
        if ( parini%init_temp_dynamics==0.d0) then
            atoms%vat(:,:)=0.d0
        else
            call set_velocities(atoms, ekin_target)
        endif
    endif

    call ekin_temprature(atoms,temp,vcm,rcm,totmass) 
    temp=temp*(1.5d0*atoms%nat)/(0.5*nof)
    rcm_init = rcm 
    !____________________________________________________________________
    do imd=1+rmd,nmd+rmd
        parini%time_dynamics = (imd-1)*dt
        epotold=atoms%epot

        !if(parini%vflip_dynamics) then
        !    call update_ratp(atoms)
        !    do iat=1,atoms%nat
        !        if (atoms%cellvec(3,3)-atoms%ratp(3,iat) < 5.5d0 .and. atoms%vat(3,iat) > 0.d0)then
        !            atoms%vat(3,iat) = -atoms%vat(3,iat)
        !        endif
        !    enddo
        !    call update_ratp(atoms)
        !    do iat=1,atoms%nat
        !        if (atoms%ratp(3,iat) < 5.5d0 .and. atoms%vat(3,iat) < 0.d0)then
        !            atoms%vat(3,iat) = -atoms%vat(3,iat)
        !        endif
        !    enddo
        !endif
        call ekin_temprature(atoms,temp,vcm,rcm,totmass) 
        temp=temp*(1.5d0*atoms%nat)/(0.5*nof)
        call cal_potential_forces(parini,atoms)

        if(parini%fix_cm_dynamics) then
            sumf1 = sum(atoms%fat(1,:))/atoms%nat
            sumf2 = sum(atoms%fat(2,:))/atoms%nat
            sumf3 = sum(atoms%fat(3,:))/atoms%nat
            write(23,'(i5,3es15.4)') imd, sumf1*atoms%nat, sumf2*atoms%nat, sumf3*atoms%nat
            atoms%fat(1,:)=atoms%fat(1,:)-sumf1
            atoms%fat(2,:)=atoms%fat(2,:)-sumf2
           ! atoms%fat(3,:)=atoms%fat(3,:)-sumf3
        endif

        if(parini%wall_repulsion_dynamics) then
            call plane_repulsion(atoms)
        endif

        etot=atoms%epot+atoms%ekin
        etotold=etot
        enhc=atoms%epot+atoms%ekin+0.5*sum(dzeta**2*mass_q)+nof*kt*zeta(1)+sum(zeta*kt)

        write(*,'(a,2e20.10)') 'epotold,epot',epotold,atoms%epot
        write(21,'(i9,5es25.15)') imd,etot-atoms%ebattery,atoms%epot-atoms%ebattery,atoms%ekin,temp,enhc
       ! write(22,'(i9,6es20.10)') imd,rcm(1:3),vcm(1:3)

    !___________________  some steps temperature rescaling for pre_equilibrium  __________________

        if (imd<20 .and. (.not. parini%restart_dynamics)) then
            tt=(ekin_target/atoms%ekin)/(3.d0*atoms%nat)*nof
            atoms%vat =  atoms%vat*sqrt(tt)
        endif
! ___________________  chain _____________________________________________

        azeta(1) = (2.0*atoms%ekin-nof*kt)/mass_q(1);
        do ith = 2, ntherm
            azeta(ith) = (mass_q(ith-1)*dzeta(ith-1)*dzeta(ith-1)-kt);
        enddo
        dzeta(ntherm) = dzeta(ntherm) +azeta(ntherm) *dt4;
        do ith = 1, ntherm-1
            tt=exp(-dzeta(ntherm+1-ith)*dt8)
            dzeta(ntherm-ith) = (dzeta(ntherm-ith) *tt+azeta(ntherm-ith) *dt4) *tt
        enddo
        do ith = 1, ntherm
            zeta(ith)  = zeta(ith)  +dzeta(ith)*dt2;
        enddo
        tt = exp(-dzeta(1)*dt2);
        atoms%vat  =atoms%vat  * tt;
        atoms%ekin =atoms%ekin * tt*tt;
        azeta(1) = (2.0*atoms%ekin-nof*kt)/mass_q(1);
        do ith = 2, ntherm
            azeta(ith) = (mass_q(ith-1)*dzeta(ith-1)*dzeta(ith-1)-kt)/mass_q(ith);
        enddo
        do ith = 1, ntherm-1
            tt = exp(-dzeta(ith+1)*dt8)
            dzeta(ith) =(dzeta(ith) * tt + azeta(ith) *dt4)*  tt
            azeta(ith+1) = (mass_q(ith)*dzeta(ith)*dzeta(ith)-kt)/mass_q(2);
        enddo
        dzeta(ntherm) =dzeta(ntherm) + azeta(ntherm) *dt4;
! ___________________ end of  chain _________________________________________    

        call update_ratp(atoms)
        do iat=1,atoms%nat
            do j=1,3
                if (atoms%bemoved(j,iat)) then
                    atoms%ratp(j,iat) = atoms%ratp(j,iat) + atoms%vat(j,iat)*dt2
                else
                    atoms%vat(j,iat) = 0.d0
                    atoms%ratp(j,iat) = rat_init(j,iat) 
                endif
            enddo
        enddo
        call update_rat(atoms)

        call cal_potential_forces(parini,atoms)
        call update_ratp(atoms)
        do iat=1,atoms%nat
            do j=1,3
                if (atoms%bemoved(j,iat)) then
                    atoms%vat(j,iat) = atoms%vat(j,iat) + atoms%fat(j,iat) / atoms%amass(iat)*dt
                    atoms%ratp(j,iat) = atoms%ratp(j,iat) + atoms%vat(j,iat)*dt2
                else
                    atoms%vat(j,iat) = 0.d0
                    atoms%ratp(j,iat) = rat_init(j,iat) 
                endif
            enddo
        enddo
        call update_rat(atoms)
        call ekin_temprature(atoms,temp,vcm,rcm,totmass) 
        temp=temp*(1.5d0*atoms%nat)/(0.5*nof)

! ____________________  chain ______________________________________________
        azeta(1) = (2.0*atoms%ekin-nof*kt)/mass_q(1);
        do ith = 2, ntherm
            azeta(ith) = (mass_q(ith-1)*dzeta(ith-1)*dzeta(ith-1)-kt);
        enddo
        dzeta(ntherm) = dzeta(ntherm) +azeta(ntherm) *dt4;
        do ith = 1, ntherm-1
            tt=exp(-dzeta(ntherm+1-ith)*dt8)
            dzeta(ntherm-ith) = (dzeta(ntherm-ith) *tt+azeta(ntherm-ith) *dt4) *tt
        enddo
        do ith = 1, ntherm
            zeta(ith)  = zeta(ith)  +dzeta(ith)*dt2;
        enddo
        tt = exp(-dzeta(1)*dt2);
        atoms%vat  =atoms%vat  * tt;
        atoms%ekin =atoms%ekin * tt*tt;
        azeta(1) = (2.0*atoms%ekin-nof*kt)/mass_q(1);
        do ith = 2, ntherm
            azeta(ith) = (mass_q(ith-1)*dzeta(ith-1)*dzeta(ith-1)-kt)/mass_q(ith);
        enddo
        do ith = 1, ntherm-1
            tt = exp(-dzeta(ith+1)*dt8)
            dzeta(ith) =(dzeta(ith) * tt + azeta(ith) *dt4)*  tt
            azeta(ith+1) = (mass_q(ith)*dzeta(ith)*dzeta(ith)-kt)/mass_q(2);
        enddo
        dzeta(ntherm) =dzeta(ntherm) + azeta(ntherm) *dt4;

!____________________ write restart ________________________________________________   
     if(mod(imd,nfreq)==0) then
         file_info%file_position='append'
         call acf_write(file_info,atoms=atoms,strkey='trajectory')
         call cal_potential_forces(parini,atoms)

         open(unit=1003,file="restart.dat",status='replace')
         write(1003,*) "MD_STEP:" , imd
         write(1003,*) "ATOM:" 
         do iat=1,atoms%nat
             write(1003,'(3es25.17)') atoms%ratp(1,iat),atoms%ratp(2,iat),atoms%ratp(3,iat)
         enddo
         write(1003,*) "VELOCITY:" 
         do iat=1,atoms%nat
             write(1003,'(3es25.17)') atoms%vat(1,iat),atoms%vat(2,iat),atoms%vat(3,iat)
         enddo
         write(1003,*) "THERMOSTAT:" 
         do ith=1,ntherm
             write(1003,'(2es25.17)') zeta(ith),dzeta(ith)
         enddo
         write(1003,*) "ENER:" ,  atoms%epot
         close(1003)

         write(1111,*) '#'
         write(1111,*) '#    imd = ',imd
         write(1111,'(a5,2a5,4a25)') '#','','','dx','dy','dz','dr' 
         msd1= 0.d0
         msd2= 0.d0
         msd3= 0.d0
         call update_ratp(atoms)
         do iat=1,atoms%nat
             dx(1:3)=atoms%ratp(:,iat)-rat_init(:,iat)
             rsq=(dx(1)**2+dx(2)**2+dx(3)**2)
             r=sqrt(dx(1)**2+dx(2)**2+dx(3)**2)
             msd1 = msd1 + rsq                !all directions
             msd2 = msd2 + dx(1)**2+dx(2)**2  !x,y directions
             msd3 = msd3 + dx(3)**2           !z   direction
             write(1111,'(i5,2a5,4es25.17)')iat," ",atoms%sat(iat), dx, r
         enddo
         write(1112,'(i15,3es25.14)') imd-1 , msd1/atoms%nat , msd2/atoms%nat, msd3/atoms%nat 
     endif

        etotold=etot
    enddo !end of loop over imd
    close(21)

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
    use mod_atoms, only: typ_atoms, update_ratp, update_rat
    implicit none
    type(typ_atoms):: atoms
    integer :: iat

        if (trim(atoms%boundcond)=='free') then
            continue
        else if (trim(atoms%boundcond)=='bulk') then
            call update_ratp(atoms)
            do iat=1,atoms%nat
                atoms%ratp(1,iat)=modulo(modulo(atoms%ratp(1,iat),atoms%cellvec(1,1)),atoms%cellvec(1,1))
                atoms%ratp(2,iat)=modulo(modulo(atoms%ratp(2,iat),atoms%cellvec(2,2)),atoms%cellvec(2,2))
                atoms%ratp(3,iat)=modulo(modulo(atoms%ratp(2,iat),atoms%cellvec(3,3)),atoms%cellvec(3,3))
            enddo
            call update_rat(atoms,upall=.true.)
        else if (trim(atoms%boundcond)=='slab') then
            call update_ratp(atoms)
            do iat=1,atoms%nat
                atoms%ratp(1,iat)=modulo(modulo(atoms%ratp(1,iat),atoms%cellvec(1,1)),atoms%cellvec(1,1))
                atoms%ratp(2,iat)=modulo(modulo(atoms%ratp(2,iat),atoms%cellvec(2,2)),atoms%cellvec(2,2))
            enddo
            call update_rat(atoms,upall=.true.)
        else
            write(*,*)"periodic BC is just modified for slab "
            stop
        endif

end subroutine back_to_cell
!***********************************************************************************************
!It works just for repolsive of wall in z direction.
subroutine plane_repulsion(atoms)
    use mod_atoms, only: typ_atoms, update_ratp
    implicit none
    integer :: iat
    real(8) rl, ru,t1,t2
    type(typ_atoms):: atoms
    call update_ratp(atoms)
    do iat=1,atoms%nat
        rl=atoms%ratp(3,iat)
        if (rl<5.5d0) then
            t1=200.d0*exp(-1.5d0*rl)
            atoms%fat(3,iat)=atoms%fat(3,iat)+t1
            !atoms%epot = atoms%epot+t1/3.d0
        endif
        ru=atoms%cellvec(3,3)-atoms%ratp(3,iat)
        if (ru<5.5d0) then
            t2=200.d0*exp(-1.5d0*ru)
            atoms%fat(3,iat)=atoms%fat(3,iat)-t2
            !atoms%epot = atoms%epot+t2/3.d0
        endif
    enddo
end subroutine plane_repulsion
!*************************************************************************************************************
subroutine thermostat_evolution(atoms,zeta_next,zeta,zeta_prev,dzeta,mass_q,kt,ntherm,imd)
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
subroutine get_atomic_mass(atoms,totmass)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms):: atoms
    real(8):: totmass,mass_conv = 1822.888484264545
    integer:: iat
    totmass=0.d0
    do iat=1,atoms%nat
        if     (trim(atoms%sat(iat))=='H') then
            atoms%amass(iat)= mass_conv* 1.0079
        else if(trim(atoms%sat(iat))=='He') then   
            atoms%amass(iat)= mass_conv* 4.0026 
        else if(trim(atoms%sat(iat))=='Li') then   
            atoms%amass(iat)= mass_conv* 6.941
        else if(trim(atoms%sat(iat))=='Be') then   
            atoms%amass(iat)= mass_conv* 9.0122
        else if(trim(atoms%sat(iat))=='B') then
            atoms%amass(iat)= mass_conv* 10.811 
        else if(trim(atoms%sat(iat))=='C') then
            atoms%amass(iat)= mass_conv* 12.0107 
        else if(trim(atoms%sat(iat))=='N') then
            atoms%amass(iat)= mass_conv* 14.0067 
        else if(trim(atoms%sat(iat))=='O') then
            atoms%amass(iat)= mass_conv* 15.9994 
        else if(trim(atoms%sat(iat))=='F') then
            atoms%amass(iat)= mass_conv* 18.9984 
        else if(trim(atoms%sat(iat))=='Ne') then   
            atoms%amass(iat)= mass_conv* 20.1797     
        else if(trim(atoms%sat(iat))=='Na') then   
            atoms%amass(iat)= mass_conv* 22.9897     
        else if(trim(atoms%sat(iat))=='Mg') then   
            atoms%amass(iat)= mass_conv* 24.305
        else if(trim(atoms%sat(iat))=='Al') then   
            atoms%amass(iat)= mass_conv* 26.9815     
        else if(trim(atoms%sat(iat))=='Si') then   
            atoms%amass(iat)= mass_conv* 28.0855     
        else if(trim(atoms%sat(iat))=='P') then
            atoms%amass(iat)= mass_conv* 30.9738 
        else if(trim(atoms%sat(iat))=='S') then
            atoms%amass(iat)= mass_conv* 32.065
        else if(trim(atoms%sat(iat))=='Cl') then   
            atoms%amass(iat)= mass_conv* 35.453 
        else if(trim(atoms%sat(iat))=='K') then
            atoms%amass(iat)= mass_conv* 39.0983 
        else if(trim(atoms%sat(iat))=='Ar') then   
            atoms%amass(iat)= mass_conv* 39.948
        else if(trim(atoms%sat(iat))=='Ca') then   
            atoms%amass(iat)= mass_conv* 40.078
        else if(trim(atoms%sat(iat))=='Sc') then   
            atoms%amass(iat)= mass_conv* 44.9559     
        else if(trim(atoms%sat(iat))=='Ti') then   
            atoms%amass(iat)= mass_conv* 47.867 
        else if(trim(atoms%sat(iat))=='V') then
            atoms%amass(iat)= mass_conv* 50.9415 
        else if(trim(atoms%sat(iat))=='Cr') then   
            atoms%amass(iat)= mass_conv* 51.9961     
        else if(trim(atoms%sat(iat))=='Mn') then   
            atoms%amass(iat)= mass_conv* 54.938
        else if(trim(atoms%sat(iat))=='Fe') then   
            atoms%amass(iat)= mass_conv* 55.845 
        else if(trim(atoms%sat(iat))=='Ni') then   
            atoms%amass(iat)= mass_conv* 58.6934     
        else if(trim(atoms%sat(iat))=='Co') then   
            atoms%amass(iat)= mass_conv* 58.9332     
        else if(trim(atoms%sat(iat))=='Cu') then   
            atoms%amass(iat)= mass_conv* 63.546
        else if(trim(atoms%sat(iat))=='Zn') then   
            atoms%amass(iat)= mass_conv* 65.39 
        else if(trim(atoms%sat(iat))=='Ga') then   
            atoms%amass(iat)= mass_conv* 69.723 
        else if(trim(atoms%sat(iat))=='Ge') then   
            atoms%amass(iat)= mass_conv* 72.64 
        else if(trim(atoms%sat(iat))=='As') then   
            atoms%amass(iat)= mass_conv* 74.9216     
        else if(trim(atoms%sat(iat))=='Se') then   
            atoms%amass(iat)= mass_conv* 78.96 
        else if(trim(atoms%sat(iat))=='Br') then   
            atoms%amass(iat)= mass_conv* 79.904 
        else if(trim(atoms%sat(iat))=='Kr') then   
            atoms%amass(iat)= mass_conv* 83.8 
        else if(trim(atoms%sat(iat))=='Rb') then   
            atoms%amass(iat)= mass_conv* 85.4678     
        else if(trim(atoms%sat(iat))=='Sr') then   
            atoms%amass(iat)= mass_conv* 87.62 
        else if(trim(atoms%sat(iat))=='Y') then
            atoms%amass(iat)= mass_conv* 88.9059 
        else if(trim(atoms%sat(iat))=='Zr') then   
            atoms%amass(iat)= mass_conv* 91.224
        else if(trim(atoms%sat(iat))=='Nb') then   
            atoms%amass(iat)= mass_conv* 92.9064     
        else if(trim(atoms%sat(iat))=='Mo') then   
            atoms%amass(iat)= mass_conv* 95.94 
        else if(trim(atoms%sat(iat))=='Tc') then   
            atoms%amass(iat)= mass_conv* 98 
        else if(trim(atoms%sat(iat))=='Ru') then   
            atoms%amass(iat)= mass_conv* 101.07
        else if(trim(atoms%sat(iat))=='Rh') then   
            atoms%amass(iat)= mass_conv* 102.9055    
        else if(trim(atoms%sat(iat))=='Pd') then   
            atoms%amass(iat)= mass_conv* 106.42
        else if(trim(atoms%sat(iat))=='Ag') then   
            atoms%amass(iat)= mass_conv* 107.8682    
        else if(trim(atoms%sat(iat))=='Cd') then   
            atoms%amass(iat)= mass_conv* 112.411     
        else if(trim(atoms%sat(iat))=='In') then   
            atoms%amass(iat)= mass_conv* 114.818     
        else if(trim(atoms%sat(iat))=='Sn') then   
            atoms%amass(iat)= mass_conv* 118.71
        else if(trim(atoms%sat(iat))=='Sb') then   
            atoms%amass(iat)= mass_conv* 121.76 
        else if(trim(atoms%sat(iat))=='I') then
            atoms%amass(iat)= mass_conv* 126.9045
        else if(trim(atoms%sat(iat))=='Te') then   
            atoms%amass(iat)= mass_conv* 127.6 
        else if(trim(atoms%sat(iat))=='Xe') then   
            atoms%amass(iat)= mass_conv* 131.293     
        else if(trim(atoms%sat(iat))=='Cs') then   
            atoms%amass(iat)= mass_conv* 132.9055    
        else if(trim(atoms%sat(iat))=='Ba') then   
            atoms%amass(iat)= mass_conv* 137.327     
        else if(trim(atoms%sat(iat))=='La') then   
            atoms%amass(iat)= mass_conv* 138.9055    
        else if(trim(atoms%sat(iat))=='Ce') then   
            atoms%amass(iat)= mass_conv* 140.116     
        else if(trim(atoms%sat(iat))=='Pr') then   
            atoms%amass(iat)= mass_conv* 140.9077    
        else if(trim(atoms%sat(iat))=='Nd') then   
            atoms%amass(iat)= mass_conv* 144.24 
        else if(trim(atoms%sat(iat))=='Pm') then   
            atoms%amass(iat)= mass_conv* 145
        else if(trim(atoms%sat(iat))=='Sm') then   
            atoms%amass(iat)= mass_conv* 150.36
        else if(trim(atoms%sat(iat))=='Eu') then   
            atoms%amass(iat)= mass_conv* 151.964     
        else if(trim(atoms%sat(iat))=='Gd') then   
            atoms%amass(iat)= mass_conv* 157.25
        else if(trim(atoms%sat(iat))=='Tb') then   
            atoms%amass(iat)= mass_conv* 158.9253    
        else if(trim(atoms%sat(iat))=='Dy') then   
            atoms%amass(iat)= mass_conv* 162.5
        else if(trim(atoms%sat(iat))=='Ho') then   
            atoms%amass(iat)= mass_conv* 164.9303    
        else if(trim(atoms%sat(iat))=='Er') then   
            atoms%amass(iat)= mass_conv* 167.259     
        else if(trim(atoms%sat(iat))=='Tm') then   
            atoms%amass(iat)= mass_conv* 168.9342    
        else if(trim(atoms%sat(iat))=='Yb') then   
            atoms%amass(iat)= mass_conv* 173.04
        else if(trim(atoms%sat(iat))=='Lu') then   
            atoms%amass(iat)= mass_conv* 174.967     
        else if(trim(atoms%sat(iat))=='Hf') then   
            atoms%amass(iat)= mass_conv* 178.49
        else if(trim(atoms%sat(iat))=='Ta') then   
            atoms%amass(iat)= mass_conv* 180.9479    
        else if(trim(atoms%sat(iat))=='W') then
            atoms%amass(iat)= mass_conv* 183.84
        else if(trim(atoms%sat(iat))=='Re') then   
            atoms%amass(iat)= mass_conv* 186.207     
        else if(trim(atoms%sat(iat))=='Os') then   
            atoms%amass(iat)= mass_conv* 190.23
        else if(trim(atoms%sat(iat))=='Ir') then   
            atoms%amass(iat)= mass_conv* 192.217     
        else if(trim(atoms%sat(iat))=='Pt') then   
            atoms%amass(iat)= mass_conv* 195.078     
        else if(trim(atoms%sat(iat))=='Au') then   
            atoms%amass(iat)= mass_conv* 196.9665    
        else if(trim(atoms%sat(iat))=='Hg') then   
            atoms%amass(iat)= mass_conv* 200.59
        else if(trim(atoms%sat(iat))=='Tl') then   
            atoms%amass(iat)= mass_conv* 204.3833    
        else if(trim(atoms%sat(iat))=='Pb') then   
            atoms%amass(iat)= mass_conv* 207.2 
        else if(trim(atoms%sat(iat))=='Bi') then 
            atoms%amass(iat)= mass_conv* 208.9804    
        else
            write(*,*)"unknown atomic type"
            stop
        endif
        totmass=totmass+atoms%amass(iat)
    enddo
end subroutine get_atomic_mass
!*****************************************************************************************
subroutine write_trajectory_velocity(parini,atoms,file_info,rat_init,imd,ntherm,zeta,dzeta)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_file_info, update_ratp
    use mod_acf, only: acf_write
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variables
    type(typ_atoms):: atoms
    type(typ_file_info):: file_info
    integer:: iat
    integer:: imd, ntherm, ith
    real(8):: zeta(ntherm), dzeta(ntherm)
    real(8):: rat_init(3,atoms%nat)
    real(8):: r, dx(3) , rsq, msd1, msd2, msd3
    integer::  vfile
    real(8):: sumf1, sumf2, sumf3
    logical:: lfist= .true.
    lfist= .true.

    if(mod(imd-1,100)==0) then
        write(1111,*) '#'
        write(1111,*) '#    imd = ',imd, parini%time_dynamics
        write(1111,*) '#'
    endif
    msd1= 0.d0
    msd2= 0.d0
    msd3= 0.d0
    call update_ratp(atoms)
    do iat=1,atoms%nat

        dx(1:3)=atoms%ratp(:,iat)-rat_init(:,iat)
        rsq=(dx(1)**2+dx(2)**2+dx(3)**2)
        r=sqrt(dx(1)**2+dx(2)**2+dx(3)**2)
        msd1 = msd1 + rsq                !all directions
        msd2 = msd2 + dx(1)**2+dx(2)**2  !x,y directions
        msd3 = msd3 + dx(3)**2           !z   direction
        if(mod(imd,1000)==0) then
            write(1111,'(i5,2a5,4es25.17)')iat," ",atoms%sat(iat), dx, r
        endif
    enddo
    write(1112,'(i15,3es25.14)') imd-1 , msd1/atoms%nat , msd2/atoms%nat, msd3/atoms%nat 
    if (lfist) then
        open(unit=1003,file="velocity0",status='replace')
        vfile=1003
        lfist = .false.
    else
        open(unit=1002,file="velocity1",status='replace')
        vfile = 1002
        lfist = .true.
    endif

    file_info%file_position='append'
    call acf_write(file_info,atoms=atoms,strkey='trajectory')

    if(mod(imd-1,500)==0) then
        write(vfile,*) '#'
        write(vfile,*) '#    imd = ', imd
        write(vfile,*) '#'
        do iat=1,atoms%nat
            write(vfile,'(3es25.17)') atoms%vat(1,iat),atoms%vat(2,iat),atoms%vat(3,iat)
        enddo
        do ith=1,ntherm
            write(vfile,'(2es25.17)') zeta(ith),dzeta(ith)
        enddo
    endif
    if (lfist) then
        close(1002)
    else
        close(1003)
    endif
end subroutine write_trajectory_velocity
!*****************************************************************************************
