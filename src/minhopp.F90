!*****************************************************************************************
subroutine minimahopping(parini)
    use mod_parini, only: typ_parini
    use mod_task, only: time_exceeded
    use mod_minhopp, only: nstep, nlmin, nlminx, ekin, istep, ihopp, kerathopp, ediff, etoler, re_sm, &
        nlmin_old, minter, eref, nbuf, earr, dt, count_md, count_opt, escaped, accepted
    use mod_processors, only: parallel, nproc, iproc, imaster, mpi_comm_abz
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_copy
    use mod_opt, only: typ_paropt
    use yaml_output
    !minima hopping program with restart option.
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    integer:: ierr,iconf
    type(typ_paropt):: paropt_prec, paropt
    type(typ_atoms):: atoms_curr, atoms_hopp
    type(typ_atoms_arr):: atoms_allproc
    type(typ_atoms_arr):: atoms_locmin
#if defined(MPI)
    include 'mpif.h'
#endif
    call init_minimahopping(parini,atoms_curr,atoms_hopp,atoms_allproc,atoms_locmin,paropt,paropt_prec)
    !if(iproc==0) call sleep(20)
    call yaml_sequence_open('MH steps')
    call yaml_sequence(advance='no')
    call relax_minhopp(parini,atoms_curr,paropt_prec,paropt)
    !stop
    call yaml_mapping_open('input structure relaxed',flow=.true.)
    call yaml_map('iproc',iproc)
    call yaml_map('epot',atoms_curr%epot,fmt='(es24.15)')
    call yaml_mapping_close()
    !write(*,'(a,i4,e24.15)') 'input(relaxed): iproc,atoms_curr%epot',iproc,atoms_curr%epot
    if(nlmin>0) then
        write(*,'(a,1x,i4,2e24.15)') 'new/old energy for input configuration: iproc,new,old', &
            !iproc,atoms_curr%epot,atoms_allproc%epotall(iproc+1)
            iproc,atoms_curr%epot,atoms_allproc%atoms(iproc+1)%epot
        !if(abs(atoms_curr%epot-atoms_allproc%epotall(iproc+1))>etoler) write(*,'(a)') 'WARNING: new/old energies for input file differ'
        if(abs(atoms_curr%epot-atoms_allproc%atoms(iproc+1)%epot)>etoler) write(*,'(a)') 'WARNING: new/old energies for input file differ'
    endif
    if(nlmin==0) call minhopp_newrun_initialization(atoms_curr,atoms_locmin)
    re_sm=min(atoms_curr%epot,earr(1))
    call yaml_map('re_sm',re_sm,fmt='(es20.12)')
    !call yaml_mapping_open('?????')
    !call yaml_map('iproc',iproc,fmt='(i8)')
    !call yaml_map('istep',istep,fmt='(i8)')
    !call yaml_map('ihopp',ihopp,fmt='(i8)')
    !call yaml_map('epot',atoms_curr%epot,fmt='(es20.12)')
    !call yaml_map('ediff',,fmt='()')
    !call yaml_map('',,fmt='()')
    !call yaml_mapping_close()
    !write(*,'(a,i4,e24.15)') 'iproc,initial re_sm',iproc,re_sm
    write(2,'(a1,i3.3,5x,a,3x,a,10x,a,12x,a,2(6x,a))') '#',iproc,'istep','ihopp','epot','ediff','ekin','dt' !,'',''
    write(2,'(i4.3,i10,i6,e24.15,3f10.4,15x,a,2f9.1)') iproc,istep,ihopp,atoms_curr%epot,ediff,ekin,dt,'  --',count_md,count_opt
    !rathopp(1:3,1:nat)=atoms_curr%rat(1:3,1:nat)
    !atoms_hopp%rat(1:3,1:atoms_curr%nat)=atoms_curr%rat(1:3,1:atoms_curr%nat)
    call atom_copy(atoms_curr,atoms_hopp,'atoms_curr->atoms_hopp')
    if(parallel) call send_minimum_to_all(atoms_curr)
    nlmin_old=nlmin
    call request_receive(atoms_allproc)
    escaped=.true.
    accepted=.true.
    call report_minhopp_iteration_info(atoms_curr)
    do istep=1,nstep+1 !outer loop for hopping
        call yaml_sequence(advance='no')
        !collecting local minima data if send check by other processes.
        call collect_data_from_all_processors(5,atoms_curr,atoms_allproc,atoms_locmin)
        !write(*,'(a,i4,3i7,e24.15)') 'iproc,istep,ihopp,nlmin,erat ',iproc,istep-1,ihopp,nlmin,atoms_curr%epot
        !Energy has reached taregt eref and global minimum is presumably found
        if(re_sm<=eref) then
            call yaml_map('success, eref is achieved by iproc',iproc)
            !write(*,'(a,i4)') 'success: eref is achieved. iproc= ',iproc
            exit
        endif
        if(escaped) then
            !if(nlmin>=nlminx) exit
            !if(nlmin-nlmin_old>=minter) then !this prevents writing intermediate
            !even if minter=1, for the moment let's send in every cycle but I need
            !to find a better way control this action.
            if(1>0) then 
                call yaml_mapping_open('escaped',flow=.true.)
                call yaml_map('iproc',iproc)
                call yaml_map('nlmin',nlmin)
                call yaml_map('nlmin_old',nlmin_old)
                call yaml_map('minter',minter)
                call yaml_mapping_close()
                !write(*,'(a,i4,2i7,i6)') 'iproc,nlmin,nlmin_old,minter',iproc,nlmin,nlmin_old,minter
                call send_minhopp_parameters_to_all(atoms_curr)
            endif
        endif
        if(nlmin-nlmin_old>=minter) call write_minhopp(atoms_allproc,atoms_locmin)
        if(istep>nstep) exit
        !check whether time exceeded
        call check_whether_time_exceeded
        if(time_exceeded) exit
        !-------------------------------------------------------------
        !atoms_hopp%rat(1:3,1:atoms_curr%nat)=atoms_curr%rat(1:3,1:atoms_curr%nat)
        call atom_copy(atoms_curr,atoms_hopp,'atoms_curr->atoms_hopp')
        call mdescape(parini,atoms_hopp)
        call relax_minhopp(parini,atoms_hopp,paropt_prec,paropt)
        if(abs(atoms_curr%epot-atoms_hopp%epot)<etoler) then !failed to escape.
            call escape_failed(parini,atoms_curr%epot,atoms_hopp%epot)
            escaped=.false.
        else
            escaped=.true.
        endif
        !-------------------------------------------------------------
        if(escaped) then
            !check whether new minimum
            call hunt2(min(nlmin,nlminx+nbuf),earr(1),atoms_hopp%epot,kerathopp)
            if(abs(atoms_hopp%epot-earr(kerathopp))<etoler) then 
                call already_visited_minimum(parini)
            else
                call new_minimum(atoms_hopp)
            endif
            ihopp=ihopp+1
            !Monte Carlo step for local minima hopping
            if(atoms_hopp%epot-atoms_curr%epot<ediff) then !local minimum accepted.
                accepted=.true.
                call local_minimum_accepted(atoms_hopp,atoms_curr,atoms_locmin)
            else !local minima rejected.
                accepted=.false.
                call local_minimum_rejected(atoms_hopp)
            endif
        endif
        call report_minhopp_iteration_info(atoms_curr)
    enddo !end of outer loop for hopping
    call yaml_sequence_close()
    call send_minhopp_parameters_to_all(atoms_curr)
#if defined(MPI)
    if(parallel) call MPI_BARRIER(mpi_comm_abz,ierr)
#endif
    call yaml_map('collecting accepted minima sent by other procs, final, iproc',iproc)
    !write(*,'(a,i4)') 'final collect of accepted minima sent by all other processors: iproc',iproc
    call collect_data_from_all_processors(50,atoms_curr,atoms_allproc,atoms_locmin)
    call cancel_excessive_irecv
    !if(parallel) call MPI_BARRIER(mpi_comm_abz,ierr)
    call final_minimahopping(parini,atoms_curr,atoms_hopp,atoms_allproc,atoms_locmin,paropt,paropt_prec)
end subroutine minimahopping
!*****************************************************************************************
subroutine init_minimahopping(parini,atoms_curr,atoms_hopp,atoms_allproc,atoms_locmin,paropt,paropt_prec)
    use mod_parini, only: typ_parini
    use mod_minhopp, only: nlmin, alpha_soften, npminx, nstep, alpha1, alpha2, &
        beta1, beta2, beta3, eref, etoler, minter, mdmin, nsoften, nrandoff
    use mod_processors, only: nproc, iproc, imaster, mpi_comm_abz, parallel
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_copy, set_ndof
    use mod_yaml_conf, only: read_yaml_conf
    use mod_opt, only: typ_paropt
    use mod_potential, only: potential
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_paropt), intent(inout):: paropt_prec, paropt
    type(typ_atoms), intent(inout):: atoms_curr, atoms_hopp
    type(typ_atoms_arr), intent(inout):: atoms_allproc
    type(typ_atoms_arr), intent(inout):: atoms_locmin
    !local variables
    integer:: ios, ierr, iconf, iat, jproc
    character(50):: filename
    character(3):: fn
    nstep=parini%nstep_minhopp
    nsoften=parini%nsoften_minhopp
    mdmin=parini%mdmin_minhopp
    minter=parini%minter_minhopp
    nrandoff=parini%nrandoff_minhopp
    npminx=parini%npminx_minhopp
    etoler=parini%etoler_minhopp
    eref=parini%eref_minhopp
    alpha1=parini%alpha1_minhopp
    alpha2=parini%alpha2_minhopp
    beta1=parini%beta1_minhopp
    beta2=parini%beta2_minhopp
    beta3=parini%beta3_minhopp
    nstep=max(1,nstep/nproc) !nstep in input.ini is considered to be for all processors.
    potential=trim(parini%potential_potential)
    paropt=parini%paropt_geopt
    call initminimize(paropt)
    if(parini%two_level_geopt) then
        paropt_prec=parini%paropt_geopt_prec
        call initminimize(paropt_prec)
    endif
    call print_minhopp_parameters
    !reading nlmin,nlminx,ediffarr,ekinarr,dtarr and !distribute among all proc.
    call read_minhopp_parameters(parini)
    if(iproc==imaster) then
        call yaml_comment('reading cellvec sizes and atomic positions from posinp.acf')
        !write(*,'(a)') 'reading cellvec sizes and atomic positions from posinp.acf'
        !call acf_read(parini,'posinp.acf',nproc,atoms_all=atoms_allproc)
        !call acf_read_new(parini,'posinp.acf',nproc,atoms_arr=atoms_allproc)
        !inquire(file='posinp.yaml',exist=yaml_exists)
        !inquire(file='posinp.acf',exist=acf_exists)
        !if(yaml_exists) then
            call read_yaml_conf(parini,'posinp.yaml',nproc,atoms_arr=atoms_allproc)
        !elseif(acf_exists) then
        !    call acf_read_new(parini,'posinp.acf',10000,atoms_arr=atoms_allproc)
        !else
        !endif
        !write(*,*) atoms_allproc%atoms%nat,atoms_curr%nat
        !write(*,*) atoms_allproc%atoms%bemoved(:,:)
        !write(*,*) atoms_curr%bemoved(:,:)
        !call set_ndof(atoms_allproc%atoms)
        call set_ndof(atoms_allproc%atoms(1))
        !write(*,'(a,i6)') 'ndof= ',atoms_allproc%atoms%ndof
        call yaml_map('ndof',atoms_allproc%atoms(1)%ndof)
        !write(*,'(a,i6)') 'ndof= ',atoms_allproc%atoms(1)%ndof
        !call atom_copy(atoms_allproc%atoms,atoms_curr,'atoms_allproc%atoms->atoms_curr')
        call atom_copy(atoms_allproc%atoms(1),atoms_curr,'atoms_allproc%atoms->atoms_curr')
        !do jproc=1,nproc !HERE: what is loop doing?
        !    call atom_copy(atoms_allproc%atoms(jproc),atoms_curr,'atoms_allproc%atoms->atoms_curr')
        !    if(.not. allocated(atoms_curr%rat)) then
        !        allocate(atoms_curr%rat(3,atoms_allproc%atoms(jproc)%nat))
        !    endif
        !enddo
    endif
    call readnat(atoms_curr) !read number of atoms.
    call read_poscur_alborz(atoms_curr,atoms_allproc) !Read initial positions.
    call allocate_minhopp_arrays2(atoms_curr%nat,nproc)
    if(nlmin>0) call read_earr
    !call allocateatomsarrays(nproc)
    call setpot_init(parini,atoms_curr,paropt,paropt_prec)
    !-----------------------------------------------------------------
    !atoms_hopp%ndof=atoms_curr%ndof
    !call atom_allocate(atoms_hopp,atoms_curr%nat,0,0,sat=.true.,vat=.true.,amass=.true.,fat=.true.,bemoved=.true.)
    call atom_copy(atoms_curr,atoms_hopp,'atoms_curr->atoms_hopp')
    !atoms_hopp%sat(1:atoms_hopp%nat)=atoms_curr%sat(1:atoms_hopp%nat)
    !atoms_hopp%bemoved(1:3,1:atoms_hopp%nat)=atoms_curr%bemoved(1:3,1:atoms_hopp%nat)
    !atoms_hopp%boundcond=atoms_curr%boundcond
    !atoms_hopp%cellvec(1:3,1:3)=atoms_curr%cellvec(1:3,1:3)
    !-----------------------------------------------------------------
    call set_amass(atoms_curr)
    call set_amass(atoms_hopp)
    !atoms_hopp%amass(1:atoms_hopp%nat)=atoms_hopp%amass(1:atoms_hopp%nat)
    !-----------------------------------------------------------------
    !if restart, read low energy configuration previously found.
    atoms_locmin%nconfmax=npminx
    if(nlmin==0 .or. iproc/=imaster) then
        !call atom_copy(atoms_curr,atoms_locmin%atoms,'atoms_curr->atoms_locmin%atoms')
        allocate(atoms_locmin%atoms(atoms_locmin%nconfmax))
        !call atom_copy(atoms_curr,atoms_locmin%atoms(1),'atoms_curr->atoms_locmin%atoms')
    endif
    if(nlmin>0) then
        if(iproc/=imaster) then
        do iconf=1,atoms_locmin%nconfmax
            call atom_copy(atoms_curr,atoms_locmin%atoms(iconf),'atoms_curr->atoms_locmin%atoms')
        enddo
        endif
        if(iproc==imaster) then
            !call acf_read(parini,'poslow.acf',npminx,atoms_arr=atoms_locmin)
            !call acf_read_new(parini,'poslow.acf',npminx,atoms_arr=atoms_locmin)
            call read_yaml_conf(parini,'poslow.yaml',npminx,atoms_arr=atoms_locmin)
        endif
        call read_poslow(atoms_locmin)
    endif
    !call execute_command_line("mkdir -p monminhopp")
    call system("mkdir -p monminhopp")
    write(fn,'(i3.3)') iproc;filename='monminhopp/monitoring.'//fn
    open(unit=2,file=filename,status='replace',iostat=ios)
    if(ios/=0) then;write(*,'(a,1x,a)') 'ERROR: failure openning',filename;stop;endif
    !nsoften=7
    !call readparopt(paropt_prec,paropt)
    !write(*,*) 'HERE ',trim(paropt%approach)
    alpha_soften=paropt%alphax*0.5d0 !*0.8d0
end subroutine init_minimahopping
!*****************************************************************************************
subroutine final_minimahopping(parini,atoms_curr,atoms_hopp,atoms_allproc,atoms_locmin,paropt,paropt_prec)
    use mod_parini, only: typ_parini
    use mod_minhopp, only: count_md_tot, count_opt_tot, count_soften_tot, &
        fcall_tot_all, fcall_tot_all_md, fcall_tot_all_opt, fcall_tot_all_soften
    use mod_processors, only: iproc, imaster, mpi_comm_abz, parallel
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_deallocate
    use mod_opt, only: typ_paropt
    use mod_potential, only: fcalls
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_paropt), intent(inout):: paropt_prec, paropt
    type(typ_atoms), intent(inout):: atoms_curr, atoms_hopp
    type(typ_atoms_arr), intent(inout):: atoms_allproc
    type(typ_atoms_arr), intent(inout):: atoms_locmin
    !local variables
    integer:: ierr, iconf
#if defined(MPI)
    include 'mpif.h'
    call MPI_REDUCE(fcalls        ,fcall_tot_all       ,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_abz,ierr)
    call MPI_REDUCE(count_md_tot    ,fcall_tot_all_md    ,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_abz,ierr)
    call MPI_REDUCE(count_opt_tot   ,fcall_tot_all_opt   ,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_abz,ierr)
    call MPI_REDUCE(count_soften_tot,fcall_tot_all_soften,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_abz,ierr)
    if(parallel) call MPI_BARRIER(mpi_comm_abz,ierr)
#else
    fcall_tot_all       =fcalls          
    fcall_tot_all_md    =count_md_tot    
    fcall_tot_all_opt   =count_opt_tot   
    fcall_tot_all_soften=count_soften_tot
#endif
    call write_minhopp(atoms_allproc,atoms_locmin)
    call print_final_statistics
    call deallocate_minhopp_arrays
    call atom_deallocate(atoms_curr)
    call atom_deallocate(atoms_hopp)
    !call atom_deallocate_old(atoms_locmin%atoms)
    !call atom_deallocate_old(atoms_allproc%atoms)
    !call atom_all_deallocate(atoms_locmin,ratall=.true.,epotall=.true.)
    !call atom_all_deallocate(atoms_allproc,ratall=.true.,epotall=.true.)
    do iconf=1,atoms_locmin%nconfmax
        call atom_deallocate(atoms_locmin%atoms(iconf))
    enddo
    deallocate(atoms_locmin%atoms)
    do iconf=1,atoms_allproc%nconfmax
        call atom_deallocate(atoms_allproc%atoms(iconf))
    enddo
    deallocate(atoms_allproc%atoms)
    close(2)
    if(parini%two_level_geopt) then
        call finalminimize(paropt_prec)
    endif
    call finalminimize(paropt)
    call setpot_final(parini,atoms_curr)
end subroutine final_minimahopping
!*****************************************************************************************
subroutine set_amass(atoms_hopp)
    use mod_atoms, only: typ_atoms
    use mod_processors, only: iproc, imaster
    implicit none
    type(typ_atoms), intent(inout):: atoms_hopp
    !local variables
    integer:: ios, iat
    logical:: amass_exists
    if(.not. allocated(atoms_hopp%amass)) stop 'ERROR: atoms_hopp%amass not allocated in set_amass'
    inquire(file='amass.dat',exist=amass_exists)
    if(amass_exists) then
        if(iproc==imaster) then
            write(*,'(a)') 'amass.dat exists, minhopp will try to use its contents.'
        endif
        open(unit=21,file='amass.dat',status='old',iostat=ios)
        if(ios/=0) then;write(*,'(a)') 'ERROR: failure openning amass.dat';stop;endif
        do iat=1,atoms_hopp%nat
            read(21,*) atoms_hopp%amass(iat)
        enddo
        close(21)
    else
        atoms_hopp%amass(1:atoms_hopp%nat)=1.d0 !16608.045513137928d0
    endif
end subroutine set_amass
!*****************************************************************************************
subroutine relax_minhopp(parini,atoms,paropt_prec,paropt)
    use mod_parini, only: typ_parini
    use mod_minhopp, only: die, count_opt, count_opt_tot, istep
    use mod_opt, only: typ_paropt
    use mod_atoms, only: typ_atoms, get_rat, set_rat
    use mod_potential, only: fcalls
    use mod_processors, only: iproc
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_paropt), intent(inout):: paropt_prec, paropt
    !local variables
    integer:: istat
    real(8):: count_tt
    real(8), allocatable:: rat_tmp(:,:) !a temporary array
    type(typ_paropt):: paropt_backup
    character(32):: ttfn
    count_tt=fcalls
    allocate(rat_tmp(3,atoms%nat),stat=istat)
    if(istat/=0) write(*,*) 'ERROR: failure allocating array rat_tmp'
    if(parini%two_level_geopt) then
        if(parini%trajectory_minhopp) then
            paropt_prec%trajectory=parini%trajectory_minhopp
            paropt_prec%print_force=parini%print_force_minhopp
            write(ttfn,'(a9,i3.3,a5,i7.7,a8)') 'traj_proc',iproc,'_step',istep,'_pgo.bin'
            paropt_prec%filename=ttfn
        endif
        call setpot_geopt_prec
        call minimize(parini,iproc,atoms,paropt_prec)
    endif
    call setpot_geopt
    call get_rat(atoms,rat_tmp)
    !paropt%nitfire=10
    if(parini%trajectory_minhopp) then
        paropt%trajectory=parini%trajectory_minhopp
        paropt%print_force=parini%print_force_minhopp
        write(ttfn,'(a9,i3.3,a5,i7.7,a8)') 'traj_proc',iproc,'_step',istep,'_mgo.bin'
        paropt%filename=ttfn
    endif
    call minimize(parini,iproc,atoms,paropt)
    if(.not. paropt%converged) then
        call set_rat(atoms,rat_tmp,setall=.true.)
        paropt_backup=paropt
        if(trim(paropt%approach)=='BFGS'  ) paropt_backup%approach='FIRE'
        if(trim(paropt%approach)=='SDCG'  ) paropt_backup%approach='FIRE'
        if(trim(paropt%approach)=='SDDIIS') paropt_backup%approach='FIRE'
        if(trim(paropt%approach)=='NLBFGS') paropt_backup%approach='FIRE'
        if(trim(paropt%approach)=='FIRE'  ) paropt_backup%approach='SD'
        call get_rat(atoms,rat_tmp)
        if(parini%trajectory_minhopp) then
            paropt_backup%trajectory=parini%trajectory_minhopp
            write(ttfn,'(a9,i3.3,a5,i7.7,a8)') 'traj_proc',iproc,'_step',istep,'_bgo.acf'
            paropt_backup%filename=ttfn
        endif
        call minimize(parini,iproc,atoms,paropt_backup)
        if(.not. paropt_backup%converged) then
            call set_rat(atoms,rat_tmp,setall=.true.)
            die=.true.
        endif
    endif
    count_opt=fcalls-count_tt
    count_opt_tot=count_opt_tot+count_opt
    deallocate(rat_tmp,stat=istat)
    if(istat/=0) write(*,*) 'ERROR: failure deallocating array rat_tmp'
end subroutine relax_minhopp
!*****************************************************************************************
subroutine print_minhopp_parameters
    use mod_task, only: time_start
    use mod_minhopp, only: beta1, beta2, beta3, alpha1, alpha2, mdmin, minter, npminx
    !use mod_atoms, only:
    use mod_processors, only: iproc, imaster
    use yaml_output
    implicit none
    real(8)::ratio,tt,p1,p2
    ratio=-log(alpha2)/log(alpha1)
    tt=beta2**ratio*beta3
    if(iproc==imaster) then
        call yaml_mapping_open('minhopp parameters',flow=.true.)
        call yaml_map('beta1',beta1,fmt='(f20.5)')
        call yaml_map('beta2',beta2,fmt='(f20.5)')
        call yaml_map('beta3',beta3,fmt='(f20.5)')
        call yaml_map('alpha1',alpha1,fmt='(f20.5)')
        call yaml_map('alpha2',alpha2,fmt='(f20.5)')
        p1=ratio/(1.d0+ratio);p2=1.d0/(1.d0+ratio)
        call yaml_map('p1',p1,fmt='(es11.3)')
        call yaml_map('p2',p2,fmt='(es11.3)')
        call yaml_map('critical ratio (.ge. 1)',tt,fmt='(es12.4)')
        call yaml_map('mdmin',mdmin,fmt='(i8)')
        call yaml_map('npminx',npminx,fmt='(i8)')
        call yaml_map('minter',minter,fmt='(i8)')
        call yaml_mapping_close()
        !write(*,'(a,3f20.5)') 'beta1,beta2,beta3',beta1,beta2,beta3
        !write(*,'(a,2f20.5)') 'alpha1,alpha2',alpha1,alpha2
        !p1=ratio/(1.d0+ratio);p2=1.d0/(1.d0+ratio)
        !write(*,'(a,2e10.3)') 'predicted fraction accepted, rejected',p1,p2
        !write(*,'(a,e11.4)') 'critical ratio (.ge. 1)',tt
        !write(*,'(a,i4,2i6)') 'mdmin,npminx,minter',mdmin,npminx,minter
    endif
    !if(tt<0.999999999999999d0) stop 'incompatible alpha s, beta s'
end subroutine print_minhopp_parameters
!*****************************************************************************************
subroutine read_earr
    use mod_minhopp, only: nlmin, nlminx, nbuf, earr, nvisit
    use mod_processors, only: iproc, parallel, imaster, mpi_comm_abz
    implicit none
    !local variables
    integer:: ios, k, ierr
#if defined(MPI)
    include 'mpif.h'
#endif
    earr(0)=-1.d100
    if(iproc==imaster) then 
        open(unit=12,file='earr.dat',status='old',iostat=ios)
        if(ios/=0) stop 'ERROR: failure openning earr.dat'
        do k=1,nlmin
            read(12,*) earr(k),nvisit(k)
            if(earr(k)<earr(k-1)) stop 'wrong ordering in earr.dat'
        enddo
        write(*,'(a)') 'master processor read earr.dat'
        close(12)
    endif
#if defined(MPI)
    if(parallel) then
        call MPI_BCAST(earr,(nlminx+nbuf+1),MPI_DOUBLE_PRECISION,imaster,mpi_comm_abz,ierr)
        call MPI_BCAST(nvisit,(nlminx+nbuf+1),MPI_INTEGER,imaster,mpi_comm_abz,ierr)
    endif
#endif
end subroutine read_earr
!*****************************************************************************************
subroutine readnat(atoms_curr)
    use mod_minhopp, only: nlmin
    use mod_processors, only: iproc, parallel, imaster, mpi_comm_abz
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms):: atoms_curr
    !local variables
    integer:: ios, ierr
#if defined(MPI)
    include 'mpif.h'
#endif
    !if(nlmin==-1) stop 'ERROR: nlmin=-1, is not set.'
    if(iproc==imaster) then
    !    open(unit=9,file='posinp.xyz',status='old',iostat=ios)
    !    if(ios/=0) stop 'ERROR: failure openning posinp.xyz'
    !    read(9,*) atoms_curr%nat
    !    close(9)
    !    write(*,'(a,i7)') 'nat ',atoms_curr%nat
    endif
#if defined(MPI)
    if(parallel) then
        call MPI_BCAST(atoms_curr%nat,1,MPI_INTEGER,imaster,mpi_comm_abz,ierr)
    endif
#endif
end subroutine readnat
!*****************************************************************************************
subroutine read_poscur_alborz(atoms_curr,atoms_allproc)
    !use mod_minhopp, only: eratallproc
    use mod_processors, only: iproc, nproc, parallel, imaster, mpi_comm_abz
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_copy
    implicit none
    type(typ_atoms):: atoms_curr
    type(typ_atoms_arr):: atoms_allproc
    !local variables
    integer:: iat, ios, k, j, ierr, nconf, jproc, ii1, iconf
    character(20):: tatomnames
    !integer::status_mpi(MPI_STATUS_SIZE)
    !logical:: jobmustdie
    character(3):: str_motion
    real(8):: cv(3,3), x, y, z
#if defined(MPI)
    include 'mpif.h'
#endif
    !atoms_curr%nat=atoms_allproc%atoms%nat
    !write(*,*) atoms_curr%nat
    !call atom_allocate_new(atoms_curr,atoms_curr%nat,0,0) !this is commented because atoms_allproc%atoms->atoms_curr
    !call atom_copy(atoms_allproc%atoms,atoms_curr,'atoms_allproc%atoms->atoms_curr')
    !atoms_allproc%atoms%nat=atoms_curr%nat
    !atoms_allproc%nconfmax=nproc
    !atoms_allproc%nconf=atoms_allproc%nconfmax
    !call atom_all_allocate(atoms_allproc,ratall=.true.,epotall=.true.)
    !atoms_allproc%epotall(1:nproc)=1000.d0
    !jobmustdie=.false.
    if(iproc==imaster) then
        !open(unit=9,file='posinp.xyz',status='old',iostat=ios)
        !if(ios/=0) stop 'ERROR: failure openning posinp.xyz'
        !write(*,'(a)') 'reading cellvec sizes and atomic positions from posinp.xyz'
        !!read(55,*,iostat=k) cpulimit
        !nconf=1
        !do !nconf=0,nproc-1
        !    read(9,*,iostat=k) ii1,atoms_allproc%epotall(nconf)
        !    if(k<0) then;write(*,*) 'LINE-N: k=',k;exit;endif
        !    read(9,*,iostat=k) atoms_curr%boundcond,cv(1,1),cv(2,2),cv(3,3),cv(1,2),cv(1,3),cv(2,3)
        !    cv(2,1)=0.d0 ; cv(3,1)=0.d0 ; cv(3,2)=0.d0
        !    atoms_curr%cellvec(1:3,1:3)=cv(1:3,1:3)
        !    if(k<0) then;write(*,*) 'LINE-C: k=',k;exit;endif
        !    if(trim(atoms_curr%boundcond)/='free' .and. trim(atoms_curr%boundcond)/='wire' .and. &
        !       trim(atoms_curr%boundcond)/='slab' .and. trim(atoms_curr%boundcond)/='bulk') then
        !       stop 'ERROR: boundary condition in posinp.xyz is unknown.'
        !   endif
        !   if(nconf==1) then
                !cv(1:3,1:3)=atoms_allproc%atoms%cellvec(1:3,1:3)
                cv(1:3,1:3)=atoms_allproc%atoms(1)%cellvec(1:3,1:3)
                !write(*,'(a,3f20.10)') 'cellvec ',cv(1,1),cv(2,1),cv(3,1)
                !write(*,'(a,3f20.10)') 'cellvec ',cv(1,2),cv(2,2),cv(3,2)
                !write(*,'(a,3f20.10)') 'cellvec ',cv(1,3),cv(2,3),cv(3,3)
        !    endif
        !    atoms_curr%ndof=0
        !    do iat=1,atoms_curr%nat
        !        !read(9,*,iostat=k) tatomnames,(ratallproc(j,iat,nconf-1),j=1,3),str_motion
        !        read(9,*,iostat=k) tatomnames,x,y,z,str_motion
        !        atoms_allproc%ratall(1,iat,nconf)=x
        !        atoms_allproc%ratall(2,iat,nconf)=y
        !        atoms_allproc%ratall(3,iat,nconf)=z
        !        call string2bemoved(str_motion,atoms_curr%bemoved(1,iat))
        !        if(nconf==1) then
        !            atoms_curr%sat(iat)=trim(tatomnames)
        !            !write(*,*) iat,trim(tatomnames)
        !        else
        !            if(atoms_curr%sat(iat)/=trim(tatomnames)) jobmustdie=.true.
        !        endif
        !        !write(*,'(a,3e24.15,i6)') trim(tatomnames),(ratallproc(j,iat,nconf-1),j=1,3),k
        !        if(k<0) then;write(*,*) 'LINE-A: k=',k,iat;exit;endif
        !    enddo !end of loop over iat
        !    call set_ndof(atoms_curr)
        !    if(k<0) then;write(*,*) 'LINE-E: k=',k;exit;endif
        !    if(nconf==nproc) exit
        !    nconf=nconf+1
        !enddo !end of loop over nconf
        !if(k<0) nconf=nconf-1
        !write(*,'(a,i4)') 'number of configuration in posinp.xyz = ',nconf
        !close(9)
        !if(nconf<nproc) then
        !    write(*,'(a,2i4)') 'WARNING: number of configurations in posinp.xyz fewer than nproc.',nconf,nproc
        !endif
        nconf=atoms_allproc%nconf
        !write(*,*) 'NCONF= ',nconf
        do jproc=nconf-1,nproc-1
            !ratallproc(1:3,1:atoms_curr%nat,jproc)=ratallproc(1:3,1:atoms_curr%nat,mod(jproc,nconf))
            !!atoms_allproc%ratall(1:3,1:atoms_allproc%atoms%nat,jproc+1)= &
            !!    atoms_allproc%ratall(1:3,1:atoms_allproc%atoms%nat,mod(jproc,nconf)+1)
            !!atoms_allproc%epotall(jproc+1)=atoms_allproc%epotall(mod(jproc,nconf)+1)
            call atom_copy(atoms_allproc%atoms(mod(jproc,nconf)+1),atoms_allproc%atoms(jproc+1),'atoms_allproc%atoms->atoms_allproc%atoms')
        enddo
        !call set_ndof(atoms_allproc%atoms)
        !write(*,'(a,i6)') 'ndof= ',atoms_allproc%atoms%ndof
    endif !end of if for proc==imaster
    !if(jobmustdie) stop 'ERROR: something wrong with atoms symbol in posinp.xyz'
#if defined(MPI)
    if(parallel) then
        call MPI_BCAST(atoms_curr%cellvec,9,MPI_DOUBLE_PRECISION,imaster,mpi_comm_abz,ierr)
        call MPI_BCAST(atoms_curr%bemoved,3*atoms_curr%nat,MPI_LOGICAL,imaster,mpi_comm_abz,ierr)
        call MPI_BCAST(atoms_curr%ndof,1,MPI_INTEGER,imaster,mpi_comm_abz,ierr)
        call MPI_BCAST(atoms_curr%sat,5*atoms_curr%nat,MPI_CHARACTER,imaster,mpi_comm_abz,ierr)
        if(iproc/=imaster) atoms_curr%boundcond=''
        call MPI_BCAST(atoms_curr%boundcond,4,MPI_CHARACTER,imaster,mpi_comm_abz,ierr)
    endif
#endif
    if(trim(atoms_curr%boundcond)=='unkn' .or. trim(atoms_curr%boundcond)=='unknown') then
        write(*,'(a)') 'ERROR: boundary conditions is unknow in read_poscur in minhopp'
        stop
    endif
    atoms_allproc%nconf=nproc
    atoms_allproc%nconfmax=nproc
    if(iproc/=imaster) then
        !assigning nat and nproc and allocating arrays for iproc==imaster is
        !done while reading from file.
        !atoms_allproc%atoms%nat=atoms_curr%nat 
        !call atom_all_allocate(atoms_allproc,ratall=.true.,epotall=.true.)
        allocate(atoms_allproc%atoms(atoms_allproc%nconfmax))
        do iconf=1,atoms_allproc%nconfmax
            call atom_copy(atoms_curr,atoms_allproc%atoms(iconf),'atoms_curr->atoms_allproc%atoms')
        enddo
    endif
#if defined(MPI)
    if(parallel) then
        !call MPI_BCAST(atoms_allproc%ratall,3*atoms_curr%nat*nproc,MPI_DOUBLE_PRECISION,imaster,mpi_comm_abz,ierr)
        !call MPI_BCAST(atoms_allproc%epotall,nproc,MPI_DOUBLE_PRECISION,imaster,mpi_comm_abz,ierr)
        call MPI_atom_arr_copy(atoms_curr%nat,atoms_allproc)
    endif
#endif
    !atoms_curr%rat(1:3,1:atoms_curr%nat)=atoms_allproc%ratall(1:3,1:atoms_curr%nat,iproc+1)
    call atom_copy(atoms_allproc%atoms(iproc+1),atoms_curr,'allproc%atoms(iproc+1)=>atoms_curr')
end subroutine read_poscur_alborz
!*****************************************************************************************
subroutine read_minhopp_parameters(parini)
    use mod_parini, only: typ_parini
    use mod_processors, only: iproc, parallel, nproc, imaster, mpi_comm_abz
    use mod_atoms, only: atom_copy
    use mod_minhopp, only: nrandoff, ediff, ekin, dt, nlmin, nlminx, eref, etoler, ekinarr, &
        dtarr, ediffarr, nstep
    use mod_utils
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    integer:: ios, i, mproc, k, jproc, ierr
    real(8):: tts
#if defined(MPI)
    include 'mpif.h'
    !integer::status_mpi(MPI_STATUS_SIZE)
#endif
    if(iproc==imaster) then
        open(unit=11,file='input.minhopp',status='old',iostat=ios)
        if(ios/=0) stop 'ERROR: failure openning input.minhopp'
        call yaml_comment('master proc reads input parameters from input.minhopp and sends to all.')
        !write(*,'(a)') 'master proc reads input parameters from input.minhopp and sends to all.'
        read(11,*) nlmin
        !read(11,*) nlminx
        nlminx=nlmin+nstep+1 !+1 is need if nlmin=0 (and nbuf=0)
        !read(11,*) eref
        call yaml_mapping_open('input minima info',flow=.true.)
        call yaml_map('nlmin',nlmin)
        call yaml_map('nlminx',nlminx)
        call yaml_map('eref',eref)
        call yaml_mapping_close()
        !write(*,'(a,1x,2i8,e24.15)') 'nlmin,nlminx,eref ',nlmin,nlminx,eref
        if(nlmin==0) then
            call yaml_comment('This is a new run.')
            !write(*,'(a)') 'This is a new run.'
        else
            call yaml_comment('This is a restart run.')
            !write(*,'(a)') 'This is a restart run.'
        endif
    endif
    !eref=-0.100000000000000E+05
#if defined(MPI)
    if(parallel) then
        call MPI_BCAST(nlmin,1,MPI_INTEGER,imaster,mpi_comm_abz,ierr)
        call MPI_BCAST(nlminx,1,MPI_INTEGER,imaster,mpi_comm_abz,ierr)
    endif
#endif
    !etoler=accur
    if(nlmin>nlminx) stop 'nlmin>nlminx'
    call allocate_minhopp_arrays1(nproc)
    if(iproc==imaster) then
        do mproc=0,nproc-1
            read(11,*,iostat=k) ediffarr(mproc),ekinarr(mproc),dtarr(mproc)
            dtarr(mproc)=dtarr(mproc) !*41.34137333656d0
            if(k<0) then;write(*,*) 'k=',k;exit;endif
        enddo
        call yaml_map('number of sets for input parameters in input.minhopp',mproc)
        !write(*,'(a,i4,a)') 'There were ',mproc,' sets of input parameters in input.minhopp'
        close(11)
        if(mproc<nproc) then
            write(*,'(a)') 'WARNING: number of sets of input parameters fewer than nproc.'
        endif
        do jproc=mproc,nproc-1
            ediffarr(jproc)=ediffarr(mod(jproc,mproc))
            ekinarr(jproc)=ekinarr(mod(jproc,mproc))
            dtarr(jproc)=dtarr(mod(jproc,mproc))
        enddo
    endif
#if defined(MPI)
    if(parallel) then
        call MPI_BCAST(ediffarr,nproc,MPI_DOUBLE_PRECISION,imaster,mpi_comm_abz,ierr)
        call MPI_BCAST(ekinarr,nproc,MPI_DOUBLE_PRECISION,imaster,mpi_comm_abz,ierr)
        call MPI_BCAST(dtarr,nproc,MPI_DOUBLE_PRECISION,imaster,mpi_comm_abz,ierr)
    endif
#endif
    ediff=ediffarr(iproc);ekin=ekinarr(iproc);dt=dtarr(iproc)
    do  i=1,nrandoff+iproc*1000
        if(trim(parini%rng_type)=='only_for_tests') then
            call random_number_generator_simple(tts)
        else
            call random_number(tts)
        endif
    enddo
    call yaml_mapping_open('other minhopp parameters',flow=.true.)
    call yaml_map('iproc',iproc)
    call yaml_map('ediff',ediff)
    call yaml_map('ekin',ekin)
    call yaml_map('dt',dt)
    call yaml_mapping_close()
    !write(*,'(a,i4,3e15.3)') 'In: iproc,ediff,ekin,dt',iproc,ediff,ekin,dt
end subroutine read_minhopp_parameters
!*****************************************************************************************
subroutine minhopp_newrun_initialization(atoms_curr,atoms_locmin)
    use mod_minhopp, only: nlmin, nlmin_l, earr, nvisit
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_copy
    implicit none
    type(typ_atoms), intent(in):: atoms_curr
    type(typ_atoms_arr), intent(inout):: atoms_locmin
    nlmin=1
    nlmin_l=1
    earr(1)=atoms_curr%epot
    nvisit(1)=1
    !elocmin(1)=atoms_curr%epot
    !atoms_locmin%atoms%units=trim(atoms_curr%units)
    !atoms_locmin%epotall(1)=atoms_curr%epot
    !poslocmin(1:3,1:atoms_curr%nat,1)=atoms_curr%rat(1:3,1:atoms_curr%nat)
    !atoms_locmin%ratall(1:3,1:atoms_curr%nat,1)=atoms_curr%rat(1:3,1:atoms_curr%nat)
    call atom_copy(atoms_curr,atoms_locmin%atoms(1),'atoms_curr->atoms_locmin%atoms')
    !npmin=1
    atoms_locmin%nconf=1
    call send_minhopp_parameters_to_all(atoms_curr)
end subroutine minhopp_newrun_initialization
!*****************************************************************************************
subroutine read_poslow(atoms_locmin)
    use mod_processors, only: iproc, nproc, imaster, parallel, mpi_comm_abz
    use mod_minhopp, only: nlmin
    use mod_atoms, only: typ_atoms_arr
    implicit none
    type(typ_atoms_arr), intent(inout):: atoms_locmin
    !local variables
    integer:: natposlow, k, iat, ios, ierr
    character(20):: tatomnames
    integer:: nconf
#if defined(MPI)
    include 'mpif.h'
#endif
    !if(iproc==imaster) then
    !    open(unit=9,file='poslow.xyz',status='old',iostat=ios)
    !    if(ios/=0) then;write(*,'(a)') 'master processor could not read poslow.xyz';return;endif
    !    do nconf=1,min(atoms_locmin%nconfmax,nlmin)
    !        read(9,*,iostat=k) natposlow,atoms_locmin%epotall(nconf)
    !        if(k<0) then;write(*,*) 'k=',k;exit;endif
    !        if(atoms_locmin%atoms%nat/=natposlow) &
    !            write(*,'(a)') 'WARNING: number of atoms in poslow is incorrect.'
    !        read(9,*,iostat=k)
    !        if(k<0) then;write(*,*) 'k=',k;exit;endif
    !        do iat=1,atoms_locmin%atoms%nat
    !            read(9,*,iostat=k) tatomnames,atoms_locmin%ratall(1:3,iat,nconf)
    !            if(k<0) then;write(*,*) 'k=',k;exit;endif
    !        enddo
    !        if(k<0) then;write(*,*) 'k=',k;exit;endif
    !    enddo
    !    nconf=nconf-1
    !    if(nlmin>=atoms_locmin%nconfmax .and. nconf<atoms_locmin%nconfmax) then
    !        write(*,'(a)') 'WARNING: in this restart run npminx is increased.'
    !    endif
    !    close(9)
    !    atoms_locmin%nconf=nconf
    !endif
#if defined(MPI)
    if(parallel) then
        call MPI_BARRIER(mpi_comm_abz,ierr)
        call MPI_BCAST(atoms_locmin%nconf,1,MPI_INTEGER,imaster,mpi_comm_abz,ierr)
        !call MPI_BCAST(atoms_locmin%ratall,3*atoms_locmin%atoms%nat*atoms_locmin%nconf, &
        !    MPI_DOUBLE_PRECISION,imaster,mpi_comm_abz,ierr)
        !call MPI_BCAST(atoms_locmin%epotall,atoms_locmin%nconf,MPI_DOUBLE_PRECISION,imaster, &
        !    mpi_comm_abz,ierr)
        call MPI_atom_arr_copy(atoms_locmin%atoms(1)%nat,atoms_locmin)
    endif
#endif
end subroutine read_poslow
!*****************************************************************************************
subroutine send_minimum_to_all(atoms_curr)
    use mod_minhopp, only: nbuf, abufall, mtagarr1
    use mod_processors, only: nproc, iproc, mpi_comm_abz
    use mod_atoms, only: typ_atoms, get_rat
    implicit none
    type(typ_atoms), intent(in):: atoms_curr
    !local variables
    integer::ireq,ierr,jproc,iat,ibuf,itag
#if defined(MPI)
    include 'mpif.h'
#endif
    ibuf=mod(mtagarr1(iproc),nbuf)
    call get_rat(atoms_curr,abufall(1,ibuf))
    abufall(3*atoms_curr%nat+1,ibuf)=atoms_curr%epot
    itag=mtagarr1(iproc)
#if defined(MPI)
    do jproc=0,nproc-1
        if(jproc/=iproc) then
            call MPI_ISEND(abufall(1,ibuf),3*atoms_curr%nat+1,MPI_DOUBLE_PRECISION,jproc,itag,mpi_comm_abz, &
                ireq,ierr)
            write(*,'(a,2i4,i7,e24.15)') 'minimum sent:     iproc,jproc  ,itag,epossent    ', &
                iproc,jproc,itag,atoms_curr%epot
        endif
    enddo
#endif
    mtagarr1(iproc)=mtagarr1(iproc)+2
end subroutine send_minimum_to_all
!*****************************************************************************************
subroutine send_minhopp_parameters_to_all(atoms_curr)
    use mod_minhopp, only: nbuf, abufall, mtagarr2, ediff, ekin, dt
    use mod_processors, only: nproc, iproc, mpi_comm_abz
    use mod_atoms, only: typ_atoms, get_rat
    implicit none
    type(typ_atoms), intent(in):: atoms_curr
    !local variables
    integer::ireq,ierr,jproc,iat,itag,ibuf
#if defined(MPI)
    include 'mpif.h'
#endif
    ibuf=mod(mtagarr2(iproc),nbuf)
    call get_rat(atoms_curr,abufall(1,ibuf))
    abufall(3*atoms_curr%nat+1,ibuf)=atoms_curr%epot
    abufall(3*atoms_curr%nat+2,ibuf)=ediff
    abufall(3*atoms_curr%nat+3,ibuf)=ekin
    abufall(3*atoms_curr%nat+4,ibuf)=dt
#if defined(MPI)
    do jproc=0,nproc-1
        if(jproc/=iproc) then
            itag=mtagarr2(iproc)
            call MPI_ISEND(abufall(1,ibuf),3*atoms_curr%nat+4,MPI_DOUBLE_PRECISION,jproc,itag,mpi_comm_abz,ireq,ierr)
            !write(*,'(a,2i4,i7,e24.15)') 'data sent:     iproc,jproc  ,itag,erat              ', &
            !    iproc,jproc,itag,atoms_curr%epot
        endif
    enddo
#endif
    mtagarr2(iproc)=mtagarr2(iproc)+2
end subroutine send_minhopp_parameters_to_all
!*****************************************************************************************
subroutine mdescape(parini,atoms_hopp)
    use mod_parini, only: typ_parini
    use mod_minhopp, only: av_ekinetic, mdmin, dt, ekin, istep, count_md, count_md_tot, istep
    use mod_processors, only: iproc, nproc
    use mod_bin, only: write_bin_conf
    use mod_atoms, only: typ_atoms, typ_file_info, update_rat, update_ratp
    use mod_potential, only: fcalls
    use yaml_output
    !Does a MD run with the atomic positions in atoms_hopp
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms_hopp
    !local variables
    type(typ_file_info):: file_info
    integer:: nummax, nummin, iat, imd, istat, itt, imd_ideal, imd_maxdev
    real(8):: en0000, econs_max, econs_min, rkin, at1, at2, at3, devcon, enmin1, enmin2
    real(8):: tt1, tt2, tt
    logical:: mdpresumablyescaped
    real(8):: count_tt
    character(20):: filename
    character(100):: comment
    real(8), allocatable:: fatwa(:,:) !work array for atomic forces
    real(8):: epot0
    character(32):: ttfn
    !write(*,*) 'enter mdescape. iproc,imd',iproc,imd
    allocate(fatwa(3,atoms_hopp%nat),stat=istat)
    if(istat/=0) write(*,*) 'ERROR: allocating array fatwa'
    call velopt(parini,atoms_hopp) !initialize velocities
    !vat(1:3,1:atoms_hopp%nat)=atoms_hopp%vat(1:3,1:atoms_hopp%nat)
    fatwa(1:3,1:atoms_hopp%nat)=0.d0
    nummax=0;nummin=0;enmin1=0.d0;en0000=0.d0;econs_max=-1.d100;econs_min=1.d100
    mdpresumablyescaped=.false.
    count_tt=fcalls
    call setpot_mdescape
    !tt1=0.5d0*dt;tt2=dt*tt1
    if(parini%trajectory_minhopp) then
        write(ttfn,'(a9,i3.3,a5,i7.7,a8)') 'traj_proc',iproc,'_step',istep,'_mde.bin'
        file_info%filename_positions=ttfn
        file_info%file_position='new'
        file_info%print_force=parini%print_force_minhopp
        !call acf_write(file_info,atoms=atoms_hopp,strkey='mdescape')
        call write_bin_conf(file_info,atoms_hopp,strkey='mdescape')
    endif
    call cal_potential_forces(parini,atoms_hopp)
    epot0=atoms_hopp%epot
    call yaml_sequence_open('MD escape')
    do imd=1,10000
        !Evolution of the system according to 'VELOCITY VERLET' algorithm
        rkin=0.d0
        call update_ratp(atoms_hopp)
        do iat=1,atoms_hopp%nat
            tt2=0.5d0*dt**2/atoms_hopp%amass(iat)
            atoms_hopp%ratp(1,iat)=atoms_hopp%ratp(1,iat)+dt*atoms_hopp%vat(1,iat)+tt2*fatwa(1,iat)
            atoms_hopp%ratp(2,iat)=atoms_hopp%ratp(2,iat)+dt*atoms_hopp%vat(2,iat)+tt2*fatwa(2,iat)
            atoms_hopp%ratp(3,iat)=atoms_hopp%ratp(3,iat)+dt*atoms_hopp%vat(3,iat)+tt2*fatwa(3,iat)
            rkin=rkin+atoms_hopp%amass(iat)*(atoms_hopp%vat(1,iat)**2+atoms_hopp%vat(2,iat)**2+atoms_hopp%vat(3,iat)**2)
            !if(int(fcalls)==27) write(*,'(a,f6.2,3f10.5)') 'MASS ',amass(iat),vat(1:3,iat)
        enddo
        call update_rat(atoms_hopp)
        rkin=rkin*.5d0
        enmin2=enmin1
        enmin1=en0000
        if(parini%trajectory_minhopp) then
            file_info%file_position='append'
            !call acf_write(file_info,atoms=atoms_hopp,strkey='mdescape')
            call write_bin_conf(file_info,atoms_hopp,strkey='mdescape')
        endif
        call cal_potential_forces(parini,atoms_hopp)
        en0000=atoms_hopp%epot-epot0
        if(enmin1>enmin2 .and. enmin1>en0000)  nummax=nummax+1
        if(enmin1<enmin2 .and. enmin1<en0000)  nummin=nummin+1
        econs_max=max(econs_max,rkin+atoms_hopp%epot)
        econs_min=min(econs_min,rkin+atoms_hopp%epot)
        if(parini%iverbose>=0) then
            !write(*,'(a3,i3.3,a,i10,2e15.6,2i3)') 'MD:',iproc,' imd,en0000,rkin,nummax,nummin', &
            !    imd,en0000,rkin,nummax,nummin
            call yaml_sequence(advance='no')
            call yaml_mapping_open('MDEs',flow=.true.)
            !call yaml_map('iproc',iproc,fmt='(i3.3)')
            call yaml_map('iter',imd,fmt='(i3.3)')
            call yaml_map('en0000',en0000,fmt='(e14.5)')
            call yaml_map('rkin',rkin,fmt='(e14.5)')
            call yaml_map('nmax',nummax,fmt='(i3)')
            call yaml_map('nmin',nummin,fmt='(i3)')
            call yaml_mapping_close()
        endif
        if(nummin>=mdmin) then
            if(nummax/=nummin) write(*,'(a,3i4)') 'WARNING: iproc,nummin,nummax',iproc,nummin,nummax
            mdpresumablyescaped=.true.
            exit
        endif
        do iat=1,atoms_hopp%nat
            tt1=0.5d0*dt/atoms_hopp%amass(iat)
            at1=atoms_hopp%fat(1,iat)
            at2=atoms_hopp%fat(2,iat)
            at3=atoms_hopp%fat(3,iat)
            if(atoms_hopp%bemoved(1,iat)) atoms_hopp%vat(1,iat)=atoms_hopp%vat(1,iat)+tt1*(at1+fatwa(1,iat))
            if(atoms_hopp%bemoved(2,iat)) atoms_hopp%vat(2,iat)=atoms_hopp%vat(2,iat)+tt1*(at2+fatwa(2,iat))
            if(atoms_hopp%bemoved(3,iat)) atoms_hopp%vat(3,iat)=atoms_hopp%vat(3,iat)+tt1*(at3+fatwa(3,iat))
            fatwa(1,iat)=at1
            fatwa(2,iat)=at2
            fatwa(3,iat)=at3
        enddo
    enddo !end of loop over imd
    call yaml_sequence_close()
    devcon=econs_max-econs_min
    !adjust time step to meet precision criterion
    imd_ideal=30+(mdmin-1)*30
    imd_maxdev=mdmin*5
    itt=min(max(imd-imd_ideal,-imd_maxdev),imd_maxdev)
    !tt=1.d-2*(real(itt,8)+50.d0)+0.5d0
    tt=real(itt,8)/real(imd_maxdev,8)
    if(mdpresumablyescaped) then
        devcon=devcon/real(atoms_hopp%ndof,8)
        if(devcon/ekin<5.d-3) then
            dt=dt*1.20d0
        elseif(devcon/ekin>5.d-1) then
            dt=dt/2.d0
        else
            dt=dt*(1.d0+tt*0.1d0)
        endif
    else
        write(*,'(a)') 'TOO MANY MD STEPS'
        dt=2.d0*dt
    endif
    av_ekinetic=av_ekinetic+ekin
    !istep=istep+1
    count_md=fcalls-count_tt
    count_md_tot=count_md_tot+count_md
    deallocate(fatwa,stat=istat)
    if(istat/=0) write(*,*) 'ERROR: deallocating array fatwa'
end subroutine mdescape
!*****************************************************************************************
!MPI_IRECV(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, REQUEST, IERROR)
!MPI_RECV (BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, STATUS, IERROR) 
subroutine collect_data_from_all_processors(ntry,atoms_curr,atoms_allproc,atoms_locmin)
    use mod_processors, only: parallel, iproc
    use mod_minhopp, only: ekinarr, ediffarr, dtarr, dt, ekin, ediff
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_copy
    implicit none
    integer, intent(in):: ntry
    type(typ_atoms), intent(in):: atoms_curr
    type(typ_atoms_arr), intent(inout):: atoms_allproc, atoms_locmin
    !local variables
    integer:: itry
    !integer, save:: icall=0
    !icall=icall+1
    if(.not. parallel) return
    do itry=1,ntry
        call request_receive(atoms_allproc)
        call test_receive(atoms_allproc,atoms_locmin)
    enddo
    !ratallproc(1:3,1:atoms_curr%nat,iproc)=atoms_curr%rat(1:3,1:atoms_curr%nat)
    !atoms_allproc%ratall(1:3,1:atoms_curr%nat,iproc+1)=atoms_curr%rat(1:3,1:atoms_curr%nat)
    call atom_copy(atoms_curr,atoms_allproc%atoms(iproc+1),'atoms_curr->atoms_allproc%atoms')
    !if(iproc==0) then
    !    do iat=1,nat
    !    write(100+icall,'(3e24.15)') ratallproc(1:3,iat,0)
    !    enddo
    !endif
    !atoms_allproc%epotall(iproc+1)=atoms_curr%epot
    ediffarr(iproc)=ediff
    ekinarr(iproc)=ekin
    dtarr(iproc)=dt
end subroutine collect_data_from_all_processors
!*****************************************************************************************
subroutine request_receive(atoms_allproc)
    use mod_processors, only: iproc, nproc, mpi_comm_abz
    use mod_minhopp, only: mtagarr1, mtagarr2, abuf1, abuf2, do_req1, do_req2, ireqarr1, ireqarr2
    use mod_atoms, only: typ_atoms_arr
    implicit none
    type(typ_atoms_arr), intent(in):: atoms_allproc
    !local variables
    integer:: jproc, ierr
#if defined(MPI)
    include 'mpif.h'
    do jproc=0,nproc-1
        if(jproc==iproc) cycle
        if(do_req1(jproc)) then
        !call MPI_IRECV(abuf1(1,jproc),3*atoms_allproc%atoms%nat+1,MPI_DOUBLE_PRECISION,jproc,mtagarr1(jproc), &
        !    mpi_comm_abz, ireqarr1(jproc),ierr)
        call MPI_IRECV(abuf1(1,jproc),3*atoms_allproc%atoms(jproc)%nat+1,MPI_DOUBLE_PRECISION,jproc,mtagarr1(jproc), &
            mpi_comm_abz, ireqarr1(jproc),ierr)
        do_req1(jproc)=.false.
        endif
    enddo
    do jproc=0,nproc-1
        if(jproc==iproc) cycle
        if(do_req2(jproc)) then
        !call MPI_IRECV(abuf2(1,jproc),3*atoms_allproc%atoms%nat+4,MPI_DOUBLE_PRECISION,jproc,mtagarr2(jproc), &
        !    mpi_comm_abz, ireqarr2(jproc),ierr)
        call MPI_IRECV(abuf2(1,jproc),3*atoms_allproc%atoms(jproc)%nat+4,MPI_DOUBLE_PRECISION,jproc,mtagarr2(jproc), &
            mpi_comm_abz, ireqarr2(jproc),ierr)
        do_req2(jproc)=.false.
        endif
    enddo
#endif
end subroutine request_receive
!*****************************************************************************************
subroutine test_receive(atoms_allproc,atoms_locmin)
    use mod_processors, only: iproc, nproc
    use mod_minhopp, only: nlmin, nlminx, abuf1, abuf2, etoler, earr, re_sm, abuf, mtagarr1, mtagarr2, &
        nvisit, do_req1, do_req2, ireqarr1, ireqarr2, nbuf, ekinarr, ediffarr, dtarr
    use mod_atoms, only: typ_atoms_arr, set_rat
    implicit none
    type(typ_atoms_arr), intent(inout):: atoms_allproc, atoms_locmin
    !local variables
    integer:: ierr, jproc, itag, iat
    logical:: flag
    real(8):: eposreceived   
    integer:: ieposreceived
#if defined(MPI)
    include 'mpif.h'
    integer:: status_mpi(MPI_STATUS_SIZE)
    do jproc=0,nproc-1
        if(jproc==iproc) cycle
        !if(do_req1(jproc)) then
        call MPI_TEST(ireqarr1(jproc),flag,status_mpi,ierr)
        if(.not. flag) cycle
        itag=status_mpi(MPI_TAG)
        if(itag/=mtagarr1(jproc) .or. jproc/=status_mpi(MPI_SOURCE)) then
            write(*,*) 'ERROR: iproc,jproc ',iproc,jproc
        endif
        !eposreceived=abuf1(3*atoms_allproc%atoms%nat+1,jproc)
        eposreceived=abuf1(3*atoms_allproc%atoms(jproc+1)%nat+1,jproc)
        re_sm=min(re_sm,eposreceived)
        write(*,'(a,2i4,i7,e24.15)') 'minimum recEIVED: iproc,jproc,itag,eposreceived', &
            iproc,jproc,itag,eposreceived
        call hunt2(min(nlmin,nlminx+nbuf),earr(1),eposreceived,ieposreceived)
        if(abs(eposreceived-earr(ieposreceived))>etoler) then
            !write(*,'(a,i4)') 'new minimum: iproc',iproc
            nlmin=nlmin+1
            if(nlmin>=nlminx+nbuf) return
            call insert_alborz(ieposreceived,eposreceived)
            !call save_low_conf_alborz(atoms_allproc%atoms%nat,abuf1(1,jproc),eposreceived,atoms_locmin)
            stop 'ERROR: fix this line before parallel runs' !call save_low_conf_alborz(atoms_allproc%atoms(jproc+1)%nat,abuf1(1,jproc),eposreceived,atoms_locmin)
        else
            nvisit(ieposreceived)=nvisit(ieposreceived)+1
        endif
        mtagarr1(jproc)=mtagarr1(jproc)+2
        do_req1(jproc)=.true.
        !endif
    enddo
    do jproc=0,nproc-1
        if(jproc==iproc) cycle
        !if(do_req2(jproc)) then
        call MPI_TEST(ireqarr2(jproc),flag,status_mpi,ierr)
        !if(.not. flag .and. nproc>1) cycle
        if(.not. flag) cycle
        itag=status_mpi(MPI_TAG)
        if(itag/=mtagarr2(jproc) .or. jproc/=status_mpi(MPI_SOURCE)) then
            write(*,*) 'ERROR: iproc,jproc ',iproc,jproc
        endif
        !if(nproc>1) then
        !do iat=1,atoms_allproc%atoms%nat
        call set_rat(atoms_allproc%atoms(jproc+1),abuf2(1,jproc),setall=.true.)
        !atoms_allproc%epotall(jproc+1)=abuf2(3*atoms_allproc%atoms%nat+1,jproc)
        atoms_allproc%atoms(jproc+1)%epot=abuf2(3*atoms_allproc%atoms(jproc+1)%nat+1,jproc)
        !ediffarr(jproc)=abuf2(3*atoms_allproc%atoms%nat+2,jproc)
        !ekinarr(jproc)=abuf2(3*atoms_allproc%atoms%nat+3,jproc)
        !dtarr(jproc)=abuf2(3*atoms_allproc%atoms%nat+4,jproc)
        ediffarr(jproc)=abuf2(3*atoms_allproc%atoms(jproc+1)%nat+2,jproc)
        ekinarr(jproc)=abuf2(3*atoms_allproc%atoms(jproc+1)%nat+3,jproc)
        dtarr(jproc)=abuf2(3*atoms_allproc%atoms(jproc+1)%nat+4,jproc)
        !write(*,'(a,2i4,i7,e24.15)') 'data received: iproc,jproc,itag,atoms_allproc%epotall(jproc+1)', &
        !    iproc,jproc,itag,atoms_allproc%epotall(jproc+1)
        !mtagarr2(jproc)=itag
        mtagarr2(jproc)=mtagarr2(jproc)+2
        do_req2(jproc)=.true.
        !endif
    enddo
#endif
end subroutine test_receive
!*****************************************************************************************
!MPI_CANCEL(REQUEST, IERROR)
subroutine cancel_excessive_irecv
    use mod_processors, only: iproc, nproc
    use mod_minhopp, only: mtagarr1, mtagarr2, do_req1, do_req2, ireqarr1, ireqarr2
    implicit none
    !local variables
    integer:: ierr, jproc
    logical:: flag
#if defined(MPI)
    include 'mpif.h'
    integer:: status_mpi(MPI_STATUS_SIZE)
    do jproc=0,nproc-1
        if(jproc==iproc) cycle
        if(do_req1(jproc)) then
            flag=.true.
        else
            call MPI_TEST(ireqarr1(jproc),flag,status_mpi,ierr)
        endif
        if(.not. flag) then
            call MPI_CANCEL(ireqarr1(jproc),ierr)
        endif
        if(do_req2(jproc)) then
            flag=.true.
        else
            call MPI_TEST(ireqarr2(jproc),flag,status_mpi,ierr)
        endif
        if(.not. flag) then
            call MPI_CANCEL(ireqarr2(jproc),ierr)
        endif
    enddo
#endif
end subroutine cancel_excessive_irecv
!*****************************************************************************************
subroutine insert_alborz(kepos,epos)
    !inserts the energy epos at position kepos and shifts up all other energies
    use mod_minhopp, only: nbuf, nlmin, earr, nvisit
    implicit none
    integer, intent(in):: kepos
    real(8), intent(in):: epos
    !local variables
    integer:: k
    do k=nlmin-1,kepos+1,-1
        earr(k+1)=earr(k)
        nvisit(k+1)=nvisit(k)
    enddo
    earr(kepos+1)=epos
    nvisit(kepos+1)=1
end subroutine insert_alborz
!*****************************************************************************************
subroutine save_low_conf_alborz(atoms,atoms_locmin)
    !save configuration if it is among the lowest ones in energy
    use mod_minhopp, only: etoler
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_copy
    use mod_processors, only: iproc
    implicit none
    type(typ_atoms), intent(in):: atoms
    type(typ_atoms_arr), intent(inout):: atoms_locmin
    !local variables
    integer:: k, j
    if(atoms_locmin%nconf<atoms_locmin%nconfmax) then
        do k=1,atoms_locmin%nconf
            !if(epos<atoms_locmin%epotall(k)-etoler) then
            if(atoms%epot<atoms_locmin%atoms(k)%epot-etoler) then
                do j=atoms_locmin%nconf,k,-1
                    !atoms_locmin%ratall(1:3,1:nat,j+1)=atoms_locmin%ratall(1:3,1:nat,j)
                    !atoms_locmin%epotall(j+1)=atoms_locmin%epotall(j)
                    call atom_copy(atoms_locmin%atoms(j),atoms_locmin%atoms(j+1),'atoms_locmin%atoms(j)->atoms_locmin%atoms(j+1)')
                enddo
                !atoms_locmin%ratall(1:3,1:nat,k)=pos(1:3,1:nat)
                !atoms_locmin%epotall(k)=epos
                !atoms_locmin%atoms(k)%rat(1:3,1:nat)=pos(1:3,1:nat)
                !atoms_locmin%atoms(k)%epot=epos
                call atom_copy(atoms,atoms_locmin%atoms(k),'atoms->atoms_locmin%atoms(k)')
                atoms_locmin%nconf=atoms_locmin%nconf+1
                exit
            endif
        enddo
        if(k==atoms_locmin%nconf+1) then
            atoms_locmin%nconf=atoms_locmin%nconf+1
            !atoms_locmin%epotall(atoms_locmin%nconf)=epos
            !atoms_locmin%ratall(1:3,1:nat,atoms_locmin%nconf)=pos(1:3,1:nat)
            !atoms_locmin%atoms(atoms_locmin%nconf)%epot=epos
            !atoms_locmin%atoms(atoms_locmin%nconf)%rat(1:3,1:nat)=pos(1:3,1:nat)
            call atom_copy(atoms,atoms_locmin%atoms(k),'atoms->atoms_locmin%atoms(k)')
        endif
    else
        do k=1,min(atoms_locmin%nconf+1,atoms_locmin%nconfmax)
            !if(epos<atoms_locmin%epotall(k)-etoler) then
            if(atoms%epot<atoms_locmin%atoms(k)%epot-etoler) then
                do j=atoms_locmin%nconfmax-1,k,-1
                    !atoms_locmin%ratall(1:3,1:nat,j+1)=atoms_locmin%ratall(1:3,1:nat,j)
                    !atoms_locmin%epotall(j+1)=atoms_locmin%epotall(j)
                    !atoms_locmin%atoms(j+1)%rat(1:3,1:nat)=atoms_locmin%atoms(j)%rat(1:3,1:nat)
                    !atoms_locmin%atoms(j+1)%epot=atoms_locmin%atoms(j)%epot
                    call atom_copy(atoms_locmin%atoms(j),atoms_locmin%atoms(j+1),'atoms_locmin%atoms(j)->atoms_locmin%atoms(j+1)')
                enddo
                !atoms_locmin%ratall(1:3,1:nat,k)=pos(1:3,1:nat)
                !atoms_locmin%epotall(k)=epos
                !atoms_locmin%atoms(k)%rat(1:3,1:nat)=pos(1:3,1:nat)
                !atoms_locmin%atoms(k)%epot=epos
                call atom_copy(atoms,atoms_locmin%atoms(k),'atoms->atoms_locmin%atoms(k)')
                exit
            endif
        enddo
    endif
end subroutine save_low_conf_alborz
!*****************************************************************************************
subroutine velopt(parini,atoms)
    use mod_parini, only: typ_parini
    !assigns velocities according to boltzmann distribution
    use mod_minhopp, only: ekin, ekin, dt, count_soften, count_soften_tot, nsoften !, vdamp
    use mod_atoms, only: typ_atoms, update_ratp
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    !local variables
    real(8):: sumvx, sumvy, sumvz, rkinsum, sclvel, rkin, totmass
    !real(8):: dx,dy, dz, r
    integer:: iat
    !call gausdist_alborz(3*nat,vat)
    if(trim(parini%rng_type)=='only_for_tests') then
        call randdist_simple(0.15d0,3*atoms%nat,atoms%vat)
    else
        call randdist(0.15d0,3*atoms%nat,atoms%vat)
    endif
    !do iat=1,atoms%nat
    !    dx=atoms%rat(1,iat)-4.45819d0/0.529177d0
    !    dy=atoms%rat(2,iat)-4.45819d0/0.529177d0
    !    dz=atoms%rat(3,iat)-4.45819d0/0.529177d0
    !    r=sqrt(dx**2+dy**2+dz**2)
    !    if(r>11.d0) atoms%vat(1:3,iat)=0.d0
    !enddo
    if(atoms%ndof==3*atoms%nat) then
        call elim_moment_alborz(atoms%nat,atoms%vat)
    else
        do iat=1,atoms%nat
            if(.not. atoms%bemoved(1,iat)) atoms%vat(1,iat)=0.d0
            if(.not. atoms%bemoved(2,iat)) atoms%vat(2,iat)=0.d0
            if(.not. atoms%bemoved(3,iat)) atoms%vat(3,iat)=0.d0
        enddo
    endif
    if(trim(atoms%boundcond)=='free') then
        call update_ratp(atoms)
        call elim_torque_reza_alborz(atoms%nat,atoms%ratp,atoms%vat)
    endif
    !do iat=1,nat
    !    vat(1,iat)=vat(1,iat)*vdamp(3*iat-2)
    !    vat(2,iat)=vat(2,iat)*vdamp(3*iat-1)
    !    vat(3,iat)=vat(3,iat)*vdamp(3*iat-0)
    !enddo
    call soften(parini,nsoften,atoms,count_soften,count_soften_tot)
    !do iat=1,nat
    !    vat(1,iat)=vat(1,iat)*vdamp(3*iat-2)
    !    vat(2,iat)=vat(2,iat)*vdamp(3*iat-1)
    !    vat(3,iat)=vat(3,iat)*vdamp(3*iat-0)
    !enddo
    !Computation of total momentum
    sumvx=0.d0;sumvy=0.d0;sumvz=0.d0
    totmass=0.d0
    do iat=1,atoms%nat
        sumvx=sumvx+atoms%vat(1,iat)*atoms%amass(iat)
        sumvy=sumvy+atoms%vat(2,iat)*atoms%amass(iat)
        sumvz=sumvz+atoms%vat(3,iat)*atoms%amass(iat)
        totmass=totmass+atoms%amass(iat)
    enddo
    sumvx=sumvx/totmass
    sumvy=sumvy/totmass
    sumvz=sumvz/totmass
    !Correction to initial velocities in order to have zero total initial momentum .
    do iat=1,atoms%nat
        atoms%vat(1,iat)=atoms%vat(1,iat)-sumvx
        atoms%vat(2,iat)=atoms%vat(2,iat)-sumvy
        atoms%vat(3,iat)=atoms%vat(3,iat)-sumvz
        if(.not. atoms%bemoved(1,iat)) atoms%vat(1,iat)=0.d0
        if(.not. atoms%bemoved(2,iat)) atoms%vat(2,iat)=0.d0
        if(.not. atoms%bemoved(3,iat)) atoms%vat(3,iat)=0.d0
    enddo
    !Kinetic energy of the random velocities
    rkinsum= 0.d0      
    do iat=1,atoms%nat
        !Here it is assumed if a degree of freedom is fixed, its velocity is zero.
        rkinsum= rkinsum+atoms%amass(iat)*(atoms%vat(1,iat)**2+atoms%vat(2,iat)**2+atoms%vat(3,iat)**2)
    enddo
    rkin=.5d0*rkinsum/(atoms%ndof-3)
    !Rescaling of velocities to get reference kinetic energy
    sclvel= dsqrt(ekin/rkin)
    do iat=1,atoms%nat
        atoms%vat(1,iat)=atoms%vat(1,iat)*sclvel
        atoms%vat(2,iat)=atoms%vat(2,iat)*sclvel
        atoms%vat(3,iat)=atoms%vat(3,iat)*sclvel
    enddo
end subroutine velopt
!*****************************************************************************************
subroutine soften(parini,nstep,atoms0,count_soften,count_soften_tot)
    use mod_parini, only: typ_parini
    use mod_processors, only: iproc
    use mod_potential, only: fcalls
    use mod_minhopp, only: lprint, alpha_soften, istep
    use mod_atoms, only: typ_atoms, typ_file_info, atom_deallocate, atom_copy
    use mod_bin, only: write_bin_conf
    use mod_atoms, only: update_ratp, update_rat
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: nstep
    type(typ_atoms), intent(inout):: atoms0
    real(8), intent(inout):: count_soften, count_soften_tot
    !local variables
    type(typ_file_info):: file_info
    type(typ_atoms):: atoms
    real(8):: alpha, svxyz, fd2, sdf, curv, curv0, res, etot0
    real(8):: eps_vxyz, eps_vxyz_init, count_tt
    !real(8), save:: eps_vxyz=2.d-2 !used for SOFC
    !real(8), save:: eps_vxyz=4.d-2 !used for Li clusters
    integer:: iter, i, iat, ixyz
    character(32):: ttfn
    count_tt=fcalls
    !eps_vxyz_init=1.d-1*atoms0%nat
    !alpha=2.d-3 !inputs_md%betax
    alpha=alpha_soften !0.8d0*alphax !2.d-2 !inputs_md%betax
    !if(lprint) write(*,'(a,i3)') 'softening begins ',iproc
    call setpot_soften
    !write(*,*) 'HERE ',trim(atoms0%boundcond),atoms0%ndof,atoms0%bemoved
    !perfstatus='normal'
    if(parini%trajectory_minhopp) then
        write(ttfn,'(a9,i3.3,a5,i7.7,a8)') 'traj_proc',iproc,'_step',istep,'_sof.bin'
        file_info%filename_positions=ttfn
        file_info%file_position='new'
        !call acf_write(file_info,atoms=atoms0,strkey='soften')
        call write_bin_conf(file_info,atoms0,strkey='soften')
    endif
    call cal_potential_forces(parini,atoms0)
    etot0=atoms0%epot
    !stop
    !call call_bigdft(nproc,iproc,atoms0,rxyz,inputs_md,etot0,atoms0%fat,strten,fnoise,rst,infocode)
    !scale velocity to generate dimer 
    !call atomic_dot(atoms0,vxyz,vxyz,svxyz)
    call atom_copy(atoms0,atoms,'atoms0->atoms')
    svxyz=0.d0
    do iat=1,atoms0%nat
        svxyz=svxyz+atoms0%vat(1,iat)**2+atoms0%vat(2,iat)**2+atoms0%vat(3,iat)**2
    enddo
    !eps_vxyz_init=sqrt(svxyz)
    !eps_vxyz_init=1.d-3 !LJ38
    !if(lprint) write(*,'(a,f12.5)') '#eps_vxyz= ',eps_vxyz
    !eps_vxyz_init=2.d-2 !used for SOFC with perfstatus='accurate'
    eps_vxyz_init=2.d-2 !used for VASP
    !eps_vxyz_init=10.d-2 !used for SIESTA
    eps_vxyz=eps_vxyz_init
    atoms0%vat(1:3,1:atoms0%nat)=atoms0%vat(1:3,1:atoms0%nat)/sqrt(svxyz)*eps_vxyz
    call yaml_sequence_open('SOFTEN')
    do iter=1,nstep
        call update_ratp(atoms0)
        write(*,*) atoms0%nat,atoms%nat
        do iat=1,atoms0%nat
            atoms%ratp(1,iat)=atoms0%ratp(1,iat)+atoms0%vat(1,iat)
            atoms%ratp(2,iat)=atoms0%ratp(2,iat)+atoms0%vat(2,iat)
            atoms%ratp(3,iat)=atoms0%ratp(3,iat)+atoms0%vat(3,iat)
        enddo
        call update_rat(atoms)
        if(parini%trajectory_minhopp) then
            file_info%file_position='append'
            !call acf_write(file_info,atoms=atoms,strkey='soften')
            call write_bin_conf(file_info,atoms,strkey='soften')
        endif
        call cal_potential_forces(parini,atoms)
        fd2=2.d0*(atoms%epot-etot0)/eps_vxyz**2
        sdf=0.d0
        svxyz=0.d0
        do iat=1,atoms0%nat
            sdf=sdf+atoms0%vat(1,iat)*atoms%fat(1,iat)+atoms0%vat(2,iat)*atoms%fat(2,iat)+atoms0%vat(3,iat)*atoms%fat(3,iat)
            svxyz=svxyz+atoms0%vat(1,iat)**2+atoms0%vat(2,iat)**2+atoms0%vat(3,iat)**2
        enddo
        curv=-sdf/svxyz
        if(iter==1) curv0=curv
        res=0.d0
        do i=1,3*atoms0%nat
            iat=(i-1)/3+1
            ixyz=mod(i-1,3)+1
            atoms%fat(ixyz,iat)=atoms%fat(ixyz,iat)+curv*atoms0%vat(ixyz,iat)
            res=res+atoms%fat(ixyz,iat)**2
        enddo
        res=sqrt(res)
        if(parini%iverbose>=0) then
            !write(*,'(a,i3.3,a,i3,5(f12.5))') 'SOFTEN:',iproc,' iter,curv,fd2,de,res,eps_vxyz: ', &
            !    iter,curv,fd2,atoms%epot-etot0,res,eps_vxyz
            call yaml_sequence(advance='no')
            call yaml_mapping_open('SOFT',flow=.true.)
            !call yaml_map('iproc',iproc,fmt='(i3.3)')
            call yaml_map('iter',iter,fmt='(i3)')
            call yaml_map('curv',curv,fmt='(f10.3)')
            call yaml_map('fd2',fd2,fmt='(f10.3)')
            call yaml_map('de',atoms%epot-etot0,fmt='(f10.3)')
            call yaml_map('res',res,fmt='(f10.3)')
            call yaml_map('eps_v',eps_vxyz,fmt='(f10.3)')
            call yaml_mapping_close()
        endif
        if(curv<0.d0 .or. fd2<0.d0) then
            if(lprint) write(*,'(a,i3)') 'WARNING: negative curvature. iproc=',iproc
            exit
        endif
        if(atoms%epot-etot0<1.d-3) eps_vxyz=min(1.2d0*eps_vxyz,2.d0*eps_vxyz_init) !used for VASP
        !if(atoms%epot-etot0<1.d-2) eps_vxyz=eps_vxyz*1.1d0 !used for SIESTA
        call update_ratp(atoms)
        call update_ratp(atoms0)
        do i=1,3*atoms0%nat
            iat=(i-1)/3+1
            ixyz=mod(i-1,3)+1
            if(atoms0%bemoved(ixyz,iat)) then
                atoms%ratp(ixyz,iat)=atoms%ratp(ixyz,iat)+alpha*atoms%fat(ixyz,iat)
                atoms0%vat(ixyz,iat)=atoms%ratp(ixyz,iat)-atoms0%ratp(ixyz,iat)
            endif
        enddo
        call update_rat(atoms)
        if(atoms0%ndof==3*atoms0%nat) then
            call elim_moment_alborz(atoms0%nat,atoms0%vat)
        else
            !I am not sure what would be the best to do here.
        endif
        if(trim(atoms0%boundcond)=='free') then
            call elim_torque_reza_alborz(atoms0%nat,atoms0%ratp,atoms0%vat)
        endif
        svxyz=0.d0
        do iat=1,atoms0%nat
            svxyz=svxyz+atoms0%vat(1,iat)**2+atoms0%vat(2,iat)**2+atoms0%vat(3,iat)**2
        enddo
        if(res<=curv*eps_vxyz*5.d-2) exit
        svxyz=eps_vxyz/sqrt(svxyz)
        do i=1,3*atoms0%nat
            iat=(i-1)/3+1
            ixyz=mod(i-1,3)+1
            atoms0%vat(ixyz,iat)=atoms0%vat(ixyz,iat)*svxyz
        enddo
    enddo !end of loop over iter
    call yaml_sequence_close()
    count_soften=fcalls-count_tt
    count_soften_tot=count_soften_tot+count_soften
    call atom_deallocate(atoms)
end subroutine soften
!*****************************************************************************************
subroutine write_minhopp(atoms_allproc,atoms_locmin)
    use mod_atoms, only: typ_atoms_arr, typ_file_info
    use mod_processors, only: iproc, imaster
    use mod_minhopp, only: nlmin_old, nlmin, earr, etoler
    use mod_yaml_conf, only: write_yaml_conf
    use yaml_output
    implicit none
    type(typ_atoms_arr), intent(inout):: atoms_allproc, atoms_locmin
    !local variables
    type(typ_file_info):: file_info
    integer:: ipmin
    integer:: iconf, iat
    nlmin_old=nlmin
    if(iproc/=imaster) return
    call yaml_mapping_open('writing minhopp results',flow=.true.)
    call yaml_map('iproc',iproc)
    call yaml_map('nlmin',nlmin)
    call yaml_mapping_close()
    !write(*,'(a,i4,i7)') 'writing minhopp results: iproc,nlmin',iproc,nlmin
    !call write_poscur(atoms_allproc)
    !file_info%filename_positions='posout.acf'
    !file_info%file_position='new'
    !call acf_write(file_info,atoms_arr=atoms_allproc,strkey='poscur')
    !call acf_write_new(file_info,atoms_arr=atoms_allproc,strkey='poscur')
    file_info%filename_positions='posout.yaml'
    do iconf=1,atoms_allproc%nconf
        if(iconf==1) then
            file_info%file_position='new'
        else
            file_info%file_position='append'
        endif
        call write_yaml_conf(file_info,atoms_allproc%atoms(iconf),strkey='poscur')
    enddo
    !call write_poslow(atoms_locmin)
    do ipmin=1,atoms_locmin%nconf
        !if(abs(earr(ipmin)-atoms_locmin%epotall(ipmin))>etoler) then
        if(abs(earr(ipmin)-atoms_locmin%atoms(ipmin)%epot)>etoler) then
            write(*,'(a,i4)') 'ERROR: earr(ipmin) and elocmin(ipmin) are not equal: ipmin',ipmin
        endif
    enddo
    do ipmin=2,atoms_locmin%nconf
        !if(abs(atoms_locmin%epotall(ipmin-1)-atoms_locmin%epotall(ipmin))<etoler) then
        if(abs(atoms_locmin%atoms(ipmin-1)%epot-atoms_locmin%atoms(ipmin)%epot)<etoler) then
            write(*,'(a,2i4)') 'ERROR: duplication in elocmin array: ipmin-1,ipmin',ipmin-1,ipmin
        endif
    enddo
    do ipmin=2,atoms_locmin%nconf
        !if(.not. (atoms_locmin%epotall(ipmin-1)<atoms_locmin%epotall(ipmin)-etoler)) then
        if(.not. (atoms_locmin%atoms(ipmin-1)%epot<atoms_locmin%atoms(ipmin)%epot-etoler)) then
            write(*,'(a,2i4)') 'ERROR: ranking error in elocmin array: ipmin-1,ipmin',ipmin-1,ipmin
        endif
    enddo
    !file_info%filename_positions='poslow.acf'
    !file_info%file_position='new'
    !call acf_write_new(file_info,atoms_arr=atoms_locmin,strkey='poslow')
    file_info%filename_positions='poslow.yaml'
    do iconf=1,atoms_locmin%nconf
        if(iconf==1) then
            file_info%file_position='new'
        else
            file_info%file_position='append'
        endif
        call write_yaml_conf(file_info,atoms_locmin%atoms(iconf),strkey='poslow')
    enddo
    call yaml_comment('master proc wrote poslow.xyz')
    !write(*,*) 'master proc wrote poslow.xyz'
    call write_earr
    call write_minhopp_parameters
end subroutine write_minhopp
!*****************************************************************************************
subroutine write_minhopp_parameters
    use mod_processors, only: iproc, nproc, imaster
    use mod_minhopp, only: nlmin, ekinarr, dtarr, ediffarr, eref
    use yaml_output
    implicit none
    !local variables
    integer::ios, jproc
    real(8):: tt
    if(iproc==imaster) then
        open(unit=11,file='input.minhopp',status='replace',iostat=ios)
        if(ios/=0) stop 'ERROR: failure openning input.minhopp'
        write(11,'(i8,1x,a30)') nlmin,'number of minima already found'
        !write(11,'(i,1x,a)') 2*nlmin,'number of minima to be found in consecutive run'
        do jproc=0,nproc-1
            tt=dtarr(jproc) !/41.34137333656d0
            write(11,'(3(e15.6,1x),a)') ediffarr(jproc),ekinarr(jproc),tt,'  ediff,ekin,dt'
        enddo
        close(11)
        call yaml_comment('master proc wrote input/output parameters in input.minhopp.')
        !write(*,'(a)') 'master proc wrote input/output parameters in input.minhopp.'
    endif
end subroutine write_minhopp_parameters
!*****************************************************************************************
subroutine write_earr
    use mod_minhopp, only: earr, nvisit, nlmin, nlminx, nbuf, eref
    use yaml_output
    implicit none
    integer:: mm, k, ios
    real(8):: tt
    open(unit=12,file='earr.dat',status='replace',iostat=ios)
    if(ios/=0) stop 'ERROR: failure openning earr.dat'
    mm=min(nlmin,nlminx+nbuf)
    write(12,'(e24.15,i5,2f10.5,e14.5)') earr(1),nvisit(1),0.d0,0.d0,earr(1)-eref
    do k=2,mm
        tt=earr(k)
        write(12,'(e24.15,i5,2f10.5,e14.5)') tt,nvisit(k),tt-earr(1),tt-earr(k-1),tt-eref
    enddo
    call yaml_comment('master proc wrote energies in earr.dat, can be used in a restart run.')
    !write(*,'(a)') 'master proc wrote energies in earr.dat, can be used in a restart run.'
    close(12)
end subroutine write_earr
!*****************************************************************************************
subroutine escape_failed(parini,erat,erathopp)
    use mod_parini, only: typ_parini
    use mod_minhopp, only: istep_sam, esep, ekin, beta1, ihopp, istep, ediff, istep_old, &
        istep_new, dt, count_md, count_opt
    use mod_processors, only: iproc
    implicit none
    type(typ_parini), intent(in):: parini
    real(8), intent(in):: erat, erathopp
    !local variables
    real(8):: t1, t2, t3
    istep_sam=istep_sam+1
    esep=esep+(erat-erathopp)**2
    ekin=min(ekin*beta1,parini%ekinmax_minhopp)
    dt=dt/1.10d0
    t1=real(istep_sam,8)/real(istep,8)
    t2=real(istep_old,8)/real(istep,8)
    t3=real(istep_new,8)/real(istep,8)
    write(2,'(i4.3,i10,i6,e24.15,3f10.4,3f5.2,a,2f9.1)') iproc,istep,ihopp,erathopp,ediff,ekin,dt,t1,t2,t3,'  S-',count_md,count_opt
end subroutine escape_failed
!*****************************************************************************************
subroutine local_minimum_accepted(atoms_hopp,atoms_curr,atoms_locmin)
    use mod_minhopp, only: istep_sam, ekin, ihopp, istep, ediff, istep_old, istep_new, dt, &
        nlmin, nlmin_l, ihopp_acc, alpha1, nvisit, nbuf, kerathopp, newmin, av_ediff, count_md, count_opt
    use mod_processors, only: parallel, iproc
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_copy
    implicit none
    type(typ_atoms), intent(in):: atoms_hopp
    type(typ_atoms), intent(inout):: atoms_curr
    type(typ_atoms_arr), intent(inout):: atoms_locmin
    !local variables
    real(8):: t1, t2, t3
    t1=real(istep_sam,8)/real(istep,8)
    t2=real(istep_old,8)/real(istep,8)
    t3=real(istep_new,8)/real(istep,8)
    write(2,'(i4.3,i10,i6,e24.15,3f10.4,3f5.2,l3,a,2f9.1)') iproc,istep,ihopp,atoms_hopp%epot,ediff,ekin,dt,t1,t2,t3,newmin,'A',count_md,count_opt
    ihopp_acc=ihopp_acc+1
    ediff=ediff*alpha1 !standard feedback
    av_ediff=av_ediff+ediff
    !atoms_curr%epot=atoms_hopp%epot
    !atoms_curr%rat(1:3,1:nat)=rathopp(1:3,1:nat)
    !atoms_curr%rat(1:3,1:atoms_curr%nat)=atoms_hopp%rat(1:3,1:atoms_curr%nat)
    call atom_copy(atoms_hopp,atoms_curr,'atoms_hopp->atoms_curr')
    if(parallel) call send_minimum_to_all(atoms_curr)
    if(newmin) then !if new local minimum.
        nlmin=nlmin+1
        nlmin_l=nlmin_l+1
        !add minimum to history list
        call insert_alborz(kerathopp,atoms_hopp%epot)
        !save configuration if it is among the lowest ones in energy
        call save_low_conf_alborz(atoms_curr,atoms_locmin)  
    else !old minimum revisited.
        nvisit(kerathopp)=nvisit(kerathopp)+1
    endif
end subroutine local_minimum_accepted
!*****************************************************************************************
subroutine local_minimum_rejected(atoms_hopp)
    use mod_minhopp, only: istep_sam, ekin, ihopp, istep, ediff, istep_old, istep_new, &
        alpha2, nvisit, newmin, ihopp_rej, dt, count_md, count_opt
    use mod_processors, only: iproc
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(in):: atoms_hopp
    !local variables
    real(8)::t1,t2,t3
    t1=real(istep_sam,8)/real(istep,8)
    t2=real(istep_old,8)/real(istep,8)
    t3=real(istep_new,8)/real(istep,8)
    write(2,'(i4.3,i10,i6,e24.15,3f10.4,3f5.2,l3,a,2f9.1)') iproc,istep,ihopp,atoms_hopp%epot,ediff,ekin,dt,t1,t2,t3,newmin,'R',count_md,count_opt
    ihopp_rej=ihopp_rej+1
    ediff=ediff*alpha2 !standard feedback
    !ediff=ediff*alpha2**(1.d0+.1d0*log(real(nvisit(kerathopp),8))) !enhanced feedback
end subroutine local_minimum_rejected
!*****************************************************************************************
subroutine report_minhopp_iteration_info(atoms_curr)
    use mod_minhopp, only: istep, ihopp, nlmin, escaped, accepted
    use mod_processors, only: iproc
    use mod_atoms, only: typ_atoms, typ_atoms_arr
    use yaml_output
    !minima hopping program with restart option.
    implicit none
    type(typ_atoms), intent(in):: atoms_curr
    !local variables

    call yaml_mapping_open('minhopp iteration info',flow=.true.)
    call yaml_map('iproc',iproc,fmt='(i8)')
    call yaml_map('istep',istep,fmt='(i8)')
    call yaml_map('ihopp',ihopp,fmt='(i8)')
    call yaml_map('nlmin',nlmin,fmt='(i8)')
    call yaml_map('escaped',escaped)
    call yaml_map('accepted',accepted)
    call yaml_map('epot',atoms_curr%epot,fmt='(es20.12)')
    call yaml_mapping_close()
end subroutine report_minhopp_iteration_info
!*****************************************************************************************
subroutine already_visited_minimum(parini)
    use mod_parini, only: typ_parini
    use mod_minhopp, only: newmin, istep_old, ekin, beta2, nvisit, kerathopp , dt, ediff, alpha3
    implicit none
    type(typ_parini), intent(in):: parini
    newmin=.false.
    istep_old=istep_old+1
    ekin=min(ekin*beta2,parini%ekinmax_minhopp) !standard feedback
    !ekin=ekin*beta2*(1.d0+.33d0*log10(real(nvisit(kerathopp),8))) !enhanced feedback, Stefan
    !ekin=ekin*beta2*(1.d0+.10d0*log10(real(nvisit(kerathopp),8))) !enhanced feedback, Sandro
    !ekin=ekin*(1.d0+(beta2-1.d0)*(1.d0+1.d-1*real(nvisit(kerathopp),8))) !enhanced feedback
    dt=dt/1.15d0
    ediff=ediff*alpha3 !standard feedback added by Alireza
end subroutine already_visited_minimum
!*****************************************************************************************
subroutine new_minimum(atoms_hopp)
    use mod_minhopp, only: newmin, istep_new, ekin, beta3, kerathopp, nlminx, earr, egap, &
        re_sm,nlmin, etoler, dt
    use mod_processors, only:iproc
    use mod_bin, only: write_bin_conf
    use mod_atoms, only: typ_atoms, typ_file_info
    use mod_minhopp, only: ibest
    use yaml_output
    implicit none
    type(typ_atoms), intent(in):: atoms_hopp
    !local variables
    type(typ_file_info):: file_info
    real(8)::tt
    newmin=.true.
    istep_new=istep_new+1
    ekin=ekin*beta3
    dt=dt*1.10d0
    !determine energy separation between distinct minima
    if(kerathopp+1<=nlminx) then
        tt=min(atoms_hopp%epot-earr(kerathopp),earr(kerathopp+1)-atoms_hopp%epot)
        if(tt>etoler*(1.1d0)) egap=min(egap,tt)
        !write(*,*) 'tt,egap',tt,egap
    endif
    if(kerathopp==0) then
        re_sm=min(re_sm,atoms_hopp%epot)
        call yaml_mapping_open('new lowest',flow=.true.)
        call yaml_map('iproc',iproc)
        call yaml_map('nlmin',nlmin)
        call yaml_map('erathopp',atoms_hopp%epot,fmt='(es20.12)')
        call yaml_map('de',atoms_hopp%epot-earr(1),fmt='(f10.5)')
        call yaml_mapping_close()
        !write(*,'(a,i4,i7,e24.15,f10.5)') 'new lowest: iproc,nlmin,erathopp,dE', &
        !    iproc,nlmin,atoms_hopp%epot,atoms_hopp%epot-earr(1)
        file_info%filename_positions='posbest.bin'
        if(ibest==1) then
            file_info%file_position='new'
        else
            file_info%file_position='append'
        endif
        !call acf_write(file_info,atoms=atoms_hopp,strkey='posbest')
        call write_bin_conf(file_info,atoms_hopp,strkey='posbest')
        ibest=ibest+1
    endif
end subroutine new_minimum
!*****************************************************************************************
subroutine print_final_statistics
    use mod_processors, only: iproc
    use mod_task, only: time_start
    use mod_minhopp, only: nlmin_l, nlmin, istep_sam, istep_old, istep_new, &
        istep, count_md_tot, count_opt_tot, count_soften_tot, av_ediff, earr, eref, &
        ihopp_acc, ihopp_rej, ihopp, ekin, av_ekinetic, nvisit,egap, esep, dt, ediff, &
        fcall_tot_all, fcall_tot_all_md, fcall_tot_all_opt, fcall_tot_all_soften
    use mod_potential, only: fcalls
    use yaml_output
    implicit none
    !local variables
    integer:: i, itt, iss
    real(8):: tt1, tt2, t1, t2, t3
    logical:: reached
    call yaml_mapping_open('number of minima found by each processors',flow=.true.)
    call yaml_map('iproc',iproc)
    call yaml_map('nlmin_l',nlmin_l)
    call yaml_map('nlmin',nlmin)
    call yaml_mapping_close()
    !write(*,'(a,i4,2i7)') 'number of minima found by each processors: iproc,nlmin_l,nlmin',iproc,nlmin_l,nlmin
    !ratios from all the global counters
    if(istep>0 .and. ihopp>0) then
        t1=real(istep_sam,8)/real(istep,8)
        t2=real(istep_old,8)/real(istep,8)
        t3=real(istep_new,8)/real(istep,8)
        call yaml_mapping_open('ratio stuck',flow=.true.)
        call yaml_map('same',t1,fmt='(e11.3)')
        call yaml_map('old',t2,fmt='(e11.3)')
        call yaml_map('new',t3,fmt='(e11.3)')
        call yaml_mapping_close()
        !write(*,'(a,i4,3e15.3)') 'ratio stuck: iproc,same,old,new',iproc,t1,t2,t3
        tt1=real(ihopp_acc,8)/real(ihopp,8)
        tt2=real(ihopp_rej,8)/real(ihopp,8)
        call yaml_mapping_open('accepting/rejecting ratio',flow=.true.)
        call yaml_map('iproc',iproc)
        call yaml_map('accept',tt1,fmt='(e11.3)')
        call yaml_map('reject',tt2,fmt='(e11.3)')
        call yaml_mapping_close()
        !write(*,'(a,i4,2e15.3)') 'accepting/rejecting ratio: iproc,acc,rej',iproc,tt1,tt2 
    endif
    call yaml_mapping_open('force calls',flow=.true.)
    call yaml_map('total',fcalls)
    call yaml_map('MD',count_md_tot)
    call yaml_map('GEOPT',count_opt_tot)
    call yaml_map('SOFT',count_soften_tot)
    call yaml_mapping_close()
    !write(*,'(a,i4,4f15.1)') 'icount: iproc,total,md,geopt,soften ',&
    !    iproc,fcalls,count_md_tot,count_opt_tot,count_soften_tot
    if(iproc==0) then
        if(earr(1)<eref) then
            reached=.true.
        else
            reached=.false.
        endif
        call yaml_mapping_open('force calls all processes')
        call yaml_map('fcall_tot_all',fcall_tot_all)
        call yaml_map('fcall_tot_all_md',fcall_tot_all_md)
        call yaml_map('fcall_tot_all_opt',fcall_tot_all_opt)
        call yaml_map('fcall_tot_all_soften',fcall_tot_all_soften)
        call yaml_map('reached',reached)
        call yaml_mapping_close()
        !write(*,'(a,4f15.1,l3)') 'ICOUNT: TOTAL,MD,GEOPT,SOFTEN ',&
        !    fcall_tot_all,fcall_tot_all_md,fcall_tot_all_opt,fcall_tot_all_soften,reached
    endif
    if(istep>0 .and. ihopp>0) then
        av_ediff=av_ediff/real(ihopp,8)
        av_ekinetic=av_ekinetic/real(istep,8)
        call yaml_mapping_open('average minima hopping parameters',flow=.true.)
        call yaml_map('iproc',iproc)
        call yaml_map('av_ediff',av_ediff)
        call yaml_map('av_ekinetic',av_ekinetic)
        call yaml_mapping_close()
        !write(*,'(a,i4,2e15.3)') 'average minima hopping parameters: iproc,av_ediff,av_ekinetic', &
        !    iproc,av_ediff,av_ekinetic
    endif
    itt=0.d0;iss=0.d0
    do i=1,nlmin
        itt=max(itt,nvisit(i))
        iss=iss+nvisit(i)
    enddo
    call yaml_mapping_open('some minhopp results')
    call yaml_map('iproc',iproc)
    call yaml_map('most frequent visits',itt)
    call yaml_map('avg. numb. visits per minimum',real(iss,8)/real(nlmin,8),fmt='(e10.3)')
    call yaml_map('minimum energy separation',egap,fmt='(e10.2)')
    !write(*,'(a,i4,i8)') 'most frequent visits: iproc,most',iproc,itt
    !write(*,'(a,i4,e10.3)') 'avg. numb. visits per minimum: iproc,average',iproc,real(iss,8)/real(nlmin,8)
    !write(*,'(a,i4,e10.2)') 'minimum energy separation: iproc,egap',iproc,egap
    if(istep_sam>0 .and. nlmin>1) then
        esep=sqrt(esep/real(istep_sam,8))
        call yaml_map('esep',esep,fmt='(e10.2)')
        !write(*,'(a,i4,e10.2)') 'average energy separation: iproc,esep',iproc,esep
    endif
    call yaml_map('ediff',ediff,fmt='(e15.3)')
    call yaml_map('ekin',ekin,fmt='(e15.3)')
    call yaml_map('dt',dt,fmt='(e15.3)')
    !write(*,'(a,i4,3e15.3)') 'Out: iproc,ediff,ekin,dt',iproc,ediff,ekin,dt
    call yaml_mapping_close()
end subroutine print_final_statistics
!*****************************************************************************************
subroutine MPI_atom_arr_copy(nat,atoms_arr)
    use mod_processors, only: iproc, nproc, parallel, imaster, mpi_comm_abz
    use mod_atoms, only: typ_atoms, typ_atoms_arr, get_rat, set_rat
    implicit none
    integer,intent(in)::nat
    type(typ_atoms_arr),intent(inout):: atoms_arr
    real(8),allocatable::ratall(:,:,:),epotall(:,:),cvall(:,:,:)
    integer:: iconf,ierr
#if defined(MPI)
    include 'mpif.h'
#endif
    allocate(ratall(3,nat,atoms_arr%nconf),epotall(1,atoms_arr%nconf),cvall(3,3,atoms_arr%nconf))
    do iconf=1,atoms_arr%nconf
        call get_rat(atoms_arr%atoms(iconf),ratall(1,1,iconf))
        epotall(1,iconf)=atoms_arr%atoms(iconf)%epot
        cvall(1:3,1:3,iconf)=atoms_arr%atoms(iconf)%cellvec
    enddo
#if defined(MPI)
    call MPI_BCAST(ratall,3*nat*atoms_arr%nconf,MPI_DOUBLE_PRECISION,imaster,mpi_comm_abz,ierr)
    call MPI_BCAST(epotall,atoms_arr%nconf,MPI_DOUBLE_PRECISION,imaster,mpi_comm_abz,ierr)
    call MPI_BCAST(cvall,3*3*atoms_arr%nconf,MPI_DOUBLE_PRECISION,imaster,mpi_comm_abz,ierr)
#endif
    do iconf=1,atoms_arr%nconf
        !call atoms_copy(atoms_curr,atoms_arr%atoms(1),'atoms_arr%atoms->atoms_curr')
       ! call atoms_copy(atoms_curr,atoms_arr%atoms(1),'atoms_curr->atoms_arr%atoms')
        call set_rat(atoms_arr%atoms(iconf),ratall(1,1,iconf),setall=.true.)
        atoms_arr%atoms(iconf)%epot=epotall(1,iconf)
        atoms_arr%atoms(iconf)%cellvec=cvall(1:3,1:3,iconf)
    enddo
    deallocate(ratall,epotall,cvall)
end subroutine MPI_atom_arr_copy
!*****************************************************************************************
subroutine hunt2(n,x,p,ip)
    !p is in interval [x(ip),x(ip+1)[ ; x(0)=-Infinity ; x(n+1) = Infinity
    use mod_minhopp, only: etoler
    implicit none
    integer, intent(in):: n
    real(8), intent(in):: x(n), p
    integer, intent(out):: ip
    !local variables
    integer:: i
    do i=1,n
        if(p<x(i)) exit
    enddo
    i=i-1
    if(i/=n) then
        if(abs(p-x(i+1))<etoler) i=i+1
    endif
    ip=i
end subroutine hunt2
!*****************************************************************************************
