!*****************************************************************************************
subroutine ann_train(parini)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr, typ_ekf
    use mod_atoms, only: typ_atoms_arr
    use dynamic_memory
    !use ifport
    use mod_processors, only: iproc, mpi_comm_abz
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    type(typ_ann_arr):: ann_arr
    type(typ_ekf):: ekf
    type(typ_atoms_arr):: atoms_train
    type(typ_atoms_arr):: atoms_valid
    type(typ_symfunc_arr):: symfunc_train
    type(typ_symfunc_arr):: symfunc_valid
    integer:: iat, jat, ialpha, i, ig, iconf, ng, ios, nat, ib
    real(8):: time1, time2, time3
    real(8):: tt, epot
    !character(6):: path1
    !character(2):: path2
    character(15):: fnout
    real(8), allocatable:: ratred(:,:)
    call f_routine(id='ann_train')
    !write(path1,'(a3,i3.3)') 'set',iproc
    !path2='..'
    !istat=chdir(path1)
    !if(istat/=0) stop 'ERROR: could not change directory to path1'
    !write(*,*) trim(parini%stypat_ann)
    !call count_words(parini%stypat_ann,ann_arr%n)
    !read(parini%stypat_ann,*) ann_arr%stypat(1:ann_arr%n)
    ann_arr%n=parini%ntypat
    if(parini%bondbased_ann) then
        ann_arr%n=4
    endif
    if(ann_arr%n==0) stop 'ERROR: number of type of atoms zero in ann_train'
    !if(parini%bondbased_ann) then
    !    do i=1,ann_arr%n
    !        ann_arr%ltypat(i)=1
    !    enddo
    !    do i=1,ann_arr%n
    !        ann_arr%stypat(i)=ann_arr%stypat(1)
    !        write(*,*) i,ann_arr%stypat(i)
    !    enddo
    !else
    !    do i=1,ann_arr%n
    !        ann_arr%ltypat(i)=i
    !        write(*,*) i,ann_arr%stypat(i)
    !    enddo
    !endif
    write(*,*) 'Here', ann_arr%n
    allocate(ann_arr%ann(ann_arr%n))
    ann_arr%approach=trim(parini%approach_ann)
    call read_input_ann(parini,iproc,ann_arr)
    !---------------------------------------------
    ekf%num(1:10)=0
    ekf%n=0
    do i=1,ann_arr%n
        do ialpha=1,ann_arr%ann(i)%nl
            ekf%num(i)=ekf%num(i)+(ann_arr%ann(i)%nn(ialpha-1)+1)*ann_arr%ann(i)%nn(ialpha)
        enddo
        ekf%loc(i)=ekf%n+1
        ekf%n=ekf%n+ekf%num(i)
        write(*,'(a,3i5)') 'EKF: ',ekf%loc(i),ekf%num(i),ekf%n
    enddo
    call ann_allocate(ekf,ann_arr)
    call read_data(parini,'list_posinp_train',atoms_train)
    call read_data(parini,'list_posinp_valid',atoms_valid)
    if(iproc==0) then
        write(*,'(a,i)') 'number of ANN wights:             ',ekf%n
        write(*,'(a,i)') 'number of training data points:   ',atoms_train%nconf
        write(*,'(a,i)') 'number of validating data points: ',atoms_valid%nconf
    endif
    do iconf=1,atoms_train%nconf
        if(trim(atoms_train%atoms(iconf)%boundcond)=='bulk') then
        ratred=f_malloc([1.to.3,1.to.atoms_train%atoms(iconf)%nat],id='ratred')
        call rxyz_cart2int_alborz(atoms_train%atoms(iconf)%nat,atoms_train%atoms(iconf)%cellvec,atoms_train%atoms(iconf)%rat,ratred)
        call backtocell_alborz(atoms_train%atoms(iconf)%nat,atoms_train%atoms(iconf)%cellvec,ratred)
        call rxyz_int2cart_alborz(atoms_train%atoms(iconf)%nat,atoms_train%atoms(iconf)%cellvec,ratred,atoms_train%atoms(iconf)%rat)
        call f_free(ratred)
        endif
        do iat=1,atoms_train%atoms(iconf)%nat
            do i=1,ann_arr%n
                if(trim(atoms_train%atoms(iconf)%sat(iat))==trim(parini%stypat(i))) then
                    atoms_train%atoms(iconf)%itypat(iat)=parini%ltypat(i)
                    exit
                endif
            enddo
        enddo
    enddo
    do iconf=1,atoms_valid%nconf
        if(trim(atoms_valid%atoms(iconf)%boundcond)=='bulk') then
        ratred=f_malloc([1.to.3,1.to.atoms_valid%atoms(iconf)%nat],id='ratred')
        call rxyz_cart2int_alborz(atoms_valid%atoms(iconf)%nat,atoms_valid%atoms(iconf)%cellvec,atoms_valid%atoms(iconf)%rat,ratred)
        call backtocell_alborz(atoms_valid%atoms(iconf)%nat,atoms_valid%atoms(iconf)%cellvec,ratred)
        call rxyz_int2cart_alborz(atoms_valid%atoms(iconf)%nat,atoms_valid%atoms(iconf)%cellvec,ratred,atoms_valid%atoms(iconf)%rat)
        call f_free(ratred)
        endif
        do iat=1,atoms_valid%atoms(iconf)%nat
            do i=1,ann_arr%n
                if(trim(atoms_valid%atoms(iconf)%sat(iat))==trim(parini%stypat(i))) then
                    atoms_valid%atoms(iconf)%itypat(iat)=parini%ltypat(i)
                    exit
                endif
            enddo
        enddo
    enddo
    !allocate(atoms_train%inclusion(atoms_train%nconf),source=0)
    if(iproc==0) then
        !write(fnout,'(a12,i3.3)') 'err_train',iproc
        fnout='err_train'
        open(unit=11,file=fnout,status='replace',iostat=ios)
        if(ios/=0) then;write(*,'(a)') 'ERROR: failure openning output file';stop;endif
        !write(fnout,'(a12,i3.3)') 'err_valid',iproc
        fnout='err_valid'
        open(unit=12,file=fnout,status='replace',iostat=ios)
        if(ios/=0) then;write(*,'(a)') 'ERROR: failure openning output file';stop;endif
    endif
    !-------------------------------------------------------
    call cpu_time(time1)
    call set_gbounds(parini,ann_arr,atoms_valid,'bounds_valid',symfunc_valid)
    call cpu_time(time2)
    call set_gbounds(parini,ann_arr,atoms_train,'bounds_train',symfunc_train)
    call cpu_time(time3)
    write(*,'(a,2f10.1)') 'TIMING: evaluation symmetry functions: ',time2-time1,time3-time2
    !do i=1,ann_arr%n !HERE
    !    do ig=1,symfunc_valid%symfunc(1)%ng !There must ng of ann not typ_symfunc
    !        ann_arr%ann(i)%gbounds(1,ig)=ann_arr%ann(1)%gbounds(1,ig)
    !        ann_arr%ann(i)%gbounds(2,ig)=ann_arr%ann(1)%gbounds(2,ig)
    !        ann_arr%ann(i)%two_over_gdiff(ig)=ann_arr%ann(1)%two_over_gdiff(ig)
    !    enddo
    !enddo
    !The following must be done after set_gbounds is called for training set.
    if (parini%bondbased_ann) then
    !-------------------------------- bond symmetry function --------------------------------
        do iconf=1,atoms_valid%nconf
            if(atoms_valid%atoms(iconf)%ntypat>1) stop 'ERROR: this part not ready for ntypat>1'
            if(symfunc_valid%symfunc(iconf)%linked_lists%maxbound_rad/=2) stop 'ERROR: correct next line'
            do ib=1,1
                !iat=symfunc%linked_lists%bound_rad(1,ib)
                !jat=symfunc%linked_lists%bound_rad(2,ib)
                ng=symfunc_valid%symfunc(iconf)%ng
                do ig=1,ng
                    tt=symfunc_valid%symfunc(iconf)%y(ig,ib)
                    tt=(tt-ann_arr%ann(1)%gbounds(1,ig))*ann_arr%ann(1)%two_over_gdiff(ig)-1.d0
                    symfunc_valid%symfunc(iconf)%y(ig,ib)=tt
                enddo
            enddo
        enddo
        do iconf=1,atoms_train%nconf
            if(symfunc_train%symfunc(iconf)%linked_lists%maxbound_rad/=2) stop 'ERROR: correct next line'
            do ib=1,1 !symfunc_train%symfunc(iconf)%linked_lists%maxbound_rad
                !iat=symfunc%linked_lists%bound_rad(1,ib)
                !jat=symfunc%linked_lists%bound_rad(2,ib)
                ng=symfunc_train%symfunc(iconf)%ng
                do ig=1,ng
                    tt=symfunc_train%symfunc(iconf)%y(ig,ib)
                    tt=(tt-ann_arr%ann(1)%gbounds(1,ig))*ann_arr%ann(1)%two_over_gdiff(ig)-1.d0
                    symfunc_train%symfunc(iconf)%y(ig,ib)=tt
                enddo
            enddo
        enddo
    !---------------------------------------------------------------------------------------
    else
    do iconf=1,atoms_valid%nconf
        do iat=1,atoms_valid%atoms(iconf)%nat
            ng=symfunc_valid%symfunc(iconf)%ng
            i=atoms_valid%atoms(iconf)%itypat(iat)
            do ig=1,ng
                tt=symfunc_valid%symfunc(iconf)%y(ig,iat)
                tt=(tt-ann_arr%ann(i)%gbounds(1,ig))*ann_arr%ann(i)%two_over_gdiff(ig)-1.d0
                symfunc_valid%symfunc(iconf)%y(ig,iat)=tt
            enddo
        enddo
    enddo
    do iconf=1,atoms_train%nconf
        do iat=1,atoms_train%atoms(iconf)%nat
            ng=symfunc_train%symfunc(iconf)%ng
            i=atoms_train%atoms(iconf)%itypat(iat)
            do ig=1,ng
                tt=symfunc_train%symfunc(iconf)%y(ig,iat)
                tt=(tt-ann_arr%ann(i)%gbounds(1,ig))*ann_arr%ann(i)%two_over_gdiff(ig)-1.d0
                symfunc_train%symfunc(iconf)%y(ig,iat)=tt
            enddo
        enddo
    enddo
    endif
    ann_arr%ann(1)%ebounds(1)= 1.d20
    ann_arr%ann(1)%ebounds(2)=-1.d20
    do iconf=1,atoms_train%nconf
        tt=atoms_train%atoms(iconf)%epot/atoms_train%atoms(iconf)%nat
        ann_arr%ann(1)%ebounds(1)=min(ann_arr%ann(1)%ebounds(1),tt)
        ann_arr%ann(1)%ebounds(2)=max(ann_arr%ann(1)%ebounds(2),tt)
    enddo
    if(atoms_train%nconf==1) then
        tt=atoms_train%atoms(1)%epot !/atoms_train%atoms(1)%nat
        ann_arr%ann(1)%ebounds(1)=(tt-1.d0)/atoms_train%atoms(1)%nat
        ann_arr%ann(1)%ebounds(2)=(tt+1.d0)/atoms_train%atoms(1)%nat
    endif
    ann_arr%ann(1)%ebounds(1)=-1.d0
    ann_arr%ann(1)%ebounds(2)= 1.d0
    write(*,'(a,2es14.5)') 'ebounds: ',ann_arr%ann(1)%ebounds(1),ann_arr%ann(1)%ebounds(2)
    tt=2.d0/(ann_arr%ann(1)%ebounds(2)-ann_arr%ann(1)%ebounds(1))
    do iconf=1,atoms_train%nconf
        nat=1 !atoms_train%atoms(iconf)%nat
        epot=atoms_train%atoms(iconf)%epot
        symfunc_train%symfunc(iconf)%epot=(epot/nat-ann_arr%ann(1)%ebounds(1))*tt-1.d0
        !write(*,'(a,i6,f8.3)') 'train: ',iconf,symfunc_train%symfunc(iconf)%epot
    enddo
    do iconf=1,atoms_valid%nconf
        nat=1 !atoms_valid%atoms(iconf)%nat
        epot=atoms_valid%atoms(iconf)%epot
        symfunc_valid%symfunc(iconf)%epot=(epot/nat-ann_arr%ann(1)%ebounds(1))*tt-1.d0
        !write(*,'(a,i6,f8.3)') 'valid: ',iconf,symfunc_valid%symfunc(iconf)%epot
    enddo
    !-------------------------------------------------------
    ekf%x=f_malloc([1.to.ekf%n],id='ekf%x')
    call set_annweights(parini,ekf)
    ann_arr%compute_symfunc=.false.
    if(trim(parini%optimizer_ann)=='behler') then
        call ekf_behler(parini,ann_arr,symfunc_train,symfunc_valid,atoms_train,atoms_valid,ekf)
    else if(trim(parini%optimizer_ann)=='rivals') then
        call ekf_rivals(parini,ann_arr,symfunc_train,symfunc_valid,atoms_train,atoms_valid,ekf)
    else if(trim(parini%optimizer_ann)=='rivals_tmp') then
        call ekf_rivals_tmp(parini,ann_arr,symfunc_train,symfunc_valid,atoms_train,atoms_valid,ekf)
    else if(trim(parini%optimizer_ann)=='lm') then
        call ann_lm(parini,ann_arr,atoms_train,atoms_valid,symfunc_train,symfunc_valid,ekf)
    else
        write(*,*) 'ERROR: unknown optimzer in ANN training'
    endif

    !call convert_x_ann(ekf%n,ekf%x,ann_arr) !HERE
    if(iproc==0) then
        do i=1,ann_arr%n
            call write_ann(parini,trim(parini%stypat(i))//'.ann.param',ann_arr%ann(i))
        enddo
    endif
    call f_free(ekf%x)

    if(iproc==0) then
        close(11)
        close(12)
    endif

    call ann_deallocate(ann_arr)

    do iconf=1,atoms_train%nconf
        call f_free(symfunc_train%symfunc(iconf)%linked_lists%prime_bound)
        call f_free(symfunc_train%symfunc(iconf)%linked_lists%bound_rad)
        call f_free(symfunc_train%symfunc(iconf)%linked_lists%bound_ang)
    enddo
    do iconf=1,atoms_valid%nconf
        call f_free(symfunc_valid%symfunc(iconf)%linked_lists%prime_bound)
        call f_free(symfunc_valid%symfunc(iconf)%linked_lists%bound_rad)
        call f_free(symfunc_valid%symfunc(iconf)%linked_lists%bound_ang)
    enddo

    !do iconf=1,atoms_train%nconf
    !    call atom_deallocate(atoms_train%atoms(iconf))
    !enddo
    !do iconf=1,atoms_valid%nconf
    !    call atom_deallocate(atoms_valid%atoms(iconf))
    !enddo


    !istat=chdir(path2)
    !if(istat/=0) stop 'ERROR: could not change directory to path2'
    !deallocate(atoms_train%inclusion)
    call f_release_routine()
end subroutine ann_train
!*****************************************************************************************
subroutine ann_evaluate(parini,iter,ann_arr,symfunc_arr,atoms_arr,ifile,partb)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr, typ_ekf
    use mod_atoms, only: typ_atoms, typ_atoms_arr
    use mod_processors, only: iproc
    use mod_tightbinding, only: typ_partb
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iter
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc_arr), intent(inout):: symfunc_arr
    type(typ_atoms_arr), intent(inout):: atoms_arr
    integer, intent(in):: ifile
    type(typ_partb), optional, intent(inout):: partb
    !local variables
    type(typ_atoms):: atoms
    real(8):: rmse, errmax, tt, pi
    real(8):: frmse, ttx, tty, ttz, ppx, ppy, ppz, tt1, tt2, tt3, ttn, tta, ttmax
    integer:: iconf, ierrmax, iat, nat_tot, nconf_force
    real(8):: time1=0.d0
    real(8):: time2=0.d0
    real(8), save:: time_p=0.d0
    real(8):: dtime1, dtime2
    integer:: ilarge1, ilarge2, ilarge3
    character(28):: frmt1='(i6,5f10.3,i7,i5,3i6,a40,i6)'
    character(28):: frmt2='(i6,5e10.1,i7,i5,3i6,a40,i6)'
    character(28):: frmt
    call cpu_time(time1)
    pi=4.d0*atan(1.d0)
    rmse=0.d0
    frmse=0.d0
    ttn=0.d0
    tta=0.d0
    nat_tot=0
    nconf_force=0
    errmax=0.d0
    ierrmax=0
    ilarge1=0
    ilarge2=0
    ilarge3=0
    ann_arr%event='evalu'
    configuration: do iconf=1,atoms_arr%nconf
        call atom_copy_old(atoms_arr%atoms(iconf),atoms,'atoms_arr%atoms(iconf)->atoms')
        call eval_cal_ann_main(parini,atoms,symfunc_arr%symfunc(iconf),ann_arr)
        if(iter==parini%nstep_ekf) then
            write(40+ifile,'(2i6,2es24.15,es14.5)') iconf,atoms_arr%atoms(iconf)%nat, &
                atoms_arr%atoms(iconf)%epot/atoms_arr%atoms(iconf)%nat,atoms%epot/atoms_arr%atoms(iconf)%nat, &
                (atoms%epot-atoms_arr%atoms(iconf)%epot)/atoms_arr%atoms(iconf)%nat
        endif
        tt=abs(atoms%epot-atoms_arr%atoms(iconf)%epot)/atoms_arr%atoms(iconf)%nat
        if(tt>1.d-2) ilarge1=ilarge1+1
        if(tt>1.d-3) ilarge2=ilarge2+1
        if(tt>1.d-4) ilarge3=ilarge3+1
        if(tt>errmax) then
            errmax=tt
            ierrmax=iconf
        endif
        rmse=rmse+tt**2
        !if(tt>1.d-2) then
        !    atoms_arr%inclusion(iconf)=0
        !else
        !    atoms_arr%inclusion(iconf)=1
        !endif
        !write(22,'(a,i5.5)') 'configuration ',iconf
        if(atoms%nat<=parini%nat_force) then
        nat_tot=nat_tot+atoms%nat
        nconf_force=nconf_force+1
        do iat=1,atoms%nat
            ttx=atoms_arr%atoms(iconf)%fat(1,iat)
            tty=atoms_arr%atoms(iconf)%fat(2,iat)
            ttz=atoms_arr%atoms(iconf)%fat(3,iat)
            ppx=atoms%fat(1,iat)
            ppy=atoms%fat(2,iat)
            ppz=atoms%fat(3,iat)
            tt1=sqrt(ttx**2+tty**2+ttz**2)
            tt2=sqrt(ppx**2+ppy**2+ppz**2)
            tt3=(ppx-ttx)**2+(ppy-tty)**2+(ppz-ttz)**2
            frmse=frmse+tt3
            tt3=sqrt(tt3)
            !write(25+ifile,'(i4,3es14.5)') iat,tt1,tt2,tt3
            !write(22,'(3es19.10)') ttx,tty,ttz
        enddo
        ttn=ttn+ann_arr%fchi_norm
        tta=tta+ann_arr%fchi_angle
        !write(44,'(2i7,4es14.5)') iter,iconf,ann_arr%fchi_norm,ann_arr%fchi_angle,ttn/nconf_force,tta/nconf_force
        endif
    enddo configuration
    !stop 'HERe'
    rmse=sqrt(rmse/real(atoms_arr%nconf,8))
    if(nconf_force==0) nconf_force=1
    ttn=ttn/real(nconf_force,8)
    tta=tta/real(nconf_force,8)
    if(nat_tot==0) nat_tot=1
    frmse=sqrt(frmse/real(3*nat_tot,8))
    if(iproc==0) then
        rmse=rmse*1.d3
        errmax=errmax*1.d3
        if(rmse>99999.d0) then
            frmt=frmt2
        else
            frmt=frmt1
        endif
        write(ifile,frmt) iter,rmse,ttn,tta,frmse,errmax,ierrmax,atoms_arr%atoms(ierrmax)%nat, &
            ilarge1,ilarge2,ilarge3,trim(atoms_arr%fn(ierrmax)),atoms_arr%lconf(ierrmax)
        !if(ifile==11) then
        !tt1=1.d0/(ann_arr%ann(1)%gausswidth*sqrt(2.d0*pi))+ann_arr%ann(1)%hardness
        !tt2=1.d0/(ann_arr%ann(2)%gausswidth*sqrt(2.d0*pi))+ann_arr%ann(2)%hardness
        !ttmax=max(tt1,tt2)
        !write(*,'(a,i6,5f10.3,i7,3f10.3)') 'ETA ',iter,rmse,ttn,tta,frmse,errmax,ierrmax,tt1,tt2,ttmax
        !endif
    endif
    call cpu_time(time2)
    dtime1=time1-time_p
    dtime2=time2-time1
    write(*,'(a,2f20.2)') 'TIME ',dtime1,dtime2
    time_p=time2
end subroutine ann_evaluate
!*****************************************************************************************
!Subroutine cal_ann_main is called during training process.
!It does the same job as subroutine eval_cal_ann_main with a few more
!commands related to ANN weights.
!These two subroutine are supposed to be merged.
!eval_cal_ann_main is called for task_training but only for
!obtaining energy/forces when information on errors on
!training and validation data points are expected.
!eval_cal_ann_main is called by potential_ANN
subroutine eval_cal_ann_main(parini,atoms,symfunc,ann_arr)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_tightbinding, only: typ_partb
    use mod_ann, only: typ_ann_arr, typ_symfunc, typ_ekf
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc), intent(inout):: symfunc
    !local variables
    type(typ_partb):: partb
    type(typ_ekf):: ekf
    if(trim(ann_arr%event)=='potential') then
        !allocate(symfunc%y(ann_arr%ann(1)%nn(0),atoms%nat))
    endif
    if(trim(ann_arr%approach)=='atombased') then
        call cal_ann_atombased(parini,atoms,symfunc,ann_arr,ekf)
    elseif(trim(ann_arr%approach)=='eem1') then
        call cal_ann_eem1(parini,atoms,symfunc,ann_arr,ekf)
    elseif(trim(ann_arr%approach)=='eem2') then
        call cal_ann_eem2(parini,atoms,symfunc,ann_arr,ekf)
    elseif(trim(ann_arr%approach)=='tb') then
        call cal_ann_tb(parini,partb,atoms,ann_arr,symfunc,ekf)
    else
        write(*,'(2a)') 'ERROR: unknown approach in ANN, ',trim(ann_arr%approach)
        stop
    endif
    if(trim(ann_arr%event)=='potential') then
        !deallocate(symfunc%y)
        !call ann_deallocate(ann_arr)
    endif
end subroutine eval_cal_ann_main
!*****************************************************************************************
subroutine set_gbounds(parini,ann_arr,atoms_arr,strmess,symfunc_arr)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr
    use mod_atoms, only: typ_atoms_arr
    use mod_linked_lists, only: typ_pia_arr
    use mod_processors, only: iproc, nproc, mpi_comm_abz
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms_arr), intent(inout):: atoms_arr
    character(*), intent(in):: strmess
    type(typ_symfunc_arr), intent(inout):: symfunc_arr
    !local variables 
    character(5):: sat
    real(8), allocatable:: gminarr(:,:), gmaxarr(:,:) !, poll_period
    real(8):: eps=epsilon(1.d0)
    real(8):: ttx, tty, ttz
    real(8), allocatable:: wa(:)
    integer:: i, ig, iconf, iat, jat, i0, ios, nat_t, ng_t, nb_t, n, ib
    integer, allocatable:: iatmin(:,:), iatmax(:,:), iconfmin(:,:), iconfmax(:,:)
    integer:: ibmin(100), ibmax(100)
    integer:: jproc, ierr, ireq, ireq_tmp, nreq, mreq, ii, nwa !, itry
    type(typ_pia_arr):: pia_arr_tmp
    integer, allocatable:: ireqarr(:)
    logical:: flag
    character(30):: filename
    character(100):: smsg
    integer:: ngmax
#if defined(MPI)
    include 'mpif.h'
    integer:: status_mpi(MPI_STATUS_SIZE)
#endif
    call f_routine(id='set_gbounds')
#if defined(MPI)
    associate(MPI_DP=>MPI_DOUBLE_PRECISION)
    nreq=atoms_arr%nconf/nproc
    if(mod(atoms_arr%nconf,nproc)>iproc) nreq=nreq+1
    nreq=atoms_arr%nconf-nreq
    ireqarr=f_malloc([1.to.nreq],id='ireqarr')
#endif
    ngmax=200
    gminarr=f_malloc([1.to.ngmax,1.to.parini%ntypat],id='gminarr')
    gminarr=huge(1.d20)
    gmaxarr=f_malloc([1.to.ngmax,1.to.parini%ntypat],id='gmaxarr')
    gmaxarr=-huge(1.d20)
    iatmin=f_malloc0([1.to.ngmax,1.to.parini%ntypat],id='iatmin')
    iatmax=f_malloc0([1.to.ngmax,1.to.parini%ntypat],id='iatmax')
    iconfmin=f_malloc0([1.to.ngmax,1.to.parini%ntypat],id='iconfmin')
    iconfmax=f_malloc0([1.to.ngmax,1.to.parini%ntypat],id='iconfmax')
    ibmin(1:100)=0 ; ibmax(1:100)=0
    if(.not. allocated(symfunc_arr%symfunc)) then
        symfunc_arr%nconf=atoms_arr%nconf
        allocate(symfunc_arr%symfunc(symfunc_arr%nconf))
    endif
    !write(*,'(a,i3,i6)') 'iproc,nconf ',iproc,atoms_arr%nconf
    do iconf=1,atoms_arr%nconf
        symfunc_arr%symfunc(iconf)%ng=ann_arr%ann(1)%nn(0) !HERE
        symfunc_arr%symfunc(iconf)%nat=atoms_arr%atoms(iconf)%nat
    enddo
    configuration: do iconf=1+iproc,atoms_arr%nconf,nproc
        if(trim(parini%symfunc)=='write' .or. trim(parini%symfunc)=='only_calculate') then
            call symmetry_functions(parini,ann_arr,atoms_arr%atoms(iconf),symfunc_arr%symfunc(iconf),.false.)
            deallocate(symfunc_arr%symfunc(iconf)%y0dr)
            if(parini%bondbased_ann .and. symfunc_arr%symfunc(iconf)%linked_lists%maxbound_rad/=2) then
                stop 'ERROR: correct next line'
            endif
            !call ann_deallocate(ann_arr)
            if(parini%symfunc_type_ann=='behler') then
            !call f_free(symfunc_arr%symfunc(iconf)%linked_lists%prime_bound)
            if(.not. parini%bondbased_ann) then
            !call f_free(symfunc_arr%symfunc(iconf)%linked_lists%bound_rad)
            !call f_free(symfunc_arr%symfunc(iconf)%linked_lists%bound_ang)
            endif
            endif
        elseif(trim(parini%symfunc)/='read') then
            stop 'ERROR: arini%symfunc contains none of the three acceptable possibilies'
        endif
        if(trim(parini%symfunc)=='write') then
            !Symmetry functions are written into files to be used for
            !subsequent training runs.
            if(trim(strmess)=='bounds_train') then
                write(filename,'(a25,i5.5)') '../symfunc/train.symfunc.',iconf
            elseif(trim(strmess)=='bounds_valid') then
                write(filename,'(a25,i5.5)') '../symfunc/valid.symfunc.',iconf
            else
                stop 'ERROR: invalid content in strmess in gset_bounds '
            endif
            open(unit=311,file=trim(filename),status='replace',form='unformatted', &
                access='stream',iostat=ios)
            if(ios/=0) then
                write(*,'(2a)') 'ERROR: failure openning file ',trim(filename)
                stop
            endif
            associate(nb=>symfunc_arr%symfunc(iconf)%linked_lists%maxbound_rad)
            associate(ng=>symfunc_arr%symfunc(iconf)%ng)
            associate(nat=>atoms_arr%atoms(iconf)%nat)
            if(nat<=parini%nat_force) then
                nwa=3+nat*(3+ng)+ng*3*nb
            else
                nwa=3+nat*(3+ng)
            endif
            wa=f_malloc([1.to.nwa],id='wa')
            wa(1)=real(nat,8)
            wa(2)=real(ng,8)
            wa(3)=real(nb,8)
            do iat=1,atoms_arr%atoms(iconf)%nat
                wa(3+iat*3-2)=atoms_arr%atoms(iconf)%rat(1,iat)
                wa(3+iat*3-1)=atoms_arr%atoms(iconf)%rat(2,iat)
                wa(3+iat*3-0)=atoms_arr%atoms(iconf)%rat(3,iat)
            enddo
            n=3+3*atoms_arr%atoms(iconf)%nat
            bondbased: if(parini%bondbased_ann)then
            !-------------------------- bond symmetry functions ---------------------------
            write(*,*) 'ERROR: writing/reading symmetry function values from files not'
            write(*,*) 'working yet, for several reasons for example allocation of wa'
            stop
            !----------------------------------------------------------------------------
            else
            do iat=1,atoms_arr%atoms(iconf)%nat
                do ig=1,symfunc_arr%symfunc(iconf)%ng
                    n=n+1
                    wa(n)=symfunc_arr%symfunc(iconf)%y(ig,iat)
                enddo
            enddo
            if(nat<=parini%nat_force) then
            do ib=1,nb
                do i=1,3
                    do ig=1,symfunc_arr%symfunc(iconf)%ng
                        n=n+1
                        wa(n)=symfunc_arr%symfunc(iconf)%y0d(ig,i,ib)
                    enddo
                enddo
            enddo
            else
                deallocate(symfunc_arr%symfunc(iconf)%y0d)
            endif
            endif bondbased
            end associate
            end associate
            end associate
            
            write(311,iostat=ios) wa
            if(ios/=0) then
                write(*,'(2a)') 'ERROR: failure writing to file ',trim(filename)
                stop
            endif
            !write(311,'(2i6)') atoms_arr%atoms(iconf)%nat,symfunc_arr%symfunc(iconf)%ng
            !do iat=1,atoms_arr%atoms(iconf)%nat
            !    do ig=1,symfunc_arr%symfunc(iconf)%ng
            !        write(311,'(es25.16)') symfunc_arr%symfunc(iconf)%y(ig,iat)
            !    enddo
            !enddo
            close(311)
            call f_free(wa)
        elseif(trim(parini%symfunc)=='read') then
            !Symmetry functions which are previously calculated and written by
            !some other run is going to be read from files
            if(trim(strmess)=='bounds_train') then
                write(filename,'(a25,i5.5)') '../symfunc/train.symfunc.',iconf
            elseif(trim(strmess)=='bounds_valid') then
                write(filename,'(a25,i5.5)') '../symfunc/valid.symfunc.',iconf
            else
                stop 'ERROR: invalid content in strmess in gset_bounds '
            endif
            !open(unit=311,file=trim(filename),status='old',iostat=ios)
            open(unit=311,file=trim(filename),status='old',form='unformatted', &
                access='stream',iostat=ios)
            if(ios/=0) then
                write(*,'(2a)') 'ERROR: failure openning file ',trim(filename)
                stop
            endif
            !read(311,*) nat_t,ng_t
            associate(nb=>symfunc_arr%symfunc(iconf)%linked_lists%maxbound_rad)
            associate(ng=>symfunc_arr%symfunc(iconf)%ng)
            associate(nat=>atoms_arr%atoms(iconf)%nat)
            symfunc_arr%symfunc(iconf)%linked_lists%rcut=ann_arr%rcut
            symfunc_arr%symfunc(iconf)%linked_lists%triplex=.true.
            call call_linkedlist(parini,atoms_arr%atoms(iconf),symfunc_arr%symfunc(iconf)%linked_lists,pia_arr_tmp)
            deallocate(pia_arr_tmp%pia)
            allocate(symfunc_arr%symfunc(iconf)%y(ng,nat))
            if(nat<=parini%nat_force) then
                allocate(symfunc_arr%symfunc(iconf)%y0d(ng,3,nb))
                nwa=3+nat*(3+ng)+ng*3*nb
            else
                nwa=3+nat*(3+ng)
            endif
            wa=f_malloc([1.to.nwa],id='wa')
            read(311,iostat=ios) wa
            if(ios/=0) then
                write(*,'(2a)') 'ERROR: failure reading from file ',trim(filename)
                stop
            endif
            nat_t=nint(wa(1))
            ng_t=nint(wa(2))
            nb_t=nint(wa(3))
            if(nat_t/=nat .or. ng_t/=ng .or. nb_t/=nb) then
                write(*,'(a,7i6)') 'ERROR: inconsistent nat or ng ',iconf,nat_t,nat,ng_t,ng,nb_t,nb
                stop
            endif
            do iat=1,nat
                ttx=abs(wa(3+iat*3-2)-atoms_arr%atoms(iconf)%rat(1,iat))
                tty=abs(wa(3+iat*3-1)-atoms_arr%atoms(iconf)%rat(2,iat))
                ttz=abs(wa(3+iat*3-0)-atoms_arr%atoms(iconf)%rat(3,iat))
                if(ttx>eps .or. tty>eps .or. ttz>eps) then
                    smsg='ERROR: inconsistency of configuration in symmetry functions file. '
                    write(*,'(a,3es14.5,i7,a,i5)') trim(smsg),ttx,tty,ttz,iconf,trim(atoms_arr%fn(iconf)),atoms_arr%lconf(iconf)
                    stop
                endif
            enddo
            n=3+3*atoms_arr%atoms(iconf)%nat
            bondbased: if(parini%bondbased_ann)then
            !--------------------- bond symmetry functions -------------------------------
            write(*,*) 'ERROR: writing/reading symmetry function values from files not'
            write(*,*) 'working yet, for several reasons for example allocation of wa'
            stop
            !do iat=1,atoms_arr%atoms(iconf)%nat
            !do jat=1,atoms_arr%atoms(iconf)%nat
            !    do ig=1,symfunc_arr%symfunc(iconf)%ng
            !        n=n+1
            !        symfunc_arr%symfunc(iconf)%y_bond(ig,iat,jat)=wa(n)
            !    enddo
            !enddo
            !enddo
            !-----------------------------------------------------------------------------
            else
            do iat=1,nat
                do ig=1,ng
                    n=n+1
                    symfunc_arr%symfunc(iconf)%y(ig,iat)=wa(n)
                enddo
            enddo
            if(nat<=parini%nat_force) then
            do ib=1,nb
                do i=1,3
                    do ig=1,ng
                        n=n+1
                        symfunc_arr%symfunc(iconf)%y0d(ig,i,ib)=wa(n)
                    enddo
                enddo
            enddo
            endif
            endif bondbased            
            end associate
            end associate
            end associate
            !do iat=1,atoms_arr%atoms(iconf)%nat
            !    do ig=1,symfunc_arr%symfunc(iconf)%ng
            !        read(311,*) symfunc_arr%symfunc(iconf)%y(ig,iat)
            !    enddo
            !enddo
            close(311)
            call f_free(wa)
            endif
#if defined(MPI)
        do jproc=0,nproc-1
            if(jproc==iproc) cycle
            associate(ntot=>symfunc_arr%symfunc(iconf)%ng*symfunc_arr%symfunc(iconf)%nat)
            call MPI_ISEND(symfunc_arr%symfunc(iconf)%y(1,1),ntot,MPI_DP,jproc,iconf,mpi_comm_abz,ireq_tmp,ierr)
            end associate
        enddo
#endif
    enddo configuration
#if defined(MPI)
    if(nproc>1) then
    ireq=0
    do iconf=1,atoms_arr%nconf
        jproc=mod(iconf-1,nproc)
        if(jproc==iproc) cycle
        ireq=ireq+1
        associate(ntot=>symfunc_arr%symfunc(iconf)%ng*symfunc_arr%symfunc(iconf)%nat)
        call MPI_IRECV(symfunc_arr%symfunc(iconf)%y(1,1),ntot,MPI_DP,jproc,iconf,mpi_comm_abz,ireqarr(ireq),ierr)
        end associate
    enddo
    call MPI_BARRIER(mpi_comm_abz,ierr)

    !itry=0
    ireq=0
    mreq=nreq
    do
        !itry=itry+1
        ireq=ireq+1
        flag=.false.
        call MPI_TEST(ireqarr(ireq),flag,status_mpi,ierr)
        if(flag) then
            !write(41,'(i3,i6,l2,4i6)') status_mpi(MPI_SOURCE),status_mpi(MPI_TAG),flag,mreq,nreq,ireq,itry
            ii=ireqarr(ireq)
            ireqarr(ireq)=ireqarr(mreq)
            ireqarr(mreq)=ii
            mreq=mreq-1
        endif
        if(mreq==0) exit
        if(ireq>=mreq) ireq=0
    enddo
    call MPI_BARRIER(mpi_comm_abz,ierr)
    endif
#endif
    do iconf=1,atoms_arr%nconf
        !if(mod(iconf-1,nproc)==iproc) cycle
        !write(41,'(i6,i3)',advance='no') mod(iconf-1,nproc),iproc
        !associate(n=>symfunc_arr%symfunc(iconf)%ng*symfunc_arr%symfunc(iconf)%nat)
        !call MPI_RECV(symfunc_arr%symfunc(iconf)%y(1,1),n,MPI_DP,mod(iconf-1,nproc),iconf,mpi_comm_abz,status_mpi,ierr)
        !end associate
        !write(41,'(i6,i3)') status_mpi(MPI_TAG),status_mpi(MPI_SOURCE)
        if(parini%bondbased_ann)then
        !-------------------------------- bond symmetry function --------------------------------
        if(symfunc_arr%symfunc(iconf)%linked_lists%maxbound_rad/=2) stop 'ERROR: correct next line'
        do ib=1,1 !symfunc%linked_lists%maxbound_rad
            do ig=1,symfunc_arr%symfunc(iconf)%ng
                if(symfunc_arr%symfunc(iconf)%y(ig,ib)<gminarr(ig,1)) then
                    gminarr(ig,1)=symfunc_arr%symfunc(iconf)%y(ig,ib)
                    ibmin(ig)=ib
                    iconfmin(ig,1)=iconf
                endif
                if(symfunc_arr%symfunc(iconf)%y(ig,ib)>gmaxarr(ig,1)) then
                    gmaxarr(ig,1)=symfunc_arr%symfunc(iconf)%y(ig,ib)
                    ibmax(ig)=ib
                    iconfmax(ig,1)=iconf
                endif
            enddo
        enddo
        !-----------------------------------------------------------------------------------------
        else
        do iat=1,atoms_arr%atoms(iconf)%nat
            i=atoms_arr%atoms(iconf)%itypat(iat)
            do ig=1,symfunc_arr%symfunc(iconf)%ng
                if(symfunc_arr%symfunc(iconf)%y(ig,iat)<gminarr(ig,i)) then
                    gminarr(ig,i)=symfunc_arr%symfunc(iconf)%y(ig,iat)
                    iatmin(ig,i)=iat
                    iconfmin(ig,i)=iconf
                endif
                if(symfunc_arr%symfunc(iconf)%y(ig,iat)>gmaxarr(ig,i)) then
                    gmaxarr(ig,i)=symfunc_arr%symfunc(iconf)%y(ig,iat)
                    iatmax(ig,i)=iat
                    iconfmax(ig,i)=iconf
                endif
            enddo
        enddo
        endif
    enddo
    do i=1,ann_arr%n
    do ig=1,ann_arr%ann(1)%nn(0) !HERE
        if(iproc==0) then
        if(parini%bondbased_ann)then
            write(*,'(2(i7,2i4,es20.10),1x,a)') &
                iconfmin(ig,1),atoms_arr%atoms(iconfmin(ig,1))%nat,ibmin(ig),gminarr(ig,1), &
                iconfmax(ig,1),atoms_arr%atoms(iconfmax(ig,1))%nat,ibmax(ig),gmaxarr(ig,1),trim(strmess)
            write(*,'(2(a50,i6,1x))') trim(atoms_arr%fn(iconfmin(ig,1))),atoms_arr%lconf(iconfmin(ig,1)), &
                trim(atoms_arr%fn(iconfmax(ig,1))),atoms_arr%lconf(iconfmax(ig,1))
        else
            write(*,'(2(i7,2i4,es20.10),1x,a)') &
                iconfmin(ig,i),atoms_arr%atoms(iconfmin(ig,i))%nat,iatmin(ig,i),gminarr(ig,i), &
                iconfmax(ig,i),atoms_arr%atoms(iconfmax(ig,i))%nat,iatmax(ig,i),gmaxarr(ig,i),trim(strmess)
            write(*,'(2(a50,i6,1x))') trim(atoms_arr%fn(iconfmin(ig,i))),atoms_arr%lconf(iconfmin(ig,i)), &
                trim(atoms_arr%fn(iconfmax(ig,i))),atoms_arr%lconf(iconfmax(ig,i))
        endif
        endif
    enddo
    enddo

    if(parini%bondbased_ann)then
    !------------------------------ bond symmetry function ----------------------------------------
    !IMPORTANT: the upper bound of the loop over i needs to be corrected for multicomponent systems.
    if(atoms_arr%atoms(1)%ntypat>1) stop 'ERROR: this part not ready for ntypat>1'
    do i=1,1
    do i0=1,ann_arr%ann(i)%nn(0)
        if(gminarr(i0,1)==0.d0) then
            ann_arr%ann(i)%gbounds(1,i0)=-epsilon(1.d0)
        else
            ann_arr%ann(i)%gbounds(1,i0)=gminarr(i0,1)
        endif
        if(gmaxarr(i0,1)==0.d0) then
            ann_arr%ann(i)%gbounds(2,i0)=epsilon(1.d0)
        else
            ann_arr%ann(i)%gbounds(2,i0)=gmaxarr(i0,1)
        endif
        !write(*,*) ann_arr%ann(i)%gbounds(2,i0),ann_arr%ann(i)%gbounds(1,i0)
        ann_arr%ann(i)%two_over_gdiff(i0)=2.d0/(ann_arr%ann(i)%gbounds(2,i0)-ann_arr%ann(i)%gbounds(1,i0))
    enddo
    enddo
    ! ----------------------------------------------------------------------------------------------
    else
    do i=1,ann_arr%n
    do i0=1,ann_arr%ann(i)%nn(0)
        !if(abs(gminarr(i0))<epsilon(1.d0)
        if(gminarr(i0,i)==0.d0) then
            ann_arr%ann(i)%gbounds(1,i0)=-epsilon(1.d0)
        else
            ann_arr%ann(i)%gbounds(1,i0)=gminarr(i0,i)
        endif
        if(gmaxarr(i0,i)==0.d0) then
            ann_arr%ann(i)%gbounds(2,i0)=epsilon(1.d0)
        else
            ann_arr%ann(i)%gbounds(2,i0)=gmaxarr(i0,i)
        endif
        ann_arr%ann(i)%two_over_gdiff(i0)=2.d0/(ann_arr%ann(i)%gbounds(2,i0)-ann_arr%ann(i)%gbounds(1,i0))
    enddo
    enddo
    endif
#if defined(MPI)
    call f_free(ireqarr)
    end associate
#endif
    call f_free(gminarr)
    call f_free(gmaxarr)
    call f_free(iatmin)
    call f_free(iatmax)
    call f_free(iconfmin)
    call f_free(iconfmax)
    call f_release_routine()
end subroutine set_gbounds
!*****************************************************************************************
subroutine convert_x_ann(n,x,ann)
    use mod_interface
    use mod_ann, only: typ_ann
    implicit none
    integer, intent(in):: n
    real(8), intent(in):: x(n)
    type(typ_ann), intent(inout):: ann
    !local variables
    integer:: i, j, l, ialpha
    l=0
    do ialpha=1,ann%nl
        do j=1,ann%nn(ialpha)
            do i=1,ann%nn(ialpha-1)
                l=l+1
                ann%a(i,j,ialpha)=x(l)
            enddo
        enddo
        do i=1,ann%nn(ialpha)
            l=l+1
            ann%b(i,ialpha)=x(l)
        enddo
    enddo
    if(l/=n) stop 'ERROR: l/=n'
end subroutine convert_x_ann
!*****************************************************************************************
subroutine convert_ann_epotd(ann,n,epotd)
    use mod_interface
    use mod_ann, only: typ_ann
    implicit none
    type(typ_ann), intent(in):: ann
    integer, intent(in):: n
    real(8), intent(inout):: epotd(n)
    !local variables
    integer:: i, ij, l, ialpha
    l=0
    do ialpha=1,ann%nl
        do ij=1,ann%nn(ialpha)*ann%nn(ialpha-1)
            l=l+1
            epotd(l)=ann%ad(ij,ialpha)
        enddo
        do i=1,ann%nn(ialpha)
            l=l+1
            epotd(l)=ann%bd(i,ialpha)
        enddo
    enddo
    if(l/=n) stop 'ERROR: l/=n'
end subroutine convert_ann_epotd
!*****************************************************************************************
subroutine randomize_data_order(atoms_arr)
    use mod_interface
    use mod_atoms, only: typ_atoms_arr, typ_atoms
    implicit none
    type(typ_atoms_arr), intent(inout):: atoms_arr
    !local variables 
    integer:: i1, i2, iconf !, nat_t
    real(8):: t1, t2 !, rat_t(3,1000), epot_t
    type(typ_atoms):: atoms_t
    !if(natmax>1000) stop 'ERROR: natmax>1000'
    iconf=0
    do
        call random_number(t1)
        call random_number(t2)
        i1=int(t1*real(atoms_arr%nconf,8))+1
        i2=int(t2*real(atoms_arr%nconf,8))+1
        if(i1<1 .or. i1>atoms_arr%nconf) stop 'ERROR: something wrong with i1'
        if(i2<1 .or. i2>atoms_arr%nconf) stop 'ERROR: something wrong with i2'
        if(i1==i2) cycle
        iconf=iconf+1
        call atom_copy_old(atoms_arr%atoms(i1),atoms_t,'atoms_arr%atoms(i1)->atoms_t')
        call atom_copy_old(atoms_arr%atoms(i2),atoms_arr%atoms(i1),'atoms_arr%atoms(i2)->atoms_arr%atoms(i1)')
        call atom_copy_old(atoms_t,atoms_arr%atoms(i2),'atoms_t->atoms_arr%atoms(i2)')

        !nat_t=natarr(i1)
        !rat_t(1:3,1:nat_t)=ratall(1:3,1:nat_t,i1)
        !epot_t=epotall(i1)

        !natarr(i1)=natarr(i2)
        !ratall(1:3,1:natarr(i2),i1)=ratall(1:3,1:natarr(i2),i2)
        !epotall(i1)=epotall(i2)

        !natarr(i2)=nat_t
        !ratall(1:3,1:nat_t,i2)=rat_t(1:3,1:nat_t)
        !epotall(i2)=epot_t

        if(iconf>1*atoms_arr%nconf) exit
    enddo
    call atom_deallocate_old(atoms_t)
end subroutine randomize_data_order
!*****************************************************************************************
