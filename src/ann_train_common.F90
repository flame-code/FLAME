!*****************************************************************************************
subroutine ann_evaluate(parini,iter,ann_arr,symfunc_arr,atoms_arr,data_set)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc
    use mod_symfunc, only: typ_symfunc_arr
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_copy_old, atom_deallocate_old
    use mod_opt_ann, only: typ_opt_ann
    use mod_processors, only: iproc
    use futile
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iter
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc_arr), intent(inout):: symfunc_arr
    type(typ_atoms_arr), intent(inout):: atoms_arr
    character(*), intent(in):: data_set
    !local variables
    type(typ_atoms):: atoms
    type(typ_symfunc):: symfunc
    type(typ_opt_ann):: opt_ann
    real(8):: dpmrmse, rmse, errmax, tt, pi
    real(8):: frmse, ttx, tty, ttz, ppx, ppy, ppz, tt1, tt2, tt3, ttn, tta
    real(8):: dpx,dpy,dpz
    real(8):: dpx_ref,dpy_ref,dpz_ref
    integer:: iconf, ierrmax, iat, nat_tot, nconf_force
    real(8):: time1=0.d0
    real(8):: time2=0.d0
    real(8), save:: time_p=0.d0
    real(8):: dtime1, dtime2
    integer:: ilarge1, ilarge2, ilarge3, iunit, ios
    !character(28):: frmt1='(i6,5f10.3,i7,i5,3i6,a40,i6)'
    !character(28):: frmt2='(i6,5e10.1,i7,i5,3i6,a40,i6)'
    !character(28):: frmt
    character(28):: fmt_main
    character(28):: fmt1='(f10.3)'
    character(28):: fmt2='(e10.1)'
    character(15):: filename
    character(len=8):: str_key
    call cpu_time(time1)
    pi=4.d0*atan(1.d0)
    rmse=0.d0
    dpmrmse=0.d0
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
    if(parini%save_symfunc_behnam) then
        ann_arr%compute_symfunc=.false.
    else
        ann_arr%compute_symfunc=.true.
    endif
    if(parini%print_energy) then
        write(filename,'(a12,i3.3)') 'detailed_err',iter
        iunit=f_get_free_unit(10**5)
        if(trim(data_set)=="train") then
            open(unit=iunit,file=trim(filename),status='unknown',iostat=ios)
        elseif(trim(data_set)=="valid") then
            open(unit=iunit,file=trim(filename),status='old',access='append',iostat=ios)
        endif
        if(ios/=0) then
            write(*,'(a,a)') 'ERROR: failure openning file: ',trim(filename)
            stop
        endif
        !write(iunit,'(a2,a44,4a23)') "#", " ","E_dft","E_ann","E_dft-E_ann/atom (Ha)","E_dft-E_ann (eV)"
    endif
    configuration: do iconf=1,atoms_arr%nconf
        if(.not. atoms_arr%conf_inc(iconf)) cycle
        call atom_copy_old(atoms_arr%atoms(iconf),atoms,'atoms_arr%atoms(iconf)->atoms')
        if(parini%save_symfunc_behnam) then
            call cal_ann_main(parini,atoms,symfunc_arr%symfunc(iconf),ann_arr,opt_ann)
        else
            call cal_ann_main(parini,atoms,symfunc,ann_arr,opt_ann)
        endif
        !if(iter==parini%nstep_opt_ann) then
        !    write(40+ifile,'(2i6,2es24.15,es14.5)') iconf,atoms_arr%atoms(iconf)%nat, &
        !        atoms_arr%atoms(iconf)%epot/atoms_arr%atoms(iconf)%nat,atoms%epot/atoms_arr%atoms(iconf)%nat, &
        !        (atoms%epot-atoms_arr%atoms(iconf)%epot)/atoms_arr%atoms(iconf)%nat
        !endif
        tt=abs(atoms%epot-atoms_arr%atoms(iconf)%epot)/atoms_arr%atoms(iconf)%nat
        tt1=atoms%epot/atoms_arr%atoms(iconf)%nat
        tt2=atoms_arr%atoms(iconf)%epot/atoms_arr%atoms(iconf)%nat
        !HERE
        if(parini%print_energy) then
          !write(iunit,'(i7,3es14.5,a35,i6,a)') iconf,tt,tt1,tt2,atoms_arr%lconf(iconf),trim(data_set)
          write(iunit,'(i7,es16.8,es16.8,a40,i6,a)') iconf,tt1,tt2,trim(atoms_arr%fn(iconf)),atoms_arr%lconf(iconf),trim(data_set)
        endif
        if(tt>1.d-2) ilarge1=ilarge1+1
        if(tt>1.d-3) ilarge2=ilarge2+1
        if(tt>1.d-4) ilarge3=ilarge3+1
        if(tt>errmax) then
            errmax=tt
            ierrmax=iconf
        endif
        rmse=rmse+tt**2
        dpx_ref=atoms_arr%atoms(iconf)%dpm(1)
        dpy_ref=atoms_arr%atoms(iconf)%dpm(2)
        dpz_ref=atoms_arr%atoms(iconf)%dpm(3)
        dpx=atoms%dpm(1)
        dpy=atoms%dpm(2)
        dpz=atoms%dpm(3)
        dpmrmse=dpmrmse+(dpx_ref-dpx)**2+(dpy_ref-dpy)**2+(dpz_ref-dpz)**2
        write(1391,'(3es18.8,a3,3es18.8,a3,es18.8)')dpx,dpy,dpz,' | ',dpx_ref,dpy_ref,dpz_ref,' | ',dpmrmse
        !write(22,'(a,i5.5)') 'configuration ',iconf
        !if(atoms%nat<=parini%nat_force) then
        nat_tot=nat_tot+atoms%nat
        nconf_force=nconf_force+1
        do iat=1,atoms%nat
            ttx=atoms_arr%atoms(iconf)%fat(1,iat)
            tty=atoms_arr%atoms(iconf)%fat(2,iat)
            ttz=atoms_arr%atoms(iconf)%fat(3,iat)
            ppx=atoms%fat(1,iat)
            ppy=atoms%fat(2,iat)
            ppz=atoms%fat(3,iat)
            !write(41,'(2i6,6f7.3)') iconf,iat,ttx,tty,ttz,ppx,ppy,ppz
            tt1=sqrt(ttx**2+tty**2+ttz**2)
            tt2=sqrt(ppx**2+ppy**2+ppz**2)
            tt3=(ppx-ttx)**2+(ppy-tty)**2+(ppz-ttz)**2
            frmse=frmse+tt3
            tt3=sqrt(tt3)
        enddo
        ttn=ttn+ann_arr%fchi_norm
        tta=tta+ann_arr%fchi_angle
        call atom_deallocate_old(atoms)
        !write(44,'(2i7,4es14.5)') iter,iconf,ann_arr%fchi_norm,ann_arr%fchi_angle,ttn/nconf_force,tta/nconf_force
        !endif
    enddo configuration
    rmse=sqrt(rmse/real(atoms_arr%nconf_inc,8))
    dpmrmse=sqrt(dpmrmse/real(3*atoms_arr%nconf_inc,8))
    if(nconf_force==0) nconf_force=1
    ttn=ttn/real(nconf_force,8)
    tta=tta/real(nconf_force,8)
    if(nat_tot==0) nat_tot=1
    frmse=sqrt(frmse/real(3*nat_tot,8))
    if(iproc==0) then
        rmse=rmse*1.d3
        errmax=errmax*1.d3
        if(rmse>99999.d0) then
            !frmt=frmt2
            fmt_main=fmt2
        else
            !frmt=frmt1
            fmt_main=fmt1
        endif

        !write(ifile,frmt) iter,rmse,ttn,tta,frmse,errmax,ierrmax,atoms_arr%atoms(ierrmax)%nat, &
        !    ilarge1,ilarge2,ilarge3,trim(atoms_arr%fn(ierrmax)),atoms_arr%lconf(ierrmax)
        write(str_key,'(a)') trim(data_set)
        call yaml_mapping_open(trim(str_key),flow=.true.,unit=ann_arr%iunit)
        call yaml_map('iter',iter,unit=ann_arr%iunit)
        call yaml_map('rmse',rmse,fmt=trim(fmt_main),unit=ann_arr%iunit)
        call yaml_map('ttn',ttn,fmt=trim(fmt_main),unit=ann_arr%iunit)
        call yaml_map('tta',tta,fmt=trim(fmt_main),unit=ann_arr%iunit)
        call yaml_map('frmse',frmse,fmt=trim(fmt_main),unit=ann_arr%iunit)
        call yaml_map('dpmrmse',dpmrmse,fmt=trim(fmt_main),unit=ann_arr%iunit)
        call yaml_map('errmax',errmax,fmt=trim(fmt_main),unit=ann_arr%iunit)
        call yaml_map('ierrmax',ierrmax,unit=ann_arr%iunit)
        call yaml_map('nat',atoms_arr%atoms(ierrmax)%nat,unit=ann_arr%iunit)
        call yaml_map('ilarge1',ilarge1,unit=ann_arr%iunit)
        call yaml_map('ilarge2',ilarge2,unit=ann_arr%iunit)
        call yaml_map('ilarge3',ilarge3,unit=ann_arr%iunit)
        call yaml_map('fn_ierrmax',trim(atoms_arr%fn(ierrmax)),unit=ann_arr%iunit)
        call yaml_map('lconf_ierrmax',atoms_arr%lconf(ierrmax),unit=ann_arr%iunit)
        call yaml_mapping_close(unit=ann_arr%iunit)
    endif
    call cpu_time(time2)
    dtime1=time1-time_p
    dtime2=time2-time1
    !write(*,'(a,2f20.2)') 'TIME ',dtime1,dtime2
    time_p=time2
    if(parini%print_energy) then
        close(iunit)
    endif
    ann_arr%compute_symfunc=.false.
end subroutine ann_evaluate
!*****************************************************************************************
subroutine set_opt_ann_grad(ann_arr,grad,n,g)
    use mod_ann, only: typ_ann_arr
    implicit none
    type(typ_ann_arr), intent(in):: ann_arr
    real(8), intent(in):: grad(ann_arr%nweight_max,ann_arr%nann)
    integer, intent(in):: n
    real(8), intent(out):: g(n)
    !local variables
    integer:: iann, iw, i
    i=0
    do iann=1,ann_arr%nann
        do iw=1,ann_arr%num(iann)
            i=i+1
            g(ann_arr%loc(iann)+iw-1)=grad(iw,iann)
        enddo
    enddo
    if(i/=n) then
        write(*,'(a)') 'ERROR: i/=n in set_opt_ann_grad:'
        write(*,'(a,i8)') 'i= ',i
        write(*,'(a,i8)') 'n= ',n
        stop
    endif
end subroutine set_opt_ann_grad
!*****************************************************************************************
subroutine set_annweights(parini,opt_ann,ann_arr)
    use mod_parini, only: typ_parini
    use mod_processors, only: iproc
    use mod_opt_ann, only: typ_opt_ann
    use mod_ann, only: typ_ann_arr, convert_ann_x
    use mod_flm_futile
    use mod_utils
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_opt_ann), intent(inout):: opt_ann
    type(typ_ann_arr), intent(in):: ann_arr
    !local variables
    integer:: ierr, i, ia
    real(8):: tt
    character(50):: approach
    real(8), allocatable:: x(:)
#if defined(MPI)
    include 'mpif.h'
#endif
    x=f_malloc([1.to.opt_ann%n],id='x')
    if (parini%restart_param) then
        do ia=1,ann_arr%nann
            call convert_ann_x(ann_arr%num(ia),x(ann_arr%loc(ia)),ann_arr%ann(ia))
        enddo
        return
    endif
    !approach='uniform'
    approach='pure_electrostatic'
    !approach='type_dependent'
    if(iproc==0) then
        if(trim(approach)=='uniform') then
            if(trim(parini%rng_type)=='only_for_tests') then
                call random_number_generator_simple(opt_ann%n,x)
            else
                call random_number(x)
            endif
            x(1:opt_ann%n)=(x(1:opt_ann%n)-0.5d0)*2.d0*parini%ampl_rand
        elseif(trim(approach)=='pure_electrostatic') then
            do i=1,opt_ann%n
                if(i<=opt_ann%n/2) then
                if(trim(parini%rng_type)=='only_for_tests') then
                    call random_number_generator_simple(x(i))
                else
                    call random_number(x(i))
                endif
                    tt=-2.d0*parini%ampl_rand
                    x(i)=(x(i)-0.5d0)*tt
                else
                    x(i)=-x(i-opt_ann%n/2)
                endif
            enddo
        elseif(trim(approach)=='type_dependent') then
            do i=1,opt_ann%n
                if(i<=opt_ann%n/4) then
                    if(trim(parini%rng_type)=='only_for_tests') then
                        call random_number_generator_simple(x(i))
                    else
                        call random_number(x(i))
                    endif
                    tt=2.d0*parini%ampl_rand
                    x(i)=(x(i)-0.5d0)*tt
                elseif(i<=opt_ann%n/2) then
                    x(i)=-x(i-opt_ann%n/4)
                elseif(i<=3*opt_ann%n/4) then
                    if(trim(parini%rng_type)=='only_for_tests') then
                        call random_number_generator_simple(x(i))
                    else
                        call random_number(x(i))
                    endif
                    tt=2.d0*parini%ampl_rand
                    x(i)=-(x(i)-0.5d0)*tt
                else
                    x(i)=-x(i-opt_ann%n/4)
                endif
            enddo
        else
            stop 'ERROR: unknown approach in set_annweights'
        endif

    endif
#if defined(MPI)
    call MPI_BCAST(x,opt_ann%n,MPI_DOUBLE_PRECISION,0,parini%mpi_env%mpi_comm,ierr)
#endif
    call opt_ann%init_x_opt_ann(opt_ann%n,x)
    call f_free(x)
end subroutine set_annweights
!*****************************************************************************************
subroutine analyze_epoch_init(parini,ann_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    !local variables
    if(.not. (trim(ann_arr%approach)=='cent3' .or. trim(ann_arr%approach)=='cent2')) return
    ann_arr%natsum(1:10)=0
    ann_arr%qmin(1:10)=huge(1.d0)
    ann_arr%qmax(1:10)=-huge(1.d0)
    ann_arr%qsum(1:10)=0.d0
    ann_arr%chi_min(1:10)=huge(1.d0)
    ann_arr%chi_max(1:10)=-huge(1.d0)
    ann_arr%chi_sum(1:10)=0.d0
    ann_arr%chi_delta(1:10)=0.d0
end subroutine analyze_epoch_init
!*****************************************************************************************
subroutine analyze_epoch_print(parini,iter,ann_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_flm_futile
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iter
    type(typ_ann_arr), intent(inout):: ann_arr
    !local variables
    integer:: i, ios
    real(8):: ttavg, ttmin, ttmax, ssavg, ssmin, ssmax
    character(50):: fn_charge, fn_chi
    character(20):: str_key
    if(.not. (trim(ann_arr%approach)=='cent3' .or. trim(ann_arr%approach)=='cent2' )) return
    do i=1,parini%ntypat
        !fn_charge='charge.'//trim(parini%stypat(i))
        !fn_chi='chi.'//trim(parini%stypat(i))
        !if(iter==0) then
        !    open(unit=61,file=trim(fn_charge),status='replace',iostat=ios)
        !    if(ios/=0) then
        !        write(*,'(2a)') 'ERROR: failure openning ',trim(fn_charge)
        !        stop
        !    endif
        !    open(unit=71,file=trim(fn_chi),status='replace',iostat=ios)
        !    if(ios/=0) then
        !        write(*,'(2a)') 'ERROR: failure openning ',trim(fn_chi)
        !        stop
        !    endif
        !else
        !    open(unit=61,file=trim(fn_charge),status='old',position='append',iostat=ios)
        !    if(ios/=0) then
        !        write(*,'(2a)') 'ERROR: failure openning ',trim(fn_charge)
        !        stop
        !    endif
        !    open(unit=71,file=trim(fn_chi),status='old',position='append',iostat=ios)
        !    if(ios/=0) then
        !        write(*,'(2a)') 'ERROR: failure openning ',trim(fn_chi)
        !        stop
        !    endif
        !endif
        ttavg=ann_arr%qsum(i)/real(ann_arr%natsum(i),8)
        ttmin=ann_arr%qmin(i)
        ttmax=ann_arr%qmax(i)
        !write(61,'(i6,4f8.3)') iter,ttavg,ttmin,ttmax,ttmax-ttmin
        
        write(str_key,'(2a)') 'charge_',trim(parini%stypat(i))
        call yaml_mapping_open(trim(str_key),flow=.true.,unit=ann_arr%iunit)
        call yaml_map('iter',iter,unit=ann_arr%iunit)
        call yaml_map('qavg',ttavg,fmt='(f8.3)',unit=ann_arr%iunit)
        call yaml_map('qmin',ttmin,fmt='(f8.3)',unit=ann_arr%iunit)
        call yaml_map('qmax',ttmax,fmt='(f8.3)',unit=ann_arr%iunit)
        call yaml_map('qvar',ttmax-ttmin,fmt='(f8.3)',unit=ann_arr%iunit)
        call yaml_mapping_close(unit=ann_arr%iunit)

        !write(61,'(i6,4es14.5)') iter,ttavg,ttmin,ttmax,ttmax-ttmin
        ssavg=ann_arr%chi_sum(i)/real(ann_arr%natsum(i),8)
        ssmin=ann_arr%chi_min(i)
        ssmax=ann_arr%chi_max(i)
        !write(71,'(i6,5f8.3)') iter,ssavg,ssmin,ssmax,ssmax-ssmin,ann_arr%chi_delta(i)
        write(str_key,'(2a)') 'chi_',trim(parini%stypat(i))
        call yaml_mapping_open(trim(str_key),flow=.true.,unit=ann_arr%iunit)
        call yaml_map('iter',iter,unit=ann_arr%iunit)
        call yaml_map('chiavg',ssavg,fmt='(f8.3)',unit=ann_arr%iunit)
        call yaml_map('chimin',ssmin,fmt='(f8.3)',unit=ann_arr%iunit)
        call yaml_map('chimax',ssmax,fmt='(f8.3)',unit=ann_arr%iunit)
        call yaml_map('chivar',ssmax-ssmin,fmt='(f8.3)',unit=ann_arr%iunit)
        call yaml_mapping_close(unit=ann_arr%iunit)
        !write(71,'(i6,4es14.5)') iter,ssavg,ssmin,ssmax,ssmax-ssmin
        !if (trim(parini%stypat(i))=='O' .and. ssmax-ssmin> 0.01) stop
        !close(61)
        !close(71)
    enddo
end subroutine analyze_epoch_print
!*****************************************************************************************
subroutine convert_opt_x_ann_arr(opt_ann,ann_arr)
    use mod_ann, only: typ_ann_arr, convert_x_ann_arr
    use mod_opt_ann, only: typ_opt_ann
    use mod_flm_futile
    implicit none
    type(typ_opt_ann), intent(in):: opt_ann
    type(typ_ann_arr), intent(inout):: ann_arr
    !local variables
    real(8), allocatable:: x(:)
    x=f_malloc([1.to.opt_ann%n],id='x')
    call opt_ann%get_opt_ann_x(opt_ann%n,x)
    call convert_x_ann_arr(opt_ann%n,x,ann_arr)
    call f_free(x)
end subroutine convert_opt_x_ann_arr
!*****************************************************************************************
