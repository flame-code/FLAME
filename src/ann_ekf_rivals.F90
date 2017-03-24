!*****************************************************************************************
subroutine ekf_rivals(parini,ann_arr,symfunc_train,symfunc_valid,atoms_train,atoms_valid,ekf)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr, typ_ekf
    use mod_atoms, only: typ_atoms, typ_atoms_arr
    use mod_processors, only: iproc, mpi_comm_abz
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc_arr), intent(inout):: symfunc_train, symfunc_valid
    type(typ_atoms_arr), intent(inout):: atoms_train
    type(typ_atoms_arr), intent(inout):: atoms_valid
    type(typ_ekf), intent(inout):: ekf
    !local variables
    real(8), allocatable:: f(:) !Kalman gain matrix
    real(8), allocatable:: p(:,:) !covariance matrix
    real(8), allocatable:: v1(:) !work array
    type(typ_atoms):: atoms
    integer:: i, j, iter, iconf, ios, ia
    real(8):: DDOT, tt, den
    real(8):: r, rinv, r0, rf, alpha
    real(8):: time_s, time_e, time1, time2, time3 !, time4
    real(8):: dtime, dtime1, dtime2, dtime3, dtime4, dtime5, dtime6
    real(8):: tt1, tt2, tt3, tt4, tt5, tt6
    allocate(f(ekf%n),p(ekf%n,ekf%n),ekf%g(ekf%n),v1(ekf%n),ekf%epotd(ekf%num(1)))
    p(1:ekf%n,1:ekf%n)=0.d0
    do i=1,ekf%n
        p(i,i)=1.d-2
    enddo
    if(trim(parini%approach_ann)=='eem1') then
        r0=100000.d0
        alpha=120.d-2
        rf=1.d-6
    elseif(trim(parini%approach_ann)=='eem2') then
        r0=100000.d0
        alpha=120.d-2
        rf=1.d-6
    elseif(trim(parini%approach_ann)=='tb') then
        r0=100000.d0
        alpha=120.d-2
        rf=1.d-6
    else
        r0=1.d0
        alpha=5.d-1
        rf=1.d-8
    endif
    do iter=0,parini%nstep_ekf
        call cpu_time(time_s)
        !call randomize_data_order(atoms_train)
        do ia=1,ann_arr%n
            call convert_x_ann(ekf%num(ia),ekf%x(ekf%loc(ia)),ann_arr%ann(ia))
        enddo
        if(iproc==0) then
            call write_ann_all(parini,ann_arr,iter)
        endif
        if(mod(iter,1)==0) then
            call analyze_epoch_init(parini,atoms_train,ann_arr)
            call ann_evaluate(parini,iter,ann_arr,symfunc_train,atoms_train,11)
            call ann_evaluate(parini,iter,ann_arr,symfunc_valid,atoms_valid,12)
            call analyze_epoch_print(parini,iter,atoms_train,ann_arr)
        endif
        if(iter==parini%nstep_ekf) exit
        call cpu_time(time1)
        dtime1=time1-time_s !time to get ANN energies for atoms_train and atoms_valid
        dtime2=0.d0 !time to convert ANN 1D array to ANN typ_ann
        dtime3=0.d0 !time to calculate ANN and its derivatives w.r.t. weights
        dtime4=0.d0 !time to convert derivative of ANN in typ_ann to 1D array
        dtime5=0.d0 !time to matrix-vector multiplication in Kalman filter
        dtime6=0.d0 !time of the rest of Kalman filter algorithm
        r=(r0-rf)*exp(-alpha*(iter))+rf
        rinv=1.d0/r
        write(31,'(i6,es14.5)') iter,r
        do iconf=1,atoms_train%nconf
            ann_arr%event='train'
            call atom_copy_old(atoms_train%atoms(iconf),atoms,'atoms_train%atoms(iconf)->atoms')
            call cal_ann_main(parini,atoms,symfunc_train%symfunc(iconf),ann_arr,ekf)
            call cpu_time(time1)
            call DGEMV('T',ekf%n,ekf%n,1.d0,p,ekf%n,ekf%g,1,0.d0,v1,1)
            !call cal_matvec_mpi(ekf%n,p,ekf%g,v1)
            call cpu_time(time2)
            tt=DDOT(ekf%n,ekf%g,1,v1,1)
            den=1.d0/(tt+r)
            do j=1,ekf%n
                do i=1,ekf%n
                    p(i,j)=p(i,j)-v1(i)*v1(j)*den
                    !if(i==j) p(i,j)=p(i,j)+1.d-1
                    !write(21,'(2i5,es20.10)') i,j,p(i,j)
                enddo
            enddo
                    !write(21,'(a)') '----------------------------------------'
            !write(*,'(a,i7,i6,2es15.5)') 'forgetting ',iter,iconf,r,den
            call DGEMV('N',ekf%n,ekf%n,rinv,p,ekf%n,ekf%g,1,0.d0,f,1)
            do i=1,ekf%n
                !write(81,*) iter,iconf,i,f(i)*(epotall(iconf)-ekf%epot)
                ekf%x(i)=ekf%x(i)+f(i)*(symfunc_train%symfunc(iconf)%epot-atoms%epot)
            enddo
            call cpu_time(time3)
            dtime5=dtime5+time2-time1
            dtime6=dtime6+time3-time2
        enddo
        call cpu_time(time_e)
        dtime=time_e-time_s
        !dtime=dtime1+dtime2+dtime3+dtime4+dtime5+dtime6
        tt1=dtime1/dtime*100.d0
        tt2=dtime2/dtime*100.d0
        tt3=dtime3/dtime*100.d0
        tt4=dtime4/dtime*100.d0
        tt5=dtime5/dtime*100.d0
        tt6=dtime6/dtime*100.d0
        tt=tt1+tt2+tt3+tt4+tt5+tt6
        if(iter==0) then
            open(unit=41,file='time_real.prc',status='replace',iostat=ios)
            if(ios/=0) then;write(*,'(a)') 'ERROR: failure openning output file';stop;endif
            write(41,'(a1,4x,7a10)') '#',' dtime1 ',' dtime2 ',' dtime3 ',' dtime4 ',' dtime5 ',' dtime6 ',' sum '
            open(unit=42,file='time_frac.prc',status='replace',iostat=ios)
            if(ios/=0) then;write(*,'(a)') 'ERROR: failure openning output file';stop;endif
            write(42,'(a1,4x,7a10)') '#','    tt1 ','    tt2 ','    tt3 ','    tt4 ','    tt5 ','    tt6 ', 'sum '
        endif
        write(41,'(i5,7f10.2)') iter,dtime1,dtime2,dtime3,dtime4,dtime5,dtime6,dtime
        write(42,'(i5,7f10.1)') iter,tt1,tt2,tt3,tt4,tt5,tt6,tt
    enddo
    close(41)
    close(42)
    deallocate(f,p,ekf%g,v1,ekf%epotd)
end subroutine ekf_rivals
!*****************************************************************************************
subroutine ekf_rivals_tmp(parini,ann_arr,symfunc_train,symfunc_valid,atoms_train,atoms_valid,ekf)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr, typ_ekf
    use mod_atoms, only: typ_atoms, typ_atoms_arr
    use mod_processors, only: iproc, mpi_comm_abz
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc_arr), intent(inout):: symfunc_train, symfunc_valid
    type(typ_atoms_arr), intent(inout):: atoms_train
    type(typ_atoms_arr), intent(inout):: atoms_valid
    type(typ_ekf), intent(inout):: ekf
    !local variables
    real(8), allocatable:: f(:) !Kalman gain matrix
    real(8), allocatable:: p(:,:) !covariance matrix
    real(8), allocatable:: v1(:) !work array
    type(typ_atoms):: atoms, atoms_ann_ref
    integer:: i, j, iter, iconf, ios, ia
    real(8):: DDOT, tt, den
    real(8):: r, rinv, r0, rf, alpha
    character(16):: fn
    character(50):: filename
    real(8):: time_s, time_e, time1, time2, time3 !, time4
    real(8):: dtime, dtime1, dtime2, dtime3, dtime4, dtime5, dtime6
    real(8):: tt1, tt2, tt3, tt4, tt5, tt6
    type(typ_atoms_arr):: atoms_ref
    integer:: ind(200), mat
    real(8):: de_ann, de_ref, de
    allocate(f(ekf%n),p(ekf%n,ekf%n),ekf%g(ekf%n),v1(ekf%n),ekf%epotd(ekf%num(1)))
    p(1:ekf%n,1:ekf%n)=0.d0
    do i=1,ekf%n
        p(i,i)=1.d-2
    enddo
    if(trim(parini%approach_ann)=='eem1') then
        r0=300.d0
        alpha=20.d-2
        rf=1.d-8
    elseif(trim(parini%approach_ann)=='eem2') then
        r0=100.d0
        alpha=80.d-2
        rf=1.d-8
    else
        r0=1.d0
        alpha=5.d-1
        rf=1.d-8
    endif
    call set_ref_energy(parini,atoms_train,atoms_ref,ind)
    do iter=0,parini%nstep_ekf
        call cpu_time(time_s)
        !call randomize_data_order(atoms_train)
        do ia=1,ann_arr%n
            call convert_x_ann(ekf%num(ia),ekf%x(ekf%loc(ia)),ann_arr%ann(ia))
        enddo
        if(iproc==0) then
            write(fn,'(a11,i5.5)') '.ann.param.',iter
            do i=1,ann_arr%n
                filename=trim(parini%stypat(i))//trim(fn)
                write(*,'(a)') trim(filename)
                call write_ann(parini,filename,ann_arr%ann(i))
            enddo
        endif
        if(mod(iter,1)==0) then
            call ann_evaluate(parini,iter,ann_arr,symfunc_train,atoms_train,11)
            call ann_evaluate(parini,iter,ann_arr,symfunc_valid,atoms_valid,12)
        endif
        if(iter==parini%nstep_ekf) exit
        call cpu_time(time1)
        dtime1=time1-time_s !time to get ANN energies for atoms_train and atoms_valid
        dtime2=0.d0 !time to convert ANN 1D array to ANN typ_ann
        dtime3=0.d0 !time to calculate ANN and its derivatives w.r.t. weights
        dtime4=0.d0 !time to convert derivative of ANN in typ_ann to 1D array
        dtime5=0.d0 !time to matrix-vector multiplication in Kalman filter
        dtime6=0.d0 !time of the rest of Kalman filter algorithm
        r=(r0-rf)*exp(-alpha*(iter))+rf
        rinv=1.d0/r
        !write(31,'(i6,es14.5)') iter,r
        do iconf=1,atoms_train%nconf
            ann_arr%event='train'

            mat=atoms_train%atoms(iconf)%nat
            call atom_copy_old(atoms_ref%atoms(mat),atoms_ann_ref,'atoms_train%atoms(iconf)->atoms')
            call cal_ann_main(parini,atoms_ann_ref,symfunc_train%symfunc(ind(mat)),ann_arr,ekf)

            call atom_copy_old(atoms_train%atoms(iconf),atoms,'atoms_train%atoms(iconf)->atoms')
            call cal_ann_main(parini,atoms,symfunc_train%symfunc(iconf),ann_arr,ekf)
            call cpu_time(time1)
            call DGEMV('T',ekf%n,ekf%n,1.d0,p,ekf%n,ekf%g,1,0.d0,v1,1)
            !call cal_matvec_mpi(ekf%n,p,ekf%g,v1)
            call cpu_time(time2)
            tt=DDOT(ekf%n,ekf%g,1,v1,1)
            den=1.d0/(tt+r)
            do j=1,ekf%n
                do i=1,ekf%n
                    p(i,j)=p(i,j)-v1(i)*v1(j)*den
                    !if(i==j) p(i,j)=p(i,j)+1.d-1
                    !write(21,'(2i5,es20.10)') i,j,p(i,j)
                enddo
            enddo
                    !write(21,'(a)') '----------------------------------------'
            !write(*,'(a,i7,i6,2es15.5)') 'forgetting ',iter,iconf,r,den
            call DGEMV('N',ekf%n,ekf%n,rinv,p,ekf%n,ekf%g,1,0.d0,f,1)
            if(iter<0) then
                de_ann=atoms%epot
                de_ref=symfunc_train%symfunc(iconf)%epot
            else
                de_ann=atoms%epot-atoms_ann_ref%epot
                de_ref=symfunc_train%symfunc(iconf)%epot-atoms_ref%atoms(mat)%epot
            endif
            write(89,'(i2,i4,2es14.5)') iter,mat,abs(de_ref-de_ann),abs(symfunc_train%symfunc(iconf)%epot-atoms%epot)
            de=de_ref-de_ann
            de=sign(min(2.d-2/(iter*0+1),abs(de)),de)
            do i=1,ekf%n
                !write(81,*) iter,iconf,i,f(i)*(epotall(iconf)-ekf%epot)
                ekf%x(i)=ekf%x(i)+f(i)*(min((iter+1)*1.d0,1.d0)*de) !+max(1.d0-iter*2.d-2,0.d0)*(symfunc_train%symfunc(iconf)%epot-atoms%epot))
            enddo
            call cpu_time(time3)
            dtime5=dtime5+time2-time1
            dtime6=dtime6+time3-time2
        enddo
        call cpu_time(time_e)
        dtime=time_e-time_s
        !dtime=dtime1+dtime2+dtime3+dtime4+dtime5+dtime6
        tt1=dtime1/dtime*100.d0
        tt2=dtime2/dtime*100.d0
        tt3=dtime3/dtime*100.d0
        tt4=dtime4/dtime*100.d0
        tt5=dtime5/dtime*100.d0
        tt6=dtime6/dtime*100.d0
        tt=tt1+tt2+tt3+tt4+tt5+tt6
        if(iter==0) then
            open(unit=41,file='time_real.prc',status='replace',iostat=ios)
            if(ios/=0) then;write(*,'(a)') 'ERROR: failure openning output file';stop;endif
            write(41,'(a1,4x,7a10)') '#',' dtime1 ',' dtime2 ',' dtime3 ',' dtime4 ',' dtime5 ',' dtime6 ',' sum '
            open(unit=42,file='time_frac.prc',status='replace',iostat=ios)
            if(ios/=0) then;write(*,'(a)') 'ERROR: failure openning output file';stop;endif
            write(42,'(a1,4x,7a10)') '#','    tt1 ','    tt2 ','    tt3 ','    tt4 ','    tt5 ','    tt6 ', 'sum '
        endif
        write(41,'(i5,7f10.2)') iter,dtime1,dtime2,dtime3,dtime4,dtime5,dtime6,dtime
        write(42,'(i5,7f10.1)') iter,tt1,tt2,tt3,tt4,tt5,tt6,tt
    enddo
    close(41)
    close(42)
    deallocate(f,p,ekf%g,v1,ekf%epotd)
end subroutine ekf_rivals_tmp
!*****************************************************************************************
subroutine set_ref_energy(parini,atoms_train,atoms_ref,ind)
    use mod_interface
    use mod_parini, only: typ_parini
    !use mod_ann, only: typ_symfunc_arr
    use mod_atoms, only: typ_atoms, typ_atoms_arr
    implicit none
    type(typ_parini), intent(in):: parini
    !type(typ_symfunc_arr), intent(in):: symfunc_train, symfunc_valid
    type(typ_atoms_arr), intent(in):: atoms_train
    !type(typ_atoms_arr), intent(in):: atoms_valid
    type(typ_atoms_arr), intent(inout):: atoms_ref
    integer, intent(out):: ind(200)
    !local variables
    !real(8), allocatable:: 
    integer:: iconf, mat
    atoms_ref%nconf=200
    allocate(atoms_ref%atoms(atoms_ref%nconf))
    do mat=1,200
        atoms_ref%atoms(mat)%epot=huge(1.d0)
    enddo
    write(*,*) 'atoms_train%nconf ',atoms_train%nconf
    ind(1:200)=0
    do iconf=1,atoms_train%nconf
        write(*,*) 'nat ',atoms_train%atoms(iconf)%nat
        mat=atoms_train%atoms(iconf)%nat
        if(atoms_train%atoms(iconf)%epot<atoms_ref%atoms(mat)%epot) then
            ind(mat)=iconf
            atoms_ref%atoms(mat)%epot=atoms_train%atoms(iconf)%epot
        endif
    enddo
    do mat=1,200
        if(ind(mat)==0) cycle
        iconf=ind(mat)
        call atom_copy_old(atoms_train%atoms(iconf),atoms_ref%atoms(mat),'copy to atoms_ref')
    enddo
end subroutine set_ref_energy
!*****************************************************************************************
subroutine analyze_epoch_init(parini,atoms_train,ann_arr)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr
    use mod_ann, only: typ_ann_arr
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms_arr), intent(in):: atoms_train
    type(typ_ann_arr), intent(inout):: ann_arr
    !local variables
    if(.not. (trim(ann_arr%approach)=='eem1' .or. trim(ann_arr%approach)=='eem2')) return
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
subroutine analyze_epoch_print(parini,iter,atoms_train,ann_arr)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr
    use mod_ann, only: typ_ann_arr
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iter
    type(typ_atoms_arr), intent(in):: atoms_train
    type(typ_ann_arr), intent(inout):: ann_arr
    !local variables
    integer:: i, ios
    real(8):: ttavg, ttmin, ttmax, ssavg, ssmin, ssmax
    character(50):: fn_charge, fn_chi
    if(.not. (trim(ann_arr%approach)=='eem1' .or. trim(ann_arr%approach)=='eem2')) return
    do i=1,parini%ntypat
        fn_charge='charge.'//trim(parini%stypat(i))
        fn_chi='chi.'//trim(parini%stypat(i))
        if(iter==0) then
            open(unit=61,file=trim(fn_charge),status='replace',iostat=ios)
            if(ios/=0) then
                write(*,'(2a)') 'ERROR: failure openning ',trim(fn_charge)
                stop
            endif
            open(unit=71,file=trim(fn_chi),status='replace',iostat=ios)
            if(ios/=0) then
                write(*,'(2a)') 'ERROR: failure openning ',trim(fn_chi)
                stop
            endif
        else
            open(unit=61,file=trim(fn_charge),status='old',position='append',iostat=ios)
            if(ios/=0) then
                write(*,'(2a)') 'ERROR: failure openning ',trim(fn_charge)
                stop
            endif
            open(unit=71,file=trim(fn_chi),status='old',position='append',iostat=ios)
            if(ios/=0) then
                write(*,'(2a)') 'ERROR: failure openning ',trim(fn_chi)
                stop
            endif
        endif
        ttavg=ann_arr%qsum(i)/real(ann_arr%natsum(i),8)
        ttmin=ann_arr%qmin(i)
        ttmax=ann_arr%qmax(i)
        write(61,'(i6,4f8.3)') iter,ttavg,ttmin,ttmax,ttmax-ttmin
        !write(61,'(i6,4es14.5)') iter,ttavg,ttmin,ttmax,ttmax-ttmin
        ssavg=ann_arr%chi_sum(i)/real(ann_arr%natsum(i),8)
        ssmin=ann_arr%chi_min(i)
        ssmax=ann_arr%chi_max(i)
        write(71,'(i6,5f8.3)') iter,ssavg,ssmin,ssmax,ssmax-ssmin,ann_arr%chi_delta(i)
        !write(71,'(i6,4es14.5)') iter,ssavg,ssmin,ssmax,ssmax-ssmin
        close(61)
        close(71)
    enddo
end subroutine analyze_epoch_print
!*****************************************************************************************
