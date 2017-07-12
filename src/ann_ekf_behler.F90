!*****************************************************************************************
subroutine ekf_behler(parini,ann_arr,symfunc_train,symfunc_valid,atoms_train,atoms_valid,ekf)
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
    real(8):: DDOT, tt, den, alambda, alambdainv, alambda0
    real(8):: time_s, time_e, time1, time2, time3 !, time4
    real(8):: dtime, dtime1, dtime2, dtime3, dtime4, dtime5, dtime6
    real(8):: tt1, tt2, tt3, tt4, tt5, tt6
    allocate(f(ekf%n),p(ekf%n,ekf%n),ekf%g(ekf%n),v1(ekf%n),ekf%epotd(ekf%num(1)))
    p(1:ekf%n,1:ekf%n)=0.d0
    do i=1,ekf%n
        p(i,i)=1.d-2
    enddo
    !alambda0=0.9994d0
    !alambda=0.99d0
    alambda0=0.9998d0
    alambda=0.997d0
    do iter=0,parini%nstep_ekf
        call cpu_time(time_s)
        !call randomize_data_order(atoms_train)
        do ia=1,ann_arr%n
            call convert_x_ann(ekf%num(ia),ekf%x(ekf%loc(ia)),ann_arr%ann(ia))
        enddo
        if(iproc==0) then
            if( ann_arr%exists_yaml_file) then
                call write_ann_all_yaml(parini,ann_arr,iter)
            else
                call write_ann_all(parini,ann_arr,iter) 
            endif
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
        do iconf=1,atoms_train%nconf
            alambda=min(alambda0*alambda+1.d0-alambda0,0.999999d0)
            alambdainv=1.d0/alambda
            !write(31,'(2i6,f14.10)') iter,iconf,alambda
            ann_arr%event='train'
            call atom_copy_old(atoms_train%atoms(iconf),atoms,'atoms_train%atoms(iconf)->atoms')
            call cal_ann_main(parini,atoms,symfunc_train%symfunc(iconf),ann_arr,ekf)
            call cpu_time(time1)
            !ekf%g(1:ekf%n)=ekf%g(1:ekf%n)*(ekf%epot-epotall(iconf))
            call DGEMV('T',ekf%n,ekf%n,1.d0,p,ekf%n,ekf%g,1,0.d0,v1,1)
            !call cal_matvec_mpi(ekf%n,p,ekf%g,v1)
            call cpu_time(time2)
            tt=DDOT(ekf%n,ekf%g,1,v1,1)
            den=alambdainv/(1.d0+tt*alambdainv)
            !write(*,'(i7,i5,f14.10,es11.2)') iter,iconf,alambda,den
            !call DGEMV('T',ekf%n,ekf%n,den,p,ekf%n,ekf%g,1,0.d0,f,1)
            f(1:ekf%n)=v1(1:ekf%n)*den
            do j=1,ekf%n
                do i=1,ekf%n
                    p(i,j)=(p(i,j)-f(i)*v1(j))*alambdainv
                    !if(i==j) p(i,j)=p(i,j)+1.d-1
                enddo
            enddo
                    !write(21,'(a)') '----------------------------------------'
            !write(*,'(a,i7,i6,2es15.5)') 'forgetting ',iter,iconf,alambda,den
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
end subroutine ekf_behler
!*****************************************************************************************
