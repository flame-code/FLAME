!*****************************************************************************************
!Subroutine cal_ann_main is called during training process.
!It does the same job as subroutine eval_cal_ann_main with a few more
!commands related to ANN weights.
!These two subroutine are supposed to be merged.
!eval_cal_ann_main is called for task_training but only for
!obtaining energy/forces when information on errors on
!training and validation data points are expected.
!eval_cal_ann_main is called by potential_ANN
subroutine cal_ann_main(parini,atoms,symfunc,ann_arr,ekf)
    use mod_interface
    use mod_tightbinding, only: typ_partb
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr, typ_symfunc, typ_ekf
    !use mod_tightbinding, only: typ_partb
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc), intent(inout):: symfunc
    type(typ_ekf), intent(inout):: ekf
    !local variables
    !real(8):: g, g_tb, dis, E0, E1 
    integer:: i, j, iat
    type(typ_partb):: partb
    if(trim(ann_arr%approach)=='atombased') then
        allocate(ekf%gs(ekf%num(1),atoms%nat)) !HERE
        call convert_x_ann(ekf%num(1),ekf%x(ekf%loc(1)),ann_arr%ann(1))
        call cal_ann_atombased(parini,atoms,symfunc,ann_arr,ekf)
        ekf%g(1:ekf%n)=0.d0
        do iat=1,atoms%nat
            i=atoms%itypat(iat)
            do j=1,ekf%num(1)
                ekf%g(ekf%loc(i)+j-1)=ekf%g(ekf%loc(i)+j-1)+ekf%gs(j,iat)
            enddo
        enddo
        deallocate(ekf%gs) !HERE
    elseif(trim(ann_arr%approach)=='eem1' .or. trim(ann_arr%approach)=='cent1') then
        call cal_ann_eem1(parini,atoms,symfunc,ann_arr,ekf)
    elseif(trim(ann_arr%approach)=='cent2') then
        call cal_ann_eem2(parini,atoms,symfunc,ann_arr,ekf)
    elseif(trim(ann_arr%approach)=='tb') then
        call cal_ann_tb(parini,partb,atoms,ann_arr,symfunc,ekf)
       ! do i=1,ekf%n
       ! dis=1.d-4
       ! g_tb=ekf%g(i)
       ! !!Finite difference 
       !     E0=atoms%epot
       !     ekf%x(i)=ekf%x(i)+dis
       !     call cal_ann_tb(parini,partb,atoms,ann_arr,symfunc,ekf)
       !     E1=atoms%epot
       !     g=(E1-E0)/dis
       !     write(*,'(a,2es19.10,es14.5,2es19.10)') 'FD-TEST',g_tb,g,g-g_tb,E0,E1
       ! !endif
       ! enddo
       ! stop 'TTTTTT'
    else
        write(*,'(2a)') 'ERROR: unknown approach in ANN, ',trim(ann_arr%approach)
        stop
    endif
end subroutine cal_ann_main
!*****************************************************************************************
subroutine pre_train(parini,ann_arr,symfunc_train,symfunc_valid,atoms_train,atoms_valid,ekf)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr, typ_ekf
    use mod_atoms, only: typ_atoms, typ_atoms_arr
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc_arr), intent(inout):: symfunc_train, symfunc_valid
    type(typ_atoms_arr), intent(inout):: atoms_train
    type(typ_atoms_arr), intent(inout):: atoms_valid
    type(typ_ekf), intent(inout):: ekf
    !local variables
    type(typ_atoms):: atoms
    integer:: iconf, istep, iat, ia
    real(8):: anat(100), g(100), rmse
    !return
    ann_arr%event='evalu'
    do ia=1,ann_arr%n
        call convert_x_ann(ekf%num(ia),ekf%x(ekf%loc(ia)),ann_arr%ann(ia))
    enddo
    do istep=1,50
        rmse=0.d0
        g=0.d0
        do iconf=1,atoms_train%nconf
            call atom_copy_old(atoms_train%atoms(iconf),atoms,'atoms_train%atoms(iconf)->atoms')
            call cal_ann_main(parini,atoms,symfunc_train%symfunc(iconf),ann_arr,ekf)
            rmse=rmse+((symfunc_train%symfunc(iconf)%epot-atoms%epot)/atoms%nat)**2
            anat=0.d0
            do iat=1,atoms%nat
                anat(atoms%itypat(iat))=anat(atoms%itypat(iat))+1.d0
            enddo
            do ia=1,ann_arr%n
                g(ia)=g(ia)+2.d0*anat(ia)*(atoms%epot-symfunc_train%symfunc(iconf)%epot)/atoms%nat**2
            enddo
        enddo
        rmse=sqrt(rmse/atoms_train%nconf)
        write(*,'(a,i4,3es19.10)') 'pretrain: ',istep,rmse*1.d3, &
            ann_arr%ann(1)%ener_ref,ann_arr%ann(2)%ener_ref
        do ia=1,ann_arr%n
            ann_arr%ann(ia)%ener_ref=ann_arr%ann(ia)%ener_ref-2.d-2*g(ia)
        enddo
    enddo
    !stop
end subroutine pre_train
!*****************************************************************************************
