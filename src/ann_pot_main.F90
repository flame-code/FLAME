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
    !if(trim(ann_arr%event)=='potential') then
    !    allocate(symfunc%y(ann_arr%ann(1)%nn(0),atoms%nat))
    !endif
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
    elseif(trim(ann_arr%approach)=='eem1') then
        call cal_ann_eem1(parini,atoms,symfunc,ann_arr,ekf)
    elseif(trim(ann_arr%approach)=='eem2') then
        allocate(ekf%gc(ekf%num(1),atoms%nat)) !HERE
        allocate(ekf%gs(ekf%num(1),atoms%nat)) !HERE
        call convert_x_ann(ekf%num(1),ekf%x(ekf%loc(1)),ann_arr%ann(1))
        call convert_x_ann(ekf%num(2),ekf%x(ekf%loc(2)),ann_arr%ann(2))
        !call convert_x_ann(ekf%num(3),ekf%x(ekf%loc(3)),ann_arr%ann(3))
        !call convert_x_ann(ekf%num(4),ekf%x(ekf%loc(4)),ann_arr%ann(4))
        call cal_ann_eem2(parini,atoms,symfunc,ann_arr,ekf)
        ekf%g(1:ekf%n)=0.d0
        do iat=1,atoms%nat
            i=atoms%itypat(iat)
            do j=1,ekf%num(1)
                !ekf%g(ekf%loc(i)+j-1)=ekf%g(ekf%loc(i)+j-1)+atoms%qat(iat)*ekf%gc(j,iat)
                ekf%g(ekf%loc(i)+j-1)=ekf%g(ekf%loc(i)+j-1)+(atoms%qat(iat)+atoms%zat(iat))*ekf%gc(j,iat)
            enddo
!            do j=1,ekf%num(1)
!                ekf%g(ekf%loc(i)+j-1)=ekf%g(ekf%loc(i)+j-1)+ekf%gs(j,iat)
!            enddo
        enddo
        deallocate(ekf%gc) !HERE
        deallocate(ekf%gs) !HERE
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
    !if(trim(ann_arr%event)=='potential') then
    !    deallocate(symfunc%y)
    !endif
end subroutine cal_ann_main
!*****************************************************************************************
