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
    !real(8), allocatable:: xt(:), gt(:)
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
        call cal_ann_cent1(parini,atoms,symfunc,ann_arr,ekf)
    elseif(trim(ann_arr%approach)=='cent2') then
        call cal_ann_eem2(parini,atoms,symfunc,ann_arr,ekf)
    elseif(trim(ann_arr%approach)=='tb') then
        call cal_ann_tb(parini,partb,atoms,ann_arr,symfunc,ekf)
        ! E0=atoms%epot
        ! allocate(xt(ekf%n),gt(ekf%n))
        ! xt(1:ekf%n)=ekf%x(1:ekf%n)
        ! gt(1:ekf%n)=ekf%g(1:ekf%n)
        ! do i=1,ekf%n
        !     g_tb=gt(i)
        !     !!Finite difference 
        !     dis=1.d-4 !*abs(ekf%x(i))
        !     ekf%x(i)=ekf%x(i)+dis
        !     call cal_ann_tb(parini,partb,atoms,ann_arr,symfunc,ekf)
        !     E1=atoms%epot
        !     g=(E1-E0)/dis
        !     write(*,'(a,2es19.10,es14.5,2es19.10)') 'FD-TEST',g_tb,g,g-g_tb,E0,E1
        !     ekf%x(i)=xt(i)
        ! enddo
        ! stop 'TTTTTTTTTTTTTTTT'
    else
        write(*,'(2a)') 'ERROR: unknown approach in ANN, ',trim(ann_arr%approach)
        stop
    endif
end subroutine cal_ann_main
!*****************************************************************************************
subroutine prefit_cent_ener_ref(parini,ann_arr,symfunc_train,symfunc_valid,atoms_train,atoms_valid,ekf)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr, typ_ekf
    use mod_atoms, only: typ_atoms, typ_atoms_arr
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc_arr), intent(inout):: symfunc_train, symfunc_valid
    type(typ_atoms_arr), intent(inout):: atoms_train
    type(typ_atoms_arr), intent(inout):: atoms_valid
    type(typ_ekf), intent(inout):: ekf
    !local variables
    type(typ_atoms):: atoms
    integer:: iconf, istep, iat, ia, isatur, nsatur
    real(8):: anat(100), g(100), rmse, rmse_old, de0, alpha, tt
    real(8), allocatable:: epotall(:), eref_all(:)
    !return
    ann_arr%event='evalu'
    do ia=1,ann_arr%n
        call convert_x_ann(ekf%num(ia),ekf%x(ekf%loc(ia)),ann_arr%ann(ia))
    enddo
    nsatur=3
    isatur=0
    epotall=f_malloc([1.to.atoms_train%nconf],id='epotall')
    eref_all=f_malloc([1.to.atoms_train%nconf],id='eref_all')
    do iconf=1,atoms_train%nconf
        tt=0.d0
        do iat=1,atoms_train%atoms(iconf)%nat
            tt=tt+ann_arr%ann(atoms_train%atoms(iconf)%itypat(iat))%ener_ref
        enddo
        eref_all(iconf)=tt
    enddo
    alpha=1.d0/real(atoms_train%nconf,8)
    do istep=0,50
        rmse=0.d0
        g=0.d0
        do iconf=1,atoms_train%nconf
            call atom_copy_old(atoms_train%atoms(iconf),atoms,'atoms_train%atoms(iconf)->atoms')
            if(istep==0) then
                call cal_ann_main(parini,atoms,symfunc_train%symfunc(iconf),ann_arr,ekf)
                epotall(iconf)=atoms%epot
            else
                tt=0.d0
                do iat=1,atoms%nat
                    tt=tt+ann_arr%ann(atoms%itypat(iat))%ener_ref
                enddo
                atoms%epot=epotall(iconf)-eref_all(iconf)+tt
            endif
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
        if(istep==0) rmse_old=rmse
        if(istep>0 .and. rmse<rmse_old .and. abs(rmse_old-rmse)<1.d-3) then
            isatur=isatur+1
        else
            isatur=0
        endif
        write(*,'(a,i4,3es19.10,i3)') 'pretrain: ',istep,rmse*1.d3, &
            ann_arr%ann(1)%ener_ref,ann_arr%ann(2)%ener_ref,isatur
        if(rmse*1.d3<2.d1) exit
        if(isatur>nsatur) exit
        do ia=1,ann_arr%n
            de0=alpha*g(ia)
            de0=sign(min(abs(de0),1.d0),de0)
            ann_arr%ann(ia)%ener_ref=ann_arr%ann(ia)%ener_ref-de0
        enddo
        rmse_old=rmse
    enddo
    call f_free(epotall)
    call f_free(eref_all)
    !stop
end subroutine prefit_cent_ener_ref
!*****************************************************************************************
subroutine prefit_cent(parini,ann_arr,symfunc_train,symfunc_valid,atoms_train,atoms_valid,ekf)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr, typ_ekf
    use mod_atoms, only: typ_atoms, typ_atoms_arr
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc_arr), intent(inout):: symfunc_train, symfunc_valid
    type(typ_atoms_arr), intent(inout):: atoms_train
    type(typ_atoms_arr), intent(inout):: atoms_valid
    type(typ_ekf), intent(inout):: ekf
    !local variables
    type(typ_atoms):: atoms
    integer:: iconf, istep, iat, ia, isatur, nsatur
    real(8):: anat(100), g(100), rmse, rmse_old, de0, alpha, tt
    ann_arr%event='evalu'
    do ia=1,ann_arr%n
        call convert_x_ann(ekf%num(ia),ekf%x(ekf%loc(ia)),ann_arr%ann(ia))
    enddo
    nsatur=3
    isatur=0
    alpha=2.d0/real(atoms_train%nconf,8)
    do istep=0,50
        rmse=0.d0
        g=0.d0
        do iconf=1,atoms_train%nconf
            call atom_copy_old(atoms_train%atoms(iconf),atoms,'atoms_train%atoms(iconf)->atoms')
            call cal_ann_main(parini,atoms,symfunc_train%symfunc(iconf),ann_arr,ekf)
            rmse=rmse+((symfunc_train%symfunc(iconf)%epot-atoms%epot)/atoms%nat)**2
            anat=0.d0
            do iat=1,atoms%nat
                anat(atoms%itypat(iat))=anat(atoms%itypat(iat))+atoms%qat(iat)
            enddo
            do ia=1,ann_arr%n
                g(ia)=g(ia)+2.d0*anat(ia)*(atoms%epot-symfunc_train%symfunc(iconf)%epot)/atoms%nat**2
            enddo
        enddo
        rmse=sqrt(rmse/atoms_train%nconf)
        if(istep==0) rmse_old=rmse
        if(istep>0 .and. rmse<rmse_old .and. abs(rmse_old-rmse)<1.d-3) then
            isatur=isatur+1
        else
            isatur=0
        endif
        write(*,'(a,i4,3es19.10,i3)') 'prefit: ',istep,rmse*1.d3, &
            ann_arr%ann(1)%chi0,ann_arr%ann(2)%chi0,isatur
        if(rmse*1.d3<3.d0) exit
        if(isatur>nsatur) exit
        do ia=1,ann_arr%n
            de0=alpha*g(ia)
            de0=sign(min(abs(de0),1.d0),de0)
            ann_arr%ann(ia)%chi0=ann_arr%ann(ia)%chi0-de0
        enddo
        rmse_old=rmse
    enddo
end subroutine prefit_cent
!*****************************************************************************************
