!*****************************************************************************************
subroutine get_fcn_ann(parini,idp,str_dataset,ann_arr,opt_ann,fcn_ann,fcn_ref)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_opt_ann, only: typ_opt_ann, set_opt_ann_grad
    use mod_atoms, only: typ_atoms, atom_copy_old, update_ratp
    use mod_callback_ann, only: atoms_train=>atoms_train_t
    use mod_callback_ann, only: symfunc_train=>symfunc_train_t
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: idp
    character(*), intent(in):: str_dataset
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_opt_ann), intent(inout):: opt_ann
    real(8), intent(out):: fcn_ann
    real(8), intent(out):: fcn_ref
    !local variables
    type(typ_atoms):: atoms
    real(8), allocatable:: ann_grad(:,:)
    integer:: iat, i, j, iconf, ixyz
    !-----------------------------------------------------------------
    if(trim(ann_arr%approach)=='atombased') then
        iconf=idp
    elseif(trim(ann_arr%approach)=='cent1') then
        iconf=idp
    elseif(trim(ann_arr%approach)=='cent2') then
        iconf=idp
    elseif(trim(ann_arr%approach)=='centt') then
        iconf=idp
    elseif(trim(ann_arr%approach)=='cent3') then
        iconf=int((idp-1)/3)+1
        ixyz=mod(idp-1,3)+1
    elseif(trim(ann_arr%approach)=='tb') then
        iconf=idp
    endif
    !-----------------------------------------------------------------
    call atom_copy_old(atoms_train%atoms(iconf),atoms,'atoms_train%atoms(iconf)->atoms')
    if(trim(ann_arr%approach)=='cent2') then
        if (allocated(ann_arr%a)) deallocate(ann_arr%a)
        if(.not. allocated(ann_arr%ann_amat_train(iconf)%amat)) then 
            allocate(ann_arr%a(1:(atoms%nat+1)*(atoms%nat+1)))
            allocate(ann_arr%ann_amat_train(iconf)%amat(1:(atoms%nat+1)*(atoms%nat+1)))
            ann_arr%ann_amat_train(iconf)%amat=0.d0
            ann_arr%a=0.d0
            ann_arr%amat_initiated=.false.
        else
            allocate(ann_arr%a(1:(atoms%nat+1)*(atoms%nat+1)))
            ann_arr%a(1:(atoms%nat+1)*(atoms%nat+1))=ann_arr%ann_amat_train(iconf)%amat(1:(atoms%nat+1)*(atoms%nat+1))
            ann_arr%amat_initiated=.true.
        end if 
        if (allocated(ann_arr%Xq)) deallocate(ann_arr%Xq)
        if(.not. allocated(ann_arr%ann_chiQPar_train(iconf)%chiQPar)) then 
            allocate(ann_arr%Xq((atoms%nat),(atoms%nat)))
            allocate(ann_arr%ann_chiQPar_train(iconf)%chiQPar(1:(atoms%nat),1:(atoms%nat)))
            ann_arr%ann_chiQPar_train(iconf)%chiQPar=0.d0
            ann_arr%Xq=0.d0
            ann_arr%chiQPar_initiated=.false.
        else
            allocate(ann_arr%Xq(1:(atoms%nat),1:(atoms%nat)))
            ann_arr%Xq(1:(atoms%nat),1:(atoms%nat))=ann_arr%ann_chiQPar_train(iconf)%chiQPar(1:(atoms%nat),1:(atoms%nat))
            ann_arr%chiQPar_initiated=.true.
        end if 
        if (allocated(ann_arr%EP)) deallocate(ann_arr%EP)
        if(.not. allocated(ann_arr%ann_EPar_train(iconf)%EPar)) then 
            allocate(ann_arr%EP(1:atoms%nat))
            allocate(ann_arr%ann_EPar_train(iconf)%EPar(1:(atoms%nat)))
            ann_arr%ann_EPar_train(iconf)%EPar=0.d0
            ann_arr%EP=0.d0
            ann_arr%EPar_initiated=.false.
        else
            allocate(ann_arr%EP(1:(atoms%nat)))
            ann_arr%EP(1:(atoms%nat))=ann_arr%ann_EPar_train(iconf)%EPar(1:(atoms%nat))
            ann_arr%EPar_initiated=.true.
        end if 
    end if 
    call cal_ann_main(parini,atoms,symfunc_train%symfunc(iconf),ann_arr,opt_ann)
    if(trim(ann_arr%approach)=='cent2') then
        if(.not. ann_arr%amat_initiated) then
            ann_arr%ann_amat_train(iconf)%amat=ann_arr%a
        end if
        deallocate(ann_arr%a)
        if(.not. ann_arr%chiQPar_initiated) then
            ann_arr%ann_chiQPar_train(iconf)%chiQPar=ann_arr%Xq
        end if
        deallocate(ann_arr%Xq)
        if(.not. ann_arr%EPar_initiated) then
            ann_arr%ann_EPar_train(iconf)%EPar=ann_arr%EP
        end if
        deallocate(ann_arr%EP)
    end if
    !-----------------------------------------------------------------
    allocate(ann_grad(ann_arr%nweight_max,ann_arr%nann),source=0.d0)
    if(trim(ann_arr%approach)=='atombased') then
        do iat=1,atoms%nat
            i=atoms%itypat(iat)
            do j=1,ann_arr%nweight_max
                ann_grad(j,i)=ann_grad(j,i)+ann_arr%g_per_atom(j,iat)
            enddo
        enddo
        call set_opt_ann_grad(ann_arr,ann_grad,opt_ann)
    elseif(trim(ann_arr%approach)=='cent1') then
        do iat=1,atoms%nat
            i=atoms%itypat(iat)
            do j=1,ann_arr%nweight_max
                ann_grad(j,i)=ann_grad(j,i)+atoms%qat(iat)*ann_arr%g_per_atom(j,iat)
            enddo
        enddo
        call set_opt_ann_grad(ann_arr,ann_grad,opt_ann)
    elseif(trim(ann_arr%approach)=='cent2') then
        do iat=1,atoms%nat
            i=atoms%itypat(iat)
            do j=1,ann_arr%nweight_max
                !ann_grad(j,i)=ann_grad(j,i)+atoms%qat(iat)*ann_arr%g_per_atom(j,iat)
                ann_grad(j,i)=ann_grad(j,i)+ann_arr%g_per_atom(j,iat)
            enddo
        enddo
        call set_opt_ann_grad(ann_arr,ann_grad,opt_ann)
    elseif(trim(ann_arr%approach)=='centt') then
        do iat=1,atoms%nat
            i=atoms%itypat(iat)
            do j=1,ann_arr%nweight_max
                ann_grad(j,i)=ann_grad(j,i)+(atoms%zat(iat)+atoms%qat(iat))*ann_arr%g_per_atom(j,iat)
            enddo
        enddo
        call set_opt_ann_grad(ann_arr,ann_grad,opt_ann)
    elseif(trim(ann_arr%approach)=='cent3') then
        call update_ratp(atoms)
        do iat=1,atoms%nat
            i=atoms%itypat(iat)
            do j=1,ann_arr%nweight_max
                ann_grad(j,i)=ann_grad(j,i)+(atoms%ratp(ixyz,iat))*ann_arr%dqat_weights(j,iat)
            enddo
        enddo
        call set_opt_ann_grad(ann_arr,ann_grad,opt_ann)
    endif
    deallocate(ann_grad)
    !-----------------------------------------------------------------
    if(trim(ann_arr%approach)=='cent3') then
        fcn_ann=atoms%dpm(ixyz)
        fcn_ref=atoms_train%atoms(iconf)%dpm(ixyz)
        write(*,'(a,2f10.3)') 'fcn_ann,fcn_ref ',fcn_ann,fcn_ref
    else
        fcn_ann=atoms%epot
        fcn_ref=atoms_train%atoms(iconf)%epot
    endif
end subroutine get_fcn_ann
!*****************************************************************************************
subroutine cal_ann_main(parini,atoms,symfunc,ann_arr,opt_ann)
    use mod_tightbinding, only: typ_partb
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc
    use mod_opt_ann, only: typ_opt_ann
    !use mod_tightbinding, only: typ_partb
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc), intent(inout):: symfunc
    type(typ_opt_ann), intent(inout):: opt_ann
    !local variables
    !real(8):: g, g_tb, dis, E0, E1 
    !real(8), allocatable:: xt(:), gt(:)
    integer:: i, j, iat
    type(typ_partb):: partb
    if(trim(ann_arr%approach)=='atombased') then
        call cal_ann_atombased(parini,atoms,symfunc,ann_arr)
    elseif(trim(ann_arr%approach)=='eem1' .or. trim(ann_arr%approach)=='cent1') then
        call cal_ann_cent1(parini,atoms,symfunc,ann_arr)
    elseif(trim(ann_arr%approach)=='cent2') then
        call cal_ann_cent2(parini,atoms,symfunc,ann_arr)
    elseif(trim(ann_arr%approach)=='centt') then
        call cal_ann_centt(parini,atoms,symfunc,ann_arr)
    elseif(trim(ann_arr%approach)=='cent3') then
        call cal_ann_cent3(parini,atoms,symfunc,ann_arr)
    elseif(trim(ann_arr%approach)=='tb') then
        call cal_ann_tb(parini,partb,atoms,ann_arr,symfunc,opt_ann)
        !if(trim(ann_arr%event)=='train') then
        ! E0=atoms%epot
        ! allocate(xt(opt_ann%n),gt(opt_ann%n))
        ! xt(1:opt_ann%n)=opt_ann%x(1:opt_ann%n)
        ! gt(1:opt_ann%n)=ann_grad(1:opt_ann%n)
        ! do i=1,opt_ann%n
        !     g_tb=gt(i)
        !     !!Finite difference 
        !     dis=1.d-4 !*abs(opt_ann%x(i))
        !     opt_ann%x(i)=opt_ann%x(i)+dis
        !     call cal_ann_tb(parini,partb,atoms,ann_arr,symfunc,opt_ann)
        !     E1=atoms%epot
        !     g=(E1-E0)/dis
        !     write(*,'(a,2es19.10,es14.5,2es19.10)') 'FD-TEST',g_tb,g,g-g_tb,E0,E1
        !     opt_ann%x(i)=xt(i)
        ! enddo
        ! stop 'TTTTTTTTTTTTTTTT'
        !endif
    else
        write(*,'(2a)') 'ERROR: unknown approach in ANN, ',trim(ann_arr%approach)
        stop
    endif
end subroutine cal_ann_main
!*****************************************************************************************
subroutine prefit_cent_ener_ref(parini,ann_arr,symfunc_train,symfunc_valid,atoms_train,atoms_valid,opt_ann)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc_arr
    use mod_opt_ann, only: typ_opt_ann, convert_opt_x_ann_arr
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_copy_old
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc_arr), intent(inout):: symfunc_train, symfunc_valid
    type(typ_atoms_arr), intent(inout):: atoms_train
    type(typ_atoms_arr), intent(inout):: atoms_valid
    type(typ_opt_ann), intent(inout):: opt_ann
    !local variables
    type(typ_atoms):: atoms
    integer:: iconf, istep, iat, ia, isatur, nsatur
    real(8):: anat(100), g(100), rmse, rmse_old, de0, alpha, tt
    real(8), allocatable:: epotall(:), eref_all(:)
    !return
    ann_arr%event='train'
    call convert_opt_x_ann_arr(opt_ann,ann_arr)
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
                call cal_ann_main(parini,atoms,symfunc_train%symfunc(iconf),ann_arr,opt_ann)
                epotall(iconf)=atoms%epot
            else
                tt=0.d0
                do iat=1,atoms%nat
                    tt=tt+ann_arr%ann(atoms%itypat(iat))%ener_ref
                enddo
                atoms%epot=epotall(iconf)-eref_all(iconf)+tt
            endif
            rmse=rmse+((atoms_train%atoms(iconf)%epot-atoms%epot)/atoms%nat)**2
            anat=0.d0
            do iat=1,atoms%nat
                anat(atoms%itypat(iat))=anat(atoms%itypat(iat))+1.d0
            enddo
            do ia=1,ann_arr%nann
                g(ia)=g(ia)+2.d0*anat(ia)*(atoms%epot-atoms_train%atoms(iconf)%epot)/atoms%nat**2
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
        do ia=1,ann_arr%nann
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
subroutine prefit_cent(parini,ann_arr,symfunc_train,symfunc_valid,atoms_train,atoms_valid,opt_ann)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc_arr
    use mod_opt_ann, only: typ_opt_ann, convert_opt_x_ann_arr
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_copy_old
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc_arr), intent(inout):: symfunc_train, symfunc_valid
    type(typ_atoms_arr), intent(inout):: atoms_train
    type(typ_atoms_arr), intent(inout):: atoms_valid
    type(typ_opt_ann), intent(inout):: opt_ann
    !local variables
    type(typ_atoms):: atoms
    integer:: iconf, istep, iat, ia, isatur, nsatur
    real(8):: anat1(100), g1(100), rmse, rmse_old, dchi0, dhardness, alpha1, alpha2, tt
    real(8):: anat2(100), g2(100), qnet
    ann_arr%event='train'
    call convert_opt_x_ann_arr(opt_ann,ann_arr)
    nsatur=3
    isatur=0
    alpha1=0.2d0/real(atoms_train%nconf,8)
    alpha2=0.02d0/real(atoms_train%nconf,8)
    do istep=0,50
        rmse=0.d0
        g1=0.d0
        g2=0.d0
        do iconf=1,atoms_train%nconf
            call atom_copy_old(atoms_train%atoms(iconf),atoms,'atoms_train%atoms(iconf)->atoms')
            call cal_ann_main(parini,atoms,symfunc_train%symfunc(iconf),ann_arr,opt_ann)
            rmse=rmse+((atoms_train%atoms(iconf)%epot-atoms%epot)/atoms%nat)**2
            anat1=0.d0
            anat2=0.d0
            do iat=1,atoms%nat
                if(trim(ann_arr%approach)=='eem1' .or. trim(ann_arr%approach)=='cent1' &
                    .or. trim(ann_arr%approach)=='cent2') then
                    qnet=atoms%qat(iat)
                elseif(trim(ann_arr%approach)=='centt') then
                    qnet=atoms%zat(iat)+atoms%qat(iat)
                else
                    write(*,'(2a)') 'ERROR: unknown approach in ANN, ',trim(ann_arr%approach)
                    stop
                endif
                anat1(atoms%itypat(iat))=anat1(atoms%itypat(iat))+qnet
                anat2(atoms%itypat(iat))=anat2(atoms%itypat(iat))+0.5d0*qnet**2
            enddo
            do ia=1,ann_arr%nann
                g1(ia)=g1(ia)+2.d0*anat1(ia)*(atoms%epot-atoms_train%atoms(iconf)%epot)/atoms%nat**2
                g2(ia)=g2(ia)+2.d0*anat2(ia)*(atoms%epot-atoms_train%atoms(iconf)%epot)/atoms%nat**2
            enddo
        enddo
        rmse=sqrt(rmse/atoms_train%nconf)
        if(istep==0) rmse_old=rmse
        if(istep>0 .and. rmse<rmse_old .and. abs(rmse_old-rmse)<1.d-4) then
            isatur=isatur+1
        else
            isatur=0
        endif
        write(*,'(a,i4,5es19.10,i3)') 'prefit: ',istep,rmse*1.d3, &
            ann_arr%ann(1)%chi0,ann_arr%ann(2)%chi0, &
            ann_arr%ann(1)%hardness,ann_arr%ann(2)%hardness,isatur
        if(rmse*1.d3<1.d0) exit
        if(isatur>nsatur) exit
        do ia=1,ann_arr%nann
            dchi0=alpha1*g1(ia)
            dchi0=sign(min(abs(dchi0),1.d-1),dchi0)
            ann_arr%ann(ia)%chi0=ann_arr%ann(ia)%chi0-dchi0
            dhardness=alpha2*g2(ia)
            dhardness=sign(min(abs(dhardness),1.d-1),dhardness)
            ann_arr%ann(ia)%hardness=ann_arr%ann(ia)%hardness-dhardness
        enddo
        rmse_old=rmse
    enddo
end subroutine prefit_cent
!*****************************************************************************************
