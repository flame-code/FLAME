!*****************************************************************************************
subroutine cal_ann_tb(parini,partb,atoms,ann_arr,symfunc,opt_ann)
    use mod_parini, only: typ_parini
    use mod_tightbinding, only: typ_partb
    use mod_potl, only: potl_typ
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_ann, only: typ_ann_arr, convert_ann_epotd
    use mod_symfunc, only: typ_symfunc
    use mod_opt_ann, only: typ_opt_ann, set_opt_ann_grad
    use mod_linked_lists, only: typ_pia_arr, typ_linked_lists
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_opt_ann), intent(inout):: opt_ann
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc), intent(inout):: symfunc
    type(typ_partb), intent(inout):: partb
    !local variables
    type(typ_linked_lists):: linked_lists
    type(typ_pia_arr):: pia_arr
    type(potl_typ):: pplocal
    real(8), allocatable:: hgen(:,:), dhgen(:,:)
    real(8), allocatable:: ann_grad(:,:)
    integer:: iat, jat, ng, i, j, k, isat, ib, nb, ixyz
    real(8):: hgen_der(4,1:atoms%nat,1:atoms%nat)  , ttxyz !derivative of 
    real(8):: epotn, tt, epotdh, c, dx, dy, dz, r, rsq, hbar, fc, dfc, tt1
    real(8):: rc, rs, pi
    atoms%fat=0.d0
    partb%paircut=ann_arr%rcut
    allocate(partb%dedh(4,atoms%nat,atoms%nat),source=0.d0)
    allocate(ann_grad(ann_arr%nweight_max,ann_arr%nann))
    linked_lists%rcut=partb%paircut !ann_arr%rcut
    linked_lists%triplex=.true.
    call call_linkedlist(parini,atoms,.true.,linked_lists,pia_arr)
    allocate(partb%hgenall0(linked_lists%maxbound_rad),source=0.d0)
    allocate(partb%hgenall1(linked_lists%maxbound_rad),source=0.d0)
    allocate(partb%hgenall2(linked_lists%maxbound_rad),source=0.d0)
    allocate(partb%hgenall3(linked_lists%maxbound_rad),source=0.d0)
    allocate(partb%dhgenall0(linked_lists%maxbound_rad),source=0.d0)
    allocate(partb%dhgenall1(linked_lists%maxbound_rad),source=0.d0)
    allocate(partb%dhgenall2(linked_lists%maxbound_rad),source=0.d0)
    allocate(partb%dhgenall3(linked_lists%maxbound_rad),source=0.d0)
    if(ann_arr%compute_symfunc) then
        call symmetry_functions(parini,ann_arr,atoms,symfunc,.true.)
    endif
    !if(symfunc%linked_lists%maxbound_rad/=2) stop 'ERROR: correct next line'
    nb=symfunc%linked_lists%maxbound_rad!/2
    allocate(hgen(4,nb))
    allocate(dhgen(4,nb))
    if(trim(ann_arr%event)=='train') then
        allocate(ann_arr%g_per_bond(ann_arr%num(1),4,nb))
    endif
    over_i: do i=1,4
        over_ib: do ib=1,nb
            ng=ann_arr%ann(i)%nn(0)
            ann_arr%ann(i)%y(1:ng,0)=symfunc%y(1:ng,ib)
            if(trim(ann_arr%event)=='train') then
            !write(*,*) "symfunc", ng, nb, symfunc%y(1,ib)
                call cal_architecture_der(ann_arr%ann(i),hgen(i,ib))
                !write(*,*) "hopping", hgen(i,ib)
                call convert_ann_epotd(ann_arr%ann(i),ann_arr%num(i),ann_arr%g_per_bond(1,i,ib))
                !write(*,*) "dhda", ann_arr%g_per_bond(1,i,ib)
            elseif(trim(ann_arr%event)=='potential' .or. trim(ann_arr%event)=='evalu') then
                call cal_architecture(ann_arr%ann(i),hgen(i,ib))          
                dhgen(i,ib)=0.d0
                do j=1,ann_arr%ann(i)%nn(0)
                    dhgen(i,ib)=dhgen(i,ib) + ann_arr%ann(i)%d(j)*symfunc%y0d_bond(j,ib)
                enddo
            else
                stop 'ERROR: in cal_ann_tb undefined content for ann_arr%event'
            endif
        enddo over_ib
    enddo over_i
    pi=4.d0*atan(1.d0)
    call update_ratp(atoms)
    do ib=1,nb
        iat=symfunc%linked_lists%bound_rad(1,ib)
        jat=symfunc%linked_lists%bound_rad(2,ib)
        if(trim(ann_arr%event)=='potential') then
        dx=atoms%ratp(1,iat)-atoms%ratp(1,jat)
        dy=atoms%ratp(2,iat)-atoms%ratp(2,jat)
        dz=atoms%ratp(3,iat)-atoms%ratp(3,jat)
        rsq=dx**2+dy**2+dz**2
        r=sqrt(rsq)
        rc=partb%paircut
        rs=0.9d0*rc
        if(r<rs) then
            fc=1.d0
            dfc=0.d0
        elseif(r<rc) then
            tt=(r-rs)/(rc-rs)
            fc= cos((1.d0-(1.d0-tt**2)**3)*pi*0.5d0)
            dfc=sin((1.d0-(1.d0-tt**2)**3)*pi*0.5d0)*(-3.d0*pi*(1.d0-tt**2)**2*(r-rs)/(rc-rs)**2)
        else
            fc=0.d0
            dfc=0.d0
        endif
        do i=1,4
            hbar=hgen(i,ib)
            hgen(i,ib)=hbar*fc
            dhgen(i,ib)=dhgen(i,ib)*fc+hbar*dfc
            !do j=1, ann_arr%num(1)
            !    ann_arr%g_per_bond(j,i,ib)=fc*ann_arr%g_per_bond(j,i,ib)
            !enddo
        enddo
        endif
        
        partb%hgenall0(ib)=hgen(1,ib)
        partb%hgenall1(ib)=hgen(2,ib)
        partb%hgenall2(ib)=hgen(3,ib)
        partb%hgenall3(ib)=hgen(4,ib)
        partb%dhgenall0(ib)=dhgen(1,ib)
        partb%dhgenall1(ib)=dhgen(2,ib)
        partb%dhgenall2(ib)=dhgen(3,ib)
        partb%dhgenall3(ib)=dhgen(4,ib)
    enddo
        partb%event=ann_arr%event
        call lenoskytb_ann(parini,ann_arr,pia_arr,linked_lists,partb,atoms,atoms%nat,c)
        if(trim(ann_arr%event)=='train') then
            ann_grad=0.d0
            do i=1,4
                do j=1,ann_arr%num(1)
                    tt1=0.d0
                    do ib=1,nb
                    iat=symfunc%linked_lists%bound_rad(1,ib)
                    jat=symfunc%linked_lists%bound_rad(2,ib)
                    !write(*,'(a,5i4,2es14.5,i4)') 'DEDH',i,j,ib,iat,jat,partb%dedh(i,iat,jat),ann_arr%g_per_bond(j,i,ib),opt_ann%loc(i)+j-1
                    tt1=tt1+partb%dedh(i,iat,jat)*ann_arr%g_per_bond(j,i,ib)
                    enddo
                    ann_grad(j,i)=tt1
                enddo
            enddo
            call set_opt_ann_grad(ann_arr,ann_grad,opt_ann)
        endif
    deallocate(hgen)
    deallocate(dhgen)
    if(trim(ann_arr%event)=='train') then
        deallocate(ann_arr%g_per_bond)
        deallocate(partb%dedh)
    endif
    !if(trim(ann_arr%event)=='potential' .or. trim(parini%symfunc)=='do_not_save') then
    if(trim(ann_arr%event)=='potential') then
        deallocate(symfunc%y)
        deallocate(symfunc%y0d_bond)
        deallocate(symfunc%y0dr)
    endif
    !if(.not. (trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train' .and. \
    !    trim(parini%symfunc)/='do_not_save')) then
    if(.not. (trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train')) then
        deallocate(symfunc%linked_lists%prime_bound)
        deallocate(symfunc%linked_lists%bound_rad)
        deallocate(symfunc%linked_lists%bound_ang)
    endif
    deallocate(linked_lists%prime_bound)
    deallocate(linked_lists%bound_rad)
    deallocate(linked_lists%bound_ang)
    deallocate(ann_grad)
end subroutine cal_ann_tb
!*****************************************************************************************
subroutine lenoskytb_ann(parini,ann_arr,pia_arr,linked_lists,partb,atoms,natsi,count_md)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    use mod_potl, only: potl_typ
    use mod_linked_lists, only: typ_pia_arr, typ_linked_lists
    use mod_tightbinding, only: typ_partb, lenosky
    use mod_const, only: ha2ev, bohr2ang
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_linked_lists), intent(in):: linked_lists
    type(typ_pia_arr), intent(in):: pia_arr
    type(typ_partb), intent(inout):: partb
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in):: natsi
    real(8), intent(inout):: count_md
    !local variables
    integer, save:: icall=0, firstcall=1
    type(potl_typ), save:: pplocal
    type(typ_partb), save:: partb_pair
    real(8), allocatable:: fat(:,:)
    integer:: iat
    if(firstcall==1) then
        write(*,'(a)') 'Reading spline potential coeff.cls'
        call prmst38c(partb_pair,pplocal) !Reads potential 
        firstcall=0
    endif
    partb%usepairpot=partb_pair%usepairpot
    partb%paircut=partb_pair%paircut 

    icall=icall+1
    !write(*,'(a,f)') 'paircut= ',partb%paircut
    call lenoskytb_init(partb,atoms,natsi,linked_lists)
    count_md=count_md+1.d0
    !PRINT SOME WARNINGS
    if(natsi>atoms%nat) write(*,'(a)') 'WARNING natsi = ',natsi,' is greater than number of atoms = ',atoms%nat
    if(atoms%nat == 0) write(*,'(a)') 'WARNING lenoskytb called with zero atoms'
    if(atoms%nat < 0) write(*,'(a)') 'WARNING lenoskytb called with negative number of atoms'
    !write(*,*) 'natsi ',natsi
    call totalenergy(pia_arr,linked_lists,parini,partb,atoms,natsi,pplocal) 

    allocate(fat(3,atoms%nat))
    do iat=1,atoms%nat
        fat(1,iat)=atoms%fat(1,iat)
        fat(2,iat)=atoms%fat(2,iat)
        fat(3,iat)=atoms%fat(3,iat)
        atoms%fat(1,iat)=0.d0
        atoms%fat(2,iat)=0.d0
        atoms%fat(3,iat)=0.d0
    enddo

    call pairenergy(parini,partb,atoms,pplocal,natsi)
    atoms%epot=atoms%epot+partb%pairen !+atoms%nat*-0.789592525650303d+04/ha2ev
    ann_arr%ener_ref=0.d0
    do iat=1,atoms%nat
        ann_arr%ener_ref=ann_arr%ener_ref+ann_arr%ann(atoms%itypat(iat))%ener_ref
    enddo
    atoms%epot=atoms%epot+ann_arr%ener_ref

    !write(61,*) atoms%epot
    !stop 'BBBBBBBBB'

    do iat=1,atoms%nat
        atoms%fat(1,iat)=fat(1,iat)+atoms%fat(1,iat)
        atoms%fat(2,iat)=fat(2,iat)+atoms%fat(2,iat)
        atoms%fat(3,iat)=fat(3,iat)+atoms%fat(3,iat)
    enddo
    deallocate(fat)

    call lenoskytb_final(partb)
end subroutine lenoskytb_ann
!*****************************************************************************************
subroutine fit_hgen(parini,ann_arr,opt_ann)
    use mod_parini, only: typ_parini
    !use mod_tightbinding, only: typ_partb
    !use mod_potl, only: potl_typ
    use mod_atoms, only: typ_atoms, atom_allocate_old, set_rat_iat, update_ratp
    use mod_ann, only: typ_ann_arr, convert_x_ann
    use mod_symfunc, only: typ_symfunc
    use mod_opt_ann, only: typ_opt_ann, get_opt_ann_x, set_opt_ann_x
    use mod_parlm, only: typ_parlm
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    !type(typ_atoms), intent(inout):: atoms
    type(typ_opt_ann), intent(inout):: opt_ann
    type(typ_ann_arr), intent(inout):: ann_arr
    !type(typ_symfunc), intent(inout):: symfunc
    !type(typ_partb), intent(inout):: partb
    !local variables
    type(typ_atoms):: atoms
    type(typ_symfunc):: symfunc
    integer:: nb, i, ib, ng, iann, ios, ii, m
    real(8):: fc, dfc, rsq, r, rmse_hgen(4), alpha_arr(4), gnrm_arr(4)
    type(typ_parlm):: parlm
    real(8), allocatable:: yall(:,:,:)
    real(8), allocatable:: grad(:,:)
    real(8):: hgen_ltb(4,325), dis_ltb(325)
    real(8), allocatable:: x_arr(:,:)
    character(10):: str_tmp
    call atom_allocate_old(atoms,2,0,0)
    atoms%boundcond='free'
    atoms%cellvec=0.d0
    atoms%cellvec(1,1)=30.d0
    atoms%cellvec(2,2)=30.d0
    atoms%cellvec(3,3)=30.d0
    call set_rat_iat(atoms,1,(/5.d0,5.d0,5.d0/))
    call set_rat_iat(atoms,2,(/7.d0,5.d0,5.d0/))

    open(unit=101,file='hgen.ltb',status='old',iostat=ios)
    if(ios/=0) then;write(*,'(a)') 'ERROR: failure openning hgen.ltb ';stop;endif
    do i=1,325
        read(101,*) str_tmp,hgen_ltb(1,i),hgen_ltb(2,i),hgen_ltb(3,i),hgen_ltb(4,i),dis_ltb(i)
        !write(44,'(a,5es14.5)') 'hgen-L',hgen_ltb(1,i),hgen_ltb(2,i),hgen_ltb(3,i),hgen_ltb(4,i),dis_ltb(i)
    enddo
    close(101)
    allocate(grad(ann_arr%num(1),4))
    nb=1
    allocate(ann_arr%g_per_bond(ann_arr%num(1),4,nb))
    allocate(yall(ann_arr%ann(1)%nn(0),1000,325))
    call update_ratp(atoms)
    do i=1,325
        !call atom_copy_old(atoms_train%atoms(iconf),atoms,'atoms_train%atoms(iconf)->atoms')
        atoms%ratp(1,2)=atoms%ratp(1,1)+dis_ltb(i)
        call symmetry_functions(parini,ann_arr,atoms,symfunc,.true.)
        ng=ann_arr%ann(1)%nn(0)
        nb=symfunc%linked_lists%maxbound_rad
        yall(1:ng,1:nb,i)=symfunc%y(1:ng,1:nb)
        deallocate(symfunc%linked_lists%prime_bound)
        deallocate(symfunc%linked_lists%bound_rad)
        deallocate(symfunc%linked_lists%bound_ang)
        deallocate(symfunc%y)
        deallocate(symfunc%y0d)
        deallocate(symfunc%y0d_bond)
        deallocate(symfunc%y0dr)
    enddo
    !-------------------------------------------------------
    parlm%xtol=5.d-3
    parlm%ftol=5.d-3
    parlm%gtol=5.d-3
    m=325
    parlm%n=opt_ann%n/4
    x_arr=f_malloc([1.to.ann_arr%nweight_max,1.to.ann_arr%nann],id='x_arr')
    do iann=1,4
    call init_lmder_modified(parlm,m,m)
    call get_opt_ann_x(ann_arr,opt_ann,x_arr)
    !parlm%x(1:parlm%n)=opt_ann%x(opt_ann%loc(iann):opt_ann%loc(iann)+ann_arr%num(1)-1)
    parlm%x(1:parlm%n)=x_arr(1:parlm%n,iann)
    do
        call lmder_modified(parlm,m,m)
        write(*,*) 'nfev,njev: ',parlm%nfev,parlm%njev
        if(parlm%finish) exit
        !do iann=1,ann_arr%nann
        !    call convert_x_ann(ann_arr%num(iann),parlm%x,ann_arr%ann(iann))
        !enddo
        if(parlm%icontinue==700) then
            call convert_x_ann(ann_arr%num(iann),parlm%wa2,ann_arr%ann(iann))
            call fcn_hgen(m,parlm%n,parlm%wa2,parlm%wa4,parlm%fjac,m,parlm%iflag,iann,ann_arr,hgen_ltb,yall)
        else
            call convert_x_ann(ann_arr%num(iann),parlm%x,ann_arr%ann(iann))
            call fcn_hgen(m,parlm%n,parlm%x,parlm%fvec,parlm%fjac,m,parlm%iflag,iann,ann_arr,hgen_ltb,yall)
        endif
    enddo
    !opt_ann%x(opt_ann%loc(iann):opt_ann%loc(iann)+ann_arr%num(1)-1)=parlm%x(1:parlm%n)
    x_arr(1:parlm%n,iann)=parlm%x(1:parlm%n)
    call set_opt_ann_x(ann_arr,x_arr,opt_ann)
    write(*,*) 'info= ',parlm%info
    call final_lmder_modified(parlm)
    enddo
    !-------------------------------------------------------
    call f_free(x_arr)
    deallocate(grad)
    deallocate(ann_arr%g_per_bond)
end subroutine fit_hgen
!*****************************************************************************************
subroutine fcn_hgen(m,n,x,fvec,fjac,ldfjac,iflag,iann,ann_arr,hgen_ltb,yall)
    use mod_ann, only: typ_ann_arr, convert_ann_epotd
    use mod_symfunc, only: typ_symfunc
    implicit none
    integer, intent(in):: m, n, ldfjac, iflag, iann
    type(typ_ann_arr), intent(inout):: ann_arr
    real(8), intent(in):: x(n), hgen_ltb(4,325), yall(ann_arr%ann(1)%nn(0),1000,325)
    real(8), intent(inout):: fvec(m), fjac(m,n)
    !local variables
    integer:: i, ib, ii, nb, ng
    real(8):: res, rmse
    real(8), allocatable:: hgen(:,:), dhgen(:,:)
    nb=1
    allocate(hgen(4,nb))
    allocate(dhgen(4,nb))
    if(iflag==1) then
        do i=1,m
            nb=1 !symfunc%linked_lists%maxbound_rad!/2
            do ib=1,nb
                ng=ann_arr%ann(iann)%nn(0)
                ann_arr%ann(iann)%y(1:ng,0)=yall(1:ng,ib,i) !symfunc%y(1:ng,ib)
                call cal_architecture_der(ann_arr%ann(iann),hgen(iann,ib))
                write(*,*) 'HGEN ',hgen(iann,ib)
                call convert_ann_epotd(ann_arr%ann(iann),n,ann_arr%g_per_bond(1,iann,ib))
                res=hgen(iann,ib)-hgen_ltb(iann,i)
                fvec(i)=res**2
            enddo
        enddo
        !do i=1,m
        !    call exppp(n,x,pnt(i),f,g)
        !    res=f-f_p(i)
        !    fvec(i)=res**2
        !enddo
    elseif(iflag==2) then
        do i=1,m
            nb=1 !symfunc%linked_lists%maxbound_rad!/2
            do ib=1,nb
                ng=ann_arr%ann(iann)%nn(0)
                ann_arr%ann(iann)%y(1:ng,0)=yall(1:ng,ib,i) !symfunc%y(1:ng,ib)
                call cal_architecture_der(ann_arr%ann(iann),hgen(iann,ib))
                call convert_ann_epotd(ann_arr%ann(iann),n,ann_arr%g_per_bond(1,iann,ib))
                !rmse_hgen(iann)=rmse_hgen(iann)+(hgen(iann,ib)-hgen_ltb(iann,i))**2
                res=hgen(iann,ib)-hgen_ltb(iann,i)
                do ii=1,n
                    fjac(i,ii)=2.d0*res*ann_arr%g_per_bond(ii,iann,ib)
                enddo
            enddo
        enddo
        !do i=1,m
        !    call exppp(n,x,pnt(i),f,g)
        !    res=f-f_p(i)
        !    fjac(i,1:3)=g(1:3)*res*2.d0
        !enddo
    elseif(iflag==0) then
        rmse=0.d0
        do i=1,m
            nb=1 !symfunc%linked_lists%maxbound_rad!/2
            do ib=1,nb
                ng=ann_arr%ann(iann)%nn(0)
                ann_arr%ann(iann)%y(1:ng,0)=yall(1:ng,ib,i) !symfunc%y(1:ng,ib)
                call cal_architecture_der(ann_arr%ann(iann),hgen(iann,ib))
                call convert_ann_epotd(ann_arr%ann(iann),n,ann_arr%g_per_bond(1,iann,ib))
                res=hgen(iann,ib)-hgen_ltb(iann,i)
                rmse=rmse+res**2
            enddo
        enddo
        rmse=sqrt(rmse/m)
        write(*,*) 'RMSE: ',rmse
    endif
    deallocate(hgen)
    deallocate(dhgen)
end subroutine fcn_hgen
!*****************************************************************************************
