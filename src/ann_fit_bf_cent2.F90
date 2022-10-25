!*****************************************************************************************
module mod_fit_bf_cent2
    implicit none
    private
    public:: get_basis_functions_cent2
contains
!*****************************************************************************************
subroutine get_basis_functions_cent2(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp, typ_atoms_arr, atom_deallocate_old
    !use mod_ann, only: typ_ann_arr
    use mod_ann_io_yaml, only: read_data_yaml
    use mod_electrostatics, only: typ_poisson
    use mod_linkedlists, only: typ_linkedlists
    use mod_linked_lists, only: typ_pia_arr
    use mod_flm_futile
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    type(typ_poisson):: poisson, poisson_ref
    type(typ_atoms):: atoms
    type(typ_atoms_arr):: atoms_arr
    !type(typ_ann_arr):: ann_arr
    integer:: iat, iconf, istep
    integer:: ix, iy, iz, itypat
    !real(8):: ehartree_scn_excl
    !real(8):: dpx, dpy, dpz
    !real(8):: time1, time2, time3, time4, time5, time6, time7, time8, time9
    real(8):: tt0, tt1, tt2, tt3, tt4, ts1, ts2, tp1, tp2, cost
    real(8):: gausswidth(1), q_tmp(1), p_tmp(3), voxel
    real(8):: x, y, z !, r
    real(8):: errmax, rmse, alpha
    real(8):: gg1, gg2, gg3, gg4, gg5, gg6, alpha_gg
    real(8):: bc_s1, bc_s2, bc_p1, bc_p2
    real(8):: bcg_s1(2), bcg_s2(2), bcg_p1(2), bcg_p2(2), gradtot
    real(8), allocatable:: c_s(:)
    real(8), allocatable:: c_p(:)
    real(8), allocatable:: cgrad_s(:)
    real(8), allocatable:: cgrad_p(:)
    real(8), allocatable:: gw_s1(:)
    real(8), allocatable:: gw_s2(:)
    real(8), allocatable:: gw_p1(:)
    real(8), allocatable:: gw_p2(:)
    real(8), allocatable:: gwgrad_s1(:)
    real(8), allocatable:: gwgrad_s2(:)
    real(8), allocatable:: gwgrad_p1(:)
    real(8), allocatable:: gwgrad_p2(:)
    real(8), allocatable:: vgrad_s1_t(:,:,:)
    real(8), allocatable:: vgrad_s2_t(:,:,:)
    real(8), allocatable:: vgrad_p1_t(:,:,:)
    real(8), allocatable:: vgrad_p2_t(:,:,:)
    real(8), allocatable:: vgrad_s1(:,:,:,:)
    real(8), allocatable:: vgrad_s2(:,:,:,:)
    real(8), allocatable:: vgrad_p1(:,:,:,:)
    real(8), allocatable:: vgrad_p2(:,:,:,:)
    real(8), allocatable:: pot_s1(:,:,:)
    real(8), allocatable:: pot_s2(:,:,:)
    real(8), allocatable:: pot_p1(:,:,:)
    real(8), allocatable:: pot_p2(:,:,:)
    real(8), allocatable:: pot_iat_s(:,:,:,:)
    real(8), allocatable:: pot_iat_p(:,:,:,:)
    real(8), allocatable:: pot_ion(:,:,:)
    call f_routine(id='get_basis_functions_cent2')
    !pi=4.d0*atan(1.d0)
    call read_data_yaml(parini,'list_posinp_cent2.yaml',atoms_arr)
    call get_poisson_ref(parini,poisson_ref,atoms)
    do iat=1,atoms%nat
        write(*,*) atoms%ratp(1,iat),atoms%ratp(2,iat),atoms%ratp(3,iat)
    enddo
    write(*,*) poisson_ref%xyz111(1),poisson_ref%xyz111(2),poisson_ref%xyz111(3)




    !-------------------------------------------------------
    poisson%ngpx=poisson_ref%ngpx
    poisson%ngpy=poisson_ref%ngpy
    poisson%ngpz=poisson_ref%ngpz
    poisson%hgrid(1:3,1:3)=0.d0
    poisson%hgrid(1,1)=poisson_ref%hgrid(1,1)
    poisson%hgrid(2,2)=poisson_ref%hgrid(2,2)
    poisson%hgrid(3,3)=poisson_ref%hgrid(3,3)
    poisson%xyz111=poisson_ref%xyz111
    poisson%cal_scn=poisson_ref%cal_scn
    poisson%screening_factor=poisson_ref%screening_factor
    poisson%bc=poisson_ref%bc
    poisson%task_finit="alloc_rho"
    call init_hartree(parini,atoms,poisson,gausswidth)

    allocate(c_s(atoms%nat))
    allocate(c_p(atoms%nat))
    allocate(cgrad_s(atoms%nat))
    allocate(cgrad_p(atoms%nat))
    allocate(gw_s1(parini%ntypat))
    allocate(gw_s2(parini%ntypat))
    allocate(gw_p1(parini%ntypat))
    allocate(gw_p2(parini%ntypat))
    allocate(gwgrad_s1(parini%ntypat))
    allocate(gwgrad_s2(parini%ntypat))
    allocate(gwgrad_p1(parini%ntypat))
    allocate(gwgrad_p2(parini%ntypat))
    allocate(vgrad_s1_t(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    allocate(vgrad_s2_t(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    allocate(vgrad_p1_t(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    allocate(vgrad_p2_t(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    allocate(vgrad_s1(poisson%ngpx,poisson%ngpy,poisson%ngpz,atoms%nat))
    allocate(vgrad_s2(poisson%ngpx,poisson%ngpy,poisson%ngpz,atoms%nat))
    allocate(vgrad_p1(poisson%ngpx,poisson%ngpy,poisson%ngpz,atoms%nat))
    allocate(vgrad_p2(poisson%ngpx,poisson%ngpy,poisson%ngpz,atoms%nat))
    allocate(pot_s1(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    allocate(pot_s2(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    allocate(pot_p1(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    allocate(pot_p2(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    allocate(pot_iat_s(poisson%ngpx,poisson%ngpy,poisson%ngpz,parini%ntypat))
    allocate(pot_iat_p(poisson%ngpx,poisson%ngpy,poisson%ngpz,parini%ntypat))
    allocate(pot_ion(poisson%ngpx,poisson%ngpy,poisson%ngpz))

    do iat=1,atoms%nat
        do itypat=1,parini%ntypat
            if(trim(atoms%sat(iat))==trim(parini%stypat(itypat))) then
                atoms%itypat(iat)=parini%ltypat(itypat)
                exit
            endif
        enddo
    enddo

    gausswidth=0.5d0

    !q_tmp=2.d0
    !call put_gto_sym_ortho(parini,poisson%bc,.true.,1,atoms%ratp(1,1),q_tmp,gausswidth, &
    !    6.d0*maxval(gausswidth),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
    !    poisson%hgrid,poisson%rho)

    open(unit=2134,file="gw.inp",status='old')
    read(2134,*) gg1,gg2,gg3,alpha_gg
    read(2134,*) gg4,gg5,gg6,alpha
    close(2134)
    !gausswidth=0.5d0
    pot_ion=0.d0
    do iat=1,atoms%nat
    itypat=atoms%itypat(iat)
!    if(trim(atoms%sat(iat))=='Mg') then
!        q_tmp(1)=2.d0
!    elseif(trim(atoms%sat(iat))=='O') then
!        q_tmp(1)=6.d0
!    else
!        stop 'ERROR: unknown type in get_basis_functions_cent2'
!    endif
!    call cal_pot_gauss_s(parini,poisson,pot_ion,.false.,atoms%ratp(1,iat),gausswidth(1),q_tmp(1),vgrad_s1_t)
    if(trim(atoms%sat(iat))=='O') then
    gw_s1(itypat)=gg1 !1.193d0
    gw_s2(itypat)=gg2 !0.613d0
    bc_s1= gw_s1(itypat)**3/(gw_s1(itypat)**3-gw_s2(itypat)**3)
    bc_s2=-gw_s2(itypat)**3/(gw_s1(itypat)**3-gw_s2(itypat)**3)
    q_tmp(1)=gg3*bc_s1
    call cal_pot_gauss_s(parini,poisson,pot_ion,.false.,atoms%ratp(1,iat),gw_s1(itypat),q_tmp(1),vgrad_s1_t)
    q_tmp(1)=gg3*bc_s2
    call cal_pot_gauss_s(parini,poisson,pot_ion,.false.,atoms%ratp(1,iat),gw_s2(itypat),q_tmp(1),vgrad_s1_t)
    endif
    enddo

    do itypat=1,parini%ntypat
        if(trim(parini%stypat(itypat))=='Mg') then
            gw_s1(itypat)=2.598d0
            gw_s2(itypat)=1.254d0
            gw_p1(itypat)=2.719d0
            gw_p2(itypat)=1.143d0
            !2.683   1.092   3.012   1.117
            !2.631   1.130   2.746   1.110
            !2.598   1.254   2.719   1.143
        elseif(trim(parini%stypat(itypat))=='O') then
            gw_s1(itypat)=gg4
            gw_s2(itypat)=gg5
            gw_p1(itypat)=1.826d0
            gw_p2(itypat)=0.982d0
            !1.707   1.217   1.878   0.995
            !1.552   1.418   1.834   0.984
            !1.459   1.273   1.826   0.982
        else
            stop 'ERROR: unknown type in get_basis_functions_cent2'
        endif
    enddo

    !q_tmp=-2.d0*bc_s1
    !call put_gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,1),q_tmp,gw_s1, &
    !    6.d0*gw_s1,poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
    !    poisson%hgrid,poisson%rho)
    !q_tmp=-2.d0*bc_s2
    !call put_gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,1),q_tmp,gw_s2, &
    !    6.d0*gw_s2,poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
    !    poisson%hgrid,poisson%rho)
    !q_tmp=-2.d0*b3
    !call put_gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,1),q_tmp,a3, &
    !    6.d0*a3,poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
    !    poisson%hgrid,poisson%rho)

    !p_tmp=0.d0
    !p_tmp(1)=0.3634384d0 !1.d0
    !call put_gto_p_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,1),p_tmp,gausswidth, &
    !    6.d0*maxval(gausswidth),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)

    !call get_hartree(parini,poisson,atoms,gausswidth,ehartree_scn_excl)
    
!    q_tmp(1)=1.d0
!    call cal_pot_gauss_s(parini,poisson,.true.,atoms%ratp(1,1),gw_s1,q_tmp(1),vgrad_s1)
!    pot_s1=poisson%pot
!    gw_s1=gw_s1+1.d-2
!    call cal_pot_gauss_s(parini,poisson,.true.,atoms%ratp(1,1),gw_s1,q_tmp(1),vgrad_s2)
!    pot_s2=poisson%pot
!    !do iz=1,poisson%ngpz
!    do iz=poisson%ngpz/2,poisson%ngpz/2
!    do iy=poisson%ngpy/2,poisson%ngpy/2
!    !do ix=poisson%ngpx/2,poisson%ngpx/2
!    do ix=1,poisson%ngpx
!        x=poisson%xyz111(1)+(ix-1)*poisson%hgrid(1,1)
!        y=poisson%xyz111(2)+(iy-1)*poisson%hgrid(2,2)
!        z=poisson%xyz111(3)+(iz-1)*poisson%hgrid(3,3)
!        write(83,'(3f8.3,2es19.10)') x,y,z,(pot_s2(ix,iy,iz)-pot_s1(ix,iy,iz))/1.d-2,vgrad_s1(ix,iy,iz)
!    enddo
!    enddo
!    enddo
!    stop 'TTTTTTTTTTTTTTTTTTTTT'
    
!    p_tmp=0.d0
!    p_tmp(1)=1.d0
!    call cal_pot_gauss_p(parini,poisson,.true.,atoms%ratp(1,1),gw_s1,p_tmp(1),vgrad_p1)
!    pot_p1=poisson%pot
!    gw_s1=gw_s1+1.d-2
!    call cal_pot_gauss_p(parini,poisson,.true.,atoms%ratp(1,1),gw_s1,p_tmp(1),vgrad_p2)
!    pot_p2=poisson%pot
!    !do iz=1,poisson%ngpz
!    do iz=poisson%ngpz/2,poisson%ngpz/2
!    do iy=poisson%ngpy/2,poisson%ngpy/2
!    !do ix=poisson%ngpx/2,poisson%ngpx/2
!    do ix=1,poisson%ngpx
!        x=poisson%xyz111(1)+(ix-1)*poisson%hgrid(1,1)
!        y=poisson%xyz111(2)+(iy-1)*poisson%hgrid(2,2)
!        z=poisson%xyz111(3)+(iz-1)*poisson%hgrid(3,3)
!        write(83,'(3f8.3,2es19.10)') x,y,z,(pot_p2(ix,iy,iz)-pot_p1(ix,iy,iz))/1.d-2,vgrad_p1(ix,iy,iz)
!    enddo
!    enddo
!    enddo
!    stop 'TTTTTTTTTTTTTTTTTTTTT'




!    write(*,*) 'atoms%nat= ',atoms%nat
!    write(*,*) 'parini%ntypat= ',parini%ntypat
!    write(*,*) 'atoms%itypat= ',atoms%itypat
    do iat=1,atoms%nat
        if(trim(atoms%sat(iat))=='Mg') then
        c_s(iat)=-2.d0+1.025d0
        elseif(trim(atoms%sat(iat))=='O') then
        c_s(iat)=gg6-1.025d0 !-4.d0
        else
            stop 'ERROR: unknown type in get_basis_functions_cent2'
        endif
        !c_s(iat)=gg6 !-4.d0
        c_p(iat)=0.d0 !0.3634384d0
    enddo
!    stop 'MMMMMMMMMMMMMMMMMMMMMMM'

    do istep=0,500

    poisson%pot=pot_ion

    voxel=poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3)

    vgrad_s1=0.d0
    vgrad_s2=0.d0
    vgrad_p1=0.d0
    vgrad_p2=0.d0
    do iat=1,atoms%nat
    itypat=atoms%itypat(iat)
    bc_s1= gw_s1(itypat)**3/(gw_s1(itypat)**3-gw_s2(itypat)**3)
    bc_s2=-gw_s2(itypat)**3/(gw_s1(itypat)**3-gw_s2(itypat)**3)
    bc_p1= gw_p1(itypat)**5/(gw_p1(itypat)**5-gw_p2(itypat)**5)
    bc_p2=-gw_p2(itypat)**5/(gw_p1(itypat)**5-gw_p2(itypat)**5)
    q_tmp(1)=1.d0
    call cal_pot_gauss_s(parini,poisson,pot_s1,.true.,atoms%ratp(1,iat),gw_s1(itypat),q_tmp(1),vgrad_s1_t)
    call cal_pot_gauss_s(parini,poisson,pot_s2,.true.,atoms%ratp(1,iat),gw_s2(itypat),q_tmp(1),vgrad_s2_t)
    p_tmp=0.d0
    p_tmp(1)=1.d0
    call cal_pot_gauss_p(parini,poisson,pot_p1,.true.,atoms%ratp(1,iat),gw_p1(itypat),p_tmp(1),vgrad_p1_t)
    call cal_pot_gauss_p(parini,poisson,pot_p2,.true.,atoms%ratp(1,iat),gw_p2(itypat),p_tmp(1),vgrad_p2_t)

    poisson%pot=poisson%pot+ &
        c_s(iat)*(bc_s1*pot_s1+bc_s2*pot_s2)+ &
        c_p(iat)*(bc_p1*pot_p1+bc_p2*pot_p2)

    bcg_s1(1)=-3.d0*gw_s1(itypat)**2*gw_s2(itypat)**3/(gw_s1(itypat)**3-gw_s2(itypat)**3)**2
    bcg_s1(2)= 3.d0*gw_s1(itypat)**3*gw_s2(itypat)**2/(gw_s1(itypat)**3-gw_s2(itypat)**3)**2
    bcg_s2(1)= 3.d0*gw_s2(itypat)**3*gw_s1(itypat)**2/(gw_s2(itypat)**3-gw_s1(itypat)**3)**2
    bcg_s2(2)=-3.d0*gw_s2(itypat)**2*gw_s1(itypat)**3/(gw_s2(itypat)**3-gw_s1(itypat)**3)**2
    bcg_p1(1)=-5.d0*gw_p1(itypat)**4*gw_p2(itypat)**5/(gw_p1(itypat)**5-gw_p2(itypat)**5)**2
    bcg_p1(2)= 5.d0*gw_p1(itypat)**5*gw_p2(itypat)**4/(gw_p1(itypat)**5-gw_p2(itypat)**5)**2
    bcg_p2(1)= 5.d0*gw_p2(itypat)**5*gw_p1(itypat)**4/(gw_p2(itypat)**5-gw_p1(itypat)**5)**2
    bcg_p2(2)=-5.d0*gw_p2(itypat)**4*gw_p1(itypat)**5/(gw_p2(itypat)**5-gw_p1(itypat)**5)**2
    do iz=1,poisson%ngpz
    do iy=1,poisson%ngpy
    do ix=1,poisson%ngpx
        tt1=bcg_s1(1)*pot_s1(ix,iy,iz)+bcg_s2(1)*pot_s2(ix,iy,iz)+bc_s1*vgrad_s1_t(ix,iy,iz)
        tt2=bcg_s1(2)*pot_s1(ix,iy,iz)+bcg_s2(2)*pot_s2(ix,iy,iz)+bc_s2*vgrad_s2_t(ix,iy,iz)
        tt3=bcg_p1(1)*pot_p1(ix,iy,iz)+bcg_p2(1)*pot_p2(ix,iy,iz)+bc_p1*vgrad_p1_t(ix,iy,iz)
        tt4=bcg_p1(2)*pot_p1(ix,iy,iz)+bcg_p2(2)*pot_p2(ix,iy,iz)+bc_p2*vgrad_p2_t(ix,iy,iz)
        vgrad_s1(ix,iy,iz,itypat)=vgrad_s1(ix,iy,iz,itypat)+c_s(iat)*tt1
        vgrad_s2(ix,iy,iz,itypat)=vgrad_s2(ix,iy,iz,itypat)+c_s(iat)*tt2
        vgrad_p1(ix,iy,iz,itypat)=vgrad_p1(ix,iy,iz,itypat)+c_p(iat)*tt3
        vgrad_p2(ix,iy,iz,itypat)=vgrad_p2(ix,iy,iz,itypat)+c_p(iat)*tt4
        pot_iat_s(ix,iy,iz,iat)=bc_s1*pot_s1(ix,iy,iz)+bc_s2*pot_s2(ix,iy,iz)
        pot_iat_p(ix,iy,iz,iat)=bc_p1*pot_p1(ix,iy,iz)+bc_p2*pot_p2(ix,iy,iz)
    enddo
    enddo
    enddo
    enddo !end of loop over iat
    !----------------------------------------------------------
    cost=0.d0
    do iz=1,poisson%ngpz
    do iy=1,poisson%ngpy
    do ix=1,poisson%ngpx
        tt0=poisson%pot(ix,iy,iz)-poisson_ref%pot(ix,iy,iz)
        cost=cost+tt0**2
    enddo
    enddo
    enddo
    cost=cost*voxel
    do itypat=1,parini%ntypat
        ts1=0.d0
        ts2=0.d0
        tp1=0.d0
        tp2=0.d0
        do iz=1,poisson%ngpz
        do iy=1,poisson%ngpy
        do ix=1,poisson%ngpx
            tt0=poisson%pot(ix,iy,iz)-poisson_ref%pot(ix,iy,iz)
            ts1=ts1+vgrad_s1(ix,iy,iz,itypat)*tt0
            ts2=ts2+vgrad_s2(ix,iy,iz,itypat)*tt0
            tp1=tp1+vgrad_p1(ix,iy,iz,itypat)*tt0
            tp2=tp2+vgrad_p2(ix,iy,iz,itypat)*tt0
        enddo
        enddo
        enddo
        gwgrad_s1(itypat)=2.d0*ts1*voxel
        gwgrad_s2(itypat)=2.d0*ts2*voxel
        gwgrad_p1(itypat)=2.d0*tp1*voxel
        gwgrad_p2(itypat)=2.d0*tp2*voxel
    enddo
    do iat=1,atoms%nat
        tt1=0.d0
        tt2=0.d0
        do iz=1,poisson%ngpz
        do iy=1,poisson%ngpy
        do ix=1,poisson%ngpx
            tt0=poisson%pot(ix,iy,iz)-poisson_ref%pot(ix,iy,iz)
            tt1=tt1+pot_iat_s(ix,iy,iz,iat)*tt0
            tt2=tt2+pot_iat_p(ix,iy,iz,iat)*tt0
        enddo
        enddo
        enddo
        cgrad_s(iat)=2.d0*tt1*voxel
        cgrad_p(iat)=2.d0*tt2*voxel
    enddo
    gradtot=0.d0
    do iat=1,atoms%nat
        gradtot=gradtot+cgrad_s(iat)
    enddo
    gradtot=gradtot/real(atoms%nat,kind=8)
    do iat=1,atoms%nat
        cgrad_s(iat)=cgrad_s(iat)-gradtot
    enddo
    write(*,'(a,i6,es14.5,6es10.1,7f8.3,es10.1)') 'cost ',istep,cost,gwgrad_s1(1),gwgrad_s2(1),gwgrad_p1(1),gwgrad_p2(1),cgrad_s(1),cgrad_p(1),gw_s1(1),gw_s2(1),gw_p1(1),gw_p2(1),c_s(1),c_s(2),sum(c_p(:)),alpha

    !alpha=4.d-4
    !if(abs(gwgrad_s1(1))<1.d-1 .and. abs(gwgrad_s2(1))<1.d-1) alpha=4.d-2
    gw_s1=gw_s1- 1.d0*alpha*gwgrad_s1
    gw_s2=gw_s2- 1.d0*alpha*gwgrad_s2
    gw_p1=gw_p1-30.d0*alpha*gwgrad_p1
    gw_p2=gw_p2-20.d0*alpha*gwgrad_p2
    c_s=c_s-1.d-1*cgrad_s
    c_p=c_p-1.d-1*cgrad_p
    !write(*,*) 'qtot= ',sum(c_s)

    !bc_s1= gw_s1**3/(gw_s1**3-gw_s2**3)
    !bc_s2=-gw_s2**3/(gw_s1**3-gw_s2**3)
    !bc_p1= gw_p1**5/(gw_p1**5-gw_p2**5)
    !bc_p2=-gw_p2**5/(gw_p1**5-gw_p2**5)


    enddo !end of loop over istep


    errmax=0.d0
    rmse=0.d0
    do iz=1,poisson%ngpz
    do iy=1,poisson%ngpy
    do ix=1,poisson%ngpx
        errmax=max(errmax,abs(poisson%pot(ix,iy,iz)-poisson_ref%pot(ix,iy,iz)))
        rmse=rmse+(poisson%pot(ix,iy,iz)-poisson_ref%pot(ix,iy,iz))**2
    enddo
    enddo
    enddo
    rmse=sqrt(rmse/real(poisson%ngpx*poisson%ngpy*poisson%ngpz,kind=8))
    write(*,'(a,2es14.5)') 'errmax,rmse= ',errmax,rmse
    open(unit=2134,file="gw.out",status='replace')
    write(2134,'(4f15.10)') gg1,gg2,gg3,alpha_gg
    write(2134,'(4f15.10)') gw_s1(1),gw_s2(1),gg6,alpha
    close(2134)
    !do iz=1,poisson%ngpz
    do iz=poisson%ngpz/2,poisson%ngpz/2
    do iy=poisson%ngpy/2,poisson%ngpy/2
    !do ix=poisson%ngpx/2,poisson%ngpx/2
    do ix=1,poisson%ngpx
        x=poisson%xyz111(1)+(ix-1)*poisson%hgrid(1,1)
        y=poisson%xyz111(2)+(iy-1)*poisson%hgrid(2,2)
        z=poisson%xyz111(3)+(iz-1)*poisson%hgrid(3,3)
        write(83,'(3f8.3,2es19.10)') x,y,z,poisson_ref%pot(ix,iy,iz),poisson%pot(ix,iy,iz)
    enddo
    enddo
    enddo

    do iconf=1,atoms_arr%nconf
        call atom_deallocate_old(atoms_arr%atoms(iconf))
    enddo
    deallocate(atoms_arr%atoms)
    call fini_hartree(parini,atoms,poisson_ref)
    call fini_hartree(parini,atoms,poisson)

    call f_release_routine()
end subroutine get_basis_functions_cent2
!*****************************************************************************************
subroutine cal_pot_gauss_p(parini,poisson,pot,reset,xyz,gw,p,vgrad)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp, typ_atoms_arr, atom_deallocate_old
    use mod_ann_io_yaml, only: read_data_yaml
    use mod_electrostatics, only: typ_poisson
    use mod_linkedlists, only: typ_linkedlists
    use mod_linked_lists, only: typ_pia_arr
    use mod_flm_futile
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_poisson), intent(in):: poisson
    real(8), intent(inout):: pot(poisson%ngpx,poisson%ngpy,poisson%ngpz)
    logical, intent(in):: reset
    real(8), intent(in):: xyz(3), gw, p(3)
    real(8), intent(inout):: vgrad(poisson%ngpx,poisson%ngpy,poisson%ngpz)
    !local variables
    !type(typ_atoms):: atoms
    real(8):: pi, sf, ss0, ss1, ss2, ss3, ss4
    real(8):: x, y, z, r, tt1, p_dot_r
    real(8):: uu1, vv1, qq1, sft
    real(8):: yy1, ff1, hh1
    real(8):: tg0, tg1, tg2, tg3, tg4
    integer:: ix, iy, iz
    pi=4.d0*atan(1.d0)
    if(reset) then
        pot=0.d0
    endif
    sf=parini%screening_factor
    do iz=1,poisson%ngpz
    do iy=1,poisson%ngpy
    do ix=1,poisson%ngpx
        x=poisson%xyz111(1)+(ix-1)*poisson%hgrid(1,1)-xyz(1)
        y=poisson%xyz111(2)+(iy-1)*poisson%hgrid(2,2)-xyz(2)
        z=poisson%xyz111(3)+(iz-1)*poisson%hgrid(3,3)-xyz(3)
        r=sqrt(x**2+y**2+z**2)
        p_dot_r=x*p(1)+y*p(2)+z*p(3)
        if(r<1.d-10) then
            ss0=4.d0*p_dot_r/(3.d0*sqrt(pi)*gw**3)
            tg0=-3.d0*4.d0*p_dot_r/(3.d0*sqrt(pi)*gw**4)
            sft=sf
            tt1=3.d0*sqrt(pi)*gw**3*(1.d0+sft**2*gw**2)**2
            ss1=4.d0*p_dot_r*(1.d0+3.d0*sft**2*gw**2)/tt1
            yy1=3.d0+10.d0*gw**2*sft**2+15.d0*gw**4*sft**4
            tg1=(-4.d0*p_dot_r*yy1)/(3.d0*sqrt(pi)*gw**4*tt1**3)
            sft=2.d0*sf
            tt1=3.d0*sqrt(pi)*gw**3*(1.d0+sft**2*gw**2)**2
            ss2=4.d0*p_dot_r*(1.d0+3.d0*sft**2*gw**2)/tt1
            yy1=3.d0+10.d0*gw**2*sft**2+15.d0*gw**4*sft**4
            tg2=(-4.d0*p_dot_r*yy1)/(3.d0*sqrt(pi)*gw**4*tt1**3)
            sft=3.d0*sf
            tt1=3.d0*sqrt(pi)*gw**3*(1.d0+sft**2*gw**2)**2
            ss3=4.d0*p_dot_r*(1.d0+3.d0*sft**2*gw**2)/tt1
            yy1=3.d0+10.d0*gw**2*sft**2+15.d0*gw**4*sft**4
            tg3=(-4.d0*p_dot_r*yy1)/(3.d0*sqrt(pi)*gw**4*tt1**3)
            sft=4.d0*sf
            tt1=3.d0*sqrt(pi)*gw**3*(1.d0+sft**2*gw**2)**2
            ss4=4.d0*p_dot_r*(1.d0+3.d0*sft**2*gw**2)/tt1
            yy1=3.d0+10.d0*gw**2*sft**2+15.d0*gw**4*sft**4
            tg4=(-4.d0*p_dot_r*yy1)/(3.d0*sqrt(pi)*gw**4*tt1**3)
        else
            ss0=p_dot_r*(-2.d0*r*exp(-(r/gw)**2)/(sqrt(pi)*gw)+erf(r/gw))/r**3
            tg0=-4.d0*p_dot_r*exp(-(r/gw)**2)/(sqrt(pi)*gw**4)
            sft=sf
            tt1=1.d0+sft**2*gw**2
            uu1=exp(-(sft*r)**2/tt1)*sqrt(pi)*gw
            vv1=(1.d0+2.d0*r**2*sft**2+gw**2*sft**2)*erf(r/(gw*sqrt(tt1)))
            qq1=sqrt(pi)*r**3*gw*(tt1)**1.5d0
            ss1=p_dot_r*(-2.d0*r*sqrt(tt1)*exp(-(r/gw)**2)+uu1*vv1)/qq1
            hh1=2.d0*r*tt1*(sft**2*gw**4*tt1-2.d0*r**2*(1.d0+3.d0*sft**2*gw**2*tt1))
            yy1=sft**2*gw**5*sqrt(pi*tt1)*erf(r/(gw*sqrt(tt1)))
            ff1=(-4.d0*r**4*sft**4+tt1**2+4.d0*r**2*sft**2*tt1)*exp(-(sft*r)**2/tt1)
            tg1=p_dot_r*(hh1*exp(-(r/gw)**2)-yy1*ff1)/(sqrt(pi)*r**3*gw**4*tt1**4)
            sft=2.d0*sf
            tt1=1.d0+sft**2*gw**2
            uu1=exp(-(sft*r)**2/tt1)*sqrt(pi)*gw
            vv1=(1.d0+2.d0*r**2*sft**2+gw**2*sft**2)*erf(r/(gw*sqrt(tt1)))
            qq1=sqrt(pi)*r**3*gw*(tt1)**1.5d0
            ss2=p_dot_r*(-2.d0*r*sqrt(tt1)*exp(-(r/gw)**2)+uu1*vv1)/qq1
            hh1=2.d0*r*tt1*(sft**2*gw**4*tt1-2.d0*r**2*(1.d0+3.d0*sft**2*gw**2*tt1))
            yy1=sft**2*gw**5*sqrt(pi*tt1)*erf(r/(gw*sqrt(tt1)))
            ff1=(-4.d0*r**4*sft**4+tt1**2+4.d0*r**2*sft**2*tt1)*exp(-(sft*r)**2/tt1)
            tg2=p_dot_r*(hh1*exp(-(r/gw)**2)-yy1*ff1)/(sqrt(pi)*r**3*gw**4*tt1**4)
            sft=3.d0*sf
            tt1=1.d0+sft**2*gw**2
            uu1=exp(-(sft*r)**2/tt1)*sqrt(pi)*gw
            vv1=(1.d0+2.d0*r**2*sft**2+gw**2*sft**2)*erf(r/(gw*sqrt(tt1)))
            qq1=sqrt(pi)*r**3*gw*(tt1)**1.5d0
            ss3=p_dot_r*(-2.d0*r*sqrt(tt1)*exp(-(r/gw)**2)+uu1*vv1)/qq1
            hh1=2.d0*r*tt1*(sft**2*gw**4*tt1-2.d0*r**2*(1.d0+3.d0*sft**2*gw**2*tt1))
            yy1=sft**2*gw**5*sqrt(pi*tt1)*erf(r/(gw*sqrt(tt1)))
            ff1=(-4.d0*r**4*sft**4+tt1**2+4.d0*r**2*sft**2*tt1)*exp(-(sft*r)**2/tt1)
            tg3=p_dot_r*(hh1*exp(-(r/gw)**2)-yy1*ff1)/(sqrt(pi)*r**3*gw**4*tt1**4)
            sft=4.d0*sf
            tt1=1.d0+sft**2*gw**2
            uu1=exp(-(sft*r)**2/tt1)*sqrt(pi)*gw
            vv1=(1.d0+2.d0*r**2*sft**2+gw**2*sft**2)*erf(r/(gw*sqrt(tt1)))
            qq1=sqrt(pi)*r**3*gw*(tt1)**1.5d0
            ss4=p_dot_r*(-2.d0*r*sqrt(tt1)*exp(-(r/gw)**2)+uu1*vv1)/qq1
            hh1=2.d0*r*tt1*(sft**2*gw**4*tt1-2.d0*r**2*(1.d0+3.d0*sft**2*gw**2*tt1))
            yy1=sft**2*gw**5*sqrt(pi*tt1)*erf(r/(gw*sqrt(tt1)))
            ff1=(-4.d0*r**4*sft**4+tt1**2+4.d0*r**2*sft**2*tt1)*exp(-(sft*r)**2/tt1)
            tg4=p_dot_r*(hh1*exp(-(r/gw)**2)-yy1*ff1)/(sqrt(pi)*r**3*gw**4*tt1**4)
        endif
        !write(*,'(3i4,5es14.5)') ix,iy,iz,ss0,ss1,ss2,ss3,ss4
        pot(ix,iy,iz)=pot(ix,iy,iz)+ss0-4.d0*ss1+6.d0*ss2-4.d0*ss3+ss4
        vgrad(ix,iy,iz)=tg0-4.d0*tg1+6.d0*tg2-4.d0*tg3+tg4
        !poisson%pot(ix,iy,iz)=ss0-2.d0*ss1+ss2
    enddo
    enddo
    enddo
end subroutine cal_pot_gauss_p
!*****************************************************************************************
subroutine cal_pot_gauss_s(parini,poisson,pot,reset,xyz,gw,q,vgrad)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp, typ_atoms_arr, atom_deallocate_old
    use mod_ann_io_yaml, only: read_data_yaml
    use mod_electrostatics, only: typ_poisson
    use mod_linkedlists, only: typ_linkedlists
    use mod_linked_lists, only: typ_pia_arr
    use mod_flm_futile
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_poisson), intent(in):: poisson
    real(8), intent(inout):: pot(poisson%ngpx,poisson%ngpy,poisson%ngpz)
    logical, intent(in):: reset
    real(8), intent(in):: xyz(3), gw, q
    real(8), intent(inout):: vgrad(poisson%ngpx,poisson%ngpy,poisson%ngpz)
    !local variables
    !type(typ_atoms):: atoms
    real(8):: pi, sf, ss0, ss1, ss2, ss3, ss4
    real(8):: x, y, z, r, tt1, uu1, ww1, qq1, vv1, sft
    real(8):: tg0, tg1, tg2, tg3, tg4
    integer:: ix, iy, iz
    pi=4.d0*atan(1.d0)
    if(reset) then
        pot=0.d0
    endif
    sf=parini%screening_factor
    do iz=1,poisson%ngpz
    do iy=1,poisson%ngpy
    do ix=1,poisson%ngpx
        x=poisson%xyz111(1)+(ix-1)*poisson%hgrid(1,1)-xyz(1)
        y=poisson%xyz111(2)+(iy-1)*poisson%hgrid(2,2)-xyz(2)
        z=poisson%xyz111(3)+(iz-1)*poisson%hgrid(3,3)-xyz(3)
        r=sqrt(x**2+y**2+z**2)
        if(r<1.d-10) then
            ss0=2.d0/(sqrt(pi)*gw)
            ss1=2.d0/(sqrt(pi)*gw*(1.d0+sf**2*gw**2))
            ss2=2.d0/(sqrt(pi)*gw*(1.d0+(2.d0*sf)**2*gw**2))
            ss3=2.d0/(sqrt(pi)*gw*(1.d0+(3.d0*sf)**2*gw**2))
            ss4=2.d0/(sqrt(pi)*gw*(1.d0+(4.d0*sf)**2*gw**2))
            tg0=-2.d0/(sqrt(pi)*gw**2)
            tg1=-2.d0*(1.d0+3.d0*sf**2*gw**2)/(sqrt(pi)*(gw*(1.d0+sf**2*gw**2))**2)
            tg2=-2.d0*(1.d0+3.d0*(2.d0*sf)**2*gw**2)/(sqrt(pi)*(gw*(1.d0+(2.d0*sf)**2*gw**2))**2)
            tg3=-2.d0*(1.d0+3.d0*(3.d0*sf)**2*gw**2)/(sqrt(pi)*(gw*(1.d0+(3.d0*sf)**2*gw**2))**2)
            tg4=-2.d0*(1.d0+3.d0*(4.d0*sf)**2*gw**2)/(sqrt(pi)*(gw*(1.d0+(4.d0*sf)**2*gw**2))**2)
        else
            ss0=erf(r/gw)/r
            tt1=1.d0+sf**2*gw**2
            ss1=erf(r/(gw*sqrt(tt1)))*exp(-(sf*r)**2/tt1)/(r*sqrt(tt1))
            tt1=1.d0+(2.d0*sf)**2*gw**2
            ss2=erf(r/(gw*sqrt(tt1)))*exp(-(2.d0*sf*r)**2/tt1)/(r*sqrt(tt1))
            tt1=1.d0+(3.d0*sf)**2*gw**2
            ss3=erf(r/(gw*sqrt(tt1)))*exp(-(3.d0*sf*r)**2/tt1)/(r*sqrt(tt1))
            tt1=1.d0+(4.d0*sf)**2*gw**2
            ss4=erf(r/(gw*sqrt(tt1)))*exp(-(4.d0*sf*r)**2/tt1)/(r*sqrt(tt1))
            tg0=-2.d0*exp(-(r/gw)**2)/(sqrt(pi)*gw**2)
            sft=sf
            tt1=1.d0+sft**2*gw**2
            qq1=1.d0+2.d0*sft**2*gw**2
            ww1=(1.d0-2.d0*r**2*sft**2+gw**2*sft**2)*erf(r/(gw*sqrt(tt1)))
            uu1=sqrt(pi)*r*gw**2*tt1**2.5d0
            vv1=exp(-(sft*r)**2/tt1)*sqrt(pi)*gw**3*sft**2
            tg1=(-2.d0*r*sqrt(tt1)*qq1*exp(-(r/gw)**2)-vv1*ww1)/uu1
            sft=2.d0*sf
            tt1=1.d0+sft**2*gw**2
            qq1=1.d0+2.d0*sft**2*gw**2
            ww1=(1.d0-2.d0*r**2*sft**2+gw**2*sft**2)*erf(r/(gw*sqrt(tt1)))
            uu1=sqrt(pi)*r*gw**2*tt1**2.5d0
            vv1=exp(-(sft*r)**2/tt1)*sqrt(pi)*gw**3*sft**2
            tg2=(-2.d0*r*sqrt(tt1)*qq1*exp(-(r/gw)**2)-vv1*ww1)/uu1
            sft=3.d0*sf
            tt1=1.d0+sft**2*gw**2
            qq1=1.d0+2.d0*sft**2*gw**2
            ww1=(1.d0-2.d0*r**2*sft**2+gw**2*sft**2)*erf(r/(gw*sqrt(tt1)))
            uu1=sqrt(pi)*r*gw**2*tt1**2.5d0
            vv1=exp(-(sft*r)**2/tt1)*sqrt(pi)*gw**3*sft**2
            tg3=(-2.d0*r*sqrt(tt1)*qq1*exp(-(r/gw)**2)-vv1*ww1)/uu1
            sft=4.d0*sf
            tt1=1.d0+sft**2*gw**2
            qq1=1.d0+2.d0*sft**2*gw**2
            ww1=(1.d0-2.d0*r**2*sft**2+gw**2*sft**2)*erf(r/(gw*sqrt(tt1)))
            uu1=sqrt(pi)*r*gw**2*tt1**2.5d0
            vv1=exp(-(sft*r)**2/tt1)*sqrt(pi)*gw**3*sft**2
            tg4=(-2.d0*r*sqrt(tt1)*qq1*exp(-(r/gw)**2)-vv1*ww1)/uu1
        endif
        pot(ix,iy,iz)=pot(ix,iy,iz)+q*(ss0-4.d0*ss1+6.d0*ss2-4.d0*ss3+ss4)
        vgrad(ix,iy,iz)=q*(tg0-4.d0*tg1+6.d0*tg2-4.d0*tg3+tg4)
    enddo
    enddo
    enddo
end subroutine cal_pot_gauss_s
!*****************************************************************************************
subroutine get_dpm_fit_bf_cent2(atoms,dpx,dpy,dpz,dpm_err)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(in):: atoms
    real(8), intent(out):: dpm_err 
    real(8), intent(out):: dpx, dpy, dpz
    !local variables
    integer:: iat
    real(8):: centroid_x, centroid_y, centroid_z
    dpx=0.d0
    dpy=0.d0
    dpz=0.d0
    centroid_x=sum(atoms%ratp(1,:))/atoms%nat
    centroid_y=sum(atoms%ratp(2,:))/atoms%nat
    centroid_z=sum(atoms%ratp(3,:))/atoms%nat
    write(*,'(a,5f8.2)') 'DPM ',atoms%qat(1),atoms%qat(2),atoms%ratp(1,2)-atoms%ratp(1,1),atoms%zat(1),atoms%zat(2)
    do iat=1, atoms%nat
        if(trim(atoms%sat(iat))=='Mg') then
            dpx=dpx+(atoms%ratp(1,iat)-centroid_x)*(atoms%qat(iat)+2.d0)
            dpy=dpy+(atoms%ratp(2,iat)-centroid_y)*(atoms%qat(iat)+2.d0)
            dpz=dpz+(atoms%ratp(3,iat)-centroid_z)*(atoms%qat(iat)+2.d0)
        elseif(trim(atoms%sat(iat))=='O') then
            dpx=dpx+(atoms%ratp(1,iat)-centroid_x)*(atoms%qat(iat)+4.d0)
            dpy=dpy+(atoms%ratp(2,iat)-centroid_y)*(atoms%qat(iat)+4.d0)
            dpz=dpz+(atoms%ratp(3,iat)-centroid_z)*(atoms%qat(iat)+4.d0)
        endif
    enddo
    dpm_err=((dpx-atoms%dpm(1))**2+(dpy-atoms%dpm(2))**2+(dpz-atoms%dpm(3))**2)
end subroutine get_dpm_fit_bf_cent2
!*****************************************************************************************
subroutine get_poisson_ref(parini,poisson,atoms)
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms, get_rat, update_ratp, update_rat, atom_deallocate_old
    use mod_flm_futile
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_poisson), intent(inout):: poisson
    type(typ_atoms), intent(inout):: atoms
    !local variables
    type(typ_poisson):: poisson_ion
    integer:: igpx, igpy, igpz, iat
    integer:: jgpx, jgpy, jgpz
    real(8):: rgcut_a, qtot, pi
    real(8):: ehartree_scn_excl, tt1
    real(8):: dx, dy, dz, r2, coeff, gwt
    !real(8):: xmin, ymin, zmin, xmax, ymax, zmax
    real(8), allocatable::  gausswidth(:)
    real(8), allocatable::  rho(:,:,:)
    integer:: nbgpx, nbgpy, nbgpz
    pi=4.d0*atan(1.d0)
    if(parini%ewald) then
        write(*,*) 'ERROR: ewald=True is wrong when reading from cube file.'
        stop
    endif
    if(trim(parini%psolver)=='kwald') then
        write(*,*) 'ERROR: psolver=kwald is wrong for grid base charge density.'
        stop
    endif
    if(parini%cal_scn .and. parini%screening_factor==0.d0) then
        stop 'ERROR: cal_scn is TRUE and screening_factor is 0, MEANINGLESS!'
    endif
    !-------------------------------------------------------
    call cube_read('rho.cube',atoms,poisson)
    allocate(rho(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    do igpz=1,poisson%ngpz
    do igpy=1,poisson%ngpy
    do igpx=1,poisson%ngpx
        rho(igpx,igpy,igpz)=poisson%rho(igpx,igpy,igpz)
    enddo
    enddo
    enddo
    write(*,*) allocated(poisson%rho)
    call f_free(poisson%rho)
    poisson%rho=f_malloc0([1.to.(poisson%ngpx+80),1.to.(poisson%ngpy+80),1.to.(poisson%ngpz+80)],id='poisson%rho')
    do igpz=1,poisson%ngpz
    do igpy=1,poisson%ngpy
    do igpx=1,poisson%ngpx
        poisson%rho(igpx+40,igpy+40,igpz+40)=rho(igpx,igpy,igpz)
    enddo
    enddo
    enddo
    call update_ratp(atoms)
    do iat=1,atoms%nat
        atoms%ratp(1,iat)=atoms%ratp(1,iat)+real(40,kind=8)*poisson%hgrid(1,1)
        atoms%ratp(2,iat)=atoms%ratp(2,iat)+real(40,kind=8)*poisson%hgrid(2,2)
        atoms%ratp(3,iat)=atoms%ratp(3,iat)+real(40,kind=8)*poisson%hgrid(3,3)
    enddo
    call update_rat(atoms)
    deallocate(rho)
    poisson%ngpx=poisson%ngpx+80
    poisson%ngpy=poisson%ngpy+80
    poisson%ngpz=poisson%ngpz+80
    !-------------------------------------------------------
    atoms%boundcond='free'
    poisson%bc=atoms%boundcond
    allocate(gausswidth(atoms%nat))
    gausswidth=0.5d0 !parini%gaussian_width !TO_BE_CORRECTED
    poisson%cal_scn=parini%cal_scn
    poisson%screening_factor=parini%screening_factor
    poisson%task_finit=""
    call init_hartree(parini,atoms,poisson,gausswidth)
    poisson%cell(1)=poisson%hgrid(1,1)*poisson%ngpx
    poisson%cell(2)=poisson%hgrid(2,2)*poisson%ngpy
    poisson%cell(3)=poisson%hgrid(3,3)*poisson%ngpz
    !-------------------------------------------------------
    poisson_ion%alpha=0.5d0 !parini%gaussian_width !TO_BE_CORRECTED
    poisson_ion%rgcut=parini%rgcut_ewald*poisson_ion%alpha
    poisson_ion%ngpx=poisson%ngpx
    poisson_ion%ngpy=poisson%ngpy
    poisson_ion%ngpz=poisson%ngpz
    poisson_ion%hgrid(1:3,1:3)=0.d0
    poisson_ion%hgrid(1,1)=poisson%hgrid(1,1)
    poisson_ion%hgrid(2,2)=poisson%hgrid(2,2)
    poisson_ion%hgrid(3,3)=poisson%hgrid(3,3)
    poisson_ion%xyz111=poisson%xyz111
    rgcut_a=parini%rgcut_ewald*maxval(gausswidth) !parini%gaussian_width !3.d0
    nbgpx=int(rgcut_a/poisson_ion%hgrid(1,1))+2
    nbgpy=int(rgcut_a/poisson_ion%hgrid(2,2))+2
    nbgpz=int(rgcut_a/poisson_ion%hgrid(3,3))+2
    poisson_ion%task_finit="alloc_rho"
    call init_hartree(parini,atoms,poisson_ion,gausswidth)
    poisson_ion%reset_rho=.true.
    poisson_ion%nat=atoms%nat
    poisson_ion%cv=atoms%cellvec
    poisson_ion%bc=atoms%boundcond
    poisson_ion%q(1:poisson_ion%nat)=atoms%zat(1:atoms%nat)
    poisson_ion%gw(1:poisson_ion%nat)=gausswidth(1:atoms%nat)
    call get_rat(atoms,poisson_ion%rcart)
    !call put_charge_density(parini,poisson_ion)
    nbgpx=int(7.d0*maxval(gausswidth(1:atoms%nat))/poisson_ion%hgrid(1,1))+2
    nbgpy=int(7.d0*maxval(gausswidth(1:atoms%nat))/poisson_ion%hgrid(2,2))+2
    nbgpz=int(7.d0*maxval(gausswidth(1:atoms%nat))/poisson_ion%hgrid(3,3))+2
    poisson_ion%rho=0.d0
    do iat=1,atoms%nat
    gwt=gausswidth(iat)
    coeff=atoms%zat(iat)/(gwt**3*pi**1.5d0)
    jgpx=int(poisson_ion%rcart(1,iat)/poisson%hgrid(1,1))
    jgpy=int(poisson_ion%rcart(2,iat)/poisson%hgrid(2,2))
    jgpz=int(poisson_ion%rcart(3,iat)/poisson%hgrid(3,3))
    do igpz=jgpz-nbgpz,jgpz+nbgpz
    do igpy=jgpy-nbgpy,jgpy+nbgpy
    do igpx=jgpx-nbgpx,jgpx+nbgpx
        dx=(igpx-1)*poisson%hgrid(1,1)-poisson_ion%rcart(1,iat)
        dy=(igpy-1)*poisson%hgrid(2,2)-poisson_ion%rcart(2,iat)
        dz=(igpz-1)*poisson%hgrid(3,3)-poisson_ion%rcart(3,iat)
        r2=dx**2+dy**2+dz**2
        if(r2<10.d0**2*gwt**2) then
            poisson_ion%rho(igpx,igpy,igpz)=poisson_ion%rho(igpx,igpy,igpz)+coeff*exp(-r2/gwt**2)
        endif
    enddo
    enddo
    enddo
    enddo
    tt1=0.d0
    do igpz=1,poisson%ngpz
    do igpy=1,poisson%ngpy
    do igpx=1,poisson%ngpx
        tt1=tt1+poisson_ion%rho(igpx,igpy,igpz)
    enddo
    enddo
    enddo
    tt1=tt1*(poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3))
    if(parini%mpi_env%iproc==0) then
    write(*,*) 'TT1 ',tt1
    endif
    !-------------------------------------------------------
    qtot=0.d0
    do igpz=1,poisson%ngpz
        do igpy=1,poisson%ngpy
            do igpx=1,poisson%ngpx
                !poisson%rho(igpx,igpy,igpz)=poisson_ion%rho(igpx,igpy,igpz)-poisson%rho(igpx,igpy,igpz)
                poisson%rho(igpx,igpy,igpz)=-poisson%rho(igpx,igpy,igpz)
                qtot=qtot+poisson%rho(igpx,igpy,igpz)
            enddo
        enddo
    enddo
    qtot=qtot*poisson_ion%hgrid(1,1)*poisson_ion%hgrid(2,2)*poisson_ion%hgrid(3,3)
    if(parini%mpi_env%iproc==0) then
        write(*,*) 'qtot= ',qtot
    endif
    !-------------------------------------------------------
    call update_ratp(atoms)
    call get_hartree(parini,poisson,atoms,gausswidth,ehartree_scn_excl)
    !if(parini%mpi_env%iproc==0) then
    !write(*,'(a,es24.15,es14.5)') 'ehartree_scn_excl ',ehartree_scn_excl,poisson%screening_factor
    !endif
    !atoms%fat=0.d0
    !call force_gto_sym_ortho(parini,poisson_ion%bc,atoms%nat,poisson_ion%rcart, &
    !    poisson_ion%q,gausswidth,6.d0,poisson_ion%xyz111, &
    !    poisson_ion%ngpx,poisson_ion%ngpx,poisson_ion%ngpy,poisson_ion%ngpz, &
    !    poisson_ion%hgrid,poisson%pot,atoms%fat)
    !if(parini%mpi_env%iproc==0) then
    !do iat=1,atoms%nat
    !    write(*,'(a,i4,3es19.10)') 'FAT ',iat,atoms%fat(1,iat),atoms%fat(2,iat),atoms%fat(3,iat)
    !enddo
    !endif
    !-------------------------------------------------------
    call fini_hartree(parini,atoms,poisson_ion)
end subroutine get_poisson_ref
!*****************************************************************************************
subroutine get_proc_stake(mpi_env,n,is,ie)
    use mod_flm_futile
    implicit none
    type(mpi_environment), intent(in):: mpi_env
    integer, intent(in):: n
    integer, intent(out):: is, ie
    !local variables
    integer:: m, mproc
    if(mpi_env%nproc>1) then
        m=n/mpi_env%nproc
        is=mpi_env%iproc*m+1
        mproc=mod(n,mpi_env%nproc)
        is=is+max(0,mpi_env%iproc-mpi_env%nproc+mproc)
        if(mpi_env%iproc>mpi_env%nproc-mproc-1) m=m+1
        ie=is+m-1
    else
        is=1
        ie=n
    endif
end subroutine get_proc_stake
!*****************************************************************************************
end module mod_fit_bf_cent2
!*****************************************************************************************
