!*****************************************************************************************
subroutine init_potential_ann(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_potential, only: ann_arr, ann_boundcheck
    use mod_ann, only: set_number_of_ann
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    !local variables
    integer:: i, iat
    character(10):: fn
    character (50)::fname
    !write(*,*) trim(parini%stypat_ann)
    !call count_words(parini%stypat_ann,ann_arr%nann)
    ann_arr%approach=trim(parini%approach_ann)
    !write (*,*) 'parini         ', ann_arr%approach
    call set_number_of_ann(parini,ann_arr)
    if(ann_arr%nann==0) stop 'ERROR: number of type of atoms zero in init_potential_ann'
    call yaml_map('number of ANN',ann_arr%nann)
    !write(*,*) 'ann_arr%nann= ',ann_arr%nann
    !read(parini%stypat_ann,*) ann_arr%stypat(1:ann_arr%nann)
    !do i=1,ann_arr%nann
    !    ann_arr%ltypat(i)=i
    !    write(*,*) i,ann_arr%stypat(i)
    !enddo
    if(ann_arr%nann==0) stop 'ERROR: number of type of atoms zero in init_potential_ann'
    allocate(ann_arr%ann(ann_arr%nann))
    do iat=1,atoms%nat
        do i=1,ann_arr%nann !this should be changed for bond based ANN
            if(trim(adjustl(atoms%sat(iat)))==trim(parini%stypat(i))) then
                atoms%itypat(iat)=parini%ltypat(i)
                exit
            endif
        enddo
    enddo
    if(parini%bondbased_ann .and. trim(ann_arr%approach)=='tb') then
        if(parini%ntypat>1) then
            stop 'ERROR: writing ANN parameters for tb available only ntypat=1'
        endif
        !write(fn_tt,'(i1)') iann
        fname=trim(parini%stypat(1))//'1'//'.ann.param.yaml'
    else
        fname = trim(parini%stypat(1))//'.ann.param.yaml'
    endif
    inquire(file=trim(fname),exist=ann_arr%exists_yaml_file)
    if( ann_arr%exists_yaml_file) then
        call read_ann_yaml(parini,ann_arr)
    else
        call read_ann(parini,ann_arr)
    endif
    ann_boundcheck=trim(parini%potential_ann_boundcheck)
    ann_arr%event='potential'
end subroutine init_potential_ann
!*****************************************************************************************
subroutine cal_potential_ann(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, atom_deallocate_old, get_rat
    use mod_potential, only: ann_arr, fcalls, fcalls_sec, potential, potential_sec, ann_boundcheck
    use mod_symfunc, only: typ_symfunc
    use mod_opt_ann, only: typ_opt_ann
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    !local variables
    integer:: iat, jat, i0, i
    real(8):: epoti, fcalls_t
    type(typ_symfunc):: symfunc
    type(typ_opt_ann):: opt_ann
    if(trim(potential)=='ann') then
        fcalls_t=fcalls
    elseif(trim(potential_sec)=='ann') then
        fcalls_t=fcalls_sec
    else
        stop 'ERROR: why is cal_potential_ann called?'
    endif
    if(trim(atoms%boundcond)=='free') then
        atoms%natim=atoms%nat
        if(.not. allocated(atoms%ratim)) then
            allocate(atoms%ratim(3,atoms%natim),source=0.d0)
        endif
        call get_rat(atoms,atoms%ratim)
    !elseif(trim(atoms%boundcond)=='bulk' .or. trim(atoms%boundcond)=='slab' .or. &
    !       trim(atoms%boundcond)=='wire') then
    !    call atom_build_periodic_images(atoms,8.d0)
    !else
    !    write(*,'(2a)') 'ERROR: unknown boundary conditions, boundcond=',trim(atoms%boundcond)
    !    stop
    endif
    !atoms%epot=0.d0
    !atoms%fat(1:3,1:atoms%nat)=0.d0
    !call cal_ann_cent1(atoms,symfunc,ann_arr,opt_ann)
    if(trim(ann_arr%event)/='potential') then
        write(*,*) 'ERROR: ann_arr%event different from potential'
        stop
    endif
    call cal_ann_main(parini,atoms,symfunc,ann_arr,opt_ann)
!    do iat=1,atoms%nat
!        call symmetry_functions(ann_arr%ann(i),iat,atoms,.true.)
!        !if(trim(ann_boundcheck)=='weak' .or. trim(ann_boundcheck)=='strong') then
!        !do i0=1,ann%nn(0)
!        !    if(ann%y(i0,0)<-1.0d0) then
!        !        write(*,'(a,i3,f8.2,f9.0)') 'WARNING: ann%y(i0,0)<-1.d0), i0=',i0,ann%y(i0,0),fcalls_t
!        !    endif
!        !    if(ann%y(i0,0)> 1.0d0) then
!        !        write(*,'(a,i3,f8.2,f9.0)') 'WARNING: ann%y(i0,0)> 1.d0), i0=',i0,ann%y(i0,0),fcalls_t
!        !    endif
!        !enddo
!        !endif
!        !call cal_architecture(ann,epoti)
!        !write(*,*) iat,epoti
!        atoms%epot=atoms%epot+epoti
!    enddo
    if(parini%add_repulsive) then
        call add_repulsive_potential(parini,atoms)
    endif
    if(allocated(atoms%ratim)) then
        call atom_deallocate_old(atoms,ratim=.true.)
    endif
    atoms%natim=0
end subroutine cal_potential_ann
!*****************************************************************************************
subroutine final_potential_ann
    use mod_potential, only: ann_arr
    implicit none
    !local variables
    deallocate(ann_arr%ann)
end subroutine final_potential_ann
!*****************************************************************************************
subroutine add_repulsive_potential(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, set_rcov, update_ratp
    use mod_linked_lists, only: typ_linked_lists
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    !local variables
    integer:: iat, jat, maincell_iat, maincell
    integer:: ip, jp, jpt, il, jl, iz, iy, ix, jx, jy, jz, iatp
    real(8):: rcovmax, cell(3), epot_rep, fx, fy, fz, ttt, a !, b, c, d
    real(8):: rc, rcsq, dx, dy, dz, xiat, yiat, ziat, r, rsq, tt1, tt2, tt4
    real(8):: frac1, frac12, ratio1, ratio2, t1, t2, t4, t6, t12, t14, sigma
    type(typ_linked_lists):: linked_lists
    !integer, save:: icall=0
    !icall=icall+1
    associate(rcmax=>linked_lists%rcut)
    call set_rcov(atoms)
    rcovmax=maxval(atoms%rcov(1:atoms%nat))
    rcmax=2.d0*rcovmax
    call linkedlists_init(parini,atoms,cell,linked_lists)
    frac1=0.72d0
    frac12=frac1**12
    !-------------------------------------------------------
    epot_rep=0.d0
    linked_lists%fat=0.d0
    call update_ratp(linked_lists%typ_atoms)
    include '../src/act1_cell_linkedlist.inc'
    do  iat=ip,il
        !qiat=linked_lists%qat(iat)
        xiat=linked_lists%ratp(1,iat)
        yiat=linked_lists%ratp(2,iat)
        ziat=linked_lists%ratp(3,iat)
        jpt=linked_lists%prime(ix+linked_lists%limnbx(1,jy-iy,jz-iz),jy,jz)
        jp=(iat-ip+1)*((isign(1,ip-jpt)+1)/2)+jpt
        jl=linked_lists%last(ix+linked_lists%limnbx(2,jy-iy,jz-iz),jy,jz)
        maincell_iat=linked_lists%maincell(iat)
        iatp=linked_lists%perm(iat)
        do  jat=jp,jl
            dx=xiat-linked_lists%ratp(1,jat)
            dy=yiat-linked_lists%ratp(2,jat)
            dz=ziat-linked_lists%ratp(3,jat)
            rsq= dx*dx+dy*dy+dz*dz
            maincell=maincell_iat+linked_lists%maincell(jat)
            rc=(atoms%rcov(iatp)+atoms%rcov(linked_lists%perm(jat)))*0.7d0
            rcsq=rc**2
            if(rsq<rcsq .and. maincell >-1) then
                r=sqrt(rsq)
                !---------------------------------
                !c=0.1d0
                !b=-2.d0*c/rc
                !a=c/rcsq
                !epot_rep=epot_rep+(c+r*(b+r*a))
                !ttt=-(2.d0*a*r+b)/r
                !---------------------------------
                !d=0.2d0
                !c=-3.d0*d/rc
                !b=3.d0*d/rcsq
                !a=-d/(rc*rcsq)
                !epot_rep=epot_rep+(d+r*(c+r*(b+r*a)))
                !ttt=-(c+r*(2.d0*b+r*3.d0*a))/r
                !---------------------------------
                !a=2.d0
                !tt1=(1.d0-rsq/rcsq)
                !tt2=tt1*tt1
                !tt4=tt2*tt2
                !epot_rep=epot_rep+a*tt4*tt2
                !ttt=12.d0*a*tt4*tt1/rcsq
                !---------------------------------
                sigma=rc*frac1
                t1=sigma/r
                t2=t1**2
                t4=t2**2
                t6=t2*t4
                t12=t6**2
                t14=t2*t12
                ratio1=r/rc
                ratio2=ratio1**2
                epot_rep=epot_rep+(t12+6.d0*frac12*ratio2-7.d0*frac12)
                ttt=12.d0*t14/sigma**2-12.d0*frac12/rcsq
                !---------------------------------
                fx=ttt*dx;fy=ttt*dy;fz=ttt*dz
                !write(*,'(a,i6,9f10.5)') 'HERE ',icall,r,rc,a,b,c,c+r*(b+r*a),fx,fy,fz
                !-----------------------------------
                linked_lists%fat(1,iat)=linked_lists%fat(1,iat)+fx
                linked_lists%fat(2,iat)=linked_lists%fat(2,iat)+fy
                linked_lists%fat(3,iat)=linked_lists%fat(3,iat)+fz
                linked_lists%fat(1,jat)=linked_lists%fat(1,jat)-fx
                linked_lists%fat(2,jat)=linked_lists%fat(2,jat)-fy
                linked_lists%fat(3,jat)=linked_lists%fat(3,jat)-fz
            endif
        enddo
    enddo
    include '../src/act2_cell_linkedlist.inc'
    !-------------------------------------------------------
    atoms%epot=atoms%epot+epot_rep
    do iat=1,linked_lists%natim
        iatp=linked_lists%perm(iat)
        atoms%fat(1,iatp)=atoms%fat(1,iatp)+linked_lists%fat(1,iat)
        atoms%fat(2,iatp)=atoms%fat(2,iatp)+linked_lists%fat(2,iat)
        atoms%fat(3,iatp)=atoms%fat(3,iatp)+linked_lists%fat(3,iat)
    enddo
    !-------------------------------------------------------
    call linkedlists_final(linked_lists)
    end associate
end subroutine add_repulsive_potential
!*****************************************************************************************
