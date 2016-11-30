!*****************************************************************************************
subroutine cal_force_chi_part1(parini,symfunc,iat,atoms,out_ann,ann_arr)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_symfunc), intent(in):: symfunc
    integer, intent(in):: iat
    type(typ_atoms), intent(in):: atoms
    real(8), intent(in):: out_ann
    type(typ_ann_arr), intent(inout):: ann_arr
    !local variables
    real(8):: sxx, sxy, sxz, syx, syy, syz, szx, szy, szz
    real(8):: ttx, tty, ttz, tt1, tt2
    integer:: ib, i, j
    i=atoms%itypat(iat)
    ann_arr%chi_i(iat)=out_ann
    tt1=tanh(ann_arr%ann(i)%prefactor_chi*out_ann)
    ann_arr%chi_o(iat)=ann_arr%ann(i)%ampl_chi*tt1+ann_arr%ann(i)%chi0
    if(.not. (trim(ann_arr%event)=='evalu' .and. atoms%nat>parini%nat_force)) then
        tt2=ann_arr%ann(i)%ampl_chi*ann_arr%ann(i)%prefactor_chi*(1.d0-tt1**2)
        do ib=symfunc%linked_lists%prime_bound(iat),symfunc%linked_lists%prime_bound(iat+1)-1
            ttx=0.d0 ; tty=0.d0 ; ttz=0.d0
            do j=1,ann_arr%ann(i)%nn(0)
                ttx=ttx+ann_arr%ann(i)%d(j)*symfunc%y0d(j,1,ib)
                tty=tty+ann_arr%ann(i)%d(j)*symfunc%y0d(j,2,ib)
                ttz=ttz+ann_arr%ann(i)%d(j)*symfunc%y0d(j,3,ib)
            enddo
            ann_arr%fatpq(1,ib)=ttx*tt2
            ann_arr%fatpq(2,ib)=tty*tt2
            ann_arr%fatpq(3,ib)=ttz*tt2
        enddo
    endif
    if(trim(ann_arr%event)=='potential') then
        do ib=symfunc%linked_lists%prime_bound(iat),symfunc%linked_lists%prime_bound(iat+1)-1
            sxx=0.d0 ; sxy=0.d0 ; sxz=0.d0
            syx=0.d0 ; syy=0.d0 ; syz=0.d0
            szx=0.d0 ; szy=0.d0 ; szz=0.d0
            do j=1,ann_arr%ann(i)%nn(0)
                sxx=sxx+ann_arr%ann(i)%d(j)*symfunc%y0dr(j,1,ib)
                sxy=sxy+ann_arr%ann(i)%d(j)*symfunc%y0dr(j,2,ib)
                sxz=sxz+ann_arr%ann(i)%d(j)*symfunc%y0dr(j,3,ib)
                syx=syx+ann_arr%ann(i)%d(j)*symfunc%y0dr(j,4,ib)
                syy=syy+ann_arr%ann(i)%d(j)*symfunc%y0dr(j,5,ib)
                syz=syz+ann_arr%ann(i)%d(j)*symfunc%y0dr(j,6,ib)
                szx=szx+ann_arr%ann(i)%d(j)*symfunc%y0dr(j,7,ib)
                szy=szy+ann_arr%ann(i)%d(j)*symfunc%y0dr(j,8,ib)
                szz=szz+ann_arr%ann(i)%d(j)*symfunc%y0dr(j,9,ib)
            enddo
            ann_arr%stresspq(1,1,ib)=-sxx*tt2
            ann_arr%stresspq(2,1,ib)=-syx*tt2
            ann_arr%stresspq(3,1,ib)=-szx*tt2
            ann_arr%stresspq(1,2,ib)=-sxy*tt2
            ann_arr%stresspq(2,2,ib)=-syy*tt2
            ann_arr%stresspq(3,2,ib)=-szy*tt2
            ann_arr%stresspq(1,3,ib)=-sxz*tt2
            ann_arr%stresspq(2,3,ib)=-syz*tt2
            ann_arr%stresspq(3,3,ib)=-szz*tt2
        enddo
    endif
end subroutine cal_force_chi_part1
!*****************************************************************************************
subroutine cal_force_chi_part2(parini,symfunc,atoms,ann_arr)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_symfunc), intent(in):: symfunc
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    !local variables
    real(8):: ttx, tty, ttz, qnet, hinv(3,3), vol
    integer:: ib, i, j, iat, jat
    call getvol_alborz(atoms%cellvec,vol)
    atoms%stress(1:3,1:3)=atoms%stress(1:3,1:3)*vol !*atoms%nat !not certain if this is needed!!!
    if(.not. (trim(ann_arr%event)=='evalu' .and. atoms%nat>parini%nat_force)) then
        do ib=1,symfunc%linked_lists%maxbound_rad
            iat=symfunc%linked_lists%bound_rad(1,ib)
            jat=symfunc%linked_lists%bound_rad(2,ib)
            qnet=atoms%qat(iat)
            ttx=ann_arr%fatpq(1,ib)*qnet
            tty=ann_arr%fatpq(2,ib)*qnet
            ttz=ann_arr%fatpq(3,ib)*qnet
            ann_arr%fat_chi(1,iat)=ann_arr%fat_chi(1,iat)+ttx
            ann_arr%fat_chi(2,iat)=ann_arr%fat_chi(2,iat)+tty
            ann_arr%fat_chi(3,iat)=ann_arr%fat_chi(3,iat)+ttz
            ann_arr%fat_chi(1,jat)=ann_arr%fat_chi(1,jat)-ttx
            ann_arr%fat_chi(2,jat)=ann_arr%fat_chi(2,jat)-tty
            ann_arr%fat_chi(3,jat)=ann_arr%fat_chi(3,jat)-ttz
            atoms%stress(1,1)=atoms%stress(1,1)+ann_arr%stresspq(1,1,ib)*qnet
            atoms%stress(2,1)=atoms%stress(2,1)+ann_arr%stresspq(2,1,ib)*qnet
            atoms%stress(3,1)=atoms%stress(3,1)+ann_arr%stresspq(3,1,ib)*qnet
            atoms%stress(1,2)=atoms%stress(1,2)+ann_arr%stresspq(1,2,ib)*qnet
            atoms%stress(2,2)=atoms%stress(2,2)+ann_arr%stresspq(2,2,ib)*qnet
            atoms%stress(3,2)=atoms%stress(3,2)+ann_arr%stresspq(3,2,ib)*qnet
            atoms%stress(1,3)=atoms%stress(1,3)+ann_arr%stresspq(1,3,ib)*qnet
            atoms%stress(2,3)=atoms%stress(2,3)+ann_arr%stresspq(2,3,ib)*qnet
            atoms%stress(3,3)=atoms%stress(3,3)+ann_arr%stresspq(3,3,ib)*qnet
        enddo
        do iat=1,atoms%nat
            atoms%fat(1,iat)=atoms%fat(1,iat)+ann_arr%fat_chi(1,iat)
            atoms%fat(2,iat)=atoms%fat(2,iat)+ann_arr%fat_chi(2,iat)
            atoms%fat(3,iat)=atoms%fat(3,iat)+ann_arr%fat_chi(3,iat)
        enddo
    endif
    call invertmat_alborz(atoms%cellvec,hinv)
    do i=1,3
    do j=1,3
        atoms%celldv(i,j)=vol*(atoms%stress(i,1)*hinv(j,1)+atoms%stress(i,2)*hinv(j,2)+atoms%stress(i,3)*hinv(j,3))
    enddo
    enddo
end subroutine cal_force_chi_part2
!*****************************************************************************************
subroutine repulsive_potential_cent(parini,atoms)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_linked_lists, only: typ_linked_lists
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    !local variables
    integer:: iat, jat, maincell_iat, maincell
    integer:: ip, jp, jpt, il, jl, iz, iy, ix, jx, jy, jz, iatp
    real(8):: cell(3), epot_rep, fx, fy, fz, ttt, a , b, c, d, g, h
    real(8):: rc, rcsq, dx, dy, dz, xiat, yiat, ziat, r, rsq
    real(8):: rt2, rtinv2, rtinv4, rtinv12, rcsqinv
    type(typ_linked_lists):: linked_lists
    !integer, save:: icall=0
    !icall=icall+1
    associate(rcmax=>linked_lists%rcut)
    !-------------------------------------------------------
    !{{B -> -126 A, C -> 560 A, D -> -945 A, G -> 720 A, H -> -210 A}}
    a=1/2.d2;
    b=-126.d0*a
    c=560.d0*a
    d=-945.d0*a
    g=720.d0*a
    h=-210.d0*a
    rc=5.d0
    !-------------------------------------------------------
    linked_lists%rcut=rc
    call linkedlists_init(parini,atoms,cell,linked_lists)
    epot_rep=0.d0
    rcsq=rc**2
    rcsqinv=1.d0/rcsq
    linked_lists%fat=0.d0
    include '../src/act1_cell_linkedlist.inc'
    do  iat=ip,il
        !qiat=linked_lists%qat(iat)
        xiat=linked_lists%rat(1,iat)
        yiat=linked_lists%rat(2,iat)
        ziat=linked_lists%rat(3,iat)
        jpt=linked_lists%prime(ix+linked_lists%limnbx(1,jy-iy,jz-iz),jy,jz)
        jp=(iat-ip+1)*((isign(1,ip-jpt)+1)/2)+jpt
        jl=linked_lists%last(ix+linked_lists%limnbx(2,jy-iy,jz-iz),jy,jz)
        maincell_iat=linked_lists%maincell(iat)
        iatp=linked_lists%perm(iat)
        do  jat=jp,jl
            dx=xiat-linked_lists%rat(1,jat)
            dy=yiat-linked_lists%rat(2,jat)
            dz=ziat-linked_lists%rat(3,jat)
            rsq=dx*dx+dy*dy+dz*dz
            maincell=maincell_iat+linked_lists%maincell(jat)
            if(rsq<rcsq .and. maincell >-1) then
                !---------------------------------
                rt2=rsq/rcsq
                rtinv2=1.d0/rt2
                rtinv4=rtinv2**2
                rtinv12=rtinv4**3
                epot_rep=epot_rep+a*rtinv12+(((b*rt2+c)*rt2+d)*rt2+g)*rt2+h
                ttt=-13.d0*a*rcsqinv*rtinv12*rtinv2+((8.d0*b*rcsqinv*rt2+6.d0*rcsqinv*c)*rt2+4.d0*rcsqinv*d)*rt2+2.d0*g*rcsqinv
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
end subroutine repulsive_potential_cent
!*****************************************************************************************
