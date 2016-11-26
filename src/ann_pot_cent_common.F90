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
