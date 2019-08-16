!*****************************************************************************************
subroutine symmetry_functions_driver_bond(parini,ann_arr,atoms,symfunc)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc
    use mod_atoms, only: typ_atoms
    use mod_linked_lists, only: typ_pia_arr !,typ_linked_lists
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(in):: atoms
    type(typ_symfunc), intent(inout):: symfunc
    !local variables
    integer:: ig
    integer:: iat !, jat, kat
    type(typ_pia_arr):: pia_arr
    integer:: isat, jsat, ksat, ib, ia, ibij, ibik, istat
    integer:: jat
    real(8):: rij, drij(3), fcij, fcdij
    real(8):: cutoff_function, cutoff_function_der
    external cutoff_function, cutoff_function_der
    associate(rc=>symfunc%linked_lists%rcut)
    symfunc%linked_lists%rcut=ann_arr%rcut
    symfunc%linked_lists%triplex=.true.
    call call_linkedlist(parini,atoms,.true.,symfunc%linked_lists,pia_arr)
    !write(*,*) 'HERE ',symfunc%linked_lists%maxbound_rad
    !stop
    !if(symfunc%linked_lists%maxbound_rad/=2) stop 'ERROR: correct next line'
    allocate(symfunc%y(ann_arr%ann(1)%nn(0),symfunc%linked_lists%maxbound_rad),stat=istat,source=0.d0)
    if(istat/=0) stop 'ERROR: unable to allocate array symfunc%y'
    !-------------------------------------------------------------------------------------
    associate(ng=>ann_arr%ann(1)%nn(0))
    !write(*,*) ng,atoms%nat,atoms%maxbound_rad,allocated(ann_arr%y0d)
    allocate(symfunc%y0d(ng,3,symfunc%linked_lists%maxbound_rad),stat=istat,source=0.d0)
    allocate(symfunc%y0d_bond(ng,symfunc%linked_lists%maxbound_rad),stat=istat,source=0.d0)
    !write(*,*) ng,atoms%nat,atoms%maxbound_rad,allocated(symfunc%y0d)
    if(istat/=0) stop 'ERROR: unable to allocate array symfunc%y0d.'
    allocate(symfunc%y0dr(ng,9,symfunc%linked_lists%maxbound_rad),stat=istat,source=0.d0)
    if(istat/=0) stop 'ERROR: unable to allocate array symfunc%y0dr.'
    !stop
    do ib=1,symfunc%linked_lists%maxbound_rad
        iat=symfunc%linked_lists%bound_rad(1,ib)
        jat=symfunc%linked_lists%bound_rad(2,ib)
        !write(*,*) 'BEFORE ',ib,iat,jat
        if(iat>jat) cycle !TO_BE_CORRECTED
        !write(*,*) 'AFTER ',ib,iat,jat
        isat=atoms%itypat(iat)
        jsat=atoms%itypat(symfunc%linked_lists%bound_rad(2,ib))
        pia_arr%pia(ib)%fc=cutoff_function(pia_arr%pia(ib)%r,rc)
        pia_arr%pia(ib)%fcd=cutoff_function_der(pia_arr%pia(ib)%r,rc)
        !rij=pia_arr%pia(ib)%r
        !drij=pia_arr%pia(ib)%dr
        !fcij=pia_arr%pia(ib)%fc
        !fcdij=pia_arr%pia(ib)%fcd
        !call symmetry_functions_g02_atom(ann_arr,pia_arr%pia(ib),ib,iat,isat,jsat)
        call symmetry_functions_g01_bond(ann_arr,ib,pia_arr%pia(ib),symfunc)
    enddo
    !do ia=1,symfunc%linked_lists%maxbound_ang
    !    ibij=symfunc%linked_lists%bound_ang(1,ia)
    !    ibik=symfunc%linked_lists%bound_ang(2,ia)
    !    if(symfunc%linked_lists%bound_rad(1,ibij)/=symfunc%linked_lists%bound_rad(1,ibik)) then
    !        stop 'ERROR: the centered atoms for two bonds is different.'
    !    endif
    !    iat=symfunc%linked_lists%bound_rad(1,ibij)
    !    isat=atoms%itypat(iat)
    !    jsat=atoms%itypat(symfunc%linked_lists%bound_rad(2,ibij))
    !    ksat=atoms%itypat(symfunc%linked_lists%bound_rad(2,ibik))
    !    call symmetry_functions_g05_atom(ann_arr,pia_arr%pia(ibij),pia_arr%pia(ibik),ibij,ibik,iat,isat,jsat,ksat)
    !    rij=pia_arr%pia(ibij)%r
    !    drij=pia_arr%pia(ibij)%dr
    !    fcij=

    !    call symmetry_functions_g02_bond(ann_arr,atoms,iat,jat,rij,drij,fcij,fcdij,rik,drik,fcik,fcdik,rjk,drjk,fcjk,fcdjk,isat)
    !enddo
    !-------------------------------------------------------------------------------------
    !if(parini%iverbose>=2) then
    !do iat=1,atoms%nat
    !do ig=1,ann_arr%ann(1)%nn(0)
    !    write(77,'(2i4,es24.15)') ig,iat,ann_arr%yall(ig,iat)
    !enddo
    !enddo
    !endif
   ! do iat=1,atoms%nat
   ! do jat=1,atoms%nat
   ! do ig=1,ann_arr%ann(1)%nn(0)
   !     write(88,'(2i4,3es18.9)') ig,jat,ann_arr%y0d(ig,1,jat,iat),ann_arr%y0d(ig,2,jat,iat),ann_arr%y0d(ig,3,jat,iat)
   !     write(99,'(2i4,9es18.9)') ig,jat,ann_arr%y0dr(ig,1,jat,iat),ann_arr%y0dr(ig,2,jat,iat),ann_arr%y0dr(ig,3,jat,iat) &
   !                                     ,ann_arr%y0dr(ig,4,jat,iat),ann_arr%y0dr(ig,5,jat,iat),ann_arr%y0dr(ig,6,jat,iat) &
   !                                     ,ann_arr%y0dr(ig,7,jat,iat),ann_arr%y0dr(ig,8,jat,iat),ann_arr%y0dr(ig,9,jat,iat)
   ! enddo
   ! enddo
   ! enddo
    end associate
    end associate
    !stop
end subroutine symmetry_functions_driver_bond
!**********************************************************************************************
subroutine symmetry_functions_driver_bond_tmp(ann_arr,atoms)
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms, get_rat
    implicit none
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(in):: atoms
    !local variables
    real(8):: drij(3), drik(3), drjk(3)
    integer:: i0, ig, isat
    integer:: iat, jat, kat, kat_maincell
    real(8):: fcij, fcik, fcdij, fcdik, fcjk, fcdjk
    real(8):: rij, rik, rjk, rc, en, eta1, eta2, zeta, vij, rs
    real(8):: eval(3)
    real(8):: cutoff_function, cutoff_function_der
    external cutoff_function, cutoff_function_der
    integer, parameter::lwork=100
    integer :: info
    real(8), dimension(lwork) :: work
    real(8), allocatable:: rat(:,:)
    do iat=1,atoms%nat
    isat=atoms%itypat(iat)
    rc=ann_arr%rcut
    rs=ann_arr%ann(isat)%g2rs(1)
    eta1=ann_arr%ann(isat)%g2eta(1)
    eta2=ann_arr%ann(isat)%g4eta(1)
    zeta=ann_arr%ann(isat)%g4zeta(1)
    allocate(rat(3,atoms%nat))
    call get_rat(atoms,rat)
    do jat=1,atoms%nat
        if(jat==iat) cycle
        !distance between bond atoms
        drij(1)=rat(1,jat)-rat(1,iat)
        drij(2)=rat(2,jat)-rat(2,iat)
        drij(3)=rat(3,jat)-rat(3,iat)
        rij=sqrt(drij(1)**2+drij(2)**2+drij(3)**2)
        !-----------------------------------------------------------------
        fcij=cutoff_function(rij,rc)
        fcdij=cutoff_function_der(rij,rc)
        vij=exp(-eta1*(rij-rs)**2)
        do kat=1,atoms%natim
            if(kat==iat .or. kat==jat) cycle
            drik(1)=atoms%ratim(1,kat)-rat(1,iat)
            drik(2)=atoms%ratim(2,kat)-rat(2,iat)
            drik(3)=atoms%ratim(3,kat)-rat(3,iat)
            drjk(1)=atoms%ratim(1,kat)-rat(1,jat)
            drjk(2)=atoms%ratim(2,kat)-rat(2,jat)
            drjk(3)=atoms%ratim(3,kat)-rat(3,jat)
            rik=sqrt(drik(1)**2+drik(2)**2+drik(3)**2)
            rjk=sqrt(drjk(1)**2+drjk(2)**2+drjk(3)**2)
            fcik=cutoff_function(rik,rc)
            fcdik=cutoff_function_der(rik,rc)
            fcjk=cutoff_function(rjk,rc)
            fcdjk=cutoff_function_der(rjk,rc)
            kat_maincell=mod(kat-1,atoms%nat)+1
            call symmetry_functions_g02_bond(ann_arr,iat,jat,rij,drij,fcij,fcdij,rik,drik,fcik,fcdik,rjk,drjk,fcjk,fcdjk)
            if(ann_arr%ann(isat)%ng4>0) then
            call symmetry_functions_g04_bond(ann_arr,iat,jat,rij,drij,fcij,fcdij,rik,drik,fcik,fcdik,rjk,drjk,fcjk,fcdjk)
            endif
        enddo
    enddo
    i0=ann_arr%ann(isat)%ng1+ann_arr%ann(isat)%ng2+ann_arr%ann(isat)%ng3+ann_arr%ann(isat)%ng4
    enddo
    deallocate(rat)
end subroutine symmetry_functions_driver_bond_tmp
!**********************************************************************************************
subroutine symmetry_functions_g01_bond(ann_arr,ib,pia,symfunc)
    use mod_linked_lists, only: typ_pia
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc
    implicit none
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_pia), intent(in):: pia
    integer, intent(in):: ib
    type(typ_symfunc), intent(inout):: symfunc
    real(8):: fcdij
    real(8):: rij
    real(8):: drij(3)
    !local variables
    integer:: kat, ig, i0
    real(8):: rs, rc, vij, eta
    real(8):: ttei, tt1i
    real(8):: ttjx, ttjy, ttjz, kat_maincell
    i0=0
    do ig=1,ann_arr%ann(1)%ng1
        i0=i0+1
        rs=ann_arr%ann(1)%g1rs(ig)
        eta=ann_arr%ann(1)%g1eta(ig) 
        vij=exp(-eta*(pia%r-rs)**2)
        ttei=vij*pia%fc
        tt1i=(-2.d0*eta*(pia%r-rs)*pia%fc+pia%fcd)*vij/pia%r
        symfunc%y0d_bond(i0,ib)=tt1i*pia%r
        ttjx=tt1i*pia%dr(1)
        ttjy=tt1i*pia%dr(2)
        ttjz=tt1i*pia%dr(3)
        !write(*,*) 'PIA', pia%r, pia%dr(3)
        symfunc%y0d(i0,1,ib)=symfunc%y0d(i0,1,ib)+ttjx
        symfunc%y0d(i0,2,ib)=symfunc%y0d(i0,2,ib)+ttjy
        symfunc%y0d(i0,3,ib)=symfunc%y0d(i0,3,ib)+ttjz
        symfunc%y(i0,ib)=symfunc%y(i0,ib)+ttei
    enddo
    !ann_arr%yall_bond(i0,iat,jat)=ann_arr%yall_bond(i0,iat,jat)*fcij*exp(-eta*rij**2)
end subroutine symmetry_functions_g01_bond
!******************************************************************************************
subroutine symmetry_functions_g02_bond(ann_arr,iat,jat,rij,drij,fcij,fcdij,rik,drik,fcik,fcdik,rjk,drjk,fcjk,fcdjk)
    use mod_ann, only: typ_ann_arr
    implicit none
    type(typ_ann_arr), intent(inout):: ann_arr
    integer, intent(in):: iat, jat
    real(8), intent(in):: fcij, fcik, fcjk, fcdij, fcdik, fcdjk
    real(8), intent(in):: rij, rik, rjk
    real(8), intent(in):: drij(3), drik(3), drjk(3)
    !local variables
    integer:: kat, ig, i0
    real(8):: rs, rc, vij, vik, vjk, eta
    real(8):: tte, tt1, tt2, tt1_x, tt1_y, tt1_z, tt2_x, tt2_y, tt2_z
    real(8):: ttx, tty, ttz, kat_maincell
    i0=ann_arr%ann(1)%ng1+ann_arr%ann(1)%ng2
    do ig=1,ann_arr%ann(1)%ng2
        i0=i0+1
        rs=ann_arr%ann(1)%g2rs(ig)
        eta=ann_arr%ann(1)%g2eta(ig) 
            vik=exp(-eta*(rik-rs)**2)
            vjk=exp(-eta*(rjk-rs)**2)
            tte=vik*fcik+vjk*fcjk        
            tt1=(-2.d0*eta*(rik-rs)*fcik+fcdik)*vik/rik
            tt1_x=tt1*drik(1)
            tt1_y=tt1*drik(2)
            tt1_z=tt1*drik(3)
            tt2=(-2.d0*eta*(rjk-rs)*fcjk+fcdjk)*vjk/rjk
            tt2_x=tt2*drjk(1)
            tt2_y=tt2*drjk(2)
            tt2_z=tt2*drjk(3)
            ttx=tt1_x+tt2_x
            tty=tt1_y+tt2_y
            ttz=tt1_z+tt2_z
            ann_arr%y0d_bond(i0,1,iat,jat)=ann_arr%y0d_bond(i0,1,iat,jat)+ttx
            ann_arr%y0d_bond(i0,2,iat,jat)=ann_arr%y0d_bond(i0,2,iat,jat)+tty
            ann_arr%y0d_bond(i0,3,iat,jat)=ann_arr%y0d_bond(i0,3,iat,jat)+ttz
            ann_arr%y0d_bond(i0,1,jat,iat)=ann_arr%y0d_bond(i0,1,jat,iat)-ttx
            ann_arr%y0d_bond(i0,2,jat,iat)=ann_arr%y0d_bond(i0,2,jat,iat)-tty
            ann_arr%y0d_bond(i0,3,jat,iat)=ann_arr%y0d_bond(i0,3,jat,iat)-ttz
            ann_arr%yall_bond(i0,iat,jat)=ann_arr%yall_bond(i0,iat,jat)+tte
    enddo
    ann_arr%yall_bond(i0,iat,jat)=ann_arr%yall_bond(i0,iat,jat)*fcij*exp(-eta*rij**2)
end subroutine symmetry_functions_g02_bond
!*****************************************************************************************
subroutine symmetry_functions_g04_bond(ann_arr,iat,jat,rij,drij,fcij,fcdij,rik,drik,fcik,fcdik,rjk,drjk,fcjk,fcdjk)
    use mod_ann, only: typ_ann_arr
    implicit none
    type(typ_ann_arr), intent(inout):: ann_arr
    integer, intent(in):: iat, jat
    real(8), intent(in):: fcij, fcik, fcjk, fcdij, fcdik, fcdjk
    real(8), intent(in):: rij, rik, rjk
    real(8), intent(in):: drij(3), drik(3), drjk(3)
    !local variables
    integer:: ig, i0
    real(8):: cos_theta, tte, ttx, tty, ttz, tt1, tt2
    i0=ann_arr%ann(1)%ng1+ann_arr%ann(1)%ng2+ann_arr%ann(1)%ng3
    do ig=1,ann_arr%ann(1)%ng4
        i0=i0+1
        !if(isat==jsat) then
            cos_theta=(drik(1)*drjk(1)+drik(2)*drjk(2)+drik(3)*drjk(3))/(rik*rjk)
        !else
        !    cos_theta=(drij(1)*drik(1)+drij(2)*drik(2)+drij(3)*drik(3))/(rij*rik)
        !endif
        tte=(1.d0+ann_arr%ann(1)%g4lambda(ig)*cos_theta)**ann_arr%ann(1)%g4zeta(ig)*exp(-ann_arr%ann(1)%g4eta(ig)*(rik**2+rjk**2))*fcik*fcjk
        ann_arr%yall_bond(i0,iat,jat)=ann_arr%yall_bond(i0,iat,jat)+tte
    enddo
     ann_arr%yall_bond(i0,iat,jat)=ann_arr%yall_bond(i0,iat,jat)*2.d0**(1.d0-ann_arr%ann(1)%g4zeta(1))*fcij*exp(-ann_arr%ann(1)%g4eta(1)*rij**2)
end subroutine symmetry_functions_g04_bond
!*****************************************************************************************
