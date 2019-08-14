!*****************************************************************************************
subroutine shortrange_init(atoms,shortrange,linked_lists,spline)
    use mod_atoms, only: typ_atoms
    use mod_shortrange, only: typ_shortrange
    use mod_linked_lists, only: typ_linked_lists
    use mod_spline, only: typ_spline
    implicit none
    type(typ_atoms), intent(in):: atoms
    type(typ_shortrange), intent(inout):: shortrange
    type(typ_linked_lists), intent(in):: linked_lists
    type(typ_spline), intent(inout):: spline
    !local variables
    integer:: istat
    shortrange%ntypinter=(atoms%ntypat**2+atoms%ntypat)/2
    allocate(spline%fsp(0:4,0:spline%nsp-1,shortrange%ntypinter),stat=istat)
    if(istat/=0) stop 'ERROR: failure allocating array fsp.'
    allocate(spline%fdsp(0:3,0:spline%nsp,shortrange%ntypinter),stat=istat)
    if(istat/=0) stop 'ERROR: failure allocating array fdsp.'
    call build_shortrange_spline(shortrange,spline,linked_lists%rcut,shortrange%alpha)
end subroutine shortrange_init
!*****************************************************************************************
subroutine shortrange_final(linked_lists,spline)
    use mod_linked_lists, only: typ_linked_lists
    use mod_spline, only: typ_spline
    implicit none
    type(typ_linked_lists), intent(in):: linked_lists
    type(typ_spline), intent(inout):: spline
    !local variables
    integer:: istat
    deallocate(spline%fsp,stat=istat)
    if(istat/=0) then
        write(*,*) 'ERROR: failure deallocating array fsp.'
    endif
    deallocate(spline%fdsp,stat=istat)
    if(istat/=0) then
        write(*,*) 'ERROR: failure deallocating array fdsp.'
    endif
end subroutine shortrange_final
!*****************************************************************************************
subroutine set_interaction(atoms,shortrange)
    use mod_atoms, only: typ_atoms
    use mod_shortrange, only: typ_shortrange
    implicit none
    type(typ_atoms), intent(inout):: atoms
    type(typ_shortrange), intent(inout):: shortrange
    !local variables
    integer:: ntypinter
    integer:: itypinter, itypat, jtypat !, nall, iall
    integer:: leni, lenj
    !character(6):: strint(10) !string of interaction
    !real(8):: parameters(5,10)
    !character(20):: strtmpi, strtmpj
    !character(4):: namatnamat1, namatnamat2
    shortrange%ntypinter=(atoms%ntypat**2+atoms%ntypat)/2
    allocate (shortrange%interaction(10,10),shortrange%qq(10))
    !nall=10
    !write(*,*) 'nall=',nall
    !call tosifumi_parameters(strint,parameters)
    itypinter=0
    do itypat=1,atoms%ntypat
        !strtmpi=atoms%stypat(itypat)
        !leni=len_trim(strtmpi)
        do jtypat=itypat,atoms%ntypat
            itypinter=itypinter+1
            shortrange%interaction(itypat,jtypat)=itypinter
            shortrange%interaction(jtypat,itypat)=itypinter
            !strtmpj=atoms%stypat(jtypat)
            !lenj=len_trim(strtmpj)
            !namatnamat1=strtmpi(1:leni)//strtmpj(1:lenj)
            !namatnamat2=strtmpj(1:lenj)//strtmpi(1:leni)
            !do iall=1,nall
            !    if(namatnamat1==strint(iall) .or. namatnamat2==strint(iall)) then
                shortrange%qq(itypinter)=atoms%qtypat(itypat)*atoms%qtypat(jtypat)
            !    endif
            !enddo
            !write(*,'(i3,f10.4,2(2x,a))') itypinter,tosifumi%ccc(itypinter),namatnamat1,namatnamat2
        enddo
    enddo
    !write(*,*) shortrange%interaction(1,1),shortrange%interaction(1,2),shortrange%interaction(1,3)
    !write(*,*) shortrange%interaction(2,1),shortrange%interaction(2,2),shortrange%interaction(2,3)
    !write(*,*) shortrange%interaction(3,1),shortrange%interaction(3,2),shortrange%interaction(3,3)
    !write(*,*)
    !write(*,*) shortrange%qq(1),shortrange%qq(2),shortrange%qq(3)
    !stop
end subroutine set_interaction
!*****************************************************************************************
subroutine cal_shortenergy(parini,shortrange,atoms,linked_lists,spline,alpha,cell,epot_short)
    use mod_parini, only: typ_parini
    use mod_shortrange, only: typ_shortrange
    use mod_linked_lists, only: typ_linked_lists
    use mod_spline, only: typ_spline
    use mod_atoms, only: typ_atoms, update_ratp
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_shortrange), intent(in):: shortrange
    type(typ_atoms), intent(inout):: atoms
    type(typ_linked_lists), intent(inout):: linked_lists
    type(typ_spline), intent(in):: spline
    real(8), intent(in):: alpha
    real(8), intent(out):: cell(3)
    real(8), intent(out):: epot_short !short range electrostatic energy
    !local variables
    real(8):: dx, dy, dz, r, rsq, xiat, yiat, ziat, alphainv, twosqrtinv
    real(8):: t, tt, tt1, tt2, tt3, ttt
    real(8):: rcutsq, fx, fy, fz, pi, hspinv, rhspinv, rinv, qiat, qiatjat, spf, spfd
    integer:: ip, jp, jpt, jl, il
    integer:: ix, iy, iz, jy, jz, iat, jat, ipat, isp, maincell, maincell_iat
    integer::itypinter
    pi=4.d0*atan(1.d0)
    hspinv=1.d0/spline%hsp
    call yaml_map('hsp in shortenergy',hspinv,fmt='(es22.14)')
    !write(*,*) 'inside shortenergy  hsp=',spline%hsp
    call linkedlists_init(parini,atoms,cell,linked_lists)
    !-------------------------------------------------------
    epot_short=0.d0
    alphainv=1.d0/alpha
    twosqrtinv=1.d0/sqrt(2.d0)
    linked_lists%fat=0.d0
    rcutsq=linked_lists%rcut**2
    call update_ratp(linked_lists%typ_atoms)
    include 'act1_cell_linkedlist.inc'
    do  iat=ip,il
        qiat=linked_lists%qat(iat)
        xiat=linked_lists%ratp(1,iat)
        yiat=linked_lists%ratp(2,iat)
        ziat=linked_lists%ratp(3,iat)
        jpt=linked_lists%prime(ix+linked_lists%limnbx(1,jy-iy,jz-iz),jy,jz)
        jp=(iat-ip+1)*((isign(1,ip-jpt)+1)/2)+jpt
        jl=linked_lists%last(ix+linked_lists%limnbx(2,jy-iy,jz-iz),jy,jz)
        maincell_iat=linked_lists%maincell(iat)
        do  jat=jp,jl
            dx=xiat-linked_lists%ratp(1,jat)
            dy=yiat-linked_lists%ratp(2,jat)
            dz=ziat-linked_lists%ratp(3,jat)
            rsq= dx*dx+dy*dy+dz*dz
            maincell=maincell_iat+linked_lists%maincell(jat)
            if(rsq<rcutsq .and. maincell >-1) then
                r=sqrt(rsq)
                itypinter=shortrange%interaction(linked_lists%itypat(iat),linked_lists%itypat(jat))
                qiatjat=linked_lists%qat(jat)*qiat
                !-----------------------------------
                !t=erfc(r*alphainv*twosqrtinv)
                !tt=linked_lists%qat(iat)*linked_lists%qat(jat)*t/r
                !epot_short=epot_short+tt !qat(iat)*qat(jat)*t/r
                !tt1=linked_lists%qat(iat)*linked_lists%qat(jat)
                !tt2=sqrt(2.d0/pi)*exp(-0.5d0*rsq*alphainv**2)*alphainv/rsq
                !tt3=t/(r*rsq)
                !ttt=tt1*(tt2+tt3)
                !fx=ttt*dx;fy=ttt*dy;fz=ttt*dz
                !-----------------------------------
                rinv=1.d0/r
                rhspinv=r*hspinv
                isp=floor(rhspinv)
                t=rhspinv-isp
                spf=spline%fsp(0,isp,itypinter)+ &
                   (spline%fsp(1,isp,itypinter)+ &
                   (spline%fsp(2,isp,itypinter)+ &
                   (spline%fsp(3,isp,itypinter)+ &
                    spline%fsp(4,isp,itypinter)*t)*t)*t)*t
                spfd=spline%fdsp(0,isp,itypinter)+(spline%fdsp(1,isp,itypinter)+ &
                    (spline%fdsp(2,isp,itypinter)+spline%fdsp(3,isp,itypinter)*t)*t)*t
                epot_short=epot_short+qiatjat*rinv-spf
                !write(*,*) linked_lists%itypat(iat),linked_lists%itypat(jat)
                !write(*,*) itypinter,shortrange%qq(itypinter)
                !write(51,'(f20.10)') spf !*qiatjat
                !stop
                ttt=(spfd+qiatjat*rinv*rinv)*rinv
                fx=ttt*dx;fy=ttt*dy;fz=ttt*dz
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
    include 'act2_cell_linkedlist.inc'
    !-------------------------------------------------------
    do iat=1,linked_lists%natim
        ipat=linked_lists%perm(iat)
        atoms%fat(1,ipat)=atoms%fat(1,ipat)+linked_lists%fat(1,iat)
        atoms%fat(2,ipat)=atoms%fat(2,ipat)+linked_lists%fat(2,iat)
        atoms%fat(3,ipat)=atoms%fat(3,ipat)+linked_lists%fat(3,iat)
    enddo
    !-------------------------------------------------------
    call linkedlists_final(linked_lists)
end subroutine cal_shortenergy
!*****************************************************************************************
