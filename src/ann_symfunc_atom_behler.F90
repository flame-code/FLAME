!*****************************************************************************************
subroutine symmetry_functions_driver(parini,ann_arr,atoms,symfunc)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc
    use mod_atoms, only: typ_atoms
    use mod_linked_lists, only: typ_pia_arr !,typ_linked_lists
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    !integer, intent(in):: iat
    type(typ_atoms), intent(in):: atoms
    type(typ_symfunc), intent(inout):: symfunc
    !local variables
    integer:: ig, i
    integer:: iat !, jat, kat
    type(typ_pia_arr):: pia_arr
    !real(8):: rc, en
    !real(8):: eval(3)
    !integer, parameter::lwork=100
    !integer :: info
    !real(8), dimension(lwork) :: work
    integer:: isat, jsat, ksat, ib, ia, ibij, ibik, istat
    !type(typ_linked_lists):: linked_lists
    call f_routine(id='symmetry_functions_driver')
    associate(rc=>symfunc%linked_lists%rcut)
    symfunc%linked_lists%rcut=ann_arr%rcut
    symfunc%linked_lists%triplex=.true.
    call call_linkedlist(parini,atoms,.true.,symfunc%linked_lists,pia_arr)
    !-------------------------------------------------------------------------------------
    associate(ng=>ann_arr%ann(1)%nn(0))
    !write(*,*) ng,atoms%nat,atoms%maxbound_rad,allocated(ann_arr%y0d)
    !allocate(ann_arr%y0d(ng,3,symfunc%linked_lists%maxbound_rad),stat=istat,source=0.d0)
    !write(*,*) ng,atoms%nat,atoms%maxbound_rad,allocated(ann_arr%y0d)
    !if(istat/=0) stop 'ERROR: unable to allocate array ann_arr%y0d.'
    !allocate(ann_arr%y0dr(ng,9,symfunc%linked_lists%maxbound_rad),stat=istat,source=0.d0)
    !if(istat/=0) stop 'ERROR: unable to allocate array ann_arr%y0dr.'
    symfunc%y=f_malloc0((/1.to.ng,1.to.atoms%nat/),id='symfunc%y')
    symfunc%y0d=f_malloc0((/1.to.ng,1.to.3,1.to.symfunc%linked_lists%maxbound_rad/),id='symfunc%y0d')
    symfunc%y0dr=f_malloc0((/1.to.ng,1.to.9,1.to.symfunc%linked_lists%maxbound_rad/),id='symfunc%y0dr')

    !stop
    do ib=1,symfunc%linked_lists%maxbound_rad
        iat=symfunc%linked_lists%bound_rad(1,ib)
        isat=atoms%itypat(iat)
        jsat=atoms%itypat(symfunc%linked_lists%bound_rad(2,ib))
        pia_arr%pia(ib)%fc=cutoff_function(pia_arr%pia(ib)%r,rc)
        pia_arr%pia(ib)%fcd=cutoff_function_der(pia_arr%pia(ib)%r,rc)
        call symmetry_functions_g02_atom(ann_arr,pia_arr%pia(ib),ib,iat,isat,jsat,symfunc)
    enddo
    do ia=1,symfunc%linked_lists%maxbound_ang
        ibij=symfunc%linked_lists%bound_ang(1,ia)
        ibik=symfunc%linked_lists%bound_ang(2,ia)
        if(symfunc%linked_lists%bound_rad(1,ibij)/=symfunc%linked_lists%bound_rad(1,ibik)) then
            stop 'ERROR: the centered atoms for two bonds is different.'
        endif
        iat=symfunc%linked_lists%bound_rad(1,ibij)
        isat=atoms%itypat(iat)
        jsat=atoms%itypat(symfunc%linked_lists%bound_rad(2,ibij))
        ksat=atoms%itypat(symfunc%linked_lists%bound_rad(2,ibik))
        call symmetry_functions_g05_atom(ann_arr,pia_arr%pia(ibij),pia_arr%pia(ibik),ibij,ibik,iat,isat,jsat,ksat,symfunc)
    enddo
    !-------------------------------------------------------------------------------------
    if(parini%iverbose>2) then
    do iat=1,atoms%nat
    do ig=1,ann_arr%ann(1)%nn(0)
        write(77,'(2i4,es24.15)') ig,iat,symfunc%y(ig,iat)
    enddo
    enddo
    endif
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
    deallocate(pia_arr%pia)
    call f_release_routine()
    !stop
end subroutine symmetry_functions_driver
!*****************************************************************************************
subroutine symmetry_functions_g02_atom(ann_arr,pia,ib,iat,isat,jsat,symfunc)
    use mod_interface
    use mod_ann, only: typ_ann_arr, typ_symfunc
    use mod_linked_lists, only: typ_pia
    implicit none
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_pia), intent(in):: pia
    integer, intent(in):: ib, iat, isat, jsat
    type(typ_symfunc), intent(inout):: symfunc
    !local variables
    integer:: i0, ig
    real(8):: ttei,ttej, ttix,ttiy,ttiz,ttjx,ttjy,ttjz,tt1i,tt1j
    real(8):: etai,etaj, rs, vi,vj
    i0=ann_arr%ann(isat)%ng1
    do ig=1,ann_arr%ann(isat)%ng2
        i0=i0+1
        if((ann_arr%ann(isat)%g2i(ig)/=0).and.(.not.(jsat==ann_arr%ann(isat)%g2i(ig)))) cycle
        rs=ann_arr%ann(jsat)%g2rs(ig)
        !The central atom is i:
        etaj=ann_arr%ann(jsat)%g2eta(ig)
        vi=exp(-etaj*(pia%r-rs)**2)
        ttei=vi*pia%fc
        tt1i=(-2.d0*etaj*(pia%r-rs)*pia%fc+pia%fcd)*vi/pia%r
        ttjx=tt1i*pia%dr(1)
        ttjy=tt1i*pia%dr(2)
        ttjz=tt1i*pia%dr(3)
        symfunc%y0d(i0,1,ib)=ttjx
        symfunc%y0d(i0,2,ib)=ttjy
        symfunc%y0d(i0,3,ib)=ttjz
        symfunc%y(i0,iat)=symfunc%y(i0,iat)+ttei
        
        !The lines we need to calculate the radial part of the stress tensor(without volume):
        symfunc%y0dr(i0,1,ib)=symfunc%y0dr(i0,1,ib)+ttjx*pia%dr(1) !sigma(1,1) 
        symfunc%y0dr(i0,2,ib)=symfunc%y0dr(i0,2,ib)+ttjx*pia%dr(2) !sigma(1,2)
        symfunc%y0dr(i0,3,ib)=symfunc%y0dr(i0,3,ib)+ttjx*pia%dr(3) !sigma(1,3)
                                            
        symfunc%y0dr(i0,4,ib)=symfunc%y0dr(i0,4,ib)+ttjy*pia%dr(1) !sigma(2,1) 
        symfunc%y0dr(i0,5,ib)=symfunc%y0dr(i0,5,ib)+ttjy*pia%dr(2) !sigma(2,2)
        symfunc%y0dr(i0,6,ib)=symfunc%y0dr(i0,6,ib)+ttjy*pia%dr(3) !sigma(2,3)
                                            
        symfunc%y0dr(i0,7,ib)=symfunc%y0dr(i0,7,ib)+ttjz*pia%dr(1) !sigma(3,1) 
        symfunc%y0dr(i0,8,ib)=symfunc%y0dr(i0,8,ib)+ttjz*pia%dr(2) !sigma(3,2)
        symfunc%y0dr(i0,9,ib)=symfunc%y0dr(i0,9,ib)+ttjz*pia%dr(3) !sigma(3,3)
    enddo
end subroutine symmetry_functions_g02_atom
!*****************************************************************************************
subroutine symmetry_functions_g04_atom(ann_arr,isat,iat,jsat,jat_maincell,ksat,kat_maincell,rij,rik,rjk,drij,drik,drjk,fcij,fcdij,fcik,fcdik,fcjk,fcdjk,symfunc)
    use mod_interface
    use mod_ann, only: typ_ann_arr, typ_symfunc
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_ann_arr), intent(inout):: ann_arr
    integer, intent(in):: isat, iat, jsat, jat_maincell, ksat, kat_maincell
    real(8), intent(in):: rij, rik, rjk, drij(3), drik(3), drjk(3), fcij, fcdij, fcik, fcdik, fcjk, fcdjk
    type(typ_symfunc), intent(inout):: symfunc
    !local variables
    integer:: i0, ig
    real(8):: tte, ttjx, ttjy, ttjz, ttkx, ttky, ttkz, tt1, tt2, ttj, ttk, tt4, tt5, tt6, tt7
    real(8):: ss1, ss2, ss3, ss4, ss5
    real(8):: cos_theta, uijk, vijk
    real(8):: zeta, alam, etaj, etak, etai, zzz
    !I don't understand the following line
    i0=ann_arr%ann(isat)%ng1+ann_arr%ann(isat)%ng2+ann_arr%ann(isat)%ng3
    do ig=1,ann_arr%ann(isat)%ng4
        i0=i0+1
        zeta=ann_arr%ann(ksat)%g4zeta(ig)
        alam=ann_arr%ann(ksat)%g4lambda(ig)
        !etaj=ann_arr%ann(jsat)%g4eta(ig)
        etak=ann_arr%ann(ksat)%g4eta(ig)
        zzz=2.d0**(1.d0-zeta)
        cos_theta=(drij(1)*drik(1)+drij(2)*drik(2)+drij(3)*drik(3))/(rij*rik)
        uijk=(1.d0+alam*cos_theta)**zeta
        vijk=exp(-(etak*(rij**2+rik**2+rjk**2)))
        tte=zzz*uijk*vijk*fcij*fcik*fcjk
        ss1=vijk*fcij*fcik*fcjk
        ss2=uijk*fcij*fcik*fcjk
        ss3=uijk*vijk*fcik*fcjk
        ss4=uijk*vijk*fcij*fcjk
        ss5=uijk*vijk*fcij*fcik

        tt1=zeta*alam*(1.d0+alam*cos_theta)**(zeta-1)*ss1
        tt2=tt1/(rij*rik)
        ttj=0.d0               !HERE etaj must be checked
        ttk=4.d0*etak*vijk*ss2 !HERE etak must be checked
        tt4=-tt1*cos_theta/rij**2-ttj+ss3*fcdij/rij
        tt6=-ss3*fcdjk/rjk
        tt5=-tt1*cos_theta/rik**2-ttk+ss4*fcdik/rik
        tt7=-ss5*fcdjk/rjk


        ttjx=zzz*(tt2*drik(1)+tt4*drij(1)+tt6*drjk(1))
        ttjy=zzz*(tt2*drik(2)+tt4*drij(2)+tt6*drjk(2))
        ttjz=zzz*(tt2*drik(3)+tt4*drij(3)+tt6*drjk(3))

        ttkx=zzz*(tt2*drij(1)+tt5*drik(1)+tt7*drjk(1))
        ttky=zzz*(tt2*drij(2)+tt5*drik(2)+tt7*drjk(2))
        ttkz=zzz*(tt2*drij(3)+tt5*drik(3)+tt7*drjk(3))

        !ann_arr%y0d(i0,1,iat         )=ann_arr%y0d(i0,1,iat         )-ttjx-ttkx
        !ann_arr%y0d(i0,2,iat         )=ann_arr%y0d(i0,2,iat         )-ttjy-ttky
        !ann_arr%y0d(i0,3,iat         )=ann_arr%y0d(i0,3,iat         )-ttjz-ttkz
        !ann_arr%y0d(i0,1,jat_maincell)=ann_arr%y0d(i0,1,jat_maincell)+ttjx
        !ann_arr%y0d(i0,2,jat_maincell)=ann_arr%y0d(i0,2,jat_maincell)+ttjy
        !ann_arr%y0d(i0,3,jat_maincell)=ann_arr%y0d(i0,3,jat_maincell)+ttjz
        !ann_arr%y0d(i0,1,kat_maincell)=ann_arr%y0d(i0,1,kat_maincell)+ttkx
        !ann_arr%y0d(i0,2,kat_maincell)=ann_arr%y0d(i0,2,kat_maincell)+ttky
        !ann_arr%y0d(i0,3,kat_maincell)=ann_arr%y0d(i0,3,kat_maincell)+ttkz
        symfunc%y(i0,iat)=symfunc%y(i0,iat)+tte
    enddo
end subroutine symmetry_functions_g04_atom
!*****************************************************************************************
subroutine symmetry_functions_g05_atom(ann_arr,piaij,piaik,ibij,ibik,iat,isat,jsat,ksat,symfunc)
    use mod_interface
    use mod_ann, only: typ_ann_arr, typ_symfunc
    use mod_linked_lists, only: typ_pia
    implicit none
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_pia), intent(in):: piaij, piaik
    integer, intent(in):: ibij, ibik, isat, iat, jsat, ksat
    type(typ_symfunc), intent(inout):: symfunc
    !local variables
    integer:: i0, ig, ii1, ii2
    real(8):: ttei, ttek, ttix, ttiy, ttiz, ttej, ttjx, ttjy, ttjz, ttkx, ttky, ttkz, tt1, tt2, tti, ttj, ttk, tt4, tt5
    real(8):: ss1, ss2, ss3, ss4
    real(8):: cos_theta_i, cos_theta_k, cos_theta_j, ui, uk, uj, vi, vk, vj
    real(8):: zeta, alam, etai, etaj, etak, zzz
    i0=ann_arr%ann(isat)%ng1+ann_arr%ann(isat)%ng2+ann_arr%ann(isat)%ng3+ann_arr%ann(isat)%ng4
    do ig=1,ann_arr%ann(isat)%ng5
        i0=i0+1
        ii1=ann_arr%ann(isat)%g5i(1,ig)+ann_arr%ann(isat)%g5i(2,ig)
        ii2=abs(ann_arr%ann(isat)%g5i(1,ig)-ann_arr%ann(isat)%g5i(2,ig))
        if((ann_arr%ann(isat)%g5i(1,ig)/=0).and. (.not.((jsat+ksat)==ii1 .and. abs(jsat-ksat)==ii2))) cycle
        zeta=ann_arr%ann(ksat)%g5zeta(ig)
        alam=ann_arr%ann(ksat)%g5lambda(ig)
        etai=ann_arr%ann(isat)%g5eta(ig)
        etaj=ann_arr%ann(jsat)%g5eta(ig)
        etak=ann_arr%ann(ksat)%g5eta(ig)
        zzz=2.d0**(1.d0-zeta)

        cos_theta_i=(piaij%dr(1)*piaik%dr(1)+piaij%dr(2)*piaik%dr(2)+piaij%dr(3)*piaik%dr(3))/(piaij%r*piaik%r)
        ui=(1.d0+alam*cos_theta_i)**zeta
        vi=exp(-(etaj*piaij%r**2+etak*piaik%r**2))
        ttei=zzz*ui*vi*piaij%fc*piaik%fc
        ss1=vi*piaij%fc*piaik%fc
        ss2=ui*piaij%fc*piaik%fc
        ss3=ui*vi*piaik%fc
        ss4=ui*vi*piaij%fc
        tt1=zeta*alam*(1.d0+alam*cos_theta_i)**(zeta-1)*ss1
        tt2=tt1/(piaij%r*piaik%r)
        ttj=2.d0*etaj*vi*ss2 !HERE etaj must be checked
        ttk=2.d0*etak*vi*ss2 !HERE etaj must be checked
        tt4=-tt1*cos_theta_i/piaij%r**2-ttj+ss3*piaij%fcd/piaij%r
        tt5=-tt1*cos_theta_i/piaik%r**2-ttk+ss4*piaik%fcd/piaik%r
        ttjx=zzz*(tt2*piaik%dr(1)+tt4*piaij%dr(1))
        ttjy=zzz*(tt2*piaik%dr(2)+tt4*piaij%dr(2))
        ttjz=zzz*(tt2*piaik%dr(3)+tt4*piaij%dr(3))
        ttkx=zzz*(tt2*piaij%dr(1)+tt5*piaik%dr(1))
        ttky=zzz*(tt2*piaij%dr(2)+tt5*piaik%dr(2))
        ttkz=zzz*(tt2*piaij%dr(3)+tt5*piaik%dr(3))
        symfunc%y0d(i0,1,ibij)=symfunc%y0d(i0,1,ibij)+ttjx
        symfunc%y0d(i0,2,ibij)=symfunc%y0d(i0,2,ibij)+ttjy
        symfunc%y0d(i0,3,ibij)=symfunc%y0d(i0,3,ibij)+ttjz
        symfunc%y0d(i0,1,ibik)=symfunc%y0d(i0,1,ibik)+ttkx
        symfunc%y0d(i0,2,ibik)=symfunc%y0d(i0,2,ibik)+ttky
        symfunc%y0d(i0,3,ibik)=symfunc%y0d(i0,3,ibik)+ttkz
        symfunc%y(i0,iat)=symfunc%y(i0,iat)+ttei

        !The lines we need to calculate the angular part of the stress tensor(without volume):
        symfunc%y0dr(i0,1,ibij)=symfunc%y0dr(i0,1,ibij)+(1.d0*ttjx*piaij%dr(1))  !sigma(1,1)
        symfunc%y0dr(i0,2,ibij)=symfunc%y0dr(i0,2,ibij)+(1.d0*ttjx*piaij%dr(2))  !sigma(1,2)
        symfunc%y0dr(i0,3,ibij)=symfunc%y0dr(i0,3,ibij)+(1.d0*ttjx*piaij%dr(3))  !sigma(1,3)
        
        symfunc%y0dr(i0,4,ibij)=symfunc%y0dr(i0,4,ibij)+(1.d0*ttjy*piaij%dr(1))  !sigma(2,1)
        symfunc%y0dr(i0,5,ibij)=symfunc%y0dr(i0,5,ibij)+(1.d0*ttjy*piaij%dr(2))  !sigma(2,2)
        symfunc%y0dr(i0,6,ibij)=symfunc%y0dr(i0,6,ibij)+(1.d0*ttjy*piaij%dr(3))  !sigma(2,3)
        
        symfunc%y0dr(i0,7,ibij)=symfunc%y0dr(i0,7,ibij)+(1.d0*ttjz*piaij%dr(1))  !sigma(3,1)
        symfunc%y0dr(i0,8,ibij)=symfunc%y0dr(i0,8,ibij)+(1.d0*ttjz*piaij%dr(2))  !sigma(3,2)
        symfunc%y0dr(i0,9,ibij)=symfunc%y0dr(i0,9,ibij)+(1.d0*ttjz*piaij%dr(3))  !sigma(3,3)
        
        symfunc%y0dr(i0,1,ibik)=symfunc%y0dr(i0,1,ibik)+(1.d0*ttkx*piaik%dr(1))  !sigma(1,1)
        symfunc%y0dr(i0,2,ibik)=symfunc%y0dr(i0,2,ibik)+(1.d0*ttkx*piaik%dr(2))  !sigma(1,2)
        symfunc%y0dr(i0,3,ibik)=symfunc%y0dr(i0,3,ibik)+(1.d0*ttkx*piaik%dr(3))  !sigma(1,3)
        
        symfunc%y0dr(i0,4,ibik)=symfunc%y0dr(i0,4,ibik)+(1.d0*ttky*piaik%dr(1))  !sigma(2,1)
        symfunc%y0dr(i0,5,ibik)=symfunc%y0dr(i0,5,ibik)+(1.d0*ttky*piaik%dr(2))  !sigma(2,2)
        symfunc%y0dr(i0,6,ibik)=symfunc%y0dr(i0,6,ibik)+(1.d0*ttky*piaik%dr(3))  !sigma(2,3)
        
        symfunc%y0dr(i0,7,ibik)=symfunc%y0dr(i0,7,ibik)+(1.d0*ttkz*piaik%dr(1))  !sigma(3,1)
        symfunc%y0dr(i0,8,ibik)=symfunc%y0dr(i0,8,ibik)+(1.d0*ttkz*piaik%dr(2))  !sigma(3,2)
        symfunc%y0dr(i0,9,ibik)=symfunc%y0dr(i0,9,ibik)+(1.d0*ttkz*piaik%dr(3))  !sigma(3,3)
    enddo
end subroutine symmetry_functions_g05_atom
!*****************************************************************************************
subroutine symmetry_functions_g06_atom(ann,iat,jat_maincell,r,dr,fc,fcd)
    use mod_interface
    use mod_ann, only: typ_ann
    implicit none
    type(typ_ann), intent(inout):: ann
    integer, intent(in):: iat, jat_maincell
    real(8), intent(in):: r, dr(3), fc, fcd
    !local variables
    integer:: ig, igt
    real(8) :: rcov, tt
    stop 'ERROR: symmetry functions G6: not ready yet'
    do ig=1,ann%ng6,3
        igt=(ig-1)/3+1
        rcov=1.d0
        !rcov=1.36d0
        !rcov=ann%g6eta(igt)
        !fc=exp(-r/(4.d0*rcov))
        !fc=amass(jat)*exp(-r/(4.d0*rcov))
        !fc=amass(jat)*exp(-r**2/(6.d0*rcov)**2)
        !tt=fc*exp(-ann%g6eta(igt)*(r)**2)
        tt=fc*exp(-1.d0*ann%g6eta(igt)*r**1)
        ann%teneria(1,1,igt)=ann%teneria(1,1,igt)+tt*(dr(2)*dr(2)+dr(3)*dr(3)+2.d0*(rcov*0.5d0)**2)
        ann%teneria(2,2,igt)=ann%teneria(2,2,igt)+tt*(dr(1)*dr(1)+dr(3)*dr(3)+2.d0*(rcov*0.5d0)**2)
        ann%teneria(3,3,igt)=ann%teneria(3,3,igt)+tt*(dr(1)*dr(1)+dr(2)*dr(2)+2.d0*(rcov*0.5d0)**2)
        ann%teneria(2,1,igt)=ann%teneria(2,1,igt)-tt*(dr(1)*dr(2))
        ann%teneria(3,1,igt)=ann%teneria(3,1,igt)-tt*(dr(1)*dr(3))
        ann%teneria(3,2,igt)=ann%teneria(3,2,igt)-tt*(dr(2)*dr(3))
        !write(*,*) fc
    enddo
end subroutine symmetry_functions_g06_atom
!*****************************************************************************************
subroutine symmetry_functions_g02(ann,iat,atoms,i0)
    use mod_interface
    use mod_ann, only: typ_ann
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_ann), intent(inout):: ann
    integer, intent(in):: iat
    type(typ_atoms), intent(in):: atoms
    integer, intent(inout):: i0
    !local variables
    integer:: jat, ig, jat_maincell
    real(8):: tte, ttx, tty, ttz, tt1
    real(8):: fc, fcd, r, dx, dy, dz
    real(8):: eta, rc, rs, vij
    stop 'ERROR: this routine does not work!'
    do ig=1,ann%ng2
        i0=i0+1
        rs=ann%g2rs(ig)
        eta=ann%g2eta(ig)
!HERE        rc=ann%g2rc(ig)
        tte=0.d0
        do jat=1,atoms%natim
            if(jat==iat) cycle
            dx=atoms%ratim(1,jat)-atoms%rat(1,iat)
            dy=atoms%ratim(2,jat)-atoms%rat(2,iat)
            dz=atoms%ratim(3,jat)-atoms%rat(3,iat)
            r=sqrt(dx**2+dy**2+dz**2)
            fc=cutoff_function(r,rc)
            fcd=cutoff_function_der(r,rc)
            vij=exp(-eta*(r-rs)**2)
            tte=tte+vij*fc
            tt1=(-2.d0*eta*(r-rs)*fc+fcd)*vij/r
            ttx=tt1*dx
            tty=tt1*dy
            ttz=tt1*dz
            jat_maincell=mod(jat-1,atoms%nat)+1
            !following lines commented due to moving yall
            !ann%y0d(i0,1,jat_maincell)=ann%y0d(i0,1,jat_maincell)+ttx
            !ann%y0d(i0,2,jat_maincell)=ann%y0d(i0,2,jat_maincell)+tty
            !ann%y0d(i0,3,jat_maincell)=ann%y0d(i0,3,jat_maincell)+ttz
            !ann%y0d(i0,1,iat)=ann%y0d(i0,1,iat)-ttx
            !ann%y0d(i0,2,iat)=ann%y0d(i0,2,iat)-tty
            !ann%y0d(i0,3,iat)=ann%y0d(i0,3,iat)-ttz
        enddo
        !ann%y0d(i0,1:3,1:atoms%nat)=ann%y0d(i0,1:3,1:atoms%nat)*two_over_gdiff/(log(10.d0)*(tte+1.d0))
        !tte=log10(tte+1.d0)
        !ann%yall(i0,iat)=tte !commented due to moving yall
    enddo
end subroutine symmetry_functions_g02
!*****************************************************************************************
subroutine symmetry_functions_g04(ann,iat,atoms,i0)
    use mod_interface
    use mod_ann, only: typ_ann
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_ann), intent(inout):: ann
    integer, intent(in):: iat
    type(typ_atoms), intent(in):: atoms
    integer, intent(inout):: i0
    !local variables
    integer:: jat, kat, ig
    real(8):: cos_theta, tte, ttx, tty, ttz, tt1, tt2
    real(8):: fcij, fcik, fcjk
    real(8):: rij, rik, rjk
    real(8):: dxij, dyij, dzij, dxik, dyik, dzik, dxjk, dyjk, dzjk
    stop 'ERROR: this routine does not work!'
    do ig=1,ann%ng4
        i0=i0+1
        tte=0.d0
        do jat=1,atoms%natim
            if(jat==iat) cycle
            dxij=atoms%ratim(1,jat)-atoms%rat(1,iat)
            dyij=atoms%ratim(2,jat)-atoms%rat(2,iat)
            dzij=atoms%ratim(3,jat)-atoms%rat(3,iat)
            rij=sqrt(dxij**2+dyij**2+dzij**2)
!HERE            fcij=cutoff_function(rij,ann%g4rc(ig))
            do kat=jat+1,atoms%natim
                if(kat==iat) cycle
                dxik=atoms%ratim(1,kat)-atoms%rat(1,iat)
                dyik=atoms%ratim(2,kat)-atoms%rat(2,iat)
                dzik=atoms%ratim(3,kat)-atoms%rat(3,iat)
                rik=sqrt(dxik**2+dyik**2+dzik**2)
!HERE                fcik=cutoff_function(rik,ann%g4rc(ig))
                dxjk=atoms%ratim(1,kat)-atoms%ratim(1,jat)
                dyjk=atoms%ratim(2,kat)-atoms%ratim(2,jat)
                dzjk=atoms%ratim(3,kat)-atoms%ratim(3,jat)
                rjk=sqrt(dxjk**2+dyjk**2+dzjk**2)
!HERE                fcjk=cutoff_function(rjk,ann%g4rc(ig))
                cos_theta=(dxij*dxik+dyij*dyik+dzij*dzik)/(rij*rik)
                tte=tte+(1.d0+ann%g4lambda(ig)*cos_theta)**ann%g4zeta(ig)*exp(-ann%g4eta(ig)*(rij**2+rik**2+rjk**2))*fcij*fcik*fcjk
            enddo
        enddo
        tte=tte*2.d0**(1.d0-ann%g4zeta(ig))
        !tte=log10(1.d0+tte)
        !ann%yall(i0,iat)=tte !commented due to moving yall
    enddo
end subroutine symmetry_functions_g04
!*****************************************************************************************
subroutine symmetry_functions_g05(ann,iat,atoms,i0)
    use mod_interface
    use mod_ann, only: typ_ann
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_ann), intent(inout):: ann
    integer, intent(in):: iat
    type(typ_atoms), intent(in):: atoms
    integer, intent(inout):: i0
    !local variables
    integer:: jat, kat, ig, jat_maincell, kat_maincell
    real(8):: tte, ttjx, ttjy, ttjz, ttkx, ttky, ttkz, tt1, tt2, tt3, tt4, tt5
    real(8):: ss1, ss2, ss3, ss4
    real(8):: cos_theta, fcij, fcik, fcdij, fcdik, uijk, vijk
    real(8):: rij, rik, rjk
    real(8):: dxij, dyij, dzij, dxik, dyik, dzik
    real(8):: zeta, alam, eta, rc, zzz
    stop 'ERROR: this routine does not work!'
    do ig=1,ann%ng5
        i0=i0+1
        zeta=ann%g5zeta(ig)
        alam=ann%g5lambda(ig)
        eta=ann%g5eta(ig)
!HERE        rc=ann%g5rc(ig)
        zzz=2.d0**(1.d0-zeta)
        tte=0.d0
        do jat=1,atoms%natim
            if(jat==iat) cycle
            dxij=atoms%ratim(1,jat)-atoms%rat(1,iat)
            dyij=atoms%ratim(2,jat)-atoms%rat(2,iat)
            dzij=atoms%ratim(3,jat)-atoms%rat(3,iat)
            rij=sqrt(dxij**2+dyij**2+dzij**2)
            fcij=cutoff_function(rij,rc)
            fcdij=cutoff_function_der(rij,rc)
            do kat=jat+1,atoms%natim
                if(kat==iat) cycle
                dxik=atoms%ratim(1,kat)-atoms%rat(1,iat)
                dyik=atoms%ratim(2,kat)-atoms%rat(2,iat)
                dzik=atoms%ratim(3,kat)-atoms%rat(3,iat)
                rik=sqrt(dxik**2+dyik**2+dzik**2)
                fcik=cutoff_function(rik,rc)
                fcdik=cutoff_function_der(rik,rc)
                cos_theta=(dxij*dxik+dyij*dyik+dzij*dzik)/(rij*rik)
                uijk=(1.d0+alam*cos_theta)**zeta
                vijk=exp(-eta*(rij**2+rik**2))
                tte=tte+uijk*vijk*fcij*fcik
                ss1=vijk*fcij*fcik
                ss2=uijk*fcij*fcik
                ss3=uijk*vijk*fcik
                ss4=uijk*vijk*fcij

                tt1=zeta*alam*(1.d0+alam*cos_theta)**(zeta-1)*ss1
                tt2=tt1/(rij*rik)
                tt3=2.d0*eta*vijk*ss2
                tt4=-tt1*cos_theta/rij**2-tt3+ss3*fcdij/rij
                tt5=-tt1*cos_theta/rik**2-tt3+ss4*fcdik/rik

                ttjx=tt2*dxik+tt4*dxij
                ttjy=tt2*dyik+tt4*dyij
                ttjz=tt2*dzik+tt4*dzij

                ttkx=tt2*dxij+tt5*dxik
                ttky=tt2*dyij+tt5*dyik
                ttkz=tt2*dzij+tt5*dzik

                jat_maincell=mod(jat-1,atoms%nat)+1
                kat_maincell=mod(kat-1,atoms%nat)+1
                !following lines commented due to moving yall
                !ann%y0d(i0,1,iat)=ann%y0d(i0,1,iat)-ttjx-ttkx
                !ann%y0d(i0,2,iat)=ann%y0d(i0,2,iat)-ttjy-ttky
                !ann%y0d(i0,3,iat)=ann%y0d(i0,3,iat)-ttjz-ttkz
                !ann%y0d(i0,1,jat_maincell)=ann%y0d(i0,1,jat_maincell)+ttjx
                !ann%y0d(i0,2,jat_maincell)=ann%y0d(i0,2,jat_maincell)+ttjy
                !ann%y0d(i0,3,jat_maincell)=ann%y0d(i0,3,jat_maincell)+ttjz
                !ann%y0d(i0,1,kat_maincell)=ann%y0d(i0,1,kat_maincell)+ttkx
                !ann%y0d(i0,2,kat_maincell)=ann%y0d(i0,2,kat_maincell)+ttky
                !ann%y0d(i0,3,kat_maincell)=ann%y0d(i0,3,kat_maincell)+ttkz
            enddo
        enddo
        !following line commented due to moving yall
        !ann%y0d(i0,1:3,1:atoms%nat)=zzz*ann%y0d(i0,1:3,1:atoms%nat)
        tte=tte*zzz
        !ann%yall(i0,iat)=tte !commented due to moving yall
    enddo
end subroutine symmetry_functions_g05
!*****************************************************************************************
subroutine symmetry_functions_g06(ann,iat,atoms,i0)
    use mod_interface
    use mod_ann, only: typ_ann
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_ann), intent(inout):: ann
    integer, intent(in):: iat
    type(typ_atoms), intent(in):: atoms
    integer, intent(inout):: i0
    !local variables
    integer:: jat, ig
    real(8):: fc, r, dx, dy, dz
    real(8):: teneria(1:3,1:3), eval(3)
    integer, parameter::lwork=100
    integer :: info
    real(8) :: rcov
    real(8), dimension(lwork) :: work
    stop 'ERROR: this routine does not work!'
    do ig=1,ann%ng6,3
        teneria(1:3,1:3)=0.d0
        rcov=ann%g6eta((ig-1)/3+1)
        do jat=1,atoms%natim
            dx=atoms%ratim(1,jat)-atoms%rat(1,iat)
            dy=atoms%ratim(2,jat)-atoms%rat(2,iat)
            dz=atoms%ratim(3,jat)-atoms%rat(3,iat)
            r=sqrt(dx**2+dy**2+dz**2)
!HERE            fc=cutoff_function(r,ann%g6rc((ig-1)/3+1))
            !fc=exp(-r/(4.d0*rcov))
            !fc=amass(jat)*exp(-r/(4.d0*rcov))
            !fc=amass(jat)*exp(-r**2/(6.d0*rcov)**2)
            teneria(1,1)=teneria(1,1)+fc*(dy*dy+dz*dz+2.d0*(rcov*0.5d0)**2)
            teneria(2,2)=teneria(2,2)+fc*(dx*dx+dz*dz+2.d0*(rcov*0.5d0)**2)
            teneria(3,3)=teneria(3,3)+fc*(dx*dx+dy*dy+2.d0*(rcov*0.5d0)**2)
            teneria(1,2)=teneria(1,2)-fc*(dx*dy)
            teneria(1,3)=teneria(1,3)-fc*(dx*dz)
            teneria(2,3)=teneria(2,3)-fc*(dy*dz)
        enddo
        teneria(2,1)=teneria(1,2)
        teneria(3,1)=teneria(1,3)
        teneria(3,2)=teneria(2,3)
        !diagonalize inertia tensor
        call DSYEV('V','L',3,teneria,3,eval,work,lwork,info)
        i0=i0+1
        !ann%yall(i0,iat)=eval(1) !commented due to moving yall
        i0=i0+1
        !ann%yall(i0,iat)=eval(2) !commented due to moving yall
        i0=i0+1
        !ann%yall(i0,iat)=eval(3) !commented due to moving yall
    enddo
end subroutine symmetry_functions_g06
!*****************************************************************************************
function cutoff_function(r, rc) result(fc)
    use mod_interface
    implicit none
    real(8), intent(in):: r, rc
    !local variables
    real(8):: fc, pi
    if(r<rc) then
        fc=(1.d0-(r/rc)**2)**3
    else
        fc=0.d0
    endif
    !pi=4.d0*atan(1.d0)
    !if(r<rc*.5d0) then
    !    fc=1.d0
    !elseif(r<rc) then
    !    fc=cos((1.d0-((1.d0-((r-rc*0.5d0)/(rc*0.5d0))**2)**3))*pi*0.5d0)
    !else
    !    fc=0.d0
    !endif
end function cutoff_function
!*****************************************************************************************
function cutoff_function_der(r, rc) result(fcd)
    use mod_interface
    implicit none
    real(8), intent(in):: r, rc
    !local variables
    real(8):: fcd, pi
    if(r<rc) then
        fcd=-6.d0*(r/rc**2)*(1.d0-(r/rc)**2)**2
    else
        fcd=0.d0
    endif
end function cutoff_function_der
!*****************************************************************************************
