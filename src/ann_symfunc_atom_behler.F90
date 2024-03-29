!*****************************************************************************************
subroutine symmetry_functions_driver(parini,ann_arr,atoms,mpi_env,symfunc)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_symfunc_data, only: typ_symfunc_data
    use mod_atoms, only: typ_atoms
    use mod_linked_lists, only: typ_pia_arr
    use wrapper_MPI, only: mpi_environment
    use wrapper_MPI, only: fmpi_allreduce, FMPI_SUM
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(in):: atoms
    type(mpi_environment), intent(in):: mpi_env
    type(typ_symfunc_data), intent(inout):: symfunc
    !local variables
    integer:: ig, i
    integer:: iats, iate, mat, mproc
    integer:: iat, ibs, ibe
    type(typ_pia_arr):: pia_arr
    real(8):: cutoff_function, cutoff_function_der
    !real(8):: time0, time1, time2, time3
    external cutoff_function, cutoff_function_der
    integer:: isat, jsat, ksat, ib, ia, ibij, ibik, istat
    !real(8), allocatable:: y(:,:), y0d(:,:,:), y0dr(:,:,:)
    !type(typ_linked_lists):: linked_lists
    call f_routine(id='symmetry_functions_driver')
    associate(rc=>symfunc%linked_lists%rcut)
    symfunc%linked_lists%rcut=ann_arr%rcut
    symfunc%linked_lists%triplex=.true.
    !call cpu_time(time0)
    call call_linkedlist(parini,atoms,.true.,symfunc%linked_lists,pia_arr)
    !call cpu_time(time1)
    !-------------------------------------------------------------------------------------
    associate(ng=>ann_arr%ann(1)%nn(0))
    if(mpi_env%nproc>1) then
        mat=atoms%nat/mpi_env%nproc
        iats=mpi_env%iproc*mat+1
        mproc=mod(atoms%nat,mpi_env%nproc)
        iats=iats+max(0,mpi_env%iproc-mpi_env%nproc+mproc)
        if(mpi_env%iproc>mpi_env%nproc-mproc-1) mat=mat+1
        iate=iats+mat-1
    else
        iats=1
        iate=atoms%nat
        !mat=atoms%nat
    endif
    ibs=symfunc%linked_lists%prime_bound(iats)
    ibe=symfunc%linked_lists%prime_bound(iate+1)-1
    symfunc%y=f_malloc0((/1.to.ng,1.to.atoms%nat/),id='symfunc%y')
    symfunc%y0d=f_malloc0((/1.to.ng,1.to.3,ibs.to.ibe/),id='symfunc%y0d')
    symfunc%y0dr=f_malloc0((/1.to.ng,1.to.9,ibs.to.ibe/),id='symfunc%y0dr')
    !call flm_zero(ng*atoms%nat,symfunc%y)
    !call flm_zero(ng*3*(ibe-ibs+1),symfunc%y0d)
    !call flm_zero(ng*9*(ibe-ibs+1),symfunc%y0dr)
    !symfunc%y=0.d0
    !symfunc%y0d=0.d0
    !symfunc%y0dr=0.d0
    !call cpu_time(time2)
    do ib=1,symfunc%linked_lists%maxbound_rad
        pia_arr%pia(ib)%fc=cutoff_function(pia_arr%pia(ib)%r,rc)
        pia_arr%pia(ib)%fcd=cutoff_function_der(pia_arr%pia(ib)%r,rc)
        iat=symfunc%linked_lists%bound_rad(1,ib)
        if(iat<iats .or. iat>iate) cycle
        isat=atoms%itypat(iat)
        jsat=atoms%itypat(symfunc%linked_lists%bound_rad(2,ib))
        call symmetry_functions_g02_atom(ann_arr,pia_arr%pia(ib),ib,iat,isat,jsat,symfunc)
    enddo
    do ia=1,symfunc%linked_lists%maxbound_ang
        ibij=symfunc%linked_lists%bound_ang(1,ia)
        iat=symfunc%linked_lists%bound_rad(1,ibij)
        if(iat<iats .or. iat>iate) cycle
        ibik=symfunc%linked_lists%bound_ang(2,ia)
        if(symfunc%linked_lists%bound_rad(1,ibij)/=symfunc%linked_lists%bound_rad(1,ibik)) then
            stop 'ERROR: the centered atoms for two bonds is different.'
        endif
        isat=atoms%itypat(iat)
        jsat=atoms%itypat(symfunc%linked_lists%bound_rad(2,ibij))
        ksat=atoms%itypat(symfunc%linked_lists%bound_rad(2,ibik))
        call symmetry_functions_g05_atom(ann_arr,pia_arr%pia(ibij),pia_arr%pia(ibik),ibij,ibik,iat,isat,jsat,ksat,symfunc)
    enddo
    !call cpu_time(time3)
    !write(*,*) 'SF time ',time2-time1
    !write(*,*) 'SS1 time ',time1-time0
    !write(*,*) 'SS2 time ',time2-time1
    !write(*,*) 'SS3 time ',time3-time2
    !write(*,*) 'SSt time ',time3-time0
    !-------------------------------------------------------------------------------------
    !if(parini%iverbose>2) then
    !do iat=1,atoms%nat
    !do ig=1,ann_arr%ann(1)%nn(0)
    !    write(77,'(2i4,es24.15)') ig,iat,symfunc%y(ig,iat)
    !enddo
    !enddo
    !endif
    end associate
    end associate
    deallocate(pia_arr%pia)
    call f_release_routine()
    !stop
end subroutine symmetry_functions_driver
!*****************************************************************************************
subroutine symmetry_functions_g02_atom(ann_arr,pia,ib,iat,isat,jsat,symfunc)
    use mod_ann, only: typ_ann_arr
    use mod_symfunc_data, only: typ_symfunc_data
    use mod_linked_lists, only: typ_pia
    implicit none
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_pia), intent(in):: pia
    integer, intent(in):: ib, iat, isat, jsat
    type(typ_symfunc_data), intent(inout):: symfunc
    !local variables
    integer:: i0, ig
    real(8):: ttei,ttej, ttix,ttiy,ttiz,ttjx,ttjy,ttjz,tt1i,tt1j
    real(8):: etai,etaj, rs, vi,vj,sign_chi0
    real(8):: factor
    i0=ann_arr%ann(isat)%ng1
    sign_chi0=sign(1.d0,ann_arr%ann(jsat)%chi0)
    if (trim(ann_arr%ann(isat)%method)=="angle2") then
        factor= sign_chi0
    else
        factor= 1.d0
    endif
    do ig=1,ann_arr%ann(isat)%ng2
        i0=i0+1
        if((ann_arr%ann(isat)%g2i(ig)/=0).and.(.not.(jsat==ann_arr%ann(isat)%g2i(ig)))) cycle
        rs=ann_arr%ann(jsat)%g2rs(ig)
        !The central atom is i:
        etaj=ann_arr%ann(jsat)%g2eta(ig)

        vi=exp(-etaj*(pia%r-rs)**2) * factor

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
    use mod_ann, only: typ_ann_arr
    use mod_symfunc_data, only: typ_symfunc_data
    implicit none
    type(typ_ann_arr), intent(inout):: ann_arr
    integer, intent(in):: isat, iat, jsat, jat_maincell, ksat, kat_maincell
    real(8), intent(in):: rij, rik, rjk, drij(3), drik(3), drjk(3), fcij, fcdij, fcik, fcdik, fcjk, fcdjk
    type(typ_symfunc_data), intent(inout):: symfunc
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
    use mod_ann, only: typ_ann_arr
    use mod_symfunc_data, only: typ_symfunc_data
    use mod_linked_lists, only: typ_pia
    implicit none
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_pia), intent(in):: piaij, piaik
    integer, intent(in):: ibij, ibik, isat, iat, jsat, ksat
    type(typ_symfunc_data), intent(inout):: symfunc
    !local variables
    integer:: i0, ig, ii1, ii2
    real(8):: ttei, ttek, ttix, ttiy, ttiz, ttej, ttjx, ttjy, ttjz, ttkx, ttky, ttkz, tt1, tt2, tti, ttj, ttk, tt4, tt5
    real(8):: ss1, ss2, ss3, ss4
    real(8):: cos_theta_i, cos_theta_k, cos_theta_j, ui, uk, uj, vi, vk, vj
    real(8):: zeta, alam, etai, etaj, etak, zzz
    real(8):: sign_chi0j, sign_chi0k, sign_chi0i 
    real(8):: factor, fcj_fck, one_rij, one_rik, one_rijk, rijsq, riksq
    real(8):: cos_one_rijsq, cos_one_riksq, fcd_one_rij, fcd_one_rik 
    logical:: cal_stress 
    i0=ann_arr%ann(isat)%ng1+ann_arr%ann(isat)%ng2+ann_arr%ann(isat)%ng3+ann_arr%ann(isat)%ng4
    sign_chi0i=sign(1.d0,ann_arr%ann(isat)%chi0)
    sign_chi0j=sign(1.d0,ann_arr%ann(jsat)%chi0)
    sign_chi0k=sign(1.d0,ann_arr%ann(ksat)%chi0)

    if (trim(ann_arr%ann(isat)%method)=="angle2") then
        factor=(sign_chi0i*(sign_chi0k+sign_chi0j)-1)
    else
        factor=1.d0
    endif
    cos_theta_i=(piaij%dr(1)*piaik%dr(1)+piaij%dr(2)*piaik%dr(2)+piaij%dr(3)*piaik%dr(3))/(piaij%r*piaik%r)
    fcj_fck=piaij%fc*piaik%fc
    one_rij=1.d0/piaij%r
    one_rik=1.d0/piaik%r
    one_rijk=one_rij*one_rik
    cos_one_rijsq=cos_theta_i/piaij%r**2
    cos_one_riksq=cos_theta_i/piaik%r**2
    fcd_one_rij=piaij%fcd/piaij%r
    fcd_one_rik=piaik%fcd/piaik%r
    rijsq=piaij%r**2
    riksq=piaik%r**2
    do ig=1,ann_arr%ann(isat)%ng5
        i0=i0+1
        ii1=ann_arr%ann(isat)%g5i(1,ig)+ann_arr%ann(isat)%g5i(2,ig)
        ii2=abs(ann_arr%ann(isat)%g5i(1,ig)-ann_arr%ann(isat)%g5i(2,ig))
        if((ann_arr%ann(isat)%g5i(1,ig)/=0).and. (.not.((jsat+ksat)==ii1 .and. abs(jsat-ksat)==ii2))) cycle
        zeta=ann_arr%ann(isat)%g5zeta(ig)
        alam=ann_arr%ann(isat)%g5lambda(ig)
        etai=ann_arr%ann(isat)%g5eta(ig)
        etaj=ann_arr%ann(jsat)%g5eta(ig)
        etak=ann_arr%ann(ksat)%g5eta(ig)
        zzz=2.d0**(1.d0-zeta)
        ui=(1.0000000000001d0+alam*cos_theta_i)**zeta
        vi=exp(-(etaj*rijsq+etak*riksq))*factor
        ttei=zzz*ui*vi*fcj_fck
        ss1=vi*fcj_fck
        ss2=ui*fcj_fck
        ss3=ui*vi*piaik%fc
        ss4=ui*vi*piaij%fc
        tt1=zeta*alam*(1.0000000000001d0+alam*cos_theta_i)**(zeta-1)*ss1
        tt2=tt1*one_rijk
        ttj=2.d0*etaj*vi*ss2 !HERE etaj must be checked
        ttk=2.d0*etak*vi*ss2 !HERE etaj must be checked
        tt4=-tt1*cos_one_rijsq-ttj+ss3*fcd_one_rij
        tt5=-tt1*cos_one_riksq-ttk+ss4*fcd_one_rik
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
function cutoff_function(r, rc) result(fc)
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
    !Compare it with: PHYSICAL REVIEW B 93, 155203 (2016) -> comparison is done:
    !the second derivative in the cutoff function in PRB paper does not vanish.
end function cutoff_function
!*****************************************************************************************
function cutoff_function_der(r, rc) result(fcd)
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
