!*****************************************************************************************
subroutine get_hartree(parini,ewald_p3d,atoms,gausswidth,ehartree,g)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_ewald_p3d
    use time_profiling
    use mod_timing , only: TCAT_PSOLVER
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ewald_p3d),intent(inout):: ewald_p3d
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(in):: gausswidth(atoms%nat)
    real(8), intent(out):: ehartree, g(atoms%nat)
    !local variables
    real(8):: dpm, pi, gtot, ecut, epotreal, alphasq
    integer:: iat, igpx, igpy, igpz
    real(8), allocatable:: gwsq(:), ratred(:,:), gg(:) 
    real(8), allocatable::  ewaldwidth(:)
    real(8):: stress(3,3), kmax, c, vol, talpha
    call f_routine(id='get_hartree')
    call f_timing(TCAT_PSOLVER,'ON')
    pi=4.d0*atan(1.d0)
! !   if (parini%ewald .and. parini%alpha_ewald<0.d0) then
!        call getvol_alborz(atoms%cellvec,vol)
!        c=2
!        write(*,*)"alpha optimize", 1.d0/(sqrt(pi)*(atoms%nat/vol**2)**(1.d0/6.d0))
!         
! !   end if

    dpm=0.d0
    do iat=1,atoms%nat
        dpm=dpm+atoms%qat(iat)*atoms%rat(3,iat)
    enddo
    dpm=dpm*2.d0*pi*ewald_p3d%poisson_p3d%ngpx*ewald_p3d%poisson_p3d%ngpy/(ewald_p3d%cell(1)*ewald_p3d%cell(2))
    !do iat=1,atoms%nat
    !    write(33,'(2i4,3es14.5)') iter,iat,atoms%qat(iat),atoms%rat(3,iat),dpm
    !enddo
    epotreal=0.d0
    if(parini%ewald) then
        gg=f_malloc([1.to.atoms%nat],id='gg')
        ewaldwidth=f_malloc([1.to.atoms%nat],id='ewaldwidth')
       ewaldwidth(:)=ewald_p3d%alpha
    end if

    if(trim(parini%psolver_ann)=='kwald') then
        if (parini%ewald) then
            !kmax=2.d0/ewald_p3d%alpha*sqrt(-log(1.d-3))
            kmax=2.d0/ewald_p3d%alpha*sqrt(-log(1.d3*parini%tolerance_ewald))
            ecut=kmax**2/2.d0
            write(*,*)"ecut", ecut, "alpha",ewald_p3d%alpha
        else
            ecut=parini%ecut_ewald
        !!!!!  else
        !!!!!      talpha=minval(gausswidth)
        !!!!!      kmax=2.d0/talpha*sqrt(-log(1.d3*parini%tolerance_ewald))
        !!!!!      !kmax=2.d0/talpha*sqrt(-log(1.d-3))
        !!!!!      ecut=kmax**2/2.d0*(4.d0*pi**2)
        !!!!!      write(*,*)"ecut", ecut, "alpha",ewald_p3d%alpha
        endif
        gwsq=f_malloc([1.to.atoms%nat],id='gwsq')
        ratred=f_malloc([1.to.3,1.to.atoms%nat],id='ratred')
        if(parini%ewald) then
     !       gwsq(1:atoms%nat)=ewaldwidth(1:atoms%nat)**2
     !       call kwald(atoms%nat,atoms%rat,ratred,atoms%qat,atoms%cellvec,gwsq,ecut,ehartree,atoms%fat,g,atoms%stress,atoms%celldv)
            alphasq=ewald_p3d%alpha**2
            call kwald_samare(parini%iverbose,atoms%nat,atoms%rat,ratred,atoms%qat,atoms%cellvec,alphasq,ecut,ehartree,atoms%fat,g,atoms%stress,atoms%celldv)
         else
            gwsq(1:atoms%nat)=gausswidth(1:atoms%nat)**2
            call kwald(parini%iverbose,atoms%nat,atoms%rat,ratred,atoms%qat,atoms%cellvec,gwsq,ecut,ehartree,atoms%fat,g,atoms%stress,atoms%celldv)
        end if
        if(parini%ewald) then
            call real_part(parini,atoms,gausswidth,ewald_p3d%alpha,epotreal,gg,stress)
            atoms%stress=atoms%stress+stress
            ehartree=ehartree+epotreal
            g=g+gg
        end if
        
        call f_free(gwsq)
        call f_free(ratred)
        if(parini%ewald) then
            call f_free(gg)
            call f_free(ewaldwidth)
        endif
        call f_timing(TCAT_PSOLVER,'OF')
        call f_release_routine()
        return
    endif
    if(parini%ewald) then
        call putgaussgrid(parini,atoms%boundcond,.true.,atoms%nat,atoms%rat,atoms%qat,ewaldwidth,ewald_p3d)
    else
        call putgaussgrid(parini,atoms%boundcond,.true.,atoms%nat,atoms%rat,atoms%qat,gausswidth,ewald_p3d)
    end if
    if(trim(atoms%boundcond)=='bulk') then
        if(trim(parini%psolver_ann)=='bigdft') then
            call cal_hartree_pot_bps(ewald_p3d,atoms,ehartree)
            ehartree=ehartree+epotreal
        else
            stop 'ERROR: unknown psolver_ann'
        endif
    elseif(trim(atoms%boundcond)=='slab') then
        call calculate_potener_pot(parini,ewald_p3d%poisson_p3d,ewald_p3d%cell,ewald_p3d%hgx,ewald_p3d%hgy,ewald_p3d%hgz,ehartree,dpm)
        ehartree=ehartree+epotreal
        do igpz=1,ewald_p3d%poisson_p3d%ngpz
        do igpy=1,ewald_p3d%poisson_p3d%ngpy
        do igpx=1,ewald_p3d%poisson_p3d%ngpx
            ewald_p3d%poisson_p3d%rho(igpx,igpy,igpz)=ewald_p3d%poisson_p3d%pot(igpx,igpy,igpz)
        enddo
        enddo
        enddo
    elseif(trim(atoms%boundcond)=='wire') then
        stop 'ERROR: wire BCs is not complete yet.'
    elseif(trim(atoms%boundcond)=='free') then
        stop 'ERROR: free BCs is not complete yet.'
    else
        write(*,'(2a)') 'ERROR: unknown BC in calparam ',trim(atoms%boundcond)
        stop
    endif

    if(trim(parini%psolver_ann)/='kwald') then
        g(1:atoms%nat)=0.d0
        if(parini%ewald) then
            atoms%fat=0.d0
            call get_g_from_pot(parini,atoms,ewald_p3d,ewaldwidth,g)
            call real_part(parini,atoms,gausswidth,ewald_p3d%alpha,epotreal,gg,stress)
            ehartree=ehartree+epotreal
            g=g+gg
        else
            call get_g_from_pot(parini,atoms,ewald_p3d,gausswidth,g)
        end if
    endif
    call apply_external_field(parini,atoms,ewald_p3d,ehartree,g)
    if(parini%ewald) then
        call f_free(gg)
        call f_free(ewaldwidth)
    endif
    call f_timing(TCAT_PSOLVER,'OF')
    call f_release_routine()
end subroutine get_hartree
!*****************************************************************************************
subroutine apply_external_field(parini,atoms,ewald_p3d,ehartree,g)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_ewald_p3d
    !use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ewald_p3d),intent(inout):: ewald_p3d
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(inout):: ehartree, g(atoms%nat)
    !local variables
    !real(8):: dpm, pi, gtot, ecut, epotreal, alphasq
    integer:: iat, igpx, igpy, igpz
    !real(8), allocatable:: gwsq(:), ratred(:,:), gg(:) 
    !real(8), allocatable::  ewaldwidth(:)
    real(8):: ext_pot, dipole
    real(8):: dv, pi, beta

    pi=4.d0*atan(1.d0)
    if((.not. ewald_p3d%poisson_p3d%point_particle) .and. trim(parini%bias_type)=='fixed_efield') then
        do igpz=1,ewald_p3d%poisson_p3d%ngpz
            !igpz=1 is not necessarily z=0, now in this way external potential is not
            !zero at z=0 but since shift the potential, such that external potential
            !to be zero at z=0, has no effect on energy, we keep it this way.
            !ext_pot=parini%efield*igpz*ewald_p3d%hgz
            !do igpy=1,ewald_p3d%poisson_p3d%ngpy
            !    do igpx=1,ewald_p3d%poisson_p3d%ngpx
            !        ewald_p3d%poisson_p3d%pot(igpx,igpy,igpz)=ewald_p3d%poisson_p3d%pot(igpx,igpy,igpz)+ext_pot
            !    enddo
            !enddo
        enddo
        dipole=0.d0
        do iat=1,atoms%nat
            dipole=dipole+atoms%qat(iat)*atoms%rat(3,iat)
        enddo
        ehartree=ehartree+parini%efield*0.5d0*dipole
        do iat=1,atoms%nat
            !atoms%fat(3,iat)=atoms%fat(3,iat)-parini%efield*0.5d0*atoms%qat(iat)
            g(iat)=g(iat)+parini%efield*0.5d0*atoms%rat(3,iat)
        enddo
    elseif((.not. ewald_p3d%poisson_p3d%point_particle) .and. trim(parini%bias_type)=='fixed_potdiff') then
        dipole=0.d0
        do iat=1,atoms%nat
            dipole=dipole+atoms%qat(iat)*atoms%rat(3,iat)
        enddo
        beta=-dipole*2.d0*pi/(ewald_p3d%cell(1)*ewald_p3d%cell(2))
        dv=parini%vu_ewald-parini%vl_ewald
        write(*,*)'real pot = vu , vl ', parini%vu_ewald+dipole, parini%vl_ewald-dipole 

        ehartree=ehartree+(dv+2*beta)*0.5d0*dipole/ewald_p3d%cell(3)
        do iat=1,atoms%nat
            !atoms%fat(3,iat)=atoms%fat(3,iat)-parini%efield*0.5d0*atoms%qat(iat)
            g(iat)=g(iat)+(dv+2*beta)*0.5d0*atoms%rat(3,iat)/ewald_p3d%cell(3)
            g(iat)=g(iat)+(beta)*atoms%rat(3,iat)/ewald_p3d%cell(3)
        enddo
        !efield=0.d0 !to be corrected by Samare
        !do igpz=1,ewald_p3d%poisson_p3d%ngpz
        !    !igpz=1 is not necessarily z=0, now in this way external potential is not
        !    !zero at z=0 but since shift the potential, such that external potential
        !    !to be zero at z=0, has no effect on energy, we keep it this way.
        !    ext_pot=efield*igpz*hz
        !    do igpy=1,ewald_p3d%poisson_p3d%ngpy
        !        do igpx=1,ewald_p3d%poisson_p3d%ngpx
        !            ewald_p3d%poisson_p3d%pot(igpx,igpy,igpz)=ewald_p3d%poisson_p3d%pot(igpx,igpy,igpz)+ext_pot
        !        enddo
        !    enddo
        !enddo
    endif
end subroutine apply_external_field
!*****************************************************************************************
subroutine get_g_from_pot(parini,atoms,ewald_p3d,gausswidth,g)
    use mod_interface
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_ewald_p3d
    use mod_parini, only: typ_parini
    use dynamic_memory
    implicit none
    type(typ_atoms), intent(in):: atoms
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
    type(typ_parini), intent(in):: parini
    !work array that is bigger than rho array, big enough to include of 
    !grid points that are outside of box.
    !local variables
    real(8), allocatable:: wx(:), wy(:), wz(:) !values of one dimensional Gaussian functions
    real(8):: rhoz, rhoyz, pi
    real(8):: g(atoms%nat) 
    real(8):: hgxinv, hgyinv, hgzinv, hgxhgyhgz
    real(8):: width_inv, width_inv_xat, width_inv_yat, width_inv_zat
    real(8):: width_inv_hgx, width_inv_hgy, width_inv_hgz,width
    real(8):: xat, yat, zat, facqiat, fac
    integer:: iat, iw, ix, iy, iz, iatox, iatoy, iatoz, jx, jy, jz, iyt, izt
    real(8):: fx, fy, fz, dx, dy, dz
    real(8), intent(in):: gausswidth(atoms%nat)
    integer:: ngpx, ngpy, ngpz, iii
    real(8), allocatable:: wa(:,:,:)
    associate(nagpx=>ewald_p3d%nagpx,nagpy=>ewald_p3d%nagpy,nagpz=>ewald_p3d%nagpz)
    ngpx=ewald_p3d%poisson_p3d%ngpx
    ngpy=ewald_p3d%poisson_p3d%ngpy
    ngpz=ewald_p3d%poisson_p3d%ngpz
    !wa(1-nagpx:ngpx+nagpx,1-nagpy:ngpy+nagpy,1-nagpz:ngpz+nagpz)=>ewald_p3d%poisson_p3d%pot
    wa=f_malloc([1-nagpx.to.ngpx+nagpx,1-nagpy.to.ngpy+nagpy,1-nagpz.to.ngpz+nagpz],id='wa')
    wx=f_malloc([-ewald_p3d%nbgpx.to.ewald_p3d%nbgpx],id='wx')
    wy=f_malloc([-ewald_p3d%nbgpy.to.ewald_p3d%nbgpy],id='wy')
    wz=f_malloc([-ewald_p3d%nbgpz.to.ewald_p3d%nbgpz],id='wz')
    pi=4.d0*atan(1.d0)
    hgxhgyhgz=ewald_p3d%hgx*ewald_p3d%hgy*ewald_p3d%hgz
    hgxinv=1.d0/ewald_p3d%hgx
    hgyinv=1.d0/ewald_p3d%hgy
    hgzinv=1.d0/ewald_p3d%hgz
    !---------------------------------------------------------------------------
    do iz=1-nagpz,ngpz+nagpz
        izt=iz+(sign(ngpz,-iz)+sign(ngpz,ngpz-iz))/2
        do iy=1-ewald_p3d%nagpy,ngpy+ewald_p3d%nagpy
            iyt=iy+(sign(ngpy,-iy)+sign(ngpy,ngpy-iy))/2
            do ix=1-ewald_p3d%nagpx,0
                wa(ix,iy,iz)=ewald_p3d%poisson_p3d%rho(ix+ngpx,iyt,izt)
            enddo
            do ix=1,ngpx
                wa(ix,iy,iz)=ewald_p3d%poisson_p3d%rho(ix,iyt,izt)
            enddo
            do ix=ngpx+1,ngpx+ewald_p3d%nagpx
                wa(ix,iy,iz)=ewald_p3d%poisson_p3d%rho(ix-ngpx,iyt,izt)
            enddo
        enddo
    enddo
    !---------------------------------------------------------------------------
        if(trim(atoms%boundcond)=='bulk') then
            iii=0
        elseif(trim(atoms%boundcond)=='slab') then
            iii=1
        endif
    !initialize the density 
    if(parini%ewald) then
        width= ewald_p3d%alpha
        width_inv=1.d0/width
        fac=1.d0/(width*sqrt(pi))**3
        width_inv_hgx=width_inv*ewald_p3d%hgx
        width_inv_hgy=width_inv*ewald_p3d%hgy
        width_inv_hgz=width_inv*ewald_p3d%hgz
    endif


    do iat=1,atoms%nat  
        !shift the Gaussian centers
        iatox=nint(atoms%rat(1,iat)*hgxinv)+1
        iatoy=nint(atoms%rat(2,iat)*hgyinv)+1
        iatoz=nint(atoms%rat(3,iat)*hgzinv)+1+ewald_p3d%nbgpz*iii
        xat=atoms%rat(1,iat)-(iatox-1)*ewald_p3d%hgx
        yat=atoms%rat(2,iat)-(iatoy-1)*ewald_p3d%hgy
        zat=atoms%rat(3,iat)-(iatoz-1-ewald_p3d%nbgpz*iii)*ewald_p3d%hgz
        !construct the one-dimensional gaussians

        if(.not.parini%ewald) then
            width=gausswidth(iat)
            width_inv=1.d0/width
            fac=1.d0/(width*sqrt(pi))**3
            width_inv_hgx=width_inv*ewald_p3d%hgx
            width_inv_hgy=width_inv*ewald_p3d%hgy
            width_inv_hgz=width_inv*ewald_p3d%hgz
        endif

        width_inv_xat=width_inv*xat
        width_inv_yat=width_inv*yat
        width_inv_zat=width_inv*zat
        do iw=-ewald_p3d%nbgpx,ewald_p3d%nbgpx
            wx(iw)=exp(-(width_inv_hgx*iw-width_inv_xat)**2)
        enddo
        do iw=-ewald_p3d%nbgpy,ewald_p3d%nbgpy
            wy(iw)=exp(-(width_inv_hgy*iw-width_inv_yat)**2)
        enddo
        do iw=-ewald_p3d%nbgpz,ewald_p3d%nbgpz
            wz(iw)=exp(-(width_inv_hgz*iw-width_inv_zat)**2)
        enddo

        do iz=-ewald_p3d%nbgpz,ewald_p3d%nbgpz
        jz=iatoz+iz
        rhoz=fac*wz(iz)
            do iy=-ewald_p3d%nbgpy,ewald_p3d%nbgpy
                rhoyz=rhoz*wy(iy)
                jy=iatoy+iy
                do ix=ewald_p3d%mboundg(1,iy,iz),ewald_p3d%mboundg(2,iy,iz)
                    jx=iatox+ix
                    g(iat)=g(iat)+rhoyz*wx(ix)*wa(jx,jy,jz)
                enddo
            enddo
        enddo
        g(iat)=g(iat)*hgxhgyhgz
    enddo
    call f_free(wx)
    call f_free(wy)
    call f_free(wz)
    call f_free(wa)
    end associate
    ! call f_release_routine()
end subroutine get_g_from_pot
!*****************************************************************************************
subroutine real_part(parini,atoms,gausswidth,alpha,epotreal,gg,stress)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_linked_lists, only: typ_linked_lists
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_linked_lists):: linked_lists
    real(8):: dx, dy, dz, r, rsq, xiat, yiat, ziat, alphainv, twosqrtinv
    real(8):: t, tt, tt1, tt2, tt3, ttt
    real(8):: rcutsq, fx, fy, fz, pi, hspinv, rhspinv, rinv, qiat, qiatjat, spf, spfd
    integer:: ip, jp, jpt, jl, il, iat_maincell, jat_maincell
    integer:: ix, iy, iz, jy, jz, iat, jat, ipat, isp, maincell, maincell_iat
    real(8), intent(in):: gausswidth(atoms%nat),alpha
    real(8):: gg(atoms%nat),rr
    real(8):: cell(3) , vol, stress(3,3)
    real(8)::epotreal,alphatwoinv,ralphasq,rbetasq,rbetainv,alphasq,betainv
    pi = 4.d0 *atan(1.d0)
    call getvol_alborz(atoms%cellvec,vol)
    rr=linked_lists%rcut
    linked_lists%rcut =sqrt(2.d0)*max(maxval(gausswidth(:)),alpha)*sqrt(-log(parini%tolerance_ewald))
    call linkedlists_init(parini,atoms,cell,linked_lists)
    rcutsq=linked_lists%rcut**2
    alphatwoinv =1.d0/(sqrt(2.d0)*alpha)
    alphasq=(alphatwoinv)**2
    epotreal=0.d0
    gg=0.d0
    tt2=2.d0/sqrt(pi)
    stress = 0.d0
    include 'act1_cell_linkedlist.inc'
    do  iat=ip,il
        qiat=linked_lists%qat(iat)
        xiat=linked_lists%rat(1,iat)
        yiat=linked_lists%rat(2,iat)
        ziat=linked_lists%rat(3,iat)
        jpt=linked_lists%prime(ix+linked_lists%limnbx(1,jy-iy,jz-iz),jy,jz)
        jp=(iat-ip+1)*((isign(1,ip-jpt)+1)/2)+jpt
        jl=linked_lists%last(ix+linked_lists%limnbx(2,jy-iy,jz-iz),jy,jz)
        maincell_iat=linked_lists%maincell(iat)
        do  jat=jp,jl
            dx=xiat-linked_lists%rat(1,jat)
            dy=yiat-linked_lists%rat(2,jat)
            dz=ziat-linked_lists%rat(3,jat)
            rsq= dx*dx+dy*dy+dz*dz
            maincell=maincell_iat+linked_lists%maincell(jat)
            if(rsq<rcutsq .and. maincell >-1) then
                iat_maincell=linked_lists%perm(iat)
                jat_maincell=linked_lists%perm(jat)
                qiatjat=linked_lists%qat(jat)*qiat
                r=sqrt(rsq)
                betainv=1.d0/sqrt(gausswidth(iat_maincell)**2+gausswidth(jat_maincell)**2)
                rbetainv=r*betainv
                tt1= 1.d0/r *(-erfc(rbetainv)+erfc(r*alphatwoinv))
                epotreal = epotreal + tt1*qiatjat
                gg(iat_maincell) = gg(iat_maincell)+tt1*linked_lists%qat(jat)
                gg(jat_maincell) = gg(jat_maincell)+tt1*qiat

                tt3=tt2*(exp(-rsq*alphasq)*alphatwoinv-exp(-rbetainv**2)*betainv)
                ttt=(tt1+tt3)/rsq*qiatjat
                fx=ttt*dx;fy=ttt*dy;fz=ttt*dz
                !-----------------------------------
                atoms%fat(1,iat_maincell)=atoms%fat(1,iat_maincell)+fx
                atoms%fat(2,iat_maincell)=atoms%fat(2,iat_maincell)+fy
                atoms%fat(3,iat_maincell)=atoms%fat(3,iat_maincell)+fz
                atoms%fat(1,jat_maincell)=atoms%fat(1,jat_maincell)-fx
                atoms%fat(2,jat_maincell)=atoms%fat(2,jat_maincell)-fy
                atoms%fat(3,jat_maincell)=atoms%fat(3,jat_maincell)-fz
                stress(1,1)=stress(1,1)+fx*dx
                stress(1,2)=stress(1,2)+fx*dy
                stress(1,3)=stress(1,3)+fx*dz
                stress(2,1)=stress(2,1)+fy*dx
                stress(2,2)=stress(2,2)+fy*dy
                stress(2,3)=stress(2,3)+fy*dz
                stress(3,1)=stress(3,1)+fz*dx
                stress(3,2)=stress(3,2)+fz*dy
                stress(3,3)=stress(3,3)+fz*dz
            endif
        enddo
    enddo
    include 'act2_cell_linkedlist.inc'
    stress=stress/vol
    tt2=0.d0
    t=sqrt(2.d0/pi)
    do iat=1,atoms%nat
        tt2=tt2+(1/gausswidth(iat)-1/alpha)*atoms%qat(iat)**2
        gg(iat)=gg(iat)+(1/gausswidth(iat)-1/alpha)*atoms%qat(iat)*t
    end do
    tt2=tt2/sqrt(2*pi)
    epotreal = epotreal+tt2
    linked_lists%rcut=rr

end subroutine real_part
!*****************************************************************************************
