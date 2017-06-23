!*****************************************************************************************
subroutine erfc_surface_zero(parini,atoms,ewald_p3d,nlayer)
    use mod_interface
    use mod_electrostatics, only: typ_ewald_p3d
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_linked_lists
    use mod_parini, only: typ_parini
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
    type(typ_atoms), intent(inout):: atoms
    type(typ_linked_lists):: linked_lists

    !local variables
    real(8):: dx, dy, dz, sclinv, r, rsq, xiat, yiat, ziat, alphainv, twosqrtinv
    real(8):: t, tt, tt1, tt2, tt3, ttt,xgps,ygps,zgps,dzysq,dzsq
    real(8):: fx,fy, fz, pi, hspinv, rhspinv, rinv, qiatp, qiatpjatp, spf, spfd
    real(8):: hxinv, hyinv, hzinv ,zgpu,zgpl,hx,hy,x,y
    real(8):: dzlsq,dzusq,rcutsq
    real(8):: cell(3)
    real(8):: dnlayer, hgxinv , hgyinv, hgzinv 
    real(8):: xat, yat, zat

    integer:: nimat !number of image atoms.
    integer:: nlayer, igpx,igpy,igpz,mx,my,mz, mlimnlayer
    integer:: nkat, njat, jat, kat, jl, npl, npu
    integer:: ix, iy, iz, jx, jy, jz, kx, ky, kz
    integer:: nsclx,nscly,nx,ny,nz , jat1, jat2, iat, isp
    integer:: ip, il, iatox, iatoy, iatoz
    integer:: ngpx, ngpy, ngpz, nbgpx, nbgpy, nbgpz 
    real(8), allocatable:: ratp(:,:,:),qatp(:,:)
    integer, allocatable::  mboundg(:,:,:),mboundgy(:,:)

    alphainv=1.d0/ewald_p3d%alpha;twosqrtinv=1.d0!/sqrt(2.d0)
    ewald_p3d%poisson_p3d%pots=0.d0
    linked_lists%rcut = ewald_p3d%linked_lists%rcut/sqrt(2.d0)
    rcutsq= ewald_p3d%linked_lists%rcut**2

    npl=ewald_p3d%poisson_p3d%npl
    npu=ewald_p3d%poisson_p3d%npu
    zgpl=(npl-1-ewald_p3d%nbgpz)*ewald_p3d%hgz
    zgpu=(npu-1-ewald_p3d%nbgpz)*ewald_p3d%hgz

    ngpx=ewald_p3d%poisson_p3d%ngpx
    ngpy=ewald_p3d%poisson_p3d%ngpy
    ngpz=ewald_p3d%poisson_p3d%ngpz

    nbgpx = int(linked_lists%rcut/ewald_p3d%hgx)+1
    nbgpy = int(linked_lists%rcut/ewald_p3d%hgy)+1
    nbgpz = int(linked_lists%rcut/ewald_p3d%hgz)+1
    
    allocate(mboundg(1:2,-nbgpy:nbgpy,-nbgpz:nbgpz))
    allocate(mboundgy(1:2,-nbgpz:nbgpz))
    call determine_limitsphere(ewald_p3d,mboundg,mboundgy,nbgpx,nbgpy,nbgpz)

    call linkedlists_init(parini,atoms,cell,linked_lists)
    mx = linked_lists%mx
    my = linked_lists%my
    mz = linked_lists%mz
    hzinv=real(mz,8)/cell(3)

    dnlayer=(nlayer-1)*ewald_p3d%hgz
    mlimnlayer=floor(dnlayer*hzinv)+1
    hgxinv=1.d0/ewald_p3d%hgx
    hgyinv=1.d0/ewald_p3d%hgy
    hgzinv=1.d0/ewald_p3d%hgz

    do kz=linked_lists%mz,max(linked_lists%mz-mlimnlayer-linked_lists%mlimnb,1),-1
    do ky=1,linked_lists%my
    do kx=1,linked_lists%mx
    ip=linked_lists%prime(kx,ky,kz)
    il=linked_lists%last(kx,ky,kz)
    do  iat=ip,il
        if (zgpu-linked_lists%rat(3,iat)>linked_lists%rcut+dnlayer) cycle
        iatox=nint(linked_lists%rat(1,iat)*hgxinv)+1
        iatoy=nint(linked_lists%rat(2,iat)*hgyinv)+1
        iatoz=nint(linked_lists%rat(3,iat)*hgzinv)+1+ewald_p3d%nbgpz
        xat=linked_lists%rat(1,iat)-(iatox-1)*ewald_p3d%hgx
        yat=linked_lists%rat(2,iat)-(iatoy-1)*ewald_p3d%hgy
        zat=linked_lists%rat(3,iat)-(iatoz-1-ewald_p3d%nbgpz)*ewald_p3d%hgz
        do iz=-nbgpz,nbgpz
            jz=iatoz+iz
            if (.not. (jz<=npu .and. jz>npu-nlayer)) cycle
            dzsq= (iz*ewald_p3d%hgz-zat)**2
            do iy=mboundgy(1,iz),mboundgy(2,iz)
                dzysq= (iy*ewald_p3d%hgy-yat)**2+dzsq
                jy=modulo(iatoy+iy-1,ngpy)+1
                do ix=mboundg(1,iy,iz),mboundg(2,iy,iz)
                    dx= ix*ewald_p3d%hgx-xat
                    rsq= dx*dx+dzysq
                    if(rsq.lt.rcutsq )then
                        jx=modulo(iatox+ix-1,ngpx)+1
                        r=sqrt(rsq)
                        t=erfc(r*alphainv)
                        tt=linked_lists%qat(iat)*t/r
                        ewald_p3d%poisson_p3d%pots(jx,jy,jz)=ewald_p3d%poisson_p3d%pots(jx,jy,jz)+tt
                    endif
                enddo
            enddo
        enddo
    enddo
    enddo
    enddo
    enddo

    do kz=1,min(mlimnlayer+linked_lists%mlimnb+1,linked_lists%mz)
    do ky=1,linked_lists%my
    do kx=1,linked_lists%mx
    ip=linked_lists%prime(kx,ky,kz)
    il=linked_lists%last(kx,ky,kz)
    do  iat=ip,il
        if (-zgpl+linked_lists%rat(3,iat)>linked_lists%rcut+dnlayer) cycle
        iatox=nint(linked_lists%rat(1,iat)*hgxinv)+1
        iatoy=nint(linked_lists%rat(2,iat)*hgyinv)+1
        iatoz=nint(linked_lists%rat(3,iat)*hgzinv)+1+ewald_p3d%nbgpz
        xat=linked_lists%rat(1,iat)-(iatox-1)*ewald_p3d%hgx
        yat=linked_lists%rat(2,iat)-(iatoy-1)*ewald_p3d%hgy
        zat=linked_lists%rat(3,iat)-(iatoz-1-ewald_p3d%nbgpz)*ewald_p3d%hgz
        do iz=-nbgpz,nbgpz
            jz=iatoz+iz
            if (.not. (jz>=npl .and. jz<npl+nlayer .and. jz<=npu-nlayer)) cycle
            dzsq= (iz*ewald_p3d%hgz-zat)**2
            do iy=mboundgy(1,iz),mboundgy(2,iz)
                jy=modulo(iatoy+iy-1,ngpy)+1
                dzysq= (iy*ewald_p3d%hgy-yat)**2+dzsq
                do ix=mboundg(1,iy,iz),mboundg(2,iy,iz)
                    dx= ix*ewald_p3d%hgx-xat
                    rsq= dx*dx+dzysq
                    if(rsq.lt.rcutsq )then
                        jx=modulo(iatox+ix-1,ngpx)+1
                        r=sqrt(rsq)
                        t=erfc(r*alphainv*twosqrtinv)
                        tt=linked_lists%qat(iat)*t/r
                        ewald_p3d%poisson_p3d%pots(jx,jy,jz)=ewald_p3d%poisson_p3d%pots(jx,jy,jz)+tt
                    endif
                enddo
            enddo
        enddo
    enddo
    enddo
    enddo
    enddo
    call linkedlists_final(linked_lists)
    deallocate(mboundg)
    deallocate(mboundgy)
end subroutine erfc_surface_zero
!*****************************************************************************************
subroutine sollaplaceq(poisson_p3d,hz,cell,vl,vu)
    use mod_interface
    use mod_electrostatics, only: typ_poisson_p3d
    implicit none
    include 'fftw3.f'
    type(typ_poisson_p3d), intent(inout):: poisson_p3d
    !local variables
    real(8):: vl, vu , zlmzu , sinhzlmzu, zlmzuinv
    real(8):: cell(3)
    real(8):: hz , vlmvu, vlzumvuzl 
    real(8):: t,tt,ttt, ghz,ciz,potl,potu
    real(8):: pi, fourpi, fourpisq, gsq, gsqx, gsqy, g, hzsq
    real(8):: fourpisqcellxinvsq, fourpisqcellyinvsq,valuengpxyinv
    real(8):: hzlu,hzliz,hzuiz 
    real(8):: tel,teu,telu 
    integer::ix,iy,iz,ixt,iyt,npu,npl
    integer(8), allocatable:: plan_bs(:),plan_fs(:)
    npl=poisson_p3d%npl
    npu=poisson_p3d%npu
    pi=4.d0*atan(1.d0)
    hzsq=hz**2
    fourpi=4.d0*pi 
    fourpisq=fourpi*pi 
    fourpisqcellxinvsq=fourpisq/cell(1)**2 
    fourpisqcellyinvsq=fourpisq/cell(2)**2 
    allocate(plan_bs(npl:npu))
    allocate(plan_fs(1:2))
    do iz=npl,npu
        call dfftw_plan_dft_c2r_2d(plan_bs(iz),poisson_p3d%ngpx, &
            poisson_p3d%ngpy,poisson_p3d%pots(1,1,iz),poisson_p3d%pots(1,1,iz),fftw_estimate)
    enddo
    call dfftw_plan_dft_r2c_2d(plan_fs(1),poisson_p3d%ngpx, &
            poisson_p3d%ngpy,poisson_p3d%pots(1,1,npl), &
            poisson_p3d%pots(1,1,npl),fftw_estimate)
    call dfftw_plan_dft_r2c_2d(plan_fs(2),poisson_p3d%ngpx, &
            poisson_p3d%ngpy,poisson_p3d%pots(1,1,npu), &
            poisson_p3d%pots(1,1,npu),fftw_estimate)

    do iy=1,poisson_p3d%ngpy
    do ix=1,poisson_p3d%ngpx
        poisson_p3d%pots(ix,iy,npl)=-(poisson_p3d%pots(ix,iy,npl)+poisson_p3d%rho(ix,iy,npl))+vl
        poisson_p3d%pots(ix,iy,npu)=-(poisson_p3d%pots(ix,iy,npu)+poisson_p3d%rho(ix,iy,npu))+vu    
    enddo 
    enddo 
    call dfftw_execute(plan_fs(1))
    call dfftw_execute(plan_fs(2))   

   ! k=0 , l=0
     potl=poisson_p3d%pots(1,1,npl)
     potu=poisson_p3d%pots(1,1,npu)

    !-----------------------------------------
    ix=1;iy=1
     do iz=npl+1,npu-1
         ciz=((potl-potu)*iz-(potl*npu-potu*npl))/real((npl-npu),8)
         poisson_p3d%pots(ix,iy,iz)=ciz
     enddo 
    !-----------------------------------------
     hzlu=-hz*(npl-npu)
     if (hzlu<=30.d0) then
         do iz=npl+1,npu-1
             hzliz=-hz*(npl-iz)
             hzuiz=-hz*(iz-npu)
             !-----------------------------------------
                 ix=poisson_p3d%ngpx+1;iy=1
                 gsq=fourpisqcellxinvsq*((ix-1)/2)**2 + fourpisqcellyinvsq*(iy-1)**2
                 g=sqrt(gsq)
                 t= poisson_p3d%pots(ix,iy,npu)*sinh(g*hzliz)
                 tt=poisson_p3d%pots(ix,iy,npl)*sinh(g*hzuiz)
                 ciz=(t+tt)/sinh(g*hzlu)
                 poisson_p3d%pots(ix,iy,iz)=ciz
             !-----------------------------------------
                 ix=1;iy=poisson_p3d%ngpy/2+1
                 gsq=fourpisqcellyinvsq*(iy-1)**2
                 g=sqrt(gsq)
                 t= poisson_p3d%pots(ix,iy,npu)*sinh(g*hzliz)
                 tt=poisson_p3d%pots(ix,iy,npl)*sinh(g*hzuiz)
                 ciz=(t+tt)/sinh(g*hzlu)
                 poisson_p3d%pots(ix,iy,iz)=ciz
             !-----------------------------------------
                 ix=poisson_p3d%ngpx+1;iy=poisson_p3d%ngpy/2+1
                 gsq=fourpisqcellxinvsq*((ix-1)/2)**2 + fourpisqcellyinvsq*(iy-1)**2
                 g=sqrt(gsq)
                 t= poisson_p3d%pots(ix,iy,npu)*sinh(g*hzliz)
                 tt=poisson_p3d%pots(ix,iy,npl)*sinh(g*hzuiz)
                 ciz=(t+tt)/sinh(g*hzlu)
                 poisson_p3d%pots(ix,iy,iz)=ciz
             !-----------------------------------------
                 ix=1
                 do iy=2,poisson_p3d%ngpy/2
                     gsq=fourpisqcellyinvsq*(iy-1)**2
                     g=sqrt(gsq)
                     t= poisson_p3d%pots(ix,iy,npu)*sinh(g*hzliz)
                     tt=poisson_p3d%pots(ix,iy,npl)*sinh(g*hzuiz)
                     ciz=(t+tt)/sinh(g*hzlu)
                     poisson_p3d%pots(ix,iy,iz)=ciz
                     poisson_p3d%pots(ix,poisson_p3d%ngpy-iy+2,iz)=ciz

                     ixt=ix+1
                     t= poisson_p3d%pots(ixt,iy,npu)*sinh(g*hzliz)
                     tt=poisson_p3d%pots(ixt,iy,npl)*sinh(g*hzuiz)
                     ciz=(t+tt)/sinh(g*hzlu)
                     poisson_p3d%pots(ixt,iy,iz)=ciz
                     poisson_p3d%pots(ixt,poisson_p3d%ngpy-iy+2,iz)=-ciz
                 enddo
             !-----------------------------------------
                 ix=poisson_p3d%ngpx+1
                 gsqx=fourpisqcellxinvsq*(poisson_p3d%ngpx/2)**2
                 do iy=2,poisson_p3d%ngpy/2
                     gsq=gsqx+fourpisqcellyinvsq*(iy-1)**2
                     g=sqrt(gsq)
                     t= poisson_p3d%pots(ix,iy,npu)*sinh(g*hzliz)
                     tt=poisson_p3d%pots(ix,iy,npl)*sinh(g*hzuiz)
                     ciz=(t+tt)/sinh(g*hzlu)
                     poisson_p3d%pots(ix,iy,iz)=ciz
                     poisson_p3d%pots(ix,poisson_p3d%ngpy-iy+2,iz)=ciz

                     ixt=ix+1
                     t= poisson_p3d%pots(ixt,iy,npu)*sinh(g*hzliz)
                     tt=poisson_p3d%pots(ixt,iy,npl)*sinh(g*hzuiz)
                     ciz=(t+tt)/sinh(g*hzlu)
                     poisson_p3d%pots(ixt,iy,iz)=ciz
                     poisson_p3d%pots(ixt,poisson_p3d%ngpy-iy+2,iz)=-ciz
                 enddo
             !-----------------------------------------
                 iy=1
                 gsqy=fourpisqcellyinvsq*(iy-1)**2
                 do ix=3,poisson_p3d%ngpx,2
                     gsq=gsqy+fourpisqcellxinvsq*((ix-1)/2)**2
                     g=sqrt(gsq)
                     t= poisson_p3d%pots(ix,iy,npu)*sinh(g*hzliz)
                     tt=poisson_p3d%pots(ix,iy,npl)*sinh(g*hzuiz)
                     ciz=(t+tt)/sinh(g*hzlu)
                     poisson_p3d%pots(ix,iy,iz)=ciz

                     ixt=ix+1
                     t= poisson_p3d%pots(ixt,iy,npu)*sinh(g*hzliz)
                     tt=poisson_p3d%pots(ixt,iy,npl)*sinh(g*hzuiz)
                     ciz=(t+tt)/sinh(g*hzlu)
                     poisson_p3d%pots(ixt,iy,iz)=ciz
                enddo
            !-----------------------------------------
                iy=poisson_p3d%ngpy/2+1
                gsqy=fourpisqcellyinvsq*(iy-1)**2
                do ix=3,poisson_p3d%ngpx,2
                    gsq=gsqy+fourpisqcellxinvsq*((ix-1)/2)**2
                    g=sqrt(gsq)
                    t= poisson_p3d%pots(ix,iy,npu)*sinh(g*hzliz)
                    tt=poisson_p3d%pots(ix,iy,npl)*sinh(g*hzuiz)
                    ciz=(t+tt)/sinh(g*hzlu)
                    poisson_p3d%pots(ix,iy,iz)=ciz

                    ixt=ix+1
                    t= poisson_p3d%pots(ixt,iy,npu)*sinh(g*hzliz)
                    tt=poisson_p3d%pots(ixt,iy,npl)*sinh(g*hzuiz)
                    ciz=(t+tt)/sinh(g*hzlu)
                    poisson_p3d%pots(ixt,iy,iz)=ciz
                enddo
            !-----------------------------------------
                do iy=2,poisson_p3d%ngpy/2
                    iyt=poisson_p3d%ngpy-iy+2
                    gsqy=fourpisqcellyinvsq*(iy-1)**2
                    do ix=3,poisson_p3d%ngpx,2
                        gsq=gsqy+fourpisqcellxinvsq*((ix-1)/2)**2
                        g=sqrt(gsq)
                        t= poisson_p3d%pots(ix,iy,npu)*sinh(g*hzliz)
                        tt=poisson_p3d%pots(ix,iy,npl)*sinh(g*hzuiz)
                        ciz=(t+tt)/sinh(g*hzlu)
                        poisson_p3d%pots(ix,iy,iz)=ciz

                        ixt=ix+1
                        t= poisson_p3d%pots(ixt,iy,npu)*sinh(g*hzliz)
                        tt=poisson_p3d%pots(ixt,iy,npl)*sinh(g*hzuiz)
                        ciz=(t+tt)/sinh(g*hzlu)
                        poisson_p3d%pots(ixt,iy,iz)=ciz

                        t= poisson_p3d%pots(ix,iyt,npu)*sinh(g*hzliz)
                        tt=poisson_p3d%pots(ix,iyt,npl)*sinh(g*hzuiz)
                        ciz=(t+tt)/sinh(g*hzlu)
                        poisson_p3d%pots(ix,iyt,iz)=ciz

                        ixt=ix+1
                        t= poisson_p3d%pots(ixt,iyt,npu)*sinh(g*hzliz)
                        tt=poisson_p3d%pots(ixt,iyt,npl)*sinh(g*hzuiz)
                        ciz=(t+tt)/sinh(g*hzlu)
                        poisson_p3d%pots(ixt,iyt,iz)=ciz
                    enddo
                enddo
        enddo
            !-----------------------------------------
            write(*,*)"part sinh finished"

    else
         do iz=npl+1,npu-1
             hzliz=-hz*(npl-iz)
             hzuiz=-hz*(iz-npu)
             !-----------------------------------------
                 ix=poisson_p3d%ngpx+1;iy=1
                 gsq=fourpisqcellxinvsq*((ix-1)/2)**2 + fourpisqcellyinvsq*(iy-1)**2
                 g=sqrt(gsq)
                 teu =exp(-g*hzuiz)
                 tel =exp(-g*hzliz)
                 telu=teu*tel!exp(-g*hzlu)

                 t= poisson_p3d%pots(ix,iy,npu)*(teu-tel*telu)
                 tt=poisson_p3d%pots(ix,iy,npl)*(tel-teu*telu)
                 ciz=(t+tt)/(1.d0-telu**2)
                 poisson_p3d%pots(ix,iy,iz)=ciz
             !-----------------------------------------
                 ix=1;iy=poisson_p3d%ngpy/2+1
                 gsq=fourpisqcellyinvsq*(iy-1)**2
                 g=sqrt(gsq)
                 teu =exp(-g*hzuiz)
                 tel =exp(-g*hzliz)
                 telu=teu*tel!exp(-g*hzlu)

                 t= poisson_p3d%pots(ix,iy,npu)*(teu-tel*telu)
                 tt=poisson_p3d%pots(ix,iy,npl)*(tel-teu*telu)
                 ciz=(t+tt)/(1.d0-telu**2)
                 poisson_p3d%pots(ix,iy,iz)=ciz
             !-----------------------------------------
                 ix=poisson_p3d%ngpx+1;iy=poisson_p3d%ngpy/2+1
                 gsq=fourpisqcellxinvsq*((ix-1)/2)**2 + fourpisqcellyinvsq*(iy-1)**2
                 g=sqrt(gsq)
                 teu =exp(-g*hzuiz)
                 tel =exp(-g*hzliz)
                 telu=teu*tel!exp(-g*hzlu)

                 t= poisson_p3d%pots(ix,iy,npu)*(teu-tel*telu)
                 tt=poisson_p3d%pots(ix,iy,npl)*(tel-teu*telu)
                 ciz=(t+tt)/(1.d0-telu**2)
                 poisson_p3d%pots(ix,iy,iz)=ciz
             !-----------------------------------------
                 ix=1
                 do iy=2,poisson_p3d%ngpy/2
                     gsq=fourpisqcellyinvsq*(iy-1)**2
                     g=sqrt(gsq)
                     teu =exp(-g*hzuiz)
                     tel =exp(-g*hzliz)
                     telu=teu*tel!exp(-g*hzlu)

                     t= poisson_p3d%pots(ix,iy,npu)*(teu-tel*telu)
                     tt=poisson_p3d%pots(ix,iy,npl)*(tel-teu*telu)
                     ciz=(t+tt)/(1.d0-telu**2)
                     poisson_p3d%pots(ix,iy,iz)=ciz
                     poisson_p3d%pots(ix,poisson_p3d%ngpy-iy+2,iz)=ciz

                     ixt=ix+1
                     t= poisson_p3d%pots(ixt,iy,npu)*(teu-tel*telu)
                     tt=poisson_p3d%pots(ixt,iy,npl)*(tel-teu*telu)
                     ciz=(t+tt)/(1.d0-telu**2)
                     poisson_p3d%pots(ixt,iy,iz)=ciz
                     poisson_p3d%pots(ixt,poisson_p3d%ngpy-iy+2,iz)=-ciz
                 enddo
             !-----------------------------------------
                 ix=poisson_p3d%ngpx+1
                 gsqx=fourpisqcellxinvsq*(poisson_p3d%ngpx/2)**2
                 do iy=2,poisson_p3d%ngpy/2
                     gsq=gsqx+fourpisqcellyinvsq*(iy-1)**2
                     g=sqrt(gsq)
                     teu =exp(-g*hzuiz)
                     tel =exp(-g*hzliz)
                     telu=teu*tel!exp(-g*hzlu)

                     t= poisson_p3d%pots(ix,iy,npu)*(teu-tel*telu)
                     tt=poisson_p3d%pots(ix,iy,npl)*(tel-teu*telu)
                     ciz=(t+tt)/(1.d0-telu**2)
                     poisson_p3d%pots(ix,iy,iz)=ciz
                     poisson_p3d%pots(ix,poisson_p3d%ngpy-iy+2,iz)=ciz

                     ixt=ix+1
                     t= poisson_p3d%pots(ixt,iy,npu)*(teu-tel*telu)
                     tt=poisson_p3d%pots(ixt,iy,npl)*(tel-teu*telu)
                     ciz=(t+tt)/(1.d0-telu**2)
                     poisson_p3d%pots(ixt,iy,iz)=ciz
                     poisson_p3d%pots(ixt,poisson_p3d%ngpy-iy+2,iz)=-ciz
                 enddo
             !-----------------------------------------
                 iy=1
                 gsqy=fourpisqcellyinvsq*(iy-1)**2
                 do ix=3,poisson_p3d%ngpx,2
                     gsq=gsqy+fourpisqcellxinvsq*((ix-1)/2)**2
                     g=sqrt(gsq)
                     teu =exp(-g*hzuiz)
                     tel =exp(-g*hzliz)
                     telu=teu*tel!exp(-g*hzlu)

                     t= poisson_p3d%pots(ix,iy,npu)*(teu-tel*telu)
                     tt=poisson_p3d%pots(ix,iy,npl)*(tel-teu*telu)
                     ciz=(t+tt)/(1.d0-telu**2)
                     poisson_p3d%pots(ix,iy,iz)=ciz

                     ixt=ix+1
                     t= poisson_p3d%pots(ixt,iy,npu)*(teu-tel*telu)
                     tt=poisson_p3d%pots(ixt,iy,npl)*(tel-teu*telu)
                     ciz=(t+tt)/(1.d0-telu**2)
                     poisson_p3d%pots(ixt,iy,iz)=ciz
                enddo
            !-----------------------------------------
                iy=poisson_p3d%ngpy/2+1
                gsqy=fourpisqcellyinvsq*(iy-1)**2
                do ix=3,poisson_p3d%ngpx,2
                    gsq=gsqy+fourpisqcellxinvsq*((ix-1)/2)**2
                    g=sqrt(gsq)
                    teu =exp(-g*hzuiz)
                    tel =exp(-g*hzliz)
                    telu=teu*tel!exp(-g*hzlu)

                    t= poisson_p3d%pots(ix,iy,npu)*(teu-tel*telu)
                    tt=poisson_p3d%pots(ix,iy,npl)*(tel-teu*telu)
                    ciz=(t+tt)/(1.d0-telu**2)
                    poisson_p3d%pots(ix,iy,iz)=ciz

                    ixt=ix+1
                    t= poisson_p3d%pots(ixt,iy,npu)*(teu-tel*telu)
                    tt=poisson_p3d%pots(ixt,iy,npl)*(tel-teu*telu)
                    ciz=(t+tt)/(1.d0-telu**2)
                    poisson_p3d%pots(ixt,iy,iz)=ciz
                enddo
            !-----------------------------------------
                do iy=2,poisson_p3d%ngpy/2
                    iyt=poisson_p3d%ngpy-iy+2
                    gsqy=fourpisqcellyinvsq*(iy-1)**2
                    do ix=3,poisson_p3d%ngpx,2
                        gsq=gsqy+fourpisqcellxinvsq*((ix-1)/2)**2
                        g=sqrt(gsq)
                        teu =exp(-g*hzuiz)
                        tel =exp(-g*hzliz)
                        telu=teu*tel!exp(-g*hzlu)

                        t= poisson_p3d%pots(ix,iy,npu)*(teu-tel*telu)
                        tt=poisson_p3d%pots(ix,iy,npl)*(tel-teu*telu)
                        ciz=(t+tt)/(1.d0-telu**2)
                        poisson_p3d%pots(ix,iy,iz)=ciz

                        ixt=ix+1
                        t= poisson_p3d%pots(ixt,iy,npu)*(teu-tel*telu)
                        tt=poisson_p3d%pots(ixt,iy,npl)*(tel-teu*telu)
                        ciz=(t+tt)/(1.d0-telu**2)
                        poisson_p3d%pots(ixt,iy,iz)=ciz

                        t= poisson_p3d%pots(ix,iyt,npu)*(teu-tel*telu)
                        tt=poisson_p3d%pots(ix,iyt,npl)*(tel-teu*telu)
                        ciz=(t+tt)/(1.d0-telu**2)
                        poisson_p3d%pots(ix,iyt,iz)=ciz

                        ixt=ix+1
                        t= poisson_p3d%pots(ixt,iyt,npu)*(teu-tel*telu)
                        tt=poisson_p3d%pots(ixt,iyt,npl)*(tel-teu*telu)
                        ciz=(t+tt)/(1.d0-telu**2)
                        poisson_p3d%pots(ixt,iyt,iz)=ciz
                    enddo
                enddo
        enddo
            !-----------------------------------------
            write(*,*)"part exp finished"
    endif
  
    do iz=npl,npu
        call dfftw_execute(plan_bs(iz))
    enddo

    valuengpxyinv=1.d0/real(poisson_p3d%ngpx*poisson_p3d%ngpy,8)
    do iz=npl,npu
    do iy=1,poisson_p3d%ngpy
    do ix=1,poisson_p3d%ngpx
        poisson_p3d%pots(ix,iy,iz)=poisson_p3d%pots(ix,iy,iz)*valuengpxyinv
    enddo
    enddo
    enddo
 
    do iz=npl,npu
        call dfftw_destroy_plan(plan_bs(iz))
    enddo
    call dfftw_destroy_plan(plan_fs)
    deallocate(plan_bs)
    deallocate(plan_fs)
end subroutine sollaplaceq
!*****************************************************************************************
 subroutine calculate_force_ener_plane(atoms,ewald_p3d,epot)
    use mod_interface
    use mod_electrostatics, only: typ_ewald_p3d
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
    type(typ_atoms), intent(inout):: atoms

    !local variables
    integer:: ix,iy,iz ,jx,jy,jz, kx,ky,kz, iat
    integer:: nlgx,nlgy,nlgz 
    integer:: ngpx,ngpy,ngpz 
    integer:: ngpx2,ngpy2,ngpz2 
    integer:: npl,npu 
    real(8), allocatable:: wx(:), wy(:), wz(:) 
    real(8), allocatable:: LGx(:), LGy(:), LGz(:) 
    real(8), allocatable:: DLGx(:), DLGy(:), DLGz(:) 
    real(8), allocatable::pots_local(:,:,:),potg(:,:,:),pot_atom(:)
    integer:: ix1,iy1,iz1,ii 
    real(8):: time1,time2
    real(8):: x,y,z ,t,tl ,epot ,t1,t2
    real(8):: tfx,tfy,tfz ,tmp
    real(8):: fatp(3,atoms%nat) 
    nlgx=9;  nlgy=9; nlgz=8
    ngpx=ewald_p3d%poisson_p3d%ngpx
    ngpy=ewald_p3d%poisson_p3d%ngpy
    ngpz=ewald_p3d%poisson_p3d%ngpz
    npl=ewald_p3d%poisson_p3d%npl
    npu=ewald_p3d%poisson_p3d%npu
    do iat=1,atoms%nat
        if (atoms%rat(3,iat) < (int((nlgz/2.d0))*ewald_p3d%hgz)) then
             write(*,*) 'ERROR:atoms are too close to lower plane ',atoms%rat(3,iat),npl,ceiling(nlgz/2.d0)*ewald_p3d%hgz 
             stop
        else if (ewald_p3d%cell(3)-atoms%rat(3,iat) < (int((nlgz/2.d0))*ewald_p3d%hgz)) then
             write(*,*) 'ERROR:atoms are too close to upper plane' 
             stop
        endif
    enddo
    allocate(pots_local(-nlgx+1:ewald_p3d%poisson_p3d%ngpx+nlgx,-nlgy+1:ewald_p3d%poisson_p3d%ngpy+nlgy,npl:npu))
    allocate(wx(1:nlgx))
    allocate(wy(1:nlgy))
    allocate(wz(1:nlgz))
    allocate(LGx(1:nlgx))
    allocate(LGy(1:nlgy))
    allocate(LGz(1:nlgz))
    allocate(DLGx(1:nlgx))
    allocate(DLGy(1:nlgy))
    allocate(DLGz(1:nlgz))
    allocate(pot_atom(1:atoms%nat))
    do iz=npl,npu
    do ix=1-nlgx,ewald_p3d%poisson_p3d%ngpx+nlgx
    do iy=-nlgy+1,ewald_p3d%poisson_p3d%ngpy+nlgy
        pots_local(ix,iy,iz)=ewald_p3d%poisson_p3d%pots(modulo(ix-1,ngpx)+1,modulo(iy-1,ngpy)+1,iz)
    enddo
    enddo
    enddo

    call LG_weight(nlgx,nlgy,nlgz,ewald_p3d%hgx,ewald_p3d%hgy,ewald_p3d%hgz,wx,wy,wz)
    pot_atom=0.d0
    t1=0.d0
    t2=0.d0
    fatp=0.d0
    call cpu_time (time1)
    do iat=1,atoms%nat
        x=atoms%rat(1,iat)
        y=atoms%rat(2,iat)
        z=atoms%rat(3,iat)
        !call LGW(nlgx, wx,ewald_p3d%hgx, x, LGx, DLGx, ix1, 0)
        !call LGW(nlgy, wy,ewald_p3d%hgy, y, LGy, DLGy, iy1, 0)
        !call LGW(nlgz, wz,ewald_p3d%hgz, z, LGz, DLGz, iz1, ewald_p3d%nbgpz )
        
        call LGW4(nlgx, wx,ewald_p3d%hgx, x, LGx, DLGx, ix1, 0)
        call LGW4(nlgy, wy,ewald_p3d%hgy, y, LGy, DLGy, iy1, 0)
        call LGW4(nlgz, wz,ewald_p3d%hgz, z, LGz, DLGz, iz1, ewald_p3d%nbgpz )
        do iz=1,nlgz 
            jz=iz1+iz-1
            do iy=1,nlgy 
                jy=iy1+iy-1
                tl=LGy(iy)*LGz(iz)
                tfy=DLGy(iy)*LGz(iz)
                tfz=LGy(iy)*DLGz(iz)
                t=0.d0 ; tfx=0.d0
                do ix=1,nlgx
                    tmp=pots_local(ix1+ix-1,jy,jz)
                    t=t+tmp*LGx(ix)
                    tfx=tfx+tmp*DLGx(ix)
                enddo 
                pot_atom(iat)=pot_atom(iat)+t*tl
                fatp(1,iat)=fatp(1,iat)+tfx*tl
                fatp(2,iat)=fatp(2,iat)+tfy*t
                fatp(3,iat)=fatp(3,iat)+tfz*t
            enddo 
        enddo 
    enddo 
    call cpu_time (time2)
    write(*,*)"------------------------------------------"
    write(*,*)"time for interpolation =",time2-time1
    epot=0.d0
    do iat=1,atoms%nat
        epot=epot+pot_atom(iat)*atoms%qat(iat)
        atoms%fat(1,iat)=atoms%fat(1,iat)-fatp(1,iat)*atoms%qat(iat)
        atoms%fat(2,iat)=atoms%fat(2,iat)-fatp(2,iat)*atoms%qat(iat)
        atoms%fat(3,iat)=atoms%fat(3,iat)-fatp(3,iat)*atoms%qat(iat)
    enddo
    epot=0.5*epot    

    write(*,*)"epotplane=" ,epot
 !   write(*,*)"sum force",sum(fatp(1,:)),sum(fatp(2,:)),sum(fatp(3,:))
   ! write(*,*)
   ! write(*,*)"------------------------------------------"
    deallocate(pots_local)
    deallocate(wx)
    deallocate(wy)
    deallocate(wz)
    deallocate(LGx)
    deallocate(LGy)
    deallocate(LGz)
    deallocate(DLGx)
    deallocate(DLGy)
    deallocate(DLGz)
    deallocate(pot_atom)
end subroutine calculate_force_ener_plane
!*****************************************************************************************
subroutine LG_weight(nlx,nly,nlz,hx,hy,hz,wx,wy,wz)
    implicit none
    integer:: nlx ,nly, nlz !number of point for Lagrange interpolation
    integer:: i  ,maxnl
    real(8):: hx ,hy ,hz , hxp , hyp, hzp
    real(8), allocatable::factorial(:)
    real(8):: wx(nlx), wy(nly), wz(nlz) 
    maxnl=max(nlx,nly,nlz)
    allocate(factorial(0:maxnl-1))
    factorial(0)=1
    do i=1,maxnl-1
        factorial(i)=factorial(i-1)*i
    enddo
    do i=1,nlx
        wx(i)=(-1)**i*factorial(nlx-i)*factorial(i-1)
    enddo
    do i=1,nly
        wy(i)=(-1)**i*factorial(nly-i)*factorial(i-1)
    enddo
    do i=1,nlz
        wz(i)=(-1)**i*factorial(nlz-i)*factorial(i-1)
    enddo
    if (mod(nlx,2)/=0) wx=-wx
    if (mod(nly,2)/=0) wy=-wy
    if (mod(nlz,2)/=0) wz=-wz
    hxp=hx**(nlx-1)          ! how to optimize it ?
    wx(1:nlx) =1.d0/(wx(1:nlx)*hxp)
    hyp=hy**(nly-1)
    wy(1:nly) =1.d0/(wy(1:nly)*hyp)
    hzp=hz**(nlz-1)
    wz(1:nlz) =1.d0/(wz(1:nlz)*hzp)
    deallocate(factorial)
end subroutine LG_weight
!*********************************************************************************************
subroutine LGW(n, w, h, x, LGx, DLGx, ix1, nbgp)
implicit none
    integer:: ix, jx,kx, ix1,ixo ,n ,n2 ,nbgp
    real(8):: w(n), q(n),LGx(n),DLGx(n),h ,x ,diffx,x1,protot,t
    n2=(n+1)/2
    !ixo=int(real(x/h),4)+1
    ixo=nint(x/h)+1+nbgp
    ix1=ixo-(n2-1)
    x1=(ix1-1-nbgp)*h
    diffx=x-x1
    do ix=0,n-1
        q(ix+1)=diffx-(ix)*h
    enddo
    LGx=1.d0 ; DLGx=0.d0
    do ix=1,n
        do jx=1,n
            if (ix==jx) cycle
            LGx(ix)=LGx(ix)*q(jx)
            t=1.d0
            do kx=1,n
                if (kx==jx .or. kx==ix) cycle
                    t=t*q(kx)
            enddo
            DLGx(ix)=DLGx(ix)+t
        enddo
        LGx(ix)=LGx(ix)*w(ix)
        DLGx(ix)=DLGx(ix)*w(ix)
    enddo
end subroutine LGW
!*****************************************************************************
subroutine LGW4(n, w, h, x, LGx, DLGx, ix1, nbgp)
implicit none
    
    integer:: ix, jx, ix1,ixo ,n ,n2 ,nbgp
    real(8):: w(n), q(n), qinv(n), LGx(n), DLGx(n), h ,x ,diffx,x1,protot
    real(8):: t,tt 
    n2=(n+1)/2
    !ixo=int(real((x/h)),4)+1+nbgp
    ixo=nint(x/h)+1+nbgp
    ix1=ixo-(n2-1)
    x1=(ix1-1-nbgp)*h
    diffx=x-x1
    protot=1.d0
    DLGx=0.d0
    if( (ixo-1-nbgp)*h/=x)then
        do ix=0,n-1
            q(ix+1)=diffx-ix*h
            protot=protot*q(ix+1)
            qinv(ix+1)=1.d0/q(ix+1)
        enddo

        do ix=1,n
            LGx(ix)=protot*qinv(ix)
            LGx(ix)=LGx(ix)*w(ix)
            t=0.d0            
            do jx=1,n
               if (ix==jx) cycle
               t=t+qinv(jx)
            enddo
            DLGx(ix)=LGx(ix)*t
        enddo
    else
        LGx=0.d0
        LGx(n2)=1.d0
        do ix=0,n-1
            if (ix==n2-1) cycle
            q(ix+1)=diffx-ix*h
            protot=protot*q(ix+1)
        enddo
        do ix=1,n
            if (ix==n2) cycle
            DLGx(ix)=protot/q(ix)
        enddo
        tt=sum(DLGx)
        DLGx(n2)=tt
        DLGx(1:n)=DLGx(1:n)*w(1:n)
    endif
end subroutine LGW4
!*******************************************************************************************
subroutine surface_charge(parini,ewald_p3d,pot_short,vl,vu)
    use mod_electrostatics, only: typ_ewald_p3d
    use mod_parini, only: typ_parini
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
    integer::ix, iy, iz, npl,npu
    real(8):: t, tt ,density(ewald_p3d%poisson_p3d%ngpx,ewald_p3d%poisson_p3d%ngpy,2),vl,vu
    real(8)::hgzinv,pi,pot_layerl,pot_layeru,pot_short(ewald_p3d%poisson_p3d%ngpx,ewald_p3d%poisson_p3d%ngpy,2,5)
    real(8)::pot_layerl2,pot_layeru2
    real(8)::pot_layerl3,pot_layeru3
    real(8)::pot_layerl4,pot_layeru4
    real(8)::pot_layerl0,pot_layeru0
    real(8):: E,d
    pi=4*atan(1.d0)
    npl=ewald_p3d%poisson_p3d%npl
    npu=ewald_p3d%poisson_p3d%npu
    hgzinv=1.d0/(ewald_p3d%hgz*4.d0*pi)
    t=0.d0
    tt=0.d0
    d = ewald_p3d%cell(3)
            E =- (vu-vl)/d

    do iy=1,ewald_p3d%poisson_p3d%ngpy
    do ix=1,ewald_p3d%poisson_p3d%ngpx
            pot_layerl4=ewald_p3d%poisson_p3d%pots(ix,iy,npl+4)+ewald_p3d%poisson_p3d%rho(ix,iy,npl+4)+pot_short(ix,iy,1,5)
            pot_layerl3=ewald_p3d%poisson_p3d%pots(ix,iy,npl+3)+ewald_p3d%poisson_p3d%rho(ix,iy,npl+3)+pot_short(ix,iy,1,4)
            pot_layerl2=ewald_p3d%poisson_p3d%pots(ix,iy,npl+2)+ewald_p3d%poisson_p3d%rho(ix,iy,npl+2)+pot_short(ix,iy,1,3)
            pot_layerl =ewald_p3d%poisson_p3d%pots(ix,iy,npl+1)+ewald_p3d%poisson_p3d%rho(ix,iy,npl+1)+pot_short(ix,iy,1,2)
              vl       =ewald_p3d%poisson_p3d%pots(ix,iy,npl  )+ewald_p3d%poisson_p3d%rho(ix,iy,npl  )+pot_short(ix,iy,1,1)
            !density(ix,iy,1)=-0.5d0*(-3.d0*vl+4*pot_layerl-pot_layerl2)* hgzinv
            density(ix,iy,1)=-(-25.d0/12.d0*vl+4.d0*pot_layerl-3.d0*pot_layerl2+4.d0/3.d0*pot_layerl3-0.25d0*pot_layerl4)* hgzinv
            t=t+ density(ix,iy,1)
            pot_layeru4=ewald_p3d%poisson_p3d%pots(ix,iy,npu-4)+ewald_p3d%poisson_p3d%rho(ix,iy,npu-4)+pot_short(ix,iy,2,5)
            pot_layeru3=ewald_p3d%poisson_p3d%pots(ix,iy,npu-3)+ewald_p3d%poisson_p3d%rho(ix,iy,npu-3)+pot_short(ix,iy,2,4)
            pot_layeru2=ewald_p3d%poisson_p3d%pots(ix,iy,npu-2)+ewald_p3d%poisson_p3d%rho(ix,iy,npu-2)+pot_short(ix,iy,2,3)
            pot_layeru =ewald_p3d%poisson_p3d%pots(ix,iy,npu-1)+ewald_p3d%poisson_p3d%rho(ix,iy,npu-1)+pot_short(ix,iy,2,2)
                    vu =ewald_p3d%poisson_p3d%pots(ix,iy,npu  )+ewald_p3d%poisson_p3d%rho(ix,iy,npu  )+pot_short(ix,iy,2,1)
            !density(ix,iy,2)=0.5d0*(3.d0*vu-4.d0*pot_layeru+pot_layeru2)* hgzinv
            density(ix,iy,2)=-(-25.d0/12.d0*vu+4.d0*pot_layeru-3.d0*pot_layeru2+4.d0/3.d0*pot_layeru3-0.25d0*pot_layeru4)* hgzinv
            tt=tt+ density(ix,iy,2)
    enddo
    enddo
  !  do iy=1,ewald_p3d%poisson_p3d%ngpy
  !  do ix=1,ewald_p3d%poisson_p3d%ngpx
  !      write(1000,*)density(ix,iy,1)
  !      write(1001,*)density(ix,iy,2)
  !  enddo
  !  write(1000,*)
  !  write(1001,*)
  !  enddo
    t=t*ewald_p3d%hgx*ewald_p3d%hgy
    tt=tt*ewald_p3d%hgx*ewald_p3d%hgy
    if (trim(parini%bias_field)=='yes') then
        t =t -E/(4*pi)*ewald_p3d%cell(1)*ewald_p3d%cell(2)
        tt=tt+E/(4*pi)*ewald_p3d%cell(1)*ewald_p3d%cell(2)
    endif
    write(*,'(a,es25.13)')'charge on lower plane' ,t
    write(*,'(a,es25.13)')'charge on upper plane',tt
    write(77,'(3es25.13)')vu-vl,t,tt
    vu=ewald_p3d%vu
    vl=ewald_p3d%vl
end subroutine surface_charge
!*****************************************************************************************
!This subroutine determines the limits of grids in a sphere.
subroutine determine_limitsphere(ewald_p3d,mboundg,mboundgy,nbgpx,nbgpy,nbgpz)
    use mod_interface
    use mod_electrostatics, only: typ_ewald_p3d
    implicit none
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
    !local variables
    integer:: ix, iy, iz, mboundg(2,-nbgpy:nbgpy,-nbgpz:nbgpz), mboundgy(2,-nbgpz:nbgpz)
    integer:: nbgpx, nbgpy, nbgpz

    real(8):: rgcut, rgcutsq
    rgcut=max(ewald_p3d%hgx*nbgpx,ewald_p3d%hgy*nbgpy,ewald_p3d%hgz*nbgpz)
    rgcutsq=rgcut**2
    do iz=-nbgpz,nbgpz
        do iy=-nbgpy,nbgpy
            mboundg(1,iy,iz)=1
            mboundg(2,iy,iz)=0
        enddo
    enddo
    do iz=0,nbgpz
    do iy=-nbgpy,nbgpy
    do ix=0,nbgpx
        if(ix**2*ewald_p3d%hgx**2+iy**2*ewald_p3d%hgy**2+iz**2*ewald_p3d%hgz**2<=rgcutsq) then
            mboundg(1,iy,iz)=-ix
            mboundg(2,iy,iz)=ix
        endif
    enddo
    enddo
    enddo
    do iz=-nbgpz,-1
        do iy=-nbgpy,nbgpy
            mboundg(1,iy,iz)=mboundg(1,iy,-iz)
            mboundg(2,iy,iz)=mboundg(2,iy,-iz)
        enddo
    enddo
    do iz=-nbgpz,nbgpz
        mboundgy(1,iz)=1 !-mlimnb
        mboundgy(2,iz)=0 !mlimnb
    enddo
    do iz=-nbgpz,nbgpz
    do iy=1,nbgpy
    if(mboundg(1,iy,iz)<=mboundg(2,iy,iz)) then
        mboundgy(1,iz)=-iy
        mboundgy(2,iz)=iy
    endif
    enddo
    enddo


end subroutine determine_limitsphere
!***********************************************************************************************************************
