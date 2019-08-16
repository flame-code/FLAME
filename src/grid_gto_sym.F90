!*****************************************************************************************
subroutine put_gto_sym(parini,bc,reset,nat,rxyz,qat,gw,rgcut,ngx,ngy,ngz,hgrid,rho)
    use mod_atoms, only: typ_atoms
    use mod_parini, only: typ_parini
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: bc
    logical, intent(in):: reset
    integer, intent(in):: nat
    real(8), intent(in):: rxyz(3,nat)
    real(8), intent(in):: qat(nat)
    real(8), intent(in):: gw(nat)
    real(8), intent(in):: rgcut
    integer, intent(in):: ngx, ngy, ngz
    real(8), intent(in):: hgrid(3,3)
    real(8), intent(inout):: rho(ngx,ngy,ngz)
    !local variables
    !work arrays to save the values of one dimensional gaussian function.
    real(8):: pi
    real(8):: facqiat, fac
    real(8):: vol
    real(8):: hxx, hxy, hxz, hyx, hyy, hyz, hzx, hzy, hzz
    real(8):: hyx_g, hyy_g, hyz_g, hzx_g, hzy_g, hzz_g
    real(8):: hrxinv, hryinv, hrzinv
    real(8):: dmx, dmy, dmz, dmsq, gwsq_inv
    real(8):: ximg, yimg, zimg , cv(3,3)
    integer:: imgx, imgy, imgz
    integer:: iat, igx, igy, igz, jgx, jgy, jgz, igzysq
    integer:: nlimsq, iix
    integer:: nbgx, nbgy, nbgz, nagx, nagy, nagz
    real(8), allocatable:: wa(:,:,:)
    real(8), allocatable:: ratred(:,:)
    real(8), allocatable:: exponentval(:), expval(:)
    real(8):: sqpi, tt
    integer:: ibcz

    allocate(ratred(3,nat))
    hxx=hgrid(1,1) ; hxy=hgrid(2,1) ; hxz=hgrid(3,1)
    hyx=hgrid(1,2) ; hyy=hgrid(2,2) ; hyz=hgrid(3,2)
    hzx=hgrid(1,3) ; hzy=hgrid(2,3) ; hzz=hgrid(3,3)
    cv(1,1)=hxx*ngx ; cv(2,1)=hxy*ngx ; cv(3,1)=hxz*ngx
    cv(1,2)=hyx*ngy ; cv(2,2)=hyy*ngy ; cv(3,2)=hyz*ngy
    cv(1,3)=hzx*ngz ; cv(2,3)=hzy*ngz ; cv(3,3)=hzz*ngz
    call init_grid_param(nat,rxyz,cv,rgcut,ngx,ngy,ngz,ratred,vol,nlimsq,nagx,nagy,nagz,nbgx,nbgy,nbgz)
    !-------------------------------------------------------
    ibcz=0
    if(trim(bc)=='bulk') then
        !variables already are set to correct value
    elseif(trim(bc)=='slab') then
        ibcz=1
        nagz=0
        if(hgrid(1,3)/=0.d0 .or. hgrid(2,3)/=0.d0) then
            write(*,'(a,2es14.5)') &
                'ERROR: for slab BC, third cell vector must be along z direction', &
                hgrid(1,3),hgrid(2,3)
            stop
        endif
    else
        write(*,'(2a)') 'ERROR: using this routine is meaningful only for bulk and slab BCs, bc= ',trim(bc)
        stop
    endif
    !-------------------------------------------------------
    wa=f_malloc0([1-nagx.to.ngx+nagx,1-nagy.to.ngy+nagy,1-nagz.to.ngz+nagz],id='wa')
    allocate(exponentval(-nbgx:nbgx),expval(-nbgx:nbgx))
   !  allocate(wa(1-nagx:ngx+nagx,1-nagy:ngy+nagy,1-nagz:ngz+nagz))

    hrxinv=real(ngx,8) !inverse of grid spacing in reduced coordinates
    hryinv=real(ngy,8) !inverse of grid spacing in reduced coordinates
    hrzinv=real(ngz,8) !inverse of grid spacing in reduced coordinates
    !-------------------------------------------------------
    pi=4.d0*atan(1.d0)
    sqpi=sqrt(pi)
    !-------------------------------------------------------
    do iat=1,nat
        gwsq_inv=1.d0/gw(iat)**2
        fac=1.d0/(gw(iat)*sqpi)**3
        imgx=nint(ratred(1,iat)*hrxinv)+1
        imgy=nint(ratred(2,iat)*hryinv)+1
        imgz=nint(ratred(3,iat)*hrzinv)+1+nbgz*ibcz
        facqiat=fac*qat(iat)
        do igz=-nbgz,nbgz
            jgz=igz+imgz
            hzx_g=(jgz-1)*hzx
            hzy_g=(jgz-1)*hzy
            hzz_g=(jgz-1-nbgz*ibcz)*hzz
            do igy=-nbgy,nbgy
                igzysq=igz**2+igy**2
                if(igzysq>nlimsq) cycle
                jgy=igy+imgy
                hyx_g=(jgy-1)*hyx
                hyy_g=(jgy-1)*hyy
                hyz_g=(jgy-1)*hyz
                tt=nlimsq-igzysq
                iix=min(floor(sqrt(tt)),nbgx)
                do igx=-iix,iix
                    jgx=igx+imgx
                    ximg=(jgx-1)*hxx+hyx_g+hzx_g    
                    yimg=(jgx-1)*hxy+hyy_g+hzy_g
                    zimg=(jgx-1)*hxz+hyz_g+hzz_g
                    dmx=ximg-rxyz(1,iat)
                    dmy=yimg-rxyz(2,iat)
                    dmz=zimg-rxyz(3,iat)
                    dmsq=dmx**2+dmy**2+dmz**2
                    exponentval(igx)=-dmsq*gwsq_inv
                enddo
#ifdef HAVE_MKL
                call vdexp(2*iix+1,exponentval(-iix),expval(-iix))
#else
                do igx=-iix,iix
                    expval(igx)=exp(exponentval(igx))
                enddo
#endif
                wa(imgx-iix:imgx+iix,jgy,jgz)=wa(imgx-iix:imgx+iix,jgy,jgz)+facqiat*expval(-iix:iix)
            enddo
        enddo
    enddo
    !---------------------------------------------------------------------------
    if(reset) then
        !if the input array of charge density does not contain any previous value
        !wanted to be preserved.
        rho = 0.d0
    endif
    call charge_back_to_cell(ngx,ngy,ngz,nagx,nagy,nagz,0,wa,rho)
    !---------------------------------------------------------------------------
    deallocate(ratred)
    deallocate(exponentval,expval)
    call f_free(wa)
end subroutine put_gto_sym
!*****************************************************************************************
subroutine rqgrad_gto_sym(parini,bc,nat,rxyz,qat,gw,rgcut,lda,ngx,ngy,ngz,hgrid,pot,rgrad,qgrad)
    use mod_atoms, only: typ_atoms
    use mod_parini, only: typ_parini
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: bc
    integer, intent(in):: nat
    real(8), intent(in):: rxyz(3,nat)
    real(8), intent(in):: qat(nat)
    real(8), intent(in):: gw(nat)
    real(8), intent(in):: rgcut
    integer, intent(in):: lda, ngx, ngy, ngz
    real(8), intent(in):: hgrid(3,3)
    real(8), intent(in):: pot(lda,ngy,ngz)
    real(8), intent(out):: rgrad(3,nat), qgrad(nat)
    !local variables
    !work arrays to save the values of one dimensional gaussian function.
    real(8):: pi
    real(8):: facqiat, fac
    real(8):: vol
    real(8):: hxx, hxy, hxz, hyx, hyy, hyz, hzx, hzy, hzz
    real(8):: hyx_g, hyy_g, hyz_g, hzx_g, hzy_g, hzz_g
    real(8):: hrxinv, hryinv, hrzinv
    real(8):: vol_voxel, ttx, tty, ttz, ttq, tt1
    real(8):: dmsq, gwsq_inv
    real(8):: ximg, yimg, zimg, cv(3,3)
    integer:: imgx, imgy, imgz
    integer:: iat, igx, igy, igz, jgx, jgy, jgz, igzysq
    integer:: nlimsq, iix
    integer:: nbgx, nbgy, nbgz, nagx, nagy, nagz
    real(8), allocatable:: wa(:,:,:)
    real(8), allocatable:: ratred(:,:)
    real(8), allocatable:: exponentval(:), expval(:)
    real(8), allocatable:: dmxarr(:), dmyarr(:), dmzarr(:)
    real(8):: sqpi, tt
    integer:: ibcz

    allocate(ratred(3,nat))
    hxx=hgrid(1,1) ; hxy=hgrid(2,1) ; hxz=hgrid(3,1)
    hyx=hgrid(1,2) ; hyy=hgrid(2,2) ; hyz=hgrid(3,2)
    hzx=hgrid(1,3) ; hzy=hgrid(2,3) ; hzz=hgrid(3,3)
    cv(1,1)=hxx*ngx ; cv(2,1)=hxy*ngx ; cv(3,1)=hxz*ngx
    cv(1,2)=hyx*ngy ; cv(2,2)=hyy*ngy ; cv(3,2)=hyz*ngy
    cv(1,3)=hzx*ngz ; cv(2,3)=hzy*ngz ; cv(3,3)=hzz*ngz
    call init_grid_param(nat,rxyz,cv,rgcut,ngx,ngy,ngz,ratred,vol,nlimsq,nagx,nagy,nagz,nbgx,nbgy,nbgz)
    !-------------------------------------------------------
    ibcz=0
    if(trim(bc)=='bulk') then
        !variables already are set to correct value
    elseif(trim(bc)=='slab') then
        ibcz=1
        nagz=0
        if(hgrid(1,3)/=0.d0 .or. hgrid(2,3)/=0.d0) then
            write(*,'(a,2es14.5)') &
                'ERROR: for slab BC, third cell vector must be along z direction', &
                hgrid(1,3),hgrid(2,3)
            stop
        endif
    else
        write(*,'(2a)') 'ERROR: using this routine is meaningful only for bulk and slab BCs, bc= ',trim(bc)
        stop
    endif
    !-------------------------------------------------------
    wa=f_malloc0([1-nagx.to.ngx+nagx,1-nagy.to.ngy+nagy,1-nagz.to.ngz+nagz],id='wa')
    allocate(exponentval(-nbgx:nbgx),expval(-nbgx:nbgx))
    allocate(dmxarr(-nbgx:nbgx),dmyarr(-nbgx:nbgx),dmzarr(-nbgx:nbgx))
   !  allocate(wa(1-nagx:ngx+nagx,1-nagy:ngy+nagy,1-nagz:ngz+nagz))

    hrxinv=real(ngx,8) !inverse of grid spacing in reduced coordinates
    hryinv=real(ngy,8) !inverse of grid spacing in reduced coordinates
    hrzinv=real(ngz,8) !inverse of grid spacing in reduced coordinates
    !-------------------------------------------------------
    pi=4.d0*atan(1.d0)
    sqpi=sqrt(pi)
    !---------------------------------------------------------------------------
    call potential_on_extended_grid(lda,ngx,ngy,ngz,nagx,nagy,nagz,0,pot,wa)
    !-------------------------------------------------------
    vol_voxel=vol/(ngx*ngy*ngz)
    do iat=1,nat
        gwsq_inv=1.d0/gw(iat)**2
        fac=1.d0/(gw(iat)*sqpi)**3
        imgx=nint(ratred(1,iat)*hrxinv)+1
        imgy=nint(ratred(2,iat)*hryinv)+1
        imgz=nint(ratred(3,iat)*hrzinv)+1+nbgz*ibcz
        facqiat=fac*qat(iat)
        ttq=0.d0
        ttx=0.d0
        tty=0.d0
        ttz=0.d0
        do igz=-nbgz,nbgz
            jgz=igz+imgz
            hzx_g=(jgz-1)*hzx
            hzy_g=(jgz-1)*hzy
            hzz_g=(jgz-1-nbgz*ibcz)*hzz
            do igy=-nbgy,nbgy
                igzysq=igz**2+igy**2
                if(igzysq>nlimsq) cycle
                jgy=igy+imgy
                hyx_g=(jgy-1)*hyx
                hyy_g=(jgy-1)*hyy
                hyz_g=(jgy-1)*hyz
                tt=nlimsq-igzysq
                iix=min(floor(sqrt(tt)),nbgx)
                do igx=-iix,iix
                    jgx=igx+imgx
                    ximg=(jgx-1)*hxx+hyx_g+hzx_g    
                    yimg=(jgx-1)*hxy+hyy_g+hzy_g
                    zimg=(jgx-1)*hxz+hyz_g+hzz_g
                    dmxarr(igx)=ximg-rxyz(1,iat)
                    dmyarr(igx)=yimg-rxyz(2,iat)
                    dmzarr(igx)=zimg-rxyz(3,iat)
                    dmsq=dmxarr(igx)**2+dmyarr(igx)**2+dmzarr(igx)**2
                    exponentval(igx)=-dmsq*gwsq_inv
                enddo
#ifdef HAVE_MKL
                call vdexp(2*iix+1,exponentval(-iix),expval(-iix))
#else
                do igx=-iix,iix
                    expval(igx)=exp(exponentval(igx))
                enddo
#endif
                do igx=-iix,iix
                    ttq=ttq+fac*expval(igx)*wa(igx+imgx,jgy,jgz)
                    tt1=facqiat*expval(igx)*wa(igx+imgx,jgy,jgz)*(2.d0*gwsq_inv)
                    ttx=ttx+tt1*dmxarr(igx)
                    tty=tty+tt1*dmyarr(igx)
                    ttz=ttz+tt1*dmzarr(igx)
                enddo
            enddo
        enddo
        qgrad(iat)=qgrad(iat)+ttq*vol_voxel
        rgrad(1,iat)=rgrad(1,iat)+ttx*vol_voxel
        rgrad(2,iat)=rgrad(2,iat)+tty*vol_voxel
        rgrad(3,iat)=rgrad(3,iat)+ttz*vol_voxel
    enddo
    !---------------------------------------------------------------------------
    deallocate(ratred)
    deallocate(exponentval,expval)
    deallocate(dmxarr,dmyarr,dmzarr)
    call f_free(wa)
end subroutine rqgrad_gto_sym
!*****************************************************************************************
subroutine force_gto_sym(parini,bc,nat,rxyz,qat,gw,rgcut,lda,ngx,ngy,ngz,hgrid,pot,fat)
    use mod_atoms, only: typ_atoms
    use mod_parini, only: typ_parini
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: bc
    integer, intent(in):: nat
    real(8), intent(in):: rxyz(3,nat)
    real(8), intent(in):: qat(nat)
    real(8), intent(in):: gw(nat)
    real(8), intent(in):: rgcut
    integer, intent(in):: lda, ngx, ngy, ngz
    real(8), intent(in):: hgrid(3,3)
    real(8), intent(in):: pot(lda,ngy,ngz)
    real(8), intent(out):: fat(3,nat)
    !local variables
    !work arrays to save the values of one dimensional gaussian function.
    real(8):: pi
    real(8):: facqiat, fac
    real(8):: vol
    real(8):: hxx, hxy, hxz, hyx, hyy, hyz, hzx, hzy, hzz
    real(8):: hyx_g, hyy_g, hyz_g, hzx_g, hzy_g, hzz_g
    real(8):: hrxinv, hryinv, hrzinv
    real(8):: vol_voxel, ttx, tty, ttz, tt1
    real(8):: dmsq, gwsq_inv
    real(8):: ximg, yimg, zimg, cv(3,3)
    integer:: imgx, imgy, imgz
    integer:: iat, igx, igy, igz, jgx, jgy, jgz, igzysq
    integer:: nlimsq, iix
    integer:: nbgx, nbgy, nbgz, nagx, nagy, nagz
    real(8), allocatable:: wa(:,:,:)
    real(8), allocatable:: ratred(:,:)
    real(8), allocatable:: exponentval(:), expval(:)
    real(8), allocatable:: dmxarr(:), dmyarr(:), dmzarr(:)
    real(8):: sqpi, tt
    integer:: ibcz

    allocate(ratred(3,nat))
    hxx=hgrid(1,1) ; hxy=hgrid(2,1) ; hxz=hgrid(3,1)
    hyx=hgrid(1,2) ; hyy=hgrid(2,2) ; hyz=hgrid(3,2)
    hzx=hgrid(1,3) ; hzy=hgrid(2,3) ; hzz=hgrid(3,3)
    cv(1,1)=hxx*ngx ; cv(2,1)=hxy*ngx ; cv(3,1)=hxz*ngx
    cv(1,2)=hyx*ngy ; cv(2,2)=hyy*ngy ; cv(3,2)=hyz*ngy
    cv(1,3)=hzx*ngz ; cv(2,3)=hzy*ngz ; cv(3,3)=hzz*ngz
    call init_grid_param(nat,rxyz,cv,rgcut,ngx,ngy,ngz,ratred,vol,nlimsq,nagx,nagy,nagz,nbgx,nbgy,nbgz)
    !-------------------------------------------------------
    ibcz=0
    if(trim(bc)=='bulk') then
        !variables already are set to correct value
    elseif(trim(bc)=='slab') then
        ibcz=1
        nagz=0
        if(hgrid(1,3)/=0.d0 .or. hgrid(2,3)/=0.d0) then
            write(*,'(a,2es14.5)') &
                'ERROR: for slab BC, third cell vector must be along z direction', &
                hgrid(1,3),hgrid(2,3)
            stop
        endif
    else
        write(*,'(2a)') 'ERROR: using this routine is meaningful only for bulk and slab BCs, bc= ',trim(bc)
        stop
    endif
    !-------------------------------------------------------
    wa=f_malloc0([1-nagx.to.ngx+nagx,1-nagy.to.ngy+nagy,1-nagz.to.ngz+nagz],id='wa')
    allocate(exponentval(-nbgx:nbgx),expval(-nbgx:nbgx))
    allocate(dmxarr(-nbgx:nbgx),dmyarr(-nbgx:nbgx),dmzarr(-nbgx:nbgx))
   !  allocate(wa(1-nagx:ngx+nagx,1-nagy:ngy+nagy,1-nagz:ngz+nagz))

    hrxinv=real(ngx,8) !inverse of grid spacing in reduced coordinates
    hryinv=real(ngy,8) !inverse of grid spacing in reduced coordinates
    hrzinv=real(ngz,8) !inverse of grid spacing in reduced coordinates
    !-------------------------------------------------------
    pi=4.d0*atan(1.d0)
    sqpi=sqrt(pi)
    !---------------------------------------------------------------------------
    call potential_on_extended_grid(lda,ngx,ngy,ngz,nagx,nagy,nagz,0,pot,wa)
    !-------------------------------------------------------
    vol_voxel=vol/(ngx*ngy*ngz)
    do iat=1,nat
        gwsq_inv=1.d0/gw(iat)**2
        fac=1.d0/(gw(iat)*sqpi)**3
        imgx=nint(ratred(1,iat)*hrxinv)+1
        imgy=nint(ratred(2,iat)*hryinv)+1
        imgz=nint(ratred(3,iat)*hrzinv)+1+nbgz*ibcz
        facqiat=fac*qat(iat)
        ttx=0.d0
        tty=0.d0
        ttz=0.d0
        do igz=-nbgz,nbgz
            jgz=igz+imgz
            hzx_g=(jgz-1)*hzx
            hzy_g=(jgz-1)*hzy
            hzz_g=(jgz-1-nbgz*ibcz)*hzz
            do igy=-nbgy,nbgy
                igzysq=igz**2+igy**2
                if(igzysq>nlimsq) cycle
                jgy=igy+imgy
                hyx_g=(jgy-1)*hyx
                hyy_g=(jgy-1)*hyy
                hyz_g=(jgy-1)*hyz
                tt=nlimsq-igzysq
                iix=min(floor(sqrt(tt)),nbgx)
                do igx=-iix,iix
                    jgx=igx+imgx
                    ximg=(jgx-1)*hxx+hyx_g+hzx_g    
                    yimg=(jgx-1)*hxy+hyy_g+hzy_g
                    zimg=(jgx-1)*hxz+hyz_g+hzz_g
                    dmxarr(igx)=ximg-rxyz(1,iat)
                    dmyarr(igx)=yimg-rxyz(2,iat)
                    dmzarr(igx)=zimg-rxyz(3,iat)
                    dmsq=dmxarr(igx)**2+dmyarr(igx)**2+dmzarr(igx)**2
                    exponentval(igx)=-dmsq*gwsq_inv
                enddo
#ifdef HAVE_MKL
                call vdexp(2*iix+1,exponentval(-iix),expval(-iix))
#else
                do igx=-iix,iix
                    expval(igx)=exp(exponentval(igx))
                enddo
#endif
                do igx=-iix,iix
                    tt1=facqiat*expval(igx)*wa(igx+imgx,jgy,jgz)*(2.d0*gwsq_inv)
                    ttx=ttx+tt1*dmxarr(igx)
                    tty=tty+tt1*dmyarr(igx)
                    ttz=ttz+tt1*dmzarr(igx)
                enddo
            enddo
        enddo
        fat(1,iat)=fat(1,iat)-ttx*vol_voxel
        fat(2,iat)=fat(2,iat)-tty*vol_voxel
        fat(3,iat)=fat(3,iat)-ttz*vol_voxel
    enddo
    !---------------------------------------------------------------------------
    deallocate(ratred)
    deallocate(exponentval,expval)
    deallocate(dmxarr,dmyarr,dmzarr)
    call f_free(wa)
end subroutine force_gto_sym
!*****************************************************************************************
subroutine gwrqgrad_gto_sym(parini,bc,nat,rxyz,qat,gw,rgcut,lda,ngx,ngy,ngz,hgrid,pot,rgrad,qgrad,agrad)
    use mod_atoms, only: typ_atoms
    use mod_parini, only: typ_parini
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: bc
    integer, intent(in):: nat
    real(8), intent(in):: rxyz(3,nat)
    real(8), intent(in):: qat(nat) 
    real(8), intent(in):: gw(nat)
    real(8), intent(in):: rgcut
    integer, intent(in):: lda, ngx, ngy, ngz
    real(8), intent(in):: hgrid(3,3)
    real(8), intent(in):: pot(lda,ngy,ngz)
    real(8), intent(out):: rgrad(3,nat), qgrad(nat), agrad(nat)
    !local variables
    !work arrays to save the values of one dimensional gaussian function.
    real(8):: pi
    real(8):: facqiat, fac
    real(8):: vol
    real(8):: hxx, hxy, hxz, hyx, hyy, hyz, hzx, hzy, hzz
    real(8):: hyx_g, hyy_g, hyz_g, hzx_g, hzy_g, hzz_g
    real(8):: hrxinv, hryinv, hrzinv
    real(8):: vol_voxel, ttx, tty, ttz, ttq, tt1, tta
    real(8):: dmsq, gwsq_inv, gw_inv
    real(8):: ximg, yimg, zimg, cv(3,3)
    integer:: imgx, imgy, imgz
    integer:: iat, igx, igy, igz, jgx, jgy, jgz, igzysq
    integer:: nlimsq, iix
    integer:: nbgx, nbgy, nbgz, nagx, nagy, nagz
    real(8), allocatable:: wa(:,:,:)
    real(8), allocatable:: ratred(:,:)
    real(8), allocatable:: exponentval(:), expval(:)
    real(8), allocatable:: dmxarr(:), dmyarr(:), dmzarr(:)
    real(8):: sqpi, tt
    integer:: ibcz

    allocate(ratred(3,nat))
    hxx=hgrid(1,1) ; hxy=hgrid(2,1) ; hxz=hgrid(3,1)
    hyx=hgrid(1,2) ; hyy=hgrid(2,2) ; hyz=hgrid(3,2)
    hzx=hgrid(1,3) ; hzy=hgrid(2,3) ; hzz=hgrid(3,3)
    cv(1,1)=hxx*ngx ; cv(2,1)=hxy*ngx ; cv(3,1)=hxz*ngx
    cv(1,2)=hyx*ngy ; cv(2,2)=hyy*ngy ; cv(3,2)=hyz*ngy
    cv(1,3)=hzx*ngz ; cv(2,3)=hzy*ngz ; cv(3,3)=hzz*ngz
    call init_grid_param(nat,rxyz,cv,rgcut,ngx,ngy,ngz,ratred,vol,nlimsq,nagx,nagy,nagz,nbgx,nbgy,nbgz)
    !-------------------------------------------------------
    ibcz=0
    if(trim(bc)=='bulk') then
        !variables already are set to correct value
    elseif(trim(bc)=='slab') then
        ibcz=1
        nagz=0
        if(hgrid(1,3)/=0.d0 .or. hgrid(2,3)/=0.d0) then
            write(*,'(a,2es14.5)') &
                'ERROR: for slab BC, third cell vector must be along z direction', &
                hgrid(1,3),hgrid(2,3)
            stop
        endif
    else
        write(*,'(2a)') 'ERROR: using this routine is meaningful only for bulk and slab BCs, bc= ',trim(bc)
        stop
    endif
    !-------------------------------------------------------
    wa=f_malloc0([1-nagx.to.ngx+nagx,1-nagy.to.ngy+nagy,1-nagz.to.ngz+nagz],id='wa')
    allocate(exponentval(-nbgx:nbgx),expval(-nbgx:nbgx))
    allocate(dmxarr(-nbgx:nbgx),dmyarr(-nbgx:nbgx),dmzarr(-nbgx:nbgx))
   !  allocate(wa(1-nagx:ngx+nagx,1-nagy:ngy+nagy,1-nagz:ngz+nagz))

    hrxinv=real(ngx,8) !inverse of grid spacing in reduced coordinates
    hryinv=real(ngy,8) !inverse of grid spacing in reduced coordinates
    hrzinv=real(ngz,8) !inverse of grid spacing in reduced coordinates
    !-------------------------------------------------------
    pi=4.d0*atan(1.d0)
    sqpi=sqrt(pi)
    !---------------------------------------------------------------------------
    call potential_on_extended_grid(lda,ngx,ngy,ngz,nagx,nagy,nagz,0,pot,wa)
    !-------------------------------------------------------
    vol_voxel=vol/(ngx*ngy*ngz)
    do iat=1,nat
        gw_inv = 1.d0/gw(iat)
        gwsq_inv=1.d0/gw(iat)**2
        fac=1.d0/(gw(iat)*sqpi)**3
        imgx=nint(ratred(1,iat)*hrxinv)+1
        imgy=nint(ratred(2,iat)*hryinv)+1
        imgz=nint(ratred(3,iat)*hrzinv)+1+nbgz*ibcz
        facqiat=fac*qat(iat)
        ttq=0.d0
        ttx=0.d0
        tty=0.d0
        ttz=0.d0
        !------------------------
        tta=0.d0
        do igz=-nbgz,nbgz
            jgz=igz+imgz
            hzx_g=(jgz-1)*hzx
            hzy_g=(jgz-1)*hzy
            hzz_g=(jgz-1-nbgz*ibcz)*hzz
            do igy=-nbgy,nbgy
                igzysq=igz**2+igy**2
                if(igzysq>nlimsq) cycle
                jgy=igy+imgy
                hyx_g=(jgy-1)*hyx
                hyy_g=(jgy-1)*hyy
                hyz_g=(jgy-1)*hyz
                tt=nlimsq-igzysq
                iix=min(floor(sqrt(tt)),nbgx)
                do igx=-iix,iix
                    jgx=igx+imgx
                    ximg=(jgx-1)*hxx+hyx_g+hzx_g    
                    yimg=(jgx-1)*hxy+hyy_g+hzy_g
                    zimg=(jgx-1)*hxz+hyz_g+hzz_g
                    dmxarr(igx)=ximg-rxyz(1,iat)
                    dmyarr(igx)=yimg-rxyz(2,iat)
                    dmzarr(igx)=zimg-rxyz(3,iat)
                    dmsq=dmxarr(igx)**2+dmyarr(igx)**2+dmzarr(igx)**2
                    exponentval(igx)=-dmsq*gwsq_inv
                enddo
#ifdef HAVE_MKL
                call vdexp(2*iix+1,exponentval(-iix),expval(-iix))
#else
                do igx=-iix,iix
                    expval(igx)=exp(exponentval(igx))
                enddo
#endif
                do igx=-iix,iix
                    ttq=ttq+fac*expval(igx)*wa(igx+imgx,jgy,jgz)
                    tt1=facqiat*expval(igx)*wa(igx+imgx,jgy,jgz)*(2.d0*gwsq_inv)
                    ttx=ttx+tt1*dmxarr(igx)
                    tty=tty+tt1*dmyarr(igx)
                    ttz=ttz+tt1*dmzarr(igx)
                    !-----------------------
                    tta=tta+(-3.d0 - (2.d0*exponentval(igx)))*facqiat*gw_inv*expval(igx)*wa(igx+imgx,jgy,jgz)

                    !_______________________
                enddo
            enddo
        enddo
        qgrad(iat)=qgrad(iat)+ttq*vol_voxel
        rgrad(1,iat)=rgrad(1,iat)+ttx*vol_voxel
        rgrad(2,iat)=rgrad(2,iat)+tty*vol_voxel
        rgrad(3,iat)=rgrad(3,iat)+ttz*vol_voxel
        !-----------------------
        agrad(iat) = agrad(iat)+tta*vol_voxel
        !_______________________
    enddo
    !---------------------------------------------------------------------------
    deallocate(ratred)
    deallocate(exponentval,expval)
    deallocate(dmxarr,dmyarr,dmzarr)
    call f_free(wa)
end subroutine gwrqgrad_gto_sym
!*****************************************************************************************
subroutine rhograd_gto_sym(parini,bc,reset,nat,rxyz,cv,qat,gw,rgcut,ngx,ngy,ngz,rho,rho_q_par,rho_a_par)
    use mod_atoms, only: typ_atoms
    use mod_parini, only: typ_parini
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: bc
    logical, intent(in):: reset
    integer, intent(in):: nat
    real(8), intent(in):: rxyz(3,nat)
    real(8), intent(in):: cv(3,3)
    real(8), intent(in):: qat(nat)
    real(8), intent(in):: gw(nat)
    real(8), intent(in):: rgcut
    integer, intent(in):: ngx, ngy, ngz
    real(8), intent(inout):: rho(ngx,ngy,ngz),rho_a_par(ngx,ngy,ngz),rho_q_par(ngx,ngy,ngz)
    !local variables
    !work arrays to save the values of one dimensional gaussian function.
    real(8):: pi
    real(8):: facqiat, fac
    real(8):: vol
    real(8):: hxx, hxy, hxz, hyx, hyy, hyz, hzx, hzy, hzz
    real(8):: hyx_g, hyy_g, hyz_g, hzx_g, hzy_g, hzz_g
    real(8):: hrxinv, hryinv, hrzinv
    real(8):: dmx, dmy, dmz, dmsq, gwsq_inv, gw_inv, gwcub_inv
    real(8):: ximg, yimg, zimg
    integer:: imgx, imgy, imgz
    integer:: iat, igx, igy, igz, jgx, jgy, jgz, igzysq
    integer:: nlimsq, iix, iiy, iiz
    integer:: nbgx, nbgy, nbgz, nagx, nagy, nagz
    integer:: ifinalx, igxs, igxf 
    real(8), allocatable:: wa(:,:,:),wb(:,:,:),wc(:,:,:)
    real(8), allocatable:: ratred(:,:)
    real(8), allocatable:: exponentval(:), expval(:)
    real(8):: sqpi, tt
    integer:: istartx, istarty, istartz

    allocate(ratred(3,nat))
    call init_grid_param(nat,rxyz,cv,rgcut,ngx,ngy,ngz,ratred,vol,nlimsq,nagx,nagy,nagz,nbgx,nbgy,nbgz)
    hxx=cv(1,1)/ngx ; hxy=cv(2,1)/ngx ; hxz=cv(3,1)/ngx
    hyx=cv(1,2)/ngy ; hyy=cv(2,2)/ngy ; hyz=cv(3,2)/ngy
    hzx=cv(1,3)/ngz ; hzy=cv(2,3)/ngz ; hzz=cv(3,3)/ngz
    if(trim(bc)/='bulk') stop 'ERROR: rhograd_gto_sym works only for bulk BC'
    wa=f_malloc0([1-nagx.to.ngx+nagx,1-nagy.to.ngy+nagy,1-nagz.to.ngz+nagz],id='wa')
    wb=f_malloc0([1-nagx.to.ngx+nagx,1-nagy.to.ngy+nagy,1-nagz.to.ngz+nagz],id='wb')
    wc=f_malloc0([1-nagx.to.ngx+nagx,1-nagy.to.ngy+nagy,1-nagz.to.ngz+nagz],id='wc')
    allocate(exponentval(-nbgx:nbgx),expval(-nbgx:nbgx))
   !  allocate(wa(1-nagx:ngx+nagx,1-nagy:ngy+nagy,1-nagz:ngz+nagz))

    hrxinv=real(ngx,8) !inverse of grid spacing in reduced coordinates
    hryinv=real(ngy,8) !inverse of grid spacing in reduced coordinates
    hrzinv=real(ngz,8) !inverse of grid spacing in reduced coordinates
    !-------------------------------------------------------
    pi=4.d0*atan(1.d0)
    sqpi=sqrt(pi)
    !-------------------------------------------------------
    do iat=1,nat
        gw_inv=1.d0/gw(iat)
        gwsq_inv=gw_inv*gw_inv
        gwcub_inv=gw_inv*gwsq_inv
        fac=1.d0/(gw(iat)*sqpi)**3
        facqiat=fac*qat(iat)
        imgx=nint(ratred(1,iat)*hrxinv)+1
        imgy=nint(ratred(2,iat)*hryinv)+1
        imgz=nint(ratred(3,iat)*hrzinv)+1
        do igz=-nbgz,nbgz
            jgz=igz+imgz
            hzx_g=(jgz-1)*hzx
            hzy_g=(jgz-1)*hzy
            hzz_g=(jgz-1)*hzz
            do igy=-nbgy,nbgy
                igzysq=igz**2+igy**2
                if(igzysq>nlimsq) cycle
                jgy=igy+imgy
                hyx_g=(jgy-1)*hyx
                hyy_g=(jgy-1)*hyy
                hyz_g=(jgy-1)*hyz
                tt=nlimsq-igzysq
                iix=min(floor(sqrt(tt)),nbgx)
                do igx=-iix,iix
                    jgx=igx+imgx
                    ximg=(jgx-1)*hxx+hyx_g+hzx_g    
                    yimg=(jgx-1)*hxy+hyy_g+hzy_g
                    zimg=(jgx-1)*hxz+hyz_g+hzz_g
                    dmx=ximg-rxyz(1,iat)
                    dmy=yimg-rxyz(2,iat)
                    dmz=zimg-rxyz(3,iat)
                    dmsq=dmx**2+dmy**2+dmz**2
                    exponentval(igx)=-dmsq*gwsq_inv
                enddo
#ifdef HAVE_MKL
                call vdexp(2*iix+1,exponentval(-iix),expval(-iix))
#else
                do igx=-iix,iix
                    expval(igx)=exp(exponentval(igx))
                enddo
#endif
                wa(imgx-iix:imgx+iix,jgy,jgz)=wa(imgx-iix:imgx+iix,jgy,jgz)+facqiat*expval(-iix:iix)
                wb(imgx-iix:imgx+iix,jgy,jgz)=wb(imgx-iix:imgx+iix,jgy,jgz)+fac*expval(-iix:iix)
                wc(imgx-iix:imgx+iix,jgy,jgz)=wc(imgx-iix:imgx+iix,jgy,jgz)+(-2.d0*exponentval(-iix:iix) - 3.d0)*gw_inv*facqiat*expval(-iix:iix)

            enddo
        enddo
    enddo
    !---------------------------------------------------------------------------
    if(reset) then
        !if the input array of charge density does not contain any previous value
        !wanted to be preserved.
        rho = 0.d0
        rho_a_par = 0.d0
        rho_q_par = 0.d0
    endif
        !if the input array of charge density already some value that must be preserved.
    istartx = modulo((1-nagx-1),ngx)+1
    ifinalx = modulo((ngx+nagx-1),ngx)+1
    igxs = 1-nagx+(ngx-istartx+1)
    igxf = ngx+nagx-ifinalx+1
    istarty = modulo((-nagy),ngy)+1
    istartz = modulo((-nagz),ngz)+1
    iiz=istartz-1

    do igz=1-nagz,ngz+nagz
        iiy=istarty-1
        iiz=iiz+1
        if (iiz==ngz+1) iiz=1
        do igy=1-nagy,ngy+nagy
            iiy=iiy+1
            if (iiy==ngy+1) iiy=1
            iix=istartx-1
            do igx=1-nagx,igxs-1
                iix=iix+1
                rho(iix,iiy,iiz)=rho(iix,iiy,iiz)+wa(igx,igy,igz)
                rho_q_par(iix,iiy,iiz)=rho_q_par(iix,iiy,iiz)+wb(igx,igy,igz)
                rho_a_par(iix,iiy,iiz)=rho_a_par(iix,iiy,iiz)+wc(igx,igy,igz)
            enddo
            do igx=igxs,igxf-1,ngx
                do iix=1,ngx
                    jgx=igx+iix-1
                    rho(iix,iiy,iiz)=rho(iix,iiy,iiz)+wa(jgx,igy,igz)
                    rho_q_par(iix,iiy,iiz)=rho_q_par(iix,iiy,iiz)+wb(jgx,igy,igz)
                    rho_a_par(iix,iiy,iiz)=rho_a_par(iix,iiy,iiz)+wc(jgx,igy,igz)
                enddo
            enddo
            iix=0
            do igx=igxf,ngx+nagx
                iix=iix+1
                rho(iix,iiy,iiz)=rho(iix,iiy,iiz)+wa(igx,igy,igz)
                rho_q_par(iix,iiy,iiz)=rho_q_par(iix,iiy,iiz)+wb(igx,igy,igz)
                rho_a_par(iix,iiy,iiz)=rho_a_par(iix,iiy,iiz)+wc(igx,igy,igz)
            enddo
        enddo
    enddo
    !---------------------------------------------------------------------------
    deallocate(ratred)
    deallocate(exponentval,expval)
    call f_free(wa)
end subroutine rhograd_gto_sym
!*****************************************************************************************
