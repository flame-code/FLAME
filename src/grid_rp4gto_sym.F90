!*****************************************************************************************
subroutine put_rp4gto_sym(parini,bc,reset,nat,rxyz,cv,qat,gw,rgcut,ngx,ngy,ngz,rho,rho_q_par,rho_a_par)
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
    real(8):: hxx_g, hxy_g, hxz_g, hyx_g, hyy_g, hyz_g, hzx_g, hzy_g, hzz_g
    real(8):: hrxinv, hryinv, hrzinv
    real(8):: dmx, dmy, dmz, dmsq, gwsq_inv, gw_inv, gwcub_inv
    real(8):: ximg, yimg, zimg
    integer:: imgx, imgy, imgz
    integer:: ncellx, ncelly, ncellz
    integer:: iat, igx, igy, igz, jgx, jgy, jgz, igzysq
    integer:: iii, nlimsq, iix, iiy, iiz
    integer:: nbgx, nbgy, nbgz, nagx, nagy, nagz, nex, ney, nez
    integer:: finalx, igxs, igxf 
    real(8), allocatable:: wa(:,:,:),wb(:,:,:),wc(:,:,:)
    real(8), allocatable:: ratred(:,:)
    real(8), allocatable:: exponentval(:), expval(:),rp4(:)
    real(8):: sqpi, tt
    integer:: istartx, istarty, istartz

    allocate(ratred(3,nat))
    call init_grid_param(nat,rxyz,cv,rgcut,ngx,ngy,ngz,ratred,vol,nlimsq,nagx,nagy,nagz,nbgx,nbgy,nbgz)
    hxx=cv(1,1)/ngx ; hxy=cv(2,1)/ngx ; hxz=cv(3,1)/ngx
    hyx=cv(1,2)/ngy ; hyy=cv(2,2)/ngy ; hyz=cv(3,2)/ngy
    hzx=cv(1,3)/ngz ; hzy=cv(2,3)/ngz ; hzz=cv(3,3)/ngz

    wa=f_malloc0([1-nagx.to.ngx+nagx,1-nagy.to.ngy+nagy,1-nagz.to.ngz+nagz],id='wa')
    wb=f_malloc0([1-nagx.to.ngx+nagx,1-nagy.to.ngy+nagy,1-nagz.to.ngz+nagz],id='wb')
    wc=f_malloc0([1-nagx.to.ngx+nagx,1-nagy.to.ngy+nagy,1-nagz.to.ngz+nagz],id='wc')
    allocate(exponentval(-nbgx:nbgx),expval(-nbgx:nbgx),rp4(-nbgx:nbgx))
   !  allocate(wa(1-nagx:ngx+nagx,1-nagy:ngy+nagy,1-nagz:ngz+nagz))

    hrxinv=real(ngx,8) !inverse of grid spacing in reduced coordinates
    hryinv=real(ngy,8) !inverse of grid spacing in reduced coordinates
    hrzinv=real(ngz,8) !inverse of grid spacing in reduced coordinates
    !if(trim(bc)=='bulk') then
    !    iii=0
    !elseif(trim(bc)=='slab') then
    !    iii=1
    !endif
    !-------------------------------------------------------
    pi=4.d0*atan(1.d0)
    sqpi=sqrt(pi)
    !-------------------------------------------------------
    do iat=1,nat
        gw_inv=1.d0/gw(iat)
        gwsq_inv=gw_inv*gw_inv
        gwcub_inv=gw_inv*gwsq_inv
        fac=4.d0/(15.d0*(gw(iat)**7)*(sqpi**3))
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
                    rp4(igx)=dmsq**2
                enddo
                call vdexp(2*iix+1,exponentval(-iix),expval(-iix))
                wa(imgx-iix:imgx+iix,jgy,jgz)=wa(imgx-iix:imgx+iix,jgy,jgz)+facqiat*rp4(-iix:iix)*expval(-iix:iix)
                wb(imgx-iix:imgx+iix,jgy,jgz)=wb(imgx-iix:imgx+iix,jgy,jgz)+fac*rp4(-iix:iix)*expval(-iix:iix)
                wc(imgx-iix:imgx+iix,jgy,jgz)=wc(imgx-iix:imgx+iix,jgy,jgz)+facqiat*rp4(-iix:iix)*expval(-iix:iix)*(-7.d0 - 2.d0*exponentval(-iix:iix))*gw_inv

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
    finalx = modulo((ngx+nagx-1),ngx)+1
    igxs = 1-nagx+(ngx-istartx+1)
    igxf = ngx+nagx-finalx+1
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
end subroutine put_rp4gto_sym
!*****************************************************************************************
