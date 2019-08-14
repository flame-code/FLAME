!*****************************************************************************************
subroutine put_gto_sym_ortho(parini,bc,reset,nat,rxyz,qat,gw,rgcut,ngx,ngy,ngz,hgrid,rho)
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
    !work array that is bigger than rho array, big enough to include of 
    !grid points that are outside of box.
    !local variables
    !work arrays to save the values of one dimensional gaussian function.
    real(8), allocatable:: wx(:), wy(:), wz(:)
    real(8):: rhoz, rhoyz, pi
    real(8):: hgxinv, hgyinv, hgzinv
    real(8):: width_inv, width_inv_xat, width_inv_yat, width_inv_zat
    real(8):: width_inv_hgx, width_inv_hgy, width_inv_hgz
    real(8):: xat, yat, zat, facqiat, fac, width
    real(8):: hx, hy, hz
    integer:: iat, iw, ix, iy, iz, iatox, iatoy, iatoz, jx, jy, jz
    integer:: nbgx, nbgy, nbgz, nagx, nagy, nagz
    integer:: ibcx, ibcy, ibcz
    real(8), allocatable:: wa(:,:,:)
    integer, allocatable:: mboundg(:,:,:)
    call f_routine(id='put_gto_sym_ortho')
    if(hgrid(2,1)/=0.d0 .or. hgrid(3,1)/=0.d0 .or. &
       hgrid(1,2)/=0.d0 .or. hgrid(3,2)/=0.d0 .or. &
       hgrid(1,3)/=0.d0 .or. hgrid(2,3)/=0.d0) then
        write(*,'(a)') 'ERROR: this routine is only for orthogonal cell:'
        write(*,'(3es14.5)') hgrid(1,1),hgrid(2,1),hgrid(3,1)
        write(*,'(3es14.5)') hgrid(1,2),hgrid(2,2),hgrid(3,2)
        write(*,'(3es14.5)') hgrid(1,3),hgrid(2,3),hgrid(3,3)
        stop
    endif
    hx=hgrid(1,1)
    hy=hgrid(2,2)
    hz=hgrid(3,3)
    nbgx=int(rgcut/hx)+2
    nbgy=int(rgcut/hy)+2
    nbgz=int(rgcut/hz)+2
    nagx=nbgx+1
    nagy=nbgy+1
    nagz=nbgz+1
    !-------------------------------------------------------
    ibcx=0 !zero means periodic
    ibcy=0 !zero means periodic
    ibcz=0 !zero means periodic
    if(trim(bc)=='bulk') then
        !variables already are set to correct value
    elseif(trim(bc)=='slab') then
        ibcz=1
        nagz=0
    elseif(trim(bc)=='wire') then
        ibcy=1
        ibcz=1
        nagy=0
        nagz=0
    elseif(trim(bc)=='free') then
        ibcx=1
        ibcy=1
        ibcz=1
        nagx=0
        nagy=0
        nagz=0
    else
        write(*,'(2a)') 'ERROR: unknown BC, bc= ',trim(bc)
        stop
    endif
    !-------------------------------------------------------
    mboundg=f_malloc([1.to.2,-nbgy.to.nbgy,-nbgz.to.nbgz],id='mboundg')
    call get_glimitsphere(hx,hy,hz,nbgx,nbgy,nbgz,mboundg)
    wa=f_malloc([1-nagx.to.ngx+nagx,1-nagy.to.ngy+nagy,1-nagz.to.ngz+nagz],id='wa')
    wx=f_malloc([-nbgx.to.nbgx],id='wx')
    wy=f_malloc([-nbgy.to.nbgy],id='wy')
    wz=f_malloc([-nbgz.to.nbgz],id='wz')
    pi=4.d0*atan(1.d0)
    hgxinv=1.d0/hx
    hgyinv=1.d0/hy
    hgzinv=1.d0/hz
    !1+nbg?*ibc? gives the index of grid point whose ? coordinate is zero (?=x,y,z)
    wa=0.d0
    do iat=1,nat
        !shift the gaussian centers
        iatox=nint(rxyz(1,iat)*hgxinv)+1+nbgx*ibcx
        iatoy=nint(rxyz(2,iat)*hgyinv)+1+nbgy*ibcy
        iatoz=nint(rxyz(3,iat)*hgzinv)+1+nbgz*ibcz
        xat=rxyz(1,iat)-(iatox-1-nbgx*ibcx)*hx
        yat=rxyz(2,iat)-(iatoy-1-nbgy*ibcy)*hy
        zat=rxyz(3,iat)-(iatoz-1-nbgz*ibcz)*hz
        !construct the one-dimensional gaussians

        width=gw(iat)
        width_inv=1.d0/width
        fac=1.d0/(width*sqrt(pi))**3
        width_inv_hgx=width_inv*hx
        width_inv_hgy=width_inv*hy
        width_inv_hgz=width_inv*hz

        width_inv_xat=width_inv*xat
        width_inv_yat=width_inv*yat
        width_inv_zat=width_inv*zat
        do iw=-nbgx,nbgx
            wx(iw)=exp(-(width_inv_hgx*iw-width_inv_xat)**2)
        enddo
        do iw=-nbgy,nbgy
            wy(iw)=exp(-(width_inv_hgy*iw-width_inv_yat)**2)
        enddo
        do iw=-nbgz,nbgz
            wz(iw)=exp(-(width_inv_hgz*iw-width_inv_zat)**2)
        enddo
        facqiat=fac*qat(iat)
        do iz=-nbgz,nbgz
            rhoz=facqiat*wz(iz)
            jz=iatoz+iz
            do iy=-nbgy,nbgy
                rhoyz=rhoz*wy(iy)
                jy=iatoy+iy
                do ix=mboundg(1,iy,iz),mboundg(2,iy,iz)
                    jx=iatox+ix
                    !write(*,'(5i5)') iat,iatox,ix,iatoy,iy
                    wa(jx,jy,jz)=wa(jx,jy,jz)+rhoyz*wx(ix)
                enddo
            enddo
        enddo
    enddo
    !---------------------------------------------------------------------------
    if(reset) then
        !if the input array of charge density does not contain any previous value
        !wanted to be preserved.
        rho = 0.d0
    endif
    call charge_back_to_cell(ngx,ngy,ngz,nagx,nagy,nagz,ibcx,wa,rho)
    !do iz=1,ngz
    !    do iy=1,ngy
    !        do ix=1,ngx
    !            write(61,'(3i4,es20.10)') ix,iy,iz,rho(ix,iy,iz)
    !        enddo
    !    enddo
    !enddo
    !stop
    call f_free(wx)
    call f_free(wy)
    call f_free(wz)
    call f_free(wa)
    call f_free(mboundg)
    call f_release_routine()
end subroutine put_gto_sym_ortho
!*****************************************************************************************
subroutine qgrad_gto_sym_ortho(parini,bc,nat,rxyz,qat,gw,rgcut,lda,ngx,ngy,ngz,hgrid,pot,g)
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
    real(8), intent(inout):: g(nat)
    !local variables
    real(8), allocatable:: wx(:), wy(:), wz(:) !values of one dimensional Gaussian functions
    real(8):: rhoz, rhoyz, pi
    real(8):: hgxinv, hgyinv, hgzinv, hgxhgyhgz
    real(8):: width_inv, width_inv_xat, width_inv_yat, width_inv_zat
    real(8):: width_inv_hgx, width_inv_hgy, width_inv_hgz,width
    real(8):: xat, yat, zat, fac
    integer:: iat, iw, ix, iy, iz, iatox, iatoy, iatoz, jx, jy, jz, iyt, izt
    real(8), allocatable:: wa(:,:,:)
    integer, allocatable:: mboundg(:,:,:)
    integer:: nbgx, nbgy, nbgz, nagx, nagy, nagz
    integer:: ibcx, ibcy, ibcz
    real(8):: hx, hy, hz
    real(8):: ttg
    if(hgrid(2,1)/=0.d0 .or. hgrid(3,1)/=0.d0 .or. &
       hgrid(1,2)/=0.d0 .or. hgrid(3,2)/=0.d0 .or. &
       hgrid(1,3)/=0.d0 .or. hgrid(2,3)/=0.d0) then
        write(*,'(a)') 'ERROR: this routine is only for orthogonal cell:'
        write(*,'(3es14.5)') hgrid(1,1),hgrid(2,1),hgrid(3,1)
        write(*,'(3es14.5)') hgrid(1,2),hgrid(2,2),hgrid(3,2)
        write(*,'(3es14.5)') hgrid(1,3),hgrid(2,3),hgrid(3,3)
        stop
    endif
    hx=hgrid(1,1)
    hy=hgrid(2,2)
    hz=hgrid(3,3)
    nbgx=int(rgcut/hx)+2
    nbgy=int(rgcut/hy)+2
    nbgz=int(rgcut/hz)+2
    nagx=nbgx+1
    nagy=nbgy+1
    nagz=nbgz+1
    !-------------------------------------------------------
    ibcx=0 !zero means periodic
    ibcy=0 !zero means periodic
    ibcz=0 !zero means periodic
    if(trim(bc)=='bulk') then
        !variables already are set to correct value
    elseif(trim(bc)=='slab') then
        ibcz=1
        nagz=0
    elseif(trim(bc)=='wire') then
        ibcy=1
        ibcz=1
        nagy=0
        nagz=0
    elseif(trim(bc)=='free') then
        ibcx=1
        ibcy=1
        ibcz=1
        nagx=0
        nagy=0
        nagz=0
    else
        write(*,'(2a)') 'ERROR: unknown BC, bc= ',trim(bc)
        stop
    endif
    !-------------------------------------------------------
    mboundg=f_malloc([1.to.2,-nbgy.to.nbgy,-nbgz.to.nbgz],id='mboundg')
    call get_glimitsphere(hx,hy,hz,nbgx,nbgy,nbgz,mboundg)
    wa=f_malloc([1-nagx.to.ngx+nagx,1-nagy.to.ngy+nagy,1-nagz.to.ngz+nagz],id='wa')
    wx=f_malloc([-nbgx.to.nbgx],id='wx')
    wy=f_malloc([-nbgy.to.nbgy],id='wy')
    wz=f_malloc([-nbgz.to.nbgz],id='wz')
    pi=4.d0*atan(1.d0)
    hgxhgyhgz=hx*hy*hz
    hgxinv=1.d0/hx
    hgyinv=1.d0/hy
    hgzinv=1.d0/hz
    !---------------------------------------------------------------------------
    call potential_on_extended_grid(lda,ngx,ngy,ngz,nagx,nagy,nagz,ibcx,pot,wa)
    !initialize the density 
    do iat=1,nat  
        !shift the Gaussian centers
        iatox=nint(rxyz(1,iat)*hgxinv)+1+nbgx*ibcx
        iatoy=nint(rxyz(2,iat)*hgyinv)+1+nbgy*ibcy
        iatoz=nint(rxyz(3,iat)*hgzinv)+1+nbgz*ibcz
        xat=rxyz(1,iat)-(iatox-1-nbgx*ibcx)*hx
        yat=rxyz(2,iat)-(iatoy-1-nbgy*ibcy)*hy
        zat=rxyz(3,iat)-(iatoz-1-nbgz*ibcz)*hz
        !construct the one-dimensional gaussians

        width=gw(iat)
        width_inv=1.d0/width
        fac=1.d0/(width*sqrt(pi))**3
        width_inv_hgx=width_inv*hx
        width_inv_hgy=width_inv*hy
        width_inv_hgz=width_inv*hz

        width_inv_xat=width_inv*xat
        width_inv_yat=width_inv*yat
        width_inv_zat=width_inv*zat
        do iw=-nbgx,nbgx
            wx(iw)=exp(-(width_inv_hgx*iw-width_inv_xat)**2)
        enddo
        do iw=-nbgy,nbgy
            wy(iw)=exp(-(width_inv_hgy*iw-width_inv_yat)**2)
        enddo
        do iw=-nbgz,nbgz
            wz(iw)=exp(-(width_inv_hgz*iw-width_inv_zat)**2)
        enddo

        ttg=0.d0
        do iz=-nbgz,nbgz
        jz=iatoz+iz
        rhoz=fac*wz(iz)
            do iy=-nbgy,nbgy
                rhoyz=rhoz*wy(iy)
                jy=iatoy+iy
                do ix=mboundg(1,iy,iz),mboundg(2,iy,iz)
                    jx=iatox+ix
                    ttg=ttg+rhoyz*wx(ix)*wa(jx,jy,jz)
                enddo
            enddo
        enddo
        g(iat)=g(iat)+ttg*hgxhgyhgz
    enddo
    call f_free(wx)
    call f_free(wy)
    call f_free(wz)
    call f_free(wa)
    call f_free(mboundg)
    ! call f_release_routine()
end subroutine qgrad_gto_sym_ortho
!*****************************************************************************************
subroutine force_gto_sym_ortho(parini,bc,nat,rxyz,qat,gw,rgcut,lda,ngx,ngy,ngz,hgrid,pot,fat)
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
    real(8), intent(inout):: pot(lda,ngy,ngz)
    real(8), intent(out):: fat(3,nat)
    !local variables
    real(8), allocatable:: wx(:), wy(:), wz(:) !values of one dimensional Gaussian functions
    real(8), allocatable:: vx(:), vy(:), vz(:) !derivatives of wx,wy,wz arrays
    real(8):: rhoz, rhoyz, pi
    real(8):: hgxinv, hgyinv, hgzinv, hgxhgyhgz
    real(8):: width, width_inv, width_inv_xat, width_inv_yat, width_inv_zat
    real(8):: width_inv_hgx, width_inv_hgy, width_inv_hgz
    real(8):: xat, yat, zat, facqiat, fac
    integer:: iat, ivw, ix, iy, iz, iatox, iatoy, iatoz, jx, jy, jz, iyt, izt
    integer:: ibcx, ibcy, ibcz
    real(8):: fx, fy, fz, dx, dy, dz, derrhoz, derrhoyz, derrhozy
    integer:: nbgx, nbgy, nbgz, nagx, nagy, nagz
    real(8):: hx, hy, hz
    real(8), allocatable:: wa(:,:,:)
    integer, allocatable:: mboundg(:,:,:)
    call f_routine(id='force_gto_sym_ortho')
    if(hgrid(2,1)/=0.d0 .or. hgrid(3,1)/=0.d0 .or. &
       hgrid(1,2)/=0.d0 .or. hgrid(3,2)/=0.d0 .or. &
       hgrid(1,3)/=0.d0 .or. hgrid(2,3)/=0.d0) then
        write(*,'(a)') 'ERROR: this routine is only for orthogonal cell:'
        write(*,'(3es14.5)') hgrid(1,1),hgrid(2,1),hgrid(3,1)
        write(*,'(3es14.5)') hgrid(1,2),hgrid(2,2),hgrid(3,2)
        write(*,'(3es14.5)') hgrid(1,3),hgrid(2,3),hgrid(3,3)
        stop
    endif
    hx=hgrid(1,1)
    hy=hgrid(2,2)
    hz=hgrid(3,3)
    nbgx=int(rgcut/hx)+2
    nbgy=int(rgcut/hy)+2
    nbgz=int(rgcut/hz)+2
    nagx=nbgx+1
    nagy=nbgy+1
    nagz=nbgz+1
    !-------------------------------------------------------
    ibcx=0 !zero means periodic
    ibcy=0 !zero means periodic
    ibcz=0 !zero means periodic
    if(trim(bc)=='bulk') then
        !variables already are set to correct value
    elseif(trim(bc)=='slab') then
        ibcz=1
        nagz=0
    elseif(trim(bc)=='wire') then
        ibcy=1
        ibcz=1
        nagy=0
        nagz=0
    elseif(trim(bc)=='free') then
        ibcx=1
        ibcy=1
        ibcz=1
        nagx=0
        nagy=0
        nagz=0
    else
        write(*,'(2a)') 'ERROR: unknown BC, bc= ',trim(bc)
        stop
    endif
    !-------------------------------------------------------
    mboundg=f_malloc([1.to.2,-nbgy.to.nbgy,-nbgz.to.nbgz],id='mboundg')
    call get_glimitsphere(hx,hy,hz,nbgx,nbgy,nbgz,mboundg)
    wa=f_malloc([1-nagx.to.ngx+nagx,1-nagy.to.ngy+nagy,1-nagz.to.ngz+nagz],id='wa')
    wx=f_malloc([-nbgx.to.nbgx],id='wx')
    wy=f_malloc([-nbgy.to.nbgy],id='wy')
    wz=f_malloc([-nbgz.to.nbgz],id='wz')
    vx=f_malloc([-nbgx.to.nbgx],id='vx')
    vy=f_malloc([-nbgy.to.nbgy],id='vy')
    vz=f_malloc([-nbgz.to.nbgz],id='vz')
    pi=4.d0*atan(1.d0)
    hgxhgyhgz=hx*hy*hz
    hgxinv=1.d0/hx
    hgyinv=1.d0/hy
    hgzinv=1.d0/hz
    !---------------------------------------------------------------------------
    call potential_on_extended_grid(lda,ngx,ngy,ngz,nagx,nagy,nagz,ibcx,pot,wa)
    !initialize the density 
    do iat=1,nat  
        !shift the Gaussian centers
        iatox=nint(rxyz(1,iat)*hgxinv)+1+nbgx*ibcx
        iatoy=nint(rxyz(2,iat)*hgyinv)+1+nbgy*ibcy
        iatoz=nint(rxyz(3,iat)*hgzinv)+1+nbgz*ibcz
        xat=rxyz(1,iat)-(iatox-1-nbgx*ibcx)*hx
        yat=rxyz(2,iat)-(iatoy-1-nbgy*ibcy)*hy
        zat=rxyz(3,iat)-(iatoz-1-nbgz*ibcz)*hz
        width=gw(iat)
        width_inv=1.d0/width
        fac=2.d0/(width*(width*sqrt(pi))**3)
        width_inv_hgx=width_inv*hx
        width_inv_hgy=width_inv*hy
        width_inv_hgz=width_inv*hz
        width_inv_xat=width_inv*xat
        width_inv_yat=width_inv*yat
        width_inv_zat=width_inv*zat
        !construct the one-dimensional gaussians
        do ivw=-nbgx,nbgx
            dx=width_inv_hgx*ivw-width_inv_xat
            wx(ivw)=exp(-dx*dx)
            vx(ivw)=-wx(ivw)*dx
        enddo
        do ivw=-nbgy,nbgy
            dy=width_inv_hgy*ivw-width_inv_yat
            wy(ivw)=exp(-dy*dy)
            vy(ivw)=-wy(ivw)*dy
        enddo
        do ivw=-nbgz,nbgz
            dz=width_inv_hgz*ivw-width_inv_zat
            wz(ivw)=exp(-dz*dz)
            vz(ivw)=-wz(ivw)*dz
        enddo
        fx=0.d0;fy=0.d0;fz=0.d0
        facqiat=fac*qat(iat)
        do iz=-nbgz,nbgz
            rhoz=facqiat*wz(iz)
            derrhoz=facqiat*vz(iz)
            jz=iatoz+iz
            do iy=-nbgy,nbgy
                rhoyz=rhoz*wy(iy)
                derrhozy=rhoz*vy(iy)
                derrhoyz=derrhoz*wy(iy)
                jy=iatoy+iy
                do ix=mboundg(1,iy,iz),mboundg(2,iy,iz)
                    jx=iatox+ix
                    fx=fx+rhoyz*vx(ix)*wa(jx,jy,jz)
                    fy=fy+derrhozy*wx(ix)*wa(jx,jy,jz)
                    fz=fz+derrhoyz*wx(ix)*wa(jx,jy,jz)
                enddo
            enddo
        enddo
        fat(1,iat)=fat(1,iat)+fx*hgxhgyhgz
        fat(2,iat)=fat(2,iat)+fy*hgxhgyhgz
        fat(3,iat)=fat(3,iat)+fz*hgxhgyhgz
    enddo
    call f_free(wa)
    call f_free(wx)
    call f_free(wy)
    call f_free(wz)
    call f_free(vx)
    call f_free(vy)
    call f_free(vz)
    call f_free(mboundg)
    call f_release_routine()
end subroutine force_gto_sym_ortho
!*****************************************************************************************
