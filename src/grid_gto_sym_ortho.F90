!*****************************************************************************************
subroutine put_gto_sym_ortho(parini,bc,reset,nat,rxyz,cv,qat,gw,rgcut,ngx,ngy,ngz,hgrid,rho)
    use mod_interface
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
    hx=hgrid(1,1)
    hy=hgrid(2,2)
    hz=hgrid(3,3)
    nbgx=int(rgcut/hx)+2
    nbgy=int(rgcut/hy)+2
    nbgz=int(rgcut/hz)+2
    nagx=nbgx+1
    nagy=nbgy+1
    nagz=0
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
!    write(*,*) hgxinv,hgyinv,hgzinv
    !assaigning the gaussian charge density on grid points, saving in the 
    !work array which is bigger than rho array.
    !write(*,*) '********************************** ',associated(wa)
    !if(trim(bc)=='bulk') then
    !    iii=0
    !elseif(trim(bc)=='slab') then
    !    iii=1
    !endif
    ibcz=1
    wa=0.d0
    do iat=1,nat
        !shift the gaussian centers
        iatox=nint(rxyz(1,iat)*hgxinv)+1
        iatoy=nint(rxyz(2,iat)*hgyinv)+1
        iatoz=nint(rxyz(3,iat)*hgzinv)+1+nbgz*ibcz
        xat=rxyz(1,iat)-(iatox-1)*hx
        yat=rxyz(2,iat)-(iatoy-1)*hy
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
    do iz=1-nagz,ngz+nagz
        do iy=1-nagy,0
            do ix=1-nagx,0
                wa(ix+ngx,iy+ngy,iz)=wa(ix+ngx,iy+ngy,iz)+wa(ix,iy,iz)
            enddo
            do ix=1,ngx
                wa(ix,iy+ngy,iz)=wa(ix,iy+ngy,iz)+wa(ix,iy,iz)
            enddo
            do ix=ngx+1,ngx+nagx
                wa(ix-ngx,iy+ngy,iz)=wa(ix-ngx,iy+ngy,iz)+wa(ix,iy,iz)
            enddo
        enddo
        do iy=1,ngy
            do ix=1-nagx,0
                wa(ix+ngx,iy,iz)=wa(ix+ngx,iy,iz)+wa(ix,iy,iz)
            enddo
            do ix=ngx+1,ngx+nagx
                wa(ix-ngx,iy,iz)=wa(ix-ngx,iy,iz)+wa(ix,iy,iz)
            enddo
        enddo
        do iy=ngy+1,ngy+nagy
            do ix=1-nagx,0
                wa(ix+ngx,iy-ngy,iz)=wa(ix+ngx,iy-ngy,iz)+wa(ix,iy,iz)
            enddo
            do ix=1,ngx
                wa(ix,iy-ngy,iz)=wa(ix,iy-ngy,iz)+wa(ix,iy,iz)
            enddo
            do ix=ngx+1,ngx+nagx
                wa(ix-ngx,iy-ngy,iz)=wa(ix-ngx,iy-ngy,iz)+wa(ix,iy,iz)
            enddo
        enddo
    enddo
    do iz=1-nagz,0
        do iy=1,ngy
            do ix=1,ngx
                wa(ix,iy,iz+ngz)=wa(ix,iy,iz+ngz)+wa(ix,iy,iz)
            enddo
        enddo
    enddo
    do iz=ngz+1,ngz+nagz
        do iy=1,ngy
            do ix=1,ngx
                wa(ix,iy,iz-ngz)=wa(ix,iy,iz-ngz)+wa(ix,iy,iz)
            enddo
        enddo
    enddo
    if(reset) then
        do iz=1,ngz
            do iy=1,ngy
                do ix=1,ngx
                    rho(ix,iy,iz)=wa(ix,iy,iz)
                enddo
            enddo
        enddo
    else
        do iz=1,ngz
            do iy=1,ngy
                do ix=1,ngx
                    rho(ix,iy,iz)=rho(ix,iy,iz)+wa(ix,iy,iz)
                enddo
            enddo
        enddo
    endif
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
subroutine force_gto_sym_ortho(parini,atoms,gw,rgcut,lda,ngx,ngy,ngz,hgrid,pot)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(in):: gw(atoms%nat)
    real(8), intent(in):: rgcut
    integer, intent(in):: lda, ngx, ngy, ngz
    real(8), intent(in):: hgrid(3,3)
    real(8), intent(inout):: pot(lda,ngy,ngz)
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
    hx=hgrid(1,1)
    hy=hgrid(2,2)
    hz=hgrid(3,3)
    nbgx=int(rgcut/hx)+2
    nbgy=int(rgcut/hy)+2
    nbgz=int(rgcut/hz)+2
    nagx=nbgx+1
    nagy=nbgy+1
    nagz=0
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
    do iz=1-nagz,ngz+nagz
        izt=iz+(sign(ngz,-iz)+sign(ngz,ngz-iz))/2
        do iy=1-nagy,ngy+nagy
            iyt=iy+(sign(ngy,-iy)+sign(ngy,ngy-iy))/2
            do ix=1-nagx,0
                wa(ix,iy,iz)=pot(ix+ngx,iyt,izt)
            enddo
            do ix=1,ngx
                wa(ix,iy,iz)=pot(ix,iyt,izt)
            enddo
            do ix=ngx+1,ngx+nagx
                wa(ix,iy,iz)=pot(ix-ngx,iyt,izt)
            enddo
        enddo
    enddo
    !---------------------------------------------------------------------------
    !if(trim(atoms%boundcond)=='bulk') then
    !    iii=0
    !elseif(trim(atoms%boundcond)=='slab') then
    !    iii=1
    !endif
    ibcz=1
    !initialize the density 
    do iat=1,atoms%nat  
        !shift the Gaussian centers
        iatox=nint(atoms%rat(1,iat)*hgxinv)+1
        iatoy=nint(atoms%rat(2,iat)*hgyinv)+1
        iatoz=nint(atoms%rat(3,iat)*hgzinv)+1+nbgz*ibcz
        xat=atoms%rat(1,iat)-(iatox-1)*hx
        yat=atoms%rat(2,iat)-(iatoy-1)*hy
        zat=atoms%rat(3,iat)-(iatoz-1-nbgz*ibcz)*hz
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
        facqiat=fac*atoms%qat(iat)
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
        atoms%fat(1,iat)=atoms%fat(1,iat)+fx*hgxhgyhgz
        atoms%fat(2,iat)=atoms%fat(2,iat)+fy*hgxhgyhgz
        atoms%fat(3,iat)=atoms%fat(3,iat)+fz*hgxhgyhgz
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
subroutine qgrad_gto_sym_ortho(parini,atoms,gw,rgcut,lda,ngx,ngy,ngz,hgrid,pot,g)
    use mod_interface
    use mod_atoms, only: typ_atoms
    use mod_parini, only: typ_parini
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    real(8), intent(in):: gw(atoms%nat)
    real(8), intent(in):: rgcut
    integer, intent(in):: lda, ngx, ngy, ngz
    real(8), intent(in):: hgrid(3,3)
    real(8), intent(in):: pot(lda,ngy,ngz)
    real(8), intent(inout):: g(atoms%nat)
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
    hx=hgrid(1,1)
    hy=hgrid(2,2)
    hz=hgrid(3,3)
    nbgx=int(rgcut/hx)+2
    nbgy=int(rgcut/hy)+2
    nbgz=int(rgcut/hz)+2
    nagx=nbgx+1
    nagy=nbgy+1
    nagz=0
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
    do iz=1-nagz,ngz+nagz
        izt=iz+(sign(ngz,-iz)+sign(ngz,ngz-iz))/2
        do iy=1-nagy,ngy+nagy
            iyt=iy+(sign(ngy,-iy)+sign(ngy,ngy-iy))/2
            do ix=1-nagx,0
                wa(ix,iy,iz)=pot(ix+ngx,iyt,izt)
            enddo
            do ix=1,ngx
                wa(ix,iy,iz)=pot(ix,iyt,izt)
            enddo
            do ix=ngx+1,ngx+nagx
                wa(ix,iy,iz)=pot(ix-ngx,iyt,izt)
            enddo
        enddo
    enddo
    !---------------------------------------------------------------------------
    !    if(trim(atoms%boundcond)=='bulk') then
    !        iii=0
    !    elseif(trim(atoms%boundcond)=='slab') then
    !        iii=1
    !    endif
    ibcz=1
    !initialize the density 
    do iat=1,atoms%nat  
        !shift the Gaussian centers
        iatox=nint(atoms%rat(1,iat)*hgxinv)+1
        iatoy=nint(atoms%rat(2,iat)*hgyinv)+1
        iatoz=nint(atoms%rat(3,iat)*hgzinv)+1+nbgz*ibcz
        xat=atoms%rat(1,iat)-(iatox-1)*hx
        yat=atoms%rat(2,iat)-(iatoy-1)*hy
        zat=atoms%rat(3,iat)-(iatoz-1-nbgz*ibcz)*hz
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
