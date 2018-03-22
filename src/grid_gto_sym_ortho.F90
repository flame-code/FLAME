!*****************************************************************************************
subroutine put_gto_sym_ortho(parini,bc,reset,nat,rxyz,qat,gausswidth,poisson)
    use mod_interface
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    use mod_parini, only: typ_parini
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: bc
    logical, intent(in):: reset
    integer, intent(in):: nat
    real(8), intent(in):: rxyz(3,nat)
    real(8), intent(in):: qat(nat)
    real(8), intent(in):: gausswidth(nat)
    type(typ_poisson), intent(inout):: poisson
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
    integer:: iat, iw, ix, iy, iz, iatox, iatoy, iatoz, jx, jy, jz
    integer:: nbgpx, nbgpy, nbgpz, nagpx, nagpy, nagpz
    integer:: ngpx, ngpy, ngpz, iii
    real(8), allocatable:: wa(:,:,:)
    integer, allocatable:: mboundg(:,:,:)
    call f_routine(id='put_gto_sym_ortho')
    nbgpx=int(poisson%rgcut/poisson%hx)+2
    nbgpy=int(poisson%rgcut/poisson%hy)+2
    nbgpz=int(poisson%rgcut/poisson%hz)+2
    nagpx=nbgpx+1
    nagpy=nbgpy+1
    nagpz=0
    mboundg=f_malloc([1.to.2,-nbgpy.to.nbgpy,-nbgpz.to.nbgpz],id='mboundg')
    call get_glimitsphere(poisson,nbgpx,nbgpy,nbgpz,mboundg)
    ngpx=poisson%ngpx
    ngpy=poisson%ngpy
    ngpz=poisson%ngpz
    !wa(1-nagpx:ngpx+nagpx,1-nagpy:ngpy+nagpy,1-nagpz:ngpz+nagpz)=>poisson%pot
    wa=f_malloc([1-nagpx.to.ngpx+nagpx,1-nagpy.to.ngpy+nagpy,1-nagpz.to.ngpz+nagpz],id='wa')
    wx=f_malloc([-nbgpx.to.nbgpx],id='wx')
    wy=f_malloc([-nbgpy.to.nbgpy],id='wy')
    wz=f_malloc([-nbgpz.to.nbgpz],id='wz')
    pi=4.d0*atan(1.d0)
    hgxinv=1.d0/poisson%hx
    hgyinv=1.d0/poisson%hy
    hgzinv=1.d0/poisson%hz
!    write(*,*) hgxinv,hgyinv,hgzinv
    !assaigning the gaussian charge density on grid points, saving in the 
    !work array which is bigger than rho array.
    !write(*,*) '********************************** ',associated(wa)
    if(trim(bc)=='bulk') then
        iii=0
    elseif(trim(bc)=='slab') then
        iii=1
    endif
    wa=0.d0
    if (parini%ewald) then
        width= poisson%alpha
        width_inv=1.d0/width
        fac=1.d0/(width*sqrt(pi))**3
        width_inv_hgx=width_inv*poisson%hx
        width_inv_hgy=width_inv*poisson%hy
        width_inv_hgz=width_inv*poisson%hz
    endif
    do iat=1,nat
        !shift the gaussian centers
        iatox=nint(rxyz(1,iat)*hgxinv)+1
        iatoy=nint(rxyz(2,iat)*hgyinv)+1
        iatoz=nint(rxyz(3,iat)*hgzinv)+1+nbgpz*iii
        xat=rxyz(1,iat)-(iatox-1)*poisson%hx
        yat=rxyz(2,iat)-(iatoy-1)*poisson%hy
        zat=rxyz(3,iat)-(iatoz-1-nbgpz*iii)*poisson%hz
        !construct the one-dimensional gaussians

        if (.not.parini%ewald) then
            width=gausswidth(iat)
            width_inv=1.d0/width
            fac=1.d0/(width*sqrt(pi))**3
            width_inv_hgx=width_inv*poisson%hx
            width_inv_hgy=width_inv*poisson%hy
            width_inv_hgz=width_inv*poisson%hz
        endif

        width_inv_xat=width_inv*xat
        width_inv_yat=width_inv*yat
        width_inv_zat=width_inv*zat
        do iw=-nbgpx,nbgpx
            wx(iw)=exp(-(width_inv_hgx*iw-width_inv_xat)**2)
        enddo
        do iw=-nbgpy,nbgpy
            wy(iw)=exp(-(width_inv_hgy*iw-width_inv_yat)**2)
        enddo
        do iw=-nbgpz,nbgpz
            wz(iw)=exp(-(width_inv_hgz*iw-width_inv_zat)**2)
        enddo
        facqiat=fac*qat(iat)
        do iz=-nbgpz,nbgpz
            rhoz=facqiat*wz(iz)
            jz=iatoz+iz
            do iy=-nbgpy,nbgpy
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
    do iz=1-nagpz,ngpz+nagpz
        do iy=1-nagpy,0
            do ix=1-nagpx,0
                wa(ix+ngpx,iy+ngpy,iz)=wa(ix+ngpx,iy+ngpy,iz)+wa(ix,iy,iz)
            enddo
            do ix=1,ngpx
                wa(ix,iy+ngpy,iz)=wa(ix,iy+ngpy,iz)+wa(ix,iy,iz)
            enddo
            do ix=ngpx+1,ngpx+nagpx
                wa(ix-ngpx,iy+ngpy,iz)=wa(ix-ngpx,iy+ngpy,iz)+wa(ix,iy,iz)
            enddo
        enddo
        do iy=1,ngpy
            do ix=1-nagpx,0
                wa(ix+ngpx,iy,iz)=wa(ix+ngpx,iy,iz)+wa(ix,iy,iz)
            enddo
            do ix=ngpx+1,ngpx+nagpx
                wa(ix-ngpx,iy,iz)=wa(ix-ngpx,iy,iz)+wa(ix,iy,iz)
            enddo
        enddo
        do iy=ngpy+1,ngpy+nagpy
            do ix=1-nagpx,0
                wa(ix+ngpx,iy-ngpy,iz)=wa(ix+ngpx,iy-ngpy,iz)+wa(ix,iy,iz)
            enddo
            do ix=1,ngpx
                wa(ix,iy-ngpy,iz)=wa(ix,iy-ngpy,iz)+wa(ix,iy,iz)
            enddo
            do ix=ngpx+1,ngpx+nagpx
                wa(ix-ngpx,iy-ngpy,iz)=wa(ix-ngpx,iy-ngpy,iz)+wa(ix,iy,iz)
            enddo
        enddo
    enddo
    do iz=1-nagpz,0
        do iy=1,ngpy
            do ix=1,ngpx
                wa(ix,iy,iz+ngpz)=wa(ix,iy,iz+ngpz)+wa(ix,iy,iz)
            enddo
        enddo
    enddo
    do iz=ngpz+1,ngpz+nagpz
        do iy=1,ngpy
            do ix=1,ngpx
                wa(ix,iy,iz-ngpz)=wa(ix,iy,iz-ngpz)+wa(ix,iy,iz)
            enddo
        enddo
    enddo
    if(reset) then
        do iz=1,ngpz
            do iy=1,ngpy
                do ix=1,ngpx
                    poisson%rho(ix,iy,iz)=wa(ix,iy,iz)
                enddo
            enddo
        enddo
    else
        do iz=1,ngpz
            do iy=1,ngpy
                do ix=1,ngpx
                    poisson%rho(ix,iy,iz)=poisson%rho(ix,iy,iz)+wa(ix,iy,iz)
                enddo
            enddo
        enddo
    endif
    !do iz=1,ngpz
    !    do iy=1,ngpy
    !        do ix=1,ngpx
    !            write(61,'(3i4,es20.10)') ix,iy,iz,poisson%rho(ix,iy,iz)
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
subroutine force_gto_sym_ortho(parini,atoms,poisson,gausswidth)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_poisson), intent(inout):: poisson
    real(8), intent(in):: gausswidth(atoms%nat)
    !local variables
    real(8), allocatable:: wx(:), wy(:), wz(:) !values of one dimensional Gaussian functions
    real(8), allocatable:: vx(:), vy(:), vz(:) !derivatives of wx,wy,wz arrays
    real(8):: rhoz, rhoyz, pi
    real(8):: hgxinv, hgyinv, hgzinv, hgxhgyhgz
    real(8):: width, width_inv, width_inv_xat, width_inv_yat, width_inv_zat
    real(8):: width_inv_hgx, width_inv_hgy, width_inv_hgz
    real(8):: xat, yat, zat, facqiat, fac
    integer:: iat, ivw, ix, iy, iz, iatox, iatoy, iatoz, jx, jy, jz, iyt, izt, iii
    real(8):: fx, fy, fz, dx, dy, dz, derrhoz, derrhoyz, derrhozy
    integer:: nbgpx, nbgpy, nbgpz, nagpx, nagpy, nagpz
    real(8), allocatable:: wa(:,:,:)
    integer, allocatable:: mboundg(:,:,:)
    real(8):: dipole !,ext_pot
    real(8):: dv, beta                                           
    call f_routine(id='force_gto_sym_ortho')
    associate(ngpx=>poisson%ngpx)
    associate(ngpy=>poisson%ngpy)
    associate(ngpz=>poisson%ngpz)
    nbgpx=int(poisson%rgcut/poisson%hx)+2
    nbgpy=int(poisson%rgcut/poisson%hy)+2
    nbgpz=int(poisson%rgcut/poisson%hz)+2
    nagpx=nbgpx+1
    nagpy=nbgpy+1
    nagpz=0
    mboundg=f_malloc([1.to.2,-nbgpy.to.nbgpy,-nbgpz.to.nbgpz],id='mboundg')
    call get_glimitsphere(poisson,nbgpx,nbgpy,nbgpz,mboundg)
    !wa(1-nagpx:ngpx+nagpx,1-nagpy:ngpy+nagpy,1-nagpz:ngpz+nagpz)=>poisson%pot
    wa=f_malloc([1-nagpx.to.ngpx+nagpx,1-nagpy.to.ngpy+nagpy,1-nagpz.to.ngpz+nagpz],id='wa')
    wx=f_malloc([-nbgpx.to.nbgpx],id='wx')
    wy=f_malloc([-nbgpy.to.nbgpy],id='wy')
    wz=f_malloc([-nbgpz.to.nbgpz],id='wz')
    vx=f_malloc([-nbgpx.to.nbgpx],id='vx')
    vy=f_malloc([-nbgpy.to.nbgpy],id='vy')
    vz=f_malloc([-nbgpz.to.nbgpz],id='vz')
    pi=4.d0*atan(1.d0)
    hgxhgyhgz=poisson%hx*poisson%hy*poisson%hz
    hgxinv=1.d0/poisson%hx
    hgyinv=1.d0/poisson%hy
    hgzinv=1.d0/poisson%hz
    !---------------------------------------------------------------------------
    do iz=1-nagpz,ngpz+nagpz
        izt=iz+(sign(ngpz,-iz)+sign(ngpz,ngpz-iz))/2
        do iy=1-nagpy,ngpy+nagpy
            iyt=iy+(sign(ngpy,-iy)+sign(ngpy,ngpy-iy))/2
            do ix=1-nagpx,0
                wa(ix,iy,iz)=poisson%pot(ix+ngpx,iyt,izt)
            enddo
            do ix=1,ngpx
                wa(ix,iy,iz)=poisson%pot(ix,iyt,izt)
            enddo
            do ix=ngpx+1,ngpx+nagpx
                wa(ix,iy,iz)=poisson%pot(ix-ngpx,iyt,izt)
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
    do iat=1,atoms%nat  
        !shift the Gaussian centers
        iatox=nint(atoms%rat(1,iat)*hgxinv)+1
        iatoy=nint(atoms%rat(2,iat)*hgyinv)+1
        iatoz=nint(atoms%rat(3,iat)*hgzinv)+1+nbgpz*iii
        xat=atoms%rat(1,iat)-(iatox-1)*poisson%hx
        yat=atoms%rat(2,iat)-(iatoy-1)*poisson%hy
        zat=atoms%rat(3,iat)-(iatoz-1-nbgpz*iii)*poisson%hz
        width=gausswidth(iat)
        width_inv=1.d0/width
        fac=2.d0/(width*(width*sqrt(pi))**3)
        width_inv_hgx=width_inv*poisson%hx
        width_inv_hgy=width_inv*poisson%hy
        width_inv_hgz=width_inv*poisson%hz
        width_inv_xat=width_inv*xat
        width_inv_yat=width_inv*yat
        width_inv_zat=width_inv*zat
        !construct the one-dimensional gaussians
        do ivw=-nbgpx,nbgpx
            dx=width_inv_hgx*ivw-width_inv_xat
            wx(ivw)=exp(-dx*dx)
            vx(ivw)=-wx(ivw)*dx
        enddo
        do ivw=-nbgpy,nbgpy
            dy=width_inv_hgy*ivw-width_inv_yat
            wy(ivw)=exp(-dy*dy)
            vy(ivw)=-wy(ivw)*dy
        enddo
        do ivw=-nbgpz,nbgpz
            dz=width_inv_hgz*ivw-width_inv_zat
            wz(ivw)=exp(-dz*dz)
            vz(ivw)=-wz(ivw)*dz
        enddo
        fx=0.d0;fy=0.d0;fz=0.d0
        facqiat=fac*atoms%qat(iat)
        do iz=-nbgpz,nbgpz
            rhoz=facqiat*wz(iz)
            derrhoz=facqiat*vz(iz)
            jz=iatoz+iz
            do iy=-nbgpy,nbgpy
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
    if((.not. poisson%point_particle) .and. trim(parini%bias_type)=='fixed_efield') then
        do iat=1,atoms%nat
            atoms%fat(3,iat)=atoms%fat(3,iat)+parini%efield*0.5d0*atoms%qat(iat)
        enddo
    elseif((.not. poisson%point_particle) .and. trim(parini%bias_type)=='fixed_potdiff') then
        dipole=0.d0
        do iat=1,atoms%nat
            dipole=dipole+atoms%qat(iat)*atoms%rat(3,iat)
        enddo
        beta=-dipole*2.d0*pi/(poisson%cell(1)*poisson%cell(2)) 
        dv=parini%vu_ewald-parini%vl_ewald 

        do iat=1,atoms%nat
            atoms%fat(3,iat)=atoms%fat(3,iat)-(dv+2*beta)*0.5d0*atoms%qat(iat)/poisson%cell(3)
            atoms%fat(3,iat)=atoms%fat(3,iat)-(beta)*atoms%qat(iat)/poisson%cell(3)
            atoms%fat(3,iat)=atoms%fat(3,iat)+0.5d0*atoms%qat(iat)/poisson%cell(3)*dv
        enddo
    endif
    call f_free(wa)
    call f_free(wx)
    call f_free(wy)
    call f_free(wz)
    call f_free(vx)
    call f_free(vy)
    call f_free(vz)
    end associate
    end associate
    end associate
    call f_free(mboundg)
    call f_release_routine()
end subroutine force_gto_sym_ortho
!*****************************************************************************************
subroutine qgrad_gto_sym_ortho(parini,atoms,poisson,gausswidth,g)
    use mod_interface
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    use mod_parini, only: typ_parini
    use dynamic_memory
    implicit none
    type(typ_atoms), intent(in):: atoms
    type(typ_poisson), intent(inout):: poisson
    type(typ_parini), intent(in):: parini
    !local variables
    real(8), allocatable:: wx(:), wy(:), wz(:) !values of one dimensional Gaussian functions
    real(8):: rhoz, rhoyz, pi
    real(8):: g(atoms%nat) 
    real(8):: hgxinv, hgyinv, hgzinv, hgxhgyhgz
    real(8):: width_inv, width_inv_xat, width_inv_yat, width_inv_zat
    real(8):: width_inv_hgx, width_inv_hgy, width_inv_hgz,width
    real(8):: xat, yat, zat, fac
    integer:: iat, iw, ix, iy, iz, iatox, iatoy, iatoz, jx, jy, jz, iyt, izt
    real(8), intent(in):: gausswidth(atoms%nat)
    integer:: ngpx, ngpy, ngpz, iii
    real(8), allocatable:: wa(:,:,:)
    integer, allocatable:: mboundg(:,:,:)
    integer:: nbgpx, nbgpy, nbgpz, nagpx, nagpy, nagpz
    real(8):: ttg
    nbgpx=int(poisson%rgcut/poisson%hx)+2
    nbgpy=int(poisson%rgcut/poisson%hy)+2
    nbgpz=int(poisson%rgcut/poisson%hz)+2
    nagpx=nbgpx+1
    nagpy=nbgpy+1
    nagpz=0
    mboundg=f_malloc([1.to.2,-nbgpy.to.nbgpy,-nbgpz.to.nbgpz],id='mboundg')
    call get_glimitsphere(poisson,nbgpx,nbgpy,nbgpz,mboundg)
    ngpx=poisson%ngpx
    ngpy=poisson%ngpy
    ngpz=poisson%ngpz
    !wa(1-nagpx:ngpx+nagpx,1-nagpy:ngpy+nagpy,1-nagpz:ngpz+nagpz)=>poisson%pot
    wa=f_malloc([1-nagpx.to.ngpx+nagpx,1-nagpy.to.ngpy+nagpy,1-nagpz.to.ngpz+nagpz],id='wa')
    wx=f_malloc([-nbgpx.to.nbgpx],id='wx')
    wy=f_malloc([-nbgpy.to.nbgpy],id='wy')
    wz=f_malloc([-nbgpz.to.nbgpz],id='wz')
    pi=4.d0*atan(1.d0)
    hgxhgyhgz=poisson%hx*poisson%hy*poisson%hz
    hgxinv=1.d0/poisson%hx
    hgyinv=1.d0/poisson%hy
    hgzinv=1.d0/poisson%hz
    !---------------------------------------------------------------------------
    do iz=1-nagpz,ngpz+nagpz
        izt=iz+(sign(ngpz,-iz)+sign(ngpz,ngpz-iz))/2
        do iy=1-nagpy,ngpy+nagpy
            iyt=iy+(sign(ngpy,-iy)+sign(ngpy,ngpy-iy))/2
            do ix=1-nagpx,0
                wa(ix,iy,iz)=poisson%pot(ix+ngpx,iyt,izt)
            enddo
            do ix=1,ngpx
                wa(ix,iy,iz)=poisson%pot(ix,iyt,izt)
            enddo
            do ix=ngpx+1,ngpx+nagpx
                wa(ix,iy,iz)=poisson%pot(ix-ngpx,iyt,izt)
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
        width= poisson%alpha
        width_inv=1.d0/width
        fac=1.d0/(width*sqrt(pi))**3
        width_inv_hgx=width_inv*poisson%hx
        width_inv_hgy=width_inv*poisson%hy
        width_inv_hgz=width_inv*poisson%hz
    endif


    do iat=1,atoms%nat  
        !shift the Gaussian centers
        iatox=nint(atoms%rat(1,iat)*hgxinv)+1
        iatoy=nint(atoms%rat(2,iat)*hgyinv)+1
        iatoz=nint(atoms%rat(3,iat)*hgzinv)+1+nbgpz*iii
        xat=atoms%rat(1,iat)-(iatox-1)*poisson%hx
        yat=atoms%rat(2,iat)-(iatoy-1)*poisson%hy
        zat=atoms%rat(3,iat)-(iatoz-1-nbgpz*iii)*poisson%hz
        !construct the one-dimensional gaussians

        if(.not.parini%ewald) then
            width=gausswidth(iat)
            width_inv=1.d0/width
            fac=1.d0/(width*sqrt(pi))**3
            width_inv_hgx=width_inv*poisson%hx
            width_inv_hgy=width_inv*poisson%hy
            width_inv_hgz=width_inv*poisson%hz
        endif

        width_inv_xat=width_inv*xat
        width_inv_yat=width_inv*yat
        width_inv_zat=width_inv*zat
        do iw=-nbgpx,nbgpx
            wx(iw)=exp(-(width_inv_hgx*iw-width_inv_xat)**2)
        enddo
        do iw=-nbgpy,nbgpy
            wy(iw)=exp(-(width_inv_hgy*iw-width_inv_yat)**2)
        enddo
        do iw=-nbgpz,nbgpz
            wz(iw)=exp(-(width_inv_hgz*iw-width_inv_zat)**2)
        enddo

        ttg=0.d0
        do iz=-nbgpz,nbgpz
        jz=iatoz+iz
        rhoz=fac*wz(iz)
            do iy=-nbgpy,nbgpy
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
