!!--------------------------------------------------------------------------------------------------
subroutine best_charge_density(parini)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson_p3d, typ_ewald_p3d,typ_poisson
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    type(typ_poisson_p3d):: poisson_p3d
    type(typ_ewald_p3d):: ewald_p3d ,ewald_p3d_rough
    type(typ_atoms):: atoms
    integer:: istat, igpx, igpy, igpz, iat
   real(8):: cell(3), epot, rgcut_a, t1, t2, t3, t4, pi
    real(8):: ehartree 
    real(8),allocatable::  gausswidth(:)
    pi=4.d0*atan(1.d0)
    write(*,*) "OK num 1"
    call cube_read('electronic_density.cube',atoms,ewald_p3d%poisson_p3d%typ_poisson)
    ewald_p3d%poisson_p3d%ngpx=ewald_p3d%poisson_p3d%ngpx
    ewald_p3d%poisson_p3d%ngpy=ewald_p3d%poisson_p3d%ngpy
    ewald_p3d%poisson_p3d%ngpz=ewald_p3d%poisson_p3d%ngpz
    ewald_p3d%hgx=ewald_p3d%poisson_p3d%hx
    ewald_p3d%hgy=ewald_p3d%poisson_p3d%hy
    ewald_p3d%hgz=ewald_p3d%poisson_p3d%hz
    write(*,*) "hx,hy,hz",ewald_p3d%poisson_p3d%typ_poisson%hx,ewald_p3d%poisson_p3d%hy,ewald_p3d%poisson_p3d%hz
    write(*,*) "nx,ny,nz",ewald_p3d%poisson_p3d%ngpx,ewald_p3d%poisson_p3d%ngpy,ewald_p3d%poisson_p3d%ngpz
    write(*,*) "OK num 2"   
    atoms%boundcond='bulk'
    allocate(ewald_p3d%poisson_p3d%pot(ewald_p3d%poisson_p3d%ngpx+2,ewald_p3d%poisson_p3d%ngpy,ewald_p3d%poisson_p3d%ngpz),stat=istat)
    call construct_ewald_bps(parini,atoms,ewald_p3d)
    write(*,*) "OK num 3"
   ! call putgaussgrid(parini,atoms%boundcond,.true.,atoms%nat,atoms%rat,atoms%qat,gausswidth,ewald_p3d)
    write(*,*) "OK num 4"
    write(*,*)ewald_p3d%poisson_p3d%ngpx
    call cal_hartree_pot_bps(ewald_p3d,atoms,ehartree)
    write(*,*) "OK num 5"
    write(*,*) "ehartree",ehartree
    call destruct_ewald_bps(ewald_p3d)
    write(*,*) "OK num 6"
end subroutine best_charge_density
!!*****************************************************************************************
!subroutine best_charge_density(parini,ewald_p3d,atoms,symfunc,ann_arr,ekf)
!    use mod_interface
!    use mod_parini, only: typ_parini
!    use mod_atoms, only: typ_atoms
!    use mod_electrostatics, only: typ_poisson_p3d, typ_ewald_p3d
!    use mod_ann, only: typ_ann_arr, typ_symfunc, typ_ekf, typ_cent
!    use dynamic_memory
!    implicit none
!    type(typ_parini), intent(in):: parini
!    type(typ_atoms), intent(inout):: atoms
!    type(typ_poisson_p3d):: poisson_p3d
!    type(typ_ann_arr), intent(inout):: ann_arr
!    type(typ_symfunc), intent(inout):: symfunc
!    type(typ_ekf), intent(inout):: ekf
!    type(typ_ewald_p3d),intent(inout):: ewald_p3d
!    !local variables
!    type(typ_cent):: cent
!    integer:: iat, i, j, ng, ia
!    real(8):: epot_c, out_ann, ehartree, dx, dy, dz ,hardness, spring_const
!    real(8):: time1, time2, time3, time4, time5, time6, time7, time8
!    real(8):: tt1, tt2, tt3, fx_es, fy_es, fz_es, hinv(3,3), vol, fnet(3)
!    write(*,*) "OK num 1"
!    call cube_read('electronic_density.cube',atoms,ewald_p3d%poisson_p3d%typ_poisson)
!    write(*,*) "OK num 2"
!   ! call cal_hartree_pot_bps(ewald_p3d,atoms,ehartree)
!     write(*,*) "OK num 3"
!     write(*,*) "epot_es:", ehartree
!end subroutine best_charge_density
!*****************************************************************************************
subroutine gauss_gradient_rzx(parini,bc,nat,rxyz,cv,qat,gw,rgcut,ngx,ngy,ngz,pot,rgrad,qgrad)
    use mod_interface
    use mod_atoms, only: typ_atoms
    use mod_parini, only: typ_parini
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: bc
    integer, intent(in):: nat
    real(8), intent(in):: rxyz(3,nat)
    real(8), intent(in):: cv(3,3)
    real(8), intent(in):: qat(nat) 
    real(8), intent(in):: gw(nat)
    real(8), intent(in):: rgcut
    integer, intent(in):: ngx, ngy, ngz
    real(8), intent(inout):: pot(ngx,ngy,ngz)
    real(8), intent(out):: rgrad(3,nat), qgrad(nat)
    !local variables
    !work arrays to save the values of one dimensional gaussian function.
    real(8):: pi
    real(8):: facqiat, fac
    real(8):: cell(3) !dimensions of a smaller orthogonal cell for replication
    real(8):: vol
    real(8):: cvinv(3) !cell vectors of inverse coordinate, actual one at a time
    real(8):: htx, hty, htz
    real(8):: hxx, hxy, hxz, hyx, hyy, hyz, hzx, hzy, hzz
    real(8):: hxx_g, hxy_g, hxz_g, hyx_g, hyy_g, hyz_g, hzx_g, hzy_g, hzz_g
    real(8):: hrxinv, hryinv, hrzinv
    real(8):: cvinv_norm, vol_voxel, ttx, tty, ttz, ttq, tt1
    real(8):: dmsq, gwsq_inv
    real(8):: xred, yred, zred
    real(8):: ximg, yimg, zimg
    integer:: imgx, imgy, imgz
    integer:: ncellx, ncelly, ncellz
    integer:: iat, igx, igy, igz, jgx, jgy, jgz, igyt, igzt, igzysq
    integer:: iii, nlim, nlimsq, iix, iiy, iiz
    integer:: nbgx, nbgy, nbgz, nagx, nagy, nagz, nex, ney, nez
    integer:: finalx, igxs, igxf 
    real(8), allocatable:: wa(:,:,:)
    real(8), allocatable:: ratred(:,:)
    real(8), allocatable:: exponentval(:), expval(:)
    real(8), allocatable:: dmxarr(:), dmyarr(:), dmzarr(:)
    real(8):: sqpi, tt
    integer:: istartx, istarty, istartz

    allocate(ratred(3,nat))
    call rxyz_cart2int_alborz(nat,cv,rxyz,ratred)
    do iat=1,nat
        xred=ratred(1,iat)
        yred=ratred(2,iat)
        zred=ratred(3,iat)
        if(xred<0.d0) write(*,*) 'ATOM OUTSIDE: iat,sx ',iat,xred
        if(yred<0.d0) write(*,*) 'ATOM OUTSIDE: iat,sy ',iat,yred
        if(zred<0.d0) write(*,*) 'ATOM OUTSIDE: iat,sz ',iat,zred
        if(.not. (xred<1.d0)) write(*,*) 'ATOM OUTSIDE: iat,sx ',iat,xred
        if(.not. (yred<1.d0)) write(*,*) 'ATOM OUTSIDE: iat,sy ',iat,yred
        if(.not. (zred<1.d0)) write(*,*) 'ATOM OUTSIDE: iat,sz ',iat,zred
    enddo

    !reciprocal lattice to be used to determine the distance of corners of
    !the parallelepiped to its facets. Then those distances are used to
    !determine the number of grid points in each direction that are within
    !the cutoff of Gaussian function.
    call cell_vol(nat,cv,vol)
    vol=abs(vol)*nat
    call cross_product_alborz(cv(1,1),cv(1,2),cvinv)
    cvinv_norm=sqrt(cvinv(1)**2+cvinv(2)**2+cvinv(3)**2)
    cell(3)=vol/cvinv_norm
    call cross_product_alborz(cv(1,2),cv(1,3),cvinv)
    cvinv_norm=sqrt(cvinv(1)**2+cvinv(2)**2+cvinv(3)**2)
    cell(1)=vol/cvinv_norm
    call cross_product_alborz(cv(1,1),cv(1,3),cvinv)
    cvinv_norm=sqrt(cvinv(1)**2+cvinv(2)**2+cvinv(3)**2)
    cell(2)=vol/cvinv_norm
    if(parini%iverbose>1) then
        write(*,*) 'cell  ', cell(1),cell(2),cell(3)
    endif
    htx=cell(1)/real(ngx,8)
    hty=cell(2)/real(ngy,8)
    htz=cell(3)/real(ngz,8)
    nbgx=int(rgcut/htx)+2
    nbgy=int(rgcut/hty)+2
    nbgz=int(rgcut/htz)+2
    nagx=nbgx+1
    nagy=nbgy+1
    nagz=nbgz+1
    nlim = max(nagx,nagy,nagz)
    nlimsq = nlim**2
    !detemining the largest dimension for the pseudogrid.
    hxx=cv(1,1)/ngx ; hxy=cv(2,1)/ngx ; hxz=cv(3,1)/ngx
    hyx=cv(1,2)/ngy ; hyy=cv(2,2)/ngy ; hyz=cv(3,2)/ngy
    hzx=cv(1,3)/ngz ; hzy=cv(2,3)/ngz ; hzz=cv(3,3)/ngz

    wa=f_malloc0([1-nagx.to.ngx+nagx,1-nagy.to.ngy+nagy,1-nagz.to.ngz+nagz],id='wa')
    allocate(exponentval(-nbgx:nbgx),expval(-nbgx:nbgx))
    allocate(dmxarr(-nbgx:nbgx),dmyarr(-nbgx:nbgx),dmzarr(-nbgx:nbgx))
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
    !---------------------------------------------------------------------------
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
                wa(igx,igy,igz)=pot(iix,iiy,iiz)
            enddo
            do igx=igxs,igxf-1,ngx
                do iix=1,ngx
                    jgx=igx+iix-1
                    wa(jgx,igy,igz)=pot(iix,iiy,iiz)
                enddo
            enddo
            iix=0
            do igx=igxf,ngx+nagx
                iix=iix+1
                wa(igx,igy,igz)=pot(iix,iiy,iiz)
            enddo
        enddo
    enddo
    !-------------------------------------------------------
    vol_voxel=vol/(ngx*ngy*ngz)
    do iat=1,nat
        gwsq_inv=1.d0/gw(iat)**2
        fac=1.d0/(gw(iat)*sqrt(pi))**3
        imgx=nint(ratred(1,iat)*hrxinv)+1
        imgy=nint(ratred(2,iat)*hryinv)+1
        imgz=nint(ratred(3,iat)*hrzinv)+1
        facqiat=fac*qat(iat)
        ttq=0.d0
        ttx=0.d0
        tty=0.d0
        ttz=0.d0
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
                    dmxarr(igx)=ximg-rxyz(1,iat)
                    dmyarr(igx)=yimg-rxyz(2,iat)
                    dmzarr(igx)=zimg-rxyz(3,iat)
                    dmsq=dmxarr(igx)**2+dmyarr(igx)**2+dmzarr(igx)**2
                    exponentval(igx)=-dmsq*gwsq_inv
                enddo
                call vdexp( 2*iix+1, exponentval(-iix), expval(-iix) )
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
end subroutine gauss_gradient_rzx
!*****************************************************************************************

