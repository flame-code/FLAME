!*****************************************************************************************
!This subroutine determines the limits of grids in a sphere.
subroutine determine_glimitsphere(poisson,nbgpx,nbgpy,nbgpz,mboundg)
    use mod_interface
    use mod_electrostatics, only: typ_poisson
    implicit none
    type(typ_poisson), intent(inout):: poisson
    integer, intent(in):: nbgpx, nbgpy, nbgpz
    integer, intent(out):: mboundg(1:2,-nbgpy:nbgpy,-nbgpz:nbgpz)
    !local variables
    integer:: ix, iy, iz
    real(8):: rgcut, rgcutsq
    rgcut=max(poisson%hx*nbgpx,poisson%hy*nbgpy,poisson%hz*nbgpz)
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
        if(ix**2*poisson%hx**2+iy**2*poisson%hy**2+iz**2*poisson%hz**2<=rgcutsq) then
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
end subroutine determine_glimitsphere
!*****************************************************************************************
subroutine gauss_init(nat,rxyz,cv,rgcut,ngx,ngy,ngz,ratred,vol,nlimsq,nagx,nagy,nagz,nbgx,nbgy,nbgz)
    use mod_interface
    implicit none
    integer, intent(in):: nat
    real(8), intent(in):: rxyz(3,nat)
    real(8), intent(in):: cv(3,3)
    real(8), intent(in):: rgcut
    integer, intent(in):: ngx, ngy, ngz
    real(8), intent(out):: ratred(3,nat), vol
    integer, intent(out):: nlimsq, nagx, nagy, nagz, nbgx, nbgy, nbgz
    !local variables
    real(8):: htx, hty, htz, xred, yred, zred, cvinv_norm
    integer:: iat, nlim
    real(8):: cvinv(3) !cell vectors of inverse coordinate, actual one at a time
    real(8):: cell(3) !dimensions of a smaller orthogonal cell for replication
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
    !if(parini%iverbose>1) then
    !    write(*,*) 'cell  ', cell(1),cell(2),cell(3)
    !endif
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
end subroutine gauss_init
!*****************************************************************************************
