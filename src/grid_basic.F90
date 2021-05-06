!*****************************************************************************************
!This subroutine determines the limits of grids in a sphere.
subroutine get_glimitsphere(hx,hy,hz,nbgpx,nbgpy,nbgpz,mboundg)
    use mod_electrostatics, only: typ_poisson
    implicit none
    real(8), intent(in):: hx, hy, hz
    integer, intent(in):: nbgpx, nbgpy, nbgpz
    integer, intent(out):: mboundg(1:2,-nbgpy:nbgpy,-nbgpz:nbgpz)
    !local variables
    integer:: ix, iy, iz
    real(8):: rgcut, rgcutsq
    rgcut=max(hx*nbgpx,hy*nbgpy,hz*nbgpz)
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
        if(ix**2*hx**2+iy**2*hy**2+iz**2*hz**2<=rgcutsq) then
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
end subroutine get_glimitsphere
!*****************************************************************************************
subroutine init_grid_param(nat,rxyz,cv,rgcut,ngx,ngy,ngz,ratred,vol,nlimsq,nagx,nagy,nagz,nbgx,nbgy,nbgz)
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
end subroutine init_grid_param
!*****************************************************************************************
subroutine charge_back_to_cell(ngx,ngy,ngz,nagx,nagy,nagz,ibcx,wa,rho)
    use mod_atoms, only: typ_atoms
    use mod_parini, only: typ_parini
    use dynamic_memory
    implicit none
    integer, intent(in):: ngx, ngy, ngz, nagx, nagy, nagz, ibcx
    real(8), intent(in):: wa(1-nagx:ngx+nagx,1-nagy:ngy+nagy,1-nagz:ngz+nagz)
    real(8), intent(inout):: rho(ngx,ngy,ngz)
    !local variables
    !work arrays to save the values of one dimensional gaussian function.
    integer:: igx, igy, igz, jgx
    integer:: iix, iiy, iiz
    integer:: ifinalx, igxs, igxf
    integer:: istartx, istarty, istartz
    istartx = modulo((1-nagx-1),ngx)+1
    ifinalx = modulo((ngx+nagx-1),ngx)+1
    igxs = 1-nagx+(ngx-istartx+1)
    igxf = ngx+nagx-ifinalx+1+ibcx*ngx
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
            enddo
            do igx=igxs,igxf-1,ngx
                do iix=1,ngx
                    jgx=igx+iix-1
                    rho(iix,iiy,iiz)=rho(iix,iiy,iiz)+wa(jgx,igy,igz)
                enddo
            enddo
            iix=0
            do igx=igxf,ngx+nagx
                iix=iix+1
                rho(iix,iiy,iiz)=rho(iix,iiy,iiz)+wa(igx,igy,igz)
            enddo
        enddo
    enddo
end subroutine charge_back_to_cell
!*****************************************************************************************
subroutine potential_on_extended_grid(lda,ngx,ngy,ngz,nagx,nagy,nagz,ibcx,pot,wa)
    use mod_atoms, only: typ_atoms
    use mod_parini, only: typ_parini
    use dynamic_memory
    implicit none
    integer, intent(in):: lda, ngx, ngy, ngz, nagx, nagy, nagz, ibcx
    real(8), intent(in):: pot(lda,ngy,ngz)
    real(8), intent(out):: wa(1-nagx:ngx+nagx,1-nagy:ngy+nagy,1-nagz:ngz+nagz)
    !local variables
    !work arrays to save the values of one dimensional gaussian function.
    integer:: igx, igy, igz, jgx
    integer:: iix, iiy, iiz
    integer:: ifinalx, igxs, igxf
    integer:: istartx, istarty, istartz
    istartx = modulo((1-nagx-1),ngx)+1
    ifinalx = modulo((ngx+nagx-1),ngx)+1
    igxs = 1-nagx+(ngx-istartx+1)
    igxf = ngx+nagx-ifinalx+1+ibcx*ngx
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
end subroutine potential_on_extended_grid
!*****************************************************************************************
