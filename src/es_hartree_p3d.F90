!*****************************************************************************************
subroutine init_psolver_p3d(poisson)
    use mod_electrostatics, only: typ_poisson
    use dynamic_memory
    implicit none
    include 'fftw3.f'
    type(typ_poisson), intent(inout):: poisson
    !local variables
    integer:: iz
    call f_routine(id='init_psolver_p3d')
    !poisson%plan_f=f_malloc((/1.to.poisson%ngpz/),id='poisson%plan_f')
    !poisson%plan_b=f_malloc((/1.to.poisson%ngpz/),id='poisson%plan_b')
    allocate(poisson%plan_f(1:poisson%ngpz))
    allocate(poisson%plan_b(1:poisson%ngpz))
    do iz=1,poisson%ngpz
        call dfftw_plan_dft_r2c_2d(poisson%plan_f(iz),poisson%ngpx, &
            poisson%ngpy,poisson%rho(1,1,iz),poisson%pot(1,1,iz),fftw_estimate)
        call dfftw_plan_dft_c2r_2d(poisson%plan_b(iz),poisson%ngpx, &
            poisson%ngpy,poisson%pot(1,1,iz),poisson%pot(1,1,iz),fftw_estimate)
    enddo
    call f_release_routine()
end subroutine init_psolver_p3d
!*****************************************************************************************
subroutine fini_psolver_p3d(poisson)
    use mod_electrostatics, only: typ_poisson
    implicit none
    type(typ_poisson), intent(inout):: poisson
    !local variables
    integer:: istat, iz
    do iz=1,poisson%ngpz
        call dfftw_destroy_plan(poisson%plan_f(iz))
        call dfftw_destroy_plan(poisson%plan_b(iz))
    enddo
    deallocate(poisson%plan_f,stat=istat)
    if(istat/=0) then
        stop 'ERROR: failure deallocating plan_f of type typ_poisson'
    endif
    deallocate(poisson%plan_b,stat=istat)
    if(istat/=0) then
        stop 'ERROR: failure deallocating plan_b of type typ_poisson'
    endif
end subroutine fini_psolver_p3d
!*****************************************************************************************
subroutine get_psolver_p3d(parini,poisson,cell,hx,hy,hz,epot)
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_poisson), intent(inout):: poisson
    real(8):: cell(3) !cell array contains size of the simulation box.
    real(8):: hx, hy, hz
    real(8):: epot
    !local variables
    real(8):: valuengpxyinv
    integer:: igpx, igpy, igpz
    do igpz=1,poisson%ngpz
        call dfftw_execute(poisson%plan_f(igpz))
    enddo
    call solve_syslinequ_p3d(poisson,hz,cell)
    do igpz=1,poisson%ngpz
        call dfftw_execute(poisson%plan_b(igpz))
    enddo
    valuengpxyinv=1.d0/real(poisson%ngpx*poisson%ngpy,8)
    !DSCAL cannot be used due the first dimension of pot which is ngpx+2 NOT ngpx
    !call DSCAL(poisson%ngpx*poisson%ngpy*poisson%ngpz,valuengpxyinv,poisson%pot,1)
    epot=0.d0
    do igpz=1,poisson%ngpz
        do igpy=1,poisson%ngpy
            do igpx=1,poisson%ngpx
                poisson%pot(igpx,igpy,igpz)=poisson%pot(igpx,igpy,igpz)*valuengpxyinv
                epot=epot+poisson%pot(igpx,igpy,igpz)*poisson%rho(igpx,igpy,igpz)
            enddo
        enddo
    enddo
    epot=epot*0.5d0*hx*hy*hz
end subroutine get_psolver_p3d
!*****************************************************************************************
!This subroutine calculates the systems of linear equations using LAPACK routines.
subroutine solve_syslinequ_p3d(poisson,hz,cell)
    use mod_electrostatics, only: typ_poisson
    implicit none
    type(typ_poisson), intent(inout):: poisson
    !nem-1 is the degree of the polynomial used to calculate the 
    !integrals needed to prepare the right-hand side of the systems of 
    !linear equations.
    real(8):: hz, cell(3)
    !local variables
    integer, parameter:: nem=8 
    real(8):: pi, fourpi, fourpisq, gsq, gsqx, gsqy, g, hzsq
    real(8):: fourpisqcellxinvsq, fourpisqcellyinvsq, ciz
    real(8):: d(poisson%ngpz+2*8) !nem was replaced by 8 to be able to compile interface_mod.F90
    real(8):: e1(poisson%ngpz), e2(poisson%ngpz-1), c(poisson%ngpz)
    integer::ix,iy,iz,info,ixt,iyt
    !local variables
    !beta_grid is proportion to dipole moment and it is related to beta in the
    !P3D paper (i.e. poisson%beta) but multiplied by factor (-ngpx*ngpy)
    real(8):: beta_grid 
    pi=4.d0*atan(1.d0)
    hzsq=hz**2
    fourpi=4.d0*pi 
    fourpisq=fourpi*pi 
    fourpisqcellxinvsq=fourpisq/cell(1)**2 
    fourpisqcellyinvsq=fourpisq/cell(2)**2 
    d=0.d0
    !-----------------------------------------
    ix=1;iy=1
    gsq=0.d0
    g=0.d0
    e1=2.d0;e2=-1.d0
    call dpttrf(poisson%ngpz,e1,e2,info)
    if(info/= 0) write(*,*) 'Factorization failed',info
    do iz=1,poisson%ngpz
        d(nem+iz)=fourpi*poisson%pot(ix,iy,iz)
    enddo
    call get_beta_grid(hzsq,poisson%ngpz,d(1+nem),beta_grid)
    poisson%beta = beta_grid /(-poisson%ngpx*poisson%ngpy)
    call prepare00(poisson%ngpz,nem,d,c,hz)
    c(1)=c(1)-beta_grid
    c(poisson%ngpz)=c(poisson%ngpz)+beta_grid
    call dpttrs(poisson%ngpz,1,e1,e2,c,poisson%ngpz,info)
    if (info/= 0) write(*,*) 'Solution failed',info
    do iz=1,poisson%ngpz
        poisson%pot(ix,iy,iz)=c(iz)
    enddo
    !-----------------------------------------
    ix=poisson%ngpx+1;iy=1
    gsq=fourpisqcellxinvsq*((ix-1)/2)**2 + fourpisqcellyinvsq*(iy-1)**2
    g=sqrt(gsq)
    call fdcoeff(poisson%ngpz,e1,e2,g,gsq,hz,hzsq)
    call dpttrf(poisson%ngpz,e1,e2,info)
    if(info/= 0) write(*,*) 'Factorization failed',info
    do iz=1,poisson%ngpz
        d(nem+iz)=fourpi*poisson%pot(ix,iy,iz)
    enddo
    call prepare(poisson%ngpz,nem,d,c,gsq,hz,hzsq)
    call dpttrs(poisson%ngpz,1,e1,e2,c,poisson%ngpz,info)
    if (info/= 0) write(*,*) 'Solution failed',info
    do iz=1,poisson%ngpz
        poisson%pot(ix,iy,iz)=c(iz)
    enddo
    !-----------------------------------------
    ix=1;iy=poisson%ngpy/2+1
    gsq=fourpisqcellyinvsq*(iy-1)**2
    g=sqrt(gsq)
    call fdcoeff(poisson%ngpz,e1,e2,g,gsq,hz,hzsq)
    call dpttrf(poisson%ngpz,e1,e2,info)
    if(info/= 0) write(*,*) 'Factorization failed',info
    do iz=1,poisson%ngpz
        d(nem+iz)=fourpi*poisson%pot(ix,iy,iz)
    enddo

    call prepare(poisson%ngpz,nem,d,c,gsq,hz,hzsq)
    call dpttrs(poisson%ngpz,1,e1,e2,c,poisson%ngpz,info)
    if (info/= 0) write(*,*) 'Solution failed',info
    do iz=1,poisson%ngpz
        poisson%pot(ix,iy,iz)=c(iz)
    enddo
    !-----------------------------------------
    ix=poisson%ngpx+1;iy=poisson%ngpy/2+1
    gsq=fourpisqcellxinvsq*((ix-1)/2)**2 + fourpisqcellyinvsq*(iy-1)**2
    g=sqrt(gsq)
    call fdcoeff(poisson%ngpz,e1,e2,g,gsq,hz,hzsq)
    call dpttrf(poisson%ngpz,e1,e2,info)
    if(info/= 0) write(*,*) 'Factorization failed',info
    do iz=1,poisson%ngpz
        d(nem+iz)=fourpi*poisson%pot(ix,iy,iz)
    enddo
    call prepare(poisson%ngpz,nem,d,c,gsq,hz,hzsq)
    call dpttrs(poisson%ngpz,1,e1,e2,c,poisson%ngpz,info)
    if (info/= 0) write(*,*) 'Solution failed',info
    do iz=1,poisson%ngpz
        poisson%pot(ix,iy,iz)=c(iz)
    enddo
    !-----------------------------------------
    ix=1
    do iy=2,poisson%ngpy/2
        gsq=fourpisqcellyinvsq*(iy-1)**2
        g=sqrt(gsq)
        call fdcoeff(poisson%ngpz,e1,e2,g,gsq,hz,hzsq)
        call dpttrf(poisson%ngpz,e1,e2,info)
        if(info/= 0) write(*,*) 'Factorization failed',info
        do iz=1,poisson%ngpz
            d(nem+iz)=fourpi*poisson%pot(ix,iy,iz)
        enddo
        call prepare(poisson%ngpz,nem,d,c,gsq,hz,hzsq)
        call dpttrs(poisson%ngpz,1,e1,e2,c,poisson%ngpz,info)
        if (info/= 0) write(*,*) 'Solution failed',info
        do iz=1,poisson%ngpz
            ciz=c(iz)
            poisson%pot(ix,iy,iz)=ciz
            poisson%pot(ix,poisson%ngpy-iy+2,iz)=ciz
        enddo
        ixt=ix+1
        do iz=1,poisson%ngpz
            d(nem+iz)=fourpi*poisson%pot(ixt,iy,iz)
        enddo
        call prepare(poisson%ngpz,nem,d,c,gsq,hz,hzsq)
        call dpttrs(poisson%ngpz,1,e1,e2,c,poisson%ngpz,info)
        if (info/= 0) write(*,*) 'Solution failed',info
        do iz=1,poisson%ngpz
            ciz=c(iz)
            poisson%pot(ixt,iy,iz)=ciz
            poisson%pot(ixt,poisson%ngpy-iy+2,iz)=-ciz
        enddo
    enddo
    !-----------------------------------------
    ix=poisson%ngpx+1
    gsqx=fourpisqcellxinvsq*(poisson%ngpx/2)**2
    do iy=2,poisson%ngpy/2
        gsq=gsqx+fourpisqcellyinvsq*(iy-1)**2
        g=sqrt(gsq)
        call fdcoeff(poisson%ngpz,e1,e2,g,gsq,hz,hzsq)
        call dpttrf(poisson%ngpz,e1,e2,info)
        if(info/= 0) write(*,*) 'Factorization failed',info
        do iz=1,poisson%ngpz
            d(nem+iz)=fourpi*poisson%pot(ix,iy,iz)
        enddo
        call prepare(poisson%ngpz,nem,d,c,gsq,hz,hzsq)
        call dpttrs(poisson%ngpz,1,e1,e2,c,poisson%ngpz,info)
        if (info/= 0) write(*,*) 'Solution failed',info
        do iz=1,poisson%ngpz
            ciz=c(iz)
            poisson%pot(ix,iy,iz)=ciz
            poisson%pot(ix,poisson%ngpy-iy+2,iz)=ciz
        enddo
        ixt=ix+1
        do iz=1,poisson%ngpz
            d(nem+iz)=fourpi*poisson%pot(ixt,iy,iz)
        enddo
        call prepare(poisson%ngpz,nem,d,c,gsq,hz,hzsq)
        call dpttrs(poisson%ngpz,1,e1,e2,c,poisson%ngpz,info)
        if (info/= 0) write(*,*) 'Solution failed',info
        do iz=1,poisson%ngpz
            ciz=c(iz)
            poisson%pot(ixt,iy,iz)=ciz
            poisson%pot(ixt,poisson%ngpy-iy+2,iz)=-ciz
        enddo
    enddo
    !-----------------------------------------
    iy=1
    gsqy=fourpisqcellyinvsq*(iy-1)**2
    do ix=3,poisson%ngpx,2
        gsq=gsqy+fourpisqcellxinvsq*((ix-1)/2)**2
        g=sqrt(gsq)
        call fdcoeff(poisson%ngpz,e1,e2,g,gsq,hz,hzsq)
        call dpttrf(poisson%ngpz,e1,e2,info)
        if(info/= 0) write(*,*) 'Factorization failed',info
        do iz=1,poisson%ngpz
            d(nem+iz)=fourpi*poisson%pot(ix,iy,iz)
        enddo
        call prepare(poisson%ngpz,nem,d,c,gsq,hz,hzsq)
        call dpttrs(poisson%ngpz,1,e1,e2,c,poisson%ngpz,info)
        if (info/= 0) write(*,*) 'Solution failed',info
        do iz=1,poisson%ngpz
            poisson%pot(ix,iy,iz)=c(iz)
        enddo
        ixt=ix+1
        do iz=1,poisson%ngpz
            d(nem+iz)=fourpi*poisson%pot(ixt,iy,iz)
        enddo
        call prepare(poisson%ngpz,nem,d,c,gsq,hz,hzsq)
        call dpttrs(poisson%ngpz,1,e1,e2,c,poisson%ngpz,info)
        if (info/= 0) write(*,*) 'Solution failed',info
        do iz=1,poisson%ngpz
            poisson%pot(ixt,iy,iz)=c(iz)
        enddo
    enddo
    !-----------------------------------------
    iy=poisson%ngpy/2+1
    gsqy=fourpisqcellyinvsq*(iy-1)**2
    do ix=3,poisson%ngpx,2
        gsq=gsqy+fourpisqcellxinvsq*((ix-1)/2)**2
        g=sqrt(gsq)
        call fdcoeff(poisson%ngpz,e1,e2,g,gsq,hz,hzsq)
        call dpttrf(poisson%ngpz,e1,e2,info)
        if(info/= 0) write(*,*) 'Factorization failed',info
        do iz=1,poisson%ngpz
            d(nem+iz)=fourpi*poisson%pot(ix,iy,iz)
        enddo
        call prepare(poisson%ngpz,nem,d,c,gsq,hz,hzsq)
        call dpttrs(poisson%ngpz,1,e1,e2,c,poisson%ngpz,info)
        if (info/= 0) write(*,*) 'Solution failed',info
        do iz=1,poisson%ngpz
            poisson%pot(ix,iy,iz)=c(iz)
        enddo
        ixt=ix+1
        do iz=1,poisson%ngpz
            d(nem+iz)=fourpi*poisson%pot(ixt,iy,iz)
        enddo
        call prepare(poisson%ngpz,nem,d,c,gsq,hz,hzsq)
        call dpttrs(poisson%ngpz,1,e1,e2,c,poisson%ngpz,info)
        if (info/= 0) write(*,*) 'Solution failed',info
        do iz=1,poisson%ngpz
            poisson%pot(ixt,iy,iz)=c(iz)
        enddo
    enddo
    !-----------------------------------------
    do iy=2,poisson%ngpy/2
        iyt=poisson%ngpy-iy+2
        gsqy=fourpisqcellyinvsq*(iy-1)**2
        do ix=3,poisson%ngpx,2
            gsq=gsqy+fourpisqcellxinvsq*((ix-1)/2)**2
            g=sqrt(gsq)
            call fdcoeff(poisson%ngpz,e1,e2,g,gsq,hz,hzsq)
            call dpttrf(poisson%ngpz,e1,e2,info)
            if(info/= 0) write(*,*) 'Factorization failed',info
            do iz=1,poisson%ngpz
                d(nem+iz)=fourpi*poisson%pot(ix,iy,iz)
            enddo
            call prepare(poisson%ngpz,nem,d,c,gsq,hz,hzsq)
            call dpttrs(poisson%ngpz,1,e1,e2,c,poisson%ngpz,info)
            if (info/= 0) write(*,*) 'Solution failed',info
            do iz=1,poisson%ngpz
                poisson%pot(ix,iy,iz)=c(iz)
            enddo
            ixt=ix+1
            do iz=1,poisson%ngpz
                d(nem+iz)=fourpi*poisson%pot(ixt,iy,iz)
            enddo
            call prepare(poisson%ngpz,nem,d,c,gsq,hz,hzsq)
            call dpttrs(poisson%ngpz,1,e1,e2,c,poisson%ngpz,info)
            if (info/= 0) write(*,*) 'Solution failed',info
            do iz=1,poisson%ngpz
                poisson%pot(ixt,iy,iz)=c(iz)
            enddo
            do iz=1,poisson%ngpz
                d(nem+iz)=fourpi*poisson%pot(ix,iyt,iz)
            enddo
            call prepare(poisson%ngpz,nem,d,c,gsq,hz,hzsq)
            call dpttrs(poisson%ngpz,1,e1,e2,c,poisson%ngpz,info)
            if (info/= 0) write(*,*) 'Solution failed',info
            do iz=1,poisson%ngpz
                poisson%pot(ix,iyt,iz)=c(iz)
            enddo
            ixt=ix+1
            do iz=1,poisson%ngpz
                d(nem+iz)=fourpi*poisson%pot(ixt,iyt,iz)
            enddo
            call prepare(poisson%ngpz,nem,d,c,gsq,hz,hzsq)
            call dpttrs(poisson%ngpz,1,e1,e2,c,poisson%ngpz,info)
            if (info/= 0) write(*,*) 'Solution failed',info
            do iz=1,poisson%ngpz
                poisson%pot(ixt,iyt,iz)=c(iz)
            enddo
        enddo
    enddo
    !-----------------------------------------
end subroutine solve_syslinequ_p3d
!*****************************************************************************************
subroutine fdcoeff(ngpz,e1,e2,g,gsq,hz,hzsq)
    implicit none
    integer::ngpz,ngpzm1
    real(8)::e1(ngpz) !Diagonal elements of the matrix
    real(8)::e2(ngpz-1) !Offdiagonal elements of the matrix
    real(8)::gsqhzsq,a,hz,hzsq,gsq,g,diagonal_fl,diagonal,offdiagonal
    real(8)::hz4,hz6,asq,acb,t1,t2,t3,t4,t5,t6,hz4asq,hzsqa,hz6acb
    gsqhzsq=gsq*hzsq;hz4=hzsq**2;hz6=hz4*hzsq;a=-gsq;asq=a**2
    acb=asq*a;hz4asq=hz4*asq;hzsqa=hzsq*a;hz6acb=hz6*acb
    t1=(9.d0*hz4asq*(80080.d0 -1232.d0*hzsqa +3.d0*hz4asq))
    t2=-308880.d0+35640.d0*hzsqa-450.d0*hz4asq +hz6acb
    t3=11.d0*hz4asq*(65520.d0 - 624.d0*hzsqa + hz4asq)
    t4=-3603600.d0 + 120120.d0*hzsqa - 770.d0*hz4asq + hz6acb
    t5=t1/(56.d0*t2) 
    t6=t3/(72.d0*t4)
    diagonal=-2.d0 + (2.d0*hzsqa)/3.d0 - t5 - t6
    diagonal_fl=0.5d0*diagonal - g*hz
    offdiagonal=1.d0 + (hzsqa)/6.d0 - 0.5d0*(t5-t6)
    ngpzm1=ngpz-1
    e1(1)=-diagonal_fl   
    e1(2:ngpzm1)=-diagonal 
    e1(ngpz)=-diagonal_fl 
    e2(1:ngpzm1)=-offdiagonal   
end subroutine fdcoeff
!*****************************************************************************************
subroutine prepare00(ngpz,nem,f,c,hz)
    implicit none
    integer::ngpz,i,nem,n
    real(8)::f(ngpz+2*nem),c(ngpz),eta(6),hz 
    real(8)::coefftot1(17),coefftoti(17),coefftotn(17)
    eta=0.d0
    call prepcoeff(hz,eta,coefftot1,coefftoti,coefftotn)
    do i=1,ngpz
        c(i)=f(i   )*coefftoti( 1)+f(i+ 1)*coefftoti( 2)+ &
             f(i+ 2)*coefftoti( 3)+f(i+ 3)*coefftoti( 4)+ &
             f(i+ 4)*coefftoti( 5)+f(i+ 5)*coefftoti( 6)+ &
             f(i+ 6)*coefftoti( 7)+f(i+ 7)*coefftoti( 8)+ &
             f(i+ 8)*coefftoti( 9)+f(i+ 9)*coefftoti(10)+ &
             f(i+10)*coefftoti(11)+f(i+11)*coefftoti(12)+ &
             f(i+12)*coefftoti(13)+f(i+13)*coefftoti(14)+ &
             f(i+14)*coefftoti(15)+f(i+15)*coefftoti(16)+ &
             f(i+16)*coefftoti(17)
    enddo
end subroutine prepare00
!*****************************************************************************************
subroutine prepare(ngpz,nem,f,c,gsq,hz,hzsq)
    implicit none
    integer::ngpz,i,nem,n
    real(8)::f(ngpz+2*nem),c(ngpz),eta(6),gsq,hz,hzsq,a
    real(8)::coefftot1(17),coefftoti(17),coefftotn(17)
    real(8)::t1,t2,t3,t4
    real(8)::asq,acb,hzcb,hz4,hz5,hz6,hza,hzsqa,hzcba,hz4asq,hz5asq,hz6acb
    a=-gsq;asq=a**2;acb=asq*a;hzcb=hzsq*hz;hz4=hzsq**2;hz5=hz4*hz;hz6=hz4*hzsq
    hza=hz*a;hzsqa=hzsq*a;hzcba=hzcb*a;hz4asq=hz4*asq;hz5asq=hz5*asq
    hz6acb=hz6*acb
    t1=1.d0/(-308880.d0  + 35640.d0 *hzsqa - 450.d0*hz4asq + hz6acb)
    t2=1.d0/(-3603600.d0 + 120120.d0*hzsqa - 770.d0*hz4asq + hz6acb)
    t3=(- 9.d0*sqrt(1.5d0)*hza*(80080.d0*hz-1232.d0*hzcba+3.d0*hz5asq))
    t4=(-11.d0*sqrt(2.5d0)*hza*(65520.d0*hz- 624.d0*hzcba+     hz5asq))
    eta(1)=-hz*t3/28.d0*t1
    eta(2)=-hz*t4/12.d0*t2
    eta(3)=-hz*(-11.d0*hz4asq*(-234.d0+hzsqa))/(2.d0*sqrt(14.d0))*t1
    eta(4)=-hz*(-13.d0*hz4asq*(-330.d0+hzsqa))/(6.d0*sqrt( 2.d0))*t2
    eta(5)=-hz*(-13.d0*sqrt(5.5d0)*hz6acb)/28.d0*t1
    eta(6)=-hz*(- 5.d0*sqrt(6.5d0)*hz6acb)/12.d0*t2
    call prepcoeff(hz,eta,coefftot1,coefftoti,coefftotn)
    c(1)=f( 1)*coefftot1( 1)+f( 2)*coefftot1( 2)+f( 3)*coefftot1( 3)+ &
         f( 4)*coefftot1( 4)+f( 5)*coefftot1( 5)+f( 6)*coefftot1( 6)+ &
         f( 7)*coefftot1( 7)+f( 8)*coefftot1( 8)+f( 9)*coefftot1( 9)+ &
         f(10)*coefftot1(10)+f(11)*coefftot1(11)+f(12)*coefftot1(12)+ &
         f(13)*coefftot1(13)+f(14)*coefftot1(14)+f(15)*coefftot1(15)+ &
         f(16)*coefftot1(16)+f(17)*coefftot1(17)
    do i=2,ngpz-1
        c(i)=f(i   )*coefftoti( 1)+f(i+ 1)*coefftoti( 2)+ &
             f(i+ 2)*coefftoti( 3)+f(i+ 3)*coefftoti( 4)+ &
             f(i+ 4)*coefftoti( 5)+f(i+ 5)*coefftoti( 6)+ &
             f(i+ 6)*coefftoti( 7)+f(i+ 7)*coefftoti( 8)+ &
             f(i+ 8)*coefftoti( 9)+f(i+ 9)*coefftoti(10)+ &
             f(i+10)*coefftoti(11)+f(i+11)*coefftoti(12)+ &
             f(i+12)*coefftoti(13)+f(i+13)*coefftoti(14)+ &
             f(i+14)*coefftoti(15)+f(i+15)*coefftoti(16)+ &
             f(i+16)*coefftoti(17)
    enddo
    c(ngpz)=f(ngpz   )*coefftotn( 1)+f(ngpz+ 1)*coefftotn( 2)+ &
            f(ngpz+ 2)*coefftotn( 3)+f(ngpz+ 3)*coefftotn( 4)+ &
            f(ngpz+ 4)*coefftotn( 5)+f(ngpz+ 5)*coefftotn( 6)+ &
            f(ngpz+ 6)*coefftotn( 7)+f(ngpz+ 7)*coefftotn( 8)+ &
            f(ngpz+ 8)*coefftotn( 9)+f(ngpz+ 9)*coefftotn(10)+ &
            f(ngpz+10)*coefftotn(11)+f(ngpz+11)*coefftotn(12)+ &
            f(ngpz+12)*coefftotn(13)+f(ngpz+13)*coefftotn(14)+ &
            f(ngpz+14)*coefftotn(15)+f(ngpz+15)*coefftotn(16)+ &
            f(ngpz+16)*coefftotn(17)
end subroutine prepare
!*****************************************************************************************
subroutine prepcoeff(hz,eta,coefftot1,coefftoti,coefftotn)
    implicit none
    real(8)::hz,coeff(16,8),coefftot1(17),coefftoti(17),coefftotn(17),eta(6)
    real(8)::hzsq,hzeta(6)

    hzsq=hz**2
    hzeta(1:6)=hz*eta(1:6)

    coeff( 1,1)=-2.07475957726675297654306541812d-6*hzsq        
    coeff( 2,1)=0.0000359890672760720953609220026898d0*hzsq        
    coeff( 3,1)=-0.000298648715170959195883980612687d0*hzsq        
    coeff( 4,1)=0.00158901097031623523663240784657d0*hzsq        
    coeff( 5,1)=-0.00617538689758571634933290706711d0*hzsq        
    coeff( 6,1)=0.0192996081122292884486187733413d0*hzsq        
    coeff( 7,1)=-0.0557823926274722850104545992778d0*hzsq        
    coeff( 8,1)=0.37940939599653687502368641491d0*hzsq        
    coeff( 9,1)=0.197580778610955699152088819067d0*hzsq        
    coeff(10,1)=-0.0489457731460799627205850188871d0*hzsq        
    coeff(11,1)=0.0178683823757735831245790785162d0*hzsq        
    coeff(12,1)=-0.0058462283895129089051565449532d0*hzsq        
    coeff(13,1)=0.00152286284046759745526906712369d0*hzsq        
    coeff(14,1)=-0.000288447535444826114566232835763d0*hzsq        
    coeff(15,1)=0.0000349466931209943808501750196947d0*hzsq        
    coeff(16,1)=-2.02259583241986812983112737582d-6*hzsq        

    coeff( 1,2)=-2.02259583241986812983112737582d-6*hzsq        
    coeff( 2,2)=0.0000349466931209943808501750196947d0*hzsq        
    coeff( 3,2)=-0.000288447535444826114566232835763d0*hzsq        
    coeff( 4,2)=0.00152286284046759745526906712369d0*hzsq        
    coeff( 5,2)=-0.0058462283895129089051565449532d0*hzsq        
    coeff( 6,2)=0.0178683823757735831245790785162d0*hzsq        
    coeff( 7,2)=-0.0489457731460799627205850188871d0*hzsq        
    coeff( 8,2)=0.197580778610955699152088819067d0*hzsq        
    coeff( 9,2)=0.37940939599653687502368641491d0*hzsq        
    coeff(10,2)=-0.0557823926274722850104545992778d0*hzsq        
    coeff(11,2)=0.0192996081122292884486187733413d0*hzsq        
    coeff(12,2)=-0.00617538689758571634933290706711d0*hzsq        
    coeff(13,2)=0.00158901097031623523663240784657d0*hzsq        
    coeff(14,2)=-0.000298648715170959195883980612687d0*hzsq        
    coeff(15,2)=0.0000359890672760720953609220026898d0*hzsq        
    coeff(16,2)=-2.07475957726675297654306541812d-6*hzsq        

    coeff( 1,3)=2.02995291983762385617462574996d-6         *hzeta(1)
    coeff( 2,3)=-0.0000351409287690181621486755726107d0    *hzeta(1)     
    coeff( 3,3)=0.000290805426942967652041646808273d0      *hzeta(1)   
    coeff( 4,3)=-0.0015410601225584615873251953534d0       *hzeta(1)  
    coeff( 5,3)=0.00595073265312861074872981931172d0       *hzeta(1)  
    coeff( 6,3)=-0.0183784369137370031736020870667d0       *hzeta(1)  
    coeff( 7,3)=0.0515729699894665248765145948659d0        *hzeta(1) 
    coeff( 8,3)=-0.241986045289324966161173283844d0        *hzeta(1) 
    coeff( 9,3)=-0.241986045289324966161173283844d0        *hzeta(1) 
    coeff(10,3)=0.0515729699894665248765145948659d0        *hzeta(1) 
    coeff(11,3)=-0.0183784369137370031736020870667d0       *hzeta(1)  
    coeff(12,3)=0.00595073265312861074872981931172d0       *hzeta(1)  
    coeff(13,3)=-0.0015410601225584615873251953534d0       *hzeta(1)  
    coeff(14,3)=0.000290805426942967652041646808273d0      *hzeta(1)   
    coeff(15,3)=-0.0000351409287690181621486755726107d0    *hzeta(1)     
    coeff(16,3)=2.02995291983762385617462574996d-6         *hzeta(1)

    coeff( 1,4)=-2.40225990404051644648415508219d-8        *hzeta(2)
    coeff( 2,4)=4.79971611784916615735237313764d-7         *hzeta(2)
    coeff( 3,4)=-4.6962108400065933158894431998d-6         *hzeta(2)
    coeff( 4,4)=0.0000304403395187933758612590850302d0     *hzeta(2)   
    coeff( 5,4)=-0.000151359054696339561892657097308d0     *hzeta(2)   
    coeff( 6,4)=0.000656902635771458952244692183703d0      *hzeta(2)  
    coeff( 7,4)=-0.00311559134264715472402190126528d0      *hzeta(2)  
    coeff( 8,4)=0.0595982178730348695159249349471d0        *hzeta(2)
    coeff( 9,4)=-0.0595982178730348695159249349471d0       *hzeta(2) 
    coeff(10,4)=0.00311559134264715472402190126528d0       *hzeta(2) 
    coeff(11,4)=-0.000656902635771458952244692183703d0     *hzeta(2)   
    coeff(12,4)=0.000151359054696339561892657097308d0      *hzeta(2)  
    coeff(13,4)=-0.0000304403395187933758612590850302d0    *hzeta(2)    
    coeff(14,4)=4.6962108400065933158894431998d-6          *hzeta(2)
    coeff(15,4)=-4.79971611784916615735237313764d-7        *hzeta(2)
    coeff(16,4)=2.40225990404051644648415508219d-8         *hzeta(2)

    coeff( 1,5)=-2.42252529077693364853270248854d-7        *hzeta(3)
    coeff( 2,5)=4.1911252055266878343545750128d-6          *hzeta(3)
    coeff( 3,5)=-0.0000346494722679002303083813075042d0    *hzeta(3)    
    coeff( 4,5)=0.000183305231435073822450800905397d0      *hzeta(3)  
    coeff( 5,5)=-0.000705432142530625000849885262722d0     *hzeta(3)   
    coeff( 6,5)=0.00216052102232372719680853232496d0       *hzeta(3) 
    coeff( 7,5)=-0.00587068801655842125788540835559d0      *hzeta(3)  
    coeff( 8,5)=0.0042629945049216964753148403907d0        *hzeta(3)
    coeff( 9,5)=0.0042629945049216964753148403907d0        *hzeta(3)
    coeff(10,5)=-0.00587068801655842125788540835559d0      *hzeta(3)  
    coeff(11,5)=0.00216052102232372719680853232496d0       *hzeta(3) 
    coeff(12,5)=-0.000705432142530625000849885262722d0     *hzeta(3)   
    coeff(13,5)=0.000183305231435073822450800905397d0      *hzeta(3)  
    coeff(14,5)=-0.0000346494722679002303083813075042d0    *hzeta(3)    
    coeff(15,5)=4.1911252055266878343545750128d-6          *hzeta(3)
    coeff(16,5)=-2.42252529077693364853270248854d-7        *hzeta(3)
           
    coeff( 1,6)=5.83127337270514931512673488815d-9          *hzeta(4)
    coeff( 2,6)=-1.16442996091987472838683859737d-7         *hzeta(4) 
    coeff( 3,6)=1.13829029007821803399163899338d-6          *hzeta(4)
    coeff( 4,6)=-7.36664860376531573168334444592d-6         *hzeta(4) 
    coeff( 5,6)=0.0000365144904971176539462378534455d0      *hzeta(4)    
    coeff( 6,6)=-0.000157247997860582845557676582966d0      *hzeta(4)    
    coeff( 7,6)=0.000723808477916397455397344358808d0       *hzeta(4)   
    coeff( 8,6)=-0.00158558194383446901290889707547d0       *hzeta(4)   
    coeff( 9,6)=0.00158558194383446901290889707547d0        *hzeta(4)  
    coeff(10,6)=-0.000723808477916397455397344358808d0      *hzeta(4)    
    coeff(11,6)=0.000157247997860582845557676582966d0       *hzeta(4)   
    coeff(12,6)=-0.0000365144904971176539462378534455d0     *hzeta(4)     
    coeff(13,6)=7.36664860376531573168334444592d-6          *hzeta(4)
    coeff(14,6)=-1.13829029007821803399163899338d-6         *hzeta(4) 
    coeff(15,6)=1.16442996091987472838683859737d-7          *hzeta(4)
    coeff(16,6)=-5.83127337270514931512673488815d-9         *hzeta(4) 
           
    coeff( 1,7)=6.7887145393459615466705721742d-9          *hzeta(5)
    coeff( 2,7)=-1.16641515637242695500627487408d-7        *hzeta(5)  
    coeff( 3,7)=9.53656214378404625298477483164d-7         *hzeta(5) 
    coeff( 4,7)=-4.9467919183633031318985593389d-6         *hzeta(5) 
    coeff( 5,7)=0.0000182859695229315051184003611123d0     *hzeta(5)     
    coeff( 6,7)=-0.0000503466965190755747908449801666d0    *hzeta(5)      
    coeff( 7,7)=0.0000787467362088750652838476061824d0     *hzeta(5)     
    coeff( 8,7)=-0.0000425830207076482003708489483571d0    *hzeta(5)      
    coeff( 9,7)=-0.0000425830207076482003708489483571d0    *hzeta(5)      
    coeff(10,7)=0.0000787467362088750652838476061824d0     *hzeta(5)     
    coeff(11,7)=-0.0000503466965190755747908449801666d0    *hzeta(5)      
    coeff(12,7)=0.0000182859695229315051184003611123d0     *hzeta(5)     
    coeff(13,7)=-4.9467919183633031318985593389d-6         *hzeta(5) 
    coeff(14,7)=9.53656214378404625298477483164d-7         *hzeta(5) 
    coeff(15,7)=-1.16641515637242695500627487408d-7        *hzeta(5)  
    coeff(16,7)=6.7887145393459615466705721742d-9          *hzeta(5)
    
    coeff( 1,8)=-1.86362149128123761093706596947d-10       *hzeta(6)
    coeff( 2,8)=3.69550692262139012695831864162d-9         *hzeta(6)
    coeff( 3,8)=-3.57211805024307870448999045143d-8        *hzeta(6)
    coeff( 4,8)=2.26613525160693991429104581318d-7         *hzeta(6)
    coeff( 5,8)=-1.07834084132734977617139382693d-6        *hzeta(6)
    coeff( 6,8)=4.16778896192710562551682611264d-6         *hzeta(6)
    coeff( 7,8)=-0.108861142873203668915883477502d-4       *hzeta(6)
    coeff( 8,8)=0.176759490429403574997785746502d-4        *hzeta(6)
    coeff( 9,8)=-0.176759490429403574997785746502d-4       *hzeta(6)
    coeff(10,8)=0.108861142873203668915883477502d-4        *hzeta(6)
    coeff(11,8)=-4.16778896192710562551682611264d-6        *hzeta(6)
    coeff(12,8)=1.07834084132734977617139382693d-6         *hzeta(6)
    coeff(13,8)=-2.26613525160693991429104581318d-7        *hzeta(6)
    coeff(14,8)=3.57211805024307870448999045143d-8         *hzeta(6)
    coeff(15,8)=-3.69550692262139012695831864162d-9        *hzeta(6)
    coeff(16,8)=1.86362149128123761093706596947d-10        *hzeta(6)
    !---------------------------------------------------------------------------
    coefftot1(1)=coeff(1,2)+coeff(1,3)- &
                 coeff(1,4)+coeff(1,5)- &
                 coeff(1,6)+coeff(1,7)- &
                 coeff(1,8)
    coefftot1(2:16)=(coeff(2:16,2)+coeff(1:15,1))+ &
                   ( coeff(2:16,3))+ & 
                   (-coeff(2:16,4))+ &
                   ( coeff(2:16,5))+ &
                   (-coeff(2:16,6))+ &
                   ( coeff(2:16,7))+ &
                   (-coeff(2:16,8))
    coefftot1(17)=coeff(16,1)
    !---------------------------------------------------------------------------
    coefftoti(1)=coeff(1,2)+coeff(1,3)+ &
                 coeff(1,4)+coeff(1,5)+ &
                 coeff(1,6)+coeff(1,7)+ &
                 coeff(1,8)
    coefftoti(2:16)=(coeff(2:16,2)+coeff(1:15,1))+ &
                    (coeff(2:16,3)+coeff(1:15,3))+ & 
                    (coeff(2:16,4)-coeff(1:15,4))+ &
                    (coeff(2:16,5)+coeff(1:15,5))+ &
                    (coeff(2:16,6)-coeff(1:15,6))+ &
                    (coeff(2:16,7)+coeff(1:15,7))+ &
                    (coeff(2:16,8)-coeff(1:15,8))
    coefftoti(17)=coeff(16,1)+coeff(16,3)- &
                  coeff(16,4)+coeff(16,5)- &
                  coeff(16,6)+coeff(16,7)- &
                  coeff(16,8)
    !---------------------------------------------------------------------------
    coefftotn(1)=coeff(1,2)
    coefftotn(2:16)=(coeff(2:16,2)+coeff(1:15,1))+ &
                   (coeff(1:15,3))+ & 
                   (coeff(1:15,4))+ &
                   (coeff(1:15,5))+ &
                   (coeff(1:15,6))+ &
                   (coeff(1:15,7))+ &
                   (coeff(1:15,8))
    coefftotn(17)=coeff(16,1)+coeff(16,3)+ &
                  coeff(16,4)+coeff(16,5)+ &
                  coeff(16,6)+coeff(16,7)+ &
                  coeff(16,8)
    !---------------------------------------------------------------------------
end subroutine prepcoeff
!*****************************************************************************************
!This subroutine calculates the dipole moment which defines the boundary 
!condition for the single ODE in which gamma=0
subroutine get_beta_grid(hzsq,ngpz,analc00,beta_grid)
    use yaml_output
    implicit none
    real(8), intent(in):: hzsq, analc00(ngpz)
    integer, intent(in):: ngpz
    real(8), intent(out):: beta_grid
    !local variables
    integer:: iz
    beta_grid=0.d0
    do iz=1,ngpz-1
        beta_grid=beta_grid+iz*analc00(iz)
    enddo
    beta_grid=beta_grid*0.5d0*hzsq
    call yaml_map('inside get_beta_grid',beta_grid,fmt='(e30.17)')
    !write(*,'(a22,e30.17)') 'inside get_beta_grid beta_grid=',beta_grid
end subroutine get_beta_grid
!*****************************************************************************************

