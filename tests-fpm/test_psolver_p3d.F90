!*****************************************************************************************
subroutine test_psolver_p3d()
    use iso_fortran_env, only: error_unit, output_unit
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson
    use mod_colors, only: green_passed, red_failed
    implicit none
    real(8):: errmax, err_epot
    type(typ_parini):: parini
    type(typ_poisson):: poisson
    real(8), allocatable:: pot(:,:,:)
    integer:: ix, iy, iz, ios, ii
    real(8):: cell(3), x, y, z, a, b, c, hx, hy, hz, epot, epot_ref, tt
    real(8):: time1, time2
    a=1.05d0
    b=0.95d0
    c=0.97d0
    cell(1)=1.d0
    cell(2)=1.d0
    cell(3)=12.d0
    poisson%ngpx=24
    poisson%ngpy=24
    poisson%ngpz=120
    poisson%lda=poisson%ngpx+2
    poisson%hgrid(1,1)=cell(1)/real(poisson%ngpx,kind=8)
    poisson%hgrid(2,2)=cell(2)/real(poisson%ngpy,kind=8)
    poisson%hgrid(3,3)=cell(3)/real(poisson%ngpz,kind=8)
    hx=poisson%hgrid(1,1)
    hy=poisson%hgrid(2,2)
    hz=poisson%hgrid(3,3)
    allocate(poisson%rho(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    allocate(poisson%pot(poisson%lda,poisson%ngpy,poisson%ngpz))
    allocate(pot(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    call init_psolver_p3d(poisson)
    call get_density(cell,a,b,c,poisson)
    !call cpu_time(time1)
    call get_psolver_p3d(parini,poisson,cell,hx,hy,hz,epot)
    !call cpu_time(time2)
    !write(*,*) 'time=',time2-time1
    !open(unit=761,file='/home/ghasemi/pot_p3d.txt',status='replace',iostat=ios)
    !ii=0
    !do ix=1,poisson%ngpx
    !do iy=1,poisson%ngpy
    !do iz=1,poisson%ngpz
    !    !write(761,'(4es14.5)') (ix-1)*hx,(iy-1)*hy,(iz-1)*hz,poisson%rho(ix,iy,iz)
    !    ii=ii+1
    !    tt=poisson%pot(ix,iy,iz)
    !    !if(abs(tt)<0.0000005d0) tt=0.d0
    !    !write(761,'(f13.6)',advance='no') tt
    !    write(761,'(es24.15)',advance='no') tt
    !    if(mod(ii,6)==0) write(761,*)
    !enddo
    !enddo
    !enddo
    !if(mod(ii,6)/=0) write(761,*)
    !close(761)
    call get_potential(cell,a,b,c,poisson,pot)
    !open(unit=761,file='/home/ghasemi/pot_analytic.txt',status='replace',iostat=ios)
    !ii=0
    !do ix=1,poisson%ngpx
    !do iy=1,poisson%ngpy
    !do iz=1,poisson%ngpz
    !    ii=ii+1
    !    tt=pot(ix,iy,iz)
    !    !if(abs(tt)<0.0000005d0) tt=0.d0
    !    !write(761,'(f13.6)',advance='no') tt
    !    write(761,'(es24.15)',advance='no') tt
    !    if(mod(ii,6)==0) write(761,*)
    !enddo
    !enddo
    !enddo
    !if(mod(ii,6)/=0) write(761,*)
    !close(761)

    errmax=0.d0
    epot_ref=0.d0
    do ix=1,poisson%ngpx
    do iy=1,poisson%ngpy
    do iz=1,poisson%ngpz
        errmax=max(errmax,abs(poisson%pot(ix,iy,iz)-pot(ix,iy,iz)))
        epot_ref=epot_ref+poisson%rho(ix,iy,iz)*pot(ix,iy,iz)
    enddo
    enddo
    enddo
    epot_ref=epot_ref*0.5d0*hx*hy*hz
    if(errmax<1.d-12) then
        write(output_unit,'(2a)') green_passed,' in test_psolver_p3d: errmax'
    else
        write(error_unit,'(2a,es14.5)') red_failed,' in test_psolver_p3d: errmax=  ',errmax
        call exit(1)
    end if
    err_epot=abs(epot-epot_ref)
    if(err_epot<1.d-13) then
        write(output_unit,'(2a)') green_passed,' in test_psolver_p3d: err_epot'
    else
        write(error_unit,'(2a,es14.5)') red_failed,' in test_psolver_p3d: err_epot=  ',err_epot
        call exit(1)
    end if

    call fini_psolver_p3d(poisson)
    deallocate(poisson%rho)
    deallocate(poisson%pot)
    deallocate(pot)
end subroutine test_psolver_p3d
!*****************************************************************************************
subroutine get_potential(cell,a,b,c,poisson,pot)
    use mod_electrostatics, only: typ_poisson
    implicit none
    real(8), intent(in):: cell(3), a, b, c
    type(typ_poisson), intent(in):: poisson
    real(8), intent(out):: pot(poisson%ngpx,poisson%ngpy,poisson%ngpz)
    !local variables
    integer:: ix, iy, iz
    real(8):: hx, hy, hz, pi, x, y, z
    real(8):: z0, tz, tyz
    pi=4.d0*atan(1.d0)
    hx=poisson%hgrid(1,1)
    hy=poisson%hgrid(2,2)
    hz=poisson%hgrid(3,3)
    z0=cell(3)/2.d0
    do iz=1,poisson%ngpz
        z=(iz-1)*hz
        tz=exp(-(z-z0)**2/c**2)
        do iy=1,poisson%ngpy
            y=(iy-1)*hy
            tyz=tz*sin(b*sin(2.d0*pi*y/cell(2)))
            do ix=1,poisson%ngpx
                x=(ix-1)*hx
                pot(ix,iy,iz)=tyz*sin(a*sin(2.d0*pi*x/cell(1)))
            enddo
        enddo
    enddo
end subroutine get_potential
!*****************************************************************************************
subroutine get_density(cell,a,b,c,poisson)
    use mod_electrostatics, only: typ_poisson
    implicit none
    real(8), intent(in):: cell(3), a, b, c
    type(typ_poisson), intent(inout):: poisson
    !local variables
    integer:: ix, iy, iz
    real(8):: hx, hy, hz, pi, x, y, z
    real(8):: z0, t1, t2, t3
    pi=4.d0*atan(1.d0)
    hx=poisson%hgrid(1,1)
    hy=poisson%hgrid(2,2)
    hz=poisson%hgrid(3,3)
    z0=cell(3)/2.0
    t1=1.d0/(2.d0*c**4*cell(1)**2*cell(2)**2*pi)
    do iz=1,poisson%ngpz
    do iy=1,poisson%ngpy
    do ix=1,poisson%ngpx
        z=(iz-1)*hz
        y=(iy-1)*hy
        x=(ix-1)*hx
        t2=exp(-(z-z0)**2/c**2)
        t3=(2.d0*b*c**4*cell(1)**2*pi**2*cos(b*sin((2.d0*pi*y)/cell(2)))* &
            sin((2.d0*pi*y)/cell(2))*sin(a*sin((2.d0*pi*x)/cell(1)))+ &
            (2.d0*a*c**4*cell(2)**2*pi**2*cos(a*sin((2.d0*pi*x)/cell(1)))* &
            sin((2.d0*pi*x)/cell(1))+(cell(2)**2*(cell(1)**2*(c**2-2.d0*(z-z0)**2)+ &
            a**2*c**4*pi**2*(1.d0+cos((4.d0*pi*x)/cell(1))))+ &
            b**2*c**4*cell(1)**2*pi**2*(1.d0+cos((4.d0*pi*y)/cell(2))))* &
            sin(a*sin((2.d0*pi*x)/cell(1))))*sin(b*sin((2.d0*pi*y)/cell(2))))
        poisson%rho(ix,iy,iz)=t1*t2*t3
    enddo
    enddo
    enddo
end subroutine get_density
!*****************************************************************************************
