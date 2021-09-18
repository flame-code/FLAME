!*****************************************************************************************
subroutine get_psolver_kspace_exprnscreening(ngpx,ngpy,ngpz,hgrid,rho,sf,npow,pot)
    use mod_greenf_kspace, only: typ_greenf_kspace
    implicit none
    include 'fftw3.f'
    integer, intent(in):: ngpx, ngpy, ngpz, npow
    real(8), intent(in):: rho(ngpx,ngpy,ngpz), hgrid(3,3), sf
    real(8), intent(out):: pot(ngpx,ngpy,ngpz)
    !local variables
    integer:: igpx, igpy, igpz
    integer:: igpxt, igpyt, igpzt
    integer:: ngpx_ext, ngpy_ext, ngpz_ext, nadd
    real(8):: tt1, pi, fourpi, hmin, daw, sf2
    real(8):: akx, aky, akz, aknorm2
    real(8), allocatable:: rho_ext(:,:,:), pot_ext(:,:,:)
    integer(8):: planf, planb
    type(typ_greenf_kspace):: greenf_kspace
    sf2=sf**2
    pi=4.d0*atan(1.d0)
    fourpi=4.d0*pi
    hmin=min(hgrid(1,1),hgrid(2,2),hgrid(3,3))
    nadd=floor((-log(1.d-4))**(1.d0/real(npow,kind=8))/(sf*hmin))
    !if(mod(nadd,2)/=0) nadd=nadd+1
    write(*,*) 'nadd= ',nadd
    ngpx_ext=ngpx+2*nadd
    ngpy_ext=ngpy+2*nadd
    ngpz_ext=ngpz+2*nadd
    if(mod(ngpx_ext,2)/=0) ngpx_ext=ngpx_ext+1
    if(mod(ngpy_ext,2)/=0) ngpy_ext=ngpy_ext+1
    if(mod(ngpz_ext,2)/=0) ngpz_ext=ngpz_ext+1
    call greenf_kspace%init_greenf_kspace(npow,'numeric')
    allocate(pot_ext(ngpx_ext+2,ngpy_ext,ngpz_ext),rho_ext(ngpx_ext,ngpy_ext,ngpz_ext))
    call dfftw_plan_dft_r2c_3d(planf,ngpx_ext,ngpy_ext,ngpz_ext,rho_ext(1,1,1),pot_ext(1,1,1),fftw_estimate)
    call dfftw_plan_dft_c2r_3d(planb,ngpx_ext,ngpy_ext,ngpz_ext,pot_ext(1,1,1),pot_ext(1,1,1),fftw_patient)
    rho_ext=0.d0
    do igpz=1,ngpz
    do igpy=1,ngpy
    do igpx=1,ngpx
        rho_ext(nadd+igpx,nadd+igpy,nadd+igpz)=rho(igpx,igpy,igpz)
    enddo
    enddo
    enddo
    call dfftw_execute(planf)
    do igpz=1,ngpz_ext
    do igpy=1,ngpy_ext
    do igpx=1,ngpx_ext
        igpxt=igpx
        igpyt=igpy
        igpzt=igpz
        if(igpy>ngpy_ext/2) igpyt=ngpy_ext-igpy+2
        if(igpz>ngpz_ext/2) igpzt=ngpz_ext-igpz+2
        akx=2.d0*pi*((igpxt-1)/2)/(real(ngpx_ext,kind=8)*hgrid(1,1))
        aky=2.d0*pi*( igpyt-1   )/(real(ngpy_ext,kind=8)*hgrid(2,2))
        akz=2.d0*pi*( igpzt-1   )/(real(ngpz_ext,kind=8)*hgrid(3,3))
        aknorm2=akx**2+aky**2+akz**2
        call greenf_kspace%get_greenf_kspace_single(sqrt(aknorm2),sf,tt1)
        pot_ext(igpx,igpy,igpz)=pot_ext(igpx,igpy,igpz)*tt1
    enddo
    enddo
    enddo
    call greenf_kspace%fini_greenf_kspace()
    call dfftw_execute(planb)
    tt1=1.d0/(real(ngpx_ext,kind=8)*real(ngpy_ext,kind=8)*real(ngpz_ext,kind=8))
    do igpz=1,ngpz
    do igpy=1,ngpy
    do igpx=1,ngpx
        pot(igpx,igpy,igpz)=pot_ext(nadd+igpx,nadd+igpy,nadd+igpz)*tt1
    enddo
    enddo
    enddo
    !-----------------------------------------------------------------
    call dfftw_destroy_plan(planf)
    call dfftw_destroy_plan(planb)
    deallocate(pot_ext,rho_ext)
end subroutine get_psolver_kspace_exprnscreening
!*****************************************************************************************
subroutine get_psolver_kspace_gaussscreening(ngpx,ngpy,ngpz,hgrid,rho,sf,pot)
    implicit none
    include 'fftw3.f'
    integer, intent(in):: ngpx, ngpy, ngpz
    real(8), intent(in):: rho(ngpx,ngpy,ngpz), hgrid(3,3), sf
    real(8), intent(out):: pot(ngpx,ngpy,ngpz)
    !local variables
    integer:: igpx, igpy, igpz
    integer:: igpxt, igpyt, igpzt
    integer:: ngpx_ext, ngpy_ext, ngpz_ext, nadd
    real(8):: tt1, pi, fourpi, hmin, daw, sf2
    real(8):: akx, aky, akz, aknorm2
    real(8), allocatable:: rho_ext(:,:,:), pot_ext(:,:,:)
    integer(8):: planf, planb
    sf2=sf**2
    pi=4.d0*atan(1.d0)
    fourpi=4.d0*pi
    hmin=min(hgrid(1,1),hgrid(2,2),hgrid(3,3))
    nadd=floor(sqrt(-log(1.d-10)/sf2)/hmin)
    if(mod(nadd,2)/=0) nadd=nadd+1
    write(*,*) 'nadd= ',nadd
    ngpx_ext=ngpx+nadd
    ngpy_ext=ngpy+nadd
    ngpz_ext=ngpz+nadd
    allocate(pot_ext(ngpx_ext+2,ngpy_ext,ngpz_ext),rho_ext(ngpx_ext,ngpy_ext,ngpz_ext))
    call dfftw_plan_dft_r2c_3d(planf,ngpx_ext,ngpy_ext,ngpz_ext,rho_ext(1,1,1),pot_ext(1,1,1),fftw_estimate)
    call dfftw_plan_dft_c2r_3d(planb,ngpx_ext,ngpy_ext,ngpz_ext,pot_ext(1,1,1),pot_ext(1,1,1),fftw_patient)
    rho_ext=0.d0
    do igpz=1,ngpz
    do igpy=1,ngpy
    do igpx=1,ngpx
        rho_ext(igpx,igpy,igpz)=rho(igpx,igpy,igpz)
    enddo
    enddo
    enddo
    call dfftw_execute(planf)
    do igpz=1,ngpz_ext
    do igpy=1,ngpy_ext
    do igpx=1,ngpx_ext
        igpxt=igpx
        igpyt=igpy
        igpzt=igpz
        if(igpy>ngpy_ext/2) igpyt=ngpy_ext-igpy+2
        if(igpz>ngpz_ext/2) igpzt=ngpz_ext-igpz+2
        akx=2.d0*pi*((igpxt-1)/2)/(real(ngpx_ext,kind=8)*hgrid(1,1))
        aky=2.d0*pi*( igpyt-1   )/(real(ngpy_ext,kind=8)*hgrid(2,2))
        akz=2.d0*pi*( igpzt-1   )/(real(ngpz_ext,kind=8)*hgrid(3,3))
        aknorm2=akx**2+aky**2+akz**2
        if(abs(aknorm2)<1.d-20) then
            tt1=2.d0*pi/sf2
        else
            tt1=fourpi*daw(sqrt(aknorm2/sf2)*0.5d0)/(sqrt(aknorm2*sf2))
        endif
        pot_ext(igpx,igpy,igpz)=pot_ext(igpx,igpy,igpz)*tt1
    enddo
    enddo
    enddo
    call dfftw_execute(planb)
    tt1=1.d0/(real(ngpx_ext,kind=8)*real(ngpy_ext,kind=8)*real(ngpz_ext,kind=8))
    do igpz=1,ngpz
    do igpy=1,ngpy
    do igpx=1,ngpx
        pot(igpx,igpy,igpz)=pot_ext(igpx,igpy,igpz)*tt1
    enddo
    enddo
    enddo
    !-----------------------------------------------------------------
    call dfftw_destroy_plan(planf)
    call dfftw_destroy_plan(planb)
    deallocate(pot_ext,rho_ext)
end subroutine get_psolver_kspace_gaussscreening
!*****************************************************************************************
subroutine get_psolver_kspace_expscreening(ngpx,ngpy,ngpz,hgrid,rho,sf,pot)
    implicit none
    include 'fftw3.f'
    integer, intent(in):: ngpx, ngpy, ngpz
    real(8), intent(in):: rho(ngpx,ngpy,ngpz), hgrid(3,3), sf
    real(8), intent(out):: pot(ngpx,ngpy,ngpz)
    !local variables
    integer:: igpx, igpy, igpz
    integer:: igpxt, igpyt, igpzt
    integer:: ngpx_ext, ngpy_ext, ngpz_ext, nadd
    real(8):: tt1, pi, fourpi, hmin
    real(8):: akx, aky, akz, aknorm2
    real(8), allocatable:: rho_ext(:,:,:), pot_ext(:,:,:)
    integer(8):: planf, planb
    pi=4.d0*atan(1.d0)
    fourpi=4.d0*pi
    hmin=min(hgrid(1,1),hgrid(2,2),hgrid(3,3))
    nadd=floor(-log(1.d-10)/(sf*hmin))
    if(mod(nadd,2)/=0) nadd=nadd+1
    write(*,*) 'nadd= ',nadd
    ngpx_ext=ngpx+nadd
    ngpy_ext=ngpy+nadd
    ngpz_ext=ngpz+nadd
    allocate(pot_ext(ngpx_ext+2,ngpy_ext,ngpz_ext),rho_ext(ngpx_ext,ngpy_ext,ngpz_ext))
    call dfftw_plan_dft_r2c_3d(planf,ngpx_ext,ngpy_ext,ngpz_ext,rho_ext(1,1,1),pot_ext(1,1,1),fftw_estimate)
    call dfftw_plan_dft_c2r_3d(planb,ngpx_ext,ngpy_ext,ngpz_ext,pot_ext(1,1,1),pot_ext(1,1,1),fftw_patient)
    rho_ext=0.d0
    do igpz=1,ngpz
    do igpy=1,ngpy
    do igpx=1,ngpx
        rho_ext(igpx,igpy,igpz)=rho(igpx,igpy,igpz)
    enddo
    enddo
    enddo
    call dfftw_execute(planf)
    do igpz=1,ngpz_ext
    do igpy=1,ngpy_ext
    do igpx=1,ngpx_ext
        igpxt=igpx
        igpyt=igpy
        igpzt=igpz
        if(igpy>ngpy_ext/2) igpyt=ngpy_ext-igpy+2
        if(igpz>ngpz_ext/2) igpzt=ngpz_ext-igpz+2
        akx=2.d0*pi*((igpxt-1)/2)/(real(ngpx_ext,kind=8)*hgrid(1,1))
        aky=2.d0*pi*( igpyt-1   )/(real(ngpy_ext,kind=8)*hgrid(2,2))
        akz=2.d0*pi*( igpzt-1   )/(real(ngpz_ext,kind=8)*hgrid(3,3))
        aknorm2=akx**2+aky**2+akz**2
        tt1=fourpi/(aknorm2+sf**2)
        pot_ext(igpx,igpy,igpz)=pot_ext(igpx,igpy,igpz)*tt1
        !if(abs(pot_kspace(igpx,igpy,igpz))>1.d-10) then
        !write(*,'(3i4,f10.1)') igpx,igpy,igpz,pot_kspace(igpx,igpy,igpz)
        !endif
    enddo
    enddo
    enddo
    call dfftw_execute(planb)
    tt1=1.d0/(real(ngpx_ext,kind=8)*real(ngpy_ext,kind=8)*real(ngpz_ext,kind=8))
    do igpz=1,ngpz
    do igpy=1,ngpy
    do igpx=1,ngpx
        pot(igpx,igpy,igpz)=pot_ext(igpx,igpy,igpz)*tt1
    enddo
    enddo
    enddo
    call dfftw_destroy_plan(planf)
    call dfftw_destroy_plan(planb)
end subroutine get_psolver_kspace_expscreening
!*****************************************************************************************
