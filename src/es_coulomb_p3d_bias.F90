subroutine bias_potener_forces(parini,poisson,atoms,epotplane)
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms
    use mod_potential, only: potential 
    use mod_parini, only: typ_parini
    use dynamic_memory
    use yaml_output
    implicit none
    type(typ_poisson), intent(inout):: poisson
    type(typ_atoms), intent(inout):: atoms
    type(typ_parini), intent(in):: parini
    !local variables
    real(8):: epotlong, epotplane 
    real(8):: time(10)
    !The following dummy varables definition to be deleted later.
    !real(8):: totrho
    real(8):: beta, pi, charge,c,charge0,E,tmp
    integer:: igpx, igpy, igpz, iat
    integer:: ix, iy, iz, jx, jy, jz, kx, ky, kz
    integer:: npl, npu, nlayer, ngpx, ngpy, ngpz, nbgpz
    real(8),allocatable :: pots_layer(:,:,:,:)
    real(8):: vl, vu, A, d, rl, ru, dipole_correction, dipole,dv
    real(8):: tt, tt1, tt2

    pi=4.d0*atan(1.d0)
    epotplane=0.d0
    beta = poisson%beta*(-poisson%ngpx*poisson%ngpy)
    ngpx= poisson%ngpx
    ngpy= poisson%ngpy
    nbgpz=int(poisson%rgcut/poisson%hz)+2
    if(trim(parini%bias_type)=='p3dbias') then
        vl=parini%vl_ewald
        vu=parini%vu_ewald+parini%vu_ac_ewald*sin(parini%frequency_ewald*parini%time_dynamics)
        vu=vu 
        dv= (vu-vl)
        d= poisson%cell(3)
        !write(*,*)'--------------------------------------------------------'
!        write(*,*)"distance between planes",d
        A= poisson%ngpx*poisson%ngpy*poisson%hx*poisson%hy
        c= A/(4.d0*pi*d)
        dipole = beta*(poisson%hx*poisson%hy)
        charge0= -dipole/(2*pi*d)
        charge = -dipole/(2*pi*d)+c*(dv)
        call yaml_mapping_open('p3dbias info') !,flow=.true.)
        call yaml_map('dipole moment',dipole/(2*pi))
        call yaml_map('real pot (lower)',beta/(poisson%ngpx*poisson%ngpy)+vl)
        call yaml_map('real pot (upper)',-beta/(poisson%ngpx*poisson%ngpy)+vu)
        call yaml_map('charge on upper plane',charge)
        call yaml_map('C_0',c)
        call yaml_map('K',charge/dv/c)
        call yaml_mapping_close()
        !write(*,*)'dipole = ', dipole/(2*pi)
        !write(*,*)'real pot = ', beta/(poisson%ngpx*poisson%ngpy)+vl,&
        !                        -beta/(poisson%ngpx*poisson%ngpy)+vu
        !write(*,*)'charge on upper  plate  ', charge
        !write(*,*)"C_0 = ",c, " K =" , charge/dv/c


        !********************************************************************
        ! Esperesso energy 
        !dipole_correction = 3/(4*pi)*dipole**2/(poisson%cell(3)*poisson%cell(2)*poisson%cell(1))
        !dipole_correction = 0.d0
        !dipole_correction =dipole_correction +charge0*(vu-vl)!0.5*c*(vu-vl)**2
        !********************************************************************
        dipole_correction = 0.d0
        dipole_correction =dipole_correction -0.5*charge0*(dv)!+0.5*c*(-dv)**2
        atoms%ebattery=-charge0*dv
        poisson%npu=poisson%ngpz-nbgpz
        poisson%npl=1+nbgpz  
        npl=poisson%npl
        npu=poisson%npu

        allocate(poisson%pots(1:poisson%ngpx+2,1:poisson%ngpy,npl:npu))
        poisson%pots=0.d0
        nlayer=1
        if (parini%cal_charge) then 
            nlayer=5
            allocate(pots_layer(1:poisson%ngpx,1:poisson%ngpy,1:2,1:nlayer))
        endif
        pots_layer = 0.d0
        if(.not. (trim(potential)=='ann')) then
            call erfc_surface_zero(parini,atoms,poisson,nlayer)
            if (parini%cal_charge) then
                pots_layer(1:ngpx,1:ngpy,1,:)=poisson%pots(1:ngpx,1:ngpy,npl:npl+nlayer-1)
                pots_layer(1:ngpx,1:ngpy,2,:)=poisson%pots(1:ngpx,1:ngpy,npu:npu-(nlayer-1):-1)
                poisson%pots(:,:,npl+1:npl+nlayer-1)  =0.d0
                poisson%pots(:,:,npu-(nlayer-1):npu-1)=0.d0
            endif
        endif
        call sollaplaceq(poisson,poisson%hz,poisson%cell,vl,vu)
        call calculate_force_ener_plane(atoms,poisson,epotplane,nbgpz)
        

        if (parini%cal_charge) then 
            call surface_charge(parini,poisson,pots_layer,vl,vu)
            deallocate(pots_layer)
        endif
     !           open(unit=55, file="pots.txt" )
     !              do iy=1,poisson%ngpy!min(9,poisson%ngpy)
     !              do ix=1,poisson%ngpx!min(9,poisson%ngpx)
     !              do iz=poisson%npl,poisson%npu
     !                  !write(55,*)  (iz-1-nbgpz)*poisson%hz, -poisson%pots(ix,iy,iz)+&
     !                  !             (-2*beta/(poisson%ngpx*poisson%ngpy)+vu-vl)/d*(iz-1-nbgpz)*poisson%hz+beta/(poisson%ngpx*poisson%ngpy)+vl
     !                  write(55,'(4es20.8)')  (iz-1-nbgpz)*poisson%hz, poisson%pots(ix,iy,iz),&
     !                               (vu-vl)/d*(iz-1-nbgpz)*poisson%hz+vl,&
     !                               (-2*beta/(poisson%ngpx*poisson%ngpy)+vu-vl)/d*(iz-1-nbgpz)*poisson%hz+beta/(poisson%ngpx*poisson%ngpy)+vl

     !              enddo 
     !                  write(55,*)   
     !              enddo 
     !              enddo 
     !           close(55)
     !           open(unit=56, file="efield.txt" )
     !              do iy=1,poisson%ngpy,10!min(9,poisson%ngpy)
     !              do ix=1,poisson%ngpx,10!min(9,poisson%ngpx)
     !              do iz=poisson%npl+1,poisson%npu-1,2
     !                  !write(55,*)  (iz-1-nbgpz)*poisson%hz, -poisson%pots(ix,iy,iz)+&
     !                  !             (-2*beta/(poisson%ngpx*poisson%ngpy)+vu-vl)/d*(iz-1-nbgpz)*poisson%hz+beta/(poisson%ngpx*poisson%ngpy)+vl
     !                  tt= -(poisson%pots(ix,iy,iz-1)-poisson%pots(ix,iy,iz+1))/(2*poisson%hz)/0.529177210d0*27.211385d0
     !                  tt1= (vu-vl)/d/0.529177210d0*27.211385d0
     !                  tt2= (-2*beta/(poisson%ngpx*poisson%ngpy)+vu-vl)/d/0.529177210d0*27.211385d0
     !                  write(56,'(4es20.8)')  (iz-1-nbgpz)*poisson%hz*0.529177210d0, tt,tt1,tt2
     !              enddo 
     !                  write(56,*)   
     !              enddo 
     !              enddo 
     !           close(56)
     !           open(unit=66, file="pot_up.txt" )
     !              do iy=1,poisson%ngpy
     !              do ix=1,poisson%ngpx
     !                  iz=poisson%npu
     !                  write(66,*)  ix*poisson%hx , iy*poisson%hy, poisson%pots(ix,iy,iz)-(-beta/(poisson%ngpx*poisson%ngpy)+vu)
     !              enddo 
     !                  write(66,*)  
     !              enddo 
     !           close(66)
        epotplane = epotplane+dipole_correction
        deallocate(poisson%pots)
    end if
end subroutine bias_potener_forces
!*****************************************************************************************
subroutine erfc_surface_zero(parini,atoms,poisson,nlayer)
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_electrostatics, only: typ_linked_lists
    use mod_parini, only: typ_parini
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_poisson), intent(inout):: poisson
    type(typ_atoms), intent(inout):: atoms
    type(typ_linked_lists):: linked_lists

    !local variables
    real(8):: dx, dy, dz, sclinv, r, rsq, xiat, yiat, ziat, alphainv, twosqrtinv
    real(8):: t, tt, tt1, tt2, tt3, ttt,xgps,ygps,zgps,dzysq,dzsq
    real(8):: fx,fy, fz, pi, hspinv, rhspinv, rinv, qiatp, qiatpjatp, spf, spfd
    real(8):: hxinv, hyinv, hzinv ,zgpu,zgpl,hx,hy,x,y
    real(8):: dzlsq,dzusq,rcutsq
    real(8):: cell(3)
    real(8):: dnlayer, hgxinv , hgyinv, hgzinv 
    real(8):: xat, yat, zat

    integer:: nimat !number of image atoms.
    integer:: nlayer, igpx,igpy,igpz,mx,my,mz, mlimnlayer
    integer:: nkat, njat, jat, kat, jl, npl, npu
    integer:: ix, iy, iz, jx, jy, jz, kx, ky, kz
    integer:: nsclx,nscly,nx,ny,nz , jat1, jat2, iat, isp
    integer:: ip, il, iatox, iatoy, iatoz
    integer:: ngpx, ngpy, ngpz, nbgpx, nbgpy, nbgpz, nbgpz_poisson
    real(8), allocatable:: ratp(:,:,:),qatp(:,:)
    integer, allocatable::  mboundg(:,:,:),mboundgy(:,:)

    alphainv=1.d0/poisson%alpha;twosqrtinv=1.d0!/sqrt(2.d0)
    poisson%pots=0.d0
    linked_lists%rcut = poisson%linked_lists%rcut/sqrt(2.d0)
    rcutsq= poisson%linked_lists%rcut**2

    nbgpz_poisson=int(poisson%rgcut/poisson%hz)+2
    npl=poisson%npl
    npu=poisson%npu
    zgpl=(npl-1-nbgpz_poisson)*poisson%hz
    zgpu=(npu-1-nbgpz_poisson)*poisson%hz

    ngpx=poisson%ngpx
    ngpy=poisson%ngpy
    ngpz=poisson%ngpz

    nbgpx = int(linked_lists%rcut/poisson%hx)+1
    nbgpy = int(linked_lists%rcut/poisson%hy)+1
    nbgpz = int(linked_lists%rcut/poisson%hz)+1
    
    allocate(mboundg(1:2,-nbgpy:nbgpy,-nbgpz:nbgpz))
    allocate(mboundgy(1:2,-nbgpz:nbgpz))
    call determine_limitsphere(poisson,mboundg,mboundgy,nbgpx,nbgpy,nbgpz)

    call linkedlists_init(parini,atoms,cell,linked_lists)
    mx = linked_lists%mx
    my = linked_lists%my
    mz = linked_lists%mz
    hzinv=real(mz,8)/cell(3)

    dnlayer=(nlayer-1)*poisson%hz
    mlimnlayer=floor(dnlayer*hzinv)+1
    hgxinv=1.d0/poisson%hx
    hgyinv=1.d0/poisson%hy
    hgzinv=1.d0/poisson%hz

    call update_ratp(linked_lists%typ_atoms)
    do kz=linked_lists%mz,max(linked_lists%mz-mlimnlayer-linked_lists%mlimnb,1),-1
    do ky=1,linked_lists%my
    do kx=1,linked_lists%mx
    ip=linked_lists%prime(kx,ky,kz)
    il=linked_lists%last(kx,ky,kz)
    do  iat=ip,il
        if (zgpu-linked_lists%ratp(3,iat)>linked_lists%rcut+dnlayer) cycle
        iatox=nint(linked_lists%ratp(1,iat)*hgxinv)+1
        iatoy=nint(linked_lists%ratp(2,iat)*hgyinv)+1
        iatoz=nint(linked_lists%ratp(3,iat)*hgzinv)+1+nbgpz_poisson
        xat=linked_lists%ratp(1,iat)-(iatox-1)*poisson%hx
        yat=linked_lists%ratp(2,iat)-(iatoy-1)*poisson%hy
        zat=linked_lists%ratp(3,iat)-(iatoz-1-nbgpz_poisson)*poisson%hz
        do iz=-nbgpz,nbgpz
            jz=iatoz+iz
            if (.not. (jz<=npu .and. jz>npu-nlayer)) cycle
            dzsq= (iz*poisson%hz-zat)**2
            do iy=mboundgy(1,iz),mboundgy(2,iz)
                dzysq= (iy*poisson%hy-yat)**2+dzsq
                jy=modulo(iatoy+iy-1,ngpy)+1
                do ix=mboundg(1,iy,iz),mboundg(2,iy,iz)
                    dx= ix*poisson%hx-xat
                    rsq= dx*dx+dzysq
                    if(rsq.lt.rcutsq )then
                        jx=modulo(iatox+ix-1,ngpx)+1
                        r=sqrt(rsq)
                        t=erfc(r*alphainv)
                        tt=linked_lists%qat(iat)*t/r
                        poisson%pots(jx,jy,jz)=poisson%pots(jx,jy,jz)+tt
                    endif
                enddo
            enddo
        enddo
    enddo
    enddo
    enddo
    enddo

    do kz=1,min(mlimnlayer+linked_lists%mlimnb+1,linked_lists%mz)
    do ky=1,linked_lists%my
    do kx=1,linked_lists%mx
    ip=linked_lists%prime(kx,ky,kz)
    il=linked_lists%last(kx,ky,kz)
    do  iat=ip,il
        if (-zgpl+linked_lists%ratp(3,iat)>linked_lists%rcut+dnlayer) cycle
        iatox=nint(linked_lists%ratp(1,iat)*hgxinv)+1
        iatoy=nint(linked_lists%ratp(2,iat)*hgyinv)+1
        iatoz=nint(linked_lists%ratp(3,iat)*hgzinv)+1+nbgpz_poisson
        xat=linked_lists%ratp(1,iat)-(iatox-1)*poisson%hx
        yat=linked_lists%ratp(2,iat)-(iatoy-1)*poisson%hy
        zat=linked_lists%ratp(3,iat)-(iatoz-1-nbgpz_poisson)*poisson%hz
        do iz=-nbgpz,nbgpz
            jz=iatoz+iz
            if (.not. (jz>=npl .and. jz<npl+nlayer .and. jz<=npu-nlayer)) cycle
            dzsq= (iz*poisson%hz-zat)**2
            do iy=mboundgy(1,iz),mboundgy(2,iz)
                jy=modulo(iatoy+iy-1,ngpy)+1
                dzysq= (iy*poisson%hy-yat)**2+dzsq
                do ix=mboundg(1,iy,iz),mboundg(2,iy,iz)
                    dx= ix*poisson%hx-xat
                    rsq= dx*dx+dzysq
                    if(rsq.lt.rcutsq )then
                        jx=modulo(iatox+ix-1,ngpx)+1
                        r=sqrt(rsq)
                        t=erfc(r*alphainv*twosqrtinv)
                        tt=linked_lists%qat(iat)*t/r
                        poisson%pots(jx,jy,jz)=poisson%pots(jx,jy,jz)+tt
                    endif
                enddo
            enddo
        enddo
    enddo
    enddo
    enddo
    enddo
    call linkedlists_final(linked_lists)
    deallocate(mboundg)
    deallocate(mboundgy)
end subroutine erfc_surface_zero
!*****************************************************************************************
subroutine sollaplaceq(poisson,hz,cell,vl,vu)
    use mod_electrostatics, only: typ_poisson
    implicit none
    include 'fftw3.f'
    type(typ_poisson), intent(inout):: poisson
    !local variables
    real(8):: vl, vu , zlmzu , sinhzlmzu, zlmzuinv
    real(8):: cell(3)
    real(8):: hz , vlmvu, vlzumvuzl 
    real(8):: t,tt,ttt, ghz,ciz,potl,potu
    real(8):: pi, fourpi, fourpisq, gsq, gsqx, gsqy, g, hzsq
    real(8):: fourpisqcellxinvsq, fourpisqcellyinvsq,valuengpxyinv
    real(8):: hzlu,hzliz,hzuiz, zz 
    real(8):: tel,teu,telu
    real(8) :: temp_exp_2 ,temp_exp_l , temp_exp_zl1, temp_exp_zl2, temp_exp_z  , temp_exp_zl  
    integer::ix,iy,iz,ixt,iyt,npu,npl,izz
    integer(8), allocatable:: plan_bs(:),plan_fs(:)
    npl=poisson%npl
    npu=poisson%npu
    pi=4.d0*atan(1.d0)
    hzsq=hz**2
    fourpi=4.d0*pi 
    fourpisq=fourpi*pi 
    fourpisqcellxinvsq=fourpisq/cell(1)**2 
    fourpisqcellyinvsq=fourpisq/cell(2)**2 
    allocate(plan_bs(npl:npu))
    allocate(plan_fs(1:2))
    do iz=npl,npu
        call dfftw_plan_dft_c2r_2d(plan_bs(iz),poisson%ngpx, &
            poisson%ngpy,poisson%pots(1,1,iz),poisson%pots(1,1,iz),fftw_estimate)
    enddo
    call dfftw_plan_dft_r2c_2d(plan_fs(1),poisson%ngpx, &
            poisson%ngpy,poisson%pots(1,1,npl), &
            poisson%pots(1,1,npl),fftw_estimate)
    call dfftw_plan_dft_r2c_2d(plan_fs(2),poisson%ngpx, &
            poisson%ngpy,poisson%pots(1,1,npu), &
            poisson%pots(1,1,npu),fftw_estimate)

    do iy=1,poisson%ngpy
    do ix=1,poisson%ngpx
        poisson%pots(ix,iy,npl)=-(poisson%pots(ix,iy,npl)+poisson%pot(ix,iy,npl))+vl
        poisson%pots(ix,iy,npu)=-(poisson%pots(ix,iy,npu)+poisson%pot(ix,iy,npu))+vu    
    enddo 
    enddo 
    call dfftw_execute(plan_fs(1))
    call dfftw_execute(plan_fs(2))   

   ! k=0 , l=0
     potl=poisson%pots(1,1,npl)
     potu=poisson%pots(1,1,npu)

    !-----------------------------------------
    ix=1;iy=1
    do iz=npl+1,npu-1
        ciz=((potl-potu)*iz-(potl*npu-potu*npl))/real((npl-npu),8)
        poisson%pots(ix,iy,iz)=ciz
    enddo 
    !-----------------------------------------
    do iz=npl+1,npu-1
        izz = iz - npl
        zz = izz*hz
        !-----------------------------------------
        ix=poisson%ngpx+1;iy=1
        gsq=fourpisqcellxinvsq*((ix-1)/2)**2 + fourpisqcellyinvsq*(iy-1)**2
        g=sqrt(gsq)
        !****************************************************  
        temp_exp_l  = exp(-g*cell(3))
        temp_exp_zl1= exp(g*(zz-cell(3)))
        temp_exp_z  = exp(-g*(zz))
        temp_exp_2  = temp_exp_l**2 
        temp_exp_zl2= temp_exp_zl1*temp_exp_l     !exp(g*(zz-2.d0*cell(3)))
        temp_exp_zl = temp_exp_z  *temp_exp_l     ! exp(-g*(zz+cell(3)))
        t = temp_exp_zl1 - temp_exp_zl
        tt= temp_exp_z - temp_exp_zl2
        ttt = 1.d0/(1.d0-temp_exp_2)

        ciz= (poisson%pots(ix,iy,npu)*t+poisson%pots(ix,iy,npl)*tt)*ttt
        poisson%pots(ix,iy,iz)=ciz
        !-----------------------------------------
        ix=1;iy=poisson%ngpy/2+1
        gsq=fourpisqcellyinvsq*(iy-1)**2
        g=sqrt(gsq)
        !****************************************************  
        temp_exp_l  = exp(-g*cell(3))
        temp_exp_zl1= exp(g*(zz-cell(3)))
        temp_exp_z  = exp(-g*(zz))
        temp_exp_2  = temp_exp_l**2 
        temp_exp_zl2= temp_exp_zl1*temp_exp_l     !exp(g*(zz-2.d0*cell(3)))
        temp_exp_zl = temp_exp_z  *temp_exp_l     ! exp(-g*(zz+cell(3)))
        t = temp_exp_zl1 - temp_exp_zl
        tt= temp_exp_z - temp_exp_zl2
        ttt = 1.d0/(1.d0-temp_exp_2)

        ciz= (poisson%pots(ix,iy,npu)*t+poisson%pots(ix,iy,npl)*tt)*ttt
        poisson%pots(ix,iy,iz)=ciz
        !-----------------------------------------
        ix=poisson%ngpx+1;iy=poisson%ngpy/2+1
        gsq=fourpisqcellxinvsq*((ix-1)/2)**2 + fourpisqcellyinvsq*(iy-1)**2
        g=sqrt(gsq)
        !****************************************************  
        temp_exp_l  = exp(-g*cell(3))
        temp_exp_zl1= exp(g*(zz-cell(3)))
        temp_exp_z  = exp(-g*(zz))
        temp_exp_2  = temp_exp_l**2 
        temp_exp_zl2= temp_exp_zl1*temp_exp_l     !exp(g*(zz-2.d0*cell(3)))
        temp_exp_zl = temp_exp_z  *temp_exp_l     ! exp(-g*(zz+cell(3)))
        t = temp_exp_zl1 - temp_exp_zl
        tt= temp_exp_z - temp_exp_zl2
        ttt = 1.d0/(1.d0-temp_exp_2)

        ciz= (poisson%pots(ix,iy,npu)*t+poisson%pots(ix,iy,npl)*tt)*ttt
        !****************************************************  
        poisson%pots(ix,iy,iz)=ciz
        !-----------------------------------------
        ix=1
        do iy=2,poisson%ngpy/2
            gsq=fourpisqcellyinvsq*(iy-1)**2
            g=sqrt(gsq)
            !****************************************************  
            temp_exp_l  = exp(-g*cell(3))
            temp_exp_zl1= exp(g*(zz-cell(3)))
            temp_exp_z  = exp(-g*(zz))
            temp_exp_2  = temp_exp_l**2 
            temp_exp_zl2= temp_exp_zl1*temp_exp_l     !exp(g*(zz-2.d0*cell(3)))
            temp_exp_zl = temp_exp_z  *temp_exp_l     ! exp(-g*(zz+cell(3)))
            t = temp_exp_zl1 - temp_exp_zl
            tt= temp_exp_z - temp_exp_zl2
            ttt = 1.d0/(1.d0-temp_exp_2)

            ciz= (poisson%pots(ix,iy,npu)*t+poisson%pots(ix,iy,npl)*tt)*ttt
            poisson%pots(ix,iy,iz)=ciz
            poisson%pots(ix,poisson%ngpy-iy+2,iz)=ciz

            ixt=ix+1
            ciz= (poisson%pots(ixt,iy,npu)*t+poisson%pots(ixt,iy,npl)*tt)*ttt
            poisson%pots(ixt,iy,iz)=ciz
            poisson%pots(ixt,poisson%ngpy-iy+2,iz)=-ciz
        enddo
        !-----------------------------------------
        ix=poisson%ngpx+1
        gsqx=fourpisqcellxinvsq*(poisson%ngpx/2)**2
        do iy=2,poisson%ngpy/2
            gsq=gsqx+fourpisqcellyinvsq*(iy-1)**2
            g=sqrt(gsq)
            !****************************************************  
            temp_exp_l  = exp(-g*cell(3))
            temp_exp_zl1= exp(g*(zz-cell(3)))
            temp_exp_z  = exp(-g*(zz))
            temp_exp_2  = temp_exp_l**2 
            temp_exp_zl2= temp_exp_zl1*temp_exp_l     !exp(g*(zz-2.d0*cell(3)))
            temp_exp_zl = temp_exp_z  *temp_exp_l     ! exp(-g*(zz+cell(3)))
            t = temp_exp_zl1 - temp_exp_zl
            tt= temp_exp_z - temp_exp_zl2
            ttt = 1.d0/(1.d0-temp_exp_2)

            ciz= (poisson%pots(ix,iy,npu)*t+poisson%pots(ix,iy,npl)*tt)*ttt
            poisson%pots(ix,iy,iz)=ciz
            poisson%pots(ix,poisson%ngpy-iy+2,iz)=ciz

            ixt=ix+1
            ciz= (poisson%pots(ixt,iy,npu)*t+poisson%pots(ixt,iy,npl)*tt)*ttt
            poisson%pots(ixt,iy,iz)=ciz
            poisson%pots(ixt,poisson%ngpy-iy+2,iz)=-ciz
        enddo
        !-----------------------------------------
        iy=1
        gsqy=fourpisqcellyinvsq*(iy-1)**2
        do ix=3,poisson%ngpx,2
            gsq=gsqy+fourpisqcellxinvsq*((ix-1)/2)**2
            g=sqrt(gsq)
            !****************************************************  
            temp_exp_l  = exp(-g*cell(3))
            temp_exp_zl1= exp(g*(zz-cell(3)))
            temp_exp_z  = exp(-g*(zz))
            temp_exp_2  = temp_exp_l**2 
            temp_exp_zl2= temp_exp_zl1*temp_exp_l     !exp(g*(zz-2.d0*cell(3)))
            temp_exp_zl = temp_exp_z  *temp_exp_l     ! exp(-g*(zz+cell(3)))
            t = temp_exp_zl1 - temp_exp_zl
            tt= temp_exp_z - temp_exp_zl2
            ttt = 1.d0/(1.d0-temp_exp_2)

            ciz= (poisson%pots(ix,iy,npu)*t+poisson%pots(ix,iy,npl)*tt)*ttt
            poisson%pots(ix,iy,iz)=ciz

            ixt=ix+1
            ciz= (poisson%pots(ixt,iy,npu)*t+poisson%pots(ixt,iy,npl)*tt)*ttt
            poisson%pots(ixt,iy,iz)=ciz
        enddo
        !-----------------------------------------
        iy=poisson%ngpy/2+1
        gsqy=fourpisqcellyinvsq*(iy-1)**2
        do ix=3,poisson%ngpx,2
            gsq=gsqy+fourpisqcellxinvsq*((ix-1)/2)**2
            g=sqrt(gsq)
            !****************************************************  
            temp_exp_l  = exp(-g*cell(3))
            temp_exp_zl1= exp(g*(zz-cell(3)))
            temp_exp_z  = exp(-g*(zz))
            temp_exp_2  = temp_exp_l**2 
            temp_exp_zl2= temp_exp_zl1*temp_exp_l     !exp(g*(zz-2.d0*cell(3)))
            temp_exp_zl = temp_exp_z  *temp_exp_l     ! exp(-g*(zz+cell(3)))
            t = temp_exp_zl1 - temp_exp_zl
            tt= temp_exp_z - temp_exp_zl2
            ttt = 1.d0/(1.d0-temp_exp_2)

            ciz= (poisson%pots(ix,iy,npu)*t+poisson%pots(ix,iy,npl)*tt)*ttt
            poisson%pots(ix,iy,iz)=ciz

            ixt=ix+1
            ciz= (poisson%pots(ixt,iy,npu)*t+poisson%pots(ixt,iy,npl)*tt)*ttt
            poisson%pots(ixt,iy,iz)=ciz
        enddo
        !-----------------------------------------
        do iy=2,poisson%ngpy/2
            iyt=poisson%ngpy-iy+2
            gsqy=fourpisqcellyinvsq*(iy-1)**2
            do ix=3,poisson%ngpx,2
                gsq=gsqy+fourpisqcellxinvsq*((ix-1)/2)**2
                g=sqrt(gsq)
                !****************************************************  
                temp_exp_l  = exp(-g*cell(3))
                temp_exp_zl1= exp(g*(zz-cell(3)))
                temp_exp_z  = exp(-g*(zz))
                temp_exp_2  = temp_exp_l**2 
                temp_exp_zl2= temp_exp_zl1*temp_exp_l     !exp(g*(zz-2.d0*cell(3)))
                temp_exp_zl = temp_exp_z  *temp_exp_l     ! exp(-g*(zz+cell(3)))
                t = temp_exp_zl1 - temp_exp_zl
                tt= temp_exp_z - temp_exp_zl2
                ttt = 1.d0/(1.d0-temp_exp_2)

                ciz= (poisson%pots(ix,iy,npu)*t+poisson%pots(ix,iy,npl)*tt)*ttt
                poisson%pots(ix,iy,iz)=ciz

                ixt=ix+1
                ciz= (poisson%pots(ixt,iy,npu)*t+poisson%pots(ixt,iy,npl)*tt)*ttt
                poisson%pots(ixt,iy,iz)=ciz

                ciz= (poisson%pots(ix,iyt,npu)*t+poisson%pots(ix,iyt,npl)*tt)*ttt
                poisson%pots(ix,iyt,iz)=ciz

                ixt=ix+1
                ciz= (poisson%pots(ixt,iyt,npu)*t+poisson%pots(ixt,iyt,npl)*tt)*ttt
                poisson%pots(ixt,iyt,iz)=ciz
            enddo
        enddo
    enddo
  
    do iz=npl,npu
        call dfftw_execute(plan_bs(iz))
    enddo

    valuengpxyinv=1.d0/real(poisson%ngpx*poisson%ngpy,8)
    do iz=npl,npu
    do iy=1,poisson%ngpy
    do ix=1,poisson%ngpx
        poisson%pots(ix,iy,iz)=poisson%pots(ix,iy,iz)*valuengpxyinv
    enddo
    enddo
    enddo
 
    do iz=npl,npu
        call dfftw_destroy_plan(plan_bs(iz))
    enddo
    call dfftw_destroy_plan(plan_fs(1))
    call dfftw_destroy_plan(plan_fs(2))
    deallocate(plan_bs)
    deallocate(plan_fs)
end subroutine sollaplaceq
!*****************************************************************************************
 subroutine calculate_force_ener_plane(atoms,poisson,epot,nbgpz)
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms, update_ratp
    use yaml_output
    implicit none
    type(typ_poisson), intent(inout):: poisson
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in):: nbgpz
    !local variables
    integer:: ix,iy,iz ,jx,jy,jz, kx,ky,kz, iat
    integer:: nlgx,nlgy,nlgz 
    integer:: ngpx,ngpy,ngpz 
    integer:: ngpx2,ngpy2,ngpz2 
    integer:: npl,npu 
    real(8), allocatable:: wx(:), wy(:), wz(:) 
    real(8), allocatable:: LGx(:), LGy(:), LGz(:) 
    real(8), allocatable:: DLGx(:), DLGy(:), DLGz(:) 
    real(8), allocatable::pots_local(:,:,:),potg(:,:,:),pot_atom(:)
    integer:: ix1,iy1,iz1,ii 
    real(8):: time1,time2
    real(8):: x,y,z ,t,tl ,epot ,t1,t2
    real(8):: tfx,tfy,tfz ,tmp
    real(8):: fatp(3,atoms%nat) 
    nlgx=9;  nlgy=9; nlgz=9
    ngpx=poisson%ngpx
    ngpy=poisson%ngpy
    ngpz=poisson%ngpz
    npl=poisson%npl
    npu=poisson%npu
    call update_ratp(atoms)
    do iat=1,atoms%nat
        if (atoms%ratp(3,iat) < (int((nlgz/2.d0))*poisson%hz)) then
             write(*,*) 'ERROR:atoms are too close to lower plane ',atoms%ratp(3,iat),npl,ceiling(nlgz/2.d0)*poisson%hz 
             stop
        else if (poisson%cell(3)-atoms%ratp(3,iat) < (int((nlgz/2.d0))*poisson%hz)) then
             write(*,*) 'ERROR:atoms are too close to upper plane' 
             stop
        endif
    enddo
    allocate(pots_local(-nlgx+1:poisson%ngpx+nlgx,-nlgy+1:poisson%ngpy+nlgy,npl:npu))
    allocate(wx(1:nlgx))
    allocate(wy(1:nlgy))
    allocate(wz(1:nlgz))
    allocate(LGx(1:nlgx))
    allocate(LGy(1:nlgy))
    allocate(LGz(1:nlgz))
    allocate(DLGx(1:nlgx))
    allocate(DLGy(1:nlgy))
    allocate(DLGz(1:nlgz))
    allocate(pot_atom(1:atoms%nat))
    do iz=npl,npu
    do ix=1-nlgx,poisson%ngpx+nlgx
    do iy=-nlgy+1,poisson%ngpy+nlgy
        pots_local(ix,iy,iz)=poisson%pots(modulo(ix-1,ngpx)+1,modulo(iy-1,ngpy)+1,iz)
    enddo
    enddo
    enddo

    call LG_weight(nlgx,nlgy,nlgz,poisson%hx,poisson%hy,poisson%hz,wx,wy,wz)
    pot_atom=0.d0
    t1=0.d0
    t2=0.d0
    fatp=0.d0
    call cpu_time (time1)
    call update_ratp(atoms)
    do iat=1,atoms%nat
        x=atoms%ratp(1,iat)
        y=atoms%ratp(2,iat)
        z=atoms%ratp(3,iat)
        !call LGW(nlgx, wx,poisson%hx, x, LGx, DLGx, ix1, 0)
        !call LGW(nlgy, wy,poisson%hy, y, LGy, DLGy, iy1, 0)
        !call LGW(nlgz, wz,poisson%hz, z, LGz, DLGz, iz1, nbgpz )
        
        call LGW4(nlgx, wx,poisson%hx, x, LGx, DLGx, ix1, 0)
        call LGW4(nlgy, wy,poisson%hy, y, LGy, DLGy, iy1, 0)
        call LGW4(nlgz, wz,poisson%hz, z, LGz, DLGz, iz1, nbgpz )
        do iz=1,nlgz 
            jz=iz1+iz-1
            do iy=1,nlgy 
                jy=iy1+iy-1
                tl=LGy(iy)*LGz(iz)
                tfy=DLGy(iy)*LGz(iz)
                tfz=LGy(iy)*DLGz(iz)
                t=0.d0 ; tfx=0.d0
                do ix=1,nlgx
                    tmp=pots_local(ix1+ix-1,jy,jz)
                    t=t+tmp*LGx(ix)
                    tfx=tfx+tmp*DLGx(ix)
                enddo 
                pot_atom(iat)=pot_atom(iat)+t*tl
                fatp(1,iat)=fatp(1,iat)+tfx*tl
                fatp(2,iat)=fatp(2,iat)+tfy*t
                fatp(3,iat)=fatp(3,iat)+tfz*t
            enddo 
        enddo 
    enddo 
    call cpu_time (time2)
!    write(*,*)"------------------------------------------"
!    write(*,*)"time for interpolation =",time2-time1
    epot=0.d0
    do iat=1,atoms%nat
        epot=epot+pot_atom(iat)*atoms%qat(iat)
        atoms%fat(1,iat)=atoms%fat(1,iat)-fatp(1,iat)*atoms%qat(iat)
        atoms%fat(2,iat)=atoms%fat(2,iat)-fatp(2,iat)*atoms%qat(iat)
        atoms%fat(3,iat)=atoms%fat(3,iat)-fatp(3,iat)*atoms%qat(iat)
    enddo
    epot=0.5*epot    

    !call yaml_map('epotplane',epot)
    !write(*,*)"epotplane=" ,epot
 !   write(*,*)"sum force",sum(fatp(1,:)),sum(fatp(2,:)),sum(fatp(3,:))
   ! write(*,*)
   ! write(*,*)"------------------------------------------"
    deallocate(pots_local)
    deallocate(wx)
    deallocate(wy)
    deallocate(wz)
    deallocate(LGx)
    deallocate(LGy)
    deallocate(LGz)
    deallocate(DLGx)
    deallocate(DLGy)
    deallocate(DLGz)
    deallocate(pot_atom)
end subroutine calculate_force_ener_plane
!*****************************************************************************************
subroutine LG_weight(nlx,nly,nlz,hx,hy,hz,wx,wy,wz)
    implicit none
    integer:: nlx ,nly, nlz !number of point for Lagrange interpolation
    integer:: i  ,maxnl
    real(8):: hx ,hy ,hz , hxp , hyp, hzp
    real(8), allocatable::factorial(:)
    real(8):: wx(nlx), wy(nly), wz(nlz) 
    maxnl=max(nlx,nly,nlz)
    allocate(factorial(0:maxnl-1))
    factorial(0)=1
    do i=1,maxnl-1
        factorial(i)=factorial(i-1)*i
    enddo
    do i=1,nlx
        wx(i)=(-1)**i*factorial(nlx-i)*factorial(i-1)
    enddo
    do i=1,nly
        wy(i)=(-1)**i*factorial(nly-i)*factorial(i-1)
    enddo
    do i=1,nlz
        wz(i)=(-1)**i*factorial(nlz-i)*factorial(i-1)
    enddo
    if (mod(nlx,2)/=0) wx=-wx
    if (mod(nly,2)/=0) wy=-wy
    if (mod(nlz,2)/=0) wz=-wz
    hxp=hx**(nlx-1)          ! how to optimize it ?
    wx(1:nlx) =1.d0/(wx(1:nlx)*hxp)
    hyp=hy**(nly-1)
    wy(1:nly) =1.d0/(wy(1:nly)*hyp)
    hzp=hz**(nlz-1)
    wz(1:nlz) =1.d0/(wz(1:nlz)*hzp)
    deallocate(factorial)
end subroutine LG_weight
!*********************************************************************************************
subroutine LGW(n, w, h, x, LGx, DLGx, ix1, nbgp)
implicit none
    integer:: ix, jx,kx, ix1,ixo ,n ,n2 ,nbgp
    real(8):: w(n), q(n),LGx(n),DLGx(n),h ,x ,diffx,x1,protot,t
    n2=(n+1)/2
    !ixo=int(real(x/h),4)+1
    ixo=nint(x/h)+1+nbgp
    ix1=ixo-(n2-1)
    x1=(ix1-1-nbgp)*h
    diffx=x-x1
    do ix=0,n-1
        q(ix+1)=diffx-(ix)*h
    enddo
    LGx=1.d0 ; DLGx=0.d0
    do ix=1,n
        do jx=1,n
            if (ix==jx) cycle
            LGx(ix)=LGx(ix)*q(jx)
            t=1.d0
            do kx=1,n
                if (kx==jx .or. kx==ix) cycle
                    t=t*q(kx)
            enddo
            DLGx(ix)=DLGx(ix)+t
        enddo
        LGx(ix)=LGx(ix)*w(ix)
        DLGx(ix)=DLGx(ix)*w(ix)
    enddo
end subroutine LGW
!*****************************************************************************
subroutine LGW4(n, w, h, x, LGx, DLGx, ix1, nbgp)
implicit none
    
    integer:: ix, jx, ix1,ixo ,n ,n2 ,nbgp
    real(8):: w(n), q(n), qinv(n), LGx(n), DLGx(n), h ,x ,diffx,x1,protot
    real(8):: t,tt 
    n2=(n+1)/2
    !ixo=int(real((x/h)),4)+1+nbgp
    ixo=nint(x/h)+1+nbgp
    ix1=ixo-(n2-1)
    x1=(ix1-1-nbgp)*h
    diffx=x-x1
    protot=1.d0
    DLGx=0.d0
    if( (ixo-1-nbgp)*h/=x)then
        do ix=0,n-1
            q(ix+1)=diffx-ix*h
            protot=protot*q(ix+1)
            qinv(ix+1)=1.d0/q(ix+1)
        enddo

        do ix=1,n
            LGx(ix)=protot*qinv(ix)
            LGx(ix)=LGx(ix)*w(ix)
            t=0.d0            
            do jx=1,n
               if (ix==jx) cycle
               t=t+qinv(jx)
            enddo
            DLGx(ix)=LGx(ix)*t
        enddo
    else
        LGx=0.d0
        LGx(n2)=1.d0
        do ix=0,n-1
            if (ix==n2-1) cycle
            q(ix+1)=diffx-ix*h
            protot=protot*q(ix+1)
        enddo
        do ix=1,n
            if (ix==n2) cycle
            DLGx(ix)=protot/q(ix)
        enddo
        tt=sum(DLGx)
        DLGx(n2)=tt
        DLGx(1:n)=DLGx(1:n)*w(1:n)
    endif
end subroutine LGW4
!*******************************************************************************************
subroutine surface_charge(parini,poisson,pot_short,vl,vu)
    use mod_electrostatics, only: typ_poisson
    use mod_parini, only: typ_parini
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_poisson), intent(inout):: poisson
    integer::ix, iy, iz, npl,npu
    real(8):: t, tt ,density(poisson%ngpx,poisson%ngpy,2),vl,vu
    real(8)::hgzinv,pi,pot_layerl,pot_layeru,pot_short(poisson%ngpx,poisson%ngpy,2,5)
    real(8):: E,d
    real(8), parameter, dimension(5) :: cf5 = [-25.d0/12.d0,4.d0,-3.d0,+4.d0/3.d0,-0.25d0]
    pi=4*atan(1.d0)
    npl=poisson%npl
    npu=poisson%npu
    hgzinv=1.d0/(poisson%hz*4.d0*pi)
    t=0.d0
    tt=0.d0
    d = poisson%cell(3)
    E =- (vu-vl)/d

    do iy=1,poisson%ngpy
    do ix=1,poisson%ngpx
        density(ix,iy,1) =-dot_product(poisson%pots(ix,iy,npl:npl+4)+poisson%pot(ix,iy,npl:npl+4)+pot_short(ix,iy,1,1:5),cf5) * hgzinv
        t=t+ density(ix,iy,1)
        density(ix,iy,2) =-dot_product(poisson%pots(ix,iy,npu:npu-4:-1)+poisson%pot(ix,iy,npu:npu-4:-1)+pot_short(ix,iy,2,1:5),cf5) * hgzinv
        tt=tt+ density(ix,iy,2)
    enddo
    enddo

    t=t*poisson%hx*poisson%hy
    tt=tt*poisson%hx*poisson%hy
    !if(trim(parini%bias_type)=='fixed_efield' .or. trim(parini%bias_type)=='fixed_potdiff') then
    !    t =t -E/(4*pi)*poisson%cell(1)*poisson%cell(2)
    !    tt=tt+E/(4*pi)*poisson%cell(1)*poisson%cell(2)
    !endif
    call yaml_mapping_open('charge on planes',flow=.true.)
    call yaml_map('lower',t,fmt='(es25.13)')
    call yaml_map('upper',tt,fmt='(es25.13)')
    call yaml_mapping_close()
    !write(*,'(a,es25.13)')'charge on lower plane' ,t
    !write(*,'(a,es25.13)')'charge on upper plane',tt
    vl=poisson%vl
    vu=poisson%vu
end subroutine surface_charge
!*****************************************************************************************
!This subroutine determines the limits of grids in a sphere.
subroutine determine_limitsphere(poisson,mboundg,mboundgy,nbgpx,nbgpy,nbgpz)
    use mod_electrostatics, only: typ_poisson
    implicit none
    type(typ_poisson), intent(inout):: poisson
    !local variables
    integer:: ix, iy, iz, mboundg(2,-nbgpy:nbgpy,-nbgpz:nbgpz), mboundgy(2,-nbgpz:nbgpz)
    integer:: nbgpx, nbgpy, nbgpz

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
    do iz=-nbgpz,nbgpz
        mboundgy(1,iz)=1 !-mlimnb
        mboundgy(2,iz)=0 !mlimnb
    enddo
    do iz=-nbgpz,nbgpz
    do iy=1,nbgpy
    if(mboundg(1,iy,iz)<=mboundg(2,iy,iz)) then
        mboundgy(1,iz)=-iy
        mboundgy(2,iz)=iy
    endif
    enddo
    enddo


end subroutine determine_limitsphere
!***********************************************************************************************************************
subroutine bias_field_potener_forces(parini,poisson,atoms,epotplane)
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_parini, only: typ_parini
    use dynamic_memory
    implicit none
    type(typ_poisson), intent(inout):: poisson
    type(typ_atoms), intent(inout):: atoms
    type(typ_parini), intent(in):: parini
    !local variables
    real(8):: epotlong, epotplane !, epotshort
    real(8):: time(10)
    !The following dummy varables definition to be deleted later.
    !real(8):: totrho
    real(8):: beta, pi, charge,c,charge0,E,tmp
    integer:: igpx, igpy, igpz, iat
    integer:: ix, iy, iz, jx, jy, jz, kx, ky, kz
    integer:: npl, npu, nlayer, ngpx, ngpy, ngpz
    real(8),allocatable :: pots_layer(:,:,:,:)
    real(8):: vl, vu, A, d, rl, ru, dipole_correction
    real(8):: pot_correction ,dipole,dv,beta2 

    beta = poisson%beta*(-poisson%ngpx*poisson%ngpy)
    pi=4.d0*atan(1.d0)
    ngpz = poisson%ngpz
    ngpx = poisson%ngpx
    ngpy = poisson%ngpy
    epotplane=0.d0
    d = poisson%cell(3)
    if(poisson%point_particle .and. trim(parini%bias_type)=='fixed_efield') then
        !E = parini%efield
    elseif(poisson%point_particle .and. trim(parini%bias_type)=='fixed_potdiff') then
        dipole=0.d0
        call update_ratp(atoms)
        do iat=1,atoms%nat
            dipole= dipole + atoms%ratp(3,iat)*atoms%qat(iat)
        enddo
        beta2 = - 2.d0*pi/(poisson%cell(2)*poisson%cell(1))*dipole
        dv = parini%vu_ewald-parini%vl_ewald 
        write(*,*)'real pot = vu , vl ',parini%vu_ewald+ beta2,parini%vl_ewald- beta2 
        E =- (dv+2.d0*beta2)/d
    endif
       A= poisson%ngpx*poisson%ngpy*poisson%hx*poisson%hy
       c= A/(4.d0*pi*d)
       charge0= -dipole/(d)
       charge = -dipole/(d)+c*(dv)
        epotplane = 0.d0
        epotplane = epotplane !+0.5*c*dv**2
        epotplane = epotplane + (dv+beta2)/d*dipole
        do iat=1,atoms%nat
            atoms%fat(3,iat)=atoms%fat(3,iat)-((dv+2.d0*beta2)/poisson%cell(3))*atoms%qat(iat)
        enddo

!!*****************************************************************************
!        poisson%npu=poisson%ngpz-poisson%nbgpz
!        poisson%npl=1+poisson%nbgpz  
!        npl=poisson%npl
!        npu=poisson%npu
!
!        allocate(poisson%pots(1:poisson%ngpx+2,1:poisson%ngpy,npl:npu))
!!        write(*,*)"npu,npl",poisson%npu,poisson%npl
!        nlayer=1
!        if (parini%cal_charge) then 
!            nlayer=5
!            allocate(pots_layer(1:poisson%ngpx,1:poisson%ngpy,1:2,1:nlayer))
!        endif
!        poisson%pots=0.d0
!        call erfc_surface_zero(parini,atoms,poisson,nlayer)
!        if (parini%cal_charge) then
!            pots_layer(1:ngpx,1:ngpy,1,:)=poisson%pots(1:ngpx,1:ngpy,npl:npl+nlayer-1)
!            pots_layer(1:ngpx,1:ngpy,2,:)=poisson%pots(1:ngpx,1:ngpy,npu:npu-(nlayer-1):-1)
!        endif
!        poisson%pots=0.d0
!
!        if (parini%cal_charge) then 
!            call surface_charge(parini,poisson,pots_layer,vl,vu)
!            deallocate(pots_layer)
!        endif
!        deallocate(poisson%pots)
!!*****************************************************************************
end subroutine bias_field_potener_forces
