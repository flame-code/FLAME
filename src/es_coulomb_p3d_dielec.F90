subroutine dielec_potener_forces(parini,poisson,atoms,epot_dielec)
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms
    use mod_potential, only: potential 
    use mod_parini, only: typ_parini
    use dynamic_memory
    implicit none
    type(typ_poisson), intent(inout):: poisson
    type(typ_atoms), intent(inout):: atoms
    type(typ_parini), intent(in):: parini
    real(8), intent(out):: epot_dielec
    !local variables
    !Do be added by Farhad.
    REAL(8), PARAMETER :: PI = 3.14159265358979312
    real(8):: time(10)
    !The following dummy varables definition to be deleted later.
    real(8):: beta, charge,c,charge0,E,tmp
    integer:: igpx, igpy, igpz, iat
    integer:: ix, iy, iz, jx, jy, jz, kx, ky, kz
    integer:: npl, npu, nlayer, ngpx, ngpy, ngpz, nbgpz, nlayer2
    real(8),allocatable :: pots_layer(:,:,:,:)
    real(8):: vl, vu, A, d, rl, ru, dipole_correction, dipole

    epot_dielec=0.d0
    beta = poisson%beta*(-poisson%ngpx*poisson%ngpy)  
    ngpx= poisson%ngpx
    ngpy= poisson%ngpy
    nbgpz=int(poisson%rgcut/poisson%hz)+2
    poisson%npu=poisson%ngpz-nbgpz 
    poisson%npl=1+nbgpz  
    npl=poisson%npl
    npu=poisson%npu
    nlayer=7 
    nlayer2=(nlayer-1)/2

    allocate(poisson%pots(1:poisson%ngpx+2,1:poisson%ngpy,npl-nlayer2:npu+nlayer2))
    allocate(poisson%dpot(1:poisson%ngpx+2,1:poisson%ngpy,2))
    poisson%pots=0.d0
    poisson%dpot=0.d0
    
    !**********************************************************************************
    poisson%npl=poisson%npl-nlayer2
    poisson%npu=poisson%npu+nlayer2
    allocate(pots_layer(1:poisson%ngpx,1:poisson%ngpy,1:2,1:nlayer)) 
    pots_layer = 0.d0
    if(.not. (trim(potential)=='ann')) then
        call erfc_surface_zero(parini,atoms,poisson,nlayer) 
        pots_layer(1:ngpx,1:ngpy,1,:)=poisson%pots(1:ngpx,1:ngpy,npl-nlayer2:npl+nlayer2) 
        pots_layer(1:ngpx,1:ngpy,2,:)=poisson%pots(1:ngpx,1:ngpy,npu+nlayer2:npu-nlayer2:-1)
        poisson%pots  =0.d0
    endif
    poisson%npl=poisson%npl+nlayer2
    poisson%npu=poisson%npu-nlayer2
    !**********************************************************************************

    call diff_pot_pp(parini,poisson,pots_layer,vl,vu,nlayer)
    call sollaplaceq_dielctric(parini,poisson,poisson%hz,poisson%cell)
    call calculate_force_ener_plane(atoms,poisson,epot_dielec,nbgpz)

    deallocate(pots_layer)
    deallocate(poisson%pots)
    deallocate(poisson%dpot)
end subroutine dielec_potener_forces
!*******************************************************************************************
subroutine sollaplaceq_dielctric(parini,poisson,hz,cell)
    use mod_electrostatics, only: typ_poisson
    use mod_parini, only: typ_parini
    implicit none
    include 'fftw3.f'
    type(typ_parini), intent(in):: parini

    type(typ_poisson), intent(inout):: poisson
    !local variables
    real(8):: vl, vu , zlmzu , sinhzlmzu, zlmzuinv
    real(8):: cell(3)
    real(8):: hz , vlmvu, vlzumvuzl 
    real(8):: t,tt,ttt, ghz,fkl,potl,potu
    real(8):: pi, fourpi, fourpisq, gsq, gsqx, gsqy, g, hzsq
    real(8):: fourpisqcellxinvsq, fourpisqcellyinvsq,valuengpxyinv
    real(8):: hzlu,hzliz,hzuiz 
    real(8):: tel,teu,telu, beta, beta_1, beta_2, r_1, r_2 
    integer::ix,iy,iz,ixt,iyt,npu,npl, izz
    integer(8), allocatable:: plan_bs(:),plan_fs(:)
    real(8) :: ep_1_p, ep_1_m, ep_2_p, ep_2_m
    real(8) :: eps0, eps1, eps2, temp_sinh, temp_cosh
    real(8) :: temp_exp_2  , temp_exp_zl1, temp_exp_zl2, temp_exp_z  , temp_exp_zl ,temp_exp 
    real(8) :: tt1, tt2, tt3, tt4, tt5 
    pi=4.d0*atan(1.d0)
    fourpi=4.d0*pi 
    npl=poisson%npl
    npu=poisson%npu
    hzsq=hz**2
    fourpisq=fourpi*pi 
    fourpisqcellxinvsq=fourpisq/cell(1)**2 
    fourpisqcellyinvsq=fourpisq/cell(2)**2 
    
    eps0 = parini%dielec_const/fourpi
    eps1 = parini%dielec_const1/fourpi
    eps2 = parini%dielec_const2/fourpi

    ep_1_p = eps0 + eps1
    ep_1_m = eps0 - eps1
    ep_2_p = eps0 + eps2
    ep_2_m = eps0 - eps2
    r_1 = ep_1_m / ep_1_p
    r_2 = ep_2_m / ep_2_p
    allocate(plan_bs(npl:npu))
    allocate(plan_fs(1:2))
    do iz=npl,npu
        call dfftw_plan_dft_c2r_2d(plan_bs(iz),poisson%ngpx, &
            poisson%ngpy,poisson%pots(1,1,iz),poisson%pots(1,1,iz),fftw_estimate)
    enddo
    call dfftw_plan_dft_r2c_2d(plan_fs(1),poisson%ngpx, &
            poisson%ngpy,poisson%dpot(1,1,1), &
            poisson%dpot(1,1,1),fftw_estimate)
    call dfftw_plan_dft_r2c_2d(plan_fs(2),poisson%ngpx, &
            poisson%ngpy,poisson%dpot(1,1,2), &
            poisson%dpot(1,1,2),fftw_estimate)

      do iy=1,poisson%ngpy
      do ix=1,poisson%ngpx
          poisson%dpot(ix,iy,2)=poisson%dpot(ix,iy,2)*(eps0 - eps2)
          poisson%dpot(ix,iy,1)=poisson%dpot(ix,iy,1)*(eps0 - eps1)
      enddo 
      enddo
      call dfftw_execute(plan_fs(1))
      call dfftw_execute(plan_fs(2))   

   ! k=0 , l=0
    beta = poisson%beta
    beta_1 = - beta * 2.0d0 / ep_1_p * (1.0d0 - r_2) / (1.0d0 - r_1*r_2)
    beta_2 = - beta * 2.0d0 / ep_2_p * (r_1 - 1.0d0) / (1.0d0 - r_1*r_2)

    !-----------------------------------------
    ix=1;iy=1
         poisson%pots(ix,iy,npu)= beta_1 + beta / eps0
         poisson%pots(ix,iy,npl)= beta_2 - beta / eps0
     do iz=npl+1,npu-1
         izz = iz - npl
         fkl= (beta_1 - beta_2 + 2.0d0 * beta / eps0) *(izz * hz) / cell(3) + (beta_2 - beta / eps0) 
         poisson%pots(ix,iy,iz)=fkl
     enddo 
    !-----------------------------------------
     hzlu=-hz*(npl-npu)
     do iz=npl,npu
         izz = iz - npl
         hzliz=-hz*(npl-iz)
         hzuiz=-hz*(iz-npu)
         !-----------------------------------------
         ix=poisson%ngpx+1;iy=1
         gsq=fourpisqcellxinvsq*((ix-1)/2)**2 + fourpisqcellyinvsq*(iy-1)**2
         g=sqrt(gsq)
         !************************************************
         temp_exp    = exp(-g*cell(3))
         temp_exp_2  = temp_exp**2
         temp_exp_zl1= exp(g*(izz*hz-cell(3)))
         temp_exp_zl2= temp_exp_zl1*temp_exp 
         temp_exp_z  = exp(-g*(izz*hz))
         temp_exp_zl = temp_exp_z*temp_exp 

         tt1= g *(ep_1_p*ep_2_p-ep_1_m*ep_2_m*temp_exp_2)
         tt2= - ep_2_p * temp_exp_zl1
         tt3= + ep_1_m * temp_exp_zl2 
         tt4= - ep_2_m * temp_exp_zl
         tt5= + ep_1_p * temp_exp_z

         t  = (tt2* poisson%dpot(ix,iy,1) + tt3 *  poisson%dpot(ix,iy,2))/tt1
         tt = (tt4* poisson%dpot(ix,iy,1) + tt5 *  poisson%dpot(ix,iy,2))/tt1
         fkl = t+ tt 
         !************************************************
         poisson%pots(ix,iy,iz)=fkl
         !-----------------------------------------
         ix=1;iy=poisson%ngpy/2+1
         gsq=fourpisqcellyinvsq*(iy-1)**2
         g=sqrt(gsq)
         !************************************************
         temp_exp    = exp(-g*cell(3))
         temp_exp_2  = temp_exp**2
         temp_exp_zl1= exp(g*(izz*hz-cell(3)))
         temp_exp_zl2= temp_exp_zl1*temp_exp 
         temp_exp_z  = exp(-g*(izz*hz))
         temp_exp_zl = temp_exp_z*temp_exp 

         tt1= g *(ep_1_p*ep_2_p-ep_1_m*ep_2_m*temp_exp_2)
         tt2= - ep_2_p * temp_exp_zl1
         tt3= + ep_1_m * temp_exp_zl2 
         tt4= - ep_2_m * temp_exp_zl
         tt5= + ep_1_p * temp_exp_z

         t  = (tt2* poisson%dpot(ix,iy,1) + tt3 *  poisson%dpot(ix,iy,2))/tt1
         tt = (tt4* poisson%dpot(ix,iy,1) + tt5 *  poisson%dpot(ix,iy,2))/tt1
         fkl = t+ tt 
         !************************************************
         poisson%pots(ix,iy,iz)=fkl
         !-----------------------------------------
         ix=poisson%ngpx+1;iy=poisson%ngpy/2+1
         gsq=fourpisqcellxinvsq*((ix-1)/2)**2 + fourpisqcellyinvsq*(iy-1)**2
         g=sqrt(gsq)
         !************************************************
         temp_exp    = exp(-g*cell(3))
         temp_exp_2  = temp_exp**2
         temp_exp_zl1= exp(g*(izz*hz-cell(3)))
         temp_exp_zl2= temp_exp_zl1*temp_exp 
         temp_exp_z  = exp(-g*(izz*hz))
         temp_exp_zl = temp_exp_z*temp_exp 

         tt1= g *(ep_1_p*ep_2_p-ep_1_m*ep_2_m*temp_exp_2)
         tt2= - ep_2_p * temp_exp_zl1
         tt3= + ep_1_m * temp_exp_zl2 
         tt4= - ep_2_m * temp_exp_zl
         tt5= + ep_1_p * temp_exp_z

         t  = (tt2* poisson%dpot(ix,iy,1) + tt3 *  poisson%dpot(ix,iy,2))/tt1
         tt = (tt4* poisson%dpot(ix,iy,1) + tt5 *  poisson%dpot(ix,iy,2))/tt1
         fkl = t+ tt 
         !************************************************
         poisson%pots(ix,iy,iz)=fkl
         !-----------------------------------------
         ix=1
         do iy=2,poisson%ngpy/2
             gsq=fourpisqcellyinvsq*(iy-1)**2
             g=sqrt(gsq)
             !************************************************
             temp_exp    = exp(-g*cell(3))
             temp_exp_2  = temp_exp**2
             temp_exp_zl1= exp(g*(izz*hz-cell(3)))
             temp_exp_zl2= temp_exp_zl1*temp_exp 
             temp_exp_z  = exp(-g*(izz*hz))
             temp_exp_zl = temp_exp_z*temp_exp 

             tt1= g *(ep_1_p*ep_2_p-ep_1_m*ep_2_m*temp_exp_2)
             tt2= - ep_2_p * temp_exp_zl1
             tt3= + ep_1_m * temp_exp_zl2 
             tt4= - ep_2_m * temp_exp_zl
             tt5= + ep_1_p * temp_exp_z

             t  = (tt2* poisson%dpot(ix,iy,1) + tt3 *  poisson%dpot(ix,iy,2))/tt1
             tt = (tt4* poisson%dpot(ix,iy,1) + tt5 *  poisson%dpot(ix,iy,2))/tt1
             fkl = t+ tt 
             !************************************************
             poisson%pots(ix,iy,iz)=fkl
             poisson%pots(ix,poisson%ngpy-iy+2,iz)=fkl

             ixt=ix+1
             !************************************************
             t  = (tt2* poisson%dpot(ixt,iy,1) + tt3 *  poisson%dpot(ixt,iy,2))/tt1
             tt = (tt4* poisson%dpot(ixt,iy,1) + tt5 *  poisson%dpot(ixt,iy,2))/tt1
             fkl = t+ tt 
             !************************************************
             poisson%pots(ixt,iy,iz)=fkl
             poisson%pots(ixt,poisson%ngpy-iy+2,iz)=-fkl
         enddo
         !-----------------------------------------
         ix=poisson%ngpx+1
         gsqx=fourpisqcellxinvsq*(poisson%ngpx/2)**2
         do iy=2,poisson%ngpy/2
             gsq=gsqx+fourpisqcellyinvsq*(iy-1)**2
             g=sqrt(gsq)
             !************************************************
             temp_exp    = exp(-g*cell(3))
             temp_exp_2  = temp_exp**2
             temp_exp_zl1= exp(g*(izz*hz-cell(3)))
             temp_exp_zl2= temp_exp_zl1*temp_exp 
             temp_exp_z  = exp(-g*(izz*hz))
             temp_exp_zl = temp_exp_z*temp_exp 

             tt1= g *(ep_1_p*ep_2_p-ep_1_m*ep_2_m*temp_exp_2)
             tt2= - ep_2_p * temp_exp_zl1
             tt3= + ep_1_m * temp_exp_zl2 
             tt4= - ep_2_m * temp_exp_zl
             tt5= + ep_1_p * temp_exp_z

             t  = (tt2* poisson%dpot(ix,iy,1) + tt3 *  poisson%dpot(ix,iy,2))/tt1
             tt = (tt4* poisson%dpot(ix,iy,1) + tt5 *  poisson%dpot(ix,iy,2))/tt1
             fkl = t+ tt 
             !************************************************
             poisson%pots(ix,iy,iz)=fkl
             poisson%pots(ix,poisson%ngpy-iy+2,iz)=fkl

             ixt=ix+1
             !************************************************
             t  = (tt2* poisson%dpot(ixt,iy,1) + tt3 *  poisson%dpot(ixt,iy,2))/tt1
             tt = (tt4* poisson%dpot(ixt,iy,1) + tt5 *  poisson%dpot(ixt,iy,2))/tt1
             fkl = t+ tt 
             !************************************************
             poisson%pots(ixt,iy,iz)=fkl
             poisson%pots(ixt,poisson%ngpy-iy+2,iz)=-fkl
         enddo
         !-----------------------------------------
         iy=1
         gsqy=fourpisqcellyinvsq*(iy-1)**2
         do ix=3,poisson%ngpx,2
             gsq=gsqy+fourpisqcellxinvsq*((ix-1)/2)**2
             g=sqrt(gsq)
             !************************************************
             temp_exp    = exp(-g*cell(3))
             temp_exp_2  = temp_exp**2
             temp_exp_zl1= exp(g*(izz*hz-cell(3)))
             temp_exp_zl2= temp_exp_zl1*temp_exp 
             temp_exp_z  = exp(-g*(izz*hz))
             temp_exp_zl = temp_exp_z*temp_exp 

             tt1= g *(ep_1_p*ep_2_p-ep_1_m*ep_2_m*temp_exp_2)
             tt2= - ep_2_p * temp_exp_zl1
             tt3= + ep_1_m * temp_exp_zl2 
             tt4= - ep_2_m * temp_exp_zl
             tt5= + ep_1_p * temp_exp_z

             t  = (tt2* poisson%dpot(ix,iy,1) + tt3 *  poisson%dpot(ix,iy,2))/tt1
             tt = (tt4* poisson%dpot(ix,iy,1) + tt5 *  poisson%dpot(ix,iy,2))/tt1
             fkl = t+ tt 
             !************************************************
             poisson%pots(ix,iy,iz)=fkl

             ixt=ix+1
             !************************************************
             t  = (tt2* poisson%dpot(ixt,iy,1) + tt3 *  poisson%dpot(ixt,iy,2))/tt1
             tt = (tt4* poisson%dpot(ixt,iy,1) + tt5 *  poisson%dpot(ixt,iy,2))/tt1
             fkl = t+ tt 
             !************************************************
             poisson%pots(ixt,iy,iz)=fkl
         enddo
         !-----------------------------------------
         iy=poisson%ngpy/2+1
         gsqy=fourpisqcellyinvsq*(iy-1)**2
         do ix=3,poisson%ngpx,2
             gsq=gsqy+fourpisqcellxinvsq*((ix-1)/2)**2
             g=sqrt(gsq)
             !************************************************
             temp_exp    = exp(-g*cell(3))
             temp_exp_2  = temp_exp**2
             temp_exp_zl1= exp(g*(izz*hz-cell(3)))
             temp_exp_zl2= temp_exp_zl1*temp_exp 
             temp_exp_z  = exp(-g*(izz*hz))
             temp_exp_zl = temp_exp_z*temp_exp 

             tt1= g *(ep_1_p*ep_2_p-ep_1_m*ep_2_m*temp_exp_2)
             tt2= - ep_2_p * temp_exp_zl1
             tt3= + ep_1_m * temp_exp_zl2 
             tt4= - ep_2_m * temp_exp_zl
             tt5= + ep_1_p * temp_exp_z

             t  = (tt2* poisson%dpot(ix,iy,1) + tt3 *  poisson%dpot(ix,iy,2))/tt1
             tt = (tt4* poisson%dpot(ix,iy,1) + tt5 *  poisson%dpot(ix,iy,2))/tt1
             fkl = t+ tt 
             !************************************************
             poisson%pots(ix,iy,iz)=fkl

             ixt=ix+1
             !************************************************
             t  = (tt2* poisson%dpot(ixt,iy,1) + tt3 *  poisson%dpot(ixt,iy,2))/tt1
             tt = (tt4* poisson%dpot(ixt,iy,1) + tt5 *  poisson%dpot(ixt,iy,2))/tt1
             fkl = t+ tt 
             !************************************************
             poisson%pots(ixt,iy,iz)=fkl
         enddo
         !-----------------------------------------
         do iy=2,poisson%ngpy/2
             iyt=poisson%ngpy-iy+2
             gsqy=fourpisqcellyinvsq*(iy-1)**2
             do ix=3,poisson%ngpx,2
                 gsq=gsqy+fourpisqcellxinvsq*((ix-1)/2)**2
                 g=sqrt(gsq)
                 !************************************************
                 temp_exp    = exp(-g*cell(3))
                 temp_exp_2  = temp_exp**2
                 temp_exp_zl1= exp(g*(izz*hz-cell(3)))
                 temp_exp_zl2= temp_exp_zl1*temp_exp 
                 temp_exp_z  = exp(-g*(izz*hz))
                 temp_exp_zl = temp_exp_z*temp_exp 

                 tt1= g *(ep_1_p*ep_2_p-ep_1_m*ep_2_m*temp_exp_2)
                 tt2= - ep_2_p * temp_exp_zl1
                 tt3= + ep_1_m * temp_exp_zl2 
                 tt4= - ep_2_m * temp_exp_zl
                 tt5= + ep_1_p * temp_exp_z

                 t  = (tt2* poisson%dpot(ix,iy,1) + tt3 *  poisson%dpot(ix,iy,2))/tt1
                 tt = (tt4* poisson%dpot(ix,iy,1) + tt5 *  poisson%dpot(ix,iy,2))/tt1
                 fkl = t+ tt 
                 !************************************************
                 poisson%pots(ix,iy,iz)=fkl

                 ixt=ix+1
                 !************************************************
                 t  = (tt2* poisson%dpot(ixt,iy,1) + tt3 *  poisson%dpot(ixt,iy,2))/tt1
                 tt = (tt4* poisson%dpot(ixt,iy,1) + tt5 *  poisson%dpot(ixt,iy,2))/tt1
                 fkl = t+ tt 
                 !************************************************
                 poisson%pots(ixt,iy,iz)=fkl

                 !************************************************
                 t  = (tt2* poisson%dpot(ix,iyt,1) + tt3 *  poisson%dpot(ix,iyt,2))/tt1
                 tt = (tt4* poisson%dpot(ix,iyt,1) + tt5 *  poisson%dpot(ix,iyt,2))/tt1
                 fkl = t+ tt 
                 !************************************************
                 poisson%pots(ix,iyt,iz)=fkl

                 ixt=ix+1
                 !************************************************
                 t  = (tt2* poisson%dpot(ixt,iyt,1) + tt3 *  poisson%dpot(ixt,iyt,2))/tt1
                 tt = (tt4* poisson%dpot(ixt,iyt,1) + tt5 *  poisson%dpot(ixt,iyt,2))/tt1
                 fkl = t+ tt 
                 !************************************************
                 poisson%pots(ixt,iyt,iz)=fkl
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
    call dfftw_destroy_plan(plan_fs)
    deallocate(plan_bs)
    deallocate(plan_fs)
end subroutine sollaplaceq_dielctric

!*****************************************************************************************
!subroutine diff_pot_pp2(parini,poisson,pot_short,vl,vu,nlayer)
!    use mod_electrostatics, only: typ_poisson
!    use mod_parini, only: typ_parini
!    implicit none
!
!!FARHAD:: define a real kind with at least 15-digit precision and range 300 ::
!    integer, parameter :: DP = selected_real_kind(P=15,R=300)
!    
!    type(typ_parini), intent(in):: parini
!    type(typ_poisson), intent(inout):: poisson
!    integer::ix, iy, iz, npl,npu,nlayer
!    real(8):: vl,vu
!    real(8)::hgzinv,pi
!    real(8):: E
!    real(8) :: pot_short(poisson%ngpx,poisson%ngpy,2,nlayer)
!!FARHAD:: cf8 :: coefficients for Forward and Backward Num-Diff, order :: 8
!    real(kind=dp), parameter, dimension(9) :: cf8 = [-2.717857142857143_dp,8.0_dp,-14.0_dp,&
!    &18.666666666666668_dp,-17.5_dp,11.2_dp,-4.666666666666667_dp,1.142857142857143_dp,-0.125_dp]
!!FARHAD:: unused vars ::
!!    real(8):: d
!!    real(8):: t, tt ,density(poisson%ngpx,poisson%ngpy,2)
!!    real(kind=dp) :: pot_layerl , pot_layeru
!!    real(kind=dp) :: pot_layerl2, pot_layeru2
!!    real(kind=dp) :: pot_layerl3, pot_layeru3
!!    real(kind=dp) :: pot_layerl4, pot_layeru4
!!    real(kind=dp) :: pot_layerl5, pot_layeru5
!!    real(kind=dp) :: pot_layerl6, pot_layeru6
!!    real(kind=dp) :: pot_layerl7, pot_layeru7
!!    real(kind=dp) :: pot_layerl8, pot_layeru8
!!    real(8), parameter :: pi = 3.14159265358979312
!
!    npl=poisson%npl
!    npu=poisson%npu
!    hgzinv=1.d0/(poisson%hz)
!!    t=0.d0
!!    tt=0.d0
!!    d = poisson%cell(3)
!
!    do iy=1,poisson%ngpy
!        do ix=1,poisson%ngpx
!!FARHAD:: maybe "nlayer" can be set as global and the order of index of %pot could be inverse for better Perf.
!          poisson%dpot(ix,iy,2) = dot_product( poisson%pot(ix,iy,npl:npl+8) + pot_short(ix,iy,1,1:9),cf8) * hgzinv
!          poisson%dpot(ix,iy,1) = - dot_product( poisson%pot(ix,iy,npu:npu-8:-1)+pot_short(ix,iy,2,1:9) , cf8) * hgzinv
!!                pot_layerl4 = poisson%pot(ix,iy,npl+4)+pot_short(ix,iy,1,5)
!!                pot_layerl3 = poisson%pot(ix,iy,npl+3)+pot_short(ix,iy,1,4)
!!                pot_layerl2 = poisson%pot(ix,iy,npl+2)+pot_short(ix,iy,1,3)
!!                pot_layerl  = poisson%pot(ix,iy,npl+1)+pot_short(ix,iy,1,2)
!!                vl          = poisson%pot(ix,iy,npl  )+pot_short(ix,iy,1,1)
!!                !density(ix,iy,1)=-0.5d0*(-3.d0*vl+4*pot_layerl-pot_layerl2)* hgzinv
!!                poisson%dpot(ix,iy,2)=(-25.d0/12.d0*vl+4.d0*pot_layerl-3.d0*pot_layerl2+4.d0/3.d0*pot_layerl3-0.25d0*pot_layerl4)* hgzinv
!!                
!!                pot_layeru4 = poisson%pot(ix,iy,npu-4)+pot_short(ix,iy,2,5)
!!                pot_layeru3 = poisson%pot(ix,iy,npu-3)+pot_short(ix,iy,2,4)
!!                pot_layeru2 = poisson%pot(ix,iy,npu-2)+pot_short(ix,iy,2,3)
!!                pot_layeru  = poisson%pot(ix,iy,npu-1)+pot_short(ix,iy,2,2)
!!                vu          = poisson%pot(ix,iy,npu  )+pot_short(ix,iy,2,1)
!!                !density(ix,iy,2)=0.5d0*(3.d0*vu-4.d0*pot_layeru+pot_layeru2)* hgzinv
!!                poisson%dpot(ix,iy,1)=-(-25.d0/12.d0*vu+4.d0*pot_layeru-3.d0*pot_layeru2+4.d0/3.d0*pot_layeru3-0.25d0*pot_layeru4)* hgzinv
!        enddo
!    enddo
!end subroutine diff_pot_pp2
!*****************************************************************************************
subroutine diff_pot_pp(parini,poisson,pot_short,vl,vu,nlayer)
    use mod_electrostatics, only: typ_poisson
    use mod_parini, only: typ_parini
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_poisson), intent(inout):: poisson
    integer::ix, iy, iz, npl,npu,nlayer
    real(8):: vl,vu
    real(8)::hgzinv,pi
    real(8):: E
    real(8) :: pot_short(poisson%ngpx,poisson%ngpy,2,nlayer)
    real(8), parameter, dimension(7) :: cf7 = [-1.d0/60.d0,3.d0/20.d0,-3.d0/4.d0,0.d0,3.d0/4.d0,-3.d0/20.d0,1.d0/60.d0]

    npl=poisson%npl
    npu=poisson%npu
    hgzinv=1.d0/(poisson%hz)

    do iy=1,poisson%ngpy
        do ix=1,poisson%ngpx
          poisson%dpot(ix,iy,2) =  dot_product( poisson%pot(ix,iy,npl-3:npl+3   )+pot_short(ix,iy,1,1:7) , cf7) * hgzinv
          poisson%dpot(ix,iy,1) = -dot_product( poisson%pot(ix,iy,npu+3:npu-3:-1)+pot_short(ix,iy,2,1:7) , cf7) * hgzinv
        enddo
    enddo
end subroutine diff_pot_pp
!*****************************************************************************************
