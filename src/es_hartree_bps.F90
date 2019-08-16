!*****************************************************************************************
subroutine get_psolver_bps(poisson,atoms,ehartree)
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
#if defined(HAVE_BPS)
    use Poisson_Solver, only: H_POTENTIAL
#endif
    use dynamic_memory
    implicit none
    type(typ_poisson),intent(inout):: poisson
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(out):: ehartree
    !local variables
    integer:: igpx, igpy, igpz
#if defined(HAVE_BPS)
    real(8), allocatable:: pot_ion(:)
    real(8):: stress_t(6)
    !type(coulomb_operator):: pkernel
    if(trim(atoms%boundcond)/='bulk' .and. trim(atoms%boundcond)/='free') then
        write(*,*) 'ERROR: BPS for BCs other than bulk or free not implemented yet in FLAME.'
        stop
    endif
    pot_ion=f_malloc([1.to.1],id='pot_ion') !an array which must be provided but will not be used.
#endif

#if defined(HAVE_BPS)
    !The reason I call H_POTENTIAL with rho and after that pot=rho,
    !is because I was worried wether pot is allocated with (ngpx+2,ngpy,ngpz),
    !however, I can check or make it to be allocated as (ngpx,ngpy,ngpz), if
    !BPS is used.
    !stress_t(1:6)=0.d0
    call H_POTENTIAL('G',poisson%pkernel,poisson%rho, &
        pot_ion,ehartree,0.d0,.false.,stress_tensor=stress_t) !,quiet='yes')
    !write(*,'(a,6es14.5)') 'STRESS ',stress_t(1:6)
    !ordering from BigDFT ---> (11,22,33,23,13,12)
    atoms%stress(1,1)=stress_t(1)
    atoms%stress(2,2)=stress_t(2)
    atoms%stress(3,3)=stress_t(3)
    atoms%stress(2,3)=stress_t(4)
    atoms%stress(1,3)=stress_t(5)
    atoms%stress(1,2)=stress_t(6)
    atoms%stress(3,2)=stress_t(4)
    atoms%stress(3,1)=stress_t(5)
    atoms%stress(2,1)=stress_t(6)
    do igpz=1,poisson%ngpz
    do igpy=1,poisson%ngpy
    do igpx=1,poisson%ngpx
        poisson%pot(igpx,igpy,igpz)=poisson%rho(igpx,igpy,igpz)
    enddo
    enddo
    enddo
    call f_free(pot_ion)
#else
    stop 'ERROR: Alborz is not linked with Poisson solvers in BigDFT.'
    ehartree=0.d0 !this is just to be able to compile.
#endif
end subroutine get_psolver_bps
!*****************************************************************************************
subroutine init_psolver_bps(parini,atoms,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    use dynamic_memory
    use dictionaries, dict_set => set
    !use wrapper_mpi, only: mpi_environment, MPI_COMM_WORLD
#if defined(HAVE_BPS)
    use Poisson_Solver, only: pkernel_init, pkernel_set
#endif
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    type(typ_poisson), intent(inout):: poisson
    !local variables
    character(len=1):: geocode
    integer:: n01, n02, n03, itype_scf, iproc=0, nproc=1
    integer:: nxyz(3), ndims(3)
    real(kind=8):: hgrids(3), hx, hy, hz
    real(kind=8):: cv1(3), cv2(3), cv3(3), ang_bc, ang_ac, ang_ab
    real(kind=8):: alpha_bc,beta_ac, gamma_ab ,pi=4.d0*atan(1.d0)
    type(dictionary), pointer :: dict_input=>null()
    !type(mpi_environment):: bigdft_mpi
#if defined(HAVE_BPS)
    write(*,*) 'REZA-1'
    !call f_lib_initialize() 
    write(*,*) 'REZA-2'
    !bigdft_mpi%mpi_comm=MPI_COMM_WORLD !workaround to be removed
    nxyz=(/64,64,64/)
    if(trim(atoms%boundcond)=='bulk') then
        geocode='P'
    elseif(trim(atoms%boundcond)=='slab') then
        geocode='S'
    elseif(trim(atoms%boundcond)=='free') then
        geocode='F'
    else
        write(*,*) 'ERROR: unknown atoms%boundcond in init_psolver_bps',trim(atoms%boundcond)
    endif
    !nxyz=options//'ndim'
    !geocode=options//'geocode'
    !call dict_free(options)
    n01=poisson%ngpx !nxyz(1)
    n02=poisson%ngpy !nxyz(2)
    n03=poisson%ngpz !nxyz(3)
    !Step size
    !order of the scaling functions chosen
    itype_scf=16
    !calculate the kernel in parallel for each processor
    ndims=(/n01,n02,n03/)
    hx=sqrt(sum(poisson%hgrid(1:3,1)**2))
    hy=sqrt(sum(poisson%hgrid(1:3,2)**2))
    hz=sqrt(sum(poisson%hgrid(1:3,3)**2))
    write(*,'(a,3f20.10)') 'norm: hx,hy,hz ',hx,hy,hz
    hgrids=(/hx,hy,hz/)
    cv1(1:3)=atoms%cellvec(1:3,1)
    cv2(1:3)=atoms%cellvec(1:3,2)
    cv3(1:3)=atoms%cellvec(1:3,3)
    ang_bc=acos(dot_product(cv2,cv3)/sqrt(dot_product(cv2,cv2)*dot_product(cv3,cv3)))
    ang_ac=acos(dot_product(cv1,cv3)/sqrt(dot_product(cv1,cv1)*dot_product(cv3,cv3)))
    ang_ab=acos(dot_product(cv1,cv2)/sqrt(dot_product(cv1,cv1)*dot_product(cv2,cv2)))
    !write(*,'(a,3f15.5)') 'alpha,beta,gamma ',ang_bc,ang_ac,ang_ab
    write(*,*) 'REZA-3'
    write(*,*) 'iproc,nproc', iproc, nproc
    write(*,*) 'geocode : ',geocode
    dict_input=>dict_new('kernel' .is. dict_new('isf_order' .is. itype_scf))
    alpha_bc = abs(ang_bc)!+pi/2.d0
    beta_ac = abs(ang_ac)!+pi/2.d0
    gamma_ab = abs(ang_ab)!+pi/2.d0
    !write(*,*) iproc,nproc,geocode,ndims, hgrids,alpha_bc,beta_ac,gamma_ab
    poisson%pkernel=pkernel_init(iproc,nproc,dict_input,geocode,ndims,hgrids,alpha_bc,beta_ac,gamma_ab)
    call dict_free(dict_input)
    write(*,*) 'REZA-4'
    call pkernel_set(poisson%pkernel,verbose=.true.)
    write(*,*) 'REZA-5'
    !write(*,'(a,2es20.12)') 'Hartree ',ehartree,ehartree-ehartree_ref
#else
    stop 'ERROR: Alborz is not linked with Poisson solvers in BigDFT.'
#endif
end subroutine init_psolver_bps
!*****************************************************************************************
subroutine fini_psolver_bps(poisson)
    use mod_electrostatics, only: typ_poisson
#if defined(HAVE_BPS)
    use wrapper_mpi, only: mpi_environment, MPI_COMM_WORLD
    use Poisson_Solver, only: pkernel_free
#endif
    implicit none
    type(typ_poisson), intent(inout):: poisson
#if defined(HAVE_BPS)
    call pkernel_free(poisson%pkernel)
    !call f_lib_finalize()
#else
    stop 'ERROR: Alborz is not linked with Poisson solvers in BigDFT.'
#endif
end subroutine fini_psolver_bps
!*****************************************************************************************
subroutine set_ngp_bps(parini,atoms,poisson_rough,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
#if defined(HAVE_BPS)
    use module_fft_sg, only: i_data, ndata
#endif
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    type(typ_poisson), intent(in):: poisson_rough
    type(typ_poisson), intent(inout):: poisson
    !local variables
    real(8):: dh1, dh2, harr(3)
    real(8):: cell(3), vol, cvinv(3), cvinv_norm(3)
    integer:: i, ndim(3), id
    call cell_vol(atoms%nat,atoms%cellvec,vol)
    vol=abs(vol)*atoms%nat
    call cross_product_alborz(atoms%cellvec(1,1),atoms%cellvec(1,2),cvinv)
    cvinv_norm(3)=sqrt(cvinv(1)**2+cvinv(2)**2+cvinv(3)**2)
    cell(3)=vol/cvinv_norm(3)
    call cross_product_alborz(atoms%cellvec(1,2),atoms%cellvec(1,3),cvinv)
    cvinv_norm(1)=sqrt(cvinv(1)**2+cvinv(2)**2+cvinv(3)**2)
    cell(1)=vol/cvinv_norm(1)
    call cross_product_alborz(atoms%cellvec(1,1),atoms%cellvec(1,3),cvinv)
    cvinv_norm(2)=sqrt(cvinv(1)**2+cvinv(2)**2+cvinv(3)**2)
    cell(2)=vol/cvinv_norm(2)
    !write(*,*) cell(1:3)
    poisson%ngpx=int(cell(1)/poisson_rough%hx)+1
    poisson%ngpy=int(cell(2)/poisson_rough%hy)+1
    poisson%ngpz=int(cell(3)/poisson_rough%hz)+1
    !write(*,*) poisson_rough%hx,poisson_rough%hy,poisson_rough%hz
    !write(*,*) poisson%ngpx,poisson%ngpy,poisson%ngpz

#if defined(HAVE_BPS)

    ndim(1)=poisson%ngpx
    ndim(2)=poisson%ngpy
    ndim(3)=poisson%ngpz
    harr(1)=poisson_rough%hx
    harr(2)=poisson_rough%hy
    harr(3)=poisson_rough%hz
    do id=1,3
        do i=1,ndata
            if(ndim(id)<i_data(i)) exit
        enddo
        if(i==ndata .or. i==ndata+1) then 
            write(*,*) 'ERROR: invalid number of grid points in calparam: ',id,ndim(id)
            stop
        endif
        if(ndim(id)/=i_data(i-1)) then
            dh1=cell(id)/real(i_data(i-1),8)
            dh2=cell(id)/real(i_data(i),8)
            if(abs(dh1-harr(id))<abs(dh2-harr(id))) then
                ndim(id)=i_data(i-1)
            else
                ndim(id)=i_data(i)
            endif
        endif
    enddo
    poisson%ngpx=max(16,ndim(1))
    poisson%ngpy=max(16,ndim(2))
    poisson%ngpz=max(16,ndim(3))
    !write(*,*) ndim(:)
    !stop
#else
    stop 'ERROR: FLAME is not linked with Poisson solvers in BigDFT.'
#endif
end subroutine set_ngp_bps
!*****************************************************************************************
