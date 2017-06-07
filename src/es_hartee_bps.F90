!*****************************************************************************************
subroutine cal_hartree_pot_bps(ewald_p3d,atoms,ehartree)
    use mod_interface
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_ewald_p3d
#if defined(HAVE_BPS)
    use Poisson_Solver, only: H_POTENTIAL
#endif
    use dynamic_memory
    implicit none
    type(typ_ewald_p3d),intent(inout):: ewald_p3d
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(out):: ehartree
    !local variables
    integer:: igpx, igpy, igpz
#if defined(HAVE_BPS)
    real(8), allocatable:: pot_ion(:)
    real(8):: stress_t(6)
    !type(coulomb_operator):: pkernel
    pot_ion=f_malloc([1.to.1],id='pot_ion') !an array which must be provided but will not be used.
#endif

#if defined(HAVE_BPS)
    !stress_t(1:6)=0.d0
    call H_POTENTIAL('G',ewald_p3d%poisson_p3d%pkernel,ewald_p3d%poisson_p3d%rho, &
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
    do igpz=1,ewald_p3d%poisson_p3d%ngpz
    do igpy=1,ewald_p3d%poisson_p3d%ngpy
    do igpx=1,ewald_p3d%poisson_p3d%ngpx
        ewald_p3d%poisson_p3d%pot(igpx,igpy,igpz)=ewald_p3d%poisson_p3d%rho(igpx,igpy,igpz)
    enddo
    enddo
    enddo
    call f_free(pot_ion)
#else
    stop 'ERROR: Alborz is not linked with Poisson solvers in BigDFT.'
    ehartree=0.d0 !this is just to be able to compile.
#endif
end subroutine cal_hartree_pot_bps
!*****************************************************************************************
subroutine construct_ewald_bps(parini,atoms,ewald_p3d)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_ewald_p3d
    use dynamic_memory
    use dictionaries, dict_set => set
    !use wrapper_mpi, only: mpi_environment, MPI_COMM_WORLD
#if defined(HAVE_BPS)
    use Poisson_Solver, only: pkernel_init, pkernel_set
#endif
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
    !local variables
    character(len=1):: geocode
    integer:: n01, n02, n03, itype_scf, iproc=0, nproc=1
    integer:: nxyz(3), ndims(3)
    real(kind=8):: hgrids(3)
    real(kind=8):: cv1(3), cv2(3), cv3(3), ang_bc, ang_ac, ang_ab
    type(dictionary), pointer :: dict_input
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
        write(*,*) 'ERROR: unknown atoms%boundcond in construct_ewald_bps',trim(atoms%boundcond)
    endif
    !nxyz=options//'ndim'
    !geocode=options//'geocode'
    !call dict_free(options)
    n01=ewald_p3d%poisson_p3d%ngpx !nxyz(1)
    n02=ewald_p3d%poisson_p3d%ngpy !nxyz(2)
    n03=ewald_p3d%poisson_p3d%ngpz !nxyz(3)
    !Step size
    !ewald_p3d%hgx=atoms%cellvec(1,1)/real(n01,kind=8)
    !ewald_p3d%hgy=atoms%cellvec(2,2)/real(n02,kind=8)
    !ewald_p3d%hgz=atoms%cellvec(3,3)/real(n03,kind=8)
    write(*,'(a,3f15.10)') 'hx,hy,hz ',ewald_p3d%hgx,ewald_p3d%hgy,ewald_p3d%hgz
    !order of the scaling functions chosen
    itype_scf=16
    !calculate the kernel in parallel for each processor
    ndims=(/n01,n02,n03/)
    hgrids=(/ewald_p3d%hgx,ewald_p3d%hgy,ewald_p3d%hgz/)
    cv1(1:3)=atoms%cellvec(1:3,1)
    cv2(1:3)=atoms%cellvec(1:3,2)
    cv3(1:3)=atoms%cellvec(1:3,3)
    ang_bc=acos(dot_product(cv2,cv3)/sqrt(dot_product(cv2,cv2)*dot_product(cv3,cv3)))
    ang_ac=acos(dot_product(cv1,cv3)/sqrt(dot_product(cv1,cv1)*dot_product(cv3,cv3)))
    ang_ab=acos(dot_product(cv1,cv2)/sqrt(dot_product(cv1,cv1)*dot_product(cv2,cv2)))
    !write(*,'(a,3f15.5)') 'alpha,beta,gamma ',ang_bc,ang_ac,ang_ab
    write(*,*) 'REZA-3'
    write(*,*) iproc, nproc
    write(*,*) geocode
    dict_input=>dict_new('kernel' .is. dict_new('isf_order' .is. itype_scf))
    ewald_p3d%poisson_p3d%pkernel=pkernel_init(iproc,nproc,dict_input,geocode,ndims, &
        hgrids,alpha_bc=ang_bc,beta_ac=ang_ac,gamma_ab=ang_ab)
    call dict_free(dict_input)
    write(*,*) 'REZA-4'
    call pkernel_set(ewald_p3d%poisson_p3d%pkernel,verbose=.true.)
    write(*,*) 'REZA-5'
    !write(*,'(a,2es20.12)') 'Hartree ',ehartree,ehartree-ehartree_ref
#else
    stop 'ERROR: Alborz is not linked with Poisson solvers in BigDFT.'
#endif
end subroutine construct_ewald_bps
!*****************************************************************************************
subroutine destruct_ewald_bps(ewald_p3d)
    use mod_interface
    use mod_electrostatics, only: typ_ewald_p3d
#if defined(HAVE_BPS)
    use wrapper_mpi, only: mpi_environment, MPI_COMM_WORLD
    use Poisson_Solver, only: pkernel_free
#endif
    implicit none
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
#if defined(HAVE_BPS)
    call pkernel_free(ewald_p3d%poisson_p3d%pkernel)
    !call f_lib_finalize()
#else
    stop 'ERROR: Alborz is not linked with Poisson solvers in BigDFT.'
#endif
end subroutine destruct_ewald_bps
!*****************************************************************************************
subroutine set_ngp_bps(atoms,ewald_p3d_rough,ewald_p3d)
    use mod_interface
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_ewald_p3d
#if defined(HAVE_BPS)
    use module_fft_sg, only: i_data, ndata
#endif
    implicit none
    type(typ_atoms), intent(in):: atoms
    type(typ_ewald_p3d), intent(in):: ewald_p3d_rough
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
    !local variables
    real(8):: dh1, dh2, harr(3)
    real(8):: cell(3), vol, cvinv(3), cvinv_norm
    integer:: i, ndim(3), id
    call cell_vol(atoms%nat,atoms%cellvec,vol)
    vol=abs(vol)*atoms%nat
    call cross_product_alborz(atoms%cellvec(1,1),atoms%cellvec(1,2),cvinv)
    cvinv_norm=sqrt(cvinv(1)**2+cvinv(2)**2+cvinv(3)**2)
    cell(3)=vol/cvinv_norm
    call cross_product_alborz(atoms%cellvec(1,2),atoms%cellvec(1,3),cvinv)
    cvinv_norm=sqrt(cvinv(1)**2+cvinv(2)**2+cvinv(3)**2)
    cell(1)=vol/cvinv_norm
    call cross_product_alborz(atoms%cellvec(1,1),atoms%cellvec(1,3),cvinv)
    cvinv_norm=sqrt(cvinv(1)**2+cvinv(2)**2+cvinv(3)**2)
    cell(2)=vol/cvinv_norm
    !write(*,*) cell(1:3)
    ewald_p3d%poisson_p3d%ngpx=int(cell(1)/ewald_p3d_rough%hgx)+1
    ewald_p3d%poisson_p3d%ngpy=int(cell(2)/ewald_p3d_rough%hgy)+1
    ewald_p3d%poisson_p3d%ngpz=int(cell(3)/ewald_p3d_rough%hgz)+1

#if defined(HAVE_BPS)
    ndim(1)=ewald_p3d%poisson_p3d%ngpx
    ndim(2)=ewald_p3d%poisson_p3d%ngpy
    ndim(3)=ewald_p3d%poisson_p3d%ngpz
    harr(1)=ewald_p3d_rough%hgx
    harr(2)=ewald_p3d_rough%hgy
    harr(3)=ewald_p3d_rough%hgz
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
    ewald_p3d%poisson_p3d%ngpx=max(20,ndim(1))
    ewald_p3d%poisson_p3d%ngpy=max(20,ndim(2))
    ewald_p3d%poisson_p3d%ngpz=max(20,ndim(3))
    !write(*,*) ndim(:)
    !stop
    ewald_p3d%hgx=sqrt(sum(atoms%cellvec(1:3,1)**2))/real(ewald_p3d%poisson_p3d%ngpx,8)
    ewald_p3d%hgy=sqrt(sum(atoms%cellvec(1:3,2)**2))/real(ewald_p3d%poisson_p3d%ngpy,8)
    ewald_p3d%hgz=sqrt(sum(atoms%cellvec(1:3,3)**2))/real(ewald_p3d%poisson_p3d%ngpz,8)
    !ewald_p3d%hgx=cell(1)/real(ewald_p3d%poisson_p3d%ngpx,8)
    !ewald_p3d%hgy=cell(2)/real(ewald_p3d%poisson_p3d%ngpy,8)
    !ewald_p3d%hgz=cell(3)/real(ewald_p3d%poisson_p3d%ngpz,8)
#else
    stop 'ERROR: Alborz is not linked with Poisson solvers in BigDFT.'
#endif
end subroutine set_ngp_bps
!*****************************************************************************************
