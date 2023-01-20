!*****************************************************************************************
module mod_trial_energy
    implicit none
    private
    public:: trial_energy_allocate, trial_energy_deallocate
    public:: get_trial_energy, get_rmse
    type, public:: typ_trial_energy
        private
        integer, public:: ntrial=0
        real(8), public:: ehartree_scn_excl
        real(8), allocatable, public:: energy(:)
        real(8), allocatable, public:: disp(:,:)
        integer, allocatable, public:: iat_list(:)
        real(8), allocatable, public:: EP(:,:)
        real(8), allocatable, public:: E_all(:)
        real(8), allocatable, public:: EP_n(:)

    end type typ_trial_energy
contains
!*****************************************************************************************
subroutine trial_energy_allocate(ntrial,trial_energy,nbf)
    use dynamic_memory
    implicit none
    integer, intent(in):: ntrial
    type(typ_trial_energy), intent(out):: trial_energy
    integer, intent(in):: nbf
    !local variables
    integer:: istat
    if(allocated(trial_energy%energy)) stop 'ERROR: trial_energy%energy is already allocated'
    if(allocated(trial_energy%disp)) stop 'ERROR: trial_energy%disp is already allocated'
    if(allocated(trial_energy%iat_list)) stop 'ERROR: trial_energy%iat_list is already allocated'
    trial_energy%ntrial=ntrial
    trial_energy%energy=f_malloc0([1.to.ntrial],id='trial_energy%energy')
    trial_energy%iat_list=f_malloc0([1.to.ntrial],id='trial_energy%iat_list')
    trial_energy%disp=f_malloc0([1.to.3,1.to.ntrial],id='trial_energy%disp')
    trial_energy%EP_n=f_malloc0([1.to.ntrial],id='trial_energy%EP_n')
    trial_energy%E_all=f_malloc0([1.to.ntrial],id='trial_energy%E_all')
    trial_energy%EP=f_malloc0([1.to.nbf,1.to.ntrial],id='trial_energy%EP')
end subroutine trial_energy_allocate
!*****************************************************************************************
subroutine trial_energy_deallocate(trial_energy)
    use dynamic_memory
    implicit none
    type(typ_trial_energy), intent(inout):: trial_energy
    !local variables
    integer:: istat
    if(allocated(trial_energy%energy)) call f_free(trial_energy%energy)
    if(allocated(trial_energy%disp)) call f_free(trial_energy%disp)
    if(allocated(trial_energy%iat_list)) call f_free(trial_energy%iat_list)
    if(allocated(trial_energy%EP_n)) call f_free(trial_energy%EP_n)
    if(allocated(trial_energy%E_all)) call f_free(trial_energy%E_all)
    if(allocated(trial_energy%EP)) call f_free(trial_energy%EP)
end subroutine trial_energy_deallocate
!*****************************************************************************************
subroutine get_trial_energy(parini,atoms,poisson,nbf,bz,gwz,trial_energy,qtot,dpm)
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms, get_rat, update_ratp, atom_deallocate_old, update_rat
    use mod_processors, only: get_proc_stake
    use mod_flm_futile
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_poisson), intent(inout):: poisson
    integer, intent(in):: nbf
    real(8), intent(in):: bz(atoms%nat), gwz(atoms%nat)
    type(typ_trial_energy), intent(inout):: trial_energy
    real(8), intent(out):: qtot, dpm(3)
    !local variables
    type(typ_poisson):: poisson_ion
    integer:: igpx, igpy, igpz, iat, ntrial, itrial, nsegx, nsegy, nsegz
    integer:: jgpx, jgpy, jgpz
    real(8):: rgcut_a, pi!, qtot_e, qtot_i
    real(8):: ehartree_scn_excl, tt1, tt2
    real(8):: xyz(3), dxyz(3), epot_trial, gwt, q_tmp(1), center(3)
    real(8):: dx, dy, dz, r2, coeff, x, y, z !, rloc, c1, c2
    real(8):: xmin, ymin, zmin, xmax, ymax, zmax
    real(8):: q_one(1), gw_one(1)
    real(8), allocatable::  gausswidth(:)
    real(8), allocatable::  rat_trial(:,:)
    real(8), allocatable::  rho(:,:,:)
    real(8), allocatable::  pot(:,:,:)
    integer:: nbgpx, nbgpy, nbgpz, ix, iy, iz
    integer:: nex, ney, nez, nd
    integer:: itrials, itriale, ierr
#if defined(MPI)
    include 'mpif.h'
#endif
    pi=4.d0*atan(1.d0)
    if(parini%ewald) then
        write(*,*) 'ERROR: ewald=True is wrong when reading from cube file.'
        stop
    endif
    if(trim(parini%psolver)=='kwald') then
        write(*,*) 'ERROR: psolver=kwald is wrong for grid base charge density.'
        stop
    endif
    if(parini%cal_scn .and. parini%screening_factor==0.d0) then
        stop 'ERROR: cal_scn is TRUE and screening_factor is 0, MEANINGLESS!'
    endif
    !-------------------------------------------------------
    do iat=1,atoms%nat
        if(trim(atoms%sat(iat))=='Li') atoms%zat(iat)=3.d0
        if(trim(atoms%sat(iat))=='S' ) atoms%zat(iat)=6.d0
    enddo
    !-------------------------------------------------------
    if(.true.) then
    nd=2
    allocate(rho(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    do igpz=1,poisson%ngpz
    do igpy=1,poisson%ngpy
    do igpx=1,poisson%ngpx
        rho(igpx,igpy,igpz)=poisson%rho(igpx,igpy,igpz)
    enddo
    enddo
    enddo
    call f_free(poisson%rho)
    poisson%ngpz=poisson%ngpz/nd
    poisson%ngpy=poisson%ngpy/nd
    poisson%ngpx=poisson%ngpx/nd
    poisson%rho=f_malloc0([1.to.(poisson%ngpx),1.to.(poisson%ngpy),1.to.(poisson%ngpz)],id='poisson%rho')
    do igpz=1,poisson%ngpz
    do igpy=1,poisson%ngpy
    do igpx=1,poisson%ngpx
        poisson%rho(igpx,igpy,igpz)=rho(nd*igpx-nd+1,nd*igpy-nd+1,nd*igpz-nd+1)
    enddo
    enddo
    enddo
    deallocate(rho)
    poisson%hgrid=poisson%hgrid*real(nd,kind=8)
    endif
    !-------------------------------------------------------
    allocate(rho(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    do igpz=1,poisson%ngpz
    do igpy=1,poisson%ngpy
    do igpx=1,poisson%ngpx
        rho(igpx,igpy,igpz)=poisson%rho(igpx,igpy,igpz)
    enddo
    enddo
    enddo
    !write(*,*) allocated(poisson%rho)
    call f_free(poisson%rho)
    nex=20
    ney=20
    nez=20
    poisson%rho=f_malloc0([1.to.(poisson%ngpx+2*nex),1.to.(poisson%ngpy+2*ney),1.to.(poisson%ngpz+2*nez)],id='poisson%rho')
    do igpz=1,poisson%ngpz
    do igpy=1,poisson%ngpy
    do igpx=1,poisson%ngpx
        poisson%rho(igpx+nex,igpy+ney,igpz+nez)=rho(igpx,igpy,igpz)
    enddo
    enddo
    enddo
    call update_ratp(atoms)
    do iat=1,atoms%nat
        atoms%ratp(1,iat)=atoms%ratp(1,iat)+real(nex,kind=8)*poisson%hgrid(1,1)
        atoms%ratp(2,iat)=atoms%ratp(2,iat)+real(ney,kind=8)*poisson%hgrid(2,2)
        atoms%ratp(3,iat)=atoms%ratp(3,iat)+real(nez,kind=8)*poisson%hgrid(3,3)
    enddo
    call update_rat(atoms)
    deallocate(rho)
    poisson%ngpx=poisson%ngpx+2*nex
    poisson%ngpy=poisson%ngpy+2*ney
    poisson%ngpz=poisson%ngpz+2*nez
    !-------------------------------------------------------
    atoms%boundcond='free'
    poisson%bc=atoms%boundcond
    allocate(gausswidth(atoms%nat))
    gausswidth=0.5d0 !parini%gaussian_width !TO_BE_CORRECTED
    poisson%cal_scn=parini%cal_scn
    poisson%screening_factor=parini%screening_factor
    poisson%task_finit=""
    call init_hartree(parini,atoms,poisson,gausswidth)
    poisson%cell(1)=poisson%hgrid(1,1)*poisson%ngpx
    poisson%cell(2)=poisson%hgrid(2,2)*poisson%ngpy
    poisson%cell(3)=poisson%hgrid(3,3)*poisson%ngpz
    !-------------------------------------------------------
    poisson_ion%alpha=0.5d0 !parini%gaussian_width !TO_BE_CORRECTED
    poisson_ion%rgcut=parini%rgcut_ewald*poisson_ion%alpha
    poisson_ion%ngpx=poisson%ngpx
    poisson_ion%ngpy=poisson%ngpy
    poisson_ion%ngpz=poisson%ngpz
    poisson_ion%hgrid(1:3,1:3)=0.d0
    poisson_ion%hgrid(1,1)=poisson%hgrid(1,1)
    poisson_ion%hgrid(2,2)=poisson%hgrid(2,2)
    poisson_ion%hgrid(3,3)=poisson%hgrid(3,3)
    poisson_ion%xyz111=poisson%xyz111
    rgcut_a=parini%rgcut_ewald*maxval(gausswidth) !parini%gaussian_width !3.d0
    nbgpx=int(rgcut_a/poisson_ion%hgrid(1,1))+2
    nbgpy=int(rgcut_a/poisson_ion%hgrid(2,2))+2
    nbgpz=int(rgcut_a/poisson_ion%hgrid(3,3))+2
    poisson_ion%task_finit="alloc_rho"
    call init_hartree(parini,atoms,poisson_ion,gausswidth)
    poisson_ion%reset_rho=.true.
    poisson_ion%nat=atoms%nat
    poisson_ion%cv=atoms%cellvec
    poisson_ion%bc=atoms%boundcond
    poisson_ion%q(1:poisson_ion%nat)=atoms%zat(1:atoms%nat)
    poisson_ion%gw(1:poisson_ion%nat)=gausswidth(1:atoms%nat)
    call get_rat(atoms,poisson_ion%rcart)
    !call put_charge_density(parini,poisson_ion)
!    nbgpx=int(7.d0*maxval(gausswidth(1:atoms%nat))/poisson_ion%hgrid(1,1))+2
!    nbgpy=int(7.d0*maxval(gausswidth(1:atoms%nat))/poisson_ion%hgrid(2,2))+2
!    nbgpz=int(7.d0*maxval(gausswidth(1:atoms%nat))/poisson_ion%hgrid(3,3))+2
    poisson_ion%rho=0.d0
    do iat=1,atoms%nat
    !itypat=atoms%itypat(iat)
    q_tmp(1)=atoms%zat(iat)*bz(iat)
    call put_gto_sym_ortho(parini,poisson_ion%bc,.false.,1,atoms%ratp(1,iat),q_tmp,gwz(iat), &
        5.d0*gwz(iat),poisson_ion%xyz111,poisson_ion%ngpx,poisson_ion%ngpy,poisson_ion%ngpz, &
        poisson_ion%hgrid,poisson_ion%rho)
    q_tmp(1)=atoms%zat(iat)*(1.d0-bz(iat))
    call put_r2gto_sym_ortho(parini,poisson_ion%bc,.false.,1,atoms%ratp(1,iat),q_tmp,gwz(iat), &
        5.d0*gwz(iat),poisson_ion%xyz111,poisson_ion%ngpx,poisson_ion%ngpy,poisson_ion%ngpz, &
        poisson_ion%hgrid,poisson_ion%rho)
    !if(trim(atoms%sat(iat))=='Mg') gwt=0.6d0
    !if(trim(atoms%sat(iat))=='O' ) gwt=0.3d0
!    gwt=gausswidth(iat)
!    coeff=atoms%zat(iat)/(gwt**3*pi**1.5d0)
    !coeff=atoms%zat(iat)*2.d0/(3.d0*gwt**5*pi**1.5d0)
    !coeff=atoms%zat(iat)*4.d0/(15.d0*gwt**7*pi**1.5d0)
    !if(trim(atoms%sat(iat))=='Mg') then
    !    !0.65406138674  2 -5.223929095  0.913704167481045 rloc nloc c1 .. cnloc
    !    rloc=0.65406138674d0
    !    c1=-5.223929095d0
    !    c2=0.913704167481045d0
    !endif
    !if(trim(atoms%sat(iat))=='O') then
    !    !0.3454999999    2 -11.7435870154  1.90653967947 rloc nloc c1 .. cnloc
    !    rloc=0.3454999999d0
    !    c1=-11.7435870154d0
    !    c2=1.90653967947d0
    !endif
    !c1=-c1
    !c2=-c2
!    jgpx=int(poisson_ion%rcart(1,iat)/poisson%hgrid(1,1))
!    jgpy=int(poisson_ion%rcart(2,iat)/poisson%hgrid(2,2))
!    jgpz=int(poisson_ion%rcart(3,iat)/poisson%hgrid(3,3))
!    do igpz=jgpz-nbgpz,jgpz+nbgpz
!    do igpy=jgpy-nbgpy,jgpy+nbgpy
!    do igpx=jgpx-nbgpx,jgpx+nbgpx
!        dx=(igpx-1)*poisson%hgrid(1,1)-poisson_ion%rcart(1,iat)
!        dy=(igpy-1)*poisson%hgrid(2,2)-poisson_ion%rcart(2,iat)
!        dz=(igpz-1)*poisson%hgrid(3,3)-poisson_ion%rcart(3,iat)
!        r2=dx**2+dy**2+dz**2
!        if(r2<10.d0**2*gwt**2) then
!        !if(r2<10.d0**2*rloc**2) then
!            poisson_ion%rho(igpx,igpy,igpz)=poisson_ion%rho(igpx,igpy,igpz)+coeff*exp(-r2/gwt**2)
!            !poisson_ion%rho(igpx,igpy,igpz)=poisson_ion%rho(igpx,igpy,igpz)+coeff*r2*exp(-r2/gwt**2)
!            !poisson_ion%rho(igpx,igpy,igpz)=poisson_ion%rho(igpx,igpy,igpz)+coeff*r2**2*exp(-r2/gwt**2)
!            !tt1=exp(-r2/(2.d0*rloc**2))
!            !tt2=-3.d0*rloc**4*(c1-2.d0*c2)+rloc**2*(c1-7.d0*c2)*r2+c2*r2**2+rloc**3*sqrt(2.d0/pi)*atoms%zat(iat)
!            !poisson_ion%rho(igpx,igpy,igpz)=poisson_ion%rho(igpx,igpy,igpz)+tt1*tt2/(4.d0*rloc**6*pi)
!        endif
!    enddo
!    enddo
!    enddo
    enddo
    tt1=0.d0
    do igpz=1,poisson%ngpz
    do igpy=1,poisson%ngpy
    do igpx=1,poisson%ngpx
        tt1=tt1+poisson_ion%rho(igpx,igpy,igpz)
    enddo
    enddo
    enddo
    tt1=tt1*(poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3))
    if(parini%mpi_env%iproc==0) then
    write(*,*) 'TT1 ',tt1
    endif
    !-------------------------------------------------------
    qtot=0.d0
    do igpz=1,poisson%ngpz
        do igpy=1,poisson%ngpy
            do igpx=1,poisson%ngpx
                poisson%rho(igpx,igpy,igpz)=poisson_ion%rho(igpx,igpy,igpz)-poisson%rho(igpx,igpy,igpz)
                qtot=qtot+poisson%rho(igpx,igpy,igpz)
            enddo
        enddo
    enddo
    qtot=qtot*poisson_ion%hgrid(1,1)*poisson_ion%hgrid(2,2)*poisson_ion%hgrid(3,3)
    if(parini%mpi_env%iproc==0) then
        write(*,'(a,f20.12)') 'qtot= ',qtot
    endif
    !-------------------------------------------------------
    center(1)=0.d0
    center(2)=0.d0
    center(3)=0.d0
    do iat=1,atoms%nat
        center(1)=center(1)+atoms%ratp(1,iat)
        center(2)=center(2)+atoms%ratp(2,iat)
        center(3)=center(3)+atoms%ratp(3,iat)
    enddo
    center(1)=center(1)/atoms%nat
    center(2)=center(2)/atoms%nat
    center(3)=center(3)/atoms%nat
    dpm(1)=0.d0
    dpm(2)=0.d0
    dpm(3)=0.d0
    do igpz=1,poisson%ngpz
    do igpy=1,poisson%ngpy
    do igpx=1,poisson%ngpx
        x=poisson%xyz111(1)+(igpx-1)*poisson%hgrid(1,1)-center(1)
        y=poisson%xyz111(2)+(igpy-1)*poisson%hgrid(2,2)-center(2)
        z=poisson%xyz111(3)+(igpz-1)*poisson%hgrid(3,3)-center(3)
        dpm(1)=dpm(1)+x*poisson%rho(igpx,igpy,igpz)
        dpm(2)=dpm(2)+y*poisson%rho(igpx,igpy,igpz)
        dpm(3)=dpm(3)+z*poisson%rho(igpx,igpy,igpz)
    enddo
    enddo
    enddo
    dpm(1)=dpm(1)*(poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3))
    dpm(2)=dpm(2)*(poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3))
    dpm(3)=dpm(3)*(poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3))
    if(parini%mpi_env%iproc==0) then
        write(*,'(a,3f10.5)') 'DPM= ',dpm(1),dpm(2),dpm(3)
    endif
    !-------------------------------------------------------
    call update_ratp(atoms)
    if(parini%mpi_env%iproc==0) then
        call get_hartree(parini,poisson,atoms,gausswidth,ehartree_scn_excl)
    endif
    !write(*,*) 'ALLOCATED= ',parini%mpi_env%iproc,allocated(poisson%pot)
    if(parini%mpi_env%nproc>1) then
    call MPI_BARRIER(parini%mpi_env%mpi_comm,ierr)
    call MPI_BCAST(poisson%pot,poisson%ngpx*poisson%ngpy*poisson%ngpz,MPI_DOUBLE_PRECISION,0,parini%mpi_env%mpi_comm,ierr)
    endif
    if(parini%mpi_env%iproc==0) then
    write(*,'(a,es24.15,es14.5)') 'ehartree_scn_excl ',ehartree_scn_excl,poisson%screening_factor
    endif
    !-------------------------------------------------------
    if(.true.) then
    nd=1
    allocate(rho(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    allocate(pot(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    do igpz=1,poisson%ngpz
    do igpy=1,poisson%ngpy
    do igpx=1,poisson%ngpx
        rho(igpx,igpy,igpz)=poisson%rho(igpx,igpy,igpz)
        pot(igpx,igpy,igpz)=poisson%pot(igpx,igpy,igpz)
    enddo
    enddo
    enddo
    call f_free(poisson%rho)
    call f_free(poisson%pot)
    poisson%ngpz=poisson%ngpz/nd
    poisson%ngpy=poisson%ngpy/nd
    poisson%ngpx=poisson%ngpx/nd
    poisson%rho=f_malloc0([1.to.(poisson%ngpx),1.to.(poisson%ngpy),1.to.(poisson%ngpz)],id='poisson%rho')
    poisson%pot=f_malloc0([1.to.(poisson%ngpx),1.to.(poisson%ngpy),1.to.(poisson%ngpz)],id='poisson%pot')
    do igpz=1,poisson%ngpz
    do igpy=1,poisson%ngpy
    do igpx=1,poisson%ngpx
        poisson%rho(igpx,igpy,igpz)=rho(nd*igpx-nd+1,nd*igpy-nd+1,nd*igpz-nd+1)
        poisson%pot(igpx,igpy,igpz)=pot(nd*igpx-nd+1,nd*igpy-nd+1,nd*igpz-nd+1)
    enddo
    enddo
    enddo
    deallocate(rho)
    deallocate(pot)
    poisson%hgrid=poisson%hgrid*real(nd,kind=8)
    endif
    !-------------------------------------------------------
    atoms%fat=0.d0
    call force_gto_sym_ortho(parini,poisson_ion%bc,atoms%nat,poisson_ion%rcart, &
        poisson_ion%q,gausswidth,6.d0,poisson_ion%xyz111, &
        poisson_ion%ngpx,poisson_ion%ngpx,poisson_ion%ngpy,poisson_ion%ngpz, &
        poisson_ion%hgrid,poisson%pot,atoms%fat)
    if(parini%mpi_env%iproc==0) then
    do iat=1,atoms%nat
        write(*,'(a,i4,3es19.10)') 'FAT ',iat,atoms%fat(1,iat),atoms%fat(2,iat),atoms%fat(3,iat)
    enddo
    endif
    !-------------------------------------------------------
    poisson_ion%gw(1:poisson_ion%nat)=1.d0
    xmin= huge(1.d0)
    ymin= huge(1.d0)
    zmin= huge(1.d0)
    xmax=-huge(1.d0)
    ymax=-huge(1.d0)
    zmax=-huge(1.d0)
    do iat=1,atoms%nat
        if(poisson_ion%rcart(1,iat)<xmin) xmin=poisson_ion%rcart(1,iat)
        if(poisson_ion%rcart(2,iat)<ymin) ymin=poisson_ion%rcart(2,iat)
        if(poisson_ion%rcart(3,iat)<zmin) zmin=poisson_ion%rcart(3,iat)
        if(poisson_ion%rcart(1,iat)>xmax) xmax=poisson_ion%rcart(1,iat)
        if(poisson_ion%rcart(2,iat)>ymax) ymax=poisson_ion%rcart(2,iat)
        if(poisson_ion%rcart(3,iat)>zmax) zmax=poisson_ion%rcart(3,iat)
    enddo
    dxyz(1)=1.5d0
    dxyz(2)=1.5d0
    dxyz(3)=1.5d0
    nsegx=int((poisson%hgrid(1,1)*poisson%ngpx-16.d0)/dxyz(1))+1
    nsegy=int((poisson%hgrid(2,2)*poisson%ngpy-16.d0)/dxyz(2))+1
    nsegz=int((poisson%hgrid(3,3)*poisson%ngpz-16.d0)/dxyz(3))+1
    dxyz(1)=(poisson%hgrid(1,1)*poisson%ngpx-16.d0)/real(nsegx,kind=8)
    dxyz(2)=(poisson%hgrid(2,2)*poisson%ngpy-16.d0)/real(nsegy,kind=8)
    dxyz(3)=(poisson%hgrid(3,3)*poisson%ngpz-16.d0)/real(nsegz,kind=8)
    ntrial=(nsegx+1)*(nsegy+1)*(nsegz+1)
    call trial_energy_allocate(ntrial,trial_energy,nbf)
    trial_energy%ehartree_scn_excl=ehartree_scn_excl
    if(parini%mpi_env%iproc==0) then
    write(*,'(a,3i3,i6)') 'nsegx,nsegy,nsegz,ntrial ',nsegx,nsegy,nsegz,ntrial
    endif
    allocate(rat_trial(3,ntrial))
    itrial=0
    do iz=0,nsegz
    do iy=0,nsegy
    do ix=0,nsegx
        itrial=itrial+1
        rat_trial(1,itrial)=8.d0+dxyz(1)*ix
        rat_trial(2,itrial)=8.d0+dxyz(2)*iy
        rat_trial(3,itrial)=8.d0+dxyz(3)*iz
    enddo
    enddo
    enddo
    call get_proc_stake(parini%mpi_env,ntrial,itrials,itriale)
    !write(*,'(a,4i8)') 'iproc,itrials,itriale,ntrial ',parini%mpi_env%iproc,itrials,itriale,ntrial
    trial_energy%energy=0.d0
    trial_energy%disp=0.d0
    trial_energy%iat_list=0
    do itrial=itrials,itriale
        xyz(1)=rat_trial(1,itrial)-poisson_ion%rcart(1,1)
        xyz(2)=rat_trial(2,itrial)-poisson_ion%rcart(2,1)
        xyz(3)=rat_trial(3,itrial)-poisson_ion%rcart(3,1)
        q_one(1)=1.d0
        gw_one(1)=1.d0
        call put_gto_sym_ortho(parini,poisson%bc,.true.,1,rat_trial(1,itrial),q_one,gw_one, &
            5.d0,poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
        epot_trial=0.d0
        do igpz=1,poisson%ngpz
        do igpy=1,poisson%ngpy
        do igpx=1,poisson%ngpx
            epot_trial=epot_trial+poisson%rho(igpx,igpy,igpz)*poisson%pot(igpx,igpy,igpz)
        enddo
        enddo
        enddo
        epot_trial=epot_trial*(poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3))
        trial_energy%energy(itrial)=epot_trial
        trial_energy%disp(1,itrial)=xyz(1)
        trial_energy%disp(2,itrial)=xyz(2)
        trial_energy%disp(3,itrial)=xyz(3)
        trial_energy%iat_list(itrial)=1
    enddo
    !-------------------------------------------------------
    do iat=1,atoms%nat
        if(trim(atoms%sat(iat))=='Li') atoms%zat(iat)=0.d0
        if(trim(atoms%sat(iat))=='S' ) atoms%zat(iat)=0.d0
    enddo
    !-------------------------------------------------------
    if(parini%mpi_env%nproc>1) then
    call fmpi_allreduce(trial_energy%energy(1),ntrial,op=FMPI_SUM,comm=parini%mpi_env%mpi_comm)
    call fmpi_allreduce(trial_energy%disp(1,1),3*ntrial,op=FMPI_SUM,comm=parini%mpi_env%mpi_comm)
    call fmpi_allreduce(trial_energy%iat_list(1),ntrial,op=FMPI_SUM,comm=parini%mpi_env%mpi_comm)
    endif
    if(parini%mpi_env%iproc==0) then
    write(*,'(a,6f8.1)') 'MINMAX ',xmin,ymin,zmin,xmax,ymax,zmax
    write(*,'(a,1f8.1)') 'BBBBBB ',poisson%hgrid(1,1)*poisson%ngpx
    write(*,'(a,1f8.1)') 'BBBBBB ',poisson%hgrid(2,2)*poisson%ngpy
    write(*,'(a,1f8.1)') 'BBBBBB ',poisson%hgrid(3,3)*poisson%ngpz
    endif
    !-------------------------------------------------------
    !call fini_hartree(parini,atoms,poisson)
    call fini_hartree(parini,atoms,poisson_ion)
    !call atom_deallocate_old(atoms)
end subroutine get_trial_energy
!*****************************************************************************************
subroutine get_rmse(trial_energy,nbf,qq,rmse)
    implicit none
    type(typ_trial_energy), intent(inout):: trial_energy
    integer, intent(in):: nbf
    real(8), intent(in):: qq(nbf)
    real(8), intent(out):: rmse
    !local variables
    integer:: itrial, ibf
    real(8):: tt
    rmse=0.d0
    do itrial=1,trial_energy%ntrial
        tt=0.d0
        do ibf=1,nbf
            tt=tt+qq(ibf)*trial_energy%EP(ibf,itrial)
        enddo
        trial_energy%E_all(itrial)=trial_energy%EP_n(itrial)+tt
        rmse=rmse+(trial_energy%E_all(itrial)-trial_energy%energy(itrial))**2
    enddo
    rmse=sqrt(rmse/real(trial_energy%ntrial,kind=8))
end subroutine get_rmse
!*****************************************************************************************
end module mod_trial_energy
!*****************************************************************************************
