!*****************************************************************************************
module mod_atoms
    implicit none
    private
    public:: get_rat, get_rat_iat, swap_rat
    public:: set_rat, set_rat_iat, set_rat_atoms
    public:: update_rat, update_ratp
    public:: atom_allocate, atom_deallocate, atom_allocate_old, atom_deallocate_old
    public:: atom_all_allocate, atom_all_deallocate
    public:: atom_copy, atom_copy_old
    public:: atom_build_periodic_images
    public:: set_typat, atom_ddot
    public:: set_ndof
    public:: bemoved2string, string2bemoved
    public:: atoms_all_writexyz
    public:: atom_normalizevector
    public:: atom_calnorm, atom_calmaxforcecomponent
    public:: checkallatomstobeincell
    public:: determinexyzminmax
    public:: set_rcov, set_qat, set_atomic_mass, sat_to_iatom, iatom_to_sat
    type, public:: typ_atoms
        private
        integer, public:: nat=-1 !number of atoms
        integer, public:: natim=0 !number of atoms of all periodic images including itself 
        integer, public:: ndof=-1 !number of degrees of freedom
        real(8), public:: cellvec(3,3)=-1.d0 !cell vectors
        real(8), public:: celldv(3,3)=0.d0 !
        real(8), public:: stress(3,3)=0.d0 !
        real(8), public:: epot=0.d0 !potential energy
        real(8), public:: ekin=0.d0 !kinetic energy
        real(8), public:: etot=0.d0 !total energy
        real(8), public:: enth=0.d0 !enthalpy
        real(8), public:: ebattery=0.d0 !energy of external work of battery in p3d_bias
        real(8), public:: qtot=0.d0 !total charge
        real(8), public:: ztot=0.d0 !total ionic charge
        real(8), public:: dpm(3)=0.d0 !total electric dipole moment
        real(8), public:: elecfield(3)=0.d0 !external electric field
        real(8), public:: pressure=0.d0 !external pressure
        real(8), public:: tol
        real(8), public:: qtypat(20)=1.d20
        integer, public:: ntypat=-1
        integer, public:: ltypat(20)=-1
        integer, public:: nfp=-1
        character(5), public:: stypat(20)='none'
        character(20), public:: boundcond='unknown'
        character(10), public:: units='angstrom'
        character(10), public:: units_length_io='atomic'
        !coordinates type only at time of reading from file and writing to file.
        character(10), public:: coordinates_type='cartesian'
        character(50), public:: alloclist='all'
        character(5), allocatable, public:: sat(:) !symbol of atoms
        real(8), allocatable:: rat(:,:) !atomic positions
        real(8), allocatable, public:: ratp(:,:) !public atomic positions
        real(8), allocatable, public:: ratim(:,:) !atomic positions of periodic images
        real(8), allocatable, public:: vat(:,:) !atomic velocities
        real(8), allocatable, public:: amass(:) !atomic mass
        real(8), allocatable, public:: fat(:,:) !atomic forces
        logical, allocatable, public:: bemoved(:,:) !status to be moved or not
        real(8), allocatable, public:: qat(:) !atomic charges
        real(8), allocatable, public:: zat(:) !ionic charges
        real(8), allocatable, public:: rcov(:) !covalent radii
        real(8), allocatable, public:: fp(:) !fingerprint
        real(8), allocatable, public:: trial_ref_energy(:) 
        real(8), allocatable, public:: trial_ref_disp(:,:) 
        integer, allocatable, public:: trial_ref_nat(:) 
        integer, public:: ntrial=0
        integer, allocatable, public:: itypat(:) !The type of each atom is set in this array
        !contains
        !procedure:: atoms_assign
        !generic:: assignment(=) => atoms_assign
    end type typ_atoms
    type, public:: typ_atoms_all
        type(typ_atoms):: atoms
        integer:: nconfmax=-1 !maximum number of configurations
        integer:: nconf=-1 !number of configurations
        real(8), allocatable:: epotall(:) !potential energy
        real(8), allocatable:: qtotall(:) !total charge
        real(8), allocatable:: ratall(:,:,:) !atomic positions
        real(8), allocatable:: fatall(:,:,:) !atomic positions
        real(8), allocatable:: fpall(:,:) !fingerprint
    end type typ_atoms_all
    type, public:: typ_atoms_arr
        type(typ_atoms), allocatable:: atoms(:)
        integer:: nconfmax=-1 !maximum number of configurations
        integer:: nconf=-1 !number of configurations
        !integer, allocatable:: inclusion(:)
        character(50), allocatable:: fn(:)
        integer, allocatable:: lconf(:)
        logical, allocatable:: conf_inc(:)
        integer:: nconf_inc=-1
    end type typ_atoms_arr
    !type typ_atoms_arr
    !    integer:: natmax=10000
    !    integer, allocatable:: natarr(:)
    !    real(8), allocatable:: cellvecarr(:,:,:)
    !    real(8), allocatable:: epotarr(:)
    !    character(20), allocatable:: boundcondarr(:)
    !    character(5), allocatable:: sat(:,:)
    !    real(8), allocatable:: ratarr(:,:,:)
    !    logical, allocatable:: bemoved(:,:,:)
    !end type typ_atoms_arr
    type, public:: typ_file_info
        character(256):: filename_positions='unknown'
        character(256):: filename_forces='unknown'
        character(50):: file_position='unknown'
        character(50):: frmt='(a5,2x,3es24.15,2x,3l1)'
        logical:: print_force
        integer:: nconf=0
    end type typ_file_info
    !contains
    !subroutine atoms_assign(a,b)
    !    !import typ_atoms
    !    class(typ_atoms), intent(inout):: a
    !    class(typ_atoms), intent(in):: b
    !    write(*,'(a)') 'ERROR: assignment is not accept for type typ_atoms,'
    !    write(*,'(a)') '       please use subroutine atom_copy.'
    !    stop
    !end subroutine atoms_assign
    type, public:: type_pairs
        integer ,allocatable:: posat2nd(:,:)
    end type type_pairs
contains
!*****************************************************************************************
subroutine get_rat_iat(atoms,iat,xyz)
    implicit none
    type(typ_atoms), intent(in):: atoms
    integer, intent(in):: iat
    real(8), intent(out):: xyz(3)
    !local variables
    xyz(1)=atoms%rat(1,iat)
    xyz(2)=atoms%rat(2,iat)
    xyz(3)=atoms%rat(3,iat)
end subroutine get_rat_iat
!*****************************************************************************************
subroutine get_rat(atoms,rat)
    implicit none
    type(typ_atoms), intent(in):: atoms
    real(8), intent(out):: rat(3,atoms%nat)
    !local variables
    integer:: iat
    do iat=1,atoms%nat
        rat(1,iat)=atoms%rat(1,iat)
        rat(2,iat)=atoms%rat(2,iat)
        rat(3,iat)=atoms%rat(3,iat)
    enddo
end subroutine get_rat
!*****************************************************************************************
subroutine update_ratp(atoms)
    implicit none
    type(typ_atoms), intent(inout):: atoms
    !local variables
    integer:: iat
    do iat=1,atoms%nat
        atoms%ratp(1,iat)=atoms%rat(1,iat)
        atoms%ratp(2,iat)=atoms%rat(2,iat)
        atoms%ratp(3,iat)=atoms%rat(3,iat)
    enddo
end subroutine update_ratp
!*****************************************************************************************
subroutine update_rat(atoms,upall)
    implicit none
    type(typ_atoms), intent(inout):: atoms
    logical, intent(in), optional:: upall
    !local variables
    integer:: iat
    logical:: upall_in
    upall_in=.false.
    if(present(upall)) then
        upall_in=upall
    endif
    if(upall_in) then
        do iat=1,atoms%nat
            atoms%rat(1,iat)=atoms%ratp(1,iat)
            atoms%rat(2,iat)=atoms%ratp(2,iat)
            atoms%rat(3,iat)=atoms%ratp(3,iat)
        enddo
    else
        do iat=1,atoms%nat
            if(atoms%bemoved(1,iat)) atoms%rat(1,iat)=atoms%ratp(1,iat)
            if(atoms%bemoved(2,iat)) atoms%rat(2,iat)=atoms%ratp(2,iat)
            if(atoms%bemoved(3,iat)) atoms%rat(3,iat)=atoms%ratp(3,iat)
        enddo
    endif
end subroutine update_rat
!*****************************************************************************************
subroutine set_rat(atoms,rat,setall)
    implicit none
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(in):: rat(3,atoms%nat)
    logical, intent(in), optional:: setall
    !local variables
    integer:: iat
    logical:: setall_in
    setall_in=.false.
    if(present(setall)) then
        setall_in=setall
    endif
    if(setall_in) then
        do iat=1,atoms%nat
            atoms%rat(1,iat)=rat(1,iat)
            atoms%rat(2,iat)=rat(2,iat)
            atoms%rat(3,iat)=rat(3,iat)
        enddo
    else
        do iat=1,atoms%nat
            if(atoms%bemoved(1,iat)) atoms%rat(1,iat)=rat(1,iat)
            if(atoms%bemoved(2,iat)) atoms%rat(2,iat)=rat(2,iat)
            if(atoms%bemoved(3,iat)) atoms%rat(3,iat)=rat(3,iat)
        enddo
    endif
end subroutine set_rat
!*****************************************************************************************
subroutine set_rat_atoms(atoms,atoms_src,setall)
    implicit none
    type(typ_atoms), intent(inout):: atoms
    type(typ_atoms), intent(in):: atoms_src
    logical, intent(in), optional:: setall
    !local variables
    integer:: iat
    logical:: setall_in
    setall_in=.false.
    if(present(setall)) then
        setall_in=setall
    endif
    if(setall_in) then
        do iat=1,atoms%nat
            atoms%rat(1,iat)=atoms_src%rat(1,iat)
            atoms%rat(2,iat)=atoms_src%rat(2,iat)
            atoms%rat(3,iat)=atoms_src%rat(3,iat)
        enddo
    else
        do iat=1,atoms%nat
            if(atoms%bemoved(1,iat)) atoms%rat(1,iat)=atoms_src%rat(1,iat)
            if(atoms%bemoved(2,iat)) atoms%rat(2,iat)=atoms_src%rat(2,iat)
            if(atoms%bemoved(3,iat)) atoms%rat(3,iat)=atoms_src%rat(3,iat)
        enddo
    endif
end subroutine set_rat_atoms
!*****************************************************************************************
subroutine set_rat_iat(atoms,iat,xyz)
    implicit none
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in):: iat
    real(8), intent(in):: xyz(3)
    !local variables
    if(atoms%bemoved(1,iat)) atoms%rat(1,iat)=xyz(1)
    if(atoms%bemoved(2,iat)) atoms%rat(2,iat)=xyz(2)
    if(atoms%bemoved(3,iat)) atoms%rat(3,iat)=xyz(3)
end subroutine set_rat_iat
!*****************************************************************************************
subroutine swap_rat(atoms,iat,jat)
    implicit none
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in):: iat
    integer, intent(in):: jat
    !local variables
    real(8):: xyz(3)
    if(atoms%bemoved(1,iat)) stop 'ERROR: swap fixed atom'
    if(atoms%bemoved(2,iat)) stop 'ERROR: swap fixed atom'
    if(atoms%bemoved(3,iat)) stop 'ERROR: swap fixed atom'
    if(atoms%bemoved(1,jat)) stop 'ERROR: swap fixed atom'
    if(atoms%bemoved(2,jat)) stop 'ERROR: swap fixed atom'
    if(atoms%bemoved(3,jat)) stop 'ERROR: swap fixed atom'
    xyz(1)=atoms%rat(1,iat)
    xyz(2)=atoms%rat(2,iat)
    xyz(3)=atoms%rat(3,iat)
    atoms%rat(1,iat)=atoms%rat(1,jat)
    atoms%rat(2,iat)=atoms%rat(2,jat)
    atoms%rat(3,iat)=atoms%rat(3,jat)
    atoms%rat(1,jat)=xyz(1)
    atoms%rat(2,jat)=xyz(2)
    atoms%rat(3,jat)=xyz(3)
end subroutine swap_rat
!*****************************************************************************************
subroutine atom_allocate(atoms,nat,natim,nfp,ntrial)
    use dynamic_memory
    implicit none
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in):: nat, natim, nfp
    integer,optional,intent(in) :: ntrial
    !local variables
    integer:: iat !, ifp
    integer:: ind, ind_all
    !write(*,*) 'in atom_allocate: HERE'
    !if(atoms%natim>0 .and. trim(atoms%boundcond)=='free') then
    !    write(*,'(a)') 'WARNING: Do you really need atoms of periodic images with free BC'
    !endif
    if(present(ntrial)) then
    end if
    if(nat<1) stop 'ERROR: in atom_allocate: nat must be larger than zero'
    ind_all=index(atoms%alloclist,'all')
    !-----------------------------------------------------------------
    if(natim>0) then
        if(allocated(atoms%ratim)) stop 'ERROR: ratim is already allocated'
        !if(atoms%natim/=natim) stop 'ERROR: atoms%natim/=natim'
        atoms%ratim=f_malloc0([1.to.3,1.to.natim],id='atoms%ratim')
        atoms%alloclist=atoms%alloclist//':ratim'
    endif
    if(nfp>0) then
        if(allocated(atoms%fp)) stop 'ERROR: fp is already allocated'
        atoms%fp=f_malloc0([1.to.nfp],id='atoms%fp')
        atoms%alloclist=atoms%alloclist//':fp'
    endif
    !-----------------------------------------------------------------
    !allocation of components based on alloclist
    ind=index(atoms%alloclist,'sat')
    if(ind_all>0 .or. ind>0) then
        if(allocated(atoms%sat)) stop 'ERROR: sat is already allocated'
        atoms%sat=f_malloc_str(5,[1.to.nat],id='atoms%sat')
        do iat=1,nat
            atoms%sat(iat)='     '
        enddo
    endif
    ind=index(atoms%alloclist,'rat')
    if(ind_all>0 .or. ind>0) then
        if(allocated(atoms%rat)) stop 'ERROR: rat is already allocated'
        atoms%rat=f_malloc0([1.to.3,1.to.nat],id='atoms%rat')
    endif
    ind=index(atoms%alloclist,'ratp')
    if(ind_all>0 .or. ind>0) then
        if(allocated(atoms%ratp)) stop 'ERROR: ratp is already allocated'
        atoms%ratp=f_malloc0([1.to.3,1.to.nat],id='atoms%ratp')
    endif
    ind=index(atoms%alloclist,'vat')
    if(ind_all>0 .or. ind>0) then
        if(allocated(atoms%vat)) stop 'ERROR: vat is already allocated'
        atoms%vat=f_malloc0([1.to.3,1.to.nat],id='atoms%vat')
    endif
    ind=index(atoms%alloclist,'amass')
    if(ind_all>0 .or. ind>0) then
        if(allocated(atoms%amass)) stop 'ERROR: amass is already allocated'
        atoms%amass=f_malloc0([1.to.nat],id='atoms%amass')
    endif
    ind=index(atoms%alloclist,'fat')
    if(ind_all>0 .or. ind>0) then
        if(allocated(atoms%fat)) stop 'ERROR: fat is already allocated'
        atoms%fat=f_malloc0([1.to.3,1.to.nat],id='atoms%fat')
    endif
    ind=index(atoms%alloclist,'bemoved')
    if(ind_all>0 .or. ind>0) then
        if(allocated(atoms%bemoved)) stop 'ERROR: bemoved is already allocated'
        atoms%bemoved=f_malloc([1.to.3,1.to.nat],id='atoms%bemoved')
        atoms%bemoved=.true.
    endif
    ind=index(atoms%alloclist,'qat')
    if(ind_all>0 .or. ind>0) then
        if(allocated(atoms%qat)) stop 'ERROR: qat is already allocated'
        atoms%qat=f_malloc0([1.to.nat],id='atoms%qat')
    endif
    ind=index(atoms%alloclist,'zat')
    if(ind_all>0 .or. ind>0) then
        if(allocated(atoms%zat)) stop 'ERROR: zat is already allocated'
        atoms%zat=f_malloc0([1.to.nat],id='atoms%zat')
    endif
    ind=index(atoms%alloclist,'rcov')
    if(ind_all>0 .or. ind>0) then
        if(allocated(atoms%rcov)) stop 'ERROR: rcov is already allocated'
        atoms%rcov=f_malloc0([1.to.nat],id='atoms%rcov')
    endif
    ind=index(atoms%alloclist,'itypat')
    if(ind_all>0 .or. ind>0) then
        if(allocated(atoms%itypat)) stop 'ERROR: qat is already allocated'
        atoms%itypat=f_malloc0([1.to.nat],id='atoms%itypat')
    endif
    ind=index(atoms%alloclist,'trial_ref_energy')
    if((ind_all>0 .or. ind>0) .and. present(ntrial)) then
        if(allocated(atoms%trial_ref_energy)) stop 'ERROR: trial_ref_energy is already allocated'
        if(allocated(atoms%trial_ref_disp)) stop 'ERROR: trial_ref_disp is already allocated'
        if(allocated(atoms%trial_ref_nat)) stop 'ERROR: trial_ref_nat is already allocated'
        atoms%trial_ref_energy=f_malloc0([1.to.ntrial],id='atoms%trial_ref_energy')
        atoms%trial_ref_nat=f_malloc0([1.to.ntrial],id='atoms%trial_ref_nat')
        atoms%trial_ref_disp=f_malloc0([1.to.3,1.to.ntrial],id='atoms%trial_ref_disp')
    endif
    if(present(ntrial)) atoms%ntrial=ntrial
    atoms%nat=nat
    atoms%natim=natim
    atoms%nfp=nfp
end subroutine atom_allocate
!*****************************************************************************************
subroutine atom_deallocate(atoms)
    use dynamic_memory
    implicit none
    type(typ_atoms), intent(inout):: atoms
    !local variables
    integer:: ind_all, ind
    ind_all=index(atoms%alloclist,'all')
    !-----------------------------------------------------------------
    ind=index(atoms%alloclist,'ratim')
    if(ind>0) then
        if(.not. allocated(atoms%ratim)) then
            stop 'ERROR: ratim is not allocated'
        else
            call f_free(atoms%ratim)
        endif
    endif
    ind=index(atoms%alloclist,'fp')
    if(ind>0) then
        if(.not. allocated(atoms%fp)) then
            stop 'ERROR: fp is not allocated'
        else
            call f_free(atoms%fp)
        endif
    endif
    !-----------------------------------------------------------------
    ind=index(atoms%alloclist,'sat')
    if((ind_all>0 .and. allocated(atoms%sat)) .or. ind>0) then
        if(.not. allocated(atoms%sat)) then
            stop 'ERROR: sat is not allocated'
        else
            call f_free_str(5,atoms%sat)
        endif
    endif
    ind=index(atoms%alloclist,'rat')
    if((ind_all>0 .and. allocated(atoms%rat)) .or. ind>0) then
        if(.not. allocated(atoms%rat)) then
            stop 'ERROR: rat is not allocated'
        else
            call f_free(atoms%rat)
        endif
    endif
    ind=index(atoms%alloclist,'ratp')
    if((ind_all>0 .and. allocated(atoms%ratp)) .or. ind>0) then
        if(.not. allocated(atoms%ratp)) then
            stop 'ERROR: ratp is not allocated'
        else
            call f_free(atoms%ratp)
        endif
    endif
    ind=index(atoms%alloclist,'vat')
    if((ind_all>0 .and. allocated(atoms%vat)) .or. ind>0) then
        if(.not. allocated(atoms%vat)) then
            stop 'ERROR: vat is not allocated'
        else
            call f_free(atoms%vat)
        endif
    endif
    ind=index(atoms%alloclist,'amass')
    if((ind_all>0 .and. allocated(atoms%amass)) .or. ind>0) then
        if(.not. allocated(atoms%amass)) then
            stop 'ERROR: amass is not allocated'
        else
            call f_free(atoms%amass)
        endif
    endif
    ind=index(atoms%alloclist,'fat')
    if((ind_all>0 .and. allocated(atoms%fat)) .or. ind>0) then
        if(.not. allocated(atoms%fat)) then
            stop 'ERROR: fat is not allocated'
        else
            call f_free(atoms%fat)
        endif
    endif
    ind=index(atoms%alloclist,'bemoved')
    if((ind_all>0 .and. allocated(atoms%bemoved)) .or. ind>0) then
        if(.not. allocated(atoms%bemoved)) then
            stop 'ERROR: bemoved is not allocated'
        else
            call f_free(atoms%bemoved)
        endif
    endif
    ind=index(atoms%alloclist,'qat')
    if((ind_all>0 .and. allocated(atoms%qat)) .or. ind>0) then
        if(.not. allocated(atoms%qat)) then
            stop 'ERROR: qat is not allocated'
        else
            call f_free(atoms%qat)
        endif
    endif
    ind=index(atoms%alloclist,'zat')
    if((ind_all>0 .and. allocated(atoms%zat)) .or. ind>0) then
        if(.not. allocated(atoms%zat)) then
            stop 'ERROR: zat is not allocated'
        else
            call f_free(atoms%zat)
        endif
    endif
    ind=index(atoms%alloclist,'rcov')
    if((ind_all>0 .and. allocated(atoms%rcov)) .or. ind>0) then
        if(.not. allocated(atoms%rcov)) then
            stop 'ERROR: rcov is not allocated'
        else
            call f_free(atoms%rcov)
        endif
    endif
    ind=index(atoms%alloclist,'itypat')
    if((ind_all>0 .and. allocated(atoms%itypat)) .or. ind>0) then
        if(.not. allocated(atoms%itypat)) then
            stop 'ERROR: itypat is not allocated'
        else
            call f_free(atoms%itypat)
        endif
    endif
    ind=index(atoms%alloclist,'trial_ref_energy')
    if((ind_all>0 .and. allocated(atoms%trial_ref_energy)) .or. ind>0) then
        if(.not. allocated(atoms%trial_ref_energy)) then
            stop 'ERROR: trial_ref_energy is not allocated'
        else
            call f_free(atoms%trial_ref_energy)
        endif
    endif
    ind=index(atoms%alloclist,'trial_ref_disp')
    if((ind_all>0 .and. allocated(atoms%trial_ref_disp)) .or. ind>0) then
        if(.not. allocated(atoms%trial_ref_disp)) then
            stop 'ERROR: trial_ref_disp is not allocated'
        else
            call f_free(atoms%trial_ref_disp)
        endif
    endif
    ind=index(atoms%alloclist,'trial_ref_nat')
    if((ind_all>0 .and. allocated(atoms%trial_ref_nat)) .or. ind>0) then
        if(.not. allocated(atoms%trial_ref_nat)) then
            stop 'ERROR: trial_ref_disp is not allocated'
        else
            call f_free(atoms%trial_ref_nat)
        endif
    endif
end subroutine atom_deallocate
!*****************************************************************************************
subroutine atom_allocate_old(atoms,nat,natim,nfp,sat,vat,amass,fat,bemoved,qat,zat,rcov,typat&
                            ,ntrial,trial_ref_energy,trial_ref_nat,trial_ref_disp)
    implicit none
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in):: nat, natim, nfp
    integer, optional, intent(in):: ntrial
    logical, optional, intent(in):: sat, vat, amass
    logical, optional, intent(in):: fat, bemoved, qat, zat, rcov, typat
    logical, optional, intent(in):: trial_ref_energy,trial_ref_nat,trial_ref_disp
    !local variables
    logical:: all_of_them, l_arg(12)
    logical:: sat_t, vat_t, amass_t, fat_t, bemoved_t, qat_t, zat_t, rcov_t
    logical:: typat_t, trial_ref_energy_t, trial_ref_nat_t, trial_ref_disp_t
    sat_t=.false. ; vat_t=.false. ; amass_t=.false. ; fat_t=.false. ; bemoved_t=.false.
    qat_t=.false. ; zat_t=.false. ; rcov_t=.false.
    typat_t=.false. ; trial_ref_energy_t=.false. ; trial_ref_nat_t=.false.
    trial_ref_disp_t=.false.
    !integer:: iat, ifp
    !write(*,*) 'in atom_allocate_old: HERE'
    atoms%nat=nat
    atoms%natim=natim
    atoms%nfp=nfp
    if (present(ntrial)) then
        atoms%ntrial=ntrial
    end if
    !if(atoms%natim>0 .and. trim(atoms%boundcond)=='free') then
    !    write(*,'(a)') 'WARNING: Do you really need atoms of periodic images with free BC'
    !endif
    if(atoms%nat<1) stop 'ERROR: in atom_allocate_old: nat must be larger than zero'
    l_arg=(/present(sat),present(vat),present(amass),present(fat), &
        present(bemoved),present(qat),present(zat),present(rcov),present(typat), &
        present(trial_ref_energy),present(trial_ref_nat),present(trial_ref_disp)/)
    if(any(l_arg)) then
        all_of_them=.false.
    else
        all_of_them=.true.
    endif
    if(all_of_them) sat_t=.true.
    if(present(sat)) then;if(sat) sat_t=.true.;endif
    if(sat_t .and. .not. allocated(atoms%sat)) then
        allocate(atoms%sat(atoms%nat),source='     ')
    endif
    if(atoms%nat>0 .and. .not. allocated(atoms%rat)) then
        allocate(atoms%rat(3,atoms%nat),source=0.d0)
    endif
    if(atoms%nat>0 .and. .not. allocated(atoms%ratp)) then
        allocate(atoms%ratp(3,atoms%nat),source=0.d0)
    endif
    if(atoms%natim>0 .and. .not. allocated(atoms%ratim)) then
        allocate(atoms%ratim(3,atoms%natim),source=0.d0)
    endif
    if(all_of_them) vat_t=.true.
    if(present(vat)) then;if(vat) vat_t=.true.;endif
    if(vat_t .and. .not. allocated(atoms%vat)) then
        allocate(atoms%vat(3,atoms%nat),source=0.d0)
    endif
    if(all_of_them) amass_t=.true.
    if(present(amass)) then;if(amass) amass_t=.true.;endif
    if(amass_t .and. .not. allocated(atoms%amass)) then
        allocate(atoms%amass(atoms%nat),source=0.d0)
    endif
    if(all_of_them) fat_t=.true.
    if(present(fat)) then;if(fat) fat_t=.true.;endif
    if(fat_t .and. .not. allocated(atoms%fat)) then
        allocate(atoms%fat(3,atoms%nat),source=0.d0)
    endif
    if(all_of_them) bemoved_t=.true.
    if(present(bemoved)) then;if(bemoved) bemoved_t=.true.;endif
    if(bemoved_t .and. .not. allocated(atoms%bemoved)) then
        allocate(atoms%bemoved(3,atoms%nat),source=.true.)
    endif
    if(all_of_them) qat_t=.true.
    if(present(qat)) then;if(qat) qat_t=.true.;endif
    if(qat_t .and. .not. allocated(atoms%qat)) then
        allocate(atoms%qat(atoms%nat),source=0.d0)
    endif
    if(all_of_them) zat_t=.true.
    if(present(zat)) then;if(zat) zat_t=.true.;endif
    if(zat_t .and. .not. allocated(atoms%zat)) then
        allocate(atoms%zat(atoms%nat),source=0.d0)
    endif
    if(all_of_them) rcov_t=.true.
    if(present(rcov)) then;if(rcov) rcov_t=.true.;endif
    if(rcov_t .and. .not. allocated(atoms%rcov)) then
        allocate(atoms%rcov(atoms%nat),source=0.d0)
    endif
    if(nfp>0 .and. .not. allocated(atoms%fp)) then
        allocate(atoms%fp(atoms%nfp),source=0.d0)
    endif
    !if(atoms%nat>0 .and. .not. allocated(atoms%trial_ref_energy)) then
    !    allocate(atoms%trial_ref_energy(1:atoms%nat),source=0.d0)
    !endif
    if(all_of_them) typat_t=.true.
    if(present(typat)) then;if(typat) typat_t=.true.;endif
    if(typat_t .and. .not. allocated(atoms%itypat)) then
        allocate(atoms%itypat(atoms%nat),source=0)
    endif
    if(all_of_them) trial_ref_energy_t=.true.
    if(present(trial_ref_energy)) then;if(trial_ref_energy) trial_ref_energy_t=.true.;endif
    if(trial_ref_energy_t .and. .not. allocated(atoms%trial_ref_energy)) then
        allocate(atoms%trial_ref_energy(atoms%ntrial),source=0.d0)
    endif
    if(all_of_them) trial_ref_disp_t=.true.
    if(present(trial_ref_disp)) then;if(trial_ref_disp) trial_ref_disp_t=.true.;endif
    if(trial_ref_disp_t .and. .not. allocated(atoms%trial_ref_disp)) then
        allocate(atoms%trial_ref_disp(3,atoms%ntrial),source=0.d0)
    endif
    if(all_of_them) trial_ref_nat_t=.true.
    if(present(trial_ref_nat)) then;if(trial_ref_nat) trial_ref_nat_t=.true.;endif
    if(trial_ref_nat_t .and. .not. allocated(atoms%trial_ref_nat)) then
        allocate(atoms%trial_ref_nat(atoms%ntrial),source=0)
    endif
end subroutine atom_allocate_old
!*****************************************************************************************
subroutine atom_deallocate_old(atoms,sat,rat,ratim,vat,amass,fat,bemoved,qat,zat,rcov,fp,typat,&
                               trial_ref_energy,trial_ref_nat, trial_ref_disp)
    implicit none
    type(typ_atoms), intent(inout):: atoms
    logical, optional, intent(in):: sat, rat, ratim, vat, amass
    logical, optional, intent(in):: fat, bemoved, qat, zat, rcov, fp, typat
    logical, optional, intent(in):: trial_ref_energy, trial_ref_nat, trial_ref_disp
    !local variables
    logical:: all_of_them, l_arg(15)
    !integer::
    logical:: sat_t, vat_t, amass_t, fat_t, bemoved_t, qat_t, zat_t, rcov_t, fp_t, rat_t
    logical:: typat_t, trial_ref_energy_t, trial_ref_nat_t, trial_ref_disp_t, ratim_t
    sat_t=.false. ; vat_t=.false. ; amass_t=.false. ; fat_t=.false. ; bemoved_t=.false.
    qat_t=.false. ; zat_t=.false. ; rcov_t=.false. ; fp_t=.false. ; ratim_t=.false.
    typat_t=.false. ; trial_ref_energy_t=.false. ; trial_ref_nat_t=.false.
    trial_ref_disp_t=.false. ; rat_t=.false.
    l_arg=(/present(sat),present(rat),present(vat),present(amass),present(fat), &
        present(ratim),present(bemoved),present(qat),present(zat),present(rcov),&
        present(fp), present(typat), &
        present(trial_ref_energy),present(trial_ref_nat), present(trial_ref_disp)/)
    if(any(l_arg)) then
        all_of_them=.false.
    else
        all_of_them=.true.
    endif
    if(all_of_them) sat_t=.true.
    if(present(sat)) then;if(sat) sat_t=.true.;endif
    if(sat_t) then
        if(.not. allocated(atoms%sat)) then
            if(.not. all_of_them) then
                stop 'ERROR: sat is not allocated'
            endif
        else
            deallocate(atoms%sat)
        endif
    endif
    if(all_of_them) rat_t=.true.
    if(present(rat)) then;if(rat) rat_t=.true.;endif
    if(rat_t) then
        if(.not. allocated(atoms%rat)) then
            if(.not. all_of_them) then
                stop 'ERROR: rat is not allocated'
            endif
        else
            deallocate(atoms%rat)
        endif
    endif
    if(rat_t) then
        if(.not. allocated(atoms%ratp)) then
            if(.not. all_of_them) then
                stop 'ERROR: ratp is not allocated'
            endif
        else
            deallocate(atoms%ratp)
        endif
    endif
    if(all_of_them) ratim_t=.true.
    if(present(ratim)) then;if(ratim) ratim_t=.true.;endif
    if(ratim_t) then
        if(.not. allocated(atoms%ratim)) then
            if(.not. all_of_them) then
                stop 'ERROR: ratim is not allocated'
            endif
        else
            deallocate(atoms%ratim)
        endif
    endif
    if(all_of_them) vat_t=.true.
    if(present(vat)) then;if(vat) vat_t=.true.;endif
    if(vat_t) then
        if(.not. allocated(atoms%vat)) then
            if(.not. all_of_them) then
                stop 'ERROR: vat is not allocated'
            endif
        else
            deallocate(atoms%vat)
        endif
    endif
    if(all_of_them) amass_t=.true.
    if(present(amass)) then;if(amass) amass_t=.true.;endif
    if(amass_t) then
        if(.not. allocated(atoms%amass)) then
            if(.not. all_of_them) then
                stop 'ERROR: amass is not allocated'
            endif
        else
            deallocate(atoms%amass)
        endif
    endif
    if(all_of_them) fat_t=.true.
    if(present(fat)) then;if(fat) fat_t=.true.;endif
    if(fat_t) then
        if(.not. allocated(atoms%fat)) then
            if(.not. all_of_them) then
                stop 'ERROR: fat is not allocated'
            endif
        else
            deallocate(atoms%fat)
        endif
    endif
    if(all_of_them) bemoved_t=.true.
    if(present(bemoved)) then;if(bemoved) bemoved_t=.true.;endif
    if(bemoved_t) then
        if(.not. allocated(atoms%bemoved)) then
            if(.not. all_of_them) then
                stop 'ERROR: bemoved is not allocated'
            endif
        else
            deallocate(atoms%bemoved)
        endif
    endif
    if(all_of_them) qat_t=.true.
    if(present(qat)) then;if(qat) qat_t=.true.;endif
    if(qat_t) then
        if(.not. allocated(atoms%qat)) then
            if(.not. all_of_them) then
                stop 'ERROR: qat is not allocated'
            endif
        else
            deallocate(atoms%qat)
        endif
    endif
    if(all_of_them) zat_t=.true.
    if(present(zat)) then;if(zat) zat_t=.true.;endif
    if(zat_t) then
        if(.not. allocated(atoms%zat)) then
            if(.not. all_of_them) then
                stop 'ERROR: zat is not allocated'
            endif
        else
            deallocate(atoms%zat)
        endif
    endif
    if(all_of_them) rcov_t=.true.
    if(present(rcov)) then;if(rcov) rcov_t=.true.;endif
    if(rcov_t) then
        if(.not. allocated(atoms%rcov)) then
            if(.not. all_of_them) then
                stop 'ERROR: rcov is not allocated'
            endif
        else
            deallocate(atoms%rcov)
        endif
    endif
    if(all_of_them) fp_t=.true.
    if(present(fp)) then;if(fp) fp_t=.true.;endif
    if(fp_t) then
        if(.not. allocated(atoms%fp)) then
            if(.not. all_of_them) then
                stop 'ERROR: fp is not allocated'
            endif
        else
            deallocate(atoms%fp)
        endif
    endif
    if(all_of_them) typat_t=.true.
    if(present(typat)) then;if(typat) typat_t=.true.;endif
    if(typat_t) then
        if(.not. allocated(atoms%itypat)) then
            if(.not. all_of_them) then
                stop 'ERROR: itypat is not allocated'
            endif
        else
            deallocate(atoms%itypat)
        endif
    endif
    if(all_of_them) trial_ref_energy_t=.true.
    if(present(trial_ref_energy)) then;if(trial_ref_energy) trial_ref_energy_t=.true.;endif
    if(trial_ref_energy_t) then
        if(.not. allocated(atoms%trial_ref_energy)) then
            if(.not. all_of_them) then
                stop 'ERROR: trial_ref_energy is not allocated'
            endif
        else
            deallocate(atoms%trial_ref_energy)
        endif
    endif
    if(all_of_them) trial_ref_nat_t=.true.
    if(present(trial_ref_nat)) then;if(trial_ref_nat) trial_ref_nat_t=.true.;endif
    if(trial_ref_nat_t) then
        if(.not. allocated(atoms%trial_ref_nat)) then
            if(.not. all_of_them) then
                stop 'ERROR: trial_ref_nat is not allocated'
            endif
        else
            deallocate(atoms%trial_ref_nat)
        endif
    endif
    if(all_of_them) trial_ref_disp_t=.true.
    if(present(trial_ref_disp)) then;if(trial_ref_disp) trial_ref_disp_t=.true.;endif
    if(trial_ref_disp_t) then
        if(.not. allocated(atoms%trial_ref_disp)) then
            if(.not. all_of_them) then
                stop 'ERROR: trial_ref_disp is not allocated'
            endif
        else
            deallocate(atoms%trial_ref_disp)
        endif
    endif
end subroutine atom_deallocate_old
!*****************************************************************************************
subroutine atom_all_allocate(atoms_all,ratall,fatall,epotall,fpall,qtotall)
    implicit none
    type(typ_atoms_all), intent(inout):: atoms_all
    logical, optional, intent(in):: ratall, fatall, epotall, fpall, qtotall
    !local variables
    !integer:: iat, ifp
    if(atoms_all%atoms%nat<1) stop 'ERROR: atoms%nat must be larger than zero'
    if(atoms_all%nconfmax<1) stop 'ERROR: atoms_all%nconfmax must be larger than zero'
    if(present(ratall) .and. ratall) then
        if(allocated(atoms_all%ratall)) stop 'ERROR: ratall is already allocated'
        allocate(atoms_all%ratall(3,atoms_all%atoms%nat,atoms_all%nconfmax),source=0.d0)
    endif
    if(present(fatall) .and. fatall) then
        if(allocated(atoms_all%fatall)) stop 'ERROR: fatall is already allocated'
        allocate(atoms_all%fatall(3,atoms_all%atoms%nat,atoms_all%nconfmax),source=0.d0)
    endif
    if(present(epotall) .and. epotall) then
        if(allocated(atoms_all%epotall)) stop 'ERROR: epotall is already allocated'
        allocate(atoms_all%epotall(atoms_all%nconfmax),source=0.d0)
    endif
    if(present(fpall) .and. fpall) then
        if(allocated(atoms_all%fpall)) stop 'ERROR: fpall is already allocated'
        if(atoms_all%atoms%nfp<0) stop 'ERROR: nfp is not set.'
        allocate(atoms_all%fpall(atoms_all%atoms%nfp,atoms_all%nconfmax),source=0.d0)
    endif
    if(present(qtotall) .and. qtotall) then
        if(allocated(atoms_all%qtotall)) stop 'ERROR: qtotall is already allocated'
        allocate(atoms_all%qtotall(atoms_all%nconfmax),source=0.d0)
    endif
end subroutine atom_all_allocate
!*****************************************************************************************
subroutine atom_all_deallocate(atoms_all,ratall,fatall,epotall,fpall,qtotall)
    implicit none
    type(typ_atoms_all), intent(inout):: atoms_all
    logical, optional, intent(in):: ratall, fatall, epotall, fpall, qtotall
    !local variables
    if(present(ratall) .and. ratall) then
        if(.not. allocated(atoms_all%ratall)) stop 'ERROR: ratall is not allocated'
        deallocate(atoms_all%ratall)
    endif
    if(present(fatall) .and. fatall) then
        if(.not. allocated(atoms_all%fatall)) stop 'ERROR: fatall is not allocated'
        deallocate(atoms_all%fatall)
    endif
    if(present(epotall) .and. epotall) then
        if(.not. allocated(atoms_all%epotall)) stop 'ERROR: epotall is not allocated'
        deallocate(atoms_all%epotall)
    endif
    if(present(fpall) .and. fpall) then
        if(.not. allocated(atoms_all%fpall)) stop 'ERROR: fpall is not allocated'
        deallocate(atoms_all%fpall)
    endif
    if(present(qtotall) .and. qtotall) then
        if(.not. allocated(atoms_all%qtotall)) stop 'ERROR: qtotall is not allocated'
        deallocate(atoms_all%qtotall)
    endif
end subroutine atom_all_deallocate
!*****************************************************************************************
subroutine atom_copy(at_inp,at_out,str_message)
    use dynamic_memory
    implicit none
    type(typ_atoms), intent(in):: at_inp
    type(typ_atoms), intent(inout):: at_out
    character(*):: str_message
    !local variables
    integer:: iat, ishape(2), ifp, ind_all, ind
    character(100):: err_mess
    if(at_inp%nat<1) stop 'ERROR: at_inp%nat must be larger than zero'
    ind_all=index(at_inp%alloclist,'all')
    at_out%nat=at_inp%nat
    !if(ind_all>0) then
        at_out%ndof=at_inp%ndof
        at_out%boundcond=at_inp%boundcond
        at_out%units=at_inp%units
        at_out%coordinates_type=at_inp%coordinates_type
        at_out%cellvec(1,1)=at_inp%cellvec(1,1)
        at_out%cellvec(2,1)=at_inp%cellvec(2,1)
        at_out%cellvec(3,1)=at_inp%cellvec(3,1)
        at_out%cellvec(1,2)=at_inp%cellvec(1,2)
        at_out%cellvec(2,2)=at_inp%cellvec(2,2)
        at_out%cellvec(3,2)=at_inp%cellvec(3,2)
        at_out%cellvec(1,3)=at_inp%cellvec(1,3)
        at_out%cellvec(2,3)=at_inp%cellvec(2,3)
        at_out%cellvec(3,3)=at_inp%cellvec(3,3)
        at_out%epot=at_inp%epot
        at_out%ekin=at_inp%ekin
        at_out%etot=at_inp%etot
        at_out%nfp=at_inp%nfp
        at_out%tol=at_inp%tol
        at_out%qtot=at_inp%qtot
        at_out%units_length_io=at_inp%units_length_io
        at_out%alloclist=at_inp%alloclist
        at_out%dpm(1)=at_inp%dpm(1)
        at_out%dpm(2)=at_inp%dpm(2)
        at_out%dpm(3)=at_inp%dpm(3)
        at_out%elecfield(1)=at_inp%elecfield(1)
        at_out%elecfield(2)=at_inp%elecfield(2)
        at_out%elecfield(3)=at_inp%elecfield(3)

    !endif
    !copying array at_inp%sat to at_out%sat
    ind=index(at_inp%alloclist,'sat')
    if(ind_all>0 .or. ind>0) then
        if(allocated(at_inp%sat)) then
            if(allocated(at_out%sat)) then
                ishape(1:1)=shape(at_out%sat)
                if(at_inp%nat/=ishape(1)) then
                    call f_free_str(5,at_out%sat)
                endif
            endif
            if(.not. allocated(at_out%sat)) then
                at_out%sat=f_malloc_str(5,[1.to.at_out%nat],id='at_out%sat')
            endif
            do iat=1,at_inp%nat
                at_out%sat(iat)=at_inp%sat(iat)
            enddo
        else
            err_mess='ERROR: sat in at_inp%alloclist but at_inp%sat not allocated:'
            write(*,'(a,1x,a)') trim(err_mess),trim(str_message)
            stop
        endif
    endif
    !copying array at_inp%rat to at_out%rat
    ind=index(at_inp%alloclist,'rat')
    if(ind_all>0 .or. ind>0) then
        if(allocated(at_inp%rat)) then
            if(allocated(at_out%rat)) then
                ishape(1:2)=shape(at_out%rat)
                if(at_inp%nat/=ishape(2)) then
                    call f_free(at_out%rat)
                endif
            endif
            if(.not. allocated(at_out%rat)) then
                at_out%rat=f_malloc([1.to.3,1.to.at_out%nat],id='at_out%rat')
            endif
            do iat=1,at_inp%nat
                at_out%rat(1,iat)=at_inp%rat(1,iat)
                at_out%rat(2,iat)=at_inp%rat(2,iat)
                at_out%rat(3,iat)=at_inp%rat(3,iat)
            enddo
        else
            err_mess='ERROR: rat in at_inp%alloclist but at_inp%rat not allocated:'
            write(*,'(a,1x,a)') trim(err_mess),trim(str_message)
            stop
        endif
        if(allocated(at_inp%ratp)) then
            if(allocated(at_out%ratp)) then
                ishape(1:2)=shape(at_out%ratp)
                if(at_inp%nat/=ishape(2)) then
                    call f_free(at_out%ratp)
                endif
            endif
            if(.not. allocated(at_out%ratp)) then
                at_out%ratp=f_malloc([1.to.3,1.to.at_out%nat],id='at_out%ratp')
            endif
            do iat=1,at_inp%nat
                at_out%ratp(1,iat)=at_inp%ratp(1,iat)
                at_out%ratp(2,iat)=at_inp%ratp(2,iat)
                at_out%ratp(3,iat)=at_inp%ratp(3,iat)
            enddo
        else
            err_mess='ERROR: ratp in at_inp%alloclist but at_inp%ratp not allocated:'
            write(*,'(a,1x,a)') trim(err_mess),trim(str_message)
            stop
        endif
    endif
    !copying array at_inp%ratim to at_out%ratim
    ind=index(at_inp%alloclist,'ratim')
    if(ind>0) then
        if(allocated(at_inp%ratim)) then
            if(allocated(at_out%ratim)) then
                ishape(1:2)=shape(at_out%ratim)
                if(at_inp%natim/=ishape(2)) then
                    call f_free(at_out%ratim)
                endif
            endif
            if(.not. allocated(at_out%ratim)) then
                at_out%ratim=f_malloc([1.to.3,1.to.at_out%natim],id='at_out%ratim')
            endif
            do iat=1,at_inp%natim
                at_out%ratim(1,iat)=at_inp%ratim(1,iat)
                at_out%ratim(2,iat)=at_inp%ratim(2,iat)
                at_out%ratim(3,iat)=at_inp%ratim(3,iat)
            enddo
        else
            err_mess='ERROR: ratim in at_inp%alloclist but at_inp%ratim not allocated:'
            write(*,'(a,1x,a)') trim(err_mess),trim(str_message)
            stop
        endif
    endif
    !copying array at_inp%vat to at_out%vat
    ind=index(at_inp%alloclist,'vat')
    if(ind_all>0 .or. ind>0) then
        if(allocated(at_inp%vat)) then
            if(allocated(at_out%vat)) then
                ishape(1:2)=shape(at_out%vat)
                if(at_inp%nat/=ishape(2)) then
                    call f_free(at_out%vat)
                endif
            endif
            if(.not. allocated(at_out%vat)) then
                at_out%vat=f_malloc([1.to.3,1.to.at_out%nat],id='at_out%vat')
            endif
            do iat=1,at_inp%nat
                at_out%vat(1,iat)=at_inp%vat(1,iat)
                at_out%vat(2,iat)=at_inp%vat(2,iat)
                at_out%vat(3,iat)=at_inp%vat(3,iat)
            enddo
        else
            err_mess='ERROR: vat in at_inp%alloclist but at_inp%vat not allocated:'
            write(*,'(a,1x,a)') trim(err_mess),trim(str_message)
            stop
        endif
    endif
    !copying array at_inp%amass to at_out%amass
    ind=index(at_inp%alloclist,'amass')
    if(ind_all>0 .or. ind>0) then
        if(allocated(at_inp%amass)) then
            if(allocated(at_out%amass)) then
                ishape(1:1)=shape(at_out%amass)
                if(at_inp%nat/=ishape(1)) then
                    call f_free(at_out%amass)
                endif
            endif
            if(.not. allocated(at_out%amass)) then
                at_out%amass=f_malloc([1.to.at_out%nat],id='at_out%amass')
            endif
            do iat=1,at_inp%nat
                at_out%amass(iat)=at_inp%amass(iat)
            enddo
        else
            err_mess='ERROR: amass in at_inp%alloclist but at_inp%amass not allocated:'
            write(*,'(a,1x,a)') trim(err_mess),trim(str_message)
            stop
        endif
    endif
    !copying array at_inp%fat to at_out%fat
    ind=index(at_inp%alloclist,'fat')
    if(ind_all>0 .or. ind>0) then
        if(allocated(at_inp%fat)) then
            if(allocated(at_out%fat)) then
                ishape(1:2)=shape(at_out%fat)
                if(at_inp%nat/=ishape(2)) then
                    call f_free(at_out%fat)
                endif
            endif
            if(.not. allocated(at_out%fat)) then
                at_out%fat=f_malloc([1.to.3,1.to.at_out%nat],id='at_out%fat')
            endif
            do iat=1,at_inp%nat
                at_out%fat(1,iat)=at_inp%fat(1,iat)
                at_out%fat(2,iat)=at_inp%fat(2,iat)
                at_out%fat(3,iat)=at_inp%fat(3,iat)
            enddo
        else
            err_mess='ERROR: fat in at_inp%alloclist but at_inp%fat not allocated:'
            write(*,'(a,1x,a)') trim(err_mess),trim(str_message)
            stop
        endif
    endif
    !copying array at_inp%bemoved to at_out%bemoved
    ind=index(at_inp%alloclist,'bemoved')
    if(ind_all>0 .or. ind>0) then
        if(allocated(at_inp%bemoved)) then
            if(allocated(at_out%bemoved)) then
                ishape(1:2)=shape(at_out%bemoved)
                if(at_inp%nat/=ishape(2)) then
                    call f_free(at_out%bemoved)
                endif
            endif
            if(.not. allocated(at_out%bemoved)) then
                at_out%bemoved=f_malloc([1.to.3,1.to.at_out%nat],id='at_out%bemoved')
            endif
            do iat=1,at_inp%nat
                at_out%bemoved(1,iat)=at_inp%bemoved(1,iat)
                at_out%bemoved(2,iat)=at_inp%bemoved(2,iat)
                at_out%bemoved(3,iat)=at_inp%bemoved(3,iat)
            enddo
        else
            err_mess='ERROR: bemoved in at_inp%alloclist but at_inp%bemoved not allocated:'
            write(*,'(a,1x,a)') trim(err_mess),trim(str_message)
            stop
        endif
    endif
    !copying array at_inp%qat to at_out%qat
    ind=index(at_inp%alloclist,'qat')
    if(ind_all>0 .or. ind>0) then
        if(allocated(at_inp%qat)) then
            if(allocated(at_out%qat)) then
                ishape(1:1)=shape(at_out%qat)
                if(at_inp%nat/=ishape(1)) then
                    call f_free(at_out%qat)
                endif
            endif
            if(.not. allocated(at_out%qat)) then
                at_out%qat=f_malloc([1.to.at_out%nat],id='at_out%qat')
            endif
            do iat=1,at_inp%nat
                at_out%qat(iat)=at_inp%qat(iat)
            enddo
        else
            err_mess='ERROR: qat in at_inp%alloclist but at_inp%qat not allocated:'
            write(*,'(a,1x,a)') trim(err_mess),trim(str_message)
            stop
        endif
    endif
    !copying array at_inp%zat to at_out%zat
    ind=index(at_inp%alloclist,'zat')
    if(ind_all>0 .or. ind>0) then
        if(allocated(at_inp%zat)) then
            if(allocated(at_out%zat)) then
                ishape(1:1)=shape(at_out%zat)
                if(at_inp%nat/=ishape(1)) then
                    call f_free(at_out%zat)
                endif
            endif
            if(.not. allocated(at_out%zat)) then
                at_out%zat=f_malloc([1.to.at_out%nat],id='at_out%zat')
            endif
            do iat=1,at_inp%nat
                at_out%zat(iat)=at_inp%zat(iat)
            enddo
        else
            err_mess='ERROR: zat in at_inp%alloclist but at_inp%zat not allocated:'
            write(*,'(a,1x,a)') trim(err_mess),trim(str_message)
            stop
        endif
    endif
    !copying array at_inp%rcov to at_out%rcov
    ind=index(at_inp%alloclist,'rcov')
    if(ind_all>0 .or. ind>0) then
        if(allocated(at_inp%rcov)) then
            if(allocated(at_out%rcov)) then
                ishape(1:1)=shape(at_out%rcov)
                if(at_inp%nat/=ishape(1)) then
                    call f_free(at_out%rcov)
                endif
            endif
            if(.not. allocated(at_out%rcov)) then
                at_out%rcov=f_malloc([1.to.at_out%nat],id='at_out%rcov')
            endif
            do iat=1,at_inp%nat
                at_out%rcov(iat)=at_inp%rcov(iat)
            enddo
        else
            err_mess='ERROR: rcov in at_inp%alloclist but at_inp%rcov not allocated:'
            write(*,'(a,1x,a)') trim(err_mess),trim(str_message)
            stop
        endif
    endif
    !copying array at_inp%fp to at_out%fp
    ind=index(at_inp%alloclist,'fp')
    if(ind>0) then
        if(allocated(at_inp%fp)) then
            if(allocated(at_out%fp)) then
                ishape(1:1)=shape(at_out%fp)
                if(at_inp%nfp/=ishape(1)) then
                    call f_free(at_out%fp)
                endif
            endif
            if(.not. allocated(at_out%fp)) then
                at_out%fp=f_malloc([1.to.at_out%nfp],id='at_out%fp')
            endif
            do iat=1,at_inp%nfp
                at_out%fp(iat)=at_inp%fp(iat)
            enddo
        else
            err_mess='ERROR: fp in at_inp%alloclist but at_inp%fp not allocated:'
            write(*,'(a,1x,a)') trim(err_mess),trim(str_message)
            stop
        endif
    endif
    !copying array at_inp%itypat to at_out%itypat
    ind=index(at_inp%alloclist,'itypat')
    if(ind_all>0 .or. ind>0) then
        if(allocated(at_inp%itypat)) then
            if(allocated(at_out%itypat)) then
                ishape(1:1)=shape(at_out%itypat)
                if(at_inp%nat/=ishape(1)) then
                    call f_free(at_out%itypat)
                endif
            endif
            if(.not. allocated(at_out%itypat)) then
                at_out%itypat=f_malloc([1.to.at_out%nat],id='at_out%itypat')
            endif
            do iat=1,at_inp%nat
                at_out%itypat(iat)=at_inp%itypat(iat)
            enddo
        else
            err_mess='ERROR: itypat in at_inp%alloclist but at_inp%itypat not allocated:'
            write(*,'(a,1x,a)') trim(err_mess),trim(str_message)
            stop
        endif
    endif
end subroutine atom_copy
!*****************************************************************************************
subroutine atom_copy_old(at_inp,at_out,str_message,sat,rat,ratim,vat,amass,fat,bemoved,qat,zat,rcov,fp,typat, &
                         trial_ref_energy,trial_ref_nat,trial_ref_disp)
    implicit none
    type(typ_atoms), intent(in):: at_inp
    type(typ_atoms), intent(inout):: at_out
    character(*):: str_message
    logical, optional, intent(in):: sat, rat, ratim, vat, amass
    logical, optional, intent(in):: fat, bemoved, qat, zat, rcov, fp, typat
    logical, optional, intent(in):: trial_ref_energy, trial_ref_nat, trial_ref_disp
    !local variables
    integer:: iat, ishape(2), ifp
    logical:: prsnt
    logical:: all_of_them, l_arg(15)
    character(100):: err_mess
    !write(*,*) 'in atom_copy_old: HERE'
    if(at_inp%nat<1) stop 'ERROR: atoms%nat must be larger than zero'
    !err_mess='ERROR: cannot copy typ_atoms contents with different nat:'
    !if(at_inp%nat/=at_out%nat) then
    !    write(*,'(a,2i6,a)') trim(err_mess),at_inp%nat,at_out%nat,trim(str_message)
    !    stop
    !endif
    l_arg=(/present(sat),present(rat),present(vat),present(amass),present(fat), &
        present(ratim),present(bemoved),present(qat),present(zat),present(rcov),&
        present(fp), present(typat), &
        present(trial_ref_energy), present(trial_ref_nat), present(trial_ref_disp)/)
    if(any(l_arg)) then
        all_of_them=.false.
    else
        all_of_them=.true.
    endif
    !if(at_inp%nat/=at_out%nat) all_of_them=.true.
    at_out%nat=at_inp%nat
    !if(all_of_them) then
        at_out%ndof=at_inp%ndof
        at_out%boundcond=at_inp%boundcond
        at_out%units=at_inp%units
        at_out%coordinates_type=at_inp%coordinates_type
        at_out%cellvec(1,1)=at_inp%cellvec(1,1)
        at_out%cellvec(2,1)=at_inp%cellvec(2,1)
        at_out%cellvec(3,1)=at_inp%cellvec(3,1)
        at_out%cellvec(1,2)=at_inp%cellvec(1,2)
        at_out%cellvec(2,2)=at_inp%cellvec(2,2)
        at_out%cellvec(3,2)=at_inp%cellvec(3,2)
        at_out%cellvec(1,3)=at_inp%cellvec(1,3)
        at_out%cellvec(2,3)=at_inp%cellvec(2,3)
        at_out%cellvec(3,3)=at_inp%cellvec(3,3)
        at_out%epot=at_inp%epot
        at_out%ekin=at_inp%ekin
        at_out%etot=at_inp%etot
        at_out%nfp=at_inp%nfp
        at_out%tol=at_inp%tol
        at_out%qtot=at_inp%qtot
        at_out%dpm(1)=at_inp%dpm(1)
        at_out%dpm(2)=at_inp%dpm(2)
        at_out%dpm(3)=at_inp%dpm(3)
        at_out%elecfield(1)=at_inp%elecfield(1)
        at_out%elecfield(2)=at_inp%elecfield(2)
        at_out%elecfield(3)=at_inp%elecfield(3)
    !endif
    !copying array at_inp%sat to at_out%sat
    if(present(sat)) then ; prsnt=sat ; else ; prsnt=.false. ;  endif
    if(allocated(at_inp%sat)) then
        if(all_of_them .or. prsnt) then
            if(allocated(at_out%sat)) then
                ishape(1:1)=shape(at_out%sat)
                if(at_inp%nat/=ishape(1)) then
                    call atom_deallocate_old(at_out,sat=.true.)
                endif
            endif
            if(.not. allocated(at_out%sat)) then
                call atom_allocate_old(at_out,at_inp%nat,at_inp%natim,at_inp%nfp,sat=.true.)
            endif
            do iat=1,at_inp%nat
                at_out%sat(iat)=at_inp%sat(iat)
            enddo
        endif
    else
        err_mess='ERROR: cannot copy typ_atoms%sat when source is not allocated:'
        if(prsnt) then
            write(*,'(a,1x,a)') trim(err_mess),trim(str_message)
            stop
        endif
    endif
    !copying array at_inp%rat to at_out%rat
    if(present(rat)) then ; prsnt=rat ; else ; prsnt=.false. ;  endif
    if(allocated(at_inp%rat)) then
        if(all_of_them .or. prsnt) then
            ishape(1:2)=shape(at_out%rat)
            if(allocated(at_out%rat) .and. at_inp%nat/=ishape(2)) then
                call atom_deallocate_old(at_out,rat=.true.)
            endif
            if(.not. allocated(at_out%rat)) then
                call atom_allocate_old(at_out,at_inp%nat,at_inp%natim,at_inp%nfp)
            endif
            do iat=1,at_inp%nat
                at_out%rat(1,iat)=at_inp%rat(1,iat)
                at_out%rat(2,iat)=at_inp%rat(2,iat)
                at_out%rat(3,iat)=at_inp%rat(3,iat)
            enddo
        endif
    else
        err_mess='ERROR: cannot copy typ_atoms%rat when source is not allocated:'
        if(prsnt) then
            write(*,'(a,1x,a)') trim(err_mess),trim(str_message)
            stop
        endif
    endif
    if(allocated(at_inp%ratp)) then
        if(all_of_them .or. prsnt) then
            ishape(1:2)=shape(at_out%ratp)
            if(allocated(at_out%ratp) .and. at_inp%nat/=ishape(2)) then
                call atom_deallocate_old(at_out,rat=.true.)
            endif
            if(.not. allocated(at_out%ratp)) then
                call atom_allocate_old(at_out,at_inp%nat,at_inp%natim,at_inp%nfp)
            endif
            do iat=1,at_inp%nat
                at_out%ratp(1,iat)=at_inp%ratp(1,iat)
                at_out%ratp(2,iat)=at_inp%ratp(2,iat)
                at_out%ratp(3,iat)=at_inp%ratp(3,iat)
            enddo
        endif
    else
        err_mess='ERROR: cannot copy typ_atoms%ratp when source is not allocated:'
        if(prsnt) then
            write(*,'(a,1x,a)') trim(err_mess),trim(str_message)
            stop
        endif
    endif
    !copying array at_inp%ratim to at_out%ratim
    if(present(ratim)) then ; prsnt=ratim ; else ; prsnt=.false. ;  endif
    if(allocated(at_inp%ratim)) then
        if(all_of_them .or. prsnt) then
            if(allocated(at_out%ratim)) then
                ishape(1:2)=shape(at_out%ratim)
                if(at_inp%natim/=ishape(2)) then
                    call atom_deallocate_old(at_out,ratim=.true.)
                endif
            endif
            if(.not. allocated(at_out%ratim)) then
                call atom_allocate_old(at_out,at_inp%nat,at_inp%natim,at_inp%nfp)
            endif
            do iat=1,at_inp%natim
                at_out%ratim(1,iat)=at_inp%ratim(1,iat)
                at_out%ratim(2,iat)=at_inp%ratim(2,iat)
                at_out%ratim(3,iat)=at_inp%ratim(3,iat)
            enddo
        endif
    else
        err_mess='ERROR: cannot copy typ_atoms%ratim when source is not allocated:'
        if(prsnt) then
            write(*,'(a,1x,a)') trim(err_mess),trim(str_message)
            stop
        endif
    endif
    !copying array at_inp%vat to at_out%vat
    if(present(vat)) then ; prsnt=vat ; else ; prsnt=.false. ;  endif
    if(allocated(at_inp%vat)) then
        if(all_of_them .or. prsnt) then
            if(allocated(at_out%vat)) then
                ishape(1:2)=shape(at_out%vat)
                if(at_inp%nat/=ishape(2)) then
                    call atom_deallocate_old(at_out,vat=.true.)
                endif
            endif
            if(.not. allocated(at_out%vat)) then
                call atom_allocate_old(at_out,at_inp%nat,at_inp%natim,at_inp%nfp,vat=.true.)
            endif
            do iat=1,at_inp%nat
                at_out%vat(1,iat)=at_inp%vat(1,iat)
                at_out%vat(2,iat)=at_inp%vat(2,iat)
                at_out%vat(3,iat)=at_inp%vat(3,iat)
            enddo
        endif
    else
        err_mess='ERROR: cannot copy typ_atoms%vat when source is not allocated:'
        if(prsnt) then
            write(*,'(a,1x,a)') trim(err_mess),trim(str_message)
            stop
        endif
    endif
    !copying array at_inp%amass to at_out%amass
    if(present(amass)) then ; prsnt=amass ; else ; prsnt=.false. ;  endif
    if(allocated(at_inp%amass)) then
        if(all_of_them .or. prsnt) then
            if(allocated(at_out%amass)) then
                ishape(1:1)=shape(at_out%amass)
                if(at_inp%nat/=ishape(1)) then
                    call atom_deallocate_old(at_out,amass=.true.)
                endif
            endif
            if(.not. allocated(at_out%amass)) then
                call atom_allocate_old(at_out,at_inp%nat,at_inp%natim,at_inp%nfp,amass=.true.)
            endif
            do iat=1,at_inp%nat
                at_out%amass(iat)=at_inp%amass(iat)
            enddo
        endif
    else
        err_mess='ERROR: cannot copy typ_atoms%amass when source is not allocated:'
        if(prsnt) then
            write(*,'(a,1x,a)') trim(err_mess),trim(str_message)
            stop
        endif
    endif
    !copying array at_inp%fat to at_out%fat
    if(present(fat)) then ; prsnt=fat ; else ; prsnt=.false. ;  endif
    if(allocated(at_inp%fat)) then
        if(all_of_them .or. prsnt) then
            if(allocated(at_out%fat)) then
                ishape(1:2)=shape(at_out%fat)
                if(at_inp%nat/=ishape(2)) then
                    call atom_deallocate_old(at_out,fat=.true.)
                endif
            endif
            if(.not. allocated(at_out%fat)) then
                call atom_allocate_old(at_out,at_inp%nat,at_inp%natim,at_inp%nfp,fat=.true.)
            endif
            do iat=1,at_inp%nat
                at_out%fat(1,iat)=at_inp%fat(1,iat)
                at_out%fat(2,iat)=at_inp%fat(2,iat)
                at_out%fat(3,iat)=at_inp%fat(3,iat)
            enddo
        endif
    else
        err_mess='ERROR: cannot copy typ_atoms%fat when source is not allocated:'
        if(prsnt) then
            write(*,'(a,1x,a)') trim(err_mess),trim(str_message)
            stop
        endif
    endif
    !copying array at_inp%bemoved to at_out%bemoved
    if(present(bemoved)) then ; prsnt=bemoved ; else ; prsnt=.false. ;  endif
    if(allocated(at_inp%bemoved)) then
        if(all_of_them .or. prsnt) then
            if(allocated(at_out%bemoved)) then
                ishape(1:2)=shape(at_out%bemoved)
                if(at_inp%nat/=ishape(2)) then
                    call atom_deallocate_old(at_out,bemoved=.true.)
                endif
            endif
            if(.not. allocated(at_out%bemoved)) then
                call atom_allocate_old(at_out,at_inp%nat,at_inp%natim,at_inp%nfp,bemoved=.true.)
            endif
            do iat=1,at_inp%nat
                at_out%bemoved(1,iat)=at_inp%bemoved(1,iat)
                at_out%bemoved(2,iat)=at_inp%bemoved(2,iat)
                at_out%bemoved(3,iat)=at_inp%bemoved(3,iat)
            enddo
        endif
    else
        err_mess='ERROR: cannot copy typ_atoms%bemoved when source is not allocated:'
        if(prsnt) then
            write(*,'(a,1x,a)') trim(err_mess),trim(str_message)
            stop
        endif
    endif
    !copying array at_inp%qat to at_out%qat
    if(present(qat)) then ; prsnt=qat ; else ; prsnt=.false. ;  endif
    if(allocated(at_inp%qat)) then
        if(all_of_them .or. prsnt) then
            if(allocated(at_out%qat)) then
                ishape(1:1)=shape(at_out%qat)
                if(at_inp%nat/=ishape(1)) then
                    call atom_deallocate_old(at_out,qat=.true.)
                endif
            endif
            if(.not. allocated(at_out%qat)) then
                call atom_allocate_old(at_out,at_inp%nat,at_inp%natim,at_inp%nfp,qat=.true.)
            endif
            do iat=1,at_inp%nat
                at_out%qat(iat)=at_inp%qat(iat)
            enddo
        endif
    else
        err_mess='ERROR: cannot copy typ_atoms%qat when source is not allocated:'
        if(prsnt) then
            write(*,'(a,1x,a)') trim(err_mess),trim(str_message)
            stop
        endif
    endif
    !copying array at_inp%zat to at_out%zat
    if(present(zat)) then ; prsnt=zat ; else ; prsnt=.false. ;  endif
    if(allocated(at_inp%zat)) then
        if(all_of_them .or. prsnt) then
            if(allocated(at_out%zat)) then
                ishape(1:1)=shape(at_out%zat)
                if(at_inp%nat/=ishape(1)) then
                    call atom_deallocate_old(at_out,zat=.true.)
                endif
            endif
            if(.not. allocated(at_out%zat)) then
                call atom_allocate_old(at_out,at_inp%nat,at_inp%natim,at_inp%nfp,zat=.true.)
            endif
            do iat=1,at_inp%nat
                at_out%zat(iat)=at_inp%zat(iat)
            enddo
        endif
    else
        err_mess='ERROR: cannot copy typ_atoms%zat when source is not allocated:'
        if(prsnt) then
            write(*,'(a,1x,a)') trim(err_mess),trim(str_message)
            stop
        endif
    endif
    !copying array at_inp%rcov to at_out%rcov
    if(present(rcov)) then ; prsnt=rcov ; else ; prsnt=.false. ;  endif
    if(allocated(at_inp%rcov)) then
        if(all_of_them .or. prsnt) then
            if(allocated(at_out%rcov)) then
                ishape(1:1)=shape(at_out%rcov)
                if(at_inp%nat/=ishape(1)) then
                    call atom_deallocate_old(at_out,rcov=.true.)
                endif
            endif
            if(.not. allocated(at_out%rcov)) then
                call atom_allocate_old(at_out,at_inp%nat,at_inp%natim,at_inp%nfp,rcov=.true.)
            endif
            do iat=1,at_inp%nat
                at_out%rcov(iat)=at_inp%rcov(iat)
            enddo
        endif
    else
        err_mess='ERROR: cannot copy typ_atoms%rcov when source is not allocated:'
        if(prsnt) then
            write(*,'(a,1x,a)') trim(err_mess),trim(str_message)
            stop
        endif
    endif
    !copying array at_inp%fp to at_out%fp
    if(present(fp)) then ; prsnt=fp ; else ; prsnt=.false. ;  endif
    if(allocated(at_inp%fp)) then
        if(all_of_them .or. prsnt) then
            if(allocated(at_out%fp)) then
                ishape(1:1)=shape(at_out%fp)
                if(at_inp%nat/=ishape(1)) then
                    call atom_deallocate_old(at_out,fp=.true.)
                endif
            endif
            if(.not. allocated(at_out%fp)) then
                call atom_allocate_old(at_out,at_inp%nat,at_inp%natim,at_inp%nfp)
            endif
            do ifp=1,at_inp%nfp
                at_out%fp(ifp)=at_inp%fp(ifp)
            enddo
        endif
    else
        err_mess='ERROR: cannot copy typ_atoms%fp when source is not allocated:'
        if(prsnt) then
            write(*,'(a,1x,a)') trim(err_mess),trim(str_message)
            stop
        endif
    endif
    !copying array at_inp%itypat to at_out%itypat
    if(present(typat)) then ; prsnt=typat ; else ; prsnt=.false. ;  endif
    if(allocated(at_inp%itypat)) then
        if(all_of_them .or. prsnt) then
            if(allocated(at_out%itypat)) then
                ishape(1:1)=shape(at_out%itypat)
                if(at_inp%nat/=ishape(1)) then
                    call atom_deallocate_old(at_out,typat=.true.)
                endif
            endif
            if(.not. allocated(at_out%itypat)) then
                call atom_allocate_old(at_out,at_inp%nat,at_inp%natim,at_inp%nfp,typat=.true.)
            endif
            do iat=1,at_inp%nat
                at_out%itypat(iat)=at_inp%itypat(iat)
            enddo
        endif
    else
        err_mess='ERROR: cannot copy typ_atoms%itypat when source is not allocated:'
        if(prsnt) then
            write(*,'(a,1x,a)') trim(err_mess),trim(str_message)
            stop
        endif
    endif
    !copying array at_inp%trial_ref_energy to at_out%trial_ref_energy
    if(present(trial_ref_energy)) then ; prsnt=trial_ref_energy ; else ; prsnt=.false. ;  endif
    if(allocated(at_inp%trial_ref_energy)) then
        if(all_of_them .or. prsnt) then
            if(allocated(at_out%trial_ref_energy)) then
                ishape(1:1)=shape(at_out%trial_ref_energy)
                if(at_inp%ntrial/=ishape(1)) then
                    call atom_deallocate_old(at_out,trial_ref_energy=.true.)
                endif
            endif
            if(.not. allocated(at_out%trial_ref_energy)) then
                call atom_allocate_old(at_out,at_inp%nat,at_inp%natim,at_inp%nfp,ntrial=at_inp%ntrial,trial_ref_energy=.true.)
            endif
            do iat=1,at_inp%ntrial
                at_out%trial_ref_energy(iat)=at_inp%trial_ref_energy(iat)
            enddo
        endif
    else
        err_mess='ERROR: cannot copy typ_atoms%trial_ref_energy when source is not allocated:'
        if(prsnt) then
            write(*,'(a,1x,a)') trim(err_mess),trim(str_message)
            stop
        endif
    endif
    !copying array at_inp%trial_ref_disp to at_out%trial_ref_disp
    if(present(trial_ref_disp)) then ; prsnt=trial_ref_disp ; else ; prsnt=.false. ;  endif
    if(allocated(at_inp%trial_ref_disp)) then
        if(all_of_them .or. prsnt) then
            if(allocated(at_out%trial_ref_disp)) then
                ishape(1:2)=shape(at_out%trial_ref_disp)
                if(at_inp%ntrial/=ishape(1)) then
                    call atom_deallocate_old(at_out,trial_ref_disp=.true.)
                endif
            endif
            if(.not. allocated(at_out%trial_ref_disp)) then
                call atom_allocate_old(at_out,at_inp%nat,at_inp%natim,at_inp%nfp,ntrial=at_inp%ntrial,trial_ref_disp=.true.)
            endif
            do iat=1,at_inp%ntrial
                at_out%trial_ref_disp(1,iat)=at_inp%trial_ref_disp(1,iat)
                at_out%trial_ref_disp(2,iat)=at_inp%trial_ref_disp(2,iat)
                at_out%trial_ref_disp(3,iat)=at_inp%trial_ref_disp(3,iat)
            enddo
        endif
    else
        err_mess='ERROR: cannot copy typ_atoms%trial_ref_disp when source is not allocated:'
        if(prsnt) then
            write(*,'(a,1x,a)') trim(err_mess),trim(str_message)
            stop
        endif
    endif
    !copying array at_inp%trial_ref_nat to at_out%trial_ref_nat
    if(present(trial_ref_nat)) then ; prsnt=trial_ref_nat ; else ; prsnt=.false. ;  endif
    if(allocated(at_inp%trial_ref_nat)) then
        if(all_of_them .or. prsnt) then
            if(allocated(at_out%trial_ref_nat)) then
                ishape(1:1)=shape(at_out%trial_ref_nat)
                if(at_inp%ntrial/=ishape(1)) then
                    call atom_deallocate_old(at_out,trial_ref_nat=.true.)
                endif
            endif
            if(.not. allocated(at_out%trial_ref_nat)) then
                call atom_allocate_old(at_out,at_inp%nat,at_inp%natim,at_inp%nfp,ntrial=at_inp%ntrial,trial_ref_nat=.true.)
            endif
            do iat=1,at_inp%ntrial
                at_out%trial_ref_nat(iat)=at_inp%trial_ref_nat(iat)
            enddo
        endif
    else
        err_mess='ERROR: cannot copy typ_atoms%trial_ref_nat when source is not allocated:'
        if(prsnt) then
            write(*,'(a,1x,a)') trim(err_mess),trim(str_message)
            stop
        endif
    endif
end subroutine atom_copy_old
!*****************************************************************************************
subroutine atom_build_periodic_images(atoms,rcut)
    implicit none
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(in):: rcut
    !local variables
    integer:: nx, ny, nz, ix, iy, iz, iat, iatim, natim
    if(trim(atoms%boundcond)/='bulk') then
        write(*,'(a)') 'ERROR: atom_build_periodic_images works only for bulk BC'
        stop
    endif
    call n_rep_dim_alborz(atoms%cellvec,rcut,nx,ny,nz)
    natim=atoms%nat*(2*nx+1)*(2*ny+1)*(2*nz+1)
    !write(*,*) 'nx,ny,nz,natim',nx,ny,nz,natim
    call atom_allocate_old(atoms,atoms%nat,natim,atoms%nfp)
    iatim=0
    do iat=1,atoms%nat
        iatim=iatim+1
        atoms%ratim(1,iatim)=atoms%rat(1,iat)
        atoms%ratim(2,iatim)=atoms%rat(2,iat)
        atoms%ratim(3,iatim)=atoms%rat(3,iat)
    enddo
    do iz=-nz,nz
    do iy=-ny,ny
    do ix=-nx,nx
        if(ix==0 .and. iy==0 .and. iz==0) cycle
        do iat=1,atoms%nat
            iatim=iatim+1
            atoms%ratim(1,iatim)=atoms%rat(1,iat)+ix*atoms%cellvec(1,1)+iy*atoms%cellvec(1,2)+iz*atoms%cellvec(1,3)
            atoms%ratim(2,iatim)=atoms%rat(2,iat)+ix*atoms%cellvec(2,1)+iy*atoms%cellvec(2,2)+iz*atoms%cellvec(2,3)
            atoms%ratim(3,iatim)=atoms%rat(3,iat)+ix*atoms%cellvec(3,1)+iy*atoms%cellvec(3,2)+iz*atoms%cellvec(3,3)
        enddo
    enddo
    enddo
    enddo
    if(iatim/=atoms%natim) then
        write(*,'(a,2i7)') 'ERROR: iatim/=atoms%natim in atom_build_periodic_images',iatim,atoms%natim
        stop
    endif
end subroutine atom_build_periodic_images
!*****************************************************************************************
subroutine set_typat(atoms)
    implicit none
    type(typ_atoms), intent(inout):: atoms
    !local variables
    integer:: iat, itypat
    logical:: type_is_new
    if(atoms%nat<0) then
        stop 'ERROR: nat<0 in set_typat'
    endif
    atoms%ntypat=1
    atoms%stypat(1)=trim(atoms%sat(1))
    atoms%ltypat(1)=1
    atoms%itypat(1)=1
    do iat=2,atoms%nat
        !write(*,*) 'TYPE ',trim(atoms%sat(iat)),trim(atoms%stypat(atoms%ntypat))
        type_is_new=.true.
        do itypat=1,atoms%ntypat
            if(trim(atoms%sat(iat))==trim(atoms%stypat(itypat))) then
                type_is_new=.false.
                exit
            endif
        enddo
        if(type_is_new) then
            atoms%ntypat=atoms%ntypat+1
            atoms%ltypat(atoms%ntypat)=1
            atoms%stypat(atoms%ntypat)=trim(atoms%sat(iat))
        else
            atoms%ltypat(itypat)=atoms%ltypat(itypat)+1
        endif
        atoms%itypat(iat)=itypat
    enddo
end subroutine set_typat
!*****************************************************************************************
subroutine set_ndof(atoms)
    implicit none
    type(typ_atoms), intent(inout):: atoms
    !local variables
    integer:: iat
    atoms%ndof=0
    do iat=1,atoms%nat
        if(atoms%bemoved(1,iat)) atoms%ndof=atoms%ndof+1
        if(atoms%bemoved(2,iat)) atoms%ndof=atoms%ndof+1
        if(atoms%bemoved(3,iat)) atoms%ndof=atoms%ndof+1
    enddo
end subroutine set_ndof
!*****************************************************************************************
subroutine bemoved2string(bemoved,str_motion)
    implicit none
    logical:: bemoved(3)
    character(3):: str_motion
    if(bemoved(1)) then
        str_motion='T'
    else
        str_motion='F'
    endif
    if(bemoved(2)) then
        str_motion=trim(str_motion)//'T'
    else
        str_motion=trim(str_motion)//'F'
    endif
    if(bemoved(3)) then
        str_motion=trim(str_motion)//'T'
    else
        str_motion=trim(str_motion)//'F'
    endif
end subroutine bemoved2string
!*****************************************************************************************
subroutine string2bemoved(str_motion,bemoved)
    implicit none
    character(3):: str_motion
    logical:: bemoved(3)
    if(str_motion(1:1)=='T') then
        bemoved(1)=.true.
    elseif(str_motion(1:1)=='F') then
        bemoved(1)=.false.
    else
        stop 'ERROR: incorrect symbol for atomic motion'
    endif
    if(str_motion(2:2)=='T') then
        bemoved(2)=.true.
    elseif(str_motion(2:2)=='F') then
        bemoved(2)=.false.
    else
        stop 'ERROR: incorrect symbol for atomic motion'
    endif
    if(str_motion(3:3)=='T') then
        bemoved(3)=.true.
    elseif(str_motion(3:3)=='F') then
        bemoved(3)=.false.
    else
        stop 'ERROR: incorrect symbol for atomic motion'
    endif
end subroutine string2bemoved
!*****************************************************************************************
subroutine atoms_all_writexyz(filename,fn_position,atoms_all,strkey)
    implicit none
    character(*), intent(in):: filename, fn_position, strkey
    type(typ_atoms_all), intent(in):: atoms_all
    !local variables
    integer:: iconf, iat, ios
    real(8):: cv(3,3), x, y, z
    character(3):: str_motion
    character(4):: tch1
    character(40):: tch2
    if(len_trim(strkey)>36) stop 'ERROR: strkey too in atoms_all_writexyz'
    if(trim(fn_position)=='new') then
        open(unit=1358, file=filename,status='replace',iostat=ios)
    elseif(trim(fn_position)=='append') then
        open(unit=1358, file=filename,status='old',position='append',iostat=ios)
    else
        write(*,'(2a)') 'ERROR: fn_position is unknown, ',trim(fn_position)
    endif
    if(ios/=0) then;write(*,'(2a)') 'ERROR: failure openning ',trim(filename);stop;endif
    !write(1358,'(i6,2x,a)')   nat,trim(comment)
    !write(1358,'(a,2x,6(e24.15))') trim(boundcond),cellvec(1,1),cellvec(2,2),cellvec(3,3),cellvec(1,2),cellvec(1,3),cellvec(2,3)
    !do iat=1,nat
    !    x=rat(1,iat) !*0.529177d0
    !    y=rat(2,iat) !*0.529177d0
    !    z=rat(3,iat) !*0.529177d0
    !    write(1358,'(a5,2x,3es24.15,2x,3l1)') trim(sat(iat)),x,y,z, &
    !        bemoved(1,iat),bemoved(2,iat),bemoved(3,iat)
    !enddo
    do iconf=1,atoms_all%nconf
        write(tch1,'(i4.4)') iconf 
        tch2=trim(strkey)//tch1
        write(1358,'(i6,1x,e24.15,1x,a)') atoms_all%atoms%nat,atoms_all%epotall(iconf),trim(tch2)
        cv(1,1)=atoms_all%atoms%cellvec(1,1)
        cv(2,2)=atoms_all%atoms%cellvec(2,2)
        cv(3,3)=atoms_all%atoms%cellvec(3,3)
        cv(1,2)=atoms_all%atoms%cellvec(1,2)
        cv(1,3)=atoms_all%atoms%cellvec(1,3)
        cv(2,3)=atoms_all%atoms%cellvec(2,3)
        write(1358,'(a,2x,6e22.13)') trim(atoms_all%atoms%boundcond),cv(1,1),cv(2,2),cv(3,3),cv(1,2),cv(1,3),cv(2,3)
        do iat=1,atoms_all%atoms%nat
            call bemoved2string(atoms_all%atoms%bemoved(1,iat),str_motion)
            x=atoms_all%ratall(1,iat,iconf)
            y=atoms_all%ratall(2,iat,iconf)
            z=atoms_all%ratall(3,iat,iconf)
            write(1358,'(a5,2x,3e22.13,2x,a3)') trim(atoms_all%atoms%sat(iat)),x,y,z,str_motion
        enddo
    enddo
    close(1358)
end subroutine atoms_all_writexyz
!*****************************************************************************************
subroutine atom_normalizevector(nat,bemoved,vxyz)
    implicit none
    integer, intent(in):: nat
    logical, intent(in):: bemoved(3,nat)
    real(8), intent(inout):: vxyz(3,nat)
    !local variables
    integer:: iat, ixyz, i
    real(8):: vnrm
    !integer, save:: icall=0
    !icall=icall+1
    vnrm=0.d0
    do iat=1,nat
        if(bemoved(1,iat)) vnrm=vnrm+vxyz(1,iat)**2
        if(bemoved(2,iat)) vnrm=vnrm+vxyz(2,iat)**2
        if(bemoved(3,iat)) vnrm=vnrm+vxyz(3,iat)**2
    enddo
    !do i=1,3*nat
    !    iat=(i-1)/3+1
    !    ixyz=mod(i-1,3)+1
    !    if(bemoved(ixyz,iat)) vnrm=vnrm+vxyz(ixyz,iat)**2
    !enddo
    vnrm=sqrt(vnrm)
    !write(21,'(i5,es24.15)') icall,vnrm
    if(vnrm==0.d0) stop 'ERROR: vnrm=0 in atom_normalizevector'
    do iat=1,nat
        if(bemoved(1,iat)) vxyz(1,iat)=vxyz(1,iat)/vnrm
        if(bemoved(2,iat)) vxyz(2,iat)=vxyz(2,iat)/vnrm
        if(bemoved(3,iat)) vxyz(3,iat)=vxyz(3,iat)/vnrm
    enddo
    !do i=1,3*nat
    !    iat=(i-1)/3+1
    !    ixyz=mod(i-1,3)+1
    !    if(bemoved(ixyz,iat)) vxyz(ixyz,iat)=vxyz(ixyz,iat)/vnrm
    !enddo
end subroutine atom_normalizevector
!*****************************************************************************************
function atom_ddot(nat,vec1,vec2,bemoved) result(res)
    implicit none
    integer, intent(in):: nat
    real(8), intent(in):: vec1(3,nat), vec2(3,nat)
    logical, intent(in):: bemoved(3,nat)
    real(8):: res
    !local variables
    integer:: iat
    res=0.d0
    do iat=1,nat
        if(bemoved(1,iat)) res=res+vec1(1,iat)*vec2(1,iat)
        if(bemoved(2,iat)) res=res+vec1(2,iat)*vec2(2,iat)
        if(bemoved(3,iat)) res=res+vec1(3,iat)*vec2(3,iat)
    enddo
end function atom_ddot
!*****************************************************************************************
subroutine atom_calnorm(nat,bemoved,vec,vnrm)
    implicit none
    integer, intent(in):: nat
    logical, intent(in):: bemoved(3,nat)
    real(8), intent(in):: vec(3,nat)
    real(8), intent(out):: vnrm
    !local variables
    integer:: iat
    vnrm=0.d0
    do iat=1,nat
        if(bemoved(1,iat)) vnrm=vnrm+vec(1,iat)**2
        if(bemoved(2,iat)) vnrm=vnrm+vec(2,iat)**2
        if(bemoved(3,iat)) vnrm=vnrm+vec(3,iat)**2
    enddo
    vnrm=sqrt(vnrm)
end subroutine atom_calnorm
!*****************************************************************************************
subroutine atom_calmaxforcecomponent(nat,bemoved,vec,vmax)
    implicit none
    integer, intent(in):: nat
    logical, intent(in):: bemoved(3,nat)
    real(8), intent(in):: vec(3,nat)
    real(8), intent(out):: vmax
    !local variables
    integer:: iat
    vmax=0.d0
    do iat=1,nat
        if(bemoved(1,iat)) vmax=max(vmax,abs(vec(1,iat)))
        if(bemoved(2,iat)) vmax=max(vmax,abs(vec(2,iat)))
        if(bemoved(3,iat)) vmax=max(vmax,abs(vec(3,iat)))
    enddo
end subroutine atom_calmaxforcecomponent
!*****************************************************************************************
!subroutine readatomicpositions
!    use mod_atoms, only:nat,rat,itypat,indatinp,unfrozenat,ntypat,atomnames,unitsfnposinp,cellvec,fnposinp,maxntypat
!    use mod_unitsconversion, only:angstroemtobohr
!    implicit none
!    integer:: iat,ios,ittt,ityp,iatf,iatnf,istat, natnf
!    character(20):: tatomnames
!    logical:: condition
!    real(8), allocatable:: rat_t(:,:)
!    integer, allocatable:: itypat_t(:)
!    !integer:: ntypat
!    if(.not. allocated(rat)) stop 'ERROR: reading rat because is not allocated.'
!    !read atomic positions and determine the type of it 
!    open(unit=1,file=fnposinp,status='old',iostat=ios)
!    if(ios/=0) stop 'ERROR: openning fnposinp.'
!    write(*,*) 'reading atomic positions from ',fnposinp
!    read(1,*)  ittt,natnf,unitsfnposinp 
!    write(*,*) 'natnf,unitsfnposinp',natnf,unitsfnposinp
!    stop 'CORRECT reading cellvec'
!    !read(1,*)  cellvec(1,1),cellvec(1,2),cellvec(1,3)
!    !write(*,'(a10,3f30.16)') 'cell',cell(1:3)
!    !if Angstroem convert to Bohr.
!    if(unitsfnposinp=='angstroem') cellvec(1:3,1:3)=cellvec(1:3,1:3)*angstroemtobohr
!    ntypat=0
!    allocate(rat_t(3,nat),stat=istat)
!    if(istat/=0) write(*,*) 'ERROR: allocating array rat_t.'
!    allocate(itypat_t(nat),stat=istat)
!    if(istat/=0) write(*,*) 'ERROR: allocating array itypat_t.'
!    do iat=1,nat
!        condition=.true.
!        read(1,*) tatomnames,rat_t(1,iat),rat_t(2,iat),rat_t(3,iat),unfrozenat(iat)
!        write(*,'(a8,3e30.12,l4)')  tatomnames,rat_t(1:3,iat),unfrozenat(iat)
!        do ityp=1,ntypat
!            if(tatomnames==atomnames(ityp)) then
!                itypat_t(iat)=ityp
!                condition=.false. 
!                exit
!            endif
!        enddo
!        if(condition) then
!            ntypat=ntypat+1
!            if(ntypat>maxntypat) stop 'more than maxntypat atomnames not permitted'
!            atomnames(ntypat)=tatomnames
!            itypat_t(iat)=ntypat
!        endif
!        if (unitsfnposinp=='angstroem') then   !if Angstroem convert to Bohr
!            rat_t(1:3,iat)=rat_t(1:3,iat)*angstroemtobohr
!        elseif(unitsfnposinp=='atomic' .or. unitsfnposinp=='bohr') then
!        else
!            write(*,*) 'length unitsfnposinp in input file unrecognized'
!            write(*,*) 'recognizable unitsfnposinp are angstroem or atomic = bohr'
!            stop
!        endif
!    enddo
!    close(1)
!    write(*,*) 'ntypat=',ntypat
!    iatnf=0 
!    do iat=1,nat
!        if(unfrozenat(iat)) iatnf=iatnf+1
!    enddo
!    if(iatnf/=natnf) then
!        write(*,*) 'iatnf,natnf',iatnf,natnf
!        stop 'ERROR: iatnf is not equal to natnf.'
!    endif
!    iatnf=0;iatf=natnf
!    if(natnf<nat) then
!        do iat=1,nat
!            if(unfrozenat(iat)) then
!                iatnf=iatnf+1
!                indatinp(iatnf)=iat
!                rat(1:3,iatnf)=rat_t(1:3,iat)
!                itypat(iatnf)=itypat_t(iat)
!            else
!                iatf=iatf+1
!                indatinp(iatf)=iat
!                rat(1:3,iatf)=rat_t(1:3,iat)
!                itypat(iatf)=itypat_t(iat)
!            endif
!        enddo
!    endif
!    unfrozenat(1:natnf)=.true.
!    unfrozenat(natnf+1:nat)=.false.
!    deallocate(rat_t,stat=istat)
!    if(istat/=0) write(*,*) 'ERROR: deallocating array rat_t.'
!    deallocate(itypat_t,stat=istat)
!    if(istat/=0) write(*,*) 'ERROR: deallocating array itypat_t.'
!end subroutine readatomicpositions
!*****************************************************************************************
!subroutine initqat
!    use mod_atoms, only:nat,qat,itypat,ntypat,qtypat,atomnames
!    implicit none
!    integer:: ityp,iat
!    character(8):: tatnam
!    do ityp=1,ntypat
!        tatnam=atomnames(ityp)
!        if(tatnam=='Li') then
!            qtypat(ityp)=1.d0
!        elseif(tatnam=='Be') then
!            qtypat(ityp)=2.d0
!        elseif(tatnam=='O') then
!            qtypat(ityp)=-2.d0
!        elseif(tatnam=='F') then
!            qtypat(ityp)=-1.d0
!        elseif(tatnam=='Na') then
!            qtypat(ityp)=1.d0
!        elseif(tatnam=='Mg') then
!            qtypat(ityp)=2.d0
!        elseif(tatnam=='S') then
!            qtypat(ityp)=-2.d0
!        elseif(tatnam=='Cl') then
!            qtypat(ityp)=-1.d0
!        elseif(tatnam=='K') then
!            qtypat(ityp)=1.d0
!        elseif(tatnam=='Ca') then
!            qtypat(ityp)=2.d0
!        elseif(tatnam=='Se') then
!            qtypat(ityp)=-2.d0
!        elseif(tatnam=='Br') then
!            qtypat(ityp)=-1.d0
!        elseif(tatnam=='Rb') then
!            qtypat(ityp)=1.d0
!        elseif(tatnam=='Sr') then
!            qtypat(ityp)=2.d0
!        elseif(tatnam=='Te') then
!            qtypat(ityp)=-2.d0
!        elseif(tatnam=='I') then
!            qtypat(ityp)=-1.d0
!        elseif(tatnam=='Cs') then
!            qtypat(ityp)=1.d0
!        elseif(tatnam=='Ba') then
!            qtypat(ityp)=2.d0
!        elseif(tatnam=='At') then
!            qtypat(ityp)=-1.d0
!        else 
!            write(*,*) 'To assign charge of atom, ',trim(tatnam),'is not recognized.'
!        endif
!    enddo
!    do iat=1,nat
!        qat(iat)=qtypat(itypat(iat))
!    enddo
!end subroutine initqat
!*****************************************************************************************
!subroutine writeposout
!    use mod_atoms, only:nat,cellvec,indatinp,unitsfnposinp,rat,unfrozenat,atomnames,itypat,fnposout
!    use mod_unitsconversion, only:bohrtoangstroem
!    implicit none
!    integer:: iat,ios,istat,natnf
!    character(8):: tt
!    real(8), allocatable:: rat_t(:,:)
!    integer, allocatable:: itypat_t(:)
!    logical, allocatable:: unfrozenat_t(:)
!    allocate(rat_t(3,nat),stat=istat)
!    if(istat/=0) write(*,*) 'ERROR: allocating array rat_t.'
!    allocate(itypat_t(nat),stat=istat)
!    if(istat/=0) write(*,*) 'ERROR: allocating array itypat_t.'
!    allocate(unfrozenat_t(nat),stat=istat)
!    if(istat/=0) write(*,*) 'ERROR: allocating array unfrozenat_t.'
!
!    !write(fn,'(i5.5)') f
!    !fnposout='posout'//fn//'.xyz'
!    if(fnposout=='none') stop 'ERROR: reading nat because fnposinp is not inialized.'
!    write(*,'(2a)') 'Openning atomic position output file, fnposout=',fnposout
!    open(unit=1,file=fnposout,status='replace',iostat=ios)
!    if(ios/=0) stop 'ERROR: openning output file.'
!    write(1,*) nat,natnf,unitsfnposinp
!    stop 'CORRECT writing cell'
!    !write(1,*) cell(1:3)*bohrtoangstroem
!    do iat=1,nat
!        !rat_t(1:3,iat)=rat(1:3,indatinp(iat))
!        !itypat_t(iat)=itypat(indatinp(iat))
!        !unfrozenat_t(iat)=unfrozenat(indatinp(iat))
!        rat_t(1:3,indatinp(iat))=rat(1:3,iat)
!        itypat_t(indatinp(iat))=itypat(iat)
!        unfrozenat_t(indatinp(iat))=unfrozenat(iat)
!    enddo
!    do iat=1,nat
!        tt='        '
!        tt=tt(1:8-len(trim(atomnames(itypat_t(iat)))))//trim(atomnames(itypat_t(iat)))
!        write(1,'(a8,3e30.12,l4)') tt,rat_t(1:3,iat)*bohrtoangstroem,unfrozenat_t(iat)
!    enddo
!    close(1)
!    deallocate(rat_t,stat=istat)
!    if(istat/=0) write(*,*) 'ERROR: deallocating array rat_t.'
!    deallocate(itypat_t,stat=istat)
!    if(istat/=0) write(*,*) 'ERROR: deallocating array itypat_t.'
!    deallocate(unfrozenat_t,stat=istat)
!    if(istat/=0) write(*,*) 'ERROR: deallocating array unfrozenat_t.'
!end subroutine writeposout
!*****************************************************************************************
subroutine checkallatomstobeincell(nat,rat,cellvec,allatomsincell)
    implicit none
    integer:: iat,nat
    real(8):: rat(3,nat),cellvec(3,3)
    logical:: allatomsincell
    stop 'CORRECT working with cellvec'
    allatomsincell=.true.
    do iat=1,nat
    if(rat(1,iat)<0.d0) then;allatomsincell=.false.;exit;endif
    if(rat(2,iat)<0.d0) then;allatomsincell=.false.;exit;endif
    if(rat(3,iat)<0.d0) then;allatomsincell=.false.;exit;endif
    if(rat(1,iat)>cellvec(1,1)) then;allatomsincell=.false.;exit;endif
    if(rat(2,iat)>cellvec(2,2)) then;allatomsincell=.false.;exit;endif
    if(rat(3,iat)>cellvec(3,3)) then;allatomsincell=.false.;exit;endif
    enddo
    if(.not. allatomsincell) write(*,*) 'Not all atoms are inside the cell.'
end subroutine checkallatomstobeincell
!*****************************************************************************************
subroutine determinexyzminmax(nat,rat,cellvec,xmin,ymin,zmin,xmax,ymax,zmax)
    implicit none
    integer:: iat,nat
    real(8):: rat(3,nat),cellvec(3,3)
    real(8):: xmin,ymin,zmin,xmax,ymax,zmax,x,y,z
    !stop 'CORRECT working with cellvec'
    xmin=1.d20;xmax=-1.d20;ymin=1.d20;ymax=-1.d20;zmin=1.d20;zmax=-1.d20
    do iat=1,nat
        x=rat(1,iat);y=rat(2,iat);z=rat(3,iat)
        xmin=min(x,xmin);ymin=min(y,ymin);zmin=min(z,zmin)
        xmax=max(x,xmax);ymax=max(y,ymax);zmax=max(z,zmax)
    enddo
    if(xmin<0.d0) write(*,*) 'xmin<0.d0'
    if(ymin<0.d0) write(*,*) 'ymin<0.d0'
    if(zmin<0.d0) write(*,*) 'zmin<0.d0'
    if(xmax>cellvec(1,1)) write(*,*) 'xmax>cell(1)'
    if(ymax>cellvec(2,2)) write(*,*) 'ymax>cell(2)'
    if(zmax>cellvec(3,3)) write(*,*) 'zmax>cell(3)'
end subroutine determinexyzminmax
!*****************************************************************************************
subroutine set_rcov(atoms)
    use mod_processors, only: iproc
    use mod_const, only: bohr2ang
    implicit none
    type(typ_atoms), intent(inout) :: atoms
    !local variables
    integer:: iat
    if(.not. allocated(atoms%sat)) stop 'ERROR: in set_rcov: sat not allocated.'
    if(.not. allocated(atoms%rcov)) stop 'ERROR: in set_rcov: rcov not allocated.'
    do iat=1,atoms%nat
        if(trim(atoms%sat(iat))=='LJ') then
            atoms%rcov(iat)=2.d0**(1.d0/6.d0)/2.d0
        else if(trim(atoms%sat(iat))=='H') then
            atoms%rcov(iat)=0.31d0
        else if(trim(atoms%sat(iat))=='He') then
            atoms%rcov(iat)=0.28d0
        else if(trim(atoms%sat(iat))=='Li') then
            atoms%rcov(iat)=1.28d0
        else if(trim(atoms%sat(iat))=='Be') then
            atoms%rcov(iat)=0.96d0
        else if(trim(atoms%sat(iat))=='B') then
            atoms%rcov(iat)=0.84d0
        else if(trim(atoms%sat(iat))=='C') then
            atoms%rcov(iat)=0.76d0
        else if(trim(atoms%sat(iat))=='N') then
            atoms%rcov(iat)=0.71d0
        else if(trim(atoms%sat(iat))=='O') then
            atoms%rcov(iat)=0.66d0
        else if(trim(atoms%sat(iat))=='F') then
            atoms%rcov(iat)=0.57d0
        else if(trim(atoms%sat(iat))=='Ne') then
            atoms%rcov(iat)=0.58d0
        else if(trim(atoms%sat(iat))=='Na') then
            atoms%rcov(iat)=1.66d0
        else if(trim(atoms%sat(iat))=='Mg') then
            atoms%rcov(iat)=1.41d0
        else if(trim(atoms%sat(iat))=='Al') then
            atoms%rcov(iat)=1.21d0
        else if(trim(atoms%sat(iat))=='Si') then
            atoms%rcov(iat)=1.11d0
        else if(trim(atoms%sat(iat))=='P') then
            atoms%rcov(iat)=1.07d0
        else if(trim(atoms%sat(iat))=='S') then
            atoms%rcov(iat)=1.05d0
        else if(trim(atoms%sat(iat))=='Cl') then
            atoms%rcov(iat)=1.02d0
        else if(trim(atoms%sat(iat))=='Ar') then
            atoms%rcov(iat)=1.06d0
        else if(trim(atoms%sat(iat))=='K') then
            atoms%rcov(iat)=2.03d0
        else if(trim(atoms%sat(iat))=='Ca') then
            atoms%rcov(iat)=1.76d0
        else if(trim(atoms%sat(iat))=='Sc') then
            atoms%rcov(iat)=1.70d0
        else if(trim(atoms%sat(iat))=='Ti') then
            atoms%rcov(iat)=1.60d0
        else if(trim(atoms%sat(iat))=='V') then
            atoms%rcov(iat)=1.53d0
        else if(trim(atoms%sat(iat))=='Cr') then
            atoms%rcov(iat)=1.39d0
        else if(trim(atoms%sat(iat))=='Mn') then
            atoms%rcov(iat)=1.39d0
        else if(trim(atoms%sat(iat))=='Fe') then
            atoms%rcov(iat)=1.32d0
        else if(trim(atoms%sat(iat))=='Co') then
            atoms%rcov(iat)=1.26d0
        else if(trim(atoms%sat(iat))=='Ni') then
            atoms%rcov(iat)=1.24d0
        else if(trim(atoms%sat(iat))=='Cu') then
            atoms%rcov(iat)=1.32d0
        else if(trim(atoms%sat(iat))=='Zn') then
            atoms%rcov(iat)=1.22d0
        else if(trim(atoms%sat(iat))=='Ga') then
            atoms%rcov(iat)=1.22d0
        else if(trim(atoms%sat(iat))=='Ge') then
            atoms%rcov(iat)=1.20d0
        else if(trim(atoms%sat(iat))=='As') then
            atoms%rcov(iat)=1.19d0
        else if(trim(atoms%sat(iat))=='Se') then
            atoms%rcov(iat)=1.20d0
        else if(trim(atoms%sat(iat))=='Br') then
            atoms%rcov(iat)=1.20d0
        else if(trim(atoms%sat(iat))=='Kr') then
            atoms%rcov(iat)=1.16d0
        else if(trim(atoms%sat(iat))=='Rb') then
            atoms%rcov(iat)=2.20d0
        else if(trim(atoms%sat(iat))=='Sr') then
            atoms%rcov(iat)=1.95d0
        else if(trim(atoms%sat(iat))=='Y') then
            atoms%rcov(iat)=1.90d0
        else if(trim(atoms%sat(iat))=='Zr') then
            atoms%rcov(iat)=1.75d0
        else if(trim(atoms%sat(iat))=='Nb') then
            atoms%rcov(iat)=1.64d0
        else if(trim(atoms%sat(iat))=='Mo') then
            atoms%rcov(iat)=1.54d0
        else if(trim(atoms%sat(iat))=='Tc') then
            atoms%rcov(iat)=1.47d0
        else if(trim(atoms%sat(iat))=='Ru') then
            atoms%rcov(iat)=1.46d0
        else if(trim(atoms%sat(iat))=='Rh') then
            atoms%rcov(iat)=1.42d0
        else if(trim(atoms%sat(iat))=='Pd') then
            atoms%rcov(iat)=1.39d0
        else if(trim(atoms%sat(iat))=='Ag') then
            atoms%rcov(iat)=1.45d0
        else if(trim(atoms%sat(iat))=='Cd') then
            atoms%rcov(iat)=1.44d0
        else if(trim(atoms%sat(iat))=='In') then
            atoms%rcov(iat)=1.42d0
        else if(trim(atoms%sat(iat))=='Sn') then
            atoms%rcov(iat)=1.39d0
        else if(trim(atoms%sat(iat))=='Sb') then
            atoms%rcov(iat)=1.39d0
        else if(trim(atoms%sat(iat))=='Te') then
            atoms%rcov(iat)=1.38d0
        else if(trim(atoms%sat(iat))=='I') then
            atoms%rcov(iat)=1.39d0
        else if(trim(atoms%sat(iat))=='Xe') then
            atoms%rcov(iat)=1.40d0
        else if(trim(atoms%sat(iat))=='Cs') then
            atoms%rcov(iat)=2.44d0
        else if(trim(atoms%sat(iat))=='Ba') then
            atoms%rcov(iat)=2.15d0
        else if(trim(atoms%sat(iat))=='La') then
            atoms%rcov(iat)=2.07d0
        else if(trim(atoms%sat(iat))=='Ce') then
            atoms%rcov(iat)=2.04d0
        else if(trim(atoms%sat(iat))=='Pr') then
            atoms%rcov(iat)=2.03d0
        else if(trim(atoms%sat(iat))=='Nd') then
            atoms%rcov(iat)=2.01d0
        else if(trim(atoms%sat(iat))=='Pm') then
            atoms%rcov(iat)=1.99d0
        else if(trim(atoms%sat(iat))=='Sm') then
            atoms%rcov(iat)=1.98d0
        else if(trim(atoms%sat(iat))=='Eu') then
            atoms%rcov(iat)=1.98d0
        else if(trim(atoms%sat(iat))=='Gd') then
            atoms%rcov(iat)=1.96d0
        else if(trim(atoms%sat(iat))=='Tb') then
            atoms%rcov(iat)=1.94d0
        else if(trim(atoms%sat(iat))=='Dy') then
            atoms%rcov(iat)=1.92d0
        else if(trim(atoms%sat(iat))=='Ho') then
            atoms%rcov(iat)=1.92d0
        else if(trim(atoms%sat(iat))=='Er') then
            atoms%rcov(iat)=1.89d0
        else if(trim(atoms%sat(iat))=='Tm') then
            atoms%rcov(iat)=1.90d0
        else if(trim(atoms%sat(iat))=='Yb') then
            atoms%rcov(iat)=1.87d0
        else if(trim(atoms%sat(iat))=='Lu') then
            atoms%rcov(iat)=1.87d0
        else if(trim(atoms%sat(iat))=='Hf') then
            atoms%rcov(iat)=1.75d0
        else if(trim(atoms%sat(iat))=='Ta') then
            atoms%rcov(iat)=1.70d0
        else if(trim(atoms%sat(iat))=='W') then
            atoms%rcov(iat)=1.62d0
        else if(trim(atoms%sat(iat))=='Re') then
            atoms%rcov(iat)=1.51d0
        else if(trim(atoms%sat(iat))=='Os') then
            atoms%rcov(iat)=1.44d0
        else if(trim(atoms%sat(iat))=='Ir') then
            atoms%rcov(iat)=1.41d0
        else if(trim(atoms%sat(iat))=='Pt') then
            atoms%rcov(iat)=1.36d0
        else if(trim(atoms%sat(iat))=='Au') then
            atoms%rcov(iat)=1.36d0
        else if(trim(atoms%sat(iat))=='Hg') then
            atoms%rcov(iat)=1.32d0
        else if(trim(atoms%sat(iat))=='Tl') then
            atoms%rcov(iat)=1.45d0
        else if(trim(atoms%sat(iat))=='Pb') then
            atoms%rcov(iat)=1.46d0
        else if(trim(atoms%sat(iat))=='Bi') then
            atoms%rcov(iat)=1.48d0
        else if(trim(atoms%sat(iat))=='Po') then
            atoms%rcov(iat)=1.40d0
        else if(trim(atoms%sat(iat))=='At') then
            atoms%rcov(iat)=1.50d0
        else if(trim(atoms%sat(iat))=='Rn') then
            atoms%rcov(iat)=1.50d0
        else if(trim(atoms%sat(iat))=='Fr') then
            atoms%rcov(iat)=2.60d0
        else if(trim(atoms%sat(iat))=='Ra') then
            atoms%rcov(iat)=2.21d0
        else if(trim(atoms%sat(iat))=='Ac') then
            atoms%rcov(iat)=2.15d0
        else if(trim(atoms%sat(iat))=='Th') then
            atoms%rcov(iat)=2.06d0
        else if(trim(atoms%sat(iat))=='Pa') then
            atoms%rcov(iat)=2.00d0
        else if(trim(atoms%sat(iat))=='U') then
            atoms%rcov(iat)=1.96d0
        else if(trim(atoms%sat(iat))=='Np') then
            atoms%rcov(iat)=1.90d0
        else if(trim(atoms%sat(iat))=='Pu') then
            atoms%rcov(iat)=1.87d0
        else if(trim(atoms%sat(iat))=='Am') then
            atoms%rcov(iat)=1.80d0
        else if(trim(atoms%sat(iat))=='Cm') then
            atoms%rcov(iat)=1.69d0
        else
            write(*,*) 'ERROR: no covalent radius stored for atomtype=',trim(atoms%sat(iat))
            stop
        endif
        !if(iproc==0) then
        !    write(*,'(a,a5,f10.5)') 'RCOV:',trim(atoms%sat(iat)),atoms%rcov(iat)
        !endif
        atoms%rcov(iat)=atoms%rcov(iat)/bohr2ang
    enddo
end subroutine set_rcov
!*****************************************************************************************
subroutine set_qat(atoms)
    use mod_processors, only: iproc
    implicit none
    type(typ_atoms), intent(inout) :: atoms
    !local variables
    integer:: iat, itypat
    if(atoms%nat<0) then
        stop 'ERROR: nat<0 in set_qat'
    endif
    do itypat=1,atoms%ntypat
        if(trim(atoms%stypat(itypat))=='Na') then
            atoms%qtypat(itypat)=1.d0
        elseif(trim(atoms%stypat(itypat))=='S') then
            atoms%qtypat(itypat)=-0.9d0
        elseif(trim(atoms%stypat(itypat))=='Cl') then
            atoms%qtypat(itypat)=-1.d0
        elseif(trim(atoms%stypat(itypat))=='Li') then
            atoms%qtypat(itypat)=1.0d0
        elseif(trim(atoms%stypat(itypat))=='Zn') then
            atoms%qtypat(itypat)=1.d0
        elseif(trim(atoms%stypat(itypat))=='Ti') then
            atoms%qtypat(itypat)=2.0d0
        elseif(trim(atoms%stypat(itypat))=='Zr') then
            atoms%qtypat(itypat)=4.0d0
        elseif(trim(atoms%stypat(itypat))=='Y') then
            atoms%qtypat(itypat)=3.0d0
        elseif(trim(atoms%stypat(itypat))=='Sn') then
            atoms%qtypat(itypat)=2.0d0
        elseif(trim(atoms%stypat(itypat))=='O') then
            atoms%qtypat(itypat)=-1.2d0
        elseif(trim(atoms%stypat(itypat))=='Si') then
            atoms%qtypat(itypat)=1.1d0
        elseif(trim(atoms%stypat(itypat))=='Pb') then
            atoms%qtypat(itypat)=0.9d0
        elseif(trim(atoms%stypat(itypat))=='Te') then
            atoms%qtypat(itypat)=-0.9d0
        elseif(trim(atoms%stypat(itypat))=='Sn') then
            atoms%qtypat(itypat)=0.9d0
        elseif(trim(atoms%stypat(itypat))=='Se') then
            atoms%qtypat(itypat)=-0.9d0
        elseif(trim(atoms%stypat(itypat))=='Ca') then
            atoms%qtypat(itypat)=1.3d0
        elseif(trim(atoms%stypat(itypat))=='F') then
            atoms%qtypat(itypat)=-0.65d0
        elseif(trim(atoms%stypat(itypat))=='Al') then
            atoms%qtypat(itypat)=3.0d0
        elseif(trim(atoms%stypat(itypat))=='W') then
            atoms%qtypat(itypat)=0.8d0
        elseif(trim(atoms%stypat(itypat))=='S') then
            atoms%qtypat(itypat)=-0.4d0
        elseif(trim(atoms%stypat(itypat))=='K') then
            atoms%qtypat(itypat)=1.d0
        elseif(trim(atoms%stypat(itypat))=='Br') then
            atoms%qtypat(itypat)=-1.d0
        elseif(trim(atoms%stypat(itypat))=='Sr') then
            atoms%qtypat(itypat)=2.d0
        elseif(trim(atoms%stypat(itypat))=='Al') then
            atoms%qtypat(itypat)=3.d0
        elseif(trim(atoms%stypat(itypat))=='He') then
            atoms%qtypat(itypat)=0.d0
        elseif(trim(atoms%stypat(itypat))=='Mg') then
            atoms%qtypat(itypat)=1.2d0
        else
            write(*,*) 'ERROR: no atomic charge stored for atoms%stypat=',trim(atoms%stypat(itypat))
            stop
        endif
    enddo
    if(.not. allocated(atoms%sat)) stop 'ERROR: in set_qat: sat not allocated.'
    do iat=1,atoms%nat
        if(trim(atoms%sat(iat))=='Na') then
            atoms%qat(iat)=1.d0
        else if(trim(atoms%sat(iat))=='Cl') then
            atoms%qat(iat)=-1.d0
        else if(trim(atoms%sat(iat))=='Li') then
            atoms%qat(iat)=1.d0
        else if(trim(atoms%sat(iat))=='Zn') then
            atoms%qat(iat)=1.d0
        else if(trim(atoms%sat(iat))=='Ti') then
            atoms%qat(iat)=2.0d0
        else if(trim(atoms%sat(iat))=='Zr') then
            atoms%qat(iat)=4.0d0
        else if(trim(atoms%sat(iat))=='Y') then
            atoms%qat(iat)=3.0d0
        else if(trim(atoms%sat(iat))=='Sn') then
            atoms%qat(iat)=2.0d0
        else if(trim(atoms%sat(iat))=='Si') then
            atoms%qat(iat)=1.1d0
        else if(trim(atoms%sat(iat))=='O') then
            atoms%qat(iat)=-1.2d0
        else if(trim(atoms%sat(iat))=='Pb') then
            atoms%qat(iat)=0.9d0
        else if(trim(atoms%sat(iat))=='Te') then
            atoms%qat(iat)=-0.9d0
        else if(trim(atoms%sat(iat))=='Sn') then
            atoms%qat(iat)=0.9d0
        else if(trim(atoms%sat(iat))=='Se') then
            atoms%qat(iat)=-0.9d0
        else if(trim(atoms%sat(iat))=='Ca') then
            atoms%qat(iat)=1.3d0
        else if(trim(atoms%sat(iat))=='F') then
            atoms%qat(iat)=-0.65d0
        else if(trim(atoms%sat(iat))=='Al') then
            atoms%qat(iat)=3.0d0
        elseif(trim(atoms%sat(iat))=='W') then
            atoms%qat(itypat)=0.8d0
        elseif(trim(atoms%sat(iat))=='S') then
            atoms%qat(itypat)=-0.4d0
        else if(trim(atoms%sat(iat))=='K') then
            atoms%qat(iat)=1.d0
        else if(trim(atoms%sat(iat))=='Br') then
            atoms%qat(iat)=-1.d0
        else if(trim(atoms%sat(iat))=='Sr') then
            atoms%qat(iat)=2.d0
        else if(trim(atoms%sat(iat))=='Al') then
            atoms%qat(iat)=3.d0
        else if(trim(atoms%sat(iat))=='He') then
            atoms%qat(iat)=0.d0
        else if(trim(atoms%sat(iat))=='Mg') then
            atoms%qat(iat)=1.2d0
        else
            write(*,*) 'ERROR: no atomic charge stored for atoms%sat=',trim(atoms%sat(iat))
            stop
        endif
        !if(iproc==0) then
        !    write(*,'(a,a5,f10.5)') 'charge of atom:',trim(atoms%sat(iat)),atoms%qat(iat)
        !endif
    enddo
end subroutine set_qat
!*****************************************************************************************
subroutine set_atomic_mass(atoms)
    use mod_processors, only: iproc
    use mod_const, only: bohr2ang
    implicit none
    type(typ_atoms), intent(inout) :: atoms
    !local variables
    integer:: iat
    if(.not. allocated(atoms%sat)) stop 'ERROR: in set_amass:sat not allocated.'
    if(.not. allocated(atoms%amass)) stop 'ERROR: in set_amass: amass not allocated.'
    do iat=1,atoms%nat
        if(trim(atoms%sat(iat))=='LJ') then
            atoms%amass(iat)=1.d0
        else if(trim(atoms%sat(iat))=='H') then
            atoms%amass(iat)=1.00797d0
        else if(trim(atoms%sat(iat))=='He') then
            atoms%amass(iat)=4.00260d0
        else if(trim(atoms%sat(iat))=='Li') then
            atoms%amass(iat)=6.941d0
        else if(trim(atoms%sat(iat))=='Be') then
            atoms%amass(iat)=9.012180
        else if(trim(atoms%sat(iat))=='B') then
            atoms%amass(iat)=10.81d0
        else if(trim(atoms%sat(iat))=='C') then
            atoms%amass(iat)=12.011d0
        else if(trim(atoms%sat(iat))=='N') then
            atoms%amass(iat)=14.0067d0
        else if(trim(atoms%sat(iat))=='O') then
            atoms%amass(iat)=15.9994d0
        else if(trim(atoms%sat(iat))=='F') then
            atoms%amass(iat)=18.998403d0
        else if(trim(atoms%sat(iat))=='Ne') then
            atoms%amass(iat)=20.179d0
        else if(trim(atoms%sat(iat))=='Na') then
            atoms%amass(iat)=22.98977d0
        else if(trim(atoms%sat(iat))=='Mg') then
            atoms%amass(iat)=24.305d0
        else if(trim(atoms%sat(iat))=='Al') then
            atoms%amass(iat)=26.98154d0
        else if(trim(atoms%sat(iat))=='Si') then
            atoms%amass(iat)=28.0855d0
        else if(trim(atoms%sat(iat))=='P') then
            atoms%amass(iat)=30.97376d0
        else if(trim(atoms%sat(iat))=='S') then
            atoms%amass(iat)=32.06d0
        else if(trim(atoms%sat(iat))=='Cl') then
            atoms%amass(iat)=35.453d0
        else if(trim(atoms%sat(iat))=='Ar') then
            atoms%amass(iat)=39.948d0
        else if(trim(atoms%sat(iat))=='K') then
            atoms%amass(iat)=39.0983d0
        else if(trim(atoms%sat(iat))=='Ca') then
            atoms%amass(iat)=40.08d0
        else if(trim(atoms%sat(iat))=='Sc') then
            atoms%amass(iat)=44.9559d0
        else if(trim(atoms%sat(iat))=='Ti') then
            atoms%amass(iat)=47.90d0
        else if(trim(atoms%sat(iat))=='V') then
            atoms%amass(iat)=50.9415d0
        else if(trim(atoms%sat(iat))=='Cr') then
            atoms%amass(iat)=51.996d0
        else if(trim(atoms%sat(iat))=='Mn') then
            atoms%amass(iat)=54.9380d0
        else if(trim(atoms%sat(iat))=='Fe') then
            atoms%amass(iat)=55.847d0
        else if(trim(atoms%sat(iat))=='Co') then
            atoms%amass(iat)=58.9332d0
        else if(trim(atoms%sat(iat))=='Ni') then
            atoms%amass(iat)=58.70d0
        else if(trim(atoms%sat(iat))=='Cu') then
            atoms%amass(iat)=63.546d0
        else if(trim(atoms%sat(iat))=='Zn') then
            atoms%amass(iat)=65.38d0
        else if(trim(atoms%sat(iat))=='Ga') then
            atoms%amass(iat)=69.72d0
        else if(trim(atoms%sat(iat))=='Ge') then
            atoms%amass(iat)=72.59d0
        else if(trim(atoms%sat(iat))=='As') then
            atoms%amass(iat)=74.9216d0
        else if(trim(atoms%sat(iat))=='Se') then
            atoms%amass(iat)=78.96d0
        else if(trim(atoms%sat(iat))=='Br') then
            atoms%amass(iat)=79.904d0
        else if(trim(atoms%sat(iat))=='Kr') then
            atoms%amass(iat)=83.80d0
        else if(trim(atoms%sat(iat))=='Rb') then
            atoms%amass(iat)=85.4678d0
        else if(trim(atoms%sat(iat))=='Sr') then
            atoms%amass(iat)=87.62d0
        else if(trim(atoms%sat(iat))=='Y') then
            atoms%amass(iat)=88.9059d0
        else if(trim(atoms%sat(iat))=='Zr') then
            atoms%amass(iat)=91.22d0
        else if(trim(atoms%sat(iat))=='Nb') then
            atoms%amass(iat)=92.9064d0
        else if(trim(atoms%sat(iat))=='Mo') then
            atoms%amass(iat)=95.94d0
        else if(trim(atoms%sat(iat))=='Tc') then
            atoms%amass(iat)=98d0
        else if(trim(atoms%sat(iat))=='Ru') then
            atoms%amass(iat)=101.07d0
        else if(trim(atoms%sat(iat))=='Rh') then
            atoms%amass(iat)=102.9055d0
        else if(trim(atoms%sat(iat))=='Pd') then
            atoms%amass(iat)=106.4d0
        else if(trim(atoms%sat(iat))=='Ag') then
            atoms%amass(iat)=107.868d0
        else if(trim(atoms%sat(iat))=='Cd') then
            atoms%amass(iat)=112.41d0
        else if(trim(atoms%sat(iat))=='In') then
            atoms%amass(iat)=114.82d0
        else if(trim(atoms%sat(iat))=='Sn') then
            atoms%amass(iat)=118.69d0
        else if(trim(atoms%sat(iat))=='Sb') then
            atoms%amass(iat)=121.75d0
        else if(trim(atoms%sat(iat))=='Te') then
            atoms%amass(iat)=127.60d0
        else if(trim(atoms%sat(iat))=='I') then
            atoms%amass(iat)=126.9045d0
        else if(trim(atoms%sat(iat))=='Xe') then
            atoms%amass(iat)=131.30d0
        else if(trim(atoms%sat(iat))=='Cs') then
            atoms%amass(iat)=132.9054d0
        else if(trim(atoms%sat(iat))=='Ba') then
            atoms%amass(iat)=137.33d0
        else if(trim(atoms%sat(iat))=='La') then
            atoms%amass(iat)=138.9055d0
        else if(trim(atoms%sat(iat))=='Ce') then
            atoms%amass(iat)=140.12d0
        else if(trim(atoms%sat(iat))=='Pr') then
            atoms%amass(iat)=140.9077d0
        else if(trim(atoms%sat(iat))=='Nd') then
            atoms%amass(iat)=144.24d0
        else if(trim(atoms%sat(iat))=='Pm') then
            atoms%amass(iat)=145d0
        else if(trim(atoms%sat(iat))=='Sm') then
            atoms%amass(iat)=150.4d0
        else if(trim(atoms%sat(iat))=='Eu') then
            atoms%amass(iat)=151.96d0
        else if(trim(atoms%sat(iat))=='Gd') then
            atoms%amass(iat)=157.25d0
        else if(trim(atoms%sat(iat))=='Tb') then
            atoms%amass(iat)=158.9254d0
        else if(trim(atoms%sat(iat))=='Dy') then
            atoms%amass(iat)=162.50d0
        else if(trim(atoms%sat(iat))=='Ho') then
            atoms%amass(iat)=164.9304d0
        else if(trim(atoms%sat(iat))=='Er') then
            atoms%amass(iat)=167.26d0
        else if(trim(atoms%sat(iat))=='Tm') then
            atoms%amass(iat)=168.9342d0
        else if(trim(atoms%sat(iat))=='Yb') then
            atoms%amass(iat)=173.04d0
        else if(trim(atoms%sat(iat))=='Lu') then
            atoms%amass(iat)=174.967d0
        else if(trim(atoms%sat(iat))=='Hf') then
            atoms%amass(iat)=178.49d0
        else if(trim(atoms%sat(iat))=='Ta') then
            atoms%amass(iat)=180.9479d0
        else if(trim(atoms%sat(iat))=='W') then
            atoms%amass(iat)=183.85d0
        else if(trim(atoms%sat(iat))=='Re') then
            atoms%amass(iat)=186.207d0
        else if(trim(atoms%sat(iat))=='Os') then
            atoms%amass(iat)=190.2d0
        else if(trim(atoms%sat(iat))=='Ir') then
            atoms%amass(iat)=192.22d0
        else if(trim(atoms%sat(iat))=='Pt') then
            atoms%amass(iat)=195.09d0
        else if(trim(atoms%sat(iat))=='Au') then
            atoms%amass(iat)=196.9665d0
        else if(trim(atoms%sat(iat))=='Hg') then
            atoms%amass(iat)=200.59d0
        else if(trim(atoms%sat(iat))=='Tl') then
            atoms%amass(iat)=204.37d0
        else if(trim(atoms%sat(iat))=='Pb') then
            atoms%amass(iat)=207.2d0
        else if(trim(atoms%sat(iat))=='Bi') then
            atoms%amass(iat)=208.9804d0
        else if(trim(atoms%sat(iat))=='Po') then
            atoms%amass(iat)=209d0
        else if(trim(atoms%sat(iat))=='At') then
            atoms%amass(iat)=210d0
        else if(trim(atoms%sat(iat))=='Rn') then
            atoms%amass(iat)=222d0
        else if(trim(atoms%sat(iat))=='Fr') then
            atoms%amass(iat)=223d0
        else if(trim(atoms%sat(iat))=='Ra') then
            atoms%amass(iat)=226.0254d0
        else if(trim(atoms%sat(iat))=='Ac') then
            atoms%amass(iat)=227.0278d0
        else if(trim(atoms%sat(iat))=='Th') then
            atoms%amass(iat)=232.0381d0
        else if(trim(atoms%sat(iat))=='Pa') then
            atoms%amass(iat)=231.0359d0
        else if(trim(atoms%sat(iat))=='U') then
            atoms%amass(iat)=238.029d0
        else if(trim(atoms%sat(iat))=='Np') then
            atoms%amass(iat)=237.0482d0
        else if(trim(atoms%sat(iat))=='Pu') then
            atoms%amass(iat)=242d0
        else if(trim(atoms%sat(iat))=='Am') then
            atoms%amass(iat)=243d0
        else if(trim(atoms%sat(iat))=='Cm') then
            atoms%amass(iat)=247d0
        else
            write(*,*) 'ERROR: no atomic mass stored for atomtype=',trim(atoms%sat(iat))
            stop
        endif
        !if(iproc==0) then
        !    write(*,'(a,a5,f10.5)') 'RCOV:',trim(atoms%sat(iat)),atoms%amass(iat)
        !endif
        atoms%amass(iat)=atoms%amass(iat)/0.00054858d0
    enddo
end subroutine set_atomic_mass
!*****************************************************************************************
subroutine sat_to_iatom(sat,iatom)
    use mod_processors, only: iproc
    implicit none
    character(*), intent(in) :: sat
    integer, intent(out) :: iatom
    !local variables
    if(trim(sat)=='H') then
        iatom=1
    else if(trim(sat)=='He') then
        iatom=2
    else if(trim(sat)=='Li') then
        iatom=3
    else if(trim(sat)=='Be') then
        iatom=4
    else if(trim(sat)=='B') then
        iatom=5
    else if(trim(sat)=='C') then
        iatom=6
    else if(trim(sat)=='N') then
        iatom=7
    else if(trim(sat)=='O') then
        iatom=8
    else if(trim(sat)=='F') then
        iatom=9
    else if(trim(sat)=='Ne') then
        iatom=10
    else if(trim(sat)=='Na') then
        iatom=11
    else if(trim(sat)=='Mg') then
        iatom=12
    else if(trim(sat)=='Al') then
        iatom=13
    else if(trim(sat)=='Si') then
        iatom=14
    else if(trim(sat)=='P') then
        iatom=15
    else if(trim(sat)=='S') then
        iatom=16
    else if(trim(sat)=='Cl') then
        iatom=17
    else if(trim(sat)=='Ar') then
        iatom=18
    else if(trim(sat)=='K') then
        iatom=19
    else if(trim(sat)=='Ca') then
        iatom=20
    else if(trim(sat)=='Sc') then
        iatom=21
    else if(trim(sat)=='Ti') then
        iatom=22
    else if(trim(sat)=='V') then
        iatom=23
    else if(trim(sat)=='Cr') then
        iatom=24
    else if(trim(sat)=='Mn') then
        iatom=25
    else if(trim(sat)=='Fe') then
        iatom=26
    else if(trim(sat)=='Co') then
        iatom=27
    else if(trim(sat)=='Ni') then
        iatom=28
    else if(trim(sat)=='Cu') then
        iatom=29
    else if(trim(sat)=='Zn') then
        iatom=30
    else if(trim(sat)=='Ga') then
        iatom=31
    else if(trim(sat)=='Ge') then
        iatom=32
    else if(trim(sat)=='As') then
        iatom=33
    else if(trim(sat)=='Se') then
        iatom=34
    else if(trim(sat)=='Br') then
        iatom=35
    else if(trim(sat)=='Kr') then
        iatom=36
    else if(trim(sat)=='Rb') then
        iatom=37
    else if(trim(sat)=='Sr') then
        iatom=38
    else if(trim(sat)=='Y') then
        iatom=39
    else if(trim(sat)=='Zr') then
        iatom=40
    else if(trim(sat)=='Nb') then
        iatom=41
    else if(trim(sat)=='Mo') then
        iatom=42
    else if(trim(sat)=='Tc') then
        iatom=43
    else if(trim(sat)=='Ru') then
        iatom=44
    else if(trim(sat)=='Rh') then
        iatom=45
    else if(trim(sat)=='Pd') then
        iatom=46
    else if(trim(sat)=='Ag') then
        iatom=47
    else if(trim(sat)=='Cd') then
        iatom=48
    else if(trim(sat)=='In') then
        iatom=49
    else if(trim(sat)=='Sn') then
        iatom=50
    else if(trim(sat)=='Sb') then
        iatom=51
    else if(trim(sat)=='Te') then
        iatom=52
    else if(trim(sat)=='I') then
        iatom=53
    else if(trim(sat)=='Xe') then
        iatom=54
    else if(trim(sat)=='Cs') then
        iatom=55
    else if(trim(sat)=='Ba') then
        iatom=56
    else if(trim(sat)=='La') then
        iatom=57
    else if(trim(sat)=='Ce') then
        iatom=58
    else if(trim(sat)=='Pr') then
        iatom=59
    else if(trim(sat)=='Nd') then
        iatom=60
    else if(trim(sat)=='Pm') then
        iatom=61
    else if(trim(sat)=='Sm') then
        iatom=62
    else if(trim(sat)=='Eu') then
        iatom=63
    else if(trim(sat)=='Gd') then
        iatom=64
    else if(trim(sat)=='Tb') then
        iatom=65
    else if(trim(sat)=='Dy') then
        iatom=66
    else if(trim(sat)=='Ho') then
        iatom=67
    else if(trim(sat)=='Er') then
        iatom=68
    else if(trim(sat)=='Tm') then
        iatom=69
    else if(trim(sat)=='Yb') then
        iatom=70
    else if(trim(sat)=='Lu') then
        iatom=71
    else if(trim(sat)=='Hf') then
        iatom=72
    else if(trim(sat)=='Ta') then
        iatom=73
    else if(trim(sat)=='W') then
        iatom=74
    else if(trim(sat)=='Re') then
        iatom=75
    else if(trim(sat)=='Os') then
        iatom=76
    else if(trim(sat)=='Ir') then
        iatom=77
    else if(trim(sat)=='Pt') then
        iatom=78
    else if(trim(sat)=='Au') then
        iatom=79
    else if(trim(sat)=='Hg') then
        iatom=80
    else if(trim(sat)=='Tl') then
        iatom=81
    else if(trim(sat)=='Pb') then
        iatom=82
    else if(trim(sat)=='Bi') then
        iatom=83
    else if(trim(sat)=='Po') then
        iatom=84
    else if(trim(sat)=='At') then
        iatom=85
    else if(trim(sat)=='Rn') then
        iatom=86
    else if(trim(sat)=='Fr') then
        iatom=87
    else if(trim(sat)=='Ra') then
        iatom=88
    else if(trim(sat)=='Ac') then
        iatom=89
    else if(trim(sat)=='Th') then
        iatom=90
    else if(trim(sat)=='Pa') then
        iatom=91
    else if(trim(sat)=='U') then
        iatom=92
    else if(trim(sat)=='Np') then
        iatom=93
    else if(trim(sat)=='Pu') then
        iatom=94
    else if(trim(sat)=='Am') then
        iatom=95
    else if(trim(sat)=='Cm') then
        iatom=96
    else if(trim(sat)=='Bk') then
        iatom=97
    else if(trim(sat)=='Cf') then
        iatom=98
    else if(trim(sat)=='Es') then
        iatom=99
    else if(trim(sat)=='Fm') then
        iatom=100
    else if(trim(sat)=='Md') then
        iatom=101
    else if(trim(sat)=='No') then
        iatom=102
    else if(trim(sat)=='Lr') then
        iatom=103
    else if(trim(sat)=='Rf') then
        iatom=104
    else if(trim(sat)=='Db') then
        iatom=105
    else if(trim(sat)=='Sg') then
        iatom=106
    else if(trim(sat)=='Bh') then
        iatom=107
    else if(trim(sat)=='Hs') then
        iatom=108
    else if(trim(sat)=='Mt') then
        iatom=109
    else if(trim(sat)=='Ds') then
        iatom=110
    else if(trim(sat)=='LJ') then
        iatom=201
    else
        write(*,*) 'ERROR: no atomic number stored for atomtype=',trim(sat)
        stop
    endif
end subroutine sat_to_iatom
!*****************************************************************************************
subroutine iatom_to_sat(iatom,sat)
    use mod_processors, only: iproc
    implicit none
    integer, intent(in) :: iatom
    character(*), intent(out) :: sat
    !local variables
    character(5):: elements(1:109)=(/ &
      ' H','He','Li','Be',' B',' C',' N',' O',' F','Ne','Na',       &
      'Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca','Sc','Ti',    &
      ' V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As',  &
      'Se','Br','Kr','Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru',  &
      'Rh','Pd','Ag','Cd','In','Sn','Sb','Te',' I','Xe','Cs',  &
      'Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy', &
      'Ho','Er','Tm','Yb','Lu','Hf','Ta',' W','Re','Os','Ir',  &
      'Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra', &
      'Ac','Th','Pa',' U','Np','Pu','Am','Cm','Bk','Cf','Es',  &
      'Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt'      /)

    if(iatom<1 .or. iatom>109) then
        write(*,*) 'ERROR: no symbol stored for atomic number=',iatom
        stop
    else
        sat=adjustl(elements(iatom))
    endif
end subroutine iatom_to_sat
!*****************************************************************************************
end module mod_atoms
!*****************************************************************************************
