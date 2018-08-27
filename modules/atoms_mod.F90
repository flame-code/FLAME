!*****************************************************************************************
module mod_atoms
    implicit none
    type typ_atoms
        integer:: nat=-1 !number of atoms
        integer:: natim=0 !number of atoms of all periodic images including itself 
        integer:: ndof=-1 !number of degrees of freedom
        real(8):: cellvec(3,3)=-1.d0 !cell vectors
        real(8):: celldv(3,3)=0.d0 !
        real(8):: stress(3,3)=0.d0 !
        real(8):: epot=0.d0 !potential energy
        real(8):: ekin=0.d0 !kinetic energy
        real(8):: etot=0.d0 !total energy
        real(8):: enth=0.d0 !enthalpy
        real(8):: ebattery=0.d0 !energy of external work of battery in p3d_bias
        real(8):: qtot=0.d0 !total charge of the configuration
        real(8):: ztot=0.d0 !total ionic charge of the configuration
        real(8):: pressure=0.d0 !external pressure
        real(8):: tol
        real(8):: qtypat(20)=1.d20
        integer:: ntypat=-1
        integer:: ltypat(20)=-1
        integer:: nfp=-1
        character(5):: stypat(20)='none'
        character(20):: boundcond='unknown'
        character(10):: units='angstrom'
        !coordinates type only at time of reading from file and writing to file.
        character(10):: coordinates_type='cartesian'
        character(50):: alloclist='all'
        character(5), allocatable:: sat(:) !symbol of atoms
        real(8), allocatable:: rat(:,:) !atomic positions
        real(8), allocatable:: ratim(:,:) !atomic positions of periodic images
        real(8), allocatable:: vat(:,:) !atomic velocities
        real(8), allocatable:: amass(:) !atomic mass
        real(8), allocatable:: fat(:,:) !atomic forces
        logical, allocatable:: bemoved(:,:) !status to be moved or not
        real(8), allocatable:: qat(:) !atomic charges
        real(8), allocatable:: zat(:) !ionic charges
        real(8), allocatable:: rcov(:) !covalent radii
        real(8), allocatable:: fp(:) !fingerprint
        integer, allocatable:: itypat(:) !The type of each atom is set in this array
        !contains
        !procedure:: atoms_assign
        !generic:: assignment(=) => atoms_assign
    end type typ_atoms
    type typ_atoms_all
        type(typ_atoms):: atoms
        integer:: nconfmax=-1 !maximum number of configurations
        integer:: nconf=-1 !number of configurations
        real(8), allocatable:: epotall(:) !potential energy
        real(8), allocatable:: qtotall(:) !total charge
        real(8), allocatable:: ratall(:,:,:) !atomic positions
        real(8), allocatable:: fatall(:,:,:) !atomic positions
        real(8), allocatable:: fpall(:,:) !fingerprint
    end type typ_atoms_all
    type typ_atoms_arr
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
    type typ_file_info
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
    type type_pairs
        integer ,allocatable:: posat2nd(:,:)
    end type type_pairs
end module mod_atoms
!*****************************************************************************************
