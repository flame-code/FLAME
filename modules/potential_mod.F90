!*****************************************************************************************
module mod_potential
    use mod_ann, only: typ_ann_arr
    use mod_electrostatics, only: typ_poisson
    use mod_shortrange, only: typ_tosifumi, typ_shortrange
    implicit none
    integer, parameter:: natmax=1000
    !all potentials
    real(8):: fcalls
    real(8):: fcalls_sec
    character(50):: potential='unknown'
    character(50):: potential_sec='unknown'
    !selected potentials
    real(8):: cell(3)
    integer:: ntypat
    integer:: natarr(128)
    character(5):: stypat(128)
    real(8):: cellvec(3,3)
    integer:: l_atomic_num(natmax)
    character(20):: perfstatus='unknown'
    character(20):: perfstatus_old='unknown'
    character(5), allocatable :: sat(:)
    !Lenosky TB potential
    integer:: natsi
    !ANN potential
    type(typ_ann_arr):: ann_arr
    character(50):: ann_boundcheck='none'
    character(50):: approach_ann='unknown'
    !MPMD potential
    !BigDFT potential
    logical:: bigdft_restart=.false.
    !VASP potential
    character(256):: comment, comment2
    logical:: vasp5
    logical:: single_point_calculation=.true.
    logical:: vasp_restart=.false.
    logical:: atom_motion(3,natmax)
    !VCBLJ potential
    real(8):: sigmalj(2,2) !sigma(i,j) is the sigma between particle i and particle j
    real(8):: epslj(2,2) !eps(i,j) is the epsilon between particle i and particle j
    real(8):: rcut(2,2) !rcut(i,j) is the cutoff distance between particle i and j
    !real(8):: atmass(2)
    !real(8):: atmassinv(2)
    integer, allocatable:: kinds(:)
    !coulomb potential
    type(typ_poisson):: poisson
    !force potentials
    type(typ_shortrange):: shortrange
    !socket
    INTEGER:: sock_socket, sock_inet, sock_port        ! socket ID & address of the socket
    CHARACTER(LEN=1024) :: sock_host
    INTEGER, PARAMETER  :: MSGLEN=12   ! length of the headers of the driver/wrapper communication protocol
    CHARACTER(LEN=60)   :: sock_extra_string="                                                            "
    real(8)             :: sock_ecutwf(2)
    logical             :: reset=.true.

end module mod_potential
!*****************************************************************************************
