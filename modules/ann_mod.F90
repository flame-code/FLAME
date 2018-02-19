!*****************************************************************************************
module mod_ann
    use dictionaries
    use mod_linked_lists, only: typ_linked_lists
    use mod_electrostatics, only: typ_ewald_p3d
    implicit none
    type typ_ann
        type(dictionary), pointer :: dict
        integer:: nl !number of hidden layer plus one
        integer:: nn(0:10)
        !integer:: n0=-1
        !integer:: n1=-1
        !integer:: n2=-1
        !integer:: n_all=n0*n1+n1*n2+n2+n1+n2+1
        integer:: ng1=-1
        integer:: ng2=-1
        integer:: ng3=-1
        integer:: ng4=-1
        integer:: ng5=-1
        integer:: ng6=-1
        real(8):: a(140,140,10), b(140,10), x(140,10), y(140,0:10), yd(140,10), ad(140*140,10), bd(140,10)
        real(8):: d(140)
        real(8):: gbounds(2,140)
        real(8):: two_over_gdiff(140)
        real(8):: ebounds(2)=(/-1.d0,1.d0/)
        !real(8):: rc1(140)
        real(8):: gausswidth
        real(8):: gausswidth_ion
        real(8):: chi0
        real(8):: hardness
        real(8):: spring_const
        real(8):: zion
        real(8):: qinit
        real(8):: rionic
        real(8):: ener_ref
        real(8):: ampl_chi=-1.d0
        real(8):: prefactor_chi=-1.d0
        character(20):: method

        !The 1st type of symmetry functions introduced by Behler
        real(8):: g1eta(140)
        real(8):: g1rs(140)

        !The 2nd type of symmetry functions introduced by Behler
        real(8):: g2eta(140)
        real(8):: g2rs(140)
        integer:: g2i(140)

        !The 3rd type of symmetry functions introduced by Behler
        real(8):: g3kappa(140)

        !The 4th type of symmetry functions introduced by Behler
        real(8):: g4eta(140)
        real(8):: g4zeta(140)
        real(8):: g4lambda(140)

        !The 5th type of symmetry functions introduced by Behler
        real(8):: g5eta(140)
        real(8):: g5zeta(140)
        real(8):: g5lambda(140)
        integer:: g5i(2,140)

        !The 6th type of symmetry functions introduced by Ghasemi
        real(8):: g6eta(140)
        real(8):: teneria(3,3,30)

        !some other variables
        real(8):: his(1000,140)

        character(256):: hlines(10)
        
    end type typ_ann
    type typ_ann_arr
        logical:: exists_yaml_file = .false.
        integer:: n=-1
        integer:: natmax=1000
        logical:: compute_symfunc=.true.
        character(30):: event='unknown'
        character(50):: approach='unknown'
        real(8):: rcut=-1.d0
        real(8):: ener_ref
        real(8):: epot_es
        real(8):: fchi_angle
        real(8):: fchi_norm
        !real(8), allocatable:: yall(:,:)
        !real(8), allocatable:: y0d(:,:,:)
        integer:: natsum(10)
        !real(8):: repfac(10,10)
        real(8):: reprcut(10,10)
        real(8):: qmax(10)
        real(8):: qmin(10)
        real(8):: qsum(10)
        real(8):: chi_max(10)
        real(8):: chi_min(10)
        real(8):: chi_sum(10)
        real(8):: chi_delta(10)
        real(8):: yall_bond(100,100,100)
        real(8):: y0d_bond(100,3,100,100)
        !real(8), allocatable:: y0dr(:,:,:)
        real(8), allocatable:: a(:)
        real(8), allocatable:: chi_i(:)
        real(8), allocatable:: chi_o(:)
        real(8), allocatable:: chi_d(:)
        real(8), allocatable:: fat_chi(:,:)
        real(8), allocatable:: g_per_atom(:,:)
        real(8), allocatable:: g_per_bond(:,:,:)
        real(8), allocatable:: fatpq(:,:)
        real(8), allocatable:: stresspq(:,:,:)
        integer, allocatable:: ipiv(:)
        real(8), allocatable:: qq(:)
        type(typ_ann), allocatable:: ann(:)
    end type typ_ann_arr
    type typ_symfunc
        integer:: ng=-1
        integer:: nat=-1
        real(8):: epot
        real(8), allocatable:: y(:,:)
        real(8), allocatable:: y0d_bond(:,:)
        real(8), allocatable:: y0d(:,:,:)
        real(8), allocatable:: y0dr(:,:,:)
        type(typ_linked_lists):: linked_lists
    end type typ_symfunc
    type typ_symfunc_arr
        integer:: nconf=-1
        type(typ_symfunc), allocatable:: symfunc(:)
    end type typ_symfunc_arr
    type typ_ekf
        integer:: n=-1
        integer:: loc(10)
        integer:: num(10)
        real(8), allocatable:: x(:)
        real(8), allocatable:: epotd(:)
        real(8), allocatable:: g(:) !gradient of neural artificial neural network output
        real(8), allocatable:: gc(:,:)
        real(8), allocatable:: gs(:,:)
    end type typ_ekf
    type typ_cent
        real(8), allocatable:: gwi(:)
        real(8), allocatable:: gwe(:)
        real(8), allocatable:: gwit(:)
        real(8), allocatable:: rel(:,:)
        real(8), allocatable:: qgrad(:)
        real(8), allocatable:: rgrad(:,:)
        type(typ_ewald_p3d):: ewald_p3d
    end type typ_cent
end module mod_ann
!*****************************************************************************************
!module data_point
!    implicit none
!    integer:: m
!    integer:: natmax=-1
!    integer, allocatable:: natarr(:)
!    real(8), allocatable:: pnt(:), f_p(:)
!    real(8), allocatable:: ratall(:,:,:), epotall(:)
!end module data_point
!*****************************************************************************************
module mod_parlm
    implicit none
    type typ_parlm
        real(8):: ftol=1.d-8
        real(8):: xtol=1.d-8
        real(8):: gtol=1.d-8
        real(8):: factor=100.d0
        integer:: maxfev=1000
        integer:: nprint=1
        integer:: n=0
        integer:: mode, info, nfev, njev
        integer:: iter
        integer:: icontinue
        integer:: iflag
        logical:: finish
        real(8):: epsmch
        real(8):: fnorm
        !real(8):: fnorm1
        real(8):: xnorm
        real(8):: gnorm
        real(8):: pnorm
        real(8):: par
        real(8):: delta
        real(8):: actred
        real(8):: prered
        real(8):: ratio
        real(8), allocatable:: wa1(:), wa2(:), wa3(:), wa4(:), qtf(:)
        real(8), allocatable:: x(:)
        real(8), allocatable:: fvec(:)
        real(8), allocatable:: fjac(:,:)
        real(8), allocatable:: diag(:)
        integer, allocatable:: ipvt(:)
    end type typ_parlm
end module mod_parlm
!*****************************************************************************************
