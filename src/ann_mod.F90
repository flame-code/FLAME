!*****************************************************************************************
module mod_ann
    use dictionaries
    use mod_linked_lists, only: typ_linked_lists
    use mod_electrostatics, only: typ_poisson
    implicit none
    private
    public:: ann_allocate, ann_deallocate, set_number_of_ann
    type, public:: typ_ann
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
        real(8):: a(350,350,10), b(350,10), x(350,10), y(350,0:10), yd(350,10), ad(350*350,10), bd(350,10)
        real(8):: d(350)
        real(8):: gbounds(2,350)
        real(8):: two_over_gdiff(350)
        real(8):: ebounds(2)=(/-1.d0,1.d0/)
        !real(8):: rc1(350)
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
        real(8):: g1eta(350)
        real(8):: g1rs(350)

        !The 2nd type of symmetry functions introduced by Behler
        real(8):: g2eta(350)
        real(8):: g2rs(350)
        integer:: g2i(350)

        !The 3rd type of symmetry functions introduced by Behler
        real(8):: g3kappa(350)

        !The 4th type of symmetry functions introduced by Behler
        real(8):: g4eta(350)
        real(8):: g4zeta(350)
        real(8):: g4lambda(350)

        !The 5th type of symmetry functions introduced by Behler
        real(8):: g5eta(350)
        real(8):: g5zeta(350)
        real(8):: g5lambda(350)
        integer:: g5i(2,350)

        !The 6th type of symmetry functions introduced by Ghasemi
        real(8):: g6eta(350)
        real(8):: teneria(3,3,30)

        !some other variables
        real(8):: his(1000,350)

        character(256):: hlines(10)
        
    end type typ_ann
    type, public:: typ_ann_arr
        logical:: exists_yaml_file = .false.
        integer:: iunit
        integer:: nann=-1
        integer:: natmax=1000
        logical:: compute_symfunc=.true.
        logical:: cal_force=.true.
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
    type, public:: typ_cent
        real(8), allocatable:: gwi(:)
        real(8), allocatable:: gwe(:)
        real(8), allocatable:: gwit(:)
        real(8), allocatable:: rel(:,:)
        real(8), allocatable:: qgrad(:)
        real(8), allocatable:: rgrad(:,:)
        type(typ_poisson):: poisson
    end type typ_cent
contains
!*****************************************************************************************
subroutine set_number_of_ann(parini,ann_arr)
    use mod_parini, only: typ_parini
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    !local variables
    ann_arr%nann=parini%ntypat
    if(parini%bondbased_ann) then
        ann_arr%nann=4
    endif
end subroutine set_number_of_ann
!*****************************************************************************************
subroutine ann_allocate(n,num,ann_arr)
    !use mod_ann, only: typ_ann_arr
    !use mod_opt_ann, only: typ_opt_ann
    use dynamic_memory
    implicit none
    integer, intent(in):: n
    integer, intent(in):: num(n)
    type(typ_ann_arr), intent(inout):: ann_arr
    !local variables
    integer:: istat, ng
    ng=ann_arr%ann(1)%nn(0)
    !allocate(ann_arr%yall(ng,natmax),stat=istat)
    !if(istat/=0) stop 'ERROR: unable to allocate array ann_arr%yall.'
    !allocate(ann_arr%y0d(ng,3,natmax,natmax),stat=istat)
    !if(istat/=0) stop 'ERROR: unable to allocate array ann_arr%y0d.'
    !allocate(ann_arr%y0dr(ng,9,natmax,natmax),stat=istat)
    !if(istat/=0) stop 'ERROR: unable to allocate array ann_arr%y0dr.'
    allocate(ann_arr%fat_chi(1:3,1:ann_arr%natmax))
    allocate(ann_arr%chi_i(1:ann_arr%natmax))
    allocate(ann_arr%chi_o(1:ann_arr%natmax))
    allocate(ann_arr%chi_d(1:ann_arr%natmax))
    allocate(ann_arr%a(1:(ann_arr%natmax+1)*(ann_arr%natmax+1)))
    ann_arr%fat_chi=0.d0
    ann_arr%chi_i=0.d0
    ann_arr%chi_o=0.d0
    ann_arr%chi_d=0.d0
    ann_arr%a=0.d0
    allocate(ann_arr%g_per_atom(1:num(1),1:ann_arr%natmax))
    !symfunc%linked_lists%maxbound_rad is assumed 10000
    allocate(ann_arr%fatpq(1:3,1:10000))
    allocate(ann_arr%stresspq(1:3,1:3,1:10000))
    allocate(ann_arr%ipiv(1:ann_arr%natmax+1))
    allocate(ann_arr%qq(1:ann_arr%natmax+1))
end subroutine ann_allocate
!*****************************************************************************************
subroutine ann_deallocate(ann_arr)
    !use mod_ann, only: typ_ann_arr
    use dynamic_memory
    implicit none
    type(typ_ann_arr), intent(inout):: ann_arr
    !local variables
    integer:: istat
    !deallocate(ann_arr%yall,stat=istat)
    !if(istat/=0) stop 'ERROR: unable to deallocate array ann_arr%yall.'
    !deallocate(ann_arr%y0d,stat=istat)
    !if(istat/=0) stop 'ERROR: unable to deallocate array ann_arr%y0d.'
    !deallocate(ann_arr%y0dr,stat=istat)
    !if(istat/=0) stop 'ERROR: unable to deallocate array ann_arr%y0dr.'
    deallocate(ann_arr%chi_i)
    deallocate(ann_arr%chi_o)
    deallocate(ann_arr%chi_d)
    deallocate(ann_arr%a)
    deallocate(ann_arr%fat_chi)
    deallocate(ann_arr%g_per_atom)
    deallocate(ann_arr%fatpq)
    deallocate(ann_arr%stresspq)
    deallocate(ann_arr%ipiv)
    deallocate(ann_arr%qq)
end subroutine ann_deallocate
!*****************************************************************************************
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
