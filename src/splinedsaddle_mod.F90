module minimization_sp
    implicit none
    type parameterminimization_sp
        !general parameters for all methods
        integer::ifile=6
        character(10)::approach='unknown'
        real(8)::fmaxtol=-1.d0
        real(8)::eps=1.d-5
        integer::iter=0
        integer::iflag=0
        logical::converged
        real(8)::epotitm1
        real(8)::avgalpha
        real(8)::avgnum
        integer::iprint=1

        !parameters for SD and SDDIIS
        real(8)::alpha=-1.d0
        real(8)::alphax=-1.d0
        integer::maxforcecall=10000
        real(8)::anoise=-1.d0

        !parameters for SD
        integer::itsd
        integer::nitsd=-1
        integer::nsatur=-1
        integer::isatur
        real(8)::alphamin=-1.d0
        real(8)::alphamax=-1.d0
        real(8)::fnrmtolsatur=-1.d0
        real(8)::epotitm2
        real(8)::fnrmitm1
        real(8)::fnrmitm2
        logical::sdsaturated
        logical::sdminimum
        logical::care

        !parameters for FIRE
        real(8)::dt=-1.d0
        real(8)::dtmax=-1.d0
        real(8)::finc=-1.d0
        real(8)::fdec=-1.d0
        real(8)::falpha=-1.d0
        real(8)::alphastart=-1.d0
        integer::ndowntol=-1
        integer::itfire=0

        !parameters for DIIS
        integer::idsx
        integer, allocatable::ipiv(:)
        real(8), allocatable::a(:,:,:),b(:)
        real(8)::emin
        real(8)::fnrmlowest
        logical::diisminimum
        logical::diisdivergence
        integer::itdiis
        integer::ld
        integer::nd

        !parameters for print information
        integer::mp=6
        integer::lp=6

        !parameters for line search routine
        integer::maxfev=20
        real(8)::ftol=1.d-4
        real(8)::gtol=9.d-1
        real(8)::stpmin=1.d-20
        real(8)::stpmax=1.d+20
        !real(8)::xtol=epsilon(xtol)
        integer::info
    end type parameterminimization_sp
end module minimization_sp


!> Module used by the program splined_saddle
module modulesplinedsaddle
    implicit none
    type parametersplinedsaddle
        !integer, parameter::npmax=20
        !integer::napmax=50
        real(8)::s(0:200)
        real(8)::h(200)
        real(8)::c(0:200)
        real(8)::y(0:200)
        real(8)::e1(200-1)
        real(8)::e2(200-2)
        real(8)::cv(0:50)
        real(8)::ex(0:50)
        real(8)::exd(0:50)
        real(8)::sv(0:50)
        real(8)::e1v(49)
        real(8)::e2v(48)
        real(8)::hv(50)
        real(8)::a(50)
        real(8)::b(50)
        real(8)::tmax
        real(8)::htol=1.d-6  !1.d-6
        real(8)::vdtol=1.d-5  !1.d-10
        real(8)::exends(2)
        real(8)::exends_b(2)
        integer::ifile=6
        integer::ns=2
        integer::ns2=0  !45
        integer::npv
        integer::ncount_ll
        integer::ncount
        real(8)::time_ll
        real(8)::time
        real(8)::epotci
        logical::granot
        character(20)::hybrid
        character(20)::doneb
        character(20)::docineb
        character(20)::pickbestanchorpoints
        logical::do_fill_ex_exd
        character(20)::runstat
        character(10)::typintpol
    end type parametersplinedsaddle
end module modulesplinedsaddle
!*****************************************************************************************
