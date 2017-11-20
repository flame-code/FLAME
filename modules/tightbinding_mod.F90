!**************************************************************************************
module mod_tightbinding
    implicit none
    logical:: lenosky=.false.
    type typ_partb
        integer:: norb !number of orbitals (= number of electrons)
        integer:: norbcut
        real(8):: temp_fermi=100.d0 !Fermi temperature in eV, is set by code also
        integer:: usetb
        real(8):: paircut
        integer:: usepairpot
        integer:: nstride=4
        real(8):: pairen
        real(8):: eband !total energy for bands
        real(8):: frac=0.5d0
        integer:: extra=20
        character(30):: event='unknown'
        integer, allocatable:: indorb(:)
        integer, allocatable:: indat(:)
        integer, allocatable:: norbat(:)
        real(8), allocatable:: dedh(:,:,:)
        real(8), allocatable:: hgenall0(:)
        real(8), allocatable:: hgenall1(:)
        real(8), allocatable:: hgenall2(:)
        real(8), allocatable:: hgenall3(:)
        real(8), allocatable:: dhgenall0(:)
        real(8), allocatable:: dhgenall1(:)
        real(8), allocatable:: dhgenall2(:)
        real(8), allocatable:: dhgenall3(:)
        real(8), allocatable:: eval(:)
        real(8), allocatable:: evec(:,:)
        real(8), allocatable:: focc(:)
        real(8), allocatable:: tbmat(:,:)
    end type typ_partb
end module mod_tightbinding
!**************************************************************************************
module mod_splinetb
    implicit none 
        integer, parameter :: NSPMAX=20
    type spline_typ
        real(8):: x(0:NSPMAX)
        real(8):: y(0:NSPMAX)
        real(8):: yp1
        real(8):: ypn
        real(8):: y2(0:NSPMAX)
        integer:: npt
        integer:: iset1,isetn,ider1,idern
    end type spline_typ
end module mod_splinetb
!***************************************************************************************
module mod_potl
    use mod_splinetb, only: NSPMAX, spline_typ
    implicit none
        integer, parameter :: MRHOMAX=1
        integer, parameter :: MHMAX=11
        integer, parameter :: MEPSMAX=10
        integer, parameter :: MRHOMAX2=1
    type potl_typ
        type(spline_typ) :: phi(0:2),f,g,u(0:MRHOMAX-1),rho(0:MRHOMAX-1),h(0:MHMAX-1)
        real(8) :: eps(0:MEPSMAX-1)
        type(spline_typ) :: xf(0:MRHOMAX2-1),xg(0:MRHOMAX2-1),xu(0:MRHOMAX2-1),xrho(0:MRHOMAX2-1)
    end type potl_typ
end module mod_potl
!***************************************************************************************
