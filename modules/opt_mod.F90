!*****************************************************************************************
module mod_opt
    character(36), parameter:: frmt_base="a4,i3.3,1x,i5,es20.12,es11.3,2es12.5"
    type typ_paropt
        !parameters for some or all methods
        character(10):: approach='unknown'
        character(10):: approach_current='unknown'
        character(20):: precaution='normal'
        logical:: print_force=.false.
        real(8):: strfact
        !character(36):: frmt_base !="a4,i3.3,1x,i5,es23.15,es11.3,2es12.5"
        !parameter (frmt_base="a4,i3.3,1x,i5,es23.15,es11.3,2es12.5")
        integer:: iter=0
        integer:: iflag=0
        integer:: maxforcecall=-1
        integer:: ifail=0
        integer:: nfail=-1
        integer:: nsatur=-1
        integer:: isatur
        integer:: feedback=-1
        integer:: nit=-1
        logical:: converged
        logical:: lprint=.false.
        logical:: param_reported=.false.
        logical:: cellrelax=.false.
        real(8):: dxmax=-1.d0
        real(8):: fmax
        real(8):: fmaxtol=-1.d0
        real(8):: eps=1.d-5
        real(8):: funits=-1.d0
        real(8):: alphax=-1.d0
        real(8):: alpha=-1.d0
        real(8):: condnum=10.d0
        real(8):: anoise=-1.d0
        real(8):: epotold
        real(8):: epotitm1
        real(8):: epotitm2
        real(8):: fnrmitm1
        real(8):: fnrmitm2
        logical:: trajectory=.false.
        character(256):: filename='unknown'
        !parameters for SD
        integer:: itsd
        !integer:: nitsd=-1
        logical:: sdsaturated
        logical:: care
        logical:: xmoved
        logical:: optional_control_on_saturation
        real(8):: alphamin=-1.d0
        real(8):: alphamax=-1.d0
        real(8):: fnrmtolsatur=-1.d0
        character(10):: sdsaturation
        !parameters for CG
        integer:: itcg=-1
        !integer:: nitcg=-1
        logical:: dolinesearch
        real(8):: avgalpha=-1.d0
        real(8):: avgnum=-1.d0
        real(8):: alpha0=-1.d0
        !parameters for DIIS
        integer:: idsx
        integer, allocatable::ipiv(:)
        real(8), allocatable::a(:,:,:),b(:)
        !parameters for L-BFGS
        logical:: diagco=.false.
        !parameters for BFGS
        integer:: ifnrminc
        logical:: reset_hess=.false.
        real(8):: evalmin
        real(8):: evalmax
        real(8):: zeta
        real(8):: prefactor
        integer:: increment
        !parameters for FIRE
        integer:: ndown
        integer:: ndowntol=-1
        integer:: itfire=0
        !integer:: nitfire=-1
        real(8):: dt=-1.d0
        real(8):: dt_start=-1.d0
        real(8):: dtmin=-1.d0
        real(8):: dtmax=-1.d0
        real(8):: finc=-1.d0
        real(8):: fdec=-1.d0
        real(8):: falpha=-1.d0
        real(8):: alphastart=-1.d0
        !parameters for SQNM
        integer:: nhist=10
        real(8):: beta_lat
        real(8):: beta_at
        real(8):: maxrise
        real(8):: cutoffRatio
        real(8):: steepthresh
        real(8):: trustr
        !parameters for print information
        integer:: mp=6
        integer:: lp=6
        integer:: iprint(2)=(/1,0/)
        !parameters for line search routine
        integer:: maxfev=20
        real(8):: ftol=1.d-4
        real(8):: gtol=9.d-1
        real(8):: stpmin=1.d-20
        real(8):: stpmax=1.d+20
        real(8):: xtol=epsilon(xtol)
    end type typ_paropt
end module mod_opt
!*****************************************************************************************
