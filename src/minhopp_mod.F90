!*****************************************************************************************
module mod_minhopp
    implicit none
    integer:: mdmin=3
    !---------------------------------------------
    !original values for parameters of minima hopping feedback
    real(8):: beta1=1.05d0
    real(8):: beta2=1.05d0
    real(8):: beta3=1.d0/1.05d0
    real(8):: alpha1=1.d0/1.02d0
    real(8):: alpha2=1.02d0
    real(8):: alpha3=1.00d0
    !---------------------------------------------
    !LJ good values for parameters of minima hopping feedback
    !real(8), parameter:: beta1=1.05d0
    !real(8), parameter:: beta2=1.02d0 !/1.05d0
    !real(8), parameter:: beta3=1.d0/1.05d0
    !real(8), parameter:: alpha1=1.d0/1.02d0
    !real(8), parameter:: alpha2=1.05d0
    !real(8), parameter:: alpha3=1.02d0
    !---------------------------------------------
    !LJ38: good values for parameters of minima hopping feedback
    !real(8), parameter:: beta1=1.10d0
    !real(8), parameter:: beta2=1.05d0 !/1.05d0
    !real(8), parameter:: beta3=1.d0/1.10d0
    !real(8), parameter:: alpha1=1.d0/1.10d0
    !real(8), parameter:: alpha2=1.10d0
    !real(8), parameter:: alpha3=1.05d0
    !---------------------------------------------
    !aggressive values for parameters of minima hopping feedback
    !real(8), parameter:: beta1=1.06d0
    !real(8), parameter:: beta2=1.04d0
    !real(8), parameter:: beta3=1.d0/1.1d0
    !real(8), parameter:: alpha1=1.d0/1.1d0
    !real(8), parameter:: alpha2=1.1d0
    !---------------------------------------------
    !conservative values for parameters of minima hopping feedback
    !real(8), parameter:: beta1=1.05d0
    !real(8), parameter:: beta2=1.05d0
    !real(8), parameter:: beta3=1.d0/1.1d0
    !real(8), parameter:: alpha1=1.d0/1.05d0
    !real(8), parameter:: alpha2=1.02d0
    !---------------------------------------------
    integer:: nlmin=-1 !number of local minima already found.
    integer:: nlminx=-1 !number of local minima to be found including nlmin. 
    integer:: nlmin_l=0
    integer:: nstep=0
    integer:: nlmin_old
    integer:: minter=1 !interval for writing intermediate results
    integer, parameter:: nbuf=1000 
    integer:: npminx=1000 !maximum number of minima whose configurations will be stored
    real(8):: eref=-1.d50 !reference energy.
    !real(8):: accur
    !real(8):: accur !accuracy for rounding.
    real(8):: ekin
    real(8):: av_ekinetic=0.d0
    real(8):: av_ediff=0.d0
    integer:: istep=0
    integer:: istep_sam=0
    integer:: istep_old=0
    integer:: istep_new=0
    integer:: ihopp=0
    integer:: ihopp_acc=0
    integer:: ihopp_rej=0
    integer:: ibest=1
    real(8):: egap=1.d100
    real(8):: esep=0.d0
    !real(8):: erat
    real(8):: etoler=1.d-2
    real(8):: re_sm
    !real(8):: erathopp
    real(8):: dt !time step in MD.
    real(8):: ediff !upper bound for energy difference to accept the hopping.
    real(8):: count_md=0.d0
    real(8):: count_opt=0.d0
    real(8):: count_soften=0.d0
    real(8):: count_md_tot=0.d0
    real(8):: count_opt_tot=0.d0
    real(8):: count_soften_tot=0.d0
    real(8):: fcall_tot_all=0.d0
    real(8):: fcall_tot_all_md=0.d0
    real(8):: fcall_tot_all_opt=0.d0
    real(8):: fcall_tot_all_soften=0.d0
    logical:: newmin
    logical:: escaped
    logical:: accepted
    !logical:: write_results_now=.false.
    logical:: lprint=.true.
    logical:: die=.false.
    integer:: nrandoff=0
    integer:: nsoften=7
    real(8):: alpha_soften
    integer:: kerathopp=1 
    integer, allocatable:: itagintermediate(:)
    integer, allocatable:: mtagarr1(:)
    integer, allocatable:: mtagarr2(:)
    real(8), allocatable:: earr(:) !energy array.
    integer, allocatable:: nvisit(:)
    real(8), allocatable:: abuf(:)
    real(8), allocatable:: abufall(:,:)
    real(8), allocatable:: dtarr(:)
    real(8), allocatable:: ediffarr(:)
    real(8), allocatable:: ekinarr(:)
    real(8), allocatable:: abuf1(:,:)
    real(8), allocatable:: abuf2(:,:)
    logical, allocatable:: do_req1(:)
    logical, allocatable:: do_req2(:)
    integer, allocatable:: ireqarr1(:)
    integer, allocatable:: ireqarr2(:)
end module mod_minhopp
!*****************************************************************************************
