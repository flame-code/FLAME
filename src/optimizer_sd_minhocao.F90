!> @file
!!  Routines for Stefan's new minimization method
!! @author Bastian Schaefer
!! @section LICENCE
!!    Copyright (C) 2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/8l.txt .
!!    For the list of contributors, see ~/AUTHORS
!subroutine geopt(nat,wpos,etot,fout,fnrmtol,count,count_sd,displr)
!subroutine sqnm(nproc,iproc,verbosity,ncount_bigdft,fail,nat)
subroutine GEOPT_SD(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,counter,folder)
 use global, only: units,max_kpt,ka1,kb1,kc1,confine
 use steepest_descent
 use defs_basis
 use interface_code
!subroutine sqnm(runObj,outsIO,nproc,iproc,verbosity,ncount_bigdft,fail)
!call_bigdft has to be run once on runObj and outs !before calling this routine
!sqnm will return to caller the energies and coordinates used/obtained from the last accepted iteration step
!   use module_base
!   use bigdft_run!module_types
!   use yaml_output
   use module_sqn, only: modify_gradient_minhocao, getSubSpaceEvecEval_minhocao !, findbonds
   use mod_parini, only: typ_parini
   implicit none
   type(typ_parini), intent(in):: parini
   type(typ_parini), intent(inout):: parres
   !parameter
!   integer, intent(in)                    :: nproc
!   integer, intent(in)                    :: iproc
!   integer, intent(in)                    :: verbosity
!   type(run_objects), intent(inout)       :: runObj
!   type(DFT_global_output), intent(inout) :: outsIO
!   integer, intent(inout)                 :: ncount_bigdft
!   logical, intent(out)                   :: fail
   !local variables
   character(len=*), parameter :: subname='sqnm'
   integer :: infocode,info !< variables containing state codes
   integer :: nhistx !< maximum history length
   integer :: nhist  !< actual history length
   integer :: ndim   !< dimension of significant subspace
   integer :: nit    !< maximum number of iterations
!   integer :: nat    !< number of atoms
   integer :: istat,iall
   integer :: lwork
   integer :: it,i,iat,l,j,idim,jdim,ihist,icheck !<counter variables
   integer :: itswitch
   integer :: imode=1
   integer :: nbond=1
!   type(DFT_global_output) :: outs
   logical :: success=.false.
   logical :: debug !< set .true. for debug output to fort.100
   logical :: steep !< steepest descent flag
   real(8) :: displr !< (non-physical) integrated path length,
                      !< includes displacements from rejctions
                      !< (denoted as dsplr in geopt.mon)
   real(8) :: displp !< (physical) integrated path length,
                      !< includes NO displacements from rejections
                      !< (denoted as dsplp in geopt.mon)
   real(8) :: etot 
!   real(8) :: fnrm
   real(8) :: fmax
   real(8) :: fluct
   real(8) :: fnoise
   real(8) :: betax !< initial step size (gets not changed)
   real(8) :: beta_stretchx
   real(8) :: beta  !< current step size
   real(8) :: beta_stretch  !< current step size in bond-stretching directions
   real(8) :: cosangle
   real(8) :: tt
   real(8) :: maxd !< maximum displacement of single atom
   real(8) :: scl
   real(8) :: etotold
   real(8) :: detot
   real(8) :: etotp
   real(8) :: dt
   real(8) :: cutoffRatio !< if fraction of eigenvalues of overlapmatrix of
                           !< displacements is smaller that cutoffRatio,
                           !< then those displacements are regarded
                           !< as linear dependent and are not taken into account
                           !< for building the significant subspace
   real(8) :: maxrise !< energy ist not allowed to rise more than maxrise in single step
   real(8) :: steepthresh !< if fnrm is larger that steepthresh, steepest descent is used
   real(8) :: trustr !< a single atoms is not allowed to be dsiplaced more than by trustr
   integer,  allocatable, dimension(:,:)   :: iconnect
   real(8), allocatable, dimension(:,:,:) :: rxyz
   real(8), allocatable, dimension(:,:,:) :: rxyzraw
   real(8), allocatable, dimension(:,:,:) :: fxyz
   real(8), allocatable, dimension(:,:,:) :: fxyzraw
   real(8), allocatable, dimension(:,:,:) :: fstretch
   real(8), allocatable, dimension(:,:)   :: rxyzOld
   real(8), allocatable, dimension(:,:)   :: delta
   real(8), allocatable, dimension(:,:,:) :: ff
   real(8), allocatable, dimension(:,:,:) :: rr
   real(8), allocatable, dimension(:,:,:) :: rrr
   real(8), allocatable, dimension(:,:,:) :: fff
   real(8), allocatable, dimension(:,:)   :: aa
   real(8), allocatable, dimension(:,:)   :: dd
   real(8), allocatable, dimension(:)     :: eval
   real(8), allocatable, dimension(:)     :: work
   real(8), allocatable, dimension(:)     :: res
   real(8), allocatable, dimension(:)     :: scpr
   real(8), allocatable, dimension(:)     :: rnorm
   real(8), allocatable, dimension(:)     :: wold
   real(8), allocatable, dimension(:)     :: rcov
   character(len=4)                        :: fn4
   character(len=40)                       :: comment
   character(len=12)                        :: cdmy12_1
   character(len=12)                        :: cdmy12_2
   character(len=9)                        :: cdmy9
   character(len=8)                        :: cdmy8
    !functions
    real(8) :: ddot,dnrm2
!my own
   logical:: biomode,fail,getwfk
   real(8):: energy,forcemax,frac_fluct,fmax_at,fmax_lat,enthalpy,beta_scale
   integer:: ncount_cluster_x,iexit,iprec,lattdeg
   real(8):: latvec_in(3,3),xred_in(3,parini%nat),fcart_in(3,parini%nat),strten_in(6),etot_in,counter,pressure,latvec0(3,3),enthalpy_old
   character(40)::filename,folder
   frac_fluct=0.d0
   fluct=0.d0







   pressure=parini%target_pressure_habohr
   biomode=.false.
   !set parameters
!   nit=runObj%inputs%ncount_cluster_x
!   nat=runObj%nat
!   betax=runObj%inputs%betax
!   nhistx=runObj%inputs%nhistx
!   maxrise=runObj%inputs%maxrise
!   cutoffRatio=runObj%inputs%cutoffratio
!   steepthresh=runObj%inputs%steepthresh
!   trustr=runObj%inputs%trustr
   nit=        parini%paropt_geopt%nit
!   nat=        nat
   betax=       sd_beta_at
   beta_scale=  sd_beta_lat/sd_beta_at
   nhistx=     10!nhistx
   maxrise=    1.d-6!maxrise
   cutoffRatio=1.d-4!cutoffratio
   steepthresh=0.0000002d0!steepthresh
   trustr=     0.2d0!trustr
   if(biomode)imode=2

!!   if (iproc==0.and.verbosity > 0) then
!!      call yaml_mapping_open('Geometry parameters')
!!         call yaml_map('Geometry Method','GEOPT_SQNM')
!!         call yaml_map('nhistx',nhistx)
!!         call yaml_map('biomode',runObj%inputs%biomode)
!!         call yaml_map('betax', betax,fmt='(1pe21.14)')
!!         call yaml_map('beta_stretchx', runObj%inputs%beta_stretchx,fmt='(1pe21.14)')
!!         call yaml_map('maxrise', maxrise,fmt='(1pe21.14)')
!!         call yaml_map('cutoffRatio', cutoffRatio,fmt='(1pe21.14)')
!!         call yaml_map('steepthresh', steepthresh,fmt='(1pe21.14)')
!!         call yaml_map('trustr', trustr,fmt='(1pe21.14)')
!!      call yaml_mapping_close()
!!   end if

   !init variables
   debug=.false.
   fail=.true.
   displr=0.0d0
   displp=0.0d0
   fluct=0.0d0
   icheck=0
   detot=0.0d0
   itswitch=0
   ndim=0
   nhist=0
   beta=betax
!   beta_stretch=runObj%inputs%beta_stretchx
   beta_stretch=1.d0!beta_stretchx
   maxd=1.0d0

   ! allocate arrays
!   lwork=1000+10*nat**2
   lwork=1000+10*(parini%nat+3)**2
!   rxyz = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='rxyz')
!   rxyzraw = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='rxyzraw')
!   fxyz = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='fxyz')
!   fxyzraw = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='fxyzraw')
!   fstretch = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='fstretch')
!   rxyzOld = f_malloc((/ 1.to.3, 1.to.nat/),id='rxyzOld')
!   delta = f_malloc((/ 1.to.3, 1.to.nat/),id='delta')
!   aa = f_malloc((/ nhistx, nhistx /),id='aa')
!   eval = f_malloc(nhistx,id='eval')
!   res = f_malloc(nhistx,id='res')
!   rnorm = f_malloc(nhistx,id='rnorm')
!   work = f_malloc(lwork,id='work')
!   ff = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='ff')
!   rr = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='rr')
!   dd = f_malloc((/ 3, nat /),id='dd')
!   fff = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='fff')
!   rrr = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='rrr')
!   scpr = f_malloc(nhistx,id='scpr')
!   rcov     = f_malloc((/ 1.to.nat/),id='rcov')
!   iconnect = f_malloc((/ 1.to.2, 1.to.1000/),id='iconnect')
   allocate(rxyz(1:3,1:parini%nat+3,0:nhistx))
   allocate(rxyzraw(1:3,1:parini%nat+3,0:nhistx))
   allocate(fxyz(1:3,1:parini%nat+3,0:nhistx)) 
   allocate(fxyzraw(1:3,1:parini%nat+3,0:nhistx)) 
   allocate(fstretch(1:3,1:parini%nat+3,0:nhistx)) 
   allocate(rxyzOld(1:3,1:parini%nat+3)) 
   allocate(delta(1:3,1:parini%nat+3)) 
   allocate(aa(nhistx,nhistx)) 
   allocate(eval(nhistx)) 
   allocate(res(nhistx)) 
   allocate(rnorm(nhistx)) 
   allocate(work(lwork)) 
   allocate(ff(1:3,1:parini%nat+3,0:nhistx)) 
   allocate(rr(1:3,1:parini%nat+3,0:nhistx)) 
   allocate(dd(3,parini%nat+3)) 
   allocate(fff(1:3,1:parini%nat+3,0:nhistx)) 
   allocate(rrr(1:3,1:parini%nat+3,0:nhistx)) 
   allocate(scpr(nhistx)) 
!   allocate(rcov(1:nat))     
!   allocate(iconnect(1:2,1:1000)) 
!   if(runObj%inputs%biomode)then
!        call give_rcov_sqnm(iproc,runObj%atoms,runObj%nat,rcov)
!        call findbonds('(SQNM)',iproc,10,runObj%nat,&
!             rcov,runObj%rxyz,nbond,iconnect)
!   endif 
!   wold = f_malloc((/ 1.to.nbond/),id='wold')
!   allocate(wold(1:nbond))
!   wold =0.0d0


!   call init_global_output(outs, runObj%nat)

   !copy outs_datatype
!   call copy_global_output(outsIO,outs)



!!!!!!   call energyandforces(nat,rxyz(1,1,0),fxyz(1,1,0),etot)
!!  not necessary, call_bigdft allready called outside
!   call call_bigdft(runObj,outs,nproc,iproc,infocode)
!   ncount_bigdft=ncount_bigdft+1

!! copy to internal variables
!   call vcopy(3*runObj%nat, runObj%rxyz(1,1), 1,rxyz(1,1,0), 1)
   rxyz(:,1:parini%nat,0)=xred_in(:,:);rxyz(:,parini%nat+1:parini%nat+3,0)=latvec_in(:,:)
!   call vcopy(3*runObj%nat, runObj%rxyz(1,1), 1,rxyzOld(1,1), 1)
!   call vcopy(3*outs%fdim, outs%fxyz(1,1), 1, fxyz(1,1,0), 1)
!   fxyz(:,1:3*nat,0)=_in(:,:);rxyz(:,3*nat+1:3*nat+3,0)=latvec_in(:,:)
!   fxyz(:,:,0)=fxyz0(:,:)
!   etot=outs%energy

rxyz(:,:,1)=rxyz(:,:,0) 
fxyz(:,:,1)=0.d0
enthalpy_old=10.d0
counter=0.d0
do it=1,nit
      lattdeg=1
      call get_BFGS_forces_strainlatt(parini,parres,rxyz(:,:,0),fxyz(:,:,0),enthalpy,getwfk,iprec,latvec0,&
             &lattdeg,latvec_in,xred_in,etot_in,fcart_in,strten_in)
      call get_enthalpy(latvec_in,etot_in,pressure,enthalpy)
      call convcheck(parini,parini%nat,latvec_in,fcart_in,strten_in,parini%target_pressure_habohr,parini%paropt_geopt%strfact,fmax,fmax_at,fmax_lat,parini%paropt_geopt%fmaxtol,iexit)
      counter=real(it,8)
      write(*,'(i8,a10,5es15.7)') it," GEOPT_SD ",enthalpy,fmax,fmax_at,fmax_lat,beta
      if(iexit==1) then
            write(*,*) " SD converged in iteration",it
            return
      endif 
      if(enthalpy.lt.enthalpy_old) then
         beta=beta*1.05d0
         enthalpy_old=enthalpy
         fxyz(:,:,1)=fxyz(:,:,0)
         rxyz(:,:,1)=rxyz(:,:,0)
      else
         beta=max(beta*0.5d0,1.d-4)
         if(beta.gt.1.d-4) then
         fxyz(:,:,0)=fxyz(:,:,1)
         rxyz(:,:,0)=rxyz(:,:,1)
         else
         enthalpy_old=enthalpy
         endif
      endif
      rxyz(:,parini%nat+1:parini%nat+3,0)=rxyz(:,parini%nat+1:parini%nat+3,0)+beta*beta_scale*fxyz(:,parini%nat+1:parini%nat+3,0)
      rxyz(:,1:parini%nat,0)=rxyz(:,1:parini%nat,0)+beta*fxyz(:,1:parini%nat,0)
enddo

deallocate(rxyz)
deallocate(rxyzOld)
deallocate(delta)
deallocate(rxyzraw)
deallocate(fxyz)
deallocate(fxyzraw)
deallocate(fstretch)
deallocate(aa)
deallocate(eval)
deallocate(res)
deallocate(rnorm)
deallocate(work)
deallocate(ff)
deallocate(rr)
deallocate(dd)
deallocate(fff)
deallocate(rrr)
deallocate(scpr)
!deallocate(wold)
!deallocate(rcov )   
!deallocate(iconnect)
!   call deallocate_global_output(outs)
end subroutine
