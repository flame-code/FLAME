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
subroutine GEOPT_sqnm(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,counter,folder)
 use global, only: units,max_kpt,ka1,kb1,kc1,confine
 use defs_basis
 use interface_code
 use modsocket, only: sock_extra_string
!subroutine sqnm(runObj,outsIO,nproc,iproc,verbosity,ncount_bigdft,fail)
!call_bigdft has to be run once on runObj and outs !before calling this routine
!sqnm will return to caller the energies and coordinates used/obtained from the last accepted iteration step
!   use module_base
!   use bigdft_run!module_types
!   use yaml_output
   use module_sqn, only: modify_gradient_minhocao, getSubSpaceEvecEval_minhocao!, findbonds
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
   real(8) :: fnoise
   real(8) :: betax !< initial step size (gets not changed)
   real(8) :: betalatx !< initial step size (gets not changed)
   real(8) :: betalat_scale
   real(8) :: beta_stretchx
   real(8) :: beta  !< current step size
   real(8) :: betalat  !< current step size
   real(8) :: beta_stretch  !< current step size in bond-stretching directions
   real(8) :: cosangle
   real(8) :: cosangle_lat
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
   real(8):: energy,forcemax,fmax_at,fmax_lat,enthalpy
   integer:: ncount_cluster_x,iexit,iprec
   real(8):: latvec_in(3,3),xred_in(3,parini%nat),fcart_in(3,parini%nat),strten_in(6),etot_in,counter,pressure
   character(40)::filename,folder
   logical:: multiprec,cellfix_done
   real(8)::tolmxf_switch,cellfix_switch
   real(8),allocatable:: hessinv(:,:),metric(:,:)
!For running on cartesian forces
   logical:: cart_forces
   real(8):: pos_tmp(3,parini%nat),latvec_old(3,3)
!Latvec correction io
   integer:: latvec_io
latvec_io=0
!multiprec is hardcoded and, if true, starts a geopt with iprec==2, and then switches 
!to iprec==1 when the fmax==tolmxf_switch. The switch only occurs once
 max_kpt=.false.
 multiprec=.true.
 tolmxf_switch=1.d-3
 cellfix_switch=1.d-3
 cellfix_done=.false.
 cart_forces=.false.


   debug=.false.


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
   betax=      parini%paropt_geopt%beta_at !betax
   betalatx=   parini%paropt_geopt%beta_lat !betax
   betalat_scale=betalatx/betax
   nhistx=     parini%paropt_geopt%nhist!nhistx
   maxrise=    parini%paropt_geopt%maxrise!1.d-6!maxrise
   cutoffRatio=parini%paropt_geopt%cutoffRatio!2.d-4!cutoffratio
   steepthresh=parini%paropt_geopt%steepthresh!1000.d0!steepthresh
   trustr=     parini%paropt_geopt%trustr!0.2d0!trustr
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
   icheck=0
   detot=0.0d0
   itswitch=0
   ndim=0
   nhist=0
   beta=betax
   betalat=betalatx
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
   allocate(hessinv(3*(parini%nat+3),3*(parini%nat+3)),metric(3*(parini%nat+3),3*(parini%nat+3)))
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
if(cart_forces) then
   call rxyz_int2cart(latvec_in,xred_in,rxyz(:,1:parini%nat,0),parini%nat);rxyz(:,parini%nat+1:parini%nat+3,0)=latvec_in(:,:)
   call rxyz_int2cart(latvec_in,xred_in,rxyzOld(:,1:parini%nat),parini%nat);rxyzOld(:,parini%nat+1:parini%nat+3)=latvec_in(:,:)
else
!   call vcopy(3*runObj%nat, runObj%rxyz(1,1), 1,rxyz(1,1,0), 1)
   rxyz(:,1:parini%nat,0)=xred_in(:,:);rxyz(:,parini%nat+1:parini%nat+3,0)=latvec_in(:,:)
!   call vcopy(3*runObj%nat, runObj%rxyz(1,1), 1,rxyzOld(1,1), 1)
   rxyzOld(:,1:parini%nat)=xred_in(:,:);rxyzOld(:,parini%nat+1:parini%nat+3)=latvec_in(:,:)
endif
!   call vcopy(3*outs%fdim, outs%fxyz(1,1), 1, fxyz(1,1,0), 1)
!   fxyz(:,1:3*nat,0)=_in(:,:);rxyz(:,3*nat+1:3*nat+3,0)=latvec_in(:,:)
!   fxyz(:,:,0)=fxyz0(:,:)
!   etot=outs%energy
!   etot=energy
   
call sqnm_invhess(parini%nat,latvec_in,metric,hessinv)


!   call minenergyandforces(parini,parres,iproc,nproc,.false.,imode,runObj,outs,nat,rxyz(1,1,0),&
!       rxyzraw(1,1,0),fxyz(1,1,0),fstretch(1,1,0),fxyzraw(1,1,0),&
!       etot,iconnect,nbond,wold,beta_stretchx,beta_stretch)
!   call minenergyandforces(parini,parres,iproc,nproc,.true.,imode,nat,rxyz(1,1,0),&
!       rxyzraw(1,1,0),fxyz(1,1,0),fstretch(1,1,0),fxyzraw(1,1,0),&
!       etot,iconnect,nbond,wold,beta_stretchx,beta_stretch)
   write(fn4,'(i4.4)') 0
   sock_extra_string="SQNM"//trim(fn4)
   if(cart_forces) then
          call rxyz_cart2int(rxyz(:,parini%nat+1:parini%nat+3,0),pos_tmp,rxyz(:,1:parini%nat,0),parini%nat)
          rxyz(:,1:parini%nat,0)=pos_tmp(:,:)
   endif
   call minenergyandforces(parini,parres,.true.,imode,parini%nat,rxyz(1,1,0),&
       rxyzraw(1,1,0),fxyz(1,1,0),fstretch(1,1,0),fxyzraw(1,1,0),&
       etot,beta_stretchx,beta_stretch,&
       latvec_in,xred_in,etot_in,fcart_in,strten_in,iprec)
   if(cart_forces) then
       call rxyz_int2cart(rxyz(:,parini%nat+1:parini%nat+3,0),xred_in,rxyz(:,1:parini%nat,0),parini%nat) 
       fxyz(:,1:parini%nat,0)=fcart_in(:,:)
   endif
   if(imode==2)rxyz(:,:,0)=rxyz(:,:,0)+beta_stretch*fstretch(:,:,0)

!   call fnrmandforcemax(fxyzraw(1,1,0),fnrm,fmax,nat)
   call convcheck(parini,parini%nat,latvec_in,fcart_in,strten_in,parini%target_pressure_habohr,parini%paropt_geopt%strfact,fmax,fmax_at,fmax_lat,parini%paropt_geopt%fmaxtol,iexit)
!   fnrm=sqrt(fnrm)
!   if (fmax < 3.d-1) call updatefluctsum(outs%fnoise,fluct)

   etotold=1.d100!etot
   etotp=etot

!   if (iproc==0.and.verbosity > 0) then
   it=0
if(parini%verb>0) then
       !avoid space for leading sign (numbers are positive, anyway)
       write(cdmy8,'(es8.1)')abs(maxd)
       write(cdmy12_1,'(es12.5)')abs(displr)
       write(cdmy12_2,'(es12.5)')abs(displp)
       write(cdmy9,'(es9.2)')abs(beta)


!       write(16,'(i5,1x,i5,2x,a10,2x,1es21.14,2x,es9.2,es11.3,3es10.2,2x,a6,a8,1x,a4,i3.3,1x,a5,a7,2(1x,a6,a11))') &
!       ncount_bigdft,0,'GEOPT_SQNM',etotp,detot,fmax,fnrm,fluct*runObj%inputs%frac_fluct,fluct, &
!       'beta=',trim(adjustl(cdmy9)),'dim=',ndim,'maxd=',&
!       trim(adjustl(cdmy8)),'dsplr=',trim(adjustl(cdmy12_1)),&
!       'dsplp=',trim(adjustl(cdmy12_2))
       write(*,'(a,i4,2x,1es21.14,2x,es9.2,3es10.3,2x,a6,a8,1x,a,es9.2,1x,a4,i3.3,1x,a5,a7,2(1x,a6,a11),i4)') &
       " # GEOPT SQNM    ",it,etotp,detot,fmax,fmax_at,fmax_lat, &
       'beta=',trim(adjustl(cdmy9)),"betalat=",betalat,'dim=',ndim,'maxd=',&
       trim(adjustl(cdmy8)),'dsplr=',trim(adjustl(cdmy12_1)),&
       'dsplp=',trim(adjustl(cdmy12_2)),iprec
!MHM: Write output to file in every step***********************************
       counter=real(it,8)
       write(*,*) "Pressure, Energy",pressure,etot_in
       call get_enthalpy(latvec_in,etot_in,pressure,enthalpy)
       write(fn4,'(i4.4)') it
       filename=trim(folder)//"posgeopt."//fn4//".ascii"
       units=units
       write(*,*) "# Writing the positions in SQNM:",filename
       call write_atomic_file_ascii(parini,filename,parini%nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &parini%char_type(1:parini%ntypat_global),parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,etot_in,pressure,enthalpy,etotp)
       if(parini%verb.ge.3) then
       filename=trim(folder)//"posgeopt."//fn4//".vasp"
       call write_atomic_file_poscar(parini,filename,parini%nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &parini%char_type(1:parini%ntypat_global),parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,etot_in,pressure,enthalpy,etotp)
       endif
endif


!Initial iprec after running the first force call
 if(multiprec) iprec=2
   etotold=1.d100!etot
   etotp=etot
       if(iprec==2.and.multiprec.and.fmax.lt.tolmxf_switch) then
          iprec=1
          etotold=etot
       endif

   if(iexit==1) then
     write(*,'(a,i4,2(1x,es25.15))') " # SQNM converged before entering iterations", it,enthalpy,fmax
     max_kpt=.false.
     return 
   endif
   do it=1,nit!start main loop
!  do it=1,nit-1!start main loop (nit-1 if first bigdft call is NOT done outside, but inside this subroutine)
!      if (debug.and.iproc==0) write(100,*) 'it:',it,etot,fnrm,itswitch
      nhist=nhist+1

!      if (fnrm.gt.steepthresh .or. it.le.itswitch ) then
      if (fmax.gt.steepthresh .or. it.le.itswitch ) then
         ndim=0
         steep=.true.
         if (it.gt.itswitch) itswitch=it+nhistx
!         if (debug.and.iproc==0) write(100,*) "STEEP"
      else
         steep=.false.
      endif

      ! make space in the history list
      if (nhist.gt.nhistx) then
         nhist=nhistx
         do ihist=0,nhist-1
!            do iat=1,nat
            do iat=1,parini%nat+3
               do l=1,3
                  rxyz(l,iat,ihist)=rxyz(l,iat,ihist+1)
                  rxyzraw(l,iat,ihist)=rxyzraw(l,iat,ihist+1)
                  fxyz(l,iat,ihist)=fxyz(l,iat,ihist+1)
                  fxyzraw(l,iat,ihist)=fxyzraw(l,iat,ihist+1)
                  fstretch(l,iat,ihist)=fstretch(l,iat,ihist+1)
               enddo
            enddo
         enddo
      endif
   
      ! decompose gradient
500 continue
    call modify_gradient_minhocao(parini%nat,ndim,rrr(1,1,1),eval(1),res(1),fxyz(1,1,nhist-1),beta,betalat_scale,dd(1,1))
   
      tt=0.0d0
      dt=0.0d0
      maxd=-huge(1.0d0)
      do iat=1,parini%nat
         dt=dd(1,iat)**2+dd(2,iat)**2+dd(3,iat)**2
         tt=tt+dt
         maxd=max(maxd,dt)
      enddo
!Take care of the lattice maximal displacement
      do iat=parini%nat+1,parini%nat+3
         dt=(dd(1,iat)**2+dd(2,iat)**2+dd(3,iat)**2)/betalat_scale**2
         tt=tt+dt
         maxd=max(maxd,dt)
      enddo
      tt=sqrt(tt)
      maxd=sqrt(maxd)
   
      !trust radius approach: avoids too large steps due to large forces
      !only used when in steepest decent mode
      if(maxd>trustr .and. steep)then
!         if(debug.and.iproc==0)write(100,'(a,1x,es24.17,1x,i0)')'step too large',maxd,it
         write(*,'(a,1x,es24.17,1x,i0)')'step too large',maxd,it
!         if(iproc==0)write(16,'(a,2(1x,es9.2))')'WARNING GEOPT_SQNM: step too large: maxd, trustradius ',maxd,trustr
         scl=0.50d0*trustr/maxd
         write(*,'(a,3(1x,es9.2))')'WARNING: step too large: maxd, trustradius ',maxd,trustr,scl
         dd=dd*scl
         tt=tt*scl
         maxd=maxd*scl
!Also reset beta
         beta=max(0.5d0*beta,0.1d0*betax)
         betalat=max(0.5d0*betalat,0.1d0*betalatx)
      endif
!      displr=displr+tt
   
      !update positions
!      dd(:,nat+1:nat+3)=0.d0
      if(cart_forces) then
        latvec_old(:,:)=rxyz(:,parini%nat+1:parini%nat+3,nhist-1)
        do iat=1,parini%nat+3
           rxyz(1,iat,nhist)=rxyz(1,iat,nhist-1)-dd(1,iat)
           rxyz(2,iat,nhist)=rxyz(2,iat,nhist-1)-dd(2,iat)
           rxyz(3,iat,nhist)=rxyz(3,iat,nhist-1)-dd(3,iat)
        enddo
        call updaterxyz(latvec_old,rxyz(:,parini%nat+1:parini%nat+3,nhist),rxyz(:,1:parini%nat,nhist),parini%nat)
      else
        do iat=1,parini%nat+3
           rxyz(1,iat,nhist)=rxyz(1,iat,nhist-1)-dd(1,iat)
           rxyz(2,iat,nhist)=rxyz(2,iat,nhist-1)-dd(2,iat)
           rxyz(3,iat,nhist)=rxyz(3,iat,nhist-1)-dd(3,iat)
        enddo
      endif
   
!      call energyandforces(nat,rxyz(1,1,nhist),fxyz(1,1,nhist),etotp)
!      call vcopy(3 * runObj%nat, rxyz(1,1,nhist), 1,runObj%rxyz(1,1), 1)
!      runObj%inputs%inputPsiId=1
!      call call_bigdft(runObj,outs,nproc,iproc,infocode)
!      ncount_bigdft=ncount_bigdft+1
!      call vcopy(3 * outs%fdim, outs%fxyz(1,1), 1, fxyz(1,1,nhist), 1)
!      etotp=outs%energy
!      detot=etotp-etotold

      delta=rxyz(:,:,nhist)-rxyzOld
      displr=displr+dnrm2(3*parini%nat+9,delta(1,1),1)
!      runObj%inputs%inputPsiId=1
!      call minenergyandforces(parini,parres,iproc,nproc,.true.,imode,runObj,outs,nat,rxyz(1,1,nhist),rxyzraw(1,1,nhist),&
!                             fxyz(1,1,nhist),fstretch(1,1,nhist),fxyzraw(1,1,nhist),&
!                             etotp,iconnect,nbond,wold,beta_stretchx,beta_stretch)
       write(fn4,'(i4.4)') it
       sock_extra_string="SQNM"//trim(fn4)
   if(cart_forces) then
          call rxyz_cart2int(rxyz(:,parini%nat+1:parini%nat+3,nhist),pos_tmp,rxyz(:,1:parini%nat,nhist),parini%nat)
          rxyz(:,1:parini%nat,nhist)=pos_tmp(:,:)
   endif
       call minenergyandforces(parini,parres,.true.,imode,parini%nat,rxyz(1,1,nhist),&
           rxyzraw(1,1,nhist),fxyz(1,1,nhist),fstretch(1,1,nhist),fxyzraw(1,1,nhist),&
           etotp,beta_stretchx,beta_stretch,&
           latvec_in,xred_in,etot_in,fcart_in,strten_in,iprec)
   if(cart_forces) then
       call rxyz_int2cart(rxyz(:,parini%nat+1:parini%nat+3,nhist),xred_in,rxyz(:,1:parini%nat,nhist),parini%nat) 
       fxyz(:,1:parini%nat,nhist)=fcart_in(:,:)
   endif
      detot=etotp-etotold
!      ncount_bigdft=ncount_bigdft+1
      counter=counter+1


!      call fnrmandforcemax(fxyzraw(1,1,nhist),fnrm,fmax,nat)
      call convcheck(parini,parini%nat,latvec_in,fcart_in,strten_in,parini%target_pressure_habohr,parini%paropt_geopt%strfact,fmax,fmax_at,fmax_lat,parini%paropt_geopt%fmaxtol,iexit)
!MHM: Write output to file in every step***********************************
       counter=real(it,8)
       write(*,*) "Pressure, Energy",pressure,etot_in
       call get_enthalpy(latvec_in,etot_in,pressure,enthalpy)
if(parini%verb.gt.0) then
       write(fn4,'(i4.4)') it
       filename=trim(folder)//"posgeopt."//fn4//".ascii"
       units=units
       write(*,*) "# Writing the positions in SQNM:",filename
       call write_atomic_file_ascii(parini,filename,parini%nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &parini%char_type(1:parini%ntypat_global),parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,etot_in,pressure,enthalpy,etotp)
       if(parini%verb.ge.3) then
       filename=trim(folder)//"posgeopt."//fn4//".vasp"
       call write_atomic_file_poscar(parini,filename,parini%nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &parini%char_type(1:parini%ntypat_global),parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,etot_in,pressure,enthalpy,etotp)
       endif
endif
!*********************************************************************
         write(*,'(a,i4,2x,1es21.14,2x,es9.2,3es10.3,2x,a6,a8,1x,a,es9.2,1x,a4,i3.3,1x,a5,a7,2(1x,a6,a11),i4)') &
              " # GEOPT SQNM    ",it,etotp,detot,fmax,fmax_at,fmax_lat, &
          'beta=',trim(adjustl(cdmy9)),"betalat=",betalat,'dim=',ndim,'maxd=',&
          trim(adjustl(cdmy8)),'dsplr=',trim(adjustl(cdmy12_1)),&
          'dsplp=',trim(adjustl(cdmy12_2)),iprec
!*********************************************************************
!      fnrm=sqrt(fnrm)

!      if (iproc == 0) then
!         write(fn4,'(i4.4)') ncount_bigdft
!         write(fn4,'(i4.4)') it
!         write(comment,'(a,1pe10.3)')'SQNM:fnrm= ',fnrm
         !call bigdft_write_atomic_file(runObj,outs,'posout_'//fn4,&
         !    trim(comment))
!!$
!!$         call write_atomic_file(trim(runObj%inputs%dir_output)//'posout_'//fn4, &
!!$              outs%energy,runObj%rxyz,runObj%ixyz_int,&
!!$              runObj%atoms,trim(comment),forces=outs%fxyz)
!      endif
!!MHM: Write output to file in every step***********************************
!       write(*,*) "Pressure, Energy",pressure,etot_in
!       call get_enthalpy(latvec_in,etot_in,pressure,enthalpy)
!if(parini%verb.gt.0) then
!       write(fn4,'(i4.4)') it
!       filename=trim(folder)//"posgeoptP."//fn4//".ascii"
!       units=units
!       write(*,*) "# Writing the positions in SQNM:",filename
!       call write_atomic_file_ascii(parini,filename,nat,units,xred_in,latvec_in,fcart_in,strten_in,&
!            &char_type(1:ntypat),ntypat,typat,fixat,fixlat,etot_in,pressure,enthalpy,etotp)
!       if(parini%verb.ge.3) then
!       filename=trim(folder)//"posgeoptP."//fn4//".vasp"
!       call write_atomic_file_poscar(parini,filename,nat,units,xred_in,latvec_in,fcart_in,strten_in,&
!            &char_type(1:ntypat),ntypat,typat,fixat,fixlat,etot_in,pressure,enthalpy,etotp)
!       endif
!endif
!!*********************************************************************

!      if (fmax < 3.d-1) call updatefluctsum(outs%fnoise,fluct)
!      cosangle=-dot_double(3*nat,fxyz(1,1,nhist),1,dd(1,1),1)/&
!              sqrt(dot_double(3*nat,fxyz(1,1,nhist),1,fxyz(1,1,nhist),1)*&
!              dot_double(3*nat,dd(1,1),1,dd(1,1),1))
!      cosangle=-dot_product(reshape(fxyz(:,:,nhist),(/3*nat+9/)),reshape(dd(:,:),(/3*nat+9/)))/&
!              sqrt(dot_product(reshape(fxyz(:,:,nhist),(/3*nat+9/)),reshape(fxyz(:,:,nhist),(/3*nat+9/))*&
!              dot_product(reshape(dd(:,:),(/3*nat+9/)),reshape(dd(:,:),(/3*nat+9/)))))
      cosangle=-dot_product(reshape(fxyz(:,1:parini%nat,nhist),(/3*parini%nat/)),reshape(dd(:,1:parini%nat),(/3*parini%nat/)))/&
              sqrt(dot_product(reshape(fxyz(:,1:parini%nat,nhist),(/3*parini%nat/)),reshape(fxyz(:,1:parini%nat,nhist),(/3*parini%nat/))*&
              dot_product(reshape(dd(:,1:parini%nat),(/3*parini%nat/)),reshape(dd(:,1:parini%nat),(/3*parini%nat/)))))
      cosangle_lat=-dot_product(reshape(fxyz(:,parini%nat+1:parini%nat+3,nhist),(/9/)),reshape(dd(:,parini%nat+1:parini%nat+3),(/9/)))/&
              sqrt(dot_product(reshape(fxyz(:,parini%nat+1:parini%nat+3,nhist),(/9/)),reshape(fxyz(:,parini%nat+1:parini%nat+3,nhist),(/9/))*&
              dot_product(reshape(dd(:,parini%nat+1:parini%nat+3),(/9/)),reshape(dd(:,parini%nat+1:parini%nat+3),(/9/)))))

      if (detot.gt.maxrise .and. beta > 1.d-1*betax .and. betalat > 1.d-1*betax) then !
!         if (debug.and.iproc==0) write(100,'(a,i0,1x,e9.2)') "WARN: it,detot", it,detot
         write(*,'(a,i0,1x,e9.2)') "WARN: it,detot", it,detot
!         if (debug.and.iproc==0) write(16,'(a,i0,4(1x,e9.2))') &
         write(*,'(a,i0,4(1x,e9.2))') &
             "WARNING: Prevent energy to rise by more than maxrise: it,maxrise,detot,beta,1.e-1*betax ",&
             it,maxrise,detot,beta,1.d-1*betax
!         if (iproc==0.and.verbosity > 0) then
         if (parini%verb > 0) then
            !avoid space for leading sign (numbers are positive, anyway)
            write(cdmy8,'(es8.1)')abs(maxd)
            write(cdmy12_1,'(es12.5)')abs(displr)
            write(cdmy12_2,'(es12.5)')abs(displp)
            write(cdmy9,'(es9.2)')abs(beta)

   
!            write(16,'(i5,1x,i5,2x,a10,2x,1es21.14,2x,es9.2,es11.3,3es10.2,2x,a6,a8,1x,a4,i3.3,1x,a5,a7,2(1x,a6,a11))') &
!             ncount_bigdft,it,'GEOPT_SQNM',etotp,detot,fmax,fnrm,fluct*runObj%inputs%frac_fluct,fluct, &
!             'beta=',trim(adjustl(cdmy9)),'dim=',ndim,&
!             'maxd=',trim(adjustl(cdmy8)),'dsplr=',trim(adjustl(cdmy12_1)),&
!             'dsplp=',trim(adjustl(cdmy12_2))
!!!!!            write(*,'(a,i4,2x,1es21.14,2x,es9.2,3es10.3,2x,a6,a8,1x,a4,i3.3,1x,a5,a7,2(1x,a6,a11),i4)') &
!!!!!             " # GEOPT SQNM    ",it,etotp,detot,fmax,fmax_at,fmax_lat,&
!!!!!             'beta=',trim(adjustl(cdmy9)),'dim=',ndim,&
!!!!!             'maxd=',trim(adjustl(cdmy8)),'dsplr=',trim(adjustl(cdmy12_1)),&
!!!!!             'dsplp=',trim(adjustl(cdmy12_2)),iprec
!            call yaml_mapping_open('Geometry')
!               call yaml_map('Ncount_BigDFT',ncount_bigdft)
!               call yaml_map('Geometry step',it)
!               call yaml_map('Geometry Method','GEOPT_SQNM')
!               call yaml_map('ndim',ndim)
!               call yaml_map('etot', etotp,fmt='(1pe21.14)')
!               call yaml_map('detot',detot,fmt='(1pe21.14)')
!               call yaml_map('fmax',fmax,fmt='(1pe21.14)')
!               call yaml_map('fnrm',fnrm,fmt='(1pe21.14)')
!               call yaml_map('beta',beta,fmt='(1pe21.14)')
!               call yaml_map('beta_stretch',beta_stretch,fmt='(1pe21.14)')
!               call geometry_output(fmax,fnrm,fluct)
!            call yaml_mapping_close()
         end if
    
!         if(ncount_bigdft >= nit)then!no convergence within ncount_cluster_x energy evaluations
!!         if(it >= nit)then!no convergence within ncount_cluster_x energy evaluations
!!            !following copy of rxyz(1,1,nhist-1) to runObj is necessary for returning to the caller
!!            !the energies and coordinates used/obtained from/in the last ACCEPTED iteration step
!!            !(otherwise coordinates of last call to call_bigdft would be returned)
!!!            call vcopy(3 * runObj%nat, rxyz(1,1,nhist-1), 1,runObj%rxyz(1,1), 1)
!!            latvec_in=rxyz(:,nat+1:nat+3,nhist-1)
!!            xred_in  =rxyz(:,1:nat,nhist-1)
!!            goto 900  !sqnm will return to caller the energies and coordinates used/obtained from the last ACCEPTED iteration step
!!         endif

!         beta=min(.50d0*beta,betax)
         beta=.50d0*beta
         betalat=.50d0*betalat
!         if (debug.and.iproc==0) write(100,'(a,1x,e9.2)') 'WARNING GEOPT_SQNM: beta reset ',beta
         write(*,'(a,1x,e9.2)') 'WARNING: beta reset ',beta
         ndim=0
!         wold=0.0d0
         if(.not.steep)then
            do iat=1,parini%nat+3
               rxyz(1,iat,0)=rxyz(1,iat,nhist-1)
               rxyz(2,iat,0)=rxyz(2,iat,nhist-1)
               rxyz(3,iat,0)=rxyz(3,iat,nhist-1)
               rxyzraw(1,iat,0)=rxyzraw(1,iat,nhist-1)
               rxyzraw(2,iat,0)=rxyzraw(2,iat,nhist-1)
               rxyzraw(3,iat,0)=rxyzraw(3,iat,nhist-1)
   
               fxyz(1,iat,0)=fxyz(1,iat,nhist-1)
               fxyz(2,iat,0)=fxyz(2,iat,nhist-1)
               fxyz(3,iat,0)=fxyz(3,iat,nhist-1)
               fxyzraw(1,iat,0)=fxyzraw(1,iat,nhist-1)
               fxyzraw(2,iat,0)=fxyzraw(2,iat,nhist-1)
               fxyzraw(3,iat,0)=fxyzraw(3,iat,nhist-1)
            enddo
            nhist=1
         endif
         goto  500
      endif

!      if (iproc == 0) then
!         write(fn4,'(i4.4)') it
!         write(comment,'(a,1pe10.3)')'SQNM:fnrm= ',fnrm
!         !call bigdft_write_atomic_file(runObj,outs,'posoutP_'//fn4,&
!         !     trim(comment))
!      endif

      delta=rxyz(:,:,nhist)-rxyzOld
      displp=displp+dnrm2(3*parini%nat+9,delta(1,1),1)
      rxyzOld=rxyz(:,:,nhist)
!      displp=displp+tt
!      if (iproc==0.and.verbosity > 0) then
      if (parini%verb > 0) then
         !avoid space for leading sign (numbers are positive, anyway)
         write(cdmy8,'(es8.1)')abs(maxd)
         write(cdmy12_1,'(es12.5)')abs(displr)
         write(cdmy12_2,'(es12.5)')abs(displp)
         write(cdmy9,'(es9.2)')abs(beta)


!         write(16,'(i5,1x,i5,2x,a10,2x,1es21.14,2x,es9.2,es11.3,3es10.2,2x,a6,a8,1x,a4,i3.3,1x,a5,a7,2(1x,a6,a11))') &
!          ncount_bigdft,it,'GEOPT_SQNM',etotp,detot,fmax,fnrm,fluct*runObj%inputs%frac_fluct,fluct, &
!          'beta=',trim(adjustl(cdmy9)),'dim=',ndim,'maxd=',&
!          trim(adjustl(cdmy8)),'dsplr=',trim(adjustl(cdmy12_1)),&
!          'dsplp=',trim(adjustl(cdmy12_2))
!         write(16,'(i5,1x,i5,2x,a10,2x,1es21.14,2x,es9.2,es11.3,3es10.2,2x,a6,a8,1x,a4,i3.3,1x,a5,a7,2(1x,a6,a11))') &
!          ncount_bigdft,it,'GEOPT_SQNM',etotp,detot,fmax,fnrm,fluct*frac_fluct,fluct, &
!         call yaml_mapping_open('Geometry')
!            call yaml_map('Ncount_BigDFT',ncount_bigdft)
!            call yaml_map('Geometry step',it)
!            call yaml_map('Geometry Method','GEOPT_SQNM')
!            call yaml_map('ndim',ndim)
!            call yaml_map('etot', etotp,fmt='(1pe21.14)')
!            call yaml_map('detot',detot,fmt='(1pe21.14)')
!            call yaml_map('fmax',fmax,fmt='(1pe21.14)')
!            call yaml_map('fnrm',fnrm,fmt='(1pe21.14)')
!            call yaml_map('beta',beta,fmt='(1pe21.14)')
!            call yaml_map('beta_stretch',beta_stretch,fmt='(1pe21.14)')
!            call geometry_output(fmax,fnrm,fluct)
!         call yaml_mapping_close()
      end if

      etot    = etotp
      etotold = etot
      !copy outs_datatype
!      call copy_global_output(outs,outsIO)

      if(detot .gt. maxrise)then
!         if (iproc==0) write(16,'(a,i0,4(1x,e9.2))') &
         write(*,'(a,i0,4(1x,e9.2))') &
             "WARNING : Allowed energy to rise by more than maxrise: it,maxrise,detot,beta,1.d-1*betax ",&
             it,maxrise,detot,beta,1.d-1*betax
      endif


!      if (fnrm.le.fnrmtol) goto 1000
!      call convcheck(parini,fmax,fluct*runObj%inputs%frac_fluct,runObj%inputs%forcemax,icheck)
      call convcheck(parini,parini%nat,latvec_in,fcart_in,strten_in,parini%target_pressure_habohr,parini%paropt_geopt%strfact,fmax,fmax_at,fmax_lat,parini%paropt_geopt%fmaxtol,iexit)
!      if(icheck>5)then
      if(iexit==1)then
         goto 1000
      endif
!     if(imode==2)rxyz(:,:,nhist)=rxyz(:,:,nhist)+beta_stretch*fstretch(:,:,nhist) !has to be after convergence check,
                                                                       !otherwise energy will not match
                                                                       !the true energy of rxyz(:,:,nhist)

!      if(ncount_bigdft >= nit)then!no convergence within ncount_cluster_x energy evaluations
      if(it >= nit)then!no convergence within ncount_cluster_x energy evaluations
            goto 900  !sqnm will return to caller the energies and coordinates used/obtained from the last accepted iteration step
      endif
   
      if (cosangle.gt..200d0) then
         beta=beta*1.100d0
      else
         beta=max(beta*.850d0,0.1d0*betax)
      endif
      if (cosangle_lat.gt..200d0) then
         betalat=betalat*1.100d0
      else
         betalat=max(betalat*.850d0,1.d0*betalatx)
      endif
!      betalat_scale=betalat/beta
      betalat=betalat_scale*beta
   
!      if (debug.and.iproc==0) write(100,*) 'cosangle ',cosangle,beta
       write(*,*) 'cosangle ',cosangle,beta

!      call getSubSpaceEvecEval_minhocao('(SQNM)',iproc,verbosity,nat,nhist,nhistx,ndim,cutoffratio,lwork,work,rxyz,&
!                   &fxyz,aa,rr,ff,rrr,fff,eval,res,success)
      call getSubSpaceEvecEval_minhocao('(SQNM)',parini%verb,parini%nat,nhist,nhistx,ndim,cutoffratio,lwork,work,rxyz,&
                   &fxyz,aa,rr,ff,rrr,fff,eval,res,success)
      if(.not.success)stop 'subroutine minimizer_sqnm: no success in getSubSpaceEvecEval_minhocao.'


!Set precision if necessary
       if(parini%usewf_geopt) then
           getwfk=.true.
       else
           getwfk=.false.
       endif
       if(iprec==2.and.multiprec.and.fmax.lt.tolmxf_switch) then
           getwfk=.false.
           iprec=1
       endif
!Reset everything, recompute cell and stuff
         if((multiprec.and.it.ge.parini%paropt_geopt%nit/2).or.&
          &(fmax.lt.1.0d0*tolmxf_switch)) max_kpt=.true.
         if(fmax.lt.cellfix_switch.and..not.cellfix_done.and.(.not.(any(parini%fixlat).or.any(parini%fixat).or.confine.ge.1))) then
!Only perform the cell correction once, presumably close to the end of the optimization run
             if(cart_forces) then
                    call rxyz_cart2int(rxyz(:,parini%nat+1:parini%nat+3,nhist),pos_tmp,rxyz(:,1:parini%nat,nhist),parini%nat)
             call correct_latvec(rxyz(:,parini%nat+1:parini%nat+3,nhist),pos_tmp(:,:),parini%nat,parini%correctalg,latvec_io)
                    call rxyz_int2cart(rxyz(:,parini%nat+1:parini%nat+3,nhist),pos_tmp,rxyz(:,1:parini%nat,nhist),parini%nat)
             else
             call correct_latvec(rxyz(:,parini%nat+1:parini%nat+3,nhist),rxyz(:,1:parini%nat,nhist),parini%nat,parini%correctalg,latvec_io)
             endif
             cellfix_done=.true.
             if(latvec_io.ne.0) then
                max_kpt=.false.
                getwfk=.false.
                ka1=0;kb1=0;kc1=0
   !Reset all
                maxrise=    parini%paropt_geopt%maxrise!1.d-6!maxrise
                trustr=     parini%paropt_geopt%trustr!0.2d0!trustr
                displr=0.0d0
                displp=0.0d0
                icheck=0
                detot=0.0d0
                itswitch=0
                ndim=0
                nhist=0
                beta=betax
                maxd=1.0d0
             endif
         endif 
   enddo!end main loop

900 continue

   !if code gets here, it failed
!   if(debug.and.iproc==0) write(100,*) it,etot,fnrm
!   if(iproc==0) write(16,'(a,3(1x,i0))') &
   write(*,'(a,2(1x,i0))') &
       "WARNING: SQNM not converged: it,nit: ", &
!       "WARNING: SQNM not converged: it,ncount_bigdft,ncount_cluster_x: ", &
!       it,ncount_bigdft,runObj%inputs%ncount_cluster_x
       it,nit
!   stop "No convergence "
   fail=.true.
   goto 2000

1000 continue!converged successfully
   
!   if(iproc==0) write(16,'(2(a,1x,i0))') "SQNM converged at iteration ",it,". Needed bigdft calls: ",ncount_bigdft
   write(*,'((a,1x,i5))') "SQNM converged at iteration ",it
!   if(iproc==0)  call yaml_map('Iterations when SQNM converged',it)
!   write(*,*) 'Iterations when SQNM converged',it
   fail=.false.
   
!   etot=etotp
!   do iat=1,nat
!      do l=1,3
!         wpos(l,iat)= rxyz(l,iat,nhist)
!         fout(l,iat)= fxyz(l,iat,nhist)
!      enddo
!   enddo
2000 continue
!deallocations
!   call f_free(rxyz)
!   call f_free(rxyzOld)
!   call f_free(delta)
!   call f_free(rxyzraw)
!   call f_free(fxyz)
!   call f_free(fxyzraw)
!   call f_free(fstretch)
!   call f_free(aa)
!   call f_free(eval)
!   call f_free(res)
!   call f_free(rnorm)
!   call f_free(work)
!   call f_free(ff)
!   call f_free(rr)
!   call f_free(dd)
!   call f_free(fff)
!   call f_free(rrr)
!   call f_free(scpr)
!   call f_free(wold)
!   call f_free(rcov )   
!   call f_free(iconnect)
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
max_kpt=.false.
end subroutine
!subroutine minenergyandforces(parini,parres,iproc,nproc,eeval,imode,runObj,outs,nat,rat,rxyzraw,fat,fstretch,&
!           fxyzraw,epot,iconnect,nbond_,wold,alpha_stretch0,alpha_stretch)
subroutine minenergyandforces(parini,parres,eeval,imode,nat,rat,rxyzraw,fat,fstretch,&
           fxyzraw,epot,alpha_stretch0,alpha_stretch,&
           latvec_in,xred_in,etot_in,fcart_in,strten_in,iprec)
!    use module_base
!    use bigdft_run!module_types
    use module_sqn
!    use module_interfaces
    use mod_parini, only: typ_parini
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_parini), intent(inout):: parres
    !parameter
    integer, intent(in)           :: imode
    integer, intent(in)           :: nat
!    type(run_objects), intent(inout)       :: runObj
!    type(DFT_global_output), intent(inout) :: outs
!    integer, intent(in)           :: nbond_
!    integer, intent(in)           :: iconnect(2,nbond_)
    real(8),intent(inout)        :: rat(3,nat+3)
    real(8),intent(out)          :: rxyzraw(3,nat+3)
    real(8),intent(out)          :: fxyzraw(3,nat+3)
    real(8),intent(inout)        :: fat(3,nat+3)
    real(8),intent(out)          :: fstretch(3,nat+3)
!    real(8), intent(inout)       :: wold(nbond_)
    real(8), intent(in)          :: alpha_stretch0
    real(8), intent(inout)       :: alpha_stretch
    real(8), intent(inout)       :: epot
    logical, intent(in)          :: eeval
    !internal
    integer :: infocode
!my own
    real(8):: rxyz(3,nat+3)
    real(8):: fxyz(3,nat+3)
    real(8):: enthalpy
    real(8):: force_all(3,nat+3)
    real(8):: latvec0(3,3),latvec_in(3,3),xred_in(3,nat),etot_in,fcart_in(3,nat),strten_in(6) 
    logical:: getwfk
    integer:: lattdeg=1,iprec
    
    latvec0=0.d0
!    rxyzraw=rat
!    if(eeval)call energyandforces(nat,alat,rat,fat,fnoise,epot)
!    fxyzraw=fat
!    fstretch=0.0d0

!    call vcopy(3 * runObj%nat, rat(1,1), 1,rxyzraw(1,1), 1)
    rxyzraw(:,:)=rat(:,:)
    if(eeval)then
!        call vcopy(3 * runObj%nat, rat(1,1), 1,runObj%rxyz(1,1), 1)
        rxyz(:,:)=rat(:,:)
!        runObj%inputs%inputPsiId=1
!        call call_bigdft(runObj,outs,infocode)
          
         getwfk=.false.
!         call get_energyandforces_single(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,getwfk)
!!!!         call get_BFGS_forces_strainlatt(parini,rat,force_all,enthalpy,getwfk,iprec,latvec0,&
!!!!             &lattdeg,latvec_in,xred_in,etot_in,fcart_in,strten_in)
          call get_BFGS_forces_PR(parini,parres,rat,force_all,enthalpy,getwfk,iprec,latvec_in,xred_in,etot_in,fcart_in,strten_in)
    endif
!    call vcopy(3 * outs%fdim, outs%fxyz(1,1), 1, fat(1,1), 1)
    fat=force_all
!    call vcopy(3 * outs%fdim, fat(1,1), 1,fxyzraw(1,1), 1)
    fxyzraw(:,:)=fat(:,:)
!    epot=outs%energy
    epot=enthalpy
    fstretch=0.0d0
!    if(imode==2)then
!        call projectbond(nat,nbond_,rat,fat,fstretch,iconnect,&
!             wold,alpha_stretch0,alpha_stretch)
!    endif

end subroutine minenergyandforces
!subroutine give_rcov_sqnm(iproc,atoms,nat,rcov)
!!subroutine give_rcov_sqnm(iproc,atomnames,nat,ntypat,rcov)
!!!  use module_base, only: 8
!!!  use module_types
!!!  use yaml_output
!!  implicit none
!!  !Arguments
!!  integer, intent(in) :: iproc,nat
!!!  type(atoms_data), intent(in) :: atoms
!!  real(8), intent(out) :: rcov(nat)
!!  !Local variables
!!  integer :: iat
!!!my own
!!  integer:: iatype(nat),ntypat
!!  character(2):: atomnames(ntypat)
!!
!!
!!  do iat=1,nat
!!     if (trim(atomnames(iatype(iat)))=='H') then
!!        rcov(iat)=0.75d0
!!!        rcov(iat)=0.75d0*0.529177211d0
!!     else if (trim(atomnames(iatype(iat)))=='LJ')then
!!        rcov(iat)=0.56d0
!!     else if (trim(atomnames(iatype(iat)))=='He')then
!!        rcov(iat)=0.75d0
!!     else if (trim(atomnames(iatype(iat)))=='Li')then
!!        rcov(iat)=3.40d0
!!     else if (trim(atomnames(iatype(iat)))=='Be')then
!!        rcov(iat)=2.30d0
!!     else if (trim(atomnames(iatype(iat)))=='B' )then
!!        rcov(iat)=1.55d0
!!     else if (trim(atomnames(iatype(iat)))=='C' )then
!!        rcov(iat)=1.45d0
!!!        rcov(iat)=1.45d0*0.529177211d0
!!     else if (trim(atomnames(iatype(iat)))=='N' )then
!!        rcov(iat)=1.42d0
!!!        rcov(iat)=1.42d0*0.529177211d0
!!     else if (trim(atomnames(iatype(iat)))=='O' )then
!!        rcov(iat)=1.38d0
!!!        rcov(iat)=1.38d0*0.529177211d0
!!     else if (trim(atomnames(iatype(iat)))=='F' )then
!!        rcov(iat)=1.35d0
!!     else if (trim(atomnames(iatype(iat)))=='Ne')then
!!        rcov(iat)=1.35d0
!!     else if (trim(atomnames(iatype(iat)))=='Na')then
!!        rcov(iat)=3.40d0
!!     else if (trim(atomnames(iatype(iat)))=='Mg')then
!!        rcov(iat)=2.65d0
!!     else if (trim(atomnames(iatype(iat)))=='Al')then
!!        rcov(iat)=2.23d0
!!     else if (trim(atomnames(iatype(iat)))=='Si')then
!!        rcov(iat)=2.09d0
!!     else if (trim(atomnames(iatype(iat)))=='P' )then
!!        rcov(iat)=2.00d0
!!     else if (trim(atomnames(iatype(iat)))=='S' )then
!!        rcov(iat)=1.92d0
!!     else if (trim(atomnames(iatype(iat)))=='Cl')then
!!        rcov(iat)=1.87d0
!!     else if (trim(atomnames(iatype(iat)))=='Ar')then
!!        rcov(iat)=1.80d0
!!     else if (trim(atomnames(iatype(iat)))=='K' )then
!!        rcov(iat)=4.00d0
!!     else if (trim(atomnames(iatype(iat)))=='Ca')then
!!        rcov(iat)=3.00d0
!!     else if (trim(atomnames(iatype(iat)))=='Sc')then
!!        rcov(iat)=2.70d0
!!     else if (trim(atomnames(iatype(iat)))=='Ti')then
!!        rcov(iat)=2.70d0
!!     else if (trim(atomnames(iatype(iat)))=='V' )then
!!        rcov(iat)=2.60d0
!!     else if (trim(atomnames(iatype(iat)))=='Cr')then
!!        rcov(iat)=2.60d0
!!     else if (trim(atomnames(iatype(iat)))=='Mn')then
!!        rcov(iat)=2.50d0
!!     else if (trim(atomnames(iatype(iat)))=='Fe')then
!!        rcov(iat)=2.50d0
!!     else if (trim(atomnames(iatype(iat)))=='Co')then
!!        rcov(iat)=2.40d0
!!     else if(trim(atomnames(iatype(iat)))=='Ni')then
!!        rcov(iat)=2.30d0
!!     else if (trim(atomnames(iatype(iat)))=='Cu')then
!!        rcov(iat)=2.30d0
!!     else if (trim(atomnames(iatype(iat)))=='Zn')then
!!        rcov(iat)=2.30d0
!!     else if (trim(atomnames(iatype(iat)))=='Ga')then
!!        rcov(iat)=2.10d0
!!     else if (trim(atomnames(iatype(iat)))=='Ge')then
!!        rcov(iat)=2.40d0
!!     else if (trim(atomnames(iatype(iat)))=='As')then
!!        rcov(iat)=2.30d0
!!     else if (trim(atomnames(iatype(iat)))=='Se')then
!!        rcov(iat)=2.30d0
!!     else if (trim(atomnames(iatype(iat)))=='Br')then
!!        rcov(iat)=2.20d0
!!     else if (trim(atomnames(iatype(iat)))=='Kr')then
!!        rcov(iat)=2.20d0
!!     else if (trim(atomnames(iatype(iat)))=='Rb')then
!!        rcov(iat)=4.50d0
!!     else if (trim(atomnames(iatype(iat)))=='Sr')then
!!        rcov(iat)=3.30d0
!!     else if (trim(atomnames(iatype(iat)))=='Y' )then
!!        rcov(iat)=3.30d0
!!     else if (trim(atomnames(iatype(iat)))=='Zr')then
!!        rcov(iat)=3.00d0
!!     else if (trim(atomnames(iatype(iat)))=='Nb')then
!!        rcov(iat)=2.92d0
!!     else if (trim(atomnames(iatype(iat)))=='Mo')then
!!        rcov(iat)=2.83d0
!!     else if (trim(atomnames(iatype(iat)))=='Tc')then
!!        rcov(iat)=2.75d0
!!     else if (trim(atomnames(iatype(iat)))=='Ru')then
!!        rcov(iat)=2.67d0
!!     else if (trim(atomnames(iatype(iat)))=='Rh')then
!!        rcov(iat)=2.58d0
!!     else if (trim(atomnames(iatype(iat)))=='Pd')then
!!        rcov(iat)=2.50d0
!!     else if (trim(atomnames(iatype(iat)))=='Ag')then
!!        rcov(iat)=2.50d0
!!     else if (trim(atomnames(iatype(iat)))=='Cd')then
!!        rcov(iat)=2.50d0
!!     else if (trim(atomnames(iatype(iat)))=='In')then
!!        rcov(iat)=2.30d0
!!     else if (trim(atomnames(iatype(iat)))=='Sn')then
!!        rcov(iat)=2.66d0
!!     else if (trim(atomnames(iatype(iat)))=='Sb')then
!!        rcov(iat)=2.66d0
!!     else if (trim(atomnames(iatype(iat)))=='Te')then
!!        rcov(iat)=2.53d0
!!     else if (trim(atomnames(iatype(iat)))=='I' )then
!!        rcov(iat)=2.50d0
!!     else if (trim(atomnames(iatype(iat)))=='Xe')then
!!        rcov(iat)=2.50d0
!!     else if (trim(atomnames(iatype(iat)))=='Cs')then
!!        rcov(iat)=4.50d0
!!     else if (trim(atomnames(iatype(iat)))=='Ba')then
!!        rcov(iat)=4.00d0
!!     else if (trim(atomnames(iatype(iat)))=='La')then
!!        rcov(iat)=3.50d0
!!     else if (trim(atomnames(iatype(iat)))=='Ce')then
!!        rcov(iat)=3.50d0
!!     else if (trim(atomnames(iatype(iat)))=='Pr')then
!!        rcov(iat)=3.44d0
!!     else if (trim(atomnames(iatype(iat)))=='Nd')then
!!        rcov(iat)=3.38d0
!!     else if (trim(atomnames(iatype(iat)))=='Pm')then
!!        rcov(iat)=3.33d0
!!     else if (trim(atomnames(iatype(iat)))=='Sm')then
!!        rcov(iat)=3.27d0
!!     else if (trim(atomnames(iatype(iat)))=='Eu')then
!!        rcov(iat)=3.21d0
!!     else if (trim(atomnames(iatype(iat)))=='Gd')then
!!        rcov(iat)=3.15d0
!!     else if (trim(atomnames(iatype(iat)))=='Td')then
!!        rcov(iat)=3.09d0
!!     else if (trim(atomnames(iatype(iat)))=='Dy')then
!!        rcov(iat)=3.03d0
!!     else if (trim(atomnames(iatype(iat)))=='Ho')then
!!        rcov(iat)=2.97d0
!!     else if (trim(atomnames(iatype(iat)))=='Er')then
!!        rcov(iat)=2.92d0
!!     else if (trim(atomnames(iatype(iat)))=='Tm')then
!!        rcov(iat)=2.92d0
!!     else if (trim(atomnames(iatype(iat)))=='Yb')then
!!        rcov(iat)=2.80d0
!!     else if (trim(atomnames(iatype(iat)))=='Lu')then
!!        rcov(iat)=2.80d0
!!     else if (trim(atomnames(iatype(iat)))=='Hf')then
!!        rcov(iat)=2.90d0
!!     else if (trim(atomnames(iatype(iat)))=='Ta')then
!!        rcov(iat)=2.70d0
!!     else if (trim(atomnames(iatype(iat)))=='W' )then
!!        rcov(iat)=2.60d0
!!     else if (trim(atomnames(iatype(iat)))=='Re')then
!!        rcov(iat)=2.60d0
!!     else if (trim(atomnames(iatype(iat)))=='Os')then
!!        rcov(iat)=2.50d0
!!     else if (trim(atomnames(iatype(iat)))=='Ir')then
!!        rcov(iat)=2.50d0
!!     else if (trim(atomnames(iatype(iat)))=='Pt')then
!!        rcov(iat)=2.60d0
!!     else if (trim(atomnames(iatype(iat)))=='Au')then
!!        rcov(iat)=2.70d0
!!     else if (trim(atomnames(iatype(iat)))=='Hg')then
!!        rcov(iat)=2.80d0
!!     else if (trim(atomnames(iatype(iat)))=='Tl')then
!!        rcov(iat)=2.50d0
!!     else if (trim(atomnames(iatype(iat)))=='Pb')then
!!        rcov(iat)=3.30d0
!!     else if (trim(atomnames(iatype(iat)))=='Bi')then
!!        rcov(iat)=2.90d0
!!     else if (trim(atomnames(iatype(iat)))=='Po')then
!!        rcov(iat)=2.80d0
!!     else if (trim(atomnames(iatype(iat)))=='At')then
!!        rcov(iat)=2.60d0
!!     else if (trim(atomnames(iatype(iat)))=='Rn')then
!!        rcov(iat)=2.60d0
!!     else
!!        call yaml_comment('(SQNM) no covalent radius stored for this atomtype '&
!!             //trim(atomnames(iatype(iat))))
!!        stop
!!     endif
!!     if (iproc == 0) then
!!        call yaml_map('(SQNM) RCOV:'//trim(atomnames&
!!                 (iatype(iat))),rcov(iat))
!!     endif
!!  enddo
!!end subroutine give_rcov_sqnm

subroutine sqnm_invhess(nat,h,metric,hessinv)
implicit none
integer:: nat,info,i,j,k
real(8):: metric(3*(nat+3),3*(nat+3)),hessinv(3*(nat+3),3*(nat+3))
real(8):: h(3,3),hinv(3,3),g(3,3),ginv(3,3)
real(8):: garbage,omega
real(8):: hessinv_at(3*nat,3*nat),hessinv_lat(9,9)
real(8):: eval(3*(nat+3)),eval_at(3*nat),eval_lat(9)
real(8),allocatable:: work(:)
integer:: lwork
   lwork=3*(3*(nat+3))
   allocate(work(lwork))
   call invmat(3, h, hinv, omega)
      omega = abs(omega) 
! ... generate metric to work with scaled ionic coordinates
   g = MATMUL(TRANSPOSE(h),h)
   call invmat(3,g,ginv,garbage)
   metric = 0.d0
   FORALL ( k=0:nat-1,   i=1:3, j=1:3 )     metric(i+3*k,j+3*k) = g(i,j)
   FORALL ( k=nat:nat+2, i=1:3, j=1:3 )     metric(i+3*k,j+3*k) = 0.04 * omega * ginv(i,j)
   call invmat(3*(nat+3),metric,hessinv,garbage)
   call invmat(3*nat,metric(1:3*nat,1:3*nat),hessinv_at,garbage)
   call invmat(9,metric(3*(nat)+1:3*(nat+3),3*(nat)+1:3*(nat+3)),hessinv_lat,garbage)
!Diagonalize and get eigenvalues
   call dsyev('v','l', 3*(nat+3), hessinv, 3*(nat+3), eval, work, lwork, info)
        do i=1,3*(nat+3)
            write(*,*) ' All eigenvalues: ',i,eval(i)
        enddo
   call dsyev('v','l', 3*(nat), hessinv_at, 3*(nat), eval_at, work, lwork, info)
        do i=1,3*(nat)
            write(*,*) ' Atomic eigenvalues: ',i,eval_at(i)
        enddo
   call dsyev('v','l', 9, hessinv_lat, 9, eval_lat, work, lwork, info)
        do i=1,9
            write(*,*) ' Lattice eigenvalues: ',i,eval_lat(i)
        enddo
end subroutine
