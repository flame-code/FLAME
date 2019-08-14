!  type, public :: state_properties
!     real(gp) :: energy, fnoise, pressure      !< Total energy, noise over forces and pressure
!     type(energy_terms) :: energs              !< All energy terms
!     integer :: fdim                           !< Dimension of allocated forces (second dimension)
!     real(gp), dimension(:,:), pointer :: fxyz !< Atomic forces
!     real(gp), dimension(6) :: strten          !< Stress Tensor
!  end type state_properties

!*****************************************************************************************
subroutine sqnm(parini,atoms,paropt,count_sqnm,fail)
   !use module_base
   !use bigdft_run!module_types
   !use yaml_output
   !use module_sqn, only: modify_gradient, getSubSpaceEvecEval, findbonds
   use mod_parini, only: typ_parini
   use mod_atoms, only: typ_atoms, set_rat, get_rat
   use mod_opt, only: typ_paropt, frmt_base
   use mod_processors, only: iproc, nproc
   use dynamic_memory
   use yaml_output
   implicit none
   !parameter
   !integer, intent(in)                    :: nproc
   !integer, intent(in)                    :: iproc
   !integer, intent(in)                    :: verbosity
   !type(run_objects), intent(inout)       :: runObj
   !type(state_properties), intent(inout) :: outsIO
   type(typ_parini), intent(in):: parini
   type(typ_atoms), intent(inout):: atoms
   type(typ_paropt), intent(inout):: paropt
   real(8), intent(inout):: count_sqnm
   logical, intent(out):: fail
   !local variables
   character(len=*), parameter :: subname='sqnm'
   character(41), parameter:: frmt='('//frmt_base//',i5)'
   integer:: verbosity
   integer :: infocode,info !< variables containing state codes
   integer :: nhistx !< maximum history length
   integer :: nhist  !< actual history length
   integer :: ndim   !< dimension of significant subspace
   integer :: nit    !< maximum number of iterations
   integer :: nat    !< number of atoms
   integer :: istat,iall
   integer :: lwork
   integer :: it,i,iat,l,j,idim,jdim,ihist,icheck !<counter variables
   integer :: itswitch
   integer :: imode
   integer :: nbond
   !type(state_properties) :: outs
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
   real(8) :: fnrm
   real(8) :: fmax
   real(8) :: fluct
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
   integer :: idxtmp
   integer,  allocatable, dimension(:,:)   :: iconnect
   integer,  allocatable, dimension(:)     :: idx!index array for keeping track of history
   real(8), allocatable, dimension(:,:,:) :: rxyz
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
   integer :: ifail
   logical:: biomode
   !functions
   real(8) :: dnrm2, DDOT
   real(8) :: fnoise=0.d0
   real(8) :: tt1, tt2, tt3
   !type(f_tree) :: f_info
   call f_routine(id='sqnm')
   call yaml_sequence_open('SQNM optimization iterations')
   call yaml_sequence(advance='no')
   nat=atoms%nat
    call cal_potential_forces(parini,atoms)

   !f_info=f_tree_new()

   !set parameters
   if(paropt%nit<0) then
       nit=1000
   else
    nit=paropt%nit !runObj%inputs%ncount_cluster_x
   endif
   betax=paropt%alphax !runObj%inputs%betax
   nhistx=paropt%nhist !runObj%inputs%nhistx
   maxrise=1.d-4 !runObj%inputs%maxrise
   cutoffRatio=1.d-4 !runObj%inputs%cutoffratio
   steepthresh=1.d-1*sqrt(real(nat,8)) !runObj%inputs%steepthresh
   trustr=0.2d0 !runObj%inputs%trustr
   verbosity=parini%iverbose
   nbond=1
   ifail=0
   biomode=.false. !later one keyword must be added in input.ini for biomolecules, until must be false
   if(biomode) imode=2

   if (iproc==0 .and. verbosity>0) then
         !call yaml_mapping_open('Geometry parameters')
         !call yaml_map('Geometry Method','GEOPT_SQNM')
         !call yaml_map('nhistx',nhistx)
         !call yaml_map('biomode',runObj%inputs%biomode)
         !call yaml_map('betax', betax,fmt='(1pe21.14)')
         !call yaml_map('beta_stretchx', runObj%inputs%beta_stretchx,fmt='(1pe21.14)')
         !call yaml_map('maxrise', maxrise,fmt='(1pe21.14)')
         !call yaml_map('cutoffRatio', cutoffRatio,fmt='(1pe21.14)')
         !call yaml_map('steepthresh', steepthresh,fmt='(1pe21.14)')
         !call yaml_map('trustr', trustr,fmt='(1pe21.14)')
         !call yaml_mapping_close()
   end if

   !init variables
   debug=.false.
   fail=.true.
   displr=0.d0
   displp=0.d0
   fluct=paropt%anoise
   icheck=0
   detot=0.d0
   itswitch=0
   ndim=0
   nhist=0
   beta=betax*0.2d0
   beta_stretch=0.4d0 !runObj%inputs%beta_stretchx
   maxd=0.d0

   ! allocate arrays
   !lwork=1000+10*nat**2
   lwork=10*nhistx**2
   rxyz = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='rxyz')
   idx = f_malloc((/ 0.to.nhistx /),id='idx')
   fxyz = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='fxyz')
   fxyzraw = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='fxyzraw')
   fstretch = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='fstretch')
   rxyzOld = f_malloc((/ 1.to.3, 1.to.nat/),id='rxyzOld')
   delta = f_malloc((/ 1.to.3, 1.to.nat/),id='delta')
   aa = f_malloc((/ nhistx, nhistx /),id='aa')
   eval = f_malloc(nhistx,id='eval')
   res = f_malloc(nhistx,id='res')
   rnorm = f_malloc(nhistx,id='rnorm')
   work = f_malloc(lwork,id='work')
   ff = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='ff')
   rr = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='rr')
   dd = f_malloc((/ 3, nat /),id='dd')
   fff = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='fff')
   rrr = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='rrr')
   scpr = f_malloc(nhistx,id='scpr')
   rcov     = f_malloc((/ 1.to.nat/),id='rcov')
   iconnect = f_malloc((/ 1.to.2, 1.to.1000/),id='iconnect')
   if(biomode)then
        call give_rcov_sqnm(iproc,atoms,rcov)
        call findbonds('(SQNM)',iproc,10,atoms,rcov,nbond,iconnect)
   endif 
   wold = f_malloc((/ 1.to.nbond/),id='wold')
   wold =0.d0
   do i=0,nhistx
    idx(i)=i
   enddo


   !call init_state_properties(outs, runObj%atoms%astruct%nat)
   !copy outs_datatype
   !call copy_state_properties(outsIO,outs)



!!!!!!   call energyandforces(nat,rxyz(1,1,0),fxyz(1,1,0),etot)
!!  not necessary, bigdft_state allready called outside
!   call bigdft_state(runObj,outs,nproc,iproc,infocode)
!   count_sqnm=count_sqnm+1

!! copy to internal variables
   !call vcopy(3*runObj%atoms%astruct%nat, runObj%atoms%astruct%rxyz(1,1), 1,rxyz(1,1,idx(0)), 1)
   !call vcopy(3*runObj%atoms%astruct%nat, runObj%atoms%astruct%rxyz(1,1), 1,rxyzOld(1,1), 1)
   !call vcopy(3*outs%fdim, outs%fxyz(1,1), 1, fxyz(1,1,idx(0)), 1)
   !rxyz(1:3,1:nat,idx(0))=atoms%rat(1:3,1:nat)
   call get_rat(atoms,rxyz(1,1,idx(0)))
   !rxyzOld(1:3,1:nat)=atoms%rat(1:3,1:nat)
   call get_rat(atoms,rxyzOld)
   fxyz(1:3,1:nat,idx(0))=atoms%fat(1:3,1:nat)
   !----------------------------------------------------------------------------
   !The following lines added by Alireza to avoid moving fixed atoms.
    do iat=1,atoms%nat
        if(.not. atoms%bemoved(1,iat)) fxyz(1,iat,idx(0))=0.d0
        if(.not. atoms%bemoved(2,iat)) fxyz(2,iat,idx(0))=0.d0
        if(.not. atoms%bemoved(3,iat)) fxyz(3,iat,idx(0))=0.d0
    enddo
   !----------------------------------------------------------------------------

   etot=atoms%epot !outs%energy

   call minenergyandforces_alborz(parini,iproc,nproc,.false.,imode,atoms,nat,rxyz(1,1,idx(0)),&
       fxyz(1,1,idx(0)),fstretch(1,1,idx(0)),fxyzraw(1,1,idx(0)),&
       etot,iconnect,nbond,wold,beta_stretchx,beta_stretch,infocode)
   if(imode==2)rxyz(:,:,idx(0))=rxyz(:,:,idx(0))+beta_stretch*fstretch(:,:,idx(0))

   !call fnrmandforcemax(fxyzraw(1,1,idx(0)),fnrm,fmax,nat)
   call calmaxforcecomponent(3*nat,fxyzraw(1,1,idx(0)),fmax)
   call calnorm(3*nat,fxyzraw(1,1,idx(0)),fnrm)
   fnrm=sqrt(fnrm)
   if (fmax < 3.d-1) fluct=paropt%anoise !call updatefluctsum(fnoise,fluct)

   etotold=etot
   etotp=etot

   if (iproc==0.and.verbosity > 0) then
       !avoid space for leading sign (numbers are positive, anyway)
       write(cdmy8,'(es8.1)')abs(maxd)
       write(cdmy12_1,'(es12.5)')abs(displr)
       write(cdmy12_2,'(es12.5)')abs(displp)
       write(cdmy9,'(es9.2)')abs(beta)

       write(16,'(i5,1x,i5,2x,a10,2x,1es21.14,2x,es9.2,es11.3,es10.2,2x,a6,a8,1x,a4,i3.3,1x,a5,a7,2(1x,a6,a11))') &
       nint(count_sqnm),0,'GEOPT_SQNM',etotp,detot,fmax,fnrm, &
       'beta=',trim(adjustl(cdmy9)),'dim=',ndim,'maxd=',&
       trim(adjustl(cdmy8)),'dsplr=',trim(adjustl(cdmy12_1)),&
       'dsplp=',trim(adjustl(cdmy12_2))

       !write(16,'(i5,1x,i5,2x,a10,2x,1es21.14,2x,es9.2,es11.3,es10.2,2x,a6,a8,1x,a4,i3.3,1x,a5,a7,2(1x,a6,a11))') &
       !count_sqnm,0,'GEOPT_SQNM',etotp,detot,fmax,fnrm, &
       !'beta=',trim(adjustl(cdmy9)),'dim=',ndim,'maxd=',&
       !trim(adjustl(cdmy8)),'dsplr=',trim(adjustl(cdmy12_1)),&
       !'dsplp=',trim(adjustl(cdmy12_2))
       !call f_utils_flush(16)
   endif

   do it=1,nit!start main loop
!  do it=1,nit-1!start main loop (nit-1 if first bigdft call is NOT done outside, but inside this subroutine)
      if (debug.and.iproc==0) write(100,*) 'it:',it,etot,fnrm,itswitch
      nhist=nhist+1

      if (fnrm.gt.steepthresh .or. it.le.itswitch ) then
         ndim=0
         steep=.true.
         if (it.gt.itswitch) itswitch=it+nhistx
         if (debug.and.iproc==0) write(100,*) "STEEP"
      else
         steep=.false.
      endif
      write(21,*) it,steep

      ! cycle the index array
      if (nhist.gt.nhistx) then
         nhist=nhistx
         idxtmp=idx(0)
         do ihist=0,nhist-1
            idx(ihist)=idx(ihist+1)
         enddo
         idx(nhist)=idxtmp
      endif
   
      ! decompose gradient
500 continue
    call modify_gradient(nat,ndim,rrr(1,1,1),eval(1),res(1),fxyz(1,1,idx(nhist-1)),beta,dd(1,1))
   
      tt=0.d0
      dt=0.d0
      maxd=-huge(1.d0)
      do iat=1,nat
         dt=dd(1,iat)**2+dd(2,iat)**2+dd(3,iat)**2
         tt=tt+dt
         maxd=max(maxd,dt)
      enddo
      tt=sqrt(tt)
      maxd=sqrt(maxd)
   
      !trust radius approach: avoids too large steps due to large forces
      !only used when in steepest decent mode
      if(maxd>trustr .and. steep)then
         if(debug.and.iproc==0)write(100,'(a,1x,es24.17,1x,i0)')'step too large',maxd,it
         if(iproc==0)then
            write(16,'(a,2(1x,es9.2))')'WARNING GEOPT_SQNM: step too large: maxd, trustradius ',maxd,trustr
            !call f_utils_flush(16)
         endif
         scl=0.5d0*trustr/maxd
         dd=dd*scl
         tt=tt*scl
         maxd=maxd*scl
      endif
!      displr=displr+tt
   
      !update positions
      do iat=1,nat
         rxyz(1,iat,idx(nhist))=rxyz(1,iat,idx(nhist-1))-dd(1,iat)
         rxyz(2,iat,idx(nhist))=rxyz(2,iat,idx(nhist-1))-dd(2,iat)
         rxyz(3,iat,idx(nhist))=rxyz(3,iat,idx(nhist-1))-dd(3,iat)
      enddo
   
      delta=rxyz(:,:,idx(nhist))-rxyzOld
      displr=displr+dnrm2(3*nat,delta(1,1),1)
      !runObj%inputs%inputPsiId=1
      !call bigdft_set_input_policy(INPUT_POLICY_MEMORY, runObj)
      call minenergyandforces_alborz(parini,iproc,nproc,.true.,imode,atoms,nat,rxyz(1,1,idx(nhist)),&
                             fxyz(1,1,idx(nhist)),fstretch(1,1,idx(nhist)),fxyzraw(1,1,idx(nhist)),&
                             etotp,iconnect,nbond,wold,beta_stretchx,beta_stretch,infocode)
      detot=etotp-etotold
      count_sqnm=count_sqnm+1


      !call fnrmandforcemax(fxyzraw(1,1,idx(nhist)),fnrm,fmax,nat)
      call calmaxforcecomponent(3*nat,fxyzraw(1,1,idx(nhist)),fmax)
      call calnorm(3*nat,fxyzraw(1,1,idx(nhist)),fnrm)
      fnrm=sqrt(fnrm)

      if (iproc == 0) then
         write(fn4,'(i4.4)') nint(count_sqnm)
!         write(comment,'(a,1pe10.3)')'SQNM:fnrm= ',fnrm
         if (detot.gt.maxrise .and. beta > 1.d-1*betax) then !
            write(comment,'(a,1pe10.3)')'R SQNM:fnrm= ',fnrm
         else
            write(comment,'(a,1pe10.3)')'A SQNM:fnrm= ',fnrm
         endif

         !call bigdft_write_atomic_file(runObj,outs,'posout_'//fn4,trim(comment))
      endif
      if(infocode==0)then
        ifail=0
      else
        ifail=ifail+1
      endif
      if ((infocode==0 .or. ifail>20).and.(fmax < 3.d-1)) fluct=paropt%anoise !call updatefluctsum(fnoise,fluct)
      tt1=DDOT(3*nat,fxyz(1,1,idx(nhist)),1,dd(1,1),1)
      tt2=DDOT(3*nat,fxyz(1,1,idx(nhist)),1,fxyz(1,1,idx(nhist)),1)
      tt3=DDOT(3*nat,dd(1,1),1,dd(1,1),1)
      !write(*,*) tt1,tt2,tt3
      !stop
      cosangle=-tt1/sqrt(tt2*tt3)

      if (detot.gt.maxrise .and. beta > 1.d-1*betax) then !
         if (debug.and.iproc==0) write(100,'(a,i0,1x,e9.2)') "WARN: it,detot", it,detot
         if (debug.and.iproc==0) then
             write(16,'(a,i0,4(1x,e9.2))') &
             "WARNING GEOPT_SQNM: Prevent energy to rise by more than maxrise: it,maxrise,detot,beta,1.e-1*betax ",&
             it,maxrise,detot,beta,1.d-1*betax
             !call f_utils_flush(16)
         endif
         if (iproc==0.and.verbosity > 0) then
            !avoid space for leading sign (numbers are positive, anyway)
            write(cdmy8,'(es8.1)')abs(maxd)
            write(cdmy12_1,'(es12.5)')abs(displr)
            write(cdmy12_2,'(es12.5)')abs(displp)
            write(cdmy9,'(es9.2)')abs(beta)

   
            write(16,'(i5,1x,i5,2x,a10,2x,1es21.14,2x,es9.2,es11.3,es10.2,2x,a6,a8,1x,a4,i3.3,1x,a5,a7,2(1x,a6,a11))') &
             nint(count_sqnm),it,'GEOPT_SQNM',etotp,detot,fmax,fnrm, &
             'beta=',trim(adjustl(cdmy9)),'dim=',ndim,&
             'maxd=',trim(adjustl(cdmy8)),'dsplr=',trim(adjustl(cdmy12_1)),&
             'dsplp=',trim(adjustl(cdmy12_2))
            !call f_utils_flush(16)

            !call f_tree_push(f_info//'ndim'        ,yaml_toa(ndim))
            !call f_tree_push(f_info//'etot'        ,yaml_toa(etotp,fmt='(1pe21.14)'))
            !call f_tree_push(f_info//'detot'       ,yaml_toa(detot,fmt='(1pe21.14)'))
            !call f_tree_push(f_info//'beta'        ,yaml_toa(beta,fmt='(1pe21.14)'))
            !call f_tree_push(f_info//'beta_stretch',yaml_toa(beta_stretch,fmt='(1pe21.14)'))
            !call geometry_output('GEOPT_SQNM',nint(count_sqnm),it,fmax,fnrm,fluct,f_info)
            !write(*,frmt) 'MIN: ',iproc,it,etotp,detot,fnrm,fmax,ndim
            if(it>1) then
                call yaml_sequence(advance='no')
            endif
            call yaml_mapping_open('SQNM',flow=.true.)
            call yaml_map('iproc',iproc,fmt='(i3.3)')
            call yaml_map('iter',it,fmt='(i5)')
            call yaml_map('epot',etotp,fmt='(es20.12)')
            call yaml_map('de',detot,fmt='(es11.3)')
            call yaml_map('fmax',fmax,fmt='(es12.5)')
            call yaml_map('fnrm',fnrm,fmt='(es12.5)')
            call yaml_map('ndim',ndim,fmt='(i4)')
            call yaml_mapping_close()

         end if
    
         if(nint(count_sqnm) >= nit)then!no convergence within ncount_cluster_x energy evaluations
            !following copy of rxyz(1,1,nhist-1) to runObj is necessary for returning to the caller
            !the energies and coordinates used/obtained from/in the last ACCEPTED iteration step
            !(otherwise coordinates of last call to bigdft_state would be returned)
            !call vcopy(3 * runObj%atoms%astruct%nat, rxyz(1,1,idx(nhist-1)), 1,runObj%atoms%astruct%rxyz(1,1), 1)
            !atoms%rat(1:3,1:nat)=rxyz(1:3,1:nat,idx(nhist-1))
            call set_rat(atoms,rxyz(1,1,idx(nhist-1)),setall=.true.)
            goto 900  !sqnm will return to caller the energies and coordinates used/obtained from the last ACCEPTED iteration step
         endif

         beta=.5d0*beta
         if (debug.and.iproc==0) write(100,'(a,1x,e9.2)') 'WARNING GEOPT_SQNM: beta reset ',beta
         ndim=0
         wold=0.d0
         if(.not.steep)then
            do iat=1,nat
               rxyz(1,iat,0)=rxyz(1,iat,idx(nhist-1))
               rxyz(2,iat,0)=rxyz(2,iat,idx(nhist-1))
               rxyz(3,iat,0)=rxyz(3,iat,idx(nhist-1))
   
               fxyz(1,iat,0)=fxyz(1,iat,idx(nhist-1))
               fxyz(2,iat,0)=fxyz(2,iat,idx(nhist-1))
               fxyz(3,iat,0)=fxyz(3,iat,idx(nhist-1))
               fxyzraw(1,iat,0)=fxyzraw(1,iat,idx(nhist-1))
               fxyzraw(2,iat,0)=fxyzraw(2,iat,idx(nhist-1))
               fxyzraw(3,iat,0)=fxyzraw(3,iat,idx(nhist-1))
            enddo
            nhist=1
            do i=0,nhistx
                idx(i)=i
            enddo
         endif
         goto  500
      endif

      delta=rxyz(:,:,idx(nhist))-rxyzOld
      displp=displp+dnrm2(3*nat,delta(1,1),1)
      rxyzOld=rxyz(:,:,idx(nhist))
!      displp=displp+tt
      if (iproc==0.and.verbosity > 0) then
         !avoid space for leading sign (numbers are positive, anyway)
         write(cdmy8,'(es8.1)')abs(maxd)
         write(cdmy12_1,'(es12.5)')abs(displr)
         write(cdmy12_2,'(es12.5)')abs(displp)
         write(cdmy9,'(es9.2)')abs(beta)


         write(16,'(i5,1x,i5,2x,a10,2x,1es21.14,2x,es9.2,es11.3,es10.2,2x,a6,a8,1x,a4,i3.3,1x,a5,a7,2(1x,a6,a11))') &
          nint(count_sqnm),it,'GEOPT_SQNM',etotp,detot,fmax,fnrm, &
          'beta=',trim(adjustl(cdmy9)),'dim=',ndim,'maxd=',&
          trim(adjustl(cdmy8)),'dsplr=',trim(adjustl(cdmy12_1)),&
          'dsplp=',trim(adjustl(cdmy12_2))
         !call f_utils_flush(16)

         !call f_tree_push(f_info//'ndim'        ,yaml_toa(ndim))
         !call f_tree_push(f_info//'etot'        ,yaml_toa(etotp,fmt='(1pe21.14)'))
         !call f_tree_push(f_info//'detot'       ,yaml_toa(detot,fmt='(1pe21.14)'))
         !call f_tree_push(f_info//'beta'        ,yaml_toa(beta,fmt='(1pe21.14)'))
         !call f_tree_push(f_info//'beta_stretch',yaml_toa(beta_stretch,fmt='(1pe21.14)'))
         !call geometry_output('GEOPT_SQNM',nint(count_sqnm),it,fmax,fnrm,fluct,f_info)
            !write(*,frmt) 'MIN: ',iproc,it,etotp,detot,fnrm,fmax,ndim
            if(it>1) then
                call yaml_sequence(advance='no')
            endif
            call yaml_mapping_open('SQNM',flow=.true.)
            !call yaml_map('iproc',iproc,fmt='(i3.3)')
            call yaml_map('iter',it,fmt='(i5)')
            call yaml_map('epot',etotp,fmt='(es20.12)')
            call yaml_map('de',detot,fmt='(es9.1)')
            call yaml_map('fmax',fmax,fmt='(es10.3)')
            call yaml_map('fnrm',fnrm,fmt='(es10.3)')
            call yaml_map('ndim',ndim,fmt='(i4)')
            call yaml_mapping_close()

!!$         call yaml_mapping_open('Geometry')
!!$            call yaml_map('Ncount_BigDFT',ncount_bigdft)
!!$            call yaml_map('Geometry step',it)
!!$            call yaml_map('Geometry Method','GEOPT_SQNM')
!!$            call yaml_map('ndim',ndim)
!!$            call yaml_map('etot', etotp,fmt='(1pe21.14)')
!!$            call yaml_map('detot',detot,fmt='(1pe21.14)')
!!$            call yaml_map('fmax',fmax,fmt='(1pe21.14)')
!!$            call yaml_map('fnrm',fnrm,fmt='(1pe21.14)')
!!$            call yaml_map('beta',beta,fmt='(1pe21.14)')
!!$            call yaml_map('beta_stretch',beta_stretch,fmt='(1pe21.14)')
!!$            call geometry_output(fmax,fnrm,fluct)
!!$         call yaml_mapping_close()
      end if

      etot    = etotp
      etotold = etot
      !copy outs_datatype
      !call copy_state_properties(outs,outsIO)

      if(detot .gt. maxrise)then
         if (iproc==0) then
            write(16,'(a,i0,4(1x,e9.2))') &
             "WARNING GEOPT_SQNM: Allowed energy to rise by more than maxrise: it,maxrise,detot,beta,1.d-1*betax ",&
             it,maxrise,detot,beta,1.d-1*betax
            !call f_utils_flush(16)
         endif
      endif


!      if (fnrm.le.fnrmtol) goto 1000
      !call convcheck(fmax,fluct*runObj%inputs%frac_fluct,runObj%inputs%forcemax,icheck)
      if(fmax<paropt%fmaxtol) icheck=icheck+1
      !if(icheck>5)then
      if(icheck>1)then
         goto 1000
      endif
     if(imode==2)rxyz(:,:,idx(nhist))=rxyz(:,:,idx(nhist))+&
                 beta_stretch*fstretch(:,:,idx(nhist)) !has to be after convergence check,
                                                       !otherwise energy will not match
                                                       !the true energy of rxyz(:,:,idx(nhist))

      if(nint(count_sqnm) >= nit)then!no convergence within ncount_cluster_x energy evaluations
            goto 900  !sqnm will return to caller the energies and coordinates used/obtained from the last accepted iteration step
      endif
   
      if (cosangle.gt. 0.2d0) then
         beta=beta*1.1d0
      else
         beta=max(beta*0.85d0,betax)
      endif
   
      if (debug.and.iproc==0) write(100,*) 'cosangle ',cosangle,beta

      call getSubSpaceEvecEval('(SQNM)',iproc,verbosity,nat,nhist,&
           nhistx,ndim,cutoffratio,lwork,work,idx,rxyz,fxyz,aa,rr,ff,&
           rrr,fff,eval,res,success)
      if(.not.success)stop 'subroutine minimizer_sqnm: no success in getSubSpaceEvecEval.'

   enddo!end main loop

900 continue

   !if code gets here, it failed
   if(debug.and.iproc==0) write(100,*) it,etot,fnrm
   if(iproc==0) then
        write(16,'(a,2(1x,i0))') &
       "WARNING GEOPT_SQNM: SQNM not converged: it,ncount_bigdft,ncount_cluster_x: ", &
       it,nint(count_sqnm) !,runObj%inputs%ncount_cluster_x
       !call f_utils_flush(16)
    endif
!   stop "No convergence "
   fail=.true.
   goto 2000

1000 continue!converged successfully
   
   if(iproc==0)then
         write(16,'(2(a,1x,i0))') "SQNM converged at iteration ",it,". Needed bigdft calls: ",nint(count_sqnm)
         !call f_utils_flush(16)
   endif
   !if(iproc==0)  write(*,'(a,i6,es24.15,es14.5)') 'Iterations when SQNM converged',it,etotp,fmax
   if(iproc==0) then
        call yaml_sequence(advance='no')
        call yaml_mapping_open('SQNM FINISHED') !,label='id001')
        call yaml_map('success',.true.)
        call yaml_map('iter',it,fmt='(i5)')
        call yaml_map('epot',etotp,fmt='(es20.12)')
        call yaml_map('fnrm',fnrm,fmt='(es12.5)')
        call yaml_map('fmax',fmax,fmt='(es12.5)')
        call yaml_mapping_close()
   endif
   fail=.false.
   
2000 continue
!deallocations
   call f_free(rxyz)
   call f_free(idx)
   call f_free(rxyzOld)
   call f_free(delta)
   call f_free(fxyz)
   call f_free(fxyzraw)
   call f_free(fstretch)
   call f_free(aa)
   call f_free(eval)
   call f_free(res)
   call f_free(rnorm)
   call f_free(work)
   call f_free(ff)
   call f_free(rr)
   call f_free(dd)
   call f_free(fff)
   call f_free(rrr)
   call f_free(scpr)
   call f_free(wold)
   call f_free(rcov )   
   call f_free(iconnect)
   !call f_tree_free(f_info)
   !call deallocate_state_properties(outs)
   call f_release_routine()
end subroutine sqnm
!*****************************************************************************************
!subroutine minenergyandforces_alborz(iproc,nproc,eeval,imode,outs,nat,rat,fat,fstretch,&
subroutine minenergyandforces_alborz(parini,iproc,nproc,eeval,imode,atoms,nat,rat,fat,fstretch,fxyzraw,epot,iconnect,nbond_,wold,alpha_stretch0,alpha_stretch,infocode)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, set_rat
    !use bigdft_run!module_types
    !use module_sqn
    !use module_interfaces
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in)           :: iproc,nproc,imode
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in)           :: nat
    !type(run_objects), intent(inout)       :: runObj
    !type(state_properties), intent(inout) :: outs
    integer, intent(in)           :: nbond_
    integer, intent(in)           :: iconnect(2,nbond_)
    real(8),intent(inout)        :: rat(3,nat)
    real(8),intent(out)          :: fxyzraw(3,nat)
    real(8),intent(inout)        :: fat(3,nat)
    real(8),intent(out)          :: fstretch(3,nat)
    real(8), intent(inout)       :: wold(nbond_)
    real(8), intent(in)          :: alpha_stretch0
    real(8), intent(inout)       :: alpha_stretch
    real(8), intent(inout)       :: epot
    logical, intent(in)           :: eeval
    integer,intent(out) :: infocode
    !internal
    integer:: iat

    infocode=0
    if(eeval)then
        !call vcopy(3 * runObj%atoms%astruct%nat, rat(1,1), 1,runObj%atoms%astruct%rxyz(1,1), 1)
        !atoms%rat(1:3,1:atoms%nat)=rat(1:3,1:atoms%nat)
        call set_rat(atoms,rat,setall=.true.)
        !runObj%inputs%inputPsiId=1
        !call bigdft_set_input_policy(INPUT_POLICY_MEMORY, runObj)
        !call bigdft_state(runObj,outs,infocode)
        call cal_potential_forces(parini,atoms)
    endif
    !call vcopy(3 * outs%fdim, outs%fxyz(1,1), 1, fat(1,1), 1)
    !call vcopy(3 * outs%fdim, fat(1,1), 1,fxyzraw(1,1), 1)
    fxyzraw(1:3,1:nat)=atoms%fat(1:3,1:nat)
    fat(1:3,1:nat)=atoms%fat(1:3,1:nat)
   !----------------------------------------------------------------------------
   !The following lines added by Alireza to avoid moving fixed atoms.
    do iat=1,atoms%nat
        if(.not. atoms%bemoved(1,iat)) then
            fxyzraw(1,iat)=0.d0
            fat(1,iat)=0.d0
        endif
        if(.not. atoms%bemoved(2,iat)) then
            fxyzraw(2,iat)=0.d0
            fat(2,iat)=0.d0
        endif
        if(.not. atoms%bemoved(3,iat)) then
            fxyzraw(3,iat)=0.d0
            fat(3,iat)=0.d0
        endif
    enddo
   !----------------------------------------------------------------------------
    epot=atoms%epot
    fstretch=0.d0
    if(imode==2)then
        call projectbond(nat,nbond_,rat,fat,fstretch,iconnect,&
             wold,alpha_stretch0,alpha_stretch)
    endif

end subroutine minenergyandforces_alborz
!*****************************************************************************************
subroutine give_rcov_sqnm(iproc,atoms,rcov)
    use mod_atoms, only: typ_atoms
    implicit none
    integer, intent(in) :: iproc
    type(typ_atoms), intent(in) :: atoms
    real(8), intent(out) :: rcov(atoms%nat)
    !local variables
    integer :: iat
  do iat=1,atoms%nat
     if (trim(atoms%sat(iat))=='H') then
        rcov(iat)=0.75d0
!        rcov(iat)=0.75d0*0.529177211d0
     else if (trim(atoms%sat(iat))=='LJ')then
        rcov(iat)=0.56d0
     else if (trim(atoms%sat(iat))=='He')then
        rcov(iat)=0.75d0
     else if (trim(atoms%sat(iat))=='Li')then
        rcov(iat)=3.40d0
     else if (trim(atoms%sat(iat))=='Be')then
        rcov(iat)=2.30d0
     else if (trim(atoms%sat(iat))=='B' )then
        rcov(iat)=1.55d0
     else if (trim(atoms%sat(iat))=='C' )then
        rcov(iat)=1.45d0
!        rcov(iat)=1.45d0*0.529177211d0
     else if (trim(atoms%sat(iat))=='N' )then
        rcov(iat)=1.42d0
!        rcov(iat)=1.42d0*0.529177211d0
     else if (trim(atoms%sat(iat))=='O' )then
        rcov(iat)=1.38d0
!        rcov(iat)=1.38d0*0.529177211d0
     else if (trim(atoms%sat(iat))=='F' )then
        rcov(iat)=1.35d0
     else if (trim(atoms%sat(iat))=='Ne')then
        rcov(iat)=1.35d0
     else if (trim(atoms%sat(iat))=='Na')then
        rcov(iat)=3.40d0
     else if (trim(atoms%sat(iat))=='Mg')then
        rcov(iat)=2.65d0
     else if (trim(atoms%sat(iat))=='Al')then
        rcov(iat)=2.23d0
     else if (trim(atoms%sat(iat))=='Si')then
        rcov(iat)=2.09d0
     else if (trim(atoms%sat(iat))=='P' )then
        rcov(iat)=2.00d0
     else if (trim(atoms%sat(iat))=='S' )then
        rcov(iat)=1.92d0
     else if (trim(atoms%sat(iat))=='Cl')then
        rcov(iat)=1.87d0
     else if (trim(atoms%sat(iat))=='Ar')then
        rcov(iat)=1.80d0
     else if (trim(atoms%sat(iat))=='K' )then
        rcov(iat)=4.00d0
     else if (trim(atoms%sat(iat))=='Ca')then
        rcov(iat)=3.00d0
     else if (trim(atoms%sat(iat))=='Sc')then
        rcov(iat)=2.70d0
     else if (trim(atoms%sat(iat))=='Ti')then
        rcov(iat)=2.70d0
     else if (trim(atoms%sat(iat))=='V' )then
        rcov(iat)=2.60d0
     else if (trim(atoms%sat(iat))=='Cr')then
        rcov(iat)=2.60d0
     else if (trim(atoms%sat(iat))=='Mn')then
        rcov(iat)=2.50d0
     else if (trim(atoms%sat(iat))=='Fe')then
        rcov(iat)=2.50d0
     else if (trim(atoms%sat(iat))=='Co')then
        rcov(iat)=2.40d0
     else if(trim(atoms%sat(iat))=='Ni')then
        rcov(iat)=2.30d0
     else if (trim(atoms%sat(iat))=='Cu')then
        rcov(iat)=2.30d0
     else if (trim(atoms%sat(iat))=='Zn')then
        rcov(iat)=2.30d0
     else if (trim(atoms%sat(iat))=='Ga')then
        rcov(iat)=2.10d0
     else if (trim(atoms%sat(iat))=='Ge')then
        rcov(iat)=2.40d0
     else if (trim(atoms%sat(iat))=='As')then
        rcov(iat)=2.30d0
     else if (trim(atoms%sat(iat))=='Se')then
        rcov(iat)=2.30d0
     else if (trim(atoms%sat(iat))=='Br')then
        rcov(iat)=2.20d0
     else if (trim(atoms%sat(iat))=='Kr')then
        rcov(iat)=2.20d0
     else if (trim(atoms%sat(iat))=='Rb')then
        rcov(iat)=4.50d0
     else if (trim(atoms%sat(iat))=='Sr')then
        rcov(iat)=3.30d0
     else if (trim(atoms%sat(iat))=='Y' )then
        rcov(iat)=3.30d0
     else if (trim(atoms%sat(iat))=='Zr')then
        rcov(iat)=3.00d0
     else if (trim(atoms%sat(iat))=='Nb')then
        rcov(iat)=2.92d0
     else if (trim(atoms%sat(iat))=='Mo')then
        rcov(iat)=2.83d0
     else if (trim(atoms%sat(iat))=='Tc')then
        rcov(iat)=2.75d0
     else if (trim(atoms%sat(iat))=='Ru')then
        rcov(iat)=2.67d0
     else if (trim(atoms%sat(iat))=='Rh')then
        rcov(iat)=2.58d0
     else if (trim(atoms%sat(iat))=='Pd')then
        rcov(iat)=2.50d0
     else if (trim(atoms%sat(iat))=='Ag')then
        rcov(iat)=2.50d0
     else if (trim(atoms%sat(iat))=='Cd')then
        rcov(iat)=2.50d0
     else if (trim(atoms%sat(iat))=='In')then
        rcov(iat)=2.30d0
     else if (trim(atoms%sat(iat))=='Sn')then
        rcov(iat)=2.66d0
     else if (trim(atoms%sat(iat))=='Sb')then
        rcov(iat)=2.66d0
     else if (trim(atoms%sat(iat))=='Te')then
        rcov(iat)=2.53d0
     else if (trim(atoms%sat(iat))=='I' )then
        rcov(iat)=2.50d0
     else if (trim(atoms%sat(iat))=='Xe')then
        rcov(iat)=2.50d0
     else if (trim(atoms%sat(iat))=='Cs')then
        rcov(iat)=4.50d0
     else if (trim(atoms%sat(iat))=='Ba')then
        rcov(iat)=4.00d0
     else if (trim(atoms%sat(iat))=='La')then
        rcov(iat)=3.50d0
     else if (trim(atoms%sat(iat))=='Ce')then
        rcov(iat)=3.50d0
     else if (trim(atoms%sat(iat))=='Pr')then
        rcov(iat)=3.44d0
     else if (trim(atoms%sat(iat))=='Nd')then
        rcov(iat)=3.38d0
     else if (trim(atoms%sat(iat))=='Pm')then
        rcov(iat)=3.33d0
     else if (trim(atoms%sat(iat))=='Sm')then
        rcov(iat)=3.27d0
     else if (trim(atoms%sat(iat))=='Eu')then
        rcov(iat)=3.21d0
     else if (trim(atoms%sat(iat))=='Gd')then
        rcov(iat)=3.15d0
     else if (trim(atoms%sat(iat))=='Td')then
        rcov(iat)=3.09d0
     else if (trim(atoms%sat(iat))=='Dy')then
        rcov(iat)=3.03d0
     else if (trim(atoms%sat(iat))=='Ho')then
        rcov(iat)=2.97d0
     else if (trim(atoms%sat(iat))=='Er')then
        rcov(iat)=2.92d0
     else if (trim(atoms%sat(iat))=='Tm')then
        rcov(iat)=2.92d0
     else if (trim(atoms%sat(iat))=='Yb')then
        rcov(iat)=2.80d0
     else if (trim(atoms%sat(iat))=='Lu')then
        rcov(iat)=2.80d0
     else if (trim(atoms%sat(iat))=='Hf')then
        rcov(iat)=2.90d0
     else if (trim(atoms%sat(iat))=='Ta')then
        rcov(iat)=2.70d0
     else if (trim(atoms%sat(iat))=='W' )then
        rcov(iat)=2.60d0
     else if (trim(atoms%sat(iat))=='Re')then
        rcov(iat)=2.60d0
     else if (trim(atoms%sat(iat))=='Os')then
        rcov(iat)=2.50d0
     else if (trim(atoms%sat(iat))=='Ir')then
        rcov(iat)=2.50d0
     else if (trim(atoms%sat(iat))=='Pt')then
        rcov(iat)=2.60d0
     else if (trim(atoms%sat(iat))=='Au')then
        rcov(iat)=2.70d0
     else if (trim(atoms%sat(iat))=='Hg')then
        rcov(iat)=2.80d0
     else if (trim(atoms%sat(iat))=='Tl')then
        rcov(iat)=2.50d0
     else if (trim(atoms%sat(iat))=='Pb')then
        rcov(iat)=3.30d0
     else if (trim(atoms%sat(iat))=='Bi')then
        rcov(iat)=2.90d0
     else if (trim(atoms%sat(iat))=='Po')then
        rcov(iat)=2.80d0
     else if (trim(atoms%sat(iat))=='At')then
        rcov(iat)=2.60d0
     else if (trim(atoms%sat(iat))=='Rn')then
        rcov(iat)=2.60d0
     else
        !call yaml_comment('(SQNM) no covalent radius stored for this atomtype '&
        !     //trim(atoms%sat(iat))))
        stop 'ERROR: (SQNM) no covalent radius stored for this atomtype'
     endif
     if (iproc == 0) then
        !call yaml_map('(SQNM) RCOV:'//trim(atoms%astruct%atomnames&
        !        (atoms%astruct%iatype(iat))),rcov(iat))
     endif
  enddo
end subroutine give_rcov_sqnm
!*****************************************************************************************
subroutine getSubSpaceEvecEval(label,iproc,verbosity,nat,nhist,nhistx,ndim,cutoffratio,lwork,work,idx,rxyz,fxyz,aa,rr,ff,rrr,fff,eval,res,success)
    !use module_base
    !use yaml_output
    !hard-coded parameters:
    !threshold for linear dependency:
    !if (eval(idim)/eval(nhist).gt.1.e-4_gp) then
    implicit none
    !parameters
    integer, intent(in) :: iproc,verbosity,nat,nhist,nhistx,lwork
    character(len=*), intent(in) :: label
    integer, intent(out) :: ndim
    integer, intent(in) :: idx(0:nhistx)
    real(8), intent(in) :: rxyz(3,nat,0:nhistx),fxyz(3,nat,0:nhistx)
    real(8), intent(out) :: aa(nhistx,nhistx),eval(nhistx)
    real(8), intent(out) :: work(lwork)
    real(8), intent(out) :: rr(3,nat,0:nhistx), ff(3,nat,0:nhistx)
    real(8), intent(out) :: rrr(3,nat,0:nhistx), fff(3,nat,0:nhistx)
    real(8), intent(out) :: res(nhistx)
    real(8), intent(in) :: cutoffratio
    logical, intent(out) :: success
    !internal
    integer :: i,j,l,iat,info,idim,jdim
    real(8) :: tt
    real(8) :: rnorm(nhistx)

    success=.false.

    ! calculate norms
    do i=1,nhist
        rnorm(i)=0.d0
         do iat=1,nat
             do l=1,3
                rnorm(i)=rnorm(i) + (rxyz(l,iat,idx(i))-rxyz(l,iat,idx(i-1)))**2
             enddo
         enddo
         rnorm(i)=1.d0/sqrt(rnorm(i))
    enddo

    !find linear dependencies via diagonalization of overlap matrix   
    !build overlap matrix:
    do i=1,nhist
        do j=1,nhist
            aa(i,j)=0.d0
            do iat=1,nat
                do l=1,3
                aa(i,j)=aa(i,j) + (rxyz(l,iat,idx(i))-rxyz(l,iat,idx(i-1)))&
                &*(rxyz(l,iat,idx(j))-rxyz(l,iat,idx(j-1)))
                enddo
            enddo
            aa(i,j)=aa(i,j)*rnorm(i)*rnorm(j)
         enddo
    enddo

    call DSYEV('V',"L",nhist,aa,nhistx,eval,work,lwork,info)
    if (info.ne.0) then
        !call yaml_warning(trim(adjustl(label))//' 1st DSYEV '//&
        !'(Overlapmatrix) in getSupSpaceEvecEval failed with info: '//&
        !trim(yaml_toa(info))//', iproc: '//trim(yaml_toa(iproc)))
        write(*,*) trim(adjustl(label)),' 1st DSYEV ','(Overlapmatrix) in getSupSpaceEvecEval failed with info: ',info
        return
!        stop 'info'
    endif
    if(iproc==0 .and. verbosity>=3)then
        do i=1,nhist
            !call yaml_scalar(trim(adjustl(label))//' Overlap '//&
            !'eigenvalues: '//trim(yaml_toa(i))//' '//&
            !trim(yaml_toa(eval(i))))
            write(*,*) 'eigenvalues: ',i,eval(i)
        enddo
    endif

    do idim=1,nhist
        do iat=1,nat
            do l=1,3
                rr(l,iat,idim)=0.d0
                ff(l,iat,idim)=0.d0
            enddo
        enddo
    enddo

    ndim=0
    do idim=1,nhist
        !remove linear dependencies by using
        !the overlap-matrix eigenvalues:
        if (eval(idim)/eval(nhist).gt.cutoffratio) then    ! HERE
            ndim=ndim+1

            do jdim=1,nhist
                do iat=1,nat
                    do l=1,3
                         rr(l,iat,ndim)=rr(l,iat,ndim)+&
                                 aa(jdim,idim)*rnorm(jdim)*&
                                 (rxyz(l,iat,idx(jdim))-rxyz(l,iat,idx(jdim-1)))
                         ff(l,iat,ndim)=ff(l,iat,ndim)-&
                                 aa(jdim,idim)*rnorm(jdim)*&
                                 (fxyz(l,iat,idx(jdim))-fxyz(l,iat,idx(jdim-1)))

                    enddo
                enddo
            enddo

            do iat=1,nat
                do l=1,3
                    rr(l,iat,ndim)=rr(l,iat,ndim)/&
                                    sqrt(abs(eval(idim)))
                    ff(l,iat,ndim)=ff(l,iat,ndim)/&
                                    sqrt(abs(eval(idim)))
                enddo
            enddo
        endif
    enddo

    ! Hessian matrix in significant orthogonal subspace
    do i=1,ndim
        do j=1,ndim
            aa(i,j)=0.d0
            do iat=1,nat
                do l=1,3
                    aa(i,j)=aa(i,j) + 0.5d0*(rr(l,iat,i)*ff(l,iat,j)+&
                                      &rr(l,iat,j)*ff(l,iat,i))
                enddo
            enddo
        enddo
    enddo

    call DSYEV('V',"L",ndim,aa,nhistx,eval,work,lwork,info)
    if (info.ne.0) then
        !call yaml_warning(trim(adjustl(label))//' 2nd DSYEV '//&
        !'(subpsace hessian) in getSupSpaceEvecEval failed with info: '&
        ! //trim(yaml_toa(info))//', iproc:'//trim(yaml_toa(iproc)))
        write(*,*) trim(adjustl(label)),' 2nd DSYEV ','(subpsace hessian) in getSupSpaceEvecEval failed with info: ',info
        return
!        stop 'info'
    endif

    ! calculate vectors in full 3*nat-dim space
    do i=1,ndim
        do iat=1,nat
            do l=1,3
                rrr(l,iat,i)=0.d0
                fff(l,iat,i)=0.d0
            enddo
        enddo
    enddo

    do i=1,ndim
        tt=0.d0
        do j=1,ndim
            do iat=1,nat
                do l=1,3
                    rrr(l,iat,i)=rrr(l,iat,i) + aa(j,i)*rr(l,iat,j)
                    fff(l,iat,i)=fff(l,iat,i) + aa(j,i)*ff(l,iat,j)
                enddo
            enddo
        enddo
        do iat=1,nat
            do l=1,3
                tt=tt+(fff(l,iat,i)-eval(i)*rrr(l,iat,i))**2
            enddo
        enddo
        !residuue according to Weinstein criterion
        res(i)=sqrt(tt)
        if(iproc==0 .and. verbosity>=3)&
            !call yaml_scalar(trim(adjustl(label))//' i, '//&
            !'eigenvalue, residue: '//trim(yaml_toa(i))//' '//&
            !trim(yaml_toa(eval(i)))//' '//trim(yaml_toa(res(i))))
            write(*,*) trim(adjustl(label)),' i, ','eigenvalue, residue: ',i,' ',eval(i),' ',res(i)
    enddo
    success=.true.
end subroutine getSubSpaceEvecEval
!*****************************************************************************************
subroutine modify_gradient(nat,ndim,rrr,eval,res,fxyz,alpha,dd)
    !use module_base
    implicit none
    !parameters
    integer, intent(in) :: nat
    integer, intent(in) :: ndim
    real(8), intent(out) :: dd(3,nat)
    real(8), intent(in) :: fxyz(3,nat)
    real(8), intent(in) :: rrr(3,nat,ndim)
    real(8), intent(in) :: eval(ndim)
    real(8), intent(in) :: res(ndim)
    real(8), intent(in) :: alpha
    !internal
    integer :: iat,i,l
    real(8) :: scpr(ndim)
    real(8) :: tt

    ! decompose gradient

    do iat=1,nat
        do l=1,3
            dd(l,iat)=-fxyz(l,iat)
        enddo
    enddo
    do i=1,ndim
        scpr(i)=0.d0
        do iat=1,nat
            do l=1,3
                scpr(i)=scpr(i)-fxyz(l,iat)*rrr(l,iat,i)
            enddo
        enddo
        do iat=1,nat
            do l=1,3
                dd(l,iat)=dd(l,iat)-scpr(i)*rrr(l,iat,i)
            enddo
        enddo
    enddo

    !simple sd in space orthogonal to relevant subspace
    do iat=1,nat
        do l=1,3
            dd(l,iat)=dd(l,iat)*alpha
        enddo
    enddo

    do i=1,ndim
    !quasi newton in relevant subspace
        tt=scpr(i)/sqrt(eval(i)**2+res(i)**2)
        do iat=1,nat
            do l=1,3
                dd(l,iat)=dd(l,iat)+tt*rrr(l,iat,i)
            enddo
        enddo
    enddo
end subroutine modify_gradient
!*****************************************************************************************
!has to be called before findsad (if operating in biomolecule mode)
subroutine findbonds(label,iproc,verbosity,atoms,rcov,nbond,iconnect)
    use mod_atoms, only: typ_atoms, get_rat
    implicit none
    integer, intent(in) :: iproc,verbosity
    character(len=*), intent(in) :: label
    type(typ_atoms), intent(in) :: atoms
    real(8), intent(in) :: rcov(atoms%nat)
    integer, intent(out) :: nbond
    integer, intent(out) :: iconnect(2,1000)
    !local variables
    integer :: iat,jat
    real(8) :: dist2
    real(8), allocatable:: rat(:,:)
    allocate(rat(3,atoms%nat))
    call get_rat(atoms,rat)
    nbond=0
    do iat=1,atoms%nat
        do jat=1,iat-1
            dist2=(rat(1,iat)-rat(1,jat))**2+&
                  (rat(2,iat)-rat(2,jat))**2+&
                  (rat(3,iat)-rat(3,jat))**2
            if (dist2.le.(1.2d0*(rcov(iat)+rcov(jat)))**2) then
                nbond=nbond+1
                if (nbond.gt.1000) then
                    write(*,*) 'nbond>1000, increase size of iconnect in routine which calls subroutine findbonds'
                endif
                iconnect(1,nbond)=iat
                iconnect(2,nbond)=jat
            endif
        enddo
    enddo
    if(nbond==0)stop 'nbond==0'
    if(iproc==0.and.verbosity>=2) &
        write(*,*) trim(adjustl(label)),' Found',nbond,' bonds.'
        !call yaml_scalar(trim(adjustl(label))//&
        !' Found'//trim(yaml_toa(nbond))//' bonds.')
    deallocate(rat)
end subroutine findbonds
!*****************************************************************************************
subroutine projectbond(nat,nbond,rat,fat,fstretch,iconnect,wold,alpha_stretch0,alpha_stretch)
    !use module_base, only: gp
    implicit none
    integer, intent(in) :: nat
    integer, intent(in) :: nbond
    real(8), intent(in) :: rat(3,nat)
    real(8), intent(inout) :: fat(3,nat)
    real(8), intent(inout) :: fstretch(3,nat)
    integer, intent(in) :: iconnect(2,nbond)
    real(8), intent(inout) :: wold(nbond)
    real(8), intent(in) :: alpha_stretch0
    real(8), intent(inout) :: alpha_stretch
    !internal
    integer :: iat,jat,ibond,jbond,l,nsame,info
    real(8) :: ss(nbond,nbond),w(nbond),vv(3,nat,nbond)
    real(8) :: per
    !functions
    real(8) :: ddot
    
    
    fstretch=0.d0

    !|v_i> := |rat_k>-|rat_l>
    !|F>=sum_i c_i* |v_i>
    !<v_j|F> = sum_i c_i <v_j|v_i>

    ! set up positional overlap matrix
    vv=0.d0
    do ibond=1,nbond
        iat=iconnect(1,ibond)
        jat=iconnect(2,ibond)
        do l=1,3
            vv(l,iat,ibond)=rat(l,jat)-rat(l,iat)
            vv(l,jat,ibond)=rat(l,iat)-rat(l,jat)
        enddo
    enddo
    
    ss=0.d0
    w=0.d0
    do ibond=1,nbond
        do jbond=1,nbond
            ss(ibond,jbond)=&
                       ddot(3*nat,vv(1,1,ibond),1,vv(1,1,jbond),1)
        enddo
        w(ibond)=ddot(3*nat,vv(1,1,ibond),1,fat(1,1),1)
    enddo
    
    nsame=0
    do ibond=1,nbond
        if ( wold(ibond)*w(ibond).gt.0.d0) nsame=nsame+1
        wold(ibond)=w(ibond)
    enddo
    !determine feedback on streching components of force
    per=real(nsame,8)/nbond
    if (per.gt. 0.66d0) then
        alpha_stretch=alpha_stretch*1.1d0
    else
        alpha_stretch=max(1.d-2*alpha_stretch0,&
                           alpha_stretch/1.1d0)
    endif

    call DPOSV('L', nbond, 1, ss, nbond, w, nbond, info )
    if (info.ne.0) then
        write(*,*)'info',info
        stop 'info DPOSV in minenergyforces'
    endif

    ! calculate projected force
    fstretch=0.d0
    do ibond=1,nbond
        do iat=1,nat
            do l=1,3
                fstretch(l,iat)=fstretch(l,iat)+w(ibond)*&
                                vv(l,iat,ibond)
            enddo
        enddo
    enddo

    !
    do iat=1,nat
        do l=1,3
            fat(l,iat)=fat(l,iat)-fstretch(l,iat)
        enddo
    enddo
end subroutine projectbond
!*****************************************************************************************
