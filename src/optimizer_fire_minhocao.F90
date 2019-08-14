subroutine GEOPT_FIRE_MHM(parini,parres,latvec_in,xred_in,fcart_in,strten_in,vel_in,vel_lat_in,vvol_in,etot_in,iprec,counter,folder)
 use global, only: units,max_kpt,ka1,kb1,kc1,confine
 use defs_basis
 use mod_fire
 use interface_code
 use modsocket, only: sock_extra_string
 use mod_parini, only: typ_parini
 use yaml_output
implicit none
 type(typ_parini), intent(in):: parini
 type(typ_parini), intent(inout):: parres
 real(8) :: latvec_in(3,3),xred_in(3,parini%nat),vel_in(3,parini%nat),fcart_in(3,parini%nat),etot_in,strten_in(6),vel_lat_in(3,3),vvol_in
!*********************************************************************
!Variables for my MD part
! real(8), parameter :: Ha_eV=27.21138386_dp ! 1 Hartree, in eV
! real(8), parameter :: kb_HaK=8.617343d-5/Ha_eV ! Boltzmann constant in Ha/K
! real(8), parameter :: amu_emass=1.660538782d-27/9.10938215d-31 ! 1 atomic mass unit, in electronic mass
 real(8):: amass(parini%nat)
 real(8):: pressure_ener
 real(8):: pressure_md(3,3)
 real(8):: pressure
 real(8):: vol
 real(8):: vol0
 real(8):: volprev
 real(8),dimension(3,parini%nat):: xcart
 real(8),dimension(3,parini%nat):: fposcur
 real(8),dimension(3,parini%nat):: fpospred
 real(8),dimension(3,parini%nat):: vpospred
 real(8),dimension(3,parini%nat):: poscur
 real(8),dimension(3,parini%nat):: vxyz
 real(8),dimension(3,parini%nat):: vposcur
 real(8),dimension(3,parini%nat):: pospred
 real(8),dimension(3,parini%nat):: accposcur
 real(8),dimension(3,parini%nat):: accpospred
 real(8),dimension(3,parini%nat):: accposprev
 real(8),dimension(3,parini%nat):: dxred
 real(8),dimension(3,3):: dlatvec
 real(8),dimension(3,3):: latvec
 real(8),dimension(3,3):: latvec0
 real(8),dimension(3,3):: latinv
 real(8),dimension(3,3):: unitmat
 real(8),dimension(3,3):: elmatcur
 real(8),dimension(3,3):: flatcur
 real(8),dimension(3,3):: flatpred
 real(8),dimension(3,3):: velmatcur
 real(8),dimension(3,3):: tmplat
 real(8),dimension(3,3):: tmplatt
 real(8),dimension(3,3):: tpred
 real(8),dimension(3,3):: f0
 real(8),dimension(3,3):: f0inv
 real(8),dimension(3,3):: latcur
 real(8),dimension(3,3):: vlatcur
 real(8),dimension(3,3):: vlat
 real(8),dimension(3,3):: lattrans
 real(8),dimension(3,3):: latdottrans
 real(8),dimension(3,3):: a 
 real(8),dimension(3,3):: g
 real(8),dimension(3,3):: gdot
 real(8),dimension(3,3):: ginv
 real(8),dimension(3,3):: gtot
 real(8),dimension(3,3):: sigma
 real(8),dimension(3,3):: sigmatrans
 real(8),dimension(3,3):: latpred
 real(8),dimension(3,3):: tcur
 real(8),dimension(3,3):: vlatpred
 real(8),dimension(3,3):: velmatpred
 real(8),dimension(3,3):: str_matrix
 real(8),dimension(3,3):: acclatcur
 real(8),dimension(3,3):: acclatpred
 real(8),dimension(3,3):: acclatprev
 real(8):: accvolcur
 real(8):: accvolpred
 real(8):: accvolprev
 real(8):: fvolcur
 real(8):: fvolpred
 real(8):: volcur
 real(8):: volpred
 real(8):: vvolcur
 real(8):: vvolpred
 real(8):: volpred_1_3
 real(8):: vol_1_3_inv
 real(8):: vol_1_3
 real(8):: vposcurtmp(3)
 real(8):: vpospredtmp(3)
 real(8):: crossp(3)
 real(8):: velcm(3)
 real(8):: latmass
 real(8):: latmassinv
 real(8):: latmass0
 real(8):: ekinatom
 real(8):: ekinlat
 real(8):: ekinlat_prev
 real(8):: ekinatom_prev
 real(8):: rkin
 real(8):: enthalpy
 real(8):: enthalpy_min
 real(8):: enmin1
 real(8):: enmin2
 real(8):: ent_pos_0
 real(8):: en0000
 real(8):: e_rxyz
 real(8):: econs_max
 real(8):: econs_min
 real(8):: torquenrm 
 real(8):: dt 
 real(8):: ecut_tmp
 real(8):: toldff_tmp
 real(8):: dstr(6)
 real(8):: strtarget(6)
 real(8):: counter 
 integer:: i
 integer:: j
 integer:: iat
 integer:: itime
 integer:: nummax
 integer:: nummin
 integer:: options
 integer:: iexit
 integer:: iprec
 integer:: md_type
 logical:: getwfk
 character(40)::filename,folder
 character(4) ::fn4

! FIRE VARIABLES
 real(8):: alpha,P,P_at,P_lat,fmax,fmax_at,fmax_lat,fall(3,parini%nat+3),fallnorm,vall(3,parini%nat+3),vallnorm
 integer:: nstep,istr
 logical:: multiprec
 integer:: iprec_cur
 real(8):: tolmxf_switch,alpha_lat,f_lat
 logical:: cellfix_done
 real(8):: cellfix_switch
!Latvec_variable
 integer:: latvec_io
latvec_io=0
!multiprec is hardcoded and, if true, starts a geopt with iprec==2, and then switches 
!to iprec==1 when the fmax==tolmxf_switch. The switch only occurs once
 max_kpt=.false.
 multiprec=.true.
 tolmxf_switch=1.d-3
 enthalpy_min=1.d10
 cellfix_switch=1.d-3
 cellfix_done=.false.


!Lattice alpha
alpha_lat=20.d0

call yaml_sequence_open('atoms info in FIRE')
!Assign masses to each atom (for FIRE, all atoms have the same mass)
 do iat=1,parini%nat
!   amass(iat)=amu_emass*1.d0
   amass(iat)=amu_emass*parini%amu(parini%typat_global(iat))
   !if(parini%verb.gt.0) write(*,'(a,i5,2(1x,es15.7))') " # FIRE: iat, AMU, EM: ", iat, amass(iat)/amu_emass,amass(iat)
   if(parini%verb.gt.0) then
    call yaml_sequence(advance='no')
    call yaml_mapping_open('atom',flow=.true.)
    call yaml_map('iat',iat,fmt='(i5)')
    call yaml_map('AMU',amass(iat)/amu_emass,fmt='(es15.7)')
    call yaml_map('EM',amass(iat),fmt='(es15.7)')
    call yaml_mapping_close()
   endif
 enddo
call yaml_sequence_close()

!Set FIRE parameters
alpha=alphastart
!if(fixlat(7)) alpha=1.d0
nstep=0
latmass0=latmass_rel_fire*sum(amass(:))/real(parini%nat,8)
dt=parini%paropt_geopt%dt_start
vel_in=0.d0
vel_lat_in=0.d0
P=0.d0
P_at=0.d0
P_lat=0.d0

!MD-type: 1 for PR and 2 for Cleveland and 3 for Wentzcovitch 4 Andersen
  if(.not.parini%fixlat(7)) then 
      md_type=1
  else
      md_type=4
      latmass0=0.005d0*latmass_rel_fire*sum(amass(:))/real(parini%nat,8)
  endif

    call yaml_mapping_open('FIRE info',flow=.true.)
    call yaml_map('iat',iat)
    call yaml_map('AMU',latmass0/amu_emass,fmt='(es15.7)')
    call yaml_map('EM',latmass0,fmt='(es15.7)')
    call yaml_mapping_close()
    !write(*,'(a,i5,2(1x,es15.7))') " # FIRE: latmass, AMU, EM: ", iat, latmass0/amu_emass,latmass0

!******************************************************************
!NEW VERISON OF MD: Directly implemented, simplest Parrinello-Rahman MD
!Can be accessed if the ionmov of idtset2 is equal to -13

!Here we split the routine in the Abinit native part and my new implentation
!write(*,'(a)') ' # Entering FIRE based on standalone PR MD '
call yaml_comment('Entering FIRE based on standalone PR MD',hfill='~')
!Set temporary variables, initially
  vxyz(:,:)=vel_in(:,:)
  pressure=parini%target_pressure_habohr
  unitmat=0.d0
  do i=1,3
    unitmat(i,i)=1.d0
  enddo
  vlat=vel_lat_in  !The initial cell velocity
  itime=0

!Set options=1 for Velocity Verlet of cell dynamics
!Set options=2 for Normal Verlet of cell dynamics
!Set options=3 for Beeman integration scheme, corrector-predictor
  options=3
  if(options.lt.1.or.options.gt.3) stop "Wrong integrator option"


  if(md_type.lt.1.or.md_type.gt.4) stop "Wrong algo option"

if(parini%verb.gt.0) then
    call yaml_map('GEOPT Algorithm',md_type)
    call yaml_map('MD Integrator',options)
    !write(*,'(a,i3,a,i3)') " # GEOPT Algorithm: ",md_type, ", MD Integrator: ",options
endif

!Now we run my implementation of PR MD
pressure_md=0.d0
pressure_md(1,1)=1.d0;pressure_md(2,2)=1.d0;pressure_md(3,3)=1.d0
pressure_ener=0.d0;pressure_md=pressure_md*pressure  !Here the pressure is not passed to the energyandforces, so we move on the ENERGY surface

!Transform the initial velocities given in cartesian coordinates into velocities of the intenal coordinates
       latvec=latvec_in
       call getvol(latvec,vol)
if(md_type==4) then
       vol_1_3=vol**(1.d0/3.d0)
       vol_1_3_inv=1.d0/vol_1_3
       do iat=1,parini%nat
          vposcur(:,iat)=vxyz(:,iat)*vol_1_3_inv
       enddo
else
       call invertmat(latvec,latinv,3)
       do iat=1,parini%nat
          vposcur(:,iat)=matmul(latinv,vxyz(:,iat))
       enddo
endif

!Compute f0
        call cross_product(latvec(:,2),latvec(:,3),crossp); sigma(:,1)=crossp
        call cross_product(latvec(:,3),latvec(:,1),crossp); sigma(:,2)=crossp
        call cross_product(latvec(:,1),latvec(:,2),crossp); sigma(:,3)=crossp
        sigmatrans(:,1)=sigma(1,:)
        sigmatrans(:,2)=sigma(2,:)
        sigmatrans(:,3)=sigma(3,:)
        f0=matmul(sigmatrans,sigma)
        call invertmat(f0,f0inv,3)


!Here we perform the initial force call
!****************************************************************************************************************        
!****************************************************************************************************************        
       vel_in=0.d0
       getwfk=.false.
       write(fn4,'(i4.4)') 0
       sock_extra_string="FIRE"//trim(fn4)
       call get_energyandforces_single(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,getwfk)
!Single point calculation finished. Now returning to MD part. All output variables are now in *_in
!****************************************************************************************************************        
!FIRE: check for convergence
call convcheck(parini,parini%nat,latvec_in,fcart_in,strten_in,parini%target_pressure_habohr,parini%paropt_geopt%strfact,fmax,fmax_at,fmax_lat,parini%paropt_geopt%fmaxtol,iexit)
!!!!!Compute maximal component of forces, EXCLUDING any fixed components
!!! fmax=0.0d0
!!! do iat=1,nat
!!!   do i=1,3
!!!!     if (dtsets(1)%iatfix(i,iat) /= 1) then
!!!       if( abs(fcart_in(i,iat)) >= fmax ) fmax=abs(fcart_in(i,iat))
!!!!     end if
!!!   end do
!!! end do
!!! strtarget=0.d0
!!! strtarget(1:3)=-target_pressure_habohr
!!! dstr(:)=strten_in(:)-strtarget(:)
!!!!Eventually take into account the stress
!!! do istr=1,6
!!!     if(abs(dstr(istr))*strfact >= fmax ) fmax=abs(dstr(istr))*strfact
!!! end do
!!! iexit=0
!!! if(fmax.lt.tolmxf) iexit=1
!!!
!Initial iprec after running the first force call
 if(multiprec) iprec=2

 fall=0.d0
 call fpos_flat(parini,pressure_md,fall(:,1:parini%nat),fall(:,parini%nat+1:parini%nat+3),strten_in,fcart_in,latvec_in,md_type)
 if(md_type==4) f_lat=fall(1,parini%nat+1)

!Keep track of volume
        vol0=vol
        latvec0=latvec/vol**(1.d0/3.d0)
        volprev=vol

        if(md_type==1) then
!!This is for the original Rahman Parrinello. You also need to change the volume-mass factor
        latmass=latmass0
        latmassinv=1.d0/latmass
        elseif(md_type==2) then
!This is for cleveland
        latmass=latmass0/vol**(4.d0/3.d0)
        latmassinv=1.d0/latmass
        elseif(md_type==3) then
!This is for Wentzcovitch
        latmass=latmass0
        latmassinv=1.d0/latmass
        elseif(md_type==4) then
!This is for Andersen
        latmass=latmass0
        latmassinv=1.d0/latmass
        endif

!Initialize internal coordinates and lattice vectors
if(md_type==4) then
        call rxyz_int2cart(latvec,xred_in,xcart,parini%nat)
        poscur=xcart*vol_1_3_inv           !call rxyz_cart2int(latvec,poscur,rxyz,nat)
        volcur=vol
        vvolcur=vvol_in
        latcur=latvec
else
        poscur=xred_in           !call rxyz_cart2int(latvec,poscur,rxyz,nat)
        latcur=latvec
        vlatcur=vlat
endif
!!Write every step
!        call rxyz_int2cart(latcur,poscur,rxyz,nat)
!        call wtpos_inter(nat,rxyz,latcur,0)

!MHM: initialize variable to track the minima/maxima***********
    call yaml_comment('MINHOP start FIRE')
    !write(*,*) '# MINHOP start FIRE'
    nummax=0
    nummin=0
    enmin1=0.d0
    en0000=0.d0
    econs_max=-1.d100
    econs_min=1.d100
!**************************************************************

!MHM: Write output to file in every step***********************************
!INITIAL STEP, STILL THE SAME STRUCTURE AS INPUT
       call yaml_map('Pressure',pressure)
       call yaml_map('Energy',etot_in)
       !write(*,*) "Pressure, Energy",pressure,etot_in
       call get_enthalpy(latvec_in,etot_in,pressure,enthalpy)
       ent_pos_0=enthalpy
       en0000=enthalpy-ent_pos_0
if(parini%verb.gt.0) then
       write(fn4,'(i4.4)') itime
       filename=trim(folder)//"posgeopt."//fn4//".ascii"
       units=units
       !write(*,*) "# Writing the positions in FIRE:",filename
       call yaml_map('Writing the positions in FIRE',trim(filename))
       call write_atomic_file_ascii(parini,filename,parini%nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &parini%char_type(1:parini%ntypat_global),parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,etot_in,pressure,enthalpy,en0000)
       if(parini%verb.ge.3) then
       filename=trim(folder)//"posgeopt."//fn4//".vasp"
       call write_atomic_file_poscar(parini,filename,parini%nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &parini%char_type(1:parini%ntypat_global),parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,etot_in,pressure,enthalpy,en0000)
       endif
endif
!*********************************************************************

!Exit if the forces are already converged at the beginning
        if(iexit==1) then
            call yaml_mapping_open('FIRE converged before entering iterations',flow=.true.)
            call yaml_map('itime',itime)
            call yaml_map('enthalpy',enthalpy,fmt='(es20.12)')
            call yaml_map('fmax',fmax,fmt='(es12.5)')
            call yaml_mapping_close()
          !write(*,'(a,i4,2(1x,es25.15))') " # FIRE converged before entering iterations", itime,enthalpy,fmax
          max_kpt=.false.
          return 
        endif

call acceleration_fire(pressure_md,accposcur,acclatcur,accvolcur,vposcur,&
     &vlatcur,vvolcur,strten_in,fcart_in,latcur,amass,latmass,f0inv,md_type,parini%nat)
        accposprev=accposcur
        acclatprev=acclatcur
        accvolprev=accvolcur

!Eliminate the unnecessary atmic and cell components
call elim_fixed_at(parini,parini%nat,accposcur);call elim_fixed_at(parini,parini%nat,accposprev)
if(md_type.ne.4) call elim_fixed_lat(parini,latcur,acclatcur)
if(md_type.ne.4) call elim_fixed_lat(parini,latcur,acclatprev)
call elim_fixed_at(parini,parini%nat,vposcur)
if(md_type.ne.4) call elim_fixed_lat(parini,latcur,vlatcur)


    call yaml_sequence_open('FIRE optimization iterations')
!FIRE cycles
        do itime=1,parini%paropt_geopt%nit
            call yaml_sequence(advance='no')
!          if(itime.ne.1) e_rxyz=enthalpy  !e_rxyz=e_rxyz+pressure_md(1,1)*vol   
!Check the torque on the cell for rotation
          call torque_cell(latcur,vlatcur,torquenrm)
!Now perform the step
          velcm=0.d0
          rkin=0.d0
          if(options==1) then
!For velocity verlet
!             latpred(:,:)=latcur(:,:) + dt*vlatcur(:,:) + (.5d0*dt*dt*acclatcur)
             dlatvec(:,:)=dt*vlatcur(:,:) + (.5d0*dt*dt*acclatcur)
          do iat=1,parini%nat
              dxred(:,iat)=dt*vposcur(:,iat) + dt*dt*accposcur(:,iat)
!             pospred(:,iat)=poscur(:,iat) + dt*vposcur(:,iat) + dt*dt*accposcur(:,iat)
          enddo
          call propagate(parini,parini%nat,poscur,latcur,dxred,dlatvec,pospred,latpred)
          do iat=1,parini%nat
             vpospred(:,iat)=(pospred(:,iat)-poscur(:,iat))/dt   !Update dummy variable vposcur for next iteration
          enddo
            if(md_type==4) then
!               volpred=volcur+dt*vvolcur+0.5d0*dt*dt*accvolcur
               volpred=volcur+alpha_lat*fall(1,parini%nat+1)
            endif
          elseif(options==2) then
!Instead if normal verlet is used:
!              latpred(:,:)=latcur(:,:) + dt*vlatcur + dt*dt*acclatcur
              dlatvec(:,:)=dt*vlatcur + dt*dt*acclatcur
          do iat=1,parini%nat
              dxred(:,iat)=dt*vposcur(:,iat) + dt*dt*accposcur(:,iat)
!             pospred(:,iat)=poscur(:,iat) + dt*vposcur(:,iat) + dt*dt*accposcur(:,iat)
          enddo
          call propagate(parini,parini%nat,poscur,latcur,dxred,dlatvec,pospred,latpred)
          vlatpred(:,:)=(latpred-latcur)/dt
          do iat=1,parini%nat
              vpospred(:,iat)=(pospred(:,iat)-poscur(:,iat))/dt   !Update dummy variable vposcur for next iteration
          enddo
            if(md_type==4) then
!               volpred=volcur+dt*vvolcur+dt*dt*accvolcur
!               vvolpred=(volpred-volcur)/dt
               volpred=volcur+alpha_lat*fall(1,parini%nat+1)
            endif
          elseif(options==3) then
!Predictor part of Beeman for cell
!              latpred(:,:)=latcur(:,:) + dt*vlatcur + dt*dt/6.d0*(4.d0*acclatcur(:,:)-acclatprev(:,:))
              dlatvec(:,:)=dt*vlatcur + dt*dt/6.d0*(4.d0*acclatcur(:,:)-acclatprev(:,:))
          do iat=1,parini%nat
!Predictor part of Beeman
!             pospred(:,iat)=poscur(:,iat) + dt*vposcur(:,iat) + dt*dt/6.d0*(4.d0*accposcur(:,iat)-accposprev(:,iat))             
             dxred(:,iat)=dt*vposcur(:,iat) + dt*dt/6.d0*(4.d0*accposcur(:,iat)-accposprev(:,iat))
          enddo
          call propagate(parini,parini%nat,poscur,latcur,dxred,dlatvec,pospred,latpred)
          vlatpred(:,:)=vlatcur(:,:) + 0.5d0*(3.d0*acclatcur(:,:)-acclatprev(:,:))
          do iat=1,parini%nat
             vpospred(:,iat)=vposcur(:,iat)+0.5d0*(3.d0*accposcur(:,iat)-accposprev(:,iat))
          enddo
            if(md_type==4) then
!               volpred=volcur+dt*vvolcur+dt*dt/6.d0*(4.d0*accvolcur-accvolprev)
!               vvolpred=volcur+0.5d0*(3.d0*accvolcur-accvolprev)
               volpred=volcur+alpha_lat*fall(1,parini%nat+1)
            endif
          endif

if(md_type==4) then
  do iat=1,parini%nat
     vposcurtmp=vposcur(:,iat)*vol_1_3
     rkin=rkin+amass(iat)*(vposcurtmp(1)**2+vposcurtmp(2)**2+vposcurtmp(3)**2)
     velcm=velcm+vposcurtmp*amass(iat)
  enddo
else
  do iat=1,parini%nat
             vposcurtmp=matmul(latcur,vposcur(:,iat))
             rkin=rkin+amass(iat)*(vposcurtmp(1)**2+vposcurtmp(2)**2+vposcurtmp(3)**2)
             velcm=velcm+vposcurtmp*amass(iat)
  enddo
endif



!!!
!!!
!!!
!!!!Now perform the step
!!!          velcm=0.d0
!!!          rkin=0.d0
!!!          do iat=1,nat
!!!          if(options==1.or.options==2) then
!!!             pospred(:,iat)=poscur(:,iat) + dt*vposcur(:,iat) + dt*dt*accposcur(:,iat)
!!!             vpospred(:,iat)=(pospred(:,iat)-poscur(:,iat))/dt   !Update dummy variable vposcur for next iteration
!!!          elseif(options==3) then
!!!!Predictor part of Beeman
!!!             pospred(:,iat)=poscur(:,iat) + dt*vposcur(:,iat) + dt*dt/6.d0*(4.d0*accposcur(:,iat)-accposprev(:,iat))
!!!             vpospred(:,iat)=vposcur(:,iat)+0.5d0*(3.d0*accposcur(:,iat)-accposprev(:,iat))
!!!          endif
!!!             vposcurtmp=matmul(latcur,vposcur(:,iat))
!!!             rkin=rkin+amass(iat)*(vposcurtmp(1)**2+vposcurtmp(2)**2+vposcurtmp(3)**2)
!!!             velcm=velcm+vposcurtmp*amass(iat)
!!!          enddo
!!!
!!!!Now calculate steps for the lattice according to the velocity verlet algorithm
!!!          if(options==1) then
!!!!For velocity verlet
!!!              latpred(:,:)=latcur(:,:) + dt*vlatcur(:,:) + (.5d0*dt*dt*acclatcur)
!!!          elseif(options==2) then
!!!!Instead if normal verlet is used:
!!!              latpred(:,:)=latcur(:,:) + dt*vlatcur + dt*dt*acclatcur
!!!              vlatpred(:,:)=(latpred-latcur)/dt
!!!          elseif(options==3) then
!!!!Predictor part of Beeman for cell
!!!              latpred(:,:)=latcur(:,:) + dt*vlatcur + dt*dt/6.d0*(4.d0*acclatcur(:,:)-acclatprev(:,:))
!!!              vlatpred(:,:)=vlatcur(:,:) + 0.5d0*(3.d0*acclatcur(:,:)-acclatprev(:,:))
!!!          else
!!!              stop "Wrong option for the Cell part in MD"
!!!          endif

!Kinetic energy of atoms
          ekinatom=0.5d0*rkin
          !if(parini%verb.gt.0) write(*,'(a,es15.7)')"Velocity of CM: ", sqrt(velcm(1)**2+velcm(2)**2+velcm(3)**2)
          if(parres%verb.gt.0) call yaml_map('Velocity of CM',sqrt(velcm(1)**2+velcm(2)**2+velcm(3)**2),fmt='(es15.7)')
!Kinetic energy according to Parrinello Rahman
          rkin=0.d0
if(md_type==4) then
  rkin=vvolcur**2*latmass
  ekinlat=0.5d0*rkin
else
          do i=1,3
             rkin=rkin+(vlatcur(1,i)**2+vlatcur(2,i)**2+vlatcur(3,i)**2)*latmass
          enddo
          ekinlat=0.5d0*rkin
          rkin=ekinlat+ekinatom
endif
!          call ekin_at_lat(amass,latmass,latcur,vposcur,vlatcur,ekinatom,ekinlat,f0,md_type,nat)
          rkin=ekinlat+ekinatom
!if(parini%verb.gt.0) write(*,'(a,4(1x,es15.7))') " # Torquenrm, Ekin, Enthalpy: ",torquenrm, ekinlat,ekinatom,enthalpy
if(parini%verb.gt.0) then
    call yaml_mapping_open('Torquenrm info',flow=.true.)
    call yaml_map('Torquenrm',torquenrm,fmt='(es15.7)')
    call yaml_map('Ekin',ekinatom,fmt='(es15.7)')
    call yaml_map('Enthalpy',enthalpy,fmt='(es15.7)')
    call yaml_map('ekinlat',ekinlat,fmt='(es15.7)')
    call yaml_mapping_close()
endif
!Update counter
          enmin2=enmin1
          enmin1=en0000

!Here we perform the force call
!****************************************************************************************************************        
!****************************************************************************************************************        
!Get the current acell_in and rprim_in, and also the atomic positions
if(md_type==4) then
       volpred_1_3=volpred**(1.d0/3.d0)
       latvec_in=latvec0*volpred_1_3
       latpred=latvec_in
       xcart=pospred*volpred_1_3
       call rxyz_cart2int(latvec_in,xred_in,xcart,parini%nat)
else
       xred_in=pospred
       vel_in=0.d0
       latvec_in=latpred
endif
!       getwfk=.true.
       if(parini%usewf_geopt) then
           getwfk=.true.
       else
           getwfk=.false.
       endif
       if(iprec==2.and.multiprec.and.fmax.lt.tolmxf_switch) then
           getwfk=.false.
           iprec=1
           enthalpy_min=1.d10
       endif
       write(fn4,'(i4.4)') itime
       sock_extra_string="FIRE"//trim(fn4)
       call get_energyandforces_single(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,getwfk)
!Single point calculation finished. Now returning to MD part. All output variables are now in *_in
!****************************************************************************************************************        

!FIRE: check for convergence
call convcheck(parini,parini%nat,latvec_in,fcart_in,strten_in,parini%target_pressure_habohr,parini%paropt_geopt%strfact,fmax,fmax_at,fmax_lat,parini%paropt_geopt%fmaxtol,iexit)
!!!!Compute maximal component of forces, EXCLUDING any fixed components
!!! fmax=0.0d0
!!! do iat=1,nat
!!!   do i=1,3
!!!!     if (dtsets(1)%iatfix(i,iat) /= 1) then
!!!       if( abs(fcart_in(i,iat)) >= fmax ) fmax=abs(fcart_in(i,iat))
!!!!     end if
!!!   end do
!!! end do
!!! strtarget=0.d0
!!! strtarget(1:3)=-target_pressure_habohr
!!! dstr(:)=strten_in(:)-strtarget(:)
!!!!Eventually take into account the stress
!!! do istr=1,6
!!!     if(abs(dstr(istr))*strfact >= fmax ) fmax=abs(dstr(istr))*strfact
!!! end do
!!! iexit=0
!!! if(fmax.lt.tolmxf) iexit=1


!For velocity verlet of cell
        if(options==1) then
           call acceleration_fire(pressure_md,accpospred,acclatpred,accvolpred,vpospred,&
                &vlatcur,vvolcur,strten_in,fcart_in,latpred,amass,latmass,f0inv,md_type,parini%nat)
           vlatpred(:,:)=vlatcur(:,:)+0.5d0*dt*(acclatpred+acclatcur)
           if(md_type==4) vvolpred=0.d0! vvolpred=vvolcur+0.5d0*dt*(accvolpred+accvolcur)
        elseif(options==3) then
!Corrector part of Beeman. Note that a fixed number of iterations are used!!!
           call ekin_at_lat_andersen(amass,latmass,latpred,vpospred,vlatpred,vvolpred,ekinatom_prev,ekinlat_prev,f0,md_type,parini%nat)
           do i=1,5
             call acceleration_fire(pressure_md,accpospred,acclatpred,accvolpred,vpospred,&
                  &vlatpred,vvolpred,strten_in,fcart_in,latpred,amass,latmass,f0inv,md_type,parini%nat)
             do iat=1,parini%nat
             vpospred(:,iat)=vposcur(:,iat)+dt/6.d0*(2.d0*accpospred(:,iat)+5.d0*accposcur(:,iat)-accposprev(:,iat))
             enddo
             vlatpred=vlatcur+dt/6.d0*(2.d0*acclatpred+5.d0*acclatcur-acclatprev)
             if(md_type==4) vvolpred=0.d0! vvolpred=vvolcur+dt/6.d0*(2.d0*accvolpred+5.d0*accvolcur-accvolprev)
             call ekin_at_lat_andersen(amass,latmass,latpred,vpospred,vlatpred,vvolpred,ekinatom,ekinlat,f0,md_type,parini%nat)
if(parini%verb.gt.1)write(* ,'(a,i5,1x,es15.5,es15.5)') " # Beeman: Corrector relative convergence: it, atoms, lattice: ",&
                  & i,abs(ekinatom-ekinatom_prev)/ekinatom,abs(ekinlat-ekinlat_prev)/max(ekinlat,1.d-100)
!             write(67,'(a,i5,1x,es15.5,es15.5)') " # Beeman: Corrector relative convergence: it, atoms, lattice: ",&
!                  & i,abs(ekinatom-ekinatom_prev)/ekinatom,abs(ekinlat-ekinlat_prev)/max(ekinlat,1.d-100)
             ekinlat_prev=ekinlat
             ekinatom_prev=ekinatom
           enddo
        endif
!call acceleration(pressure_md,accpospred,acclatpred,accvolpred,vpospred,vlatpred,vvolpred,strten_in,fcart_in,latcur,amass,latmass,f0inv,md_type,nat)!BUG?
call acceleration_fire(pressure_md,accpospred,acclatpred,accvolpred,vpospred,&
     &vlatpred,vvolpred,strten_in,fcart_in,latpred,amass,latmass,f0inv,md_type,parini%nat)

!Eliminate the unnecessary atmic and cell components
call elim_fixed_at(parini,parini%nat,accpospred);call elim_fixed_at(parini,parini%nat,vpospred)
if(md_type.ne.4) call elim_fixed_lat(parini,latcur,acclatpred)
if(md_type.ne.4) call elim_fixed_lat(parini,latcur,vlatpred)
!Compute the "predicted" kinetic energies:
call ekin_at_lat_andersen(amass,latmass,latpred,vpospred,vlatpred,vvolpred,ekinatom,ekinlat,f0,md_type,parini%nat)


!MHM: Write output to file in every step***********************************
       counter=real(itime,8)
       call get_enthalpy(latvec_in,etot_in,pressure,enthalpy)
       call getvol(latvec_in,vol)
       en0000=enthalpy-ent_pos_0
       if (itime >= 3 .and. enmin1 > enmin2 .and. enmin1 > en0000)  nummax=nummax+1
       if (itime >= 3 .and. enmin1 < enmin2 .and. enmin1 < en0000)  nummin=nummin+1
!       if (mpi_enreg%me == 0) write(67,'(a,i5,1x,1pe17.10)') 'FIRE ' ,itime,enthalpy
!       if (mpi_enreg%me == 0) write(*,'(a,i5,1x,1pe17.10)') ' #FIRE ',itime,enthalpy
if(parini%verb.gt.0) then
       write(fn4,'(i4.4)') itime
       filename=trim(folder)//"posgeopt."//fn4//".ascii"
       units=units
       !write(*,*) "# Writing the positions in FIRE: ",filename
       call yaml_map('Writing the positions in FIRE',trim(filename))
       call write_atomic_file_ascii(parini,filename,parini%nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &parini%char_type(1:parini%ntypat_global),parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,etot_in,pressure,enthalpy,en0000)
       if(parini%verb.ge.3) then
       filename=trim(folder)//"posgeopt."//fn4//".vasp"
       call write_atomic_file_poscar(parini,filename,parini%nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &parini%char_type(1:parini%ntypat_global),parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,etot_in,pressure,enthalpy,en0000)
       endif

        call yaml_mapping_open('FIRE',flow=.true.)
        call yaml_map('iter',itime,fmt='(i4)')
        call yaml_map('enth',enthalpy,fmt='(es20.12)')
        call yaml_map('fmax',fmax,fmt='(es10.3)')
        call yaml_map('nstep',nstep,fmt='(i4)')
        call yaml_map('iprec',iprec,fmt='(i2)')
        call yaml_map('dt',dt,fmt='(es9.2)')
        call yaml_map('fmax_at',fmax_at,fmt='(es10.3)')
        call yaml_map('fmax_lat',fmax_lat,fmt='(es10.3)')
        call yaml_map('rkin',rkin,fmt='(es9.2)')
        call yaml_map('P',P,fmt='(es9.1)')
        call yaml_map('P_at',P_at,fmt='(es9.1)')
        call yaml_map('P_lat',P_lat,fmt='(es9.1)')
        call yaml_mapping_close()
       !write(*,'(a,i4,6(1x,es17.8),3(1x,es9.2),1x,i4,1x,i4)') " # GEOPT FIRE    ",&
       !       &itime,enthalpy, fmax, fmax_at,fmax_lat,rkin,P,P_at,P_lat,dt,nstep,iprec
endif
        if(iexit==1) then
            call yaml_mapping_open('FIRE converged',flow=.true.)
            call yaml_map('itime',itime)
            call yaml_map('enthalpy',enthalpy,fmt='(es20.12)')
            call yaml_map('fmax',fmax,fmt='(es12.5)')
            call yaml_mapping_close()
          !write(*,'(a,i4,2(1x,es25.15))') " # FIRE converged", itime,enthalpy,fmax
          exit
        endif
!*********************************************************************

!        econs_max=max(econs_max,rkin+e_rxyz)
!        econs_min=min(econs_min,rkin+e_rxyz)
        econs_max=max(econs_max,rkin+enthalpy)
        econs_min=min(econs_min,rkin+enthalpy)
!        write(6,*) istep,e_rxyz-e_pos,nummax,nummin

!Update the variables for next iteration
        accposprev=accposcur
        accposcur=accpospred
        acclatprev=acclatcur
        acclatcur=acclatpred
        velmatcur=velmatpred
        flatcur=flatpred
        fposcur=fpospred
        poscur=pospred
        latcur=latpred

        accvolprev=accvolcur
        accvolcur=accvolpred
if(md_type==4)   volcur=volpred

!        vvolcur=vvolpred
!        vlatcur=vlatpred
!        vposcur=vpospred
!FIRE velocity update ****************
        call fpos_flat(parini,pressure_md,fall(:,1:parini%nat),fall(:,parini%nat+1:parini%nat+3),strten_in,fcart_in,latpred,md_type)
!        do iat=1,nat
!          fall(:,iat)=accpospred(:,iat)*amass(iat)
!        enddo
!        fall(:,nat+1:nat+3)=acclatpred*latmass 
        vall=0.d0
        vall(:,1:parini%nat)=vpospred(:,1:parini%nat)         

!        if(.not.fixlat(7)) then
        vall(:,parini%nat+1:parini%nat+3)=vlatpred(:,:)     
        if(md_type==4) then
            vall(:,parini%nat+1:parini%nat+3)=0.d0
            vall(1,parini%nat+1)=vvolpred
        endif
!           else
!           vall(:,nat+1:nat+3)=0.d0
!        endif    
        P=0.d0
        do i=1,parini%nat+3
         P=P+fall(1,i)*vall(1,i)+fall(2,i)*vall(2,i)+fall(3,i)*vall(3,i) 
         if(i==parini%nat) then
            P_at=P
!            write(16,'(a,e15.7,$)') " P_at=",P_at
         elseif(i==parini%nat+3) then
            P_lat=P-P_at
!            write(16,'(a,e15.7,$)') " P_lat=",P_lat
         endif       
        enddo
!        write(16,'(a,e15.7)') " P=",P

!Slight Modification of FIRE: Here we consider P as two different contributions, one from the atoms and one from the lattice
!        P=min(P_at,P_lat)
        fallnorm=0.d0
        do i=1,parini%nat+3
         fallnorm=fallnorm+fall(1,i)**2+fall(2,i)**2+fall(3,i)**2
        enddo       
        fall=fall/sqrt(fallnorm)
        vallnorm=0.d0
        do i=1,parini%nat+3
         vallnorm=vallnorm+vall(1,i)**2+vall(2,i)**2+vall(3,i)**2
        enddo       
        vallnorm=sqrt(vallnorm)
        vposcur=(1.d0-alpha)*vpospred+alpha*fall(:,1:parini%nat)*vallnorm
        vlatcur=(1.d0-alpha)*vlatpred+alpha*fall(:,parini%nat+1:parini%nat+3)*vallnorm
        vvolcur=(1.d0-alpha)*vvolpred+alpha*fall(1,parini%nat+1)*vallnorm
!Feedback on energy if gradient fails
        if(enthalpy.gt.enthalpy_min+1.d-3) then
             P=-1.d0
             enthalpy_min=enthalpy
        endif


         if(fmax.lt.cellfix_switch.and..not.cellfix_done.and.(.not.(any(parini%fixlat).or.any(parini%fixat).or.confine.ge.1))) then
!Only perform the cell correction once, presumably close to the end of the optimization run
             cellfix_done=.true.
             call correct_latvec(latcur,poscur,parini%nat,parini%correctalg,latvec_io)
             if(latvec_io.ne.0) then
               max_kpt=.false.
               ka1=0;kb1=0;kc1=0
             endif
         endif 

!FIRE timestep update ****************
!        if(P.le.0.d0.or.(fmax.lt.cellfix_switch.and..not.cellfix_done.and.(.not.(any(fixlat).or.any(fixat).or.confine.ge.1)))) then
        if(P.le.0.d0.or.latvec_io.ne.0) then
         latvec_io=0
         nstep=0
         dt=max(dt*fdec,dtmin)
!         accposcur=0.d0
!         acclatcur=0.d0
!         accvolcur=0.d0
         vposcur=0.d0
         vlatcur=0.d0
         vvolcur=0.d0
         accposprev=0.d0
         acclatprev=0.d0
         accvolprev=0.d0
         alpha=alphastart
         if((multiprec.and.itime.ge.parini%paropt_geopt%nit/2).or.&
          &(fmax.lt.1.0d0*tolmxf_switch)) max_kpt=.true.
         elseif(P.gt.0.d0 .and. nstep.gt.Nmin) then
           dt=min(dt*finc,dtmax)
!           alpha=max(alpha*falpha,0.1d0)!alpha*falpha
           alpha=alpha*falpha
        endif 
        nstep=nstep+1
!*************************************
!Feedback on cell
        if(md_type==4) then
          if(f_lat*fall(1,parini%nat+1).gt.0) then
             alpha_lat=alpha_lat*1.05d0
          else
             alpha_lat=alpha_lat*0.7d0
          endif
        f_lat=fall(1,parini%nat+1)
        endif


        call acceleration_fire(pressure_md,accposcur,acclatcur,accvolcur,vposcur,&
             &vlatcur,vvolcur,strten_in,fcart_in,latcur,amass,latmass,f0inv,md_type,parini%nat)
!Eliminate the unnecessary atmic and cell components
        call elim_fixed_at(parini,parini%nat,accposcur);call elim_fixed_at(parini,parini%nat,vposcur)
if(md_type.ne.4) call elim_fixed_lat(parini,latcur,acclatcur)
if(md_type.ne.4) call elim_fixed_lat(parini,latcur,vlatcur)


     enthalpy_min=min(enthalpy,enthalpy_min)
     enddo 
!FIRE
    call yaml_mapping_open('EXIT FIRE',flow=.true.)
    call yaml_map('itime',itime)
    call yaml_map('enthalpy',enthalpy,fmt='(es20.12)')
    call yaml_map('etot_in',etot_in,fmt='(es20.12)')
    call yaml_mapping_close()
     !write(*,'(a,i5,es15.7,es15.7)') ' # EXIT FIRE ',itime,enthalpy,etot_in
     max_kpt=.false.
    call yaml_sequence_close()

end subroutine GEOPT_FIRE_MHM 

!**********************************************************************************************

subroutine acceleration_fire(pressure,accpos,acclat,accvol,vpos,vlat,vvol,strten,fcart,latvec,amass,latmass,f0inv,md_type,nat) 
implicit none
integer:: iat,i,j,md_type,nat
real(8),dimension(3,nat):: accpos,vpos,fcart,fpos
real(8),dimension(3,3)  :: acclat,vlat,latvec,tmplat,pressure,a,velmat,sigma,lattrans,latdottrans,gdot,g,ginv,gtot,str_matrix
real(8),dimension(3,3)  :: term1,term2,term3,term4,term5,term5_1,term5_2,sigmatrans,f0inv
real(8):: amass(nat),latmass,crossp(3),strten(6),vol,vpostmp(3),volvel,trace3
real(8):: accvol,vvol,vol_1_3
!Get volume
           a=latvec
           vol= a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)+&
                a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)
!md_type=1
!Convert cartesian forces to reduced forces on atoms
        call invertmat(latvec,tmplat,3)
        do iat=1,nat
          fpos(:,iat)=matmul(tmplat,fcart(:,iat))
        enddo
!Get full stress matrix
           str_matrix(1,1)=strten(1)
           str_matrix(2,2)=strten(2)
           str_matrix(3,3)=strten(3)
           str_matrix(1,2)=strten(6)
           str_matrix(2,1)=strten(6)
           str_matrix(1,3)=strten(5)
           str_matrix(3,1)=strten(5)
           str_matrix(2,3)=strten(4)
           str_matrix(3,2)=strten(4)
!Update velocity part of stress tensor
           velmat=0.d0
!           do iat=1,nat
!              vpostmp=matmul(latvec,vpos(:,iat))
!              do i=1,3
!                 do j=1,3
!                 velmat(i,j)=velmat(i,j)+vpostmp(i)*vpostmp(j)*amass(iat) !Added mass on Oct 23 2015
!                 enddo
!              enddo
!           enddo
if(md_type.ge.1.and.md_type.le.3) then
!Compute sigma
        call cross_product(latvec(:,2),latvec(:,3),crossp); sigma(:,1)=crossp
        call cross_product(latvec(:,3),latvec(:,1),crossp); sigma(:,2)=crossp
        call cross_product(latvec(:,1),latvec(:,2),crossp); sigma(:,3)=crossp
!Compute the atomic acceleration
!For this we first need to calculate the matrix g^-1*dg/dt=ginv*gdot=gtot
          do i=1,3
            lattrans(:,i)=latvec(i,:)
            latdottrans(:,i)=vlat(i,:)
          enddo
          gdot=matmul(latdottrans,latvec)+matmul(lattrans,vlat)
          g=matmul(lattrans,latvec)
          call invertmat(g,ginv,3)
          gtot=matmul(ginv,gdot)
!Total acceleration
          do iat=1,nat
          accpos(:,iat)=1.d0/amass(iat)*fpos(:,iat)! - matmul(gtot,vpos(:,iat))
!          if(fixlat(7))  accpos(:,iat)=1.d0/amass(iat)*fpos(:,iat)
          enddo
elseif(md_type==4) then
!Andersen
!Be careful: the positions are not the real positions but scaled with vol**(1/3)
          vol_1_3=vol**(1.d0/3.d0)
          do iat=1,nat
              !accpos(:,iat)=(fcart(:,iat)/amass(iat)-2.d0*vpos(:,iat)*vvol/(vol*3.d0))/vol_1_3
              accpos(:,iat)=(fcart(:,iat)/amass(iat))/vol_1_3!-2.d0*vpos(:,iat)*vvol/(vol*3.d0))/vol_1_3
          enddo
!Update velocity part of stress tensor
           velmat=0.d0
!           do iat=1,nat
!              vpostmp=vol_1_3*vpos(:,iat)
!              do i=1,3
!                 do j=1,3
!                 velmat(i,j)=velmat(i,j)+vpostmp(i)*vpostmp(j)*amass(iat) !Added mass on Oct 23 2015
!                 enddo
!              enddo
!           enddo
endif
velmat=0.d0
if(md_type==1) then
!Compute the acceleration of the cell
        acclat=velmat/vol-str_matrix
!Here the pressure is applied
        acclat=acclat-pressure
!Scale with lattice mass
        acclat=acclat/latmass
!Multiply with sigma from left
        acclat=matmul(acclat,sigma)

elseif(md_type==2) then
!Fist term, same as in PR******************************
!Compute the acceleration of the cell
        term1=velmat/vol-str_matrix
!Here the pressure is applied
        term1=term1-pressure
!Scale with lattice volume and mass
        term1=term1/vol/latmass
!Combine it with the cell
        term1=matmul(term1,latvec)

!Second term*********************************************
        sigmatrans(:,1)=sigma(1,:)
        sigmatrans(:,2)=sigma(2,:)
        sigmatrans(:,3)=sigma(3,:)
        term2=matmul(sigmatrans,vlat)
        volvel=term2(1,1)+term2(2,2)+term2(3,3)
        term2=-2.d0*volvel/vol*vlat

!Third term**********************************************
!This term was taken from the moldyn manual. It is slightly different in the cleveland paper!!!
!        term3=matmul(sigmatrans,sigma)
!        term3=matmul(term3,latdottrans)
!        term3=matmul(term3,vlat)
!From cleveland
        term3=matmul(vlat,sigmatrans)
        term3=matmul(term3,sigma)
        term3=matmul(term3,latdottrans)
        trace3=term3(1,1)+term3(2,2)+term3(3,3)
        term3=1.d0/vol**2*trace3*latvec 

!Fourth term*********************************************
        term4=matmul(vlat,sigmatrans)
        term4=matmul(term4,vlat)
        term4=1.d0/vol*term4
 
!Fifth term**********************************************
        term5_1=matmul(vlat,sigmatrans)
        term5_1=matmul(term5_1,sigma)
        term5_1=matmul(term5_1,latdottrans)
        term5_2=matmul(sigma,latdottrans)
        term5_2=matmul(term5_2,vlat)
        term5_2=matmul(term5_2,sigmatrans)
        term5=1.d0/vol**2*matmul((term5_1-term5_2),latvec)
!Sum
        acclat=term1+term2+term3+term4+term5
elseif(md_type==3) then
!Compute the acceleration of the cell
        acclat=velmat/vol-str_matrix
!Here the pressure is applied
        acclat=acclat-pressure
!Scale with lattice mass
        acclat=acclat/latmass*vol**(4.d0/3.d0)
!Multiply with sigma from left
        acclat=matmul(acclat,sigma)
        acclat=matmul(acclat,f0inv)
elseif(md_type==4) then
!Andersen MD
!Compute the acceleration of the cell
        acclat=velmat/vol-str_matrix
!Here the pressure is applied
        acclat=acclat-pressure
!Scale with lattice mass
        acclat=acclat/latmass
!Compute the hydrostatic pressure stuff
        accvol=(acclat(1,1)+acclat(2,2)+acclat(3,3))/3.d0
else
stop "Wrong option in MD"
endif
end subroutine
