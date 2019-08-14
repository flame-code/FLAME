!**********************************************************************************************
subroutine MD_ANDERSEN_MHM     (parini,parres,latvec_in,xred_in,fcart_in,strten_in,vel_in,vel_lat_in,vvol_in,etot_in,iprec,counter,folder)
 use global, only: units
 use defs_basis
 use interface_code
 use mod_parini, only: typ_parini
 use yaml_output
!Implements constant cell shape MD from Andersen,
!I havent checked the equations of motion myself, taken at the moment from
!http://www.grs-sim.de/cms/upload/Carloni/Presentations/Strodel4.pdf
!The coordinates of interest are the xred of atoms and the cell volume vol
implicit none
 type(typ_parini), intent(in):: parini
 type(typ_parini), intent(inout):: parres
 real(8) :: latvec_in(3,3),xred_in(3,parini%nat),vel_in(3,parini%nat),fcart_in(3,parini%nat),etot_in,strten_in(6),vel_lat_in,vvol_in
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
 real(8),dimension(3,parini%nat):: accposcur
 real(8),dimension(3,parini%nat):: accpospred
 real(8),dimension(3,parini%nat):: accposprev
 real(8),dimension(3,parini%nat):: fpospred
 real(8),dimension(3,parini%nat):: vpospred
 real(8),dimension(3,parini%nat):: poscur
 real(8),dimension(3,parini%nat):: vxyz
 real(8),dimension(3,parini%nat):: vposcur
 real(8),dimension(3,parini%nat):: pospred
 real(8),dimension(3,parini%nat):: dxred
 real(8),dimension(3,3):: latvec
 real(8),dimension(3,3):: latcur
 real(8),dimension(3,3):: latvec0
 real(8),dimension(3,3):: latinv
 real(8),dimension(3,3):: unitmat
 real(8),dimension(3,3):: elmatcur
 real(8),dimension(3,3):: acclatcur
 real(8),dimension(3,3):: acclatpred
 real(8),dimension(3,3):: acclatprev
 real(8),dimension(3,3):: dlatvec
 real(8),dimension(3,3):: flatcur
 real(8),dimension(3,3):: flatpred
 real(8):: accvolcur
 real(8):: accvolpred
 real(8):: accvolprev
 real(8):: fvolcur
 real(8):: fvolpred
 real(8),dimension(3,3):: velmatcur
 real(8),dimension(3,3):: tmplat
 real(8),dimension(3,3):: tmplatt
 real(8),dimension(3,3):: tpred
 real(8),dimension(3,3):: f0
 real(8),dimension(3,3):: f0inv
 real(8),dimension(3,3):: vlatcur
 real(8),dimension(3,3):: vlatpred
 real(8):: volcur
 real(8):: volpred
 real(8):: vvolcur
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
 real(8):: vvolpred
 real(8):: volpred_1_3
 real(8):: vol_1_3_inv
 real(8):: vol_1_3
 real(8),dimension(3,3):: velmatpred
 real(8),dimension(3,3):: vlat
 real(8),dimension(3,3):: str_matrix
 real(8):: vposcurtmp(3)
 real(8):: vpospredtmp(3)
 real(8):: crossp(3)
 real(8):: velcm(3)
 real(8):: latmass
 real(8):: latmassinv
 real(8):: latmass0
 real(8):: ekinatom
 real(8):: ekinatom_prev
 real(8):: ekinlat
 real(8):: ekinlat_prev
 real(8):: rkin
 real(8):: enthalpy
 real(8):: enmin1
 real(8):: enmin2
 real(8):: ent_pos_0
 real(8):: en0000
 real(8):: e_rxyz
 real(8):: econs_max
 real(8):: econs_min
 real(8):: torquenrm 
 real(8):: ecut_tmp
 real(8):: toldff_tmp
 real(8):: counter
 real(8):: dt
 real(8):: dt_ratio
 integer:: i
 integer:: j
 integer:: iat
 integer:: itime
 integer:: nummax
 integer:: nummin
 integer:: iprec 
 integer:: options 
 integer:: md_type
 logical:: getwfk

 character(40)::filename,folder
 character(4) ::fn4
 write(*,*) "# PARAMETERS: Ha_eV,kb_HaK,amu_emass -->",Ha_eV,kb_HaK,amu_emass
!MD-type: 1 for PR and 2 for Cleveland and 3 for Wentzcovitch
  md_type=parres%md_algo


call yaml_sequence_open('atoms info in MD')
!Assign masses to each atom (for MD)
 do iat=1,parini%nat
   amass(iat)=amu_emass*parini%amu(parini%typat_global(iat))
   !write(*,'(a,i5,2(1x,es15.7))') " # MD: iat, AMU, EM: ", iat, parini%amu(parini%typat_global(iat)),amass(iat)
    call yaml_sequence(advance='no')
    call yaml_mapping_open('atom',flow=.true.)
    call yaml_map('iat',iat,fmt='(i5)')
    call yaml_map('AMU',parini%amu(parini%typat_global(iat)),fmt='(es15.7)')
    call yaml_map('EM',amass(iat),fmt='(es15.7)')
    call yaml_mapping_close()
 enddo
call yaml_sequence_close()
!******************************************************************
!Here we split the routine in the Abinit native part and my new implentation
 write(*,'(a)') ' # Entering standalone ANDERSEN MD  '
!The "reduced" coordinates in Andersen are quite different from the ones in PR
!Set temporary variables, initially
  vxyz(:,:)=vel_in(:,:)
  pressure=parini%target_pressure_habohr
  unitmat=0.d0
  do i=1,3
    unitmat(i,i)=1.d0
  enddo
  latmass0=parini%bmass*amu_emass !This means we use the barostat mass as the lattice mass (in ELECTRON MASS)
  vlat=vel_lat_in  !The initial cell volume velocity

  itime=0
  dt=parres%dtion_md

!Set options=1 for Velocity Verlet of cell dynamics
!Set options=2 for Normal Verlet of cell dynamics
!Set options=3 for Beeman integration scheme, corrector-predictor
  options=parres%md_integrator
  if(options.lt.1.or.options.gt.3) stop "Wrong algo option"


  write(*,'(a,i3)') " # MD Algorithm: ANDERSEN, MD Integrator: ",options

!Now we run my implementation of ANDERSEN MD
pressure_md=0.d0
pressure_md(1,1)=1.d0;pressure_md(2,2)=1.d0;pressure_md(3,3)=1.d0
pressure_ener=0.d0;pressure_md=pressure_md*pressure  !Here the pressure is not passed to the energyandforces, so we move on the ENERGY surface

!Transform the initial velocities given in cartesian coordinates into velocities of the intenal coordinates
       latvec=latvec_in
       call getvol(latvec,vol)
       vol_1_3=vol**(1.d0/3.d0) 
       vol_1_3_inv=1.d0/vol_1_3 
       do iat=1,parini%nat
          vposcur(:,iat)=vxyz(:,iat)*vol_1_3_inv
       enddo

!Keep track of volume
!In ANDERSEN, latvec0 is normed to have volume 0, so computing latvec0*vol^(1/3) will always give the current lattice vector
        vol0=vol
        latvec0=latvec/vol**(1.d0/3.d0)
        volprev=vol

!!This is for the original volume-mass factor
        latmass=latmass0
        latmassinv=1.d0/latmass
        write(*,*) "Latmasse",latmass

!Initialize internal coordinates and lattice vectors
        call rxyz_int2cart(latvec,xred_in,xcart,parini%nat)
        poscur=xcart*vol_1_3_inv           !call rxyz_cart2int(latvec,poscur,rxyz,nat)
        volcur=vol
        vvolcur=vvol_in

!MHM: initialize variable to track the minima/maxima***********
    write(*,*) '# MINHOP start MD'
    nummax=0
    nummin=0
    enmin1=0.d0
    en0000=0.d0
    econs_max=-1.d100
    econs_min=1.d100
!**************************************************************

!Compute initial kinetic energies
  latcur=latvec0*vol_1_3
  rkin=0.d0
  do iat=1,parini%nat
     vposcurtmp=vposcur(:,iat)*vol_1_3
     rkin=rkin+amass(iat)*(vposcurtmp(1)**2+vposcurtmp(2)**2+vposcurtmp(3)**2)
  enddo
  ekinatom=0.5d0*rkin
  rkin=vvolcur**2*latmass
  ekinlat=0.5d0*rkin


!MHM: Write output to file in every step***********************************
!INITIAL STEP, STILL THE SAME STRUCTURE AS INPUT
       write(*,*) "Pressure, Energy",pressure,etot_in
       call get_enthalpy(latvec_in,etot_in,pressure,enthalpy)
       ent_pos_0=enthalpy
       en0000=enthalpy-ent_pos_0
       write(*,'(a,i5,1x,4(1x,1pe17.10),3(1x,i2))') ' MD: it,enthalpy,pV,ekinat,ekinlat,nmax,nmin,mdmin ',&
             &itime,enthalpy,pressure*vol,ekinatom,ekinlat,nummax,nummin,parres%mdmin
       write(fn4,'(i4.4)') itime
       filename=trim(folder)//"posmd."//fn4//".ascii"
       units=units
       !write(*,*) "# Writing the positions in MD:",filename
       call yaml_map('Writing the positions in MD',trim(filename))
       call write_atomic_file_ascii(parini,filename,parini%nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &parini%char_type(1:parini%ntypat_global),parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,etot_in,pressure,enthalpy,en0000)
!*********************************************************************




call acceleration(pressure_md,accposcur,acclatcur,accvolcur,vposcur,&
     &vlatcur,vvolcur,strten_in,fcart_in,latcur,amass,latmass,f0inv,md_type,parini%nat) 
        accposprev=accposcur
!        acclatprev=acclatcur
        accvolprev=accvolcur

!Eliminate the unnecessary atmic and cell components
call elim_fixed_at(parini,parini%nat,accposcur);call elim_fixed_at(parini,parini%nat,accposprev)!;call elim_fixed_lat(latcur,acclatcur);call elim_fixed_lat(latcur,acclatprev)
call elim_fixed_at(parini,parini%nat,vposcur)!;call elim_fixed_lat(latcur,vlatcur)


        do itime=1,parini%nmd_dynamics

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
              dxred(:,iat)=dt*vposcur(:,iat) + 0.5d0*dt*dt*accposcur(:,iat)
!             pospred(:,iat)=poscur(:,iat) + dt*vposcur(:,iat) + dt*dt*accposcur(:,iat)
          enddo
          call propagate(parini,parini%nat,poscur,latcur,dxred,dlatvec,pospred,latpred)
          do iat=1,parini%nat
             vpospred(:,iat)=(pospred(:,iat)-poscur(:,iat))/dt   !Update dummy variable vposcur for next iteration
          enddo
            if(md_type==4) then
               volpred=volcur+dt*vvolcur+0.5d0*dt*dt*accvolcur       
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
               volpred=volcur+dt*vvolcur+dt*dt*accvolcur       
               vvolpred=(volpred-volcur)/dt
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
               volpred=volcur+dt*vvolcur+dt*dt/6.d0*(4.d0*accvolcur-accvolprev)       
               vvolpred=volcur+0.5d0*(3.d0*accvolcur-accvolprev)
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


!!!!Now calculate steps for the lattice according to the velocity verlet algorithm
!!!          if(options==1) then
!!!          elseif(options==2) then
!!!          elseif(options==3) then
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
endif
          rkin=ekinlat+ekinatom
!if(parini%verb.gt.0) write(*,'(a,3(1x,es15.7))') " # Torquenrm, Ekin, Enthalpy: ",torquenrm, rkin,enthalpy
if(parini%verb.gt.0) then
    call yaml_mapping_open('Torquenrm info',flow=.true.)
    call yaml_map('Torquenrm',torquenrm,fmt='(es15.7)')
    call yaml_map('Ekin',rkin,fmt='(es15.7)')
    call yaml_map('Enthalpy',enthalpy,fmt='(es15.7)')
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
       if(itime==1)  then
           getwfk=.false.
       elseif(itime.ne.1.and.parini%usewf_md) then
           getwfk=.true.
       else
           getwfk=.false.
       endif
       call get_energyandforces_single(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,getwfk)
!Single point calculation finished. Now returning to MD part. All output variables are now in *_in
!****************************************************************************************************************        

!For velocity verlet of cell
        if(options==1) then
           call acceleration(pressure_md,accpospred,acclatpred,accvolpred,vpospred,&
                &vlatcur,vvolcur,strten_in,fcart_in,latpred,amass,latmass,f0inv,md_type,parini%nat) 
           vlatpred(:,:)=vlatcur(:,:)+0.5d0*dt*(acclatpred+acclatcur)
           if(md_type==4) vvolpred=vvolcur+0.5d0*dt*(accvolpred+accvolcur)
        elseif(options==3) then
!Corrector part of Beeman. Note that a fixed number of iterations are used!!!
           call ekin_at_lat_andersen(amass,latmass,latpred,vpospred,vlatpred,vvolpred,ekinatom_prev,ekinlat_prev,f0,md_type,parini%nat)
           do i=1,5
             call acceleration(pressure_md,accpospred,acclatpred,accvolpred,vpospred,&
                  &vlatpred,vvolpred,strten_in,fcart_in,latpred,amass,latmass,f0inv,md_type,parini%nat) 
             do iat=1,parini%nat
             vpospred(:,iat)=vposcur(:,iat)+dt/6.d0*(2.d0*accpospred(:,iat)+5.d0*accposcur(:,iat)-accposprev(:,iat))
             enddo
             vlatpred=vlatcur+dt/6.d0*(2.d0*acclatpred+5.d0*acclatcur-acclatprev)
             if(md_type==4) vvolpred=vvolcur+dt/6.d0*(2.d0*accvolpred+5.d0*accvolcur-accvolprev)
             call ekin_at_lat_andersen(amass,latmass,latpred,vpospred,vlatpred,vvolpred,ekinatom,ekinlat,f0,md_type,parini%nat)
if(parini%verb.gt.1)write(* ,'(a,i5,1x,es15.5,es15.5)') " # Beeman: Corrector relative convergence: it, atoms, lattice: ",&
                  & i,abs(ekinatom-ekinatom_prev)/ekinatom,abs(ekinlat-ekinlat_prev)/max(ekinlat,1.d-100)
        !     write(67,'(a,i5,1x,es15.5,es15.5)') " # Beeman: Corrector relative convergence: it, atoms, lattice: ",&
        !          & i,abs(ekinatom-ekinatom_prev)/ekinatom,abs(ekinlat-ekinlat_prev)/max(ekinlat,1.d-100)
             ekinlat_prev=ekinlat
             ekinatom_prev=ekinatom
           enddo 
        endif  
call acceleration(pressure_md,accpospred,acclatpred,accvolpred,vpospred,&
     &vlatpred,vvolpred,strten_in,fcart_in,latpred,amass,latmass,f0inv,md_type,parini%nat) 

!Eliminate the unnecessary atmic and cell components
call elim_fixed_at(parini,parini%nat,accpospred);call elim_fixed_at(parini,parini%nat,vpospred)!;call elim_fixed_lat(latcur,acclatpred);call elim_fixed_lat(latcur,vlatpred)

!Compute the "predicted" kinetic energies:
call ekin_at_lat_andersen(amass,latmass,latpred,vpospred,vlatpred,vvolpred,ekinatom,ekinlat,f0,md_type,parini%nat)


!MHM: Write output to file in every step***********************************
       counter=real(itime,8)
       call get_enthalpy(latvec_in,etot_in,pressure,enthalpy)
       call getvol(latvec_in,vol)
       en0000=enthalpy-ent_pos_0
       if (itime >= 3 .and. enmin1 > enmin2 .and. enmin1 > en0000)  nummax=nummax+1
       if (itime >= 3 .and. enmin1 < enmin2 .and. enmin1 < en0000)  nummin=nummin+1
!       write(67,'(a,i5,1x,4(1x,1pe17.10),3(1x,i2))') ' MD: it,enthalpy,pV,ekinat,ekinlat,nmax,nmin,mdmin ',&
!             &itime,enthalpy,pressure*vol,ekinatom,ekinlat,nummax,nummin,parres%mdmin
       write(*,'(a,i5,1x,4(1x,1pe17.10),3(1x,i2))') ' # MD: it,enthalpy,pV,ekinat,ekinlat,nmax,nmin,mdmin ',&
             &itime,enthalpy,pressure*vol,ekinatom,ekinlat,nummax,nummin,parres%mdmin
       write(fn4,'(i4.4)') itime
       filename=trim(folder)//"posmd."//fn4//".ascii"
       units=units
       !write(*,*) "# Writing the positions in MD: ",filename
       call yaml_map('Writing the positions in MD',trim(filename))
       call write_atomic_file_ascii(parini,filename,parini%nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &parini%char_type(1:parini%ntypat_global),parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,etot_in,pressure,enthalpy,en0000)
!*********************************************************************

!        econs_max=max(econs_max,rkin+e_rxyz)
!        econs_min=min(econs_min,rkin+e_rxyz)
        econs_max=max(econs_max,rkin+enthalpy)
        econs_min=min(econs_min,rkin+enthalpy)
!        write(6,*) istep,e_rxyz-e_pos,nummax,nummin
       if (nummin.ge.parres%mdmin) then
          if (nummax.ne.nummin) &
               write(*,*) '# WARNING: nummin,nummax',nummin,nummax

!          write(67,*) " MD finished: exiting!"
          write(*,*) " MD finished: exiting!"
          exit
       endif


!Update the variables for next iteration
        accposprev=accposcur
        accposcur=accpospred
        acclatprev=acclatcur
        acclatcur=acclatpred
        velmatcur=velmatpred
        flatcur=flatpred
        fposcur=fpospred
        poscur=pospred
        vlatcur=vlatpred
        vposcur=vpospred
        latcur=latpred

        accvolprev=accvolcur
        accvolcur=accvolpred
        volcur=volpred
        vvolcur=vvolpred
 
!        e_rxyz=e_rxyz_pred
     enddo 
!Adjust MD stepsize
!Minimum number of steps per crossed minimum is 15, average should be parini%nit_per_min
     if(parini%auto_dtion_md) then
!       dt_ratio=real(itime,8)/real(parres%mdmin,8) !old version version
       dt_ratio=real(itime,8)/real(nummin,8) 
       if(dt_ratio.lt.real(parini%nit_per_min,8)) then
         parres%dtion_md=parres%dtion_md*1.d0/1.1d0
       else
         parres%dtion_md=parres%dtion_md*1.1d0 
       endif
     parres%dtion_md=min(parres%dtion_md,dt*dt_ratio/15.d0)
     write(*,'(3(a,es10.2))') " # MD: steps per minium: ",dt_ratio,&
           &", parres%dtion_md set to: ",parres%dtion_md,", upper boundary: ",dt*dt_ratio/15.d0 
     endif
   

!MD stopped, now do relaxation
     write(*,'(a,i5,es15.7,es15.7)') ' # EXIT MD ',itime,enthalpy,etot_in

end subroutine

!**********************************************************************************************

subroutine ekin_at_lat_andersen(amass,latmass,latvec,vpos,vlat,vvol,ekinat,ekinlat,f0,md_type,nat)
implicit none
integer:: iat,i,md_type,nat
real(8):: latvec(3,3),vpos(3,nat),vlat(3,3),ekinat,ekinlat,rkin,vposcurtmp(3),crossp(3),f0(3,3),vol,vvol
real(8):: latmass,amass(nat),lattrans(3,3),latdottrans(3,3),ekintrace(3,3),sigma(3,3),sigmatrans(3,3),vol_1_3
rkin=0.d0
if(md_type.lt.4) then
!Compute the kinetic energies:
  do iat=1,nat
     vposcurtmp=matmul(latvec,vpos(:,iat))
     rkin=rkin+amass(iat)*(vposcurtmp(1)**2+vposcurtmp(2)**2+vposcurtmp(3)**2)
  enddo
!  do i=1,3
!     rkin=rkin+(vlat(1,i)**2+vlat(2,i)**2+vlat(3,i)**2)*latmass
!  enddo
!  ekinlat=0.5d0*rkin
!For this we first need to calculate the matrix g^-1*dg/dt=ginv*gdot=gtot
          do i=1,3
            lattrans(:,i)=latvec(i,:)
            latdottrans(:,i)=vlat(i,:)
          enddo
else
  call getvol(latvec,vol)
  vol_1_3=vol**(1.d0/3.d0)
  do iat=1,nat
     vposcurtmp=vpos(:,iat)*vol_1_3
     rkin=rkin+amass(iat)*(vposcurtmp(1)**2+vposcurtmp(2)**2+vposcurtmp(3)**2)
  enddo
endif
  ekinat=0.5d0*rkin
  rkin=0.d0

if(md_type==1) then
  ekintrace=matmul(latdottrans,vlat)
  rkin=ekintrace(1,1)+ekintrace(2,2)+ekintrace(3,3)
  ekinlat=rkin*latmass*0.5d0
elseif(md_type==2) then
!sigma 
        call cross_product(latvec(:,2),latvec(:,3),crossp); sigma(:,1)=crossp
        call cross_product(latvec(:,3),latvec(:,1),crossp); sigma(:,2)=crossp
        call cross_product(latvec(:,1),latvec(:,2),crossp); sigma(:,3)=crossp
        sigmatrans(:,1)=sigma(1,:)
        sigmatrans(:,2)=sigma(2,:)
        sigmatrans(:,3)=sigma(3,:)
        ekintrace=matmul(vlat,sigmatrans)
        ekintrace=matmul(ekintrace,sigma)
        ekintrace=matmul(ekintrace,latdottrans)
  rkin=ekintrace(1,1)+ekintrace(2,2)+ekintrace(3,3)
  ekinlat=rkin*latmass*0.5d0
elseif(md_type==3) then
  call getvol(latvec,vol)
  ekintrace=matmul(latdottrans,f0)  
  ekintrace=matmul(ekintrace,vlat)  
  rkin=ekintrace(1,1)+ekintrace(2,2)+ekintrace(3,3)
  ekinlat=rkin*latmass*0.5d0/vol**(4.d0/3.d0)  
elseif(md_type==4) then
  rkin=vvol**2*latmass
  ekinlat=0.5d0*rkin
endif
!  write(*,*) "Ekinetic",ekinlat,rkin
end subroutine 


