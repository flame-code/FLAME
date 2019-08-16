module mod_md_rbmd
implicit none
interface
function A_omega(omega)
    implicit none
    real(8),dimension(4,4):: A_omega
    real(8)::omega(3)
end function A_omega
end interface
end module mod_md_rbmd



subroutine MD_MHM_ROT(parini,parres,latvec_in,xred_in,xred_cm_in,xcart_mol,quat_in,fcart_in,strten_in,&
                      &vel_in,vel_cm_in,vel_lat_in,l_in,vvol_in,etot_in,&
                      &masstot,intens,inprin,inaxis,lhead,llist,nmol,iprec,counter,folder)
 use global, only: units
 use defs_basis
 use interface_code
 use mod_parini, only: typ_parini
 use yaml_output
implicit none
 type(typ_parini), intent(in):: parini
 type(typ_parini), intent(inout):: parres
 real(8) ::latvec_in(3,3),xred_in(3,parini%nat),vel_in(3,parini%nat),fcart_in(3,parini%nat),etot_in,strten_in(6),vel_lat_in(3,3),vvol_in
 real(8):: amass(nmol)
 real(8):: pressure_ener
 real(8):: pressure_md(3,3)
 real(8):: pressure
 real(8):: vol
 real(8):: vol0
 real(8):: volprev
!Rigid stuff
 real(8),dimension(3,nmol):: xred_cm_in
 real(8),dimension(3,nmol):: fcart_cm
 real(8),dimension(3,nmol):: torque
 real(8),dimension(3,parini%nat) :: xcart_mol
 real(8),dimension(3,nmol):: l_in
 real(8),dimension(3,nmol):: vel_cm_in
 real(8),dimension(4,nmol):: quat_in
 real(8),dimension(4,nmol):: quatcur
 real(8),dimension(4,nmol):: quatpred
 real(8),dimension(3,3,nmol):: intens
 real(8),dimension(3,nmol):: inprin
 real(8),dimension(nmol):: masstot
 real(8),dimension(3,3,nmol):: inaxis
 integer,dimension(nmol):: lhead
 integer,dimension(parini%nat):: llist
 integer:: nmol
!Rigid stuff
 real(8),dimension(3,nmol):: xcart
 real(8),dimension(3,nmol):: fposcur
 real(8),dimension(3,nmol):: accposcur
 real(8),dimension(3,nmol):: accpospred
 real(8),dimension(3,nmol):: accposprev
 real(8),dimension(3,nmol):: fpospred
 real(8),dimension(3,nmol):: vpospred
 real(8),dimension(3,nmol):: poscur
 real(8),dimension(3,nmol):: vxyz
 real(8),dimension(3,nmol):: vposcur
 real(8),dimension(3,parini%nat):: pospred
 real(8),dimension(3,parini%nat):: dxred
 real(8),dimension(3,3):: dlatvec
 real(8),dimension(3,3):: latvec
 real(8),dimension(3,3):: latvec0
 real(8),dimension(3,3):: latinv
 real(8),dimension(3,3):: unitmat
 real(8),dimension(3,3):: elmatcur
 real(8),dimension(3,3):: acclatcur
 real(8),dimension(3,3):: acclatpred
 real(8),dimension(3,3):: acclatprev
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
 if(any(parini%fixat).and.nmol.ne.parini%nat) stop "MD with fixed atoms not implemented yet!"
 write(*,*) "# PARAMETERS: Ha_eV,kb_HaK,amu_emass -->",Ha_eV,kb_HaK,amu_emass


!Assign masses to each atom (for MD)
 do iat=1,nmol
   amass(iat)=masstot(iat)
   write(*,'(a,i5,1(1x,es15.7))') " # MD: imol, AMU, EM: ", iat,amass(iat)
 enddo
!******************************************************************
!NEW VERISON OF MD: Directly implemented, simplest Parrinello-Rahman and other MD

!Here we split the routine in the Abinit native part and my new implentation
md_type=parres%md_algo
if(md_type==1) then
 write(*,'(a)') ' # Entering standalone Parrinello Rahman MD '
elseif(md_type==2) then
 write(*,'(a)') ' # Entering standalone Cleveland MD '
elseif(md_type==3) then
 write(*,'(a)') ' # Entering standalone Wentzcovitch MD '
elseif(md_type==4) then
 write(*,'(a)') ' # Entering standalone Andersen MD '
endif

!The "reduced" coordinates in Andersen are quite different from the ones in PR
!Set temporary variables, initially
  vxyz(:,:)=vel_cm_in(:,:)
  pressure=parres%target_pressure_habohr
  unitmat=0.d0
  do i=1,3
    unitmat(i,i)=1.d0
  enddo
  latmass0=parini%bmass*amu_emass !This means we use the barostat mass as the lattice mass (in ELECTRON MASS)
  vlat=vel_lat_in  !The initial cell velocity
  itime=0
  dt=parres%dtion_md

!Set options=1 for Velocity Verlet of cell dynamics
!Set options=2 for Normal Verlet of cell dynamics
!Set options=3 for Beeman integration scheme, corrector-predictor
  options=parres%md_integrator
  if(options.lt.1.or.options.gt.3) stop "Wrong algo option"

!MD-type: 1 for PR and 2 for Cleveland and 3 for Wentzcovitch and 4 for Andersen
  md_type=parres%md_algo
  if(md_type.lt.1.or.md_type.gt.4) stop "Wrong integrator option"

  write(*,'(a,i3,a,i3)') " # MD Algorithm: ",md_type, ", MD Integrator: ",options

!Now we run my implementation of MD
pressure_md=0.d0
pressure_md(1,1)=1.d0;pressure_md(2,2)=1.d0;pressure_md(3,3)=1.d0
pressure_ener=0.d0;pressure_md=pressure_md*pressure  !Here the pressure is not passed to the energyandforces, so we move on the ENERGY surface

!Transform the initial velocities given in cartesian coordinates into velocities of the intenal coordinates
       latvec=latvec_in
       call getvol(latvec,vol)
if(md_type==4) then
       vol_1_3=vol**(1.d0/3.d0)
       vol_1_3_inv=1.d0/vol_1_3
       do iat=1,nmol
          vposcur(:,iat)=vxyz(:,iat)*vol_1_3_inv
       enddo
else
       call invertmat(latvec,latinv,3)
       do iat=1,nmol
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


!!ATTENTION: THIS IS WITH TAKING INTO ACCOUNT THE STRESS TENSOR GIVEN IN INPUT

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
        call rxyz_int2cart(latvec,xred_cm_in,xcart,nmol)
        poscur=xcart*vol_1_3_inv           !call rxyz_cart2int(latvec,poscur,rxyz,nat)
        volcur=vol
        vvolcur=vvol_in
else
        poscur=xred_cm_in           !call rxyz_cart2int(latvec,poscur,rxyz,nat)
        latcur=latvec
        vlatcur=vlat
endif

!Initiallize quat
quatcur=quat_in
quatpred=quat_in

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
!Here we need to add rotational energies!!!
if(md_type==4) then
  latcur=latvec0*vol_1_3
  rkin=0.d0
  do iat=1,nmol
     vposcurtmp=vposcur(:,iat)*vol_1_3
     rkin=rkin+amass(iat)*(vposcurtmp(1)**2+vposcurtmp(2)**2+vposcurtmp(3)**2)
  enddo
  ekinatom=0.5d0*rkin
  rkin=vvolcur**2*latmass
  ekinlat=0.5d0*rkin
else
  rkin=0.d0
  do iat=1,nmol
     vposcurtmp=matmul(latcur,vposcur(:,iat))
     rkin=rkin+amass(iat)*(vposcurtmp(1)**2+vposcurtmp(2)**2+vposcurtmp(3)**2)
  enddo
  ekinatom=0.5d0*rkin
  rkin=0.d0
  do i=1,3
     rkin=rkin+(vlatcur(1,i)**2+vlatcur(2,i)**2+vlatcur(3,i)**2)*latmass
  enddo
  ekinlat=0.5d0*rkin
endif


!MHM: Write output to file in every step***********************************
!INITIAL STEP, STILL THE SAME STRUCTURE AS INPUT
       call expand_rigid(latvec_in,xred_cm_in,quat_in,xcart_mol,lhead,llist,parini%nat,nmol,xred_in)
       write(*,*) "Pressure, Energy",pressure,etot_in
       call get_enthalpy(latvec_in,etot_in,pressure,enthalpy)
       ent_pos_0=enthalpy
       en0000=enthalpy-ent_pos_0
       write(*,'(a,i5,1x,4(1x,1pe17.10),3(1x,i2))') ' # MD: it,enthalpy,pV,ekinat,ekinlat,nmax,nmin,mdmin ',&
             &itime,enthalpy,pressure*vol,ekinatom,ekinlat,nummax,nummin,parres%mdmin
       write(fn4,'(i4.4)') itime
       filename=trim(folder)//"posmd."//fn4//".ascii"
       units=units
       !write(*,*) "# Writing the positions in MD:",filename
       call yaml_map('Writing the positions in MD',trim(filename))
       call write_atomic_file_ascii(parres,filename,parini%nat,units,xred_in,latvec_in,fcart_in,strten_in,&
       &parini%char_type(1:parini%ntypat_global),parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,etot_in,pressure,enthalpy,en0000)
!*********************************************************************


call get_fcm_torque(fcart_cm,torque,fcart_in,quat_in,xcart_mol,lhead,llist,parini%nat,nmol)
call acceleration(pressure_md,accposcur,acclatcur,accvolcur,vposcur,vlatcur,vvolcur,&
                  strten_in,fcart_cm,latcur,amass,latmass,f0inv,md_type,nmol)
        accposprev=accposcur
        acclatprev=acclatcur
        accvolprev=accvolcur

!Eliminate the unnecessary atmic and cell components
call elim_fixed_at(parini,nmol,accposcur);call elim_fixed_at(parini,nmol,accposprev)
if(md_type.ne.4) call elim_fixed_lat(parini,latcur,acclatcur)
if(md_type.ne.4) call elim_fixed_lat(parini,latcur,acclatprev)
call elim_fixed_at(parini,nmol,vposcur)
if(md_type.ne.4) call elim_fixed_lat(parini,latcur,vlatcur)


        do itime=1,parres%nmd_dynamics

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
          do iat=1,nmol
              dxred(:,iat)=dt*vposcur(:,iat) + 0.5d0*dt*dt*accposcur(:,iat)  !0.5 was missing before
!             pospred(:,iat)=poscur(:,iat) + dt*vposcur(:,iat) + dt*dt*accposcur(:,iat)
          enddo
          call propagate(parini,nmol,poscur,latcur,dxred,dlatvec,pospred,latpred)
          do iat=1,nmol
             vpospred(:,iat)=(pospred(:,iat)-poscur(:,iat))/dt   !Update dummy variable vposcur for next iteration
          enddo
            if(md_type==4) then
               volpred=volcur+dt*vvolcur+0.5d0*dt*dt*accvolcur
            endif
          elseif(options==2) then
!Instead if normal verlet is used:
              dlatvec(:,:)=dt*vlatcur + dt*dt*acclatcur
          do iat=1,nmol
              dxred(:,iat)=dt*vposcur(:,iat) + dt*dt*accposcur(:,iat)
          enddo
          call propagate(parini,nmol,poscur,latcur,dxred,dlatvec,pospred,latpred)
          vlatpred(:,:)=(latpred-latcur)/dt
          do iat=1,nmol
              vpospred(:,iat)=(pospred(:,iat)-poscur(:,iat))/dt   !Update dummy variable vposcur for next iteration
          enddo
            if(md_type==4) then
               volpred=volcur+dt*vvolcur+dt*dt*accvolcur
               vvolpred=(volpred-volcur)/dt
            endif
          elseif(options==3) then
!Predictor part of Beeman for cell
              dlatvec(:,:)=dt*vlatcur + dt*dt/6.d0*(4.d0*acclatcur(:,:)-acclatprev(:,:))
          do iat=1,nmol
!Predictor part of Beeman
             dxred(:,iat)=dt*vposcur(:,iat) + dt*dt/6.d0*(4.d0*accposcur(:,iat)-accposprev(:,iat))             
          enddo
          call propagate(parini,nmol,poscur,latcur,dxred,dlatvec,pospred,latpred)
          vlatpred(:,:)=vlatcur(:,:) + 0.5d0*(3.d0*acclatcur(:,:)-acclatprev(:,:))  
          do iat=1,nmol
             vpospred(:,iat)=vposcur(:,iat)+0.5d0*(3.d0*accposcur(:,iat)-accposprev(:,iat))
          enddo
            if(md_type==4) then
               volpred=volcur+dt*vvolcur+dt*dt/6.d0*(4.d0*accvolcur-accvolprev)
               vvolpred=volcur+0.5d0*(3.d0*accvolcur-accvolprev)
            endif
          endif


!Rotational part to predict the correct quatpred
!Implement the case for assymetric top
!call rbmd_driver






if(md_type==4) then
  do iat=1,nmol
     vposcurtmp=vposcur(:,iat)*vol_1_3
     rkin=rkin+amass(iat)*(vposcurtmp(1)**2+vposcurtmp(2)**2+vposcurtmp(3)**2)
     velcm=velcm+vposcurtmp*amass(iat)
  enddo
else
  do iat=1,nmol
             vposcurtmp=matmul(latcur,vposcur(:,iat))
             rkin=rkin+amass(iat)*(vposcurtmp(1)**2+vposcurtmp(2)**2+vposcurtmp(3)**2)
             velcm=velcm+vposcurtmp*amass(iat)
  enddo
endif

!Kinetic energy of atoms
          ekinatom=0.5d0*rkin
          !if(parres%verb.gt.0) write(*,'(a,es15.7)')"Velocity of CM: ", sqrt(velcm(1)**2+velcm(2)**2+velcm(3)**2)
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
!if(parres%verb.gt.0) write(*,'(a,3(1x,es15.7))') " # Torquenrm, Ekin, Enthalpy: ",torquenrm, rkin,enthalpy
if(parres%verb.gt.0) then
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
!Get the current acell_in and rprim_in, and also the atomic positions
if(md_type==4) then
       volpred_1_3=volpred**(1.d0/3.d0)
       latvec_in=latvec0*volpred_1_3
       latpred=latvec_in
       xcart=pospred*volpred_1_3
       call rxyz_cart2int(latvec_in,xred_cm_in,xcart,nmol)
       call expand_rigid(latvec_in,xred_cm_in,quatpred,xcart_mol,lhead,llist,parini%nat,nmol,xred_in)
else
       latvec_in=latpred
       xred_cm_in=pospred
       call expand_rigid(latvec_in,pospred,quatpred,xcart_mol,lhead,llist,parini%nat,nmol,xred_in)
endif
       if(itime==1)  then
           getwfk=.false.
       elseif(itime.ne.1.and.parini%usewf_md) then
           getwfk=.true.
       else
           getwfk=.false.
       endif
       call get_energyandforces_single(parini,parres,latvec_in,xred_cm_in,fcart_in,strten_in,etot_in,iprec,getwfk)
       call get_fcm_torque(fcart_cm,torque,fcart_in,quatpred,xcart_mol,lhead,llist,parini%nat,nmol)
!Single point calculation finished. Now returning to MD part. All output variables are now in *_in
!****************************************************************************************************************        

!For velocity verlet of cell
        if(options==1) then
           call acceleration(pressure_md,accpospred,acclatpred,accvolpred,vpospred,&
                &vlatcur,vvolcur,strten_in,fcart_cm,latpred,amass,latmass,f0inv,md_type,nmol) 
           vlatpred(:,:)=vlatcur(:,:)+0.5d0*dt*(acclatpred+acclatcur)
           if(md_type==4) vvolpred=vvolcur+0.5d0*dt*(accvolpred+accvolcur)
        elseif(options==3) then
!Corrector part of Beeman. Note that a fixed number of iterations are used!!!
           call ekin_at_lat_andersen(amass,latmass,latpred,vpospred,vlatpred,vvolpred,ekinatom_prev,ekinlat_prev,f0,md_type,nmol)
           do i=1,5
             call acceleration(pressure_md,accpospred,acclatpred,accvolpred,vpospred,&
                  &vlatpred,vvolpred,strten_in,fcart_cm,latpred,amass,latmass,f0inv,md_type,nmol) 
             do iat=1,nmol
             vpospred(:,iat)=vposcur(:,iat)+dt/6.d0*(2.d0*accpospred(:,iat)+5.d0*accposcur(:,iat)-accposprev(:,iat))
             enddo
             vlatpred=vlatcur+dt/6.d0*(2.d0*acclatpred+5.d0*acclatcur-acclatprev)
             if(md_type==4) vvolpred=vvolcur+dt/6.d0*(2.d0*accvolpred+5.d0*accvolcur-accvolprev)
             call ekin_at_lat_andersen(amass,latmass,latpred,vpospred,vlatpred,vvolpred,ekinatom,ekinlat,f0,md_type,nmol)
if(parres%verb.gt.1)write(* ,'(a,i5,1x,es15.5,es15.5)') " # Beeman: Corrector relative convergence: it, atoms, lattice: ",&
                  & i,abs(ekinatom-ekinatom_prev)/ekinatom,abs(ekinlat-ekinlat_prev)/max(ekinlat,1.d-100)
!             write(67,'(a,i5,1x,es15.5,es15.5)') " # Beeman: Corrector relative convergence: it, atoms, lattice: ",&
!                  & i,abs(ekinatom-ekinatom_prev)/ekinatom,abs(ekinlat-ekinlat_prev)/max(ekinlat,1.d-100)
             ekinlat_prev=ekinlat
             ekinatom_prev=ekinatom
           enddo 
        endif  
call acceleration(pressure_md,accpospred,acclatpred,accvolpred,vpospred,&
     &vlatpred,vvolpred,strten_in,fcart_cm,latpred,amass,latmass,f0inv,md_type,nmol)

!Eliminate the unnecessary atmic and cell components
call elim_fixed_at(parini,nmol,accpospred);call elim_fixed_at(parini,nmol,vpospred)
if(md_type.ne.4) call elim_fixed_lat(parini,latcur,acclatpred)
if(md_type.ne.4) call elim_fixed_lat(parini,latcur,vlatpred)

!Compute the "predicted" kinetic energies:
call ekin_at_lat_andersen(amass,latmass,latpred,vpospred,vlatpred,vvolpred,ekinatom,ekinlat,f0,md_type,nmol)


!MHM: Write output to file in every step***********************************
       counter=real(itime,8)
       call get_enthalpy(latvec_in,etot_in,pressure,enthalpy)
       call getvol(latvec_in,vol)
       volpred=vol
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
       call write_atomic_file_ascii(parres,filename,parini%nat,units,xred_in,latvec_in,fcart_in,strten_in,&
       &parini%char_type(1:parini%ntypat_global),parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,etot_in,pressure,enthalpy,en0000)
!*********************************************************************
        econs_max=max(econs_max,rkin+enthalpy)
        econs_min=min(econs_min,rkin+enthalpy)
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

     enddo 
!Adjust MD stepsize
!Minimum number of steps per crossed minimum is 15, average should be parres%nit_per_min
     if(parres%auto_dtion_md) then
       if(real(itime,8)/real(parres%mdmin,8).lt.real(parres%nit_per_min,8)) then
         parres%dtion_md=parres%dtion_md*1.d0/1.1d0
       else
         parres%dtion_md=parres%dtion_md*1.1d0 
       endif
     parres%dtion_md=min(parres%dtion_md,real(itime,8)*dt/(real(parres%mdmin,8)*15.d0))
     write(*,'(3(a,es10.2))') " # MD: steps per minium: ",real(itime,8)/real(parres%mdmin,8),&
           &", parres%dtion_md set to: ",parres%dtion_md,", upper boundary: ",real(itime,8)*dt/(real(parres%mdmin,8)*15.d0) 
     endif
!MD stopped, now do relaxation
     write(*,'(a,i5,es15.7,es15.7)') ' # EXIT MD ',itime,enthalpy,etot_in

end subroutine


subroutine rbmd_symasym_s1(T_t,L_t,dt,L_til_t)
!J.Chem.Phys. 110,7,3291,Matubayasi
!Symmetric top: we assume Ix=Iy
!First step: equation 15
!We construct the vector of L^tilde at T=t: L_til_t
!Input is the L_t and T_t, the angular momentum and
!the torque at time T=t
implicit none
real(8),intent(in) :: T_t(3), L_t(3),dt
real(8),intent(out):: L_til_t(3)
!integer:: nmol
L_til_t=L_t+0.5d0*dt*T_t
end subroutine

!************************************************************************************

subroutine rbmd_sym_s23(Inprin,L_til_t,dt,L_til_t5,L_til_t10)
!J.Chem.Phys. 110,7,3291,Matubayasi
!Symmetric top: we assume Ix=Iy
!Second and third steps: equation 17 and 18
!We construct the vector of L^tilde at T=t+0.5dt: L_til_t5
!Then, We construct the vector of L^tilde at T=t+dt: L_til_t10
!Input is the principal inertia tensor Inprin, symmetric Ix=Iz, and L^tilde at T=t, L_til_t
implicit none
real(8),intent(in) :: Inprin(3), L_til_t(3),dt
real(8),intent(out):: L_til_t10(3),L_til_t5(3)
real(8):: Ixz,Iangle,ca,sa
!integer:: nmol,imol
!do imol=1,nmol
  Ixz=1.d0/Inprin(1)-1.d0/Inprin(2)
  Iangle=0.5*dt*Ixz*L_til_t(3)
  ca=cos(Iangle)
  sa=sin(Iangle)

!Eq.17
  L_til_t5(1)=ca*L_til_t(1)-sa*L_til_t(2)
  L_til_t5(2)=sa*L_til_t(1)-sa*L_til_t(2)
  L_til_t5(3)=   L_til_t(3)
!Eq.18
  L_til_t10(1)=ca*L_til_t5(1)-sa*L_til_t5(2)
  L_til_t10(2)=sa*L_til_t5(1)-sa*L_til_t5(2)
  L_til_t10(3)=    L_til_t(3)
!enddo
end subroutine

!************************************************************************************

subroutine rbmd_symasym_s4(Inprin,L_til_t5,quat_t,dt,quat_t10)
    use mod_md_rbmd
!J.Chem.Phys. 110,7,3291,Matubayasi
!Symmetric top: we assume Ix=Iy
!Fourth step: equations 3 and 9/10
!We construct the new quaternion at time T=t+dt: quat_t10
!We thus need to first obtain omega_tilde at T=t+0.5*dt by using L_tilde at T=t+0.5*dt
implicit none
real(8),intent(in) :: Inprin(3),L_til_t5(3),dt,quat_t(4)
real(8),intent(out):: quat_t10(4)
real(8):: omega_til_t5(3),expmat(4,4),vnrm
real(8):: ca,sa,angle
!real(8),dimension(4,4):: A_omega
!integer:: nmol,imol
!real(8):: A_omega

!do imol=1,nmol
!Construct omega according to equation 3
  omega_til_t5(1)=L_til_t5(1)/Inprin(1)
  omega_til_t5(2)=L_til_t5(2)/Inprin(2)
  omega_til_t5(3)=L_til_t5(3)/Inprin(3)
!Do the crazy stuff at equation 10
  expmat=0.d0
  vnrm=sqrt(omega_til_t5(1)+omega_til_t5(2)+omega_til_t5(3))  
  angle=0.5d0*dt*vnrm
  ca=cos(angle)
  sa=sin(angle)
  expmat(1,1)=ca
  expmat(2,2)=ca
  expmat(3,3)=ca
  expmat(4,4)=ca
  expmat=expmat+sa/vnrm*A_omega(omega_til_t5(1:3))
  quat_t10(:)=matmul(expmat,quat_t)
!enddo
end subroutine rbmd_symasym_s4

!************************************************************************************

function A_omega(omega)
real(8),dimension(4,4):: A_omega
real(8)::omega(3)
A_omega=0.d0
A_omega(2,1)=omega(1)
A_omega(3,1)=omega(2)
A_omega(4,1)=omega(3)
A_omega(3,2)=-omega(3)
A_omega(4,2)=omega(2)
A_omega(4,3)=-omega(1)

A_omega(1,2)=-A_omega(2,1)
A_omega(1,3)=-A_omega(3,1)
A_omega(1,4)=-A_omega(4,1)
A_omega(2,3)=-A_omega(3,2)
A_omega(2,4)=-A_omega(4,2)
A_omega(3,4)=-A_omega(4,3)
end function A_omega

!************************************************************************************

subroutine rbmd_symasym_s5(T_t10,L_til_t10,dt,L_t10)
!J.Chem.Phys. 110,7,3291,Matubayasi
!Symmetric top: we assume Ix=Iy
!Fifth  step, fourth step according to the manuscript: equation 19
!We construct the vector of L at T=t+dt: L_t10
!Input is the L_til_t10 (L^tilde at T=t+dt) and the torque at T=t+dt: T_t10
implicit none
real(8),intent(in) :: T_t10(3), L_til_t10(3),dt
real(8),intent(out):: L_t10(3)
!integer:: nmol
L_t10=L_til_t10+0.5d0*dt*T_t10
end subroutine

!************************************************************************************

subroutine rbmd_asym_s23(Inprin,L_til_t,dt,L_til_t5,L_til_t10)
!J.Chem.Phys. 110,7,3291,Matubayasi
!Asymmetric top: Ix.ne.Iy.ne.Iz
!Second and third steps: equation 21 and 22
!We construct the vector of L^tilde at T=t+0.5dt: L_til_t5
!Then, We construct the vector of L^tilde at T=t+dt: L_til_t10
!Input is the principal inertia tensor Inprin, symmetric Ix=Iz, and L^tilde at T=t, L_til_t
implicit none
real(8),intent(in) :: Inprin(3), L_til_t(3),dt
real(8),intent(out):: L_til_t10(3),L_til_t5(3)
real(8):: Izy,Izx,Ianglezy,Ianglezx,cazy,sazy,cazx,sazx,Lz1,Lz2
!integer:: nmol,imol
!do imol=1,nmol
  Izy=1.d0/Inprin(3)-1.d0/Inprin(2)
  Izx=1.d0/Inprin(3)-1.d0/Inprin(1)
  Ianglezy=0.5*dt*Izy*L_til_t(2)
  cazy=cos(Ianglezy)
  sazy=sin(Ianglezy)
!Eq.21
  L_til_t5(1)= cazy*L_til_t(1)-sazy*L_til_t(3)
  Lz1        =-sazy*L_til_t(1)+cazy*L_til_t(3)
  Ianglezx=0.5*dt*Izx*L_til_t5(1)
  cazx=cos(Ianglezx)
  sazx=sin(Ianglezx)
  L_til_t5(2)= cazx*L_til_t(2)-sazx*Lz1
  L_til_t5(3)= sazx*L_til_t(2)+cazx*Lz1
!Eq.22
  L_til_t10(2)=   cazx*L_til_t5(2)-sazx*L_til_t5(3)
  Lz2         =  -sazx*L_til_t5(2)+cazx*L_til_t5(3)
  Ianglezy=0.5*dt*Izy*L_til_t10(2)
  cazy=cos(Ianglezy)
  sazy=sin(Ianglezy)
  L_til_t10(1)= cazy*L_til_t5(1)+sazy*Lz2
  L_til_t10(3)=-sazy*L_til_t5(1)+cazy*Lz2
!enddo
end subroutine

!************************************************************************************

subroutine rbmd_driver(quat_t,T_t,L_t,quat_t10,T_t10,L_t10,dt,inprin,&
           &fragsize,symtop,nmol)
!Currently implemented only to propagate the assymtric tops with fragsize.ge.3
implicit none
real(8),intent(in) :: inprin(3,nmol),L_t(3,nmol),T_t(3,nmol),quat_t(4,nmol),dt,T_t10(3,nmol)
integer,intent(in) :: nmol,fragsize(nmol)
logical,intent(in) :: symtop(nmol)
real(8),intent(out):: L_t10(3,nmol),quat_t10(4,nmol)
real(8) :: L_til_t(3,nmol),L_til_t5(3,nmol),L_til_t10(3,nmol)
integer :: imol
do imol=1,nmol
    if(.not.symtop(imol).and.fragsize(imol).ge.3) then
    call rbmd_symasym_s1(T_t(:,imol),L_t(:,imol),dt,L_til_t(:,imol))
    call rbmd_asym_s23(Inprin(:,imol),L_til_t(:,imol),dt,L_til_t5(:,imol),L_til_t10(:,imol))
    call rbmd_symasym_s4(Inprin(:,imol),L_til_t5(:,imol),quat_t(:,imol),dt,quat_t10(:,imol))
    call rbmd_symasym_s5(T_t10(:,imol),L_til_t10(:,imol),dt,L_t10(:,imol))
    endif
enddo
end subroutine

!************************************************************************************

subroutine expand_rigid(latvec,xred_cm,quat,xcart_mol,lhead,llist,nat,nmol,xred_in)
!This routine will produce the desired, expanded form of xred_in which can be directly fed into 
!the single point energy routine. Input are quat and the xred_cm, the center of masses of each molecule
!in the cell
implicit none
real(8),intent(in):: latvec(3,3),xred_cm(3,nmol),quat(4,nmol),xcart_mol(3,nat)
real(8),intent(out):: xred_in(3,nat)
real(8):: rotmat(3,3),xcart_tmp(3,nat),cmass(3,nmol)
integer:: nat,nmol,iat,imol,llist(nat),lhead(nmol)
call rxyz_int2cart(latvec,xred_cm,cmass,nmol)
do imol=1,nmol
   call quat2rotmat(rotmat,quat(:,imol))
   iat=lhead(imol)
   do while (iat.ne.0)
     xcart_tmp(:,iat)=matmul(rotmat,xcart_mol(:,iat))
     xcart_tmp(:,iat)=xcart_tmp(:,iat)+cmass(:,imol)
     iat=llist(iat)
   enddo
enddo
!Produce reduced coordinates
call rxyz_cart2int(latvec,xred_in,xcart_tmp,nat)
end subroutine expand_rigid
