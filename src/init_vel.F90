subroutine init_vel(parini,parres,vel,vel_lat,vel_vol,latvec,pos_red,latmass,temp,nsoften,folder)
 use defs_basis
 use mod_parini, only: typ_parini
 use yaml_output
implicit none
 type(typ_parini), intent(in):: parini
 type(typ_parini), intent(inout):: parres
 real(8) :: acell_in(3),xred_in(3,parini%nat),vel_in(3,parini%nat),fcart_in(3,parini%nat),etot_in,strten_in(6),vel_lat_in(3,3)
!*********************************************************************************************
 real(8):: vel(3,parini%nat),temp,pos_red(3,parini%nat),vcm(3),vel_vol
 integer:: i,iat,idim,nsoften
 real(8):: amass(parini%nat),s1,s2,v2gauss,vtest,rescale_vel,vel_lat(3,3),latvec(3,3),latmass
! real(8), parameter :: Ha_eV=27.21138386d0 ! 1 Hartree, in eV
! real(8), parameter :: kb_HaK=8.617343d-5/Ha_eV ! Boltzmann constant in Ha/K
! real(8), parameter :: amu_emass=1.660538782d-27/9.10938215d-31 ! 1 atomic mass unit, in electronic mass
 real(8), parameter :: temp_fac_lat=1.d-1 !This percentage of the temperature that should be given to the lattice 
 real(8):: pressure,curv0,curv,res,count_soft
 character(40):: folder

  vel=0.d0
  vel_vol=0.d0
  vel_lat=0.d0

    call yaml_comment('Initializing velocities')
    call yaml_mapping_open('PARAMETERS',flow=.true.)
    call yaml_map('Ha_eV',Ha_eV)
    call yaml_map('kb_HaK',kb_HaK)
    call yaml_map('amu_emass',amu_emass)
    call yaml_mapping_close()
    call yaml_map('Temp',temp)

  !write(*,*) "# Initializing velocities"
  !write(*,*) "# PARAMETERS: Ha_eV,kb_HaK,amu_emass -->",Ha_eV,kb_HaK,amu_emass
  !write(*,*) "# Temp",temp

call yaml_sequence_open('atoms info in init_vel')
!Assign masses to each atom (for MD)
 do iat=1,parini%nat
 amass(iat)=amu_emass*parini%amu(parini%typat_global(iat))
  !if(parini%verb.gt.0) write(*,'(a,i5,1x,es15.7,1x,es15.7)') " # iat, AMU, EM: ",iat,parini%amu(parini%typat_global(iat)),amass(iat)
   if(parini%verb.gt.0) then
    call yaml_sequence(advance='no')
    call yaml_mapping_open('atom',flow=.true.)
    call yaml_map('iat',iat,fmt='(i5)')
    call yaml_map('AMU',parini%amu(parini%typat_global(iat)),fmt='(es15.7)')
    call yaml_map('EM',amass(iat),fmt='(es15.7)')
    call yaml_mapping_close()
   endif
 end do
call yaml_sequence_close()


pressure=parini%target_pressure_habohr

if(parini%mol_soften) then
         if(any(parini%fixat).or.any(parini%fixlat)) stop "Fixed atoms or cell not yet implemented for molecular softening"
!Init atomic velocities
         call init_rotvels(parini,parini%nat,pos_red,latvec,temp,amass,vel)
         if(parini%verb.gt.0) then
         do iat=1,parini%nat
            write(*,'(a,i5,3(1x,es15.7))') " # iat, VEL: ",iat ,vel(:,iat)
         end do
         endif
!Now get some cell velocity
        if(trim(parini%rng_type)=='only_for_tests') then
            call gausdist_cell_simple(latvec,vel_lat)
        else
            call gausdist_cell(latvec,vel_lat)
        endif
!Soften the velocities of lattice
         if(.not.parini%fixlat(7)) then
         call soften_lat(parini,parres,latvec,pos_red,vel_lat,curv0,curv,res,pressure,count_soft,amass,parini%nsoften_minhopp,folder)
         call elim_fixed_lat(parini,latvec,vel_lat)
         endif
else
!Get random Gaussian distributed atomic velocities
         if(.not.(all(parini%fixat))) then
           if(trim(parini%rng_type)=='only_for_tests') then
            call gausdist_simple(parini%nat,vel,amass)
           else
            call gausdist(parini%nat,vel,amass)
           endif
!Soften the velocities of atoms
           call soften_pos(parini,parres,latvec,pos_red,vel,curv0,curv,res,pressure,count_soft,amass,parini%nsoften_minhopp,folder)
!Get rid of center-of-mass velocity, taking also into account fixed atoms (This has already been done in gausdist for all free atoms, but lets do it again...)
           if(.not.any(parini%fixat(:))) then
             s1=sum(amass(:))
             do idim=1,3
               s2=sum(amass(:)*vel(idim,:))
               vel(idim,:)=vel(idim,:)-s2/s1
             end do
           else
             if(trim(parini%potential_potential)=="lenosky_tb_lj") then
                write(*,'(a)') " Eliminating LJ atom velocities"  
                do iat=1,parini%nat
                  if(int(parini%znucl(parini%typat_global(iat))).gt.200) vel(:,iat)=0.d0
                enddo
             endif
             
             write(*,'(a)') " Eliminating fixed atom velocities"  
             s1=0.d0
             vcm=0.d0
             do iat=1,parini%nat
                if(.not.parini%fixat(iat)) then
                  s1=s1+amass(iat)
                  vcm(:)=vcm(:)+amass(iat)*vel(:,iat)
                endif
             enddo
             if (.not. (parini%nat.gt.1 .and. (parini%nat-count(parini%fixat))==1)) then !Dont eliminate center of mass if there is only one atom to move
             do iat=1,parini%nat
                vel(:,iat)=vel(:,iat)-vcm(:)/s1
             enddo
             endif
             call elim_fixed_at(parini,parini%nat,vel)
           endif
!Recompute v2gauss
           v2gauss=0.d0
           vtest=0.d0
           do iat=1,parini%nat
             do idim=1,3
               v2gauss=v2gauss+vel(idim,iat)*vel(idim,iat)*amass(iat)
               vtest=vtest+vel(idim,iat)/(3.d0*parini%nat)
             end do
           end do
!Now rescale the velocities to give the exact temperature
         rescale_vel=sqrt(3.d0*parini%nat*kb_HaK*temp/v2gauss)
         vel(:,:)=vel(:,:)*rescale_vel
         if(parini%verb.gt.0) then
             call yaml_sequence_open('velocity of atoms')
           do iat=1,parini%nat
            call yaml_sequence(advance='no')
            call yaml_map('VEL',vel(:,iat),fmt='(es15.7)')
             !write(*,'(a,i5,3(1x,es15.7))') " # iat, VEL: ",iat ,vel(:,iat)
           end do
             call yaml_sequence_close()
         endif
         endif

!Now get some cell velocity
         if(.not.(all(parini%fixlat(1:6))).and.(.not.parini%fixlat(7)).and.parini%bc.ne.2) then
           if(trim(parini%rng_type)=='only_for_tests') then
               call gausdist_cell_simple(latvec,vel_lat)
           else
               call gausdist_cell(latvec,vel_lat)
           endif
!Soften the velocities of lattice
           call soften_lat(parini,parres,latvec,pos_red,vel_lat,curv0,curv,res,pressure,count_soft,amass,parini%nsoften_minhopp,folder)
           call elim_fixed_lat(parini,latvec,vel_lat)
         endif
endif



!Rescale to get the correct "temperature" for the cell. This is chosen by a factor "temp_fac_lat"
if(parini%fixlat(7)) then !Andersen MD
         call random_number(vel_vol)
         vel_vol=(vel_vol-0.5d0)*2.d0
         vel_vol=vel_vol*sqrt(kb_HaK*temp*temp_fac_lat)
elseif(.not.all(parini%fixlat(1:6)).and.parini%bc.ne.2) then
!        Recompute v2gauss
         v2gauss=0.d0
         vtest=0.d0
         do i=1,3
           do idim=1,3
             v2gauss=v2gauss+vel_lat(idim,i)*vel_lat(idim,i)*latmass
             vtest=vtest+vel_lat(idim,i)/(3.d0*3.d0)
           end do
         end do
!        Now rescale the velocities to give the exact temperature*temp_fac_lat
         rescale_vel=sqrt(3.d0*3.d0*kb_HaK*temp*temp_fac_lat/v2gauss)
         vel_lat(:,:)=vel_lat(:,:)*rescale_vel
         call yaml_map('lat. vel.',vel_lat,fmt='(es15.7)')
         !do i=1,3
         !  write(*,'(a,i5,3(1x,es15.7))') " # lat, VEL: ",i,vel_lat(:,i)
         !end do
endif
end subroutine init_vel
