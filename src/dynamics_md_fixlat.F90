!**********************************************************************************************
subroutine MD_fixlat(parini,parres,latvec_in,xred_in,fcart_in,strten_in,vel_in,etot_in,iprec,counter,folder)
 use mod_parini, only: typ_parini
 use global, only: units
 use defs_basis
 use interface_code
implicit none
    type(typ_parini), intent(in):: parini
    type(typ_parini), intent(inout):: parres
    integer:: iat,iprec,istep
    real(8):: latvec_in(3,3), xred_in(3,parini%nat),fcart_in(3,parini%nat),vel_in(3,parini%nat), strten_in(6), etot_in, counter
    real(8):: rxyz(3,parini%nat),fxyz(3,parini%nat),fxyz_old(3,parini%nat),vxyz(3,parini%nat),amass(parini%nat)
    character(len=4) :: fn4
    real(8) :: e0,at1,at2,at3
    character(40):: filename,folder
    integer :: nummax,nummin, istepnext
    real(8) :: enmin1, enmin2, en0000, econs_max, econs_min, devcon
    logical:: getwfk
    real(8):: pressure,int_pressure_gpa,energy,rkin,rkin_0,dt,dt_ratio
    if((all(parini%fixlat(1:6)).and..not.parini%fixlat(7)).or.parini%bc==2) then
       continue
    else
       write(*,*) "This routine only intended for fixed cell MD calculation"
       stop
    endif
    pressure=parini%target_pressure_habohr
    dt=parres%dtion_md

    !C initialize positions,velocities, forces
    call rxyz_int2cart(latvec_in,xred_in,rxyz,parini%nat)


    !C inner (escape) loop
    nummax=0
    nummin=0
    enmin1=0.d0
    en0000=0.d0
    econs_max=-1.d100
    econs_min=1.d100
    istepnext=5

!Assign masses to each atom (for MD)
    do iat=1,parini%nat
      amass(iat)=amu_emass*parini%amu(parini%typat_global(iat))
      write(*,'(a,i5,2(1x,es15.7))') " # MD: iat, AMU, EM: ", iat, parini%amu(parini%typat_global(iat)),amass(iat)
    enddo

!INITIAL STEP, STILL THE SAME STRUCTURE AS INPUT
       getwfk=.false.
       call get_energyandforces_single(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,getwfk)
       write(*,*) "Energy",etot_in
       istep=0
       energy=etot_in
       fxyz=fcart_in
       fxyz_old=fcart_in
       vxyz=vel_in
       call elim_fixed_at(parini,parini%nat,vxyz)
       call elim_fixed_at(parini,parini%nat,fxyz)
       call elim_fixed_at(parini,parini%nat,fxyz_old)
       rkin=0.d0
       do iat=1,parini%nat
          rkin=rkin+amass(iat)*(vxyz(1,iat)**2+vxyz(2,iat)**2+vxyz(3,iat)**2)
       enddo
       rkin=rkin*.5d0
       rkin_0=rkin
       int_pressure_gpa=1.d0/3.d0*(strten_in(1)+strten_in(2)+strten_in(3))*HaBohr3_GPa
       if(parini%verb.gt.0) then
       write(*,'(a,i5,1x,3(1x,1pe17.10),3(1x,i2))') ' # MD: it,energy,ekinat,pressure,nmax,nmin,mdmin ',&
             &istep,energy,rkin,int_pressure_gpa,nummax,nummin,parres%mdmin
       write(fn4,'(i4.4)') istep
       filename=trim(folder)//"posmd."//fn4//".ascii"
       units=units
       write(*,*) "# Writing the positions in MD:",filename
       call write_atomic_file_ascii(parini,filename,parini%nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &parini%char_type(1:parini%ntypat_global),parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,etot_in,pressure,etot_in,etot_in)
       endif
!*********************************************************************
    e0 = etot_in

    md_loop: do istep=1,parini%nmd_dynamics

!C      Evolution of the system according to 'VELOCITY VERLET' algorithm
       do iat=1,parini%nat
          rxyz(:,iat)=rxyz(:,iat)+dt*vxyz(:,iat)
          rxyz(:,iat)=rxyz(:,iat)+0.5d0*dt*dt/amass(iat)*fxyz(:,iat)
         ! call daxpy(3*nat,dt,vxyz(1,1),1,rxyz(1,1),1)
         ! call daxpy(3*nat,0.5d0*dt*dt,fxyz(1,1),1,rxyz(1,1),1)
       enddo


       enmin2=enmin1
       enmin1=en0000
       if(istep.ne.1.and.parini%usewf_md) then
           getwfk=.true.
       else
           getwfk=.false.
       endif
       call rxyz_cart2int(latvec_in,xred_in,rxyz,parini%nat)
       call get_energyandforces_single(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,getwfk)
       fxyz=fcart_in
       energy=etot_in
!       call call_bigdft(runObj, outs, nproc,iproc,infocode)
       en0000=energy-e0
       if (istep >= 3 .and. enmin1 > enmin2 .and. enmin1 > en0000)  nummax=nummax+1
       if (istep >= 3 .and. enmin1 < enmin2 .and. enmin1 < en0000)  nummin=nummin+1
       do iat=1,parini%nat
          at1=fxyz(1,iat)
          at2=fxyz(2,iat)
          at3=fxyz(3,iat)
          !C Evolution of the velocities of the system
!          if (.not. atoms%lfrztyp(iat)) then
             vxyz(1,iat)=vxyz(1,iat) + (.5d0*dt/amass(iat)) * (at1 + fxyz_old(1,iat))
             vxyz(2,iat)=vxyz(2,iat) + (.5d0*dt/amass(iat)) * (at2 + fxyz_old(2,iat))
             vxyz(3,iat)=vxyz(3,iat) + (.5d0*dt/amass(iat)) * (at3 + fxyz_old(3,iat))
!          end if
          !C Memorization of old forces
          fxyz_old(1,iat) = at1
          fxyz_old(2,iat) = at2
          fxyz_old(3,iat) = at3
       end do
       call elim_fixed_at(parini,parini%nat,vxyz)
       call elim_fixed_at(parini,parini%nat,fxyz)
       call elim_fixed_at(parini,parini%nat,fxyz_old)
       rkin=0.d0
       do iat=1,parini%nat
          rkin=rkin+amass(iat)*(vxyz(1,iat)**2+vxyz(2,iat)**2+vxyz(3,iat)**2)
       enddo
       rkin=rkin*.5d0

!!!   if (atoms%astruct%geocode == 'S') then 
!!!      call fixfrag_posvel_slab(iproc,atoms%astruct%nat,rcov,atoms%astruct%rxyz,vxyz,2)
!!!   else if (atoms%astruct%geocode == 'F') then
!!!     if (istep == istepnext) then 
!!!           call fixfrag_posvel(iproc,atoms%astruct%nat,rcov,atoms%astruct%rxyz,vxyz,2,occured)
!!!        if (occured) then 
!!!          istepnext=istep+4
!!!        else
!!!          istepnext=istep+1
!!!        endif
!!!     endif
!!!   endif
       econs_max=max(econs_max,rkin+energy)
       econs_min=min(econs_min,rkin+energy)
       devcon=econs_max-econs_min
       counter=real(istep,8)
       int_pressure_gpa=1.d0/3.d0*(strten_in(1)+strten_in(2)+strten_in(3))*HaBohr3_GPa
       if(parini%verb.gt.0) then
       write(*,'(a,i5,1x,3(1x,1pe17.10),3(1x,i2))') ' # MD: it,energy,ekinat,pressure,nmax,nmin,mdmin ',&
             &istep,energy,rkin,int_pressure_gpa,nummax,nummin,parres%mdmin
       write(fn4,'(i4.4)') istep
       filename=trim(folder)//"posmd."//fn4//".ascii"
       units=units
       write(*,*) "# Writing the positions in MD: ",filename
       call write_atomic_file_ascii(parini,filename,parini%nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &parini%char_type(1:parini%ntypat_global),parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,etot_in,pressure,etot_in,etot_in)
       endif
       if (nummin.ge.parres%mdmin) then
          if (nummax.ne.nummin) &
               write(*,*) '# WARNING: nummin,nummax',nummin,nummax
             write(*,*) " MD finished: exiting!"
          exit md_loop
       endif

    end do md_loop
    !C MD stopped, now do relaxation

    !  if (iproc == 0) write(67,*) 'EXIT MD',istep
    
    ! adjust time step to meet precision criterion
!Minimum number of steps per crossed minimum is 15, average should be parini%nit_per_min
 if(parini%auto_dtion_md) then
    if(parini%energy_conservation) then
        devcon=devcon/(3*parini%nat-3)
        if (devcon/rkin_0.lt.1.d-2) then
           parres%dtion_md=parres%dtion_md*1.05d0
        else
           parres%dtion_md=parres%dtion_md*(1.d0/1.05d0)
        endif
      write(*,'(a,es10.2,es10.2,a,es10.2)') " # MD: energy conservation :",devcon,rkin_0,&
      &", parres%dtion_md set to: ",parres%dtion_md
    else 
       dt_ratio=real(istep,8)/real(nummin,8) 
       if(dt_ratio.lt.real(parini%nit_per_min,8)) then
         parres%dtion_md=parres%dtion_md*1.d0/1.1d0
       else
         parres%dtion_md=parres%dtion_md*1.1d0 
       endif
       parres%dtion_md=min(parres%dtion_md,dt*dt_ratio/15.d0)
       write(*,'(3(a,es10.2))') " # MD: steps per minium: ",dt_ratio,&
      &", parres%dtion_md set to: ",parres%dtion_md,", upper boundary: ",dt*dt_ratio/15.d0 
     endif
  endif
   

!MD stopped, now do relaxation
     write(*,'(a,i5,es15.7)') ' # EXIT MD ',istep,etot_in
    
  END SUBROUTINE

