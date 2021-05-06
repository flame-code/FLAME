 subroutine soften_pos(parini,parres,latvec,pos_red0,ddcart,curv0,curv,res,pressure,count_soft,amass,nsoft,folder)
 use global, only: units
 use defs_basis
 use interface_code
 use modsocket, only: sock_extra_string
 use mod_parini, only: typ_parini
 use yaml_output
implicit none
 type(typ_parini), intent(in):: parini
 type(typ_parini), intent(inout):: parres
 real(8) :: acell_in(3),xred_in(3,parini%nat),vel_in(3,parini%nat),etot_in,strten_in(6),vel_lat_in(3,3)
!*******************************************************************
        integer:: nsoft,i,it,nit,iprec,iat
        real(8):: curv0,curv,res,pressure,count_soft,alpha
        real(8):: rxyz(3*parini%nat)
        real(8):: latvec(9),latvec_in(9)
        real(8):: ddcart(3*parini%nat)
        real(8):: rxyzcart(3*parini%nat)
        real(8):: flat(9)
        real(8):: pos_red0(3*parini%nat)
        real(8):: pos_red_in(3*parini%nat)
        real(8):: amass(parini%nat)
        real(8):: wlat(9),wlatold(9),fxyzcart(3*parini%nat)
        real(8), allocatable :: wpos(:),fxyz(:)
        real(8):: eps_dd
        real(8):: etot
        real(8):: etot0
        real(8):: fd2
        real(8):: res1 
        real(8):: res2
        real(8):: sdd
        real(8):: sdf
        real(8):: tt
        logical:: getwfk
        character(40):: filename,folder
        character(4):: fn4         
        real(8):: pos_prev(3*parini%nat),dir_prev(3*parini%nat),dir(3*parini%nat),angle,vnrm
        logical:: decrease
        call yaml_map('Entering SOFTENING routine for ATOMS, nsoft',nsoft)
        !write(*,'(a,i5)')" # Entering SOFTENING routine for ATOMS, nsoft= ",nsoft 
        !if(parini%auto_soft) write(*,'(a)')" # Automatic softening activated" 
        if(parini%auto_soft) call yaml_comment('Automatic softening activated',hfill='~')
        decrease=.false.
!        rxyzcart=rxyz        
!First transform the atomic positions from internal to external coordinates
        call rxyz_int2cart(latvec,pos_red0,rxyz,parini%nat)

        nit=nsoft
        eps_dd=5.d-1
!        alpha=3.d0    ! step size for  Si
!        alpha=1.d0    ! step size for  C 
        alpha=parres%alpha_at
        allocate(wpos(3*parini%nat),fxyz(3*parini%nat))

!         call  rxyz_int2cart(latvec,rxyz,rxyzcart,nat)
!         call  energyandforces(nat,latvec,rxyzcart,fxyzcart,flat,pressure,etot0,count_soft)                                    !
call rxyz_cart2int(latvec,pos_red_in,rxyz,parini%nat)
latvec_in=latvec
iprec=1;getwfk=.false.
       write(fn4,'(i4.4)') 0
       sock_extra_string="SOFTAT"//trim(fn4)
call get_energyandforces_single(parini,parres,latvec_in,pos_red_in,fxyz,strten_in,etot_in,iprec,getwfk)
call get_enthalpy(latvec_in,etot_in,pressure,etot0)

!         call  fxyz_cart2int(nat,fxyzcart,fxyz,latvec)

! normalize initial guess
        sdd=0.d0
        call elim_fixed_at(parini,parini%nat,ddcart)
        do i=1,3*parini%nat
        sdd=sdd+ddcart(i)**2
        enddo
        sdd=eps_dd/sqrt(sdd)
        do i=1,3*parini%nat
        ddcart(i)=sdd*ddcart(i)
        enddo

    if(parini%verb.gt.0) call yaml_sequence_open('SOFTEN atomic iterations')
    !call yaml_sequence_open('SOFTEN_LAT')
 do it=1,nit
    if(parini%verb.gt.0) call yaml_sequence(advance='no')
        do i=1,3*parini%nat
        wpos(i)=rxyz(i)+ddcart(i)
        enddo

!Here we check for direction variation during softening
        dir_prev=dir
        dir=wpos-pos_prev
        if(it.ge.3) then
           do iat=1,parini%nat  
              angle=dot_product(dir((iat-1)*3+1:(iat)*3),dir_prev((iat-1)*3+1:(iat)*3))
              vnrm=sqrt(dot_product(dir((iat-1)*3+1:(iat)*3),dir((iat-1)*3+1:(iat)*3)))+&
                   sqrt(dot_product(dir_prev((iat-1)*3+1:(iat)*3),dir_prev((iat-1)*3+1:(iat)*3)))
              angle=angle/vnrm
              if(angle.lt.0.d0) decrease=.true.
           enddo
        endif

        pos_prev=wpos


call rxyz_cart2int(latvec,pos_red_in,wpos,parini%nat)
latvec_in=latvec
iprec=1
if(parini%usewf_soften) then
    getwfk=.true.
else
    getwfk=.false.
endif
       write(fn4,'(i4.4)') it
       sock_extra_string="SOFTAT"//trim(fn4)
call get_energyandforces_single(parini,parres,latvec_in,pos_red_in,fxyz,strten_in,etot_in,iprec,getwfk)
call get_enthalpy(latvec_in,etot_in,pressure,etot)

        fd2=2.d0*(etot-etot0)/eps_dd**2

        sdf=0.d0
        sdd=0.d0
        do i=1,3*parini%nat
        sdf=sdf+ddcart(i)*fxyz(i)
        sdd=sdd+ddcart(i)*ddcart(i)
        enddo

        curv=-sdf/sdd
        if (it.eq.1) curv0=curv
        res=0.d0
        tt=0.d0

        do i=1,3*parini%nat
        tt=tt+fxyz(i)**2
        fxyz(i)=fxyz(i)+curv*ddcart(i)
        res=res+fxyz(i)**2
        enddo
        
if(parini%verb.gt.0) then
    call yaml_mapping_open('SOFTEN',flow=.true.)
    call yaml_map('res',res,fmt='(e13.5)')
    call yaml_map('it',it,fmt='(i8)')
    call yaml_mapping_close()
    !write(*,'(a,(e13.5),i5)') ' # SOFTEN: ',res,it
endif
        res=sqrt(res)
        tt=sqrt(tt)
!if(it==1) write(*,'(a,i5,4(e13.5),e18.10)')' # SOFTEN: init atomic  it,fnrm,res,curv,fd2,etot ',it,tt,res,curv,fd2,etot-etot0
if(it==1) then
    call yaml_mapping_open('SOFTEN init atomic',flow=.true.)
    call yaml_map('it',it,fmt='(i5)')
    call yaml_map('fnrm',tt,fmt='(e13.5)')
    call yaml_map('res',res,fmt='(e13.5)')
    call yaml_map('curv',curv,fmt='(e13.5)')
    call yaml_map('fd2',fd2,fmt='(e13.5)')
    call yaml_map('de',etot-etot0,fmt='(e18.10)')
    call yaml_mapping_close()
    !if(parini%verb.gt.0) call yaml_sequence_open('SOFTEN atomic iterations')
endif
if(parini%verb.gt.0) then
        !write(*,'(a,i5,4(e13.5),e18.10)') ' # SOFTEN: it,fnrm,res,curv,fd2,etot ',it,tt,res,curv,fd2,etot-etot0
        call yaml_mapping_open('SOFTEN',flow=.true.)
        call yaml_map('it',it,fmt='(i5)')
        call yaml_map('fnrm',tt,fmt='(e13.5)')
        call yaml_map('res',res,fmt='(e13.5)')
        call yaml_map('curv',curv,fmt='(e13.5)')
        call yaml_map('fd2',fd2,fmt='(e13.5)')
        call yaml_map('de',etot-etot0,fmt='(e18.10)')
        call yaml_mapping_close()
        write(fn4,'(i4.4)') it 
        filename=trim(folder)//'possoft_at'//fn4//'.ascii'
        call write_atomic_file_ascii(parini,filename,parini%nat,units,pos_red_in,latvec_in,fxyz,strten_in,parini%char_type,&
             &parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,etot_in,pressure,curv,res)
endif
        call elim_fixed_at(parini,parini%nat,fxyz)
        do i=1,3*parini%nat
        wpos(i)=wpos(i)+alpha*fxyz(i)
        enddo

        do i=1,3*parini%nat
        ddcart(i)=wpos(i)-rxyz(i)
        enddo
         
        call elim_moment(parini%nat,ddcart(1:3*parini%nat),amass)
        call elim_fixed_at(parini,parini%nat,ddcart)
                         
        sdd=0.d0
        do i=1,3*parini%nat
        sdd=sdd+ddcart(i)*ddcart(i)
        enddo

!        if (res.le.curv*eps_dd*5.d-1) goto 1000
        sdd=eps_dd/sqrt(sdd)
        do i=1,3*parini%nat
        ddcart(i)=ddcart(i)*sdd
        enddo
      enddo
      !call yaml_sequence_close()

!       write(*,*) '# No convergence in low_cur_dir',res
1000   continue
      
!write(*,'(a,i5,4(e13.5),e18.10)')' # SOFTEN: final atomic it,fnrm,res,curv,fd2,etot ',it,tt,res,curv,fd2,etot-etot0
if(parini%verb.gt.0) call yaml_sequence_close()
call yaml_mapping_open('SOFTEN final atomic',flow=.true.)
call yaml_map('it',it,fmt='(i5)')
call yaml_map('fnrm',tt,fmt='(e13.5)')
call yaml_map('res',res,fmt='(e13.5)')
call yaml_map('curv',curv,fmt='(e13.5)')
call yaml_map('fd2',fd2,fmt='(e13.5)')
call yaml_map('de',etot-etot0,fmt='(e18.10)')
call yaml_mapping_close()
!Decrease the stepsize if necessary
       if(parini%auto_soft.and.nsoft.lt.3) write(*,'(a,e18.10)') ' # SOFTEN: increase nsoft for auto adjustment'
       if(parini%auto_soft.and.nsoft.ge.3) then
          if(decrease) then 
             parres%alpha_at=parres%alpha_at/1.1d0
          else
             parres%alpha_at=parres%alpha_at*1.1d0
          endif
          call yaml_map('new alpha_at',parres%alpha_at,fmt='(es18.10)')
       !write(*,'(a,e18.10)') ' # SOFTEN: new alpha_at :  ',parres%alpha_at
       endif


       deallocate(wpos,fxyz)
    end subroutine

!************************************************************************************

 subroutine soften_lat(parini,parres,latvec,pos_red0,ddlat,curv0,curv,res,pressure,count_soft,amass,nsoft,folder)
 use global, only: units
 use defs_basis
 use interface_code
 use modsocket, only: sock_extra_string
 use mod_parini, only: typ_parini
 use yaml_output
 implicit none
 type(typ_parini), intent(in):: parini
 type(typ_parini), intent(inout):: parres
 real(8) :: acell_in(3),xred_in(3,parini%nat),vel_in(3,parini%nat),fcart_in(3,parini%nat),etot_in,strten_in(6),vel_lat_in(3,3)
!*******************************************************************
        integer:: nsoft,i,it,nit,iprec,iat
        real(8):: curv0,curv,res,pressure,count_soft,alpha,alphalat
        real(8):: rxyz(3*parini%nat)
        real(8):: latvec(9),latvec_in(9)
        real(8):: dd(3*parini%nat)
        real(8):: ddcart(3*parini%nat)
        real(8):: rxyzcart(3*parini%nat)
        real(8):: ddlat(9)
        real(8):: ddall(3*parini%nat+9)
        real(8):: flat(9)
        real(8):: pos_red0(3*parini%nat)
        real(8):: amass(parini%nat)
        real(8):: wlat(9),wlatold(9),fxyzcart(3*parini%nat)
        real(8), allocatable :: wpos(:),fxyz(:)
        real(8):: eps_dd
        real(8):: etot
        real(8):: etot0
        real(8):: fd2
        real(8):: res1 
        real(8):: res2
        real(8):: sdd
        real(8):: sdf
        real(8):: tt
        real(8):: vol
        character(40):: filename,folder
        character(4):: fn4         
        logical:: getwfk
        real(8):: lat_prev(3*3),dir_prev(3*3),dir(3*3),angle,vnrm
        logical:: decrease
        decrease=.false.

        call yaml_map('Entering SOFTENING routine for LATTICE, nsoft',nsoft)
        !write(*,'(a,i5)')" # Entering SOFTENING routine for LATTICE, nsoft= ",nsoft
        !if(parini%auto_soft) write(*,'(a)')" # Automatic softening activated"
        if(parini%auto_soft) call yaml_comment('Automatic softening activated',hfill='~')

!        ddcart=dd
        rxyz=pos_red0
!        rxyzcart=rxyz        

        nit=nsoft
!        eps_dd=1.d0
        eps_dd=1.d0
!        alphalat=1.d0 !step size for Lattice
        alphalat=parres%alpha_lat

latvec_in=latvec
iprec=1;getwfk=.false.
       write(fn4,'(i4.4)') 0
       sock_extra_string="SOFTLAT"//trim(fn4)
call get_energyandforces_single(parini,parres,latvec_in,rxyz,fcart_in,strten_in,etot_in,iprec,getwfk)
call strten2flat(strten_in,flat,latvec,pressure)
call get_enthalpy(latvec_in,etot_in,pressure,etot0)
!Multiply the force with the unit cell volume to get the correct force
!call getvol(latvec,vol)
!flat=flat*vol

!         call  fxyz_cart2int(nat,fxyzcart,fxyz,latvec)

! normalize initial guess
9909 continue
        sdd=0.d0
        do i=1,9
        sdd=sdd+ddlat(i)**2
        enddo
        sdd=eps_dd/sqrt(sdd)
        do i=1,9
        ddlat(i)=sdd*ddlat(i)
        enddo
!Check if this eps_dd will invert the cell
        do i=1,9
        wlat(i)=latvec(i)+ddlat(i)
        enddo
        call getvol(wlat,vol)
        if(vol.le.0.d0) then
              eps_dd=eps_dd*0.5d0
              write(*,*) "Epsdd set to ",eps_dd,vol
              call getvol(latvec,vol)
              
              goto 9909
        endif


if(parini%verb.gt.0) call yaml_sequence_open('SOFTEN lattice iterations')
 do it=1,nit
    if(parini%verb.gt.0) call yaml_sequence(advance='no')
        do i=1,9
        wlat(i)=latvec(i)+ddlat(i)
        enddo

!Here we check for direction variation during softening
        dir_prev=dir
        dir=wlat-lat_prev
        if(it.ge.3) then
           do iat=1,3  
              angle=dot_product(dir((iat-1)*3+1:(iat)*3),dir_prev((iat-1)*3+1:(iat)*3))
              vnrm=sqrt(dot_product(dir((iat-1)*3+1:(iat)*3),dir((iat-1)*3+1:(iat)*3)))+&
                   sqrt(dot_product(dir_prev((iat-1)*3+1:(iat)*3),dir_prev((iat-1)*3+1:(iat)*3)))
              angle=angle/vnrm
              if(angle.lt.0.d0) decrease=.true.
           enddo
        endif

        lat_prev=wlat

latvec_in=wlat
iprec=1
if(parini%usewf_soften) then
    getwfk=.true.
else
    getwfk=.false.
endif
       write(fn4,'(i4.4)') it
       sock_extra_string="SOFTLAT"//trim(fn4)
call get_energyandforces_single(parini,parres,latvec_in,rxyz,fcart_in,strten_in,etot_in,iprec,getwfk)
call strten2flat(strten_in,flat,wlat,pressure)
call get_enthalpy(latvec_in,etot_in,pressure,etot)
!Multiply the force with the unit cell volume to get the correct force
!call getvol(wlat,vol)
!flat=flat*vol

        fd2=2.d0*(etot-etot0)/eps_dd**2

        sdf=0.d0
        sdd=0.d0
        do i=1,9
        sdf=sdf+ddlat(i)*flat(i)
        sdd=sdd+ddlat(i)*ddlat(i)
        enddo

        curv=-sdf/sdd
        if (it.eq.1) curv0=curv
        res=0.d0
        res1=0.d0
        res2=0.d0
        tt=0.d0

        do i=1,9
        tt=tt+flat(i)**2
        flat(i)=flat(i)+curv*ddlat(i)
        res=res+flat(i)**2
        enddo

if(parini%verb.gt.0) then
    call yaml_mapping_open('SOFTEN',flow=.true.)
    call yaml_map('res',res,fmt='(e13.5)')
    call yaml_map('it',it,fmt='(i8)')
    call yaml_mapping_close()
    !write(*,'(a,(e13.5),i5)') ' # SOFTEN: ',res,it
endif
        res=sqrt(res)
        tt=sqrt(tt)
!if(it==1) write(*,'(a,i5,4(e13.5),e18.10)')' # SOFTEN: init lattice  it,fnrm,res,curv,fd2,etot ',it,tt,res,curv,fd2,etot-etot0
if(it==1) then
    call yaml_mapping_open('SOFTEN init lattice',flow=.true.)
    call yaml_map('it',it,fmt='(i5)')
    call yaml_map('fnrm',tt,fmt='(e13.5)')
    call yaml_map('res',res,fmt='(e13.5)')
    call yaml_map('curv',curv,fmt='(e13.5)')
    call yaml_map('fd2',fd2,fmt='(e13.5)')
    call yaml_map('de',etot-etot0,fmt='(e18.10)')
    call yaml_mapping_close()
endif
if(parini%verb.gt.0) then
        !write(*,'(a,i5,4(e13.5),e18.10)') ' # SOFTEN: it,fnrm,res,curv,fd2,etot ',it,tt,res,curv,fd2,etot-etot0
        call yaml_mapping_open('SOFTEN',flow=.true.)
        call yaml_map('it',it,fmt='(i5)')
        call yaml_map('fnrm',tt,fmt='(e13.5)')
        call yaml_map('res',res,fmt='(e13.5)')
        call yaml_map('curv',curv,fmt='(e13.5)')
        call yaml_map('fd2',fd2,fmt='(e13.5)')
        call yaml_map('de',etot-etot0,fmt='(e18.10)')
        call yaml_mapping_close()
        write(fn4,'(i4.4)') it 
        filename=trim(folder)//'possoft_lat'//fn4//'.ascii'
        call write_atomic_file_ascii(parini,filename,parini%nat,units,rxyz,latvec_in,fcart_in,strten_in,parini%char_type,parini%ntypat_global,&
        &parini%typat_global,parini%fixat,parini%fixlat,etot_in,pressure,curv,res)
endif

        wlatold=wlat
        do i=1,9
        wlat(i)=wlat(i)+alphalat*flat(i)
!        write(*,*) " # FLAT",flat(i)
        enddo

        do i=1,9
        ddlat(i)=wlat(i)-latvec(i)
        enddo
         
        call elim_torque_cell(latvec,ddlat)
                         
        sdd=0.d0
        do i=1,9
        sdd=sdd+ddlat(i)*ddlat(i)
        enddo

!        if (res.le.curv*eps_dd*1.d-2) goto 1000
        sdd=eps_dd/sqrt(sdd)
        do i=1,9
        ddlat(i)=ddlat(i)*sdd
        enddo
      enddo

!       write(*,*) '# No convergence in low_cur_dir',res
1000   continue

!write(*,'(a,i5,4(e13.5),e18.10)')' # SOFTEN: final lattice it,fnrm,res,curv,fd2,etot ',it,tt,res,curv,fd2,etot-etot0
if(parini%verb.gt.0) call yaml_sequence_close()
call yaml_mapping_open('SOFTEN final lattice',flow=.true.)
call yaml_map('it',it,fmt='(i5)')
call yaml_map('fnrm',tt,fmt='(e13.5)')
call yaml_map('res',res,fmt='(e13.5)')
call yaml_map('curv',curv,fmt='(e13.5)')
call yaml_map('fd2',fd2,fmt='(e13.5)')
call yaml_map('de',etot-etot0,fmt='(e18.10)')
call yaml_mapping_close()
!Decrease the stepsize if necessary
       if(parini%auto_soft.and.nsoft.lt.3) write(*,'(a,e18.10)') ' # SOFTEN: increase nsoft for auto adjustment'
       if(parini%auto_soft.and.nsoft.ge.3) then
          if(decrease) then
             parres%alpha_lat=parres%alpha_lat/1.1d0
          else
             parres%alpha_lat=parres%alpha_lat*1.1d0
          endif
          call yaml_map('new alpha_at',parres%alpha_at,fmt='(es18.10)')
          call yaml_map('new alpha_lat',parres%alpha_lat,fmt='(es18.10)')
       !write(*,'(a,es18.10,1x,es18.10)') ' # SOFTEN: new alpha_at, alpha_lat : ',parres%alpha_at,parres%alpha_lat
       endif


    end subroutine
