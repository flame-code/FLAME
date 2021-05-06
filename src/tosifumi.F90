!*****************************************************************************************
subroutine set_tosifumi(atoms,tosifumi)
    use mod_atoms, only: typ_atoms
    use mod_shortrange, only: typ_tosifumi
    implicit none
    type(typ_atoms), intent(inout):: atoms
    type(typ_tosifumi), intent(inout):: tosifumi
    !local variables
    integer:: ntypinter
    integer:: itypinter, itypat, jtypat, nall, iall
    integer:: leni, lenj
    character(6):: strint(10) !string of interaction
    real(8):: parameters(5,10)
    character(20):: strtmpi, strtmpj
    character(4):: namatnamat1, namatnamat2
    tosifumi%ntypinter=(atoms%ntypat**2+atoms%ntypat)/2
    nall=10
    write(*,*) 'nall=',nall
    allocate (tosifumi%interaction(nall,nall))
    allocate (tosifumi%aaa(nall))
    allocate (tosifumi%bbb(nall))
    allocate (tosifumi%ccc(nall))
    allocate (tosifumi%ddd(nall))
    allocate (tosifumi%eee(nall))
    allocate (tosifumi%fff(nall))
    call tosifumi_parameters(strint,parameters)
    !write(*,*) 'atoms%ntypat=',atoms%ntypat
    !stop
    !write(*,*) atoms%itypat(1)
    !write(*,*) atoms%itypat(2)
    !write(*,*) atoms%itypat(3)
    !write(*,*) atoms%itypat(4)
    !stop
    itypinter=0
    do itypat=1,atoms%ntypat
        strtmpi=atoms%stypat(itypat)
        leni=len_trim(strtmpi)
        do jtypat=itypat,atoms%ntypat
            itypinter=itypinter+1
            tosifumi%interaction(itypat,jtypat)=itypinter
            tosifumi%interaction(jtypat,itypat)=itypinter
            strtmpj=atoms%stypat(jtypat)
            lenj=len_trim(strtmpj)
            namatnamat1=strtmpi(1:leni)//strtmpj(1:lenj)
            namatnamat2=strtmpj(1:lenj)//strtmpi(1:leni)
            do iall=1,nall
                if(namatnamat1==strint(iall) .or. namatnamat2==strint(iall)) then
                tosifumi%aaa(itypinter)=parameters(1,iall)
                tosifumi%bbb(itypinter)=parameters(2,iall)
                tosifumi%ccc(itypinter)=parameters(3,iall)
                tosifumi%ddd(itypinter)=parameters(4,iall)
                tosifumi%eee(itypinter)=parameters(5,iall)
                endif
            enddo
            !write(*,'(i3,f10.4,2(2x,a))') itypinter,tosifumi%ccc(itypinter),namatnamat1,namatnamat2
        enddo
    enddo
end subroutine set_tosifumi
!*****************************************************************************************
subroutine coulomb_free_direct(atoms)
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_shortrange, only: typ_tosifumi
    implicit none
    type(typ_atoms), intent(inout):: atoms
    !local variables
    real(8):: rinv, rinv2, rinvqq
    real(8):: dx, dy, dz, r,rsq, xiat, yiat, ziat
    real(8):: ttt, fx, fy, fz, qq
    integer:: iat, jat
    atoms%epot=0.d0
    do iat=1,atoms%nat
        atoms%fat(1,iat)=0.d0
        atoms%fat(2,iat)=0.d0
        atoms%fat(3,iat)=0.d0
    enddo
    call update_ratp(atoms)
    do iat=1,atoms%nat
        xiat=atoms%ratp(1,iat)
        yiat=atoms%ratp(2,iat)
        ziat=atoms%ratp(3,iat)
        do jat=iat+1,atoms%nat
            dx=xiat-atoms%ratp(1,jat)
            dy=yiat-atoms%ratp(2,jat)
            dz=ziat-atoms%ratp(3,jat)
            rsq=dx*dx+dy*dy+dz*dz
            r=sqrt(rsq)
            rinv=1.d0/r
            qq=atoms%qat(iat)*atoms%qat(jat)
            !-------------------------------------------------------------------
            !write(61,'(5i3,5f12.6)') iat,jat,atoms%itypat(iat),atoms%itypat(jat),itypinter,r,a2,a3,a4,a5
            rinvqq=rinv*qq
            rinv2=rinv**2
            atoms%epot=atoms%epot+rinvqq
            ttt=rinv2*rinvqq
            fx=ttt*dx
            fy=ttt*dy
            fz=ttt*dz
            !-------------------------------------------------------------------
            atoms%fat(1,iat)=atoms%fat(1,iat)+fx
            atoms%fat(2,iat)=atoms%fat(2,iat)+fy
            atoms%fat(3,iat)=atoms%fat(3,iat)+fz
            atoms%fat(1,jat)=atoms%fat(1,jat)-fx
            atoms%fat(2,jat)=atoms%fat(2,jat)-fy
            atoms%fat(3,jat)=atoms%fat(3,jat)-fz
        enddo
    enddo
end subroutine coulomb_free_direct
!*****************************************************************************************
subroutine calenergyforces(atoms,tosifumi)
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_shortrange, only: typ_tosifumi
    implicit none
    type(typ_atoms), intent(inout):: atoms
    type(typ_tosifumi), intent(in):: tosifumi
    !local variables
    integer::itypinter
    real(8):: a2, a3, a4, a5, rinv, rinv2, rinv6, rinv8, rinv10 !, rinvqq
    real(8):: dx, dy, dz, r,rsq,xiat,yiat,ziat,alphainv
    real(8):: t2, ttt, fx, fy, fz
    integer:: iat, jat
    !write(*,*) atoms%itypat(1)
    !write(*,*) atoms%itypat(2)
    !write(*,*) atoms%itypat(3)
    !write(*,*) atoms%itypat(4)
    !stop
    call update_ratp(atoms)
    do iat=1,atoms%nat
        xiat=atoms%ratp(1,iat)
        yiat=atoms%ratp(2,iat)
        ziat=atoms%ratp(3,iat)
        do jat=iat+1,atoms%nat
            dx=xiat-atoms%ratp(1,jat)
            dy=yiat-atoms%ratp(2,jat)
            dz=ziat-atoms%ratp(3,jat)
            rsq=dx*dx+dy*dy+dz*dz
            r=sqrt(rsq)
            rinv=1.d0/r
            itypinter=tosifumi%interaction(atoms%itypat(iat),atoms%itypat(jat))
            !write(*,*) iat,jat,itypinter
            !write(*,*) iat,jat,r
            !-------------------------------------------------------------------
            a2=tosifumi%ccc(itypinter)
            a3=tosifumi%ddd(itypinter)
            a4=tosifumi%aaa(itypinter)
            a5=tosifumi%bbb(itypinter)
            !write(61,'(5i3,5f12.6)') iat,jat,atoms%itypat(iat),atoms%itypat(jat),itypinter,r,a2,a3,a4,a5
            rinv2=rinv**2
            rinv6=rinv2**3
            rinv8=rinv6*rinv2
            rinv10=rinv8*rinv2 
            !rinvqq=rinv*tosifumi%eee(itypinter)
            t2=a4*exp(-a5*r)
            atoms%epot=atoms%epot-a2*rinv6-a3*rinv8+t2
            ttt=-6.d0*a2*rinv8-8.d0*a3*rinv10+a5*t2*rinv
            fx=ttt*dx
            fy=ttt*dy
            fz=ttt*dz
            !-------------------------------------------------------------------
            atoms%fat(1,iat)=atoms%fat(1,iat)+fx
            atoms%fat(2,iat)=atoms%fat(2,iat)+fy
            atoms%fat(3,iat)=atoms%fat(3,iat)+fz
            atoms%fat(1,jat)=atoms%fat(1,jat)-fx
            atoms%fat(2,jat)=atoms%fat(2,jat)-fy
            atoms%fat(3,jat)=atoms%fat(3,jat)-fz
        enddo
    enddo
end subroutine calenergyforces
!*****************************************************************************************
subroutine tosifumi_parameters(s,p)
    implicit none
    character(6), intent(out):: s(10)
    real(8), intent(out):: p(5,10)
    !local variables
    s( 1)='NaNa' ; p(1:5, 1)=(/ 15.5667657558d0,1.6693280757d0,  1.7548585478d0,   2.9841517285d0, 1.d0/)
    s( 2)='NaCl' ; p(1:5, 2)=(/ 46.1152157465d0,1.6693280757d0, 11.6990569852d0,  51.8496362832d0,-1.d0/)
    s( 3)='ClCl' ; p(1:5, 3)=(/128.0741185691d0,1.6693280757d0,121.1688044893d0, 869.1341909350d0, 1.d0/)
    s( 4)='NaBr' ; p(1:5, 4)=(/ 37.6572882385d0,1.5564029412d0, 14.6238212315d0,  70.8736035526d0,-1.d0/)
    s( 5)='BrBr' ; p(1:5, 5)=(/140.7136247242d0,1.5564029412d0,204.7334972405d0,1678.5853472993d0, 1.d0/)
    s( 6)='KK  ' ; p(1:5, 6)=(/ 57.1641728301d0,1.5702581602d0, 25.3827754232d0,  89.5245518560d0, 1.d0/)
    s( 7)='KCl ' ; p(1:5, 7)=(/ 65.6804690041d0,1.5702581602d0, 50.1388156507d0, 272.3038452286d0,-1.d0/)
    s( 8)='KBr ' ; p(1:5, 8)=(/102.4980177035d0,1.5796328358d0, 62.6735195634d0, 369.2887764059d0,-1.d0/)
    s( 9)='NaK ' ; p(1:5, 9)=(/ 29.8305428725d0,1.6190355256d0,  6.6740677564d0,  16.3448721674d0, 1.d0/)
    s(10)='ClBr' ; p(1:5,10)=(/134.2451990099d0,1.6118768957d0,157.5033748830d0,1207.8559176244d0, 1.d0/)
end subroutine tosifumi_parameters
!*****************************************************************************************
