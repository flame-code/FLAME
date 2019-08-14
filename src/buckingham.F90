!*****************************************************************************************
subroutine set_buckingham(atoms,tosifumi)
    use mod_atoms, only: typ_atoms
    use mod_shortrange, only: typ_tosifumi
    implicit none
    type(typ_atoms), intent(inout):: atoms
    type(typ_tosifumi), intent(inout):: tosifumi
    !local variables
    integer:: max_ntypinter, stat 
    integer:: ntypinter, it, ityp1,ityp2
    integer:: itypinter, itypat, jtypat, nall, iall
    integer:: leni, lenj
    character(6):: strint(10),tmp !string of interaction
    real(8):: parameters(5,10)
    real(8), allocatable:: parameters2(:,:),qq(:)
    character(20):: strtmpi, strtmpj
    character(4):: namatnamat1, namatnamat2
    character(2), allocatable:: type1(:), type2(:) 
    integer, allocatable:: ntypint(:,:)
    real(8):: dx, dy, dz, r,rsq,xiat,yiat,ziat,alphainv
    integer:: iat, jat

    max_ntypinter=(atoms%ntypat**2+atoms%ntypat)
    allocate (parameters2(0:max_ntypinter,4),type1(max_ntypinter),type2(max_ntypinter))
    allocate (tosifumi%interaction(atoms%ntypat,atoms%ntypat),qq(0:max_ntypinter))
    parameters2(0:max_ntypinter,1)=0.0
    parameters2(0:max_ntypinter,2)=1.0
    parameters2(0:max_ntypinter,3)=0.0
    parameters2(0:max_ntypinter,4)=0.0
    qq = 0.d0
    type1(:)='no'
    type2(:)='no'
    open (unit= 121 , file ="param_short.in")
    read(121, *) 
    ntypinter=0
    do
        ntypinter=ntypinter+1
        read(121, *, IOSTAT=stat) type1(ntypinter),type2(ntypinter),parameters2(ntypinter,:)
        if(IS_IOSTAT_END(stat)) exit 
    end do
    ntypinter=ntypinter-1
    tosifumi%ntypinter=ntypinter
    allocate(tosifumi%aaa(0:max_ntypinter/2),  tosifumi%bbb(0:max_ntypinter/2))
    allocate(tosifumi%ccc(0:max_ntypinter/2),  tosifumi%ddd(0:max_ntypinter/2))
    allocate(tosifumi%eee(0:max_ntypinter/2),  tosifumi%fff(0:max_ntypinter/2))

    tosifumi%interaction(:,:)= 0
    do it=1, ntypinter
        do itypat=1,atoms%ntypat
            if (type1(it)==atoms%stypat(itypat)) ityp1=itypat 
            if (type2(it)==atoms%stypat(itypat)) ityp2=itypat
        enddo
        tosifumi%interaction(ityp1,ityp2)=it
    enddo
    do itypat=1,atoms%ntypat
    do jtypat=1,atoms%ntypat
        if(tosifumi%interaction(itypat,jtypat)== 0) then
            tosifumi%interaction(itypat,jtypat)=tosifumi%interaction(jtypat,itypat)
        endif
    enddo
    enddo
    write(*,*) "Buckingham parameters : "
    write(*,*) "---------------------------------------------------------- "
    write(*,'(2a5,a7,4a12)') "","","inter","A","B","C","D"

    do itypat=1,atoms%ntypat
    do jtypat=1,atoms%ntypat
        if (tosifumi%interaction(itypat,jtypat)/=tosifumi%interaction(jtypat,itypat)) then
            write(*,*) 'ERROR: 2 different parameters for interaction ' ,atoms%stypat(itypat) ,'and ' ,atoms%stypat(jtypat) ,' are defined'
            stop
        endif
    enddo
    enddo

    do itypat=1,atoms%ntypat
    do jtypat=itypat,atoms%ntypat
        if (tosifumi%interaction(itypat,jtypat)==0) then
            ntypinter=ntypinter+1
            tosifumi%interaction(itypat,jtypat)=ntypinter
            tosifumi%interaction(jtypat,itypat)=ntypinter
        endif
    enddo
    enddo
    write(*,*) "---------------------------------------------------------- "
    do itypat=1,atoms%ntypat
    do jtypat=1,atoms%ntypat
        write(*,'(2a5,i7,6f12.4)') atoms%stypat(itypat),atoms%stypat(jtypat) &
        & ,tosifumi%interaction(itypat,jtypat),parameters2(tosifumi%interaction(itypat,jtypat),:)&
        &, atoms%qtypat(itypat),atoms%qtypat(jtypat) 

        qq(tosifumi%interaction(itypat,jtypat))=atoms%qtypat(itypat)*atoms%qtypat(jtypat)    

    enddo
    enddo

    write(*,*) "---------------------------------------------------------- "
    tosifumi%aaa= 0.d0
    tosifumi%bbb= 0.d0
    tosifumi%ccc= 0.d0
    tosifumi%ddd= 0.d0
    tosifumi%eee= 0.d0
    
    do it=0, ntypinter
        tosifumi%aaa(it) =parameters2(it,1)
        tosifumi%bbb(it) =parameters2(it,2)
        tosifumi%ccc(it) =parameters2(it,3)
        tosifumi%ddd(it) =parameters2(it,4)
        tosifumi%eee(it) = qq(it)
    enddo
    where (tosifumi%interaction==0)
        tosifumi%interaction=max_ntypinter/2
    endwhere
!    do iat=1,atoms%nat
!        xiat=atoms%rat(1,iat)
!        yiat=atoms%rat(2,iat)
!        ziat=atoms%rat(3,iat)
!        do jat=iat+1,atoms%nat
!            dx=xiat-atoms%rat(1,jat)
!            dy=yiat-atoms%rat(2,jat)
!            dz=ziat-atoms%rat(3,jat)
!            rsq=dx*dx+dy*dy+dz*dz
!            r=sqrt(rsq)
!            write(*,*)iat,jat,r
!        enddo
!    enddo
end subroutine set_buckingham
!*****************************************************************************************
