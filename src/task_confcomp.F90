!*****************************************************************************************
subroutine conf_comp(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_all, atom_all_allocate, atom_all_deallocate, set_rcov
    use mod_acf, only: acf_read
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    !integer:: iat
    !integer, parameter:: natmax=10**4
    type(typ_atoms_all):: atoms_all
    character(256):: comment1, comment2
    integer:: iconf, jconf, ifp
    real(8):: ttnrm, ttmax, tt, ttsum
    atoms_all%atoms%tol=parini%tol_conf_comp
    call acf_read(parini,'posinp.acf',10,atoms_all=atoms_all)
    call set_rcov(atoms_all%atoms)
    !print '(a,f7.2)','RCOV',atoms%rcov(1)
    atoms_all%atoms%nfp=10*atoms_all%atoms%nat
    !call atom_copy(atoms,atoms_all%atoms,'conf_comp: atoms->atoms_all%atoms')
    atoms_all%nconfmax=100
    call atom_all_allocate(atoms_all,fpall=.true.)
    !atoms_all%ratall(1:3,1:atoms%nat,1:atoms_all%nconfmax)=atoms_all%ratall(1:3,1:atoms%nat,1:atoms_all%nconfmax)*0.529177d0
    !call set_fpall_distance(atoms_all)
    !call set_fpall_angle(atoms_all)
    call set_fpall_ann(atoms_all)
    write(*,*) atoms_all%nconf,atoms_all%atoms%nat,atoms_all%atoms%nfp
    do iconf=1,atoms_all%nconf
        !do ifp=1,atoms_all%atoms%nfp
        !    write(100+iconf,'(i6,es12.4)') ifp,atoms_all%fpall(ifp,iconf)
        !enddo
        do jconf=iconf,atoms_all%nconf
            !if(iconf==jconf) cycle
            ttnrm=0.d0
            ttmax=0.d0
            do ifp=1,atoms_all%atoms%nfp
                !write(61,'(i6,es12.4)') ifp,atoms_all%fpall(ifp,1)
                ttnrm=ttnrm+(atoms_all%fpall(ifp,jconf)-atoms_all%fpall(ifp,iconf))**2
                ttmax=max(ttmax,abs(atoms_all%fpall(ifp,jconf)-atoms_all%fpall(ifp,iconf)))
                !tt=max(1.d0,abs(atoms_all%fpall(ifp,jconf)+atoms_all%fpall(ifp,iconf))*0.5d0)
                !ttmax=max(ttmax,abs(atoms_all%fpall(ifp,jconf)-atoms_all%fpall(ifp,iconf))/tt)
            enddo
            ttnrm=sqrt(ttnrm/atoms_all%atoms%nat)
            ttsum=abs(atoms_all%fpall(1,jconf)-atoms_all%fpall(1,iconf))
            !write(*,'(es7.0)',advance='no') ttnrm
            if(iconf==jconf) cycle
            write(21,'(4es12.4,2x,i4.4,a1,i4.4)') &
                abs(atoms_all%epotall(iconf)-atoms_all%epotall(jconf)),ttnrm,ttmax,ttsum,iconf,'_',jconf
        enddo
        !write(*,*)
    enddo
    call atom_all_deallocate(atoms_all,ratall=.true.,epotall=.true.,fpall=.true.)
end subroutine conf_comp
!*****************************************************************************************
subroutine set_fpall_ann(atoms_all)
    use mod_atoms, only: typ_atoms_all, set_rat
    use dynamic_memory
    implicit none
    type(typ_atoms_all), intent(inout):: atoms_all
    !local variables
    integer:: ifp, iat, jat, iconf, natp, natpmax
    real(8):: dx, dy, dz, r, rcut, rcov
    real(8):: tt1, tt2, tt3, tt4, shift, tt, ttfp, ttfp1, ttfp2, ttfp3, ttfp4, ttfp5
    real(8), allocatable:: ratp(:,:)
    !real(8), allocatable:: fpt(:,:)
    call f_routine(id='set_fpall_ann')
    natpmax=27*atoms_all%atoms%nat
    ratp=f_malloc0([1.to.3,1.to.natpmax],id='ratp')
    !allocate(fpt(100000),source=0.d0)
    !rcov=1.1d0
    !rcut=3.d0*rcov
    do iconf=1,atoms_all%nconf
        call set_rat(atoms_all%atoms,atoms_all%ratall(1,1,iconf),setall=.true.)
        call build_images(atoms_all%atoms,natpmax,natp,ratp)
        print '(a,i5)','natp=',natp
        !write(*,*) 'iconf=',iconf
        ifp=0
        do iat=1,atoms_all%atoms%nat
            !---------------------------
            shift=2.d0
            rcut=3.d0 !*rcov
            ttfp=0.d0
            do jat=1,natp !atoms_all%atoms%nat
                if(iat==jat) cycle
                dx=ratp(1,jat)-atoms_all%ratall(1,iat,iconf)
                dy=ratp(2,jat)-atoms_all%ratall(2,iat,iconf)
                dz=ratp(3,jat)-atoms_all%ratall(3,iat,iconf)
                r=sqrt(dx*dx+dy*dy+dz*dz)
                rcov=atoms_all%atoms%rcov(iat)+atoms_all%atoms%rcov(modulo(jat,atoms_all%atoms%nat)+1)
                if(.not. r>rcut) then
                    tt1=r/rcut
                    if(r>shift) then
                        tt1=(r-shift)/(rcut-shift)
                        tt3=(1.d0-tt1**2)**2
                    else
                        tt3=1.d0
                    endif
                    ttfp=ttfp+tt3
                endif
            enddo
            ifp=ifp+1
            atoms_all%fpall(ifp,iconf)=ttfp
            !---------------------------
            shift=4.d0
            rcut=5.d0 !*rcov
            ttfp=0.d0
            do jat=1,natp !atoms_all%atoms%nat
                if(iat==jat) cycle
                dx=ratp(1,jat)-atoms_all%ratall(1,iat,iconf)
                dy=ratp(2,jat)-atoms_all%ratall(2,iat,iconf)
                dz=ratp(3,jat)-atoms_all%ratall(3,iat,iconf)
                r=sqrt(dx*dx+dy*dy+dz*dz)
                rcov=atoms_all%atoms%rcov(iat)+atoms_all%atoms%rcov(modulo(jat,atoms_all%atoms%nat)+1)
                if(.not. r>rcut) then
                    tt1=r/rcut
                    if(r>shift) then
                        tt1=(r-shift)/(rcut-shift)
                        tt3=(1.d0-tt1**2)**2
                    else
                        tt3=1.d0
                    endif
                    ttfp=ttfp+tt3
                endif
            enddo
            ifp=ifp+1
            atoms_all%fpall(ifp,iconf)=ttfp
            !---------------------------
            shift=6.d0
            rcut=7.d0 !*rcov
            ttfp=0.d0
            do jat=1,natp !atoms_all%atoms%nat
                if(iat==jat) cycle
                dx=ratp(1,jat)-atoms_all%ratall(1,iat,iconf)
                dy=ratp(2,jat)-atoms_all%ratall(2,iat,iconf)
                dz=ratp(3,jat)-atoms_all%ratall(3,iat,iconf)
                r=sqrt(dx*dx+dy*dy+dz*dz)
                rcov=atoms_all%atoms%rcov(iat)+atoms_all%atoms%rcov(modulo(jat,atoms_all%atoms%nat)+1)
                if(.not. r>rcut) then
                    tt1=r/rcut
                    if(r>shift) then
                        tt1=(r-shift)/(rcut-shift)
                        tt3=(1.d0-tt1**2)**2
                    else
                        tt3=1.d0
                    endif
                    ttfp=ttfp+tt3
                endif
            enddo
            ifp=ifp+1
            atoms_all%fpall(ifp,iconf)=ttfp
            !---------------------------
            shift=6.d0
            rcut=7.d0 !*rcov
            !ttfp=0.d0
            ttfp1=0.d0 ; ttfp2=0.d0 ; ttfp3=0.d0 ; ttfp4=0.d0 ; ttfp5=0.d0
            do jat=1,natp !atoms_all%atoms%nat
                if(iat==jat) cycle
                dx=ratp(1,jat)-atoms_all%ratall(1,iat,iconf)
                dy=ratp(2,jat)-atoms_all%ratall(2,iat,iconf)
                dz=ratp(3,jat)-atoms_all%ratall(3,iat,iconf)
                r=sqrt(dx*dx+dy*dy+dz*dz)
                rcov=atoms_all%atoms%rcov(iat)+atoms_all%atoms%rcov(modulo(jat,atoms_all%atoms%nat)+1)
                if(.not. r>rcut) then
                    tt1=r/rcut
                    if(r>shift) then
                        tt1=(r-shift)/(rcut-shift)
                        tt3=(1.d0-tt1**2)**2
                    else
                        tt3=1.d0
                    endif
                    tt4=exp(-0.5d0*(r-rcov)**2)
                    ttfp1=ttfp1+tt3*tt4
                    tt4=exp(-0.5d0*(r-1.5d0*rcov)**2)
                    ttfp2=ttfp2+tt3*tt4
                    tt4=exp(-0.5d0*(r-2.d0*rcov)**2)
                    ttfp3=ttfp3+tt3*tt4
                    tt4=exp(-0.5d0*(r-2.5d0*rcov)**2)
                    ttfp4=ttfp4+tt3*tt4
                    tt4=exp(-0.5d0*(r-3.d0*rcov)**2)
                    ttfp5=ttfp5+tt3*tt4
                endif
            enddo
            ifp=ifp+1 ; atoms_all%fpall(ifp,iconf)=ttfp1
            ifp=ifp+1 ; atoms_all%fpall(ifp,iconf)=ttfp2
            ifp=ifp+1 ; atoms_all%fpall(ifp,iconf)=ttfp3
            ifp=ifp+1 ; atoms_all%fpall(ifp,iconf)=ttfp4
            ifp=ifp+1 ; atoms_all%fpall(ifp,iconf)=ttfp5
            !---------------------------
        enddo
        call hpsort(atoms_all%atoms%nfp,atoms_all%fpall(1,iconf))
        !tt=0.d0
        !do ifp=1,atoms_all%atoms%nfp
        !    tt=tt+atoms_all%fpall(ifp,iconf)
        !enddo
        !atoms_all%fpall(1,iconf)=tt
    enddo
    call f_free(ratp) !,fpt)
    call f_release_routine()
end subroutine set_fpall_ann
!*****************************************************************************************
subroutine set_fpall_angle(atoms_all)
    use mod_atoms, only: typ_atoms_all, set_rat
    use dynamic_memory
    implicit none
    type(typ_atoms_all), intent(inout):: atoms_all
    !local variables
    integer:: ifp, iat, jat, kat, iconf, natp, natpmax
    real(8):: dxij, dyij, dzij, rij, dxik, dyik, dzik, rik, rcut, rcov, pi, angle
    real(8):: tt1ij, tt3ij, tt1ik, tt3ik, shift, tt
    real(8), allocatable:: ratp(:,:)
    call f_routine(id='set_fpall_angle')
    pi=4.d0*atan(1.d0)
    natpmax=27*atoms_all%atoms%nat
    ratp=f_malloc0([1.to.3,1.to.natpmax],id='ratp')
    !rcov=1.1d0
    !rcut=3.d0*rcov
    shift=1.5d0
    rcut=2.d0
    do iconf=1,atoms_all%nconf
        call set_rat(atoms_all%atoms,atoms_all%ratall(1,1,iconf),setall=.true.)
        call build_images(atoms_all%atoms,natpmax,natp,ratp)
        print '(a,i5)','natp=',natp
        !write(*,*) 'iconf=',iconf
        ifp=0
        do iat=1,atoms_all%atoms%nat
            do jat=1,natp !atoms_all%atoms%nat
                if(iat==jat) cycle
                dxij=ratp(1,jat)-atoms_all%ratall(1,iat,iconf)
                dyij=ratp(2,jat)-atoms_all%ratall(2,iat,iconf)
                dzij=ratp(3,jat)-atoms_all%ratall(3,iat,iconf)
                rij=sqrt(dxij*dxij+dyij*dyij+dzij*dzij)
                if(.not. rij<rcut) cycle
                do kat=1,natp !atoms_all%atoms%nat
                    if(iat==kat .or. jat==kat) cycle
                    !dx=atoms_all%ratall(1,jat,iconf)-atoms_all%ratall(1,iat,iconf)
                    !dy=atoms_all%ratall(2,jat,iconf)-atoms_all%ratall(2,iat,iconf)
                    !dz=atoms_all%ratall(3,jat,iconf)-atoms_all%ratall(3,iat,iconf)
                    dxik=ratp(1,kat)-atoms_all%ratall(1,iat,iconf)
                    dyik=ratp(2,kat)-atoms_all%ratall(2,iat,iconf)
                    dzik=ratp(3,kat)-atoms_all%ratall(3,iat,iconf)
                    rik=sqrt(dxik*dxik+dyik*dyik+dzik*dzik)
                    if(.not. rik<rcut) cycle
                    rcov=atoms_all%atoms%rcov(iat)+atoms_all%atoms%rcov(modulo(jat,atoms_all%atoms%nat)+1)
                    !rcut=4.d0*rcov
                    ifp=ifp+1
                    tt1ik=rik/rcut
                    !tt2=r/rcov
                    if(rij>shift) then
                        tt1ij=(rij-shift)/(rcut-shift)
                        tt3ij=(1.d0-tt1ij**2)**2
                    else
                        tt3ij=1.d0
                    endif
                    if(rik>shift) then
                        tt1ik=(rik-shift)/(rcut-shift)
                        tt3ik=(1.d0-tt1ik**2)**2
                    else
                        tt3ik=1.d0
                    endif
                    !tt4=tt2**1
                    angle=acos((dxij*dxik+dyij*dyik+dzij*dzik)/(rij*rik))*180.d0/pi
                    !write(31,'(3i4,f20.10)') iat,jat,kat,angle

                    !tt4=tt2*exp(-tt2)
                    atoms_all%fpall(ifp,iconf)=angle*tt3ij*tt3ik
                    !write(*,'(6es12.4)') r,tt1,tt2,tt3,tt4,atoms_all%fpall(ifp,iconf)
                    !write(31,'(i6,es12.4)') ifp,atoms_all%fpall(ifp,iconf)
                enddo
            enddo
        enddo !end of loop over iat
        !do ifp=1,atoms_all%atoms%nfp
        !write(41,'(i6,es12.4)') ifp,atoms_all%fpall(ifp,1)
        !enddo
        call hpsort(atoms_all%atoms%nfp,atoms_all%fpall(1,iconf))
        tt=0.d0
        do ifp=1,atoms_all%atoms%nfp
            tt=tt+atoms_all%fpall(ifp,iconf)
        enddo
        atoms_all%fpall(1,iconf)=tt
    enddo
    call f_free(ratp)
    call f_release_routine()
end subroutine set_fpall_angle
!*****************************************************************************************
subroutine set_fpall_distance(atoms_all)
    use mod_atoms, only: typ_atoms_all, set_rat
    use dynamic_memory
    implicit none
    type(typ_atoms_all), intent(inout):: atoms_all
    !local variables
    integer:: ifp, iat, jat, iconf, natp, natpmax
    real(8):: dx, dy, dz, r, rcut, rcov
    real(8):: tt1, tt2, tt3, tt4, shift, tt
    real(8), allocatable:: ratp(:,:)
    call f_routine(id='set_fpall_angle')
    natpmax=27*atoms_all%atoms%nat
    ratp=f_malloc0([1.to.3,1.to.natpmax],id='ratp')
    !rcov=1.1d0
    !rcut=3.d0*rcov
    shift=2.d0
    do iconf=1,atoms_all%nconf
        call set_rat(atoms_all%atoms,atoms_all%ratall(1,1,iconf),setall=.true.)
        call build_images(atoms_all%atoms,natpmax,natp,ratp)
        print '(a,i5)','natp=',natp
        !write(*,*) 'iconf=',iconf
        ifp=0
        do iat=1,atoms_all%atoms%nat
            do jat=1,natp !atoms_all%atoms%nat
                if(iat==jat) cycle
                !dx=atoms_all%ratall(1,jat,iconf)-atoms_all%ratall(1,iat,iconf)
                !dy=atoms_all%ratall(2,jat,iconf)-atoms_all%ratall(2,iat,iconf)
                !dz=atoms_all%ratall(3,jat,iconf)-atoms_all%ratall(3,iat,iconf)
                dx=ratp(1,jat)-atoms_all%ratall(1,iat,iconf)
                dy=ratp(2,jat)-atoms_all%ratall(2,iat,iconf)
                dz=ratp(3,jat)-atoms_all%ratall(3,iat,iconf)
                !write(*,*) atoms_all%ratall(1,jat,iconf),atoms_all%ratall(1,iat,iconf)
                !write(*,*) atoms_all%ratall(2,jat,iconf),atoms_all%ratall(2,iat,iconf)
                !write(*,*) atoms_all%ratall(3,jat,iconf),atoms_all%ratall(3,iat,iconf)
                !stop
                r=sqrt(dx*dx+dy*dy+dz*dz)
                rcov=atoms_all%atoms%rcov(iat)+atoms_all%atoms%rcov(modulo(jat,atoms_all%atoms%nat)+1)
                rcut=7.d0 !*rcov
                if(.not. r>rcut) then
                    ifp=ifp+1
                    tt1=r/rcut
                    tt2=r !/rcov
                    !tt3=(1.d0-tt1**2)**2
                    if(r>shift) then
                        tt1=(r-shift)/(rcut-shift)
                        tt3=(1.d0-tt1**2)**2
                    else
                        tt3=1.d0
                    endif
                    !tt4=1.d0/tt2**0.1d0
                    !tt4=1.d0/tt2**0.5d0
                    tt4=1.d0/tt2
                    !tt4=1.d0/tt2**2
                    !tt4=1.d0/tt2**3
                    !tt4=1.d0/tt2**4
                    !tt4=1.d0/tt2**6
                    !tt4=1.d0
                    !tt4=tt2**1
                    !tt4=tt2**0.5d0
                    !tt4=tt2**1
                    !tt4=tt2*exp(-tt2)
                    atoms_all%fpall(ifp,iconf)=tt3*tt4 !*rcov**10
                    !write(*,'(6es12.4)') r,tt1,tt2,tt3,tt4,atoms_all%fpall(ifp,iconf)
                    !write(31,'(i6,es12.4)') ifp,atoms_all%fpall(ifp,iconf)
                    !ifp=ifp+1
                    !tt4=tt2**2
                    !atoms_all%fpall(ifp,iconf)=tt3*tt4
                endif
            enddo
        enddo
        !do ifp=1,atoms_all%atoms%nfp
        !write(41,'(i6,es12.4)') ifp,atoms_all%fpall(ifp,1)
        !enddo
        call hpsort(atoms_all%atoms%nfp,atoms_all%fpall(1,iconf))
        tt=0.d0
        do ifp=1,atoms_all%atoms%nfp
            tt=tt+atoms_all%fpall(ifp,iconf)
        enddo
        atoms_all%fpall(1,iconf)=tt
    enddo
    call f_free(ratp)
    call f_release_routine()
end subroutine set_fpall_distance
!*****************************************************************************************
subroutine build_images(atoms,natpmax,natp,ratp)
    use mod_atoms, only: typ_atoms, get_rat
    implicit none
    type(typ_atoms), intent(in):: atoms
    integer, intent(in):: natpmax
    integer, intent(inout):: natp
    real(8), intent(inout):: ratp(3,natpmax)
    !local variables
    integer:: iat, ix, iy, iz
    real(8):: x, y, z
    ratp(1:3,1:natpmax)=0.d0
    call get_rat(atoms,ratp)
    natp=atoms%nat
    if(trim(atoms%boundcond)=='free') then
        return
    elseif(trim(atoms%boundcond)=='wire') then
        stop 'ERROR: stopped in build_images'
    elseif(trim(atoms%boundcond)=='slab') then
        stop 'ERROR: stopped in build_images'
    endif
    do iz=-1,1
    do iy=-1,1
    do ix=-1,1
        if(ix==0 .and. iy==0 .and. iz==0) cycle
        do iat=1,atoms%nat
            natp=natp+1
            x=ix*atoms%cellvec(1,1)+iy*atoms%cellvec(1,2)+iz*atoms%cellvec(1,3)
            y=ix*atoms%cellvec(2,1)+iy*atoms%cellvec(2,2)+iz*atoms%cellvec(2,3)
            z=ix*atoms%cellvec(3,1)+iy*atoms%cellvec(3,2)+iz*atoms%cellvec(3,3)
            ratp(1,natp)=ratp(1,iat)+x
            ratp(2,natp)=ratp(2,iat)+y
            ratp(3,natp)=ratp(3,iat)+z
        enddo
    enddo
    enddo
    enddo
end subroutine build_images
!*****************************************************************************************
