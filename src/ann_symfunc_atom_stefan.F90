!******************************************************************************************************
subroutine symmetry_functions_driver_stefan(parini,ann_arr,atoms,symfunc)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc
    use mod_atoms, only: typ_atoms, set_rcov, get_rat
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_symfunc), intent(inout):: symfunc
    !local variables
    integer:: ig, i, iconf
    integer:: iat , jat, kat, istat
    real(8), allocatable:: rat(:,:)
    !-------------------------------------------------------------------------------------
    associate(ng=>ann_arr%ann(1)%nn(0))
    allocate(symfunc%y(ng,atoms%nat),stat=istat,source=0.d0)
    if(istat/=0) stop 'ERROR: unable to allocate array symfunc%y'
    allocate(symfunc%y0d(ng,3,symfunc%linked_lists%maxbound_rad),stat=istat,source=0.d0)
    if(istat/=0) stop 'ERROR: unable to allocate array symfunc%y0d'
    allocate(symfunc%y0dr(ng,9,symfunc%linked_lists%maxbound_rad),stat=istat,source=0.d0)
    if(istat/=0) stop 'ERROR: unable to allocate array symfunc%y0dr'
    symfunc%y0d=0.d0
    symfunc%y0dr=0.d0

    call set_rcov(atoms) 
    allocate(rat(3,atoms%nat))
    call get_rat(atoms,rat)
    call fingerprint_power(atoms%nat, ng/4, atoms%cellvec,rat,atoms%rcov,symfunc%y)
    deallocate(rat)
    !-------------------------------------------------------------------------------------
    !do ig=1,ann_arr%ann(1)%nn(0)
    !if(parini%iverbose>=2) then
    !endif
    end associate
end subroutine symmetry_functions_driver_stefan
!**************************************************************************************************
!program test_fplib_periodic
!  implicit real*8 (a-h,o-z)
!  parameter(nat=20,natx_sphere=50,ntypat=3,lseg=4,nconf=351)
!  dimension alat(3, 3),iassign(nat),cost(nat,nat)
!  dimension rxyz(3, nat),xyzred(3,nat),rcov(nat)
!  dimension energy(nconf) 
!  dimension fpall(lseg*natx_sphere,nat,nconf)
!  dimension nat_type(ntypat),nat_type_ref(ntypat),typ_rcov(ntypat)
!  character(len=2), dimension(ntypat) :: symb, symb_ref
!  character(len=100) :: f1, f2, filetype
!
!
!  if (lseg.eq.1) then 
!      write(*,*) 'S orbitals used'
!  else if (lseg.eq.4) then 
!      write(*,*) 'S+P orbitals used'
!  else
!      write(*,*) 'wrong lseg'
!  endif
!
!  write(*,*) "It is assumed that only ",ntypat," types of atoms can be found in the system"
!  
!
!
!  filetype = 'vasp'
!
!
!
!  if (trim(filetype)=="vasp") then
!     open(11, file='en.dat', status="old")
!     do iconf= 1, nconf
!        read(11,*) energy(iconf)
!     end do
!     close(11)
!
!     do iconf = 1, nconf
!        write(f1, '(I4)') iconf
!!        write(f1, '(I4)') 102
!        f2 = 'Cell_'//trim(adjustl(f1))//'.vasp'
! !       write(f1, '(I4.2)') iconf
! !       f2 = 'P.'//trim(adjustl(f1))
!        !write(f1, '(I5.5)') iconf
!        !f2 = 'poslocm_'//trim(adjustl(f1))//'.ascii.vasp'
!        write(*,*) 'READING ',trim(f2)
!        open(10, file=trim(f2), status="old")
!        read(10,*) 
!        read(10,*) 
!        do i = 1, 3
!           read(10, *) (alat(j, i), j = 1, 3)
!        end do
!        if (iconf.eq.1) then
!          read(10, *) (symb_ref(i), i = 1, ntypat)
!          read(10, *) (nat_type_ref(i), i = 1, ntypat)
!          do i=1,ntypat
!          symb(i)=symb_ref(i)
!          nat_type(i)=nat_type_ref(i)
!          enddo
!
!        else
!          read(10, *) (symb(i), i = 1, ntypat)
!          read(10, *) (nat_type(i), i = 1, ntypat)
!          do  i = 1, ntypat
!          if (symb_ref(i).ne.symb(i)) stop 'inconsistent symbols'
!          if (nat_type_ref(i).ne.nat_type(i))  stop 'inconsistent ntype'
!          enddo
!        endif
!        read(10, *) 
!        do i = 1, nat
!           read(10, *) (xyzred(j, i), j = 1, 3)
!        end do
!        close(10)
!
!!      transform from fractional to cartesian coordinates and to atomic units
!        convert=1.d0/0.52917720859d0
!        call frac2cart(nat, alat, xyzred, rxyz, convert)
!!        write(*,*) 'ALAT'
!        do i=1,3
!        do j=1,3
!        alat(i,j) = alat(i,j)*convert
!        enddo
!!        write(*,*) (alat(i,j),j=1,3)
!        enddo
!
!!    Assign the covalent radii
!        do i = 1, ntypat
!           call sym2rcov(symb(i), typ_rcov(i))
!        end do
!        
!        k = 0
!        do i = 1, ntypat
!          do j = 1, nat_type(i)
!            k = k + 1
!            rcov(k) = typ_rcov(i)
!          end do
!        end do
!
!        call fingerprint_periodic(nat, natx_sphere, lseg, alat, rxyz, rcov, fpall(1,1,iconf))
!     end do
!
!  end if  ! VASP
!
!  do iconf=1,nconf
!  do jconf=iconf,nconf
!
!
!  tmin=1.d100
!  tmax=-1.d100
!  tav=0.d0
!  do iat=1,nat
!  do jat=1,nat
!
!!    ti=0.d0
!!    tj=0.d0
!!    tij=0.d0
!!    do l=1,lseg*natx_sphere
!!    tij=tij+fpall(l,iat,iconf)*fpall(l,jat,jconf)
!!    ti=ti+fpall(l,iat,iconf)*fpall(l,iat,jconf)
!!    tj=tj+fpall(l,jat,iconf)*fpall(l,jat,jconf)
!!    enddo
!!    tt=1.d0-tij/sqrt(ti*tj)
!
!    tt=0.d0
!    do l=1,lseg*natx_sphere
!    tt=tt+(fpall(l,iat,iconf)-fpall(l,jat,jconf))**2
!    enddo
!    tt=sqrt(tt)
!    !tt=log(tt+1.d-12)
!    tmin=min(tmin,tt)
!    tmax=max(tmax,tt)
!    tav=tav+tt
!
!!    cost(iat,jat)=-exp(-tt*1.d4)
!!    cost(iat,jat)=-1.d0/(tt+1.d-6)
!    cost(iat,jat)=tt
!  enddo
!  enddo
!  tav=tav/nat**2
!
!  call apc(nat, cost, iassign, fpd)
!
!  tt=-1.d100
!  do iat=1,nat
!  tt=max(tt,cost(iat,iassign(iat)))
!  enddo
!
!  de=abs(energy(iconf)-energy(jconf))
!  write(10,'(i4,i4,3(2x,e14.7),3(1x,e12.5))') iconf,jconf,de,fpd,tt,tmin,tav,tmax
!
!
!  enddo
!  enddo
!  
!
!end program test_fplib_periodic

subroutine fingerprint_power(nat, npl, alat, rxyz, rcov, fpall)
  implicit real*8 (a-h,o-z)
  integer, parameter:: natx_sphere=500,lseg=4
  integer, parameter:: nwork=100
  real(8):: workalat(nwork) 
  real(8):: rxyz_sphere(3, natx_sphere),rcov_sphere(natx_sphere)
  real(8):: fpall(4,npl,nat),amplitude(natx_sphere)
  real(8):: rxyz(3,nat),rcov(nat)
  real(8):: alat(3, 3),alatalat(3,3),eigalat(3),aa(4,4),aaev(4)
  real(8), allocatable   :: om(:,:,:,:) , power(:,:,:,:)!,work(;),eval(:)
  

! parameters for cutoff function
  width_cutoff=5.d0 
  nex_cutoff=3
  radius_cutoff=sqrt(2.d0*nex_cutoff)*width_cutoff
  radius_cutoff2=radius_cutoff**2
  factor_cutoff=1.d0/(2.d0*nex_cutoff*width_cutoff**2)
!  write(*,*) 'width_cutoff,radius_cutoff',width_cutoff,radius_cutoff

  do i=1,3
  do j=1,3
  alatalat(i,j)=alat(1,i)*alat(1,j)+alat(2,i)*alat(2,j)+alat(3,i)*alat(3,j)
  enddo
  enddo
  call dsyev('N', 'L', 3, alatalat, 3, eigalat, workalat, nwork, info)
!  write(*,*) 'alat EVals',eigalat
!  write(*,*) 'ixyzmax',int(sqrt(1.d0/eigalat(1))*radius_cutoff)
  ixyzmax=0 !int(sqrt(1.d0/eigalat(1))*radius_cutoff) + 1

! loop over all center atoms
  natsmax=0
  natsmin=1000000
           !write(101,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  do iat = 1, nat
           !write(101,*) iat,' -----------------'
           !write(*,*) iat,' ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
           !write(*,*) iat,' ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
           !write(*,*) iat,' ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        
     nat_sphere=0
     do jat = 1, nat   
         do ix = -ixyzmax,ixyzmax
           do iy = -ixyzmax,ixyzmax
             do iz = -ixyzmax,ixyzmax
                xj = rxyz(1, jat) + ix*alat(1,1)+iy*alat(1,2)+iz*alat(1,3)
                yj = rxyz(2, jat) + ix*alat(2,1)+iy*alat(2,2)+iz*alat(2,3)
                zj = rxyz(3, jat) + ix*alat(3,1)+iy*alat(3,2)+iz*alat(3,3)
                dist2 = (xj-rxyz(1, iat))**2+(yj-rxyz(2, iat))**2+(zj-rxyz(3, iat))**2
              
              
                if (dist2.le.radius_cutoff2) then 
!                    write(*,*) 'dist ',jat,ix,iy,iz,sqrt(dist2)
!                    write(*,*) xj,yj,zj
!                    write(*,*) rxyz(1, jat),rxyz(2, jat),rxyz(3, jat)
                    nat_sphere=nat_sphere+1
                    if (nat_sphere.gt.natx_sphere) stop 'enlarge natx_sphere'
                    amplitude(nat_sphere)=(1.d0 - dist2*factor_cutoff)**nex_cutoff
!                    write(*,*) 'amplitude', nat_sphere,amplitude(nat_sphere)
                    rxyz_sphere(1,nat_sphere)=xj
                    rxyz_sphere(2,nat_sphere)=yj
                    rxyz_sphere(3,nat_sphere)=zj
                    rcov_sphere(nat_sphere)=rcov(jat)
                endif
                if (dist2.eq.0.d0) ll=nat_sphere
             enddo
           enddo
         enddo
     enddo
     natsmin=min(natsmin,nat_sphere)
     natsmax=max(natsmax,nat_sphere)

! set up big overlap matrix
   
          nid=nat_sphere*lseg
         
          allocate(om(lseg,nat_sphere,lseg,nat_sphere))
!          allocate(eval(nid))
          allocate(power(lseg,nat_sphere,4,npl+1))
!          allocate(work(nid))
!          call  create_om(nat_sphere,rxyz_sphere,rcov_sphere,nid,om)
          if (lseg.eq.1) then
          call  create_om_1_alborz(nat_sphere,rxyz_sphere,rcov_sphere,om)
          call mltampl_1_alborz(nat_sphere,amplitude,om)
          else if (lseg.eq.4) then
          call  OLDcreate_om_4(nat_sphere,rxyz_sphere,rcov_sphere,om)

!          call Wblock(nat_sphere,lseg,om)
!     stop
          call mltampl_4_alborz(nat_sphere,amplitude,om)
          else
              stop 'wrong lseg'
          endif

          tt=0.d0
          do i=1,nat_sphere
          do l=1,lseg
          do j=i,nat_sphere
          do k=1,lseg
          tt=max(tt,(om(k,i,l,j)-om(l,j,k,i))**2)
          if ((om(k,i,l,j)-om(l,j,k,i))**2.gt.1.d-6) write(*,*) i,j,om(k,i,l,j),om(l,j,k,i)
          enddo
          enddo
          enddo
          enddo
          if (tt.gt.1.d-6) write(*,*) 'max dev symmetry',tt

!          do j=1,nid
!            sum=0.d0
!            do i=1,nid
!            sum=sum+om(i,j)
!            enddo
!            sum=1.d0/sum
!            do i=1,nid
!            om(i,j)=om(i,j)*sum
!            enddo
!          enddo




          ipl=1
            do jat=1,nat_sphere
            do l=1,lseg
            power(l,jat,1,ipl)=om(l,jat,1,ll)
            power(l,jat,2,ipl)=om(l,jat,2,ll)
            power(l,jat,3,ipl)=om(l,jat,3,ll)
            power(l,jat,4,ipl)=om(l,jat,4,ll)
            enddo
            enddo
          do ipl=2,npl+1  
            call dgemm('N','N', nid, 4, nid, 1.d0, om, nid,  power(1,1,1,ipl-1), nid, 0.d0,power(1,1,1,ipl), nid)
          enddo

        !do ipl=2,npl+1  
        !  do i=1,4
        !  do lat=1,nat_sphere
        !  do l=1,4
        !  power(l,lat,i,ipl)=0.d0
        !  do kat=1,nat_sphere
        !  do k=1,4
        !  power(l,lat,i,ipl) = power(l,lat,i,ipl) + om(l,lat,k,kat)*power(k,kat,i,ipl-1)
        !  enddo
        !  enddo
        !  enddo
        !  enddo
        !  enddo
        !enddo

          do ipl=1,npl
          !fpall(1,ipl,iat)=1.d0/power(1,ll,1,ipl+1)
          !t1=power(2,ll,1,ipl+1)**2+power(3,ll,1,ipl+1)**2+power(4,ll,1,ipl+1)**2
          !!t2=power(1,ll,2,ipl+1)**2+power(1,ll,3,ipl+1)**2+power(1,ll,4,ipl+1)**2
          !!write(*,*) 't1,t2 ',t1,t2
          !fpall(2,ipl,iat)=1.d0/(1.d0+sqrt(t1))
          !fpall(3,ipl,iat)=1.d0/(power(2,ll,2,ipl+1)+ power(3,ll,3,ipl+1) + power(4,ll,4,ipl+1))

          do j=1,4
          do i=1,4
          aa(i,j)=power(i,ll,j,ipl+1)
          enddo
          enddo
          call dsyev('N', 'L', 4, aa, 4, aaev, workalat, nwork, info)
          do i=1,4
          fpall(i,ipl,iat)=1.d0/aaev(i)
          enddo

          enddo
      !    write(*,*) iat
      !    do i=1,4
      !    write(*,'(10(1x,e9.2))') (fpall(i,ipl,iat),ipl=1,npl)
      !    enddo

    !     deallocate(work)
    !     lwork=nid**2
    !     allocate(work(lwork))
    !     call DSYEV('N','L',nid,om,nid,eval,work,lwork,info)
    !     if (info.ne.0) stop 'info OM diagonalisation'
    !     write(*,*) ' EVAL '
    !     do i=1,nid
    !     write(*,*) i,eval(i)
    !     enddo
    !     write(77,'(a,21(e8.1))') 'ratio ',(eval(i)/eval(nid),i=nid-1,nid-20,-1)
    !     stop



          deallocate(om)
          !deallocate(eval)
          deallocate(power)
          !deallocate(work)

!!          write(*,*) 'fpall,iat=',iat,nid,nat_sphere,lseg
!          do l=1,3
!          write(*,'(10(2x,e10.3))') (fpall(l,ipl,iat),ipl=1,npl)
!          enddo
!          write(*,*)     
     
  end do
           write(*,*) 'min,max number of atoms in sphere ',natsmin,natsmax
  !         stop

end subroutine fingerprint_power

subroutine fingerprint_periodic(nat, natx_sphere, lseg, alat, rxyz, rcov, fpall)
  implicit real*8 (a-h,o-z)
  parameter(nwork=100)
  dimension workalat(nwork) 
  dimension rxyz_sphere(3, natx_sphere),rcov_sphere(natx_sphere)
  dimension fpall(lseg*natx_sphere,nat),fp(lseg*natx_sphere),amplitude(natx_sphere)
  dimension rxyz(3,nat),rcov(nat)
  dimension alat(3, 3),alatalat(3,3),eigalat(3)
  allocatable   :: om(:,:) , work(:)
  

! parameters for cutoff function
  width_cutoff=3.d0
  nex_cutoff=3
  radius_cutoff=sqrt(2.d0*nex_cutoff)*width_cutoff
  radius_cutoff2=radius_cutoff**2
  factor_cutoff=1.d0/(2.d0*nex_cutoff*width_cutoff**2)
!  write(*,*) 'width_cutoff,radius_cutoff',width_cutoff,radius_cutoff

  do i=1,3
  do j=1,3
  alatalat(i,j)=alat(1,i)*alat(1,j)+alat(2,i)*alat(2,j)+alat(3,i)*alat(3,j)
  enddo
  enddo
  call dsyev('N', 'L', 3, alatalat, 3, eigalat, workalat, nwork, info)
!  write(*,*) 'alat EVals',eigalat
!  write(*,*) 'ixyzmax',int(sqrt(1.d0/eigalat(1))*radius_cutoff)
  !ixyzmax= int(sqrt(1.d0/eigalat(1))*radius_cutoff) + 1
  ixyzmax=0

! loop over all center atoms
  natsmax=0
  natsmin=1000000
  do iat = 1, nat
        
     nat_sphere=0
     do jat = 1, nat
         do ix = -ixyzmax,ixyzmax
           do iy = -ixyzmax,ixyzmax
             do iz = -ixyzmax,ixyzmax
                xj = rxyz(1, jat) + ix*alat(1,1)+iy*alat(1,2)+iz*alat(1,3)
                yj = rxyz(2, jat) + ix*alat(2,1)+iy*alat(2,2)+iz*alat(2,3)
                zj = rxyz(3, jat) + ix*alat(3,1)+iy*alat(3,2)+iz*alat(3,3)
                dist2 = (xj-rxyz(1, iat))**2+(yj-rxyz(2, iat))**2+(zj-rxyz(3, iat))**2

                if (dist2.le.radius_cutoff2) then
!                    write(*,*) 'dist ',jat,ix,iy,iz,sqrt(dist2)
!                    write(*,*) xj,yj,zj
!                    write(*,*) rxyz(1, jat),rxyz(2, jat),rxyz(3, jat)
                    nat_sphere=nat_sphere+1
                    if (nat_sphere.gt.natx_sphere) then
                        write(*,'(a,2i5)') 'enlarge natx_sphere ',nat_sphere,natx_sphere
                        stop
                    endif

                    amplitude(nat_sphere)=(1.d0 - dist2*factor_cutoff)**nex_cutoff
                    rxyz_sphere(1,nat_sphere)=xj
                    rxyz_sphere(2,nat_sphere)=yj
                    rxyz_sphere(3,nat_sphere)=zj
                    rcov_sphere(nat_sphere)=rcov(jat)
                endif
             enddo
           enddo
         enddo
     enddo
     natsmin=min(natsmin,nat_sphere)
     natsmax=max(natsmax,nat_sphere)

! set up big overlap matrix
          nid=nat_sphere*lseg
          allocate(om(nid,nid))
!          call  create_om(nat_sphere,rxyz_sphere,rcov_sphere,nid,om)
          if (lseg.eq.1) then
          call  create_om_1_alborz(nat_sphere,rxyz_sphere,rcov_sphere,om)
          call mltampl_1_alborz(nat_sphere,amplitude,om)
          else if (lseg.eq.4) then
          call  OLDcreate_om_4(nat_sphere,rxyz_sphere,rcov_sphere,om)

!          call Wblock(nat_sphere,lseg,om)
!     stop
          call mltampl_4_alborz(nat_sphere,amplitude,om)
          else
              stop 'wrong lseg'
          endif

          tt=0.d0
          do i=1,nid
          do j=i,nid
          tt=max(tt,(om(i,j)-om(j,i))**2)
          if ((om(i,j)-om(j,i))**2.gt.1.d-6) write(*,*) i,j,om(i,j),om(j,i)
          enddo
          enddo
          if (tt.gt.1.d-6) write(*,*) 'max dev symmetry',tt

          lwork=max(1,3*nid-1)
          allocate(work(lwork))
          call DSYEV('N','L',nid,om,nid,fp,work,-1,info)
          if (info.ne.0) stop 'info query'
          lwork=nint(work(1))
          deallocate(work)
        
          allocate(work(lwork))
          call DSYEV('N','L',nid,om,nid,fp,work,lwork,info)
          if (info.ne.0) stop 'info OM diagonalisation'
          deallocate(work,om)


          do i=1,nid
          fpall(i,iat)=fp(nid+1-i)
          enddo
          do i=nid+1,lseg*natx_sphere
          fpall(i,iat)=0.d0
          enddo

!          write(*,*) 'fpall,iat=',iat,nid,nat_sphere,lseg
!          write(*,'(10(2x,e10.3))') (fpall(i,iat),i=1,lseg*natx_sphere)
     
          if (fp(1).lt.-1.d-12) then 
              write(*,*) fp(1)
              stop 'error negative EV'
          endif

  end do
           write(*,*) 'min,max number of atoms in sphere ',natsmin,natsmax

end subroutine fingerprint_periodic

          subroutine Wblock(nat_sphere,lseg,om)
          implicit real*8 (a-h,o-z)
          dimension om(lseg,nat_sphere,lseg,nat_sphere)
          do i=1,nat_sphere
          do j=1,nat_sphere
          write(*,*) i,j
          do k=1,lseg
          write(*,'(4(1x,e12.5))') (om(l,i,k,j),l=1,lseg)
          enddo
          enddo
          enddo
          return
          end

    subroutine mltampl_4_alborz(nat,amplitude,om)
    implicit real*8 (a-h,o-z)
    dimension amplitude(nat),om(4,nat,4,nat)
          do i=1,nat
          do j=1,nat
          do il=1,4
          do jl=1,4
          om(il,i,jl,j)=om(il,i,jl,j)*amplitude(i)*amplitude(j)
          enddo
          enddo
          enddo
          enddo
      end subroutine mltampl_4_alborz

    subroutine mltampl_1_alborz(nat,amplitude,om)
    implicit real*8 (a-h,o-z)
    dimension amplitude(nat),om(nat,nat)
          do i=1,nat
!          write(*,'(10(1x,e9.2))') (om(i,j),j=1,nat)
          do j=1,nat
!          write(*,'(i4,i4,3(e14.7))') i,j,om(i,j),amplitude(i),amplitude(j)
          om(i,j)=om(i,j)*amplitude(i)*amplitude(j)
          enddo
          enddo
      end subroutine mltampl_1_alborz



subroutine frac2cart(nat, alat, xyzred, rxyz, convert)
  implicit real*8 (a-h,o-z)
  dimension alat(3,3), xyzred(3,nat), rxyz(3,nat)

  do iat=1,nat
    do i = 1, 3
       t = 0
       do j = 1, 3
          t = t + xyzred(j,iat) * alat(i, j)
       end do
       rxyz(i,iat) = t*convert
    end do
  enddo

end subroutine frac2cart



         subroutine create_om_1_alborz(nat,rxyz,rcov,om)
         implicit real*8 (a-h,o-z)
         dimension rxyz(3,nat),rcov(nat),om(nat,nat)

           ! Gaussian overlap
           !  <sj|si>
           do iat=1,nat
              xi=rxyz(1,iat)  
              yi=rxyz(2,iat) 
              zi=rxyz(3,iat) 

              do jat=1,nat

                 d2=(rxyz(1,jat) -xi)**2 +(rxyz(2,jat)-yi)**2+(rxyz(3,jat)-zi)**2
                 r=.5d0/(rcov(iat)**2 + rcov(jat)**2)
                 om(jat,iat)= sqrt(4.d0*r*(rcov(iat)*rcov(jat)))**3 * exp(-d2*r)

                 !xji=rxyz(1,jat) - xi
                 !yji=rxyz(2,jat) - yi 
                 !zji=rxyz(3,jat) - zi
                 !d2=xji*xji + yji*yji + zji*zji
                 !scov=rcov(iat)**2 + rcov(jat)**2
                 !factor=(sqrt(2.d0*rcov(iat)*rcov(jat)/scov))**3
                 !facp=1.d0/sqrt(factor*scov)
                 !arg=d2*.5d0/scov
                 !fexp=factor*exp(-arg)
           !  <sj|si>
                 !om(jat,iat)= fexp
              enddo
           enddo
  end subroutine create_om_1_alborz


         subroutine create_om_4_alborz(nat,rxyz,rcov,om)
         implicit real*8 (a-h,o-z)
         dimension rxyz(3,nat),rcov(nat),om(4,nat,4,nat)


           ! Gaussian overlap
           do iat=1,nat
              xi=rxyz(1,iat)  
              yi=rxyz(2,iat) 
              zi=rxyz(3,iat) 

              do jat=1,nat
                 xji=rxyz(1,jat) - xi
                 yji=rxyz(2,jat) - yi 
                 zji=rxyz(3,jat) - zi
                 d2=xji*xji + yji*yji + zji*zji
                 scov=rcov(iat)**2 + rcov(jat)**2
                 factor=(sqrt(2.d0*rcov(iat)*rcov(jat)/scov))**3
                 facp=sqrt(scov)/(rcov(iat)*rcov(jat))
                 arg=d2*.5d0/scov
                 fexp=factor*exp(-arg)
           !  <sj|si>
                 om(1,jat,1,iat)= fexp
           !  <pj|si>
                    om(2,jat,1,iat)= xji*(facp*fexp*rcov(jat)**2)/scov
                    om(3,jat,1,iat)= yji*(facp*fexp*rcov(jat)**2)/scov
                    om(4,jat,1,iat)= zji*(facp*fexp*rcov(jat)**2)/scov
                    om(1,jat,2,iat)=-xji*(facp*fexp*rcov(iat)**2)/scov
                    om(1,jat,3,iat)=-yji*(facp*fexp*rcov(iat)**2)/scov
                    om(1,jat,4,iat)=-zji*(facp*fexp*rcov(iat)**2)/scov
            ! <pj|pi> 
                    om(2,jat,2,iat)= fexp*(1.d0 - xji*xji/scov)
                    om(3,jat,2,iat)= (fexp/scov) * yji*xji
                    om(4,jat,2,iat)= (fexp/scov) * zji*xji
                    om(2,jat,3,iat)= (fexp/scov) * xji*yji
                    om(3,jat,3,iat)= fexp*(1.d0 -  yji*yji/scov)
                    om(4,jat,3,iat)= (fexp/scov) * zji*yji
                    om(2,jat,4,iat)= (fexp/scov) * xji*zji
                    om(3,jat,4,iat)= (fexp/scov) * yji*zji
                    om(4,jat,4,iat)= fexp*(1.d0 - zji*zji/scov)
              enddo
           enddo

  end subroutine create_om_4_alborz


         subroutine OLDcreate_om_4(nat,rxyz,rcov,om)
         implicit real*8 (a-h,o-z)
         dimension rxyz(3,nat),rcov(nat),om(4,nat,4,nat)


           ! Gaussian overlap
           do iat=1,nat
              xi=rxyz(1,iat)  
              yi=rxyz(2,iat) 
              zi=rxyz(3,iat) 

              do jat=1,nat
                 xji=rxyz(1,jat) - xi
                 yji=rxyz(2,jat) - yi 
                 zji=rxyz(3,jat) - zi
                 d2=xji*xji + yji*yji + zji*zji
                 r=.5d0/(rcov(iat)**2 + rcov(jat)**2)
           !  <sj|si>
                 om(1,jat,1,iat)= sqrt(4.d0*r*(rcov(iat)*rcov(jat)))**3 * exp(-d2*r)
           !  <pj|si>
                    sji= sqrt(4.d0*r*(rcov(jat)*rcov(iat)))**3 * exp(-d2*r)
                    tt= sqrt(8.d0) *rcov(jat)*r * sji
                    om(2,jat,1,iat)= tt * xji 
                    om(3,jat,1,iat)= tt * yji 
                    om(4,jat,1,iat)= tt * zji 
                    tt= sqrt(8.d0) *rcov(iat)*r * sji
                    om(1,jat,2,iat)=-tt * xji
                    om(1,jat,3,iat)=-tt * yji
                    om(1,jat,4,iat)=-tt * zji
            ! <pj|pi> 
                    tt = -8.d0*rcov(iat)*rcov(jat) * r*r * sji 
                    om(2,jat,2,iat)= tt *(xji* xji - .5d0/r) 
                    om(3,jat,2,iat)= tt *(yji* xji         ) 
                    om(4,jat,2,iat)= tt *(zji* xji         ) 
                    om(2,jat,3,iat)= tt *(xji* yji         ) 
                    om(3,jat,3,iat)= tt *(yji* yji - .5d0/r) 
                    om(4,jat,3,iat)= tt *(zji* yji         ) 
                    om(2,jat,4,iat)= tt *(xji* zji         ) 
                    om(3,jat,4,iat)= tt *(yji* zji         ) 
                    om(4,jat,4,iat)= tt *(zji* zji - .5d0/r) 
              enddo
           enddo

  end subroutine OLDcreate_om_4


         subroutine create_om(nat,rxyz,rcov,nid,om)
         implicit real*8 (a-h,o-z)
         dimension rxyz(3,nat),rcov(nat),om(nid,nid)

  if(nid .ne. nat .and. nid .ne. 4*nat) stop ' nid should be either nat or  4*nat '

  do i=1,nid
  do j=1,nid
  om(i,j)=0.d0
  enddo
  enddo

           ! Gaussian overlap
           !  <sj|si>
           do iat=1,nat
              xi=rxyz(1,iat)  
              yi=rxyz(2,iat) 
              zi=rxyz(3,iat) 

              do jat=iat,nat
                 d2=(rxyz(1,jat) -xi)**2 +(rxyz(2,jat)-yi)**2+(rxyz(3,jat)-zi)**2
                 r=.5d0/(rcov(iat)**2 + rcov(jat)**2)
                 om(jat,iat)=om(jat,iat) + sqrt(4.d0*r*(rcov(iat)*rcov(jat)))**3 * exp(-d2*r)
              enddo
           enddo

  !!  so far only s-s have been calculated  
  if(nid == 4*nat) then  ! both s and p (nid = 4nat)


              !  <s|p>
              do iat=1,nat
                 xi=rxyz(1,iat) 
                 yi=rxyz(2,iat) 
                 zi=rxyz(3,iat) 

                 do jat=1,nat   ! NOTE: do not use  jat=iat,nat becase all elements are on the same side of the diagonal

                    xji=rxyz(1,jat) - xi
                    yji=rxyz(2,jat) - yi 
                    zji=rxyz(3,jat) - zi

                    d2=xji*xji + yji*yji + zji*zji
                    r=.5d0/(rcov(jat)**2 + rcov(iat)**2)

                    sji= sqrt(4.d0*r*(rcov(jat)*rcov(iat)))**3 * exp(-d2*r)

                    !  <pj|si>
                    tt= sqrt(8.d0) *rcov(jat)*r * sji

                    om(1+nat + (jat-1)*3 ,iat )=  om(1+nat + (jat-1)*3 ,iat ) + tt * xji 
                    om(2+nat + (jat-1)*3 ,iat )=  om(2+nat + (jat-1)*3 ,iat ) + tt * yji 
                    om(3+nat + (jat-1)*3 ,iat )=  om(3+nat + (jat-1)*3 ,iat ) + tt * zji 

                    !! !  <sj|pi> no need, because they are on the other side of the diagonal of the symmetric matrix
                    !!     tt=-sqrt8 *rcov(iat)*r * sji

                    !!     om(jat, 1+nat + (iat-1)*3 )=  om(jat, 1+nat + (iat-1)*3 ) + tt * xji 
                    !!     om(jat, 2+nat + (iat-1)*3 )=  om(jat, 2+nat + (iat-1)*3 ) + tt * yji 
                    !!     om(jat, 3+nat + (iat-1)*3 )=  om(jat, 3+nat + (iat-1)*3 ) + tt * zji 

                 enddo
              enddo


              ! <pj|pi> 
              do iat=1,nat
                 xi=rxyz(1,iat)
                 yi=rxyz(2,iat)
                 zi=rxyz(3,iat)

                 do jat=iat,nat

                    xji=rxyz(1,jat) - xi
                    yji=rxyz(2,jat) - yi 
                    zji=rxyz(3,jat) - zi

                    d2=xji*xji + yji*yji + zji*zji
                    r=.5d0/(rcov(jat)**2 + rcov(iat)**2)

                    sji= sqrt(4.d0*r*(rcov(jat)*rcov(iat)))**3 * exp(-d2*r)

                    igto=nat+1 +(iat-1)*3 
                    jgto=nat+1 +(jat-1)*3

                    tt = -8.d0*rcov(iat)*rcov(jat) * r*r * sji 

                    om(jgto   , igto  )=  om(jgto   , igto  ) + tt *(xji* xji - .5d0/r) 
                    om(jgto   , igto+1)=  om(jgto   , igto+1) + tt *(yji* xji         ) 
                    om(jgto   , igto+2)=  om(jgto   , igto+2) + tt *(zji* xji         ) 
                    om(jgto+1 , igto  )=  om(jgto+1 , igto  ) + tt *(xji* yji         ) 
                    om(jgto+1 , igto+1)=  om(jgto+1 , igto+1) + tt *(yji* yji - .5d0/r) 
                    om(jgto+1 , igto+2)=  om(jgto+1 , igto+2) + tt *(zji* yji         ) 
                    om(jgto+2 , igto  )=  om(jgto+2 , igto  ) + tt *(xji* zji         ) 
                    om(jgto+2 , igto+1)=  om(jgto+2 , igto+1) + tt *(yji* zji         ) 
                    om(jgto+2 , igto+2)=  om(jgto+2 , igto+2) + tt *(zji* zji - .5d0/r) 

                 enddo
              enddo

  endif  ! both s and p 
         end subroutine create_om
