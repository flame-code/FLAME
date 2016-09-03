subroutine get_fp_gauss(nat, ntypat, natx_sphere, typat, lseg, width_cutoff, nex_cutoff, alat, rxyz, rcov, fp)
  !use fingerprint, only :&
  !     & fp_17_natx_sphere  => natx_sphere,  &
  !     & fp_17_lseg         => lseg,         &
  !     & fp_17_width_cutoff => width_cutoff, &
  !     & fp_17_nex_cutoff   => nex_cutoff

  implicit none
  integer, intent(in) :: nat, ntypat, natx_sphere, lseg
  integer, dimension(nat), intent(in) :: typat
  real(8), intent(in) :: width_cutoff, nex_cutoff
  real(8), dimension(3,3), intent(in) :: alat
  real(8), dimension(3,nat), intent(in) :: rxyz
  real(8), dimension(nat), intent(in) :: rcov

  real(8), dimension(lseg*(ntypat+1), nat), intent(out) :: fp

  ! local variables
  integer, parameter :: nwork = 100
  integer :: i, j, k, l, ixyzmax, ix, iy, iz, iat, jat, it
  integer :: natsmax, natsmin, nat_sphere, nid, info, lwork
  integer, dimension(lseg*natx_sphere) :: ind_small
  real(8) :: radius_cutoff, radius_cutoff2, factor_cutoff
  real(8) :: xj, yj, zj, dist2, tt
  real(8), dimension(natx_sphere) :: amplitude
  real(8), dimension(nwork) ::  workalat
  real(8), dimension(lseg*natx_sphere) :: fpp
  real(8), dimension(3, natx_sphere) :: rxyz_sphere
  real(8), dimension(natx_sphere) :: rcov_sphere
  real(8), dimension(3,3) :: alatalat
  real(8), dimension(3) :: eigalat
  real(8), dimension(lseg*(ntypat+1), lseg*(ntypat+1)) :: omsa, omsb, omsaa, omsbb
  real(8), dimension(:,:), allocatable :: oms, om
  real(8), dimension(:), allocatable :: work

  !print*, 'n w', nex_cutoff, width_cutoff

  radius_cutoff = sqrt(2.d0*nex_cutoff)*width_cutoff
  radius_cutoff2 = radius_cutoff**2
  factor_cutoff = 1.d0/(2.d0*nex_cutoff*width_cutoff**2)
  print*, 'FP radius_cutoff', radius_cutoff

  do i = 1, 3
     do j = 1, 3
        alatalat(i,j) = alat(1,i)*alat(1,j)+alat(2,i)*alat(2,j)+alat(3,i)*alat(3,j)
     end do
  end do

  call dsyev('N', 'L', 3, alatalat, 3, eigalat, workalat, nwork, info)
  ixyzmax= int(sqrt(1.d0/eigalat(1))*radius_cutoff) + 1


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
!!$                    write(*,*) 'dist ',jat,ix,iy,iz,sqrt(dist2)
!!$                    write(*,*) xj,yj,zj
!!$                    write(*,*) rxyz(1, jat),rxyz(2, jat),rxyz(3, jat)
                    nat_sphere=nat_sphere+1
                    if (nat_sphere.gt.natx_sphere) then 
                        stop 'enlarge natx_sphere'
                    end if
                    amplitude(nat_sphere)=(1.d0 - dist2*factor_cutoff)**nex_cutoff
                    rxyz_sphere(1,nat_sphere)=xj
                    rxyz_sphere(2,nat_sphere)=yj
                    rxyz_sphere(3,nat_sphere)=zj
                    rcov_sphere(nat_sphere)=rcov(jat)
                    if (jat.eq.iat .and. ix.eq.0 .and. iy.eq.0 .and. iz.eq.0) then
                       it=0
                    else
                       it=typat(jat)
                    end if
                    do l=1,lseg
                       ind_small(lseg*(nat_sphere-1)+l)=it*lseg+l
                    enddo
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
     allocate(oms(nid,nid))
     if (lseg.eq.1) then
        call create_om_1(nat_sphere,rxyz_sphere,rcov_sphere,om)
        call mltampl_1(nat_sphere,amplitude,om)
     else if (lseg.eq.4) then
        call  create_om_4(nat_sphere,rxyz_sphere,rcov_sphere,om)
        call mltampl_4(nat_sphere,amplitude,om)
     else
        stop 'wrong lseg'
     endif

     tt=0.d0
     do i=1,nid
        do j=1,nid
           oms(i,j)=om(i,j)
           tt=max(tt,(om(i,j)-om(j,i))**2)
           if ((om(i,j)-om(j,i))**2.gt.1.d-6) write(*,*) i,j,om(i,j),om(j,i)
        enddo
     enddo
     if (tt.gt.1.d-6) write(*,*) 'max dev symmetry',tt

     lwork=max(1,3*nid-1)
     allocate(work(lwork))
     call DSYEV('N','L',nid,om,nid,fpp,work,-1,info)
     if (info.ne.0) stop 'info query'
     lwork=nint(work(1))
     deallocate(work)

     allocate(work(lwork))
     call DSYEV('V','L',nid,om,nid,fpp,work,lwork,info)
     if (info.ne.0) stop 'info OM diagonalisation'

!!$     do i=1,nid
!!$        fpall(i,iat)=fp(nid+1-i)
!!$     enddo
!!$     do i=nid+1,lseg*natx_sphere
!!$        fpall(i,iat)=0.d0
!!$     enddo


     if (fpp(1).lt.-1.d-12) then 
        write(*,*) fpp(1)
        stop 'error big negative EV'
     endif

     ! Contract 
     do i=1,lseg*(ntypat+1)
        do j=1,lseg*(ntypat+1)
           omsa(j,i)=0.d0
        enddo
     enddo
     do i=1,nid
        do j=1,nid
           omsa(ind_small(i),ind_small(j))=omsa(ind_small(i),ind_small(j)) + om(j,nid)*oms(i,j)*om(i,nid)
        end do
     end do

     ! Contract 
     do i=1,lseg*(ntypat+1)
        do j=1,lseg*(ntypat+1)
           omsb(j,i)=0.d0
        enddo
     enddo
     do i=1,nid
        omsb(ind_small(i),ind_small(i))=omsb(ind_small(i),ind_small(i)) + om(i,nid)*om(i,nid)
     enddo

     do i=1,lseg*(ntypat+1)
        if (omsb(i,i).eq.0.d0) then 
           !write(*,*) 'ZERO diagonal' 
           omsb(i,i)=1.d0
        endif
     enddo

     do i=1,lseg*(ntypat+1)
        do j=1,lseg*(ntypat+1)
           omsaa(j,i)=omsa(j,i)
           omsbb(j,i)=omsb(j,i)
        enddo
     enddo

     
     call  DSYGV(1,'N', 'L', lseg*(ntypat+1), omsa, lseg*(ntypat+1), omsb, lseg*(ntypat+1), fpp, WORK, LWORK, INFO )
     if (info.ne.0) then
        write(*,*) 'DSYGV info ',info
        write(*,*) 'OMSA'
        do i=1,lseg*(ntypat+1)
           write(*,'(16(1x,e9.2))') ( omsaa(i,j),j=1,lseg*(ntypat+1))
        enddo
        write(*,*) 'OMSB'
        do i=1,lseg*(ntypat+1)
           write(*,'(16(1x,e9.2))') ( omsbb(i,j),j=1,lseg*(ntypat+1))
        enddo
        stop 'DSYGV'
     endif

     if (fpp(1).lt.-1.d-12) then
        write(*,*) fpp(1)
        stop 'error big negative EV'
     endif

     do i=1,lseg*(ntypat+1)
        fp(i,iat)=fpp(lseg*(ntypat+1)+1-i)
     enddo


     deallocate(work,om,oms)
  end do
  write(*,*) '# NATX_SPHERE: min,max number of atoms in sphere ',natsmin,natsmax


end subroutine get_fp_gauss

subroutine get_distance_gauss(fp1, fp2, lseg, nat, ntypat, typat, fpd)
  implicit none
  integer, intent(in)  :: lseg, nat, ntypat
  integer, dimension(nat) :: typat
  real(8), dimension(lseg*(ntypat+1), nat) :: fp1, fp2
  real(8), intent(out) :: fpd

  integer :: iat, jat, l, i, j, itypat
  integer, dimension(:), allocatable :: iassign
  real(8) :: tmin, tmax, tt, dt
  real(8), dimension(nat, nat) :: cost
  real(8), dimension(ntypat)   :: dist
  real(8), dimension(:,:), allocatable :: costred

  do itypat = 1, ntypat 
     i = 0
     do iat = 1, nat
        if (typat(iat) == itypat) then
           i = i + 1
           j = 0
           do jat = 1, nat
              if (typat(jat) == itypat) then
                 j = j + 1
                 tt = 0.d0
                 do l = 1, lseg*(ntypat+1)
                    tt = tt + (fp1(l, iat) - fp2(l, jat))**2
                 end do
                 cost(i,j) = tt
              end if
           end do
        end if
     end do
     allocate(costred(i,i))
     allocate(iassign(i))
     costred(1:i,1:i) = cost(1:i,1:i)
     call apc(i, costred, iassign, dt)
     dist(itypat) = dt / i
     deallocate(costred)
     deallocate(iassign)
  end do

  fpd = sum(dist) / ntypat

end subroutine get_distance_gauss

subroutine mltampl_4(nat,amplitude,om)
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
end subroutine mltampl_4

subroutine mltampl_1(nat,amplitude,om)
  implicit real*8 (a-h,o-z)
  dimension amplitude(nat),om(nat,nat)
  do i=1,nat
     !          write(*,'(10(1x,e9.2))') (om(i,j),j=1,nat)
     do j=1,nat
        !          write(*,'(i4,i4,3(e14.7))') i,j,om(i,j),amplitude(i),amplitude(j)
        om(i,j)=om(i,j)*amplitude(i)*amplitude(j)
     enddo
  enddo
end subroutine mltampl_1



subroutine create_om_1(nat,rxyz,rcov,om)
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
     enddo
  enddo
end subroutine create_om_1




subroutine create_om_4(nat,rxyz,rcov,om)
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

end subroutine create_om_4
