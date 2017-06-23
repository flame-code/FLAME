!subroutine fingerprint_gaussmol(nat,rxyz,rcov,fp)
!implicit real*8 (a-h,o-z)
!dimension rxyz(3,nat),fp(nat),rcov(nat)
!real*8, allocatable, dimension(:,:) :: aa,work
!allocate(aa(nat,nat),work(nat,nat))
!
!!! Gaussian interaction
!!       do iat=1,nat
!!         aa(iat,iat)=1.d0
!!         do jat=iat+1,nat
!!           d2=(rxyz(1,iat)-rxyz(1,jat))**2+(rxyz(2,iat)-rxyz(2,jat))**2+(rxyz(3,iat)-rxyz(3,jat))**2
!!           tt=exp(-.5d0*d2/(rcov(iat)+rcov(jat)**2))
!!           aa(iat,jat)=tt
!!           aa(jat,iat)=tt
!!         enddo
!!       enddo
!!       write(667,'(20(1x,e9.2))') (rcov(iat) ,iat=1,nat)
!!       do jat=1,nat
!!       write(667,'(20(1x,e9.2))') (aa(jat,iat) ,iat=1,jat)
!!       enddo
!
!
!! Gaussian overlap
!     do iat=1,nat
!      do jat=iat,nat
!        d2=(rxyz(1,iat)-rxyz(1,jat))**2 +(rxyz(2,iat)-rxyz(2,jat))**2+(rxyz(3,iat)-rxyz(3,jat))**2
!        r=.5d0/(rcov(iat)**2 + rcov(jat)**2)
!        ! with normalized GTOs:
!        aa(jat,iat)=sqrt(2.d0*r*(2.d0*rcov(iat)*rcov(jat)))**3 * exp(-d2*r)
!        enddo
!      enddo
!
!!       write(666,'(20(1x,e9.2))') (rcov(iat) ,iat=1,nat)
!!       do jat=1,nat
!!       write(666,'(20(1x,e9.2))') (aa(jat,iat) ,iat=1,jat)
!!       enddo
!
!!! Check Gaussian overlap against Ali's formula
!!     do iat=1,nat
!!         alphai=.5d0/rcov(iat)**2
!!      do jat=iat,nat
!!         alphaj=.5d0/rcov(jat)**2
!!        d2=(rxyz(1,iat)-rxyz(1,jat))**2 +(rxyz(2,iat)-rxyz(2,jat))**2+(rxyz(3,iat)-rxyz(3,jat))**2
!!        t1=alphai*alphaj 
!!        t2=alphai+alphaj 
!!        ! with normalized GTOs:
!!        write(777,*) iat,jat, aa(jat,iat),sqrt(2.d0*sqrt(t1)/t2)**3 * exp(-t1/t2*d2)
!!        enddo
!!      enddo
!
!
!       call DSYEV('N','L',nat,aa,nat,fp,work,nat**2,info)
!       if (info.ne.0) stop 'info'
!!!       if (iproc.eq.0) write(*,'(a,20(e10.3))') '(MH) fingerprint ',(fp(i),i=1,nat)
!
!       deallocate(aa,work)
!       end subroutine 

!subroutine fingerprint_gaussmol(nat,rxyz,rcov,fp)
subroutine fingerprint_gaussmol(nat,nid,rxyz,rcov,fp)
! calculates an overlap matrix for atom centered GTO of the form:
!    s-type: 1/norm_s  exp(-(1/2)*(r/rcov)**2)
!   px type: 1/norm_p exp(-(1/2)*(r/rcov)**2) x/r  and analageously for py and pz
implicit none !real*8 (a-h,o-z)
integer  nat,nid ,iproc,  info
real*8 :: rxyz(3,nat),fp(nid),rcov(nat),tau(3),alat(3)
real*8, allocatable, dimension(:,:) :: om,work

integer igto,jgto, iat, jat
integer i1,i2,i3, n1, n2, n3  
real*8  cutoff, d2, r
real*8  sji, xi,yi,zi, xji, yji, zji   ,tt 
real*8  sqrt8 ; parameter (sqrt8=sqrt(8.d0))


   ! WARNING! check convergence to ensure that the folloing cutoff is large enough
   !! exp(-0.5*cutoff^2/rcov^2) = 1E-16  ==> cutoff^2 = 2*16*log(10)*rcov^2 ==> cutoff ~=8.5 rcov 
   !cutoff=sqrt(2*16*log(10.d0)*maxval(rcov)**2)
   cutoff=9*maxval(rcov)
     !print*, cutoff; stop

   !with these settings the fingerprints have about 9 correct decimal places
!     if (geocode == 'F') then       ! free boundary conditions
         n1=0 ; n2=0 ; n3=0
!     else if (geocode == 'S') then  ! surface boundary conditions, non-periodic direction i s
!         n1=nint(cutoff/alat(1))
!         n2=0
!         n3=nint(cutoff/alat(3))
!     else if (geocode == 'P') then  ! periodic boundary conditions
!         n1=nint(cutoff/alat(1))
!         n2=nint(cutoff/alat(2))
!         n3=nint(cutoff/alat(3))
!     else
!     stop 'unrecognized BC in fingerprint'
!     endif
!     if (n1+n2+n3.gt.30) write(*,*) 'Warning n1,n2,n3 too big ',n1,n2,n3

if(nid .ne. nat .and. nid .ne. 4*nat) stop ' nid should be either nat or  4*nat '


allocate(om(nid,nid),work(nid,nid))
om(:,:)=0.d0

    do i1=-n1,n1
    do i2=-n2,n2
    do i3=-n3,n3
    
       tau(1)=alat(1)*i1
       tau(2)=alat(2)*i2
       tau(3)=alat(3)*i3
    !   if (tau(1)*tau(1) + tau(2)*tau(2)+ tau(3)*tau(3)>cutoff*cutoff) cycle  ! to speedup
    
    ! Gaussian overlap
         !  <sj|si>
          do iat=1,nat
           xi=rxyz(1,iat) + tau(1) 
           yi=rxyz(2,iat) + tau(2)
           zi=rxyz(3,iat) + tau(3)
          
           do jat=iat,nat
             d2=(rxyz(1,jat) -xi)**2 +(rxyz(2,jat)-yi)**2+(rxyz(3,jat)-zi)**2
             r=.5d0/(rcov(iat)**2 + rcov(jat)**2)
             om(jat,iat)=om(jat,iat) + sqrt(4.d0*r*(rcov(iat)*rcov(jat)))**3 * exp(-d2*r)
             enddo
           enddo
    
    enddo !i3
    enddo !i2
    enddo !i1


!!  so far only s-s have been calculated  
if(nid == 4*nat) then  ! both s and p (nid = 4nat)

    do i1=-n1,n1
    do i2=-n2,n2
    do i3=-n3,n3
 
       tau(1)=alat(1)*i1
       tau(2)=alat(2)*i2
       tau(3)=alat(3)*i3

    !  <s|p>
    do iat=1,nat
      xi=rxyz(1,iat) + tau(1)
      yi=rxyz(2,iat) + tau(2)
      zi=rxyz(3,iat) + tau(3)

      do jat=1,nat   ! NOTE: do not use  jat=iat,nat becase all elements are on the same side of the diagonal

        xji=rxyz(1,jat) - xi
        yji=rxyz(2,jat) - yi 
        zji=rxyz(3,jat) - zi

        d2=xji*xji + yji*yji + zji*zji
        r=.5d0/(rcov(jat)**2 + rcov(iat)**2)

        sji= sqrt(4.d0*r*(rcov(jat)*rcov(iat)))**3 * exp(-d2*r)

    !  <pj|si>
        tt= sqrt8 *rcov(jat)*r * sji

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
      xi=rxyz(1,iat) + tau(1)
      yi=rxyz(2,iat) + tau(2)
      zi=rxyz(3,iat) + tau(3)

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

enddo  ! i3 
enddo  ! i2
enddo  ! i1

endif  ! both s and p 



 call DSYEV('N','L',nid,om,nid,fp,work,nid**2,info)
 if (info.ne.0) stop 'info'
! if (iproc.eq.0) write(*,'(a,1000(e10.3))') '(MH) fingerprint ',(fp(i1),i1=1,nid)

deallocate(om,work)
end subroutine fingerprint_gaussmol


       subroutine fpdistance_gaussmol(nid,fp1,fp2,d)
       implicit real*8 (a-h,o-z)
       dimension fp1(nid),fp2(nid)

       d=0.d0
       do i=1,nid
       d = d + (fp1(i)-fp2(i))**2
       enddo
       d=sqrt(d/nid)

       end subroutine

