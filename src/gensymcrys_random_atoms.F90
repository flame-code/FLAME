subroutine random_atom(LATSGP,CRYSSYS,BRAV,NSYMP,NSYM,RSYM,NSYMMAX,NAT_CELL,NAT_KINDS,NKINDS,KINDS,NGUESS,KINDSDIST_MIN,RXYZ,RED_POS_PRIM,LATVEC_PRIM)
!Idea of this subroutine: suppose that sg_ops has already been called and RSYM, etc is provided
!1.) Get symmetry operations 
!2.) Get the random sampling of atomic positions 
!3.) Pick atoms from random sampling, one by one (maybe more than one atom is placed in the cell because of wyckoff points)
!4.) Expand the primitive cell to check if the new atom is acceptable
!5.) accept/reject the position and cycle

implicit none
INTEGER:: LATSGP, CRYSSYS,BRAV, NAT_CELL,NKINDS, NGUESS,KINDS(NAT_CELL),NAT_KINDS(NKINDS),NSYMP,NSYM,NSYMMAX,COUNT_LOOP
integer:: sumkinds,i,iat,k,l,m,iat_cur,iat_cur1,ikind, counter,NPOS,IPOS,iat_cur_oldsym,counter_old,errcount
real(8):: RSYM(4,4,NSYMMAX),RED_POS(3,NGUESS),RED_POS_PRIM(3,NAT_CELL),RXYZ(3,NAT_CELL),LATVEC_PRIM(3,3),KINDSDIST_MIN(NKINDS,NKINDS)
real(8):: rand,SYM_POS_PRIM(3,NSYMMAX),rxyz_tmp(3),POS(3),dist,dist_count
real(8):: transvecall(3,3,3,3)
real(8),allocatable:: rxyzout(:,:,:,:,:)
logical:: fail
!Do some checking first
COUNT_LOOP=0
sumkinds=0
do i=1,nkinds
  sumkinds=sumkinds+nat_kinds(i)
  write(*,*) "nat_kinds(i)",nat_kinds(i)
enddo
if(NAT_CELL.ne.sumkinds) stop "Number of atoms in cell not consistent with number of atom kinds"
if(NAT_CELL.gt.NGUESS) stop "Increase NGUESS for appropriate unit cell sampling"


!1.) allocate RSYM and stuff
!Already done from outside the routine

!2.) call random atom in conventional cell sampling 
1060 continue
write(*,*) "Starting (over) new generation", COUNT_LOOP
if(COUNT_LOOP.gt.20) then
write(*,*) "No structure found that fit the requirements. Exiting"
kinds=0
return
endif
iat_cur=0
kinds=0
errcount=0
write(*,*) "NKINDS",nkinds
write(*,*) "NAT_CELL",nat_cell

call random_incell_main(LATSGP,CRYSSYS,NGUESS,RED_POS)
write(*,*) "Incell called"
do ikind=1,nkinds
dist_count=0
!  write(*,*)"nat_kinds(ikind)",nat_kinds(ikind)
counter=0
  do 
1040 continue
    if(dist_count.gt.NAT_CELL*NGUESS) then
       write(*,*) "Volume probably too small, exiting"
       kinds=0
       return 
    endif

    iat_cur_oldsym=iat_cur
    counter_old=counter
1070 continue
    call random_number(rand)
    i=int(rand*real(NGUESS,8))+1
    POS(:)=RED_POS(:,i)
    call create_unique_positions(RSYM,NSYMP,NSYM,BRAV,POS,SYM_POS_PRIM,NPOS,NSYMMAX)
    
!3.) check each of the new positions that were generated and accept/reject the positions according
!    to the interatomic distance in KINDSDIST_MIN
!    write(*,*) "NPOS",npos
    if(nat_kinds(ikind)-counter-npos.lt.0) then 
!      write(*,*) "errcount", errcount
      errcount=errcount+1
      if (errcount.gt.NAT_CELL*NGUESS) then
         count_loop=count_loop+1
         write(15,*) "going to 1060"
         goto 1060
      endif
      write(15,*) "going to 1070"
      goto 1070
    endif
    do ipos=1,npos
    rxyz_tmp=matmul(LATVEC_PRIM,SYM_POS_PRIM(:,ipos)) 
    fail=.false.     
!    if(iat_cur.eq.0) write(15,*) "going to 1030",iat_cur
!    if(iat_cur.eq.0) goto 1030
    iat_cur1=iat_cur+1 
    allocate(rxyzout(3,iat_cur1,3,3,3)) 
    rxyz(:,iat_cur1)=rxyz_tmp(:)
    kinds(iat_cur1)=ikind
    call expand_gensymcrys(rxyz,rxyzout,transvecall,latvec_prim,iat_cur1)
    write(15,*) "expanded"
       do k=1,3
        do l=1,3
         do m=1,3
          do iat=1,iat_cur1
            if(iat_cur1==iat.and.k==2.and.l==2.and.m==2) goto 1090!This is the case where we check the same atom in the base cell
               dist=0.d0    
               dist=dist+(rxyz_tmp(1)-rxyzout(1,iat,k,l,m))**2
               dist=dist+(rxyz_tmp(2)-rxyzout(2,iat,k,l,m))**2
               dist=dist+(rxyz_tmp(3)-rxyzout(3,iat,k,l,m))**2
               !write(*,*) "dist",dist,KINDSDIST_MIN(ikind,kinds(iat))**2
               if(dist.lt.KINDSDIST_MIN(ikind,kinds(iat))**2) then
                 fail=.true.
                 deallocate(rxyzout)
                 write(15,*) "going to 1030"
                 goto 1030
               endif
            1090 continue
          enddo         
         enddo
        enddo
       enddo 
     deallocate(rxyzout)
1030 continue
       if(fail) then
         dist_count=dist_count+1
         iat_cur=iat_cur_oldsym
         counter=counter_old
         write(15,*) "going to 1040"
         goto 1040
       else
         iat_cur=iat_cur+1
         counter=counter+1
         RXYZ(:,iat_cur)=rxyz_tmp
         RED_POS_PRIM(:,iat_cur)=SYM_POS_PRIM(:,ipos)
         kinds(iat_cur)=ikind
!         write(*,*) "NPOS",npos,counter,iat_cur
         if(ipos==1) write(15,*)"r", RXYZ(:,iat_cur)
         if(ipos==1) write(15,*)"p", RED_POS_PRIM(:,iat_cur)
         if(ipos==1)  write(15,*) "i",iat_cur,counter
         if(counter.ge.nat_kinds(ikind))write(15,*) "going to 1050"
         if(counter.ge.nat_kinds(ikind))goto 1050 
       endif
    enddo
  write(15,*) "Accepted"
  enddo
1050 continue   
enddo

write(*,*) "Finished random atoms"
if(nat_cell.ne.iat_cur)then

write(*,*) "Something wrong with atom numbers",nat_cell,iat_cur
endif

end subroutine


subroutine create_unique_positions(RSYM,NSYMP,NSYM,BRAV,POS,SYM_POS_PRIM,NPOS,NSYMMAX)
implicit none
integer::NSYMP,NSYM,BRAV,NPOS,NSYMMAX,iat,jat
real(8)::RSYM(4,4,NSYMMAX),POS(3),OUT_POS(3),transmat(3,3),SYM_POS_PRIM(3,NSYMMAX)
real(8):: tol,dist,rx1,ry1,rz1,rx2,ry2,rz2
write(15,*) "NEW"
tol=1.d-8
call  conv2prim(BRAV,transmat)
NPOS=1
OUT_POS=POS
OUT_POS=matmul(transmat,OUT_POS)
OUT_POS(1)=modulo(modulo(OUT_POS(1),1.d0),1.d0)
OUT_POS(2)=modulo(modulo(OUT_POS(2),1.d0),1.d0)
OUT_POS(3)=modulo(modulo(OUT_POS(3),1.d0),1.d0)
SYM_POS_PRIM(:,NPOS)=OUT_POS
write(15,*)SYM_POS_PRIM(:,NPOS)

!Outer loop goes over all symmetry operations
do iat=1,NSYM
call APPLY_SYMOP(RSYM(:,:,iat),POS,OUT_POS)
OUT_POS=matmul(transmat,OUT_POS)
OUT_POS(1)=modulo(modulo(OUT_POS(1),1.d0),1.d0)
OUT_POS(2)=modulo(modulo(OUT_POS(2),1.d0),1.d0)
OUT_POS(3)=modulo(modulo(OUT_POS(3),1.d0),1.d0)
rx1=OUT_POS(1);ry1=OUT_POS(2);rz1=OUT_POS(3)
!Inner loop to check if the atoms are projected onto itself
   do jat=1,NPOS
      rx2=SYM_POS_PRIM(1,jat);ry2=SYM_POS_PRIM(2,jat);rz2=SYM_POS_PRIM(3,jat)
      dist=0.d0
      dist=dist+min((rx1-rx2)**2,(rx1-rx2+1.d0)**2,(rx1-rx2-1.d0)**2)      
      dist=dist+min((ry1-ry2)**2,(ry1-ry2+1.d0)**2,(ry1-ry2-1.d0)**2)      
      dist=dist+min((rz1-rz2)**2,(rz1-rz2+1.d0)**2,(rz1-rz2-1.d0)**2)      
      if(dist.lt.tol) goto 1020 
   enddo 
   NPOS=NPOS+1
   SYM_POS_PRIM(:,NPOS)=OUT_POS
   write(15,*)SYM_POS_PRIM(:,NPOS)
1020 continue
enddo

   write(15,*)NPOS

end subroutine


subroutine APPLY_SYMOP(SYMOP,IN_POS,OUT_POS)
implicit none
real(8):: SYMOP(4,4),IN_POS(3),OUT_POS(3),tmp1(4),tmp2(4)
tmp1(1:3)=IN_POS(1:3)
tmp1(4)=1.d0
tmp2=matmul(SYMOP,tmp1)
OUT_POS=tmp2(1:3)
end subroutine


subroutine find_ir_multiplicity(LATSGP,RSYM,NSYMP,NSYM,BRAV,NSYMMAX,NGUESS,NPOS_IRRED,NPOS_IRRED_ARR)
implicit none
integer::NSYMP,NSYM,BRAV,NPOS,NSYMMAX,iat,jat,CRYSSYS,NPOS_IRRED,NPOS_IRRED_ARR(NSYMMAX),LATSGP,NGUESS
integer:: i,j,rest
real(8)::RSYM(4,4,NSYMMAX),POS(3),OUT_POS(3),transmat(3,3),SYM_POS_PRIM(3,NSYMMAX),RED_POS(3,NGUESS)
real(8):: tol,dist,rx1,ry1,rz1,rx2,ry2,rz2
logical:: found

NPOS_IRRED=0
NPOS_IRRED_ARR=0

call random_incell_main(LATSGP,CRYSSYS,NGUESS,RED_POS)

do i=1,NGUESS
POS(:)=RED_POS(:,i)
call create_unique_positions(RSYM,NSYMP,NSYM,BRAV,POS,SYM_POS_PRIM,NPOS,NSYMMAX)
write(16,*) NPOS

if(i==1) then
  NPOS_IRRED=NPOS_IRRED+1
  NPOS_IRRED_ARR(NPOS_IRRED)=NPOS
  goto 1020
else
  found=.false.
  do j=1,NPOS_IRRED
    if(NPOS.lt.NPOS_IRRED_ARR(j)) then
        rest=modulo(NPOS_IRRED_ARR(j),NPOS)
        if(rest==0) then
          found=.true.
          goto 1020
        endif
    elseif(NPOS.gt.NPOS_IRRED_ARR(j)) then
        rest=modulo(NPOS,NPOS_IRRED_ARR(j))
        if(rest==0) then
          found=.true.
          goto 1020
        endif
    else
        goto 1020
    endif
  enddo 
endif
  if(.not.found) then
    NPOS_IRRED=NPOS_IRRED+1
    NPOS_IRRED_ARR(NPOS_IRRED)=NPOS
    call sort(NPOS_IRRED_ARR(1:NPOS_IRRED),NPOS_IRRED)
  endif
1020 continue
enddo


end subroutine


subroutine sort(ARRAY,N)
implicit none
integer:: N, ARRAY(N),TMP,i,j
do j=1,N
do i=1,N-1
   if(ARRAY(i).gt.ARRAY(i+1)) then
     TMP=ARRAY(i)
     ARRAY(i)=ARRAY(i+1)
     ARRAY(i+1)=TMP
   endif 
enddo
enddo
end subroutine


!!************************************************************************************
!
! subroutine expand(rxyz,rxyzout,transvecall,latvec,nat)
! !This subroutine will expand the unit cell into 26 periodic cells and store them in rxyzout
! implicit none
! real*8, intent(in)  :: rxyz(3,nat),latvec(3,3)
! integer, intent(in) :: nat
! real*8, intent(out) :: rxyzout(3,nat,3,3,3) !26 periodic images plus the main cell
! integer             :: iat,iplane,icorner,iedge,m,k,l
! real*8,intent(inout):: transvecall(3,3,3,3)!,(transvecp(3,6),transvecc(3,8),transvece(3,12)
!
! do m=-1,1
!    do k=-1,1
!       do l=-1,1
!       transvecall(:,l+2,k+2,m+2)=real(l,8)*latvec(:,1)+real(k,8)*latvec(:,2)+real(m,8)*latvec(:,3)
!       enddo
!    enddo
! enddo
!
! do m=1,3
!    do k=1,3
!       do l=1,3
!       do iat=1,nat
!       rxyzout(:,iat,l,k,m)=rxyz(:,iat)+transvecall(:,l,k,m)
!       enddo
!       enddo
!    enddo
! enddo
! end
!
!!************************************************************************************
!
