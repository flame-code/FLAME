subroutine fragments(parini,latvec,xred,nfrag,xcart,fragarr,fragsize)
use mod_parini, only: typ_parini
implicit none
type(typ_parini), intent(in):: parini
real(8),dimension(3,parini%nat), INTENT(IN) :: xred
real(8):: latvec(3,3),rotmat(3,3),dproj(6)
real(8),allocatable:: pos(:,:)
integer :: nfrag, nfragold
logical :: occured,niter
real(8)::  dist, mindist, angle, vec(3), cmass(3), velcm(3), bondlength, bfactor,rnrmi,scpr
real(8):: ekin,vcm1,vcm2,vcm3,ekin0,scale,xcart(3,parini%nat)
integer::iat, jat, nmax(1), imin(2),ifrag
integer, dimension(parini%nat):: fragarr,fragsize

bfactor=1.1d0
fragarr(:)=0                  !Array, which atom belongs to which fragment
nfrag=0                       !Number of fragments
fragsize=0                    !Size of each fragment

!Calculate number of fragments and fragmentlist of the atoms
do
  nfragold=nfrag
  do iat=1,parini%nat                !Check the first atom that isn't part of a cluster yet
    if(fragarr(iat)==0) then
      nfrag=nfrag+1
      fragarr(iat)=nfrag
      fragsize(nfrag)=fragsize(nfrag)+1
      xcart(:,iat)=matmul(latvec,xred(:,iat))
      exit
    endif
  enddo
  if (nfragold==nfrag) exit
7000 niter=.false.
  do iat=1,parini%nat                !Check if all the other atoms are part of the current cluster
    do jat=1,parini%nat
    bondlength=parini%rcov(parini%typat_global(iat))+parini%rcov(parini%typat_global(jat))
    if(nfrag==fragarr(iat) .AND. jat.ne.iat .AND. fragarr(jat)==0) then
!         call pbc_distance1(latvec,xred(:,iat),xred(:,jat),dist)
         call pbc_distance2(latvec,xred(:,iat),xcart(:,iat),xred(:,jat),xcart(:,jat),dist)
         if(dist<(bfactor*bondlength)**2) then
         fragarr(jat)=nfrag
         fragsize(nfrag)=fragsize(nfrag)+1
         niter=.true.
         endif
    endif
    enddo
  enddo
  if(niter) then
  goto 7000
  endif
enddo


open(unit=46,file="pos_fragment.ascii")
write(46,*) "Fragmentation in cartesian coordinates"
allocate(pos(3,parini%nat))
pos=xcart
call latvec2dproj(dproj,latvec,rotmat,pos,parini%nat)
write(46,*) dproj(1:3)
write(46,*) dproj(4:6)
do iat=1,parini%nat
  write(46,'(3(1x,es25.15),2x,a2,1x,a1,i4)')       pos(:,iat),trim(parini%char_type(parini%typat_global(iat))),"#",fragarr(iat)
enddo
deallocate(pos)
close(46)
end subroutine fragments

!************************************************************************************

subroutine get_fragsize(fragsize,lhead,llist,nat,nmol)
implicit none
integer:: nat,nmol,iat,ifrag,llist(nat),lhead(nmol),fragsize(nmol)
fragsize=0
do ifrag=1,nmol
   iat=lhead(ifrag)
   do while (iat.ne.0)
     fragsize(ifrag)=fragsize(ifrag)+1
     iat=llist(iat)
   enddo
enddo
end subroutine

!************************************************************************************

subroutine refragment(fragarr,nat)
!This routine will rearragne the integer array fragarr such that the fragment indexes are well 
!assigned in ascending order, and new array indexes are assigned to atoms with 
!values <0
implicit none
integer:: fragarr(nat),nat,iat,jat,cnt,find,fragarr_tmp(nat)
cnt=0
fragarr_tmp=0
if(all(fragarr.le.0)) then
 do iat=1,nat
    fragarr(iat)=iat
 enddo
else
!first assign new values to existing fragments
do iat=1,nat
   if(fragarr(iat).gt.0.and.fragarr_tmp(iat)==0) then
      cnt=cnt+1
      do jat=1,nat 
        if(fragarr(iat)==fragarr(jat)) fragarr_tmp(jat)=cnt
      enddo
   endif
enddo
!Now all unassigned clusters
do iat=1,nat
  if(fragarr(iat).le.0) then
    cnt=cnt+1
    fragarr_tmp(iat)=cnt
  endif 
enddo
fragarr=fragarr_tmp
endif
end subroutine

!************************************************************************************

subroutine make_linked_list(fragarr,fragsize,lhead,llist,nat,nmol)
!This subroutine will create a linked list to represent molecules.
!lhead is an array of length nat, of which only nmol will be significant. Their entries point to a
!position in array lhead, the first atom in the molecule, 
!llist_at is an array with entries linking the atoms in the molecule. Termination flag is 0
implicit none
integer:: fragarr(nat),nat,iat,nmol,ifrag
integer:: lhead(nmol),llist(nat),fragsize(nmol)
llist=0
lhead=0
do iat=1,nat
   llist(iat)=lhead(fragarr(iat)) 
   lhead(fragarr(iat))=iat
enddo
call get_fragsize(fragsize,lhead,llist,nat,nmol)
!do ifrag=1,nmol
!   iat=lhead(ifrag)
!   do while (iat.ne.0)
!     write(*,*) ifrag, iat
!     iat=llist(iat)
!   enddo
!enddo
end subroutine

!************************************************************************************

subroutine get_cmass(cmass,masstot,xcart,amass,lhead,llist,nat,nmol)
!This subroutine will compute the center of mass and the total mass of the fragment
!given the cartesian coordinates xcart and the linked cell lists and atomic masses
implicit none
integer:: nat,nmol,iat,ifrag,lhead(nmol),llist(nat)
real(8):: xcart(3,nat),amass(nat),masstot(nmol),cmass(3,nmol)
cmass=0.d0
do ifrag=1,nmol
  masstot(ifrag)=0.d0
  !Compute the center of mass of every fragment
   iat=lhead(ifrag)
   do while (iat.ne.0)
        masstot(ifrag)=masstot(ifrag)+amass(iat)
        cmass(:,ifrag)=cmass(:,ifrag)+amass(iat)*xcart(:,iat)
        iat=llist(iat)
    enddo
  cmass(:,ifrag)=cmass(:,ifrag)/masstot(ifrag)
enddo
end subroutine

!************************************************************************************

subroutine get_inertia_tensor(parini,intens,inprin,inaxis,cmass,xcart,amass,lhead,llist,nat,nmol)
use mod_parini, only: typ_parini
!This routine will compute the intertia tensors of all molecules involved
!given the center of mass and atomic masses, the principle  moments of inertia inprin,
!and the axis of the inertia tensor for each molecule
implicit none
type(typ_parini), intent(in):: parini
integer:: nat,nmol,iat,ifrag,i,j,llist(nat),lhead(nmol),LWORK,info
real(8):: xcart(3,nat),amass(nat),cmass(3,nmol),intens(3,3,nmol),dist2,xtmp(3)
real(8):: inprin(3,nmol),inaxis(3,3,nmol),diag_inert(3,3),tmp_vec(3),tmp_val
real(8),allocatable:: work(:),eval(:)
real(8):: delta_kronecker
!Compute the tensor of intertia
do ifrag=1,nmol
   iat=lhead(ifrag)
   intens(:,:,ifrag)=0.d0
   do while (iat.ne.0)
       xtmp=xcart(:,iat)-cmass(:,ifrag)
       dist2=dot_product(xtmp,xtmp)
         do i=1,3
         do j=1,i
            intens(i,j,ifrag)=intens(i,j,ifrag)+amass(iat)*(dist2*delta_kronecker(i,j)-xtmp(i)*xtmp(j))
         enddo
         enddo
       iat=llist(iat)
   enddo
   intens(1,2,ifrag)=intens(2,1,ifrag)
   intens(1,3,ifrag)=intens(3,1,ifrag)
   intens(2,3,ifrag)=intens(3,2,ifrag)
enddo

!Compute the priciple axis and the principle moments of inertia
!Get the correct array size
diag_inert(:,:)=0.d0
allocate(WORK(1),eval(3))
LWORK=-1
call DSYEV('V','L',3,diag_inert,3,eval,WORK,LWORK,INFO)
if (info.ne.0) stop 'DSYEV in get_inertia_tensor'
LWORK=WORK(1)
deallocate(WORK)
allocate(WORK(LWORK))

do ifrag=1,nmol
if(parini%fragsize(ifrag)==1) then !We have a single atom
  inprin(:,ifrag)=0.d0
  inaxis(:,1,ifrag)=(/1.d0,0.d0,0.d0/)
  inaxis(:,2,ifrag)=(/0.d0,1.d0,0.d0/)
  inaxis(:,3,ifrag)=(/0.d0,0.d0,1.d0/)
else
!Diagonalize the tensor
          diag_inert(:,:)=intens(:,:,ifrag)
          call DSYEV('V','L',3,diag_inert,3,eval,WORK,LWORK,INFO)
          if (info.ne.0) stop 'DSYEV in get_inertia_tensor'
          if(eval(1).lt.-1.d-10) stop "Negative diagonal element in inertia tensor!!!"
          inaxis(:,1,ifrag)=sign(1.d0,eval(1))*diag_inert(:,1)
          inaxis(:,2,ifrag)=sign(1.d0,eval(2))*diag_inert(:,2)
          inaxis(:,3,ifrag)=sign(1.d0,eval(3))*diag_inert(:,3)
          inprin(1,ifrag)=abs(eval(1))
          inprin(2,ifrag)=abs(eval(2))
          inprin(3,ifrag)=abs(eval(3))
endif
enddo
deallocate(work,eval)
end subroutine

!************************************************************************************

subroutine get_fcm_torque(fcm,torque,fcart,quat,xcart_mol,lhead,llist,nat,nmol)
!Computes the total force on molecules and the torques,assuming xcart_mol has molecules with CM at origin
implicit none
integer:: nat,nmol,iat,ifrag,llist(nat),lhead(nmol)
real(8),intent(in):: fcart(3,nat),quat(4,nmol),xcart_mol(3,nat)
real(8):: fcm(3,nmol),torque(3,nmol),crossp(3),xtmp(3),rotmat(3,3)
do ifrag=1,nmol
   call quat2rotmat(rotmat,quat(:,ifrag))
   iat=lhead(ifrag)
   fcm(:,ifrag)=0.d0
   torque(:,ifrag)=0.d0
   do while (iat.ne.0)
     fcm(:,ifrag)=fcm(:,ifrag)+fcart(:,iat)   
     xtmp=matmul(rotmat,xcart_mol(:,iat))
     call cross_product(xtmp,fcart(:,iat),crossp)
     torque(:,ifrag)=torque(:,ifrag)+crossp
     iat=llist(iat)
   enddo
enddo
end subroutine

!************************************************************************************

subroutine init_cm_mol(parini,latvec,xred,xcart_shifted,xred_cm,quat,amass,masstot,intens,inprin,inaxis,lhead,llist,nat,nmol)
use mod_parini, only: typ_parini
!This routine will get the cm and shift all molecular units into xcart_shifted
!and write, for each molecule, an xyz file containing the shifted molecular unit
!We will also express the center of masses in reduced coordinates (xred_cm), and of
!course intiallize the orientation vectors quat to 0
use defs_basis, only: Ha_eV,Bohr_Ang
use global, only: units
implicit none
type(typ_parini), intent(in):: parini
real(8),intent(in):: latvec(3,3),xred(3,nat)
integer:: nat,nmol,iat,llist(nat),lhead(nmol),fragsize(nmol),imol,jmol,kmol
real(8):: xcart_in(3,nat),xcart_shifted(3,nat),cmass(3,nmol),amass(nat)
real(8):: masstot(nmol),angbohr,quat(4,nmol),xred_cm(3,nmol),xcart_tmp(3,nat)
real(8):: circular(3,3),tol,rot_c(3,3),rot_all(3,3),inprin(3,nmol),intens(3,3,nmol)
real(8):: inaxis(3,3,nmol),ident(3,3),tmp(3,nmol),quat_tmp(4),tmp_real(4),tmp_mat(3,3)
character(5):: fn
logical:: symtop(nmol),tmp_logical
tol=1.d-10
!Convert to cartesian coordinates
call rxyz_int2cart(latvec,xred,xcart_in,nat)
!Get the center of masses and save it in xred_cm
call get_cmass(cmass,masstot,xcart_in,amass,lhead,llist,nat,nmol)
call rxyz_cart2int(latvec,xred_cm,cmass,nmol)
!Compute the fragment sizes
call get_fragsize(fragsize,lhead,llist,nat,nmol)
!Get inertia tensor
call get_inertia_tensor(parini,intens,inprin,inaxis,cmass,xcart_in,amass,lhead,llist,nat,nmol)
!Circular matrix that will rotate xyz to zxy
circular(1,:)=(/0.d0,0.d0,1.d0/)
circular(2,:)=(/1.d0,0.d0,0.d0/)
circular(3,:)=(/0.d0,1.d0,0.d0/)
!Identity matrix
ident=0.d0
ident(1,1)=1.d0
ident(2,2)=1.d0
ident(3,3)=1.d0
!Initiallize the orientation vector
quat=0.d0
!Compute symmetric tops and rotate accordingly
symtop=.false.
rot_c=ident
   write(*,*) inprin
do imol=1,nmol
 if(fragsize(imol).ge.2) then
   if(abs(inprin(3,imol)-inprin(1,imol)).lt.tol) then
     rot_c=circular
     symtop(imol)=.true.
   elseif(abs(inprin(2,imol)-inprin(3,imol)).lt.tol)then  
     rot_c=matmul(circular,circular)
     symtop(imol)=.true.
   elseif(abs(inprin(1,imol)-inprin(2,imol)).lt.tol)then
     symtop(imol)=.true.
   endif
   rot_all=matmul(rot_c,transpose(inaxis(:,:,imol)))
   call rotmat2quat(rot_all,quat_tmp) 
   call qinv(quat_tmp,quat(:,imol))
  !Rotate it all
  !do imol=1,nmol
   inaxis(:,:,imol)=matmul(rot_c,inaxis(:,:,imol))
   inaxis(:,:,imol)=matmul(inaxis(:,:,imol),transpose(inaxis(:,:,imol)))
   inprin(:,imol)=matmul(rot_c,inprin(:,imol))
   if(symtop(imol)) then !Average the principle axis
       inprin(1,imol)=0.5d0*(inprin(1,imol)+inprin(2,imol))
       inprin(2,imol)=inprin(1,imol)
   endif
 else
   symtop(imol)=.true.
 endif
   iat=lhead(imol)
   do while (iat.ne.0)
     xcart_shifted(:,iat)=xcart_in(:,iat)-cmass(:,imol)
     xcart_shifted(:,iat)=matmul(rot_all,xcart_shifted(:,iat)) 
     iat=llist(iat)
   enddo
!enddo
enddo
!Rearrange the indexes of lhead, such that the first molecules are arbitrary, then such with symmetric tops, then such with 2 atoms, then with 1 atom
do imol=1,nmol
do jmol=2,nmol
 if(((fragsize(jmol).gt.fragsize(jmol-1)).and.(symtop(jmol-1).and..not.symtop(jmol))).or.&      !To move the smallest molecules down the list
   &((fragsize(jmol).gt.fragsize(jmol-1)).and.(symtop(jmol-1).and.symtop(jmol)))) then !For the case that we have other symmetric molecules
   iat=              lhead(jmol-1);      lhead(jmol-1)=     lhead(jmol);      lhead(jmol)=iat
   tmp_real=        quat(:,jmol-1);     quat(:,jmol-1)=    quat(:,jmol);     quat(:,jmol)=tmp_real
   iat=           fragsize(jmol-1);   fragsize(jmol-1)=  fragsize(jmol);   fragsize(jmol)=iat
   tmp_logical=     symtop(jmol-1);     symtop(jmol-1)=    symtop(jmol);     symtop(jmol)=tmp_logical
   tmp_real(1:3)=  cmass(:,jmol-1);    cmass(:,jmol-1)=   cmass(:,jmol);    cmass(:,jmol)=tmp_real(1:3)
   tmp_mat=intens     (:,:,jmol-1); intens(:,:,jmol-1)=intens(:,:,jmol); intens(:,:,jmol)=tmp_mat
   tmp_mat=inaxis     (:,:,jmol-1); inaxis(:,:,jmol-1)=inaxis(:,:,jmol); inaxis(:,:,jmol)=tmp_mat
   tmp_real(1:3)= inprin(:,jmol-1);   inprin(:,jmol-1)=  inprin(:,jmol);   inprin(:,jmol)=tmp_real(1:3)  
   tmp_real(1)=    masstot(jmol-1);    masstot(jmol-1)=   masstot(jmol);    masstot(jmol)=tmp_real(1)
   tmp_real(1:3)=xred_cm(:,jmol-1);  xred_cm(:,jmol-1)= xred_cm(:,jmol);  xred_cm(:,jmol)=tmp_real(1:3)
 endif  
enddo
enddo
write(*,*) symtop
!Get inertia tensor
xcart_tmp=0.d0
!call get_inertia_tensor(intens,inprin,inaxis,xcart_tmp,xcart_shifted,amass,lhead,llist,nat,nmol)
       do iat=1,nmol
          write(122,*) "Fragsize",fragsize(iat)
          write(122,*) "INTENS",iat
          write(122,*) intens(:,1,iat)
          write(122,*) intens(:,2,iat)
          write(122,*) intens(:,3,iat)
          write(122,*) "INPRIN"
          write(122,*) inprin(:,iat)
          write(122,*) "INAXIS",iat
          write(122,*) inaxis(:,1,iat)
          write(122,*) inaxis(:,2,iat)
          write(122,*) inaxis(:,3,iat)
       enddo

if(trim(units)=="angstroem") then
  angbohr=Bohr_Ang
elseif(trim(units)=="bohr") then
  angbohr=1.d0
else
  stop "Wrong file unit format"
endif

!!do imol=1,nmol
!!iat=lhead(imol)
!!   call quat2rotmat(rot_all,quat(:,imol)) 
!!   do while (iat.ne.0)
!!     write(222,'(a2,2x,3(es25.15))') trim(char_type(typat(iat))),(cmass(:,imol)+matmul(rot_all,xcart_shifted(:,iat)))*angbohr
!!     iat=llist(iat)
!!   enddo
!!enddo

do imol=1,nmol
   write(fn,'(i5.5)') imol
   open(unit=22,file="posfrag"//fn//".xyz")
   write(22,*)  fragsize(imol), trim(units)
   write(22,*)  
   iat=lhead(imol)
   do while (iat.ne.0)
     write(22,'(a2,2x,3(es25.15))') trim(parini%char_type(parini%typat_global(iat))),xcart_shifted(:,iat)*angbohr
     iat=llist(iat)
   enddo
   close(22)
enddo
end subroutine
