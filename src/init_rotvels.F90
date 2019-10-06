subroutine init_rotvels(parini,nat,xred,latvec,temp,amass,vel)
use mod_parini, only: typ_parini
!This routine will first find the correct partitioning of the system into molecules, then assign
!rotational and translational velocities to these molecules according to the tempereature temp
use global, only: units
use defs_basis, only: Bohr_Ang,pi,kb_HaK
implicit none
type(typ_parini), intent(in):: parini
integer,intent(in):: nat
real(8),intent(in):: xred(3,nat),latvec(3,3),temp,amass(nat)
real(8),intent(out):: vel(3,nat)
integer:: iat,nfrag,ifrag,LWORK,INFO,idim,omega_opt
real(8):: xcart(3,nat),ekin_rot,ekin_trans,rotmat(3,3),dproj(6),angbohr,erot_tmp,ekin_tot,v2gauss,vtest,tmp(3)
real(8):: etrans_tmp,latvec_tmp(3,3),diag_inert(3,3),dir(3),weight(3),frac_rot(2)
real(8),allocatable:: cmass(:,:),vel_cmass(:,:),intens(:,:,:),fragxcart(:,:,:),fragmass(:,:),fragvel(:,:,:),masstot(:)
real(8),allocatable:: WORK(:),eval(:),frag_ekin_rot(:),omega(:,:)
integer,allocatable:: counter(:)
integer, dimension(nat):: fragarr,fragsize
!Exact total kinetic energy at given temperature
ekin_tot=3.d0*real(nat,8)*kb_HaK*temp

!Direction to bias the rotational component:
!1 Random
!2 Along lowest intertia tensor direction
!3 Weighted according to the principle inertia
!4 A weighted random distribution of the rotational direction: principle axis weighted with 1, 1/2, 1/4
!5 A weighted random distribution of the rotational direction: principle axis weighted with 1, 1/3, 1/6
!6 A weighted random distribution of the rotational direction: principle axis weighted with 1, 1/4, 1/8
!7 A weighted random distribution of the rotational direction: principle axis weighted with 1, 1/2, 1/4
!8 A weighted random distribution of the rotational direction: principle axis weighted with 1, 1/3, 1/9
!9 A weighted random distribution of the rotational direction: principle axis weighted with 1, 1/4, 1/16
omega_opt=3

weight=1.d0
SELECT CASE (omega_opt)
   CASE (2)
      weight=(/1.d0, 0.d0, 0.d0/)
   CASE (3)
      weight=(/1.d0, 1.d0, 1.d0/)
   CASE (4)
      weight=(/1.d0, 1.d0/2.d0, 1/4.d0/)
   CASE (5)
      weight=(/1.d0, 1.d0/3.d0, 1/6.d0/)
   CASE (6)
      weight=(/1.d0, 1.d0/4.d0, 1/8.d0/)
   CASE (7)
      weight=(/1.d0, 1.d0/2.d0, 1/4.d0/)
   CASE (8)
      weight=(/1.d0, 1.d0/3.d0, 1/9.d0/)
   CASE (9)
      weight=(/1.d0, 1.d0/4.d0, 1/16.d0/)
  CASE DEFAULT
      weight=1.d0
END SELECT


!Partition the kinetic energy to go into the translational and the rotational energy
frac_rot(1)=0.3d0;frac_rot(2)=0.7d0 !Lower and upper fractional boundary for rotational energy
call random_number(ekin_rot)
ekin_rot=ekin_rot*(frac_rot(2)-frac_rot(1))+frac_rot(1)
ekin_trans=(1.d0-ekin_rot)*ekin_tot
ekin_rot=ekin_rot*ekin_tot
vel=0.d0

!First identify the fragments or molecules in the cell
call fragments(parini,latvec,xred,nfrag,xcart,fragarr,fragsize)

!Allocate the arrays corresponding to molecular quantities
allocate(cmass(3,nfrag),vel_cmass(3,nfrag),intens(3,3,nfrag))
allocate(fragxcart(3,maxval(fragsize),nfrag),fragmass(maxval(fragsize),nfrag),counter(nfrag))
allocate(frag_ekin_rot(nfrag),omega(3,nfrag),fragvel(3,maxval(fragsize),nfrag),masstot(nfrag))

!Setup the xcart array for all nfrag clusters
counter=0
do iat=1,nat
   counter(fragarr(iat))=counter(fragarr(iat))+1
   fragxcart(:,counter(fragarr(iat)),fragarr(iat))=xcart(:,iat)
   fragmass(counter(fragarr(iat)),fragarr(iat))=amass(iat)
enddo

!Write the connected molecules
open(unit=2,file="xcart.molecule.ascii")
angbohr=1.d0
if(trim(units)=="angstroem") then
  angbohr=Bohr_Ang
endif
write(2,*) nat
latvec_tmp=latvec
call latvec2dproj(dproj,latvec_tmp,rotmat,xcart,nat)
write(2,*) dproj(1:3)*angbohr
write(2,*) dproj(4:6)*angbohr
  do iat=1,nat
       write(2,'(3(1x,es25.15),2x,a2)') angbohr*xcart(:,iat),trim(parini%char_type(parini%typat_global(iat)))
  enddo
close(2)

!Compute the inertia tensor stuff
cmass=0.d0
!Now we will compute for each of the fragments the inertia tensor
do ifrag=1,nfrag
  masstot(ifrag)=0.d0
  !Compute the center of mass of every fragment
  do iat=1,fragsize(ifrag)
      masstot(ifrag)=masstot(ifrag)+fragmass(iat,ifrag)
      cmass(:,ifrag)=cmass(:,ifrag)+fragmass(iat,ifrag)*fragxcart(:,iat,ifrag)
  enddo
  cmass(:,ifrag)=cmass(:,ifrag)/masstot(ifrag)
  !Compute intertia tensor
  call inertia_tensor(fragsize(ifrag),fragxcart(:,1:fragsize(ifrag),ifrag),cmass(:,ifrag),&
       &fragmass(1:fragsize(ifrag),ifrag),intens(:,:,ifrag))
!  write(*,*) "ifrag,mass",ifrag,masstot(ifrag)
!  write(*,*) "CM", cmass(:,ifrag)
!  write(*,*) intens(:,1,ifrag)
!  write(*,*) intens(:,2,ifrag)
!  write(*,*) intens(:,3,ifrag)
!LWORK=-1
!allocate(WORK(1),eval(3))
!        call DSYEV('V','L',3,intens(:,:,ifrag),3,eval,WORK,LWORK,INFO)
!        if (info.ne.0) stop 'DSYEV in correct_hessin'
!LWORK=WORK(1)
!deallocate(WORK)
!allocate(WORK(LWORK))
!        call DSYEV('V','L',3,intens(:,:,ifrag),3,eval,WORK,LWORK,INFO)
!        if (info.ne.0) stop 'DSYEV in correct_hessin'
!        write(*,*) '---   App. eigenvalues in a.u. -------------'
!        do iat=1,3
!         write(*,'(1x,es25.15)') eval(iat)
!        enddo
!deallocate(work,eval)
end do

!Assign the target kinetic rotational energies of each fragment
!And aso the rotational exis, not scaled yet omega
fragvel=0.d0
call random_number(frag_ekin_rot)
!write(*,*) "RANDOM",frag_ekin_rot
!Scale to 1
frag_ekin_rot=frag_ekin_rot/sum(frag_ekin_rot)*ekin_rot
do ifrag=1,nfrag
!Only run if molecule is larger than 1 atom
if(fragsize(ifrag).gt.2) then
!Set omega
   if(omega_opt==1) then
!Either random direction
   call rand_sphere(omega(:,ifrag),pi)
!... or along the smallest inertia tensor  
   elseif(omega_opt.le.9) then
         !Diagonalize the tensor
          diag_inert(:,:)=intens(:,:,ifrag)
          allocate(WORK(1),eval(3))
          LWORK=-1
          call DSYEV('V','L',3,diag_inert,3,eval,WORK,LWORK,INFO)
          if (info.ne.0) stop 'DSYEV in init_rotvels'
          LWORK=WORK(1)
          deallocate(WORK)
          allocate(WORK(LWORK))
          call DSYEV('V','L',3,diag_inert,3,eval,WORK,LWORK,INFO)
          if (info.ne.0) stop 'DSYEV in correct_hessin'
!Determine direction of rotation, meaning the sign of the vector, and the weight
              write(*,'(a,i5,3(es25.15))') "Weight e  : ",ifrag,eval
              write(*,'(a,i5,3(es25.15))') "Weight w  : ",ifrag,weight
              call random_number(dir)
              dir=(dir-0.5d0)*2.d0
              write(*,'(a,i5,3(es25.15))') "Weight r  : ",ifrag,dir
              dir=dir*weight
              write(*,'(a,i5,3(es25.15))') "Weight rw : ",ifrag,dir
              dir(2)=dir(2)*eval(1)/eval(2);dir(3)=dir(3)*eval(1)/eval(3)
              write(*,'(a,i5,3(es25.15))') "Weight rwe: ",ifrag,dir
              omega(:,ifrag)=diag_inert(:,1)*dir(1)+diag_inert(:,2)*dir(2)+diag_inert(:,3)*dir(3)
          deallocate(work,eval)
          else
              stop "Wrong option for omega_opt"
    endif
       
       
          
!Set the correct rotational velocity omega
   call rot_ener(omega(:,ifrag),intens(:,:,ifrag),erot_tmp)
!   write(*,*) "Init",frag_ekin_rot(ifrag),erot_tmp
   omega(:,ifrag)=omega(:,ifrag)*sqrt(frag_ekin_rot(ifrag)/erot_tmp)
   call rot_ener(omega(:,ifrag),intens(:,:,ifrag),erot_tmp)
!   write(*,*) "Scaled",frag_ekin_rot(ifrag),erot_tmp
!Assign rotational velocities to the atoms in the fragments
   call assign_vel(fragsize(ifrag),fragxcart(:,1:fragsize(ifrag),ifrag),&
        &cmass(:,ifrag),omega(:,ifrag),fragvel(:,1:fragsize(ifrag),ifrag))
!Recompute v2gauss
         v2gauss=0.d0
         vtest=0.d0
         do iat=1,fragsize(ifrag)
           do idim=1,3
             v2gauss=v2gauss+fragvel(idim,iat,ifrag)**2*amass(iat)
             vtest=vtest+fragvel(idim,iat,ifrag)/(3.d0*real(nat,8))
           end do
         end do
!   write(*,*) "TEST", 0.5d0*v2gauss,vtest
endif
enddo

!Assign translational kinetic energy
!Get random Gaussian distributed atomic velocities
call gausdist(nfrag,vel_cmass,masstot)
call elim_moment(nfrag,vel_cmass,masstot)
!Set the correct translational velocity
etrans_tmp=0.d0
do ifrag=1,nfrag
  etrans_tmp=etrans_tmp+0.5d0*dot_product(vel_cmass(:,ifrag),vel_cmass(:,ifrag))*masstot(ifrag)
enddo
!   write(*,*) "Init",ekin_trans,etrans_tmp
   vel_cmass(:,:)=vel_cmass(:,:)*sqrt(ekin_trans/etrans_tmp)
etrans_tmp=0.d0
do ifrag=1,nfrag
  etrans_tmp=etrans_tmp+0.5d0*dot_product(vel_cmass(:,ifrag),vel_cmass(:,ifrag))*masstot(ifrag)
enddo
!   write(*,*) "Scaled",ekin_trans,etrans_tmp
   !Recompute v2gauss
         v2gauss=0.d0
         vtest=0.d0
         do iat=1,nfrag
           do idim=1,3
             v2gauss=v2gauss+vel_cmass(idim,iat)**2*masstot(iat)
             vtest=vtest+vel_cmass(idim,iat)/(3.d0*real(nfrag,8))
           end do
         end do
!   write(*,*) "TEST", 0.5d0*v2gauss,vtest

!Distribute the CM velocities to each cluster
do ifrag=1,nfrag
   do iat=1,fragsize(ifrag)
   fragvel(:,iat,ifrag)=fragvel(:,iat,ifrag)+vel_cmass(:,ifrag)
   enddo
end do

!Transform back to the current velocity structure
counter=0
do iat=1,nat
   counter(fragarr(iat))=counter(fragarr(iat))+1
   vel(:,iat)=fragvel(:,counter(fragarr(iat)),fragarr(iat))
enddo

!Check total energy
etrans_tmp=0.d0
do iat=1,nat
   etrans_tmp=etrans_tmp+0.5d0*dot_product(vel(:,iat),vel(:,iat))*amass(iat)
enddo

!write(*,*) "TOTAL", etrans_tmp,ekin_tot

deallocate(cmass,vel_cmass,intens,fragxcart,fragmass,counter)
deallocate(frag_ekin_rot,omega,fragvel,masstot)
end subroutine

!************************************************************************************

subroutine assign_vel(nat,xcart,cmass,omega,vel)
!This routine will take the cartesian coordinates and the center of mass and compute the 
!velocities of every atom in the system with nat atoms
implicit none
integer:: nat,iat
real(8):: xcart(3,nat),cmass(3),xtmp(3),vel(3,nat),omega(3)
do iat=1,nat
 xtmp=xcart(:,iat)-cmass(:)
 call cross_product(omega,xtmp,vel(:,iat))
enddo
end subroutine

!************************************************************************************

subroutine rot_ener(omega,intens,erot)
!This routine computes the rotational kinetic energy of a system with known tensor of intertia and a given 
!angular velocity omega, based on its orientation and magnitude
implicit none
real(8):: omega(3),intens(3,3),erot
erot= 0.5d0*dot_product(omega,matmul(intens,omega))
end subroutine

