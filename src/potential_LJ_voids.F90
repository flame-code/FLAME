module interface_lj_voids
  use global
  use defs_basis
  use void_lj_params
  !use cell_utils

  implicit none

  private
  public :: &
    lj_void_addon
 

contains

!!!!!!!!********************************************************************************
!!!!!!!subroutine lenosky_tb_lj(latvec,xred0,iprec,ka,kb,kc,fcart,energy,strten)
!!!!!!!!All H-Atoms have to be at the end of the rxyz-array
!!!!!!!implicit none
!!!!!!!integer:: i,n_in(3),offset_in(3),do_kpt_in,cellsearch,tb_or_meam,nec1,nec2,nec3,iprec,ka,kb,kc,iat
!!!!!!!real(8), parameter:: cut=5.24d0
!!!!!!!real(8):: xred(3,nat),rxyz(3,nat),xred0(3,nat),fcart(3,nat),latvec(3,3),latvec_ang(3,3),strten(6),stress_tb(3,3),energy,count
!!!!!!!real(8):: vol,fcart_tb(3,nat),fcart_lj(3,nat),strten_tb(6),strten_lj(6),energy_tb,energy_lj
!!!!!!!
!!!!!!!!Check if the atomic types are allright
!!!!!!!!call check_lenosky_tb_lj()
!!!!!!!
!!!!!!!!Number of kpoints
!!!!!!!n_in(1)=ka
!!!!!!!n_in(2)=kb
!!!!!!!n_in(3)=kc
!!!!!!!
!!!!!!!!Enable k-points
!!!!!!!if(ka==1.and.kb==1.and.kc==1) then 
!!!!!!!   do_kpt_in=0
!!!!!!!!Offset
!!!!!!!   offset_in(:)=0
!!!!!!!else
!!!!!!!   do_kpt_in=1
!!!!!!!!Offset
!!!!!!!   offset_in(:)=1
!!!!!!!endif
!!!!!!!
!!!!!!!!TB or MEAM
!!!!!!!tb_or_meam=0
!!!!!!!
!!!!!!!latvec_ang=latvec*Bohr_Ang
!!!!!!!xred=xred0
!!!!!!!call n_rep_dim(latvec_ang,cut,nec1,nec2,nec3)
!!!!!!!cellsearch=max(nec1,max(nec2,nec3))
!!!!!!!
!!!!!!!call backtocell(nat,latvec_ang,xred)
!!!!!!!call rxyz_int2cart(latvec_ang,xred,rxyz,nat)
!!!!!!!
!!!!!!!fcart_tb=0.d0
!!!!!!!call  lenoskytb(n_tb,rxyz(:,1:n_tb),fcart_tb(:,1:n_tb),energy_tb,count,n_silicon,latvec_ang(:,1),&
!!!!!!!      &latvec_ang(:,2),latvec_ang(:,3),stress_tb(:,1),stress_tb(:,2),stress_tb(:,3),&
!!!!!!!      &tb_or_meam,cellsearch,n_in(1),n_in(2),n_in(3),offset_in(1),offset_in(2),offset_in(3),do_kpt_in)
!!!!!!!!  stress=-stress
!!!!!!!  energy_tb=energy_tb/Ha_eV
!!!!!!!  fcart_tb=fcart_tb/Ha_eV*Bohr_Ang
!!!!!!!!Get full stress matrix
!!!!!!!call getvol(latvec_ang,vol)
!!!!!!!stress_tb=matmul(stress_tb,transpose(latvec_ang))/vol
!!!!!!!!This matrix MUST be symmetric!!!
!!!!!!!!!write(*,*) stress(:,1)           
!!!!!!!!!write(*,*) stress(:,2)           
!!!!!!!!!write(*,*) stress(:,3)           
!!!!!!!  strten_tb(1) = stress_tb(1,1)
!!!!!!!  strten_tb(2) = stress_tb(2,2)
!!!!!!!  strten_tb(3) = stress_tb(3,3)
!!!!!!!  strten_tb(6) = stress_tb(2,1)
!!!!!!!  strten_tb(5) = stress_tb(3,1)
!!!!!!!  strten_tb(4) = stress_tb(3,2)
!!!!!!!!This is not very clear yet...
!!!!!!!  strten_tb=strten_tb/Ha_eV*Bohr_Ang**3
!!!!!!!
!!!!!!!
!!!!!!!!Get the LJ stuff
!!!!!!!call lenosky_lj(latvec,xred0,fcart_lj,strten_lj,energy_lj)
!!!!!!!
!!!!!!!!Combine everything
!!!!!!!energy=energy_tb+energy_lj
!!!!!!!fcart=fcart_tb+fcart_lj
!!!!!!!strten=strten_tb+strten_lj
!!!!!!!
!!!!!!!end subroutine
!********************************************************************************
subroutine n_rep_dim(latvec,cut,nec1,nec2,nec3)
!This subroutine will return how many periodic expansions for each lattice vector direction are necessary for the periodic boundary conditions
!with for the given cut. nec1,nec2,nec3 for latvec(:,1),latvec(:,2),latvec(:,3)
implicit none
real*8 :: latvec(3,3),cut,nvec(3,3),point(3),point0(3),dist(3),eps,dd
integer:: i
integer:: nec1,nec2,nec3
! eps=1.d-6
nec1=0
nec2=0
nec3=0
call nveclatvec(latvec,nvec)
point0=(/0.d0,0.d0,0.d0/)
do i=1,3
call dist2plane(latvec(:,mod(i+1,3)+1),nvec(:,i),point0,dist(i))
! write(*,*) "cut",i,cut, dist
enddo
dist=abs(dist)
nec1=int(cut/dist(2))+1
nec2=int(cut/dist(3))+1
nec3=int(cut/dist(1))+1
end subroutine

!*********************************************************************************


        subroutine lj_void_addon(parini,latvec,xred0,fxyz,strten,etot)
!This routine will simply compute the lennard-jones interaction between the
!lj particles, and also with the other particles in the system.
!The lj particles are identified with znucl larger than 200,
!and they MUST be at the end of the list of atoms! Currenty,
!the first nat-nat_lj atoms are sent to the subroutines 
!that compute the energy/forces/etc from DFT or any other calculator, and
!the rest are computed here. In fact, DFT and the rest do never ever know
!that there are LJ perticles present in the system. We do that by
!temorarily reducing nat in the global module to nat-nat_lj

!The binary parameters are defined as follows in this scheme:

!There are NO interactions for the off diagonals and are ignored!!!

!        energy and forces for the runcated Lennard Jones potential according to PRA 8, 1504 (1973)
! input: nat: number of atoms totally
!        xred: reduced positions of atoms
!        latvec: lattice vectors of the periodic cell
!        pressure: pressure
!output: etot: energy
!        enth: enthalpy
!        strten: the stress tensor
!        fxyz: forces (negative derivative of energy with respect to positions
!
!        before calling this routine be sure to initialize the variables by calling the subroutine init_parameter(nat)

        use mod_parini, only: typ_parini
        implicit none
        type(typ_parini), intent(in):: parini
        integer:: iat, jat, kk, ll, l, nec(3), i, j, k, m
        real(8):: xred(3,parini%nat),fxyz(3,parini%nat),xred0(3,parini%nat),dxyz(3),r1red(3),r2red(3),rcut2(nat_lj)
        real(8):: etot,latvec_x(3,3),rec_nec(3),enth,stressvol(3,3),vol,strten(6)
        real(8), allocatable:: rel(:,:)
        real(8):: latvec(3,3),latvec_ang(3,3),latvecinv(3,3),celldv(3,3),pressure,gradvol(3,3),stress(3,3),tmplat(3,3)
        real(8):: trans(3,3), transinv(3,3)
        real(8):: crossp(3),d,dd,dd2,dd6,dd12,dx,dy,dz,s,s2,s6,s12,rc,rc2,rc6,rc12,tt,t1,t2,t3
        real(8):: sil_sjl1,sil_sjl2,sil_sjl3, si1_sj1, si2_sj2, si3_sj3, tkmsum1,tkmsum2,tkmsum3, cutmax,epscur
        real(8):: rcut2_lj(nat_lj),eps_lj_lj_fact,sigma_lj_lj_fact,double_count
        logical:: truncated
!Treat LJ-LJ interaction more strongly
eps_lj_lj_fact=10.d0
sigma_lj_lj_fact=1.5d0
!If truncated is true we will use the shifted form according to PRA 8, 1504 (1973), and else a simple crude cutoff
        truncated=.true.
!Convert everything from "internal atomic" units to the real "angstrom-ev" units
        latvec_ang=latvec*Bohr_Ang

        xred=xred0
        etot=0.d0
        fxyz(:,:)=0.d0
        celldv(:,:)=0.d0
        rcut2=rcutvoidlj*rcutvoidlj
        rcut2_lj=sigmavoidlj*2.d0**(1.d0/6.d0)*sigma_lj_lj_fact !This cutoff is due to the truncation at the minimum of the well, since only repulsion is used!!!
        cutmax=max(maxval(rcutvoidlj),maxval(rcut2_lj))
        rcut2_lj=rcut2_lj*rcut2_lj
        call backtocell(parini%nat,latvec_ang,xred) 
        trans=latvec_ang

         call invertmat(trans,transinv,3) 
         call n_rep_dim(latvec_ang,2.d0*cutmax,nec(1),nec(2),nec(3))


!Expand cell
    do i=1,3
      latvec_x(:,i)=real(nec(i),8)*latvec_ang(:,i)
    !Adjust reduced coordinates
      rec_nec(i)=1.d0/real(nec(i),8)
    enddo
!Interactions of LJ particles with all other particles
         do iat=nat_atoms+1,parini%nat !iat are the LJ particles
            r1red(:)=xred(:,iat)*rec_nec
            do i=0,nec(1)-1
            do j=0,nec(2)-1
            do k=0,nec(3)-1
               do jat=1,nat_atoms !jat are the "normal particles treated by DFT etc
                r2red(:)=xred(:,jat)*rec_nec
                r2red(1)=r2red(1)+real(i,8)*rec_nec(1)
                r2red(2)=r2red(2)+real(j,8)*rec_nec(2)
                r2red(3)=r2red(3)+real(k,8)*rec_nec(3)
                call pbc_distance0(latvec_x,r1red,r2red,dd,dxyz)
                !if(dd.gt.rcut2(iat-n_silicon-n_h)) then
                if(dd.gt.rcut2(iat-nat_atoms)) then
                   goto 1002
!                elseif(dd.lt.1.d-12) then
!                   goto 1002
                endif
                d=sqrt(dd)
                dx=-dxyz(1);  dy=-dxyz(2);  dz=-dxyz(3)
                dd2=1.d0/dd
                dd6=dd2*dd2*dd2
                dd12=dd6*dd6
!                s=sigmavoidlj(iat-n_silicon-n_h)   !Sigma
                s=sigmavoidlj(iat-nat_atoms)   !Sigma
                s2=s*s
                s6=s2*s2*s2
                s12=s6*s6
!                rc=1.d0/rcutvoidlj(iat-n_silicon-n_h) !Cutoff
                rc=1.d0/rcutvoidlj(iat-nat_atoms) !Cutoff
                rc2=rc*rc
                rc6=rc2*rc2*rc2
                rc12=rc6*rc6
!                epscur=epsvoidlj(iat-n_silicon-n_h)
                epscur=epsvoidlj(iat-nat_atoms)
                   etot=etot+4.d0*epscur*((s12*dd12-s6*dd6)) !Simple LJ
                   if(truncated) etot=etot+4.d0*epscur*((6.d0*s12*rc12-3.d0*s6*rc6)*dd*rc2-7.d0*s12*rc12+4.d0*s6*rc6)
                   tt=24.d0*epscur*dd2*(2.d0*s12*dd12-s6*dd6) !Simple LJ
                   if(truncated) tt=tt-8.d0*epscur*(6.d0*s12*rc12-3.d0*s6*rc6)*rc2
                   t1=dx*tt ; t2=dy*tt ; t3=dz*tt
                   fxyz(1,iat)=fxyz(1,iat)+t1
                   fxyz(2,iat)=fxyz(2,iat)+t2
                   fxyz(3,iat)=fxyz(3,iat)+t3
!Actio=Reactio
                   fxyz(1,jat)=fxyz(1,jat)-t1
                   fxyz(2,jat)=fxyz(2,jat)-t2
                   fxyz(3,jat)=fxyz(3,jat)-t3

                   !Cell gradient part
                   !My own implementation
                   si1_sj1=transinv(1,1)*dx+transinv(1,2)*dy+transinv(1,3)*dz
                   si2_sj2=transinv(2,1)*dx+transinv(2,2)*dy+transinv(2,3)*dz
                   si3_sj3=transinv(3,1)*dx+transinv(3,2)*dy+transinv(3,3)*dz
                   tkmsum1=trans(1,1)*si1_sj1+trans(1,2)*si2_sj2+trans(1,3)*si3_sj3
                   tkmsum2=trans(2,1)*si1_sj1+trans(2,2)*si2_sj2+trans(2,3)*si3_sj3
                   tkmsum3=trans(3,1)*si1_sj1+trans(3,2)*si2_sj2+trans(3,3)*si3_sj3
                   sil_sjl1=transinv(1,1)*dx+transinv(1,2)*dy+transinv(1,3)*dz
                   sil_sjl2=transinv(2,1)*dx+transinv(2,2)*dy+transinv(2,3)*dz
                   sil_sjl3=transinv(3,1)*dx+transinv(3,2)*dy+transinv(3,3)*dz

                   celldv(1,1)=celldv(1,1)+tt*tkmsum1*sil_sjl1  
                   celldv(1,2)=celldv(1,2)+tt*tkmsum1*sil_sjl2  
                   celldv(1,3)=celldv(1,3)+tt*tkmsum1*sil_sjl3  
                   celldv(2,1)=celldv(2,1)+tt*tkmsum2*sil_sjl1  
                   celldv(2,2)=celldv(2,2)+tt*tkmsum2*sil_sjl2  
                   celldv(2,3)=celldv(2,3)+tt*tkmsum2*sil_sjl3  
                   celldv(3,1)=celldv(3,1)+tt*tkmsum3*sil_sjl1  
                   celldv(3,2)=celldv(3,2)+tt*tkmsum3*sil_sjl2  
                   celldv(3,3)=celldv(3,3)+tt*tkmsum3*sil_sjl3  
                 1002 continue

        enddo
        enddo
        enddo
        enddo
        enddo

!Only the repulsive interactions between the LJ particles: Here we increse the interaction between the LJ particles by a factor of 10!!!
         do iat=nat_atoms+1,parini%nat !iat are th LJ particles 
            r1red(:)=xred(:,iat)*rec_nec
            do i=0,nec(1)-1
            do j=0,nec(2)-1
            do k=0,nec(3)-1
               do jat=nat_atoms+1,parini%nat !jat are the other LJ particles
                r2red(:)=xred(:,jat)*rec_nec
                r2red(1)=r2red(1)+real(i,8)*rec_nec(1)
                r2red(2)=r2red(2)+real(j,8)*rec_nec(2)
                r2red(3)=r2red(3)+real(k,8)*rec_nec(3)
                call pbc_distance0(latvec_x,r1red,r2red,dd,dxyz)
!                if(dd.gt.rcut2_lj(iat-n_silicon-n_h)) then
                if(dd.gt.rcut2_lj(iat-nat_atoms)) then
                   goto 1003
                elseif(dd.lt.1.d-12) then
                   goto 1003
                endif
                d=sqrt(dd)
                dx=-dxyz(1);  dy=-dxyz(2);  dz=-dxyz(3)
                dd2=1.d0/dd
                dd6=dd2*dd2*dd2
                dd12=dd6*dd6
!                s=sigmavoidlj(iat-n_silicon-n_h)*sigma_lj_lj_fact  !Sigma
                s=sigmavoidlj(iat-nat_atoms)*sigma_lj_lj_fact  !Sigma
                s2=s*s
                s6=s2*s2*s2
                s12=s6*s6
!!                rc=1.d0/rcut(iat-n_silicon-n_h) !Cutoff
!!                rc2=rc*rc
!!                rc6=rc2*rc2*rc2
!!                rc12=rc6*rc6
!                epscur=epsvoidlj(iat-n_silicon-n_h)*eps_lj_lj_fact
                epscur=epsvoidlj(iat-nat_atoms)*eps_lj_lj_fact
!                   if(abs(rcut2_lj(iat-n_silicon-n_h)-rcut2_lj(jat-n_silicon-n_h)).lt.1.d-15) then
                   if(abs(rcut2_lj(iat-nat_atoms)-rcut2_lj(jat-nat_atoms)).lt.1.d-15) then
                      double_count=0.5d0
                   else
                      double_count=1.d0
                   endif
!                   write(*,'(a,6es15.7)') "epslj",epscur,double_count,s12*dd12,s6*dd6,(4.d0*epscur*((s12*dd12-s6*dd6))+epscur)*double_count,d
                   etot=etot+(4.d0*epscur*((s12*dd12-s6*dd6))+epscur)*double_count !Repulsive LJ with shift
                   tt=24.d0*epscur*dd2*(2.d0*s12*dd12-s6*dd6)*double_count !Simple LJ
                   t1=dx*tt ; t2=dy*tt ; t3=dz*tt
                   fxyz(1,iat)=fxyz(1,iat)+t1
                   fxyz(2,iat)=fxyz(2,iat)+t2
                   fxyz(3,iat)=fxyz(3,iat)+t3
!Actio=Reactio
                   fxyz(1,jat)=fxyz(1,jat)-t1
                   fxyz(2,jat)=fxyz(2,jat)-t2
                   fxyz(3,jat)=fxyz(3,jat)-t3

                   !Cell gradient part
                   !My own implementation
                   si1_sj1=transinv(1,1)*dx+transinv(1,2)*dy+transinv(1,3)*dz
                   si2_sj2=transinv(2,1)*dx+transinv(2,2)*dy+transinv(2,3)*dz
                   si3_sj3=transinv(3,1)*dx+transinv(3,2)*dy+transinv(3,3)*dz
                   tkmsum1=trans(1,1)*si1_sj1+trans(1,2)*si2_sj2+trans(1,3)*si3_sj3
                   tkmsum2=trans(2,1)*si1_sj1+trans(2,2)*si2_sj2+trans(2,3)*si3_sj3
                   tkmsum3=trans(3,1)*si1_sj1+trans(3,2)*si2_sj2+trans(3,3)*si3_sj3
                   sil_sjl1=transinv(1,1)*dx+transinv(1,2)*dy+transinv(1,3)*dz
                   sil_sjl2=transinv(2,1)*dx+transinv(2,2)*dy+transinv(2,3)*dz
                   sil_sjl3=transinv(3,1)*dx+transinv(3,2)*dy+transinv(3,3)*dz

                   celldv(1,1)=celldv(1,1)+tt*tkmsum1*sil_sjl1  
                   celldv(1,2)=celldv(1,2)+tt*tkmsum1*sil_sjl2  
                   celldv(1,3)=celldv(1,3)+tt*tkmsum1*sil_sjl3  
                   celldv(2,1)=celldv(2,1)+tt*tkmsum2*sil_sjl1  
                   celldv(2,2)=celldv(2,2)+tt*tkmsum2*sil_sjl2  
                   celldv(2,3)=celldv(2,3)+tt*tkmsum2*sil_sjl3  
                   celldv(3,1)=celldv(3,1)+tt*tkmsum3*sil_sjl1  
                   celldv(3,2)=celldv(3,2)+tt*tkmsum3*sil_sjl2  
                   celldv(3,3)=celldv(3,3)+tt*tkmsum3*sil_sjl3  
                 1003 continue

        enddo
        enddo
        enddo
        enddo
        enddo

!Rescale
        etot=etot
        etot=etot/Ha_eV
        fxyz=fxyz/Ha_eV*Bohr_Ang
        celldv=celldv
        call getvol(latvec_ang,vol)
!!        vol=vol*real(nat,8)
!!        enth=etot+pressure*vol
!!        call stress_volume(latvec_ang,vol,pressure,stressvol)
!Transform the forces on the lattice vectors into the stress tensore according to the paper Tomas Bucko and Jurg Hafner
        do i=1,3
           tmplat(:,i)=latvec_ang(i,:)
        enddo
!The real stress tensor
        stress=-matmul(celldv,tmplat)/vol
        strten(1) = stress(1,1)
        strten(2) = stress(2,2)
        strten(3) = stress(3,3)
        strten(6) = stress(2,1)
        strten(5) = stress(3,1)
        strten(4) = stress(3,2)
!This is not very clear yet...
        strten=strten/Ha_eV*Bohr_Ang**3
!celldv has all forces to minimize the enthalpy directly on the lattice
!        celldv=celldv+stressvol
        return
        end subroutine

end module interface_lj_voids

!!!********************************************************
subroutine check_voids(parini)
use mod_parini, only: typ_parini
use global
use void_lj_params
use defs_basis
implicit none
type(typ_parini), intent(in):: parini
integer:: iat, ityp
logical:: in_atoms,in_lj
!Check typat for consistentcy and provide nat_lj
in_atoms=.false.
in_lj=.false.
nat_atoms=0
ntypat_atoms=0
ntypat_lj=0
nat_lj=0
do iat=1,parini%nat
   if(int(parini%znucl(parini%typat_global(iat))).lt.200) then
       nat_atoms=nat_atoms+1
       in_atoms=.true.
   endif
   if(int(parini%znucl(parini%typat_global(iat))).lt.200.and.in_lj) stop "Void system: First the physical atoms, then LJ pseudoparticles"
   if(int(parini%znucl(parini%typat_global(iat))).gt.200) then 
       nat_lj=nat_lj+1 
       in_lj=.true.
   endif
enddo
do ityp=1,parini%ntypat_global
 if(int(parini%znucl(ityp)).lt.200) then 
   ntypat_atoms=ntypat_atoms+1
 else
   ntypat_lj=ntypat_lj+1
 endif
enddo

end subroutine

!!!********************************************************

     subroutine voids_init_parameter()
     use global
     use void_lj_params
!The LJ particles interact with physical atoms, simply called atoms from here on.
!The sigma, epsilon and cutoff are given for each 
!LJ particle in a list lj_param.in
!There are NO interactions for the off diagonals and are ignored!!!
!In general the parmeters are symmetric, i.e. sigma(A,B)=sigma(B,A)
     integer:: iat
     logical:: file_exists
     character(40):: filename
     if(.not.allocated(sigmavoidlj)) allocate(sigmavoidlj(nat_lj)) 
     if(.not.allocated(epsvoidlj))   allocate(epsvoidlj(nat_lj)) 
     if(.not.allocated(rcutvoidlj))  allocate(rcutvoidlj(nat_lj)) 
      
     sigmalj=0.d0;epslj=0.d0
     file_exists=.false.
     filename="lj_param.in"
     INQUIRE(FILE=trim(filename), EXIST=file_exists)
     if(file_exists) then
         open(unit=33,file=trim(filename))
           do iat=1,nat_lj
           read(33,*) sigmavoidlj(iat),epsvoidlj(iat),rcutvoidlj(iat)
           enddo
         close(33)
      else
         stop "Please provide the parametrization of void's LJ interactions lj_param.in"
      endif
     end subroutine voids_init_parameter

!********************************************************
