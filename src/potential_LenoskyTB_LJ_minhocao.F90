module interface_lenosky_tb_lj
  use global
  use defs_basis
  use tb_lj_params
  !use cell_utils

  implicit none

  private
  public :: &
    lenosky_tb_lj
 

contains

!********************************************************************************
subroutine lenosky_tb_lj(parini,latvec,xred0,iprec,ka,kb,kc,fcart,energy,strten)
!All H-Atoms have to be at the end of the rxyz-array
use mod_parini, only: typ_parini
implicit none
type(typ_parini), intent(in):: parini
integer:: i,n_in(3),offset_in(3),do_kpt_in,cellsearch,tb_or_meam,nec1,nec2,nec3,iprec,ka,kb,kc,iat
real(8), parameter:: cut=5.24d0
real(8):: xred(3,parini%nat),rxyz(3,parini%nat),xred0(3,parini%nat),fcart(3,parini%nat),latvec(3,3),latvec_ang(3,3),strten(6),stress_tb(3,3),energy,count
real(8):: vol,fcart_tb(3,parini%nat),fcart_lj(3,parini%nat),strten_tb(6),strten_lj(6),energy_tb,energy_lj

!Check if the atomic types are allright
!call check_lenosky_tb_lj()

!Number of kpoints
n_in(1)=ka
n_in(2)=kb
n_in(3)=kc

!Enable k-points
if(ka==1.and.kb==1.and.kc==1) then 
   do_kpt_in=0
!Offset
   offset_in(:)=0
else
   do_kpt_in=1
!Offset
   offset_in(:)=1
endif

!TB or MEAM
tb_or_meam=0

latvec_ang=latvec*Bohr_Ang
xred=xred0
call n_rep_dim(latvec_ang,cut,nec1,nec2,nec3)
cellsearch=max(nec1,max(nec2,nec3))

call backtocell(parini%nat,latvec_ang,xred)
call rxyz_int2cart(latvec_ang,xred,rxyz,parini%nat)

fcart_tb=0.d0
call  lenoskytb(n_tb,rxyz(:,1:n_tb),fcart_tb(:,1:n_tb),energy_tb,count,n_silicon,latvec_ang(:,1),&
      &latvec_ang(:,2),latvec_ang(:,3),stress_tb(:,1),stress_tb(:,2),stress_tb(:,3),&
      &tb_or_meam,cellsearch,n_in(1),n_in(2),n_in(3),offset_in(1),offset_in(2),offset_in(3),do_kpt_in)
!  stress=-stress
  energy_tb=energy_tb/Ha_eV
  fcart_tb=fcart_tb/Ha_eV*Bohr_Ang
!Get full stress matrix
call getvol(latvec_ang,vol)
stress_tb=matmul(stress_tb,transpose(latvec_ang))/vol
!This matrix MUST be symmetric!!!
!!write(*,*) stress(:,1)           
!!write(*,*) stress(:,2)           
!!write(*,*) stress(:,3)           
  strten_tb(1) = stress_tb(1,1)
  strten_tb(2) = stress_tb(2,2)
  strten_tb(3) = stress_tb(3,3)
  strten_tb(6) = stress_tb(2,1)
  strten_tb(5) = stress_tb(3,1)
  strten_tb(4) = stress_tb(3,2)
!This is not very clear yet...
  strten_tb=strten_tb/Ha_eV*Bohr_Ang**3


!Get the LJ stuff
call lenosky_lj(parini,latvec,xred0,fcart_lj,strten_lj,energy_lj)

!Combine everything
energy=energy_tb+energy_lj
fcart=fcart_tb+fcart_lj
strten=strten_tb+strten_lj



end subroutine
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


        subroutine lenosky_lj(parini,latvec,xred0,fxyz,strten,etot)
!The binary parameters are defined as follows in this scheme:
!The parameter set 1 are for the LJ-Si interactions
!The parameter set 2 are for the LJ-H interactions
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
        real(8):: xred(3,parini%nat),fxyz(3,parini%nat),xred0(3,parini%nat),dxyz(3),r1red(3),r2red(3),rcut2(n_lj)
        real(8):: etot,latvec_x(3,3),rec_nec(3),enth,stressvol(3,3),vol,strten(6)
        real(8), allocatable:: rel(:,:)
        real(8):: latvec(3,3),latvec_ang(3,3),latvecinv(3,3),celldv(3,3),pressure,gradvol(3,3),stress(3,3),tmplat(3,3)
        real(8):: trans(3,3), transinv(3,3)
        real(8):: crossp(3),d,dd,dd2,dd6,dd12,dx,dy,dz,s,s2,s6,s12,rc,rc2,rc6,rc12,tt,t1,t2,t3
        real(8):: sil_sjl1,sil_sjl2,sil_sjl3, si1_sj1, si2_sj2, si3_sj3, tkmsum1,tkmsum2,tkmsum3, cutmax,epscur
        real(8):: rcut2_lj(n_lj),eps_lj_lj_fact,sigma_lj_lj_fact,double_count
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
        rcut2=rcuttblj*rcuttblj
        rcut2_lj=sigmatblj*2.d0**(1.d0/6.d0)*sigma_lj_lj_fact
        cutmax=max(maxval(rcuttblj),maxval(rcut2_lj))
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
!Interactions of LJ particles with Si and H
         do iat=n_silicon+n_h+1,parini%nat
            r1red(:)=xred(:,iat)*rec_nec
            do i=0,nec(1)-1
            do j=0,nec(2)-1
            do k=0,nec(3)-1
               do jat=1,n_silicon+n_h
                r2red(:)=xred(:,jat)*rec_nec
                r2red(1)=r2red(1)+real(i,8)*rec_nec(1)
                r2red(2)=r2red(2)+real(j,8)*rec_nec(2)
                r2red(3)=r2red(3)+real(k,8)*rec_nec(3)
                call pbc_distance0(latvec_x,r1red,r2red,dd,dxyz)
                if(dd.gt.rcut2(iat-n_silicon-n_h)) then
                   goto 1002
!                elseif(dd.lt.1.d-12) then
!                   goto 1002
                endif
                d=sqrt(dd)
                dx=-dxyz(1);  dy=-dxyz(2);  dz=-dxyz(3)
                dd2=1.d0/dd
                dd6=dd2*dd2*dd2
                dd12=dd6*dd6
                s=sigmatblj(iat-n_silicon-n_h)   !Sigma
                s2=s*s
                s6=s2*s2*s2
                s12=s6*s6
                rc=1.d0/rcuttblj(iat-n_silicon-n_h) !Cutoff
                rc2=rc*rc
                rc6=rc2*rc2*rc2
                rc12=rc6*rc6
                epscur=epstblj(iat-n_silicon-n_h)
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
         do iat=n_silicon+n_h+1,parini%nat
            r1red(:)=xred(:,iat)*rec_nec
            do i=0,nec(1)-1
            do j=0,nec(2)-1
            do k=0,nec(3)-1
               do jat=n_silicon+n_h+1,parini%nat
                r2red(:)=xred(:,jat)*rec_nec
                r2red(1)=r2red(1)+real(i,8)*rec_nec(1)
                r2red(2)=r2red(2)+real(j,8)*rec_nec(2)
                r2red(3)=r2red(3)+real(k,8)*rec_nec(3)
                call pbc_distance0(latvec_x,r1red,r2red,dd,dxyz)
                if(dd.gt.rcut2_lj(iat-n_silicon-n_h)) then
                   goto 1003
                elseif(dd.lt.1.d-12) then
                   goto 1003
                endif
                d=sqrt(dd)
                dx=-dxyz(1);  dy=-dxyz(2);  dz=-dxyz(3)
                dd2=1.d0/dd
                dd6=dd2*dd2*dd2
                dd12=dd6*dd6
                s=sigmatblj(iat-n_silicon-n_h)*sigma_lj_lj_fact  !Sigma
                s2=s*s
                s6=s2*s2*s2
                s12=s6*s6
!!                rc=1.d0/rcut(iat-n_silicon-n_h) !Cutoff
!!                rc2=rc*rc
!!                rc6=rc2*rc2*rc2
!!                rc12=rc6*rc6
                epscur=epstblj(iat-n_silicon-n_h)*eps_lj_lj_fact
                   if(abs(rcut2_lj(iat-n_silicon-n_h)-rcut2_lj(jat-n_silicon-n_h)).lt.1.d-15) then
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

end module interface_lenosky_tb_lj

!!!********************************************************
subroutine check_lenosky_tb_lj(parini)
use mod_parini, only: typ_parini
use global
use tb_lj_params
use defs_basis
implicit none
type(typ_parini), intent(in):: parini
integer:: iat
logical:: in_h,in_lj
!Check typat for consistentcy and provide n_silicon
in_h=.false.
in_lj=.false.
n_silicon=0
n_h=0
n_lj=0
do iat=1,parini%nat
   if(int(parini%znucl(parini%typat_global(iat))).ne.1.and.int(parini%znucl(parini%typat_global(iat))).ne.14.and.&
     &int(parini%znucl(parini%typat_global(iat))).ne.201) &
     &stop "Lenosky TB and LJ only implemented for Si and H and LJ particles"
   if(int(parini%znucl(parini%typat_global(iat)))==14) n_silicon=n_silicon+1
   if(int(parini%znucl(parini%typat_global(iat)))==14.and.(in_h.or.in_lj)) stop "Lenosky TB: First Si, then H, then LJ"
   if(int(parini%znucl(parini%typat_global(iat)))==1.and.in_lj) stop "Lenosky TB: First H, then LJ"
   if(int(parini%znucl(parini%typat_global(iat)))==1) then
       n_h=n_h+1
       in_h=.true.
   endif
   if(int(parini%znucl(parini%typat_global(iat))).gt.200) then 
       n_lj=n_lj+1 
       in_h=.true.
       in_lj=.true.
   endif
enddo
n_tb=n_silicon+n_h
end subroutine

!!!********************************************************

     subroutine lenosky_tb_lj_init_parameter()
     use global
     use tb_lj_params
!The LJ particles interact with Si and H atoms. The sigma, epsilon and cutoff are given for each 
!LJ particle in a list tb_lj_param.in
!There are NO interactions for the off diagonals and are ignored!!!
     !In general the parmeters are symmetric, i.e. sigma(A,B)=sigma(B,A)
     !In the end the array "Kinds" is allocated, if it has not already been done
     integer:: iat
     logical:: file_exists
     character(40):: filename
     if(.not.allocated(sigmatblj)) allocate(sigmatblj(n_lj)) 
     if(.not.allocated(epstblj))   allocate(epstblj(n_lj)) 
     if(.not.allocated(rcuttblj))  allocate(rcuttblj(n_lj)) 
      
     sigmalj=0.d0;epslj=0.d0
     file_exists=.false.
     filename="tb_lj_param.in"
     INQUIRE(FILE=trim(filename), EXIST=file_exists)
     if(file_exists) then
         open(unit=33,file=trim(filename))
           do iat=1,n_lj
           read(33,*) sigmatblj(iat),epstblj(iat),rcuttblj(iat)
           enddo
         close(33)
      else
         stop "Please provide the parametrization of TB-LJ in tb_lj_param.in"
      endif
     end subroutine lenosky_tb_lj_init_parameter

!********************************************************
