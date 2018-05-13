module interface_core_repulsion
  use global
  use defs_basis

  implicit none

  private
  public :: &
    core_repulsion
 

contains

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


        subroutine core_repulsion(parini,latvec,xred0,fxyz,strten,etot)
!This routine will compute a lennard-jones-like core repulsion. It is
!useful for codes which have not a well defined repulsion potential when atoms
!get too close to each other, like the DFTB code or similar.
!The repulsive part behaves like 1/r^12 without an attractive part.
!The cutoff/onset of repulsion is currently defined as sigma_lj_fact times the sum of
!the covalent radii, which can be changed according to the problem desired
        use mod_parini, only: typ_parini
        implicit none
        type(typ_parini), intent(in):: parini
        integer:: iat, jat, kk, ll, l, nec(3), i, j, k, m
        real(8):: xred(3,parini%nat),fxyz(3,parini%nat),xred0(3,parini%nat),dxyz(3),r1red(3),r2red(3)
        real(8):: etot,latvec_x(3,3),rec_nec(3),enth,stressvol(3,3),vol,strten(6)
        real(8), allocatable:: rel(:,:)
        real(8):: latvec(3,3),latvec_ang(3,3),latvecinv(3,3),celldv(3,3),pressure,gradvol(3,3),stress(3,3),tmplat(3,3)
        real(8):: trans(3,3), transinv(3,3)
        real(8):: crossp(3),d,dd,dd2,dd6,dd12,dx,dy,dz,s,s2,s6,s12,rc,rc2,rc6,rc12,tt,t1,t2,t3
        real(8):: sil_sjl1,sil_sjl2,sil_sjl3, si1_sj1, si2_sj2, si3_sj3, tkmsum1,tkmsum2,tkmsum3, cutmax,epscur
        real(8):: rcut2_lj(parini%ntypat_global,parini%ntypat_global),eps_lj,sigma_lj_fact,double_count
        real(8):: sigma_lj(parini%ntypat_global,parini%ntypat_global)
        logical:: truncated
!Set up epsilon and sigma parameters
        eps_lj=1.d0
        sigma_lj_fact=0.7d0
!If truncated is true we will use the shifted form according to PRA 8, 1504 (1973), and else a simple crude cutoff
        truncated=.true.
!Convert everything from "internal atomic" units to the real "angstrom-ev" units
        latvec_ang=latvec*Bohr_Ang
!Set up sigma matrix
        sigma_lj=0.d0
        do i=1,parini%ntypat_global
          do j=1,parini%ntypat_global
!Exclude interactions between LJ particles, since they are already LJ... hahaha
             if(int(parini%znucl(i)).lt.200.and.int(parini%znucl(j)).lt.200) sigma_lj(i,j)=(parini%rcov(i)+parini%rcov(j))*Bohr_Ang*sigma_lj_fact
          enddo
        enddo
        xred=xred0
        etot=0.d0
        fxyz(:,:)=0.d0
        celldv(:,:)=0.d0
        rcut2_lj=sigma_lj*2.d0**(1.d0/6.d0) !This cutoff is due to the truncation at the minimum of the well, since only repulsion is used!!!
        cutmax=maxval(rcut2_lj)
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

!Only the repulsive interactions between the pseudo LJ particles
         do iat=1,parini%nat !iat are th LJ particles 
            r1red(:)=xred(:,iat)*rec_nec
            do i=0,nec(1)-1
            do j=0,nec(2)-1
            do k=0,nec(3)-1
               do jat=1,parini%nat !jat are the other LJ particles
                r2red(:)=xred(:,jat)*rec_nec
                r2red(1)=r2red(1)+real(i,8)*rec_nec(1)
                r2red(2)=r2red(2)+real(j,8)*rec_nec(2)
                r2red(3)=r2red(3)+real(k,8)*rec_nec(3)
                call pbc_distance0(latvec_x,r1red,r2red,dd,dxyz)
!                if(dd.gt.rcut2_lj(iat-n_silicon-n_h)) then
                if(dd.gt.rcut2_lj(parini%typat_global(iat),parini%typat_global(jat))) then
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
                s=sigma_lj(parini%typat_global(iat),parini%typat_global(jat))  !Sigma
                if(s.lt.1.d-14) goto 1003
                s2=s*s
                s6=s2*s2*s2
                s12=s6*s6
!!                rc=1.d0/rcut(iat-n_silicon-n_h) !Cutoff
!!                rc2=rc*rc
!!                rc6=rc2*rc2*rc2
!!                rc12=rc6*rc6
!                epscur=epsvoidlj(iat-n_silicon-n_h)*eps_lj_lj_fact
                epscur=eps_lj
                      double_count=1.d0
!                   write(*,'(a,6es15.7)') "epslj",epscur,double_count,s12*dd12,s6*dd6,(4.d0*epscur*((s12*dd12-s6*dd6))+epscur)*double_count,d
                   etot=etot+(4.d0*epscur*((s12*dd12-s6*dd6))+epscur)*double_count !Repulsive LJ with shift
                   tt=24.d0*epscur*dd2*(2.d0*s12*dd12-s6*dd6)*double_count !Simple LJ
                   t1=dx*tt ; t2=dy*tt ; t3=dz*tt
                   fxyz(1,iat)=fxyz(1,iat)+t1
                   fxyz(2,iat)=fxyz(2,iat)+t2
                   fxyz(3,iat)=fxyz(3,iat)+t3
!!!Actio=Reactio
!!                   fxyz(1,jat)=fxyz(1,jat)-t1
!!                   fxyz(2,jat)=fxyz(2,jat)-t2
!!                   fxyz(3,jat)=fxyz(3,jat)-t3

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
        etot=etot*0.5d0
        etot=etot/Ha_eV
        fxyz=fxyz/Ha_eV*Bohr_Ang
        celldv=celldv*0.5d0
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

end module interface_core_repulsion
