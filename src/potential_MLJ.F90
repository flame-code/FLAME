module interface_mlj
  use global
  use defs_basis
  use mlj_params
!********************************************************
!                                                       *
! BINARY LENNARD JONES                                  *
! FIRST AND SECOND DERIVATIVE ARE 0.d0 AT CUTOFF        *
!                                                       *
!********************************************************

  implicit none

  private
  public :: &
    mlj
 

contains

        subroutine mlj(parini,latvec,xred0,fxyz,strten,etot)
        use mod_parini, only: typ_parini
        use mlj_params
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

!This lennard jones routine can handle an arbitrary number of different LJ particles which ONLY interact in a pairwise
!manner with parameters that were initiallized in the beginning with the init command

        implicit none
        type(typ_parini), intent(in):: parini
        integer:: iat, jat, kk, ll, l, nec(3), i, j, k, m
        real(8):: xred(3,parini%nat),fxyz(3,parini%nat),xred0(3,parini%nat),dxyz(3),r1red(3),r2red(3),rcut2(parini%ntypat_global,parini%ntypat_global)
        real(8):: etot,latvec_x(3,3),rec_nec(3),enth,stressvol(3,3),vol,strten(6)
        real(8), allocatable:: rel(:,:)
        real(8):: latvec(3,3),latvec_ang(3,3),latvecinv(3,3),celldv(3,3),pressure,gradvol(3,3),stress(3,3),tmplat(3,3)
        real(8):: trans(3,3), transinv(3,3)
        real(8):: crossp(3),d,dd,dd2,dd6,dd12,dx,dy,dz,s,s2,s6,s12,rc,rc2,rc6,rc12,tt,t1,t2,t3
        real(8):: sil_sjl1,sil_sjl2,sil_sjl3, si1_sj1, si2_sj2, si3_sj3, tkmsum1,tkmsum2,tkmsum3, cutmax,epscur,virial(3,3)
        logical:: truncated
!If truncated is true we will use the shifted form according to PRA 8, 1504 (1973), and else a simple crude cutoff
        truncated=.true.
!Convert everything from "internal atomic" units to the real "angstrom-ev" units
        latvec_ang=latvec*Bohr_Ang

        xred=xred0
        etot=0.d0
        fxyz(:,:)=0.d0
        celldv(:,:)=0.d0
        rcut2=rcutmlj*rcutmlj
        cutmax=maxval(rcutmlj)
        call backtocell(parini%nat,latvec_ang,xred) 
        trans=latvec_ang

         call invertmat(trans,transinv,3) 
         call n_rep_dim(latvec_ang,2.d0*cutmax,nec(1),nec(2),nec(3))

    virial=0.d0
!Expand cell
    do i=1,3
      latvec_x(:,i)=real(nec(i),8)*latvec_ang(:,i)
    !Adjust reduced coordinates
      rec_nec(i)=1.d0/real(nec(i),8)
    enddo

!Get the rest
         do iat=1,parini%nat
            r1red(:)=xred(:,iat)*rec_nec
            do i=0,nec(1)-1
            do j=0,nec(2)-1
            do k=0,nec(3)-1
               do jat=1,parini%nat
                r2red(:)=xred(:,jat)*rec_nec
                r2red(1)=r2red(1)+real(i,8)*rec_nec(1)
                r2red(2)=r2red(2)+real(j,8)*rec_nec(2)
                r2red(3)=r2red(3)+real(k,8)*rec_nec(3)
                call pbc_distance0(latvec_x,r1red,r2red,dd,dxyz)
                if(dd.gt.rcut2(parini%typat_global(iat),parini%typat_global(jat))) then
                   goto 1002
                elseif(dd.lt.1.d-12) then
                   goto 1002
                endif
                d=sqrt(dd)
                dx=-dxyz(1);  dy=-dxyz(2);  dz=-dxyz(3)
                dd2=1.d0/dd
                dd6=dd2*dd2*dd2
                dd12=dd6*dd6
                s=sigmamlj(parini%typat_global(iat),parini%typat_global(jat))   !Sigma
                s2=s*s
                s6=s2*s2*s2
                s12=s6*s6
                rc=1.d0/rcutmlj(parini%typat_global(iat),parini%typat_global(jat)) !Cutoff
                rc2=rc*rc
                rc6=rc2*rc2*rc2
                rc12=rc6*rc6
                epscur=epsmlj(parini%typat_global(iat),parini%typat_global(jat))
                   etot=etot+4.d0*epscur*((s12*dd12-s6*dd6)) !Simple LJ
                   if(truncated) etot=etot+4.d0*epscur*((6.d0*s12*rc12-3.d0*s6*rc6)*dd*rc2-7.d0*s12*rc12+4.d0*s6*rc6)
                   tt=24.d0*epscur*dd2*(2.d0*s12*dd12-s6*dd6) !Simple LJ
                   if(truncated) tt=tt-8.d0*epscur*(6.d0*s12*rc12-3.d0*s6*rc6)*rc2
                   t1=dx*tt ; t2=dy*tt ; t3=dz*tt
                   fxyz(1,iat)=fxyz(1,iat)+t1
                   fxyz(2,iat)=fxyz(2,iat)+t2
                   fxyz(3,iat)=fxyz(3,iat)+t3

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

        virial=virial/2.d0
        etot=etot*0.5d0
        etot=etot/Ha_eV
        fxyz=fxyz/Ha_eV*Bohr_Ang
        celldv=celldv*0.5d0
        call getvol(latvec_ang,vol)
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
        return
        end subroutine

end module interface_mlj

     subroutine mlj_init_parameter(parini)
     use mod_parini, only: typ_parini
     use global
     use mlj_params
     type(typ_parini), intent(inout):: parini
!We consider only the epsilon and sigma for each one of the lennard jones
!elements. The rest, meaning the interaction parameters, will be constructed
!from those elements by some mixing rule
     logical:: file_exists
     character(40):: filename
!Allocate the interaction parameters     
     if(.not.allocated(sigmamlj)) allocate(sigmamlj(parini%ntypat_global,parini%ntypat_global)) !sigma(i,j) is the sigma between particle i and particle j
     if(.not.allocated(epsmlj  )) allocate(epsmlj  (parini%ntypat_global,parini%ntypat_global)) !eps(i,j) is the epsilon between particle i and particle j
     if(.not.allocated(rcutmlj )) allocate(rcutmlj (parini%ntypat_global,parini%ntypat_global)) !rcut(i,j) is the cutoff distance between particle i and j
     if(.not.allocated(alphamlj)) allocate(alphamlj(parini%ntypat_global,parini%ntypat_global)) !cutoff discatnce

     sigmamlj=0.d0;epsmlj=0.d0
!Get the LJ parameters from a lookup table
         do ityp=1,parini%ntypat_global
           call mlj_atmdata(parini%amu,sigmamlj(ityp,ityp),epsmlj(ityp,ityp),parini%rcov(ityp),parini%char_type(ityp),parini%znucl(ityp))
           alphamlj(ityp,ityp)=2.5d0
         enddo

!We can setup simple rules for mixed LJ parameters

!SIGMA:
     do ityp=1,parini%ntypat_global
       do jtyp=1,ityp-1
!Geometric average
         sigmamlj(ityp,jtyp)=sqrt(sigmamlj(ityp,ityp)*sigmamlj(jtyp,jtyp))
!Arithmetic average
         sigmamlj(ityp,jtyp)=0.4d0*(sigmamlj(ityp,ityp)+sigmamlj(jtyp,jtyp))
         sigmamlj(jtyp,ityp)=sigmamlj(ityp,jtyp)
       enddo
     enddo

!EPSILON:
     do ityp=1,parini%ntypat_global
       do jtyp=1,ityp-1
!Geometric average
         epsmlj(ityp,jtyp)=sqrt(epsmlj(ityp,ityp)*epsmlj(jtyp,jtyp))
!Arithmetic average
         epsmlj(ityp,jtyp)=0.6d0*(epsmlj(ityp,ityp)+epsmlj(jtyp,jtyp))
         epsmlj(jtyp,ityp)=epsmlj(ityp,jtyp)
       enddo
     enddo

     rcutmlj(:,:)=sigmamlj(:,:)*2.5d0

write(*,*) "Epsilon"
do ityp=1,parini%ntypat_global
 write(*,*) epsmlj(:,ityp)
enddo
write(*,*) "Sigma"
do ityp=1,parini%ntypat_global
 write(*,*) sigmamlj(:,ityp)
enddo
write(*,*) "Rcut"
do ityp=1,parini%ntypat_global
 write(*,*) rcutmlj(:,ityp)
enddo
     end subroutine mlj_init_parameter

!********************************************************
