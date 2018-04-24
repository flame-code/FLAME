

module interface_blj
  use global
  use defs_basis
  use blj_params
!********************************************************
!                                                       *
! BINARY LENNARD JONES                                  *
! FIRST AND SECOND DERIVATIVE ARE 0.d0 AT CUTOFF        *
!                                                       *
!********************************************************

  implicit none

  private
  public :: &
    blj
 

contains

!!!subroutine blj(latvec,xred0,fxyz,strten,etot)
!!!implicit none
!!!integer:: iat
!!!real(8):: latvec(3,3),xred0(3,nat),fxyz(3,nat),strten(6),etot
!!!real(8):: dist_ang(6),latvec_rot(3,3),latvec_rot_inv(3,3),rotmat(3,3),stress(3,3)
!!!call latvec2dist_ang(dist_ang,latvec,pi)
!!!call dist_ang2latvec(dist_ang,latvec_rot,pi)
!!!call blj2(latvec_rot,xred0,fxyz,strten,etot)
!!!!Rotate forces back to original cell
!!!call rotmat_fcart_stress(latvec,latvec_rot,rotmat)
!!!do iat=1,nat
!!!   fxyz(:,iat)=matmul(rotmat,fxyz(:,iat))
!!!enddo
!!!call rotate_stresstensor(strten,rotmat)
!!!        stress(1,1) =  strten(1) 
!!!        stress(2,2) =  strten(2) 
!!!        stress(3,3) =  strten(3) 
!!!        stress(2,1) =  strten(6) 
!!!        stress(3,1) =  strten(5) 
!!!        stress(3,2) =  strten(4) 
!!!        stress(1,2) =  stress(2,1)
!!!        stress(1,3) =  stress(3,1)
!!!        stress(2,3) =  stress(3,2)
!!!!        stress = matmul(rotmat,stress)
!!!        strten(1) = stress(1,1)
!!!        strten(2) = stress(2,2)
!!!        strten(3) = stress(3,3)
!!!        strten(6) = stress(2,1)
!!!        strten(5) = stress(3,1)
!!!        strten(4) = stress(3,2)
!!!        stress=stress*Ha_eV/Bohr_Ang**3
!!!        write(*,*) stress(1,1), stress(1,1)
!!!        write(*,*) stress(2,1), stress(2,1)
!!!        write(*,*) stress(3,1), stress(3,1)
!!!        write(*,*) stress(1,2), stress(1,2)
!!!        write(*,*) stress(2,2), stress(2,2)
!!!        write(*,*) stress(3,2), stress(3,2)
!!!        write(*,*) stress(1,3), stress(1,3)
!!!        write(*,*) stress(2,3), stress(2,3)
!!!        write(*,*) stress(3,3), stress(3,3)
!!!         
!!!
!!!
!!!end subroutine

        subroutine blj(parini,latvec,xred0,fxyz,strten,etot)
        use mod_parini, only: typ_parini
        use blj_params
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

        implicit none
        type(typ_parini), intent(in):: parini
        integer:: iat, jat, kk, ll, l, nec(3), i, j, k, m
        real(8):: xred(3,parini%nat),fxyz(3,parini%nat),xred0(3,parini%nat),dxyz(3),r1red(3),r2red(3),rcut2(2,2)
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
        rcut2=rcut*rcut
        cutmax=maxval(rcut)
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
                s=sigmalj(parini%typat_global(iat),parini%typat_global(jat))   !Sigma
                s2=s*s
                s6=s2*s2*s2
                s12=s6*s6
                rc=1.d0/rcut(parini%typat_global(iat),parini%typat_global(jat)) !Cutoff
                rc2=rc*rc
                rc6=rc2*rc2*rc2
                rc12=rc6*rc6
                epscur=epslj(parini%typat_global(iat),parini%typat_global(jat))
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
                 
!!                   virial(1,1)=virial(1,1)-dx*t1!!;write(*,'(a,i3,i3,es25.15)')'vxx ',iat,jat,-dx*t1
!!                   virial(2,1)=virial(2,1)-dy*t1!!;write(*,'(a,i3,i3,es25.15)')'vyx ',iat,jat,-dy*t1
!!                   virial(3,1)=virial(3,1)-dz*t1!!;write(*,'(a,i3,i3,es25.15)')'vzx ',iat,jat,-dz*t1
!!                   virial(1,2)=virial(1,2)-dx*t2!!;write(*,'(a,i3,i3,es25.15)')'vyx ',iat,jat,-dx*t2
!!                   virial(2,2)=virial(2,2)-dy*t2!!;write(*,'(a,i3,i3,es25.15)')'vyy ',iat,jat,-dy*t2
!!                   virial(3,2)=virial(3,2)-dz*t2!!;write(*,'(a,i3,i3,es25.15)')'vzy ',iat,jat,-dz*t2
!!                   virial(1,3)=virial(1,3)-dx*t3!!;write(*,'(a,i3,i3,es25.15)')'vzx ',iat,jat,-dx*t3
!!                   virial(2,3)=virial(2,3)-dy*t3!!;write(*,'(a,i3,i3,es25.15)')'vzy ',iat,jat,-dy*t3
!!                   virial(3,3)=virial(3,3)-dz*t3!!;write(*,'(a,i3,i3,es25.15)')'vzz ',iat,jat,-dz*t3
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
!!        virial=virial/vol
!!        write(*,*) virial(1,1) ,stress(1,1)
!!        write(*,*) virial(2,1) ,stress(2,1)
!!        write(*,*) virial(3,1) ,stress(3,1)
!!        write(*,*) virial(1,2) ,stress(1,2)
!!        write(*,*) virial(2,2) ,stress(2,2)
!!        write(*,*) virial(3,2) ,stress(3,2)
!!        write(*,*) virial(1,3) ,stress(1,3)
!!        write(*,*) virial(2,3) ,stress(2,3)
!!        write(*,*) virial(3,3) ,stress(3,3)

        return
        end subroutine

end module interface_blj

!!!********************************************************
!!
!!        subroutine stress_volume(latvec,vol,pressure,stressvol)
!!        !This subroutine will compute the additional components to the negative
!!        !derivative of the enthalpy with respect to the cell variables
!!        !For the derivatives with respect to hij, the lattive vector components of the
!!        !lattece-matrix h, the relation on http://en.wikipedia.org/wiki/Determinant is used:
!!        !\frac{\partial \det(h)}{\partial h_{ij}}= \det(h)(h^{-1})_{ji}. 
!!        implicit none
!!        integer:: nat,i,j
!!        real(8):: latvec(3,3),vol,stressvol(3,3),inv_latvec(3,3),pressure
!!        stressvol=0.d0
!!        call invertmat(latvec,inv_latvec,3)
!!        do i=1,3
!!           stressvol(i,1)=-vol*inv_latvec(1,i)*pressure
!!           stressvol(i,2)=-vol*inv_latvec(2,i)*pressure
!!           stressvol(i,3)=-vol*inv_latvec(3,i)*pressure
!!        enddo
!!        end subroutine
!!
!!!********************************************************
!!
!!        subroutine cell_vol(nat,latvec,vol)
!!        !This subroutine will use the lattice vectors latvec to calculate the volume per atom of the simulation cell, vol
!!        implicit none
!!        integer:: nat
!!        real(8):: latvec(3,3),vol,a(3,3)
!!        a=latvec
!!        vol= a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)+&
!!             a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)
!!        vol=vol/real(nat,8)
!!        if(vol.lt.0.d0) stop "Negative volume!"
!!        end subroutine
!!
!!!********************************************************
!!

     subroutine blj_init_parameter(parini)
     use mod_parini, only: typ_parini
     use global
     use blj_params
     implicit none
     type(typ_parini), intent(in):: parini
     !We consider particle 1 as A- and particle 2 as B-particles     
     !In general the parmeters are symmetric, i.e. sigma(A,B)=sigma(B,A)
     !In the end the array "Kinds" is allocated, if it has not already been done
!     real(8):: alphalj!=2.5d0
     logical:: file_exists
     character(40):: filename
     sigmalj=0.d0;epslj=0.d0
     file_exists=.false.
     filename="blj_param.in"
     INQUIRE(FILE=trim(filename), EXIST=file_exists)
     if(minval(parini%znucl(:)).lt.201.or.maxval(parini%znucl(:)).gt.202) stop "BLJ particles must have znucl values of 201 and 202"
     if(file_exists) then
         open(unit=33,file=trim(filename))
         if(parini%ntypat_global==1) then
           read(33,*) sigmalj(1,1)
           read(33,*) epslj(1,1)
           read(33,*) alpha_lj
         elseif(parini%ntypat_global==2) then
           read(33,*) sigmalj(1,1),sigmalj(1,2),sigmalj(2,2)
                      sigmalj(2,1)=sigmalj(1,2)
           write(*,'(a,3(es15.7))') " # BLJ parameters: sigma(1,1), sigma(1,2)=sigma(2,1), sigma(2,2)   ",&
                 &sigmalj(1,1),sigmalj(1,2),sigmalj(2,2)
           read(33,*) epslj(1,1),epslj(1,2),epslj(2,2)
                      epslj(2,1)=epslj(1,2)
           write(*,'(a,3(es15.7))') " # BLJ parameters: eps(1,1), eps(1,2)=eps(2,1), eps(2,2)           ",&
                 &epslj(1,1),epslj(1,2),epslj(2,2)
           read(33,*) alpha_lj
           write(*,'(a,3(es15.7))') " # BLJ parameters: cut_alpha                                       ",alpha_lj
         else
           stop "Wrong number of atom types in BLJ"
         endif
         close(33)
      else
         stop "Please provide the parametrization of BLJ in blj_param.in"
      endif
!Hard coded
!!     sigmalj(1,1)=1.d0
!!     sigmalj(1,2)=0.8d0
!!     sigmalj(2,1)=sigmalj(1,2)
!!     sigmalj(2,2)=0.88d0
!!     epslj(1,1)=1.d0
!!     epslj(1,2)=1.5d0
!!     epslj(2,1)=epslj(1,2)
!!     epslj(2,2)=0.5d0
     rcut(1,1)=sigmalj(1,1)*alpha_lj
     rcut(1,2)=sigmalj(1,2)*alpha_lj
     rcut(2,1)=sigmalj(2,1)*alpha_lj
     rcut(2,2)=sigmalj(2,2)*alpha_lj

!ACHTUNG:TEMPORARY
!     rcut=9.d0
     end subroutine blj_init_parameter

!********************************************************
