!********************************************************
!                                                       *
! BINARY LENNARD JONES                                  *
! FIRST AND SECOND DERIVATIVE ARE 0.d0 AT CUTOFF        *
!                                                       *
!********************************************************
!*****************************************************************************************
subroutine init_lennardjones_vc(nat,sat)
    use mod_potential, only: sigmalj, epslj, rcut, kinds
    !We consider particle 1 as A- and particle 2 as B-particles     
    !In general the parmeters are symmetric, i.e. sigma(A,B)=sigma(B,A)
    !In the end the array "kinds" is allocated, if it has not already been done
    implicit none
    integer, intent(in):: nat
    character(5), intent(in):: sat(nat)
    !local variables
    real(8), parameter:: alphalj=2.5d0
    integer:: iat
    sigmalj(1,1)=1.d0
    sigmalj(1,2)=0.8d0
    sigmalj(2,1)=sigmalj(1,2)
    sigmalj(2,2)=0.88d0
    epslj(1,1)=1.d0
    epslj(1,2)=1.5d0
    epslj(2,1)=epslj(1,2)
    epslj(2,2)=0.5d0
    rcut(1,1)=sigmalj(1,1)*alphalj
    rcut(1,2)=sigmalj(1,2)*alphalj
    rcut(2,1)=sigmalj(2,1)*alphalj
    rcut(2,2)=sigmalj(2,2)*alphalj
    !atmass=1.d0
    !atmassinv=1.d0
    allocate(kinds(nat))
    do iat=1,nat
        if(trim(sat(iat))=="A") then
            kinds(iat)=1
        elseif(trim(sat(iat))=="B") then
            kinds(iat)=2
        else
            stop "ERROR: wrong atom types: only A and B are accepted by binary LJ"
        endif
    enddo
end subroutine init_lennardjones_vc
!*****************************************************************************************
!subroutine lennardjones_vc(nat,latvec,xred0,fxyz,celldv,stress,pressure,etot,enth)
subroutine lennardjones_vc(iproc,nat,xred0,latvec,pressure,fxyz,celldv,stress,etot,enth)
    use mod_potential, only: sigmalj, epslj, rcut, kinds
    !energy/forces for truncated Lennard Jones potential according to PRA 8, 1504 (1973)
    !input: nat: number of atoms totally
    !       xred: reduced positions of atoms
    !       latvec: lattice vectors of the periodic cell
    !       pressure: pressure
    !output:etot: energy
    !       enth: enthalpy
    !       stress: the stress tensor
    !       fxyz: forces (negative derivative of energy with respect to positions
    !       celldv: the virial stress, the negative derivative of the energy
    !               with respect to the cell variables
    implicit none
    integer, intent(in):: iproc, nat
    integer:: iat, jat, kk, ll, l, ipb, ipe, nnbrx, nec(3), i, j, k, m
    integer, allocatable:: lsta(:,:), lstb(:)
    real(8):: xred(3,nat),fxyz(3,nat),xred0(3,nat),dxyz(3),r1red(3),r2red(3),rcut2(2,2)
    real(8):: etot,latvec_x(3,3),rec_nec(3),enth,stressvol(3,3),vol
    real(8), allocatable:: rel(:,:)
    real(8):: latvec(3,3),latvecinv(3,3),celldv(3,3),pressure,gradvol(3,3),stress(3,3),tmplat(3,3)
    real(8):: trans(3,3), transinv(3,3)
    real(8):: crossp(3),d,dd,dd2,dd6,dd12,dx,dy,dz,s,s2,s6,s12,rc,rc2,rc6,rc12,tt,t1,t2,t3
    real(8):: sil_sjl1,sil_sjl2,sil_sjl3, si1_sj1, si2_sj2, si3_sj3, tkmsum1,tkmsum2,tkmsum3, cutmax,epscur
    xred=xred0
    etot=0.d0
    fxyz(:,:)=0.d0
    celldv(:,:)=0.d0
    rcut2=rcut*rcut
    cutmax=maxval(rcut)
    call backtocell_alborz(nat,latvec,xred)
    !==========================================
    !do iat=1,nat
    !    write(9999999,'(i4,3es12.3)') iat, xred(1,iat),xred(2,iat),xred(3,iat)
    !enddo
    !========================================== 
    trans=latvec
    call invertmat_alborz(trans,transinv)
    call n_rep_dim_alborz(latvec,2.d0*cutmax,nec(1),nec(2),nec(3))
    !==========================================
    !write(*,*) "nec1,nec2,nec3:",nec(1),nec(2),nec(3)
    !========================================== 
    !Expand cell
    do i=1,3
    latvec_x(:,i)=real(nec(i),8)*latvec(:,i)
    !Adjust reduced coordinates
    rec_nec(i)=1.d0/real(nec(i),8)
    enddo
    !Get the rest
    do iat=1,nat
    r1red(:)=xred(:,iat)*rec_nec
    do i=0,nec(1)-1
    do j=0,nec(2)-1
    do k=0,nec(3)-1
    do jat=1,nat
    r2red(:)=xred(:,jat)*rec_nec
    r2red(1)=r2red(1)+real(i,8)*rec_nec(1)
    r2red(2)=r2red(2)+real(j,8)*rec_nec(2)
    r2red(3)=r2red(3)+real(k,8)*rec_nec(3)
    call pbc_distance1_alborz(latvec_x,r1red,r2red,dd,dxyz)
    if(dd.gt.rcut2(kinds(iat),kinds(jat))) then
    goto 1002
    elseif(dd.lt.1.d-12) then
    goto 1002
    endif
    d=sqrt(dd)
    dx=-dxyz(1);  dy=-dxyz(2);  dz=-dxyz(3)
    dd2=1.d0/dd
    dd6=dd2*dd2*dd2
    dd12=dd6*dd6
    s=sigmalj(kinds(iat),kinds(jat))   !Sigma
    s2=s*s
    s6=s2*s2*s2
    s12=s6*s6
    rc=1.d0/rcut(kinds(iat),kinds(jat)) !Cutoff
    rc2=rc*rc
    rc6=rc2*rc2*rc2
    rc12=rc6*rc6
    epscur=epslj(kinds(iat),kinds(jat))
    etot=etot+4.d0*epscur*((s12*dd12-s6*dd6)+(6.d0*s12*rc12-3.d0*s6*rc6)*dd*rc2-7.d0*s12*rc12+4.d0*s6*rc6)
    tt=24.d0*epscur*dd2*(2.d0*s12*dd12-s6*dd6)-8.d0*epscur*(6.d0*s12*rc12-3.d0*s6*rc6)*rc2
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
    !   stress(1,1)=stress(1,1)+t1*dx
    !   stress(1,2)=stress(1,2)+t1*dy
    !   stress(1,3)=stress(1,3)+t1*dz
    !   stress(2,1)=stress(2,1)+t2*dx
    !   stress(2,2)=stress(2,2)+t2*dy
    !   stress(2,3)=stress(2,3)+t2*dz
    !   stress(3,1)=stress(3,1)+t3*dx
    !   stress(3,2)=stress(3,2)+t3*dy
    !   stress(3,3)=stress(3,3)+t3*dz
    1002 continue
    enddo
    enddo
    enddo
    enddo
    enddo
    etot=etot*0.5d0
    celldv=celldv*0.5d0
    call cell_vol(nat,latvec,vol)
    vol=vol*real(nat,8)
    enth=etot+pressure*vol
    call stress_volume_alborz(latvec,vol,pressure,stressvol)
    !Transform the forces on the lattice vectors into the stress tensore according to the paper Tomas Bucko and Jurg Hafner
    do i=1,3
    tmplat(:,i)=latvec(i,:)
    enddo
    !The real stress tensor
    stress=-matmul(celldv,tmplat)/vol
    !celldv has all forces to minimize the enthalpy directly on the lattice
    celldv=celldv+stressvol
    return
end subroutine lennardjones_vc
!*****************************************************************************************
subroutine stress_volume_alborz(latvec,vol,pressure,stressvol)
    !This subroutine will compute the additional components to the negative
    !derivative of the enthalpy with respect to the cell variables
    !For the derivatives with respect to hij, the lattive vector components of the
    !lattece-matrix h, the relation on http://en.wikipedia.org/wiki/Determinant is used:
    !\frac{\partial \det(h)}{\partial h_{ij}}= \det(h)(h^{-1})_{ji}. 
    implicit none
    integer:: nat,i,j
    real(8):: latvec(3,3),vol,stressvol(3,3),inv_latvec(3,3),pressure
    stressvol=0.d0
    call invertmat_alborz(latvec,inv_latvec)
    do i=1,3
        stressvol(i,1)=-vol*inv_latvec(1,i)*pressure
        stressvol(i,2)=-vol*inv_latvec(2,i)*pressure
        stressvol(i,3)=-vol*inv_latvec(3,i)*pressure
    enddo
end subroutine stress_volume_alborz
!*****************************************************************************************
subroutine cell_vol(nat,latvec,vol)
    !This subroutine will use the lattice vectors latvec to calculate the volume per atom of the simulation cell, vol
    implicit none
    integer:: nat
    real(8):: latvec(3,3),vol,a(3,3)
    a=latvec
    vol= a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)+&
    a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)
    vol=vol/real(nat,8)
    if(vol.lt.0.d0) stop "Negative volume!"
end subroutine cell_vol
!*****************************************************************************************
