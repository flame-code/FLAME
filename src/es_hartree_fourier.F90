!*****************************************************************************************
subroutine get_psolver_fourier(parini,poisson,atoms,gausswidth,ehartree,g)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_electrostatics, only: typ_poisson
    use time_profiling
    !use mod_timing , only: TCAT_PSOLVER
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_poisson),intent(inout):: poisson
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(in):: gausswidth(atoms%nat)
    real(8), intent(out):: ehartree, g(atoms%nat)
    !local variables
    real(8):: alphasq
    integer:: iat
    real(8), allocatable:: gwsq(:), ratred(:,:), fat(:,:)
    real(8):: stress(3,3), celldv(3,3)
    gwsq=f_malloc([1.to.atoms%nat],id='gwsq')
    ratred=f_malloc([1.to.3,1.to.atoms%nat],id='ratred')
    fat=f_malloc([1.to.3,1.to.atoms%nat],id='fat')
    call update_ratp(atoms)
    if(poisson%gw_identical) then
         !gwsq(1:atoms%nat)=ewaldwidth(1:atoms%nat)**2
         !call get_psolver_fourier_various(atoms%nat,atoms%rat,ratred,atoms%qat, &
         !    atoms%cellvec,gwsq,ecut,ehartree,fat,g,stress,atoms%celldv)
        alphasq=poisson%alpha**2
        call get_psolver_fourier_identical(parini%iverbose,atoms%nat,atoms%ratp,ratred,atoms%qat, &
            atoms%cellvec,alphasq,poisson%ecut,ehartree,fat,g,stress,atoms%celldv)
     else
        gwsq(1:atoms%nat)=gausswidth(1:atoms%nat)**2
        call get_psolver_fourier_various(parini%iverbose,atoms%nat,atoms%ratp,ratred,atoms%qat,atoms%cellvec, &
            gwsq,poisson%ecut,ehartree,fat,g,stress,atoms%celldv)
    end if
    do iat=1,atoms%nat
        atoms%fat(1,iat)=atoms%fat(1,iat)+fat(1,iat)
        atoms%fat(2,iat)=atoms%fat(2,iat)+fat(2,iat)
        atoms%fat(3,iat)=atoms%fat(3,iat)+fat(3,iat)
    enddo
    atoms%stress(1,1)=atoms%stress(1,1)+stress(1,1)
    atoms%stress(2,1)=atoms%stress(2,1)+stress(2,1)
    atoms%stress(3,1)=atoms%stress(3,1)+stress(3,1)
    atoms%stress(1,2)=atoms%stress(1,2)+stress(1,2)
    atoms%stress(2,2)=atoms%stress(2,2)+stress(2,2)
    atoms%stress(3,2)=atoms%stress(3,2)+stress(3,2)
    atoms%stress(1,3)=atoms%stress(1,3)+stress(1,3)
    atoms%stress(2,3)=atoms%stress(2,3)+stress(2,3)
    atoms%stress(3,3)=atoms%stress(3,3)+stress(3,3)
    atoms%celldv(1,1)=atoms%celldv(1,1)+celldv(1,1)
    atoms%celldv(2,1)=atoms%celldv(2,1)+celldv(2,1)
    atoms%celldv(3,1)=atoms%celldv(3,1)+celldv(3,1)
    atoms%celldv(1,2)=atoms%celldv(1,2)+celldv(1,2)
    atoms%celldv(2,2)=atoms%celldv(2,2)+celldv(2,2)
    atoms%celldv(3,2)=atoms%celldv(3,2)+celldv(3,2)
    atoms%celldv(1,3)=atoms%celldv(1,3)+celldv(1,3)
    atoms%celldv(2,3)=atoms%celldv(2,3)+celldv(2,3)
    atoms%celldv(3,3)=atoms%celldv(3,3)+celldv(3,3)
    call f_free(fat)
    call f_free(gwsq)
    call f_free(ratred)
end subroutine get_psolver_fourier
!*****************************************************************************************
subroutine get_psolver_fourier_various(iverbose,nat,rat,ratred,qat,cv,gwsq,ecut,ehartree,fat,eqd,stress,celldv)
    implicit none
    integer, intent(in):: iverbose, nat
    real(8), intent(in):: rat(3,nat), qat(nat)
    real(8), intent(in):: cv(3,3), gwsq(nat), ecut
    real(8), intent(out):: ratred(3,nat), fat(3,nat), eqd(nat), ehartree, stress(3,3), celldv(3,3)
    !local variables
    integer:: m1, m2, m3, m1_max, m2_max, m3_max, mtot, msq, msq_max
    integer:: iat, i, j, k, n, l
    real(8):: kx, ky, kz, vksq, vksqinv
    real(8):: kxp1, kyp1, kzp1, kxp12, kyp12, kzp12
    real(8):: sum_en, factor
    real(8):: vol, pi, pisq, sf, bnrm1, bnrm2, bnrm3
    real(8):: tt0, tt1, tt2, tt3, ttc, tts, twopi, fourpisq, dpx, dpy, dpz
    real(8):: recvec(3,3), g(3,3), cvinv(3,3), cvinv_vol(3,3)
    real(8):: den_dginv(3,3), dk2_dginv(3,3), dginv_da(3,3,3,3), ginv(3,3)
    real(8), allocatable:: tarr1(:), tarr2(:), tarr3(:)
    !integer, save:: icall=-1
    !icall=icall+1
    allocate(tarr1(1:nat))
    allocate(tarr2(1:nat))
    allocate(tarr3(1:nat))
    pi=4.d0*atan(1.d0)
    pisq=pi**2
    twopi=8.d0*atan(1.d0)
    fourpisq=4.d0*pisq
    msq_max=10000000
    !-------- volume of unitcell -----------------------------------
    call getvol_alborz(cv,vol)
    !-------- atomic positions in reduced coordinates --------------
    call rxyz_cart2int_alborz(nat,cv,rat,ratred)
    !-------- reciprocal lattice vectors  --------------------------
    call mycross(cv(1,2),cv(1,3),recvec(1,1))
    call mycross(cv(1,3),cv(1,1),recvec(1,2))
    call mycross(cv(1,1),cv(1,2),recvec(1,3))
    recvec(1:3,1:3)=recvec(1:3,1:3)/vol
    !write(*,'(a)') 'reciprocal lattice vectors'
    !write(*,'(3es14.5)') recvec(1,1),recvec(1,2),recvec(1,3)
    !write(*,'(3es14.5)') recvec(2,1),recvec(2,2),recvec(2,3)
    !write(*,'(3es14.5)') recvec(3,1),recvec(3,2),recvec(3,3)
    !--------- setting matrix g_inverse-----------------------------
    do i=1,3
        do j=1,3
            g(i,j)=cv(1,i)*cv(1,j)+cv(2,i)*cv(2,j)+cv(3,i)*cv(3,j)
        enddo
    enddo
    call invertmat_alborz(g,ginv)
    call invertmat_alborz(cv,cvinv)
    do i=1,3
    do j=1,3
        cvinv_vol(i,j)=vol*cvinv(j,i)
    enddo
    enddo
    !write(*,*) ginv(1:3,1:3)
    !-------- setting maximum k-vector ----------------------------
    bnrm1=sqrt(recvec(1,1)**2+recvec(2,1)**2+recvec(3,1)**2)
    bnrm2=sqrt(recvec(1,2)**2+recvec(2,2)**2+recvec(3,2)**2)
    bnrm3=sqrt(recvec(1,3)**2+recvec(2,3)**2+recvec(3,3)**2)
    m1_max=ceiling(sqrt(2.d0*ecut)/(bnrm1*2.d0*pi))
    m2_max=ceiling(sqrt(2.d0*ecut)/(bnrm2*2.d0*pi))
    m3_max=ceiling(sqrt(2.d0*ecut)/(bnrm3*2.d0*pi))
    if(iverbose>=3) then
        write(*,'(a,3i4)') 'm1_max,m2_max,m3_max',m1_max,m2_max,m3_max
    endif
    
    mtot=0
    sum_en=0.d0
    fat(1:3,1:nat)=0.d0
    eqd(1:nat)=0.d0
    den_dginv=0.d0

    do m1=0,m1_max
        if(m1==0) then
            factor=1.d0
        else
            factor=2.d0
        endif
        kxp1=real(m1,8)*recvec(1,1)
        kyp1=real(m1,8)*recvec(2,1)
        kzp1=real(m1,8)*recvec(3,1)
        dk2_dginv(1,1)=twopi**2*m1**2
        do m2=-m2_max,m2_max
            kxp12=kxp1+real(m2,8)*recvec(1,2)
            kyp12=kyp1+real(m2,8)*recvec(2,2)
            kzp12=kzp1+real(m2,8)*recvec(3,2)
            dk2_dginv(1,2)=twopi**2*m1*m2
            dk2_dginv(2,1)=dk2_dginv(1,2)
            dk2_dginv(2,2)=twopi**2*m2**2
            do m3=-m3_max,m3_max
                kx=kxp12+real(m3,8)*recvec(1,3)
                ky=kyp12+real(m3,8)*recvec(2,3)
                kz=kzp12+real(m3,8)*recvec(3,3)
                dk2_dginv(1,3)=twopi**2*m1*m3
                dk2_dginv(3,1)=dk2_dginv(1,3)
                dk2_dginv(2,3)=twopi**2*m2*m3
                dk2_dginv(3,2)=dk2_dginv(2,3)
                dk2_dginv(3,3)=twopi**2*m3**2
                !if(m1==0 .and. m2==5 .and. m3==5) then
                !    write(*,'(5es14.5)') kx,ky,kz,real(m2,8)*recvec(2,2),real(m3,8)*recvec(3,3)
                !    stop
                !endif
                msq=m1*m1+m2*m2+m3*m3
                if((msq<msq_max).and.(msq/=0)) then
                    mtot=mtot+1
                    !if(mtot>m_max) exit 
                    vksq=(kx*kx+ky*ky+kz*kz)*fourpisq
                    vksqinv=1.d0/vksq
                    !calculating the structure factor and its derivative
                    do iat=1,nat
                        tt1=(m1*ratred(1,iat)+m2*ratred(2,iat)+m3*ratred(3,iat))*twopi
                        tarr1(iat)=cos(tt1)
                        tarr2(iat)=sin(tt1)
                        tarr3(iat)=exp(-gwsq(iat)*vksq*0.25d0)
                    enddo
                    ttc=0.d0
                    tts=0.d0
                    do iat=1,nat
                        ttc=ttc+qat(iat)*tarr1(iat)*tarr3(iat)
                        tts=tts+qat(iat)*tarr2(iat)*tarr3(iat)
                    enddo
                    sf=(ttc**2+tts**2)
                    !write(*,'(3i4,i7,4es14.5,es19.10)') m1,m2,m3,mtot,kx,ky,kz,sf,sum_en
                    tt1=factor*sf*vksqinv
                    sum_en=sum_en+tt1
                    tt2=tt1*vksqinv
                    do l=1,3
                        do n=1,3
                            den_dginv(n,l)=den_dginv(n,l)+tt2*dk2_dginv(n,l)
                        enddo
                    enddo
                    do iat=1,nat
                        tt0=factor*tarr3(iat)*vksqinv
                        tt1=(tarr2(iat)*ttc-tarr1(iat)*tts)*tt0*qat(iat)
                        tt2=(tarr1(iat)*ttc+tarr2(iat)*tts)*tt0
                        fat(1,iat)=fat(1,iat)+kx*tt1
                        fat(2,iat)=fat(2,iat)+ky*tt1
                        fat(3,iat)=fat(3,iat)+kz*tt1
                        eqd(iat)=eqd(iat)+tt2
                        tt3=(0.5d0*gwsq(iat)*qat(iat)*tt2)
                        do l=1,3
                            do n=1,3
                                den_dginv(n,l)=den_dginv(n,l)+tt3*dk2_dginv(n,l)
                            enddo
                        enddo
                    enddo
                endif
            enddo
        enddo
    enddo
    ehartree=sum_en*2.d0*pi/vol
    den_dginv=den_dginv*(-2.d0*pi/vol)
    dginv_da=0.d0
    do j=1,3
    do i=1,3
        do l=1,3
        do n=1,3
            tt1=-ginv(n,j)*cv(i,1)*ginv(1,l)-ginv(n,1)*cv(i,1)*ginv(j,l) &
                -ginv(n,j)*cv(i,2)*ginv(2,l)-ginv(n,2)*cv(i,2)*ginv(j,l) &
                -ginv(n,j)*cv(i,3)*ginv(3,l)-ginv(n,3)*cv(i,3)*ginv(j,l)
            dginv_da(n,l,i,j)=dginv_da(n,l,i,j)+tt1
        enddo
        enddo
    enddo
    enddo
    tt0=ehartree/vol
    celldv=0.d0
    do j=1,3
    do i=1,3
        tt1=0.d0
        do l=1,3
        do n=1,3
            tt1=tt1+den_dginv(n,l)*dginv_da(n,l,i,j)
        enddo
        enddo
        celldv(i,j)=celldv(i,j)-tt1+tt0*cvinv_vol(i,j)
    enddo
    enddo
    stress(1:3,1:3)=0.d0
    do j=1,3
    do i=1,3
        stress(i,j)=stress(i,j)-(celldv(i,1)*cv(j,1)+celldv(i,2)*cv(j,2)+celldv(i,3)*cv(j,3))
    enddo
    enddo
    do j=1,3
    do i=1,3
        stress(i,j)=stress(i,j)*(-1.d0/vol)
    enddo
    enddo
    do iat=1,nat
        fat(1,iat)=fat(1,iat)*8.d0*pisq/vol
        fat(2,iat)=fat(2,iat)*8.d0*pisq/vol
        fat(3,iat)=fat(3,iat)*8.d0*pisq/vol
        eqd(iat)=eqd(iat)*4.d0*pi/vol
    enddo
    !dpx=0.d0 ; dpy=0.d0 ; dpz=0.d0
    !do iat=1,nat
    !    dpx=dpx+qat(iat)*rat(1,iat)
    !    dpy=dpy+qat(iat)*rat(2,iat)
    !    dpz=dpz+qat(iat)*rat(3,iat)
    !enddo
    !write(*,'(a,3f10.5)') 'dipole ',dpx,dpy,dpz
    !ehartree=ehartree+twopi*(dpx**2+dpy**2+dpz**2)/(3.d0*vol)
    !do iat=1,nat
    !    fat(1,iat)=fat(1,iat)-2.d0*qat(iat)*dpx*twopi/(3.d0*vol)
    !    fat(2,iat)=fat(2,iat)-2.d0*qat(iat)*dpy*twopi/(3.d0*vol)
    !    fat(3,iat)=fat(3,iat)-2.d0*qat(iat)*dpz*twopi/(3.d0*vol)
    !    write(21,'(2i6,es14.5)') icall,iat,4.d0*pi/(3.d0*vol)*(rat(1,iat)*dpx+rat(2,iat)*dpy+rat(3,iat)*dpz)/eqd(iat)
    !    eqd(iat)=eqd(iat)+4.d0*pi/(3.d0*vol)*(rat(1,iat)*dpx+rat(2,iat)*dpy+rat(3,iat)*dpz)
    !enddo
    deallocate(tarr1)
    deallocate(tarr2)
    deallocate(tarr3)
end subroutine get_psolver_fourier_various
!*****************************************************************************************
subroutine get_psolver_fourier_identical(iverbose,nat,rat,ratred,qat,cv,alphasq,ecut,ehartree,fat,eqd,stress,celldv)
    implicit none
    integer, intent(in):: iverbose, nat
    real(8), intent(in):: rat(3,nat), qat(nat)
    real(8), intent(in):: cv(3,3), alphasq, ecut
    real(8), intent(out):: ratred(3,nat), fat(3,nat), eqd(nat), ehartree, stress(3,3), celldv(3,3)
    !local variables
    integer:: m1, m2, m3, m1_max, m2_max, m3_max, mtot, msq, msq_max
    integer:: iat, i, j, k, n, l
    real(8):: kx, ky, kz, vksq, vksqinv
    real(8):: kxp1, kyp1, kzp1, kxp12, kyp12, kzp12
    real(8):: sum_en, factor, tt
    real(8):: vol, volinv, pi, pisq, sf, bnrm1, bnrm2, bnrm3, tarr3
    real(8):: tt0, tt1, tt2, tt3, ttc, tts, twopi, fourpisq, dpx, dpy, dpz
    real(8):: recvec(3,3), g(3,3), cvinv(3,3), cvinv_vol(3,3)
    real(8):: den_dginv(3,3), dk2_dginv(3,3), dginv_da(3,3,3,3), ginv(3,3)
    real(8):: alphasq4th, alphasq2nd
    real(8), allocatable:: tarr1(:), tarr2(:)
    allocate(tarr1(1:nat))
    allocate(tarr2(1:nat))
    pi=4.d0*atan(1.d0)
    pisq=pi**2
    twopi=8.d0*atan(1.d0)
    fourpisq=4.d0*pisq
    msq_max=10000000
    !-------- volume of unitcell -----------------------------------
    call getvol_alborz(cv,vol)
    volinv = 1.d0/vol
    !-------- atomic positions in reduced coordinates --------------
    call rxyz_cart2int_alborz(nat,cv,rat,ratred)
    !-------- reciprocal lattice vectors  --------------------------
    call mycross(cv(1,2),cv(1,3),recvec(1,1))
    call mycross(cv(1,3),cv(1,1),recvec(1,2))
    call mycross(cv(1,1),cv(1,2),recvec(1,3))
    recvec(1:3,1:3)=recvec(1:3,1:3)*volinv
    !write(*,'(a)') 'reciprocal lattice vectors'
    !write(*,'(3es14.5)') recvec(1,1),recvec(1,2),recvec(1,3)
    !write(*,'(3es14.5)') recvec(2,1),recvec(2,2),recvec(2,3)
    !write(*,'(3es14.5)') recvec(3,1),recvec(3,2),recvec(3,3)
    !--------- setting matrix g_inverse-----------------------------
    do i=1,3
        do j=1,3
            g(i,j)=cv(1,i)*cv(1,j)+cv(2,i)*cv(2,j)+cv(3,i)*cv(3,j)
        enddo
    enddo
    call invertmat_alborz(g,ginv)
    call invertmat_alborz(cv,cvinv)
    do i=1,3
    do j=1,3
        cvinv_vol(i,j)=vol*cvinv(j,i)
    enddo
    enddo
    !write(*,*) ginv(1:3,1:3)
    !-------- setting maximum k-vector ----------------------------
    bnrm1=sqrt(recvec(1,1)**2+recvec(2,1)**2+recvec(3,1)**2)
    bnrm2=sqrt(recvec(1,2)**2+recvec(2,2)**2+recvec(3,2)**2)
    bnrm3=sqrt(recvec(1,3)**2+recvec(2,3)**2+recvec(3,3)**2)
    m1_max=ceiling(sqrt(2.d0*ecut)/(bnrm1*2.d0*pi))
    m2_max=ceiling(sqrt(2.d0*ecut)/(bnrm2*2.d0*pi))
    m3_max=ceiling(sqrt(2.d0*ecut)/(bnrm3*2.d0*pi))
    if(iverbose>=3) then
        write(*,'(a,3i4)') 'm1_max,m2_max,m3_max',m1_max,m2_max,m3_max
    endif
    
    mtot=0
    sum_en=0.d0
    fat(1:3,1:nat)=0.d0
    eqd(1:nat)=0.d0
    den_dginv=0.d0

    alphasq4th=alphasq*0.25d0
    alphasq2nd=alphasq*0.5d0
    do m1=0,m1_max
        if(m1==0) then
            factor=1.d0
        else
            factor=2.d0
        endif
        kxp1=real(m1,8)*recvec(1,1)
        kyp1=real(m1,8)*recvec(2,1)
        kzp1=real(m1,8)*recvec(3,1)
        dk2_dginv(1,1)=twopi**2*m1**2
        do m2=-m2_max,m2_max
            kxp12=kxp1+real(m2,8)*recvec(1,2)
            kyp12=kyp1+real(m2,8)*recvec(2,2)
            kzp12=kzp1+real(m2,8)*recvec(3,2)
            dk2_dginv(1,2)=twopi**2*m1*m2
            dk2_dginv(2,1)=dk2_dginv(1,2)
            dk2_dginv(2,2)=twopi**2*m2**2
            do m3=-m3_max,m3_max
                kx=kxp12+real(m3,8)*recvec(1,3)
                ky=kyp12+real(m3,8)*recvec(2,3)
                kz=kzp12+real(m3,8)*recvec(3,3)
                dk2_dginv(1,3)=twopi**2*m1*m3
                dk2_dginv(3,1)=dk2_dginv(1,3)
                dk2_dginv(2,3)=twopi**2*m2*m3
                dk2_dginv(3,2)=dk2_dginv(2,3)
                dk2_dginv(3,3)=twopi**2*m3**2
                !if(m1==0 .and. m2==5 .and. m3==5) then
                !    write(*,'(5es14.5)') kx,ky,kz,real(m2,8)*recvec(2,2),real(m3,8)*recvec(3,3)
                !    stop
                !endif
                vksq=(kx*kx+ky*ky+kz*kz)*fourpisq
                msq=m1*m1+m2*m2+m3*m3
                if((vksq<=2.0*ecut) .and. msq/=0) then
                    mtot=mtot+1
                    !if(mtot>m_max) exit 
                    vksqinv=1.d0/vksq
                    !calculating the structure factor and its derivative
                    do iat=1,nat
                        tt1=(m1*ratred(1,iat)+m2*ratred(2,iat)+m3*ratred(3,iat))*twopi
                        tarr1(iat)=cos(tt1)
                        tarr2(iat)=sin(tt1)
                    enddo
                    tarr3=exp(-alphasq4th*vksq)
                    ttc=0.d0
                    tts=0.d0
                    do iat=1,nat
                        ttc=ttc+qat(iat)*tarr1(iat)*tarr3
                        tts=tts+qat(iat)*tarr2(iat)*tarr3
                    enddo
                    sf=(ttc**2+tts**2)
                    !write(*,'(3i4,i7,4es14.5,es19.10)') m1,m2,m3,mtot,kx,ky,kz,sf,sum_en
                    tt1=factor*sf*vksqinv
                    sum_en=sum_en+tt1
                    tt2=tt1*vksqinv
                    do l=1,3
                        do n=1,3
                            den_dginv(n,l)=den_dginv(n,l)+tt2*dk2_dginv(n,l)
                        enddo
                    enddo
                    tt0=factor*tarr3*vksqinv
                    do iat=1,nat
                        tt1=(tarr2(iat)*ttc-tarr1(iat)*tts)*tt0*qat(iat)
                        tt2=(tarr1(iat)*ttc+tarr2(iat)*tts)*tt0
                        fat(1,iat)=fat(1,iat)+kx*tt1
                        fat(2,iat)=fat(2,iat)+ky*tt1
                        fat(3,iat)=fat(3,iat)+kz*tt1
                        eqd(iat)=eqd(iat)+tt2
                        tt3=(alphasq2nd*qat(iat)*tt2)
                        do l=1,3
                            do n=1,3
                                den_dginv(n,l)=den_dginv(n,l)+tt3*dk2_dginv(n,l)
                            enddo
                        enddo
                    enddo
                endif
            enddo
        enddo
    enddo
    ehartree=sum_en*twopi*volinv
    den_dginv=den_dginv*(-twopi*volinv)
    dginv_da=0.d0
    do j=1,3
    do i=1,3
        do l=1,3
        do n=1,3
            tt1=-ginv(n,j)*cv(i,1)*ginv(1,l)-ginv(n,1)*cv(i,1)*ginv(j,l) &
                -ginv(n,j)*cv(i,2)*ginv(2,l)-ginv(n,2)*cv(i,2)*ginv(j,l) &
                -ginv(n,j)*cv(i,3)*ginv(3,l)-ginv(n,3)*cv(i,3)*ginv(j,l)
            dginv_da(n,l,i,j)=dginv_da(n,l,i,j)+tt1
        enddo
        enddo
    enddo
    enddo
    tt0=ehartree*volinv
    celldv=0.d0
    do j=1,3
    do i=1,3
        tt1=0.d0
        do l=1,3
        do n=1,3
            tt1=tt1+den_dginv(n,l)*dginv_da(n,l,i,j)
        enddo
        enddo
        celldv(i,j)=celldv(i,j)-tt1+tt0*cvinv_vol(i,j)
    enddo
    enddo
    stress(1:3,1:3)=0.d0
    do j=1,3
    do i=1,3
        stress(i,j)=stress(i,j)-(celldv(i,1)*cv(j,1)+celldv(i,2)*cv(j,2)+celldv(i,3)*cv(j,3))
    enddo
    enddo
    do j=1,3
    do i=1,3
        stress(i,j)=-stress(i,j)*volinv
    enddo
    enddo
    tt = 8.d0*pisq*volinv
    do iat=1,nat
        fat(1,iat)=fat(1,iat)*tt
        fat(2,iat)=fat(2,iat)*tt
        fat(3,iat)=fat(3,iat)*tt
        eqd(iat)=eqd(iat)*4.d0*pi*volinv
    enddo
    !dpx=0.d0 ; dpy=0.d0 ; dpz=0.d0
    !do iat=1,nat
    !    dpx=dpx+qat(iat)*rat(1,iat)
    !    dpy=dpy+qat(iat)*rat(2,iat)
    !    dpz=dpz+qat(iat)*rat(3,iat)
    !enddo
    !write(*,'(a,3f10.5)') 'dipole ',dpx,dpy,dpz
    !ehartree=ehartree+twopi*(dpx**2+dpy**2+dpz**2)/(3.d0*vol)
    !do iat=1,nat
    !    fat(1,iat)=fat(1,iat)-2.d0*qat(iat)*dpx*twopi/(3.d0*vol)
    !    fat(2,iat)=fat(2,iat)-2.d0*qat(iat)*dpy*twopi/(3.d0*vol)
    !    fat(3,iat)=fat(3,iat)-2.d0*qat(iat)*dpz*twopi/(3.d0*vol)
    !enddo
    deallocate(tarr1)
    deallocate(tarr2)
end subroutine get_psolver_fourier_identical
!*****************************************************************************************
