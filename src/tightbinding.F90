!*****************************************************************************************
subroutine set_indorb(partb,atoms)
    use mod_atoms, only: typ_atoms
    use mod_tightbinding, only: typ_partb
    implicit none
    type(typ_partb), intent(inout):: partb
    type(typ_atoms), intent(inout):: atoms
    !local variables
    integer:: iat, iorb, mel, nel
    character(5):: sat
    if(atoms%nat==0) stop 'ERROR: nat==0 in set_indorb'
    partb%indat(1:atoms%nat)=0
    iat=1
    mel=0
    do iorb=1,partb%norb
        if(iat>atoms%nat) stop 'ERROR: iat>nat in set_indorb'
        partb%indorb(iorb)=iat
        mel=mel+1
        sat=atoms%stypat(atoms%itypat(iat))
        if(trim(sat)=='Si') then
            nel=4
        else if(trim(sat)=='H') then
            nel=1
        else if(trim(sat)=='C') then
            nel=4
        else
            stop 'ERROR: unknown symbol of atom in set_indorb'
        endif
        if(partb%indat(iat)==0) then
            partb%indat(iat)=iorb
            partb%norbat(iat)=nel
        endif
        !write(71,'(4i5)') iorb,iat,mel,nel
        if(mel==nel) then
            iat=iat+1
            mel=0
        endif
    enddo
end subroutine set_indorb
!*****************************************************************************************
!Returns the band structure contribution to the energy/atom at the gamma point
!for nat atoms with orthogonal positions in cell given by coord,
!for a coordinate system with lattice vectors latt.  Nonorthogonal case is assumed.
!ADDS to the forces in array force[], which are already present from pair term.
subroutine gammaenergy(pia_arr,linked_lists,partb,atoms,natsi,pplocal)
    use mod_linked_lists, only: typ_pia_arr, typ_linked_lists
    use mod_atoms, only: typ_atoms, set_typat
    use mod_potl, only: potl_typ
    use mod_tightbinding, only: typ_partb, lenosky
    use mod_const, only: ha2ev
    use dynamic_memory
    implicit none
    type(typ_linked_lists), intent(in):: linked_lists
    type(typ_pia_arr), intent(in):: pia_arr
    type(typ_partb), intent(inout):: partb
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in):: natsi
    type(potl_typ), intent(inout):: pplocal
    !local variables
    real(8):: tt, reg1, sterm, totaln, beta, term1, term2
    real(8), allocatable:: ggocc(:)
    real(8), allocatable:: rho(:,:)
    integer:: iorb, jorb, korb, n1, n2, iat, jat, ixyz
    !call f_routine(id='gammaenergy')
    call set_typat(atoms) !exists in basic_atoms.F90 file
    call set_indorb(partb,atoms)
    !rho=f_malloc([1.to.partb%norb,1.to.partb%norb],id='rho')
    !ggocc=f_malloc([1.to.partb%norb],id='ggocc')
    allocate(rho(partb%norb,partb%norb))
    allocate(ggocc(partb%norb))
    !Build the TB Hamiltonian and diagonalize it to obtain eigenvalues/vectors
    call gammamat(pia_arr,linked_lists,partb,atoms,natsi,0,pplocal) 
    !do iorb=1,partb%norb
    !        write(8,'(i3,8es15.6)') iorb,partb%tbmat(iorb,1:8)
    !enddo
    !stop 'UUUUUUUUUUUUU'
    !write(*,*) partb%tbmat(1:8,1:8)
    !stop
    call forcediagonalizeg(partb)
    !Fermi dirac distribution: returns total energy :: E_TB = sum_n f_n*e_n
    !Create density matrix rho to obtain forces using Hellmann-Feynman theorem:
    !F_k = sum_ij rho_(ij)*dH_(ij)/dr_k
    do iorb=1,partb%norb
        do jorb=1,partb%norb
            rho(jorb,iorb)=0.d0
        enddo
    enddo
    call yfdocclocal(partb)
    !do iorb=1,partb%norbcut
    !    write(53,'(i4,es14.5)') iorb,partb%eval(iorb)
    !enddo
    !stop 'YYYYYYYYYYYYYYYYYYYYY'
    !write(61,'(6es14.5,i5)') partb%eval(60),partb%eval(61),partb%eval(61)-partb%eval(60),partb%focc(60),partb%focc(61),partb%focc(partb%norbcut),partb%norbcut
    !write(*,'(a,8f10.3)') 'fermi',partb%focc(1),partb%focc(2),partb%focc(3),partb%focc(4),partb%focc(5),partb%focc(6),partb%focc(7),partb%focc(8)
    !write(*,'(a,10es14.5)') 'eval ',parb%eband,partb%eval(1),partb%eval(2),partb%eval(3),partb%eval(4),partb%focc(1),partb%focc(2),partb%focc(3),partb%focc(4),partb%focc(5)
    !Need to include extra terms.  We have focc[iorb-1]...focc[NC(stride*nat)-1]
    !also fermitemp=temp in Kelvin so "beta"=1 / kT=11604 / fermitemp
    !ggocc maximum is 2.0
    term1=0.d0
    term2=0.d0
    beta=11604.d0*ha2ev/partb%temp_fermi
    do iorb=1,partb%norbcut
        term1=term1-partb%focc(iorb)*(1.d0-partb%focc(iorb)*0.5d0)
        term2=term2-partb%eval(iorb)*partb%focc(iorb)*(1.d0-partb%focc(iorb)*0.5d0)
    enddo
    do iorb=1,partb%norbcut
        ggocc(iorb)=partb%focc(iorb)-beta*partb%focc(iorb)*(1.d0-partb%focc(iorb)*0.5d0)*(partb%eval(iorb)-1.d0*term2/term1)
    enddo
    !totaln=0.d0
    !do iorb=1,partb%norbcut
    !    totaln=totaln + partb%focc(iorb) !Total number of electrons
    !enddo
    !write(*,*) 'Total number of electrons ',totaln
    !loop was over occupied only
    do iorb=1,partb%norbcut
        !For ith eigenvector partb%evec[iorb][1..stride*nat], compute rho 
        !Use factors (focc[iorb] / 2.0) for partb%evec[iorb][]  
        !Rho is quadratic in eigenfunctions, so contains this factor (probability) once
        sterm=2.d0
        do jorb=1,partb%norb
            reg1=partb%evec(jorb,iorb)*sterm*ggocc(iorb)/2.d0
            do korb=jorb+1,partb%norb
                !write(*,'(a,3i4,2es14.5)') 'RHO ',iorb,jorb,korb,reg1,partb%evec(korb,iorb)
                rho(korb,jorb)=rho(korb,jorb)+reg1*partb%evec(korb,iorb)
            enddo
        enddo
    enddo
    !stop 'ZZZZZZZZZZZZZZ'
    !Take trace rho * H for Hellmann-Feynman theorem
    !Note that other triangle of rho, and we do sum with
    !factors of two to compensate.
    if(lenosky .or. trim(partb%event)/='train' )then 
        do ixyz=1,3
            call gammamat(pia_arr,linked_lists,partb,atoms,natsi,ixyz,pplocal)
            do iorb=1,partb%norb
                iat=partb%indorb(iorb)
                do jorb=iorb+1,partb%norb
                    jat=partb%indorb(jorb)
                    tt=rho(jorb,iorb)*partb%tbmat(jorb,iorb)*2.d0
                    atoms%fat(ixyz,iat)=atoms%fat(ixyz,iat)-tt
                    atoms%fat(ixyz,jat)=atoms%fat(ixyz,jat)+tt
                enddo
            enddo
            !write(*,*) "FAt", atoms%fat(1,iat)  
        enddo
    else if(trim(partb%event)=='train')then
        do ixyz=1,4
            call gammamat(pia_arr,linked_lists,partb,atoms,natsi,ixyz,pplocal)
            do iorb=1,partb%norb
                iat=partb%indorb(iorb)
                do jorb=iorb+1,partb%norb
                    jat=partb%indorb(jorb)
                    tt=rho(jorb,iorb)*partb%tbmat(jorb,iorb)*2.d0
                    !write(*,'(a,5i4,2es14.5)') 'dedh ',ixyz,iorb,jorb,iat,jat,rho(jorb,iorb),partb%tbmat(jorb,iorb)
                    partb%dedh(ixyz,iat,jat)=partb%dedh(ixyz,iat,jat)+tt
                enddo
            enddo
        enddo
        !stop 'FFFFFFFFFFFFFFF'
    endif
    !call f_free(rho)
    !call f_free(ggocc)
    deallocate(rho)
    deallocate(ggocc)
    !call f_release_routine()
end subroutine gammaenergy
!*****************************************************************************************
!Within the unit cell, this routine
!sets up the entire matrix to be diagonalized at the gamma point.
!Lattice translation vectors are also given in orthogonal coordinates.
!Entire coupling matrix is returned in partb%tbmat, must be
!dimensioned properly on call.
!Note hgen and dhgen are arrays used to store values of splines and
!their derivatives for each atom pair.
subroutine gammamat(pia_arr,linked_lists,partb,atoms,natsi,flag2,pplocal)
    use mod_linked_lists, only: typ_pia_arr, typ_linked_lists
    use mod_tightbinding, only: typ_partb
    use mod_atoms, only: typ_atoms
    use mod_potl, only: potl_typ
    use dynamic_memory
    implicit none
    type(typ_linked_lists), intent(in):: linked_lists
    type(typ_pia_arr), intent(in):: pia_arr
    type(typ_partb), intent(inout):: partb
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in):: natsi, flag2
    type(potl_typ), intent(in):: pplocal
    !local variables
    integer:: iat, jat, iorb, jorb, ib
    integer:: indexi, indexj, atomtypei, atomtypej, norbi, norbj
    real(8), allocatable:: rex(:,:)
    real(8) eself(4), es
    !call f_routine(id='gammamat')
    !rex=f_malloc([1.to.partb%nstride,1.to.partb%nstride],id='rex')
    allocate(rex(partb%nstride,partb%nstride))
    do iorb=1,partb%norb
        do jorb=1,partb%norb
            partb%tbmat(jorb,iorb)=0.d0
        enddo
    enddo
    !Include all atoms by loop over iat and jat 
    do ib=1,linked_lists%maxbound_rad
        iat=linked_lists%bound_rad(1,ib)
        jat=linked_lists%bound_rad(2,ib)
        atomtypei=0
        if(iat>natsi) atomtypei=1
        norbi=partb%norbat(iat) !number of orbitals of atom type iat
        indexi=partb%indat(iat)-1
        atomtypej=0
        if(jat>natsi) atomtypej=1
        norbj=partb%norbat(jat) !number of orbitals of atom type jat
        indexj=partb%indat(jat)-1
        !Obtain matrix of coupling between atom iat and atom jat (iat/=jat)
        call gammacoupling(pia_arr%pia(ib),ib,partb,atoms,flag2,atomtypei,atomtypej,pplocal,rex)
        do iorb=1,norbi
            do jorb=1,norbj
                !off-diagonal terms of H_TB is constructed
                partb%tbmat(indexj+jorb,indexi+iorb)=partb%tbmat(indexj+jorb,indexi+iorb)+rex(jorb,iorb)
            enddo
        enddo
    enddo
    !Compelet matrix of coupling by including on-site energies(eself for Si and es for H)
    !eself has 4 values, then it seems diagonal terms of H_TB also themselves
    !are diagonal.
    call eselfgeneral(eself)
    do iorb=1,partb%nstride*natsi
        jorb=mod(iorb-1,partb%nstride)+1
        partb%tbmat(iorb,iorb)=partb%tbmat(iorb,iorb)+eself(jorb)
    enddo
    es=pplocal%eps(1) 
    do iorb=partb%nstride*natsi+1,partb%nstride*natsi+atoms%nat-natsi
        partb%tbmat(iorb,iorb)=partb%tbmat(iorb,iorb)+es
    enddo
    !call f_free(rex)
    deallocate(rex)
    !call f_release_routine()
end subroutine gammamat
!*****************************************************************************************
!Does the eigenvalue problem, with arguments having the same format as the
!routine diagonalize.  m is H and is real, nc by nc
!matrix.  Eigen is the list of eval, and evec is the list of eigenvectors.
subroutine forcediagonalizeg(partb)
    use mod_tightbinding, only: typ_partb
    use dynamic_memory
    implicit none
    type(typ_partb), intent(inout):: partb
    !local variables
    integer:: n, i, j, ierr, lwork, liwork, m, nc
    real(8), allocatable:: a(:,:)
    real(8), allocatable:: work(:)
    integer, allocatable:: iwork(:)
    integer, allocatable:: isuppz(:)
    real(8):: abstol=0.d0 !The absolute error tolerance for the eigenvalues.
    integer, save:: errcount=0
    integer, save:: icall=0
    !real(8):: w1(1000)
    !real(8):: w2(1000)
    icall=icall+1
    !call f_routine(id='forcediagonalizeg')
    n=partb%norb
    nc=partb%norbcut
    lwork=n*n+100*n
    liwork=n*n+50*n
    !isuppz=f_malloc([1.to.2*n],id='isuppz')
    !a=f_malloc([1.to.n,1.to.n],id='a')
    !work=f_malloc([1.to.lwork],id='work')
    !iwork=f_malloc([1.to.liwork],id='iwork')
    allocate(isuppz(2*n))
    allocate(a(n,n))
    allocate(work(lwork))
    allocate(iwork(liwork))
    !partb%tbmat=1.d-10
    !do i=1,n
    !    partb%tbmat(i,i)=1.d0
    !enddo
    do i=1,n
        do j=1,n
            a(i,j)=partb%tbmat(i,j)
        enddo
    enddo
    !ierr does not need to be set on entry work is the workspace array,
    !already set eval is also an output
    !write(*,'(a,7i6)') 'EEEEEEEEEEEEEE ',icall,errcount,size(partb%eval),size(partb%tbmat),size(a),lwork,liwork
    !if(icall==7 .or. icall==8) then
    !    write(*,'(8es14.5)') (partb%tbmat(1,j),j=1,8)
    !    write(*,'(8es14.5)') (partb%tbmat(2,j),j=1,8)
    !    write(*,'(8es14.5)') (partb%tbmat(3,j),j=1,8)
    !    write(*,'(8es14.5)') (partb%tbmat(4,j),j=1,8)
    !    write(*,'(8es14.5)') (partb%tbmat(5,j),j=1,8)
    !    write(*,'(8es14.5)') (partb%tbmat(6,j),j=1,8)
    !    write(*,'(8es14.5)') (partb%tbmat(7,j),j=1,8)
    !    write(*,'(8es14.5)') (partb%tbmat(8,j),j=1,8)
    !    !stop
    !endif
    !write(*,*)
    !write(*,*)
    !write(*,*) size(partb%eval)
    !write(*,'(a,6i5)') 'FFFFF ',icall,n,size(a),nc,size(partb%eval),size(partb%evec)
!    if(nc==n) then
        call DSYEV('V','L',n,a,n,partb%eval,work,lwork,ierr)
        do i=1,n
            do j=1,n
                partb%evec(i,j)=a(i,j)
            enddo
        enddo
!    else
!        call dsyevr('V','I','L',n,a,n,0.d0,0.d0,1,nc,abstol,m,partb%eval,partb%evec,n, &
!             isuppz,work,lwork,iwork,liwork,ierr)
!    endif
    !write(*,'(8es14.5)') (partb%evec(1,j),j=1,8)
    !write(*,'(8es14.5)') (partb%evec(2,j),j=1,8)
    !write(*,'(8es14.5)') (partb%evec(3,j),j=1,8)
    !write(*,'(8es14.5)') (partb%evec(4,j),j=1,8)
    !write(*,'(8es14.5)') (partb%evec(5,j),j=1,8)
    !write(*,'(8es14.5)') (partb%evec(6,j),j=1,8)
    !write(*,'(8es14.5)') (partb%evec(7,j),j=1,8)
    !write(*,'(8es14.5)') (partb%evec(8,j),j=1,8)
    !write(*,'(a,7f10.3)') 'FFFFF ',partb%eval(1),partb%eval(2),partb%eval(3),partb%eval(4),partb%eval(5),partb%eval(6),partb%eval(5)-partb%eval(4)
    if(ierr/= 0  .and. errcount < 250) then
        write(*,'(a36,2i8,a23)') 'TBNORTH WARNING: ierr , errcount == ',ierr,errcount,' from dsygv diagonalize'
        write(*,'(a)') 'Will not print this message when errcount exceeds 250'
        errcount=errcount+1
    endif
    !call f_free(isuppz)
    !call f_free(a)
    !call f_free(work)
    !call f_free(iwork)
    deallocate(isuppz)
    deallocate(a)
    deallocate(work)
    deallocate(iwork)
    !call f_release_routine()
end subroutine forcediagonalizeg
!*****************************************************************************************
!README: There are too many arguments setting intents of variables postponed
!after putting some variables in new derived types.
subroutine gammacoupling(pia,ib,partb,atoms,flag2,atomtypei,atomtypej,pplocal,rem)
    use mod_linked_lists, only: typ_pia !, typ_linked_lists
    use mod_tightbinding, only: typ_partb, lenosky
    use mod_atoms, only: typ_atoms
    use mod_potl, only: potl_typ
    use mod_const, only: bohr2ang, ha2ev
    use dynamic_memory
    implicit none
    !type(typ_linked_lists), intent(in):: linked_lists
    type(typ_pia), intent(in):: pia
    type(typ_partb), intent(inout):: partb
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in):: ib, flag2, atomtypei, atomtypej
    type(potl_typ), intent(in):: pplocal
    real(8), intent(out):: rem(partb%nstride,partb%nstride)
    !local variables
    real(8):: diff(3), dist
    integer:: iorb, jorb, ix, iy, iz !, iat, jat
    real(8), allocatable:: hgen(:), dhgen(:)
    !call f_routine(id='gammacoupling')
    !hgen=f_malloc([1.to.4],id='hgen')
    !dhgen=f_malloc([1.to.4],id='dhgen')
    allocate(hgen(4))
    allocate(dhgen(4))
    do jorb=1,partb%nstride
        do iorb=1,partb%nstride
            rem(iorb,jorb)=0.d0
        enddo
    enddo
    !Here want section to compute ix,iy,iz for nearest image, rather than
    !running through choices
    diff(1)=-pia%dr(1)
    diff(2)=-pia%dr(2)
    diff(3)=-pia%dr(3)
    !The distance between two arbitrary atoms
    dist=pia%r
    !Intractions are considered over a fixed domain (r_low=0.0001,r_up=paircut).
    !write(*,*) 'HERE if ',iat,jat,dist,partb%paircut
    if(dist>=1.d-4 .and. .not. dist>partb%paircut) then
        !Compute unit vector between atoms as input variable in slatercoupling()
        !write(*,'(a,3es19.10)') 'diff ',diff(1),diff(2),diff(3)
        diff(1)=diff(1)/dist
        diff(2)=diff(2)/dist
        diff(3)=diff(3)/dist
        !write(*,'(a,3es19.10)') 'diff ',diff(1),diff(2),diff(3)
        !stop
        if(flag2==0) then
            !Returns h(r)s and their derivatives using cubic splines in hgen and dhgen.
            if(lenosky)then
            call radelmgeneralsp(dist,hgen,dhgen,atomtypei,atomtypej,pplocal)
            partb%hgenall0(ib)=hgen(1)
            partb%hgenall1(ib)=hgen(2)
            partb%hgenall2(ib)=hgen(3)
            partb%hgenall3(ib)=hgen(4)
            partb%dhgenall0(ib)=dhgen(1)
            partb%dhgenall1(ib)=dhgen(2)
            partb%dhgenall2(ib)=dhgen(3)
            partb%dhgenall3(ib)=dhgen(4)
            write(44,'(a,5e26.17)') 'hgen-L',hgen(1),hgen(2),hgen(3),hgen(4),dist
            endif
        endif
        hgen(1)=partb%hgenall0(ib)
        hgen(2)=partb%hgenall1(ib)
        hgen(3)=partb%hgenall2(ib)
        hgen(4)=partb%hgenall3(ib)
        !write(55,'(a,5es14.5)') 'hgen_nn ',hgen(1),hgen(2),hgen(3),hgen(4),dist
        dhgen(1)=partb%dhgenall0(ib)
        dhgen(2)=partb%dhgenall1(ib)
        dhgen(3)=partb%dhgenall2(ib)
        dhgen(4)=partb%dhgenall3(ib)
        !Returns rem (matrix of coupling) 
        if(trim(partb%event)/='train' .or. lenosky .or. flag2==0) then
            call slatercoupling(diff,dist,hgen,dhgen,flag2,rem)
        endif
        if(flag2>0 .and. trim(partb%event)=='train') then
            call Hamiltonian_der(diff,flag2,rem)
        endif
    endif
    !call f_free(hgen)
    !call f_free(dhgen)
    deallocate(hgen)
    deallocate(dhgen)
    !call f_release_routine()
end subroutine gammacoupling
!*****************************************************************************************
!Generates the matrix of couplings for the program, if u is unit vector between atoms
!r is distance, returns coupling matrix mat.
!eselfgeneral must be treated as a seperate case. Flag is set to zero for H, 1 for S
!flag2 is 0 for normal case, 1 for x derivative, 2 for y, 3 for z, where
!u = (x,y,z) / sqrt(x^2+y^2+z^2)...
subroutine slatercoupling(u,r,hgen,dhgen,flag2,mat)
    implicit none
    real(8), intent(in):: r, u(3)
    real(8), intent(in):: hgen(4), dhgen(4)
    integer, intent(in):: flag2
    real(8), intent(out):: mat(4,4)
    !local variables
    real(8):: drdx, g(3), fv(4)
    real(8):: ess, esx, esy, esz, exx, eyy, ezz, exy, eyz, exz
    real(8):: hsss, hsps, hpps, hppp
    real(8):: dhsss, dhsps, dhpps, dhppp
    hsss=hgen(1)
    hsps=hgen(2)
    hpps=hgen(3)
    hppp=hgen(4)
    if(flag2/=0) then
        dhsss=dhgen(1)
        dhsps=dhgen(2)
        dhpps=dhgen(3)
        dhppp=dhgen(4)
    endif
    if(flag2==0) then
        ess=hsss
        esx=hsps*u(1)
        esy=hsps*u(2)
        esz=hsps*u(3)
        exx=hpps*u(1)*u(1)+(1.d0-u(1)*u(1))*hppp
        eyy=hpps*u(2)*u(2)+(1.d0-u(2)*u(2))*hppp
        ezz=hpps*u(3)*u(3)+(1.d0-u(3)*u(3))*hppp
        exy=(u(1)*u(2))*(hpps-hppp)
        eyz=(u(2)*u(3))*(hpps-hppp)
        exz=(u(1)*u(3))*(hpps-hppp)
    endif
    flag2_nonzero: if(flag2/=0) then
        !Use dr/dx=x/r=u_x
        !Assume fv (type position) is unit vector along x,y, or z depending on flag2
        !flag2/=0 is used to abtain forces using derivatives of energies.
        fv(1)=0.d0 ; fv(2)=0.d0 ; fv(3)=0.d0
        fv(flag2)=1.d0
        drdx=u(1)*fv(1)+u(2)*fv(2)+u(3)*fv(3)
        g(1)=fv(1)/r-drdx*u(1)/r  !derivative of u.x in direction of fv
        g(2)=fv(2)/r-drdx*u(2)/r
        g(3)=fv(3)/r-drdx*u(3)/r
        !Calculate derivatives with respect to r dependence of matrix elements
        ess=drdx*dhsss
        esx=drdx*dhsps*u(1)
        esy=drdx*dhsps*u(2)
        esz=drdx*dhsps*u(3)
        exx=drdx*(dhpps*u(1)*u(1)+(1-u(1)*u(1))*dhppp)
        eyy=drdx*(dhpps*u(2)*u(2)+(1-u(2)*u(2))*dhppp)
        ezz=drdx*(dhpps*u(3)*u(3)+(1-u(3)*u(3))*dhppp)
        exy=drdx* (u(1)*u(2))*(dhpps-dhppp)
        eyz=drdx* (u(2)*u(3))*(dhpps-dhppp)
        exz=drdx* (u(1)*u(3))*(dhpps-dhppp)
        !Now calculate derivatives with respect to x,y,z appearing in u=(x,y,z)/Sqrt(...)
        !Here g is the derivative of u with respect to either x,y,z (for flag2=1, 2, or 3)
        esx=esx+hsps*g(1)
        esy=esy+hsps*g(2)
        esz=esz+hsps*g(3)
        exx=exx+2.d0*hpps*g(1)*u(1)+(-2.d0*g(1)*u(1))*hppp
        eyy=eyy+2.d0*hpps*g(2)*u(2)+(-2.d0*g(2)*u(2))*hppp
        ezz=ezz+2.d0*hpps*g(3)*u(3)+(-2.d0*g(3)*u(3))*hppp
        exy=exy+(g(1)*u(2)+u(1)*g(2))*(hpps-hppp)
        eyz=eyz+(g(2)*u(3)+u(2)*g(3))*(hpps-hppp)
        exz=exz+(g(1)*u(3)+u(1)*g(3))*(hpps-hppp)
    endif flag2_nonzero
    mat(1,1)=ess ; mat(2,1)=-esx ; mat(3,1)=-esy ; mat(4,1)=-esz
    mat(1,2)=esx ; mat(2,2)= exx ; mat(3,2)= exy ; mat(4,2)= exz
    mat(1,3)=esy ; mat(2,3)= exy ; mat(3,3)= eyy ; mat(4,3)= eyz
    mat(1,4)=esz ; mat(2,4)= exz ; mat(3,4)= eyz ; mat(4,4)= ezz
end subroutine slatercoupling
!*****************************************************************************************
!Natob is number of eigenvalues,
!eval is array of eigenvalues,
!focc is array of returned eigenvalue occupancies,
!norb is number of electrons (not occupied states, but electrons),
!eband is returned energy,
!temp is temperature in degrees Kelvin,
!eigenvalues must be sorted on entry.
!Number of electrons nel must be even.  <-- Is this true?
subroutine yfdocclocal(partb)
    use mod_tightbinding, only: typ_partb, lenosky
    use mod_const, only: ha2ev
    implicit none
    type(typ_partb), intent(inout):: partb
    !local variables
    integer, parameter:: fdmaxit= 20000 !Maximum number of NR iterations
    real(8), parameter:: fdepsocc=1.d-9 !Allowed error in number of electrons
    real(8):: rsid, sumocc, sumder, delta, efermi, beta, fdexp, fder
    integer:: it, i, nocc
    do i=1,partb%norbcut
        partb%focc(i)=0.d0
    enddo
    nocc=partb%norb/2
    do i=1,nocc
        partb%focc(i)=2.d0
    enddo
    if(nocc*2/=partb%norb) partb%focc(nocc+1)=1.d0
    efermi=partb%eval(nocc)
    if(nocc*2/=partb%norb) efermi=partb%eval(nocc+1)
    beta=11604.d0*ha2ev/partb%temp_fermi !Conversion of temp into Ha
    it=0
    !Each iteration corresponds to steepest descent iteration. fermi-dirac function is minimized
    !which depend on efermi(chemical potential).
    iteration: do
        it=it+1
        if(it>fdmaxit) write(*,'(a)') 'WARNING: iteration count exceeded in yfdocc'
        sumocc=0.d0
        sumder=0.d0
        i_loop: do i=1,partb%norbcut
            delta=partb%eval(i)-efermi
            if(delta>0.d0) then ! upper than fermi energy
                fdexp=exp(-delta*beta)
                partb%focc(i)=2.d0*fdexp/(1.d0+fdexp)
                sumocc=sumocc+partb%focc(i)
                fder=2.d0*beta*fdexp/(1.d0+fdexp)/(1.d0+fdexp) ! fder:=dn/dE
                sumder=sumder+fder
            else ! lower than fermi energy
                fdexp=exp(delta*beta)
                partb%focc(i)=2.d0/(1.d0+fdexp)
                sumocc=sumocc+partb%focc(i)
                fder=2.d0*beta*fdexp/(1.d0+fdexp)/(1.d0+fdexp)
                sumder=sumder+fder
            endif
        enddo i_loop
        rsid=sumocc-real(partb%norb,8) !  rsid:= dn
        delta=(-rsid/sumder) !  delta := dE = dn/(dn/dE)
        if(abs(rsid)>FDEPSOCC) efermi=efermi-rsid/sumder
        if(.not. (abs(rsid)>FDEPSOCC .and. it<=fdmaxit)) exit
    enddo iteration
    partb%eband=0.d0
    do i=1,partb%norbcut
        partb%eband=partb%eband+partb%focc(i)*partb%eval(i)
    enddo
end subroutine yfdocclocal
!*****************************************************************************************
subroutine Hamiltonian_der(u,flag2,mat)
    implicit none
    real(8), intent(in):: u(3)
    integer, intent(in):: flag2
    real(8), intent(out):: mat(4,4)
    !local variables
    real(8):: ess, esx, esy, esz, exx, eyy, ezz, exy, eyz, exz
    !Here flag2 represent the number of "hgen"s.
    if(flag2==1) then
        ess=1
        esx=0
        esy=0
        esz=0
        exx=0
        eyy=0
        ezz=0
        exy=0
        eyz=0
        exz=0
    else if(flag2==2) then
        ess=0
        esx=u(1)
        esy=u(2)
        esz=u(3)
        exx=0
        eyy=0
        ezz=0
        exy=0
        eyz=0
        exz=0
    else if(flag2==3) then
        ess=0
        esx=0
        esy=0
        esz=0
        exx=u(1)*u(1)
        eyy=u(2)*u(2)
        ezz=u(3)*u(3)
        exy=u(1)*u(2)
        eyz=u(2)*u(3)
        exz=u(1)*u(3)
    else if(flag2==4) then
        ess=0
        esx=0
        esy=0
        esz=0
        exx=(1-u(1)*u(1))
        eyy=(1-u(2)*u(2))
        ezz=(1-u(3)*u(3))
        exy=-u(1)*u(2)
        eyz=-u(2)*u(3)
        exz=-u(1)*u(3)
    endif
    mat(1,1)=ess ; mat(2,1)=-esx ; mat(3,1)=-esy ; mat(4,1)=-esz
    mat(1,2)=esx ; mat(2,2)= exx ; mat(3,2)= exy ; mat(4,2)= exz
    mat(1,3)=esy ; mat(2,3)= exy ; mat(3,3)= eyy ; mat(4,3)= eyz
    mat(1,4)=esz ; mat(2,4)= exz ; mat(3,4)= eyz ; mat(4,4)= ezz
end subroutine Hamiltonian_der
