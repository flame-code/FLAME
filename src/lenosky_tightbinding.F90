!*****************************************************************************************
!In tight-binding(TB) method, given a set of atomic positions, the total energy is
!obtained from two terms: E_tot= E_band + E_rep
!First term of above equation(E_band) is obtained by summing lowest eigenvalues of H_TB.
!H_TB is determined by 5 radial functions, here they are represented using splines.
!In this model the 5 radial fuctions are defined over a fixed domain (r_low,r_up).
!E_rep is assumed to depend upon the positions of the atoms: pair potential
!------------------------------------- properties of system ------------------------------
!Here we have a system that contains:
!atomtype 0 = Si
!atomtype 1 = hydrogen
!silicons come first in list 1...natsi
!hydrogens are natsi+1...nat
!hydrogens have 1 electron and 1 orbital.
!silicons have 4 electrons in valence layer, then have 4 orbitals.
!------------------------------------------ Routines -------------------------------------
!lenoskytb reads an input file(contains properties of system) and then calls totalenergy.
!totalenergy calls:
!    1) gammaenergy (band energy)
!    2) pairenergy (pairwise potential)
!Then calculates forces of atoms.
!gammaenergy calls:
!    1) gammamat (for setting up TB matrix)
!    2) forcediagonalizeg, (to diagonalize TB matrix)
!    3) yfdocclocal (fermi-dirac occupancy)
!    4) then gammamat (for derivs of matrix => forces are obtained)
!gammamat calls gammacoupling
!gammacoupling calls slatercoupling
!slatercopling is constructed off-diagonal elements of H_TB, using Slater-Koster table
!*****************************************************************************************
!lenoskytb_() is interface routine by which TB routines may be called from external code
!Note silicon atoms must be first in atoms%rat
subroutine lenoskytb_alborz(parini,atoms,natsi,count_md)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_potl, only: potl_typ
    use mod_tightbinding, only: typ_partb, lenosky
    use mod_linked_lists, only: typ_pia_arr, typ_linked_lists
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in):: natsi
    real(8), intent(inout):: count_md
    !local variables
    integer, save:: firstcall=1
    type(potl_typ), save:: pplocal
    type(typ_partb), save:: partb
    type(typ_parini):: parini_tmp
    type(typ_linked_lists):: linked_lists
    type(typ_pia_arr):: pia_arr
    if(firstcall==1) then
        call yaml_map('LTB GAMMA POINT only','first call')
        !write(*,'(a)') 'GAMMA POINT only tight binding code'
        if(lenosky) then
            call yaml_comment('Reading spline potential coeff.cls') !,hfill='~')
            !write(*,'(a)') 'Reading spline potential coeff.cls'
            call prmst38c(partb,pplocal) !Reads potential 
        endif
        call yaml_map('paircut',partb%paircut,fmt='(es18.10)')
        !write(*,'(a,f)') 'paircut= ',partb%paircut
        firstcall=0
    endif
    linked_lists%rcut=partb%paircut
    linked_lists%triplex=.true.
    parini_tmp=parini !TO_BE_CORRECTED
    parini_tmp%bondbased_ann=.true. !TO_BE_CORRECTED
    call call_linkedlist(parini_tmp,atoms,.true.,linked_lists,pia_arr)
    call lenoskytb_init(partb,atoms,natsi,linked_lists)
    count_md=count_md+1.d0
    !PRINT SOME WARNINGS
    if(natsi>atoms%nat) write(*,'(a)') 'WARNING natsi = ',natsi,' is greater than number of atoms = ',atoms%nat
    if(atoms%nat == 0) write(*,'(a)') 'WARNING lenoskytb_alborz called with zero atoms'
    if(atoms%nat < 0) write(*,'(a)') 'WARNING lenoskytb_alborz called with negative number of atoms'
    call totalenergy(pia_arr,linked_lists,parini,partb,atoms,natsi,pplocal)
    !write(61,*) atoms%epot
    !stop 'BBBBBBBBBBBBBBB'
    call lenoskytb_final(partb)
    deallocate(linked_lists%prime_bound)
    deallocate(linked_lists%bound_rad)
    deallocate(linked_lists%bound_ang)
end subroutine lenoskytb_alborz
!*****************************************************************************************
subroutine lenoskytb_init(partb,atoms,natsi,linked_lists)
    use mod_tightbinding, only: typ_partb, lenosky
    use mod_atoms, only: typ_atoms
    use mod_linked_lists, only: typ_linked_lists
    use dynamic_memory
    implicit none
    type(typ_partb), intent(inout):: partb
    type(typ_atoms), intent(in):: atoms
    integer, intent(in):: natsi
    type(typ_linked_lists), intent(in):: linked_lists
    !local variables
    !call f_routine(id='lenoskytb_init')
    partb%norb=partb%nstride*natsi+(atoms%nat-natsi)
    !write(*,*) partb%nstride, natsi, partb%norb
    !Number of eigenvalues to compute=MIN(matrix*frac+extra,dimension of matrix)
    if(floor(partb%norb*partb%frac+partb%extra)>partb%norb) then
        partb%norbcut=partb%norb
    else
        partb%norbcut=floor((partb%norb*partb%frac+partb%extra))
    endif
    !---------------------------------------------------------
    !partb%indorb=f_malloc([1.to.partb%norb],id='partb%indorb')
    !partb%indat=f_malloc([1.to.atoms%nat],id='partb%indat')
    !partb%norbat=f_malloc([1.to.atoms%nat],id='partb%norbat')
    !partb%tbmat=f_malloc([1.to.partb%norb,1.to.partb%norb],id='partb%tbmat')
    !partb%evec=f_malloc([1.to.partb%norb,1.to.partb%norbcut],id='partb%evec')
    allocate(partb%indorb(partb%norb))
    allocate(partb%indat(atoms%nat))
    allocate(partb%norbat(atoms%nat))
    allocate(partb%tbmat(partb%norb,partb%norb))
    allocate(partb%evec(partb%norb,partb%norb))
    if(lenosky) then
        !partb%dhgenall0=f_malloc([1.to.atoms%nat,1.to.atoms%nat],id='partb%dhgenall0')
        !partb%dhgenall1=f_malloc([1.to.atoms%nat,1.to.atoms%nat],id='partb%dhgenall1')
        !partb%dhgenall2=f_malloc([1.to.atoms%nat,1.to.atoms%nat],id='partb%dhgenall2')
        !partb%dhgenall3=f_malloc([1.to.atoms%nat,1.to.atoms%nat],id='partb%dhgenall3')
        !partb%hgenall0=f_malloc([1.to.atoms%nat,1.to.atoms%nat],id='partb%hgenall0')
        !partb%hgenall1=f_malloc([1.to.atoms%nat,1.to.atoms%nat],id='partb%hgenall1')
        !partb%hgenall2=f_malloc([1.to.atoms%nat,1.to.atoms%nat],id='partb%hgenall2')
        !partb%hgenall3=f_malloc([1.to.atoms%nat,1.to.atoms%nat],id='partb%hgenall3')
        allocate(partb%hgenall0(linked_lists%maxbound_rad),source=0.d0)
        allocate(partb%hgenall1(linked_lists%maxbound_rad),source=0.d0)
        allocate(partb%hgenall2(linked_lists%maxbound_rad),source=0.d0)
        allocate(partb%hgenall3(linked_lists%maxbound_rad),source=0.d0)
        allocate(partb%dhgenall0(linked_lists%maxbound_rad),source=0.d0)
        allocate(partb%dhgenall1(linked_lists%maxbound_rad),source=0.d0)
        allocate(partb%dhgenall2(linked_lists%maxbound_rad),source=0.d0)
        allocate(partb%dhgenall3(linked_lists%maxbound_rad),source=0.d0)
    endif
    !partb%eval=f_malloc([1.to.partb%norb],id='partb%eval')
    !partb%focc=f_malloc([1.to.partb%norb],id='partb%focc')
    allocate(partb%eval(partb%norb))
    allocate(partb%focc(partb%norb))
    !call f_release_routine()
end subroutine lenoskytb_init
!*****************************************************************************************
subroutine totalenergy(pia_arr,linked_lists,parini,partb,atoms,natsi,pplocal)
    use mod_linked_lists, only: typ_pia_arr, typ_linked_lists
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_potl, only: potl_typ
    use mod_tightbinding, only: typ_partb, lenosky
    implicit none
    type(typ_linked_lists), intent(in):: linked_lists
    type(typ_pia_arr), intent(in):: pia_arr
    type(typ_parini), intent(in):: parini
    type(typ_partb), intent(inout):: partb
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in):: natsi
    type(potl_typ), intent(inout):: pplocal
    !local variables
    !'es' and 'eself(0:3)' are on-site energy of hydrogen and Si respectively
    !eself(0)= on-site energy of orbital s.
    !eself(1:3)= on-site energies of orbital p.
    real(8):: es !, eself(0:3)
    partb%pairen=0.d0
    if(lenosky) then
    call pairenergy(parini,partb,atoms,pplocal,natsi)
    endif
    call gammaenergy(pia_arr,linked_lists,partb,atoms,natsi,pplocal)
    !partb%pairen=partb%pairen-natsi*(2*eself(0)+2*eself(1)) !the value in parantesses is 0 
    es=pplocal%eps(1) !eselfgeneral2 is replaced by this line
    partb%pairen=partb%pairen-(atoms%nat-natsi)*es
    atoms%epot=partb%eband+partb%pairen
    !write(*,*) atoms%epot
end subroutine totalenergy
!*****************************************************************************************
!Given atoms with reciprocal lattice latt, number of atoms nat, at coordinates coord,
!return energy from pair potential, with forces in array force()
subroutine pairenergy(parini,partb,atoms,pplocal,natsi)
    use mod_parini, only: typ_parini
    use mod_tightbinding, only: typ_partb
    use mod_atoms, only: typ_atoms
    use mod_potl, only: potl_typ
    use mod_linked_lists, only: typ_pia_arr, typ_linked_lists
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_partb), intent(inout):: partb
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in):: natsi
    type(potl_typ), intent(in):: pplocal
    !local variables
    type(typ_linked_lists):: linked_lists
    type(typ_pia_arr):: pia_arr
    type(typ_parini):: parini_tmp
    integer:: iat, jat, kat, mat, iat_swap, ib
    integer:: atomtypej, atomtypem
    real(8):: r_swap
    real(8):: g(3), tempp(3)
    real(8):: r, y, der !, epot_pair !, cell(3)
    !Accumulate energy for pair of atoms i<j, add to pairen; latt(1),latt(2),latt(3) are primitive
    !translation vectors for the unit cell: coord=n*latt(1)+m*latt(2)+p*latt(3);
    !n,m, p are integers.
    !Set first value of forces and energy to 0
    do iat=1,atoms%nat
        atoms%fat(1,iat)=0.d0
        atoms%fat(2,iat)=0.d0
        atoms%fat(3,iat)=0.d0
    enddo
    partb%pairen=0.d0
    ! ---------------------------------
    linked_lists%rcut=partb%paircut !ann_arr%rcut
    linked_lists%triplex=.true.
    parini_tmp=parini !TO_BE_CORRECTED
    parini_tmp%bondbased_ann=.true. !TO_BE_CORRECTED
    call call_linkedlist(parini_tmp,atoms,.true.,linked_lists,pia_arr)
    do ib=1,linked_lists%maxbound_rad
        iat=linked_lists%bound_rad(1,ib)
        jat=linked_lists%bound_rad(2,ib)
        r=pia_arr%pia(ib)%r
        g(1)=-pia_arr%pia(ib)%dr(1)/r
        g(2)=-pia_arr%pia(ib)%dr(2)/r
        g(3)=-pia_arr%pia(ib)%dr(3)/r
        if(partb%usepairpot==1) then !Here I think the condition has not completed!!!
            atomtypej=0
            if(iat>natsi) atomtypej=1
            atomtypem=0
            if(jat>natsi) atomtypem=1
            call clssplint('pairenergy',pplocal%phi(atomtypej+atomtypem),r,y,der,1)
            atoms%fat(1,iat)=atoms%fat(1,iat)-der*(g(1))
            atoms%fat(2,iat)=atoms%fat(2,iat)-der*(g(2))
            atoms%fat(3,iat)=atoms%fat(3,iat)-der*(g(3))
            atoms%fat(1,jat)=atoms%fat(1,jat)+der*(g(1))
            atoms%fat(2,jat)=atoms%fat(2,jat)+der*(g(2))
            atoms%fat(3,jat)=atoms%fat(3,jat)+der*(g(3))
            partb%pairen=partb%pairen+y
        endif
    enddo
    deallocate(linked_lists%prime_bound)
    deallocate(linked_lists%bound_rad)
    deallocate(linked_lists%bound_ang)
end subroutine pairenergy
!*****************************************************************************************
subroutine lenoskytb_final(partb)
    use mod_tightbinding, only: typ_partb
    use dynamic_memory
    implicit none
    type(typ_partb), intent(inout):: partb
    !call f_routine(id='lenoskytb_final')
    !local variables
    !call f_free(partb%indorb)
    !call f_free(partb%indat)
    !call f_free(partb%norbat)
    !call f_free(partb%tbmat)
    !call f_free(partb%evec)
    !call f_free(partb%dhgenall0)
    !call f_free(partb%dhgenall1)
    !call f_free(partb%dhgenall2)
    !call f_free(partb%dhgenall3)
    !call f_free(partb%hgenall0)
    !call f_free(partb%hgenall1)
    !call f_free(partb%hgenall2)
    !call f_free(partb%hgenall3)
    !call f_free(partb%eval)
    !call f_free(partb%focc)
    deallocate(partb%indorb)
    deallocate(partb%indat)
    deallocate(partb%norbat)
    deallocate(partb%tbmat)
    deallocate(partb%evec)
    deallocate(partb%dhgenall0)
    deallocate(partb%dhgenall1)
    deallocate(partb%dhgenall2)
    deallocate(partb%dhgenall3)
    deallocate(partb%hgenall0)
    deallocate(partb%hgenall1)
    deallocate(partb%hgenall2)
    deallocate(partb%hgenall3)
    deallocate(partb%eval)
    deallocate(partb%focc)
    !call f_release_routine()
end subroutine lenoskytb_final
!*****************************************************************************************
!This routine returns h(r)'s for tight binding link: units eV and Angstrom
! Index i denotes
! 0 for s,s,sigma
! 1 for s,p,sigma
! 2 for p,p,sigma
! 3 for p,p,pi 
! Need to return h(*r), but for now returns a dummy argument
subroutine radelmgeneralsp(r,radar,dradar,atomtypei,atomtypej,pplocal)
    use mod_potl, only: potl_typ
    !use mod_const, only: ha2ev, bohr2ang
    implicit none
    type(potl_typ), intent(in):: pplocal
    real(8), intent(in):: r
    real(8), intent(out):: radar(0:3), dradar(0:3)
    integer, intent(in):: atomtypei, atomtypej
    !local variables
    integer:: i
    real(8) :: y, der !come out from clssplint routine.
    !If two atoms which ineract together are "Si", i start from 0 to 3 (because of 4 h(r)s).
    if(atomtypei+atomtypej==0) then
        do i=0,3
            !pplocal%h(i)s build various splines for different h(r)s.
            call clssplint('hgen_r2',pplocal%h(i),r,y,der,0)
            !According to Lenosky paper, "radar"s are "h(r)s" not "g(r)"s.
            !Thus, clssplint returns "g(r)"s.
            radar(i)=y/(r)/(r)
            dradar(i)=der/(r)/(r)-2.d0 * y/(r)/(r)/(r) 
        enddo
    !If two atoms which ineract together are "H" and "Si", i start from 0 to 1  (because of 2 h(r)s).
    else if(atomtypei+atomtypej==1) then
        do i=0,1 
            call clssplint('hgen_r2',pplocal%h(4+i),r,y,der,0)
            radar(i)=y/(r)/(r)
            dradar(i)=der/(r)/(r)-2.0 * y/(r)/(r)/(r) 
        enddo     
        radar(2)=0.d0
        dradar(2)=0.d0
        radar(3)=0.d0
        dradar(3)=0.d0
    !If two atoms which interact together are "H", i=0 (because of 1 possibility for interaction)
    else if(atomtypei+atomtypej==2) then
        i=0 
        call clssplint('hgen_r2',pplocal%h(6+i),r,y,der,0)
        radar(i)=y/(r)/(r)
        dradar(i)=der/(r)/(r)-2.0 * y/(r)/(r)/(r) 
        radar(1)=0.d0
        dradar(1)=0.d0
        radar(2)=0.d0
        dradar(2)=0.d0
        radar(3)=0.d0
        dradar(3)=0.d0
    endif
end subroutine radelmgeneralsp
!*****************************************************************************************
!This routine represented radial functions using cubic splines.
!modified version of clssplint: returns f in y and df/dx in deriv
!also modified to use endpoint deriv information if out of bounds
subroutine clssplint(str_action,s,xt,yt,derivt,extype)
    use mod_splinetb, only: NSPMAX, spline_typ
    use mod_const, only: ha2ev, bohr2ang
    implicit none
    character(*), intent(in) :: str_action
    type(spline_typ), intent(in) :: s
    !Flag to control method of extrapolation on lower bound 0=linear, 1= 1/r/r and 1/r 
    integer, intent(in)::  extype 
    real(8), intent(in) :: xt
    real(8), intent(out) :: yt, derivt
    !local variables
    integer :: klo, khi, k
    real(8) :: h, b, a
    real(8) ::  a1, b1
    real(8) ::  rc, yc, derc
    real(8) ::  x, y, deriv
    !cubic spline passing through points(s%x(1),s%y(1)),...,(s%x(npt),s%y(npt)).
    x=xt*bohr2ang
    if(x>s%x(1) .and. x < s%x(s%npt)) then
        !NR part
        klo=1
        khi=s%npt
        do while (khi-klo>1) 
        k=ishft(khi+klo,-1)
        if(s%x(k)>x) then
            khi=k
        else 
            klo=k
        endif
        enddo
        h=s%x(khi)-s%x(klo)
        if(h==0.d0) then 
            write(*,*) "clssplint called with x=  ", x
        endif
        a=(s%x(khi)-x)/h
        b=(x-s%x(klo))/h
        !See Numerical Recipes in C(William H. Press & etl,). chapter.section: 3.3 
        y=a*s%y(klo)+b*s%y(khi)+((a*a*a-a)*s%y2(klo)+(b*b*b-b)*s%y2(khi))*(h*h)/6
        deriv=-(1/h)*(s%y(klo)-s%y(khi)+((h*h)/6.d0)*((3*a*a-1)*s%y2(klo)-(3*b*b-1)*s%y2(khi)))
    else if(x<=s%x(1))then
        if(extype==0) then !extype is zero when extrapolation is linear.
            deriv=s%yp1
            y=s%y(1)+deriv*(x-s%x(1))
        else
            !Fit 1/x and 1/x/x termination
            rc=s%x(1)
            yc=s%y(1)
            derc=s%yp1
            b1=-rc*rc*rc*derc-rc*rc*yc
            a1=rc*yc-b1/rc
            deriv=-a1/x/x-2.d0 * b1/x/x/x
            y=a1/x+b1/x/x
        end if
    else
        deriv=s%ypn
        y=s%y(s%npt)+deriv*(x-s%x(s%npt))
    endif
    if(trim(str_action)=='pairenergy') then
        derivt=deriv/(ha2ev/bohr2ang)
        yt=y/ha2ev
    else if(trim(str_action)=='hgen_r2') then
        derivt=deriv/(ha2ev*bohr2ang)
        yt=y/(ha2ev*bohr2ang**2)
    else
        stop 'ERROR: unknown value for str_action in clssplint'
    endif

end subroutine clssplint
!*****************************************************************************************
!This subroutine initialize array's elements. 
subroutine eselfgeneral(eself)
    use mod_tightbinding, only: lenosky
    use mod_const, only: ha2ev
    implicit none
    real(8), intent(inout):: eself(0:3)
    eself(0)=-5.670225d0/ha2ev
    eself(1)=-eself(0)
    eself(2)=-eself(0)
    eself(3)=-eself(0)
end subroutine eselfgeneral
!*****************************************************************************************
subroutine prmst38c(partb,pplocal)
    use mod_tightbinding, only: typ_partb
    use mod_potl, only: potl_typ
    use mod_const, only: ha2ev, bohr2ang
    implicit none
    type(typ_partb):: partb
    integer:: unit
    integer:: mm
    real(8):: clsrhonecut, clstricut1
    integer:: clsusetri,clsmrho,clsmhtb,clsmepstb,clsmrho2
    type(potl_typ):: pplocal
    open(unit=1,file='coeff.cls',status='old')
    partb%paircut=5.24d0/bohr2ang
    clsrhonecut=4.d0
    clstricut1=4.d0
    read(1,*) partb%usepairpot,clsusetri,clsmrho,clsmhtb,clsmepstb,clsmrho2
    if(partb%usepairpot==1) then
        call clsfread_spline(1,pplocal%phi(0))
        call clsfread_spline(1,pplocal%phi(1))
        call clsfread_spline(1,pplocal%phi(2))
    endif
    do mm=0,clsmrho-1
        call clsfread_spline(1,pplocal%rho(mm))
        call clsfread_spline(1,pplocal%u(mm))
    enddo
    if(clsusetri==1) then
        call clsfread_spline(1,pplocal%f)
        call clsfread_spline(1,pplocal%g)
    endif
    do mm=0,clsmhtb-1
        call clsfread_spline(1,pplocal%h(mm))
    enddo
    do mm=0,clsmepstb-1
        read(1,*) pplocal%eps(mm)
        !pplocal%eps(mm)=pplocal%eps(mm)/ha2ev
    enddo
    do mm=0,clsmrho2-1
        call clsfread_spline(1,pplocal%xrho(mm))
        call clsfread_spline(1,pplocal%xu(mm)) 
        call clsfread_spline(1,pplocal%xf(mm)) 
        call clsfread_spline(1,pplocal%xg(mm))
    enddo
    close(1)
end subroutine prmst38c
!*****************************************************************************************
!Read s from unit, using format common with fwrite_spline(unit,s)/
subroutine clsfread_spline(unit,s)
    use mod_splinetb, only: NSPMAX, spline_typ
    !use mod_const, only: ha2ev, bohr2ang
    implicit none
    integer:: unit
    type(spline_typ), intent(out):: s
    !local variables
    integer:: i
    read(unit,*) s%npt
    read(unit,*) s%yp1, s%ypn
    !s%yp1=s%yp1/(ha2ev/bohr2ang)
    !s%ypn=s%ypn/(ha2ev/bohr2ang)
    !QUESTION: What is the meaning of below integer variables.
    read(unit,*) s%iset1, s%isetn, s%ider1, s%idern
    do i=1,s%npt
        read(unit,*) s%x(i), s%y(i), s%y2(i)
        !s%x(i)=s%x(i)/bohr2ang
        !s%y(i)=s%y(i)/ha2ev
        !s%y2(i)=s%y2(i)/(ha2ev/bohr2ang**2)
    enddo
    call clsspline(s)
end subroutine clsfread_spline
!*****************************************************************************************
!Given arrays x(1:n) and y(1:n) containing a tabulated function and 
!given values yp1 and ypn for the first derivative of the interpolating
!function at points 1 and n , respectively, 
!this routine returns an array y2(1:n) that contains the second derivatives 
!of the interpolating function at the tabulated points x(i)
!The goal of cubic spline interpolation is to get interpolation formula that is
!smooth in first derivatives and continuous in second derivative, both whithin
!an interval and at its boundaries.
subroutine clsspline(s)
    use mod_splinetb, only: NSPMAX, spline_typ
    implicit none
    type(spline_typ), intent(inout):: s
    !local variables
    integer:: i, k
    real(8):: p, qn, sig, un, u(0:NSPMAX)
    if((s%yp1) > 0.99d30) then
        s%y2(1) = 0.d0
        u(1) = 0.d0
    else 
        s%y2(1) = -0.5d0
        u(1) = (3.d0/(s%x(2)-s%x(1)))*((s%y(2)-s%y(1))/(s%x(2)-s%x(1))-s%yp1)
    endif
    do i=2,s%npt-1
        sig=(s%x(i)-s%x(i-1))/(s%x(i+1)-s%x(i-1))
        p=sig*(s%y2(i-1))+2.d0
        s%y2(i)=(sig-1.d0)/p
        u(i)=(s%y(i+1)-s%y(i))/(s%x(i+1)-s%x(i)) - (s%y(i)-s%y(i-1))/(s%x(i)-s%x(i-1))
        u(i)=(6.d0*u(i)/(s%x(i+1)-s%x(i-1))-sig*u(i-1))/p
    enddo
    if(s%ypn > 0.99d30) then
        qn = 0.d0
        un = 0.d0
    else 
        qn = 0.5d0
        un = (3.d0/(s%x(s%npt)-s%x(s%npt-1)))*(s%ypn-(s%y(s%npt)-s%y(s%npt-1))/(s%x(s%npt)-s%x(s%npt-1)))
    endif
    s%y2(s%npt)=(un-qn*u(s%npt-1))/(qn*(s%y2(s%npt-1))+1.d0)
    do k=s%npt-1,1,-1
        s%y2(k)= (s%y2(k)*(s%y2(k+1)))+u(k)
    enddo
end subroutine clsspline
!*****************************************************************************************
