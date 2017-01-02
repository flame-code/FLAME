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
subroutine lenoskytb_alborz(atoms,natsi,count_md)
    use mod_interface
    use mod_atoms, only: typ_atoms
    use mod_frame, only: CLSMAXATOM, CLSMAXNEIGHB
    use mod_potl, only: potl_typ
    use mod_tightbinding, only: typ_partb, lenosky
    implicit none
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in):: natsi
    real(8), intent(inout):: count_md
    !local variables
    integer, save:: firstcall=1
    type(potl_typ), save:: pplocal
    type(typ_partb), save:: partb
    if(firstcall==1) then
        write(*,'(a)') 'GAMMA POINT only tight binding code'
        if(lenosky) then
        write(*,'(a)') 'Reading spline potential coeff.cls'
        call prmst38c(partb,pplocal) !Reads potential 
        endif
        write(*,'(a,f)') 'paircut= ',partb%paircut
        call lenoskytb_init(partb,atoms,natsi)
        firstcall=0
    endif
    count_md=count_md+1.d0
    !PRINT SOME WARNINGS
    if(atoms%nat>CLSMAXATOM) then
        write(*,'(a,i,a,i)') 'WARNING Number of atoms = ',atoms%nat,' is greater than CLSMAXATOM = ',CLSMAXATOM
        write(*,'(a)') 'THIS MAY CAUSE FAILURE OR UNPREDICTABLE BEHAVIOR'
        write(*,'(a)') 'INCREASE CLSMAXATOM in defines.h and RECOMPILE'
    endif
    if(atoms%nat>CLSMAXNEIGHB) then
        write(*,'(a,i,a,i)') 'WARNING Number of atoms = ',atoms%nat,' is greater than CLSMAXNEIGHB = ',CLSMAXNEIGHB
        write(*,'(a)') 'THIS MAY CAUSE FAILURE OR UNPREDICTABLE BEHAVIOR'
        write(*,'(a)') 'INCREASE CLSMAXNEIGHB in defines.h and RECOMPILE'
    endif
    if(natsi>atoms%nat) write(*,'(a)') 'WARNING natsi = ',natsi,' is greater than number of atoms = ',atoms%nat
    if(atoms%nat == 0) write(*,'(a)') 'WARNING lenoskytb_alborz called with zero atoms'
    if(atoms%nat < 0) write(*,'(a)') 'WARNING lenoskytb_alborz called with negative number of atoms'
    call totalenergy(partb,atoms,natsi,pplocal)
    !call lenoskytb_final(partb)
end subroutine lenoskytb_alborz
!*****************************************************************************************
subroutine lenoskytb_init(partb,atoms,natsi)
    use mod_interface
    use mod_tightbinding, only: typ_partb, lenosky
    use mod_atoms, only: typ_atoms
    use dynamic_memory
    implicit none
    type(typ_partb), intent(inout):: partb
    type(typ_atoms), intent(in):: atoms
    integer, intent(in):: natsi
    !local variables
    call f_routine(id='lenoskytb_init')
    partb%norb=partb%nstride*natsi+(atoms%nat-natsi)
    !write(*,*) partb%nstride, natsi, partb%norb
    !Number of eigenvalues to compute=MIN(matrix*frac+extra,dimension of matrix)
    if(floor(partb%norb*partb%frac+partb%extra)>partb%norb) then
        partb%norbcut=partb%norb
    else
        partb%norbcut=floor((partb%norb*partb%frac+partb%extra))
    endif
    !---------------------------------------------------------
    partb%indorb=f_malloc([1.to.partb%norb],id='partb%indorb')
    partb%indat=f_malloc([1.to.atoms%nat],id='partb%indat')
    partb%norbat=f_malloc([1.to.atoms%nat],id='partb%norbat')
    partb%tbmat=f_malloc([1.to.partb%norb,1.to.partb%norb],id='partb%tbmat')
    partb%evec=f_malloc([1.to.partb%norb,1.to.partb%norbcut],id='partb%evec')
    if(lenosky) then
        partb%dhgenall0=f_malloc([1.to.atoms%nat,1.to.atoms%nat],id='partb%dhgenall0')
        partb%dhgenall1=f_malloc([1.to.atoms%nat,1.to.atoms%nat],id='partb%dhgenall1')
        partb%dhgenall2=f_malloc([1.to.atoms%nat,1.to.atoms%nat],id='partb%dhgenall2')
        partb%dhgenall3=f_malloc([1.to.atoms%nat,1.to.atoms%nat],id='partb%dhgenall3')
        partb%hgenall0=f_malloc([1.to.atoms%nat,1.to.atoms%nat],id='partb%hgenall0')
        partb%hgenall1=f_malloc([1.to.atoms%nat,1.to.atoms%nat],id='partb%hgenall1')
        partb%hgenall2=f_malloc([1.to.atoms%nat,1.to.atoms%nat],id='partb%hgenall2')
        partb%hgenall3=f_malloc([1.to.atoms%nat,1.to.atoms%nat],id='partb%hgenall3')
    endif
    partb%eval=f_malloc([1.to.partb%norbcut],id='partb%eval')
    partb%focc=f_malloc([1.to.partb%norbcut],id='partb%focc')
    call f_release_routine()
end subroutine lenoskytb_init
!*****************************************************************************************
subroutine totalenergy(partb,atoms,natsi,pplocal)
    use mod_interface
    use mod_atoms, only: typ_atoms
    use mod_potl, only: potl_typ
    use mod_tightbinding, only: typ_partb, lenosky
    use mod_const, only: ha2ev, bohr2ang
    implicit none
    type(typ_partb), intent(inout):: partb
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in):: natsi
    type(potl_typ), intent(inout):: pplocal
    !local variables
    integer:: iat
    !'es' and 'eself(0:3)' are on-site energy of hydrogen and Si respectively
    !eself(0)= on-site energy of orbital s.
    !eself(1:3)= on-site energies of orbital p.
    real(8):: es !, eself(0:3)
    partb%pairen=0.d0
    !if(.not. lenosky) then
    !do iat=1,atoms%nat
    !    atoms%rat(1,iat)=atoms%rat(1,iat)*bohr2ang
    !    atoms%rat(2,iat)=atoms%rat(2,iat)*bohr2ang
    !    atoms%rat(3,iat)=atoms%rat(3,iat)*bohr2ang
    !enddo
    !atoms%cellvec(1:3,1:3)=atoms%cellvec(1:3,1:3)*bohr2ang
    !call pairenergy(partb,atoms,pplocal,natsi)
    !do iat=1,atoms%nat
    !    atoms%rat(1,iat)=atoms%rat(1,iat)/bohr2ang
    !    atoms%rat(2,iat)=atoms%rat(2,iat)/bohr2ang
    !    atoms%rat(3,iat)=atoms%rat(3,iat)/bohr2ang
    !enddo
    !do iat=1,atoms%nat
    !    atoms%fat(1,iat)=atoms%fat(1,iat)/(ha2ev/bohr2ang)
    !    atoms%fat(2,iat)=atoms%fat(2,iat)/(ha2ev/bohr2ang)
    !    atoms%fat(3,iat)=atoms%fat(3,iat)/(ha2ev/bohr2ang)
    !enddo
    !endif
    if(lenosky) then
    call pairenergy(partb,atoms,pplocal,natsi)
    endif
    call gammaenergy(partb,atoms,natsi,pplocal)
    !partb%pairen=partb%pairen-natsi*(2*eself(0)+2*eself(1)) !the value in parantesses is 0 
    es=pplocal%eps(1) !eselfgeneral2 is replaced by this line
    partb%pairen=partb%pairen-(atoms%nat-natsi)*es
    !if(lenosky) then 
        atoms%epot=partb%eband+partb%pairen
    !endif
    !atoms%epot=partb%eband+(partb%pairen/ha2ev)
    !write(*,*) atoms%epot
end subroutine totalenergy
!*****************************************************************************************
!Given atoms with reciprocal lattice latt, number of atoms nat, at coordinates coord,
!return energy from pair potential, with forces in array force()
subroutine pairenergy(partb,atoms,pplocal,natsi)
    use mod_interface
    use mod_tightbinding, only: typ_partb, lenosky
    use mod_atoms, only: typ_atoms
    use mod_potl, only: potl_typ
    use mod_frame, only: clsframepp_type
    use mod_const, only: ha2ev, bohr2ang
    implicit none
    type(typ_partb), intent(inout):: partb
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in):: natsi
    type(potl_typ), intent(in):: pplocal
    !local variables
    integer:: iat, jat, kat, mat, iat_swap
    integer:: atomtypej, atomtypem
    real(8):: r_swap
    real(8):: g(3), tempp(3)
    type(clsframepp_type):: framepp
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
    !Set number of neigbors to 0
    do iat=1,atoms%nat
        framepp%nb(iat)=0
    enddo
    !Find atoms in space, then connect them together.
    do iat=1,atoms%nat
        do jat=iat+1,atoms%nat
            !Set CELLDIST that calculates the shortest distance between two atoms.
            !Because cells are assumed periodic.
            call CELLDIST_F90(atoms%rat(1,iat),atoms%rat(1,jat),atoms%cellvec,r)
            !Impose atoms intract together from 0 to 'paircut' distance.
            if(.not. r>partb%paircut) then
                framepp%nb(iat)=framepp%nb(iat)+1
                !Set CELLGRAD that calculates derivatives, because we need them to
                !calculate forces.
                call CELLGRAD_F90(atoms%rat(1,iat),atoms%rat(1,jat),atoms%cellvec,g)
                !add atom jat to list of atom iat, because atom jat is neighbor of atom iat according to if-condition.
                framepp%nbind(iat,framepp%nb(iat))=jat
                framepp%r(iat,framepp%nb(iat))=r
                framepp%grad(1,iat,framepp%nb(iat))=g(1)
                framepp%grad(2,iat,framepp%nb(iat))=g(2)
                framepp%grad(3,iat,framepp%nb(iat))=g(3)
                !In fact the condition of r<=paircut determine the number of neigbors
            endif
        enddo
    enddo
    !Do sorting of neighbor lists in order of increasing distance*/
    do iat=1,atoms%nat
        do jat=1,framepp%nb(iat)
            do kat=jat+1,framepp%nb(iat)
                if(framepp%r(iat,jat)>framepp%r(iat,kat)) then
                    r_swap=framepp%r(iat,jat)
                    framepp%r(iat,jat)=framepp%r(iat,kat)
                    framepp%r(iat,kat)=r_swap
                    iat_swap=framepp%nbind(iat,jat)
                    framepp%nbind(iat,jat)=framepp%nbind(iat,kat)
                    framepp%nbind(iat,kat)=iat_swap
                    tempp(1)=framepp%grad(1,iat,jat)
                    tempp(2)=framepp%grad(2,iat,jat)
                    tempp(3)=framepp%grad(3,iat,jat)
                    framepp%grad(1,iat,jat)=framepp%grad(1,iat,kat)
                    framepp%grad(2,iat,jat)=framepp%grad(2,iat,kat)
                    framepp%grad(3,iat,jat)=framepp%grad(3,iat,kat)
                    framepp%grad(1,iat,kat)=tempp(1)
                    framepp%grad(2,iat,kat)=tempp(2)
                    framepp%grad(3,iat,kat)=tempp(3)
                endif
            enddo
        enddo
    enddo
    ! ---------------------------------
    do iat=1,atoms%nat
        do jat=1,framepp%nb(iat)
            mat=framepp%nbind(iat,jat)
            g(1)=framepp%grad(1,iat,jat)
            g(2)=framepp%grad(2,iat,jat)
            g(3)=framepp%grad(3,iat,jat)
            r=framepp%r(iat,jat)
            if(partb%usepairpot==1) then !Here I think the condition has not completed!!!
                atomtypej=0
                if(iat>natsi) atomtypej=1
                atomtypem=0
                if(mat>natsi) atomtypem=1
                call clssplint(pplocal%phi(atomtypej+atomtypem),r,y,der,1)
                atoms%fat(1,iat)=atoms%fat(1,iat)-der*(g(1))
                atoms%fat(2,iat)=atoms%fat(2,iat)-der*(g(2))
                atoms%fat(3,iat)=atoms%fat(3,iat)-der*(g(3))
                atoms%fat(1,mat)=atoms%fat(1,mat)+der*(g(1))
                atoms%fat(2,mat)=atoms%fat(2,mat)+der*(g(2))
                atoms%fat(3,mat)=atoms%fat(3,mat)+der*(g(3))
                if(iat<=atoms%nat) partb%pairen=partb%pairen+y/2 !QUESTION: Why y is divided by 2 ?
                if(mat<=atoms%nat) partb%pairen=partb%pairen+y/2
            endif
        enddo !End of loop over pairs
    enddo
    end subroutine pairenergy
!*****************************************************************************************
subroutine lenoskytb_final(partb)
    use mod_interface
    use mod_tightbinding, only: typ_partb
    use dynamic_memory
    implicit none
    type(typ_partb), intent(inout):: partb
    call f_routine(id='lenoskytb_final')
    !local variables
    call f_free(partb%indorb)
    call f_free(partb%indat)
    call f_free(partb%norbat)
    call f_free(partb%tbmat)
    call f_free(partb%evec)
    call f_free(partb%dhgenall0)
    call f_free(partb%dhgenall1)
    call f_free(partb%dhgenall2)
    call f_free(partb%dhgenall3)
    call f_free(partb%hgenall0)
    call f_free(partb%hgenall1)
    call f_free(partb%hgenall2)
    call f_free(partb%hgenall3)
    call f_free(partb%eval)
    call f_free(partb%focc)
    call f_release_routine()
end subroutine lenoskytb_final
!*****************************************************************************************
subroutine VECT_SUBTRACT_F90(a,b,c) 
    use mod_interface
    implicit none
    real(8), intent(in):: a(3), b(3)
    real(8), intent(out):: c(3)
    c(1)=a(1)-b(1) 
    c(2)=a(2)-b(2) 
    c(3)=a(3)-b(3)
end subroutine VECT_SUBTRACT_F90
!*****************************************************************************************
subroutine APPLY_PBC_F90(a)
    use mod_interface
    implicit none
    real(8), intent(inout) :: a(3)
    a(1)=a(1)-floor(a(1)+0.5d0)
    a(2)=a(2)-floor(a(2)+0.5d0)
    a(3)=a(3)-floor(a(3)+0.5d0)
end subroutine APPLY_PBC_F90
!*****************************************************************************************
subroutine CELLDIST_F90(p1,p2,a,d)
    use mod_interface
    implicit none
    real(8), intent(in):: p1(3), p2(3), a(3,3)
    real(8), intent(out):: d
    !local variables
    real(8):: diff(3)
    call VECT_SUBTRACT_F90(p1,p2,diff)
    diff(1)=diff(1)/a(1,1)
    diff(2)=diff(2)/a(2,2)
    diff(3)=diff(3)/a(3,3)
    call APPLY_PBC_F90(diff)
    diff(1)=diff(1)*a(1,1)
    diff(2)=diff(2)*a(2,2)
    diff(3)=diff(3)*a(3,3)
    d=sqrt(diff(1)**2+diff(2)**2+diff(3)**2)
end subroutine CELLDIST_F90
!*****************************************************************************************
subroutine CELLGRAD_F90(p1,p2,a,g)
    use mod_interface
    implicit none
    real(8), intent(in):: p1(3), p2(3), a(3,3)
    real(8), intent(out):: g(3)
    !local variables
    real(8):: r, diff(3)
    call VECT_SUBTRACT_F90(p1,p2,diff)
    diff(1)=diff(1)/a(1,1)
    diff(2)=diff(2)/a(2,2)
    diff(3)=diff(3)/a(3,3)
    call APPLY_PBC_F90(diff)
    diff(1)=diff(1)*a(1,1)
    diff(2)=diff(2)*a(2,2)
    diff(3)=diff(3)*a(3,3)
    r=sqrt(diff(1)**2+diff(2)**2+diff(3)**2)
    g(1)=diff(1)/r
    g(2)=diff(2)/r
    g(3)=diff(3)/r
end subroutine CELLGRAD_F90
!*****************************************************************************************
!This routine returns h(r)'s for tight binding link: units eV and Angstrom
! Index i denotes
! 0 for s,s,sigma
! 1 for s,p,sigma
! 2 for p,p,sigma
! 3 for p,p,pi 
! Need to return h(*r), but for now returns a dummy argument
subroutine radelmgeneralsp(r,radar,dradar,atomtypei,atomtypej,pplocal)
    use mod_interface
    use mod_potl, only: potl_typ
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
            call clssplint(pplocal%h(i),r,y,der,0)
            !According to Lenosky paper, "radar"s are "h(r)s" not "g(r)"s.
            !Thus, clssplint returns "g(r)"s.
            radar(i)=y/(r)/(r)
            dradar(i)=der/(r)/(r)-2.d0 * y/(r)/(r)/(r) 
        enddo
    !If two atoms which ineract together are "H" and "Si", i start from 0 to 1  (because of 2 h(r)s).
    else if(atomtypei+atomtypej==1) then
        do i=0,1 
            call clssplint(pplocal%h(4+i),r,y,der,0)
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
        call clssplint(pplocal%h(6+i),r,y,der,0)
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
subroutine clssplint(s,x,y,deriv,extype)
    use mod_interface
    use mod_splinetb, only: NSPMAX, spline_typ
    implicit none
    type(spline_typ), intent(in) :: s
    !Flag to control method of extrapolation on lower bound 0=linear, 1= 1/r/r and 1/r 
    integer, intent(in)::  extype 
    real(8), intent(in) :: x
    real(8), intent(out) :: y, deriv
    !local variables
    integer :: klo, khi, k
    real(8) :: h, b, a
    real(8) ::  a1, b1
    real(8) ::  rc, yc, derc
    !cubic spline passing through points(s%x(1),s%y(1)),...,(s%x(npt),s%y(npt)).
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
end subroutine clssplint
!*****************************************************************************************
!This subroutine initialize array's elements. 
subroutine eselfgeneral(eself)
    use mod_interface
    use mod_tightbinding, only: lenosky
    implicit none
    real(8), intent(inout):: eself(0:3)
    if(lenosky) then 
        eself(0)=-5.670225d0
    else
        eself(0)=-5.670225d0/27.211385d0
    endif
    eself(1)=-eself(0)
    eself(2)=-eself(0)
    eself(3)=-eself(0)
end subroutine eselfgeneral
!*****************************************************************************************
subroutine prmst38c(partb,pplocal)
    use mod_interface
    use mod_tightbinding, only: typ_partb
    use mod_potl, only: potl_typ
    implicit none
    type(typ_partb):: partb
    integer:: unit
    integer:: mm
    real(8):: clsrhonecut, clstricut1
    integer:: clsusetri,clsmrho,clsmhtb,clsmepstb,clsmrho2
    type(potl_typ):: pplocal
    open(unit=1,file='coeff.cls',status='old')
    partb%paircut=5.24d0
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
    use mod_interface
    use mod_splinetb, only: NSPMAX, spline_typ
    implicit none
    integer:: unit
    type(spline_typ), intent(out):: s
    !local variables
    integer:: i
    read(unit,*) s%npt
    read(unit,*) s%yp1, s%ypn
    !QUESTION: What is the meaning of below integer variables.
    read(unit,*) s%iset1, s%isetn, s%ider1, s%idern
    do i=1,s%npt
        read(unit,*) s%x(i), s%y(i), s%y2(i)
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
    use mod_interface
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
