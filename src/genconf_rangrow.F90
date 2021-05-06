!*****************************************************************************************
subroutine rangrow(parini,genconf)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_file_info, atom_allocate_old, atom_deallocate_old
    use mod_atoms, only: set_rcov, determinexyzminmax, set_rat_atoms
    use mod_genconf, only: typ_genconf
    use mod_const, only: ang2bohr
    use mod_acf, only: acf_read, acf_write
    implicit none
    type(typ_parini), intent(inout):: parini
    type(typ_genconf), intent(in):: genconf
    !local variables
    type(typ_atoms):: atoms, atoms_p
    type(typ_file_info):: file_info
    integer:: iat, mat, itry, nat_yet
    real(8):: dist, distmin, dx, dy, dz, ranxyz(3), amargin, amargin_xrel, amargin_yrel, amargin_zrel
    real(8):: r, theta, phi, rmin, rmax, pi, twopi, tt
    logical:: good
    real(8):: xmin, ymin, zmin, xmax, ymax, zmax
    if(parini%fbmin_genconf<0.d0) parini%fbmin_genconf=0.75d0
    write(*,*) 'fbmin_genconf=',parini%fbmin_genconf
    call acf_read(parini,'posinp.acf',1,atoms=atoms)
    write(*,'(a,i6)') 'nat ',atoms%nat
    if(atoms%nat<1) stop 'ERROR: atoms%nat<1'
    call atom_allocate_old(atoms_p,atoms%nat+genconf%nat_add,0,0,sat=.true.,bemoved=.true.,rcov=.true.)
    call set_rat_atoms(atoms_p,atoms,setall=.true.)
    atoms_p%sat(1:atoms%nat)=atoms%sat(1:atoms%nat)
    atoms_p%bemoved(1:3,1:atoms%nat)=atoms%bemoved(1:3,1:atoms%nat)
    atoms_p%cellvec(1:3,1:3)=atoms%cellvec(1:3,1:3)
    atoms_p%boundcond=atoms%boundcond
    if(trim(parini%sat_genconf)=='unknown') then
        stop 'ERROR: sat_genconf must be defined in block [genconf] input.ini'
    else
        !atoms_p%sat(atoms%nat+1:atoms_p%nat)=trim(parini%sat_genconf)
        do iat=atoms%nat+1,atoms_p%nat
            !atoms_p%sat(atoms%nat+1:atoms_p%nat)='Na'
            if(mod(iat,2)==0) then
                atoms_p%sat(iat)='Cl'
            else
                atoms_p%sat(iat)='Na'
            endif
        enddo
    endif
    call set_rcov(atoms_p)
    pi=4.d0*atan(1.d0)
    twopi=2.d0*pi
    amargin=parini%amargin_genconf
    nat_yet=atoms%nat
    itry=0
    do
        itry=itry+1
        call random_number(tt)
        mat=int(tt*real(nat_yet,8))+1
        nat_yet=nat_yet+1
        if(nat_yet>atoms_p%nat) exit
        if(trim(atoms_p%sat(nat_yet))==trim(atoms_p%sat(mat))) then
            nat_yet=nat_yet-1
            cycle
        endif
        rmin=parini%fbmin_genconf*(atoms_p%rcov(nat_yet)+atoms_p%rcov(mat))
        rmax=1.2*                 (atoms_p%rcov(nat_yet)+atoms_p%rcov(mat))
        call random_number(r)
        call random_number(theta)
        call random_number(phi)
        r=min(max(r,rmin),rmax)
        theta=theta*pi
        phi=phi*twopi
        atoms_p%ratp(1,nat_yet)=atoms_p%ratp(1,mat)+r*sin(theta)*cos(phi)
        atoms_p%ratp(2,nat_yet)=atoms_p%ratp(2,mat)+r*sin(theta)*sin(phi)
        atoms_p%ratp(3,nat_yet)=atoms_p%ratp(3,mat)+r*cos(theta)
        good=.true.
        do iat=1,nat_yet-1
            dx=atoms_p%ratp(1,iat)-atoms_p%ratp(1,nat_yet)
            dy=atoms_p%ratp(2,iat)-atoms_p%ratp(2,nat_yet)
            dz=atoms_p%ratp(3,iat)-atoms_p%ratp(3,nat_yet)
            dist=sqrt(dx**2+dy**2+dz**2)
            !distmin=parini%fbmin_genconf*(atoms_p%rcov(iat)+atoms_p%rcov(nat_yet))
            if(trim(atoms_p%sat(nat_yet))=='Na' .and. trim(atoms_p%sat(iat))=='Na') then
                distmin=2.7d0*ang2bohr*0.9d0
            elseif(trim(atoms_p%sat(nat_yet))=='Cl' .and. trim(atoms_p%sat(iat))=='Cl') then
                distmin=2.7d0*ang2bohr*0.9d0
            else
                distmin=2.25d0*ang2bohr*0.8d0
            endif
            if(dist<distmin) then
                good=.false.
                exit
            endif
        enddo
        if(.not. good) then
            nat_yet=nat_yet-1
            cycle
        endif
        atoms_p%bemoved(1:3,nat_yet)=.true.
        write(*,'(a,2i3,a5,3es24.15,a,i6)') 'new atom #',nat_yet,mat, &
            trim(atoms_p%sat(nat_yet)),atoms_p%ratp(1,nat_yet),atoms_p%ratp(2,nat_yet),atoms_p%ratp(3,nat_yet),' is added: itry= ',itry
    enddo
    atoms_p%cellvec=1.d2
    call determinexyzminmax(atoms_p%nat,atoms_p%ratp,atoms_p%cellvec,xmin,ymin,zmin,xmax,ymax,zmax)
    do iat=1,atoms_p%nat
        atoms_p%ratp(1,iat)=atoms_p%ratp(1,iat)-xmin+amargin
        atoms_p%ratp(2,iat)=atoms_p%ratp(2,iat)-ymin+amargin
        atoms_p%ratp(3,iat)=atoms_p%ratp(3,iat)-zmin+amargin
    enddo
    atoms_p%cellvec=0.d0
    atoms_p%cellvec(1,1)=xmax-xmin+2.d0*amargin
    atoms_p%cellvec(2,2)=ymax-ymin+2.d0*amargin
    atoms_p%cellvec(3,3)=zmax-zmin+2.d0*amargin
    !atoms%sat(1:atoms_p%nat)=atoms_p%sat(1:atoms_p%nat)
    !atoms%bemoved(1:3,1:atoms_p%nat)=atoms_p%bemoved(1:3,1:atoms_p%nat)
    file_info%filename_positions='posout.acf'
    file_info%file_position='new'
    call acf_write(file_info,atoms=atoms_p,strkey='posout')
    call atom_deallocate_old(atoms,sat=.true.,rat=.true.,bemoved=.true.)
    call atom_deallocate_old(atoms_p,sat=.true.,rat=.true.,bemoved=.true.)
    !call deallocateatomsarrays
end subroutine rangrow
!*****************************************************************************************
