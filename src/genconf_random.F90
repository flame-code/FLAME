!*****************************************************************************************
subroutine genrandom(parini,genconf)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_file_info, atom_allocate_old, atom_deallocate_old
    use mod_atoms, only: set_rcov, update_rat, set_rat_atoms, update_ratp
    use mod_genconf, only: typ_genconf
    use mod_acf, only: acf_read, acf_write
    implicit none
    type(typ_parini), intent(inout):: parini
    type(typ_genconf), intent(in):: genconf
    !local variables
    type(typ_atoms):: atoms, atoms_p
    type(typ_file_info):: file_info
    integer:: ix, iy, iz, iat, mat, itry, nx, ny, nz, i
    !character(256):: comment1, comment2
    real(8), allocatable:: ratall(:,:,:,:,:)
    !logical, allocatable:: atom_motion_p(:,:)
    real(8):: r, dx, dy, dz, ranxyz(3), rmin,ang2rad,tmp
    real(8):: amargin, amargin_xrel, amargin_yrel, amargin_zrel
    real(8):: rcov_max,vol_min,vol,alpha,beta,gama,tt1
    real(8):: dmin,amin,amax,bmin,bmax,a,b,c,a_cell,b_cell,c_cell,tt,zeta 
    real(8):: rnd(3)
    logical:: reject
    real(8),parameter::pi=4.d0*atan(1.d0)
    amargin=parini%amargin_genconf
    !write(*,*)parini%stypat_genconf
    write(*,*)parini%iseed,'seed'
    !call random_seed(parini%iseed)
    !call random_seed()
    !call count_words(parini%stypat_genconf,ntypat)
    !read(parini%stypat_genconf,*) stypat(1:ntypat)
    !if(parini%ntypat==0) stop 'ERROR: number of type of atoms zero in genrandom'
    if(parini%fbmin_genconf<0.d0) parini%fbmin_genconf=0.75d0
    if(parini%fbmax_genconf<0.d0) parini%fbmax_genconf=1.5d0
    write(*,'(a,2f6.2)') 'fbmin,fbmax= ',parini%fbmin_genconf,parini%fbmax_genconf
    call acf_read(parini,'posinp.acf',1,atoms=atoms)
    write(*,'(a,i6)') 'nat ',atoms%nat
    call atom_allocate_old(atoms_p,atoms%nat+genconf%nat_add,0,0,sat=.true.,bemoved=.true.,rcov=.true.)
    allocate(ratall(3,atoms_p%nat,-1:1,-1:1,-1:1))
    if (trim(atoms%boundcond)=='free') then
        nx=0
        ny=0
        nz=0
    elseif (trim(atoms%boundcond)=='wire') then
        nx=1
        ny=0
        nz=0
    elseif (trim(atoms%boundcond)=='slab') then
        nx=1
        ny=1
        nz=0
    elseif (trim(atoms%boundcond)=='bulk') then
        nx=1
        ny=1
        nz=1
    endif
    !if(trim(parini%stypat_genconf)=='unknown') then
    !    stop 'ERROR: stypat_genconf must be defined in block [genconf] input.ini'
    !else
        do i=1,parini%ntypat
            atoms_p%sat(i:atoms_p%nat:parini%ntypat)=parini%stypat(i)
        enddo
    !endif
    call set_rcov(atoms_p)
    if (parini%variable_cell_genconf) then
        atoms%cellvec=0.d0
        vol_min =0.d0
        do iat=1,atoms_p%nat
            vol_min = vol_min+(2.d0*atoms_p%rcov(iat))**3
        enddo
        !rcov_max=maxval(atoms_p%rcov)
        !vol_min=8*atoms_p%nat*(rcov_max)**3
        vol= vol_min * 1.d0
        ang2rad=pi/180.d0
        dmin = (vol)**(1.d0/3.d0)
        amin = dmin*0.8d0
        amax = dmin*1.2d0
        call random_number(a)
        a=a*(amax-amin)+amin
        bmin = dmin*0.8d0
        bmax = dmin*1.2d0
        call random_number(b)
        b=b*(bmax-bmin)+bmin
        a_cell=a
        c=vol/(a*b)
        a=a+2.d0*(1-nx)*amargin
        b=b+bmin+2.d0*(1-ny)*amargin
        c=c+2.d0*(1-nz)*amargin
        if (parini%nonorthogonal_genconf) then 
100         continue
            call gausdist_alborz(3,rnd)
            rnd(1:3)=rnd(1:3)*20.d0
            write(*,'(a,3f10.2)') 'RND ',rnd(1),rnd(2),rnd(3)
            !call random_number(alpha)
            !call random_number(beta)
            !call random_number(gama)
            alpha=90.d0+min(max(-45.d0,rnd(1)),30.d0)
            beta= 90.d0+min(max(-45.d0,rnd(2)),30.d0)
            gama= 90.d0+min(max(-45.d0,rnd(3)),30.d0)
            alpha=alpha*ang2rad
            beta=beta*ang2rad
            gama=gama*ang2rad
            tt1=1.d0+2.d0*cos(alpha)*cos(beta)*cos(gama)-cos(alpha)**2-cos(beta)**2-cos(gama)**2
            if (tt1 <= 0.d0) goto 100
            write(*,'(a,3es20.10)') 'angles are:',alpha/ang2rad,beta/ang2rad,gama/ang2rad
            b_cell=b/sin(gama)
            atoms%cellvec(1,1)=a_cell
            atoms%cellvec(1,2)=b_cell*cos(gama)
            atoms%cellvec(2,2)=b
            tt=sqrt(tt1)
            c_cell=vol/(a_cell*b_cell*tt)
            zeta=(cos(alpha)-cos(gama)*cos(beta))/sin(gama)
            atoms%cellvec(1,3)=c_cell*cos(beta)
            atoms%cellvec(2,3)=c_cell*zeta
            atoms%cellvec(3,3)=c_cell*sqrt(1-cos(beta)**2-zeta**2)
        else
            b_cell=b
            c_cell=vol/(a_cell*b_cell)
            atoms%cellvec(1,1)=a_cell
            atoms%cellvec(2,2)=b_cell
            atoms%cellvec(3,3)=c_cell
        endif
    else
        a=atoms%cellvec(1,1)
        b=sqrt(atoms%cellvec(1,2)**2+atoms%cellvec(2,2)**2)
        c=sqrt(atoms%cellvec(1,3)**2+atoms%cellvec(2,3)**2+atoms%cellvec(3,3)**2)
    endif
    write(*,*)atoms%cellvec(1,1),atoms%cellvec(1,2), atoms%cellvec(2,2)
    write(*,*)atoms%cellvec(1,3),atoms%cellvec(2,3),atoms%cellvec(3,3)
    

    call update_ratp(atoms)
    do iz=-nz,nz
    do iy=-ny,ny
    do ix=-nx,nx
        do iat=1,atoms%nat
            ratall(1,iat,ix,iy,iz)=atoms%ratp(1,iat)+ix*atoms%cellvec(1,1)+iy*atoms%cellvec(1,2)+iz*atoms%cellvec(1,3)
            ratall(2,iat,ix,iy,iz)=atoms%ratp(2,iat)+iy*atoms%cellvec(2,2)+iz*atoms%cellvec(2,3)
            ratall(3,iat,ix,iy,iz)=atoms%ratp(3,iat)+iz*atoms%cellvec(3,3)
        enddo
    enddo
    enddo
    enddo
    call set_rat_atoms(atoms_p,atoms,setall=.true.)
    atoms_p%bemoved(1:3,1:atoms%nat)=atoms%bemoved(1:3,1:atoms%nat)
    atoms_p%cellvec(1:3,1:3)=atoms%cellvec(1:3,1:3)
    atoms_p%boundcond=atoms%boundcond
    amargin_xrel=amargin/a
    amargin_yrel=amargin/b
    amargin_zrel=amargin/c
    mat=atoms%nat
    itry=0
    do
        itry=itry+1
        mat=mat+1
        
        if(mat>atoms_p%nat) exit
        call random_number(ranxyz)
        ranxyz(1)=ranxyz(1)*(1.d0-2.d0*amargin_xrel)+amargin_xrel
        ranxyz(2)=ranxyz(2)*(1.d0-2.d0*amargin_yrel)+amargin_yrel
        ranxyz(3)=ranxyz(3)*(1.d0-2.d0*amargin_zrel)+amargin_zrel
        call  rxyz_int2cart_alborz(1,atoms%cellvec,ranxyz,atoms_p%ratp(:,mat))
        call update_rat(atoms_p,upall=.true.)
        reject=.false.
        rmin=1.d20
        loop_iz: do iz=-nz,nz
        do iy=-ny,ny
        do ix=-nx,nx
            do iat=1,mat-1
                dx=ratall(1,iat,ix,iy,iz)-atoms_p%ratp(1,mat)
                dy=ratall(2,iat,ix,iy,iz)-atoms_p%ratp(2,mat)
                dz=ratall(3,iat,ix,iy,iz)-atoms_p%ratp(3,mat)
                r=sqrt(dx**2+dy**2+dz**2)
                !if((r<0.9d0 .and. trim(atoms_p%sat(iat))=='H') .or. &
                !   (r<1.9d0 .and. trim(atoms_p%sat(iat))=='Si')) then
                if (parini%ntypat>=2 .and. trim(atoms_p%sat(iat))== trim(atoms_p%sat(mat)))then
                    tmp=1.5
                else
                    tmp=1
                endif
                if(r<parini%fbmin_genconf*tmp*(atoms_p%rcov(iat)+atoms_p%rcov(mat))) then
                !if(r<parini%fbmin_genconf/0.529) then
                    reject=.true.
                    exit loop_iz
                endif
                rmin=min(rmin,r)
            enddo
        enddo
        enddo
        enddo loop_iz
            !if (ntypat>=2 .and. trim(atoms_p%sat(iat))== trim(atoms_p%sat(mat)))then
            !    tmp=1.5
            !else
                tmp=1
            !endif
        do iat=1,mat-1
            if(rmin>parini%fbmax_genconf*tmp*(atoms_p%rcov(iat)+atoms_p%rcov(mat))) reject=.true.
        enddo
        !if(rmin>parini%fbmax_genconf/0.529) reject=.true.
        if(reject) then
            mat=mat-1
            cycle
        endif
        atoms_p%bemoved(1:3,mat)=.true.
        write(*,'(a,i3,a5,3es24.15,a,i6,es15.7)') 'new atom #',mat, &
            trim(atoms_p%sat(mat)),atoms_p%ratp(1,mat),atoms_p%ratp(2,mat),atoms_p%ratp(3,mat),' is added: itry= ',itry,rmin
        do iz=-nz,nz
        do iy=-ny,ny
        do ix=-nx,nx
            ratall(1,mat,ix,iy,iz)=atoms_p%ratp(1,iat)+ix*atoms%cellvec(1,1)+iy*atoms%cellvec(1,2)+iz*atoms%cellvec(1,3)
            ratall(2,mat,ix,iy,iz)=atoms_p%ratp(2,iat)+iy*atoms%cellvec(2,2)+iz*atoms%cellvec(2,3)
            ratall(3,mat,ix,iy,iz)=atoms_p%ratp(3,iat)+iz*atoms%cellvec(3,3)
        enddo
        enddo
        enddo
    enddo

    !atoms%sat(1:atoms_p%nat)=atoms_p%sat(1:atoms_p%nat)
    !atoms%bemoved(1:3,1:atoms_p%nat)=atoms_p%bemoved(1:3,1:atoms_p%nat)
    file_info%filename_positions='posout.acf'
    file_info%file_position='new'
    call acf_write(file_info,atoms=atoms_p,strkey='posout')
    call atom_deallocate_old(atoms,sat=.true.,rat=.true.,bemoved=.true.)
    call atom_deallocate_old(atoms_p,sat=.true.,rat=.true.,bemoved=.true.)
    deallocate(ratall)
    !call deallocateatomsarrays
end subroutine genrandom
!*****************************************************************************************
