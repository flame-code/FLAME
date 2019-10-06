!*****************************************************************************************
subroutine cube_read(filename,atoms,poisson)
    use mod_atoms, only: typ_atoms, iatom_to_sat, atom_allocate_old, update_rat
    use mod_electrostatics, only: typ_poisson
    use dynamic_memory
    implicit none
    character(*), intent(in):: filename
    type(typ_atoms), intent(out):: atoms
    type(typ_poisson), intent(out):: poisson
    !integer, intent(in):: nx, ny, nz
    !real(8), intent(out)::
    !local variables
    real(8):: tt
    integer:: nat, iat, igpx, igpy, igpz, ind, ios, iline, iatom
    character(256):: str
    character(2):: separator
    integer:: istat
    separator=achar(32)//achar(9)
    !write(*,'(2a)') separator,'REZA'
    open(unit=1358,file=trim(filename),status='old',iostat=ios)
    if(ios/=0) then;write(*,'(2a)') 'ERROR: failure openning ',trim(filename);stop;endif
    read(1358,*)
    read(1358,*)
    read(1358,*) nat
    call atom_allocate_old(atoms,nat,0,0)
    read(1358,*) poisson%ngpx,poisson%hgrid(1,1),poisson%hgrid(2,1),poisson%hgrid(3,1)
    read(1358,*) poisson%ngpy,poisson%hgrid(1,2),poisson%hgrid(2,2),poisson%hgrid(3,2)
    read(1358,*) poisson%ngpz,poisson%hgrid(1,3),poisson%hgrid(2,3),poisson%hgrid(3,3)
    atoms%cellvec(1,1)=poisson%ngpx*poisson%hgrid(1,1)
    atoms%cellvec(2,1)=poisson%ngpx*poisson%hgrid(2,1)
    atoms%cellvec(3,1)=poisson%ngpx*poisson%hgrid(3,1)
    atoms%cellvec(1,2)=poisson%ngpy*poisson%hgrid(1,2)
    atoms%cellvec(2,2)=poisson%ngpy*poisson%hgrid(2,2)
    atoms%cellvec(3,2)=poisson%ngpy*poisson%hgrid(3,2)
    atoms%cellvec(1,3)=poisson%ngpz*poisson%hgrid(1,3)
    atoms%cellvec(2,3)=poisson%ngpz*poisson%hgrid(2,3)
    atoms%cellvec(3,3)=poisson%ngpz*poisson%hgrid(3,3)
    write(*,'(2a)') 'reading ',trim(filename)
    do iat=1,atoms%nat
        read(1358,*) iatom,atoms%zat(iat),atoms%ratp(1,iat),atoms%ratp(2,iat),atoms%ratp(3,iat)
        call iatom_to_sat(iatom,atoms%sat(iat))
        !if(iatom==1) atoms%sat(iat)='H'
        !if(iatom==6) atoms%sat(iat)='C'
        !if(iatom==7) atoms%sat(iat)='N'
        !if(iatom==29) atoms%sat(iat)='Cu'
    enddo
    call update_rat(atoms,upall=.true.)
    poisson%rho=f_malloc([1.to.poisson%ngpx,1.to.poisson%ngpy,1.to.poisson%ngpz], &
        id='poisson%rho')
    !allocate(poisson%rho(poisson%ngpx,poisson%ngpy,poisson%ngpz),stat=istat)
    !if(istat/=0) stop 'ERROR: allocation of rho failed.'
    !do igpx=1,poisson%ngpx
    !    do igpy=1,poisson%ngpy
    !        do igpz=1,poisson%ngpz

    !            str=''
    !            read(1358,'(a)') str
    !            read(str,*) rho
    !        enddo
    !    enddo
    !enddo
    igpx=1
    igpy=1
    igpz=0
    iline=0
    do
        iline=iline+1
        str=''
        read(1358,'(a)') str
        !str='  1.20645E-04   1.45834E-04   1.76737E-04   2.14873E-04   2.61824E-04   3.19073E-04'
        do
            str=adjustl(str)
            !write(*,'(a)') trim(str)
            read(str,*,iostat=ios) tt
            ind=scan(str,separator)
            if(ios==0) then
                igpz=igpz+1
                if(igpz>poisson%ngpz) then
                    igpz=1
                    igpy=igpy+1
                    if(igpy>poisson%ngpy) then
                        igpy=1
                        igpx=igpx+1
                    endif
                endif
                poisson%rho(igpx,igpy,igpz)=tt
            else
            endif
            !write(*,*) tt,ios,ind
            !write(*,'(es14.5)',advance='no') tt
            !write(*,*) ind
            if(ind>1) then
                str(1:ind)=''
            else
                exit
            endif
        enddo
        !write(*,*)
        !if(iline==1) exit
        if(igpx==poisson%ngpx .and. igpy==poisson%ngpy .and. igpz==poisson%ngpz) then
            exit
        endif
        !if(iline==poisson%ngpx*poisson%ngpy*poisson%ngpz) then
        !endif
    enddo
    !print*, " elec. and ionic charges: " ,  sum(rho)*hx*hy*hz , sum(q)
    close(1358)
end subroutine cube_read
!*****************************************************************************************
subroutine cube_write(filename,atoms,poisson,rho_or_pot)
    use mod_atoms, only: typ_atoms, sat_to_iatom, get_rat
    use mod_electrostatics, only: typ_poisson
    implicit none
    character(*), intent(in):: filename
    type(typ_atoms), intent(in):: atoms
    type(typ_poisson), intent(in):: poisson
    character(*), intent(in):: rho_or_pot
    !integer, intent(in):: nx, ny, nz
    !real(8), intent(out)::
    !local variables
    real(8):: tt, quantity
    integer:: iat, igpx, igpy, igpz, item, ios, iatom
    real(8), allocatable:: rat(:,:)
    open(unit=1358,file=trim(filename),status='replace',iostat=ios)
    if(ios/=0) then;write(*,'(2a)') 'ERROR: failure openning ',trim(filename);stop;endif
    write(*,'(2a)') 'writing ',trim(filename)
    write(1358,'(a)') ''
    write(1358,'(a)') ''
    write(1358,'(i5,3f13.6)') atoms%nat,0.d0,0.d0,0.d0
    write(1358,'(i5,3f13.6)') poisson%ngpx,poisson%hgrid(1,1),poisson%hgrid(2,1),poisson%hgrid(3,1)
    write(1358,'(i5,3f13.6)') poisson%ngpy,poisson%hgrid(1,2),poisson%hgrid(2,2),poisson%hgrid(3,2)
    write(1358,'(i5,3f13.6)') poisson%ngpz,poisson%hgrid(1,3),poisson%hgrid(2,3),poisson%hgrid(3,3)
    allocate(rat(3,atoms%nat))
    call get_rat(atoms,rat)
    do iat=1,atoms%nat
        !if(trim(atoms%sat(iat))=='H') iatom=1
        !if(trim(atoms%sat(iat))=='C') iatom=6
        !if(trim(atoms%sat(iat))=='N') iatom=7
        !if(trim(atoms%sat(iat))=='Cu') iatom=29
        call sat_to_iatom(atoms%sat(iat),iatom)
        write(1358,'(i5,f13.5,3f13.6)') iatom,atoms%zat(iat),rat(1,iat),rat(2,iat),rat(3,iat)
    enddo
    item=0
    do igpx=1,poisson%ngpx
        do igpy=1,poisson%ngpy
            do igpz=1,poisson%ngpz
                if(trim(rho_or_pot)=='pot') then
                    quantity=poisson%pot(igpx,igpy,igpz)
                else if(trim(rho_or_pot)=='rho') then
                    quantity=poisson%rho(igpx,igpy,igpz)
                else
                    stop 'ERROR: what should cube_write write into file?'
                endif
                if(abs(quantity)<1.d-90) quantity=0.d0
                write(1358,'(1x,es12.5,1x)',advance='no') quantity
                !write(1358,'(1x,f13.6,1x)',advance='no') quantity
                item=item+1
                if(item==6 .or. (igpx==poisson%ngpx .and. igpy==poisson%ngpy .and. igpz==poisson%ngpz)) then
                    item=0
                    write(1358,'(a1)') ' '
                endif
            enddo
        enddo
    enddo
    close(1358)
    deallocate(rat)
end subroutine cube_write
!*****************************************************************************************
