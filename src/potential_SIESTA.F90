!*****************************************************************************************
subroutine init_cal_potential_forces_siesta(atoms_t)
    use mod_processorsm only: iproc, nproc
    use mod_potential, only: path1, path2, icount, perfstatus, cellvec
    use mod_atoms, only: typ_atoms
    use m_siesta_init
    !use m_siesta_analysis
    use m_siesta_move
    use m_siesta_end
    use m_siesta_forces
    use siesta_geom, only:na_s,na_u,xa,ucell,scell,isa !,cisa
    !use atm_types, only:species
    use atmfuncs, only:labelfis
    use parallel, only:node,nodes
    !USE m_steps, only: inicoor, fincoor
    use ifport
    implicit none
    type(typ_atoms), intent(in):: atoms_t
    !real(8):: cv_t(3,3)
    !integer:: nat, iproc, nproc
    !real(8):: rat(3,nat)
    !character(5):: sat(nat)
    !integer::iproc,nproc,nat,natmax,iat
    !real(8)::rat(3,natmax)
    !local variables
    character(50):: filename
    !logical:: typeisnew
    integer:: iat, istat, ioserr, i, itypat
    !integer:: ntypat
    !integer:: natarr(128)
    !character(5):: stypat(128)
    character(100):: command
    if(atoms_t%nat==0) stop 'ERROR: nat=0 in init_cal_potential_forces_siesta'
    cellvec(1:3,1:3)=atoms_t%cellvec(1:3,1:3)
    icount=0 !probably important for SIESTA
    write(path1,'(a3,i3.3)') 'tmp',iproc
    path2='..'
    if(trim(perfstatus)=='fast') then
        command='cp -f INPUT_DEBUG.f '//path1//'/INPUT_DEBUG'
        istat=system(command)
        if(istat/=0) stop 'ERROR: could not copy INPUT_DEBUG.f'
    elseif(trim(perfstatus)=='normal') then
        command='cp -f INPUT_DEBUG.n '//path1//'/INPUT_DEBUG'
        istat=system(command)
        if(istat/=0) stop 'ERROR: could not copy INPUT_DEBUG.n'
    elseif(trim(perfstatus)=='accurate') then
        command='cp -f INPUT_DEBUG.a '//path1//'/INPUT_DEBUG'
        istat=system(command)
        if(istat/=0) stop 'ERROR: could not copy INPUT_DEBUG.a'
    else
        stop 'ERROR: perfstatus is not set properly.'
    endif
    istat=system("sleep 1")
    if(istat/=0) stop 'ERROR: could not sleep for 1 second'
    istat=chdir(path1)
    !istat=getcwd(path)
    if(istat/=0) stop 'ERROR: could not change directory'
    !istat=system('pwd >dd')
    !stop
    !ntypat=1
    !stypat(1)=trim(atoms_t%sat(1))
    !natarr(1)=1
    !do iat=2,atoms_t%nat
    !    !write(*,*) 'TYPE ',trim(atoms_t%sat(iat)),trim(stypat(ntypat))
    !    typeisnew=.true.
    !    do itypat=1,ntypat
    !        if(trim(atoms_t%sat(iat))==trim(stypat(itypat))) then
    !            typeisnew=.false.
    !            exit
    !        endif
    !    enddo
    !    if(typeisnew) then
    !        ntypat=ntypat+1
    !        natarr(ntypat)=1
    !        stypat(ntypat)=trim(atoms_t%sat(iat))
    !    else
    !        natarr(itypat)=natarr(itypat)+1
    !    endif
    !enddo
    call set_typat(atoms_t)
    filename='posinp_siesta.xyz'
    open(unit=1358,file=filename,status='replace',iostat=ioserr)
    if(ioserr/=0) then
        write(*,*) 'ERROR: failure openning ',trim(filename)
        stop
    endif
    iat=0
    do itypat=1,atoms_t%ntypat
    do i=1,atoms_t%ltypat(itypat)
        iat=iat+1
        write(1358,'(3es24.15,i4)') atoms_t%rat(1,iat),atoms_t%rat(2,iat),atoms_t%rat(3,iat),itypat
    enddo
    enddo
    if(iat/=atoms_t%nat) stop 'ERROR: iat/=nat in SIESTA init_potential_forces'
    close(1358)

    call siesta_init()
    !stop
    !nat=na_s
    if(atoms_t%nat/=na_s) then
        write(*,*) 'ERROR: inconsistency between nat and SIESTA na_s '
    endif
    !rat(1:3,1:nat)=xa(1:3,1:nat)
    !iproc=node
    !nproc=nodes
    !cell(1)=scell(1,1)
    !cell(2)=scell(2,2)
    !cell(3)=scell(3,3)
    !if(cellvec(1:3,1:3)=scell(1:3,1:3)
    !write(*,*) 'CELL1 !',atoms_t%cellvec(1,1),atoms_t%cellvec(2,1),atoms_t%cellvec(3,1)
    !write(*,*) 'CELL1 !',atoms_t%cellvec(1,2),atoms_t%cellvec(2,2),atoms_t%cellvec(3,2)
    !write(*,*) 'CELL1 !',atoms_t%cellvec(1,3),atoms_t%cellvec(2,3),atoms_t%cellvec(3,3)
    !write(*,*) 'CELL2 ',scell(1,1),scell(2,1),scell(3,1)
    !write(*,*) 'CELL2 ',scell(1,2),scell(2,2),scell(3,2)
    !write(*,*) 'CELL2 ',scell(1,3),scell(2,3),scell(3,3)
    if(abs(scell(1,1)*0.529177d0-atoms_t%cellvec(1,1))>1.d-12 .or. & 
       abs(scell(2,1)*0.529177d0-atoms_t%cellvec(2,1))>1.d-12 .or. & 
       abs(scell(3,1)*0.529177d0-atoms_t%cellvec(3,1))>1.d-12 .or. & 
       abs(scell(1,2)*0.529177d0-atoms_t%cellvec(1,2))>1.d-12 .or. & 
       abs(scell(2,2)*0.529177d0-atoms_t%cellvec(2,2))>1.d-12 .or. & 
       abs(scell(3,2)*0.529177d0-atoms_t%cellvec(3,2))>1.d-12 .or. & 
       abs(scell(1,3)*0.529177d0-atoms_t%cellvec(1,3))>1.d-12 .or. & 
       abs(scell(2,3)*0.529177d0-atoms_t%cellvec(2,3))>1.d-12 .or. & 
       abs(scell(3,3)*0.529177d0-atoms_t%cellvec(3,3))>1.d-12) then
        stop 'ERROR: inconsistency between cell vectors: cellvec and SIESTA scell'
    endif
    !write(*,*) 'RAT ',rat(1,1),xa(1,1)
    do iat=1,atoms_t%nat
        !atoms_t%sat(iat)=trim(labelfis(isa(iat)))
        if(trim(atoms_t%sat(iat))/=trim(labelfis(isa(iat)))) then
            write(*,*) 'ERROR: inconsistency between sat and SIESTA labels ',iat
            stop
        endif
        if(abs(atoms_t%rat(1,iat)-xa(1,iat)*0.529177d0)>1.d-12 .or. &
           abs(atoms_t%rat(2,iat)-xa(2,iat)*0.529177d0)>1.d-12 .or. &
           abs(atoms_t%rat(3,iat)-xa(3,iat)*0.529177d0)>1.d-12) then
            write(*,*) 'ERROR: inconsistency between rat and SIESTA xa ',iat
            stop
        endif
        !atoms_t%sat(iat)=trim(cisa(iat))
        !atoms_t%sat(iat)=species(labelfis(isa(iat)))%label
        !if(iproc==0) write(*,*) 'SAT ',iat,trim(atoms_t%sat(iat))
    enddo
    if(iproc==0) write(*,*) 'na_s,na_u ',na_s,na_u
    if(iproc==0) write(*,*) 'node,nodes ',node,nodes
    if(iproc==0) write(*,*) 'nat ',atoms_t%nat
    !if(iproc==0) write(*,*) 'ucell ',ucell(1,1:3)
    !if(iproc==0) write(*,*) 'ucell ',ucell(2,1:3)
    !if(iproc==0) write(*,*) 'ucell ',ucell(3,1:3)
    !if(iproc==0) write(*,*) 'scell ',scell(1,1:3)
    !if(iproc==0) write(*,*) 'scell ',scell(2,1:3)
    !if(iproc==0) write(*,*) 'scell ',scell(3,1:3)
    !write(*,*) 'iproc,nproc ',node,nodes
    !call system('cd ..')
    !path='..'
    istat=chdir(path2)
    if(istat/=0) stop 'ERROR: could not change directory'
end subroutine init_cal_potential_forces_siesta
!*****************************************************************************************
subroutine cal_potential_forces_siesta(atoms)
    use mod_atoms, only: typ_atoms
    use mod_potential, only: path1, path2, icount, perfstatus
    use m_siesta_init
    !use m_siesta_analysis
    use m_siesta_move
    use m_siesta_end
    use m_siesta_forces !, only: siesta_forces, iscf_reza
    use m_energies, only: Etot
    use m_forces, only: fa
    use siesta_geom, only: xa
    use siesta_options, only: dDtol, g2cut, iscf_reza
    use ifport
    implicit none
    type(typ_atoms), intent(inout):: atoms
    !local variables
    integer:: istat, iat
    real(8):: drift(3), time1, time2
    !character(256):: path
    !call system('cd tmp000')
    !path='tmp000'
    if(trim(perfstatus)=='fast') then
        dDtol=1.d-3
    elseif(trim(perfstatus)=='normal') then
        dDtol=1.d-4
    elseif(trim(perfstatus)=='accurate') then
        dDtol=5.d-5
    else
        stop 'ERROR: perfstatus is not set properly.'
    endif
    istat=chdir(path1)
    if(istat/=0) stop 'ERROR: could not change directory'
    xa(1:3,1:n/3)=atoms%rat(1:3,1:n/3)/0.529177d0
    !write(*,*) 'IOnode ',IOnode
    call cpu_time(time1)
    call siesta_forces(icount)
    call cpu_time(time2)
    atoms%fat(1:3,1:n/3)=fa(1:3,1:n/3)*(13.60580d0/0.529177d0)
    write(*,'(a4,i3.3,2x,a10,2f10.2,i4,f10.1,es10.1)') 'REZA', &
        iproc,trim(perfstatus),time2-time1,(time2-time1)/real(iscf_reza,8),iscf_reza,g2cut,dDtol
    drift(1:3)=0.d0
    do iat=1,n/3
        drift(1:3)=drift(1:3)+atoms%fat(1:3,iat)
    enddo
    drift(1:3)=drift(1:3)/real(n/3,8)
    do iat=1,n/3
        atoms%fat(1:3,iat)=atoms%fat(1:3,iat)-drift(1:3)
    enddo
    atoms%epot=Etot*13.60580d0
    icount=icount+1
    !call system('cd ..')
    !path='..'
    istat=chdir(path2)
    if(istat/=0) stop 'ERROR: could not change directory'
end subroutine cal_potential_forces_siesta
!*****************************************************************************************
subroutine final_potential_forces_siesta
    use mod_potential, only: path1, path2
    use m_siesta_init
    !use m_siesta_analysis
    use m_siesta_move
    use m_siesta_end
    use m_siesta_forces
    use basis_types, only: basis_parameters
    use atm_types, only: species, elec_corr
    !USE m_steps, only: inicoor, fincoor
    use ifport
    implicit none
    integer:: istat
    !logical, external:: fdf_started
    !common fdf_started
    !data fdf_started /.false./

    !fdf_started=.false.
    !character(256):: path
    !call system('cd tmp000')
    !path='tmp000'
    istat=chdir(path1)
    if(istat/=0) stop 'ERROR: could not change directory'
    !call siesta_analysis( relaxd )
    call siesta_end()
    call fdf_shutdown
    if(allocated(basis_parameters)) deallocate(basis_parameters)
    if(allocated(species)) deallocate(species)
    if(allocated(elec_corr)) deallocate(elec_corr)
    !call system('cd ..')
    !path='..'
    istat=chdir(path2)
    if(istat/=0) stop 'ERROR: could not change directory'
end subroutine final_potential_forces_siesta
!*****************************************************************************************
