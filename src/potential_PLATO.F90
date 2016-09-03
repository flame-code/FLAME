!*****************************************************************************************
subroutine init_potential_forces_plato(atoms_t)
    use mod_interface
    use mod_processors, only: iproc
    use mod_potential, only: cellvec, path1, path2
    use mod_atoms, only: typ_atoms
    use ifport
    implicit none
    type(typ_atoms), intent(inout):: atoms_t
    !real(8):: cv_t(3,3)
    !integer:: nat, iproc
    !real(8):: rat(3,nat)
    !character(5):: sat(nat)
    !local variables
    character(50):: filename
    !integer::iproc,nproc,nat,natmax,iat
    !real(8)::rat(3,natmax)
    logical:: typeisnew
    integer:: iat, istat, ioserr, i, itypat, nat_p
    real(8):: cv_p(3,3), rat_p(3,1000), rat_frac(3,1000)
    real(8):: tt1, tt2, tt3
    !integer:: ntypat
    !integer:: natarr(128)
    !character(5):: sat(128)
    character(100):: command
    if(atoms_t%nat==0) stop 'ERROR: nat=0 in init_potential_forces_plato'
    if(atoms_t%nat>1000) stop 'ERROR: nat>1000 in init_potential_forces_plato'
    cellvec(1:3,1:3)=atoms_t%cellvec(1:3,1:3)
    !icount=0 !probably important for SIESTA
    write(path1,'(a3,i3.3)') 'tmp',iproc
    path2='..'
    !if(trim(perfstatus)=='fast') then
    !    command='cp -f INPUT_DEBUG.f '//path1//'/INPUT_DEBUG'
    !    istat=system(command)
    !    if(istat/=0) stop 'ERROR: could not copy INPUT_DEBUG.f'
    !elseif(trim(perfstatus)=='normal') then
        command='cp -f plato.in.n '//path1//'/plato.in'
        istat=system(command)
        if(istat/=0) stop 'ERROR: could not copy plato.in.n'
    !elseif(trim(perfstatus)=='accurate') then
    !    command='cp -f INPUT_DEBUG.a '//path1//'/INPUT_DEBUG'
    !    istat=system(command)
    !    if(istat/=0) stop 'ERROR: could not copy INPUT_DEBUG.a'
    !else
    !    stop 'ERROR: perfstatus is not set properly.'
    !endif
    !istat=system("sleep 1")
    !if(istat/=0) stop 'ERROR: could not sleep for 1 second'
    istat=chdir(path1)
    !istat=getcwd(path)
    if(istat/=0) stop 'ERROR: could not change directory'
    !ntypat=1
    !sat(1)=trim(sat(1))
    !natarr(1)=1
    !do iat=2,nat
    !    !write(*,*) 'TYPE ',trim(sat(iat)),trim(sat(ntypat))
    !    typeisnew=.true.
    !    do itypat=1,ntypat
    !        if(trim(sat(iat))==trim(sat(itypat))) then
    !            typeisnew=.false.
    !            exit
    !        endif
    !    enddo
    !    if(typeisnew) then
    !        ntypat=ntypat+1
    !        natarr(ntypat)=1
    !        sat(ntypat)=trim(sat(iat))
    !    else
    !        natarr(itypat)=natarr(itypat)+1
    !    endif
    !enddo
    call rxyz_cart2int_alborz(atoms_t%nat,cellvec,atoms_t%rat,rat_frac)
    do iat=1,atoms_t%nat
    if(rat_frac(1,iat)<0.d0) rat_frac(1,iat)=rat_frac(1,iat)+1.d0
    if(rat_frac(2,iat)<0.d0) rat_frac(2,iat)=rat_frac(2,iat)+1.d0
    if(rat_frac(3,iat)<0.d0) rat_frac(3,iat)=rat_frac(3,iat)+1.d0
    if(rat_frac(1,iat)>1.d0) rat_frac(1,iat)=rat_frac(1,iat)-1.d0
    if(rat_frac(2,iat)>1.d0) rat_frac(2,iat)=rat_frac(2,iat)-1.d0
    if(rat_frac(3,iat)>1.d0) rat_frac(3,iat)=rat_frac(3,iat)-1.d0
    atoms_t%rat(1:3,iat)=matmul(cellvec(1:3,1:3),rat_frac(1:3,iat))
    enddo
    filename='posinp_plato.xyz'
    open(unit=1358,file=filename,status='replace',iostat=ioserr)
    if(ioserr/=0) then
        write(*,*) 'ERROR: failure openning ',trim(filename)
        stop
    endif
    write(1358,'(a)') 'Atoms'
    do iat=1,atoms_t%nat
        write(1358,'(3f20.15,2x,a,2x)') rat_frac(1,iat),rat_frac(2,iat),rat_frac(3,iat),trim(atoms_t%sat(iat))
    enddo
    close(1358)
    filename='CELL_info'
    open(unit=1358,file=filename,status='replace',iostat=ioserr)
    if(ioserr/=0) then
        write(*,*) 'ERROR: failure openning ',trim(filename)
        stop
    endif
    write(1358,'(a)') 'CellVec'
    write(1358,'(3f25.15)') atoms_t%cellvec(1,1),atoms_t%cellvec(2,1),atoms_t%cellvec(3,1)
    write(1358,'(3f25.15)') atoms_t%cellvec(1,2),atoms_t%cellvec(2,2),atoms_t%cellvec(3,2)
    write(1358,'(3f25.15)') atoms_t%cellvec(1,3),atoms_t%cellvec(2,3),atoms_t%cellvec(3,3)
    write(1358,*) 
    write(1358,'(a)') 'CellSize'
    tt1=sqrt(cellvec(1,1)**2+cellvec(2,1)**2+cellvec(3,1)**2)/0.52917721d0;
    tt2=sqrt(cellvec(1,2)**2+cellvec(2,2)**2+cellvec(3,2)**2)/0.52917721d0;
    tt3=sqrt(cellvec(1,3)**2+cellvec(2,3)**2+cellvec(3,3)**2)/0.52917721d0;
    write(1358,'(3f25.15)') tt1,tt2,tt3
    close(1358)
        command='sed -i "/#template1/r CELL_info" plato.in'
        istat=system(command)
        if(istat/=0) stop 'ERROR: could not append cell information into plato.in'
        command='sed -i "/#template2/r posinp_plato.xyz" plato.in'
        istat=system(command)
        if(istat/=0) stop 'ERROR: could not append atomic position into plato.in'
        !stop

    call plato_initialize(nat_p,cv_p,rat_p)
    if(atoms_t%nat/=nat_p) then
        write(*,*) 'ERROR: inconsistency between nat and PLATO nat_p '
    endif
    !write(*,*) cv_p(1,1),atoms_t%cellvec(1,1),abs(cv_p(1,1)-atoms_t%cellvec(1,1))
    !write(*,*) cv_p(2,1),atoms_t%cellvec(2,1),abs(cv_p(2,1)-atoms_t%cellvec(2,1))
    !write(*,*) cv_p(3,1),atoms_t%cellvec(3,1),abs(cv_p(3,1)-atoms_t%cellvec(3,1))
    !write(*,*) cv_p(1,2),atoms_t%cellvec(1,2),abs(cv_p(1,2)-atoms_t%cellvec(1,2))
    !write(*,*) cv_p(2,2),atoms_t%cellvec(2,2),abs(cv_p(2,2)-atoms_t%cellvec(2,2))
    !write(*,*) cv_p(3,2),atoms_t%cellvec(3,2),abs(cv_p(3,2)-atoms_t%cellvec(3,2))
    !write(*,*) cv_p(1,3),atoms_t%cellvec(1,3),abs(cv_p(1,3)-atoms_t%cellvec(1,3))
    !write(*,*) cv_p(2,3),atoms_t%cellvec(2,3),abs(cv_p(2,3)-atoms_t%cellvec(2,3))
    !write(*,*) cv_p(3,3),atoms_t%cellvec(3,3),abs(cv_p(3,3)-atoms_t%cellvec(3,3))
    if(abs(cv_p(1,1)-atoms_t%cellvec(1,1))>1.d-12 .or. & 
       abs(cv_p(2,1)-atoms_t%cellvec(2,1))>1.d-12 .or. & 
       abs(cv_p(3,1)-atoms_t%cellvec(3,1))>1.d-12 .or. & 
       abs(cv_p(1,2)-atoms_t%cellvec(1,2))>1.d-12 .or. & 
       abs(cv_p(2,2)-atoms_t%cellvec(2,2))>1.d-12 .or. & 
       abs(cv_p(3,2)-atoms_t%cellvec(3,2))>1.d-12 .or. & 
       abs(cv_p(1,3)-atoms_t%cellvec(1,3))>1.d-12 .or. & 
       abs(cv_p(2,3)-atoms_t%cellvec(2,3))>1.d-12 .or. & 
       abs(cv_p(3,3)-atoms_t%cellvec(3,3))>1.d-12) then
        stop 'ERROR: inconsistency between cell vectors: cellvec and cv_p PLATO'
    endif
    !write(*,*) 'RAT ',atoms_t%rat(1,1),rat_p(1,1)
    !stop
    do iat=1,atoms_t%nat
        !if(trim(atoms_t%sat(iat))/=trim(labelfis(isa(iat)))) then
        !    write(*,*) 'ERROR: inconsistency between sat and PLATO labels ',iat
        !    stop
        !endif
        if(abs(atoms_t%rat(1,iat)-rat_p(1,iat))>1.d-10 .or. &
           abs(atoms_t%rat(2,iat)-rat_p(2,iat))>1.d-10 .or. &
           abs(atoms_t%rat(3,iat)-rat_p(3,iat))>1.d-10) then
            write(*,*) 'ERROR: inconsistency between rat and PLATO rat_p ',iat
            write(*,*) atoms_t%rat(1,iat),rat_p(1,iat),abs(atoms_t%rat(1,iat)-rat_p(1,iat))
            write(*,*) atoms_t%rat(2,iat),rat_p(2,iat),abs(atoms_t%rat(2,iat)-rat_p(2,iat))
            write(*,*) atoms_t%rat(3,iat),rat_p(3,iat),abs(atoms_t%rat(3,iat)-rat_p(3,iat))
            stop
        endif
        !atoms_t%sat(iat)=trim(cisa(iat))
        !atoms_t%sat(iat)=species(labelfis(isa(iat)))%label
        !if(iproc==0) write(*,*) 'SAT ',iat,trim(atoms_t%sat(iat))
    enddo
    if(iproc==0) write(*,*) 'nat ',atoms_t%nat
    istat=chdir(path2)
    if(istat/=0) stop 'ERROR: could not change directory'
end subroutine init_potential_forces_plato
!*****************************************************************************************
subroutine cal_potential_forces_plato(iproc,n,rat,fat,epot)
    use mod_interface
    use mod_potential, only: perfstatus, path1, path2
    use ifport
    implicit none
    integer, intent(in):: iproc, n
    real(8), intent(inout):: rat(3,n/3), fat(3,n/3), epot
    !local variables
    integer:: istat, iat, nat
    real(8):: drift(3), time1, time2, ResidueTol
    nat=n/3
    !character(256):: path
    !call system('cd tmp000')
    !path='tmp000'
    if(trim(perfstatus)=='fast') then
        ResidueTol=1.d-2
    elseif(trim(perfstatus)=='normal') then
        ResidueTol=5.d-3
    elseif(trim(perfstatus)=='accurate') then
        ResidueTol=1.d-3
    else
        stop 'ERROR: perfstatus is not set properly.'
    endif
    call set_plato_parameters(ResidueTol)
    istat=chdir(path1)
    if(istat/=0) stop 'ERROR: could not change directory'
    !xa(1:3,1:n/3)=rat(1:3,1:n/3)/0.529177d0
    call cpu_time(time1)
    !call siesta_forces(icount)
    call plato_energy_forces(nat,rat,epot,fat)
    call cpu_time(time2)
    !fat(1:3,1:n/3)=fa(1:3,1:n/3)*(13.60580d0/0.529177d0)
    write(*,'(a6,i3.3,2x,a10,f10.2,es10.0)') 'FCTIME',iproc,trim(perfstatus),time2-time1,ResidueTol
    drift(1:3)=0.d0
    do iat=1,n/3
        drift(1:3)=drift(1:3)+fat(1:3,iat)
    enddo
    write(*,'(a5,i3.3,2x,3es20.3)') 'DRIFT',iproc,drift(1),drift(2),drift(3)
    !drift(1:3)=drift(1:3)/real(n/3,8)
    !do iat=1,n/3
    !    fat(1:3,iat)=fat(1:3,iat)-drift(1:3)
    !enddo
    !epot=Etot*13.60580d0
    !icount=icount+1
    istat=chdir(path2)
    if(istat/=0) stop 'ERROR: could not change directory'
end subroutine cal_potential_forces_plato
!*****************************************************************************************
subroutine final_potential_forces_plato
    use mod_interface
    use ifport
    implicit none
    !integer:: istat
    !istat=chdir(path1)
    !if(istat/=0) stop 'ERROR: could not change directory'
    !istat=chdir(path2)
    !if(istat/=0) stop 'ERROR: could not change directory'
end subroutine final_potential_forces_plato
!*****************************************************************************************
