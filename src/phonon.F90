!*****************************************************************************************
subroutine cal_hessian_4p(parini)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_atoms_arr, typ_file_info
    use mod_processors, only: iproc
    use mod_potential, only: potential
    use futile
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    type(typ_atoms):: atoms
    type(typ_atoms_arr):: atoms_arr
    type(typ_file_info):: file_info
    real(8), allocatable:: hess(:,:)
    real(8), allocatable:: rat_center(:), eval(:), work(:)
    real(8):: h, rlarge, twelfth, twothird, shift, dm, tt, s, alpha
    integer:: istat, lwork, info, jj !, nat_yes
    integer:: i, iat, ixyz, j, jat, jxyz, imode, ii, iconf, iunit
    character(5):: fn
    call acf_read(parini,'posinp.acf',1,atoms=atoms)
    potential=trim(parini%potential_potential)
    !logical, allocatable:: yes(:)
    allocate(hess(3*atoms%nat,3*atoms%nat),stat=istat)
    allocate(rat_center(3*atoms%nat),stat=istat)
    allocate(eval(3*atoms%nat),stat=istat)
    lwork=100*atoms%nat
    allocate(work(lwork),stat=istat)
    !h=1.d-1
    !h=7.5d-2
    !h=5.d-2
    h=4.d-2
    !h=2.d-2
    rlarge=1.d0*1.d6
    twelfth=-1.d0/(12.d0*h)
    twothird=-2.d0/(3.d0*h)
    !if(iproc==0) write(*,*) 'HESSIAN: h',h
    !-------------------------------------------------------
    !if(iproc==0) write(*,'(a)') 'HESSIAN: reading from list_atoms.dat ...'
    !open(unit=1390,file='list_atoms.dat',status='old')
    !read(1390,*) nat_yes
    !allocate(yes(atoms%nat))
    !if(nat_yes==0) then
    !    if(iproc==0) write(*,'(a)') 'HESSIAN: a full-calculation of hessian.'
    !    yes(1:atoms%nat)=.true.
    !else
    !    if(iproc==0) write(*,'(a)') 'HESSIAN: a partial-calculation of hessian, '
    !    if(iproc==0) write(*,'(a)') 'HESSIAN: reading list of atoms.'
    !    yes(1:atoms%nat)=.false.
    !    do i=1,nat_yes
    !        read(1390,*) jj
    !        yes(jj)=.true.
    !    enddo
    !endif
    !close(1390)
    !-------------------------------------------------------
    !natsi=0
    !do iat=1,atoms%nat
    !    !if(iproc==0) write(*,*) 'GHASEMI: ',iat,trim(atoms%atomnames(atoms%iatype(iat)))
    !    if(trim(atoms%atomnames(atoms%iatype(iat)))=='Si') natsi=natsi+1
    !enddo
    !if(iproc==0) write(*,*) 'HESSIAN: natsi',natsi
    !hess(1:3*atoms%nat,1:3*atoms%nat)=0.d0
    !open(unit=1388,file='restart.dat',status='old')
    !read(1388,*) ii
    !if(ii>0) then
    !    if(iproc==0) write(*,'(a)') 'HESSIAN: it is a restart run.'
    !    do j=1,3*atoms%nat
    !    do i=1,3*atoms%nat
    !        read(1388,*) hess(i,j)
    !    enddo
    !    enddo
    !    ii=ii+1
    !elseif(ii==0) then
    !    if(iproc==0) write(*,'(a)') 'HESSIAN: it is a new run, '
    !    if(iproc==0) write(*,'(a)') 'HESSIAN: so it must be a full-calculation of hessian.'
    !    if(nat_yes/=0) stop 'HESSIAN: ERROR: a new run with nat_yes=0 '
    !    ii=1
    !else
    !    if(iproc==0) write(*,'(a)') 'HESSIAN: it is a new run, '
    !    if(iproc==0) write(*,'(a)') 'HESSIAN: but partial-calculation of hessian.'
    !    if(iproc==0) write(*,'(a)') 'HESSIAN: reading from restart_onlysi.dat ...'
    !    open(unit=1389,file='restart_onlysi.dat',status='old')
    !    read(1389,*)
    !    do j=1,3*natsi
    !    do i=1,3*natsi
    !        read(1389,*) hess(i,j)
    !    enddo
    !    enddo
    !    close(1389)
    !    ii=1
    !endif
    !close(1388)
    !-------------------------------------------------------
    !if(iproc==0) then
    !    do iat=1,atoms%nat
    !        write(*,'(a,i5,l3)') 'HESSIAN: ',iat,yes(iat)
    !    enddo
    !endif
    !-------------------------------------------------------
    call init_potential_forces(parini,atoms)
    do i=1,3*atoms%nat
        iat=(i-1)/3+1
        ixyz=mod(i-1,3)+1
        rat_center(i)=atoms%rat(ixyz,iat)
    enddo
    ii=1
    do i=ii,3*atoms%nat
        iat=(i-1)/3+1
        ixyz=mod(i-1,3)+1
        !if(.not. yes(iat)) cycle
        do j=1,3*atoms%nat
            jat=(j-1)/3+1
            jxyz=mod(j-1,3)+1
            atoms%rat(jxyz,jat)=rat_center(j)
        enddo
        !-----------------------------------------
        atoms%rat(ixyz,iat)=atoms%rat(ixyz,iat)-2*h
        call cal_potential_forces(parini,atoms)
        do j=1,3*atoms%nat
            jat=(j-1)/3+1
            jxyz=mod(j-1,3)+1
            hess(j,i)=twelfth*atoms%fat(jxyz,jat)
        enddo
        !-----------------------------------------
        atoms%rat(ixyz,iat)=atoms%rat(ixyz,iat)+h
        call cal_potential_forces(parini,atoms)
        do j=1,3*atoms%nat
            jat=(j-1)/3+1
            jxyz=mod(j-1,3)+1
            hess(j,i)=hess(j,i)-twothird*atoms%fat(jxyz,jat)
        enddo
        !-----------------------------------------
        atoms%rat(ixyz,iat)=atoms%rat(ixyz,iat)+2*h
        call cal_potential_forces(parini,atoms)
        do j=1,3*atoms%nat
            jat=(j-1)/3+1
            jxyz=mod(j-1,3)+1
            hess(j,i)=hess(j,i)+twothird*atoms%fat(jxyz,jat)
        enddo
        !-----------------------------------------
        atoms%rat(ixyz,iat)=atoms%rat(ixyz,iat)+h
        call cal_potential_forces(parini,atoms)
        do j=1,3*atoms%nat
            jat=(j-1)/3+1
            jxyz=mod(j-1,3)+1
            hess(j,i)=hess(j,i)-twelfth*atoms%fat(jxyz,jat)
        enddo
        !-----------------------------------------
        !if(iproc==0) then
        !    open(unit=1388,file='restart.dat',status='replace')
        !    write(1388,*) i
        !    do j=1,3*atoms%nat
        !    do k=1,3*atoms%nat
        !        write(1388,'(1es25.16)') hess(k,j)
        !    enddo
        !    enddo
        !    close(1388)
        !endif
        !-----------------------------------------
    enddo
    do j=1,3*atoms%nat
        jat=(j-1)/3+1
        jxyz=mod(j-1,3)+1
        atoms%rat(jxyz,jat)=rat_center(j)
    enddo
    call final_potential_forces(parini,atoms)
    !-------------------------------------------------------
    !deallocate(yes)
    !check symmetry
    dm=0.d0
    do i=1,3*atoms%nat
        do j=1,i-1
            s=.5d0*(hess(i,j)+hess(j,i))
            tt=abs(hess(i,j)-hess(j,i))/(1.d0+abs(s))
            dm=max(dm,tt)
            hess(i,j)=s
            hess(j,i)=s
        enddo
    enddo
    if(dm>1.d-1) write(*,*) 'max dev from sym',dm
    !do j=1,3*atoms%nat
    !    do i=1,3*atoms%nat
    !        write(201,*) hess(i,j)
    !    enddo
    !enddo
    !-------------------------------------------------------
    !project out rotations
    if(trim(atoms%boundcond)=='free') then
        call projectout_rotation(atoms,hess,rlarge,lwork,work)
    endif
    !-------------------------------------------------------
    !project out translations
    shift=rlarge/atoms%nat
    do j=1,3*atoms%nat-2,3
        do i=1,3*atoms%nat-2,3
            hess(i+0,j+0)=hess(i+0,j+0)+shift
            hess(i+1,j+1)=hess(i+1,j+1)+shift
            hess(i+2,j+2)=hess(i+2,j+2)+shift
        enddo
    enddo
    !-------------------------------------------------------
    !do j=1,3*atoms%nat
    !    do i=1,3*atoms%nat
    !        write(202,*) hess(i,j)
    !    enddo
    !enddo
    !-------------------------------------------------------
    !check
    call DSYEV('V','L',3*atoms%nat,hess,3*atoms%nat,eval,work,lwork,info)
    if(info/=0) stop 'DSYEV'
    if(iproc==0) then
        iunit=f_get_free_unit(10**5)
        open(unit=iunit,file='phonons.dat',status='replace')
        !write(iunit,*) '---  TB eigenvalues in a.u. -------------'
        atoms_arr%nconf=2*10+1
        allocate(atoms_arr%atoms(atoms_arr%nconf))
        do imode=1,3*atoms%nat
            write(fn,'(i5.5)') imode
            file_info%filename_positions='mode_'//fn//'.acf'
            file_info%file_position='new'
            file_info%print_force=.false.
            !tt=eval(imode)/(.529d0**2/27.2114d0)
            !write(iunit,'(a,i6,2e15.5)') 'eval (eV/A^2), a.u. ',imode,tt,eval(imode)
            !eval(imode)=tt
            write(iunit,'(i6,e15.5)') imode,eval(imode)
            if(trim(atoms%boundcond)=='free' .and. imode>3*atoms%nat-6) cycle
            if(trim(atoms%boundcond)/='free' .and. imode>3*atoms%nat-3) cycle
            iconf=0
            do j=-10,10
                iconf=iconf+1
                call atom_copy_old(atoms,atoms_arr%atoms(iconf),'atoms->atoms_arr%atoms(iconf)')
                alpha=j*5.d-2 !1.d-4/eval(imode)
                do iat=1,atoms%nat
                    atoms_arr%atoms(iconf)%rat(1,iat)=atoms%rat(1,iat)+alpha*hess(3*iat-2,imode)
                    atoms_arr%atoms(iconf)%rat(2,iat)=atoms%rat(2,iat)+alpha*hess(3*iat-1,imode)
                    atoms_arr%atoms(iconf)%rat(3,iat)=atoms%rat(3,iat)+alpha*hess(3*iat-0,imode)
                enddo
            enddo
            call acf_write_new(file_info,atoms_arr=atoms_arr,strkey='mode'//fn)
        enddo
        close(iunit)
        !open(unit=1359,file='eigenvectors.dat',status='replace')
        !write(1359,*) hess
        !close(1359)
        deallocate(atoms_arr%atoms)
    endif
end subroutine cal_hessian_4p
!*****************************************************************************************
subroutine projectout_rotation(atoms,hess,rlarge,lwork,work)
    use mod_interface
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(in):: atoms
    real(8), intent(inout):: hess(3*atoms%nat,3*atoms%nat)
    real(8), intent(in):: rlarge
    integer, intent(in):: lwork
    real(8), intent(inout):: work(lwork)
    !local variables
    real(8):: cmx, cmy, cmz
    real(8):: t1, t2, t3
    integer:: ixyz, iat, i, j
    cmx=0.d0 ; cmy=0.d0 ; cmz=0.d0
    do iat=1,atoms%nat
        cmx=cmx+atoms%rat(1,iat)
        cmy=cmy+atoms%rat(2,iat)
        cmz=cmz+atoms%rat(3,iat)
    enddo
    cmx=cmx/atoms%nat ; cmy=cmy/atoms%nat ; cmz=cmz/atoms%nat
    !x-y plane
    do i=1,3*atoms%nat-2,3
        iat=(i-1)/3+1
        ixyz=mod(i-1,3)+1
        work(i+1)= (atoms%rat(ixyz+0,iat)-cmx)
        work(i+0)=-(atoms%rat(ixyz+1,iat)-cmy)
    enddo
    t1=0.d0  ; t2=0.d0
    do i=1,3*atoms%nat-2,3
        t1=t1+work(i+0)**2
        t2=t2+work(i+1)**2
    enddo
    t1=sqrt(.5d0*rlarge/t1) ; t2=sqrt(.5d0*rlarge/t2)
    do i=1,3*atoms%nat-2,3
        work(i+0)=work(i+0)*t1
        work(i+1)=work(i+1)*t2
    enddo
    do j=1,3*atoms%nat-2,3
        do i=1,3*atoms%nat-2,3
            hess(i+0,j+0)=hess(i+0,j+0)+work(i+0)*work(j+0)
            hess(i+1,j+0)=hess(i+1,j+0)+work(i+1)*work(j+0)
            hess(i+0,j+1)=hess(i+0,j+1)+work(i+0)*work(j+1)
            hess(i+1,j+1)=hess(i+1,j+1)+work(i+1)*work(j+1)
        enddo
    enddo
    !x-z plane
    do i=1,3*atoms%nat-2,3
        iat=(i-1)/3+1
        ixyz=mod(i-1,3)+1
        work(i+2)= (atoms%rat(ixyz+0,iat)-cmx)
        work(i+0)=-(atoms%rat(ixyz+2,iat)-cmz)
    enddo
    t1=0.d0  ; t3=0.d0
    do i=1,3*atoms%nat-2,3
        t1=t1+work(i+0)**2
        t3=t3+work(i+2)**2
    enddo
    t1=sqrt(.5d0*rlarge/t1) ;  t3=sqrt(.5d0*rlarge/t3)
    do i=1,3*atoms%nat-2,3
        work(i+0)=work(i+0)*t1
        work(i+2)=work(i+2)*t3
    enddo
    do j=1,3*atoms%nat-2,3
        do i=1,3*atoms%nat-2,3
            hess(i+0,j+0)=hess(i+0,j+0)+work(i+0)*work(j+0)
            hess(i+2,j+0)=hess(i+2,j+0)+work(i+2)*work(j+0)
            hess(i+0,j+2)=hess(i+0,j+2)+work(i+0)*work(j+2)
            hess(i+2,j+2)=hess(i+2,j+2)+work(i+2)*work(j+2)
        enddo
    enddo
    !y-z plane
    do i=1,3*atoms%nat-2,3
        iat=(i-1)/3+1
        ixyz=mod(i-1,3)+1
        work(i+2)= (atoms%rat(ixyz+1,iat)-cmy)
        work(i+1)=-(atoms%rat(ixyz+2,iat)-cmz)
    enddo
    t2=0.d0 ; t3=0.d0
    do i=1,3*atoms%nat-2,3
        t2=t2+work(i+1)**2
        t3=t3+work(i+2)**2
    enddo
    t2=sqrt(.5d0*rlarge/t2) ; t3=sqrt(.5d0*rlarge/t3)
    do i=1,3*atoms%nat-2,3
        work(i+1)=work(i+1)*t2
        work(i+2)=work(i+2)*t3
    enddo
    do j=1,3*atoms%nat-2,3
        do i=1,3*atoms%nat-2,3
            hess(i+1,j+1)=hess(i+1,j+1)+work(i+1)*work(j+1)
            hess(i+2,j+1)=hess(i+2,j+1)+work(i+2)*work(j+1)
            hess(i+1,j+2)=hess(i+1,j+2)+work(i+1)*work(j+2)
            hess(i+2,j+2)=hess(i+2,j+2)+work(i+2)*work(j+2)
        enddo
    enddo
end subroutine projectout_rotation
!*****************************************************************************************
