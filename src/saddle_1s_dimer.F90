!*****************************************************************************************
subroutine dimmethimproved(parini,iproc,atoms_s,nat,ndof,rat,epot,fat,curv,uvn,paropt)
    use mod_parini, only: typ_parini
    use mod_saddle, only: dimsep, nit, epotprime, &
        do_elim_trans, do_elim_rot, dmconverged, beta
    use mod_opt, only: typ_paropt
    use mod_atoms, only: typ_atoms, typ_file_info, atom_copy_old, atom_normalizevector
    use mod_atoms, only: atom_deallocate_old
    use mod_yaml_conf, only: write_yaml_conf
    use mod_atoms, only: atom_calnorm, set_rat
    use mod_utils
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iproc, nat, ndof !number coordinates
    type(typ_atoms), intent(inout):: atoms_s
    real(8), intent(inout):: rat(3,nat) !positions
    real(8), intent(out):: fat(3,nat) !forces
    real(8), intent(out):: uvn(3,nat) !unit vector along the dimer, \hat{n}.
    real(8), intent(out):: epot !potential energy 
    real(8), intent(out):: curv !curvature
    type(typ_paropt), intent(inout):: paropt
    !local variables
    type(typ_atoms):: atoms
    type(typ_file_info):: file_info
    integer:: istat, itt, iat
    real(8):: fnrm, curv0, DDOT
    !real(8), allocatable:: uvnt(:), uvnold(:)
    character(10):: tch
    character(10):: filename
    character(100):: comment
    !allocate(uvnt(ndof),stat=istat)
    !if(istat/=0) stop 'ERROR: failure allocating uvnold'
    if(iproc==0) then
        call yaml_map('dimsep',dimsep,fmt='(f15.5)')
        !write(*,'(a,f15.5)') 'dimsep ',dimsep
    endif
    nit=0
    call atom_copy_old(atoms_s,atoms,'atoms_s->atoms')
    !atoms%rat(1:3,1:atoms%nat)=rat(1:3,1:atoms%nat)
    call set_rat(atoms,rat,setall=.true.)
    call cal_potential_forces(parini,atoms)
    fat(1:3,1:atoms%nat)=atoms%fat(1:3,1:atoms%nat)
    epotprime=atoms%epot
    !-------------------------------
    if(trim(parini%rng_type)=='only_for_tests') then
        do iat=1,nat
            call random_number_generator_simple(uvn(1,iat))
            call random_number_generator_simple(uvn(2,iat))
            call random_number_generator_simple(uvn(3,iat))
        enddo
    else
        call random_number(uvn)
    endif
    uvn(1:3,1:nat)=uvn(1:3,1:nat)-0.5d0
    call elim_moment_alborz(nat,uvn)
    call elim_torque_reza_alborz(nat,rat,uvn)
    do iat=1,nat
        if(.not. atoms_s%bemoved(1,iat)) uvn(1,iat)=0.d0
        if(.not. atoms_s%bemoved(2,iat)) uvn(2,iat)=0.d0
        if(.not. atoms_s%bemoved(3,iat)) uvn(3,iat)=0.d0
    enddo
    call yaml_map('do_elim_trans',do_elim_trans)
    call yaml_map('do_elim_rot',do_elim_rot)
    !write(*,'(a,2l4)') 'do_elim_trans,do_elim_rot ',do_elim_trans,do_elim_rot
    !if(do_elim_trans) call projtransout(ndof,uvn)
    !if(do_elim_rot) call projrotout(3*nat,ndof,rat,uvn)
    !call normalizevector(ndof,uvn)
    !write(*,*) 'REZA-1 '
    call atom_normalizevector(nat,atoms_s%bemoved,uvn)
    !call normalizevector(3*nat,uvn)
    !write(*,*) 'REZA-2 '
    !uvn(1:ndof)=-uvn(1:ndof)
    call lowestcurvature(parini,iproc,atoms_s,nat,ndof,rat,uvn,fat,1.d0,1000,curv0,curv,1)
    !call calnorm(ndof,fat,fnrm)
    call atom_calnorm(nat,atoms_s%bemoved,fat,fnrm)
    if(iproc==0) then
        call yaml_mapping_open('initial iteration of lowest mode',flow=.true.)
        call yaml_map('epotprime',epotprime,fmt='(e16.8)')
        call yaml_map('fnrm',fnrm,fmt='(e16.8)')
        call yaml_mapping_close()
        !write(*,'(a7,2i5,3e16.8)') 'dmit   ',0,0,epotprime,0.d0,fnrm
    endif
    !if(curv<0.d0 .and. fnrm<fnrmtol) then
    !    dmconverged=.true.
    !    write(*,*) 'structure seems to be already at a saddle ponit.'
    !    write(*,*) 'epot,fnrm ',epotprime,fnrm
    !    stop !replace by exit later
    !endif
    file_info%filename_positions='posout.yaml'
    file_info%file_position='append'
    !call acf_write(file_info,atoms=atoms_s,strkey='posout')
    call write_yaml_conf(file_info,atoms_s,strkey='posout')

    !call escapeconvex(iproc,3*nat,ndof,rat,uvn,fat,curv0,curv)
    !beta=10.d0
    !paropt%fmaxtol=1.d-1
    !paropt%alpha=1.d-1*paropt%alphax
    !call optimizer_saddle(parini,iproc,3*nat,ndof,rat,fat,epot,paropt,uvn)
    !beta=4.d0
    !paropt%fmaxtol=1.d-3
    paropt%alpha=2.d-1*paropt%alphax
    call optimizer_saddle(parini,iproc,atoms_s,3*nat,ndof,rat,fat,epot,paropt,uvn)
    !call sdsaddle(iproc,3*nat,ndof,rat,uvn,fat,curv0,curv)
    !fnrmtol=1.d-4
    !call cgsaddle(iproc,3*nat,ndof,rat,uvn,fat,epot,curv0,curv)
    if(paropt%converged) then
        call lowestcurvature(parini,iproc,atoms_s,nat,ndof,rat,uvn,fat,5.d-1,20,curv0,curv,1)
        call yaml_map('lowestcurvatue',curv)
        !write(*,*) 'lowestcurvatue ',curv
        dmconverged=.true.
    endif
    !-------------------------------
    !call calnorm(ndof,fat,fnrm)
    !write(*,'(a15,2i5,e20.7)') 'DM finished.',nit,icount,fnrm
    call atom_deallocate_old(atoms)
end subroutine dimmethimproved
!*****************************************************************************************
subroutine lowestcurvature(parini,iproc,atoms_s,nat,ndof,rat,uvn,fat,angletol,maxitlc,curv0,curv,nw)
    use mod_parini, only: typ_parini
    use mod_saddle, only:nit,do_elim_trans,do_elim_rot
    use mod_atoms, only: typ_atoms, atom_ddot, atom_normalizevector
    use dynamic_memory
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iproc, nat, ndof, maxitlc, nw
    type(typ_atoms), intent(inout):: atoms_s
    real(8), intent(in):: rat(3,nat), fat(3,nat), angletol
    real(8), intent(inout):: uvn(3,nat), curv0, curv
    !local variables
    real(8), allocatable:: uvnold(:,:) 
    real(8), allocatable:: uvnt(:,:,:)
    integer:: itlc, istat, iw, jw
    real(8):: pi, DDOT, cosangletol, cosangle, t2, fnrm
    if(nw>10) stop 'ERROR: nw>10'
    uvnt=f_malloc([1.to.3,1.to.nat,1.to.10],id='uvnt')
    uvnold=f_malloc([1.to.3,1.to.nat],id='uvnold')
    !if(istat/=0) stop 'ERROR: failure allocating uvnold'
    !write(*,*) 'nit= ',nit
    if(nit==0) then
        uvnold=0.d0
    else
        uvnold(1:3,1:nat)=uvn(1:3,1:nat)
    endif
    !call normalizevector(nat,uvnold)
    pi=4.d0*atan(1.d0)
    cosangletol=cos(angletol*pi/180.d0)
    call yaml_sequence_open('Dimer lowest-mode optimization iterations')
    do iw=1,nw
    do itlc=1,maxitlc
        !call elim_moment_alborz(nat,uvn)
        !call elim_torque_reza_alborz(nat,rat,uvn)
        if(do_elim_trans) call elim_moment_alborz(nat,uvn)
        if(do_elim_rot) call elim_torque_reza_alborz(nat,rat,uvn)
        do jw=1,iw-1
            !t2=DDOT(nat,uvnt(1,jw),1,uvn,1)
            t2=atom_ddot(nat,uvnt(1,1,jw),uvn,atoms_s%bemoved)
            uvn(1:3,1:nat)=uvn(1:3,1:nat)-t2*uvnt(1:3,1:nat,jw)
        enddo
        call atom_normalizevector(nat,atoms_s%bemoved,uvn)
        !call normalizevector(3*nat,uvn)
        call rotatedimer(parini,iproc,atoms_s,nat,ndof,rat,uvn,fat,curv0,curv,fnrm)
        !cosangle=DDOT(3*nat,uvn,1,uvnold,1)
        cosangle=atom_ddot(nat,uvn,uvnold,atoms_s%bemoved)
        !write(*,'(a,i5,2es15.6,3e15.6)') &
        !    'DMROT,curv0,curv,cosangle,cosangletol,fnrm ', &
        !    itlc,curv0,curv,cosangle,cosangletol,fnrm
        call yaml_sequence(advance='no')
        call yaml_mapping_open('dimer',flow=.true.)
        call yaml_map('itlc',itlc,fmt='(i5)')
        call yaml_map('curv0',curv0,fmt='(es15.6)')
        call yaml_map('curv',curv,fmt='(es15.6)')
        call yaml_map('cosangle',cosangle,fmt='(e15.6)')
        call yaml_map('cosangletol',cosangletol,fmt='(e15.6)')
        call yaml_map('fnrm',fnrm,fmt='(e15.6)')
        call yaml_mapping_close()
        !cosangle=dot_product(uvn,uvnold)
        if(cosangle>cosangletol) exit
        uvnold(1:3,1:nat)=uvn(1:3,1:nat)
    enddo
    uvnt(1:3,1:nat,iw)=uvn(1:3,1:nat)
    enddo
    call yaml_sequence_close()
    call f_free(uvnold)
    call f_free(uvnt)
    nit=nit+1
end subroutine lowestcurvature
!*****************************************************************************************
subroutine rotatedimer(parini,iproc,atoms_s,nat,ndof,rat,uvn,fat,curv0,curv,fnrm)
    use mod_parini, only: typ_parini
    use mod_saddle, only: do_elim_rot, do_elim_trans, dimsep
    use mod_atoms, only: typ_atoms, atom_ddot, atom_copy_old, atom_normalizevector
    use mod_atoms, only: atom_calnorm, atom_deallocate_old, set_rat
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iproc, nat, ndof
    real(8), intent(in):: rat(3,nat), fat(3,nat)
    real(8), intent(inout):: uvn(3,nat), curv0, curv, fnrm
    type(typ_atoms), intent(inout):: atoms_s
    !local variables
    type(typ_atoms):: atoms
    integer:: istat, iat, ixyz, i
    real(8):: epotim1, phi1, phimin, a0, a1, b1, pi, dc0, cphi1, t1, DDOT
    real(8), allocatable:: xim1(:,:) !coordinates at image one.
    real(8), allocatable:: fim1(:,:) !forces at image one.
    real(8), allocatable:: fimdiff(:,:) !difference of forces between two images.
    real(8), allocatable:: uvp(:,:) !unit vector normal to uvn ie. \hat{\phi} 
    real(8), allocatable:: uvnphi1(:,:)
    call f_routine(id='rotatedimer')
    xim1=f_malloc([1.to.3,1.to.nat],id='xim1')
    fim1=f_malloc([1.to.3,1.to.nat],id='fim1')
    fimdiff=f_malloc([1.to.3,1.to.nat],id='fimdiff')
    uvp=f_malloc([1.to.3,1.to.nat],id='uvp')
    uvnphi1=f_malloc([1.to.3,1.to.nat],id='uvnphi1')
    call atom_copy_old(atoms_s,atoms,'atoms_s->atoms')
    xim1(1:3,1:nat)=rat(1:3,1:nat)
    do i=1,3*nat
        iat=(i-1)/3+1
        ixyz=mod(i-1,3)+1
        if(.not. atoms_s%bemoved(ixyz,iat)) then
            uvp(ixyz,iat)=0.d0
            uvnphi1(ixyz,iat)=0.d0
        endif
    enddo
    !write(*,*) 'dimsep',dimsep
    pi=4.d0*atan(1.d0)
    xim1(1:3,1:nat)=rat(1:3,1:nat)+dimsep*uvn(1:3,1:nat)
    !atoms%rat(1:3,1:nat)=xim1(1:3,1:nat)
    call set_rat(atoms,xim1,setall=.true.)
    call cal_potential_forces(parini,atoms)
    fim1(1:3,1:nat)=atoms%fat(1:3,1:nat)
    epotim1=atoms%epot
    !write(91+iproc,*) 'before ',epotim1
    fimdiff(1:3,1:nat)=2.d0*(fim1(1:3,1:nat)-fat(1:3,1:nat))
    if(do_elim_trans) call elim_moment_alborz(nat,fimdiff)
    if(do_elim_rot) call elim_torque_reza_alborz(nat,rat,fimdiff)
    !t1=DDOT(3*nat,fimdiff,1,uvn,1)
    t1=atom_ddot(nat,fimdiff,uvn,atoms_s%bemoved)
    !uvp(1:3,1:nat)=fimdiff(1:3,1:nat)-t1*uvn(1:3,1:nat)
    do iat=1,nat
        if(atoms_s%bemoved(1,iat)) uvp(1,iat)=fimdiff(1,iat)-t1*uvn(1,iat)
        if(atoms_s%bemoved(2,iat)) uvp(2,iat)=fimdiff(2,iat)-t1*uvn(2,iat)
        if(atoms_s%bemoved(3,iat)) uvp(3,iat)=fimdiff(3,iat)-t1*uvn(3,iat)
    enddo
    !call normalizevector(nat,uvp)
    call atom_normalizevector(nat,atoms_s%bemoved,uvp)
    !curv0=-0.5d0*dot_product(fimdiff,uvn)/dimsep
    !curv0=-0.5d0*DDOT(3*nat,fimdiff,1,uvn,1)/dimsep
    curv0=-0.5d0*atom_ddot(nat,fimdiff,uvn,atoms_s%bemoved)/dimsep
    !dc0=-dot_product(fimdiff,uvp)/dimsep
    !dc0=-DDOT(3*nat,fimdiff,1,uvp,1)/dimsep
    dc0=-atom_ddot(nat,fimdiff,uvp,atoms_s%bemoved)/dimsep
    phi1=0.25d0*pi 
    !phi1=0.5d0*atan(-0.5d0*dc0/curv0)
    !phi1=0.5d0*atan(-0.5d0*dc0/curv0)
    !phi1=max(phi1,0.25d0*pi)
    !uvnphi1(1:3,1:nat)=cos(phi1)*uvn(1:3,1:nat)+sin(phi1)*uvp(1:3,1:nat)
    do iat=1,nat
        if(atoms_s%bemoved(1,iat)) uvnphi1(1,iat)=cos(phi1)*uvn(1,iat)+sin(phi1)*uvp(1,iat)
        if(atoms_s%bemoved(2,iat)) uvnphi1(2,iat)=cos(phi1)*uvn(2,iat)+sin(phi1)*uvp(2,iat)
        if(atoms_s%bemoved(3,iat)) uvnphi1(3,iat)=cos(phi1)*uvn(3,iat)+sin(phi1)*uvp(3,iat)
    enddo
    xim1(1:3,1:nat)=rat(1:3,1:nat)+dimsep*uvnphi1(1:3,1:nat)
    !atoms%rat(1:3,1:nat)=xim1(1:3,1:nat)
    call set_rat(atoms,xim1,setall=.true.)
    call cal_potential_forces(parini,atoms)
    fim1(1:3,1:nat)=atoms%fat(1:3,1:nat)
    epotim1=atoms%epot
    !write(91+iproc,*) 'after  ',epotim1
    fimdiff(1:3,1:nat)=2.d0*(fim1(1:3,1:nat)-fat(1:3,1:nat))
    if(do_elim_trans) call elim_moment_alborz(nat,fimdiff)
    if(do_elim_rot) call elim_torque_reza_alborz(nat,rat,fimdiff)
    !cphi1=-0.5d0*DDOT(3*nat,fimdiff,1,uvnphi1,1)/dimsep
    cphi1=-0.5d0*atom_ddot(nat,fimdiff,uvnphi1,atoms_s%bemoved)/dimsep
    a1=(curv0-cphi1+0.5d0*dc0*sin(2.d0*phi1))/(1.d0-cos(2.d0*phi1))
    b1=0.5d0*dc0
    a0=2.d0*(curv0-a1)
    phimin=0.5d0*atan(b1/a1)
    !write(*,*) 'phimin',phimin
    curv=0.5d0*a0+a1*cos(2.d0*phimin)+b1*sin(2.d0*phimin)
    if(curv>curv0) phimin=phimin+0.5d0*pi
    curv=0.5d0*a0+a1*cos(2.d0*phimin)+b1*sin(2.d0*phimin)
    !------------------NEW------------------------
    t1=atom_ddot(nat,fimdiff,uvn,atoms_s%bemoved)
    do iat=1,nat
        if(atoms_s%bemoved(1,iat)) fimdiff(1,iat)=fimdiff(1,iat)-t1*uvn(1,iat)
        if(atoms_s%bemoved(2,iat)) fimdiff(2,iat)=fimdiff(2,iat)-t1*uvn(2,iat)
        if(atoms_s%bemoved(3,iat)) fimdiff(3,iat)=fimdiff(3,iat)-t1*uvn(3,iat)
    enddo
    call atom_calnorm(nat,atoms_s%bemoved,fimdiff,fnrm)
    !---------------------------------------------
    !uvn(1:3,1:nat)=cos(phimin)*uvn(1:3,1:nat)+sin(phimin)*uvp(1:3,1:nat)
    do iat=1,nat
        if(atoms_s%bemoved(1,iat)) uvn(1,iat)=cos(phimin)*uvn(1,iat)+sin(phimin)*uvp(1,iat)
        if(atoms_s%bemoved(2,iat)) uvn(2,iat)=cos(phimin)*uvn(2,iat)+sin(phimin)*uvp(2,iat)
        if(atoms_s%bemoved(3,iat)) uvn(3,iat)=cos(phimin)*uvn(3,iat)+sin(phimin)*uvp(3,iat)
    enddo
    call atom_deallocate_old(atoms)
    call f_free(xim1)
    call f_free(fim1)
    call f_free(fimdiff)
    call f_free(uvp)
    call f_free(uvnphi1)
    call f_release_routine()
end subroutine rotatedimer
!*****************************************************************************************
!subroutine escapeconvex(iproc,n,nr,x,uvn,f,curv0,curv)
!    use mod_saddle, only: nit, epotprime, ecconverged, maxitec
!    implicit none
!    integer::n,nr,istat,itec,i,j,nj,iproc
!    real(8)::x(n),uvn(nr),f(n),curv0,curv
!    real(8)::epot,fnrm,t1,dx,DDOT
!    real(8), allocatable::feff(:),xt(:),fold(:)
!    logical::goodescape
!    character(10):: filename
!    character(100):: comment
!    !write(*,*) 'n,nr',n,nr
!    allocate(feff(nr),stat=istat);if(istat/=0) stop 'error in allocating.'
!    allocate(xt(n),stat=istat);if(istat/=0) stop 'error in allocating.'
!    allocate(fold(n),stat=istat);if(istat/=0) stop 'error in allocating.'
!    !call cal_potential_forces(parini,iproc,n,x,f,epotprime)
!    xt(1:n)=x(1:n)
!    fold(1:n)=f(1:n)
!    do i=1,10
!    x(1:nr)=xt(1:nr)
!    f(1:nr)=fold(1:nr)
!    goodescape=.false.
!    do itec=1,maxitec
!        nit=nit+1
!        if(itec/=1 .or. i/=1) call lowestcurvature(iproc,n/3,nr,x,uvn,f,3.d0,10,curv0,curv,1)
!        !if(curv<-0.05d0) then
!        if(curv<0.d0) then
!            goodescape=.true.
!            exit
!        endif
!        !stop
!        !t1=dot_product(f,uvn)
!        t1=DDOT(nr,f,1,uvn,1)
!        feff(1:nr)=-t1*uvn(1:nr)
!        call normalizevector(nr,feff)
!        if(itec==1) then
!            !dx=2.d-2/abs(curv)
!            dx=4.d-2/abs(curv)
!            dx=sign(min(abs(dx),2.d0),dx)
!            !dx=sign(min(abs(dx),1.5d0),dx)
!        else
!            dx=5.d-3/abs(curv)
!            dx=sign(min(abs(dx),0.5d0),dx)
!            !dx=sign(min(abs(dx),0.25d0),dx)
!        endif
!        x(1:nr)=x(1:nr)+dx*feff(1:nr)
!        call cal_potential_forces(parini,iproc,n,x,f,epot)
!        call calnorm(nr,f,fnrm)
!        if(iproc==0) write(*,'(a7,2i5,3e16.8,3f15.8)') 'dmitec ',nit,itec,epot,(epot-epotprime),fnrm,curv0,curv,dx
!        !call wtposout_dim(n/3,nr/3,x,uvn,nit)
!        !write(filename,'(a10)') 'posout.xyz'
!        !write(comment,'(a11,es24.15,a7,i5)') ' angstroem ',0.d0,' iter= ',0
!        !call writexyz(filename,'append',n/3,x,comment)
!        !write(filename,'(a10)') 'modout.xyz'
!        !write(comment,'(a11,es24.15,a7,i5)') ' angstroem ',curv,' iter= ',0
!        !call writexyz(filename,'append',n/3,uvn,comment)
!        if(epot>epotprime+5.d0) exit
!    enddo
!    if(goodescape) then 
!            ecconverged=.true.
!        exit
!    endif
!    enddo
!    deallocate(feff,stat=istat);if(istat/=0) stop 'error in deallocating.'
!    deallocate(xt,stat=istat);if(istat/=0) stop 'error in deallocating.'
!    deallocate(fold,stat=istat);if(istat/=0) stop 'error in deallocating.'
!end subroutine escapeconvex
!*****************************************************************************************
