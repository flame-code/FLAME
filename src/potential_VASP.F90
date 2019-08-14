!*****************************************************************************************
subroutine init_potential_forces_vasp(atoms_t)
    use mod_potential, only: sat, atom_motion, cellvec, ntypat, stypat, &
        single_point_calculation, comment, comment2, vasp5, natarr
    use mod_atoms, only: typ_atoms, set_typat
    implicit none
    type(typ_atoms), intent(inout):: atoms_t
    !real(8):: cv_t(3,3)
    !integer:: nat
    !character(5):: sat_t(nat)
    !logical:: atom_motion_t(3,nat)
    !local variables
    !logical:: typeisnew
    integer:: itypat, iat
    !character(5):: stypat(128)
    if(atoms_t%nat==0) stop 'ERROR: nat=0 in init_potential_forces_vasp'
    if(atoms_t%nat>1000) stop 'ERROR: too many atoms in init_potential_forces_vasp'
    cellvec(1:3,1:3)=atoms_t%cellvec(1:3,1:3)
    atom_motion(1:3,1:atoms_t%nat)=atoms_t%bemoved(1:3,1:atoms_t%nat)
    vasp5=.true.
    if(allocated(sat)) then
        stop 'ERROR: sat already allocated in init_potential_forces_vasp'
    else
        allocate(sat(atoms_t%nat))
    endif
    !comment='Si'
    !comment2='Si'
    !comment='K Sr Si O'
    !comment2='K Sr Si O'
    do iat=1,atoms_t%nat
        sat(iat)=trim(atoms_t%sat(iat))
    enddo
    !ntypat=1
    !comment=trim(atoms_t%sat(1))
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
    !        comment=trim(comment)//" "//trim(atoms_t%sat(iat))
    !    else
    !        natarr(itypat)=natarr(itypat)+1
    !    endif
    !enddo
    call set_typat(atoms_t)
    ntypat=atoms_t%ntypat
    comment=''
    do itypat=1,atoms_t%ntypat
        comment=trim(comment)//" "//trim(atoms_t%stypat(itypat))
        natarr(itypat)=atoms_t%ltypat(itypat)
    enddo
    comment2=trim(comment)
    single_point_calculation=.true.
    !write(*,*) 'DEBUG: ntypat ',ntypat
    !do itypat=1,ntypat
    !    write(*,*) 'DEBUG: itypat,stypat(itypat)',itypat,trim(stypat(itypat)),natarr(itypat)
    !enddo
    !stop
end subroutine init_potential_forces_vasp
!*****************************************************************************************
subroutine cal_potential_forces_vasp(atoms)
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_potential, only: atom_motion, cellvec, vasp5, comment, comment2, &
        perfstatus, perfstatus_old, vasp_restart, fcalls, natarr, ntypat
    use mod_processors, only: iproc
    implicit none
    type(typ_atoms), intent(inout):: atoms
    !local variables
    character(6):: dir
    character(50):: filename
    character(50):: filename1
    character(50):: filename2
    character(100):: command
    !local varibles
    real(8):: stress(3,3)
    integer:: iat, itry
    logical:: success
    !ntypat=1
    !natarr(1)=atoms%nat
    !ntypat=4
    !natarr(1)=2
    !natarr(2)=10
    !natarr(3)=12
    !natarr(4)=36
    write(dir,'(a3,i3.3)') 'tmp',iproc
    filename=dir//'/POSCAR'
    call update_ratp(atoms)
    call write_poscar(filename,atoms%nat,atoms%ratp,cellvec,ntypat,natarr,comment,vasp5,comment2,atom_motion)
    !write(*,*) 'IPROC ',iproc
    !stop
    if(fcalls<1.d0 .or. trim(perfstatus)/=trim(perfstatus_old)) then
        vasp_restart=.false.
    else
        vasp_restart=.true.
    endif
    perfstatus_old=trim(perfstatus)
    do itry=1,2
    if(trim(perfstatus)=='fast') then
        command='cp -f INCAR.f '//dir//'/INCAR'
        call system(command)
    elseif(trim(perfstatus)=='normal') then
        command='cp -f INCAR.n '//dir//'/INCAR'
        call system(command)
    elseif(trim(perfstatus)=='accurate') then
        command='cp -f INCAR.a '//dir//'/INCAR'
        call system(command)
    elseif(trim(perfstatus)=='soften') then
        command='cp -f INCAR.s '//dir//'/INCAR'
        call system(command)
    else
        stop 'ERROR: perfstatus is not set properly.'
    endif
    if(vasp_restart) then
        command='sed -i "s/template1/ISTART=1/g" '//dir//'/INCAR'
        call system(command)
    else
        command='sed -i "s/template1/ISTART=0/g" '//dir//'/INCAR'
        call system(command)
    endif
    !command='cp KPOINTS '//dir//'/KPOINTS'
    !call system(command)
    call system("sleep 1")
    command='cd '//dir//'; ./vrun ; cd ..'
    call system(command)
    call system("sleep 1")
    filename1=dir//'/vasprun.xml'
    filename2=dir//'/CONTCAR'
    !write(*,*) filename
    call update_ratp(atoms)
    call get_output_vasp_geopt_alborz(filename1,filename2,atoms%nat,cellvec,atoms%ratp,atoms%fat,atoms%epot,stress,success)
    if(success) exit
    vasp_restart=.false.
    enddo !end of loop over itry
    if(.not. success) then
        write(*,'(a)') 'ERROR: failed to find all requested variables from '
        write(*,'(a,i3)') '       vasprun.xml by MHM iproc= ',iproc
        stop
    endif
    !call add_repulsive_wall(iproc,atoms%nat,atoms%rat,cellvec,atoms%fat,atoms%epot)
    !do iat=1,atoms%nat
    !    write(*,*) atoms%rat(1,iat),atoms%rat(2,iat),atoms%rat(3,iat)
    !enddo
end subroutine cal_potential_forces_vasp
!*****************************************************************************************
subroutine final_potential_forces_vasp
    use mod_potential, only: sat
    implicit none
    !local variables
    if(.not. allocated(sat)) then
        stop 'ERROR: sat not allocated in final_potential_forces_vasp'
    else
        deallocate(sat)
    endif
end subroutine final_potential_forces_vasp
!*****************************************************************************************
subroutine add_repulsive_wall(iproc,nat,rat,cellvec,fat,epot)
    use mod_potential, only: sat
    implicit none
    integer, intent(in):: iproc, nat
    real(8), intent(in):: rat(3,nat), cellvec(3,3)
    real(8), intent(inout):: fat(3,nat), epot
    !local variables
    integer:: iat
    real(8):: x, y, z, ttx, tty, ttz, rc, ampl
    real(8), allocatable:: rat_t(:,:), rat_i(:,:)
    allocate(rat_t(3,nat),rat_i(3,nat))
    call rxyz_cart2int_alborz(nat,cellvec,rat,rat_i)
    do iat=1,nat
        if(rat_i(1,iat)<0.d0) rat_i(1,iat)=rat_i(1,iat)+1.d0
        if(rat_i(2,iat)<0.d0) rat_i(2,iat)=rat_i(2,iat)+1.d0
        if(rat_i(3,iat)<0.d0) rat_i(3,iat)=rat_i(3,iat)+1.d0
        if(rat_i(1,iat)>1.d0) rat_i(1,iat)=rat_i(1,iat)-1.d0
        if(rat_i(2,iat)>1.d0) rat_i(2,iat)=rat_i(2,iat)-1.d0
        if(rat_i(3,iat)>1.d0) rat_i(3,iat)=rat_i(3,iat)-1.d0
        rat_t(1:3,iat)=matmul(cellvec,rat_i(1:3,iat))
    enddo
    rc=3.d0
    ampl=1.d0
    if(cellvec(1,1)<3.d0*rc) stop 'ERROR: rc in repulsive wall is too large'
    if(cellvec(2,2)<3.d0*rc) stop 'ERROR: rc in repulsive wall is too large'
    if(cellvec(3,3)<3.d0*rc) stop 'ERROR: rc in repulsive wall is too large'
    do iat=1,nat
        if(trim(sat(iat))/='Li') cycle
        x=rat_t(1,iat) ; y=rat_t(2,iat) ; z=rat_t(3,iat)
        if(x<rc .or. x>cellvec(1,1)-rc) then
            if(x>cellvec(1,1)-rc) x=x-cellvec(1,1)
            ttx=1.d0-(x/rc)**2
        else
            ttx=0.d0
        endif
        if(y<rc .or. y>cellvec(3,3)-rc) then
            if(y>cellvec(2,2)-rc) y=y-cellvec(2,2)
            tty=1.d0-(y/rc)**2
        else
            tty=0.d0
        endif
        if(z<rc .or. z>cellvec(3,3)-rc) then
            if(z>cellvec(3,3)-rc) z=z-cellvec(3,3)
            ttz=1.d0-(z/rc)**2
        else
            ttz=0.d0
        endif
        epot=epot+ampl*(ttx*3+tty**3+ttz**3)
        fat(1,iat)=fat(1,iat)+6.d0*ampl*x*ttx**2/rc**2
        fat(2,iat)=fat(2,iat)+6.d0*ampl*y*tty**2/rc**2
        fat(3,iat)=fat(3,iat)+6.d0*ampl*z*ttz**2/rc**2
    enddo
    deallocate(rat_t,rat_i)
end subroutine add_repulsive_wall
!*****************************************************************************************
!program test
!    implicit none
!    integer:: nat, iat
!    real(8):: cellvec(3,3), stress(3,3), rat(3,60), fat(3,60), epot
!    nat=60
!    call get_output_vasp_geopt_alborz(nat,cellvec,rat,fat,epot,stress)
!    write(*,'(f16.8)') epot
!    do iat=1,nat
!        write(21,'(3f20.16)') rat(1,iat),rat(2,iat),rat(3,iat)
!        write(22,'(3f17.8)') fat(1,iat),fat(2,iat),fat(3,iat)
!    enddo
!end program test
!*****************************************************************************************
!!!   subroutine vasp_geopt(latvec,xred,fcart,strten,energy,iprec,ka,kb,kc,counter)
!!!   !This routine will setup the input file for a vasp geometry optimization
!!!   !It will also call the run script and harvest the output
!!!   use global, only: nat
!!!   implicit none
!!!   real(8):: xred(3,nat),fcart(3,nat),strten(6),energy,counter,tmp
!!!   real(8):: dproj(6),acell(3),rprim(3,3),latvec(3,3)
!!!   integer:: iat,iprec,ka,kb,kc,itype
!!!   logical:: getwfk
!!!   character(4):: tmp_char
!!!   getwfk=.false.
!!!   !Set up the input file to perform geometry optimization
!!!    call make_input_vasp_geopt(latvec,xred,iprec,ka,kb,kc,getwfk)
!!!    call system("sleep 1")
!!!   !Run the job NOW!
!!!    call system("./runjob_geovasp.sh")
!!!    call system("sleep 1")
!!!   !Now harvest the structure, energy, forces, etc
!!!    call get_output_vasp_geopt_alborz(latvec,xred,fcart,energy,strten)
!!!   !Check how many iterations have been needed
!!!    call system("grep Conjugate OUTCAR_geo_a |wc -l>tmp_count")
!!!    call system("grep Conjugate OUTCAR_geo_b |wc -l>>tmp_count")
!!!    call system("grep Conjugate OUTCAR_geo_c |wc -l>>tmp_count")
!!!    call system("sleep 1")
!!!    open(unit=32,file="tmp_count")
!!!    read(32,*) tmp
!!!    counter=counter+tmp
!!!    read(32,*) tmp
!!!    counter=counter+tmp
!!!    read(32,*) tmp
!!!    counter=counter+tmp
!!!    close(32)
!!!    call system("rm -f tmp_count")
!!!   !We create a backup of the geometry optimization file
!!!    call system("cp vasprun.xml vasprun.xml.bak")
!!!   end subroutine
!!!   
!!!   subroutine make_input_vasp_geopt(latvec,xred,iprec,ka,kb,kc,getwfk)
!!!   !This routine will append some informations to a file already containing some informations about the abininit runs
!!!   !The informations appended are:
!!!   !-The atomic informations
!!!   !A file KPOINTS is generated. There are two options available:
!!!   !Automated generation if vasp_kpt_mode==1
!!!   !Monkhorst pack if vasp_kpt_mode==2
!!!   !ATTENTION:
!!!   !The meaning of dkpt1 and dkpt2 is different depending on vasp_kpt_mode:
!!!   !accuracy is given by the integer length of dkpt for vasp_kpt_mode==1 (10 for insulators, 100 for metals)
!!!   !accuracy is 2pi/bohr*dkpt for vasp_kpt_mode==2
!!!   use global, only: nat,ntypat,znucl,typat,dkpt1,dkpt2,char_type,ntime_geopt,tolmxf,target_pressure_gpa,vasp_kpt_mode
!!!   use defs_basis,only: Bohr_Ang
!!!   implicit none
!!!   real(8):: xred(3,nat)
!!!   real(8):: dproj(6),acell(3),rprim(3,3),latvec(3,3),dkpt,angbohr
!!!   integer:: iat,iprec,ka,kb,kc,itype,nat_type(ntypat)
!!!   logical:: getwfk
!!!   character(1):: fn
!!!   character(150):: command
!!!   angbohr=1.d0/Bohr_Ang
!!!   
!!!   getwfk=.false.
!!!   
!!!   if(iprec==1) then
!!!   dkpt=dkpt1
!!!   else
!!!   dkpt=dkpt2
!!!   endif
!!!   
!!!   
!!!   
!!!   command= "cp  vasprun.xml vasprun.xml.bak"
!!!   call system(command)
!!!   command= "rm -f INCAR POSCAR KPOINTS vasprun.xml OUTCAR_geo_a* OUTCAR_geo_b* OUTCAR_geo_c*"
!!!   call system(command)
!!!   write(fn,'(i1.1)') iprec
!!!   command = "cp -f INCAR."//fn//".a INCAR_geo_a"
!!!   call system(command)
!!!   command = "cp -f INCAR."//fn//".b INCAR_geo_b"
!!!   call system(command)
!!!   command = "cp -f INCAR."//fn//" INCAR_geo_c"
!!!   call system(command)
!!!   
!!!   open(unit=87,file="INCAR_geo_a",ACCESS="APPEND")
!!!   !Setup for only one force call
!!!   write(87,'(a)') ""
!!!   !Setup for only a sequence of geopt
!!!   write(87,'(a,i5)') "NSW = ",int(ntime_geopt*0.75d0)
!!!   write(87,'(a,es25.15)') "PSTRESS = ",target_pressure_gpa*10.d0
!!!   write(87,'(a,es25.15)') "EDIFFG = ",-tolmxf
!!!   write(87,'(a)') "IBRION = 2"
!!!   write(87,'(a)') "ISIF   = 3"
!!!   close(87)
!!!   
!!!   close(87)
!!!   
!!!   open(unit=87,file="INCAR_geo_b",ACCESS="APPEND")
!!!   !Setup for only one force call
!!!   write(87,'(a)') ""
!!!   !Setup for only a sequence of geopt
!!!   write(87,'(a,i5)') "NSW = ",int(ntime_geopt*0.25d0)
!!!   write(87,'(a,es25.15)') "PSTRESS = ",target_pressure_gpa*10.d0
!!!   write(87,'(a,es25.15)') "EDIFFG = ",-tolmxf
!!!   write(87,'(a)') "IBRION = 2"
!!!   write(87,'(a)') "ISIF   = 3"
!!!   close(87)
!!!   
!!!   open(unit=87,file="INCAR_geo_c",ACCESS="APPEND")
!!!   !Setup for only one force call
!!!   write(87,'(a)') ""
!!!   !Setup for only a sequence of geopt
!!!   write(87,'(a,i5)') "NSW = ",0
!!!   write(87,'(a,es25.15)') "PSTRESS = ",target_pressure_gpa*10.d0
!!!   write(87,'(a)') "IBRION = 2"
!!!   write(87,'(a)') "ISIF   = 2"
!!!   close(87)
!!!   
!!!   !Kpoint mesh
!!!   open(unit=87,file="KPOINTS")
!!!   write(87,'(a,i5)') "# Definition of the k-point mesh ",vasp_kpt_mode
!!!   write(87,'(i5)') 0
!!!   if(dkpt==0.d0) then
!!!   write(87,'(a)') "Gamma"!"Monkhorst Pack"
!!!   write(87,'(3(1x,i5),a)') ka,kb,kc,"  # Number of gridpoints in each dimension"
!!!   write(87,'(3(1x,i5),a)') 0,0,0,"  # Shifts"
!!!   elseif(vasp_kpt_mode==2) then
!!!   call find_kpt(ka,kb,kc,latvec,dkpt)
!!!   write(87,'(a)') "Gamma"!"Monkhorst Pack"
!!!   write(87,'(3(1x,i5),a)') ka,kb,kc,"  # Number of gridpoints in each dimension"
!!!   write(87,'(3(1x,i5),a)') 0,0,0,"  # Shifts"
!!!   elseif(vasp_kpt_mode==1) then
!!!   write(87,'(a)') "Auto"
!!!   write(87,'(i5,a)') dkpt," # K-mesh length"
!!!   else
!!!   stop "Wrong kpt option"
!!!   endif
!!!   write(*,'(a,3(1x,i5))') " # KPT mesh set up as follows: ",ka,kb,kc
!!!   close(87)
!!!   
!!!   
!!!   !The atomic positions and the unit cell is defined here
!!!   !The data are read from a file, called posinp.ascii
!!!   !The units are expected to be in angstroem
!!!   !This particular part is made for Silicon and Hydrogen cells
!!!   !In VASP, the atomic types need to follow in groups, and each number of
!!!   !a specific type must be given
!!!   nat_type=0
!!!   do itype=1,ntypat
!!!   do iat=1,nat
!!!   if(typat(iat)==itype) nat_type(itype)=nat_type(itype)+1
!!!   enddo
!!!   enddo
!!!   
!!!   
!!!   open(unit=87,file="POSCAR")
!!!   write(87,'(a)') "For geopt initially"
!!!   write(87,'(es15.7)') 1.d0
!!!   write(87,*)latvec(:,1)/angbohr
!!!   write(87,*)latvec(:,2)/angbohr
!!!   write(87,*)latvec(:,3)/angbohr
!!!   write(87,*) nat_type(:)
!!!   write(87,'(a)') "Direct"
!!!   do iat=1,nat
!!!   write(87,'(3(1x,es25.15),i5)') xred(:,iat)
!!!   enddo
!!!   close(87)
!!!   end

!*****************************************************************************************
subroutine get_output_vasp_geopt_alborz(filename1,filename2,nat,latvec,xred,fcart,energy,strten,success)
    !use defs_basis
    !Since its a single call, we only have forces and stresses from one configuration!
    use mod_potential, only: single_point_calculation
    implicit none
    character(*):: filename1
    character(*):: filename2
    integer:: nat
    real(8):: fcart(3,nat),energy,strten(6),value,latvec(3,3),xred(3,nat),str_matrix(3,3),vol,a(3,3),scaling
    logical:: success
    !local variables
    real(8):: xyz(3)
    integer:: ioserr,i,iat,n,k,l,m
    character(11):: ch_tmp
    character(150)::all_line
    energy=1.d10
    fcart=1.d10
    strten=1.d10
    ch_tmp="old"
    open(unit=32,file=filename1,status='old',iostat=ioserr)
    if(ioserr/=0) stop 'ERROR: failure openning vasprun.xml'
    do while(.true.)
        read(32,'(a150)',end=99)all_line
        !!write(*,*) all_line
        n = len_trim(all_line)
        !-------------------------------------------------------------
        !k = index(all_line(1:n),"rec_basis")
        !if(k.ne.0) then
        !    read(32,*)ch_tmp
        !    read(32,*)ch_tmp
        !    read(32,*)ch_tmp
        !    cycle
        !endif
        !-------------------------------------------------------------
        !k=index(all_line(1:n),"volumeweight")
        !if(k.ne.0) cycle
        !-------------------------------------------------------------
        !k=index(all_line(1:n),"basis")
        !if(k.ne.0) then
        !    !write(*,*) "Latvec found"
        !    !write(*,*) all_line(1:n)
        !    read(32,'(a150)',end=99)all_line
        !    m=len_trim(all_line)
        !    l=scan(all_line(1:m),"/",.true.)
        !    read(all_line(1:l-2),*) ch_tmp,latvec(:,1)
        !    read(32,'(a150)',end=99)all_line
        !    m=len_trim(all_line)
        !    l=scan(all_line(1:m),"/",.true.)
        !    read(all_line(1:l-2),*) ch_tmp,latvec(:,2)
        !    read(32,'(a150)',end=99)all_line
        !    m=len_trim(all_line)
        !    l=scan(all_line(1:m),"/",.true.)
        !    read(all_line(1:l-2),*) ch_tmp,latvec(:,3)
        !    !write(*,*) latvec
        !    cycle
        !endif
        !-------------------------------------------------------------
        !k = index(all_line(1:n),"positions")
        !if(k.ne.0) then
        !!write(*,*) "XRED found"
        !do iat=1,nat
        !read(32,*)ch_tmp,xred(:,iat)
        !enddo
        !!write(*,*) xred
        !endif
        !-------------------------------------------------------------
        !k=index(all_line(1:n),"volume")
        !if(k.ne.0) then
        !    l=VERIFY(all_line(1:n)," ",.false.)
        !    !write(*,*) "Volume found"
        !    !write(*,*) all_line(l:n)
        !    read(all_line(l:n),'(a17,F16.8,a)')ch_tmp,vol,ch_tmp
        !    !write(*,*) vol
        !endif
        !-------------------------------------------------------------
        !k=index(all_line(1:n),"e_fr_energy")
        k=index(all_line(1:n),"e_wo_entrp")
        if(k.ne.0) then
            l=VERIFY(all_line(1:n)," ",.false.)
            !write(*,*) "ETOT found"
              read(all_line(l:n),'(a22,F16.8,a)')ch_tmp,energy,ch_tmp
            !write(*,*) energy
        endif
        !-------------------------------------------------------------
        k=index(all_line(1:n),"forces")
        if(k.ne.0) then
            !write(*,*) "Forces found"
            do iat=1,nat
            read(32,'(a150)',end=99)all_line
            m=len_trim(all_line)
            l=scan(all_line(1:m),"/",.true.)
            read(all_line(1:l-2),*) ch_tmp,fcart(:,iat)
            enddo
            !write(*,*) xred
        endif
        !-------------------------------------------------------------
        !k=index(all_line(1:n),"stress")
        !if(k.ne.0) then
        !    !write(*,*) "stress found"
        !    read(32,'(a150)',end=99)all_line
        !    m = len_trim(all_line)
        !    l = scan(all_line(1:m),"/",.true.)
        !    read(all_line(1:l-2),*) ch_tmp,str_matrix(:,1)
        !    read(32,'(a150)',end=99)all_line
        !    m = len_trim(all_line)
        !    l = scan(all_line(1:m),"/",.true.)
        !    read(all_line(1:l-2),*) ch_tmp,str_matrix(:,2)
        !    read(32,'(a150)',end=99)all_line
        !    m = len_trim(all_line)
        !    l = scan(all_line(1:m),"/",.true.)
        !    read(all_line(1:l-2),*) ch_tmp,str_matrix(:,3)
        !    strten(1)=-str_matrix(1,1)
        !    strten(2)=-str_matrix(2,2)
        !    strten(3)=-str_matrix(3,3)
        !    strten(6)=-str_matrix(1,2)
        !    strten(5)=-str_matrix(1,3)
        !    strten(4)=-str_matrix(2,3)
        !    !write(*,*) strten
        !endif
        !-------------------------------------------------------------
    enddo
    99 continue
    close(32)
    !if(energy==1.d10.or.strten(1)==1.d10.or.fcart(1,1)==1.d10) stop "Could not find all requested variables"
    !if(energy==1.d10 .or.fcart(1,1)==1.d10) stop "Could not find all requested variables"
    if(energy==1.d10 .or.fcart(1,1)==1.d10) then
        success=.false.
    else
        success=.true.
    endif
    !if(target_pressure_gpa.ne.0.d0) then !commented by Alireza
    !    !In the vasprun.xml file at the end you will have the enthalpy instead of the total energy in the file, so
    !    !we need to transform it back, remember pressures are in kilobar in vasp
    !    !energy=energy-target_pressure_gpa*10.d0/1.60217733d-19/1.d22*vol
    !    energy=energy-target_pressure_gpa*10.d0/1.60217733d3*vol
    !    !write(*,*) "ETOT found"
    !    !write(*,*) energy
    !endif
    !-----------------------------------------------------------------
    !since in CONTCAR the cell and atomic positions are written with higher accuracy, get it from there:
    if(.not. single_point_calculation) then
        open(unit=32,file=filename2,status='old',iostat=ioserr)
        if(ioserr/=0) stop 'ERROR: failure openning CONTCAR'
        read(32,*)
        read(32,*) scaling
        read(32,*) latvec(:,1)
        read(32,*) latvec(:,2)
        read(32,*) latvec(:,3)
        latvec=latvec*scaling
        read(32,*)
        read(32,*)
        read(32,*)
        read(32,*)
        do iat=1,nat
            read(32,*) xyz(1),xyz(2),xyz(3)
            xred(1:3,iat)=matmul(latvec,xyz(1:3))
        enddo
        close(32)
    endif
    !-----------------------------------------------------------------
    !write(*,*) "LATVEC"
    !write(*,*) LATVEC
    !write(*,*) "xred"
    !write(*,*) xred
    !Transform all to bohr
    !latvec=latvec/Bohr_Ang !commented by Alireza
    !energy=energy/Ha_eV !commented by Alireza
    !strten=strten*0.1d0/HaBohr3_GPa !commented by Alireza
    !fcart=fcart/Ha_eV*Bohr_Ang !commented by Alireza
end subroutine get_output_vasp_geopt_alborz
!*****************************************************************************************


!subroutine write_kpoint(ka,kb,kc,latvec,dg)
!  !
!  implicit none
!  !--------------------------------------------------
!  character(30) :: system
!  character(30) :: car,tmp
!  integer       :: Ka,Kb,Kc,Kv,i,j
!  real(8)       :: Ra, Rb, Rc, Ga, Gb, Gc,cosine_alpha,cosine_beta,cosine_gamma
!  real(8)       :: angle_alpha,angle_beta,angle_gamma
!  real(8)       :: R(3,3),Rt(3,3),G(3,3),C(3),Volume,dg,Dv,Da,Db,Dc,latvec(3,3)
!  real(8), parameter     :: pi = 3.14159265358970
!  !--------------------------------------------------
!  !--------------------------------------------------
!  R(1,:)=latvec(:,1)
!  R(2,:)=latvec(:,2)
!  R(3,:)=latvec(:,3)
!
!
!  !- Real Lattice Parameters ------------------------
!  Ra=dsqrt(R(1,1)**2+R(1,2)**2+R(1,3)**2)
!  Rb=dsqrt(R(2,1)**2+R(2,2)**2+R(2,3)**2)
!  Rc=dsqrt(R(3,1)**2+R(3,2)**2+R(3,3)**2)
!
!  !- Cell Angles ------------------------------------
!  cosine_alpha=(R(2,1)*R(3,1)+R(2,2)*R(3,2)+R(2,3)*R(3,3))/Rb/Rc
!  cosine_beta=(R(1,1)*R(3,1)+R(1,2)*R(3,2)+R(1,3)*R(3,3))/Ra/Rc
!  cosine_gamma=(R(1,1)*R(2,1)+R(1,2)*R(2,2)+R(1,3)*R(2,3))/Ra/Rb
!  angle_alpha=(acos(cosine_alpha)/pi)*180.d0
!  angle_beta=(acos(cosine_beta)/pi)*180.d0
!  angle_gamma=(acos(cosine_gamma)/pi)*180.d0
!
!  !- Volme ------------------------------------------
!  call cross(R(2,1),R(2,2),R(2,3),R(3,1),R(3,2),R(3,3),C(1),C(2),C(3))
!  call dot(R(1,1),R(1,2),R(1,3),C(1),C(2),C(3),Volume)
!
!  !- Reciprocal Lattice -----------------------------
!  G(1,:)=2.0*pi*C(:)/Volume
!  call cross(R(3,1),R(3,2),R(3,3),R(1,1),R(1,2),R(1,3),C(1),C(2),C(3))
!  G(2,:)=2.0*pi*C(:)/Volume
!  call cross(R(1,1),R(1,2),R(1,3),R(2,1),R(2,2),R(2,3),C(1),C(2),C(3))
!  G(3,:)=2.0*pi*C(:)/Volume
!
!  !- Reciprocal Lattice Parameters ------------------
!  Ga=dsqrt(G(1,1)**2+G(1,2)**2+G(1,3)**2)
!  Gb=dsqrt(G(2,1)**2+G(2,2)**2+G(2,3)**2)
!  Gc=dsqrt(G(3,1)**2+G(3,2)**2+G(3,3)**2)
!
!  !- Actual spacing of Monkhorst-Pack grid ----------
!
!  !- Mesh parameters of Monkhorst-Pack grid ---------
!  call kmfind(dg,Ga,Da,Ka)
!  call kmfind(dg,Gb,Db,Kb)
!  call kmfind(dg,Gc,Dc,Kc)
!
!  !! debug
!  !write(*,*)"dg = ",dg
!  !write(*,*)Ka, Kb, Kc
!  !!
!
!end subroutine write_kpoint
!
!!--------------------------------------------------
!!  subroutine to perform dot-product of 2 vectors (In1 and In2) and write
!!       the results into the scalar Out
!!       written by Eric J. Simon, Cray Research, Inc. Fri Aug  2 1996
!!
!subroutine dot(Inx1, Iny1, Inz1,Inx2, Iny2, Inz2, outt)
!  implicit none
!  real(8) :: Inx1, Iny1, Inz1, Inx2, Iny2, Inz2, outt
!  outt = Inx1*Inx2 + Iny1*Iny2 + Inz1*Inz2
!  return
!end subroutine dot
!!---------------------------------------------------
!!  subroutine to perform cross-product of 2 vectors (In1 and In2) and write
!!       the results into a third vector (Out)
!!       written by Eric J. Simon, Cray Research, Inc. Fri Aug  2 1996
!!
!!         i       j       k
!!       In1x    In1y    In1z
!!       In2x    In2y    In2z
!!
!!       In1XIn2 = In1Y*In2z-In1z*In2y, In1z*In2x-In1x*In2z, In1x*In2y-In1y*In2x
!!
!subroutine cross(In1x, In1y, In1z, In2x, In2y, In2z, Outx, Outy, Outz)
!  implicit none
!  real(8) In1x, In1y, In1z, In2x, In2y, In2z, Outx
!  real(8) Outy, Outz
!
!  Outx = In1y*In2z - In2y*In1z
!  Outy = In2x*In1z - In1x*In2z
!  Outz = In1x*In2y - In2x*In1y
!
!  return
!end subroutine cross
!!----------------------------------------------------
!! Subroutine to find out the MP Kmesh accordding to one quality of k-point separation.
!!
!subroutine kmfind(di,Gi,Dd,Kd)
!  implicit none
!  integer :: Kd,i
!  real(8)  :: Gi,di,Dd
!  real(8), parameter :: pi = 3.14159265358970
!
!  Kd=int(Gi/di/2.0/pi)
!  if (Kd == 0) Kd = 1
!  Dd=Gi/Kd/2.0/pi
!  if (Dd.ge.di) then
!     do i=1,10
!        Kd=Kd+i
!        Dd=Gi/Kd/2.0/pi
!        if (Dd.le.di) exit
!     enddo
!  endif
!
!  return
!end subroutine kmfind
!!-------------------------------------------------
