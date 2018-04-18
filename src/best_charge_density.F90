subroutine best_charge_density(parini)
    use mod_interface
    use mod_parini, only: typ_parini
    implicit none
    type(typ_parini), intent(in):: parini

    call best_charge_density_rho(parini)
    !call best_charge_density_force(parini)
    !call best_charge_density_energy(parini)

end subroutine best_charge_density 
!**************************************************
!/////////////////////////////////////////////////////////////////////////////////////////
!/////////////////////////////////////////////////////////////////////////////////////////
!///////////////////////////////////// RHO_PART //////////////////////////////////////////
!/////////////////////////////////////////////////////////////////////////////////////////
!/////////////////////////////////////////////////////////////////////////////////////////
subroutine best_charge_density_rho(parini)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_cent, typ_ann_arr
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    type(typ_poisson):: poisson ,poisson_rough,poisson_dft, poisson_cent
    type(typ_atoms):: atoms
    type(typ_cent):: cent
    type(typ_ann_arr):: ann_arr
    integer:: atoms_type,istat, igpx, igpy, igpz, iat,iter,iter_max
    integer:: gweiat,i,j,k,jj,lcn,ierr,lwork, ss, Ne! lcn = linear combinarion number
    real(8):: cell(3), epot, rgcut_a, t1, t2, t3, t4, pi,h,rho_err
    real(8):: ehartree,rmse,SDA,dft_ener,cent_ener,ener_err,force_err,ener_err_old,coeff,w
    real(8),allocatable::  dft_rho(:,:,:),cent_rho(:,:,:),cent_rho_1(:,:,:),cent_rho_a_par(:,:,:),cent_rho_q_par(:,:,:),work(:),EE(:,:),eig(:)
    real(8),allocatable:: A(:,:), Q(:,:), At(:,:), Qt(:,:), qpar(:,:), apar(:,:), qpar_t(:,:),qpar_tmp(:,:), apar_t(:,:), qtemp(:), rat_t(:,:)
    real(8),allocatable:: dft_fat(:,:),cent_fat(:,:),gwit(:)
    real(8):: sd_s,err_fdm,err_cent,volinv,err
    real(8):: q_max,a_max,qnrm,dft_q,gtot
    real(8):: errmax, peak, c1, d1,temp1,peakx,peaky,peakz,errmax_old
    integer, allocatable :: t_num(:)
    logical :: esc
    real(8) :: at_zat , t
    integer :: at_range ,xg_at , yg_at , zg_at , n_at
    pi=4.d0*atan(1.d0)
    open(unit=1370,file='plot')
    open(unit=1,file='rhoerr')
    open(unit=2,file='errmax')
    open(unit=3,file='enererr')
    open(unit=4,file='forceerr')
    open(unit=5,file='peakloc')
    open(unit=1377,file='params.txt')
!//////////////////////////////////////READING INPUT PARAMETERS///////////////////////////
    call cube_read('electronic_density.cube',atoms,poisson_dft)
    call acf_read(parini,'posinp.acf',1,atoms=atoms)
    read(1377,*) !Number of STO-NG's N
    read(1377,*) lcn
    read(1377,*) !Number of atoms types
    read(1377,*) atoms_type
    allocate(t_num(0:atoms_type),Qt(lcn,atoms_type),At(lcn,atoms_type),rat_t(1:3,atoms_type),Q(1:lcn,1:atoms%nat),A(1:lcn,1:atoms%nat))
    t_num(0) = 0
    read(1377,*) !Number of each atom type 
    do i = 1 ,atoms_type 
        read(1377,*) t_num(i)
    end do
    read(1377,*) !Q's and alpha's of atom types
    do j = 1 , atoms_type
        do i = 1 ,lcn 
            read(1377,*) Qt(i,j) , At(i,j)
        end do
    end do
    read(1377,*) !Steepest descent step
    read(1377,*) sd_s
    read(1377,*) !Minimum error to be reached
    read(1377,*) err
    read(1377,*) !Number of Electrons for Single Atoms
    read(1377,*) Ne
    read(1377,*) !to be ploted atom's number
    read(1377,*) n_at
!/////////////////////////////////E.O. READING INPUT PARAMETERS///////////////////////////
    associate(nx=>poisson_dft%ngpx,ny=>poisson_dft%ngpy,nz=>poisson_dft%ngpz)
    associate(cv1=>atoms%cellvec(1,1),cv2=>atoms%cellvec(2,2),cv3=>atoms%cellvec(3,3))
    associate(hx=>poisson_dft%hx,hy=>poisson_dft%hy,hz=>poisson_dft%hz)
    associate(x_at=>atoms%rat(1,n_at),y_at=>atoms%rat(2,n_at),z_at=>atoms%rat(3,n_at))
    allocate(dft_rho(nx,ny,nz),dft_fat(1:3,1:atoms%nat),gwit(1:atoms%nat))
    ss=0
    do j = 1 , atoms_type
        ss = ss + t_num(j-1)
        rat_t(1:3,j) = atoms%rat(1:3,ss+1)
    end do
    volinv = (cv1*cv2*cv3)/(nx*ny*nz)
    xg_at = int(x_at/hx)
    yg_at = int(y_at/hy)
    zg_at = int(z_at/hz) 
    write(*,'(a,3i4)') "specified atom's grid",xg_at,yg_at,zg_at
    poisson_dft%rho=-1.d0*poisson_dft%rho
    !this is because get_psolver_bps changes poisson_dft%rho
    dft_rho = poisson_dft%rho
    peak = maxval(abs(dft_rho))
    w = 1.d0/(Ne+0.d0)
    gwit = 0.5d0
    read(1377,*)!atoms z_at
    at_range = atoms%nat/atoms_type
    do i = 1 , atoms_type
        read(1377,*) at_zat
        atoms%zat(((i-1)*at_range)+1:at_range*i) = at_zat ! this will not be OK
                                                          ! If a molecule
                                                          ! do not have same
                                                          ! number of each atom 
                                                          ! for example :
                                                          ! SrTiO_3
    end do
    write(*,*) "ZAT" , atoms%zat
    a_max = maxval(At)
!///////////////////////////////////////////DFT PART//////////////////////////////////////
    poisson_dft%rgcut=8.d0*a_max !ask Dr. Ghasemi if it's OK !
    write(*,*) "rgcut of dft" , poisson_dft%rgcut
    atoms%boundcond='bulk'   
    dft_fat = 0.d0
    dft_ener= 0.d0 
    atoms%fat = 0.d0 
    poisson_dft%task_finit=''
    call init_hartree(parini,atoms,poisson_dft,Q(1,:))
    !call init_psolver_bps(parini,atoms,cent%poisson)
    call get_hartree(parini,poisson_dft,atoms,Q(1,:),ehartree)
    !call get_psolver_bps(poisson_dft,atoms,ehartree)
    dft_ener = ehartree 
    if (atoms%nat > 1) then
        poisson_dft%bc='bulk'   
        poisson_dft%nat=atoms%nat
        poisson_dft%rcart=atoms%rat
        poisson_dft%cv=atoms%cellvec
        poisson_dft%q=atoms%zat
        poisson_dft%gw=gwit
        !call force_gto_sym(parini,poisson_dft%bc,poisson_dft%nat,poisson_dft%rcart, &
        !    poisson_dft%cv,poisson_dft%q,poisson_dft%gw,poisson_dft%rgcut,poisson_dft%ngpx, &
        !    poisson_dft%ngpy,poisson_dft%ngpz,poisson_dft%pot,atoms%fat)
        call get_hartree_force(parini,poisson_dft,atoms)
        !call force_gto_sym(parini,'bulk',atoms%nat,atoms%rat,atoms%cellvec,atoms%zat,gwit,poisson_dft%rgcut,nx,ny,nz,poisson_dft%pot,atoms%fat)
        dft_fat = atoms%fat
    end if
    write(21,*) atoms%nat 
    write(21,*) atoms%rat
    write(21,*) atoms%cellvec
    write(21,*) atoms%zat
    write(21,*) gwit
    write(21,*) poisson_dft%rgcut
    write(21,*) nx,ny,nz
    write(21,*) atoms%fat
    write(21,*) ehartree
    write(21,*) poisson_dft%pot(:,ny/2,nz/2) 
    write(21,*) "------------------------------------------------------------"
    write(*,*) "DFT1" , dft_ener , dft_fat
    dft_q = volinv*sum(dft_rho)
!//////////////////////////////////////////CENT PART//////////////////////////////////////
    poisson_cent%hx = hx
    poisson_cent%hy = hy
    poisson_cent%hz = hz
    poisson_cent%ngpx = nx
    poisson_cent%ngpy = ny
    poisson_cent%ngpz = nz
    poisson_cent%bc='bulk'   
    poisson_cent%nat=atoms%nat
    poisson_cent%rcart=atoms%rat
    poisson_cent%cv=atoms%cellvec
    allocate(cent_rho(nx,ny,nz))
    allocate(cent_rho_1(nx,ny,nz))
    allocate(cent_rho_a_par(nx,ny,nz))
    allocate(cent_rho_q_par(nx,ny,nz))
    allocate(qpar(1:lcn,1:atoms%nat),apar(1:lcn,1:atoms%nat),qpar_t(1:lcn,1:atoms_type),qpar_tmp(1:lcn,1:atoms_type),apar_t(1:lcn,1:atoms_type))
    allocate(cent_fat(1:3,1:atoms%nat))
!///////////////////////////////////TEST EIG PART/////////////////////////////////////////
  !eig_part!  !*************************eig_part********************************************
  !eig_part!  allocate(qtemp(atoms%nat))
  !eig_part!  qtemp(1:atoms%nat) = 1.d0
  !eig_part!  allocate(EE(atoms%nat,atoms%nat),eig(atoms%nat))
  !eig_part!  lwork = atoms%nat*atoms%nat + 100*atoms%nat
  !eig_part!  allocate(work(lwork))
  !eig_part!  EE(1:atoms%nat,1:atoms%nat) = 0.d0
  !eig_part!  eig(1:atoms%nat) = 0.d0
  !eig_part!  
  !eig_part!  do i = 1 , atoms%nat
  !eig_part!      call put_gto_sym(parini,atoms%boundcond,.true.,1,atoms%rat(1:3,i),atoms%cellvec,qtemp(i),A(1,i),cent%poisson%rgcut,nx,ny,nz,cent%poisson%rho)
  !eig_part!      cent_rho_1 = cent%poisson%rho
  !eig_part!      do j = i,atoms%nat
  !eig_part!          call put_gto_sym(parini,atoms%boundcond,.true.,1,atoms%rat(1:3,j),atoms%cellvec,qtemp(j),A(1,j),cent%poisson%rgcut,nx,ny,nz,cent%poisson%rho)
  !eig_part!          cent_rho_2 = cent%poisson%rho
  !eig_part!          EE(i,j) = 2.d0*invvol*sum(cent_rho_1*cent_rho_2)
  !eig_part!      end do
  !eig_part!  end do
  !eig_part!   
  !eig_part!  call DSYEV('N','U',atoms%nat,EE,atoms%nat,eig,work,lwork,ierr)
  !eig_part!  do i =1 , atoms%nat
  !eig_part!      write(*,'(a,8f13.8)')"EE ",EE(i,:)
  !eig_part!  end do
  !eig_part!  write(*,'(a,i4,i4,8f13.8)')"eig ",iter,ierr,eig(1:atoms%nat) 
    !*********************eig_part************************************************
!/////////////////////////////////E.O. TEST EIG PART//////////////////////////////////////
    allocate(qtemp(1:atoms%nat))
    qtemp = 1.d0

    do i = 1 , lcn
        ss=0
        do j = 1 , atoms_type
            ss = ss + t_num(j-1)
            Q(i,ss+1:t_num(j)+ss) = Qt(i,j)
            A(i,ss+1:t_num(j)+ss) = At(i,j)
        end do
    end do
    
    qnrm = (dft_q+0.d0)/(sum(Q(1:lcn,1:atoms%nat)))
    Q(1:lcn,1:atoms%nat) = Q(1:lcn,1:atoms%nat)*qnrm
    
    poisson_cent%task_finit='alloc_rho'
    call init_hartree(parini,atoms,poisson_cent,Q(1,:))
    poisson_cent = poisson_dft
    do i = 1 , lcn
        ss = 0
        do j = 1 , atoms_type
            ss = ss + t_num(j-1)
            qt(i,j)=q(i,1+ss)
            at(i,j)=a(i,1+ss)
        end do
    end do
    do iter=1,huge(iter_max)
        a_max = maxval(at)
        poisson_cent%rgcut=8.d0*a_max !ask Dr. Ghasemi if it's OK !
        write(*,*) "rgcut of cent" , poisson_cent%rgcut
        !poisson_cent%rho = 0.d0
        do i = 1 , lcn
            poisson_cent%q=Q(i,1:atoms%nat)
            poisson_cent%gw=A(i,1:atoms%nat)
            !poisson_cent%basis_set='rp4'
            if (i==1) then         
                !poisson_cent%reset_rho=.true.
                !call put_charge_density(parini,poisson_cent)
                call put_rp4gto_sym(parini,atoms%boundcond,.true.,atoms%nat,atoms%rat,atoms%cellvec,Q(i,1:atoms%nat),A(i,1:atoms%nat),poisson_cent%rgcut,nx,ny,nz,poisson_cent%rho,cent_rho_q_par,cent_rho_a_par)
            else
                !poisson_cent%reset_rho=.false.
                !call put_charge_density(parini,poisson_cent)
                call put_rp4gto_sym(parini,atoms%boundcond,.false.,atoms%nat,atoms%rat,atoms%cellvec,Q(i,1:atoms%nat),A(i,1:atoms%nat),poisson_cent%rgcut,nx,ny,nz,poisson_cent%rho,cent_rho_q_par,cent_rho_a_par)
            end if
        end do
        cent_rho = poisson_cent%rho
        inquire(file='EXIT',exist=esc)
        if(esc .and. iter >2) exit !to be moved after calling gauess grid
        write(*,*) "cent_rho vs. dft_rho" , volinv*sum(cent_rho) , volinv*sum(dft_rho)
        atoms%fat = 0.d0
        !cent_fat = 0.d0
        !call init_psolver_bps(parini,atoms,cent%poisson)
        !call get_psolver_bps(cent%poisson,atoms,ehartree)
        call get_psolver(parini,poisson_cent,atoms,Q(1,:),ehartree)
        cent_ener = ehartree
        
        poisson_cent%q=atoms%zat
        poisson_cent%gw=gwit
        if(atoms%nat > 1)  then
            call get_hartree_force(parini,poisson_cent,atoms)
            !call force_gto_sym(parini,'bulk',atoms%nat,atoms%rat,atoms%cellvec,atoms%zat,cent%gwit,cent%poisson%rgcut,nx,ny,nz,cent%poisson%pot,atoms%fat)
            cent_fat = atoms%fat
        end if
        !stop('AAAAAAAAAAA')
        !call fini_psolver_bps(cent%poisson)
        write(66,'(i4, a)') iter ,"  -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
        do i = 1 , atoms%nat
            do j =1 , 3 
                write(66,'(i4,3es14.6)') i,dft_fat(j,i),abs(dft_fat(j,i)-cent_fat(j,i)),cent_fat(j,i)-dft_fat(j,i)
            end do
        end do
        do i = 1 , atoms_type
            do j = 1 , lcn
            !poisson_cent%nat=1
            !poisson_cent%q=Qt(j,i)
            !poisson_cent%gw=At(j.i)
            !poisson_cent%basis_set='rp4'
            !poisson_cent%reset_rho=.true.
            !call put_charge_density(parini,poisson_cent)
                call put_rp4gto_sym(parini,atoms%boundcond,.true.,1,rat_t(1:3,i),atoms%cellvec,Qt(j,i),At(j,i),poisson_cent%rgcut,nx,ny,nz,cent_rho_1, cent_rho_q_par, cent_rho_a_par)
                cent_rho_1 = 0.d0
                qpar_t(j,i) = 2.d0*volinv*sum(cent_rho_q_par*(cent_rho-dft_rho))
                apar_t(j,i) = 2.d0*volinv*sum(cent_rho_a_par*(cent_rho-dft_rho))
                write(*,'(a,i4,i4,3es14.6)') "par_t typ lcn qpar apar",i,j,qpar_t(j,i),apar_t(j,i)!,qpar_tmp(j,i)
                write(*,*) "par_t ","___________"
            end do 
            write(*,*) "par_t ------------------------------------------------"
        end do
        do i = 1 , lcn
            ss = 0
            do j = 1 , atoms_type
                ss = ss + t_num(j-1)
                qpar(i,1+ss:t_num(j)+ss) = qpar_t(i,j)
                apar(i,1+ss:t_num(j)+ss) = apar_t(i,j)
            end do
        end do
        ener_err = abs(cent_ener-dft_ener)
        force_err = sqrt(sum((cent_fat-dft_fat)**2)/(3*atoms%nat))
        errmax = 0.d0
        rho_err = 0.d0
        write(*,*) "peak",peak
        do i = 1 , nz
            do j = 1 , ny
                do k = 1 , nx
                    c1 = cent_rho(k,j,i)
                    d1 = dft_rho(k,j,i)
                    temp1 = (abs(c1)-abs(d1))/(max(abs(d1),(1.d-2*peak)))
                    !write(21,'(a,10es14.5)')"temp ",temp1,abs(c1),abs(d1),(1.d-2*peak),k*hx,j*hy,i*hz,atoms%rat
                    errmax_old = errmax
                    errmax=max(errmax,temp1)
                    if(errmax_old /= errmax) then
                        peakx = k*hx
                        peaky = j*hy
                        peakz = i*hz
                    end if
                    rho_err = rho_err + ((c1-d1)**2)
                end do
            end do
        end do
        write(5,'(i8,3es14.6)')iter,peakx,peaky,peakz
        rho_err = sqrt(rho_err*volinv)*w
        !write(40,*)iter,rho_err,volinv,w
        rmse = sqrt(rho_err)
        !stop 'AAAAAAAAAAAA'
        write(*,'(a,i4,4es14.6)')"SD",iter,rmse,errmax,ener_err,force_err
        write(1,*)iter,rmse
        write(2,*)iter,errmax
        write(3,*)iter,ener_err
        write(4,*)iter,force_err
        
        write(*,*)"SD --------------------------------------------------------------------"
        write(*,'(a,i4,i4,es14.6)')"GTO_val_sum ",iter,i,sum(q(:,1))
        do j = 1 , atoms_type
            do i = 1 , lcn
                write(*,'(a,i4,i4,i4,2es14.6)')"GTO_val ",iter,i,j,qt(i,j),at(i,j)
                write(*,'(a,i4,i4,i4,2es14.4)') "a-q-partial",iter,i,j,qpar_t(i,j),apar_t(i,j) 
            end do
        end do
        write(*,'(a)')"GTO --------------------------------------------------------------------"
        if (abs(errmax) < err)then
            write(*,*) "max conversion reached"
            exit
        endif
        gtot = sum(qpar)
        qpar = qpar - (gtot+0.d0)/(atoms%nat*lcn+0.d0)
        Q = Q - SD_S*qpar
        A = A - SD_S*apar
        do i = 1 , lcn
            ss = 0
            do j = 1 , atoms_type
                ss = ss + t_num(j-1)
                qt(i,j)=q(i,1+ss)
                at(i,j)=a(i,1+ss)
            end do
        end do
    end Do
!/////////////////////////////////E.O. CENT PART//////////////////////////////////////////
write(*,*) "END OF CENT CALCULATIONS"
!//////////////////////////Specified Atoms charge density part////////////////////////////
    !poisson_cent%nat=1
    !poisson_cent%q=Q(i,n_at)
    !poisson_cent%gw=At(i,n_at)
    !poisson_cent%basis_set='rp4'
    do i = 1 , lcn
        if (i==1) then         
            !poisson_cent%reset_rho=.true.
            !call put_charge_density(parini,poisson_cent)
            call put_rp4gto_sym(parini,atoms%boundcond,.true.,1,atoms%rat(1:3,n_at),atoms%cellvec,Q(i,n_at),A(i,n_at),poisson_cent%rgcut,nx,ny,nz,poisson_cent%rho,cent_rho_q_par,cent_rho_a_par)
        else
            !poisson_cent%reset_rho=.false.
            !call put_charge_density(parini,poisson_cent)
            call put_rp4gto_sym(parini,atoms%boundcond,.false.,1,atoms%rat(1:3,n_at),atoms%cellvec,Q(i,n_at),A(i,n_at),poisson_cent%rgcut,nx,ny,nz,poisson_cent%rho,cent_rho_q_par,cent_rho_a_par)
            
        end if
    end do
    write(*,*) "END OF SPECIFIED ATOM'S PLOT WRITE"
!//////////////////////////E.O. Specified Atom's charge density part//////////////////////
!////Linear combination of dft and cent charge density to investigate the force_err///////
    poisson_cent%rho = 0.d0
    poisson_cent%pot = 0.d0
    do i = 1 , lcn
        if (i==1) then         
            call put_rp4gto_sym(parini,atoms%boundcond,.true.,atoms%nat,atoms%rat,atoms%cellvec,Q(i,1:atoms%nat),A(i,1:atoms%nat),poisson_cent%rgcut,nx,ny,nz,poisson_cent%rho,cent_rho_q_par,cent_rho_a_par)
        else
            call put_rp4gto_sym(parini,atoms%boundcond,.false.,atoms%nat,atoms%rat,atoms%cellvec,Q(i,1:atoms%nat),A(i,1:atoms%nat),poisson_cent%rgcut,nx,ny,nz,poisson_cent%rho,cent_rho_q_par,cent_rho_a_par)
        end if
    end do
    cent_rho = poisson_cent%rho
    do i = 1 , 9
        poisson_cent%rho = 0.d0
        poisson_cent%pot = 0.d0
        t = i*0.1 
        atoms%fat = 0.d0
        cent_fat = 0.d0
        poisson_cent%rho=cent_rho-(t*(cent_rho-dft_rho))
        !call init_psolver_bps(parini,atoms,poisson_cent)
        !call get_psolver_bps(poisson_cent,atoms,ehartree)
        call get_psolver(parini,poisson_cent,atoms,Q(1,:),ehartree)
        cent_ener = ehartree
        if(atoms%nat > 1)  then
            call get_hartree_force(parini,poisson_cent,atoms)
            !call force_gto_sym(parini,'bulk',atoms%nat,atoms%rat,atoms%cellvec,atoms%zat,gwit,poisson_cent%rgcut,nx,ny,nz,poisson_cent%pot,atoms%fat)
            cent_fat = atoms%fat
        end if
        force_err = sqrt(sum((cent_fat-dft_fat)**2)/(3*atoms%nat))
        ener_err = abs(cent_ener-dft_ener)
        !call fini_psolver_bps(cent%poisson)
        write(77,'(3es14.6)') t,force_err,ener_err
    enddo
    do j = 1 , atoms%nat
        write(88,*) j ,dft_fat(1:3,j), sqrt(sum(dft_fat(1:3,j)**2))
    enddo
!//E.O. Linear combination of dft and cent charge density to investigate the force_err////
    write(*,'(a,es14.6)') "Na_INT",volinv*sum(poisson_cent%rho)
    do i = 1,nx
        write(1370,'(4es14.6)') i*hx , cent_rho(i,yg_at,zg_at),dft_rho(i,yg_at,zg_at),poisson_cent%rho(i,yg_at,zg_at)
    enddo 
    write(*,*) "END OF LCN CHECK OF CENT VS DFT"
    end associate
    end associate
    end associate
    end associate
   ! call fini_hartree(parini,atoms,poisson_dft)
    !call fini_hartree(parini,atoms,poisson_cent)
end subroutine best_charge_density_rho
