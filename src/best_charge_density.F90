subroutine best_charge_density(parini)
    use mod_interface
    use mod_parini, only: typ_parini
    implicit none
    type(typ_parini), intent(in):: parini

    !call best_charge_density_rho(parini)
    call best_charge_density_pot(parini)

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
    type(typ_poisson):: poisson_dft, poisson_cent
    type(typ_atoms):: atoms
    logical :: esc
    integer:: iter,iter_max, at_range ,xg_at , yg_at , zg_at , n_at,i,j,k,lcn, ss, Ne,atoms_type! lcn = linear combinarion number
    integer, allocatable :: t_num(:)
    real(8):: pi,rho_err,ehartree,rmse,dft_ener,cent_ener,ener_err,force_err,w,at_zat , t
    real(8):: sd_s,volinv,err, a_max,qnrm,dft_q,gtot, errmax, peak, c1, d1,temp1,peakx,peaky,peakz,errmax_old
    real(8),allocatable::  dft_rho(:,:,:),cent_rho(:,:,:),cent_rho_1(:,:,:),cent_rho_a_par(:,:,:),cent_rho_q_par(:,:,:)
    real(8),allocatable:: A(:,:), Q(:,:), At(:,:), Qt(:,:), qpar(:,:), apar(:,:), qpar_t(:,:),qpar_tmp(:,:), apar_t(:,:), qtemp(:), rat_t(:,:)
    real(8),allocatable:: dft_fat(:,:),cent_fat(:,:),gwit(:)
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
    poisson_dft%hx=sqrt(poisson_dft%hgrid(1,1)**2+poisson_dft%hgrid(2,1)**2+poisson_dft%hgrid(3,1)**2)
    poisson_dft%hy=sqrt(poisson_dft%hgrid(1,2)**2+poisson_dft%hgrid(2,2)**2+poisson_dft%hgrid(3,2)**2)
    poisson_dft%hz=sqrt(poisson_dft%hgrid(1,3)**2+poisson_dft%hgrid(2,3)**2+poisson_dft%hgrid(3,3)**2)
    write(*,*) poisson_dft%hx , poisson_dft%hy , poisson_dft%hz
    associate(nx=>poisson_dft%ngpx,ny=>poisson_dft%ngpy,nz=>poisson_dft%ngpz)
    associate(cv1=>atoms%cellvec(1,1),cv2=>atoms%cellvec(2,2),cv3=>atoms%cellvec(3,3))
    associate(hx=>poisson_dft%hx,hy=>poisson_dft%hy,hz=>poisson_dft%hz)
    associate(x_at=>atoms%rat(1,n_at),y_at=>atoms%rat(2,n_at),z_at=>atoms%rat(3,n_at))
    poisson_dft%lda = nx
    poisson_cent%lda = nx
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
    call init_hartree(parini,atoms,poisson_dft,Q(1,1:atoms%nat))
    !call init_psolver_bps(parini,atoms,cent%poisson)
    call get_hartree(parini,poisson_dft,atoms,Q(1,1:atoms%nat),ehartree)
    !call get_psolver_bps(poisson_dft,atoms,ehartree)
    dft_ener = ehartree 
    if (atoms%nat > 1) then
        poisson_dft%bc='bulk'   
        poisson_dft%nat=atoms%nat
        poisson_dft%rcart=atoms%rat
        poisson_dft%cv=atoms%cellvec
        poisson_dft%q=atoms%zat
        !poisson_dft%gw_ewald=gwit
        poisson_dft%gw=gwit
        !call force_gto_sym(parini,poisson_dft%bc,poisson_dft%nat,poisson_dft%rcart, &
        !    poisson_dft%cv,poisson_dft%q,poisson_dft%gw_ewald,poisson_dft%rgcut,poisson_dft%ngpx, &
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
    poisson_cent%task_finit='alloc_rho'
    call init_hartree(parini,atoms,poisson_cent,Q(1,:))
    poisson_cent = poisson_dft
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
        call get_hartree(parini,poisson_cent,atoms,Q(1,:),ehartree)
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
        call get_hartree(parini,poisson_cent,atoms,Q(1,:),ehartree)
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
!**************************************************
!/////////////////////////////////////////////////////////////////////////////////////////
!/////////////////////////////////////////////////////////////////////////////////////////
!///////////////////////////////////// POT_PART //////////////////////////////////////////
!/////////////////////////////////////////////////////////////////////////////////////////
!/////////////////////////////////////////////////////////////////////////////////////////
subroutine best_charge_density_pot(parini)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_cent, typ_ann_arr
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    type(typ_poisson):: poisson_dft, poisson_cent
    type(typ_atoms):: atoms
    logical:: esc
    integer:: atoms_type, iter, iter_max, cel_vol ,cel_vol_inv, ix, iy, iz, l, nclx, ncrx,ncly, ncry, nclz, ncrz
    integer:: at_range ,xg_at , yg_at , zg_at , n_at, i, j, k, lcn, ss, Ne, istat, nstat, ng(3)! lcn = linear combinarion number
    integer, allocatable :: t_num(:)
    real(8):: rmse_old , total_charge, start, finish, total_time, mean, std, var
    real(8):: rclx, rcrx, rcly, rcry, rclz, rcrz, rminx, rmaxx, rminy, rmaxy, rminz, rmaxz, pi,pot_err  
    real(8):: ehartree, rmse, cent_ener, dft_ener, w, exc, volinv_integration, sd_s, volinv,err, a_max,qnrm, dft_q,gtot, peak, c1, d1, err_max, temp, at_zat
    real(8),allocatable::  dft_rho(:,:,:), weight(:,:,:), dft_pot(:,:,:), cent_pot(:,:,:)
    real(8),allocatable:: A(:,:), Q(:,:), At(:,:), Qt(:,:), qpar(:,:), apar(:,:), qpar_t(:,:), rat_t(:,:)
    real(8),allocatable:: dft_fat(:,:), cent_fat(:,:), gwit(:), at_rat(:,:),cent_rho(:,:,:),cv_temp(:,:), atom_charge(:),atom_charge_type(:)
    !////////////////////////E.O. TEMP for charge on grid/////////////////////////////////
    pi=4.d0*atan(1.d0)
    nstat = 5 
    rmse =0.d0
    rmse_old=0.d0
    open(unit=1,file='poterr')
    open(unit=2,file='control')
    open(unit=1377,file='params.txt')
!//////////////////////////////////////READING INPUT PARAMETERS///////////////////////////
    call cube_read('electronic_density.cube',atoms,poisson_dft)
    allocate(at_rat(3,atoms%nat),cv_temp(3,3))
    at_rat = atoms%rat ! because acf_read corrupts atoms%rat
    poisson_dft%hx=sqrt(poisson_dft%hgrid(1,1)**2+poisson_dft%hgrid(2,1)**2+poisson_dft%hgrid(3,1)**2)
    poisson_dft%hy=sqrt(poisson_dft%hgrid(1,2)**2+poisson_dft%hgrid(2,2)**2+poisson_dft%hgrid(3,2)**2)
    poisson_dft%hz=sqrt(poisson_dft%hgrid(1,3)**2+poisson_dft%hgrid(2,3)**2+poisson_dft%hgrid(3,3)**2)
    cv_temp=atoms%cellvec
    call acf_read(parini,'posinp.acf',1,atoms=atoms)
    atoms%rat = at_rat
    atoms%cellvec=cv_temp
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
    write(*,*) poisson_dft%hx , poisson_dft%hy , poisson_dft%hz
    associate(nx=>poisson_dft%ngpx,ny=>poisson_dft%ngpy,nz=>poisson_dft%ngpz)
    associate(cv1=>cv_temp(1,1),cv2=>cv_temp(2,2),cv3=>cv_temp(3,3))
    associate(hx=>poisson_dft%hx,hy=>poisson_dft%hy,hz=>poisson_dft%hz)
    associate(x_at=>atoms%rat(1,n_at),y_at=>atoms%rat(2,n_at),z_at=>atoms%rat(3,n_at))
    write(2,*) "Number of grids : " , nx , ny , nz
    write(2,'(a,3es14.6)') 'Grid spacing : ' , hx , hy , hz
    poisson_dft%lda = nx
    poisson_cent%lda = nx
    cel_vol = nx*ny*nz
    cel_vol_inv = 1.d0/(cel_vol+0.d0)
    allocate(dft_pot(nx,ny,nz),dft_rho(nx,ny,nz),dft_fat(1:3,1:atoms%nat),weight(nx,ny,nz),gwit(1:atoms%nat))
    ss=0
    do j = 1 , atoms_type
        ss = ss + t_num(j-1)
        rat_t(1:3,j) = atoms%rat(1:3,ss+1)
    end do
    volinv = (cv1*cv2*cv3)/(nx*ny*nz)
    xg_at = int(x_at/hx)
    yg_at = int(y_at/hy)
    zg_at = int(z_at/hz) 
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
    poisson_dft%rgcut=6.d0*a_max !ask Dr. Ghasemi if it's OK !
    write(*,*) "rgcut of dft" , poisson_dft%rgcut
    dft_fat = 0.d0
    dft_ener= 0.d0 
    atoms%fat = 0.d0 
    poisson_dft%task_finit=''
    call init_hartree(parini,atoms,poisson_dft,Q(1,1:atoms%nat))
    call get_hartree(parini,poisson_dft,atoms,Q(1,1:atoms%nat),ehartree)
    dft_ener = ehartree
    dft_pot = poisson_dft%pot
    !if (atoms%nat > 1) then
    !    poisson_dft%bc=atoms%boundcond   
    !    poisson_dft%nat=atoms%nat
    !    poisson_dft%rcart=atoms%rat
    !    poisson_dft%cv=atoms%cellvec
    !    poisson_dft%q=atoms%zat
    !    poisson_dft%gw=gwit
    !    call get_hartree_force(parini,poisson_dft,atoms)
    !    dft_fat = atoms%fat
    !end if
    dft_q = hx*hy*hz*sum(dft_rho)
    write(2,*)'Total charge density of DFT : ',dft_q
!//////////////////////////////////////////CENT PART//////////////////////////////////////
    allocate(cent_pot(nx,ny,nz))
    allocate(qpar(1:lcn,1:atoms%nat),apar(1:lcn,1:atoms%nat),qpar_t(1:lcn,1:atoms_type))
    allocate(cent_fat(1:3,1:atoms%nat))
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
    do i = 1 , lcn
        ss = 0
        do j = 1 , atoms_type
            ss = ss + t_num(j-1)
            qt(i,j)=q(i,1+ss)
            at(i,j)=a(i,1+ss)
        end do
    end do
    read(1377,*) ! read left and right boundaries for x, y and z.
    read(1377,*)    rclx , rcrx
    read(1377,*)    rcly , rcry
    read(1377,*)    rclz , rcrz
    weight = 1.d0
    exc = 0.d0
    if(rclx>=0 .and. rcrx>=0 .and. rcly>=0 .and. rcry>=0 .and. rclz>=0 .and. rcrz>=0) then
        rminx = minval(atoms%rat(1,:)); rminy = minval(atoms%rat(2,:)); rminz = minval(atoms%rat(3,:));
        rmaxx = maxval(atoms%rat(1,:)); rmaxy = maxval(atoms%rat(2,:)); rmaxz = maxval(atoms%rat(3,:));
        nclx=int((rminx-rclx)/hx); ncrx=int((rmaxx+rcrx)/hx); ncly=int((rminy-rcly)/hy); 
        ncry=int((rmaxy+rcry)/hy); nclz=int((rminz-rclz)/hz); ncrz=int((rmaxz+rcrz)/hz);
        write(2,'(a,6i4)') 'Cutoff grid numbers (ncl? ncr?) : ',nclx, ncrx, ncly, ncry, nclz, ncrz
        do iz = nclz , ncrz
            do iy = ncly , ncry
                do ix = nclx , ncrx
                    weight(ix,iy,iz) = 0.d0
                    exc = exc+1.d0
                enddo
            enddo
        enddo
    else
        write(2,*) 'Cutoff : whole cell is considered'
    end if
    volinv_integration = (cv1*cv2*cv3+0.d0)/(nx*ny*nz-exc+0.d0)
    allocate(atom_charge_type(1:atoms_type),atom_charge(1:atoms%nat))
    read(1377,*) ! atoms initial charge
    do i = 1 , atoms_type
        read(1377,*) atom_charge_type(i)
    enddo
    ss=0
    do j = 1 , atoms_type
        ss = ss + t_num(j-1)
        atom_charge(ss+1:t_num(j)+ss) = atom_charge_type(j)
    end do
    ng(1) = nx ; ng(2) = ny ; ng(3) = nz
    total_time = 0.d0
    do iter=1,huge(iter_max)
        start = 0.d0
        finish = 0.d0
        call cpu_time(start)
        inquire(file='EXIT',exist=esc)
        if(esc .and. iter >2) exit 
        a_max = maxval(at)
        poisson_cent%rgcut=6.d0*a_max !ask Dr. Ghasemi if it's OK !
!-------------------------------------calculating potential and its derivatives--------------------
        call put_pot_sym_rzx(atoms%rat,hx,hy,hz,atoms%nat,Q(1:lcn,1:atoms%nat),A(1:lcn,1:atoms%nat),ng(1:3),lcn,.true. &
            ,weight(1:nx,1:ny,1:nz),poisson_dft%pot(1:nx,1:ny,1:nz),cent_pot(1:nx,1:ny,1:nz),qpar(1:lcn,1:atoms%nat),apar(1:lcn,1:atoms%nat))
        gtot = sum(qpar)
        qpar = qpar - (gtot+0.d0)/(atoms%nat*lcn+0.d0)
        do i = 1 , atoms%nat
            do l = 1 , lcn
                write(2,*) 'Atom, ITER, LCN, QPAR : ',atoms%sat(i),iter,l,qpar(l,i)
            end do
            do l = 1 , lcn
                write(2,*) 'Atom, ITER, LCN, APAR : ',atoms%sat(i),iter,l,apar(l,i)
            end do
        end do
        Q = Q - SD_S*qpar
        A = A - SD_S*apar
        do i = 1 , atoms%nat
            total_charge = 0.d0
            do l = 1 , lcn
                total_charge = total_charge+Q(l,i)
                write(2,*) 'Atom, ITER, LCN, Q : ',atoms%sat(i),iter,l,Q(l,i)
            end do
            write(2,'(a,a4,2es14.6)') 'Atom, Total charge, charge changes : ',atoms%sat(i), total_charge, -abs(atom_charge(i))+abs(total_charge)
            do l = 1 , lcn
                write(2,*) 'Atom, ITER, LCN, A : ',atoms%sat(i),iter,l,A(l,i)
            end do
            write(2,*) '---------------------------'
        end do
        write(2,*) ' Total charge of cell: ',sum(Q)
        do i = 1 , lcn
            ss = 0
            do j = 1 , atoms_type
                ss = ss + t_num(j-1)
                qt(i,j)=q(i,1+ss)
                at(i,j)=a(i,1+ss)
            end do
        end do
        err_max = 0.d0
        pot_err = 0.d0
        do i = 1 , nz
            do j = 1 , ny
                do k = 1 , nx
                    c1 = cent_pot(k,j,i)
                    d1 = dft_pot(k,j,i)
                    temp = weight(k,j,i)*abs(c1-d1)/abs(d1+0.d0)
                    if (err_max < temp ) err_max = temp
                    pot_err = pot_err + (weight(k,j,i)*(c1-d1)**2)
                end do
            end do
        end do
        pot_err= hx*hy*hz*pot_err
        rmse_old=rmse
        rmse = sqrt(pot_err)
        write(2,'(a,i6,3es14.6)')'ITER RMSE POT_ERR ERR_MAX: ',iter,rmse,pot_err,err_max
        if (abs(rmse_old-rmse) < err )then
            istat=istat+1
        else
            istat=0
        endif
        call cpu_time(finish)
        total_time = total_time + finish - start
        write(2,'(a45,f6.3,a45)') "================================ ITER TIME : ",finish-start,"(sec)========================================"
        if(istat>nstat) then
            write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!MAX CONVERSION REACHED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            write(2,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!MAX CONVERSION REACHED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            exit
        endif
    end Do
    write(2,'(a36,f10.3,f6.3)') 'Total SD time, average per iter : ' , total_time , (total_time+0.d0)/(iter+0.d0)
    ss = 0
    do j = 1 , atoms_type
        ss = ss+t_num(j-1)
        do i = 1 , lcn    
            call stdval_rzx(Q(i,ss+1:ss+t_num(j)),t_num(j),mean,std,var)
            write(2,'(a,a3,i,3es14.6)') 'LCN, MEAN, VAR, STD Q of :',atoms%sat(ss+1), i,mean,var,std
            call stdval_rzx(A(i,ss+1:ss+t_num(j)),t_num(j),mean,std,var)
            write(2,'(a,a3,i,3es14.6)') 'LCN, MEAN, VAR, STD A of :',atoms%sat(ss+1), i,mean,var,std
        enddo
    enddo
    do ix = 1 , nx
        write(3,*) ix*hx , dft_pot(ix,ny/2,nz/2) ,cent_pot(ix,ny/2,nz/2)
    end do 
!/////////////////////////////////E.O. CENT PART//////////////////////////////////////////
!/////////////////////////// CENT CHARGE DENSITY PART/////////////////////////////////////
!    allocate(poisson_cent%q(1:atoms%nat), poisson_cent%gw_ewald(1:atoms%nat), poisson_cent%rcart(1:3,1:atoms%nat))
!    poisson_cent%nat=atoms%nat
!    poisson_cent%lda = nx
!    poisson_cent%hgrid=poisson_dft%hgrid
!    poisson_cent%ngpx = nx
!    poisson_cent%ngpy = ny
!    poisson_cent%ngpz = nz
!    poisson_cent%alpha= parini%alpha_ewald
!    poisson_cent%task_finit='set_ngp:alloc_rho'
!    call init_hartree_bps(parini,atoms,poisson_cent)
!    poisson_cent%rcart=atoms%rat
!    poisson_cent%cv=cv_temp
!    a_max = maxval(at)
!    poisson_cent%bc='free'!atoms%boundcond
!    do i = 1 , lcn
!        poisson_cent%rgcut=6.d0*maxval(A(i,:))
!        poisson_cent%q=Q(i,1:atoms%nat)
!        poisson_cent%gw_ewald(1:atoms%nat)=A(i,1:atoms%nat)
!        if (i==1) then         
!            poisson_cent%reset_rho=.true.
!            call put_charge_density(parini,poisson_cent)
!        else
!            poisson_cent%reset_rho=.false.
!            call put_charge_density(parini,poisson_cent)
!        end if
!    end do
!    do ix = 1 , nx
!        write(4,*) ix*hx , dft_rho(ix,ny/2,nz/2) ,poisson_cent%rho(ix,ny/2,nz/2)
!    end do 
!    atoms%fat = 0.d0
!    call get_hartree(parini,poisson_cent,atoms,Q(1,:),ehartree)
!    cent_ener = ehartree
!    poisson_cent%q=atoms%zat
!    poisson_cent%gw_ewald=gwit
!    call get_hartree_force(parini,poisson_cent,atoms)
!    write(*,*) 'cent_ener : ' , cent_ener
!    cent_fat = atoms%fat
!/////////////////////////////////////////////////////////////////////////////////////////
!    poisson_cent%nat=atoms%nat
!    poisson_cent%ngpx = nx
!    poisson_cent%ngpy = ny
!    poisson_cent%ngpz = nz
!    poisson_cent%hgrid=poisson_dft%hgrid
!    poisson_cent%hx=poisson_dft%hx
!    poisson_cent%hy=poisson_dft%hy
!    poisson_cent%hz=poisson_dft%hz
!!????????????????????????????????????????????????????????
!    poisson_cent%alpha=parini%alpha_ewald
!    write(*,*) "C_1" , poisson_cent%ngpx,poisson_cent%ngpy,poisson_cent%ngpz
!    poisson_cent%task_finit='alloc_rho'
!    call init_hartree_bps(parini,atoms,poisson_cent)
!    write(*,*) "C_2" , poisson_cent%ngpx,poisson_cent%ngpy,poisson_cent%ngpz
!    write(*,*) "C_3" , shape(poisson_cent%rho)
!    poisson_cent%lda = nx
!    poisson_cent%cv=atoms%cellvec
!    poisson_cent%bc=atoms%boundcond
!    allocate(poisson_cent%q(1:poisson_cent%nat), poisson_cent%gw_ewald(1:poisson_cent%nat), poisson_cent%rcart(1:3,1:poisson_cent%nat))
!    poisson_cent%rcart(1:3,1:poisson_cent%nat)=atoms%rat(1:3,1:atoms%nat)
!    poisson_cent%q(1:poisson_cent%nat)=Q(1,1:atoms%nat)
!    poisson_cent%gw_ewald(1:poisson_cent%nat)=A(1,1:atoms%nat)
!    poisson_cent%reset_rho=.true.
!    call put_charge_density(parini,poisson_cent)
!    poisson_cent%reset_rho=.false.
!    do i = 2 , lcn
!        poisson_cent%q(1:poisson_cent%nat)=Q(i,1:atoms%nat)
!        poisson_cent%gw_ewald(1:poisson_cent%nat)=A(i,1:atoms%nat)
!        call put_charge_density(parini,poisson_cent)
!    enddo
!    do ix = 1 , nx
!        write(4,*) ix*hx , dft_rho(ix,ny/2,nz/2) ,poisson_cent%rho(ix,ny/2,nz/2)
!    end do  
!////////////////////////E.O. CENT CHARGE DENSITY PART////////////////////////////////////
!////////////////////////charge on grid///////////////////////////////////////////////////
    allocate(cent_rho(1:nx,1:ny,1:nz))
    call put_gto_sym_ortho_rzx(atoms%rat,poisson_dft%hgrid,atoms%nat,Q,A,ng,lcn,.true.,cent_rho)
    do ix = 1 , nx
        write(4,*) ix*hx , dft_rho(ix,ny/2,nz/2) ,cent_rho(ix,ny/2,nz/2)
    end do
!////////////////////////ENER PART///////////////////////////////////////
    poisson_cent%nat=atoms%nat
    poisson_cent%lda = nx
    poisson_cent%hgrid=poisson_dft%hgrid
    poisson_cent%ngpx = nx; poisson_cent%ngpy = ny; poisson_cent%ngpz = nz;
    poisson_cent%hx=poisson_dft%hx; poisson_cent%hy=poisson_dft%hy; poisson_cent%hz=poisson_dft%hz;
    allocate(poisson_cent%rho(1:nx,1:ny,1:nz),poisson_cent%gw_ewald(1:poisson_cent%nat), poisson_cent%rcart(1:3,1:poisson_cent%nat))
    poisson_cent%rho(1:nx,1:ny,1:nz)=cent_rho(1:nx,1:ny,1:nz)
    poisson_cent%gw_ewald(1:poisson_cent%nat)=A(1,1:atoms%nat) !!! ASK DR GHASEMI
    poisson_cent%alpha= parini%alpha_ewald
    poisson_cent%task_finit=''
    call init_hartree_bps(parini,atoms,poisson_cent)
    poisson_cent%rcart=atoms%rat
    poisson_cent%cv=cv_temp
    poisson_cent%bc='free'!atoms%boundcond
    poisson_cent%rgcut=6.d0*maxval(A(:,:))
    poisson_cent%rho=-1.d0*poisson_cent%rho
    call cube_write('cent_rho.cube',atoms,poisson_cent,'rho')
    poisson_cent%rho=-1.d0*poisson_cent%rho
    call get_hartree(parini,poisson_cent,atoms,Q(1,:),ehartree)
    cent_ener = ehartree
    write(2,'(a,4es14.6)') 'CENT energy , DFT energy , Energy difference , Err Percentage : ' , cent_ener , dft_ener , cent_ener-dft_ener,(cent_ener-dft_ener)*100.d0/dft_ener
!//////////////////////E.O. ENER PART////////////////////////////////////
!////////////////////////E.O. charge on grid//////////////////////////////////////////////
    end associate
    end associate
    end associate
    end associate
    call fini_hartree(parini,atoms,poisson_dft)
end subroutine best_charge_density_pot
!/////////////////////////////////////////////////////////////////////////////////////////
subroutine put_gto_sym_ortho_rzx(rat,hgrid,nat,qat,gw,ng,lcn,reset,cent_rho)
    implicit none
    logical :: reset
    integer :: nat, lcn,ng(3)
    real(8) :: rat(3,nat), qat(lcn,nat), gw(lcn,nat),cent_rho(1:ng(1),1:ng(2),1:ng(3)),hgrid(1:3,1:3)
    ! local variables
    integer :: i, iat, iw, iz, iy, ix
    real(8) :: pi, xat, yat, zat, width, width_inv, fac,facqiat, rhoz, rhoyz, hx, hy, hz 
    real(8) :: width_inv_hgx, width_inv_hgy, width_inv_hgz, width_inv_xat, width_inv_yat, width_inv_zat
    real(8) , allocatable:: wx(:), wy(:), wz(:)
    pi = 4.d0*atan(1.d0)
    allocate(wx(1:ng(1)),wy(1:ng(2)),wz(1:ng(3)))
    hx=sqrt(hgrid(1,1)**2+hgrid(2,1)**2+hgrid(3,1)**2)
    hy=sqrt(hgrid(1,2)**2+hgrid(2,2)**2+hgrid(3,2)**2)
    hz=sqrt(hgrid(1,3)**2+hgrid(2,3)**2+hgrid(3,3)**2)
    if(reset) cent_rho=0.d0
        do i = 1 ,lcn
            do iat=1,nat
                xat=rat(1,iat)+hx
                yat=rat(2,iat)+hy
                zat=rat(3,iat)+hz
                width=gw(i,iat)
                width_inv=1.d0/width
                fac=1.d0/(width*sqrt(pi))**3
                width_inv_hgx=width_inv*hx
                width_inv_hgy=width_inv*hy
                width_inv_hgz=width_inv*hz
                width_inv_xat=width_inv*xat
                width_inv_yat=width_inv*yat
                width_inv_zat=width_inv*zat
                do iw=1,ng(1)
                    wx(iw)=exp(-(width_inv_hgx*iw-width_inv_xat)**2)
                enddo
                do iw=1,ng(2)
                    wy(iw)=exp(-(width_inv_hgy*iw-width_inv_yat)**2)
                enddo
                do iw=i,ng(3)
                    wz(iw)=exp(-(width_inv_hgz*iw-width_inv_zat)**2)
                enddo
                facqiat=fac*qat(i,iat)
                do iz=1,ng(3)
                    rhoz=facqiat*wz(iz)
                    do iy=1,ng(2)
                        rhoyz=rhoz*wy(iy)
                        do ix=1,ng(1)
                            cent_rho(ix,iy,iz)=cent_rho(ix,iy,iz)+rhoyz*wx(ix)
                        enddo
                    enddo
                enddo
            enddo
        enddo
end subroutine put_gto_sym_ortho_rzx
!/////////////////////////////////////////////////////////////////////////////////////////
subroutine put_pot_sym_rzx(rat,hx,hy,hz,nat,qat,gw,ng,lcn,reset,weight,dft_pot,cent_pot,qpar,apar)
    implicit none
    logical :: reset
    integer , intent(in):: nat, ng(1:3), lcn
    real(8) , intent(in):: rat(1:3,1:nat), hx, hy, hz, qat(1:lcn,1:nat),gw(1:lcn,1:nat),weight(1:ng(1),1:ng(2),1:ng(3))
    real(8) , intent(in):: dft_pot(1:ng(1),1:ng(2),1:ng(3))
    real(8) , intent(out):: cent_pot(1:ng(1),1:ng(2),1:ng(3)), apar(1:lcn,1:nat), qpar(1:lcn,1:nat)
    ! local variables 
    integer :: l, iat, iz, iy, ix, nx, ny, nz
    real(8) :: asqinv, pisqrtinv, pi, r, dx, dy, dz
    real(8) :: cent_pot_a_par(1:ng(1),1:ng(2),1:ng(3)), cent_pot_q_par(1:ng(1),1:ng(2),1:ng(3))
    nx=ng(1); ny=ng(2); nz=ng(3)
    pi = 4.d0*atan(1.d0)
    if (reset) cent_pot = 0.d0
        do l = 1,lcn
            do iat=1,nat
                asqinv=1.d0/(gw(l,iat)**2)
                pisqrtinv = 1.d0/sqrt(pi)
                do iz=1,nz
                    dz=(iz-1)*hz - rat(3,iat)
                    do iy=1,ny
                        dy=(iy-1)*hy - rat(2,iat)
                        do ix=1,nx
                            dx=(ix-1)*hx - rat(1,iat)
                            r=sqrt(dx**2+dy**2+dz**2)
                            if (r <= 1.d-6) then 
                                cent_pot(ix,iy,iz)=cent_pot(ix,iy,iz)+2.d0*qat(l,iat)*pisqrtinv/gw(l,iat)
                            else
                                cent_pot(ix,iy,iz)=cent_pot(ix,iy,iz)+qat(l,iat)*erf(r/gw(l,iat))/r
                            endif
                        enddo ! of ix
                    enddo ! of iy
                enddo ! of izi
            enddo ! of iat
        enddo !of lcn
        do l = 1,lcn
            do iat=1,nat
                asqinv=1.d0/(gw(l,iat)**2)
                pisqrtinv = 1.d0/sqrt(pi)
                do iz=1,nz
                    dz=(iz-1)*hz - rat(3,iat)
                    do iy=1,ny
                        dy=(iy-1)*hy - rat(2,iat)
                        do ix=1,nx
                            dx=(ix-1)*hx - rat(1,iat)
                            r=sqrt(dx**2+dy**2+dz**2)
                            if (r <= 1.d-6) then 
                                cent_pot_a_par(ix,iy,iz) =-2.d0*qat(l,iat)*asqinv*pisqrtinv
                                cent_pot_q_par(ix,iy,iz) =2.d0*pisqrtinv/gw(l,iat)
                            else
                                cent_pot_a_par(ix,iy,iz)=-2.d0*qat(l,iat)*exp(-r**2*asqinv)*asqinv*pisqrtinv
                                cent_pot_q_par(ix,iy,iz) =erf(r/gw(l,iat))/r
                            endif
                        enddo ! of ix
                    enddo ! of iy
                enddo ! of izi
                qpar(l,iat)=2.d0*hx*hy*hz*sum(cent_pot_q_par*weight*(cent_pot-dft_pot))
                apar(l,iat)=2.d0*hx*hy*hz*sum(cent_pot_a_par*weight*(cent_pot-dft_pot))
            enddo ! of iat
        enddo !of lcn
end subroutine put_pot_sym_rzx
!/////////////////////////////////////////////////////////////////////////////////////////
subroutine stdval_rzx(f,f_len,mean,std,var)
    implicit none
    integer, intent(in) :: f_len
    real(8), intent(in) :: f(f_len)
    real(8), intent(out) :: mean, std, var
    ! local variables
    integer :: i
    real(8) :: g(f_len)
    mean = sum(f(1:f_len))/(f_len+0.d0)
    do i = 1, f_len
        g(i) = (f(i)-mean)**2
    end do
    var = sum(g(1:f_len))/(f_len+0.d0)
    std = sqrt(var)
end subroutine stdval_rzx
