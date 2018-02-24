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
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
subroutine best_charge_density_rho(parini)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson, typ_poisson
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_cent, typ_ann_arr
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    type(typ_poisson):: poisson ,poisson_rough
    type(typ_atoms):: atoms
    type(typ_cent):: cent
    type(typ_ann_arr):: ann_arr
    integer:: typ_at,istat, igpx, igpy, igpz, iat, nx ,ny,nz,iter,iter_max,gweiat,i,j,k,jj,lcn,ierr,lwork, ss, Ne! lcn = linear combinarion number
    real(8):: cell(3), epot, rgcut_a, t1, t2, t3, t4, pi,h,rho_err,hx,hy,hz
    real(8):: ehartree,rmse,SDA,dft_ener,cent_ener,ener_err,force_err,ener_err_old,coeff,w
    real(8),allocatable:: dft_rho(:,:,:), cent_rho(:,:,:),cent_rho_1(:,:,:),cent_rho_a_par(:,:,:),cent_rho_q_par(:,:,:),work(:),EE(:,:),eig(:)
    real(8),allocatable:: A(:,:), Q(:,:), At(:,:), Qt(:,:), qpar(:,:), apar(:,:), qpar_t(:,:),qpar_tmp(:,:), apar_t(:,:), qtemp(:), rat_t(:,:)
    real(8),allocatable:: dft_fat(:,:),cent_fat(:,:)
    real(8):: sd_s,err_fdm,err_cent,volinv,err
    real(8):: q_max,a_max,qnrm,dft_q,gtot
    real(8):: errmax, peak, c1, d1,temp1,cv1,cv2,cv3,peakx,peaky,peakz,errmax_old
    integer, allocatable :: t_num(:)
    logical :: esc
    real(8) :: x_at , y_at , z_at , at_zat , t
    integer :: at_range ,xg_at , yg_at , zg_at , n_at
    open(unit=1370,file='plot')
    
    open(unit=1,file='rhoerr')
    open(unit=2,file='errmax')
    open(unit=3,file='enererr')
    open(unit=4,file='forceerr')
    open(unit=5,file='peakloc')

    call cube_read('electronic_density.cube',atoms,cent%poisson)
    call acf_read(parini,'posinp.acf',1,atoms=atoms)
    open(unit=1377,file='params.txt')
    pi=4.d0*atan(1.d0)
    read(1377,*) !Number of STO-NG's N
    read(1377,*) lcn
    read(1377,*) !Number of atoms types
    read(1377,*) typ_at
    allocate(t_num(0:typ_at),Qt(lcn,typ_at),At(lcn,typ_at),rat_t(1:3,typ_at),Q(1:lcn,1:atoms%nat),A(1:lcn,1:atoms%nat))
    t_num(0) = 0
    read(1377,*) !Number of each atom type Number
    do i = 1 , typ_at
        read(1377,*) t_num(i)
    end do
    read(1377,*) !Q's and alpha's of atom types
    do j = 1 , typ_at
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
    read(1377,*) !to be ploted atom number
    read(1377,*) n_at
    ss=0
    do j = 1 , typ_at
        ss = ss + t_num(j-1)
        rat_t(1:3,j) = atoms%rat(1:3,ss+1)
    end do
    nx=cent%poisson%ngpx
    ny=cent%poisson%ngpy
    nz=cent%poisson%ngpz
    cv1=atoms%cellvec(1,1)
    cv2=atoms%cellvec(2,2)
    cv3=atoms%cellvec(3,3)
    volinv = (cv1*cv2*cv3)/(nx*ny*nz)
    hx=cent%poisson%hx
    hy=cent%poisson%hy
    hz=cent%poisson%hz
    x_at = atoms%rat(1,n_at) 
    y_at = atoms%rat(2,n_at)
    z_at = atoms%rat(3,n_at)
    xg_at = int(x_at/hx)
    yg_at = int(y_at/hy)
    zg_at = int(z_at/hz) 
    write(*,'(a,3i4)') "atom_grid",xg_at,yg_at,zg_at
    allocate(cent%poisson%pot(nx,ny,nz))
    cent%poisson%rho=-1.d0*cent%poisson%rho
    allocate(cent%gwit(1:atoms%nat))
    allocate(dft_rho(nx,ny,nz),dft_fat(1:3,1:atoms%nat))
    dft_rho = cent%poisson%rho
    peak = maxval(abs(dft_rho))
    w = 1.d0/(Ne+0.d0)
    cent%gwit(:) = 0.5d0
    read(1377,*)!atoms z_at
    at_range = atoms%nat/typ_at
    do i = 1 , typ_at
        read(1377,*) at_zat
        atoms%zat(((i-1)*at_range)+1:at_range*i) = at_zat ! this will not be OK
                                                          ! If a molecule
                                                          ! do not have same
                                                          ! number of each atom 
                                                          ! for example :
                                                          ! SrTiO_3
    end do
    write(*,*) "ZAT" , atoms%zat
    a_max = maxval(at)
    cent%poisson%rgcut=8.d0*a_max !ask Dr. Ghasemi if it's OK !
    atoms%boundcond='bulk'   
    dft_fat = 0.d0
    dft_ener= 0.d0 
    cent%poisson%pot = 0.d0
    atoms%fat = 0.d0 
    call construct_ewald_bps(parini,atoms,cent%poisson)
    call psolver_bps(cent%poisson,atoms,ehartree)
    dft_ener = ehartree 
    if (atoms%nat > 1) then
        call gauss_force(parini,'bulk',atoms%nat,atoms%rat,atoms%cellvec,atoms%zat,cent%gwit,cent%poisson%rgcut,nx,ny,nz,cent%poisson%pot,atoms%fat)
        dft_fat = atoms%fat
    end if
    write(21,*) atoms%nat 
    write(21,*) atoms%rat
    write(21,*) atoms%cellvec
    write(21,*) atoms%zat
    write(21,*) cent%gwit
    write(21,*) cent%poisson%rgcut
    write(21,*) nx,ny,nz
    write(21,*) atoms%fat
    write(21,*) ehartree
    write(21,*) cent%poisson%pot(:,ny/2,nz/2) 
    write(21,*) "------------------------------------------------------------"
    call destruct_ewald_bps(cent%poisson)
    write(*,'(a,es14.6,3es14.6)') "DFT1" , dft_ener , dft_fat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!cent part!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
    allocate(cent_rho(nx,ny,nz))
    allocate(cent_rho_1(nx,ny,nz))
    allocate(cent_rho_a_par(nx,ny,nz))
    allocate(cent_rho_q_par(nx,ny,nz))
    allocate(qpar(1:lcn,1:atoms%nat),apar(1:lcn,1:atoms%nat),qpar_t(1:lcn,1:typ_at),qpar_tmp(1:lcn,1:typ_at),apar_t(1:lcn,1:typ_at))
    allocate(cent_fat(1:3,1:atoms%nat))
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
  !eig_part!      call gauss_grid(parini,atoms%boundcond,.true.,1,atoms%rat(1:3,i),atoms%cellvec,qtemp(i),A(1,i),cent%poisson%rgcut,nx,ny,nz,cent%poisson%rho)
  !eig_part!      cent_rho_1 = cent%poisson%rho
  !eig_part!      do j = i,atoms%nat
  !eig_part!          call gauss_grid(parini,atoms%boundcond,.true.,1,atoms%rat(1:3,j),atoms%cellvec,qtemp(j),A(1,j),cent%poisson%rgcut,nx,ny,nz,cent%poisson%rho)
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
    allocate(qtemp(1:atoms%nat))
    qtemp = 1.d0

    do i = 1 , lcn
        ss=0
        do j = 1 , typ_at
            ss = ss + t_num(j-1)
            Q(i,ss+1:t_num(j)+ss) = Qt(i,j)
            A(i,ss+1:t_num(j)+ss) = At(i,j)
        end do
    end do
    
    dft_q = volinv*sum(dft_rho)
    qnrm = (dft_q+0.d0)/(sum(Q(1:lcn,1:atoms%nat)))
    Q(1:lcn,1:atoms%nat) = Q(1:lcn,1:atoms%nat)*qnrm
    
    do i = 1 , lcn
        ss = 0
        do j = 1 , typ_at
            ss = ss + t_num(j-1)
            qt(i,j)=q(i,1+ss)
            at(i,j)=a(i,1+ss)
        end do
    end do
    
    do iter=1,huge(iter_max)
        a_max = maxval(at)
        cent%poisson%rgcut=8.d0*a_max !ask Dr. Ghasemi if it's OK !
        write(*,*) "rgcut" , cent%poisson%rgcut
        cent%poisson%rho = 0.d0
        cent%poisson%pot = 0.d0
        do i = 1 , lcn
            if (i==1) then         
                call rp4gauss_grid_rzx(parini,atoms%boundcond,.true.,atoms%nat,atoms%rat,atoms%cellvec,Q(i,1:atoms%nat),A(i,1:atoms%nat),cent%poisson%rgcut,nx,ny,nz,cent%poisson%rho,cent_rho_q_par,cent_rho_a_par)
            else
                call rp4gauss_grid_rzx(parini,atoms%boundcond,.false.,atoms%nat,atoms%rat,atoms%cellvec,Q(i,1:atoms%nat),A(i,1:atoms%nat),cent%poisson%rgcut,nx,ny,nz,cent%poisson%rho,cent_rho_q_par,cent_rho_a_par)
            end if
        end do
        cent_rho = cent%poisson%rho
        inquire(file='EXIT',exist=esc)
        if(esc) exit !to be moved after calling gauess grid
        write(*,*) "cent_rho" , volinv*sum(cent_rho) , volinv*sum(dft_rho)
        atoms%fat = 0.d0
        cent_fat = 0.d0
        call construct_ewald_bps(parini,atoms,cent%poisson)
        call psolver_bps(cent%poisson,atoms,ehartree)
        cent_ener = ehartree
        if(atoms%nat > 1)  then
            call gauss_force(parini,'bulk',atoms%nat,atoms%rat,atoms%cellvec,atoms%zat,cent%gwit,cent%poisson%rgcut,nx,ny,nz,cent%poisson%pot,atoms%fat)
            cent_fat = atoms%fat
        end if
        call destruct_ewald_bps(cent%poisson)
        write(66,'(i4, a)') iter ,"  -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
        do i = 1 , atoms%nat
            do j =1 , 3 
                write(66,'(i4,3es14.6)') i,dft_fat(j,i),abs(dft_fat(j,i)-cent_fat(j,i)),cent_fat(j,i)-dft_fat(j,i)
            end do
        end do
        do i = 1 , typ_at
            do j = 1 , lcn
                call rp4gauss_grid_rzx(parini,atoms%boundcond,.true.,1,rat_t(1:3,i),atoms%cellvec,Qt(j,i),At(j,i),cent%poisson%rgcut,nx,ny,nz,cent_rho_1, cent_rho_q_par, cent_rho_a_par)
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
            do j = 1 , typ_at
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
        do j = 1 , typ_at
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
            do j = 1 , typ_at
                ss = ss + t_num(j-1)
                qt(i,j)=q(i,1+ss)
                at(i,j)=a(i,1+ss)
            end do
        end do
    end Do
    !Specified Atoms charge density part!
    do i = 1 , lcn
        if (i==1) then         
            call rp4gauss_grid_rzx(parini,atoms%boundcond,.true.,1,atoms%rat(1:3,n_at),atoms%cellvec,Q(i,n_at),A(i,n_at),cent%poisson%rgcut,nx,ny,nz,cent%poisson%rho,cent_rho_q_par,cent_rho_a_par)
        else
            call rp4gauss_grid_rzx(parini,atoms%boundcond,.false.,1,atoms%rat(1:3,n_at),atoms%cellvec,Q(i,n_at),A(i,n_at),cent%poisson%rgcut,nx,ny,nz,cent%poisson%rho,cent_rho_q_par,cent_rho_a_par)
            
        end if
    end do
    !End of Specified Atom's charge density part!
    !Linear combination of dft and cent charge density to investigate the force_err!
    cent%poisson%rho = 0.d0
    cent%poisson%pot = 0.d0
    do i = 1 , lcn
        if (i==1) then         
            call rp4gauss_grid_rzx(parini,atoms%boundcond,.true.,atoms%nat,atoms%rat,atoms%cellvec,Q(i,1:atoms%nat),A(i,1:atoms%nat),cent%poisson%rgcut,nx,ny,nz,cent%poisson%rho,cent_rho_q_par,cent_rho_a_par)
        else
            call rp4gauss_grid_rzx(parini,atoms%boundcond,.false.,atoms%nat,atoms%rat,atoms%cellvec,Q(i,1:atoms%nat),A(i,1:atoms%nat),cent%poisson%rgcut,nx,ny,nz,cent%poisson%rho,cent_rho_q_par,cent_rho_a_par)
        end if
    end do
    cent_rho = cent%poisson%rho
    do i = 1 , 9
        cent%poisson%rho = 0.d0
        cent%poisson%pot = 0.d0
        t = i*0.1 
        atoms%fat = 0.d0
        cent_fat = 0.d0
        cent%poisson%rho=cent_rho-(t*(cent_rho-dft_rho))
        call construct_ewald_bps(parini,atoms,cent%poisson)
        call psolver_bps(cent%poisson,atoms,ehartree)
        cent_ener = ehartree
        if(atoms%nat > 1)  then
            call gauss_force(parini,'bulk',atoms%nat,atoms%rat,atoms%cellvec,atoms%zat,cent%gwit,cent%poisson%rgcut,nx,ny,nz,cent%poisson%pot,atoms%fat)
            cent_fat = atoms%fat
        end if
        force_err = sqrt(sum((cent_fat-dft_fat)**2)/(3*atoms%nat))
        ener_err = abs(cent_ener-dft_ener)
        call destruct_ewald_bps(cent%poisson)
        write(77,'(3es14.6)') t,force_err,ener_err
    enddo
    do j = 1 , atoms%nat
        write(88,*) j ,dft_fat(1:3,j), sqrt(sum(dft_fat(1:3,j)**2))
    enddo
    !End Of Linear combination of dft and cent charge density to investigate the force_err!
    write(*,'(a,es14.6)') "Na_INT",volinv*sum(cent%poisson%rho)
    do i = 1,nx
        write(1370,'(4es14.6)') i*hx , cent_rho(i,yg_at,zg_at),dft_rho(i,yg_at,zg_at),cent%poisson%rho(i,yg_at,zg_at)
    enddo 
  
end subroutine best_charge_density_rho
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!BEST_CHARGE_DENSITY_FORCE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!--------------------------------------------------------------------------------------------------
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
subroutine best_charge_density_force(parini)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_cent, typ_ann_arr
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    type(typ_atoms):: atoms
    type(typ_cent):: cent
    type(typ_ann_arr):: ann_arr
    integer:: istat, igpx, igpy, igpz, iat, nx ,ny,nz,iter,iter_max,gweiat,i,j,jj,lcn ! lcn = linear combinarion number
    real(8):: cell(3), epot, rgcut_a, t1, t2, t3, t4, pi,h,rho_err
    real(8):: ehartree,exit_cond,SDA,ener,dft_ener,coeff
    real(8),allocatable::gwit(:), dft_fat_1(:,:), dft_fat_2(:,:), cent_fat_1(:,:), cent_fat_2(:,:), agrad(:)
    real(8),allocatable:: dft_rho_1(:,:,:), dft_rho_2(:,:,:), cent_rho_1(:,:,:), cent_rho_2(:,:,:)
    real(8),allocatable::err_qgrad(:,:),err_agrad(:,:), A(:,:) , Q(:,:)
    real(8):: sd_s,fd_s,err_fdm,err_cent
    real(8):: q1,q2,q3,a11,a12,a21,a22,a31,a32
    open(unit=1370,file='param.txt')
    open(unit=1377,file='param_2.txt')
    pi=4.d0*atan(1.d0)
    lcn = 1 
   !///////////////////////////1st DFT////////////////////////////////////

    call cube_read('electronic_density_1.cube',atoms,cent%poisson)
    allocate(dft_rho_1(cent%poisson%ngpx,cent%poisson%ngpy,cent%poisson%ngpz))
    cent%poisson%rho=-1.d0*cent%poisson%rho
    dft_rho_1 = -1.d0*cent%poisson%rho
    call acf_read(parini,'posinp_1.acf',1,atoms=atoms)

    atoms%boundcond='bulk'
    allocate(cent%poisson%pot(cent%poisson%ngpx,cent%poisson%ngpy,cent%poisson%ngpz),stat=istat)
    call construct_ewald_bps(parini,atoms,cent%poisson)
    call psolver_bps(cent%poisson,atoms,ehartree)
    atoms%zat(:)=1.d0
    nx=cent%poisson%ngpx
    ny=cent%poisson%ngpy
    nz=cent%poisson%ngpz
    allocate(cent%gwit(atoms%nat))
    allocate(dft_fat_1(lcn,atoms%nat))
    cent%gwit = 1.3d0
    atoms%fat = 0.d0
    dft_fat_1 = 0.d0
    call gauss_force(parini,'bulk',atoms%nat,atoms%rat,atoms%cellvec,atoms%zat,cent%gwit,cent%poisson%rgcut,nx,ny,nz,cent%poisson%pot,atoms%fat)
    dft_fat_1 = atoms%fat
    call destruct_ewald_bps(cent%poisson)
    deallocate(cent%poisson%pot)
   !///////////////////////////2nd DFT////////////////////////////////////
    cent%poisson%rho = 0.d0
    call cube_read('electronic_density_2.cube',atoms,cent%poisson)
    allocate(dft_rho_2(cent%poisson%ngpx,cent%poisson%ngpy,cent%poisson%ngpz))
    cent%poisson%rho=-1.d0*cent%poisson%rho
    dft_rho_2 = -1.d0*cent%poisson%rho
    call acf_read(parini,'posinp_2.acf',1,atoms=atoms)
    
    atoms%boundcond='bulk'
    allocate(cent%poisson%pot(cent%poisson%ngpx,cent%poisson%ngpy,cent%poisson%ngpz),stat=istat)
    call construct_ewald_bps(parini,atoms,cent%poisson)
    call psolver_bps(cent%poisson,atoms,ehartree)
    atoms%zat(:)=1.d0
    nx=cent%poisson%ngpx
    ny=cent%poisson%ngpy
    nz=cent%poisson%ngpz
    allocate(dft_fat_2(lcn,atoms%nat),cent_fat_1(lcn,atoms%nat),cent_fat_2(lcn,atoms%nat))
    cent%gwit = 1.3d0
    atoms%fat = 0.d0
    dft_fat_2 = 0.d0
    call gauss_force(parini,'bulk',atoms%nat,atoms%rat,atoms%cellvec,atoms%zat,cent%gwit,cent%poisson%rgcut,nx,ny,nz,cent%poisson%pot,atoms%fat)
    dft_fat_2 = atoms%fat
    call destruct_ewald_bps(cent%poisson)
    deallocate(cent%poisson%pot)
 
    !/////////////////CENT_FORCE////////////////////////! 
  !  read(1370,*)
  !  read(1370,*) q1,q2,q3!,a11,a12,a21,a22,a31,a32,sd_s,fd_s
  !  read(1370,*)
  !  read(1370,*) a11,a12
  !  read(1370,*)
  !  read(1370,*) a21,a22
  !  read(1370,*)
  !  read(1370,*) a31,a32
  !  read(1370,*)
  !  read(1370,*) sd_s
  !  read(1370,*)
  !  read(1370,*) fd_s
  !  write(*,*) "params : ",q1,q2,q3,a11,a12,a21,a22,a31,a32,sd_s,fd_s
     read(1377,*)
     read(1377,*) q1
     read(1377,*)
     read(1377,*) a11,a12
     read(1377,*)
     read(1377,*) sd_s
     read(1377,*)
     read(1377,*) fd_s
     write(*,*) "params : ",q1,a11,a12,sd_s,fd_s
    allocate(Q(atoms%nat,lcn),A(atoms%nat,lcn))
    Q(1:atoms%nat,1)=q1
    A(1:4,1) = a11
    A(5:8,1) = a12
    allocate(cent_rho_1(cent%poisson%ngpx,cent%poisson%ngpy,cent%poisson%ngpz))
    allocate(cent_rho_2(cent%poisson%ngpx,cent%poisson%ngpy,cent%poisson%ngpz))
    allocate(cent%poisson%pot(cent%poisson%ngpx,cent%poisson%ngpy,&
        cent%poisson%ngpz),stat=istat)
    allocate(err_qgrad(atoms%nat,lcn),err_agrad(atoms%nat,lcn))
    
    do iter=1,huge(iter_max)
         
         cent_rho_1=0.d0
         cent%poisson%pot = 0.d0
         cent%poisson%rgcut=6.d0*maxval(abs(A))
         call acf_read(parini,'posinp_1.acf',1,atoms=atoms)
         call gauss_grid(parini,atoms%boundcond,.true.,atoms%nat,atoms%rat,atoms%cellvec,Q(:,1),&
                A(:,1),cent%poisson%rgcut,nx,ny,nz,cent%poisson%rho)
         cent_rho_1 = cent%poisson%rho
         atoms%boundcond='bulk'
         atoms%zat = 1.d0
         cent%gwit = 1.3d0
         call construct_ewald_bps(parini,atoms,cent%poisson)
         call psolver_bps(cent%poisson,atoms,ehartree)
         atoms%fat = 0.d0
         call gauss_force(parini,'bulk',atoms%nat,atoms%rat,atoms%cellvec,atoms%zat,cent%gwit,cent%poisson%rgcut,nx,ny,nz,cent%poisson%pot,atoms%fat)
         cent_fat_1 = atoms%fat
         call destruct_ewald_bps(cent%poisson)
         
         cent_rho_2=0.d0
         cent%poisson%pot = 0.d0
         cent%poisson%rgcut=6.d0*maxval(abs(A))
         call acf_read(parini,'posinp_2.acf',1,atoms=atoms)
         call gauss_grid(parini,atoms%boundcond,.true.,atoms%nat,atoms%rat,atoms%cellvec,Q(:,1),&
                A(:,1),cent%poisson%rgcut,nx,ny,nz,cent%poisson%rho)
         cent_rho_2 = cent%poisson%rho
         atoms%zat = 1.d0
         cent%gwit = 1.3d0
         atoms%boundcond='bulk'
         call construct_ewald_bps(parini,atoms,cent%poisson)
         call psolver_bps(cent%poisson,atoms,ehartree)
         atoms%fat = 0.d0
         call gauss_force(parini,'bulk',atoms%nat,atoms%rat,atoms%cellvec,atoms%zat,cent%gwit,cent%poisson%rgcut,nx,ny,nz,cent%poisson%pot,atoms%fat)
         cent_fat_2 = atoms%fat
         call destruct_ewald_bps(cent%poisson)
         rho_err = sqrt(sum((cent_rho_1-dft_rho_1)**2+(cent_rho_2-dft_rho_2)**2)/(nx*ny*nz))
         err_cent = dot_product((cent_fat_1(1,:)-dft_fat_1(1,:)),(cent_fat_1(1,:)-dft_fat_1(1,:)))&
                             + dot_product((cent_fat_2(1,:)-dft_fat_2(1,:)),(cent_fat_2(1,:)-dft_fat_2(1,:)))
        !//////////////////FDM_part///////////////////////////
         do i = 1 , atoms%nat
            do j = 1 , lcn
                 A(i,j) = A(i,j) + fd_s
                 cent_rho_1=0.d0
                 cent%poisson%pot = 0.d0
                 cent%poisson%rgcut=6.d0*maxval(abs(A))
                 call acf_read(parini,'posinp_1.acf',1,atoms=atoms)
                 call gauss_grid(parini,atoms%boundcond,.true.,atoms%nat,atoms%rat,atoms%cellvec,Q(:,1),&
                        A(:,1),cent%poisson%rgcut,nx,ny,nz,cent%poisson%rho)
                 cent_rho_1 = cent%poisson%rho
                 atoms%zat = 1.d0
                 cent%gwit = 1.3d0
                 atoms%boundcond='bulk'
                 call construct_ewald_bps(parini,atoms,cent%poisson)
                 call psolver_bps(cent%poisson,atoms,ehartree)
                 atoms%fat = 0.d0
                 call gauss_force(parini,'bulk',atoms%nat,atoms%rat,atoms%cellvec,atoms%zat,cent%gwit,cent%poisson%rgcut,nx,ny,nz,cent%poisson%pot,atoms%fat)
                 cent_fat_1 = atoms%fat
                 call destruct_ewald_bps(cent%poisson)
                 
                 cent_rho_2=0.d0
                 cent%poisson%pot = 0.d0
                 cent%poisson%rgcut=6.d0*maxval(abs(A))
                 call acf_read(parini,'posinp_2.acf',1,atoms=atoms)
                 call gauss_grid(parini,atoms%boundcond,.true.,atoms%nat,atoms%rat,atoms%cellvec,Q(:,1),&
                        A(:,1),cent%poisson%rgcut,nx,ny,nz,cent%poisson%rho)
                 cent_rho_2 = cent%poisson%rho
                 atoms%zat = 1.d0
                 cent%gwit = 1.3d0
                 atoms%boundcond='bulk'
                 call construct_ewald_bps(parini,atoms,cent%poisson)
                 call psolver_bps(cent%poisson,atoms,ehartree)
                 atoms%fat = 0.d0
                 call gauss_force(parini,'bulk',atoms%nat,atoms%rat,atoms%cellvec,atoms%zat,cent%gwit,cent%poisson%rgcut,nx,ny,nz,cent%poisson%pot,atoms%fat)
                 cent_fat_2 = atoms%fat
                 call destruct_ewald_bps(cent%poisson)
                         
                 err_fdm = 0.d0
                 err_fdm = err_fdm + dot_product((cent_fat_1(1,:)-dft_fat_1(1,:)),(cent_fat_1(1,:)-dft_fat_1(1,:)))&
                                   + dot_product((cent_fat_2(1,:)-dft_fat_2(1,:)),(cent_fat_2(1,:)-dft_fat_2(1,:)))
                 err_agrad(i,j) = (err_fdm - err_cent)/fd_s
                 A(i,j) = A(i,j) - fd_s
                 !/////////////////////////////////////////////////////////////////////////////
                 Q(i,j) = Q(i,j) + fd_s
                 cent_rho_1=0.d0
                 cent%poisson%pot = 0.d0
                 cent%poisson%rgcut=6.d0*maxval(abs(A))
                 call acf_read(parini,'posinp_1.acf',1,atoms=atoms)
                 call gauss_grid(parini,atoms%boundcond,.true.,atoms%nat,atoms%rat,atoms%cellvec,Q(:,1),&
                        A(:,1),cent%poisson%rgcut,nx,ny,nz,cent%poisson%rho)
                 cent_rho_1 = cent%poisson%rho
                 atoms%zat = 1.d0
                 cent%gwit = 1.3d0
                 atoms%boundcond='bulk'
                 call construct_ewald_bps(parini,atoms,cent%poisson)
                 call psolver_bps(cent%poisson,atoms,ehartree)
                 atoms%fat = 0.d0
                 call gauss_force(parini,'bulk',atoms%nat,atoms%rat,atoms%cellvec,atoms%zat,cent%gwit,cent%poisson%rgcut,nx,ny,nz,cent%poisson%pot,atoms%fat)
                 cent_fat_1 = atoms%fat
                 call destruct_ewald_bps(cent%poisson)
                 
                 cent_rho_2=0.d0
                 cent%poisson%pot = 0.d0
                 cent%poisson%rgcut=6.d0*maxval(abs(A))
                 call acf_read(parini,'posinp_2.acf',1,atoms=atoms)
                 call gauss_grid(parini,atoms%boundcond,.true.,atoms%nat,atoms%rat,atoms%cellvec,Q(:,1),&
                        A(:,1),cent%poisson%rgcut,nx,ny,nz,cent%poisson%rho)
                 cent_rho_2 = cent%poisson%rho
                 atoms%zat = 1.d0
                 cent%gwit = 1.3d0
                 atoms%boundcond='bulk'
                 call construct_ewald_bps(parini,atoms,cent%poisson)
                 call psolver_bps(cent%poisson,atoms,ehartree)
                 atoms%fat = 0.d0
                 call gauss_force(parini,'bulk',atoms%nat,atoms%rat,atoms%cellvec,atoms%zat,cent%gwit,cent%poisson%rgcut,nx,ny,nz,cent%poisson%pot,atoms%fat)
                 cent_fat_2 = atoms%fat
                 call destruct_ewald_bps(cent%poisson)
                         
                 err_fdm = 0.d0
                 err_fdm = err_fdm + dot_product((cent_fat_1(1,:)-dft_fat_1(1,:)),(cent_fat_1(1,:)-dft_fat_1(1,:)))&
                                   + dot_product((cent_fat_2(1,:)-dft_fat_2(1,:)),(cent_fat_2(1,:)-dft_fat_2(1,:)))
                 err_qgrad(i,j) = (err_fdm - err_cent)/fd_s
                 Q(i,j) = Q(i,j) - fd_s
                 !/////////////////////////////////////////////////////////////////////////////
             end do
         end do 
         exit_cond = sqrt((err_cent))
         !do jj = 1, lcn
            write(*,'(a,i4,2es14.5,2(2x,8f6.2))')"SD ",iter,exit_cond,rho_err,q,a
         !end do
         write(*,'(a,i4,2(2x,8es14.4))') "a-q-grad",iter,err_agrad,err_qgrad
         if (exit_cond < 1.d-6)then
            write(*,*) "max conversion reached"
            write(*,*) "qat : "
            write(*,*) Q
            write(*,*) "gwe : "
            write(*,*) A
            exit
         endif
         Q = Q - SD_S*err_qgrad
         A = A - SD_S*err_agrad
    enddo
 end subroutine best_charge_density_force
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine best_charge_density_energy(parini)
     use mod_interface
 !!--------------------------------------------------------------------------------------------------
     use mod_parini, only: typ_parini
     use mod_electrostatics, only: typ_poisson, typ_poisson
     use mod_atoms, only: typ_atoms
     use mod_ann, only: typ_cent, typ_ann_arr
     implicit none
     type(typ_parini), intent(in):: parini
     !local variables
     type(typ_poisson):: poisson ,poisson_rough
     type(typ_atoms):: atoms
     type(typ_cent):: cent
     type(typ_ann_arr):: ann_arr
     integer:: istat, igpx, igpy, igpz, iat, nx ,ny, nz,iter,iter_max,gweiat
     real(8):: cell(3), epot, rgcut_a, t1, t2, t3, t4, pi,h,rho_err
     real(8):: ehartree,exit_cond,SDA,ener,dft_ener,coeff
     real(8),allocatable::  gausswidth(:) , gwit(:), diff(:,:),bigforce(:,:),agrad(:)
     real(8),allocatable:: dft_rho(:,:,:) , cent_rho(:,:,:)
     pi=4.d0*atan(1.d0)
     call cube_read('electronic_density.cube',atoms,cent%poisson)
     allocate(dft_rho(cent%poisson%ngpx,cent%poisson%ngpy,cent%poisson%ngpz))
     dft_rho = cent%poisson%rho
     call acf_read(parini,'posinp.acf',1,atoms=atoms)
     atoms%boundcond='bulk'
     allocate(cent%poisson%pot(cent%poisson%ngpx,cent%poisson%ngpy,cent%poisson%ngpz),stat=istat)
     call construct_ewald_bps(parini,atoms,cent%poisson)
     call psolver_bps(cent%poisson,atoms,ehartree)
     dft_ener = ehartree
     write(*,*) "dft_Hartree : " , dft_ener
     nx=cent%poisson%ngpx
     ny=cent%poisson%ngpy
     nz=cent%poisson%ngpz
     call destruct_ewald_bps(cent%poisson)
     deallocate(cent%poisson%pot)
    !/////////////////CENT_ENERGY////////////////////////! 
     atoms%qat = 2.d0
     SDA  = 1.d-3
     allocate(cent%gwe(atoms%nat))
     cent%gwe(1:4)=1d0
     cent%gwe(1:8)=2d0
     allocate(cent_rho(cent%poisson%ngpx,cent%poisson%ngpy,cent%poisson%ngpz))
     do iter=1,huge(iter_max)
          cent%poisson%rgcut=6.d0*maxval(abs(cent%gwe))
          call gauss_grid(parini,atoms%boundcond,.true.,atoms%nat,atoms%rat,atoms%cellvec,atoms%qat,cent%gwe,cent%poisson%rgcut,nx,ny,nz,cent%poisson%rho)
          cent_rho = cent%poisson%rho
          rho_err = sqrt(sum((cent_rho-dft_rho)**2)/(cent%poisson%ngpx**3))
          atoms%boundcond='bulk'
          allocate(cent%poisson%pot(cent%poisson%ngpx,cent%poisson%ngpy,cent%poisson%ngpz),stat=istat)
          call construct_ewald_bps(parini,atoms,cent%poisson)
          call psolver_bps(cent%poisson,atoms,ehartree)
          ener = ehartree
          write(*,*) "Hartree : " , ener ;
          allocate(cent%rgrad(1:3,atoms%nat),cent%qgrad(atoms%nat),agrad(atoms%nat))
          call gauss_gradient_rzx(parini,atoms%boundcond,atoms%nat,atoms%rat,atoms%cellvec,atoms%qat,cent%gwe,cent%poisson%rgcut,nx,ny,nz,cent%poisson%pot,cent%rgrad,cent%qgrad,agrad)
          
          exit_cond =sqrt( (dft_ener-ener)**2)
          write(*,'(a,i4,2es14.5,2(2x,8f6.2))')"SD ",iter,exit_cond,rho_err,atoms%qat,cent%gwe
          if (exit_cond < 1.d-4)then
             write(*,*) "max conversion reached"
             write(*,*) "qat : "
             write(*,*) atoms%qat
             write(*,*) "gwe : "
             write(*,*) cent%gwe
             exit
          endif
          coeff = 2.d0*(ener - dft_ener)
          atoms%qat = atoms%qat - SDA*coeff*cent%qgrad
          cent%gwe = cent%gwe - SDA*coeff*agrad
          call destruct_ewald_bps(cent%poisson)
          deallocate(cent%poisson%pot,cent%rgrad,cent%qgrad,agrad)
     enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!CHECKING AGRAD    h = 1d-8
!!CHECKING AGRAD    atoms%qat = 1.d0
!!CHECKING AGRAD    gweiat = 8
!!CHECKING AGRAD    cent%gwe(gweiat) = cent%gwe(gweiat)+h
!!CHECKING AGRAD    call gauss_grid(parini,atoms%boundcond,.true.,atoms%nat,atoms%rat,atoms%cellvec,atoms%qat,cent%gwe,cent%poisson%rgcut,nx,ny,nz,cent%poisson%rho)
!!CHECKING AGRAD    allocate(cent%poisson%pot(cent%poisson%ngpx,cent%poisson%ngpy,cent%poisson%ngpz),stat=istat)
!!CHECKING AGRAD    call construct_ewald_bps(parini,atoms,cent%poisson)
!!CHECKING AGRAD    call psolver_bps(cent%poisson,atoms,ehartree)
!!CHECKING AGRAD    ener = ehartree
!!CHECKING AGRAD    write(*,*) "ener2" , ener
!!CHECKING AGRAD    call destruct_ewald_bps(cent%poisson)
!!CHECKING AGRAD    deallocate(cent%poisson%pot)
!!CHECKING AGRAD    cent%gwe(gweiat)=cent%gwe(gweiat)-(2.d0*h)
!!CHECKING AGRAD    call gauss_grid(parini,atoms%boundcond,.true.,atoms%nat,atoms%rat,atoms%cellvec,atoms%qat,cent%gwe,cent%poisson%rgcut,nx,ny,nz,cent%poisson%rho)
!!CHECKING AGRAD    !-----------------------------------------------------------------------------------
!!CHECKING AGRAD    allocate(cent%poisson%pot(cent%poisson%ngpx,cent%poisson%ngpy,cent%poisson%ngpz),stat=istat)
!!CHECKING AGRAD    call construct_ewald_bps(parini,atoms,cent%poisson)
!!CHECKING AGRAD    call psolver_bps(cent%poisson,atoms,ehartree)
!!CHECKING AGRAD    write(*,*) "ener3",ehartree
!!CHECKING AGRAD    ener = ener - ehartree
!!CHECKING AGRAD    ener = ener/(2.d0*h)
!!CHECKING AGRAD    write(*,*) "agrad(",gweiat,")",ener ;
!!CHECKING AGRAD    call destruct_ewald_bps(cent%poisson)
!!CHECKING AGRAD    deallocate(cent%poisson%pot)
!!CHECKING AGRAD    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine best_charge_density_energy

!*****************************************************************************************
