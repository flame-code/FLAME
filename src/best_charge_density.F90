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
    use mod_electrostatics, only: typ_poisson_p3d, typ_ewald_p3d,typ_poisson
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_cent, typ_ann_arr
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    type(typ_poisson_p3d):: poisson_p3d
    type(typ_ewald_p3d):: ewald_p3d ,ewald_p3d_rough
    type(typ_atoms):: atoms
    type(typ_cent):: cent
    type(typ_ann_arr):: ann_arr
    integer:: istat, igpx, igpy, igpz, iat, nx ,ny,nz,iter,iter_max,gweiat,i,j,k,jj,lcn,ierr,lwork ! lcn = linear combinarion number
    real(8):: cell(3), epot, rgcut_a, t1, t2, t3, t4, pi,h,rho_err,hx,hy,hz
    real(8):: ehartree,rmse,SDA,dft_ener,cent_ener,ener_err,force_err,ener_err_old,coeff,w
    real(8),allocatable:: dft_rho(:,:,:), cent_rho(:,:,:),cent_rho_1(:,:,:),cent_rho_2(:,:,:),weight(:,:,:),work(:),EE(:,:),eig(:)
    real(8),allocatable::A(:,:), Q(:,:), qpar(:,:), apar(:,:), qtemp(:)
    real(8),allocatable:: dft_fat(:,:),cent_fat(:,:)
    real(8):: sd_s,err_fdm,err_cent,volinv,err
    real(8):: q_max,a_max
    real(8):: errmax, peak, c1, d1, temp1,cv1,cv2,cv3
    !open(unit=1370,file='param.txt')
    call cube_read('electronic_density.cube',atoms,cent%ewald_p3d%poisson_p3d%typ_poisson)
    call acf_read(parini,'posinp.acf',1,atoms=atoms)
    
    open(unit=1377,file='param_2.txt')
    pi=4.d0*atan(1.d0)
    read(1377,*)
    read(1377,*) lcn
    allocate(Q(lcn,atoms%nat),A(lcn,atoms%nat))
    read(1377,*)
    do i = 1 , lcn
        read(1377,*) A(i,:) , Q(i,:)
    end do
    read(1377,*)
    read(1377,*) sd_s
    read(1377,*)
    read(1377,*) err


    write(*,'(a,i,3es13.6)') "params q_max and LCN : ",lcn,A(:,1)
    write(*,'(a,es13.6)') "params max gw : ",sum(Q(:,1))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!dft part!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,'(a,es13.6,i)') "params sd_s and # of atoms : ",sd_s,atoms%nat
    nx=cent%ewald_p3d%poisson_p3d%ngpx
    ny=cent%ewald_p3d%poisson_p3d%ngpy
    nz=cent%ewald_p3d%poisson_p3d%ngpz
    cv1=atoms%cellvec(1,1)
    cv2=atoms%cellvec(2,2)
    cv3=atoms%cellvec(3,3)
    !write(30,*) cv1,cv2,cv3 
    volinv = (cv1*cv2*cv3)/(nx*ny*nz)
    !invvol = 1.d0/(3.d0*atoms%nat)
    cent%ewald_p3d%hgx=cent%ewald_p3d%poisson_p3d%hx
    cent%ewald_p3d%hgy=cent%ewald_p3d%poisson_p3d%hy
    cent%ewald_p3d%hgz=cent%ewald_p3d%poisson_p3d%hz
    hx=cent%ewald_p3d%poisson_p3d%hx
    hy=cent%ewald_p3d%poisson_p3d%hy
    hz=cent%ewald_p3d%poisson_p3d%hz
    allocate(cent%ewald_p3d%poisson_p3d%pot(nx,ny,nz),stat=istat)
    allocate(dft_rho(nx,ny,nz),weight(nx,ny,nz),dft_fat(1:3,1:atoms%nat),cent%gwit(1:atoms%nat))
    cent%ewald_p3d%poisson_p3d%rho=-1.d0*cent%ewald_p3d%poisson_p3d%rho
    dft_rho = cent%ewald_p3d%poisson_p3d%rho
    w=abs(1.d0/(volinv*sum(dft_rho)))
    cent%gwit(:) = 1.3d0
    atoms%zat(:) = 1.d0
    atoms%boundcond='bulk'   
    cent%ewald_p3d%rgcut=6.d0*a_max !ask Dr. Ghasemi if its OK !
    dft_fat = 0.d0
    dft_ener= 0.d0 
    call construct_ewald_bps(parini,atoms,cent%ewald_p3d)
    call cal_hartree_pot_bps(cent%ewald_p3d,atoms,ehartree)
    dft_ener = ehartree 
    if (atoms%nat > 1) then
        call gauss_force(parini,'bulk',atoms%nat,atoms%rat,atoms%cellvec,atoms%zat,cent%gwit,cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%pot,atoms%fat)
        dft_fat = atoms%fat
    end if
    call destruct_ewald_bps(cent%ewald_p3d)
    write(*,'(a,es14.6,3es14.6)') "DFT" , dft_ener , dft_fat
    !deallocate(cent%ewald_p3d%poisson_p3d%pot)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!cent part!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
    allocate(cent_rho(nx,ny,nz))
    allocate(cent_rho_1(nx,ny,nz))
    allocate(cent_rho_2(nx,ny,nz))
    !allocate(cent%ewald_p3d%poisson_p3d%pot(cent%ewald_p3d%poisson_p3d%ngpx,cent%ewald_p3d%poisson_p3d%ngpy,cent%ewald_p3d%poisson_p3d%ngpz),stat=istat)
    allocate(qpar(1:lcn,1:atoms%nat),apar(1:lcn,1:atoms%nat))
    allocate(cent_fat(1:3,1:atoms%nat))
    !**************************************************
    
    !do i = 1 , lcn
    !    Q(i,1:atoms%nat)=q_max/(i+0.d0)
    !    A(i,1:atoms%nat)=a_max/(i+0.d0)
    !end do
    !A(atoms%nat/2:atoms%nat) = a12
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SD part!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  !eig_part!      call gauss_grid(parini,atoms%boundcond,.true.,1,atoms%rat(1:3,i),atoms%cellvec,qtemp(i),A(1,i),cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%rho)
  !eig_part!      cent_rho_1 = cent%ewald_p3d%poisson_p3d%rho
  !eig_part!      do j = i,atoms%nat
  !eig_part!          call gauss_grid(parini,atoms%boundcond,.true.,1,atoms%rat(1:3,j),atoms%cellvec,qtemp(j),A(1,j),cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%rho)
  !eig_part!          cent_rho_2 = cent%ewald_p3d%poisson_p3d%rho
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
    do iter=1,huge(iter_max)
        cent%ewald_p3d%poisson_p3d%rho = 0.d0
        cent%ewald_p3d%poisson_p3d%pot = 0.d0
        !atoms%boundcond='bulk'
        do i = 1 , lcn
            if (i==1) then         
                call gauss_grid(parini,atoms%boundcond,.true.,atoms%nat,atoms%rat,atoms%cellvec,Q(i,1:atoms%nat),A(i,1:atoms%nat),cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%rho)
            else
                call gauss_grid(parini,atoms%boundcond,.false.,atoms%nat,atoms%rat,atoms%cellvec,Q(i,1:atoms%nat),A(i,1:atoms%nat),cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%rho)
            end if
        end do
        cent_rho = cent%ewald_p3d%poisson_p3d%rho
        atoms%fat = 0.d0
        cent_fat = 0.d0
        call construct_ewald_bps(parini,atoms,cent%ewald_p3d)
        call cal_hartree_pot_bps(cent%ewald_p3d,atoms,ehartree)
        cent_ener = ehartree
        if(atoms%nat > 1)  then
            call gauss_force(parini,'bulk',atoms%nat,atoms%rat,atoms%cellvec,atoms%zat,cent%gwit,cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%pot,atoms%fat)
            cent_fat = atoms%fat
        end if
        call destruct_ewald_bps(cent%ewald_p3d)
        do i = 1 , atoms%nat
            do j = 1 , lcn
                call gauss_grid(parini,atoms%boundcond,.true.,1,atoms%rat(1:3,i),atoms%cellvec,qtemp(i),A(j,i),cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%rho)
                cent_rho_1 = cent%ewald_p3d%poisson_p3d%rho
                qpar(j,i) = 2.d0*volinv*sum(cent_rho_1*(cent_rho-dft_rho))
            end do 
        end do 
        !rho_err = sum((cent_rho-dft_rho)**2)*invvol
        !rmse = sqrt(rho_err)
        !!ener_err_old = ener_err
        ener_err = abs(cent_ener-dft_ener)
        force_err = sqrt(sum((cent_fat-dft_fat)**2)/(3*atoms%nat))
        errmax = 0.d0
        rho_err = 0.d0
        peak = maxval(abs(dft_rho))
        write(*,*) "peak",peak
        do i = 1 , nz
            do j = 1 , ny
                do k = 1 , nx
                    c1 = cent_rho(k,j,i)
                    d1 = dft_rho(k,j,i)
                    temp1 = (abs(c1)-abs(d1))/(max(abs(d1),(1.d-2*peak)))
                    !write(21,'(a,10es14.5)')"temp ",temp1,abs(c1),abs(d1),(1.d-2*peak),k*hx,j*hy,i*hz,atoms%rat
                    errmax=max(errmax,temp1)
                    rho_err = rho_err + ((c1-d1)**2)
                end do
            end do
        end do
        rho_err = sqrt(rho_err*volinv)*w
        !write(40,*)iter,rho_err,volinv,w
        rmse = sqrt(rho_err)
        !stop 'AAAAAAAAAAAA'
        write(*,'(a,i4,4es14.6)')"SD",iter,rmse,errmax,ener_err,force_err
        write(*,*)"SD --------------------------------------------------------------------"
        write(*,'(a,i4,i4,es14.6)')"GTO_val_sum ",iter,i,sum(q(:,1))
        do i = 1 , lcn
            write(*,'(a,i4,i4,2es14.6)')"GTO_val ",iter,i,q(i,1),a(i,1)
            write(*,'(a,i4,i4,es14.4)') "a-q-partial",iter,i,qpar(i,1) 
        enddo
        write(*,'(a)')"GTO --------------------------------------------------------------------"
        if (abs(errmax) < err)then
            write(*,*) "max conversion reached"
            exit
        endif
        !coeff = 2.d0*(sum(cent_rho - dft_rho))*invvol
        Q = Q - SD_S*qpar
        !!A(:,1) = A(:,1) - SD_S*coeff*agrad
       !!call destruct_ewald_bps(cent%ewald_p3d)
       !! write(*,*) "coeff , sd_s",coeff,sd_s
       !! if (ener_err > ener_err_old) then
       !!     sd_s = 0.95d0*sd_s
       !! else
       !!     sd_s = 1.05*sd_s
       !! endif
    end Do

  
end subroutine best_charge_density_rho
! !!--------------------------------------------------------------------------------------------------
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
subroutine best_charge_density_force(parini)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson_p3d, typ_ewald_p3d,typ_poisson
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_cent, typ_ann_arr
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    type(typ_poisson_p3d):: poisson_p3d
    type(typ_ewald_p3d):: ewald_p3d ,ewald_p3d_rough
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

    call cube_read('electronic_density_1.cube',atoms,cent%ewald_p3d%poisson_p3d%typ_poisson)
    allocate(dft_rho_1(cent%ewald_p3d%poisson_p3d%ngpx,cent%ewald_p3d%poisson_p3d%ngpy,cent%ewald_p3d%poisson_p3d%ngpz))
    cent%ewald_p3d%poisson_p3d%rho=-1.d0*cent%ewald_p3d%poisson_p3d%rho
    dft_rho_1 = -1.d0*cent%ewald_p3d%poisson_p3d%rho
    call acf_read(parini,'posinp_1.acf',1,atoms=atoms)

    cent%ewald_p3d%hgx=cent%ewald_p3d%poisson_p3d%hx
    cent%ewald_p3d%hgy=cent%ewald_p3d%poisson_p3d%hy
    cent%ewald_p3d%hgz=cent%ewald_p3d%poisson_p3d%hz
    atoms%boundcond='bulk'
    allocate(cent%ewald_p3d%poisson_p3d%pot(cent%ewald_p3d%poisson_p3d%ngpx,cent%ewald_p3d%poisson_p3d%ngpy,cent%ewald_p3d%poisson_p3d%ngpz),stat=istat)
    call construct_ewald_bps(parini,atoms,cent%ewald_p3d)
    call cal_hartree_pot_bps(cent%ewald_p3d,atoms,ehartree)
    atoms%zat(:)=1.d0
    nx=cent%ewald_p3d%poisson_p3d%ngpx
    ny=cent%ewald_p3d%poisson_p3d%ngpy
    nz=cent%ewald_p3d%poisson_p3d%ngpz
    allocate(cent%gwit(atoms%nat))
    allocate(dft_fat_1(lcn,atoms%nat))
    cent%gwit = 1.3d0
    atoms%fat = 0.d0
    dft_fat_1 = 0.d0
    call gauss_force(parini,'bulk',atoms%nat,atoms%rat,atoms%cellvec,atoms%zat,cent%gwit,cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%pot,atoms%fat)
    dft_fat_1 = atoms%fat
    call destruct_ewald_bps(cent%ewald_p3d)
    deallocate(cent%ewald_p3d%poisson_p3d%pot)
   !///////////////////////////2nd DFT////////////////////////////////////
    cent%ewald_p3d%poisson_p3d%rho = 0.d0
    call cube_read('electronic_density_2.cube',atoms,cent%ewald_p3d%poisson_p3d%typ_poisson)
    allocate(dft_rho_2(cent%ewald_p3d%poisson_p3d%ngpx,cent%ewald_p3d%poisson_p3d%ngpy,cent%ewald_p3d%poisson_p3d%ngpz))
    cent%ewald_p3d%poisson_p3d%rho=-1.d0*cent%ewald_p3d%poisson_p3d%rho
    dft_rho_2 = -1.d0*cent%ewald_p3d%poisson_p3d%rho
    call acf_read(parini,'posinp_2.acf',1,atoms=atoms)
    
    cent%ewald_p3d%hgx=cent%ewald_p3d%poisson_p3d%hx
    cent%ewald_p3d%hgy=cent%ewald_p3d%poisson_p3d%hy
    cent%ewald_p3d%hgz=cent%ewald_p3d%poisson_p3d%hz
    atoms%boundcond='bulk'
    allocate(cent%ewald_p3d%poisson_p3d%pot(cent%ewald_p3d%poisson_p3d%ngpx,cent%ewald_p3d%poisson_p3d%ngpy,cent%ewald_p3d%poisson_p3d%ngpz),stat=istat)
    call construct_ewald_bps(parini,atoms,cent%ewald_p3d)
    call cal_hartree_pot_bps(cent%ewald_p3d,atoms,ehartree)
    atoms%zat(:)=1.d0
    nx=cent%ewald_p3d%poisson_p3d%ngpx
    ny=cent%ewald_p3d%poisson_p3d%ngpy
    nz=cent%ewald_p3d%poisson_p3d%ngpz
    allocate(dft_fat_2(lcn,atoms%nat),cent_fat_1(lcn,atoms%nat),cent_fat_2(lcn,atoms%nat))
    cent%gwit = 1.3d0
    atoms%fat = 0.d0
    dft_fat_2 = 0.d0
    call gauss_force(parini,'bulk',atoms%nat,atoms%rat,atoms%cellvec,atoms%zat,cent%gwit,cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%pot,atoms%fat)
    dft_fat_2 = atoms%fat
    call destruct_ewald_bps(cent%ewald_p3d)
    deallocate(cent%ewald_p3d%poisson_p3d%pot)
 
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
    allocate(cent_rho_1(cent%ewald_p3d%poisson_p3d%ngpx,cent%ewald_p3d%poisson_p3d%ngpy,cent%ewald_p3d%poisson_p3d%ngpz))
    allocate(cent_rho_2(cent%ewald_p3d%poisson_p3d%ngpx,cent%ewald_p3d%poisson_p3d%ngpy,cent%ewald_p3d%poisson_p3d%ngpz))
    allocate(cent%ewald_p3d%poisson_p3d%pot(cent%ewald_p3d%poisson_p3d%ngpx,cent%ewald_p3d%poisson_p3d%ngpy,&
        cent%ewald_p3d%poisson_p3d%ngpz),stat=istat)
    allocate(err_qgrad(atoms%nat,lcn),err_agrad(atoms%nat,lcn))
    
    do iter=1,huge(iter_max)
         
         cent_rho_1=0.d0
         cent%ewald_p3d%poisson_p3d%pot = 0.d0
         cent%ewald_p3d%rgcut=6.d0*maxval(abs(A))
         call acf_read(parini,'posinp_1.acf',1,atoms=atoms)
         call gauss_grid(parini,atoms%boundcond,.true.,atoms%nat,atoms%rat,atoms%cellvec,Q(:,1),&
                A(:,1),cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%rho)
         cent_rho_1 = cent%ewald_p3d%poisson_p3d%rho
         atoms%boundcond='bulk'
         atoms%zat = 1.d0
         cent%gwit = 1.3d0
         call construct_ewald_bps(parini,atoms,cent%ewald_p3d)
         call cal_hartree_pot_bps(cent%ewald_p3d,atoms,ehartree)
         atoms%fat = 0.d0
         call gauss_force(parini,'bulk',atoms%nat,atoms%rat,atoms%cellvec,atoms%zat,cent%gwit,cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%pot,atoms%fat)
         cent_fat_1 = atoms%fat
         call destruct_ewald_bps(cent%ewald_p3d)
         
         cent_rho_2=0.d0
         cent%ewald_p3d%poisson_p3d%pot = 0.d0
         cent%ewald_p3d%rgcut=6.d0*maxval(abs(A))
         call acf_read(parini,'posinp_2.acf',1,atoms=atoms)
         call gauss_grid(parini,atoms%boundcond,.true.,atoms%nat,atoms%rat,atoms%cellvec,Q(:,1),&
                A(:,1),cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%rho)
         cent_rho_2 = cent%ewald_p3d%poisson_p3d%rho
         atoms%zat = 1.d0
         cent%gwit = 1.3d0
         atoms%boundcond='bulk'
         call construct_ewald_bps(parini,atoms,cent%ewald_p3d)
         call cal_hartree_pot_bps(cent%ewald_p3d,atoms,ehartree)
         atoms%fat = 0.d0
         call gauss_force(parini,'bulk',atoms%nat,atoms%rat,atoms%cellvec,atoms%zat,cent%gwit,cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%pot,atoms%fat)
         cent_fat_2 = atoms%fat
         call destruct_ewald_bps(cent%ewald_p3d)
         rho_err = sqrt(sum((cent_rho_1-dft_rho_1)**2+(cent_rho_2-dft_rho_2)**2)/(nx*ny*nz))
         err_cent = dot_product((cent_fat_1(1,:)-dft_fat_1(1,:)),(cent_fat_1(1,:)-dft_fat_1(1,:)))&
                             + dot_product((cent_fat_2(1,:)-dft_fat_2(1,:)),(cent_fat_2(1,:)-dft_fat_2(1,:)))
        !//////////////////FDM_part///////////////////////////
         do i = 1 , atoms%nat
            do j = 1 , lcn
                 A(i,j) = A(i,j) + fd_s
                 cent_rho_1=0.d0
                 cent%ewald_p3d%poisson_p3d%pot = 0.d0
                 cent%ewald_p3d%rgcut=6.d0*maxval(abs(A))
                 call acf_read(parini,'posinp_1.acf',1,atoms=atoms)
                 call gauss_grid(parini,atoms%boundcond,.true.,atoms%nat,atoms%rat,atoms%cellvec,Q(:,1),&
                        A(:,1),cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%rho)
                 cent_rho_1 = cent%ewald_p3d%poisson_p3d%rho
                 atoms%zat = 1.d0
                 cent%gwit = 1.3d0
                 atoms%boundcond='bulk'
                 call construct_ewald_bps(parini,atoms,cent%ewald_p3d)
                 call cal_hartree_pot_bps(cent%ewald_p3d,atoms,ehartree)
                 atoms%fat = 0.d0
                 call gauss_force(parini,'bulk',atoms%nat,atoms%rat,atoms%cellvec,atoms%zat,cent%gwit,cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%pot,atoms%fat)
                 cent_fat_1 = atoms%fat
                 call destruct_ewald_bps(cent%ewald_p3d)
                 
                 cent_rho_2=0.d0
                 cent%ewald_p3d%poisson_p3d%pot = 0.d0
                 cent%ewald_p3d%rgcut=6.d0*maxval(abs(A))
                 call acf_read(parini,'posinp_2.acf',1,atoms=atoms)
                 call gauss_grid(parini,atoms%boundcond,.true.,atoms%nat,atoms%rat,atoms%cellvec,Q(:,1),&
                        A(:,1),cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%rho)
                 cent_rho_2 = cent%ewald_p3d%poisson_p3d%rho
                 atoms%zat = 1.d0
                 cent%gwit = 1.3d0
                 atoms%boundcond='bulk'
                 call construct_ewald_bps(parini,atoms,cent%ewald_p3d)
                 call cal_hartree_pot_bps(cent%ewald_p3d,atoms,ehartree)
                 atoms%fat = 0.d0
                 call gauss_force(parini,'bulk',atoms%nat,atoms%rat,atoms%cellvec,atoms%zat,cent%gwit,cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%pot,atoms%fat)
                 cent_fat_2 = atoms%fat
                 call destruct_ewald_bps(cent%ewald_p3d)
                         
                 err_fdm = 0.d0
                 err_fdm = err_fdm + dot_product((cent_fat_1(1,:)-dft_fat_1(1,:)),(cent_fat_1(1,:)-dft_fat_1(1,:)))&
                                   + dot_product((cent_fat_2(1,:)-dft_fat_2(1,:)),(cent_fat_2(1,:)-dft_fat_2(1,:)))
                 err_agrad(i,j) = (err_fdm - err_cent)/fd_s
                 A(i,j) = A(i,j) - fd_s
                 !/////////////////////////////////////////////////////////////////////////////
                 Q(i,j) = Q(i,j) + fd_s
                 cent_rho_1=0.d0
                 cent%ewald_p3d%poisson_p3d%pot = 0.d0
                 cent%ewald_p3d%rgcut=6.d0*maxval(abs(A))
                 call acf_read(parini,'posinp_1.acf',1,atoms=atoms)
                 call gauss_grid(parini,atoms%boundcond,.true.,atoms%nat,atoms%rat,atoms%cellvec,Q(:,1),&
                        A(:,1),cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%rho)
                 cent_rho_1 = cent%ewald_p3d%poisson_p3d%rho
                 atoms%zat = 1.d0
                 cent%gwit = 1.3d0
                 atoms%boundcond='bulk'
                 call construct_ewald_bps(parini,atoms,cent%ewald_p3d)
                 call cal_hartree_pot_bps(cent%ewald_p3d,atoms,ehartree)
                 atoms%fat = 0.d0
                 call gauss_force(parini,'bulk',atoms%nat,atoms%rat,atoms%cellvec,atoms%zat,cent%gwit,cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%pot,atoms%fat)
                 cent_fat_1 = atoms%fat
                 call destruct_ewald_bps(cent%ewald_p3d)
                 
                 cent_rho_2=0.d0
                 cent%ewald_p3d%poisson_p3d%pot = 0.d0
                 cent%ewald_p3d%rgcut=6.d0*maxval(abs(A))
                 call acf_read(parini,'posinp_2.acf',1,atoms=atoms)
                 call gauss_grid(parini,atoms%boundcond,.true.,atoms%nat,atoms%rat,atoms%cellvec,Q(:,1),&
                        A(:,1),cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%rho)
                 cent_rho_2 = cent%ewald_p3d%poisson_p3d%rho
                 atoms%zat = 1.d0
                 cent%gwit = 1.3d0
                 atoms%boundcond='bulk'
                 call construct_ewald_bps(parini,atoms,cent%ewald_p3d)
                 call cal_hartree_pot_bps(cent%ewald_p3d,atoms,ehartree)
                 atoms%fat = 0.d0
                 call gauss_force(parini,'bulk',atoms%nat,atoms%rat,atoms%cellvec,atoms%zat,cent%gwit,cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%pot,atoms%fat)
                 cent_fat_2 = atoms%fat
                 call destruct_ewald_bps(cent%ewald_p3d)
                         
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
     use mod_electrostatics, only: typ_poisson_p3d, typ_ewald_p3d,typ_poisson
     use mod_atoms, only: typ_atoms
     use mod_ann, only: typ_cent, typ_ann_arr
     implicit none
     type(typ_parini), intent(in):: parini
     !local variables
     type(typ_poisson_p3d):: poisson_p3d
     type(typ_ewald_p3d):: ewald_p3d ,ewald_p3d_rough
     type(typ_atoms):: atoms
     type(typ_cent):: cent
     type(typ_ann_arr):: ann_arr
     integer:: istat, igpx, igpy, igpz, iat, nx ,ny, nz,iter,iter_max,gweiat
     real(8):: cell(3), epot, rgcut_a, t1, t2, t3, t4, pi,h,rho_err
     real(8):: ehartree,exit_cond,SDA,ener,dft_ener,coeff
     real(8),allocatable::  gausswidth(:) , gwit(:), diff(:,:),bigforce(:,:),agrad(:)
     real(8),allocatable:: dft_rho(:,:,:) , cent_rho(:,:,:)
     pi=4.d0*atan(1.d0)
     call cube_read('electronic_density.cube',atoms,cent%ewald_p3d%poisson_p3d%typ_poisson)
     allocate(dft_rho(cent%ewald_p3d%poisson_p3d%ngpx,cent%ewald_p3d%poisson_p3d%ngpy,cent%ewald_p3d%poisson_p3d%ngpz))
     dft_rho = cent%ewald_p3d%poisson_p3d%rho
     call acf_read(parini,'posinp.acf',1,atoms=atoms)
     cent%ewald_p3d%hgx=cent%ewald_p3d%poisson_p3d%hx
     cent%ewald_p3d%hgy=cent%ewald_p3d%poisson_p3d%hy
     cent%ewald_p3d%hgz=cent%ewald_p3d%poisson_p3d%hz
     atoms%boundcond='bulk'
     allocate(cent%ewald_p3d%poisson_p3d%pot(cent%ewald_p3d%poisson_p3d%ngpx,cent%ewald_p3d%poisson_p3d%ngpy,cent%ewald_p3d%poisson_p3d%ngpz),stat=istat)
     call construct_ewald_bps(parini,atoms,cent%ewald_p3d)
     call cal_hartree_pot_bps(cent%ewald_p3d,atoms,ehartree)
     dft_ener = ehartree
     write(*,*) "dft_Hartree : " , dft_ener
     nx=cent%ewald_p3d%poisson_p3d%ngpx
     ny=cent%ewald_p3d%poisson_p3d%ngpy
     nz=cent%ewald_p3d%poisson_p3d%ngpz
     call destruct_ewald_bps(cent%ewald_p3d)
     deallocate(cent%ewald_p3d%poisson_p3d%pot)
    !/////////////////CENT_ENERGY////////////////////////! 
     atoms%qat = 2.d0
     SDA  = 1.d-3
     allocate(cent%gwe(atoms%nat))
     cent%gwe(1:4)=1d0
     cent%gwe(1:8)=2d0
     allocate(cent_rho(cent%ewald_p3d%poisson_p3d%ngpx,cent%ewald_p3d%poisson_p3d%ngpy,cent%ewald_p3d%poisson_p3d%ngpz))
     do iter=1,huge(iter_max)
          cent%ewald_p3d%rgcut=6.d0*maxval(abs(cent%gwe))
          call gauss_grid(parini,atoms%boundcond,.true.,atoms%nat,atoms%rat,atoms%cellvec,atoms%qat,cent%gwe,cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%rho)
          cent_rho = cent%ewald_p3d%poisson_p3d%rho
          rho_err = sqrt(sum((cent_rho-dft_rho)**2)/(cent%ewald_p3d%poisson_p3d%ngpx**3))
          atoms%boundcond='bulk'
          allocate(cent%ewald_p3d%poisson_p3d%pot(cent%ewald_p3d%poisson_p3d%ngpx,cent%ewald_p3d%poisson_p3d%ngpy,cent%ewald_p3d%poisson_p3d%ngpz),stat=istat)
          call construct_ewald_bps(parini,atoms,cent%ewald_p3d)
          call cal_hartree_pot_bps(cent%ewald_p3d,atoms,ehartree)
          ener = ehartree
          write(*,*) "Hartree : " , ener ;
          allocate(cent%rgrad(1:3,atoms%nat),cent%qgrad(atoms%nat),agrad(atoms%nat))
          call gauss_gradient_rzx(parini,atoms%boundcond,atoms%nat,atoms%rat,atoms%cellvec,atoms%qat,cent%gwe,cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%pot,cent%rgrad,cent%qgrad,agrad)
          
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
          call destruct_ewald_bps(cent%ewald_p3d)
          deallocate(cent%ewald_p3d%poisson_p3d%pot,cent%rgrad,cent%qgrad,agrad)
     enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!CHECKING AGRAD    h = 1d-8
!!CHECKING AGRAD    atoms%qat = 1.d0
!!CHECKING AGRAD    gweiat = 8
!!CHECKING AGRAD    cent%gwe(gweiat) = cent%gwe(gweiat)+h
!!CHECKING AGRAD    call gauss_grid(parini,atoms%boundcond,.true.,atoms%nat,atoms%rat,atoms%cellvec,atoms%qat,cent%gwe,cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%rho)
!!CHECKING AGRAD    allocate(cent%ewald_p3d%poisson_p3d%pot(cent%ewald_p3d%poisson_p3d%ngpx,cent%ewald_p3d%poisson_p3d%ngpy,cent%ewald_p3d%poisson_p3d%ngpz),stat=istat)
!!CHECKING AGRAD    call construct_ewald_bps(parini,atoms,cent%ewald_p3d)
!!CHECKING AGRAD    call cal_hartree_pot_bps(cent%ewald_p3d,atoms,ehartree)
!!CHECKING AGRAD    ener = ehartree
!!CHECKING AGRAD    write(*,*) "ener2" , ener
!!CHECKING AGRAD    call destruct_ewald_bps(cent%ewald_p3d)
!!CHECKING AGRAD    deallocate(cent%ewald_p3d%poisson_p3d%pot)
!!CHECKING AGRAD    cent%gwe(gweiat)=cent%gwe(gweiat)-(2.d0*h)
!!CHECKING AGRAD    call gauss_grid(parini,atoms%boundcond,.true.,atoms%nat,atoms%rat,atoms%cellvec,atoms%qat,cent%gwe,cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%rho)
!!CHECKING AGRAD    !-----------------------------------------------------------------------------------
!!CHECKING AGRAD    allocate(cent%ewald_p3d%poisson_p3d%pot(cent%ewald_p3d%poisson_p3d%ngpx,cent%ewald_p3d%poisson_p3d%ngpy,cent%ewald_p3d%poisson_p3d%ngpz),stat=istat)
!!CHECKING AGRAD    call construct_ewald_bps(parini,atoms,cent%ewald_p3d)
!!CHECKING AGRAD    call cal_hartree_pot_bps(cent%ewald_p3d,atoms,ehartree)
!!CHECKING AGRAD    write(*,*) "ener3",ehartree
!!CHECKING AGRAD    ener = ener - ehartree
!!CHECKING AGRAD    ener = ener/(2.d0*h)
!!CHECKING AGRAD    write(*,*) "agrad(",gweiat,")",ener ;
!!CHECKING AGRAD    call destruct_ewald_bps(cent%ewald_p3d)
!!CHECKING AGRAD    deallocate(cent%ewald_p3d%poisson_p3d%pot)
!!CHECKING AGRAD    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine best_charge_density_energy

!*****************************************************************************************
subroutine gauss_gradient_rzx(parini,bc,nat,rxyz,cv,qat,gw,rgcut,ngx,ngy,ngz,pot,rgrad,qgrad,agrad)
    use mod_interface
    use mod_atoms, only: typ_atoms
    use mod_parini, only: typ_parini
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: bc
    integer, intent(in):: nat
    real(8), intent(in):: rxyz(3,nat)
    real(8), intent(in):: cv(3,3)
    real(8), intent(in):: qat(nat) 
    real(8), intent(in):: gw(nat)
    real(8), intent(in):: rgcut
    integer, intent(in):: ngx, ngy, ngz
    real(8), intent(inout):: pot(ngx,ngy,ngz)
    real(8), intent(out):: rgrad(3,nat), qgrad(nat), agrad(nat)
    !local variables
    !work arrays to save the values of one dimensional gaussian function.
    real(8):: pi
    real(8):: facqiat, fac
    real(8):: cell(3) !dimensions of a smaller orthogonal cell for replication
    real(8):: vol
    real(8):: cvinv(3) !cell vectors of inverse coordinate, actual one at a time
    real(8):: htx, hty, htz
    real(8):: hxx, hxy, hxz, hyx, hyy, hyz, hzx, hzy, hzz
    real(8):: hxx_g, hxy_g, hxz_g, hyx_g, hyy_g, hyz_g, hzx_g, hzy_g, hzz_g
    real(8):: hrxinv, hryinv, hrzinv
    real(8):: cvinv_norm, vol_voxel, ttx, tty, ttz, ttq, tt1, tta
    real(8):: dmsq, gwsq_inv, gw_inv
    real(8):: xred, yred, zred
    real(8):: ximg, yimg, zimg
    integer:: imgx, imgy, imgz
    integer:: ncellx, ncelly, ncellz
    integer:: iat, igx, igy, igz, jgx, jgy, jgz, igyt, igzt, igzysq
    integer:: iii, nlim, nlimsq, iix, iiy, iiz
    integer:: nbgx, nbgy, nbgz, nagx, nagy, nagz, nex, ney, nez
    integer:: finalx, igxs, igxf 
    real(8), allocatable:: wa(:,:,:)
    real(8), allocatable:: ratred(:,:)
    real(8), allocatable:: exponentval(:), expval(:)
    real(8), allocatable:: dmxarr(:), dmyarr(:), dmzarr(:)
    real(8):: sqpi, tt
    integer:: istartx, istarty, istartz

    allocate(ratred(3,nat))
    call rxyz_cart2int_alborz(nat,cv,rxyz,ratred)
    do iat=1,nat
        xred=ratred(1,iat)
        yred=ratred(2,iat)
        zred=ratred(3,iat)
        if(xred<0.d0) write(*,*) 'ATOM OUTSIDE: iat,sx ',iat,xred
        if(yred<0.d0) write(*,*) 'ATOM OUTSIDE: iat,sy ',iat,yred
        if(zred<0.d0) write(*,*) 'ATOM OUTSIDE: iat,sz ',iat,zred
        if(.not. (xred<1.d0)) write(*,*) 'ATOM OUTSIDE: iat,sx ',iat,xred
        if(.not. (yred<1.d0)) write(*,*) 'ATOM OUTSIDE: iat,sy ',iat,yred
        if(.not. (zred<1.d0)) write(*,*) 'ATOM OUTSIDE: iat,sz ',iat,zred
    enddo

    !reciprocal lattice to be used to determine the distance of corners of
    !the parallelepiped to its facets. Then those distances are used to
    !determine the number of grid points in each direction that are within
    !the cutoff of Gaussian function.
    call cell_vol(nat,cv,vol)
    vol=abs(vol)*nat
    call cross_product_alborz(cv(1,1),cv(1,2),cvinv)
    cvinv_norm=sqrt(cvinv(1)**2+cvinv(2)**2+cvinv(3)**2)
    cell(3)=vol/cvinv_norm
    call cross_product_alborz(cv(1,2),cv(1,3),cvinv)
    cvinv_norm=sqrt(cvinv(1)**2+cvinv(2)**2+cvinv(3)**2)
    cell(1)=vol/cvinv_norm
    call cross_product_alborz(cv(1,1),cv(1,3),cvinv)
    cvinv_norm=sqrt(cvinv(1)**2+cvinv(2)**2+cvinv(3)**2)
    cell(2)=vol/cvinv_norm
    if(parini%iverbose>1) then
        write(*,*) 'cell  ', cell(1),cell(2),cell(3)
    endif
    htx=cell(1)/real(ngx,8)
    hty=cell(2)/real(ngy,8)
    htz=cell(3)/real(ngz,8)
    nbgx=int(rgcut/htx)+2
    nbgy=int(rgcut/hty)+2
    nbgz=int(rgcut/htz)+2
    nagx=nbgx+1
    nagy=nbgy+1
    nagz=nbgz+1
    nlim = max(nagx,nagy,nagz)
    nlimsq = nlim**2
    !detemining the largest dimension for the pseudogrid.
    hxx=cv(1,1)/ngx ; hxy=cv(2,1)/ngx ; hxz=cv(3,1)/ngx
    hyx=cv(1,2)/ngy ; hyy=cv(2,2)/ngy ; hyz=cv(3,2)/ngy
    hzx=cv(1,3)/ngz ; hzy=cv(2,3)/ngz ; hzz=cv(3,3)/ngz

    wa=f_malloc0([1-nagx.to.ngx+nagx,1-nagy.to.ngy+nagy,1-nagz.to.ngz+nagz],id='wa')
    allocate(exponentval(-nbgx:nbgx),expval(-nbgx:nbgx))
    allocate(dmxarr(-nbgx:nbgx),dmyarr(-nbgx:nbgx),dmzarr(-nbgx:nbgx))
   !  allocate(wa(1-nagx:ngx+nagx,1-nagy:ngy+nagy,1-nagz:ngz+nagz))

    hrxinv=real(ngx,8) !inverse of grid spacing in reduced coordinates
    hryinv=real(ngy,8) !inverse of grid spacing in reduced coordinates
    hrzinv=real(ngz,8) !inverse of grid spacing in reduced coordinates
    !if(trim(bc)=='bulk') then
    !    iii=0
    !elseif(trim(bc)=='slab') then
    !    iii=1
    !endif
    !-------------------------------------------------------
    pi=4.d0*atan(1.d0)
    !---------------------------------------------------------------------------
    istartx = modulo((1-nagx-1),ngx)+1
    finalx = modulo((ngx+nagx-1),ngx)+1
    igxs = 1-nagx+(ngx-istartx+1)
    igxf = ngx+nagx-finalx+1
    istarty = modulo((-nagy),ngy)+1
    istartz = modulo((-nagz),ngz)+1
    iiz=istartz-1

    do igz=1-nagz,ngz+nagz
        iiy=istarty-1
        iiz=iiz+1
        if (iiz==ngz+1) iiz=1
        do igy=1-nagy,ngy+nagy
            iiy=iiy+1
            if (iiy==ngy+1) iiy=1
            iix=istartx-1
            do igx=1-nagx,igxs-1
                iix=iix+1
                wa(igx,igy,igz)=pot(iix,iiy,iiz)
            enddo
            do igx=igxs,igxf-1,ngx
                do iix=1,ngx
                    jgx=igx+iix-1
                    wa(jgx,igy,igz)=pot(iix,iiy,iiz)
                enddo
            enddo
            iix=0
            do igx=igxf,ngx+nagx
                iix=iix+1
                wa(igx,igy,igz)=pot(iix,iiy,iiz)
            enddo
        enddo
    enddo
    !-------------------------------------------------------
    vol_voxel=vol/(ngx*ngy*ngz)
    do iat=1,nat
        gw_inv = 1.d0/gw(iat)
        gwsq_inv=1.d0/gw(iat)**2
        fac=1.d0/(gw(iat)*sqrt(pi))**3
        imgx=nint(ratred(1,iat)*hrxinv)+1
        imgy=nint(ratred(2,iat)*hryinv)+1
        imgz=nint(ratred(3,iat)*hrzinv)+1
        facqiat=fac*qat(iat)
        ttq=0.d0
        ttx=0.d0
        tty=0.d0
        ttz=0.d0
        !------------------------
        tta=0.d0
        do igz=-nbgz,nbgz
            jgz=igz+imgz
            hzx_g=(jgz-1)*hzx
            hzy_g=(jgz-1)*hzy
            hzz_g=(jgz-1)*hzz
            do igy=-nbgy,nbgy
                igzysq=igz**2+igy**2
                if(igzysq>nlimsq) cycle
                jgy=igy+imgy
                hyx_g=(jgy-1)*hyx
                hyy_g=(jgy-1)*hyy
                hyz_g=(jgy-1)*hyz
                tt=nlimsq-igzysq
                iix=min(floor(sqrt(tt)),nbgx)
                do igx=-iix,iix
                    jgx=igx+imgx
                    ximg=(jgx-1)*hxx+hyx_g+hzx_g    
                    yimg=(jgx-1)*hxy+hyy_g+hzy_g
                    zimg=(jgx-1)*hxz+hyz_g+hzz_g
                    dmxarr(igx)=ximg-rxyz(1,iat)
                    dmyarr(igx)=yimg-rxyz(2,iat)
                    dmzarr(igx)=zimg-rxyz(3,iat)
                    dmsq=dmxarr(igx)**2+dmyarr(igx)**2+dmzarr(igx)**2
                    exponentval(igx)=-dmsq*gwsq_inv
                enddo
                call vdexp( 2*iix+1, exponentval(-iix), expval(-iix) )
                do igx=-iix,iix
                    ttq=ttq+fac*expval(igx)*wa(igx+imgx,jgy,jgz)
                    tt1=facqiat*expval(igx)*wa(igx+imgx,jgy,jgz)*(2.d0*gwsq_inv)
                    ttx=ttx+tt1*dmxarr(igx)
                    tty=tty+tt1*dmyarr(igx)
                    ttz=ttz+tt1*dmzarr(igx)
                    !-----------------------
                    tta=tta+(-3.d0 - (2.d0*exponentval(igx)))*facqiat*gw_inv*expval(igx)*wa(igx+imgx,jgy,jgz)

                    !_______________________
                enddo
            enddo
        enddo
        qgrad(iat)=qgrad(iat)+ttq*vol_voxel
        rgrad(1,iat)=rgrad(1,iat)+ttx*vol_voxel
        rgrad(2,iat)=rgrad(2,iat)+tty*vol_voxel
        rgrad(3,iat)=rgrad(3,iat)+ttz*vol_voxel
        !-----------------------
        agrad(iat) = agrad(iat)+tta*vol_voxel
        !_______________________
    enddo
    !---------------------------------------------------------------------------
    deallocate(ratred)
    deallocate(exponentval,expval)
    deallocate(dmxarr,dmyarr,dmzarr)
    call f_free(wa)
end subroutine gauss_gradient_rzx
!*****************************************************************************************

