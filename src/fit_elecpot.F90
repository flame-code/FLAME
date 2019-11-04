subroutine subtask_fit_elecpot(parini)
    use mod_parini, only: typ_parini
    implicit none
    type(typ_parini), intent(in):: parini

    call fit_elecpot(parini)

end subroutine subtask_fit_elecpot 
!*****************************************************************************************
!/////////////////////////////////////////////////////////////////////////////////////////
!/////////////////////////////////////////////////////////////////////////////////////////
!/////////////////////////////////////////////////////////////////////////////////////////
!///////////////////////////////////// POT_PART //////////////////////////////////////////
!/////////////////////////////////////////////////////////////////////////////////////////
!/////////////////////////////////////////////////////////////////////////////////////////
subroutine fit_elecpot(parini)
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms, typ_atoms_arr, get_rat, set_rat, update_ratp
    use mod_atoms, only: atom_copy_old, atom_deallocate, set_rcov, get_rat_iat
    use mod_ann, only: typ_cent, typ_ann_arr
    use mod_yaml_conf, only: read_yaml_conf
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    type(typ_poisson):: poisson_dft, poisson_cent
    type(typ_atoms):: atoms
    type(typ_atoms_arr):: atoms_arr
    type(typ_cent):: cent
    logical:: esc
    integer:: atoms_type, iter, iter_max, cel_vol ,cel_vol_inv, ix, iy, iz, l, nclx, ncrx,ncly, ncry, nclz, ncrz ,nbgx, nbgy, nbgz
    integer:: at_range ,xg_at , yg_at , zg_at , n_at, i, j, k, lcn, ss, Ne, isatur, nsatur, ng(3)! lcn = linear combinarion number
    integer, allocatable :: t_num(:),zat_temp(:)
    real(8):: rmse_old , total_charge, start, finish, total_time, mean, std, var, hgx, hgy, hgz, total_force_cent, total_force_dft, force_err
    real(8):: rclx, rcrx, rcly, rcry, rclz, rcrz, rminx, rmaxx, rminy, rmaxy, rminz, rmaxz, pi,pot_err , temp1 
    real(8):: ehartree, rmse, cent_ener, dft_ener, w, exc, volinv_integration, SD_S_Q, SD_S_R, SD_S_A, volinv,err, a_max,qnrm, dft_q,gtot, peak, c1, d1, err_max, temp, at_zat
    real(8),allocatable::  dft_rho(:,:,:), weight(:,:,:),weighted_rho(:,:,:), dft_pot(:,:,:), cent_pot(:,:,:)
    real(8),allocatable:: A(:,:), Q(:,:), At(:,:), Qt(:,:), qpar(:,:), apar(:,:),rpar(:,:), qpar_t(:,:), rat_t(:,:),rat(:,:)
    real(8),allocatable:: dft_fat(:,:), cent_fat(:,:), gwit(:), at_rat(:,:),cent_rho(:,:,:),cv_temp(:,:), atom_charge(:),atom_charge_type(:)
    real(8),allocatable:: mean_arr(:), std_arr(:), a_type(:,:,:),q_type(:,:,:)
    integer,allocatable:: natom_type(:)
    real(8) :: B1, B2, h , Be1 , Be2, Beta1, Beta2
    integer :: iter2, iter3, iat
    real(8) :: dx, dy, dz, dr, dxyz(3)
    real(8) :: dpm_err(3), qpm_err(3,3)
    !////////////////////////E.O. TEMP for charge on grid/////////////////////////////////
    pi=4.d0*atan(1.d0)
    nsatur = 5 
    rmse =0.d0
    rmse_old=0.d0
    open(unit=3,file='pot.plot')
    open(unit=4,file='rho.plot')
    !open(unit=2,file='control')
    !open(unit=1377,file='params.txt')
!//////////////////////////////////////READING INPUT PARAMETERS///////////////////////////
    call cube_read('electronic_density.cube',atoms,poisson_dft)
    poisson_dft%rho=-1.d0*poisson_dft%rho
    call calc_multipoles_grid_centt(parini,atoms,poisson_dft)
    call yaml_map('DFT electric dpm',poisson_dft%dpm,fmt='(f10.4)')
    call yaml_map('DFT electric qpm',poisson_dft%qpm(1:3,1:3),fmt='(f10.4)')
    allocate(at_rat(3,atoms%nat),rat(3,atoms%nat),cv_temp(3,3),zat_temp(1:atoms%nat))
    zat_temp=atoms%zat
    !at_rat = atoms%rat ! because acf_read corrupts atoms%rat
    call get_rat(atoms,at_rat)
    hgx = poisson_dft%hgrid(1,1)
    hgy = poisson_dft%hgrid(2,2)
    hgz = poisson_dft%hgrid(3,3)
    cv_temp=atoms%cellvec
    call yaml_mapping_open('grid spacing',flow=.true.)
    call yaml_map('hx',hgx,fmt='(f6.3)')
    call yaml_map('hy',hgy,fmt='(f6.3)')
    call yaml_map('hz',hgz,fmt='(f6.3)')
    call yaml_mapping_close()
    
    !call acf_read(parini,'posinp.acf',1,atoms=atoms)
    call read_yaml_conf(parini,'posinp.yaml',1,atoms_arr)
    if(atoms_arr%nconf/=1) stop 'ERROR: atoms_arr%nconf/=1 in testforces_fd'
    call atom_copy_old(atoms_arr%atoms(1),atoms,'atoms_arr%atoms(iconf)->atoms_s')
    call atom_deallocate(atoms_arr%atoms(1))
    deallocate(atoms_arr%atoms)

    atoms%zat=zat_temp
    !atoms%rat = at_rat
    call set_rat(atoms,at_rat,setall=.true.)
    !rat(1,:) = atoms%rat(1,:)
    !rat(2,:) = atoms%rat(2,:)
    !rat(3,:) = atoms%rat(3,:)
    call get_rat(atoms,rat)
    atoms%cellvec=cv_temp
    !read(1377,*) !Number of STO-NG's N
    !read(1377,*) lcn
    lcn=parini%lcn
    !read(1377,*) !Number of atoms types
    !read(1377,*) atoms_type
    atoms_type=parini%ntypat
    allocate(t_num(0:atoms_type),Qt(lcn,atoms_type),At(lcn,atoms_type),rat_t(1:3,atoms_type),Q(1:lcn,1:atoms%nat),A(1:lcn,1:atoms%nat))
    allocate(mean_arr(lcn),std_arr(lcn))
    t_num(0) = 0
    !read(1377,*) !Number of each atom type 
    t_num(1:parini%ntypat)=0
    do iat=1,atoms%nat
        do i=1,parini%ntypat
            if(trim(atoms%sat(iat))==parini%stypat(i)) then
                t_num(i)=t_num(i)+1
            endif
        enddo
    enddo
    !write(*,*) atoms%sat(:)
    !write(*,*) atoms%itypat(:)
    !stop
    !t_num(1)=1 !CORRECT_IT
    !t_num(2)=1 !CORRECT_IT
    !write(*,*) t_num
    !stop
    !do i = 1 ,atoms_type 
    !    read(1377,*) t_num(i)
    !end do
    !read(1377,*) !Q's and alpha's of atom types
    !do j = 1 , atoms_type
    !    do i = 1 ,lcn 
    !        read(1377,*) Qt(i,j) , At(i,j)
    !    end do
    !end do
    Qt=parini%qt
    At=parini%at
    !read(1377,*) !Steepest descent step
    !read(1377,*) SD_S_Q, SD_S_A, SD_S_R
    SD_S_Q=parini%alphax_q_fit_elecpot
    SD_S_A=parini%alphax_a_fit_elecpot
    SD_S_R=parini%alphax_r_fit_elecpot
    !read(1377,*) !Minimum error to be reached
    !read(1377,*) err
    err=parini%pot_rmse_tol
    !read(1377,*) !Number of Electrons for Single Atoms
    !read(1377,*) Ne
    !read(1377,*) !to be ploted atom's number
    !read(1377,*) n_at
    n_at=parini%iat_plot
!/////////////////////////////////E.O. READING INPUT PARAMETERS///////////////////////////
    associate(nx=>poisson_dft%ngpx,ny=>poisson_dft%ngpy,nz=>poisson_dft%ngpz)
    associate(cv1=>cv_temp(1,1),cv2=>cv_temp(2,2),cv3=>cv_temp(3,3))
    associate(x_at=>rat(1,n_at),y_at=>rat(2,n_at),z_at=>rat(3,n_at))
    !write(2,*) "Number of grids : " , nx , ny , nz
    !write(2,'(a,3es14.6)') 'Grid spacing : ' , hgx , hgy , hgz
    call yaml_mapping_open('number of grid points',flow=.true.)
    call yaml_map('nx',nx)
    call yaml_map('ny',ny)
    call yaml_map('nz',nz)
    call yaml_mapping_close()
    cel_vol = nx*ny*nz
    cel_vol_inv = 1.d0/(cel_vol+0.d0)
    allocate(dft_pot(nx,ny,nz),dft_rho(nx,ny,nz),weighted_rho(nx,ny,nz),dft_fat(1:3,1:atoms%nat),weight(nx,ny,nz),gwit(1:atoms%nat))
    ss=0
    do j = 1 , atoms_type
        ss = ss + t_num(j-1)
        rat_t(1:3,j) = rat(1:3,ss+1)
    end do
    volinv = (cv1*cv2*cv3)*cel_vol_inv
    !this is because get_psolver_bps changes poisson_dft%rho
    dft_rho=poisson_dft%rho
!///////////////////////////////////////////DFT PART//////////////////////////////////////
    dft_fat = 0.d0
    dft_ener= 0.d0 
    atoms%fat = 0.d0 
    poisson_dft%task_finit=''
    call init_hartree(parini,atoms,poisson_dft,Q(1,1:atoms%nat))
    call get_hartree(parini,poisson_dft,atoms,Q(1,1:atoms%nat),ehartree)
    dft_ener = ehartree
    dft_pot = poisson_dft%pot
    dft_q = hgx*hgy*hgz*sum(dft_rho)
    !write(2,*)'Total charge density of DFT : ',dft_q
    call yaml_map('Total charge of DFT',dft_q,fmt='(f10.3)')
!//////////////////////////////////////////CENT PART//////////////////////////////////////
    allocate(cent_pot(nx,ny,nz))
    allocate(qpar(1:lcn,1:atoms%nat),apar(1:lcn,1:atoms%nat),rpar(1:3,1:atoms%nat),qpar_t(1:lcn,1:atoms_type))
    allocate(cent_fat(1:3,1:atoms%nat))
    !---------------------------------------------------------------------
    do i = 1 , lcn
        do j = 1 , parini%ntypat
            do k = 1 , atoms%nat
                if(trim(atoms%sat(k))==trim(parini%stypat(j))) then
                    a(i,k)=at(i,j)
                    q(i,k)=qt(i,j)
                end if
            end do
        end do
    end do
    qnrm = (dft_q+0.d0)/(sum(Q(1:lcn,1:atoms%nat)))
    Q(1:lcn,1:atoms%nat) = Q(1:lcn,1:atoms%nat)*qnrm
    do i = 1 , lcn
        ss=1
        do j = 1 , atoms%nat
            if(trim(atoms%sat(j))==trim(parini%stypat(ss))) then
                if(ss<=parini%ntypat) then
                    qt(i,ss)=q(i,j)
                    at(i,ss)=a(i,j)
                    ss = ss+1
                endif
            endif
        enddo
    enddo 
    weight = 1.d0
    call set_rcov(atoms)
    if(parini%cutoff_fit_elecpot) then
        call yaml_comment('Cutoff is based on exp(-r^4) and covalent radius.')
        call yaml_mapping_open('rcov',flow=.true.)
        call update_ratp(atoms)
        do i = 1 , atoms%nat
            call yaml_map(trim(atoms%sat(i)),atoms%rcov(i),fmt='(f6.3)')
            do iz = 1 , nz
                dz = (hgz*(iz-1)-atoms%ratp(3,i))**2
                do iy = 1 , ny
                    dy = (hgy*(iy-1)-atoms%ratp(2,i))**2
                    do ix = 1 , nx
                        dx = (hgx*(ix-1)-atoms%ratp(1,i))**2
                        dr = sqrt(dx+dy+dz)
                        weight(ix,iy,iz) = weight(ix,iy,iz)*(1.d0-exp(-1.d0*(dr/(4.d0*atoms%rcov(i)))**4))
                    enddo
                enddo
            enddo
        enddo
        call yaml_mapping_close()
    else
        call yaml_comment('Cutoff : whole cell is considered')
    end if
    do iz = 1 , nz
        do iy = 1 , ny
            do ix = 1 ,nx
                weighted_rho(ix,iy,iz) = dft_rho(ix,iy,iz)*weight(ix,iy,iz)
            enddo
        enddo
    enddo
    allocate(atom_charge_type(1:atoms_type),atom_charge(1:atoms%nat))
    !read(1377,*) ! atoms initial charge
    !do i = 1 , atoms_type
    !    !read(1377,*) atom_charge_type(i)
    !enddo
    do iat=1,atoms%nat
        atom_charge(iat)=atoms%zat(iat)
    enddo
    !ss=0
    !do j = 1 , atoms_type
    !    ss = ss + t_num(j-1)
    !    atom_charge(ss+1:t_num(j)+ss) = atom_charge_type(j)
    !end do
    ng(1) = nx ; ng(2) = ny ; ng(3) = nz
    total_time = 0.d0
    allocate(cent%rel(1:3,1:atoms%nat))
    call yaml_sequence_open('SD iterations')
    allocate(natom_type(parini%ntypat))
    natom_type = 0
    do i = 1 , parini%ntypat
        do j = 1 , atoms%nat
            if(trim(atoms%sat(j))==trim(parini%stypat(i))) natom_type(i)=natom_type(i)+1
        end do
    end do
    allocate(a_type(lcn,parini%ntypat,maxval(natom_type)),q_type(lcn,parini%ntypat,maxval(natom_type)))
    do iter=1,huge(iter_max)
        call yaml_sequence(advance='no')
        call yaml_map('iter',iter)
        do iat=1,atoms%nat
            atoms%qat(iat)=sum(Q(1:lcn,iat))
        enddo
        call calc_multipoles_centt(parini,atoms,poisson_cent,rat)
        call yaml_map('CENT electric dpm',poisson_cent%dpm,fmt='(f10.4)')
        call yaml_map('CENT electric qpm',poisson_cent%qpm(1:3,1:3),fmt='(f10.4)')
        dpm_err(1:3)=poisson_cent%dpm(1:3)-poisson_dft%dpm(1:3)
        qpm_err(1:3,1:3)=poisson_cent%qpm(1:3,1:3)-poisson_dft%qpm(1:3,1:3)
        call yaml_map('dpm_err',dpm_err,fmt='(f10.4)')
        call yaml_map('qpm_err',qpm_err(1:3,1:3),fmt='(f10.4)')
        start = 0.d0
        finish = 0.d0
        call cpu_time(start)
        inquire(file='EXIT',exist=esc)
        if(esc .and. iter >2) exit 
        a_max = maxval(at)
        poisson_cent%rgcut=6.d0*a_max !ask Dr. Ghasemi if it's OK !
!-------------------------------------calculating potential and its derivatives--------------------
        call put_pot_sym_rzx(rat,hgx,hgy,hgz,atoms%nat,Q(1:lcn,1:atoms%nat),A(1:lcn,1:atoms%nat),ng(1:3),lcn,.true.,&
            weight(1:nx,1:ny,1:nz),poisson_dft%pot(1:nx,1:ny,1:nz),cent_pot(1:nx,1:ny,1:nz),qpar(1:lcn,1:atoms%nat),apar(1:lcn,1:atoms%nat),rpar(1:3,1:atoms%nat))
        gtot = sum(qpar)
        qpar = qpar - (gtot+0.d0)/(atoms%nat*lcn+0.d0)
        Q = Q - SD_S_Q*qpar
        A = A - SD_S_A*apar
        rat(1:3,1:atoms%nat) = rat(1:3,1:atoms%nat) - SD_S_R*rpar(1:3,1:atoms%nat)
        call yaml_map('sum of Q',sum(Q),fmt='(es15.7)')
        do i = 1 , lcn
               ss=1
            do j = 1 , atoms%nat
                if(trim(atoms%sat(j))==trim(parini%stypat(ss))) then
                    if(ss<=parini%ntypat) then
                        qt(i,ss)=q(i,j)
                        at(i,ss)=a(i,j)
                        ss = ss+1
                    end if
                endif
            enddo
        enddo 
        do i = 1 , lcn
            do j = 1 , parini%ntypat
                ss=1
                do k = 1 , atoms%nat
                    if(trim(atoms%sat(k))==trim(parini%stypat(j))) then
                        a_type(i,j,ss)=a(i,k)
                        q_type(i,j,ss)=q(i,k)
                        ss=ss+1
                    endif
                enddo
            enddo
        enddo
        do i = 1 , lcn
            do j = 1 , parini%ntypat
                do k = 1 , atoms%nat
                    if(trim(atoms%sat(k))==trim(parini%stypat(j))) a(i,k)=at(i,j)
                end do
            end do
        end do
        call yaml_sequence_open('info for each atom')
        do i = 1 , atoms%nat
            call yaml_sequence(advance='no')
            call yaml_map('iat',i)
            call yaml_map('type',trim(atoms%sat(i)))
            call yaml_map('qpar',qpar(1:lcn,i),fmt='(es15.7)')
            call yaml_map('apar',apar(1:lcn,i),fmt='(es15.7)')
            call yaml_map('rpar',rpar(1:3,i),fmt='(es15.7)')
            call get_rat_iat(atoms,i,dxyz)
            dxyz(1:3)=rat(1:3,i)-dxyz(1:3)
            call yaml_map('dr',dxyz,fmt='(es15.7)')
            call yaml_map('rel',rat(1:3,i),fmt='(es15.7)')
            call yaml_map('Q',Q(1:lcn,i),fmt='(es15.7)')
            call yaml_map('A',A(1:lcn,i),fmt='(es15.7)')
            total_charge = 0.d0
            do l = 1 , lcn
                total_charge = total_charge+Q(l,i)
            end do
            call yaml_map('total charge',total_charge,fmt='(es15.7)')
            call yaml_map('charge changes',-abs(atom_charge(i))+abs(total_charge),fmt='(es15.7)')
        end do
        call yaml_sequence_close()
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
        pot_err= hgx*hgy*hgz*pot_err
        rmse_old=rmse
        rmse = sqrt(pot_err)
        call yaml_mapping_open('SD',flow=.true.)
        call yaml_map('iter',iter)
        call yaml_map('rmse',rmse,fmt='(es14.6)')
        call yaml_map('err_max',err_max,fmt='(es14.6)')
        if (abs(rmse_old-rmse) < err )then
            isatur=isatur+1
        else
            isatur=0
        endif
        call yaml_map('isatur',isatur)
        call yaml_mapping_close()
        call cpu_time(finish)
        total_time = total_time + finish - start
        call yaml_map('time of each SD iter',finish-start,fmt='(f6.3)')
        write(2,'(a45,f6.3,a45)') "================================ ITER TIME : ",finish-start,"(sec)========================================"
        if(isatur>nsatur) then
            call yaml_comment('MAX CONVERSION REACHED',hfill='~')
            exit
        endif
    end do !end of loop over iter
    call yaml_sequence_close()
    call yaml_map('total time of SD',total_time,fmt='(f6.3)')
    call yaml_map('average time of each SD iter',(total_time)/real(iter,kind=8),fmt='(f6.3)')
    call yaml_sequence_open('final values of params')
        do i = 1 , parini%ntypat
        call yaml_sequence(advance='no')
        call yaml_map('type',trim(parini%stypat(i)))
        do j = 1 , lcn    
            call stdval_rzx(q_type(j,i,1:natom_type(i)),natom_type(i),mean,std,var)
            mean_arr(j)=mean
            std_arr(j)=std
        enddo
        call yaml_map('avg of Q',mean_arr(1:lcn),fmt='(es14.6)')
        call yaml_map('std of Q',std_arr(1:lcn),fmt='(es14.6)')
        do j = 1 , lcn    
            call stdval_rzx(a_type(j,i,1:natom_type(i)),natom_type(i),mean,std,var)
            mean_arr(j)=mean
            std_arr(j)=std
        enddo
        call yaml_map('avg of A',mean_arr(1:lcn),fmt='(es14.6)')
        call yaml_map('std of A',std_arr(1:lcn),fmt='(es14.6)')
    enddo
    call yaml_sequence_close()
    do ix = 1 , nx
        write(3,*) ix*hgx , dft_pot(ix,ny/2,nz/2) ,cent_pot(ix,ny/2,nz/2)
    end do 
!/////////////////////////////////E.O. CENT PART//////////////////////////////////////////
    poisson_cent%hgrid = poisson_dft%hgrid
    poisson_cent%ngpx = nx
    poisson_cent%ngpy = ny
    poisson_cent%ngpz = nz
    poisson_dft%rgcut  = 8.d0*maxval(A)
    poisson_cent%rgcut = poisson_dft%rgcut
    nbgx = int(poisson_cent%rgcut/hgx)+2
    nbgy = int(poisson_cent%rgcut/hgy)+2
    nbgz = int(poisson_cent%rgcut/hgz)+2
    poisson_cent%lda = nx
    poisson_dft%lda  = nx
    poisson_cent%bc = atoms%boundcond
    poisson_dft%bc  = atoms%boundcond
    poisson_cent%nat = atoms%nat
    poisson_dft%nat  = atoms%nat
    poisson_cent%cv = atoms%cellvec
    poisson_dft%cv  = atoms%cellvec
    poisson_cent%task_finit='alloc_rho'
    call init_hartree(parini,atoms,poisson_cent,Q(1,1:atoms%nat))
    poisson_cent%rcart(1,1:atoms%nat) = rat(1,1:atoms%nat) - nbgx*hgx     
    poisson_cent%rcart(2,1:atoms%nat) = rat(2,1:atoms%nat) - nbgy*hgy     
    poisson_cent%rcart(3,1:atoms%nat) = rat(3,1:atoms%nat) - nbgz*hgz     
    poisson_dft%rcart(1,1:atoms%nat)  = atoms%ratp(1,1:atoms%nat)  - nbgx*hgx     
    poisson_dft%rcart(2,1:atoms%nat)  = atoms%ratp(2,1:atoms%nat)  - nbgy*hgy     
    poisson_dft%rcart(3,1:atoms%nat)  = atoms%ratp(3,1:atoms%nat)  - nbgz*hgz     
    poisson_dft%rho = dft_rho
    poisson_cent%reset_rho = .True.
    poisson_cent%q = Q(1,1:atoms%nat)
    poisson_cent%gw_ewald = A(1,1:atoms%nat)
    call put_charge_density(parini,poisson_cent)
    do i = 2 , lcn
        poisson_cent%reset_rho = .false.
        poisson_cent%q = Q(i,1:atoms%nat)
        poisson_cent%gw_ewald = A(i,1:atoms%nat)
        call put_charge_density(parini,poisson_cent)
    
    enddo
    allocate(cent_rho(1:nx,1:ny,1:nz))
    cent_rho = poisson_cent%rho
    do ix = 1 , nx
        write(4,'(4es16.7)') ix*hgx , dft_rho(ix,ny/2,nz/2) ,cent_rho(ix,ny/2,nz/2),weighted_rho(ix,ny/2,nz/2)
    end do

    do iat=1,atoms%nat
        atoms%qat(iat)=sum(Q(1:lcn,iat))
    enddo
    call calc_multipoles_centt(parini,atoms,poisson_cent,rat)
    call yaml_map('CENT electric dpm (FINAL)',poisson_cent%dpm,fmt='(f10.4)')
    call yaml_map('CENT electric qpm (FINAL)',poisson_cent%qpm(1:3,1:3),fmt='(f10.4)')
    dpm_err(1:3)=poisson_cent%dpm(1:3)-poisson_dft%dpm(1:3)
    qpm_err(1:3,1:3)=poisson_cent%qpm(1:3,1:3)-poisson_dft%qpm(1:3,1:3)
    call yaml_map('dpm_err (FINAL)',dpm_err,fmt='(f10.4)')
    call yaml_map('qpm_err (FINAL)',qpm_err(1:3,1:3),fmt='(f10.4)')
    poisson_cent%rho=-1.d0*poisson_cent%rho
    call cube_write('cent_rho.cube',atoms,poisson_cent,'rho')
    poisson_cent%rho=-1.d0*poisson_cent%rho
    call get_hartree(parini,poisson_cent,atoms,Q(1,1:atoms%nat),ehartree)
    call get_hartree(parini,poisson_dft,atoms,Q(1,1:atoms%nat),ehartree)
    cent_ener = ehartree
    call yaml_mapping_open('energy',flow=.true.)
    call yaml_map('CENT',cent_ener,fmt='(es15.7)')
    call yaml_map('DFT',dft_ener,fmt='(es15.7)')
    call yaml_map('dE',cent_ener-dft_ener,fmt='(es15.7)')
    call yaml_map('dE percent',(cent_ener-dft_ener)*100.d0/dft_ener,fmt='(es15.7)')
    call yaml_mapping_close()

    poisson_dft%q=atoms%zat
    poisson_cent%q=atoms%zat
    poisson_dft%gw_ewald = 0.5d0
    poisson_cent%gw_ewald=poisson_dft%gw_ewald
    call get_hartree_force(parini,poisson_dft,atoms)
    dft_fat = atoms%fat
    atoms%fat = 0.d0
    call get_hartree_force(parini,poisson_cent,atoms)
    cent_fat = atoms%fat
    call yaml_sequence_open('force due to electronic charge')
    force_err=0.d0
    do i = 1 , atoms%nat
        call yaml_sequence(advance='no')
        total_force_cent = sqrt(sum(cent_fat(1:3,i)**2))
        total_force_dft = sqrt(sum(dft_fat(1:3,i)**2))
        force_err = force_err + sum((cent_fat(1:3,i)-dft_fat(1:3,i))**2)
        call yaml_map('type',trim(atoms%sat(i)))
        call yaml_map('iat',i)
        call yaml_map('dft_fat',dft_fat(1:3,i),fmt='(es15.7)')
        call yaml_map('cent_fat',cent_fat(1:3,i),fmt='(es15.7)')
    end do
    call yaml_sequence_close()
    force_err=sqrt(force_err)
    call yaml_map('force rmse',force_err,fmt='(es15.7)')
    end associate
    end associate
    end associate
    call fini_hartree(parini,atoms,poisson_dft)
end subroutine fit_elecpot
!*****************************************************************************************
subroutine put_pot_sym_rzx(rat,hgx,hgy,hgz,nat,qat,gw,ng,lcn,reset,weight,dft_pot,cent_pot,qpar,apar,rpar)
    implicit none
    logical :: reset
    integer , intent(in):: nat, ng(1:3), lcn
    real(8) , intent(in):: rat(1:3,1:nat), hgx, hgy, hgz, qat(1:lcn,1:nat),gw(1:lcn,1:nat),weight(1:ng(1),1:ng(2),1:ng(3))
    real(8) , intent(in):: dft_pot(1:ng(1),1:ng(2),1:ng(3))
    real(8) , intent(out):: cent_pot(1:ng(1),1:ng(2),1:ng(3)), apar(1:lcn,1:nat), qpar(1:lcn,1:nat), rpar(1:3,1:nat)
    ! local variables 
    integer :: l, iat, iz, iy, ix, nx, ny, nz
    real(8) :: asqinv, pisqrtinv, pi, r, dx, dy, dz
    real(8) :: cent_pot_a_par(1:ng(1),1:ng(2),1:ng(3)), cent_pot_q_par(1:ng(1),1:ng(2),1:ng(3))
    real(8) :: cent_pot_x_par(1:ng(1),1:ng(2),1:ng(3)), cent_pot_y_par(1:ng(1),1:ng(2),1:ng(3)), cent_pot_z_par(1:ng(1),1:ng(2),1:ng(3))
    !-- FDM params ----
!    real(8) :: cent_pot_t(1:ng(1),1:ng(2),1:ng(3)), ttx(1:nat), tty(1:nat), ttz(1:nat) , h
    nx=ng(1); ny=ng(2); nz=ng(3)
    pi = 4.d0*atan(1.d0)
    if (reset) cent_pot = 0.d0
        do l = 1,lcn
            do iat=1,nat
                asqinv=1.d0/(gw(l,iat)**2)
                pisqrtinv = 1.d0/sqrt(pi)
                do iz=1,nz
                    dz=(iz-1)*hgz - rat(3,iat)
                    do iy=1,ny
                        dy=(iy-1)*hgy - rat(2,iat)
                        do ix=1,nx
                            dx=(ix-1)*hgx - rat(1,iat)
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
        rpar(1:3,1:nat) = 0.d0
        do l = 1,lcn
            do iat=1,nat
                asqinv=1.d0/(gw(l,iat)**2)
                pisqrtinv = 1.d0/sqrt(pi)
                do iz=1,nz
                    dz= rat(3,iat)- (iz-1)*hgz
                    do iy=1,ny
                        dy=rat(2,iat) - (iy-1)*hgy 
                        do ix=1,nx
                            dx=rat(1,iat)-(ix-1)*hgx
                            r=sqrt(dx**2+dy**2+dz**2)
                            if (r <= 1.d-6) then 
                                cent_pot_a_par(ix,iy,iz) =-2.d0*qat(l,iat)*asqinv*pisqrtinv
                                cent_pot_q_par(ix,iy,iz) =2.d0*pisqrtinv/gw(l,iat)
                                cent_pot_x_par(ix,iy,iz) = 0.d0
                                cent_pot_y_par(ix,iy,iz) = 0.d0
                                cent_pot_z_par(ix,iy,iz) = 0.d0
                            else
                                cent_pot_a_par(ix,iy,iz)=-2.d0*qat(l,iat)*exp(-r**2*asqinv)*asqinv*pisqrtinv
                                cent_pot_q_par(ix,iy,iz) =erf(r/gw(l,iat))/r
                                cent_pot_x_par(ix,iy,iz) =(2*qat(l,iat)*dx)/(gw(l,iat)*exp((r/gw(l,iat))**2)*Sqrt(Pi)*(r**2))&
                                    -(qat(l,iat)*dx*Erf(Sqrt(r**2)/gw(l,iat)))/(r**2)**1.5
                                cent_pot_y_par(ix,iy,iz) =(2*qat(l,iat)*dy)/(gw(l,iat)*exp((r/gw(l,iat))**2)*Sqrt(Pi)*(r**2))&
                                    -(qat(l,iat)*dy*Erf(Sqrt(r**2)/gw(l,iat)))/(r**2)**1.5
                                cent_pot_z_par(ix,iy,iz) =(2*qat(l,iat)*dz)/(gw(l,iat)*exp((r/gw(l,iat))**2)*Sqrt(Pi)*(r**2))&
                                    -(qat(l,iat)*dz*Erf(Sqrt(r**2)/gw(l,iat)))/(r**2)**1.5
                                
                            endif
                        enddo ! of ix
                    enddo ! of iy
                enddo ! of iz
                qpar(l,iat)=2.d0*hgx*hgy*hgz*sum(cent_pot_q_par*weight*(cent_pot-dft_pot))
                apar(l,iat)=2.d0*hgx*hgy*hgz*sum(cent_pot_a_par*weight*(cent_pot-dft_pot))
                rpar(1,iat)=rpar(1,iat)+2.d0*hgx*hgy*hgz*sum(cent_pot_x_par*weight*(cent_pot-dft_pot))
                rpar(2,iat)=rpar(2,iat)+2.d0*hgx*hgy*hgz*sum(cent_pot_y_par*weight*(cent_pot-dft_pot))
                rpar(3,iat)=rpar(3,iat)+2.d0*hgx*hgy*hgz*sum(cent_pot_z_par*weight*(cent_pot-dft_pot))
            enddo ! of iat
        enddo !of lcn
!        !----------------------------FDM PART------------------------------
!        h = 1.d-7
!        do iat = 1 , nat
!            write(99,'(a14,i3,3es16.7)') 'rpar_ana : ', iat, rpar(:,iat) 
!        end do
!        cent_pot_t = 0.d0
!        do iat=1,nat
!            pisqrtinv = 1.d0/sqrt(pi)
!            do iz=1,nz
!                dz=(iz-1)*hgz - rat(3,iat)
!                do iy=1,ny
!                    dy=(iy-1)*hgy - rat(2,iat)
!                    do ix=1,nx
!                        dx=(ix-1)*hgx - rat(1,iat)
!                        r=sqrt(dx**2+dy**2+dz**2)
!                        if (r <= 1.d-6) then 
!                            cent_pot_t(ix,iy,iz)=cent_pot_t(ix,iy,iz)+2.d0*qat(1,iat)*pisqrtinv/gw(1,iat)
!                        else
!                            cent_pot_t(ix,iy,iz)=cent_pot_t(ix,iy,iz)+qat(1,iat)*erf(r/gw(1,iat))/r
!                        endif
!                    enddo ! of ix
!                enddo ! of iy
!            enddo ! of izi
!            ttx(iat)=hgx*hgy*hgz*sum((weight*cent_pot_t-dft_pot)**2)
!            tty(iat)=hgx*hgy*hgz*sum((weight*cent_pot_t-dft_pot)**2)
!            ttz(iat)=hgx*hgy*hgz*sum((weight*cent_pot_t-dft_pot)**2)
!        end do
!        do iat = 1, nat
!            do iz=1,nz
!                dz=rat(3,iat)- (iz-1)*hgz
!                do iy=1,ny
!                    dy=rat(2,iat)- (iy-1)*hgy
!                    do ix=1,nx
!                        dx=rat(1,iat) - (ix-1)*hgx + h
!                        r=sqrt(dx**2+dy**2+dz**2)
!                        if (r <= 1.d-6) then 
!                            cent_pot_t(ix,iy,iz)=2.d0*qat(1,iat)*pisqrtinv/gw(1,iat)
!                        else
!                            cent_pot_t(ix,iy,iz)=qat(1,iat)*erf(r/gw(1,iat))/r
!                        endif
!                    enddo ! of ix
!                enddo ! of iy
!            enddo ! of izi
!            ttx(iat)=(-1.d0*ttx(iat)+(hgx*hgy*hgz*sum((weight*cent_pot_t-dft_pot)**2)))/h
!            do iz=1,nz
!                dz=rat(3,iat)- (iz-1)*hgz
!                do iy=1,ny
!                    dy=rat(2,iat) - (iy-1)*hgy + h
!                    do ix=1,nx
!                        dx=rat(1,iat)- (ix-1)*hgx
!                        r=sqrt(dx**2+dy**2+dz**2)
!                        if (r <= 1.d-6) then 
!                            cent_pot_t(ix,iy,iz)=2.d0*qat(1,iat)*pisqrtinv/gw(1,iat)
!                        else
!                            cent_pot_t(ix,iy,iz)=qat(1,iat)*erf(r/gw(1,iat))/r
!                        endif
!                    enddo ! of ix
!                enddo ! of iy
!            enddo ! of izi
!            tty(iat)=(-1.d0*tty(iat)+(hgx*hgy*hgz*sum((weight*cent_pot_t-dft_pot)**2)))/h
!            do iz=1,nz
!                dz=rat(3,iat)- (iz-1)*hgz +h
!                do iy=1,ny
!                    dy=rat(2,iat)- (iy-1)*hgy
!                    do ix=1,nx
!                        dx=rat(1,iat) - (ix-1)*hgx
!                        r=sqrt(dx**2+dy**2+dz**2)
!                        if (r <= 1.d-6) then 
!                            cent_pot_t(ix,iy,iz)=2.d0*qat(1,iat)*pisqrtinv/gw(1,iat)
!                        else
!                            cent_pot_t(ix,iy,iz)=qat(1,iat)*erf(r/gw(1,iat))/r
!                        endif
!                    enddo ! of ix
!                enddo ! of iy
!            enddo ! of izi
!            ttz(iat)=(-1.d0*ttz(iat)+(hgx*hgy*hgz*sum((weight*cent_pot_t-dft_pot)**2)))/h
!        enddo ! of iat
!        do iat = 1 , nat
!            write(99,'(a14,i3,3es16.7)') 'rpar_fdm : ', iat, ttx(iat),tty(iat),ttz(iat) 
!            write(99,'(a65)') '-----------------------------------------------------------------' 
!        end do
!        !---------------------------- E.O. FDM PART---------------------------------------------------------------------------
end subroutine put_pot_sym_rzx
!*****************************************************************************************
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
!*****************************************************************************************
