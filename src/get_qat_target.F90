!*****************************************************************************************
module mod_qat_target
contains
!*****************************************************************************************
subroutine get_qat_target(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_copy_old, update_ratp, update_rat
    use mod_electrostatics, only: typ_poisson
    use mod_ann, only: typ_ann_arr
    use yaml_output
    implicit none
    type(typ_parini) :: parini
    type(typ_atoms)  :: atoms
    type(typ_atoms_arr), target:: atoms_arr
    type(typ_ann_arr):: ann_arr
    type(typ_poisson):: poisson
    !local variables
    real(8), allocatable :: tmp_rat(:,:),tmp_cv(:,:)
    real(8), allocatable :: gausswidth(:)
    real(8), allocatable :: rrad(:),trial_rho(:)
    real(8), allocatable :: rho_Mg(:),pot_Mg(:),pot_scn_Mg(:)
    real(8), allocatable :: rho_O(:),pot_O(:),pot_scn_O(:)
    real(8), allocatable :: grid_rho(:,:,:)
    real(8), allocatable :: E_par(:,:)!, E_par_0(:,:)
    real(8), allocatable :: trial_energy(:), rhs(:), atoms_qat(:), qat_target(:), ipiv(:)
    real(8), allocatable :: amat(:,:),eigen_val_amat(:,:)
    real(8), allocatable :: real_eigenval(:),imag_eigenval(:),work(:)
    real(8), allocatable :: right_eigenvec(:,:),left_eigenvec(:,:)
    real(8) :: den_coeff_Mg, den_coeff_O
    real(8) :: Mg_qat, O_qat, trial_den_coeff, hgp, cost_function, energy_error
    real(8) :: dx ,dy ,dz ,dr
    real(8) :: pi, vac_space
    real(8) :: time_s, time_f
    real(8) :: Q_Mg_min , Q_Mg_max , Q_Mg_mean, Q_O_min, Q_O_max, Q_O_mean,q_tot 
    real(8) :: tt
    integer :: iat,jat,kat,trial_num
    integer :: igp, igx, igy, igz
    integer :: ngp, ngpx, ngpy, ngpz
    integer :: linearGridNumber
    integer :: info
    !INITIAL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pi = 4.d0*atan(1.d0)
    call read_data_yaml(parini,'list_posinp.yaml',atoms_arr)
    call atom_copy_old(atoms_arr%atoms(1),atoms,'atoms_arr%atoms(iconf)->atoms')
    allocate(tmp_rat(3,atoms%nat),tmp_cv(3,3))
    allocate(qat_target(1:atoms%nat))
    allocate(gausswidth(atoms%nat))
    Mg_qat = parini%q_avg_target
    O_qat =  8.d0-parini%q_avg_target
    do iat=1,atoms%nat
        if(trim(atoms%sat(iat))=='Mg') then
            qat_target(iat)=Mg_qat
            gausswidth(iat)=parini%gaussian_width_Mg
        elseif(trim(atoms%sat(iat))=='O') then
            qat_target(iat)=O_qat
            gausswidth(iat)=parini%gaussian_width_O
        end if
    end do
    q_tot=sum(qat_target(1:atoms%nat))
    call update_ratp(atoms)
    tmp_rat(1:3,1:atoms%nat)=atoms%ratp(1:3,1:atoms%nat)
    tmp_cv(1:3,1:3)=atoms%cellvec(1:3,1:3)
    vac_space = parini%free_space 
    call cal_min_cv(atoms,vac_space)
    !BPS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    poisson%cal_scn=parini%cal_scn
    poisson%screening_factor=parini%screening_factor
    poisson%alpha=parini%alpha_ewald
    poisson%task_finit="alloc_rho:set_ngp"
    call init_hartree(parini,atoms,poisson,gausswidth)
    !SINGLE ATOM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    call cpu_time(time_s)
    ngp = 10**5
    allocate(pot_Mg(0:ngp),pot_O(0:ngp),trial_rho(0:ngp))
    call get_scf_pot(atoms%cellvec,ngp,poisson%rgcut,parini%gaussian_width_Mg,parini%screening_factor,pot_Mg,.false.)
    call get_scf_pot(atoms%cellvec,ngp,poisson%rgcut,parini%gaussian_width_O,parini%screening_factor,pot_O,.true.,0.5d0,trial_rho)
    !LOCAL POT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    trial_num = atoms%ntrial ! Equal to trial gaussian numbers
    allocate(amat(1:atoms%nat+1,1:atoms%nat+1))
    allocate(E_par(1:trial_num,1:atoms%nat))
    allocate(rhs(1:atoms%nat))
    call get_amat_cent2_trial(parini,atoms,poisson,pot_Mg,pot_O,trial_rho,qat_target,ngp,E_par,rhs,amat)
    !Solve Linear Equations and get EigenValues%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    allocate(ipiv(1:atoms%nat))
    allocate(atoms_qat(1:atoms%nat+1))
    allocate(eigen_val_amat(1:atoms%nat,1:atoms%nat))
    eigen_val_amat(1:atoms%nat,1:atoms%nat)=0.d0
    do iat=1,atoms%nat
        do jat=iat,atoms%nat
           eigen_val_amat(iat,jat)=amat(iat,jat) 
        end do
    end do
    call DGETRF(atoms%nat+1,atoms%nat+1,amat,atoms%nat+1,ipiv,info)
    if(info/=0) then
        write(*,'(a19,i8)') 'ERROR: DGETRF info=',info
        stop
    endif
    atoms_qat(1:atoms%nat)=rhs(1:atoms%nat)
    atoms_qat(atoms%nat+1)=q_tot
    call DGETRS('N',atoms%nat+1,1,amat,atoms%nat+1,ipiv,atoms_qat,atoms%nat+1,info)
    if(info/=0) then
        write(*,'(a19,i8)') 'ERROR: DGETRS info=',info
        stop
    endif
    !EigenValue%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    allocate(real_eigenval(1:atoms%nat),work(1:4*atoms%nat))
    !allocate(imag_eigenval(1:atoms%nat))
    !allocate(right_eigenvec(1:atoms%nat,1:atoms%nat),left_eigenvec(1:atoms%nat,1:atoms%nat))
    !call DGEEV('N','N',atoms%nat,eigen_val_amat,atoms%nat,real_eigenval,imag_eigenval,&
    !            left_eigenvec,1,right_eigenvec,1,work,4*atoms%nat,info) 
    call DSYEV('N','U',atoms%nat,eigen_val_amat,atoms%nat,real_eigenval,work,4*atoms%nat,info)
    !COST FUNCTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    allocate(trial_energy(1:trial_num))
    cost_function = 0.d0
    energy_error = 0.d0
    do kat=1,trial_num
        trial_energy(kat)=0.d0
        do iat=1,atoms%nat
            dx=atoms%ratp(1,iat)-atoms%ratp(1,atoms%trial_ref_nat(kat))-atoms%trial_ref_disp(1,kat)
            dy=atoms%ratp(2,iat)-atoms%ratp(2,atoms%trial_ref_nat(kat))-atoms%trial_ref_disp(2,kat)
            dz=atoms%ratp(3,iat)-atoms%ratp(3,atoms%trial_ref_nat(kat))-atoms%trial_ref_disp(3,kat)
            trial_energy(kat)=trial_energy(kat)+atoms_qat(iat)*E_par(kat,iat)
            write(99,'(3i3,45es11.3)') kat,atoms%trial_ref_nat(kat),iat,E_par(kat,iat), &!E_par_0(atoms%trial_ref_nat(kat),iat), &
                                     atoms_qat(iat),atoms_qat(atoms%trial_ref_nat(kat)),sqrt(dx**2+dy**2+dz**2)
        end do
        write(55,*)kat,atoms%trial_ref_nat(kat),trial_energy(kat),atoms%trial_ref_energy(kat)
        cost_function=cost_function + (trial_energy(kat)-atoms%trial_ref_energy(kat))**2 +&
                      parini%pen_coeff*((atoms_qat(atoms%trial_ref_nat(kat))-qat_target(atoms%trial_ref_nat(kat)))**2)
        energy_error=energy_error+(trial_energy(kat)-atoms%trial_ref_energy(kat))**2
    end do
    energy_error = sqrt(energy_error/atoms%ntrial)
    call cpu_time(time_f)
    !CHARGE ANALYSIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Q_Mg_min = 1.d100
    Q_Mg_max = -1.d100
    Q_Mg_mean = 0.d0
    Q_O_min = 1.d100
    Q_O_max = -1.d100
    Q_O_mean = 0.d0
    do iat=1,atoms%nat
        if(trim(atoms%sat(iat))=='Mg') then
            if(atoms_qat(iat)>Q_Mg_max) then
                Q_Mg_max = atoms_qat(iat) 
            endif
            if(atoms_qat(iat)<Q_Mg_min) then
                Q_Mg_min = atoms_qat(iat) 
            endif
            Q_Mg_mean=Q_Mg_mean+atoms_qat(iat)
        elseif(trim(atoms%sat(iat))=='O') then
            if(atoms_qat(iat)>Q_O_max) then
                Q_O_max = atoms_qat(iat) 
            endif
            if(atoms_qat(iat)<Q_O_min) then
                Q_O_min = atoms_qat(iat) 
            endif
            Q_O_mean=Q_O_mean+atoms_qat(iat)
        end if
    end do
    Q_Mg_Mean = 2.d0*Q_Mg_Mean/atoms%nat
    Q_O_Mean = 2.d0*Q_O_Mean/atoms%nat

    !PRINTING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    write(88,'(a5,es8.1,a4,es8.1,a3,es10.3,a3,2(es8.1,a3),2(es10.3,a3),2(es8.1,a3),2(es10.3,a3),es10.3,a3,es8.1,a3,es10.3,a2,i3,a2)')&
                                            " |GMg",parini%gaussian_width_Mg," |GO",parini%gaussian_width_O," | ",&
                                            energy_error," | ",&
                                            Q_Mg_min," | ",Q_Mg_max," | ",Q_Mg_max-Q_Mg_min," | ",Q_Mg_mean," | ",&
                                            Q_O_min," | ",Q_O_max," | ",Q_O_max-Q_O_min," | ",Q_O_mean," | ",&
                                            minval(real_eigenval(:))," | ",parini%q_avg_target," | ",&
                                            maxval(abs(trial_energy(:)-atoms%trial_ref_energy))," |",&
                                            atoms%trial_ref_nat(maxloc(abs(trial_energy(:)-atoms%trial_ref_energy)))," |"

    write(77,'(a,es11.3)') 'Penalty Coeffcient:',parini%pen_coeff 
    write(77,'(a,es11.3)') 'Cost Function:     ',cost_function
    write(77,'(a,es11.3)') 'Energy Error(Ha):  ',energy_error
    write(77,'(a,es11.3)') 'Max Energy Error:  ',maxval(abs(trial_energy(:)-atoms%trial_ref_energy))
    write(77,'(a,es11.3)') 'Condition Number:  ' ,maxval(real_eigenval(:))/minval(real_eigenval(:))
    write(77,'(a,es11.3)') 'Minimum EigenVal:  ' ,minval(real_eigenval(:))
    write(77,'(a,es11.3)') 'Maximum EigenVal:  ' ,maxval(real_eigenval(:))
    write(77,'(a,es11.3)') 'Calculation Time:  ', time_f-time_s
    write(77,'(a,es11.3)') 'Total Charge:      ', sum(atoms_qat(1:atoms%nat))
    write(77,'(a,es11.3)') 'Target Mg Charge:  ', parini%q_avg_target
    write(77,'(a104)')"#-------------------------------------------------------------------------------------------------------"
    write(77,'(a104)')"# Q_Mg_min  |  Q_Mg_max  |  Q_Mg_var  |  Q_Mg_mean |  Q_O_min   |   Q_O_max  |  Q_O_var   |  Q_O_mean  |"
    write(77,'(a104)')"#-------------------------------------------------------------------------------------------------------"
    write(77,'(es11.3,a2,es11.3,a2,es11.3,a2,es11.3,a2,es11.3,a2,es11.3,a2,es11.3,a2,es11.3,a2)') &
                Q_Mg_min," |",Q_Mg_max," |",Q_Mg_max-Q_Mg_min," |",Q_Mg_mean," |",&
                Q_O_min," |",Q_O_max," |",Q_O_max-Q_O_min," |",Q_O_mean," |"
    write(77,'(a104)')"#-------------------------------------------------------------------------------------------------------"
    write(77,'(a39 )') "#--------------------------------------"
    write(77,'(a39 )') "#iat | sat |     Qat     | Qat Target |"
    write(77,'(a39 )') "#--------------------------------------"
    do iat = 1 , atoms%nat
        write(77,'(i4,3a3,es11.3,a2,es11.3,a2)')iat," | ",atoms%sat(iat)," | ",atoms_qat(iat),&
                                                " |",qat_target(iat)," |"
    end do
    !Cube Write%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !poisson%reset_rho=.true.
    !poisson%nat=atoms%nat
    !poisson%cv=atoms%cellvec
    !poisson%bc=atoms%boundcond
    !poisson%q(1:poisson%nat)=atoms_qat(1:atoms%nat)
    !poisson%gw(1:poisson%nat)=gausswidth(1:atoms%nat)
    !!call update_ratp(atoms)
    !poisson%rcart(1:3,1:poisson%nat)=atoms%ratp(1:3,1:atoms%nat)

    !write(*,*) poisson%reset_rho,poisson%nat,poisson%bc
    !do iat = 1,3
    !    write(*,*) 'CV:',poisson%cv(iat,1:3)
    !end do
    !do iat = 1, atoms%nat
    !    write(*,'(5es14.6)')poisson%q(iat),poisson%gw(iat),poisson%rcart(1:3,iat)
    !end do
    !call put_charge_density(parini,poisson)
    !CHECKING RHO ON GRID%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !allocate(grid_rho(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    !do iat = 1 , atoms%nat
    !    do igx = 1 , poisson%ngpx
    !        do igy = 1 , poisson%ngpy
    !            do igz = 1 , poisson%ngpz
    !                dx = (-nbgx+igx)*poisson%hgrid(1,1)-atoms%ratp(1,iat)
    !                dy = (-nbgy+igy)*poisson%hgrid(2,2)-atoms%ratp(2,iat)
    !                dz = (-nbgz+igz)*poisson%hgrid(3,3)-atoms%ratp(3,iat)
    !                dr = sqrt(dx**2+dy**2+dz**2)
    !                linearGridNumber=floor(dr/hgp)
    !                grid_rho(igx,igy,igz)=grid_rho(igx,igy,igz)+atoms_qat(iat)*((dr/hgp-linearGridNumber)*(rho(linearGridNumber+1)&
    !                                      -rho(linearGridNumber))+rho(linearGridNumber))
    !            end do
    !        end do
    !    end do
    !end do
    !do igx = 1 , poisson%ngpx
    !    do igy = 1 , poisson%ngpy
    !        do igz = 1 , poisson%ngpz
    !            if(abs(grid_rho(igx,igy,igz))>1.d-6.or. abs(poisson%rho(igx,igy,igz))>1.d-6) then
    !                write(66,'(i,i,i,es14.6)') igx,igy,igz,grid_rho(igx,igy,igz)
    !                write(67,'(i,i,i,es14.6)') igx,igy,igz,poisson%rho(igx,igy,igz)
    !            endif
    !        end do
    !    end do
    !end do
    !CUBE WRITE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !call cube_write('rho.cube',atoms,poisson,'rho')
    !FINALIZE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    call fini_hartree(parini,atoms,poisson)
    atoms%ratp(1:3,1:atoms%nat)=tmp_rat(1:3,1:atoms%nat)
    atoms%cellvec(1:3,1:3)=tmp_cv(1:3,1:3)
    call update_rat(atoms)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end subroutine get_qat_target
!*****************************************************************************************
subroutine cal_screened_poisson(ngp,rrad,den,sf,pot)
    !Calculates Hartree potential.
    implicit none
    integer, intent(in):: ngp
    real(8), intent(in):: rrad(0:ngp), den(0:ngp), sf
    real(8), intent(out):: pot(0:ngp)
    !local variables
    integer:: igp
    real(8):: pi, tt, tt1, tt2
    real(8), allocatable:: wa1(:), wa2(:)
    pi=4.d0*atan(1.d0)
    allocate(wa1(0:ngp),wa2(0:ngp))
    !calculating the first integral: \int_0^r sinh(sf*r') rho(r') r' dr'
    wa1(0)=.0d0
    do igp=1,ngp
        !tt1=exp(-sf*(rrad(igp)-rrad(igp-1)))-exp(-sf*(rrad(igp)+rrad(igp-1)))
        !tt2=exp(-sf*(rrad(igp)-rrad(igp-0)))-exp(-sf*(rrad(igp)+rrad(igp-0)))
        tt1=sinh(sf*rrad(igp-1))
        tt2=sinh(sf*rrad(igp-0))
        !write(42,'(4es20.10)') rrad(igp),tt1,tt2,den(igp)
        tt=den(igp-1)*tt1*rrad(igp-1)+den(igp)*tt2*rrad(igp)
        wa1(igp)=wa1(igp-1)+tt*(rrad(igp)-rrad(igp-1))*0.5d0
    enddo
    !calculating the second integral: \int_r^inf exp(-sf*r') rho(r') r' dr'
    wa2(ngp)=.0d0
    do igp=ngp-1,0,-1
        !tt1=exp(-sf*(rrad(igp+0)-rrad(igp)))-exp(-sf*(rrad(igp)+rrad(igp+0)))
        !tt2=exp(-sf*(rrad(igp+1)-rrad(igp)))-exp(-sf*(rrad(igp)+rrad(igp+1)))
        tt1=exp(-sf*rrad(igp+0))
        tt2=exp(-sf*rrad(igp+1))
        tt=den(igp)*tt1*rrad(igp)+den(igp+1)*tt2*rrad(igp+1)
        wa2(igp)=wa2(igp+1)+tt*(rrad(igp+1)-rrad(igp))*0.5d0
    enddo
    !finally summing the two integrals
    !note that Hartree potential at origin is due zero not infinity,
    !this can be obtained after resolving 0/0 ambiguity
    pot(0)=(4.d0*pi)*wa2(0)
    do igp=1,ngp
        pot(igp)=4.d0*pi*(wa1(igp)*exp(-sf*rrad(igp))+wa2(igp)*sinh(sf*rrad(igp)))/(sf*rrad(igp))
        !write(41,'(4es20.10)') rrad(igp),wa1(igp),wa2(igp),wa1(igp)*exp(-sf*rrad(igp))
    enddo
    deallocate(wa1,wa2)
end subroutine cal_screened_poisson
!*****************************************************************************************
subroutine cal_gauss_screened_poisson_gaussian(ngp,rrad,q,gw,sf,pot)
    !Calculates Hartree potential.
    implicit none
    integer, intent(in):: ngp
    real(8), intent(in):: rrad(0:ngp), q, gw, sf
    real(8), intent(out):: pot(0:ngp)
    !local variables
    integer:: igp
    real(8):: pi, tt1, tt2, tt3, tt4, tt5, ss1, ss2
    pi=4.d0*atan(1.d0)
    do igp=1,ngp
        if(rrad(igp)>10.d0*gw) then
            ss1=0.d0
        else
            ss1=exp(-(rrad(igp)/gw)**2)
        endif
        if(sf*rrad(igp)**2/(1.d0+gw**2*sf)>100.d0) then !100=10**2
            ss2=0.d0
        else
            ss2=exp(-sf*rrad(igp)**2/(1.d0+gw**2*sf))
        endif
        tt1=ss1*erf(2.d0*rrad(igp)*sqrt(sf))
        tt2=gw*ss2*sqrt(sf/(1.d0+gw**2*sf))
        tt3=erf(rrad(igp)/(gw*Sqrt(1.d0+gw**2*sf)))
        tt4=erf((rrad(igp)+2.d0*gw**2*rrad(igp)*sf)/(gw*sqrt(1.d0+gw**2*sf)))
        tt5=2.d0*gw*rrad(igp)*sqrt(sf)
        pot(igp)=q*(-tt1+tt2*(tt3+tt4))/tt5 
        tt1=ss1*erf(2.d0*rrad(igp)*sqrt(sf))
        tt2=gw*ss2*sqrt(sf/(1.d0+gw**2*sf))
        tt3=erf(rrad(igp)/(gw*sqrt(1.d0+gw**2*sf)))
        tt4=erf((rrad(igp)+2.d0*gw**2*rrad(igp)*sf)/(gw*sqrt(1.d0+gw**2*sf)))
        tt5=2.d0*gw*rrad(igp)*sqrt(sf)
        pot(igp)=pot(igp)+q*(tt1+tt2*(tt3-tt4))/tt5
    enddo
    pot(0)=q*2.d0/(sqrt(pi)*gw*(1.d0+gw**2*sf))
end subroutine cal_gauss_screened_poisson_gaussian
!*****************************************************************************************
subroutine cal_powern_screened_poisson_gaussian(ngp,rrad,weight,q,gw,sf,npow,pot)
    !Calculates Hartree potential.
    implicit none
    integer, intent(in):: ngp, npow
    real(8), intent(in):: rrad(0:ngp), weight(0:ngp), q, gw, sf
    real(8), intent(out):: pot(0:ngp)
    !local variables
    integer:: igp, jgp, mgp
    real(8):: pi, tt0, tt1, tt2, ss1, ss2, res, slimit
    real(8), allocatable:: w1(:)
    pi=4.d0*atan(1.d0)
    slimit=-log(1.d-10) ! 45.d0 !exp(-slimit)~2.9E-20
    !write(*,*) 'slimit= ',slimit
    allocate(w1(0:ngp))
    w1(0:ngp)=0.d0
    do igp=0,ngp
        ss1=(sf*rrad(igp))**npow
        if(ss1>slimit) then
            mgp=igp
            exit
        endif
        w1(igp)=exp(-ss1)
    enddo
    !write(*,*) 'mgp,ngp= ',mgp,ngp
    tt0=q/(sqrt(pi)*gw)
    do igp=1,ngp
        res=0.d0
        do jgp=0,mgp
            ss1=(rrad(jgp)-rrad(igp))**2/gw**2
            if(ss1>slimit) then
                tt1=0.d0
            else
                tt1=exp(-ss1)
            endif
            ss2=4.d0*rrad(jgp)*rrad(igp)/gw**2
            if(ss2>slimit) then
                tt2=1.d0
            else
                tt2=(1.d0-exp(-ss2))
            endif
            res=res+w1(jgp)*tt1*tt2*weight(jgp)
        enddo
        pot(igp)=res*tt0/rrad(igp)
    enddo
    res=0.d0
    do jgp=0,mgp
        ss1=rrad(jgp)**2/gw**2
        if(ss1>slimit) then
            tt1=0.d0
        else
            tt1=exp(-ss1)
        endif
        res=res+w1(jgp)*tt1*rrad(jgp)*weight(jgp)
    enddo
    pot(0)=res*4.d0*q/(sqrt(pi)*gw**3)
    deallocate(w1)
end subroutine cal_powern_screened_poisson_gaussian
!*****************************************************************************************
!subroutine cal_gauss_screened_poisson_gaussian(ngp,rrad,q,gw,sf,pot)
!    !Calculates Hartree potential.
!    implicit none
!    integer, intent(in):: ngp
!    real(8), intent(in):: rrad(0:ngp), q, gw, sf
!    real(8), intent(out):: pot(0:ngp)
!    !local variables
!    integer:: igp
!    real(8):: pi, tt1, tt2, tt3, tt4, tt5
!    pi=4.d0*atan(1.d0)
!    pot=0.d0
!    do igp=1,ngp
!        tt1=erf(2.d0*rrad(igp)*sqrt(sf))
!        tt2=gw*exp(rrad(igp)**2/(gw**2+gw**4*sf))*sqrt(sf/(1.d0+gw**2*sf))
!        tt3=erf(rrad(igp)/(gw*Sqrt(1.d0+gw**2*sf)))
!        tt4=erf((rrad(igp)+2.d0*gw**2*rrad(igp)*sf)/(gw*sqrt(1.d0+gw**2*sf)))
!        tt5=2.d0*gw*exp(rrad(igp)**2/gw**2)*rrad(igp)*sqrt(sf)
!        pot(igp)=q*(-tt1+tt2*(tt3+tt4))/tt5 
!        if(abs(tt2)>10.d20 .or. abs(tt5)>10.d20) exit
!        tt1=exp((rrad(igp)**2*sf)/(1.d0+gw**2*sf))*erf(2.d0*rrad(igp)*sqrt(sf))
!        tt2=gw*exp(rrad(igp)**2/gw**2)*sqrt(sf/(1.d0+gw**2*sf))
!        tt3=erf(rrad(igp)/(gw*sqrt(1.d0+gw**2*sf)))
!        tt4=erf((rrad(igp)+2.d0*gw**2*rrad(igp)*sf)/(gw*sqrt(1.d0+gw**2*sf)))
!        tt5=2.d0*gw*exp(rrad(igp)**2*(gw**(-2)+sf/(1.d0+gw**2*sf)))*rrad(igp)*sqrt(sf)
!        pot(igp)=pot(igp)+q*(tt1+tt2*(tt3-tt4))/tt5
!        if(abs(tt2)>10.d20 .or. abs(tt5)>10.d20) exit
!    enddo
!    pot(0)=q*2.d0/(sqrt(pi)*gw*(1.d0+gw**2*sf))
!end subroutine cal_gauss_screened_poisson_gaussian
!*****************************************************************************************
!This routine calculates:
!4*pi/r*\int_0^r rho(r') r'^2 dr' + 4*pi*\int_r^inf rho(r') r' dr'
subroutine cal_pot_hartree(ngp,rrad,den,pot_hartree)
    !Calculates Hartree potential.
    implicit none
    integer, intent(in):: ngp
    real(8), intent(in):: rrad(0:ngp), den(0:ngp)
    real(8), intent(out):: pot_hartree(0:ngp)
    !local variables
    integer:: igp
    real(8):: pi, tt
    real(8), allocatable:: wa1(:), wa2(:)
    pi=4.d0*atan(1.d0)
    allocate(wa1(0:ngp),wa2(0:ngp))
    !calculating the first integral: \int_0^r rho(r') r'^2 dr'
    wa1(0)=.0d0
    do igp=1,ngp
        !wa1(igp)=wa1(igp-1)+den(igp)*rrad(igp)*rrad(igp)*weight(igp)
        tt=den(igp-1)*rrad(igp-1)**2+den(igp)*rrad(igp)**2
        wa1(igp)=wa1(igp-1)+tt*(rrad(igp)-rrad(igp-1))*0.5d0
    enddo
    !calculating the second integral: \int_r^inf rho(r') r' dr'
    wa2(ngp)=.0d0
    do igp=ngp-1,0,-1
        tt=den(igp)*rrad(igp)+den(igp+1)*rrad(igp+1)
        wa2(igp)=wa2(igp+1)+tt*(rrad(igp+1)-rrad(igp))*0.5d0
    enddo
    !finally summing the two integrals
    !note that Hartree potential at origin is due zero not infinity,
    !this can be obtained after resolving 0/0 ambiguity
    pot_hartree(0)=(4.d0*pi)*wa2(0)
    do igp=1,ngp
        pot_hartree(igp)=4.d0*pi*(wa1(igp)/rrad(igp)+wa2(igp))
    enddo
    deallocate(wa1,wa2)
end subroutine cal_pot_hartree
!*****************************************************************************************
subroutine get_scf_pot(cv,ngp,rgcut,gw,scf,pot,calc_trial,gw_trial,trial_rho)
    logical, intent(in) :: calc_trial
    integer, intent(in) :: ngp
    real(8), intent(in) :: cv(1:3,1:3)
    real(8), intent(in) :: rgcut,gw,scf
    real(8), intent(out):: pot(0:ngp)
    real(8), optional,intent(in) :: gw_trial
    real(8), optional,intent(out):: trial_rho(0:ngp)
    !Local Variables!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: igp
    real(8) :: den_coeff, trial_den_coeff, hgp
    real(8) :: rrad(0:ngp), rho(0:ngp), pot_scn(0:ngp)
    real(8) :: pi
    pi = 4.d0*atan(1.d0)
    den_coeff = 1.d0/((pi**1.5d0)*(gw**3))
    !write(*,*) 'Ehsan',present(gw_trial)
    if(calc_trial) then
        trial_den_coeff = 1.d0/((pi**1.5d0)*(gw_trial**3)) ! Trial Density Coefficient
    end if
    hgp = 2.d0*(rgcut+maxval(cv))/ngp
    do igp = 0 , ngp
        rrad(igp)= igp*hgp
        if (rrad(igp)>(20.d0*gw)) then
            rho(igp) = 0.d0
        else
            rho(igp) = den_coeff*Exp(-1.d0*(rrad(igp)**2/gw**2)) ! Charge density with Q=1
        endif
        if(calc_trial) then
            if (rrad(igp)>(20.d0*0.5d0)) then
                trial_rho(igp) = 0.d0 
            else
                trial_rho(igp) = trial_den_coeff*Exp(-1.d0*(rrad(igp)**2/gw_trial**2)) ! trial Charge density with Q=1 and alpha = 0.5
            endif
        endif
    enddo
    pot(0:ngp)=0.d0
    call cal_screened_poisson(ngp,rrad,rho,scf,pot_scn)
    call cal_pot_hartree(ngp,rrad,rho,pot)
    do igp = 0 , ngp
        pot(igp)=pot(igp)-pot_scn(igp)
    end do
end subroutine
!*****************************************************************************************
subroutine get_amat_cent2_trial(parini,atoms,poisson,pot_Mg,pot_O,trial_rho,qat_target,ngp,E_par,rhs,amat)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_copy_old, update_ratp, update_rat
    use mod_electrostatics, only: typ_poisson
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    type(typ_poisson), intent(in):: poisson
    integer :: ngp
    real(8) ,intent(in) :: pot_Mg(0:ngp),pot_O(0:ngp),trial_rho(0:ngp),qat_target(1:atoms%nat)
    real(8) ,intent(inout):: amat(1:atoms%nat+1,1:atoms%nat+1)
    real(8) ,intent(inout):: E_par(1:atoms%ntrial,1:atoms%nat)
    real(8) ,intent(inout):: rhs(1:atoms%nat)
    !LOCAL Variables!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: iat,jat,igx,igy,igz,kat
    integer :: nbgx, nbgy, nbgz, linearGridNumber
    real(8) :: dx ,dy ,dz ,dr, hgp, tt
    real(8) :: grid_pot(poisson%ngpx,poisson%ngpy,poisson%ngpz)
    real(8) :: grid_trial_rho(poisson%ngpx,poisson%ngpy,poisson%ngpz)
    hgp = 2.d0*(poisson%rgcut+maxval(atoms%cellvec))/ngp
    nbgx = int(poisson%rgcut/poisson%hgrid(1,1))+3
    nbgy = int(poisson%rgcut/poisson%hgrid(2,2))+3
    nbgz = int(poisson%rgcut/poisson%hgrid(3,3))+3
    E_par(1:atoms%ntrial,1:atoms%nat)=0.d0
    amat(1:atoms%nat+1,1:atoms%nat+1)=0.d0
    do iat = 1 , atoms%nat
        if(trim(atoms%sat(iat))=='Mg') then
            do igx = 1 , poisson%ngpx
                do igy = 1 , poisson%ngpy
                    do igz = 1 , poisson%ngpz
                        dx = (-nbgx+igx)*poisson%hgrid(1,1)-atoms%ratp(1,iat)
                        dy = (-nbgy+igy)*poisson%hgrid(2,2)-atoms%ratp(2,iat)
                        dz = (-nbgz+igz)*poisson%hgrid(3,3)-atoms%ratp(3,iat)
                        dr = sqrt(dx**2+dy**2+dz**2)
                        linearGridNumber=floor(dr/hgp)
                        grid_pot(igx,igy,igz)=(dr/hgp-linearGridNumber)*(pot_Mg(linearGridNumber+1)&
                                              -pot_Mg(linearGridNumber))+pot_Mg(linearGridNumber)
                    end do
                end do
            end do
        elseif(trim(atoms%sat(iat))=='O') then
            do igx = 1 , poisson%ngpx
                do igy = 1 , poisson%ngpy
                    do igz = 1 , poisson%ngpz
                        dx = (-nbgx+igx)*poisson%hgrid(1,1)-atoms%ratp(1,iat)
                        dy = (-nbgy+igy)*poisson%hgrid(2,2)-atoms%ratp(2,iat)
                        dz = (-nbgz+igz)*poisson%hgrid(3,3)-atoms%ratp(3,iat)
                        dr = sqrt(dx**2+dy**2+dz**2)
                        linearGridNumber=floor(dr/hgp)
                        grid_pot(igx,igy,igz)=(dr/hgp-linearGridNumber)*(pot_O(linearGridNumber+1)&
                                              -pot_O(linearGridNumber))+pot_O(linearGridNumber)
                    end do
                end do
            end do
        endif
        rhs(iat)=0.d0
        do kat = 1 , atoms%ntrial 
            tt = 0.d0
            do igx = 1 , poisson%ngpx
                do igy = 1 , poisson%ngpy
                    do igz = 1 , poisson%ngpz
                        dx = (-nbgx+igx)*poisson%hgrid(1,1)-atoms%ratp(1,atoms%trial_ref_nat(kat))-atoms%trial_ref_disp(1,kat)
                        dy = (-nbgy+igy)*poisson%hgrid(2,2)-atoms%ratp(2,atoms%trial_ref_nat(kat))-atoms%trial_ref_disp(2,kat)
                        dz = (-nbgz+igz)*poisson%hgrid(3,3)-atoms%ratp(3,atoms%trial_ref_nat(kat))-atoms%trial_ref_disp(3,kat)
                        dr = sqrt(dx**2+dy**2+dz**2)
                        linearGridNumber=floor(dr/hgp)
                        grid_trial_rho(igx,igy,igz)=(dr/hgp-linearGridNumber)*(trial_rho(linearGridNumber+1)&
                                              -trial_rho(linearGridNumber))+trial_rho(linearGridNumber)
                        tt=tt+grid_trial_rho(igx,igy,igz)*grid_pot(igx,igy,igz)
                    end do
                end do
            end do 
            E_par(kat,iat)=tt*poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3)!-E_par_0(atoms%trial_ref_nat(kat),iat)
            rhs(iat)=rhs(iat)+E_par(kat,iat)*atoms%trial_ref_energy(kat)
        end do !kat
        rhs(iat)=2.d0*(rhs(iat)+parini%pen_coeff*qat_target(iat))
    end do !iat
    do iat = 1, atoms%nat
        do jat = iat,atoms%nat
            do kat = 1 ,atoms%ntrial
                amat(iat,jat)=amat(iat,jat)+E_par(kat,iat)*E_par(kat,jat)
            end do
            amat(iat,jat)=2.d0*amat(iat,jat)
            amat(jat,iat)=amat(iat,jat)
        end do
        amat(iat,iat)=amat(iat,iat)+2.d0*parini%pen_coeff
        amat(iat,atoms%nat+1)=1.d0
        amat(atoms%nat+1,iat)=1.d0
    end do
    amat(atoms%nat+1,atoms%nat+1)=0.d0
end subroutine
!*****************************************************************************************
subroutine cal_min_cv(atoms,vac)
    use mod_atoms, only: typ_atoms, update_ratp, update_rat
    type(typ_atoms), intent(inout):: atoms
    real(8),intent(in) :: vac
    ! local variables
    integer :: iat
    real(8) :: xmax,xmin
    real(8) :: ymax,ymin
    real(8) :: zmax,zmin
    if(atoms%cellvec(2,1)/=0.d0 .or. atoms%cellvec(3,1)/=0.d0 .or. &
       atoms%cellvec(1,2)/=0.d0 .or. atoms%cellvec(3,2)/=0.d0 .or. &
       atoms%cellvec(1,3)/=0.d0 .or. atoms%cellvec(2,3)/=0.d0) then
        write(*,'(a)') 'ERROR: cal_min_cv routine is only for orthogonal cell:'
        write(*,'(3es14.5)') atoms%cellvec(1,1),atoms%cellvec(2,1),atoms%cellvec(3,1)
        write(*,'(3es14.5)') atoms%cellvec(1,2),atoms%cellvec(2,2),atoms%cellvec(3,2)
        write(*,'(3es14.5)') atoms%cellvec(1,3),atoms%cellvec(2,3),atoms%cellvec(3,3)
        stop
    endif
    call update_ratp(atoms)
    xmin= minval(atoms%ratp(1,:))
    xmax= maxval(atoms%ratp(1,:))
    ymin= minval(atoms%ratp(2,:))
    ymax= maxval(atoms%ratp(2,:))
    zmin= minval(atoms%ratp(3,:))
    zmax= maxval(atoms%ratp(3,:))
    atoms%cellvec(1,1)=(xmax-xmin)+(2.d0*vac)
    atoms%cellvec(2,2)=(ymax-ymin)+(2.d0*vac)
    atoms%cellvec(3,3)=(zmax-zmin)+(2.d0*vac)
    do iat=1,atoms%nat
        atoms%ratp(1,iat)=atoms%ratp(1,iat)-xmin+vac
        atoms%ratp(2,iat)=atoms%ratp(2,iat)-ymin+vac
        atoms%ratp(3,iat)=atoms%ratp(3,iat)-zmin+vac
    enddo
    call update_rat(atoms)
end subroutine cal_min_cv
!*****************************************************************************************
end module mod_qat_target
!*****************************************************************************************
