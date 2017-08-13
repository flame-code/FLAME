! !!--------------------------------------------------------------------------------------------------
subroutine best_charge_density(parini)
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
    integer:: istat, igpx, igpy, igpz, iat, nx ,ny,nz,iter,iter_max,gweiat,i,j,jj
    real(8):: cell(3), epot, rgcut_a, t1, t2, t3, t4, pi,h,rho_err
    real(8):: ehartree,exit_cond,SDA,ener,dft_ener,coeff
    real(8),allocatable::gwit(:),dft_fat(:,:),cent_fat(:,:),agrad(:)
    real(8),allocatable:: dft_rho(:,:,:) , cent_rho(:,:,:)
    real(8),allocatable::err_qgrad(:,:),err_agrad(:,:), A(:,:) , Q(:,:)
    real(8):: sd_s,fd_s,err_fdm,err_cent
    real(8):: q1,q2,q3,a11,a12,a21,a22,a31,a32
    open(unit=1370,file='param.txt')
    pi=4.d0*atan(1.d0)
    
    call cube_read('electronic_density.cube',atoms,cent%ewald_p3d%poisson_p3d%typ_poisson)
    allocate(dft_rho(cent%ewald_p3d%poisson_p3d%ngpx,cent%ewald_p3d%poisson_p3d%ngpy,cent%ewald_p3d%poisson_p3d%ngpz))
    cent%ewald_p3d%poisson_p3d%rho=-1.d0*cent%ewald_p3d%poisson_p3d%rho
    dft_rho = cent%ewald_p3d%poisson_p3d%rho
    call acf_read(parini,'posinp.acf',1,atoms=atoms)
    
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
    allocate(dft_fat(3,atoms%nat),cent_fat(3,atoms%nat))
    cent%gwit = 1.3d0
    atoms%fat = 0.d0
    dft_fat = 0.d0
    call gauss_force(parini,'bulk',atoms%nat,atoms%rat,atoms%cellvec,atoms%zat,cent%gwit,cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%pot,atoms%fat)
    dft_fat = atoms%fat
    call destruct_ewald_bps(cent%ewald_p3d)
    deallocate(cent%ewald_p3d%poisson_p3d%pot)
 
    !/////////////////CENT_FORCE////////////////////////! 
    read(1370,*)
    read(1370,*) q1,q2,q3
    read(1370,*)
    read(1370,*) a11,a12
    read(1370,*)
    read(1370,*) a21,a22
    read(1370,*)
    read(1370,*) a31,a32
    read(1370,*)
    read(1370,*) sd_s
    read(1370,*)
    read(1370,*) fd_s
    write(*,*) "params : ",q1,q2,q3,a11,a12,a21,a22,a31,a32,sd_s,fd_s
    allocate(Q(atoms%nat,3),A(atoms%nat,3))
    Q(1:atoms%nat,1)=q1
    Q(1:atoms%nat,2)=q2
    Q(1:atoms%nat,3)=q3
    A(1:4,1) = a11
    A(4:8,1) = a12
    A(1:4,2) = a21
    A(4:8,2) = a22
    A(1:4,3) = a31
    A(4:8,3) = a32
    !atoms%qat = -2.d0
    !SD_s  = 1.d+1
    !FD_s = 1.d-4

    !allocate(cent%gwe(atoms%nat))
    !cent%gwe = 2.d0
    allocate(cent_rho(cent%ewald_p3d%poisson_p3d%ngpx,cent%ewald_p3d%poisson_p3d%ngpy,cent%ewald_p3d%poisson_p3d%ngpz))
    allocate(cent%ewald_p3d%poisson_p3d%pot(cent%ewald_p3d%poisson_p3d%ngpx,cent%ewald_p3d%poisson_p3d%ngpy,&
        cent%ewald_p3d%poisson_p3d%ngpz),stat=istat)
    allocate(err_qgrad(atoms%nat,3),err_agrad(atoms%nat,3))
    
    atoms%zat = 1.d0
    cent%gwit = 1.3d0
    
    do iter=1,huge(iter_max)
         
         cent_rho=0.d0
         cent%ewald_p3d%poisson_p3d%pot = 0.d0
         !cent%ewald_p3d%rgcut=6.d0*maxval(abs(cent%gwe))
         cent%ewald_p3d%rgcut=6.d0*maxval(abs(A))
         do j = 1 , 3
              if (j == 1) then
                  call gauss_grid(parini,atoms%boundcond,.true.,atoms%nat,atoms%rat,atoms%cellvec,Q(:,j),&
                      A(:,j),cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%rho)
              else
                  call gauss_grid(parini,atoms%boundcond,.false.,atoms%nat,atoms%rat,atoms%cellvec,Q(:,j),&
                      A(:,j),cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%rho)
              end if
         end do
         cent_rho = cent%ewald_p3d%poisson_p3d%rho
         rho_err = sqrt(sum((cent_rho-dft_rho)**2)/(nx*ny*nz))
         
         atoms%boundcond='bulk'
         call construct_ewald_bps(parini,atoms,cent%ewald_p3d)
         call cal_hartree_pot_bps(cent%ewald_p3d,atoms,ehartree)
         atoms%fat = 0.d0
         call gauss_force(parini,'bulk',atoms%nat,atoms%rat,atoms%cellvec,atoms%zat,cent%gwit,cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%pot,atoms%fat)
         err_cent = dot_product((atoms%fat(1,:)-dft_fat(1,:)),(atoms%fat(1,:)-dft_fat(1,:)))+&
                    dot_product((atoms%fat(2,:)-dft_fat(2,:)),(atoms%fat(2,:)-dft_fat(2,:)))+&
                    dot_product((atoms%fat(3,:)-dft_fat(3,:)),(atoms%fat(3,:)-dft_fat(3,:)))
         cent_fat = atoms%fat
         call destruct_ewald_bps(cent%ewald_p3d)

         do i = 1 , atoms%nat
            do j = 1 , 3
                 A(i,j) = A(i,j) + FD_s
                 cent%ewald_p3d%poisson_p3d%pot = 0.d0
                 do jj = 1 , 3
                    if (jj==1) then
                        call gauss_grid(parini,atoms%boundcond,.true.,atoms%nat,atoms%rat,atoms%cellvec,Q(:,jj),&
                            A(:,jj),cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%rho)
                    else
                        call gauss_grid(parini,atoms%boundcond,.false.,atoms%nat,atoms%rat,atoms%cellvec,Q(:,jj),&
                            A(:,jj),cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%rho)
                    end if
                 end do
                 atoms%boundcond='bulk'
                 call construct_ewald_bps(parini,atoms,cent%ewald_p3d)
                 call cal_hartree_pot_bps(cent%ewald_p3d,atoms,ehartree)
                 atoms%fat = 0.d0
                 call gauss_force(parini,'bulk',atoms%nat,atoms%rat,atoms%cellvec,atoms%zat,cent%gwit,&
                     cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%pot,atoms%fat)
                 err_fdm = dot_product((atoms%fat(1,:)-dft_fat(1,:)),(atoms%fat(1,:)-dft_fat(1,:)))+&
                           dot_product((atoms%fat(2,:)-dft_fat(2,:)),(atoms%fat(2,:)-dft_fat(2,:)))+&
                           dot_product((atoms%fat(3,:)-dft_fat(3,:)),(atoms%fat(3,:)-dft_fat(3,:)))
                 err_agrad(i,j) = (err_fdm - err_cent)/fd_s
                 A(i,j) = A(i,j) - fd_s 
                 call destruct_ewald_bps(cent%ewald_p3d)

                 Q(i,j) = Q(i,j) + FD_s
                 cent%ewald_p3d%poisson_p3d%pot = 0.d0
                  do jj = 1 , 3
                    if (jj==1) then
                        call gauss_grid(parini,atoms%boundcond,.true.,atoms%nat,atoms%rat,atoms%cellvec,Q(:,jj),&
                            A(:,jj),cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%rho)
                    else
                        call gauss_grid(parini,atoms%boundcond,.false.,atoms%nat,atoms%rat,atoms%cellvec,Q(:,jj),&
                            A(:,jj),cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%rho)
                    end if
                 end do
                 atoms%boundcond='bulk'
                 call construct_ewald_bps(parini,atoms,cent%ewald_p3d)
                 call cal_hartree_pot_bps(cent%ewald_p3d,atoms,ehartree)
                 atoms%fat = 0.d0
                 call gauss_force(parini,'bulk',atoms%nat,atoms%rat,atoms%cellvec,atoms%zat,cent%gwit,cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%pot,atoms%fat)
                 err_fdm = dot_product((atoms%fat(1,:)-dft_fat(1,:)),(atoms%fat(1,:)-dft_fat(1,:)))+&
                           dot_product((atoms%fat(2,:)-dft_fat(2,:)),(atoms%fat(2,:)-dft_fat(2,:)))+&
                           dot_product((atoms%fat(3,:)-dft_fat(3,:)),(atoms%fat(3,:)-dft_fat(3,:)))

                 err_qgrad(i,j) = (err_fdm - err_cent)/fd_s
                 Q(i,j) = Q(i,j) - fd_s
                 call destruct_ewald_bps(cent%ewald_p3d)
             end do
         end do 
         exit_cond = sqrt((err_cent))
         write(*,'(a,i4,2es14.5,2(2x,8f6.2))')"SD ",iter,exit_cond,rho_err,atoms%qat,cent%gwe
         write(*,'(a,i4,2(2x,8es14.4))') "a-q-grad",iter,err_agrad,err_qgrad
         if (exit_cond < 1.d-3)then
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
 end subroutine best_charge_density
! 
! ! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Energy SD! ! !subroutine best_charge_density(parini)
!Energy SD! ! !     use mod_interface
!Energy SD! ! ! !!--------------------------------------------------------------------------------------------------
!Energy SD! ! !     use mod_parini, only: typ_parini
!Energy SD! ! !     use mod_electrostatics, only: typ_poisson_p3d, typ_ewald_p3d,typ_poisson
!Energy SD! ! !     use mod_atoms, only: typ_atoms
!Energy SD! ! !     use mod_ann, only: typ_cent, typ_ann_arr
!Energy SD! ! !     implicit none
!Energy SD! ! !     type(typ_parini), intent(in):: parini
!Energy SD! ! !     !local variables
!Energy SD! ! !     type(typ_poisson_p3d):: poisson_p3d
!Energy SD! ! !     type(typ_ewald_p3d):: ewald_p3d ,ewald_p3d_rough
!Energy SD! ! !     type(typ_atoms):: atoms
!Energy SD! ! !     type(typ_cent):: cent
!Energy SD! ! !     type(typ_ann_arr):: ann_arr
!Energy SD! ! !     integer:: istat, igpx, igpy, igpz, iat, nx ,ny, nz,iter,iter_max,gweiat
!Energy SD! ! !     real(8):: cell(3), epot, rgcut_a, t1, t2, t3, t4, pi,h,rho_err
!Energy SD! ! !     real(8):: ehartree,exit_cond,SDA,ener,dft_ener,coeff
!Energy SD! ! !     real(8),allocatable::  gausswidth(:) , gwit(:), diff(:,:),bigforce(:,:),agrad(:)
!Energy SD! ! !     real(8),allocatable:: dft_rho(:,:,:) , cent_rho(:,:,:)
!Energy SD! ! !     pi=4.d0*atan(1.d0)
!Energy SD! ! !     call cube_read('electronic_density.cube',atoms,cent%ewald_p3d%poisson_p3d%typ_poisson)
!Energy SD! ! !     allocate(dft_rho(cent%ewald_p3d%poisson_p3d%ngpx,cent%ewald_p3d%poisson_p3d%ngpy,cent%ewald_p3d%poisson_p3d%ngpz))
!Energy SD! ! !     dft_rho = cent%ewald_p3d%poisson_p3d%rho
!Energy SD! ! !     call acf_read(parini,'posinp.acf',1,atoms=atoms)
!Energy SD! ! !     cent%ewald_p3d%hgx=cent%ewald_p3d%poisson_p3d%hx
!Energy SD! ! !     cent%ewald_p3d%hgy=cent%ewald_p3d%poisson_p3d%hy
!Energy SD! ! !     cent%ewald_p3d%hgz=cent%ewald_p3d%poisson_p3d%hz
!Energy SD! ! !     atoms%boundcond='bulk'
!Energy SD! ! !     allocate(cent%ewald_p3d%poisson_p3d%pot(cent%ewald_p3d%poisson_p3d%ngpx,cent%ewald_p3d%poisson_p3d%ngpy,cent%ewald_p3d%poisson_p3d%ngpz),stat=istat)
!Energy SD! ! !     call construct_ewald_bps(parini,atoms,cent%ewald_p3d)
!Energy SD! ! !     call cal_hartree_pot_bps(cent%ewald_p3d,atoms,ehartree)
!Energy SD! ! !     dft_ener = ehartree
!Energy SD! ! !     write(*,*) "dft_Hartree : " , dft_ener
!Energy SD! ! ! !   atoms%zat(:)=1.d0
!Energy SD! ! !     nx=cent%ewald_p3d%poisson_p3d%ngpx
!Energy SD! ! !     ny=cent%ewald_p3d%poisson_p3d%ngpy
!Energy SD! ! !     nz=cent%ewald_p3d%poisson_p3d%ngpz
!Energy SD! ! ! !   allocate(cent%gwit(atoms%nat))
!Energy SD! ! ! !   allocate(bigforce(3,atoms%nat))
!Energy SD! ! ! !   cent%gwit = 1.3d0
!Energy SD! ! ! !   call gauss_force(parini,'bulk',atoms%nat,atoms%rat,atoms%cellvec,atoms%zat,cent%gwit,cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%pot,atoms%fat)
!Energy SD! ! ! !   bigforce = atoms%fat
!Energy SD! ! !     call destruct_ewald_bps(cent%ewald_p3d)
!Energy SD! ! !     deallocate(cent%ewald_p3d%poisson_p3d%pot)
!Energy SD! ! !    !/////////////////CENT_FORCE////////////////////////! 
!Energy SD! ! !     atoms%qat = 2.d0
!Energy SD! ! !    ! atoms%qat =[0.991363025367481,0.991363025367481,0.991363025367481,0.991363025367482,0.991363025367481,0.991363025367481,0.991363025367481,0.991363025367481]
!Energy SD! ! ! 
!Energy SD! ! !     SDA  = 1.d-3
!Energy SD! ! !     allocate(cent%gwe(atoms%nat))
!Energy SD! ! !     cent%gwe(1:4)=1d0
!Energy SD! ! !     cent%gwe(1:8)=2d0
!Energy SD! ! !     !cent%gwe = [0.179276855138007,0.179276855138007,0.179276855138007,0.179276855138007,0.179276855138008,0.179276855138007,0.179276855138007,0.179276855138007]
!Energy SD! ! !     allocate(cent_rho(cent%ewald_p3d%poisson_p3d%ngpx,cent%ewald_p3d%poisson_p3d%ngpy,cent%ewald_p3d%poisson_p3d%ngpz))
!Energy SD! ! !     do iter=1,huge(iter_max)
!Energy SD! ! !          cent%ewald_p3d%rgcut=6.d0*maxval(abs(cent%gwe))
!Energy SD! ! !          call gauss_grid(parini,atoms%boundcond,.true.,atoms%nat,atoms%rat,atoms%cellvec,atoms%qat,cent%gwe,cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%rho)
!Energy SD! ! !          cent_rho = cent%ewald_p3d%poisson_p3d%rho
!Energy SD! ! !          rho_err = sqrt(sum((cent_rho-dft_rho)**2)/(cent%ewald_p3d%poisson_p3d%ngpx**3))
!Energy SD! ! !          !call cube_write('cent_electronic_density.cube',atoms,cent%ewald_p3d%poisson_p3d%typ_poisson,'rho')
!Energy SD! ! !          !-----------------------------------------------------------------------------------
!Energy SD! ! !          !write(*,*) "calculating gauss_force of CENT2 charge density..."
!Energy SD! ! !          !write(*,*) "reading cent electronic density..."
!Energy SD! ! !          !call cube_read('cent_electronic_density.cube',atoms,cent%ewald_p3d%poisson_p3d%typ_poisson)
!Energy SD! ! ! !!!         cent%ewald_p3d%hgx=cent%ewald_p3d%poisson_p3d%hx
!Energy SD! ! ! !!!         cent%ewald_p3d%hgy=cent%ewald_p3d%poisson_p3d%hy
!Energy SD! ! ! !!!         cent%ewald_p3d%hgz=cent%ewald_p3d%poisson_p3d%hz
!Energy SD! ! !          atoms%boundcond='bulk'
!Energy SD! ! !          allocate(cent%ewald_p3d%poisson_p3d%pot(cent%ewald_p3d%poisson_p3d%ngpx,cent%ewald_p3d%poisson_p3d%ngpy,cent%ewald_p3d%poisson_p3d%ngpz),stat=istat)
!Energy SD! ! !          call construct_ewald_bps(parini,atoms,cent%ewald_p3d)
!Energy SD! ! !          call cal_hartree_pot_bps(cent%ewald_p3d,atoms,ehartree)
!Energy SD! ! !          ener = ehartree
!Energy SD! ! !          write(*,*) "Hartree : " , ener ;
!Energy SD! ! ! !         atoms%zat(:)=1.d0
!Energy SD! ! ! !         cent%gwit = 1.3d0
!Energy SD! ! ! !         atoms%fat = 0.d0
!Energy SD! ! ! !         call gauss_force(parini,'bulk',atoms%nat,atoms%rat,atoms%cellvec,atoms%zat,cent%gwit cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%pot,atoms%fat)
!Energy SD! ! !          !write(*,*) "norm",cent%gwit(1),sqrt(sum(((bigforce(1:3,1:atoms%nat)-atoms%fat(1:3,1:atoms%nat))**2)))
!Energy SD! ! !          allocate(cent%rgrad(1:3,atoms%nat),cent%qgrad(atoms%nat),agrad(atoms%nat))
!Energy SD! ! !          call gauss_gradient_rzx(parini,atoms%boundcond,atoms%nat,atoms%rat,atoms%cellvec,atoms%qat,cent%gwe,cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%pot,cent%rgrad,cent%qgrad,agrad)
!Energy SD! ! !         ! write(*,*)"rgrad",cent%rgrad
!Energy SD! ! !         ! write(*,*)"qgrad",cent%qgrad
!Energy SD! ! !         ! write(*,*)"agrad",agrad
!Energy SD! ! !         ! write(*,*)"pot" ,shape(cent%ewald_p3d%poisson_p3d%pot)
!Energy SD! ! ! !         exit_cond = sqrt(sum(((bigforce(1:3,1:atoms%nat)-atoms%fat(1:3,1:atoms%nat))**2))/atoms%nat)
!Energy SD! ! !          exit_cond =sqrt( (dft_ener-ener)**2)
!Energy SD! ! !          write(*,'(a,i4,2es14.5,2(2x,8f6.2))')"SD ",iter,exit_cond,rho_err,atoms%qat,cent%gwe
!Energy SD! ! !          if (exit_cond < 1.d-4)then
!Energy SD! ! !             write(*,*) "max conversion reached"
!Energy SD! ! !             write(*,*) "qat : "
!Energy SD! ! !             write(*,*) atoms%qat
!Energy SD! ! !             write(*,*) "gwe : "
!Energy SD! ! !             write(*,*) cent%gwe
!Energy SD! ! !             exit
!Energy SD! ! !          endif
!Energy SD! ! !          coeff = 2.d0*(ener - dft_ener)
!Energy SD! ! !          atoms%qat = atoms%qat - SDA*coeff*cent%qgrad
!Energy SD! ! !          cent%gwe = cent%gwe - SDA*coeff*agrad
!Energy SD! ! !          call destruct_ewald_bps(cent%ewald_p3d)
!Energy SD! ! !          deallocate(cent%ewald_p3d%poisson_p3d%pot,cent%rgrad,cent%qgrad,agrad)
!Energy SD! ! !     enddo
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
!!Energy SD! ! !end subroutine best_charge_density

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
