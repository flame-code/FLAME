!*****************************************************************************************
module bar_saddle_params
implicit none
integer:: ndim,n_anchor,maxit,ncontr
real(8):: dbar,alpha_bar,bar_tol,contr_dbar,contr_bar_tol,fnrm_switch,dbar_relax
logical:: hybrid,restart,bar_contract,relax_init,relax_final

contains
  subroutine read_bar_saddle_params(nat,parini)
  use mod_parini, only: typ_parini
  integer:: nat
  type(typ_parini), intent(in):: parini
  character(100):: line
!4                    #nanch
!T T                  #relax inital endpoints, relax final endpoints
!5.d-1 14.d0 5.d-3    #parameters for saddle search: dbar, alpha_bar, bar_tol
!100                  #Maxit
!F F                  #Relax endpoints, Relax from final barends,dbar_relax
!T                    #Hybrid,fnrm_switch
!F                    #Restart
!T 1.d-1 5.d-5 10     #Bar contraction, final dbar, final fmax_tol, ncontract
  !Read input parameters
  !open(unit=2,file="input.barsad")
  !read(2,*) n_anchor
  ndim=3*nat !read(2,*) ndim
  !read(2,*) dbar,alpha_bar, bar_tol
  !read(2,*) maxit
  !read(2,*) relax_init,relax_final,dbar_relax
  !read(2,*) hybrid,fnrm_switch
  !read(2,*) restart
  !read(2,'(a100)') line;read(line,*) bar_contract;if(bar_contract) read(line,*) bar_contract,contr_dbar,contr_bar_tol,ncontr
  !close(2)
  dbar=parini%dbar
  alpha_bar=parini%alphax_bs
  bar_tol=parini%fnrmtol_coarse
  maxit=parini%nstep_bs
  bar_contract=parini%bar_contract
  contr_dbar=parini%contr_dbar
  contr_bar_tol=parini%fnrmtol_contracted
  ncontr=parini%nstep_contract
  end subroutine
end module bar_saddle_params

module spline_params
implicit none
save
integer:: n_anc,n_dim
real(8),allocatable:: sp_param(:,:),coord(:,:),t(:)
integer,private:: n,m
contains
   subroutine allocate_spline_params(n,m)
   integer :: n,m
   n_anc=n
   n_dim=m
   if(.not.allocated(sp_param)) allocate(sp_param(n_dim,n_anc))
   if(.not.allocated(coord)) allocate(coord(n_dim,n_anc))
   if(.not.allocated(t)) allocate(t(n_anc))
   end subroutine
   subroutine deallocate_spline_params()
   if(allocated(sp_param)) deallocate(sp_param)
   if(allocated(coord)) deallocate(coord)
   if(allocated(t)) deallocate(t)
   end subroutine
end module

subroutine bar_saddle(parini)
    use mod_parini, only: typ_parini
    use dictionaries
    use dynamic_memory
    use yaml_output
    use mod_processors, only: nproc, iproc
    use mod_atoms, only: typ_atoms, atom_deallocate, typ_atoms_arr
    use mod_atoms, only: atom_copy, get_rat
    use mod_potential, only: potential
    use mod_yaml_conf, only: read_yaml_conf
    use bar_saddle_params
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    integer :: iconf, infocode
    type(typ_atoms) :: atoms
    type(typ_atoms) :: atoms_2
    type(typ_atoms_arr) :: atoms_arr
    real(8), dimension(:,:), allocatable :: fxyz1
    real(8), dimension(:,:), pointer :: rxyz1, rxyz2

    !Local variables
    real(8),allocatable:: bar_vec(:),bar_cm(:),bar_max(:)
    real(8),allocatable:: e_anc(:,:)
    real(8):: etot1,etot2,fnoise,fmax,fnrm,xmin,brent,tt,e_max,e_tmp
    real(8),pointer :: rxyz_tmp(:,:)

    potential=trim(parini%potential_potential)
    !The following reads a maximum of two configurations but
    !I am sending you a posinp.yaml that includes one configuration.
    call read_yaml_conf(parini,'posinp.yaml',2,atoms_arr)
    call atom_copy(atoms_arr%atoms(1),atoms,  'atoms_arr%atoms(1)=>atoms')
    if (atoms_arr%nconf == 2) then
        call atom_copy(atoms_arr%atoms(2),atoms_2,'atoms_arr%atoms(2)=>atoms')
    endif
    do iconf=1,atoms_arr%nconf
        call atom_deallocate(atoms_arr%atoms(iconf))
    enddo
    deallocate(atoms_arr%atoms)
    allocate(bar_vec(3*atoms%nat),bar_cm(3*atoms%nat),bar_max(3*atoms%nat))
    allocate(rxyz1(3,atoms%nat),rxyz2(3,atoms%nat))
    allocate(fxyz1(3,atoms%nat))
    !Let's try to implement saddle method where intiial position is a guess for
    !the saddle point. For now, we will consider rxyz1 as the CM of the bar, and
    !rxyz2 as one of the end points of the bar. However, since we are using
    !rxyz1=rxyz2, we will have to randomply displace rxyz1 and scale the bar
    !along the direction.
    !Eventually, I would like to use rxyz1 as the estimated saddle point, and
    !rxyz2 as one of the end points of the saddle bar
    call get_rat(atoms,rxyz1)
    if (atoms_arr%nconf == 2) then
        call get_rat(atoms_2,rxyz2)
    endif
    call init_potential_forces(parini,atoms)
    !-------------------------------------------------------
    !-------------------------------------------------------
    !-------------------------------------------------------
    !-------------------------------------------------------
    !-------------------------------------------------------
    !MAX, call your routine here.
    ndim=atoms%nat*3
    !allocate(e_anc(2,n_anchor))
    !e_anc=0.d0
    !e_anc(2,:)=1.d0
    
    !write(*,*) "SAD: ndim", ndim
    !write(*,*) "SAD: fmax_tol,fnrm_tol,alpha0",fmax_tol,fnrm_tol,alpha0
    !write(*,*) "SAD: dbar,alpha_bar",dbar,alpha_bar
    !write(*,*) "SAD: bar_tol",bar_tol
    
    !Allocate arrays
    !allocate(rxyz1(ndim),fxyz1(ndim),rxyz2(ndim),fxyz2(ndim),bar_vec(ndim),bar_cm(ndim))
    
    !!Read posinp1.xyz
    !filename="posinp1.xyz"
    !call read_xyz(filename,ndim/3,rxyz1)
    !Read posinp2.xyz
!    filename="posinp2.xyz"
!    call read_xyz(filename,rxyz2,e_tmp,atoms,iproc,nproc)
    
    !Relax initial endpoints if requested
    !if(relax_init) then
    !   e_anc(2,1)=0.d0
    !   call geopt(bigdft_mpi%nproc,bigdft_mpi%iproc,rxyz1,atoms,fxyz1,strten,e_anc(1,1),rst_left,inputs,ncount_bigdft)
    !   if(iproc==0) call write_atomic_file('posinp_relaxed',e_anc(1,1),rxyz1,atoms,"Relaxed posinp.xyz",forces=fxyz1)
    !   e_anc(2,n_anchor)=0.d0
    !   call geopt(bigdft_mpi%nproc,bigdft_mpi%iproc,rxyz2,atoms,fxyz2,strten,e_anc(1,n_anchor),rst_left,inputs,ncount_bigdft)
    !   if(iproc==0) call write_atomic_file('posinp2_relaxed',e_anc(1,n_anchor),rxyz2,atoms,"Relaxed posinp2.xyz",forces=fxyz2)
    !   if(iproc==0) write(*,'(a,2(es15.7))') "SAD: Relax initial endpoints with energies: ",e_anc(1,1), e_anc(1,n_anchor)
    !endif
!n_anchor=10
!call init_bar(ndim,n_anchor,rxyz1,rxyz2,dbar,bar_vec,bar_cm,e_anc,e_max,tt,nproc,iproc,atoms,rst_left,inputs,lo_inputs,ncount_bigdft)
    !Initiallize random bar, and get its energy at the center
    call read_bar_saddle_params(atoms%nat,parini)
    if (atoms_arr%nconf == 1) then
        call init_bar_random(parini,ndim,rxyz1,dbar,bar_vec,bar_cm)
    elseif (atoms_arr%nconf == 2) then
        call init_bar_dir(ndim,rxyz1,rxyz2,dbar,bar_vec,bar_cm)
    endif
    call call_bigdft(nproc,iproc,atoms,bar_cm,etot1,fxyz1,fnoise,infocode,parini)
!if(ncount_bigdft.gt.maxit) then
!   if(iproc==0) write(*,'(a)') "SAD: MAXIT reached, exiting"
!   goto 1001
!endif
    call find_saddle(ndim,bar_vec,bar_cm,fxyz1,etot1,e_max,bar_max,alpha_bar,bar_tol,atoms,parini)
!if(ncount_bigdft.gt.maxit) then
!   if(iproc==0) write(*,'(a)') "SAD: MAXIT reached, exiting"
!   goto 1001
!endif

!!!!!!!!!!!!!call lenjon(ndim/3,bar_cm,fxyz1,etot2)
!!!!!!!!!!!!!call fmaxfnrm(ndim,fxyz1,fmax,fnrm)
!!!!!!!!!!!!!call call_bigdft(nproc,iproc,atoms,bar_cm,inputs,etot2,fxyz1,strten,fnoise,rst_left,infocode)
!!!!!!!!!!!!call call_bigdft(nproc,iproc,atoms,bar_max,inputs,etot2,fxyz1,strten,fnoise,rst_left,infocode)
!!!!!!!!!!!!ncount_bigdft=ncount_bigdft+1
!!!!!!!!!!!!call fnrmandforcemax(fxyz1,fnrm,fmax,atoms%nat)
!!!!!!!!!!!!if(iproc==0) then
!!!!!!!!!!!!      filename="possaddle.final"
!!!!!!!!!!!!      call write_atomic_file(trim(filename),etot2,bar_max,atoms,"Final Saddle",forces=fxyz1)
!!!!!!!!!!!!endif
!!!!!!!!!!!!fnrm=sqrt(fnrm)
!!!!!!!!!!!!if(iproc==0) write(*,'(a,3(es25.15))') "SAD: Final values of saddle point, bar_max: etot,fmax,fnrm           ",etot2,fmax,fnrm
!!!!!!!!!!!!if(iproc==0) write(*,'(a,2(es25.15))') "SAD: Compare saddle point                 : app_etot-etot,e_max-etot ",etot1-etot2,e_max-etot2
!!!!!!!!!!!!
!!!!!!!!!!!!!Compute barrier heights
!!!!!!!!!!!!if(iproc==0) write(*,'(a)') "SAD: Comparing barrier heights, only accurate if posinp.xyz and posinp2.xyz were intially relaxed or from restart with correct energies in posanchors"
!!!!!!!!!!!!if(iproc==0) write(*,'(a,2(es25.15))') "SAD: saddle-posinp.xyz, saddle-posinp2.xyz: ",etot2-e_anc(1,1),etot2-e_anc(1,n_anchor)
!!!!!!!!!!!!
!!!!!!!!!!!!
!!!!!!!!!!!!
!!!!!!!!!!!!
!!!!!!!!!!!!!Relax final endpoints if requested
!!!!!!!!!!!!if(relax_final) then
!!!!!!!!!!!!   if(iproc==0) write(*,'(a)')  "SAD: Relaxing the bar ends"
!!!!!!!!!!!!   if(iproc==0) write(*,'(a,es15.7)')  "SAD: The bar will be extended to a value of ",dbar_relax
!!!!!!!!!!!!   rxyz1=bar_cm-bar_vec/sqrt(dot_product(bar_vec,bar_vec))*dbar_relax
!!!!!!!!!!!!   call geopt(bigdft_mpi%nproc,bigdft_mpi%iproc,rxyz1,atoms,fxyz1,strten,e_anc(1,1),rst_left,inputs,ncount_bigdft)
!!!!!!!!!!!!   if(iproc==0) call write_atomic_file('posinp_relaxed_final',e_anc(1,1),rxyz1,atoms,"Relaxed final posinp.xyz",forces=fxyz1)
!!!!!!!!!!!!   rxyz2=bar_cm+bar_vec/sqrt(dot_product(bar_vec,bar_vec))*dbar_relax
!!!!!!!!!!!!   call geopt(bigdft_mpi%nproc,bigdft_mpi%iproc,rxyz2,atoms,fxyz2,strten,e_anc(1,n_anchor),rst_left,inputs,ncount_bigdft)
!!!!!!!!!!!!   if(iproc==0) call write_atomic_file('posinp2_relaxed_final',e_anc(1,n_anchor),rxyz2,atoms,"Relaxed final posinp2.xyz",forces=fxyz2)
!!!!!!!!!!!!   if(iproc==0) write(*,'(a,2(es15.7))') "SAD: Relax final endpoints with energies:   ",e_anc(1,1), e_anc(1,n_anchor)
!!!!!!!!!!!!!Compute barrier heights
!!!!!!!!!!!!if(iproc==0) write(*,'(a)') "SAD: Comparing barrier heights, after relaxation"
!!!!!!!!!!!!if(iproc==0) write(*,'(a,2(es25.15))') "SAD: saddle-posinp.xyz, saddle-posinp2.xyz: ",etot2-e_anc(1,1),etot2-e_anc(1,n_anchor)
!!!!!!!!!!!!endif
!!!!!!!!!!!!
!!!!!!!!!!!!1001 continue
!!!!!!!!!!!!deallocate(e_anc)
    !The following line is an example to get energy and forces
    !-------------------------------------------------------
    !-------------------------------------------------------
    !-------------------------------------------------------
    !-------------------------------------------------------
    !-------------------------------------------------------
    !call call_bigdft(nproc,iproc,atoms,rxyz1,etot1,fxyz1,fnoise,infocode,parini)
    !call yaml_map('epot',etot1,fmt='(e24.15)')
    !-------------------------------------------------------
    !-------------------------------------------------------
    !-------------------------------------------------------
    !-------------------------------------------------------
    !-------------------------------------------------------
    !-------------------------------------------------------
    call final_potential_forces(parini,atoms)
    deallocate(rxyz1,rxyz2)
    deallocate(fxyz1)
    call atom_deallocate(atoms)
    call atom_deallocate(atoms_2)
end subroutine bar_saddle
!*****************************************************************************************
subroutine init_bar_random(parini,ndim,rxyz1,dbar,bar_vec,bar_cm)
!This routine will determine where and how the initial bar should be placed
!We will set bar_cm=rxyz1, and dir a random vector. Then, bar_vec=rxyz1+dir
!normalized to dbar.
!Given this bar_cm, we will create a bar at bar_cm in the direction of dir of
!length dbar 
use mod_parini, only: typ_parini
    use bar_saddle_params, only: restart
    use mod_utils
implicit none
type(typ_parini), intent(in):: parini
integer::ndim,n_anchor,i
real(8),dimension(ndim):: bar_vec,bar_cm,rxyz1,rxyz2
real(8):: dbar
logical:: file_exists
character(40):: filename
character(5):: fn
bar_cm = rxyz1
if(trim(parini%rng_type)=='only_for_tests') then
call random_number_generator_simple(ndim,bar_vec)
else
call random_number(bar_vec)
endif
bar_vec = bar_vec*2.d0 - 1.d0
!We still need to rescale the bar
!Rescale
bar_vec=bar_vec/sqrt(dot_product(bar_vec,bar_vec))*dbar
end subroutine
!*****************************************************************************************
subroutine init_bar_dir(ndim,rxyz1,rxyz2,dbar,bar_vec,bar_cm)
!This routine will determine where and how the initial bar should be placed
!We will set bar_cm=rxyz1, and dir as rxyz2-rxyz1
!normalized to dbar.
!Given this bar_cm, we will create a bar at bar_cm in the direction of dir of
!length dbar 
    use bar_saddle_params, only: restart
implicit none
integer::ndim,n_anchor,i
real(8),dimension(ndim):: bar_vec,bar_cm,rxyz1,rxyz2
real(8):: dbar
logical:: file_exists
character(40):: filename
character(5):: fn
bar_cm = rxyz1
bar_vec = rxyz2-rxyz1
!We still need to rescale the bar
!Rescale
bar_vec=bar_vec/sqrt(dot_product(bar_vec,bar_vec))*dbar
end subroutine
!*****************************************************************************************
subroutine find_saddle(n_dim,bar_vec,bar_cm,fxyz,etot,emax,bar_max,alpha_bar,bar_tol,atoms,parini)
    use mod_parini, only: typ_parini
    use dictionaries
    use dynamic_memory
    use yaml_output
    use mod_processors, only: nproc, iproc
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, atom_deallocate, typ_atoms_arr
    use mod_atoms, only: atom_copy, get_rat
    use mod_potential, only: potential
    use bar_saddle_params, only:maxit,bar_contract,contr_dbar,contr_bar_tol,ncontr,dbar,hybrid,fnrm_switch
implicit none
type(typ_parini), intent(in):: parini
type(typ_atoms):: atoms
integer:: n_dim,iter,j,contractcount
real(8),dimension(n_dim):: bar_vec,bar_cm,fxyz,bar_vec_unit
real(8),dimension(n_dim):: bar_vec_old,bar_cm_old,bar_max
real(8),dimension(n_dim,3):: fbar_para,fbar_perp,fbar_tot
real(8),dimension(n_dim,3):: fbar_para_old,fbar_perp_old,fbar_tot_old
real(8):: etot,alpha,fmax,fnrm,etot_all(3),etot_all_old(3),alpha_bar,alpha_min,alpha_max,cosinus,bar_tol,diffbar,xmax,emax
character(40):: filename
character(5):: fn
logical:: contracting,lo_acc,skip
integer:: optimizer
!Optimizer: 1-energy feedback (not stable), 2-gradient feedback (very safe!), 3-BFGS under development
call yaml_sequence_open('BARSAD optimization iterations')
optimizer=2

if(hybrid) then
  lo_acc=.true.
else
  lo_acc=.false.
endif
skip=.false.

contracting=.false.
alpha=alpha_bar
alpha_min=alpha_bar*1.d-1
alpha_max=alpha_bar*1.d2
if(iproc==0) open(unit=3,file="Saddlepath")
!write(fn,'(i1)') 3*n_dim
!write(filename,'(a)') "("//trim(trim(fn)//"(es25.15))")
!write(*,'(a)') filename
!write(3,'(6(es25.15))') bar_cm-bar_vec,bar_cm,bar_cm+bar_vec
iter=0
if(lo_acc) then
  call forcebar_decomp(n_dim,bar_vec,bar_cm,fbar_para,fbar_perp,fbar_tot,etot_all,xmax,emax,bar_max,nproc,iproc,atoms,parini,iter)
else
  call forcebar_decomp(n_dim,bar_vec,bar_cm,fbar_para,fbar_perp,fbar_tot,etot_all,xmax,emax,bar_max,nproc,iproc,atoms,parini,iter)
endif
call fnrmandforcemax(fbar_tot(:,2),fnrm,fmax,atoms%nat)
!call fmaxfnrm(n_dim,fbar_tot(:,2),fmax,fnrm)
fbar_para_old=fbar_para;fbar_perp_old=fbar_perp;fbar_tot_old=fbar_tot
bar_vec_old=bar_vec;bar_cm_old=bar_cm
etot_all_old=etot_all

if(iproc==0) then
    call yaml_sequence(advance='no')
    call yaml_mapping_open('BARSAD',flow=.true.)
    call yaml_map('iter',iter,fmt='(i5)')
    call yaml_map('epot',etot_all(2),fmt='(es15.7)')
    call yaml_map('xmax',xmax,fmt='(es10.2)')
    call yaml_map('fmax',fmax,fmt='(es15.7)')
    call yaml_map('fnrm',fnrm,fmt='(es15.7)')
    call yaml_mapping_close()
endif



!if(iproc==0) write(*,'(a,i5,es15.7,es10.2,2(es15.7))') "SAD: Bar Saddle iter,app_etot,xmax,fmax,fnrm:                       ", iter,etot_all(2),xmax,fmax,fnrm
    call fnrmandforcemax(fbar_tot(:,2),fnrm,fmax,atoms%nat)
    if(hybrid .and. lo_acc .and. (fnrm.lt.fnrm_switch)) then
      lo_acc=.false.
      skip=.true.
      if(iproc==0) call yaml_scalar("BARSAD: Switching to higher accuracy")
      !if(iproc==0) write(*,'(a)') "SAD: Switching to higher accuracy"
    endif   
    if(fnrm.lt.bar_tol) then
    if(iproc==0) call yaml_scalar("BARSAD: Saddle point search converged before entering iterations")
    !if(iproc==0) write(*,'(a)') "SAD: Saddle point search converged before entering iterations"
    goto 1002
    endif

!stating the loop over optimization
do iter=1,maxit!(maxit-ncount_bigdft)/2
!Contract dbar slowly if desired
   if(contracting.and.contractcount.lt.ncontr) then
      dbar=dbar-diffbar/real(ncontr,8)
      contractcount=contractcount+1
   endif

!With energy feedback on alpha
1001 continue
!   call move_bar(n_dim,bar_vec,bar_cm,fbar_para,fbar_perp,fbar_tot,etot_all,alpha)
select case(optimizer)
  case(1)
   call move_bar_energy(n_dim,bar_vec,bar_cm,fbar_para,fbar_perp,fbar_tot,etot_all,alpha,alpha_min,alpha_max,iter,skip)
  case(2)
   call move_bar_force(n_dim,bar_vec,bar_cm,fbar_para,fbar_perp,fbar_tot,alpha,alpha_min,alpha_max,iter,skip)
  case(3)
   call move_bar_bfgs(iproc,n_dim,bar_vec,bar_cm,fbar_para,fbar_perp,fbar_tot,etot_all,alpha,iter,20)
  case default
   stop "NOT YET IMPLEMENTED"
end select
    if(skip) skip=.false.
    if(hybrid .and. lo_acc .and. (fnrm.lt.fnrm_switch)) then
      lo_acc=.false.
      skip=.true.
      if(iproc==0)   write(*,'(a)') "BARSAD: Switching to higher accuracy"
      !if(iproc==0)   write(*,'(a)') "SAD: Switching to higher accuracy"
    endif   
if(lo_acc) then
   call forcebar_decomp(n_dim,bar_vec,bar_cm,fbar_para,fbar_perp,fbar_tot,etot_all,xmax,emax,bar_max,nproc,iproc,atoms,parini,iter)
else
   call forcebar_decomp(n_dim,bar_vec,bar_cm,fbar_para,fbar_perp,fbar_tot,etot_all,xmax,emax,bar_max,nproc,iproc,atoms,parini,iter)
endif


   call fnrmandforcemax(fbar_tot(:,2),fnrm,fmax,atoms%nat)
!   call fmaxfnrm(n_dim,fbar_tot(:,2),fmax,fnrm)
if(iproc==0) then
    call yaml_sequence(advance='no')
    call yaml_mapping_open('BARSAD',flow=.true.)
    call yaml_map('iter',iter,fmt='(i5)')
    call yaml_map('epot',etot_all(2),fmt='(es15.7)')
    call yaml_map('xmax',xmax,fmt='(es10.2)')
    call yaml_map('fmax',fmax,fmt='(es15.7)')
    call yaml_map('fnrm',fnrm,fmt='(es15.7)')
    call yaml_map('epotl',etot_all(1),fmt='(es15.7)')
    call yaml_map('epotr',etot_all(3),fmt='(es15.7)')
    call yaml_map('alpha',alpha,fmt='(es15.7)')
    call yaml_map('dbar',dbar,fmt='(es15.7)')
    call yaml_mapping_close()
endif

    !if(iproc==0)   write(*,'(a,i5,es15.7,es10.2,6(es15.7))') "SAD: Bar Saddle iter,app_etot,xmax,fmax,fnrm,left,right,alpha,dbar: ", iter,etot_all(2),xmax,fmax,fnrm,etot_all(1),etot_all(3),alpha,dbar!acos(dot_product(bar_vec,fbar_tot(:,2))/sqrt(dot_product(bar_vec,bar_vec))/sqrt(dot_product(fbar_tot(:,2),fbar_tot(:,2))))
!   write(3,trim(filename)) bar_cm-bar_vec,bar_cm,bar_cm+bar_vec


      
if(iproc==0) then
      write(fn,'(i5.5)') iter
      filename="possaddle"//fn
              call write_atomic_file(trim(filename),etot_all(2),bar_cm(:),atoms,"Approximate saddle",fbar_tot(:,2))
      filename="possaddle.left"//fn
              call write_atomic_file(trim(filename),etot_all(1),bar_cm(:)-bar_vec(:),atoms,"Left bar end",fbar_tot(:,1))
      filename="possaddle.right"//fn
              call write_atomic_file(trim(filename),etot_all(3),bar_cm(:)+bar_vec(:),atoms,"Right bar end",fbar_tot(:,3))
      filename="possaddle_cur"
              call write_atomic_file(trim(filename),etot_all(2),bar_cm(:),atoms,"Approximate saddle",fbar_tot(:,2))
      filename="possaddle.left_cur"
              call write_atomic_file(trim(filename),etot_all(1),bar_cm(:)-bar_vec(:),atoms,"Left bar end",fbar_tot(:,1))
      filename="possaddle.right_cur"
              call write_atomic_file(trim(filename),etot_all(3),bar_cm(:)+bar_vec(:),atoms,"Right bar end",fbar_tot(:,3))
!      open(unit=2,file=trim(filename))
!      write(2,*) n_dim/3
!      write(2,*) etot_all(2)
!      do j=1,n_dim/3
!       write(2,*) "A",bar_cm(1+(j-1)*3:3+(j-1)*3)
!      enddo
!      close(2)
!      write(fn,'(i5.5)') iter
!      filename="possaddle.left."//fn//".xyz"
!      open(unit=2,file=trim(filename))
!      write(2,*) n_dim/3
!      write(2,*) etot_all(2)
!      do j=1,n_dim/3
!       write(2,*) "A",bar_cm(1+(j-1)*3:3+(j-1)*3)-bar_vec(1+(j-1)*3:3+(j-1)*3)
!      enddo
!      close(2)
!      write(fn,'(i5.5)') iter
!      filename="possaddle.right."//fn//".xyz"
!      open(unit=2,file=trim(filename))
!      write(2,*) n_dim/3
!      write(2,*) etot_all(2)
!      do j=1,n_dim/3
!       write(2,*) "A",bar_cm(1+(j-1)*3:3+(j-1)*3)+bar_vec(1+(j-1)*3:3+(j-1)*3)
!      enddo
!      close(2)
endif
 
    if(fnrm.lt.bar_tol) then
       if(bar_contract.and. .not. contracting) then
       contracting=.true.
       contractcount=0
       bar_tol=contr_bar_tol
       diffbar=dbar-contr_dbar
       !if(iproc==0)   write(*,'(a)') "SAD: Starting contraction"
       if(iproc==0) call yaml_scalar("SAD: Starting contraction")
       else
         if(lo_acc) then
            lo_acc=.false.
            if(iproc==0)   write(*,'(a)') "SAD: Switching to higher accuracy"
         else
            !if(iproc==0)   write(*,'(a)') "SAD: Saddle point search converged"
            !call yaml_sequence(advance='no')
            if(iproc==0)   call yaml_scalar("BARSAD: Saddle point search converged")
            goto 1002
         endif
       endif
    endif

!    if(etot_all(1).gt.etot_all(2).or.etot_all(3).gt.etot_all(2)) stop "Energy
!    ordering bad!"
enddo !end of loop over saddle optimization

if(iproc==0)   write(*,'(a)') "SAD: Saddle point search did not converge in given number of iterations"
1002 continue
call yaml_sequence_close()
if(iproc==0) then
      filename="possaddle_cur"
              call write_atomic_file(trim(filename),etot_all(2),bar_cm(:),atoms,"Approximate saddle",fbar_tot(:,2))
      filename="possaddle.left_cur"
              call write_atomic_file(trim(filename),etot_all(1),bar_cm(:)-bar_vec(:),atoms,"Left bar end",fbar_tot(:,1))
      filename="possaddle.right_cur"
              call write_atomic_file(trim(filename),etot_all(3),bar_cm(:)+bar_vec(:),atoms,"Right bar end",fbar_tot(:,3))
endif
etot=etot_all(2)
contains

subroutine move_bar(n_dim,bar_vec,bar_cm,fbar_para,fbar_perp,fbar_tot,etot,alpha)
!This subroutine will move the bar according to the forces at both endpoints, 
!including a shift along the bar given by the parallel forces of the CM
use bar_saddle_params, only:dbar
implicit none
integer:: n_dim
real(8),dimension(n_dim):: bar_vec,bar_cm,bar_vec_unit,rot_mean
real(8),dimension(n_dim,3)::fbar_para,fbar_perp,fbar_tot,new_bar
real(8):: dbar_cur,alpha,norm,etot(3)
!compute initial bar-length
dbar_cur=sqrt(dot_product(bar_vec,bar_vec))
!Rotation mean
rot_mean=0.5d0*(fbar_perp(:,1)+fbar_perp(:,3))
!move first endpoint
new_bar(:,1)=(bar_cm(:)-bar_vec(:))+alpha*(1.d0*fbar_perp(:,1)-rot_mean(:)-3.d0*fbar_para(:,2)+fbar_tot(:,2))
!move second endpoint
new_bar(:,3)=(bar_cm(:)+bar_vec(:))+alpha*(1.d0*fbar_perp(:,3)-rot_mean(:)-3.d0*fbar_para(:,2)+fbar_tot(:,2))
!compute new bar center
new_bar(:,2)=(new_bar(:,1)+new_bar(:,3))*0.5d0
bar_vec(:)=(new_bar(:,3)-new_bar(:,1))*0.5d0
!Rescale
norm=sqrt(dot_product(bar_vec,bar_vec))
bar_vec(:)=bar_vec(:)/norm*dbar_cur
bar_vec=bar_vec/dbar_cur*dbar
!Finalize
bar_cm=new_bar(:,2)
end subroutine

subroutine move_bar_energy(n_dim,bar_vec,bar_cm,fbar_para,fbar_perp,fbar_tot,etot_all,alpha,alpha_min,alpha_max,iter,skip)
!This subroutine will move the bar according to the forces at both endpoints, 
!including a shift along the bar given by the parallel forces of the CM
use bar_saddle_params, only:dbar
implicit none
integer:: n_dim,iter
real(8),dimension(n_dim):: bar_vec,bar_cm,bar_vec_unit,grad_trans,rot_mean
real(8),dimension(n_dim,3)::fbar_para,fbar_perp,fbar_tot,new_bar
real(8):: dbar_cur,alpha,alpha_min,alpha_max,norm,etot_all(3),anoise
logical:: skip
!Feedback variables
real(8),allocatable,save:: fbar_tot_old(:,:),fbar_perp_old(:,:),fbar_para_old(:,:),bar_cm_old(:),bar_vec_old(:),grad_trans_old(:)
real(8),save:: etot_all_old(3)

anoise=5.d-3/real(n_dim,8)

!Initiallize if first call
if(iter.le.1) then
  if(allocated(fbar_tot_old)) deallocate(fbar_tot_old,bar_cm_old,fbar_perp_old,fbar_para_old,bar_vec_old,grad_trans_old)
  if(.not.allocated(fbar_tot_old)) allocate(fbar_tot_old(n_dim,3),fbar_perp_old(n_dim,3),fbar_para_old(n_dim,3),bar_cm_old(n_dim),bar_vec_old(n_dim),grad_trans_old(n_dim))
  fbar_tot_old=fbar_tot;bar_cm_old=bar_cm;fbar_perp_old=fbar_perp;fbar_para_old=fbar_para;bar_vec_old=bar_vec
  grad_trans_old=-3.d0*fbar_para(:,2)+fbar_tot(:,2)
!Eliminate torque and stuff
  call elim_moment_unit(n_dim/3,grad_trans_old)
  call elim_torque_unit(n_dim/3,bar_cm,grad_trans_old)
endif

!compute initial bar-length
dbar_cur=sqrt(dot_product(bar_vec,bar_vec))

!Construct translational vector
grad_trans=-3.d0*fbar_para(:,2)+fbar_tot(:,2)
!Eliminate torque and stuff
call elim_moment_unit(n_dim/3,grad_trans)
call elim_torque_unit(n_dim/3,bar_cm,grad_trans)

!Energy feedback
   if(etot_all(2).gt.etot_all_old(2)+anoise.and..not.skip) then
     if(alpha.gt.alpha_min) then
       alpha=max(alpha*0.7d0,alpha_min)
       bar_vec=bar_vec_old;bar_cm=bar_cm_old
       fbar_para=fbar_para_old;fbar_perp=fbar_perp_old;fbar_tot=fbar_tot_old
       bar_vec=bar_vec_old;bar_cm=bar_cm_old
       etot_all=etot_all_old
       !if(iproc==0) write(*,'(a,es15.7)') "SAD: Reducing step size to", alpha
       if(iproc==0) call yaml_map('Reducing step size to',alpha,fmt='(es15.7)')
     else
       fbar_para_old=fbar_para;fbar_perp_old=fbar_perp;fbar_tot_old=fbar_tot
       bar_vec_old=bar_vec;bar_cm_old=bar_cm
       etot_all_old=etot_all
     endif
   else 
     alpha=min(alpha*1.05d0,alpha_max)
     fbar_para_old=fbar_para;fbar_perp_old=fbar_perp;fbar_tot_old=fbar_tot
     bar_vec_old=bar_vec;bar_cm_old=bar_cm
     etot_all_old=etot_all
   endif

!Rotation mean
rot_mean=0.5d0*(fbar_perp(:,1)+fbar_perp(:,3))
!move first endpoint
new_bar(:,1)=(bar_cm(:)-bar_vec(:))+alpha*(1.d0*fbar_perp(:,1)-rot_mean(:)+grad_trans)
!move second endpoint
new_bar(:,3)=(bar_cm(:)+bar_vec(:))+alpha*(1.d0*fbar_perp(:,3)-rot_mean(:)+grad_trans)
!compute new bar center
new_bar(:,2)=(new_bar(:,1)+new_bar(:,3))*0.5d0
bar_vec(:)=(new_bar(:,3)-new_bar(:,1))*0.5d0
!Rescale
norm=sqrt(dot_product(bar_vec,bar_vec))
bar_vec(:)=bar_vec(:)/norm*dbar_cur
bar_vec=bar_vec/dbar_cur*dbar
!Finalize
bar_cm=new_bar(:,2)
end subroutine
subroutine move_bar_force(n_dim,bar_vec,bar_cm,fbar_para,fbar_perp,fbar_tot,alpha,alpha_min,alpha_max,iter,skip)
!This subroutine will move the bar according to the forces at both endpoints, 
!including a shift along the bar given by the parallel forces of the CM
use bar_saddle_params, only:dbar
implicit none
integer:: n_dim,iter
real(8),dimension(n_dim)::bar_vec,bar_cm,bar_vec_unit,grad_trans,grad_rot,rot_mean
real(8),dimension(n_dim,3)::fbar_para,fbar_perp,fbar_tot,new_bar
real(8):: dbar_cur,alpha,alpha_min,alpha_max,norm
logical:: skip
!Feedback variables
real(8),allocatable,save:: fbar_tot_old(:,:),fbar_perp_old(:,:),fbar_para_old(:,:),bar_cm_old(:),bar_vec_old(:),grad_trans_old(:)
!Initiallize if first call
if(iter.le.1) then
  if(allocated(fbar_tot_old)) deallocate(fbar_tot_old,bar_cm_old,fbar_perp_old,fbar_para_old,bar_vec_old,grad_trans_old)
  if(.not.allocated(fbar_tot_old)) allocate(fbar_tot_old(n_dim,3),fbar_perp_old(n_dim,3),fbar_para_old(n_dim,3),bar_cm_old(n_dim),bar_vec_old(n_dim),grad_trans_old(n_dim))
  fbar_tot_old=fbar_tot;bar_cm_old=bar_cm;fbar_perp_old=fbar_perp;fbar_para_old=fbar_para;bar_vec_old=bar_vec
  grad_trans_old=-3.d0*fbar_para(:,2)+fbar_tot(:,2)
!Eliminate torque and stuff
  call elim_moment_unit(n_dim/3,grad_trans_old)
  call elim_torque_unit(n_dim/3,bar_cm,grad_trans_old)
endif

!compute initial bar-length
dbar_cur=sqrt(dot_product(bar_vec,bar_vec))

!Construct translational vector
grad_trans=-3.d0*fbar_para(:,2)+fbar_tot(:,2)
!grad_trans=-3.d0/2.d0*(fbar_para(:,1)+fbar_para(:,3))+0.5*(fbar_tot(:,1)+fbar_tot(:,3))
!Eliminate torque and stuff
call elim_moment_unit(n_dim/3,grad_trans)
call elim_torque_unit(n_dim/3,bar_cm,grad_trans)


!Gradient feedback
   cosinus=dot_product(fbar_tot(:,2),fbar_tot_old(:,2))/&
   &sqrt(dot_product(fbar_tot_old(:,2),fbar_tot_old(:,2))*dot_product(fbar_tot(:,2),fbar_tot(:,2)))
!   cosinus=dot_product(grad_trans(:),grad_trans_old(:))/&
!   &sqrt(dot_product(grad_trans_old(:),grad_trans_old(:))*dot_product(grad_trans(:),grad_trans(:)))
   if(skip) then
          cosinus=1.d0
   endif
   if(cosinus.lt.0.4d0) then
     if(alpha.gt.alpha_min) then
       alpha=max(alpha*0.7d0,alpha_min)
       bar_vec=bar_vec_old;bar_cm=bar_cm_old
       fbar_para=fbar_para_old;fbar_perp=fbar_perp_old;fbar_tot=fbar_tot_old;grad_trans=grad_trans_old
       bar_vec=bar_vec_old;bar_cm=bar_cm_old
       !if(iproc==0) write(*,'(a,es15.7)') "SAD: Reducing step size to", alpha
       if(iproc==0) call yaml_map('Reducing step size to',alpha,fmt='(es15.7)')
     else
       fbar_para_old=fbar_para;fbar_perp_old=fbar_perp;fbar_tot_old=fbar_tot;grad_trans_old=grad_trans
       bar_vec_old=bar_vec;bar_cm_old=bar_cm
!       etot_all_old=etot_all
     endif
   else 
     alpha=min(alpha*1.05d0,alpha_max)
     fbar_para_old=fbar_para;fbar_perp_old=fbar_perp;fbar_tot_old=fbar_tot;grad_trans_old=grad_trans
     bar_vec_old=bar_vec;bar_cm_old=bar_cm
!     etot_all_old=etot_all
   endif

!Rotation mean
rot_mean=0.5d0*(fbar_perp(:,1)+fbar_perp(:,3))

!move first endpoint
grad_rot=(1.d0*fbar_perp(:,1)-rot_mean(:))
!Eliminate torque and stuff
call elim_moment_unit(n_dim/3,grad_rot)
!call elim_torque_unit(n_dim/3,bar_cm,grad_rot)
call elim_torque_unit(n_dim/3,bar_cm(:)-bar_vec(:),grad_rot)
new_bar(:,1)=(bar_cm(:)-bar_vec(:))+alpha*(grad_rot+grad_trans)
!move second endpoint
grad_rot=(1.d0*fbar_perp(:,3)-rot_mean(:))
!Eliminate torque and stuff
call elim_moment_unit(n_dim/3,grad_rot)
!call elim_torque_unit(n_dim/3,bar_cm,grad_rot)
call elim_torque_unit(n_dim/3,bar_cm(:)+bar_vec(:),grad_rot)
new_bar(:,3)=(bar_cm(:)+bar_vec(:))+alpha*(grad_rot+grad_trans)
!compute new bar center
new_bar(:,2)=(new_bar(:,1)+new_bar(:,3))*0.5d0
bar_vec(:)=(new_bar(:,3)-new_bar(:,1))*0.5d0
!Rescale
norm=sqrt(dot_product(bar_vec,bar_vec))
bar_vec(:)=bar_vec(:)/norm*dbar_cur
bar_vec=bar_vec/dbar_cur*dbar
!Finalize
bar_cm=new_bar(:,2)
end subroutine


subroutine move_bar_bfgs(iproc,n_dim,bar_vec,bar_cm,fbar_para,fbar_perp,fbar_tot,etot,alpha,iter,nhist)
!This subroutine will move the bar according to the forces at both endpoints, 
!including a shift along the bar given by the parallel forces of the CM
use bar_saddle_params, only:dbar
implicit none
integer:: n_dim,i,iproc
real(8),dimension(n_dim):: bar_vec,bar_cm,bar_vec_unit,rot_mean
real(8),dimension(n_dim,3)::fbar_para,fbar_perp,fbar_tot,new_bar
real(8):: dbar_cur,alpha,norm,etot(3)
!BFGS variables
real(8),allocatable,save:: bs(:,:),by(:,:),bg(:),bq(:),bz(:),brho(:),ba(:),bb(:),bh(:),grad_prev(:)
integer::nhist,iter,bound
real(8):: grad_rot(n_dim),grad(n_dim),y(n_dim),ss(n_dim),ys,yy
!compute initial bar-length
dbar_cur=sqrt(dot_product(bar_vec,bar_vec))

!Initiallize if first call
if(iter.le.1) then
  write(*,*) "Allocating", iproc, allocated(bs)
  if(allocated(bs)) deallocate(bs,by,bg,bq,bz,brho,ba,bb,bh,grad_prev)
  if(.not.allocated(bs))  allocate(bs(n_dim,nhist),by(n_dim,nhist),bg(n_dim),bq(n_dim),bz(n_dim),brho(nhist),ba(nhist),bb(nhist),bh(n_dim),grad_prev(n_dim))
  bs=0.d0;by=0.d0;bg=0.d0;bq=0.d0;bz=0.d0;brho=0.d0;ba=0.d0;bb=0.d0;grad_prev=0.d0
endif

!Initillize boundaries
if(iter.le.nhist) then
  bound=iter
else
  bound=nhist
endif

!Set gradient to q
!grad=3.d0*fbar_para(:,2)-fbar_tot(:,2)
grad=1.d0*fbar_para(:,2)-fbar_tot(:,2)

!Eliminate torque and stuff
call elim_moment_unit(n_dim/3,grad)
call elim_torque_unit(n_dim/3,bar_cm,grad)

bq=grad

!Update the y
y=grad-grad_prev
if(iter.gt.nhist+1) then
  by(:,1:nhist-1)=by(:,2:nhist)
  by(:,nhist)=y
elseif(iter.gt.1) then
  by(:,iter-1)=y
endif
yy=dot_product(y,y)

!Update rho
if(iter.gt.nhist+1) then
  brho(1:nhist-1)=brho(2:nhist)
  ys=dot_product(by(:,nhist),bs(:,nhist))
  brho(nhist)=1.d0/ys
elseif(iter.gt.1) then
  ys=dot_product(by(:,iter-1),bs(:,iter-1))
  brho(iter-1)=1.d0/ys
endif

!Diagonal elements of the initial hessian
if(iter.le.1) then
  bh=alpha
else
  bh=ys/yy
endif


if(iproc==0)write(*,*) "BRHO", iproc,brho
if(iproc==0)write(*,*) "BY", iproc,by
if(iproc==0)write(*,*) "BS", iproc,bs
if(iproc==0)write(*,*) "YS/YY",YS/YY 

!Do the first loop
do i=bound-1,1,-1
   ba(i)=brho(i)*dot_product(bs(:,i),bq(:))
   bq=bq-ba(i)*by(:,i)
enddo
bq=bh*bq
!Do the second loop
do i=1,bound-1
   bb(i)=brho(i)*dot_product(by(:,i),bq(:))
   bq=bq+bs(:,i)*(ba(i)-bb(i))
enddo   

!Check that the step is not too large
norm=dot_product(bq,bq)
if(norm.gt.0.2d0**2) then
  bq=bq/sqrt(norm)*0.2d0
  if(iproc==0) write(*,*) "NORM",iter,norm
endif

!Update s
if(iter.gt.nhist) then
  bs(:,1:nhist-1)=bs(:,2:nhist);bs(:,nhist)=-bq(:)
else
  bs(:,iter)=-bq(:)
endif
grad_prev=grad



!Update the bar translationally
new_bar(:,1)=(bar_cm(:)-bar_vec(:))-bq(:)
new_bar(:,3)=(bar_cm(:)+bar_vec(:))-bq(:)


!ADD ROTATIONS
!Rotation mean
rot_mean=0.5d0*(fbar_perp(:,1)+fbar_perp(:,3))

!move first endpoint
grad_rot=alpha*(1.d0*fbar_perp(:,1)-rot_mean(:))

!Eliminate torque and stuff
call elim_moment_unit(n_dim/3,grad_rot)
!call elim_torque_unit(n_dim/3,bar_cm,grad_rot)
call elim_torque_unit(n_dim/3,new_bar(:,1),grad_rot)
new_bar(:,1)=new_bar(:,1)+grad_rot!alpha*(1.d0*fbar_perp(:,1)-rot_mean(:))!-3.d0*fbar_para(:,2)+fbar_tot(:,2))

!move second endpoint
grad_rot=alpha*(1.d0*fbar_perp(:,3)-rot_mean(:))
!Eliminate torque and stuff
call elim_moment_unit(n_dim/3,grad_rot)
!call elim_torque_unit(n_dim/3,bar_cm,grad_rot)
call elim_torque_unit(n_dim/3,new_bar(:,3),grad_rot)
new_bar(:,3)=new_bar(:,3)+grad_rot!alpha*(1.d0*fbar_perp(:,3)-rot_mean(:))!-3.d0*fbar_para(:,2)+fbar_tot(:,2))

!compute new bar center
new_bar(:,2)=(new_bar(:,1)+new_bar(:,3))*0.5d0
bar_vec(:)=(new_bar(:,3)-new_bar(:,1))*0.5d0
!Rescale
norm=sqrt(dot_product(bar_vec,bar_vec))
bar_vec(:)=bar_vec(:)/norm*dbar_cur
bar_vec=bar_vec/dbar_cur*dbar
!Finalize
if(iproc==0) write(*,*) "DIST",dot_product(bar_cm(:)-bq(:)-new_bar(:,2),bar_cm(:)-bq(:)-new_bar(:,2))
bar_cm=new_bar(:,2)
end subroutine

!*****************************************************************************************
end subroutine find_saddle
subroutine write_atomic_file(filename,etot,pos,atoms,comment,forces)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(in) :: atoms
    real(8)::pos(3,atoms%nat),xyz(3),forces(3,atoms%nat),etot
    integer::ip,i
    character(len=10) :: name
    character(len=2) :: symbol
    character(*)::filename
    character(*)::comment
    !-------------------------------------------------------------------------------------------
    open(unit=1389,file=filename,status='replace')
    write(1389,*) atoms%nat
    write(1389,*) trim(comment)
        do i=1,atoms%nat
                name=trim(atoms%sat(i))
                if (name(3:3)=='_') then
                   symbol=name(1:2)
                   !suffix=name(4:6)
                else if (name(2:2)=='_') then
                   symbol=name(1:1)
                   !suffix=name(3:5)
                else
                   symbol=name(1:2)
                   !suffix=' '
                end if
                if(atoms%units=='angstroemd0' .or. atoms%units=='angstroem' .or. atoms%units=='angstrom') then
                    xyz(1:3)=pos(:,i)*0.5291772108d0 !non-BigDFT
                else
                    xyz(1:3)=pos(:,i)
                endif
                write(1389,'(a,1x,3e24.15)') symbol,xyz(1),xyz(2),xyz(3)
        enddo
    close(1389)
    !-------------------------------------------------------------------------------------------
end subroutine write_atomic_file

subroutine forcebar_decomp(n_dim,bar_vec,bar_cm,fbar_para,fbar_perp,fbar_tot,etot,xmax,emax,bar_max,nproc,iproc,atoms,parini,iter)
!This subroutine will compute the energy and forces at the bar_cm, and decompose the
!force on the CM in a part parallel to the bar and perpendicular to the bar.
!all arrays: etot(1)=bar_cm-bar_vec,etot(2)=bar_cm,etot(3)=bar_cm+bar_vec
use mod_parini, only: typ_parini
use mod_atoms, only: typ_atoms
implicit none
!Bigdft Variables
integer, intent(in) :: nproc,iproc
type(typ_parini):: parini
type(typ_atoms) :: atoms
integer:: n_dim,i,iter
integer:: infocode
real(8),dimension(n_dim):: bar_vec,bar_cm,bar_vec_unit,bar_max
real(8),dimension(n_dim,3)::fbar_para,fbar_perp,fbar_tot
real(8):: etot(3),fnoise
real(8):: norm,d1,d2,x,dx,xmax,emax
!do i=1,3
!    !Compute force at barends
!    call func2d(n_dim,bar_cm+real(i-2,8)*bar_vec,fbar_tot(:,i),etot(i))
!!    call lenjon(n_dim/3,bar_cm+real(i-2,8)*bar_vec,fbar_tot(:,i),etot(i))
!    !Compute the force in the direction of the bar
!    norm=sqrt(dot_product(bar_vec,bar_vec))
!    bar_vec_unit=bar_vec/norm
!    fbar_para(:,i)=dot_product(fbar_tot(:,i),bar_vec_unit)*bar_vec_unit
!    !Compute the force perpendicular to the bar
!    fbar_perp(:,i)=fbar_tot(:,i)-fbar_para(:,i)
!enddo

!if(iproc==0) write(*,*)"SAD: Forcebar, iteration", iter
do i=1,3,2
    !Compute force at barends
!    call func2d(n_dim,bar_cm+real(i-2,8)*bar_vec,fbar_tot(:,i),etot(i))
!    call lenjon(n_dim/3,bar_cm+real(i-2,8)*bar_vec,fbar_tot(:,i),etot(i))
    call call_bigdft(nproc,iproc,atoms,bar_cm+real(i-2,8)*bar_vec,       etot(i),fbar_tot(:,i),       fnoise,infocode,parini)
    !write(*,*) i,etot(i)
    !Compute the force in the direction of the bar
    norm=sqrt(dot_product(bar_vec,bar_vec))
    bar_vec_unit=bar_vec/norm
    fbar_para(:,i)=dot_product(fbar_tot(:,i),bar_vec_unit)*bar_vec_unit
    !Compute the force perpendicular to the bar
    fbar_perp(:,i)=fbar_tot(:,i)-fbar_para(:,i)
enddo
!Interpolate energy and forces at the center of the bar
!The perpendicular force
fbar_perp(:,2)=0.5d0*(fbar_perp(:,1)+fbar_perp(:,3))
!Now the energy
!Now the parallel force
d1=dot_product(-fbar_tot(:,1),2.d0*bar_vec)
d2=dot_product(-fbar_tot(:,3),2.d0*bar_vec)
x=0.5d0
call cub_interpol(etot(1),etot(3),d1,d2,x,etot(2),dx,xmax,emax)
fbar_para(:,2)=-dx/(2.d0*norm)*bar_vec_unit
!NEW+++++++++++++
bar_max(:)=bar_cm+xmax*2.d0*norm
!NEW+++++++++++++

!write(*,*)
!-dx/(2.d0*norm),0.5d0*(dot_product(fbar_tot(:,1),bar_vec_unit)+dot_product(fbar_tot(:,3),bar_vec_unit))

!The overall force
!fbar_tot(:,2)=0.5d0*(fbar_tot(:,1)+fbar_tot(:,3))
fbar_tot(:,2)=fbar_para(:,2)+fbar_perp(:,2)

end subroutine
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine elim_moment_unit(nat,fxyz)
implicit real*8 (a-h,o-z)
dimension fxyz(3,nat)

sx=0.d0 ; sy=0.d0 ; sz=0.d0
do iat=1,nat
    sx=sx+fxyz(1,iat)
    sy=sy+fxyz(2,iat)
    sz=sz+fxyz(3,iat)
enddo
    sx=sx/nat ; sy=sy/nat ; sz=sz/nat
do iat=1,nat
    fxyz(1,iat)=fxyz(1,iat)-sx
    fxyz(2,iat)=fxyz(2,iat)-sy
    fxyz(3,iat)=fxyz(3,iat)-sz
enddo

return

end subroutine elim_moment_unit

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine elim_torque_unit(nat,rxyz,fxyz)
implicit real*8 (a-h,o-z)
dimension rxyz(3,nat),fxyz(3,nat),t(3)

! center of mass
cmx=0.d0 ; cmy=0.d0 ; cmz=0.d0
do iat=1,nat
    cmx=cmx+rxyz(1,iat)
    cmy=cmy+rxyz(2,iat)
    cmz=cmz+rxyz(3,iat)
enddo
cmx=cmx/nat ; cmy=cmy/nat ; cmz=cmz/nat

do it=1,100

    ! torque and radii in planes
    t(1)=0.d0 ; t(2)=0.d0 ; t(3)=0.d0
    sx=0.d0 ; sy=0.d0 ; sz=0.d0
    do iat=1,nat
        t(1)=t(1)+(rxyz(2,iat)-cmy)*fxyz(3,iat)-(rxyz(3,iat)-cmz)*fxyz(2,iat)
        t(2)=t(2)+(rxyz(3,iat)-cmz)*fxyz(1,iat)-(rxyz(1,iat)-cmx)*fxyz(3,iat)
        t(3)=t(3)+(rxyz(1,iat)-cmx)*fxyz(2,iat)-(rxyz(2,iat)-cmy)*fxyz(1,iat)
        sx=sx+(rxyz(1,iat)-cmx)**2
        sy=sy+(rxyz(2,iat)-cmy)**2
        sz=sz+(rxyz(3,iat)-cmz)**2
    enddo

    if (t(1)**2+t(2)**2+t(3)**2.lt.1.d-20) return

    ii=0
    tmax=0.d0
    do i=1,3
        if (t(i)**2.ge.tmax**2) then
            ii=i
            tmax=t(i)
        endif
    enddo

    !         write(*,'(i4,3(1pe11.3))') ii,t

    ! modify forces
    if (ii.eq.1) then
        cx=t(1)/(sz+sy)
        do iat=1,nat
            fxyz(2,iat)=fxyz(2,iat)+cx*(rxyz(3,iat)-cmz)
            fxyz(3,iat)=fxyz(3,iat)-cx*(rxyz(2,iat)-cmy)
        enddo
    else if(ii.eq.2) then
        cy=t(2)/(sz+sx)
        do iat=1,nat
            fxyz(1,iat)=fxyz(1,iat)-cy*(rxyz(3,iat)-cmz)
            fxyz(3,iat)=fxyz(3,iat)+cy*(rxyz(1,iat)-cmx)
        enddo
    else if(ii.eq.3) then
        cz=t(3)/(sy+sx)
        do iat=1,nat
            fxyz(1,iat)=fxyz(1,iat)+cz*(rxyz(2,iat)-cmy)
            fxyz(2,iat)=fxyz(2,iat)-cz*(rxyz(1,iat)-cmx)
        enddo
    else
        !write(*,*) t,tmax
        !stop 'wrong ii'
    endif

enddo
!write(*,'(a,3(1pe11.3))') 'WARNING REMAINING TORQUE',t

return
end subroutine elim_torque_unit

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine fnrmandforcemax(fxyz,fnrm,fmax,nat)
    implicit none
    real(8):: fxyz(3,nat),fnrm,fmax
    integer:: nat
    call calmaxforcecomponent(3*nat,fxyz,fmax)
    call calnorm(3*nat,fxyz,fnrm)
end subroutine

         
subroutine cub_interpol(e1,e2,d1,d2,x,ex,dx,xmax,emax)
!This routine will perform a cubic interpolation of a function given at 2 points
!with
!its value and derivative, in the interval of [0,1]
implicit none
real(8)::e1,e2,d1,d2,x,ex,dx
real(8):: a,b,c,d,xmax,discr,x12(2),emax
d=e1
c=d1
b=-(2.d0*c+3.d0*d+d2)+3.d0*e2
a=-(b+c+d)+e2
!This might be wrong, but mathematica gives me: a=2.d0*d-2.d0*e2+d1+d2

ex=a*x**3+b*x**2+c*x+d
dx=3.d0*a*x**2+2.d0*b*x+c
!Solve the quadratic equation to get the maximum along the line
discr=4.d0*b**2-4.d0*3.d0*a*c
if(discr.lt.0.d0) then
  x12=1.d10
else
  x12(1)=(-2.d0*b+sqrt(discr))/(2.d0*3.d0*a)
  x12(2)=(-2.d0*b-sqrt(discr))/(2.d0*3.d0*a)
endif
if(abs(x12(1)-x).lt.abs(x12(2)-x)) then
!NEW++++++++++
  emax=a*x12(1)**3+b*x12(1)**2+c*x12(1)+d
!NEW++++++++++
  xmax=x12(1)-x
else
!NEW++++++++++
  emax=a*x12(2)**3+b*x12(2)**2+c*x12(2)+d
!NEW++++++++++
  xmax=x12(2)-x
endif
end subroutine


!**********************************************************
!taken from spline_splint.f90
!**********************************************************

      subroutine spline_drv(x,y,n,yp1,ypn,y2)
      implicit none
      INTEGER n,NMAX
      REAL(8):: yp1,ypn,x(n),y(n),y2(n),xt(n),yt(n),ytp1,ytpn
      PARAMETER (NMAX=500)
      INTEGER i,k
      REAL(8):: p,qn,sig,un,u(NMAX)
      if(x(1).eq.x(n)) then
        y2=0.d0
      elseif(x(1).gt.x(2)) then
        xt=-x
        yt=y
        ytp1=-yp1
        ytpn=-ypn
        call spline_barsaddle(xt,yt,n,ytp1,ytpn,y2)
      else
        call spline_barsaddle(x,y,n,yp1,ypn,y2)
      endif
      end subroutine
      
      subroutine splint_drv(xa,ya,y2a,n,x,y,dy)
      implicit none
      INTEGER n
      REAL(8):: x,y,xa(n),y2a(n),ya(n),xat(n),yat(n),xt,dy
      INTEGER k,khi,klo,i
      REAL(8):: a,b,h
      if(xa(1).eq.xa(n)) then
        y=ya(1)
        dy=0.d0
      elseif(xa(1).gt.xa(2)) then
        xat=-xa
        yat=ya
        xt=-x
        call splint(xat,yat,y2a,n,xt,y,dy)
        dy=-dy
      else
        call splint(xa,ya,y2a,n,x,y,dy)
      endif
      end subroutine

      subroutine spline_barsaddle(x,y,n,yp1,ypn,y2)  
!From numerical recipes
!Given arrays x[1..n] and y[1..n] containing a tabulated function, i.e., yi= f(xi), with
!x1 < x2 < : : :  < xN , and given values yp1 and ypn for the first derivative of the interpolating
!function at points 1 and n, respectively, this routine returns an array y2[1..n] that contains
!the second derivatives of the interpolating function at the tabulated points xi. If yp1 and/or
!ypn are equal to 10^30 or larger, the routine is signaled to set the corresponding boundary
!condition for a natural spline, with zero second derivative on that boundary.
      implicit none
      INTEGER n,NMAX  
      REAL(8):: yp1,ypn,x(n),y(n),y2(n)  
      PARAMETER (NMAX=500)  
      INTEGER i,k  
      REAL(8):: p,qn,sig,un,u(NMAX)  
      if (yp1.gt..99d30) then  
        y2(1)=0.d0  
        u(1)=0.d0  
      else  
        y2(1)=-0.5d0  
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)  
      endif  
      do 11 i=2,n-1  
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))  
        p=sig*y2(i-1)+2.d0  
        y2(i)=(sig-1.d0)/p  
        u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p  
11    continue  
      if (ypn.gt..99d30) then  
        qn=0.d0  
        un=0.d0  
      else  
        qn=0.5d0  
        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))  
      endif  
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)  
      do 12 k=n-1,1,-1  
        y2(k)=y2(k)*y2(k+1)+u(k)  
12    continue  
      return  
      END  

      subroutine splint(xa,ya,y2a,n,x,y,dy) 
!Given the arrays xa[1..n] and ya[1..n], which tabulate a function (with the xai's in order),
!and given the array y2a[1..n], which is the output from spline above, and given a value of
!x, this routine returns a cubic-spline interpolated value y 
      implicit none
      INTEGER n  
      REAL(8):: x,y,xa(n),y2a(n),ya(n),dy  
      INTEGER k,khi,klo  
      REAL(8):: a,b,h,hy  
      klo=1  
      khi=n  
1     if (khi-klo.gt.1) then  
        k=(khi+klo)/2  
        if(xa(k).gt.x)then  
          khi=k  
        else  
          klo=k  
        endif  
      goto 1  
      endif  
      h=xa(khi)-xa(klo)  
      if (h.eq.0.d0) stop 'bad xa input in splint'  
      a=(xa(khi)-x)/h  
      b=(x-xa(klo))/h  
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d0  

!Lets additionally compute the derivative at point x with respect to x
      hy=ya(khi)-ya(klo)
      dy=hy/h+(-(3.d0*a**2-1.d0)*y2a(klo)+(3.d0*b**2-1.d0)*y2a(khi))/6.d0*h
      return  
      END  
 
!**********************************************************
!taken from brent.f90
!**********************************************************

      function brent_barsaddle(ax,bx,cx,f,tol,xmin,nbrent,nproc,iproc,atoms,parini,ncount_bigdft)  
!Given a function f, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is
!between ax and cx, and f(bx) is less than both f(ax) and f(cx)), this routine isolates
!the minimum to a fractional precision of about tol using Brent's method. The abscissa of
!the minimum is returned as xmin, and the minimum function value is returned as brent, the
!returned function value.
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use bar_saddle_params, only: maxit
implicit none
!Bigdft Variables
    integer, intent(in) :: nproc,iproc
    type(typ_parini), intent(in):: parini
    type(typ_atoms):: atoms

      INTEGER ITMAX,nbrent,ncount_bigdft
      REAL(8):: brent_barsaddle,ax,bx,cx,tol,xmin,CGOLD,ZEPS  
      REAL(8),EXTERNAL:: f  
      PARAMETER (CGOLD=.3819660,ZEPS=1.0d-10)  
      INTEGER:: iter  
      REAL(8):: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm  
      nbrent=0
      ITMAX=300!maxit-ncount_bigdft
      a=min(ax,cx)  
      b=max(ax,cx)  
      v=bx  
      w=v  
      x=v  
      e=0.d0  
      !fx=f(x,nproc,iproc,atoms,rst,inputs,ncount_bigdft)  
      fx=f(x,nproc,iproc,atoms,parini,ncount_bigdft)
      nbrent=nbrent+1
      fv=fx  
      fw=fx  
      do 11 iter=1,ITMAX  
        xm=0.5d0*(a+b)  
        tol1=tol*abs(x)+ZEPS  
        tol2=2.d0*tol1  
        if(abs(x-xm).le.(tol2-.5d0*(b-a))) goto 3  
        if(abs(e).gt.tol1) then  
          r=(x-w)*(fx-fv)  
          q=(x-v)*(fx-fw)  
          p=(x-v)*q-(x-w)*r  
          q=2.d0*(q-r)  
          if(q.gt.0.d0) p=-p  
          q=abs(q)  
          etemp=e  
          e=d  
          if(abs(p).ge.abs(.5d0*q*etemp).or.p.le.q*(a-x).or.p.ge.q*(b-x)) goto 1  
          d=p/q  
          u=x+d  
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)  
          goto 2  
        endif  
1       if(x.ge.xm) then  
          e=a-x  
        else  
          e=b-x  
        endif  
        d=CGOLD*e  
2       if(abs(d).ge.tol1) then  
          u=x+d  
        else  
          u=x+sign(tol1,d)  
        endif  
!        fu=f(u,nproc,iproc,atoms,rst,inputs,ncount_bigdft)  
        fu=f(u,nproc,iproc,atoms,parini,ncount_bigdft)
        nbrent=nbrent+1
        if(fu.le.fx) then  
          if(u.ge.x) then  
            a=x  
          else  
            b=x  
          endif  
          v=w  
          fv=fw  
          w=x  
          fw=fx  
          x=u  
          fx=fu  
        else  
          if(u.lt.x) then  
            a=u  
          else  
            b=u  
          endif  
          if(fu.le.fw .or. w.eq.x) then  
            v=w  
            fv=fw  
            w=u  
            fw=fu  
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then  
            v=u  
            fv=fu  
          endif  
        endif  
11    continue  
      stop 'brent exceed maximum iterations'  
3     xmin=x  
      brent_barsaddle=fx  
      return  
      END  
 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function time_energy(tt,nproc,iproc,atoms,parini,ncount_bigdft)
!Computes the rxyz according to a splined pathway, previously defined by t,coord
!and sp_param,
!at the time parameter tt(0:1), and outputs the energy value at the interpolated
!point. 
use spline_params
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
implicit none
!Bigdft Variables
    integer, intent(in) :: nproc,iproc
    type(typ_parini), intent(in):: parini
    type(typ_atoms):: atoms
    integer, intent(inout) :: ncount_bigdft
real(8) :: strten(6), fnoise
integer  :: infocode 
        integer:: i
        real(8):: time_energy,tt
        real(8):: rxyz(n_dim),fxyz(n_dim),etot,yd
        do i=1,n_dim
           call splint_drv(t(:),coord(i,:),sp_param(i,:),n_anc,tt,rxyz(i),yd)
        enddo
!        call func2d(n_dim,rxyz,fxyz,etot)
!        call lenjon(n_dim/3,rxyz,fxyz,etot)
!        call call_bigdft(nproc,iproc,atoms,rxyz,inputs,etot,fxyz,strten,fnoise,rst,infocode)
        call call_bigdft(nproc,iproc,atoms,rxyz,       etot,fxyz,       fnoise,infocode,parini)
        ncount_bigdft = ncount_bigdft+1
        time_energy=-etot
end function
