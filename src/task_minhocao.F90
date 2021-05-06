subroutine task_minhocao(parini,parres)
 use global
 use defs_basis
 !use cell_utils
 use interface_code
 use fingerprint
 use mod_parini, only: typ_parini
 use yaml_output
 implicit none
 type(typ_parini), intent(inout):: parini
 type(typ_parini), intent(inout):: parres
 real(8),allocatable :: amass(:)
 real(8),allocatable :: fcart(:,:)
 real(8),allocatable :: vel(:,:),xcart(:,:)
 real(8),allocatable :: xred_in(:,:),vel_in(:,:),fcart_in(:,:)
 real(8):: vel_vol_in
 real(8):: etot_in,strten_in(6),vel_lat_in(3,3),latvec_in(3,3)

!******************************************************************
!Minima Hopping Variables
  real(8), parameter :: beta1=1.10d0,beta2=1.10d0,beta3=1.d0/1.10d0
  real(8), parameter :: alpha1=1.d0/1.10d0,alpha2=1.10d0
  !real(8):: beta1, beta2, beta3, alpha1, alpha2
  integer :: npminx=100        !Number of posloc files to be saved for verbosity 0
  integer :: nwrite_inter=100  !Number of structures to be found until intermediate files are written, only for verbosity 0
!  real(8),parameter :: HaBohr2GPA = 29421.033d0
!  real(8):: elocmin(npminx,5)         !elocmin(iloc,1)=enthalpy, elocmin(iloc,2)=energy, elocmin(iloc,3)=spg, elocmin(iloc,4)=spgtol, elocmin(iloc,5)=FDOS
  real(8):: ratio
  real(8):: eref!, accur
  real(8):: latvec(3,3)
  real(8):: ediff
  real(8):: ekinetic                  !In this version of MHM ekinetic=Temperature     
  real(8):: ekinetic_max              !Maximum Temperature during MD    
  real(8):: dt
  real(8):: fp_dist
  real(8):: fp_dist_min
  real(8):: fp_delta
  real(8):: ent_delta
  real(8):: fp_sep
  real(8):: tmp_acell(3),tmp_rprim(3,3), tmp_enthalpy, tmp_energy, tmp_real
!  real(8),allocatable:: earr(:,:)
  real(8),allocatable:: e_arr(:),ent_arr(:),fp_arr(:,:),pl_arr(:,:,:),f_arr(:,:,:),str_arr(:,:)
  real(8),allocatable:: lat_arr(:,:,:),dos_arr(:),spgtol_arr(:)
  integer,allocatable:: ct_arr(:),spg_arr(:)
  real(8),allocatable:: e_arr_t(:),ent_arr_t(:),fp_arr_t(:,:),pl_arr_t(:,:,:),f_arr_t(:,:,:),str_arr_t(:,:)
  real(8),allocatable:: lat_arr_t(:,:,:),dos_arr_t(:),spgtol_arr_t(:)
  integer,allocatable:: ct_arr_t(:),spg_arr_t(:)
  integer,allocatable:: lid(:)
  integer:: nid,n_arr
!  real(8),allocatable:: poslocmin(:,:,:)
!  real(8),allocatable:: latlocmin(:,:,:)
  real(8):: av_ekinetic
  real(8):: av_ediff
  real(8):: escape
  real(8):: escape_sam
  real(8):: escape_old
  real(8):: escape_new
  real(8):: egap
  real(8):: esep
  real(8):: count_geopt
  real(8):: count_soft
  real(8):: count_md
  real(8):: re_sm
  real(8):: tt,ss
  real(8):: rejected
  real(8):: accepted
  real(8):: ecut_tmp 
  real(8):: toldff_tmp 
  real(8):: counter 
  real(8):: counter_md 
  real(8):: counter_geopt
  real(8):: tolmin
  real(8):: tolmax
  real(8):: efermi
  real(8):: strtarget(6)
  real(8):: dstr(6)
  real(8):: fmax
  real(8):: fmax_at
  real(8):: fmax_lat
  real(8):: dist_ang(6)
  integer:: istr
  integer:: n_unique
  integer:: n_nonuni
  integer:: iprec
  integer:: kid
  integer:: itype 
  integer:: iexit
  integer:: nvisit
  integer:: nhop
  integer:: nputback
  integer:: nlmin,nlmin_t,nlminx,nlmin_old
  integer:: i_stat 
  integer:: i,k,kk
  integer:: k_e_wpos
  integer:: k_e_pos
  integer:: iat
!  integer:: npmin
  integer:: ntol
  integer:: spgint
  integer:: fp_len
  integer:: clock_start,clock_cur,clock_rate,clock_max,CPUlimit,elapsed
  character(40):: filename,folder
  character(4) :: fn4
  character(5) :: fn5
  character(10) :: tmp_char
  logical:: newmin
  logical:: getwfk
  logical:: file_exists
  logical:: readfix
  logical:: readfrag
  real(8),allocatable:: xcart_mol(:,:)
  real(8),allocatable:: xcart_tmp(:,:) 
  real(8),allocatable:: intens(:,:,:)
  real(8),allocatable:: inaxis(:,:,:)
  real(8),allocatable:: inprin(:,:)
  real(8),allocatable:: masstot(:)
!Set of data belonging to pos in initial minhop
  real(8),allocatable:: pos_red(:,:),pos_fcart(:,:),pos_cm(:,:),pos_quat(:,:),fp_pos(:)
  real(8):: pos_latvec(3,3)
  real(8):: pos_strten(6)
  real(8):: e_pos, ent_pos, spgtol_pos,fdos_pos
  integer:: spg_pos
!Set of data belonging to wpos in initial minhop
  real(8),allocatable:: wpos_red(:,:), wpos_fcart(:,:),fp_wpos(:)
  real(8):: wpos_latvec(3,3)
  real(8):: wpos_strten(6)
  real(8):: e_wpos, ent_wpos, spgtol_wpos,fdos_wpos
  integer:: spg_wpos
!Set of data belonging to poshop in initial minhop
  real(8),allocatable:: poshop(:,:), poshop_fcart(:,:),fp_hop(:)
  real(8):: lathop(3,3)
  real(8):: e_hop
  real(8):: ent_hop
  real(8):: poshop_strten(6)
  real(8):: fdos_hop,spgtol_hop
  integer:: spg_hop
!Latvec correction
  integer:: latvec_io
!Allocation for refined structure
  real(8),allocatable:: sym_pos(:,:),sym_fcart(:,:)
  real(8)::sym_cell(3,3)
  integer,allocatable::sym_type(:)
  integer::sym_nat,sym_nat2 
  logical,allocatable:: sym_fixat(:)

!Print the logo
!call print_logo() 

!Delete the file CPUlimit, if it exists
!  call system("rm -f CPUlimit")
call system_clock(count=clock_start)     !Start Timer
call system_clock(count_rate=clock_rate) !Find the time rate
call system_clock(count_max=clock_max)   !Find the time max

  !alpha1=parini%alpha1_minhopp !CORRECT_IT
  !alpha2=parini%alpha2_minhopp !CORRECT_IT
  !beta1=parini%beta1_minhopp !CORRECT_IT
  !beta2=parini%beta2_minhopp !CORRECT_IT
  !beta3=parini%beta3_minhopp !CORRECT_IT



!!  !Define the boundary condition: 1: periodic, 2:free, 3:surface/slab
!!    parini%bc=1
!!  
!!  !Verbosity
!!    parini%verb=1

!Define current directory
  folder=""

!!  !Initialize auto logicals to false, otherwise it will read the params file and reset alpha_lat, alpha_at, mdmin_max, etc
!!    parini%auto_soft=.false.
!!    parini%auto_mdmin=.false.
!!    parini%auto_dtion_md=.false.
!!    !parini%alpha_at=-1.d10
!!    !parini%alpha_lat=-1.d10 

!Initialize old kpt in history
  ka1=0;kb1=0;kc1=0;max_kpt=.false.;reuse_kpt=.false.

!Some parameter setup
!  code="lammps"          !Define what code to use for force evaluation
  parini%siesta_kpt_mode=2          !Mode of generating automatic k-point mesh, only siesta
  parini%vasp_kpt_mode=2            !Mode of generating automatic k-point mesh, only vasp
  parini%abinit_kpt_mode=1          !Mode of generating automatic k-point mesh, only abinit
  parini%correctalg=2               !Method to correct cell vectors when torn out of shape

!Define if the external optimizer should be used. Only available for:
!vasp
!siesta
!dftb
!!    parini%geopt_ext=.false.


!Unset fixed cell variables
  parini%fixlat=.false.

!Set units for all input/output for this run
  units="angstroem"
!  do idtset=1,ndtset_alloc
!     dtsets(idtset)%units=trim(units)
!  enddo

!Set FDOS to zero
   fdos_wpos=0.d0
   fdos_pos=0.d0
   fdos_hop=0.d0

!Set up parameters for symmetry tool
   tolmin=0.00005;tolmax=0.5d0;ntol=20

!Read input composition of the structure params.in
!  call read_params()
!Echo the parameters
!  call params_echo()
    !write(*,*) parini%params_new
    if(parini%params_new) then
        call params_read(parini)
    else
        call params_read_for_yaml(parini)
    endif
  parres=parini
!  call params_echo()

!Read poscur into pos_red, which is equivalent to what pos was in the initial minima hopping
  allocate(pos_red(3,parini%nat),pos_fcart(3,parini%nat),xcart_tmp(3,parini%nat),xcart_mol(3,parini%nat),&
          &wpos_red(3,parini%nat),wpos_fcart(3,parini%nat),poshop_fcart(3,parini%nat),vel_in(3,parini%nat))

!Finally, reading in atomic positions
  file_exists=.false.
  filename="poscur.ascii"
  readfix=.true.;readfrag=.true.
  INQUIRE(FILE=trim(filename), EXIST=file_exists)
  if(file_exists) then
    call yaml_map('Reading from',trim(filename))
  !write(*,*) "# Reading from ",trim(filename) 
  call read_atomic_file_ascii(filename,parini%nat,units,pos_red,pos_latvec,pos_fcart,pos_strten,parini%fixat,parini%fixlat,&
       &readfix,parini%fragarr,readfrag,tmp_enthalpy,tmp_energy)
  endif
  readfix=.false.;readfrag=.false.
!From a vasp file
  if(.not.file_exists) then
     file_exists=.false.
     filename="poscur.vasp"
     INQUIRE(FILE=trim(filename), EXIST=file_exists)
     if(file_exists) then
     write(*,*) "# Reading from ",trim(filename) 
     filename="poscur.vasp"
     call read_atomic_file_poscar(filename,parini%nat,units,pos_red,pos_latvec,pos_fcart,pos_strten,parini%fixat,parini%fixlat,&
          &readfix,parini%fragarr,readfrag,tmp_enthalpy,tmp_energy)
     endif
  endif
  if(.not.file_exists) then
     stop "Provide either poscur.ascii or poscur.vasp"
  endif




!Write FIXLAT
  if(any(parini%fixlat)) write(*,*) "# FIXLAT: a,b,c,alpha,beta,gamma,shape ", parini%fixlat

!Initiallize Fingerprint
!  fp_method=11 !11: Oganov FP, 12: CALYPSO FP, 13: Modified CALYPSO, 21: molecular gaussian overlap, 22: molecular sprint
  call init_fp(parini,fp_len,pos_latvec)
  allocate(fp_pos(fp_len),fp_wpos(fp_len),fp_hop(fp_len))

!Check correct assignment of FP method
  if(parini%bc==1.and.(fp_method.lt.11.or.fp_method.gt.19)) stop "Incompatible fingerprint and boundary conditions"
  if(parini%bc==2.and.(fp_method.lt.21.or.fp_method.gt.29)) stop "Incompatible fingerprint and boundary conditions"
  if(parini%bc==3.and.(fp_method.lt.31.or.fp_method.gt.39)) stop "Incompatible fingerprint and boundary conditions"


!Assign atomic masses in electron units
 allocate(amass(parini%nat))
 do iat=1,parini%nat
   amass(iat)=amu_emass*parini%amu(parini%typat_global(iat))
 enddo

!Allocate and handle molecular crystal stuff
   call refragment(parini%fragarr,parini%nat)
   nmol=maxval(parini%fragarr)
   allocate(parini%fragsize(nmol),parini%lhead(nmol),parini%llist(parini%nat))
   call make_linked_list(parini%fragarr,parini%fragsize,parini%lhead,parini%llist,parini%nat,nmol)
!Here we will decompose the atomic position into a set of cm and quat 
   allocate(pos_cm(3,nmol),pos_quat(4,nmol),intens(3,3,nmol),inaxis(3,3,nmol),inprin(3,nmol),masstot(nmol))
   if(nmol.ne.parini%nat) then 
       call init_cm_mol(parini,pos_latvec,pos_red,xcart_mol,pos_cm,pos_quat,amass,masstot,intens,inprin,inaxis,parini%lhead,parini%llist,parini%nat,nmol)
!Get the inertia tensor and the principle inertia and the principle axis
       deallocate(xcart_tmp);allocate(xcart_tmp(3,nmol))
       xcart_tmp=0.d0
       do iat=1,nmol
          write(122,*) "INTENS",iat
          write(122,*) intens(:,1,iat)
          write(122,*) intens(:,2,iat)
          write(122,*) intens(:,3,iat)
          write(122,*) "INPRIN"
          write(122,*) inprin(:,iat)
          write(122,*) "INAXIS",iat
          write(122,*) inaxis(:,1,iat)
          write(122,*) inaxis(:,2,iat)
          write(122,*) inaxis(:,3,iat)
       enddo
   endif

!Check if confinement is desired
!  file_exists=.false.
!  filename="confine.in"
!  INQUIRE(FILE=trim(filename), EXIST=file_exists)
     if(parini%use_confine) then
!        call init_confinement(nat,filename)
        confine=1!Or 3 if confinement should be always switched on
     else
        confine=0
     endif

#if defined(SPGLIB)
!Check if poscur.ascii should be symmetrized
  file_exists=.false.
  INQUIRE(FILE="checksym.in", EXIST=file_exists)
     if(file_exists) then
       write(*,*)'# checkysm.in found, performing path symmetrization of poscur.ascii'
     open(unit=456,file="checksym.in")
     read(456,*) tolmin,tolmax,ntol
     close(456)
     !Put all the atoms back into the cell
     call backtocell(parini%nat,pos_latvec,pos_red)
     allocate(sym_pos(3,4*parini%nat),sym_type(4*parini%nat),sym_fcart(3,4*parini%nat),sym_fixat(4*parini%nat))
     sym_pos(:,1:parini%nat)=pos_red
     sym_type(1:parini%nat)=parini%typat_global
     sym_fcart=0.d0
     sym_fixat=.false.
     pos_strten=0.d0
     parini%fragarr=0
     call find_symmetry(parini,parini%nat,pos_red,pos_latvec,parini%typat_global,tolmin,tolmax,ntol,spgtol_pos,spg_pos)
     call spg_cell_refine(parini%nat,sym_nat,4*parini%nat,sym_pos,pos_latvec,sym_type,spgtol_pos,spg_pos)
     filename="sym_poscur.ascii"
     call write_atomic_file_ascii(parini,filename,sym_nat,units,sym_pos,pos_latvec,sym_fcart,pos_strten,parini%char_type,parini%ntypat_global,sym_type,sym_fixat,parini%fixlat,&
     &0.d0,parini%target_pressure_habohr,0.d0,0.d0)
     call spg_cell_primitive(sym_nat,sym_nat2,4*parini%nat,sym_pos,pos_latvec,sym_type,spgtol_pos,spg_pos)
     filename="sym_prim_poscur.ascii"
     call write_atomic_file_ascii(parini,filename,sym_nat2,units,sym_pos,pos_latvec,sym_fcart,pos_strten,parini%char_type,parini%ntypat_global,sym_type,sym_fixat,parini%fixlat,&
     &0.d0,parini%target_pressure_habohr,0.d0,0.d0)
       goto 3001 
  endif 
#endif

!Check if pathintegral should be performed
  file_exists=.false.
  INQUIRE(FILE="pathparams.in", EXIST=file_exists)
     if(file_exists) then
       write(*,*)'# pathparams.in found, performing path integral starting from poscur.ascii'
     !Put all the atoms back into the cell
       call backtocell(parini%nat,pos_latvec,pos_red)
     !Run path integral
       call pathintegral(parini,parres,pos_latvec,pos_red)
       goto 3001 
  endif 


!Check if enthalpy sampling should be performed
  file_exists=.false.
  INQUIRE(FILE="enthalpyparams.in", EXIST=file_exists)
     if(file_exists) then
       write(*,*)'# enthalpyparams.in found, performing enthalpy scan starting from poscur.ascii'
     !Put all the atoms back into the cell
       call backtocell(parini%nat,pos_latvec,pos_red)
     !Run path integral
       call enthalpyrelax(parini,parres,pos_latvec,pos_red,tolmin,tolmax,ntol,parini%findsym)
       goto 3001 
  endif 


!Check if volume sampling should be performed
  file_exists=.false.
  INQUIRE(FILE="varvol.in", EXIST=file_exists)
     if(file_exists) then
       write(*,*)'# varvol.in found, performing volume scan starting from poscur.ascii'
     !Put all the atoms back into the cell
       call backtocell(parini%nat,pos_latvec,pos_red)
     !Run path integral
       call varvol(parini,parres,pos_latvec,pos_red,tolmin,tolmax,ntol,parini%findsym)
       goto 3001 
  endif 

!Check if a list of poslows should be relaxed
  file_exists=.false.
  INQUIRE(FILE="poslowrelax.in", EXIST=file_exists)
     if(file_exists) then
       write(*,*)'# poslowrelax.in found, performing relaxations of all poslows found &
                  &in current folder, but poscur.ascii must be provided for initiallization.'
     !Put all the atoms back into the cell
       call backtocell(parini%nat,pos_latvec,pos_red)
     !Run path integral
       call poslowrelax(parini,parres,pos_latvec,pos_red,tolmin,tolmax,ntol)
       goto 3001
  endif

!Check if rotations should be performed
  file_exists=.false.
  INQUIRE(FILE="rotate.in", EXIST=file_exists)
     if(file_exists) then
       write(*,*)'# rotate.in found, performing rotations of all poscur file found &
                  &in current folder, poscur.ascii must be provided for initiallization.'
     !Put all the atoms back into the cell
       call backtocell(parini%nat,pos_latvec,pos_red)
     !Run path integral
       call rotate_like_crazy(parini,parres,pos_latvec,pos_red,tolmin,tolmax,ntol)
       goto 3001
  endif


!Write initial parameters to global.out
  ratio=0.d0
!  open(unit=67,file='global.out')
  call yaml_mapping_open('initial minima hopping parameters')
  call yaml_map('beta1',beta1,fmt='(es11.4)')
  call yaml_map('beta2',beta2,fmt='(es11.4)')
  call yaml_map('beta3',beta3,fmt='(es11.4)')
  call yaml_map('alpha1',alpha1,fmt='(es11.4)')
  call yaml_map('alpha2',alpha2,fmt='(es11.4)')
  call yaml_map('predicted fraction accepted',ratio/(1.d0+ratio),fmt='(es11.4)')
  call yaml_map('predicted fraction rejected',1.d0/(1.d0+ratio),fmt='(es11.4)')
  call yaml_mapping_close()
  !write(*,'(a,3(1x,1pe11.4))') ' # beta1,beta2,beta3',beta1,beta2,beta3
  !write(*,'(a,2(1x,1pe11.4))') ' # alpha1,alpha2',alpha1,alpha2
  !write(*,'(a,2(1x,1pe10.3))') ' # predicted fraction accepted, rejected', &
  !& ratio/(1.d0+ratio), 1.d0/(1.d0+ratio)
  call yaml_map('mdmin',parres%mdmin)
  !write(*,*) '# mdmin',parres%mdmin
  

!Read earr.dat. In this version, earr contains the enthalpies, not the energies, since they are compared during MinHopp
  open(unit=12,file='earr.dat',status='unknown')
  read(12,*) nlmin,nlminx
  if(parini%verb.gt.0) npminx=nlminx
  if(parini%verb.gt.0) nwrite_inter=1 
!  read(12,*) eref  !Eref is the reference enthalpy
!  write(67,*) 'eref=',eref
!  write(*,*) 'eref=',eref
!  read(12,*) accur
!  write(67,*) 'accuracy for rounding=',accur
!  write(*,*) 'accuracy for rounding=',accur
  read(12,*) ent_delta, fp_delta
  call yaml_mapping_open('Tolerance between structures',flow=.true.)
  call yaml_map('delta enthalpy',ent_delta,fmt='(es15.7)')
  call yaml_map('delta fingeprint',fp_delta,fmt='(es15.7)')
  call yaml_mapping_close()
  !write(*,'(a,2(es15.7))') " # Tolerance between structures: delta enthalpy, delta fingeprint: ", ent_delta,fp_delta 
  if (nlmin.gt.nlminx) stop 'nlmin>nlminx'
!  allocate(earr(0:nlminx,3),stat=i_stat) !This earr is modified from the original version: earr(i,1)=enthalpy,earr(i,2)=energy,earr(i,3)=visits
  allocate(e_arr(nlminx),ent_arr(nlminx),ct_arr(nlminx),pl_arr(3,parini%nat,nlminx),f_arr(3,parini%nat,nlminx),&
       &str_arr(6,nlminx),lat_arr(3,3,nlminx),dos_arr(nlminx),spgtol_arr(nlminx),spg_arr(nlminx),&
       &fp_arr(fp_len,nlminx),lid(nlminx))
!  earr=0.d0
!  earr(0,1)=-1.d100
  e_arr=0.d0;ent_arr=0.d0;dos_arr=0.d0
  if (nlmin == 0) then 
!     write(67,*) ' New run with nlminx=',nlminx
     call yaml_map('New run with nlminx',nlminx)
     !write(*,*) '# New run with nlminx=',nlminx
  else
!     write(67,*) ' Restart run with nlmin, nlminx=',nlmin,nlminx
     write(*,*) '# Restart run with nlmin, nlminx=',nlmin,nlminx
     do k=1,nlmin
!        read(12,*) earr(k,1),earr(k,2),earr(k,3) !Modified to contain enthalpy
        read(12,*) kk, ent_arr(k),e_arr(k),ct_arr(k),spg_arr(k),spgtol_arr(k),dos_arr(k) !Modified to contain enthalpy
!        if (earr(k,1).lt.earr(k-1,1)) stop 'wrong ordering in earr.dat'
        if (ent_arr(k).lt.ent_arr(max(1,k-1))) stop 'wrong ordering in earr.dat'
     enddo
!     write(67,*) ' read earr.dat'
     write(*,*) '# read earr.dat'
  endif
  close(12)

!Currently the fixed cell shape is implemented
  if(parini%fixlat(7)) then
          parini%md_algo = 4
          write(*,'(A,i5)') " # Variable volume algorithm selected, switching to Anderson MD: parini%md_algo = ", parini%md_algo
  endif
  if(parini%md_algo==4) then
         parini%fixlat(7)=.true.
         write(*,'(A)') " # Anderson MD selected, switching FIXLAT"
         if(any(parini%fixlat)) write(*,*) "# FIXLAT: a,b,c,alpha,beta,gamma,shape ", parini%fixlat
  endif

  if(parini%fixlat(7).and.trim(parini%paropt_geopt%approach).ne."FIRE") stop "Fixed cell shape only implemented in FIRE"

!Put all the atoms back into the cell
  !call backtocell(nat,pos_latvec,pos_red)

!Allocate some other arrays
!  allocate(poslocmin(3,nat,npminx))
!  allocate(latlocmin(3,3,npminx))
  allocate(poshop(3,parini%nat))

!Read input parameters from ioput
  write(filename,'(a6,i3.3)') 'ioput'   
  open(unit=11,file='ioput',status='old')
  read(11,*) ediff,ekinetic,ekinetic_max
  close(11)

  call yaml_mapping_open('In minhocao',flow=.true.)
  call yaml_map('ediff',ediff,fmt='(es11.3)')
  call yaml_map('temperature',ekinetic,fmt='(es11.3)')
  call yaml_map('maximum temperature',ekinetic_max,fmt='(es11.3)')
  call yaml_map('nsoften',parini%nsoften_minhopp,fmt='(i8)')
  call yaml_mapping_close()
  !write(*,'(a,1x,3(1x,e10.3),1x,i4)') ' # In :ediff,temperature,maximum temperature,nsoften',&
  !      &ediff,ekinetic,ekinetic_max,parini%nsoften_minhopp

!If restart run read previously found energies and other infos
!  elocmin=0.d0
  if (nlmin > 0) then
!     open(unit=787,file="strlist.dat")
!!     npmin=0
!!     kk=1
!!     do
!!        if (npmin.ge.npminx .or. kk.gt.min(nlmin,npminx)) then
!!           exit
!!        end if
    do kk=1,nlmin
!!        npmin=npmin+1
     if(parini%verb.gt.0) then
      write(fn5,'(i5.5)') kk
      filename = 'poslow'//fn5//'.ascii'
!        call read_atomic_file_ascii(filename,nat,units,poslocmin(1:3,1:nat,npmin),latlocmin(1:3,1:3,npmin),fixat,fixlat,readfix,&
      call read_atomic_file_ascii(filename,parini%nat,units,pl_arr(1:3,1:parini%nat,kk),lat_arr(1:3,1:3,kk),f_arr(1:3,1:parini%nat,kk),str_arr(1:6,kk),&
           &parini%fixat,parini%fixlat,readfix,parini%fragarr,readfrag,tmp_enthalpy,tmp_energy)
      write(*,*) '# read file ',filename
     elseif(kk==1) then 
      n_arr=3*parini%nat*nlmin
      filename="poslow.bin"
      call bin_read(filename,pl_arr,n_arr)
      filename="fcart.bin"
      file_exists=.false.; INQUIRE(FILE=trim(filename), EXIST=file_exists)
      if(file_exists) call bin_read(filename,f_arr,n_arr)
      n_arr=3*3*nlmin
      filename="latvec.bin"
      call bin_read(filename,lat_arr,n_arr)
      n_arr=6*nlmin
      filename="strten.bin"
      file_exists=.false.; INQUIRE(FILE=trim(filename), EXIST=file_exists)
      if(file_exists) call bin_read(filename,str_arr,n_arr)
      write(*,'(a)') ' # read binary files poslow.bin, fcart.bin, strten.bin and latvec.bin'
     endif

!      write(67,*)'read file',filename 
!Read also the strlist.dat, actually it overwrites both elocmin(npmin,1) and elocmin(npmin,2)
!        read(787,*) tmp_real, elocmin(npmin,1:5)          
!Recomute symmetries if not available for some of the structures
!        if(findsym.and.elocmin(npmin,3).lt.1.d0) then
        if(parini%findsym.and.spg_arr(kk).lt.1) then
          write(*,'(a,a)') " # Recomputing symmetry for ",trim(filename)
!          call find_symmetry(parini,nat,poslocmin(1:3,1:nat,npmin),latlocmin(1:3,1:3,npmin),typat,tolmin,tolmax,ntol,&
!          &elocmin(npmin,4),spgint)
          call find_symmetry(parini,parini%nat,pl_arr(1:3,1:parini%nat,kk),lat_arr(1:3,1:3,kk),parini%typat_global,tolmin,tolmax,ntol,&
          &spgtol_arr(kk),spg_arr(kk))
!          elocmin(npmin,3)=real(spgint,8)
        endif
!Compute the fingerprints for all poslow structures
        n_arr=fp_len*nlmin
        filename="fp.bin"
        file_exists=.false.; INQUIRE(FILE=trim(filename), EXIST=file_exists)
        if(file_exists) then 
           if(kk==1) then
             call bin_read(filename,fp_arr,n_arr)
             write(*,*) "# read fingerprints from file fp.bin"
           endif
        else
           call get_fp(parini,fp_len,pl_arr(:,:,kk),lat_arr(:,:,kk),fp_arr(:,kk))
           write(*,*) "# fingerprints computed from poslow files"
        endif
!!        kk=kk+1
     end do
!     write(67,*)'read ',nlmin,'poslow files'
  write(*,*) '# read ',nlmin,'poslow files'
!     close(787)
  endif

!Check if the informations provided should be rebinned according to the fingerprint function
!Exit on success
  file_exists=.false.
  INQUIRE(FILE="rebin.in", EXIST=file_exists)
     if(file_exists) then
       write(*,'(a)') " # rebin.in found, trying to rebin poslow files. Will write new poslow and earr.dat"
       if(nlmin.le.1) then
         write(*,'(a,i5)') " # Nothing to do, nlmin = ",nlmin
         goto 3001
       endif
       write(*,'(a,i5)') " # Running over all nlmin = ",nlmin
       allocate(e_arr_t(nlminx),ent_arr_t(nlminx),fp_arr_t(fp_len,nlminx),pl_arr_t(3,parini%nat,nlminx),f_arr_t(3,parini%nat,nlminx),&
            &str_arr_t(6,nlminx),lat_arr_t(3,3,nlminx),dos_arr_t(nlminx),spgtol_arr_t(nlminx),ct_arr_t(nlminx),spg_arr_t(nlminx))
       e_arr_t(1)=e_arr(1);ent_arr_t(1)=ent_arr(1);fp_arr_t(:,1)=fp_arr(:,1);pl_arr_t(:,:,1)=pl_arr(:,:,1)
       lat_arr_t(:,:,1)=lat_arr(:,:,1);dos_arr_t(1)=dos_arr(1);spgtol_arr_t(1)=spgtol_arr(1);ct_arr_t(1)=ct_arr(1)
       spg_arr_t(1)=spg_arr(1)
       f_arr_t=0.d0;str_arr_t=0.d0
       nlmin_t=1
       kk=1
       do 
          kk=kk+1
          call identical(parini,nlminx,nlmin_t,fp_method,fp_len,ent_arr(kk),fp_arr(:,kk),ent_arr_t,fp_arr_t,&
                         &ent_delta,fp_delta,newmin,kid,fp_dist_min,k_e_pos,n_unique,n_nonuni,lid,nid)
          if (.not.newmin) then                    !The case the structure was previously found
!!             do i=1,nid
                ct_arr_t(kid)=ct_arr_t(kid)+ct_arr(kk)
!!                if(findsym) call replace(nlminx,nlmin_t,fp_len,nat,lid(1),e_arr(lid(i)),ent_arr(lid(i)),fp_arr(:,lid(i)),pl_arr(:,:,lid(i)),&
!!                         &lat_arr(:,:,lid(i)),spg_arr(kk),spgtol_arr(kk),dos_arr(kk),&
!!                         &e_arr_t,ent_arr_t,fp_arr_t,pl_arr_t,lat_arr_t,spg_arr_t,spgtol_arr_t,dos_arr_t,ct_arr_t,findsym)
!!             enddo
          else                                     !A new minimum was found
             nlmin_t=nlmin_t+1
             call insert(nlminx,nlmin_t,fp_len,parini%nat,k_e_pos,e_arr(kk),ent_arr(kk),fp_arr(:,kk),pl_arr(:,:,kk),&
                         &lat_arr(:,:,kk),f_arr(:,:,kk),str_arr(:,kk),spg_arr(kk),spgtol_arr(kk),dos_arr(kk),&
                         &e_arr_t,ent_arr_t,fp_arr_t,pl_arr_t,lat_arr_t,f_arr_t,str_arr_t,spg_arr_t,spgtol_arr_t,dos_arr_t,ct_arr_t)
             k_e_pos=k_e_pos+1
          endif
          if(kk.eq.nlmin) exit
       enddo 
     call winter(parini,parini%nat,units,ent_pos,e_pos,pos_red,pos_latvec,pos_fcart,pos_strten,nlminx,nlmin_t,npminx,& 
          &ent_arr_t,e_arr_t,ct_arr_t,spg_arr_t,spgtol_arr_t,dos_arr_t,pl_arr_t,lat_arr_t,f_arr_t,str_arr_t,fp_arr_t,fp_len,ent_delta,fp_delta,& 
          &eref,ediff,ekinetic,ekinetic_max,dt,parini%nsoften_minhopp,parini%char_type,parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,parini%target_pressure_habohr)
     write(*,'(a,i5,a)') " # Rebinning completed, new nlmin: ",nlmin_t,". rebin.in will be deteled!"
     call system("rm rebin.in")
     goto 3001 
  endif 
     
    
!Just print a fingerprint distance grid
  file_exists=.false.
  INQUIRE(FILE="plot_fpgrid.in", EXIST=file_exists)
     if(file_exists) then
       write(*,'(a)') " # plot_fpgrid.in found. Will write Howtoplot_fingerprint"
       call plot_fp_grid(parini,nlminx,nlmin,parini%nat,fp_len,fp_arr,lat_arr,pl_arr)
     goto 3001
     endif



!Initialize statistical variables
  n_unique=0
  n_nonuni=0
  av_ekinetic=0.d0
  av_ediff=0.d0
  escape=0.d0
  escape_sam=0.d0
  escape_old=0.d0
  escape_new=0.d0
  rejected=0
  accepted=0
  egap=1.d100
  esep=0.d0
  fp_sep=0.d0
  e_hop=1.d100
  ent_hop=1.d100
  count_geopt=0.d0
  count_soft=0.d0
  count_md=0.d0
  nputback=0
  eref=0.d0

!Count the number of poslocm files in the folder
  file_exists=.true.
  nhop=-1
  do while(file_exists)
     nhop=nhop+1
     write(fn5,'(i5.5)') nhop
     filename='poslocm_'//fn5//'.ascii'
     file_exists=.false.
     INQUIRE(FILE=trim(filename), EXIST=file_exists)
  enddo
  call yaml_map('Number of poslocm files found',nhop)
!write(*,'(a,i5)') " # Number of poslocm_ files found: ",nhop
     write(fn5,'(i5.5)') nhop
     if(parini%verb.ge.2) folder="data_hop_init/"
     if(parini%verb.ge.2) call system("mkdir "//trim(folder))


!A single point calculation to get the initial energy (will be converted to enthalpy)
  iprec=1;getwfk=.false.
  call get_energyandforces_single(parini,parres,pos_latvec,pos_red,pos_fcart,pos_strten,e_pos,iprec,getwfk)
!Convert to enthalpy
  call get_enthalpy(pos_latvec,e_pos,parini%target_pressure_habohr,ent_pos)
!Write to file and exit if geopt iteration is 0
  if(parini%verb.gt.0) then
       filename=trim(folder)//"posinit.ascii"
       call write_atomic_file_ascii(parini,filename,parini%nat,units,pos_red,pos_latvec,pos_fcart,pos_strten,&
            &parini%char_type(1:parini%ntypat_global),parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,e_pos,parini%target_pressure_habohr,ent_pos,e_pos)
       if(parini%verb.ge.3) then
       filename=trim(folder)//"posinit.vasp"
       call write_atomic_file_poscar(parini,filename,parini%nat,units,pos_red,pos_latvec,pos_fcart,pos_strten,&
            &parini%char_type(1:parini%ntypat_global),parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,e_pos,parini%target_pressure_habohr,ent_pos,e_pos)
       endif
  endif
  if (parini%paropt_geopt%nit.le.0) goto 3000
!Check if the structure is already relaxed
  call convcheck(parini,parini%nat,pos_latvec,pos_fcart,pos_strten,parini%target_pressure_habohr,parini%paropt_geopt%strfact,fmax,fmax_at,fmax_lat,parini%paropt_geopt%fmaxtol,iexit)
if(iexit==1) then
     call yaml_mapping_open('Input structure already relaxed',flow=.true.)
     call yaml_map('ent_pos',ent_pos,fmt='(es20.12)')
     call yaml_map('e_pos',e_pos,fmt='(es20.12)')
     call yaml_mapping_close()
  !write(*,'(a,es15.7,es15.7)') ' # Input structure already relaxed. Proceeding with: Energy, Enthalpy: ',e_pos,ent_pos
  else 
     call yaml_mapping_open('Calling the geometry optimizer for the first time here',flow=.true.)
     call yaml_map('ent_pos',ent_pos,fmt='(es20.12)')
     call yaml_map('e_pos',e_pos,fmt='(es20.12)')
     call yaml_mapping_close()
  !write(*,'(a,es15.7,es15.7)') ' # Calling the geometry optimizer for the first time here. Energy, Enthalpy: ',e_pos,ent_pos

!Call geometry optimizer for the first time
  vel_in=0.d0
  vel_lat_in=0.d0
  vel_vol_in=0.d0
   iprec=1
   if(parini%geopt_ext) then
     call geopt_external(parini,pos_latvec,pos_red,pos_fcart,pos_strten,e_pos,iprec,parres%ka,parres%kb,parres%kc,counter)
   else
     if(confine==2) confine=1
     if(parini%paropt_geopt%approach=="RBFGS")  call GEOPT_RBFGS_MHM(parini,parres,pos_latvec,pos_red,pos_fcart,pos_strten,e_pos,iprec,counter,folder)
     if(parini%paropt_geopt%approach=="MBFGS")  call GEOPT_MBFGS_MHM(parini,parres,pos_latvec,pos_red,pos_fcart,pos_strten,e_pos,iprec,counter,folder)
     if(parini%paropt_geopt%approach=="FIRE")   call GEOPT_FIRE_MHM(parini,parres,pos_latvec,pos_red,pos_fcart,pos_strten,vel_in,vel_lat_in,&
                                                    &vel_vol_in,e_pos,iprec,counter,folder)
     if(parini%paropt_geopt%approach=="SQNM")   call      GEOPT_SQNM(parini,parres,pos_latvec,pos_red,pos_fcart,pos_strten,e_pos,iprec,counter,folder)
     if(parini%paropt_geopt%approach=="QBFGS")  call     GEOPT_qbfgs(parini,parres,pos_latvec,pos_red,pos_fcart,pos_strten,e_pos,iprec,counter,folder)
     if(parini%paropt_geopt%approach=="SD")     call        GEOPT_SD(parini,parres,pos_latvec,pos_red,pos_fcart,pos_strten,e_pos,iprec,counter,folder)
   endif

!Update GEOPT counter
     count_geopt=count_geopt+counter     
     call yaml_map('Counter of GEOPT updated',int(count_geopt))
     !write(*,'(a,i7)') " # Counter of GEOPT updated: ", int(count_geopt)
endif

!Check for symmetry
   if(parini%findsym) then
      call find_symmetry(parini,parini%nat,pos_red,pos_latvec,parini%typat_global,tolmin,tolmax,ntol,spgtol_pos,spgint)
      spg_pos=spgint
   else
      spg_pos=0
      spgtol_pos=0.d0
   endif

!Accept first structure
   accepted=accepted+1.d0

!Get FDOS
  if(parini%finddos) then
     iprec=1
     if(parini%usewf_geopt) then
        getwfk=.true.
     else
        getwfk=.false.
     endif
     call get_dos(parini,parres,pos_latvec,pos_red,efermi,fdos_pos,iprec,getwfk)
     write(*,'(a,es15.7,es15.7)') " # FDOS found: ",fdos_pos,efermi
  endif

!Re-read parameters
!   write(*,'(a)')  " # Re-reading params.in"
!   call read_params()
  !call params_read(parini)

!Get the fingerprint
   call get_fp(parini,fp_len,pos_red,pos_latvec,fp_pos)

!Convert to enthalpy
  call get_enthalpy(pos_latvec,e_pos,parini%target_pressure_habohr,ent_pos)

!  nconjgr=0
  if(parini%verb.gt.0) then
     filename='poslocm_'//fn5//'.ascii'
     call write_atomic_file_ascii(parini,filename,parini%nat,units,pos_red,pos_latvec,pos_fcart,pos_strten,parini%char_type,parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,&
     &e_pos,parini%target_pressure_habohr,ent_pos-eref,e_pos-eref)
     call dist_latvec2ang(dist_ang,pos_latvec,pi)
     call yaml_mapping_open('Wrote file',flow=.true.)
     call yaml_map('filename',trim(filename))
     call yaml_map('a',dist_ang(1),fmt='(es15.7)')
     call yaml_map('b',dist_ang(2),fmt='(es15.7)')
     call yaml_map('c',dist_ang(3),fmt='(es15.7)')
     call yaml_map('alpha',dist_ang(4),fmt='(es15.7)')
     call yaml_map('beta',dist_ang(5),fmt='(es15.7)')
     call yaml_map('gamma',dist_ang(6),fmt='(es15.7)')
     call yaml_mapping_close()
     !write(*,'(a,a,a,6(1x,es15.7))') " # Wrote file ",trim(filename),", with a,b,c,alpha,beta,gamma: ",dist_ang
  endif

!Do rounding of the received enthalpy
!!  re_pos=round(e_pos-eref,accur)
!!  rent_pos=round(ent_pos-eref,accur)
!!     write(67,'(a,1x,5(1x,1pe17.10))') 'INPUT(relaxed): ent_pos,rent_pos,e_pos,re_pos,eref ',ent_pos,rent_pos,e_pos,re_pos,eref
!!     write(*,'(a,1x,5(1x,1pe17.10))') ' # INPUT(relaxed): ent_pos,rent_pos,e_pos,re_pos,eref ',ent_pos,rent_pos,e_pos,re_pos,eref
!     write(67,'(a,1x,2(1x,1pe17.10))') 'INPUT(relaxed): ent_pos,e_pos ',ent_pos,e_pos
     call yaml_mapping_open('INPUT(relaxed)',flow=.true.)
     call yaml_map('ent_pos',ent_pos,fmt='(es20.12)')
     call yaml_map('e_pos',e_pos,fmt='(es20.12)')
     call yaml_mapping_close()
     !write(*,'(a,1x,2(1x,1pe17.10))') ' # INPUT(relaxed): ent_pos,e_pos ',ent_pos,e_pos
  if (nlmin.gt.0) then
!      write(67,'(a,2(1x,1pe24.17))') 'new/old enthalpy for input file',ent_pos
      write(*,'(a,2(1x,1pe24.17))') '# new/old enthalpy for input file',ent_pos
  endif

!Place the relaxed structure in poslocmin, latlocmin, earr, elocmin
  k_e_wpos=1
  if (nlmin == 0) then
     nlmin=1
!!     npmin=1
     ct_arr(1)=1
     e_arr(1)=e_pos
     ent_arr(1)=ent_pos
     spg_arr(1)=spg_pos
     spgtol_arr(1)=spgtol_pos
     dos_arr(1)=fdos_pos
     fp_arr(:,1)=fp_pos
!!     earr(1,1)=rent_pos
!!     earr(1,2)=re_pos
!!     earr(1,3)=1.d0
!!     elocmin(1,1)=rent_pos
!!     elocmin(1,2)=re_pos
!!     elocmin(1,3)=spg_pos
!!     elocmin(1,4)=spgtol_pos
!!     elocmin(1,5)=fdos_pos
     do iat=1,parini%nat
!!        poslocmin(1,iat,1)=pos_red(1,iat)
!!        poslocmin(2,iat,1)=pos_red(2,iat)
!!        poslocmin(3,iat,1)=pos_red(3,iat)
        pl_arr(1,iat,1)=pos_red(1,iat)
        pl_arr(2,iat,1)=pos_red(2,iat)
        pl_arr(3,iat,1)=pos_red(3,iat)
     enddo
!!     latlocmin(:,:,1)=pos_latvec(:,:)
     lat_arr(:,:,1)=pos_latvec(:,:)
     !write(*,*) '#first configuration saved'
     call yaml_comment('first configuration saved',hfill='~')
  else
!           check whether new minimum
            call identical(parini,nlminx,nlmin,fp_method,fp_len,ent_pos,fp_pos,ent_arr,fp_arr,&
                 &ent_delta,fp_delta,newmin,kid,fp_dist_min,k_e_pos,n_unique,n_nonuni,lid,nid)
!!            call hunt(earr(1,1),min(nlmin,nlminx),rent_pos,k_e_pos)
       if (.not.newmin) then
         write(*,*) '#initial minimum is old '
         if (kid .gt. nlminx .or. kid .lt. 1) stop "kid out of bounds"
         nvisit=ct_arr(kid)
!!         if(findsym) call replace(nlminx,nlmin,fp_len,nat,kid,e_pos,ent_pos,fp_pos,pos_red,&
!!           &pos_latvec,spg_pos,spgtol_pos,fdos_pos,&
!!           &e_arr,ent_arr,fp_arr,pl_arr,lat_arr,spg_arr,spgtol_arr,dos_arr,ct_arr,findsym)
  
       else
         write(*,*) '#initial minimum is new '
         nlmin=nlmin+1
         if (nlmin.gt.nlminx) stop 'nlminx too small'
         !            add minimum to history list
         call insert(nlminx,nlmin,fp_len,parini%nat,k_e_pos,e_pos,ent_pos,fp_pos,pos_red,pos_latvec,pos_fcart,pos_strten,&
                     &spg_pos,spgtol_pos,fdos_pos,e_arr,ent_arr,fp_arr,pl_arr,lat_arr,f_arr,str_arr,&
                     &spg_arr,spgtol_arr,dos_arr,ct_arr)
         k_e_pos=k_e_pos+1
         if(ct_arr(k_e_pos).ne.1) stop "New min, count must be 1"
         if (k_e_pos .gt. nlminx .or. k_e_pos .lt. 1) stop "k_e_pos out of bounds"
         nvisit=ct_arr(k_e_pos)
       endif
!!!            if (rent_pos.eq.earr(k_e_pos,1)) then
!!!              write(*,*) '#initial minimum is old '
!!!            else
!!!              write(*,*) '#initial minimum is new '
!!!              nlmin=nlmin+1
!!!!             add minimum to history list
!!!              call insert(nlminx,nlmin,k_e_pos,rent_pos,re_pos,earr(0,1))
!!!!             save configuration if it is among the lowest ones in energy
!!!                npmin=npmin+1
!!!                call save_low_conf(nat,npmin,npminx,rent_pos,re_pos,pos_red,pos_latvec,&
!!!                &spg_pos,spgtol_pos,fdos_pos,elocmin,poslocmin,latlocmin)
!!!                k_e_wpos=k_e_wpos+1
!!!            endif
  endif


!!!  re_sm=min(rent_pos,earr(1,1))
  re_sm=min(ent_pos,ent_arr(1))
  open(unit=222,file='global.mon',status='unknown',position='append')
!     write(67,*) 'initial re_sm',re_sm
     call yaml_map('initial re_sm',re_sm)
     !write(*,*) '# initial re_sm',re_sm
     write(222,'(i10,(1x,f10.0),1x,1pe21.14,2(1x,1pe10.3),i5,1x,1pe10.3)')&
                &nhop,escape,ent_pos-eref,ediff,ekinetic,spg_pos,fdos_pos
    call yaml_mapping_open('global',flow=.true.)
    call yaml_map('escape',escape,fmt='(f10.0)')
    call yaml_map('de',ent_pos-eref,fmt='(es21.14)')
    call yaml_map('ediff',ediff,fmt='(es10.3)')
    call yaml_map('ekinetic',ekinetic,fmt='(es10.3)')
    call yaml_map('spg_pos',spg_pos,fmt='(i8)')
    call yaml_map('fdos_pos',fdos_pos,fmt='(es10.3)')
    call yaml_mapping_close()
     !write(*,'(a,(1x,f10.0),1x,1pe21.14,2(1x,1pe10.3),i5,1x,1pe10.3)')" #global: ",&
     !&escape,ent_pos-eref,ediff,ekinetic,spg_pos,fdos_pos
     close(222)

  nlmin_old=nlmin
!  CPUcheck=.false.

!Outer (hopping) loop
call yaml_sequence_open('Hopping steps')
!Hopping_loop: do
1000 continue
     if (nlmin >= nlminx) then
        goto 3000
     endif
     call yaml_sequence(advance='no')
!Energy has reached target eref and global minimum is presumably found
!!!     if (re_sm <= 1.d-3) then
!!!        write(*,*)'# success: relative energy < 0.001'
!!!        goto 3000
!!!     endif

 

5555 continue


!If a file named CPUlimit exists goto end
     file_exists=.false.
     INQUIRE(FILE="CPUlimit", EXIST=file_exists)
     if(file_exists) then
         CPUlimit=0
         call system_clock(count=clock_cur)     !Current Timer
         open(unit=78,file="CPUlimit")
         read(78,*,end=1212) CPUlimit
1212 continue
         do while(clock_cur-clock_start.lt.0)
                  clock_cur=clock_cur+clock_max
         enddo
         elapsed=int(real(clock_cur-clock_start)/real(clock_rate))
         write(*,'(a,1x,i7,1x,i7)') ' # CPUlimit found, elapsed, limit',elapsed,CPUlimit
         if(elapsed.ge.CPUlimit) then
         write(*,'(a)')' # CPUlimit reached, exiting!'
         goto 3000
         endif
         close(78)
     endif 


!!CPU time check to be implemented later
!!check whether CPU time exceeded
!     tleft=1.d100
!        if(iproc==0 .and. CPUcheck)then
!        open(unit=55,file='CPUlimit_global',status='unknown')
!        read(55,*,end=555) cpulimit ; cpulimit=cpulimit*3600
!        write(*,'(a,i5,i3,2(1x,e9.2))') 'iproc,nlmin,tcpu2-tcpu1,cpulimit',iproc,nlmin,tcpu2-tcpu1,cpulimit
!        close(55)
!        call cpu_time(tcpu2)
!        tleft=cpulimit-(tcpu2-tcpu1)
!     end if
555    continue
!       call MPI_BCAST(tleft,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!       if (tleft < 0.d0) then
!       write(*,*) 'CPU time exceeded',tleft
!       goto 3000
!       endif
!          CPUcheck=.true.


!Update the number of hops being performed, the next hop is nhop 
  nhop=nhop+1
  write(fn5,'(i5.5)') nhop
  if(parini%verb.ge.2) folder="data_hop_"//fn5//"/"
  if(parini%verb.ge.2) call system("mkdir "//trim(folder))

!Copy all variables to the tomporary variables wpos_...
     do iat=1,parini%nat
        wpos_red(1,iat)=pos_red(1,iat)
        wpos_red(2,iat)=pos_red(2,iat)
        wpos_red(3,iat)=pos_red(3,iat)
     enddo
        wpos_latvec=pos_latvec
        ent_wpos=ent_pos
        e_wpos=e_pos
        wpos_fcart=pos_fcart        
        wpos_strten=pos_strten


!Here we call the softening routine
     if(confine==1) confine=2
     call init_vel(parini,parres,vel_in,vel_lat_in,vel_vol_in,wpos_latvec,wpos_red,parini%bmass*amu_emass,ekinetic,parini%nsoften_minhopp,folder)

!Call molecular dynamics 
     escape=escape+1.d0
     iprec=2
     if((all(parini%fixlat(1:6)).and..not.parini%fixlat(7)).or.parini%bc==2) then
     call MD_fixlat(parini,parres,wpos_latvec,wpos_red,wpos_fcart,wpos_strten,vel_in,e_wpos,iprec,counter,folder)
     else
     call MD_MHM(parini,parres,wpos_latvec,wpos_red,wpos_fcart,wpos_strten,vel_in,vel_lat_in,vel_vol_in,e_wpos,iprec,counter,folder)
!       call MD_ANDERSEN_MHM(parini,wpos_latvec,wpos_red,wpos_fcart,wpos_strten,vel_in,vel_lat_in,0.003d0,e_wpos,iprec,counter)
     endif

!Make the cell "good" before entering GEOPT
!Put atoms back into cell after MD
if(.not.(any(parini%fixlat).or.any(parini%fixat).or.confine.ge.1)) call correct_latvec(wpos_latvec,wpos_red,parini%nat,parini%correctalg,latvec_io)

!Convert to enthalpy
     call get_enthalpy(wpos_latvec,e_wpos,parini%target_pressure_habohr,ent_wpos)

!Show the energy and enthalpies after exiting MD
    call yaml_mapping_open('After exiting MD',flow=.true.)
    call yaml_map('energy',e_wpos)
    call yaml_map('enthalpy',ent_wpos)
    call yaml_mapping_close()
    ! write(*,'(a,2(1x,es25.15))') " # After exiting MD, we have energy, enthalpy:    ",e_wpos,ent_wpos
!     write(67,'(a,2(1x,es25.15))') "After exiting MD, we have energy, enthalpy:    ",e_wpos,ent_wpos

!Update MD counter
     count_md=count_md+counter     
     call yaml_map('Counter of MD updated',int(count_md))
     !write(*,'(a,i7)') " # Counter of MD updated: ", int(count_md)

!Keep track of kinetic energy as average, not yet decided it this should not rather be the temperature. Not yet implemented
     av_ekinetic=av_ekinetic+ekinetic

!Call geometry optimizer after MD
  vel_in=0.d0
  vel_lat_in=0.d0
   iprec=1
   if(parini%geopt_ext) then
     call geopt_external(parini,wpos_latvec,wpos_red,wpos_fcart,wpos_strten,e_wpos,iprec,parres%ka,parres%kb,parres%kc,counter)
   else
      if(confine==2) confine=1
      if(parini%paropt_geopt%approach=="RBFGS") call GEOPT_RBFGS_MHM(parini,parres,wpos_latvec,wpos_red,wpos_fcart,&
                                     &wpos_strten,e_wpos,iprec,counter,folder)
      if(parini%paropt_geopt%approach=="MBFGS") call GEOPT_MBFGS_MHM(parini,parres,wpos_latvec,wpos_red,wpos_fcart,&
                                     &wpos_strten,e_wpos,iprec,counter,folder)
      if(parini%paropt_geopt%approach=="FIRE")  call GEOPT_FIRE_MHM(parini,parres,wpos_latvec,wpos_red,wpos_fcart,&
                                     &wpos_strten,vel_in,vel_lat_in,vel_vol_in,e_wpos,iprec,counter,folder)
      if(parini%paropt_geopt%approach=="SQNM")   call    GEOPT_SQNM(parini,parres,wpos_latvec,wpos_red,wpos_fcart,&
                                     &wpos_strten,e_wpos,iprec,counter,folder)
      if(parini%paropt_geopt%approach=="QBFGS")  call   GEOPT_QBFGS(parini,parres,wpos_latvec,wpos_red,wpos_fcart,&
                                     &wpos_strten,e_wpos,iprec,counter,folder)
      if(parini%paropt_geopt%approach=="SD")     call      GEOPT_SD(parini,parres,wpos_latvec,wpos_red,wpos_fcart,&
                                     &wpos_strten,e_wpos,iprec,counter,folder)
   endif

!Update GEOPT counter
     count_geopt=count_geopt+counter     
     call yaml_map('Counter of GEOPT updated',int(count_geopt))
     !write(*,'(a,i7)') " # Counter of GEOPT updated: ", int(count_geopt)

!Check for symmetry
   if(parini%findsym) then
      call find_symmetry(parini,parini%nat,wpos_red,wpos_latvec,parini%typat_global,tolmin,tolmax,ntol,spgtol_wpos,spgint)
      spg_wpos=spgint
   else
      spg_wpos=0
      spgtol_wpos=0.d0
   endif

!Get FDOS
  if(parini%finddos) then
     iprec=1;getwfk=.true.
     call get_dos(parini,parres,wpos_latvec,wpos_red,efermi,fdos_wpos,iprec,getwfk)
     write(*,'(a,es15.7,es15.7)') " # FDOS found: ",fdos_wpos,efermi
  endif


!Re-read parameters
!   write(*,'(a)')  " # Re-reading params.in"
!   call read_params()
  !call params_read(parini)

!Convert to enthalpy
  call get_enthalpy(wpos_latvec,e_wpos,parini%target_pressure_habohr,ent_wpos)

!Put atoms back into cell after MD+GEOPT
  !call backtocell(nat,wpos_latvec,wpos_red)

!Show the energy and enthalpies after exiting GEOPT 
    call yaml_mapping_open('After exiting GEOPT',flow=.true.)
    call yaml_map('e_wpos',e_wpos,fmt='(es25.15)')
    call yaml_map('ent_wpos',ent_wpos,fmt='(es25.15)')
    call yaml_mapping_close()
  !write(*,'(a,2(1x,es25.15))') " # After exiting GEOPT, we have energy, enthalpy: ",e_wpos,ent_wpos
!     write(67,'(a,2(1x,es25.15))') "After exiting GEOPT, we have energy, enthalpy: ",e_wpos,ent_wpos


!Write out the relaxed configuration
  if(parini%verb.gt.0) then
        filename='poslocm_'//fn5//'.ascii'
        call write_atomic_file_ascii(parini,filename,parini%nat,units,wpos_red,wpos_latvec,wpos_fcart,wpos_strten,parini%char_type,&
        &parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,e_wpos,parini%target_pressure_habohr,ent_wpos-eref,e_wpos-eref)
  call dist_latvec2ang(dist_ang,wpos_latvec,pi)
     call yaml_mapping_open('Wrote file',flow=.true.)
     call yaml_map('filename',trim(filename))
     call yaml_map('a',dist_ang(1),fmt='(es15.7)')
     call yaml_map('b',dist_ang(2),fmt='(es15.7)')
     call yaml_map('c',dist_ang(3),fmt='(es15.7)')
     call yaml_map('alpha',dist_ang(4),fmt='(es15.7)')
     call yaml_map('beta',dist_ang(5),fmt='(es15.7)')
     call yaml_map('gamma',dist_ang(6),fmt='(es15.7)')
     call yaml_mapping_close()
  !write(*,'(a,a,a,6(1x,es15.7))') " # Wrote file ",trim(filename),", with a,b,c,alpha,beta,gamma: ",dist_ang
  endif

!Get the fingerprint
 call get_fp(parini,fp_len,wpos_red,wpos_latvec,fp_wpos)

!Compute relative energy and enthalpy
!!  re_wpos=round(e_wpos-eref,accur)
!!  rent_wpos=round(ent_wpos-eref,accur)

!Write some informations about current state
!!  write(67,'(a,i3,4(1x,1pe14.7))')  &
!!       'nlmin,ent_wpos,ent_pos,rent_wpos,rent_pos', nlmin,ent_wpos,ent_pos,rent_wpos,rent_pos
!!  write(67,'(a,i3,4(1x,1pe14.7))')  &
!!       'nlmin,e_wpos,  e_pos,  re_wpos,  re_pos  ', nlmin,e_wpos,e_pos,re_wpos,re_pos
!!  write(*,'(a,i3,4(1x,1pe14.7))')  &
!!       ' # nlmin,ent_wpos,ent_pos,rent_wpos,rent_pos', nlmin,ent_wpos,ent_pos,rent_wpos,rent_pos
!!  write(*,'(a,i3,4(1x,1pe14.7))')  &
!!       ' # nlmin,e_wpos,  e_pos,  re_wpos,  re_pos  ', nlmin,e_wpos,e_pos,re_wpos,re_pos
  !write(67,'(a,i3,2(1x,1pe14.7))')  &
  !     'nlmin,ent_wpos,ent_pos ', nlmin,ent_wpos,ent_pos
  !write(67,'(a,i3,2(1x,1pe14.7))')  &
  !     'nlmin,e_wpos,e_pos      ', nlmin,e_wpos,e_pos
  call yaml_mapping_open('some results',flow=.true.)
  call yaml_map('nlmin',nlmin)
  call yaml_map('ent_wpos',ent_wpos)
  call yaml_map('ent_pos',ent_pos)
  call yaml_map('e_wpos',e_wpos)
  call yaml_map('e_pos',e_pos)
  call yaml_mapping_close()
  !write(*,'(a,i3,2(1x,1pe14.7))')  &
  !     ' # nlmin,ent_wpos,ent_pos', nlmin,ent_wpos,ent_pos
  !write(*,'(a,i3,2(1x,1pe14.7))')  &
  !     ' # nlmin,e_wpos,e_pos    ', nlmin,e_wpos,e_pos

!Not escaped, here it is checked to return to another MD if necessary
!NEW
 if (abs(ent_wpos-ent_pos).lt.ent_delta) then
!   call fpdistance(nid,wfp,fp,d)
   call get_fp_distance(parini,fp_len,fp_pos,fp_wpos,fp_dist)
   write(*,'(a,2(e11.4))') ' # ID: checking escape: enthalpy difference < ent_delta, fp_dist ',e_wpos-e_pos,fp_dist
   if (fp_dist.lt.fp_delta) then ! not escaped
     escape_sam=escape_sam+1.d0
     esep=esep+(ent_pos-ent_wpos)**2
     fp_sep=max(fp_sep,fp_dist)
     ekinetic=max(min(ekinetic_max,ekinetic*beta1),100.d0)       !This may still be changed into the temperature. Or we could use T=2/(3*kb)*ekinetic

!Update the temperature within dataset 2, belonging to MD
     call wtioput(ediff,ekinetic,ekinetic_max,parini%nsoften_minhopp)
     open(unit=222,file='global.mon',status='unknown',position='append')
     write(222,'(i10,(1x,f10.0),1x,1pe21.14,2(1x,1pe10.3),i5,1x,1pe10.3,3(1x,0pf5.2),a)')  &
          nhop,escape,ent_wpos-eref,ediff,ekinetic,spg_wpos,fdos_wpos, &
          escape_sam/escape,escape_old/escape,escape_new/escape,'   S'
     close(222)
    call yaml_mapping_open('global',flow=.true.)
    call yaml_map('escape',escape,fmt='(f10.0)')
    call yaml_map('de',ent_wpos-eref,fmt='(es21.14)')
    call yaml_map('ediff',ediff,fmt='(es10.3)')
    call yaml_map('ekinetic',ekinetic,fmt='(es10.3)')
    call yaml_map('spg_pos',spg_wpos,fmt='(i8)')
    call yaml_map('fdos_pos',fdos_wpos,fmt='(es10.3)')
    call yaml_map('ratio_sam',escape_sam/escape,fmt='(f5.2)')
    call yaml_map('ratio_old',escape_old/escape,fmt='(f5.2)')
    call yaml_map('ratio_new',escape_new/escape,fmt='(f5.2)')
    call yaml_map('status','S')
    call yaml_mapping_close()
     !write(*,'(a,(1x,f10.0),1x,1pe21.14,2(1x,1pe10.3),i5,1x,1pe10.3,3(1x,0pf5.2),a)') " #global: ", &
     !     escape,ent_wpos-eref,ediff,ekinetic,spg_wpos,fdos_wpos, &
     !     escape_sam/escape,escape_old/escape,escape_new/escape,'   S'
     write(*,'(a)')' # no escape from current minimum.'

     if(parini%auto_mdmin) parres%mdmin   = min(parini%mdmin_max, parres%mdmin + 1)  ! MALM
     call yaml_map('nsoften',parini%nsoften_minhopp)
     call yaml_map('mdmin',parres%mdmin)
     !write(*,*) "# nsoften, mdmin: ", parini%nsoften_minhopp, parres%mdmin ! MALM

     goto 5555
   endif
 endif
!NEW

!!!  if (rent_pos == rent_wpos) then
!!!     escape_sam=escape_sam+1.d0
!!!     esep=esep+(ent_pos-ent_wpos)**2
!!!     ekinetic=min(ekinetic_max,ekinetic*beta1)       !This may still be changed into the temperature. Or we could use T=2/(3*kb)*ekinetic
!!!
!!!!Update the temperature within dataset 2, belonging to MD
!!!     call wtioput(ediff,ekinetic,ekinetic_max,nsoften)
!!!     open(unit=222,file='global.mon',status='unknown',position='append')
!!!     write(222,'(i10,(1x,f10.0),1x,1pe21.14,2(1x,1pe10.3),i5,1x,1pe10.3,3(1x,0pf5.2),a)')  &
!!!          nconjgr,escape,ent_wpos-eref,ediff,ekinetic,int(spg_wpos),fdos_wpos, &
!!!          escape_sam/escape,escape_old/escape,escape_new/escape,'   S'
!!!     close(222)
!!!     write(*,'(a,(1x,f10.0),1x,1pe21.14,2(1x,1pe10.3),i5,1x,1pe10.3,3(1x,0pf5.2),a)') " #global: ", &
!!!          escape,ent_wpos-eref,ediff,ekinetic,int(spg_wpos),fdos_wpos, &
!!!          escape_sam/escape,escape_old/escape,escape_new/escape,'   S'
!!!     write(*,'(a)')' # no escape from current minimum.'
!!!
!!!     if(parini%auto_mdmin) parres%mdmin   = min(mdmin_max, parres%mdmin + 1)  ! MALM
!!!     write(*,*) "# nsoften, mdmin: ", nsoften, parres%mdmin ! MALM
!!!
!!!     goto 5555
!!!  endif



!Continue since escaped
!Check whether new minimum
  call identical(parini,nlminx,nlmin,fp_method,fp_len,ent_wpos,fp_wpos,ent_arr,fp_arr,&
                 &ent_delta,fp_delta,newmin,kid,fp_dist_min,k_e_wpos,n_unique,n_nonuni,lid,nid)
  call yaml_map('k_e_wpos,kid',k_e_wpos)
  call yaml_map('kid',kid)
  !write(*,*) "k_e_wpos,kid: ",k_e_wpos,kid
!!!  call hunt(earr(1,1),min(nlmin,nlminx),rent_wpos,k_e_wpos)
  if (.not.newmin) then                    !The case the structure was previously found
!!!!     if(findsym) call replace(nlminx,nlmin,fp_len,nat,kid,e_wpos,ent_wpos,fp_wpos,wpos_red,&
!!!!       &wpos_latvec,spg_wpos,spgtol_wpos,fdos_wpos,&
!!!!       &e_arr,ent_arr,fp_arr,pl_arr,lat_arr,spg_arr,spgtol_arr,dos_arr,ct_arr,findsym)
!!  if (rent_wpos == earr(k_e_wpos,1)) then                             !The case the structure was previously found
!     k_e_wpos=k_e_wpos+1
     !write(67,'(a,i3,i4,1x,1pe14.7)')  &
     !       ' Revisited: nlmin,k_e_wpos,ent_wpos=earr',nlmin,k_e_wpos,ent_wpos
     write(*,'(a,i3,i4,1x,1pe14.7)')  &
          ' # Revisited: nlmin,k_e_wpos,ent_wpos=earr',nlmin,k_e_wpos,ent_wpos
!!     newmin=.false.
     escape_old=escape_old+1.d0
     ct_arr(kid)=ct_arr(kid)+1                               !Increasing the number of visits
!!     earr(k_e_wpos,3)=earr(k_e_wpos,3)+1.d0                           !Increasing the number of visits
     ekinetic=max(min(ekinetic_max,ekinetic*beta2),100.d0)
     nvisit=ct_arr(kid)

!     if(parini%auto_mdmin) parres%mdmin   = min(mdmin_max, parres%mdmin + 1)  ! MALM
!Update the temperature within dataset 2, belonging to MD
  else                                                                !A new minimum was found
! write intermediate results
!!     newmin=.true.
     escape_new=escape_new+1.d0
     ekinetic=max(ekinetic*beta3,100.d0)

!Determine energy separation between distinct minima
     if (k_e_wpos+1.le.nlminx)  then
!!        tt=min(rent_wpos-earr(k_e_wpos,1),earr(k_e_wpos+1,1)-rent_wpos)
        tt=min(max(ent_wpos-ent_arr(max(k_e_wpos,1)),0.d0),ent_arr(k_e_wpos+1)-ent_wpos)
!        if (tt.gt. accur*(1.1d0)) egap=min(egap,tt)
        egap=min(egap,tt)
     endif
     call wtioput(ediff,ekinetic,ekinetic_max,parini%nsoften_minhopp)
     nlmin=nlmin+1
!!       call insert(nlminx,nlmin,k_e_wpos,rent_wpos,re_wpos,earr(0,1))
     call insert(nlminx,nlmin,fp_len,parini%nat,k_e_wpos,e_wpos,ent_wpos,fp_wpos,wpos_red,wpos_latvec,wpos_fcart,wpos_strten,&
          &spg_wpos,spgtol_wpos,fdos_wpos,e_arr,ent_arr,fp_arr,pl_arr,lat_arr,f_arr,str_arr,spg_arr,spgtol_arr,dos_arr,ct_arr)
!!     call winter(parini,nat,units,re_pos,rent_pos,e_pos,pos_red,pos_latvec,npminx,nlminx,nlmin,npmin,accur, & 
!!          &earr,elocmin,poslocmin,latlocmin,eref,ediff,ekinetic,ekinetic_max,dt,nsoften,char_type,ntypat,typat,&
!!          &fixat,fixlat,target_pressure_habohr)
     k_e_wpos=k_e_wpos+1
     if (k_e_wpos .gt. nlminx .or. k_e_wpos .lt. 1) stop "k_e_wpos out of bounds"
     if(modulo(nlmin,nwrite_inter)==0) then 
     call winter(parini,parini%nat,units,ent_pos,e_pos,pos_red,pos_latvec,pos_fcart,pos_strten,nlminx,nlmin,npminx,& 
          &ent_arr,e_arr,ct_arr,spg_arr,spgtol_arr,dos_arr,pl_arr,lat_arr,f_arr,str_arr,fp_arr,fp_len,ent_delta,fp_delta,& 
          &eref,ediff,ekinetic,ekinetic_max,dt,parini%nsoften_minhopp,parini%char_type,parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,parini%target_pressure_habohr)
     endif
     nvisit=1
!!     npmin=npmin+1
!!     call save_low_conf(nat,npmin,npminx,rent_wpos,re_wpos,wpos_red,wpos_latvec,spg_wpos,&
!!     &spgtol_wpos,fdos_wpos,elocmin,poslocmin,latlocmin)

!     if(parini%auto_mdmin) parres%mdmin   = max(parini%mdmin_min, parres%mdmin - 1)  ! MALM
  endif
!  write(*,*) "# nsoften, mdmin: ", nsoften, parres%mdmin ! MALM

!This part is to finally accept the configutaion for the next run, I guess
  if (ent_wpos.lt.ent_hop) then
    ent_hop=ent_wpos
    e_hop=e_wpos
!    rent_hop=rent_wpos
!    re_hop=re_wpos
    spg_hop=spg_wpos
    spgtol_hop=spgtol_wpos
    fdos_hop=fdos_wpos
    fp_hop=fp_wpos

    nvisit=ct_arr(k_e_wpos)
    do iat=1,parini%nat
      poshop(1,iat)=wpos_red(1,iat) ; poshop(2,iat)=wpos_red(2,iat) ; poshop(3,iat)=wpos_red(3,iat)
    enddo
    lathop=wpos_latvec
    poshop_fcart=wpos_fcart
    poshop_strten=wpos_strten
  endif

!MASTER: Monte Carlo step for local minima hopping
  av_ediff=av_ediff+ediff
  if (ent_hop-ent_pos.lt.ediff) then
!Local minima accepted  
   if(parini%auto_mdmin.and.newmin)      parres%mdmin   = max(parini%mdmin_min, parres%mdmin - 1)  ! MALM
   if(parini%auto_mdmin.and..not.newmin) parres%mdmin   = min(parini%mdmin_max, parres%mdmin + 1)  ! MALM
   call yaml_map('nsoften',parini%nsoften_minhopp)
   call yaml_map('mdmin',parres%mdmin)
   !write(*,*) "# nsoften, mdmin: ", parini%nsoften_minhopp, parres%mdmin ! MALM
   accepted=accepted+1.d0
   e_pos=e_hop
   ent_pos=ent_hop
!   re_pos=re_hop
!   rent_pos=rent_hop
   spg_pos=spg_hop
   spgtol_pos=spgtol_hop
   fdos_pos=fdos_hop
   fp_pos=fp_hop
     do iat=1,parini%nat
        pos_red(1,iat)=poshop(1,iat)
        pos_red(2,iat)=poshop(2,iat)
        pos_red(3,iat)=poshop(3,iat)
     enddo
     pos_latvec=lathop
     pos_fcart=poshop_fcart
     pos_strten=poshop_strten
     call wtioput(ediff,ekinetic,ekinetic_max,parini%nsoften_minhopp)
      if (abs(ent_wpos-ent_hop).lt.1.d-10) then
         open(unit=222,file='global.mon',status='unknown',position='append')
         write(222,'(i10,(1x,f10.0),1x,1pe21.14,2(1x,1pe10.3),i5,1x,1pe10.3,3(1x,0pf5.2),l3,a,i5)')  &
          nhop,escape,ent_hop-eref,ediff,ekinetic,spg_hop,fdos_hop, &
          escape_sam/escape,escape_old/escape,escape_new/escape,newmin,' A ', nvisit
         close(222)
        call yaml_mapping_open('global',flow=.true.)
        call yaml_map('escape',escape,fmt='(f10.0)')
        call yaml_map('de',ent_hop-eref,fmt='(es21.14)')
        call yaml_map('ediff',ediff,fmt='(es10.3)')
        call yaml_map('ekinetic',ekinetic,fmt='(es10.3)')
        call yaml_map('spg_pos',spg_hop,fmt='(i8)')
        call yaml_map('fdos_pos',fdos_hop,fmt='(es10.3)')
        call yaml_map('ratio_sam',escape_sam/escape,fmt='(f5.2)')
        call yaml_map('ratio_old',escape_old/escape,fmt='(f5.2)')
        call yaml_map('ratio_new',escape_new/escape,fmt='(f5.2)')
        call yaml_map('escaped',newmin)
        call yaml_map('status','S')
        call yaml_map('nvisit',nvisit)
        call yaml_mapping_close()
         !write(*,'(a,(1x,f10.0),1x,1pe21.14,2(1x,1pe10.3),i5,1x,1pe10.3,3(1x,0pf5.2),l3,a,i5)') " #global: ", &
         ! escape,ent_hop-eref,ediff,ekinetic,spg_hop,fdos_hop, &
         ! escape_sam/escape,escape_old/escape,escape_new/escape,newmin,' A ', nvisit
      else
         open(unit=222,file='global.mon',status='unknown',position='append')
         write(222,'(i10,1x,f10.0,1x,1pe21.14,2(1x,1pe10.3),i5,1x,1pe10.3,3(1x,0pf5.2),a,i5)')  &
              nhop,escape,ent_wpos-eref,ediff,ekinetic,spg_wpos,fdos_wpos,&
              escape_sam/escape,escape_old/escape,escape_new/escape,'   I  ',int(ct_arr(k_e_wpos))
         write(222,'(i10,1x,f10.0,1x,1pe21.14,2(1x,1pe10.3),i5,1x,1pe10.3,3(1x,0pf5.2),a,i5)')  &
              nhop,escape,ent_hop-eref,ediff,ekinetic,spg_hop,fdos_hop, &
              escape_sam/escape,escape_old/escape,escape_new/escape,'   A  ',nvisit
         close(222)
        call yaml_mapping_open('global',flow=.true.)
        call yaml_map('escape',escape,fmt='(f10.0)')
        call yaml_map('de',ent_wpos-eref,fmt='(es21.14)')
        call yaml_map('ediff',ediff,fmt='(es10.3)')
        call yaml_map('ekinetic',ekinetic,fmt='(es10.3)')
        call yaml_map('spg_pos',spg_wpos,fmt='(i8)')
        call yaml_map('fdos_pos',fdos_wpos,fmt='(es10.3)')
        call yaml_map('ratio_sam',escape_sam/escape,fmt='(f5.2)')
        call yaml_map('ratio_old',escape_old/escape,fmt='(f5.2)')
        call yaml_map('ratio_new',escape_new/escape,fmt='(f5.2)')
        !call yaml_map('escaped',newmin)
        call yaml_map('status','I')
        call yaml_map('nvisit',int(ct_arr(k_e_wpos)))
        call yaml_mapping_close()
        ! write(*,'(a,1x,f10.0,1x,1pe21.14,2(1x,1pe10.3),i5,1x,1pe10.3,3(1x,0pf5.2),a,i5)') " #global: ", &
        !      escape,ent_wpos-eref,ediff,ekinetic,spg_wpos,fdos_wpos, &
        !      escape_sam/escape,escape_old/escape,escape_new/escape,'   I  ',int(ct_arr(k_e_wpos))
        call yaml_mapping_open('global',flow=.true.)
        call yaml_map('escape',escape,fmt='(f10.0)')
        call yaml_map('de',ent_hop-eref,fmt='(es21.14)')
        call yaml_map('ediff',ediff,fmt='(es10.3)')
        call yaml_map('ekinetic',ekinetic,fmt='(es10.3)')
        call yaml_map('spg_pos',spg_hop,fmt='(i8)')
        call yaml_map('fdos_pos',fdos_hop,fmt='(es10.3)')
        call yaml_map('ratio_sam',escape_sam/escape,fmt='(f5.2)')
        call yaml_map('ratio_old',escape_old/escape,fmt='(f5.2)')
        call yaml_map('ratio_new',escape_new/escape,fmt='(f5.2)')
        !call yaml_map('escaped',newmin)
        call yaml_map('status','A')
        call yaml_map('nvisit',nvisit)
        call yaml_mapping_close()
        !write(*,'(a,1x,f10.0,1x,1pe21.14,2(1x,1pe10.3),i5,1x,1pe10.3,3(1x,0pf5.2),a,i5)') " #global: ", &
        !      escape,ent_hop-eref,ediff,ekinetic,spg_hop,fdos_hop, &
        !      escape_sam/escape,escape_old/escape,escape_new/escape,'   A  ',nvisit
      endif
!      endif
      ent_hop=1.d100
      ediff=max(ediff*alpha1,1.d-4)
! write intermediate results
    call yaml_comment('WINTER')
      !write(*,*) 'WINTER'
!!      call winter(parini,nat,units,re_pos,rent_pos,e_pos,pos_red,pos_latvec,npminx,nlminx,nlmin,npmin,accur, & 
!!      &earr,elocmin,poslocmin,latlocmin,eref,ediff,ekinetic,ekinetic_max,dt,nsoften,char_type,ntypat,typat,&
!!      &fixat,fixlat,target_pressure_habohr)
     if(modulo(nlmin,nwrite_inter)==0) then 
     call winter(parini,parini%nat,units,ent_pos,e_pos,pos_red,pos_latvec,pos_fcart,pos_strten,nlminx,nlmin,npminx,& 
          &ent_arr,e_arr,ct_arr,spg_arr,spgtol_arr,dos_arr,pl_arr,lat_arr,f_arr,str_arr,fp_arr,fp_len,ent_delta,fp_delta,& 
          &eref,ediff,ekinetic,ekinetic_max,dt,parini%nsoften_minhopp,parini%char_type,parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,parini%target_pressure_habohr)
     endif
     goto 1000
  else
!local minima rejected
!    if(parini%auto_mdmin) parres%mdmin   = min(mdmin_max, parres%mdmin + 1)  ! MALM
    call yaml_map('nsoften',parini%nsoften_minhopp)
    call yaml_map('mdmin',parres%mdmin)
     !write(*,*) "# nsoften, mdmin: ", parini%nsoften_minhopp, parres%mdmin ! MALM
     open(unit=222,file='global.mon',status='unknown',position='append')
     write(222,'(i10,(1x,f10.0),1x,1pe21.14,2(1x,1pe10.3),i5,1x,1pe10.3,3(1x,0pf5.2),l3,a,i5)')  &
          nhop,escape,ent_wpos-eref,ediff,ekinetic,spg_wpos,fdos_wpos,&
          escape_sam/escape,escape_old/escape,escape_new/escape,newmin,' R ', int(ct_arr(k_e_wpos))
     close(222)
     write(*,'(a,(1x,f10.0),1x,1pe21.14,2(1x,1pe10.3),i5,1x,1pe10.3,3(1x,0pf5.2),l3,a,i5)') " #global: ",  &
          escape,ent_wpos-eref,ediff,ekinetic,spg_wpos,fdos_wpos,&
          escape_sam/escape,escape_old/escape,escape_new/escape,newmin,' R ', int(ct_arr(k_e_wpos))
     write(*,'(a,1pe21.14)')' # rejected: ew-e>ediff ',ent_wpos-ent_pos

     rejected=rejected+1.d0
     ediff=max(ediff*alpha2,1.d-4)
     call wtioput(ediff,ekinetic,ekinetic_max,parini%nsoften_minhopp)
     goto 1000
!------------------------------------------------------------
  endif

!end do hopping_loop
3000 continue   !This is the really end. Collect data
call yaml_sequence_close()

!Print some informations and rewrite restart information
     !write(67,*) 'writing final results'
     !write(*,*) ' # writing final results'
     call yaml_comment('writing final results',hfill='~')
     !write(67,*) ' found in total ',nlmin,' minima'
     !write(67,*) ' Accepeted ',int(accepted),' minima'
     !write(*,*) ' # found in total ',nlmin,' minima'
     call yaml_map('total minima found',nlmin)
     call yaml_map('minima accepeted',int(accepted))
     !write(*,*) ' # Accepeted ',int(accepted),' minima'
!!     call winter(parini,nat,units,re_pos,rent_pos,e_pos,pos_red,pos_latvec,npminx,nlminx,nlmin,npmin,accur, & 
!!      earr,elocmin,poslocmin,latlocmin,eref,ediff,ekinetic,ekinetic_max,dt,nsoften,char_type,ntypat,&
!!      typat,fixat,fixlat,target_pressure_habohr)
     call winter(parini,parini%nat,units,ent_pos,e_pos,pos_red,pos_latvec,pos_fcart,pos_strten,nlminx,nlmin,npminx,& 
          &ent_arr,e_arr,ct_arr,spg_arr,spgtol_arr,dos_arr,pl_arr,lat_arr,f_arr,str_arr,fp_arr,fp_len,ent_delta,fp_delta,& 
          &eref,ediff,ekinetic,ekinetic_max,dt,parini%nsoften_minhopp,parini%char_type,parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,parini%target_pressure_habohr)


!Print ratios from all the global counters
    call yaml_mapping_open('about minima hopping params')
     if(escape.gt.0.d0) then
        call yaml_mapping_open('ratios',flow=.true.)
        call yaml_map('same',escape_sam/escape,fmt='(es10.3)')
        call yaml_map('old',escape_old/escape,fmt='(es10.3)')
        call yaml_map('new',escape_new/escape,fmt='(es10.3)')
        call yaml_mapping_close()
        !write(*,'(a,3(1x,1pe10.3))') ' # ratio stuck,same,old,new', &
        !  escape_sam/escape,escape_old/escape,escape_new/escape
     endif
     call yaml_map('ratio_acc',accepted/(accepted+rejected),fmt='(es10.3)')
     call yaml_map('ratio_rej',rejected/(accepted+rejected),fmt='(es10.3)')
     call yaml_map('count_md',count_md,fmt='(f12.1))')
     call yaml_map('count_geopt',count_geopt,fmt='(f12.1))')
     !write(*,'(a,2(1x,1pe10.3))') ' # ratio acc,rej',accepted/(accepted+rejected),rejected/(accepted+rejected)
     !write(*,'(a,2(1x,f12.1))')   ' # count_md,count_geopt',count_md,count_geopt
     if(escape.gt.0.d0) then
         call yaml_map('average ediff',av_ediff/(accepted+rejected),fmt='(es10.3)')
         call yaml_map('average ekinetic',av_ekinetic/escape,fmt='(es10.3)')
         !write(*,'(a,2(1x,1pe10.3))') &
         ! ' # average ediff, ekinetic',av_ediff/(accepted+rejected),av_ekinetic/escape
     endif
     call yaml_map('number of configurations for which atoms escaped',nputback)
     !write(*,'(a,1x,i8)') ' # number of configurations for which atoms escaped ',nputback
     call yaml_mapping_close()

     tt=0.d0
     ss=0.d0
     do i=1,nlmin
        tt=max(tt,real(ct_arr(i),8))
        ss=ss+real(ct_arr(i),8)
     enddo
     call yaml_mapping_open('about minima',flow=.true.)
     call yaml_map('most frequent visits',tt,fmt='(f8.0)')
     call yaml_map('avg. num. visits per minimum',ss/nlmin,fmt='(es10.3)')
     call yaml_map('minimum energy separation',egap,fmt='(es9.2)')
     !write(*,'(a,f8.0)') ' # most frequent visits ',tt
     !write(*,'(a,1pe10.3)') ' # av. numb. visits per minimum',ss/nlmin
     !write(*,'(a,e9.2)') ' # minimum energy separation between presumably different configurations',egap
     if (escape_sam.gt.0) then
        esep=sqrt(esep/escape_sam)
        call yaml_map('average energy separation',esep,fmt='(es9.2)')
        !write(*,'(a,e9.2)') ' # average energy separation between presumably identical configurations',esep
     endif
     call yaml_mapping_close()

!  close(67)


! call timein(tsec(1),tsec(2))
! tsec(1)=tsec(1)-cpui
! tsec(2)=tsec(2)-walli

  deallocate(parini%char_type)
  deallocate(e_arr,ent_arr,ct_arr,pl_arr,lat_arr,dos_arr,spgtol_arr,spg_arr,fp_arr)
3001 continue
!Close socket on slave side
if(trim(parini%potential_potential)=="msock") call socket_stop()
end subroutine task_minhocao
!contains

