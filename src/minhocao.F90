subroutine task_minhocao(parini,parres)
 use mod_interface
 use global
 use defs_basis
 !use cell_utils
 use interface_code
 use fingerprint
 use mod_parini, only: typ_parini
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
call print_logo() 

!Delete the file CPUlimit, if it exists
!  call system("rm -f CPUlimit")
call system_clock(count=clock_start)     !Start Timer
call system_clock(count_rate=clock_rate) !Find the time rate
call system_clock(count_max=clock_max)   !Find the time max


!Define the boundary condition: 1: periodic, 2:free, 3:surface/slab
  bc=1

!Verbosity
  parini%verb=1

!Define current directory
  folder=""

!Initialize auto logicals to false, otherwise it will read the params file and reset alpha_lat, alpha_at, mdmin_max, etc
  auto_soft=.false.
  parini%auto_mdmin=.false.
  parini%auto_dtion_md=.false.
  alpha_at=-1.d10
  alpha_lat=-1.d10 

!Initialize old kpt in history
  ka1=0;kb1=0;kc1=0;max_kpt=.false.;reuse_kpt=.false.

!Some parameter setup
!  code="lammps"          !Define what code to use for force evaluation
  siesta_kpt_mode=2          !Mode of generating automatic k-point mesh, only siesta
  vasp_kpt_mode=2            !Mode of generating automatic k-point mesh, only vasp
  abinit_kpt_mode=1          !Mode of generating automatic k-point mesh, only abinit
  correctalg=2               !Method to correct cell vectors when torn out of shape

!Define if the external optimizer should be used. Only available for:
!vasp
!siesta
!dftb
  parini%geopt_ext=.false.


!Unset fixed cell variables
  fixlat=.false.

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
  call params_read(parini)
  parres=parini
!  call params_echo()

!Read poscur into pos_red, which is equivalent to what pos was in the initial minima hopping
  allocate(pos_red(3,nat),pos_fcart(3,nat),xcart_tmp(3,nat),xcart_mol(3,nat),&
          &wpos_red(3,nat),wpos_fcart(3,nat),poshop_fcart(3,nat),vel_in(3,nat))

!Finally, reading in atomic positions
  file_exists=.false.
  filename="poscur.ascii"
  readfix=.true.;readfrag=.true.
  INQUIRE(FILE=trim(filename), EXIST=file_exists)
  if(file_exists) then
  write(*,*) "# Reading from ",trim(filename) 
  call read_atomic_file_ascii(filename,nat,units,pos_red,pos_latvec,pos_fcart,pos_strten,fixat,fixlat,&
       &readfix,fragarr,readfrag,tmp_enthalpy,tmp_energy)
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
     call read_atomic_file_poscar(filename,nat,units,pos_red,pos_latvec,pos_fcart,pos_strten,fixat,fixlat,&
          &readfix,fragarr,readfrag,tmp_enthalpy,tmp_energy)
     endif
  endif
  if(.not.file_exists) then
     stop "Provide either poscur.ascii or poscur.vasp"
  endif




!Write FIXLAT
  if(any(fixlat)) write(*,*) "# FIXLAT: a,b,c,alpha,beta,gamma,shape ", fixlat

!Initiallize Fingerprint
!  fp_method=11 !11: Oganov FP, 12: CALYPSO FP, 13: Modified CALYPSO, 21: molecular gaussian overlap, 22: molecular sprint
  call init_fp(fp_len,pos_latvec)
  allocate(fp_pos(fp_len),fp_wpos(fp_len),fp_hop(fp_len))

!Check correct assignment of FP method
  if(bc==1.and.(fp_method.lt.11.or.fp_method.gt.19)) stop "Incompatible fingerprint and boundary conditions"
  if(bc==2.and.(fp_method.lt.21.or.fp_method.gt.29)) stop "Incompatible fingerprint and boundary conditions"
  if(bc==3.and.(fp_method.lt.31.or.fp_method.gt.39)) stop "Incompatible fingerprint and boundary conditions"


!Assign atomic masses in electron units
 allocate(amass(nat))
 do iat=1,nat
   amass(iat)=amu_emass*amu(typat(iat))
 enddo

!Allocate and handle molecular crystal stuff
   call refragment(fragarr,nat)
   nmol=maxval(fragarr)
   allocate(fragsize(nmol),lhead(nmol),llist(nat))
   call make_linked_list(fragarr,fragsize,lhead,llist,nat,nmol)
!Here we will decompose the atomic position into a set of cm and quat 
   allocate(pos_cm(3,nmol),pos_quat(4,nmol),intens(3,3,nmol),inaxis(3,3,nmol),inprin(3,nmol),masstot(nmol))
   if(nmol.ne.nat) then 
       call init_cm_mol(pos_latvec,pos_red,xcart_mol,pos_cm,pos_quat,amass,masstot,intens,inprin,inaxis,lhead,llist,nat,nmol)
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
     if(use_confine) then
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
     call backtocell(nat,pos_latvec,pos_red)
     allocate(sym_pos(3,4*nat),sym_type(4*nat),sym_fcart(3,4*nat),sym_fixat(4*nat))
     sym_pos(:,1:nat)=pos_red
     sym_type(1:nat)=typat
     sym_fcart=0.d0
     sym_fixat=.false.
     pos_strten=0.d0
     fragarr=0
     call find_symmetry(parini,nat,pos_red,pos_latvec,typat,tolmin,tolmax,ntol,spgtol_pos,spg_pos)
     call spg_cell_refine(nat,sym_nat,4*nat,sym_pos,pos_latvec,sym_type,spgtol_pos,spg_pos)
     filename="sym_poscur.ascii"
     call write_atomic_file_ascii(parini,filename,sym_nat,units,sym_pos,pos_latvec,sym_fcart,pos_strten,char_type,ntypat,sym_type,sym_fixat,fixlat,&
     &0.d0,target_pressure_habohr,0.d0,0.d0)
     call spg_cell_primitive(sym_nat,sym_nat2,4*nat,sym_pos,pos_latvec,sym_type,spgtol_pos,spg_pos)
     filename="sym_prim_poscur.ascii"
     call write_atomic_file_ascii(parini,filename,sym_nat2,units,sym_pos,pos_latvec,sym_fcart,pos_strten,char_type,ntypat,sym_type,sym_fixat,fixlat,&
     &0.d0,target_pressure_habohr,0.d0,0.d0)
       goto 3001 
  endif 
#endif

!Check if pathintegral should be performed
  file_exists=.false.
  INQUIRE(FILE="pathparams.in", EXIST=file_exists)
     if(file_exists) then
       write(*,*)'# pathparams.in found, performing path integral starting from poscur.ascii'
     !Put all the atoms back into the cell
       call backtocell(nat,pos_latvec,pos_red)
     !Run path integral
       call pathintegral(parini,pos_latvec,pos_red)
       goto 3001 
  endif 


!Check if enthalpy sampling should be performed
  file_exists=.false.
  INQUIRE(FILE="enthalpyparams.in", EXIST=file_exists)
     if(file_exists) then
       write(*,*)'# enthalpyparams.in found, performing enthalpy scan starting from poscur.ascii'
     !Put all the atoms back into the cell
       call backtocell(nat,pos_latvec,pos_red)
     !Run path integral
       call enthalpyrelax(parini,pos_latvec,pos_red,tolmin,tolmax,ntol,findsym)
       goto 3001 
  endif 


!Check if volume sampling should be performed
  file_exists=.false.
  INQUIRE(FILE="varvol.in", EXIST=file_exists)
     if(file_exists) then
       write(*,*)'# varvol.in found, performing volume scan starting from poscur.ascii'
     !Put all the atoms back into the cell
       call backtocell(nat,pos_latvec,pos_red)
     !Run path integral
       call varvol(parini,pos_latvec,pos_red,tolmin,tolmax,ntol,findsym)
       goto 3001 
  endif 

!Check if a list of poslows should be relaxed
  file_exists=.false.
  INQUIRE(FILE="poslowrelax.in", EXIST=file_exists)
     if(file_exists) then
       write(*,*)'# poslowrelax.in found, performing relaxations of all poslows found &
                  &in current folder, but poscur.ascii must be provided for initiallization.'
     !Put all the atoms back into the cell
       call backtocell(nat,pos_latvec,pos_red)
     !Run path integral
       call poslowrelax(parini,pos_latvec,pos_red,tolmin,tolmax,ntol)
       goto 3001
  endif

!Check if rotations should be performed
  file_exists=.false.
  INQUIRE(FILE="rotate.in", EXIST=file_exists)
     if(file_exists) then
       write(*,*)'# rotate.in found, performing rotations of all poscur file found &
                  &in current folder, poscur.ascii must be provided for initiallization.'
     !Put all the atoms back into the cell
       call backtocell(nat,pos_latvec,pos_red)
     !Run path integral
       call rotate_like_crazy(parini,pos_latvec,pos_red,tolmin,tolmax,ntol)
       goto 3001
  endif



!Write initial parameters to global.out
  ratio=0.d0
!  open(unit=67,file='global.out')
  write(*,'(a,3(1x,1pe11.4))') ' # beta1,beta2,beta3',beta1,beta2,beta3
  write(*,'(a,2(1x,1pe11.4))') ' # alpha1,alpha2',alpha1,alpha2
  write(*,'(a,2(1x,1pe10.3))') ' # predicted fraction accepted, rejected', &
  & ratio/(1.d0+ratio), 1.d0/(1.d0+ratio)
  write(*,*) '# parres%mdmin',parres%mdmin
  

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
  write(*,'(a,2(es15.7))') " # Tolerance between structures: delta enthalpy, delta fingeprint: ", ent_delta,fp_delta 
  if (nlmin.gt.nlminx) stop 'nlmin>nlminx'
!  allocate(earr(0:nlminx,3),stat=i_stat) !This earr is modified from the original version: earr(i,1)=enthalpy,earr(i,2)=energy,earr(i,3)=visits
  allocate(e_arr(nlminx),ent_arr(nlminx),ct_arr(nlminx),pl_arr(3,nat,nlminx),f_arr(3,nat,nlminx),&
       &str_arr(6,nlminx),lat_arr(3,3,nlminx),dos_arr(nlminx),spgtol_arr(nlminx),spg_arr(nlminx),&
       &fp_arr(fp_len,nlminx),lid(nlminx))
!  earr=0.d0
!  earr(0,1)=-1.d100
  e_arr=0.d0;ent_arr=0.d0;dos_arr=0.d0
  if (nlmin == 0) then 
!     write(67,*) ' New run with nlminx=',nlminx
     write(*,*) '# New run with nlminx=',nlminx
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
  if(fixlat(7)) then
          parini%md_algo = 4
          write(*,'(A,i5)') " # Variable volume algorithm selected, switching to Anderson MD: parini%md_algo = ", parini%md_algo
  endif
  if(parini%md_algo==4) then
         fixlat(7)=.true.
         write(*,'(A)') " # Anderson MD selected, switching FIXLAT"
         if(any(fixlat)) write(*,*) "# FIXLAT: a,b,c,alpha,beta,gamma,shape ", fixlat
  endif
  if(fixlat(7).and.trim(parini%paropt_geopt%approach).ne."FIRE") stop "Fixed cell shape only implemented in FIRE"

!Put all the atoms back into the cell
  !call backtocell(nat,pos_latvec,pos_red)

!Allocate some other arrays
!  allocate(poslocmin(3,nat,npminx))
!  allocate(latlocmin(3,3,npminx))
  allocate(poshop(3,nat))

!Read input parameters from ioput
  write(filename,'(a6,i3.3)') 'ioput'   
  open(unit=11,file='ioput',status='old')
  read(11,*) ediff,ekinetic,ekinetic_max
  close(11)

  write(*,'(a,1x,3(1x,e10.3),1x,i4)') ' # In :ediff,temperature,maximum temperature,nsoften',&
        &ediff,ekinetic,ekinetic_max,nsoften

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
      call read_atomic_file_ascii(filename,nat,units,pl_arr(1:3,1:nat,kk),lat_arr(1:3,1:3,kk),f_arr(1:3,1:nat,kk),str_arr(1:6,kk),&
           &fixat,fixlat,readfix,fragarr,readfrag,tmp_enthalpy,tmp_energy)
      write(*,*) '# read file ',filename
     elseif(kk==1) then 
      n_arr=3*nat*nlmin
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
        if(findsym.and.spg_arr(kk).lt.1) then
          write(*,'(a,a)') " # Recomputing symmetry for ",trim(filename)
!          call find_symmetry(parini,nat,poslocmin(1:3,1:nat,npmin),latlocmin(1:3,1:3,npmin),typat,tolmin,tolmax,ntol,&
!          &elocmin(npmin,4),spgint)
          call find_symmetry(parini,nat,pl_arr(1:3,1:nat,kk),lat_arr(1:3,1:3,kk),typat,tolmin,tolmax,ntol,&
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
           call get_fp(fp_len,pl_arr(:,:,kk),lat_arr(:,:,kk),fp_arr(:,kk))
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
       allocate(e_arr_t(nlminx),ent_arr_t(nlminx),fp_arr_t(fp_len,nlminx),pl_arr_t(3,nat,nlminx),f_arr_t(3,nat,nlminx),&
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
             call insert(nlminx,nlmin_t,fp_len,nat,k_e_pos,e_arr(kk),ent_arr(kk),fp_arr(:,kk),pl_arr(:,:,kk),&
                         &lat_arr(:,:,kk),f_arr(:,:,kk),str_arr(:,kk),spg_arr(kk),spgtol_arr(kk),dos_arr(kk),&
                         &e_arr_t,ent_arr_t,fp_arr_t,pl_arr_t,lat_arr_t,f_arr_t,str_arr_t,spg_arr_t,spgtol_arr_t,dos_arr_t,ct_arr_t)
             k_e_pos=k_e_pos+1
          endif
          if(kk.eq.nlmin) exit
       enddo 
     call winter(parini,nat,units,ent_pos,e_pos,pos_red,pos_latvec,pos_fcart,pos_strten,nlminx,nlmin_t,npminx,& 
          &ent_arr_t,e_arr_t,ct_arr_t,spg_arr_t,spgtol_arr_t,dos_arr_t,pl_arr_t,lat_arr_t,f_arr_t,str_arr_t,fp_arr_t,fp_len,ent_delta,fp_delta,& 
          &eref,ediff,ekinetic,ekinetic_max,dt,nsoften,char_type,ntypat,typat,fixat,fixlat,target_pressure_habohr)
     write(*,'(a,i5,a)') " # Rebinning completed, new nlmin: ",nlmin_t,". rebin.in will be deteled!"
     call system("rm rebin.in")
     goto 3001 
  endif 
     
    
!Just print a fingerprint distance grid
  file_exists=.false.
  INQUIRE(FILE="plot_fpgrid.in", EXIST=file_exists)
     if(file_exists) then
       write(*,'(a)') " # plot_fpgrid.in found. Will write Howtoplot_fingerprint"
       call plot_fp_grid(parini,nlminx,nlmin,nat,fp_len,fp_arr,lat_arr,pl_arr)
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
write(*,'(a,i5)') " # Number of poslocm_ files found: ",nhop
     write(fn5,'(i5.5)') nhop
     if(parini%verb.ge.2) folder="data_hop_init/"
     if(parini%verb.ge.2) call system("mkdir "//trim(folder))


!A single point calculation to get the initial energy (will be converted to enthalpy)
  iprec=1;getwfk=.false.
  call get_energyandforces_single(parini,pos_latvec,pos_red,pos_fcart,pos_strten,e_pos,iprec,getwfk)
!Convert to enthalpy
  call get_enthalpy(pos_latvec,e_pos,target_pressure_habohr,ent_pos)
!Write to file and exit if geopt iteration is 0
  if(parini%verb.gt.0) then
       filename=trim(folder)//"posinit.ascii"
       call write_atomic_file_ascii(parini,filename,nat,units,pos_red,pos_latvec,pos_fcart,pos_strten,&
            &char_type(1:ntypat),ntypat,typat,fixat,fixlat,e_pos,target_pressure_habohr,ent_pos,e_pos)
       if(parini%verb.ge.3) then
       filename=trim(folder)//"posinit.vasp"
       call write_atomic_file_poscar(parini,filename,nat,units,pos_red,pos_latvec,pos_fcart,pos_strten,&
            &char_type(1:ntypat),ntypat,typat,fixat,fixlat,e_pos,target_pressure_habohr,ent_pos,e_pos)
       endif
  endif
  if (parini%paropt_geopt%nit.le.0) goto 3000
!Check if the structure is already relaxed
  call convcheck(nat,pos_latvec,pos_fcart,pos_strten,target_pressure_habohr,parini%paropt_geopt%strfact,fmax,fmax_at,fmax_lat,parini%paropt_geopt%fmaxtol,iexit)
if(iexit==1) then
  write(*,'(a,es15.7,es15.7)') ' # Input structure already relaxed. Proceeding with: Energy, Enthalpy: ',e_pos,ent_pos
  else 
  write(*,'(a,es15.7,es15.7)') ' # Calling the geometry optimizer for the first time here. Energy, Enthalpy: ',e_pos,ent_pos

!Call geometry optimizer for the first time
  vel_in=0.d0
  vel_lat_in=0.d0
  vel_vol_in=0.d0
   iprec=1
   if(parini%geopt_ext) then
     call geopt_external(parini,pos_latvec,pos_red,pos_fcart,pos_strten,e_pos,iprec,ka,kb,kc,counter)
   else
     if(confine==2) confine=1
     if(parini%paropt_geopt%approach=="RBFGS")  call GEOPT_RBFGS_MHM(parini,pos_latvec,pos_red,pos_fcart,pos_strten,e_pos,iprec,counter,folder)
     if(parini%paropt_geopt%approach=="MBFGS")  call GEOPT_MBFGS_MHM(parini,pos_latvec,pos_red,pos_fcart,pos_strten,e_pos,iprec,counter,folder)
     if(parini%paropt_geopt%approach=="FIRE")   call GEOPT_FIRE_MHM(parini,pos_latvec,pos_red,pos_fcart,pos_strten,vel_in,vel_lat_in,&
                                                    &vel_vol_in,e_pos,iprec,counter,folder)
     if(parini%paropt_geopt%approach=="SQNM")   call      GEOPT_SQNM(parini,pos_latvec,pos_red,pos_fcart,pos_strten,e_pos,iprec,counter,folder)
     if(parini%paropt_geopt%approach=="QBFGS")  call     GEOPT_qbfgs(parini,pos_latvec,pos_red,pos_fcart,pos_strten,e_pos,iprec,counter,folder)
     if(parini%paropt_geopt%approach=="SD")     call        GEOPT_SD(parini,pos_latvec,pos_red,pos_fcart,pos_strten,e_pos,iprec,counter,folder)
   endif

!Update GEOPT counter
     count_geopt=count_geopt+counter     
     write(*,'(a,i7)') " # Counter of GEOPT updated: ", int(count_geopt)
endif

!Check for symmetry
   if(findsym) then
      call find_symmetry(parini,nat,pos_red,pos_latvec,typat,tolmin,tolmax,ntol,spgtol_pos,spgint)
      spg_pos=spgint
   else
      spg_pos=0
      spgtol_pos=0.d0
   endif

!Accept first structure
   accepted=accepted+1.d0

!Get FDOS
  if(finddos) then
     iprec=1
     if(usewf_geopt) then
        getwfk=.true.
     else
        getwfk=.false.
     endif
     call get_dos(parini,pos_latvec,pos_red,efermi,fdos_pos,iprec,getwfk)
     write(*,'(a,es15.7,es15.7)') " # FDOS found: ",fdos_pos,efermi
  endif

!Re-read parameters
!   write(*,'(a)')  " # Re-reading params.in"
!   call read_params()
  !call params_read(parini)

!Get the fingerprint
   call get_fp(fp_len,pos_red,pos_latvec,fp_pos)

!Convert to enthalpy
  call get_enthalpy(pos_latvec,e_pos,target_pressure_habohr,ent_pos)

!  nconjgr=0
  if(parini%verb.gt.0) then
     filename='poslocm_'//fn5//'.ascii'
     call write_atomic_file_ascii(parini,filename,nat,units,pos_red,pos_latvec,pos_fcart,pos_strten,char_type,ntypat,typat,fixat,fixlat,&
     &e_pos,target_pressure_habohr,ent_pos-eref,e_pos-eref)
     call dist_latvec2ang(dist_ang,pos_latvec,pi)
     write(*,'(a,a,a,6(1x,es15.7))') " # Wrote file ",trim(filename),", with a,b,c,alpha,beta,gamma: ",dist_ang
  endif

!Do rounding of the received enthalpy
!!  re_pos=round(e_pos-eref,accur)
!!  rent_pos=round(ent_pos-eref,accur)
!!     write(67,'(a,1x,5(1x,1pe17.10))') 'INPUT(relaxed): ent_pos,rent_pos,e_pos,re_pos,eref ',ent_pos,rent_pos,e_pos,re_pos,eref
!!     write(*,'(a,1x,5(1x,1pe17.10))') ' # INPUT(relaxed): ent_pos,rent_pos,e_pos,re_pos,eref ',ent_pos,rent_pos,e_pos,re_pos,eref
!     write(67,'(a,1x,2(1x,1pe17.10))') 'INPUT(relaxed): ent_pos,e_pos ',ent_pos,e_pos
     write(*,'(a,1x,2(1x,1pe17.10))') ' # INPUT(relaxed): ent_pos,e_pos ',ent_pos,e_pos
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
     do iat=1,nat
!!        poslocmin(1,iat,1)=pos_red(1,iat)
!!        poslocmin(2,iat,1)=pos_red(2,iat)
!!        poslocmin(3,iat,1)=pos_red(3,iat)
        pl_arr(1,iat,1)=pos_red(1,iat)
        pl_arr(2,iat,1)=pos_red(2,iat)
        pl_arr(3,iat,1)=pos_red(3,iat)
     enddo
!!     latlocmin(:,:,1)=pos_latvec(:,:)
     lat_arr(:,:,1)=pos_latvec(:,:)
     write(*,*) '#first configuration saved'
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
         call insert(nlminx,nlmin,fp_len,nat,k_e_pos,e_pos,ent_pos,fp_pos,pos_red,pos_latvec,pos_fcart,pos_strten,&
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
     write(*,*) '# initial re_sm',re_sm
     write(222,'(i10,(1x,f10.0),1x,1pe21.14,2(1x,1pe10.3),i5,1x,1pe10.3)')&
                &nhop,escape,ent_pos-eref,ediff,ekinetic,spg_pos,fdos_pos
     write(*,'(a,(1x,f10.0),1x,1pe21.14,2(1x,1pe10.3),i5,1x,1pe10.3)')" #global: ",&
     &escape,ent_pos-eref,ediff,ekinetic,spg_pos,fdos_pos
     close(222)

  nlmin_old=nlmin
!  CPUcheck=.false.

!Outer (hopping) loop
!Hopping_loop: do
1000 continue
     if (nlmin >= nlminx) then
        goto 3000
     endif
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
     do iat=1,nat
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
     call init_vel(parini,vel_in,vel_lat_in,vel_vol_in,wpos_latvec,wpos_red,parini%bmass*amu_emass,ekinetic,nsoften,folder)

!Call molecular dynamics 
     escape=escape+1.d0
     iprec=2
     if((all(fixlat(1:6)).and..not.fixlat(7)).or.bc==2) then
     call MD_fixlat(parini,parres,wpos_latvec,wpos_red,wpos_fcart,wpos_strten,vel_in,e_wpos,iprec,counter,folder)
     else
     call MD_MHM(parini,parres,wpos_latvec,wpos_red,wpos_fcart,wpos_strten,vel_in,vel_lat_in,vel_vol_in,e_wpos,iprec,counter,folder)
!       call MD_ANDERSEN_MHM(parini,wpos_latvec,wpos_red,wpos_fcart,wpos_strten,vel_in,vel_lat_in,0.003d0,e_wpos,iprec,counter)
     endif

!Make the cell "good" before entering GEOPT
!Put atoms back into cell after MD
if(.not.(any(fixlat).or.any(fixat).or.confine.ge.1)) call correct_latvec(wpos_latvec,wpos_red,nat,correctalg,latvec_io)

!Convert to enthalpy
     call get_enthalpy(wpos_latvec,e_wpos,target_pressure_habohr,ent_wpos)

!Show the energy and enthalpies after exiting MD
     write(*,'(a,2(1x,es25.15))') " # After exiting MD, we have energy, enthalpy:    ",e_wpos,ent_wpos
!     write(67,'(a,2(1x,es25.15))') "After exiting MD, we have energy, enthalpy:    ",e_wpos,ent_wpos

!Update MD counter
     count_md=count_md+counter     
     write(*,'(a,i7)') " # Counter of MD updated: ", int(count_md)

!Keep track of kinetic energy as average, not yet decided it this should not rather be the temperature. Not yet implemented
     av_ekinetic=av_ekinetic+ekinetic

!Call geometry optimizer after MD
  vel_in=0.d0
  vel_lat_in=0.d0
   iprec=1
   if(parini%geopt_ext) then
     call geopt_external(parini,wpos_latvec,wpos_red,wpos_fcart,wpos_strten,e_wpos,iprec,ka,kb,kc,counter)
   else
      if(confine==2) confine=1
      if(parini%paropt_geopt%approach=="RBFGS") call GEOPT_RBFGS_MHM(parini,wpos_latvec,wpos_red,wpos_fcart,&
                                     &wpos_strten,e_wpos,iprec,counter,folder)
      if(parini%paropt_geopt%approach=="MBFGS") call GEOPT_MBFGS_MHM(parini,wpos_latvec,wpos_red,wpos_fcart,&
                                     &wpos_strten,e_wpos,iprec,counter,folder)
      if(parini%paropt_geopt%approach=="FIRE")  call GEOPT_FIRE_MHM(parini,wpos_latvec,wpos_red,wpos_fcart,&
                                     &wpos_strten,vel_in,vel_lat_in,vel_vol_in,e_wpos,iprec,counter,folder)
      if(parini%paropt_geopt%approach=="SQNM")   call    GEOPT_SQNM(parini,wpos_latvec,wpos_red,wpos_fcart,&
                                     &wpos_strten,e_wpos,iprec,counter,folder)
      if(parini%paropt_geopt%approach=="QBFGS")  call   GEOPT_QBFGS(parini,wpos_latvec,wpos_red,wpos_fcart,&
                                     &wpos_strten,e_wpos,iprec,counter,folder)
      if(parini%paropt_geopt%approach=="SD")     call      GEOPT_SD(parini,wpos_latvec,wpos_red,wpos_fcart,&
                                     &wpos_strten,e_wpos,iprec,counter,folder)
   endif

!Update GEOPT counter
     count_geopt=count_geopt+counter     
     write(*,'(a,i7)') " # Counter of GEOPT updated: ", int(count_geopt)

!Check for symmetry
   if(findsym) then
      call find_symmetry(parini,nat,wpos_red,wpos_latvec,typat,tolmin,tolmax,ntol,spgtol_wpos,spgint)
      spg_wpos=spgint
   else
      spg_wpos=0
      spgtol_wpos=0.d0
   endif

!Get FDOS
  if(finddos) then
     iprec=1;getwfk=.true.
     call get_dos(parini,wpos_latvec,wpos_red,efermi,fdos_wpos,iprec,getwfk)
     write(*,'(a,es15.7,es15.7)') " # FDOS found: ",fdos_wpos,efermi
  endif


!Re-read parameters
!   write(*,'(a)')  " # Re-reading params.in"
!   call read_params()
  !call params_read(parini)

!Convert to enthalpy
  call get_enthalpy(wpos_latvec,e_wpos,target_pressure_habohr,ent_wpos)

!Put atoms back into cell after MD+GEOPT
  !call backtocell(nat,wpos_latvec,wpos_red)

!Show the energy and enthalpies after exiting GEOPT 
  write(*,'(a,2(1x,es25.15))') " # After exiting GEOPT, we have energy, enthalpy: ",e_wpos,ent_wpos
!     write(67,'(a,2(1x,es25.15))') "After exiting GEOPT, we have energy, enthalpy: ",e_wpos,ent_wpos


!Write out the relaxed configuration
  if(parini%verb.gt.0) then
        filename='poslocm_'//fn5//'.ascii'
        call write_atomic_file_ascii(parini,filename,nat,units,wpos_red,wpos_latvec,wpos_fcart,wpos_strten,char_type,&
        &ntypat,typat,fixat,fixlat,e_wpos,target_pressure_habohr,ent_wpos-eref,e_wpos-eref)
  call dist_latvec2ang(dist_ang,wpos_latvec,pi)
  write(*,'(a,a,a,6(1x,es15.7))') " # Wrote file ",trim(filename),", with a,b,c,alpha,beta,gamma: ",dist_ang
  endif

!Get the fingerprint
 call get_fp(fp_len,wpos_red,wpos_latvec,fp_wpos)

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
  write(*,'(a,i3,2(1x,1pe14.7))')  &
       ' # nlmin,ent_wpos,ent_pos', nlmin,ent_wpos,ent_pos
  write(*,'(a,i3,2(1x,1pe14.7))')  &
       ' # nlmin,e_wpos,e_pos    ', nlmin,e_wpos,e_pos

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
     call wtioput(ediff,ekinetic,ekinetic_max,nsoften)
     open(unit=222,file='global.mon',status='unknown',position='append')
     write(222,'(i10,(1x,f10.0),1x,1pe21.14,2(1x,1pe10.3),i5,1x,1pe10.3,3(1x,0pf5.2),a)')  &
          nhop,escape,ent_wpos-eref,ediff,ekinetic,spg_wpos,fdos_wpos, &
          escape_sam/escape,escape_old/escape,escape_new/escape,'   S'
     close(222)
     write(*,'(a,(1x,f10.0),1x,1pe21.14,2(1x,1pe10.3),i5,1x,1pe10.3,3(1x,0pf5.2),a)') " #global: ", &
          escape,ent_wpos-eref,ediff,ekinetic,spg_wpos,fdos_wpos, &
          escape_sam/escape,escape_old/escape,escape_new/escape,'   S'
     write(*,'(a)')' # no escape from current minimum.'

     if(parini%auto_mdmin) parres%mdmin   = min(parini%mdmin_max, parres%mdmin + 1)  ! MALM
     write(*,*) "# nsoften, mdmin: ", nsoften, parres%mdmin ! MALM

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
  write(*,*) "k_e_wpos,kid: ",k_e_wpos,kid
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
     call wtioput(ediff,ekinetic,ekinetic_max,nsoften)
     nlmin=nlmin+1
!!       call insert(nlminx,nlmin,k_e_wpos,rent_wpos,re_wpos,earr(0,1))
     call insert(nlminx,nlmin,fp_len,nat,k_e_wpos,e_wpos,ent_wpos,fp_wpos,wpos_red,wpos_latvec,wpos_fcart,wpos_strten,&
          &spg_wpos,spgtol_wpos,fdos_wpos,e_arr,ent_arr,fp_arr,pl_arr,lat_arr,f_arr,str_arr,spg_arr,spgtol_arr,dos_arr,ct_arr)
!!     call winter(parini,nat,units,re_pos,rent_pos,e_pos,pos_red,pos_latvec,npminx,nlminx,nlmin,npmin,accur, & 
!!          &earr,elocmin,poslocmin,latlocmin,eref,ediff,ekinetic,ekinetic_max,dt,nsoften,char_type,ntypat,typat,&
!!          &fixat,fixlat,target_pressure_habohr)
     k_e_wpos=k_e_wpos+1
     if (k_e_wpos .gt. nlminx .or. k_e_wpos .lt. 1) stop "k_e_wpos out of bounds"
     if(modulo(nlmin,nwrite_inter)==0) then 
     call winter(parini,nat,units,ent_pos,e_pos,pos_red,pos_latvec,pos_fcart,pos_strten,nlminx,nlmin,npminx,& 
          &ent_arr,e_arr,ct_arr,spg_arr,spgtol_arr,dos_arr,pl_arr,lat_arr,f_arr,str_arr,fp_arr,fp_len,ent_delta,fp_delta,& 
          &eref,ediff,ekinetic,ekinetic_max,dt,nsoften,char_type,ntypat,typat,fixat,fixlat,target_pressure_habohr)
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
    do iat=1,nat
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
   write(*,*) "# nsoften, mdmin: ", nsoften, parres%mdmin ! MALM
   accepted=accepted+1.d0
   e_pos=e_hop
   ent_pos=ent_hop
!   re_pos=re_hop
!   rent_pos=rent_hop
   spg_pos=spg_hop
   spgtol_pos=spgtol_hop
   fdos_pos=fdos_hop
   fp_pos=fp_hop
     do iat=1,nat
        pos_red(1,iat)=poshop(1,iat)
        pos_red(2,iat)=poshop(2,iat)
        pos_red(3,iat)=poshop(3,iat)
     enddo
     pos_latvec=lathop
     pos_fcart=poshop_fcart
     pos_strten=poshop_strten
     call wtioput(ediff,ekinetic,ekinetic_max,nsoften)
      if (abs(ent_wpos-ent_hop).lt.1.d-10) then
         open(unit=222,file='global.mon',status='unknown',position='append')
         write(222,'(i10,(1x,f10.0),1x,1pe21.14,2(1x,1pe10.3),i5,1x,1pe10.3,3(1x,0pf5.2),l3,a,i5)')  &
          nhop,escape,ent_hop-eref,ediff,ekinetic,spg_hop,fdos_hop, &
          escape_sam/escape,escape_old/escape,escape_new/escape,newmin,' A ', nvisit
         close(222)
         write(*,'(a,(1x,f10.0),1x,1pe21.14,2(1x,1pe10.3),i5,1x,1pe10.3,3(1x,0pf5.2),l3,a,i5)') " #global: ", &
          escape,ent_hop-eref,ediff,ekinetic,spg_hop,fdos_hop, &
          escape_sam/escape,escape_old/escape,escape_new/escape,newmin,' A ', nvisit
      else
         open(unit=222,file='global.mon',status='unknown',position='append')
         write(222,'(i10,1x,f10.0,1x,1pe21.14,2(1x,1pe10.3),i5,1x,1pe10.3,3(1x,0pf5.2),a,i5)')  &
              nhop,escape,ent_wpos-eref,ediff,ekinetic,spg_wpos,fdos_wpos,&
              escape_sam/escape,escape_old/escape,escape_new/escape,'   I  ',int(ct_arr(k_e_wpos))
         write(222,'(i10,1x,f10.0,1x,1pe21.14,2(1x,1pe10.3),i5,1x,1pe10.3,3(1x,0pf5.2),a,i5)')  &
              nhop,escape,ent_hop-eref,ediff,ekinetic,spg_hop,fdos_hop, &
              escape_sam/escape,escape_old/escape,escape_new/escape,'   A  ',nvisit
         close(222)
         write(*,'(a,1x,f10.0,1x,1pe21.14,2(1x,1pe10.3),i5,1x,1pe10.3,3(1x,0pf5.2),a,i5)') " #global: ", &
              escape,ent_wpos-eref,ediff,ekinetic,spg_wpos,fdos_wpos, &
              escape_sam/escape,escape_old/escape,escape_new/escape,'   I  ',int(ct_arr(k_e_wpos))
         write(*,'(a,1x,f10.0,1x,1pe21.14,2(1x,1pe10.3),i5,1x,1pe10.3,3(1x,0pf5.2),a,i5)') " #global: ", &
              escape,ent_hop-eref,ediff,ekinetic,spg_hop,fdos_hop, &
              escape_sam/escape,escape_old/escape,escape_new/escape,'   A  ',nvisit
      endif
!      endif
      ent_hop=1.d100
      ediff=max(ediff*alpha1,1.d-4)
! write intermediate results
      write(*,*) 'WINTER'
!!      call winter(parini,nat,units,re_pos,rent_pos,e_pos,pos_red,pos_latvec,npminx,nlminx,nlmin,npmin,accur, & 
!!      &earr,elocmin,poslocmin,latlocmin,eref,ediff,ekinetic,ekinetic_max,dt,nsoften,char_type,ntypat,typat,&
!!      &fixat,fixlat,target_pressure_habohr)
     if(modulo(nlmin,nwrite_inter)==0) then 
     call winter(parini,nat,units,ent_pos,e_pos,pos_red,pos_latvec,pos_fcart,pos_strten,nlminx,nlmin,npminx,& 
          &ent_arr,e_arr,ct_arr,spg_arr,spgtol_arr,dos_arr,pl_arr,lat_arr,f_arr,str_arr,fp_arr,fp_len,ent_delta,fp_delta,& 
          &eref,ediff,ekinetic,ekinetic_max,dt,nsoften,char_type,ntypat,typat,fixat,fixlat,target_pressure_habohr)
     endif
     goto 1000
  else
!local minima rejected
!    if(parini%auto_mdmin) parres%mdmin   = min(mdmin_max, parres%mdmin + 1)  ! MALM
     write(*,*) "# nsoften, mdmin: ", nsoften, parres%mdmin ! MALM
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
     call wtioput(ediff,ekinetic,ekinetic_max,nsoften)
     goto 1000
!------------------------------------------------------------
  endif

!end do hopping_loop
3000 continue   !This is the really end. Collect data

!Print some informations and rewrite restart information
     !write(67,*) 'writing final results'
     write(*,*) ' # writing final results'
     !write(67,*) ' found in total ',nlmin,' minima'
     !write(67,*) ' Accepeted ',int(accepted),' minima'
     write(*,*) ' # found in total ',nlmin,' minima'
     write(*,*) ' # Accepeted ',int(accepted),' minima'
!!     call winter(parini,nat,units,re_pos,rent_pos,e_pos,pos_red,pos_latvec,npminx,nlminx,nlmin,npmin,accur, & 
!!      earr,elocmin,poslocmin,latlocmin,eref,ediff,ekinetic,ekinetic_max,dt,nsoften,char_type,ntypat,&
!!      typat,fixat,fixlat,target_pressure_habohr)
     call winter(parini,nat,units,ent_pos,e_pos,pos_red,pos_latvec,pos_fcart,pos_strten,nlminx,nlmin,npminx,& 
          &ent_arr,e_arr,ct_arr,spg_arr,spgtol_arr,dos_arr,pl_arr,lat_arr,f_arr,str_arr,fp_arr,fp_len,ent_delta,fp_delta,& 
          &eref,ediff,ekinetic,ekinetic_max,dt,nsoften,char_type,ntypat,typat,fixat,fixlat,target_pressure_habohr)


!Print ratios from all the global counters
     if(escape.gt.0.d0) write(*,'(a,3(1x,1pe10.3))') ' # ratio stuck,same,old,new', &
          escape_sam/escape,escape_old/escape,escape_new/escape
     write(*,'(a,2(1x,1pe10.3))') ' # ratio acc,rej',accepted/(accepted+rejected),rejected/(accepted+rejected)
     write(*,'(a,2(1x,f12.1))')   ' # count_md,count_geopt',count_md,count_geopt
     if(escape.gt.0.d0) write(*,'(a,2(1x,1pe10.3))') &
          ' # average ediff, ekinetic',av_ediff/(accepted+rejected),av_ekinetic/escape
     write(*,'(a,1x,i8)') ' # number of configurations for which atoms escaped ',nputback

     tt=0.d0
     ss=0.d0
     do i=1,nlmin
        tt=max(tt,real(ct_arr(i),8))
        ss=ss+real(ct_arr(i),8)
     enddo
     write(*,'(a,f8.0)') ' # most frequent visits ',tt
     write(*,'(a,1pe10.3)') ' # av. numb. visits per minimum',ss/nlmin
     write(*,'(a,e9.2)') ' # minimum energy separation between presumably different configurations',egap
     if (escape_sam.gt.0) then
        esep=sqrt(esep/escape_sam)
        write(*,'(a,e9.2)') ' # average energy separation between presumably identical configurations',esep
     endif

!  close(67)


! call timein(tsec(1),tsec(2))
! tsec(1)=tsec(1)-cpui
! tsec(2)=tsec(2)-walli

  deallocate(char_type)
  deallocate(e_arr,ent_arr,ct_arr,pl_arr,lat_arr,dos_arr,spgtol_arr,spg_arr,fp_arr)
3001 continue
!Close socket on slave side
if(trim(parini%potential_potential)=="msock") call socket_stop()
end subroutine task_minhocao
!contains


!**********************************************************************************************
subroutine MD_MHM   (parini,parres,latvec_in,xred_in,fcart_in,strten_in,vel_in,vel_lat_in,vvol_in,etot_in,iprec,counter,folder)
 use mod_interface
 use global, only: target_pressure_habohr,target_pressure_gpa,nat,ntypat,znucl,amu,amutmp,typat
 use global, only: char_type,units,usewf_md
 use global, only: fixat,fixlat
 use defs_basis
 use interface_code
 use modsocket, only: sock_extra_string
 use mod_parini, only: typ_parini
implicit none
 type(typ_parini), intent(in):: parini
 type(typ_parini), intent(inout):: parres
 real(8) :: latvec_in(3,3),xred_in(3,nat),vel_in(3,nat),fcart_in(3,nat),etot_in,strten_in(6),vel_lat_in(3,3),vvol_in
!*********************************************************************
!Variables for my MD part
! real(8), parameter :: Ha_eV=27.21138386_dp ! 1 Hartree, in eV
! real(8), parameter :: kb_HaK=8.617343d-5/Ha_eV ! Boltzmann constant in Ha/K
! real(8), parameter :: amu_emass=1.660538782d-27/9.10938215d-31 ! 1 atomic mass unit, in electronic mass
 real(8):: amass(nat)
 real(8):: pressure_ener
 real(8):: pressure_md(3,3)
 real(8):: pressure
 real(8):: vol
 real(8):: vol0
 real(8):: volprev
 real(8),dimension(3,nat):: xcart
 real(8),dimension(3,nat):: fposcur
 real(8),dimension(3,nat):: accposcur
 real(8),dimension(3,nat):: accpospred
 real(8),dimension(3,nat):: accposprev
 real(8),dimension(3,nat):: fpospred
 real(8),dimension(3,nat):: vpospred
 real(8),dimension(3,nat):: poscur
 real(8),dimension(3,nat):: vxyz
 real(8),dimension(3,nat):: vposcur
 real(8),dimension(3,nat):: pospred
 real(8),dimension(3,nat):: dxred
 real(8),dimension(3,3):: dlatvec
 real(8),dimension(3,3):: latvec
 real(8),dimension(3,3):: latvec0
 real(8),dimension(3,3):: latinv
 real(8),dimension(3,3):: unitmat
 real(8),dimension(3,3):: elmatcur
 real(8),dimension(3,3):: acclatcur
 real(8),dimension(3,3):: acclatpred
 real(8),dimension(3,3):: acclatprev
 real(8),dimension(3,3):: flatcur
 real(8),dimension(3,3):: flatpred
 real(8),dimension(3,3):: velmatcur
 real(8),dimension(3,3):: tmplat
 real(8),dimension(3,3):: tmplatt
 real(8),dimension(3,3):: tpred
 real(8),dimension(3,3):: f0
 real(8),dimension(3,3):: f0inv
 real(8),dimension(3,3):: latcur
 real(8),dimension(3,3):: vlatcur
 real(8),dimension(3,3):: vlat
 real(8),dimension(3,3):: lattrans
 real(8),dimension(3,3):: latdottrans
 real(8),dimension(3,3):: a 
 real(8),dimension(3,3):: g
 real(8),dimension(3,3):: gdot
 real(8),dimension(3,3):: ginv
 real(8),dimension(3,3):: gtot
 real(8),dimension(3,3):: sigma
 real(8),dimension(3,3):: sigmatrans
 real(8),dimension(3,3):: latpred
 real(8),dimension(3,3):: tcur
 real(8),dimension(3,3):: vlatpred
 real(8),dimension(3,3):: velmatpred
 real(8),dimension(3,3):: str_matrix
 real(8),dimension(1000):: ensave   !zl
 real(8),dimension(1000):: ensmoth  !zl
 real(8):: accvolcur
 real(8):: accvolpred
 real(8):: accvolprev
 real(8):: fvolcur
 real(8):: fvolpred
 real(8):: volcur
 real(8):: volpred
 real(8):: vvolcur
 real(8):: vvolpred
 real(8):: volpred_1_3
 real(8):: vol_1_3_inv
 real(8):: vol_1_3

 real(8):: vposcurtmp(3)
 real(8):: vpospredtmp(3)
 real(8):: crossp(3)
 real(8):: velcm(3)
 real(8):: latmass
 real(8):: latmassinv
 real(8):: latmass0
 real(8):: ekinatom
 real(8):: ekinatom_prev
 real(8):: ekinlat
 real(8):: ekinlat_prev
 real(8):: rkin
 real(8):: enthalpy
 real(8):: enmin1
 real(8):: enmin2
 real(8):: ent_pos_0
 real(8):: en0000
 real(8):: e_rxyz
 real(8):: econs_max
 real(8):: econs_min
 real(8):: torquenrm 
 real(8):: ecut_tmp
 real(8):: toldff_tmp
 real(8):: counter
 real(8):: dt
 real(8):: vpressure
 real(8):: dt_ratio
 integer:: i,ii
 integer:: j
 integer:: iat
 integer:: itime
 integer:: nummax
 integer:: nummin
 integer:: iprec 
 integer:: options 
 integer:: md_type
 logical:: getwfk

 character(40)::filename,folder
 character(4) ::fn4
 write(*,*) "# PARAMETERS: Ha_eV,kb_HaK,amu_emass -->",Ha_eV,kb_HaK,amu_emass


!Assign masses to each atom (for MD)
 do iat=1,nat
   amass(iat)=amu_emass*amu(typat(iat))
if(parres%verb.gt.0) write(*,'(a,i5,2(1x,es15.7))') " # MD: iat, AMU, EM: ", iat, amu(typat(iat)),amass(iat)
 enddo
!******************************************************************
!NEW VERISON OF MD: Directly implemented, simplest Parrinello-Rahman and other MD

!Here we split the routine in the Abinit native part and my new implentation
md_type=parres%md_algo
if(md_type==1) then
 write(*,'(a)') ' # Entering standalone Parrinello Rahman MD '
elseif(md_type==2) then
 write(*,'(a)') ' # Entering standalone Cleveland MD '
elseif(md_type==3) then
 write(*,'(a)') ' # Entering standalone Wentzcovitch MD '
elseif(md_type==4) then
 write(*,'(a)') ' # Entering standalone Andersen MD '
endif

!The "reduced" coordinates in Andersen are quite different from the ones in PR
!Set temporary variables, initially
  vxyz(:,:)=vel_in(:,:)
  pressure=target_pressure_habohr
  unitmat=0.d0
  do i=1,3
    unitmat(i,i)=1.d0
  enddo
  latmass0=parini%bmass*amu_emass !This means we use the barostat mass as the lattice mass (in ELECTRON MASS)
  vlat=vel_lat_in  !The initial cell velocity
  itime=0
  dt=parres%dtion_md

!Set options=1 for Velocity Verlet of cell dynamics
!Set options=2 for Normal Verlet of cell dynamics
!Set options=3 for Beeman integration scheme, corrector-predictor
  options=parres%md_integrator
  if(options.lt.1.or.options.gt.3) stop "Wrong algo option"

!MD-type: 1 for PR and 2 for Cleveland and 3 for Wentzcovitch and 4 for Andersen
  md_type=parres%md_algo
  if(md_type.lt.1.or.md_type.gt.4) stop "Wrong integrator option"

if(parres%verb.gt.0) write(*,'(a,i3,a,i3)') " # MD Algorithm: ",md_type, ", MD Integrator: ",options

!Now we run my implementation of MD
pressure_md=0.d0
pressure_md(1,1)=1.d0;pressure_md(2,2)=1.d0;pressure_md(3,3)=1.d0
pressure_ener=0.d0;pressure_md=pressure_md*pressure  !Here the pressure is not passed to the energyandforces, so we move on the ENERGY surface

!Transform the initial velocities given in cartesian coordinates into velocities of the intenal coordinates
       latvec=latvec_in
       call getvol(latvec,vol)
if(md_type==4) then
       vol_1_3=vol**(1.d0/3.d0)
       vol_1_3_inv=1.d0/vol_1_3
       do iat=1,nat
          vposcur(:,iat)=vxyz(:,iat)*vol_1_3_inv
       enddo
else
       call invertmat(latvec,latinv,3)
       do iat=1,nat
          vposcur(:,iat)=matmul(latinv,vxyz(:,iat))
       enddo
endif

!EXPERIMENTAL: add the pressure part from the initial velocities to the MD
if(md_type==1.and.parres%md_presscomp.gt.0.d0) then
  call stress_velocity(vposcur,latvec,amass,nat,vpressure)
  write(*,'(a,f10.5)') " # Internal pressure on top of external pressure: ", vpressure*HaBohr3_GPA*parres%md_presscomp
  pressure_md(1,1)=pressure_md(1,1)+vpressure*parres%md_presscomp
  pressure_md(2,2)=pressure_md(2,2)+vpressure*parres%md_presscomp
  pressure_md(3,3)=pressure_md(3,3)+vpressure*parres%md_presscomp
endif

!Compute f0
        call cross_product(latvec(:,2),latvec(:,3),crossp); sigma(:,1)=crossp
        call cross_product(latvec(:,3),latvec(:,1),crossp); sigma(:,2)=crossp
        call cross_product(latvec(:,1),latvec(:,2),crossp); sigma(:,3)=crossp
        sigmatrans(:,1)=sigma(1,:)
        sigmatrans(:,2)=sigma(2,:)
        sigmatrans(:,3)=sigma(3,:)
        f0=matmul(sigmatrans,sigma)
        call invertmat(f0,f0inv,3)


!!ATTENTION: THIS IS WITH TAKING INTO ACCOUNT THE STRESS TENSOR GIVEN IN INPUT

!Keep track of volume
        vol0=vol
        latvec0=latvec/vol**(1.d0/3.d0)
        volprev=vol

        if(md_type==1) then
!!This is for the original Rahman Parrinello. You also need to change the volume-mass factor
        latmass=latmass0
        latmassinv=1.d0/latmass
        elseif(md_type==2) then
!This is for cleveland
        latmass=latmass0/vol**(4.d0/3.d0)
        latmassinv=1.d0/latmass
        elseif(md_type==3) then
!This is for Wentzcovitch
        latmass=latmass0
        latmassinv=1.d0/latmass
        elseif(md_type==4) then
!This is for Andersen
        latmass=latmass0
        latmassinv=1.d0/latmass
        endif

!Initialize internal coordinates and lattice vectors
if(md_type==4) then
        call rxyz_int2cart(latvec,xred_in,xcart,nat)
        poscur=xcart*vol_1_3_inv           !call rxyz_cart2int(latvec,poscur,rxyz,nat)
        volcur=vol
        vvolcur=vvol_in
else
        poscur=xred_in           !call rxyz_cart2int(latvec,poscur,rxyz,nat)
        latcur=latvec
        vlatcur=vlat
endif
!Write every step
!        call rxyz_int2cart(latcur,poscur,rxyz,nat)
!        call wtpos_inter(nat,rxyz,latcur,0)

!MHM: initialize variable to track the minima/maxima***********
    write(*,*) '# MINHOP start MD'
    nummax=0
    nummin=0
    enmin1=0.d0
    en0000=0.d0
    econs_max=-1.d100
    econs_min=1.d100
!**************************************************************

!Compute initial kinetic energies
if(md_type==4) then
  latcur=latvec0*vol_1_3
  rkin=0.d0
  do iat=1,nat
     vposcurtmp=vposcur(:,iat)*vol_1_3
     rkin=rkin+amass(iat)*(vposcurtmp(1)**2+vposcurtmp(2)**2+vposcurtmp(3)**2)
  enddo
  ekinatom=0.5d0*rkin
  rkin=vvolcur**2*latmass
  ekinlat=0.5d0*rkin
else
  rkin=0.d0
  do iat=1,nat
     vposcurtmp=matmul(latcur,vposcur(:,iat))
     rkin=rkin+amass(iat)*(vposcurtmp(1)**2+vposcurtmp(2)**2+vposcurtmp(3)**2)
  enddo
  ekinatom=0.5d0*rkin
  rkin=0.d0
  do i=1,3
     rkin=rkin+(vlatcur(1,i)**2+vlatcur(2,i)**2+vlatcur(3,i)**2)*latmass
  enddo
  ekinlat=0.5d0*rkin
endif


!MHM: Write output to file in every step***********************************
!INITIAL STEP, STILL THE SAME STRUCTURE AS INPUT
       write(*,*) "Pressure, Energy",pressure,etot_in
       call get_enthalpy(latvec_in,etot_in,pressure,enthalpy)
       ent_pos_0=enthalpy
       en0000=enthalpy-ent_pos_0
if(parres%verb.gt.0) then
       write(*,'(a,i5,1x,4(1x,1pe17.10),3(1x,i2))') ' # MD: it,enthalpy,pV,ekinat,ekinlat,nmax,nmin,mdmin ',&
             &itime,enthalpy,(pressure_md(1,1)+pressure_md(2,2)+pressure_md(3,3))/3.d0*vol,ekinatom,ekinlat,nummax,nummin,parres%mdmin
!             &itime,enthalpy,pressure*vol,ekinatom,ekinlat,nummax,nummin,parres%mdmin
       write(fn4,'(i4.4)') itime
       filename=trim(folder)//"posmd."//fn4//".ascii"
       units=units
       write(*,*) "# Writing the positions in MD:",filename
       call write_atomic_file_ascii(parres,filename,nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &char_type(1:ntypat),ntypat,typat,fixat,fixlat,etot_in,pressure,enthalpy,en0000)
endif
!*********************************************************************

call acceleration(pressure_md,accposcur,acclatcur,accvolcur,vposcur,vlatcur,vvolcur,&
                  strten_in,fcart_in,latcur,amass,latmass,f0inv,md_type,nat)
        accposprev=accposcur
        acclatprev=acclatcur
        accvolprev=accvolcur
!call acceleration(pressure_md,accposcur,acclatcur,vposcur,vlatcur,strten_in,fcart_in,latcur,amass,latmass,f0inv,md_type,nat) 
!        accposprev=accposcur
!        acclatprev=acclatcur

!Eliminate the unnecessary atmic and cell components
call elim_fixed_at(nat,accposcur);call elim_fixed_at(nat,accposprev)
if(md_type.ne.4) call elim_fixed_lat(latcur,acclatcur)
if(md_type.ne.4) call elim_fixed_lat(latcur,acclatprev)
call elim_fixed_at(nat,vposcur)
if(md_type.ne.4) call elim_fixed_lat(latcur,vlatcur)


        do itime=1,parres%nmd_dynamics

!          if(itime.ne.1) e_rxyz=enthalpy  !e_rxyz=e_rxyz+pressure_md(1,1)*vol   
!Check the torque on the cell for rotation
          call torque_cell(latcur,vlatcur,torquenrm)
!Now perform the step
          velcm=0.d0
          rkin=0.d0
          if(options==1) then
!For velocity verlet
!             latpred(:,:)=latcur(:,:) + dt*vlatcur(:,:) + (.5d0*dt*dt*acclatcur)
             dlatvec(:,:)=dt*vlatcur(:,:) + (.5d0*dt*dt*acclatcur)
          do iat=1,nat
              dxred(:,iat)=dt*vposcur(:,iat) + 0.5d0*dt*dt*accposcur(:,iat)  !0.5 was missing before
!             pospred(:,iat)=poscur(:,iat) + dt*vposcur(:,iat) + dt*dt*accposcur(:,iat)
          enddo
          call propagate(nat,poscur,latcur,dxred,dlatvec,pospred,latpred)
          do iat=1,nat
             vpospred(:,iat)=(pospred(:,iat)-poscur(:,iat))/dt   !Update dummy variable vposcur for next iteration
          enddo
            if(md_type==4) then
               volpred=volcur+dt*vvolcur+0.5d0*dt*dt*accvolcur
            endif
          elseif(options==2) then
!Instead if normal verlet is used:
!              latpred(:,:)=latcur(:,:) + dt*vlatcur + dt*dt*acclatcur
              dlatvec(:,:)=dt*vlatcur + dt*dt*acclatcur
          do iat=1,nat
              dxred(:,iat)=dt*vposcur(:,iat) + dt*dt*accposcur(:,iat)
!             pospred(:,iat)=poscur(:,iat) + dt*vposcur(:,iat) + dt*dt*accposcur(:,iat)
          enddo
          call propagate(nat,poscur,latcur,dxred,dlatvec,pospred,latpred)
          vlatpred(:,:)=(latpred-latcur)/dt
          do iat=1,nat
              vpospred(:,iat)=(pospred(:,iat)-poscur(:,iat))/dt   !Update dummy variable vposcur for next iteration
          enddo
            if(md_type==4) then
               volpred=volcur+dt*vvolcur+dt*dt*accvolcur
               vvolpred=(volpred-volcur)/dt
            endif
          elseif(options==3) then
!Predictor part of Beeman for cell
!              latpred(:,:)=latcur(:,:) + dt*vlatcur + dt*dt/6.d0*(4.d0*acclatcur(:,:)-acclatprev(:,:))
              dlatvec(:,:)=dt*vlatcur + dt*dt/6.d0*(4.d0*acclatcur(:,:)-acclatprev(:,:))
          do iat=1,nat
!Predictor part of Beeman
!             pospred(:,iat)=poscur(:,iat) + dt*vposcur(:,iat) + dt*dt/6.d0*(4.d0*accposcur(:,iat)-accposprev(:,iat))             
             dxred(:,iat)=dt*vposcur(:,iat) + dt*dt/6.d0*(4.d0*accposcur(:,iat)-accposprev(:,iat))             
          enddo
          call propagate(nat,poscur,latcur,dxred,dlatvec,pospred,latpred)
          vlatpred(:,:)=vlatcur(:,:) + 0.5d0*(3.d0*acclatcur(:,:)-acclatprev(:,:))  
          do iat=1,nat
             vpospred(:,iat)=vposcur(:,iat)+0.5d0*(3.d0*accposcur(:,iat)-accposprev(:,iat))
          enddo
            if(md_type==4) then
               volpred=volcur+dt*vvolcur+dt*dt/6.d0*(4.d0*accvolcur-accvolprev)
               vvolpred=volcur+0.5d0*(3.d0*accvolcur-accvolprev)
            endif
          endif

if(md_type==4) then
  do iat=1,nat
     vposcurtmp=vposcur(:,iat)*vol_1_3
     rkin=rkin+amass(iat)*(vposcurtmp(1)**2+vposcurtmp(2)**2+vposcurtmp(3)**2)
     velcm=velcm+vposcurtmp*amass(iat)
  enddo
else
  do iat=1,nat
             vposcurtmp=matmul(latcur,vposcur(:,iat))
             rkin=rkin+amass(iat)*(vposcurtmp(1)**2+vposcurtmp(2)**2+vposcurtmp(3)**2)
             velcm=velcm+vposcurtmp*amass(iat)
  enddo
endif


!!!!Now calculate steps for the lattice according to the velocity verlet algorithm
!!!          if(options==1) then
!!!          elseif(options==2) then
!!!          elseif(options==3) then
!!!          else
!!!              stop "Wrong option for the Cell part in MD"
!!!          endif



!Kinetic energy of atoms
          ekinatom=0.5d0*rkin
          if(parres%verb.gt.0) write(*,'(a,es15.7)')"Velocity of CM: ", sqrt(velcm(1)**2+velcm(2)**2+velcm(3)**2)
!Kinetic energy according to Parrinello Rahman
          rkin=0.d0
if(md_type==4) then
  rkin=vvolcur**2*latmass
  ekinlat=0.5d0*rkin
else
          do i=1,3
             rkin=rkin+(vlatcur(1,i)**2+vlatcur(2,i)**2+vlatcur(3,i)**2)*latmass
          enddo
          ekinlat=0.5d0*rkin
endif
          rkin=ekinlat+ekinatom
if(parres%verb.gt.0) write(*,'(a,3(1x,es15.7))') " # Torquenrm, Ekin, Enthalpy: ",torquenrm, rkin,enthalpy
!Update counter
          enmin2=enmin1
          enmin1=en0000

!Here we perform the force call
!****************************************************************************************************************        
!****************************************************************************************************************        
!Get the current acell_in and rprim_in, and also the atomic positions
if(md_type==4) then
       volpred_1_3=volpred**(1.d0/3.d0)
       latvec_in=latvec0*volpred_1_3
       latpred=latvec_in
       xcart=pospred*volpred_1_3
       call rxyz_cart2int(latvec_in,xred_in,xcart,nat)
else
       xred_in=pospred 
       latvec_in=latpred
endif
       if(itime==1)  then
           getwfk=.false.
       elseif(itime.ne.1.and.usewf_md) then
           getwfk=.true.
       else
           getwfk=.false.
       endif
       write(fn4,'(i4.4)') itime
       sock_extra_string="MD"//trim(fn4)
       call get_energyandforces_single(parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,getwfk)
!Single point calculation finished. Now returning to MD part. All output variables are now in *_in
!****************************************************************************************************************        

!For velocity verlet of cell
        if(options==1) then
           call acceleration(pressure_md,accpospred,acclatpred,accvolpred,vpospred,&
                &vlatcur,vvolcur,strten_in,fcart_in,latpred,amass,latmass,f0inv,md_type,nat) 
           vlatpred(:,:)=vlatcur(:,:)+0.5d0*dt*(acclatpred+acclatcur)
           if(md_type==4) vvolpred=vvolcur+0.5d0*dt*(accvolpred+accvolcur)
        elseif(options==3) then
!Corrector part of Beeman. Note that a fixed number of iterations are used!!!
           call ekin_at_lat_andersen(amass,latmass,latpred,vpospred,vlatpred,vvolpred,ekinatom_prev,ekinlat_prev,f0,md_type,nat)
           do i=1,5
             call acceleration(pressure_md,accpospred,acclatpred,accvolpred,vpospred,&
                  &vlatpred,vvolpred,strten_in,fcart_in,latpred,amass,latmass,f0inv,md_type,nat) 
             do iat=1,nat
             vpospred(:,iat)=vposcur(:,iat)+dt/6.d0*(2.d0*accpospred(:,iat)+5.d0*accposcur(:,iat)-accposprev(:,iat))
             enddo
             vlatpred=vlatcur+dt/6.d0*(2.d0*acclatpred+5.d0*acclatcur-acclatprev)
             if(md_type==4) vvolpred=vvolcur+dt/6.d0*(2.d0*accvolpred+5.d0*accvolcur-accvolprev)
             call ekin_at_lat_andersen(amass,latmass,latpred,vpospred,vlatpred,vvolpred,ekinatom,ekinlat,f0,md_type,nat)
if(parres%verb.gt.1)write(* ,'(a,i5,1x,es15.5,es15.5)') " # Beeman: Corrector relative convergence: it, atoms, lattice: ",&
                  & i,abs(ekinatom-ekinatom_prev)/ekinatom,abs(ekinlat-ekinlat_prev)/max(ekinlat,1.d-100)
            ! write(67,'(a,i5,1x,es15.5,es15.5)') " # Beeman: Corrector relative convergence: it, atoms, lattice: ",&
            !      & i,abs(ekinatom-ekinatom_prev)/ekinatom,abs(ekinlat-ekinlat_prev)/max(ekinlat,1.d-100)
             ekinlat_prev=ekinlat
             ekinatom_prev=ekinatom
           enddo 
        endif  
call acceleration(pressure_md,accpospred,acclatpred,accvolpred,vpospred,&
     &vlatpred,vvolpred,strten_in,fcart_in,latpred,amass,latmass,f0inv,md_type,nat)

!Eliminate the unnecessary atmic and cell components
call elim_fixed_at(nat,accpospred);call elim_fixed_at(nat,vpospred)
if(md_type.ne.4) call elim_fixed_lat(latcur,acclatpred)
if(md_type.ne.4) call elim_fixed_lat(latcur,vlatpred)

!Compute the "predicted" kinetic energies:
call ekin_at_lat_andersen(amass,latmass,latpred,vpospred,vlatpred,vvolpred,ekinatom,ekinlat,f0,md_type,nat)


!MHM: Write output to file in every step***********************************
       counter=real(itime,8)
       call get_enthalpy(latvec_in,etot_in,pressure,enthalpy)
       call getvol(latvec_in,vol)
       volpred=vol
       en0000=enthalpy-ent_pos_0
!       if (itime >= 3 .and. enmin1 > enmin2 .and. enmin1 > en0000)  nummax=nummax+1
!       if (itime >= 3 .and. enmin1 < enmin2 .and. enmin1 < en0000)  nummin=nummin+1
       ensave(itime) = en0000
       ensmoth(itime) = en0000
       if (itime >= 7) then
           ensmoth(itime-3) = sum(ensmoth(itime-6:itime))/7.d0
       end if

       !if (itime >= 3 .and. enmin1 > enmin2 .and. enmin1 > en0000)  nummax=nummax+1
       !if (itime >= 3 .and. enmin1 < enmin2 .and. enmin1 < en0000)  nummin=nummin+1
       if (itime >=9) then  
        if (itime >= 9 .and. ensmoth(itime-5)<ensmoth(itime-6) .and. ensmoth(itime-5)<ensmoth(itime-7) &
        & .and. ensmoth(itime-5)<ensmoth(itime-4) .and. ensmoth(itime-5)<ensmoth(itime-3)) nummin=nummin+1
        if (itime >= 9 .and. ensmoth(itime-5)>ensmoth(itime-6) .and. ensmoth(itime-5)>ensmoth(itime-7) &
        & .and. ensmoth(itime-5)>ensmoth(itime-4) .and. ensmoth(itime-5)>ensmoth(itime-3)) nummax=nummax+1
       endif
       !write(67,'(a,i5,1x,4(1x,1pe17.10),3(1x,i2))') ' MD: it,enthalpy,pV,ekinat,ekinlat,nmax,nmin,mdmin ',&
       !      &itime,enthalpy,pressure*vol,ekinatom,ekinlat,nummax,nummin,parres%mdmin
if(parres%verb.gt.0) then
       write(*,'(a,i5,1x,4(1x,1pe17.10),3(1x,i2))') ' # MD: it,enthalpy,pV,ekinat,ekinlat,nmax,nmin,mdmin ',&
             &itime,enthalpy,pressure*vol,ekinatom,ekinlat,nummax,nummin,parres%mdmin
       write(fn4,'(i4.4)') itime
       filename=trim(folder)//"posmd."//fn4//".ascii"
       units=units
       write(*,*) "# Writing the positions in MD: ",filename
       call write_atomic_file_ascii(parres,filename,nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &char_type(1:ntypat),ntypat,typat,fixat,fixlat,etot_in,pressure,enthalpy,en0000)
endif
!*********************************************************************

!        econs_max=max(econs_max,rkin+e_rxyz)
!        econs_min=min(econs_min,rkin+e_rxyz)
        econs_max=max(econs_max,rkin+enthalpy)
        econs_min=min(econs_min,rkin+enthalpy)
!        write(6,*) istep,e_rxyz-e_pos,nummax,nummin
       if (nummin.ge.parres%mdmin) then
          if (nummax.ne.nummin) &
               write(*,*) '# WARNING: nummin,nummax',nummin,nummax

          !write(67,*) " MD finished: exiting!"
          write(*,*) " MD finished: exiting!"
          exit
       endif


!Update the variables for next iteration
        accposprev=accposcur
        accposcur=accpospred
        acclatprev=acclatcur
        acclatcur=acclatpred
        velmatcur=velmatpred
        flatcur=flatpred
        fposcur=fpospred
        poscur=pospred
        vlatcur=vlatpred
        vposcur=vpospred
        latcur=latpred

        accvolprev=accvolcur
        accvolcur=accvolpred
        volcur=volpred
        vvolcur=vvolpred

     enddo 
!Adjust MD stepsize
!Minimum number of steps per crossed minimum is 15, average should be parres%nit_per_min
     
     if(parres%auto_dtion_md) then
!       dt_ratio=real(itime,8)/real(parres%mdmin,8) !old version version
       dt_ratio=real(itime,8)/real(nummin,8) 
       if(dt_ratio.lt.real(parres%nit_per_min,8)) then
         parres%dtion_md=parres%dtion_md*1.d0/1.1d0
       else
         parres%dtion_md=parres%dtion_md*1.1d0 
       endif
     parres%dtion_md=min(parres%dtion_md,dt*dt_ratio/15.d0)
     write(*,'(3(a,es10.2))') " # MD: steps per minium: ",dt_ratio,&
           &", parres%dtion_md set to: ",parres%dtion_md,", upper boundary: ",dt*dt_ratio/15.d0 
     endif
   

!MD stopped, now do relaxation
     write(*,'(a,i5,es15.7,es15.7)') ' # EXIT MD ',itime,enthalpy,etot_in

end subroutine
!**********************************************************************************************

subroutine acceleration_fire(pressure,accpos,acclat,accvol,vpos,vlat,vvol,strten,fcart,latvec,amass,latmass,f0inv,md_type,nat) 
use mod_interface
use global, only: fixlat
implicit none
integer:: iat,i,j,md_type,nat
real(8),dimension(3,nat):: accpos,vpos,fcart,fpos
real(8),dimension(3,3)  :: acclat,vlat,latvec,tmplat,pressure,a,velmat,sigma,lattrans,latdottrans,gdot,g,ginv,gtot,str_matrix
real(8),dimension(3,3)  :: term1,term2,term3,term4,term5,term5_1,term5_2,sigmatrans,f0inv
real(8):: amass(nat),latmass,crossp(3),strten(6),vol,vpostmp(3),volvel,trace3
real(8):: accvol,vvol,vol_1_3
!Get volume
           a=latvec
           vol= a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)+&
                a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)
!md_type=1
!Convert cartesian forces to reduced forces on atoms
        call invertmat(latvec,tmplat,3)
        do iat=1,nat
          fpos(:,iat)=matmul(tmplat,fcart(:,iat))
        enddo
!Get full stress matrix
           str_matrix(1,1)=strten(1)
           str_matrix(2,2)=strten(2)
           str_matrix(3,3)=strten(3)
           str_matrix(1,2)=strten(6)
           str_matrix(2,1)=strten(6)
           str_matrix(1,3)=strten(5)
           str_matrix(3,1)=strten(5)
           str_matrix(2,3)=strten(4)
           str_matrix(3,2)=strten(4)
!Update velocity part of stress tensor
           velmat=0.d0
!           do iat=1,nat
!              vpostmp=matmul(latvec,vpos(:,iat))
!              do i=1,3
!                 do j=1,3
!                 velmat(i,j)=velmat(i,j)+vpostmp(i)*vpostmp(j)*amass(iat) !Added mass on Oct 23 2015
!                 enddo
!              enddo
!           enddo
if(md_type.ge.1.and.md_type.le.3) then
!Compute sigma
        call cross_product(latvec(:,2),latvec(:,3),crossp); sigma(:,1)=crossp
        call cross_product(latvec(:,3),latvec(:,1),crossp); sigma(:,2)=crossp
        call cross_product(latvec(:,1),latvec(:,2),crossp); sigma(:,3)=crossp
!Compute the atomic acceleration
!For this we first need to calculate the matrix g^-1*dg/dt=ginv*gdot=gtot
          do i=1,3
            lattrans(:,i)=latvec(i,:)
            latdottrans(:,i)=vlat(i,:)
          enddo
          gdot=matmul(latdottrans,latvec)+matmul(lattrans,vlat)
          g=matmul(lattrans,latvec)
          call invertmat(g,ginv,3)
          gtot=matmul(ginv,gdot)
!Total acceleration
          do iat=1,nat
          accpos(:,iat)=1.d0/amass(iat)*fpos(:,iat)! - matmul(gtot,vpos(:,iat))
!          if(fixlat(7))  accpos(:,iat)=1.d0/amass(iat)*fpos(:,iat)
          enddo
elseif(md_type==4) then
!Andersen
!Be careful: the positions are not the real positions but scaled with vol**(1/3)
          vol_1_3=vol**(1.d0/3.d0)
          do iat=1,nat
              !accpos(:,iat)=(fcart(:,iat)/amass(iat)-2.d0*vpos(:,iat)*vvol/(vol*3.d0))/vol_1_3
              accpos(:,iat)=(fcart(:,iat)/amass(iat))/vol_1_3!-2.d0*vpos(:,iat)*vvol/(vol*3.d0))/vol_1_3
          enddo
!Update velocity part of stress tensor
           velmat=0.d0
!           do iat=1,nat
!              vpostmp=vol_1_3*vpos(:,iat)
!              do i=1,3
!                 do j=1,3
!                 velmat(i,j)=velmat(i,j)+vpostmp(i)*vpostmp(j)*amass(iat) !Added mass on Oct 23 2015
!                 enddo
!              enddo
!           enddo
endif
velmat=0.d0
if(md_type==1) then
!Compute the acceleration of the cell
        acclat=velmat/vol-str_matrix
!Here the pressure is applied
        acclat=acclat-pressure
!Scale with lattice mass
        acclat=acclat/latmass
!Multiply with sigma from left
        acclat=matmul(acclat,sigma)

elseif(md_type==2) then
!Fist term, same as in PR******************************
!Compute the acceleration of the cell
        term1=velmat/vol-str_matrix
!Here the pressure is applied
        term1=term1-pressure
!Scale with lattice volume and mass
        term1=term1/vol/latmass
!Combine it with the cell
        term1=matmul(term1,latvec)

!Second term*********************************************
        sigmatrans(:,1)=sigma(1,:)
        sigmatrans(:,2)=sigma(2,:)
        sigmatrans(:,3)=sigma(3,:)
        term2=matmul(sigmatrans,vlat)
        volvel=term2(1,1)+term2(2,2)+term2(3,3)
        term2=-2.d0*volvel/vol*vlat

!Third term**********************************************
!This term was taken from the moldyn manual. It is slightly different in the cleveland paper!!!
!        term3=matmul(sigmatrans,sigma)
!        term3=matmul(term3,latdottrans)
!        term3=matmul(term3,vlat)
!From cleveland
        term3=matmul(vlat,sigmatrans)
        term3=matmul(term3,sigma)
        term3=matmul(term3,latdottrans)
        trace3=term3(1,1)+term3(2,2)+term3(3,3)
        term3=1.d0/vol**2*trace3*latvec 

!Fourth term*********************************************
        term4=matmul(vlat,sigmatrans)
        term4=matmul(term4,vlat)
        term4=1.d0/vol*term4
 
!Fifth term**********************************************
        term5_1=matmul(vlat,sigmatrans)
        term5_1=matmul(term5_1,sigma)
        term5_1=matmul(term5_1,latdottrans)
        term5_2=matmul(sigma,latdottrans)
        term5_2=matmul(term5_2,vlat)
        term5_2=matmul(term5_2,sigmatrans)
        term5=1.d0/vol**2*matmul((term5_1-term5_2),latvec)
!Sum
        acclat=term1+term2+term3+term4+term5
elseif(md_type==3) then
!Compute the acceleration of the cell
        acclat=velmat/vol-str_matrix
!Here the pressure is applied
        acclat=acclat-pressure
!Scale with lattice mass
        acclat=acclat/latmass*vol**(4.d0/3.d0)
!Multiply with sigma from left
        acclat=matmul(acclat,sigma)
        acclat=matmul(acclat,f0inv)
elseif(md_type==4) then
!Andersen MD
!Compute the acceleration of the cell
        acclat=velmat/vol-str_matrix
!Here the pressure is applied
        acclat=acclat-pressure
!Scale with lattice mass
        acclat=acclat/latmass
!Compute the hydrostatic pressure stuff
        accvol=(acclat(1,1)+acclat(2,2)+acclat(3,3))/3.d0
else
stop "Wrong option in MD"
endif
end subroutine
!**********************************************************************************************
!**********************************************************************************************

subroutine acceleration(pressure,accpos,acclat,accvol,vpos,vlat,vvol,strten,fcart,latvec,amass,latmass,f0inv,md_type,nat) 
use mod_interface
use global, only: fixlat
implicit none
integer:: iat,i,j,md_type,nat
real(8),dimension(3,nat):: accpos,vpos,fcart,fpos
real(8),dimension(3,3)  :: acclat,vlat,latvec,tmplat,pressure,a,velmat,sigma,lattrans,latdottrans,gdot,g,ginv,gtot,str_matrix
real(8),dimension(3,3)  :: term1,term2,term3,term4,term5,term5_1,term5_2,sigmatrans,f0inv
real(8):: amass(nat),latmass,crossp(3),strten(6),vol,vpostmp(3),volvel,trace3
real(8):: accvol,vvol,vol_1_3
!Get volume
           a=latvec
           vol= a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)+&
                a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)
!md_type=1
!Convert cartesian forces to reduced forces on atoms
        call invertmat(latvec,tmplat,3)
        do iat=1,nat
          fpos(:,iat)=matmul(tmplat,fcart(:,iat))
        enddo
!Get full stress matrix
           str_matrix(1,1)=strten(1)
           str_matrix(2,2)=strten(2)
           str_matrix(3,3)=strten(3)
           str_matrix(1,2)=strten(6)
           str_matrix(2,1)=strten(6)
           str_matrix(1,3)=strten(5)
           str_matrix(3,1)=strten(5)
           str_matrix(2,3)=strten(4)
           str_matrix(3,2)=strten(4)
!Update velocity part of stress tensor
           velmat=0.d0
           do iat=1,nat
              vpostmp=matmul(latvec,vpos(:,iat))
              do i=1,3
                 do j=1,3
                 velmat(i,j)=velmat(i,j)+vpostmp(i)*vpostmp(j)*amass(iat) !Added mass on Oct 23 2015
                 enddo
              enddo
           enddo
if(md_type.ge.1.and.md_type.le.3) then
!Compute sigma
        call cross_product(latvec(:,2),latvec(:,3),crossp); sigma(:,1)=crossp
        call cross_product(latvec(:,3),latvec(:,1),crossp); sigma(:,2)=crossp
        call cross_product(latvec(:,1),latvec(:,2),crossp); sigma(:,3)=crossp
!Compute the atomic acceleration
!For this we first need to calculate the matrix g^-1*dg/dt=ginv*gdot=gtot
          do i=1,3
            lattrans(:,i)=latvec(i,:)
            latdottrans(:,i)=vlat(i,:)
          enddo
          gdot=matmul(latdottrans,latvec)+matmul(lattrans,vlat)
          g=matmul(lattrans,latvec)
          call invertmat(g,ginv,3)
          gtot=matmul(ginv,gdot)
!Total acceleration
          do iat=1,nat
          accpos(:,iat)=1.d0/amass(iat)*fpos(:,iat) - matmul(gtot,vpos(:,iat))
!          if(fixlat(7))  accpos(:,iat)=1.d0/amass(iat)*fpos(:,iat)
          enddo
elseif(md_type==4) then
!Andersen
!Be careful: the positions are not the real positions but scaled with vol**(1/3)
          vol_1_3=vol**(1.d0/3.d0)
          do iat=1,nat
              accpos(:,iat)=(fcart(:,iat)/amass(iat)-2.d0*vpos(:,iat)*vvol/(vol*3.d0))/vol_1_3
          enddo
!Update velocity part of stress tensor
           velmat=0.d0
           do iat=1,nat
              vpostmp=vol_1_3*vpos(:,iat)
              do i=1,3
                 do j=1,3
                 velmat(i,j)=velmat(i,j)+vpostmp(i)*vpostmp(j)*amass(iat) !Added mass on Oct 23 2015
                 enddo
              enddo
           enddo
endif

if(md_type==1) then
!Compute the acceleration of the cell
        acclat=velmat/vol-str_matrix
!Here the pressure is applied
        acclat=acclat-pressure
!Scale with lattice mass
        acclat=acclat/latmass
!Multiply with sigma from left
        acclat=matmul(acclat,sigma)

elseif(md_type==2) then
!Fist term, same as in PR******************************
!Compute the acceleration of the cell
        term1=velmat/vol-str_matrix
!Here the pressure is applied
        term1=term1-pressure
!Scale with lattice volume and mass
        term1=term1/vol/latmass
!Combine it with the cell
        term1=matmul(term1,latvec)

!Second term*********************************************
        sigmatrans(:,1)=sigma(1,:)
        sigmatrans(:,2)=sigma(2,:)
        sigmatrans(:,3)=sigma(3,:)
        term2=matmul(sigmatrans,vlat)
        volvel=term2(1,1)+term2(2,2)+term2(3,3)
        term2=-2.d0*volvel/vol*vlat

!Third term**********************************************
!This term was taken from the moldyn manual. It is slightly different in the cleveland paper!!!
!        term3=matmul(sigmatrans,sigma)
!        term3=matmul(term3,latdottrans)
!        term3=matmul(term3,vlat)
!From cleveland
        term3=matmul(vlat,sigmatrans)
        term3=matmul(term3,sigma)
        term3=matmul(term3,latdottrans)
        trace3=term3(1,1)+term3(2,2)+term3(3,3)
        term3=1.d0/vol**2*trace3*latvec 

!Fourth term*********************************************
        term4=matmul(vlat,sigmatrans)
        term4=matmul(term4,vlat)
        term4=1.d0/vol*term4
 
!Fifth term**********************************************
        term5_1=matmul(vlat,sigmatrans)
        term5_1=matmul(term5_1,sigma)
        term5_1=matmul(term5_1,latdottrans)
        term5_2=matmul(sigma,latdottrans)
        term5_2=matmul(term5_2,vlat)
        term5_2=matmul(term5_2,sigmatrans)
        term5=1.d0/vol**2*matmul((term5_1-term5_2),latvec)
!Sum
        acclat=term1+term2+term3+term4+term5
elseif(md_type==3) then
!Compute the acceleration of the cell
        acclat=velmat/vol-str_matrix
!Here the pressure is applied
        acclat=acclat-pressure
!Scale with lattice mass
        acclat=acclat/latmass*vol**(4.d0/3.d0)
!Multiply with sigma from left
        acclat=matmul(acclat,sigma)
        acclat=matmul(acclat,f0inv)
elseif(md_type==4) then
!Andersen MD
!Compute the acceleration of the cell
        acclat=velmat/vol-str_matrix
!Here the pressure is applied
        acclat=acclat-pressure
!Scale with lattice mass
        acclat=acclat/latmass
!Compute the hydrostatic pressure stuff
        accvol=(acclat(1,1)+acclat(2,2)+acclat(3,3))/3.d0
else
stop "Wrong option in MD"
endif
end subroutine

!**********************************************************************************************

subroutine stress_velocity(vpos,latvec,amass,nat,vpressure)
use mod_interface
implicit none
real(8):: velmat(3,3),vpostmp(3),latvec(3,3),vpos(3,nat),vpressure,a(3,3),vol,amass(nat)
integer:: iat,nat,i,j
a=latvec
vol= a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)+&
     a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)
velmat=0.d0
do iat=1,nat
   vpostmp=matmul(latvec,vpos(:,iat))
   do i=1,3
      do j=1,3
      velmat(i,j)=velmat(i,j)+vpostmp(i)*vpostmp(j)*amass(iat) !Added mass on Oct 23 2015
      enddo
   enddo
enddo
velmat=velmat/vol
vpressure=velmat(1,1)+velmat(2,2)+velmat(3,3)
vpressure=vpressure/3.d0
end subroutine

!**********************************************************************************************

!subroutine acceleration(pressure,accpos,acclat,vpos,vlat,strten,fcart,latvec,amass,latmass,f0inv,md_type,nat) 
!use global, only: fixlat
!implicit none
!integer:: iat,i,j,md_type
!real(8),dimension(3,nat):: accpos,vpos,fcart,fpos
!real(8),dimension(3,3)  :: acclat,vlat,latvec,tmplat,pressure,a,velmat,sigma,lattrans,latdottrans,gdot,g,ginv,gtot,str_matrix
!real(8),dimension(3,3)  :: term1,term2,term3,term4,term5,term5_1,term5_2,sigmatrans,f0inv
!real(8):: amass(nat),latmass,crossp(3),strten(6),vol,vpostmp(3),volvel,trace3
!!md_type=1
!!Convert cartesian forces to reduced forces on atoms
!        call invertmat(latvec,tmplat,3)
!        do iat=1,nat
!          fpos(:,iat)=matmul(tmplat,fcart(:,iat))
!        enddo
!!Get full stress matrix
!           str_matrix(1,1)=strten(1)
!           str_matrix(2,2)=strten(2)
!           str_matrix(3,3)=strten(3)
!           str_matrix(1,2)=strten(6)
!           str_matrix(2,1)=strten(6)
!           str_matrix(1,3)=strten(5)
!           str_matrix(3,1)=strten(5)
!           str_matrix(2,3)=strten(4)
!           str_matrix(3,2)=strten(4)
!!Get volume
!           a=latvec
!           vol= a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)+&
!                a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)
!!Update velocity part of stress tensor
!           velmat=0.d0
!           do iat=1,nat
!              vpostmp=matmul(latvec,vpos(:,iat))
!              do i=1,3
!                 do j=1,3
!                 velmat(i,j)=velmat(i,j)+vpostmp(i)*vpostmp(j)
!                 enddo
!              enddo
!           enddo
!!Compute sigma
!        call cross_product(latvec(:,2),latvec(:,3),crossp); sigma(:,1)=crossp
!        call cross_product(latvec(:,3),latvec(:,1),crossp); sigma(:,2)=crossp
!        call cross_product(latvec(:,1),latvec(:,2),crossp); sigma(:,3)=crossp
!!Compute the atomic acceleration
!!For this we first need to calculate the matrix g^-1*dg/dt=ginv*gdot=gtot
!          do i=1,3
!            lattrans(:,i)=latvec(i,:)
!            latdottrans(:,i)=vlat(i,:)
!          enddo
!          gdot=matmul(latdottrans,latvec)+matmul(lattrans,vlat)
!          g=matmul(lattrans,latvec)
!          call invertmat(g,ginv,3)
!          gtot=matmul(ginv,gdot)
!!Total acceleration
!          do iat=1,nat
!          accpos(:,iat)=1.d0/amass(iat)*fpos(:,iat) - matmul(gtot,vpos(:,iat))
!!          if(fixlat(7))  accpos(:,iat)=1.d0/amass(iat)*fpos(:,iat)
!          enddo
!if(md_type==1) then
!!Compute the acceleration of the cell
!        acclat=velmat/vol-str_matrix
!!Here the pressure is applied
!        acclat=acclat-pressure
!!Scale with lattice mass
!        acclat=acclat/latmass
!!Multiply with sigma from left
!        acclat=matmul(acclat,sigma)
!
!elseif(md_type==2) then
!!Fist term, same as in PR******************************
!!Compute the acceleration of the cell
!        term1=velmat/vol-str_matrix
!!Here the pressure is applied
!        term1=term1-pressure
!!Scale with lattice volume and mass
!        term1=term1/vol/latmass
!!Combine it with the cell
!        term1=matmul(term1,latvec)
!
!!Second term*********************************************
!        sigmatrans(:,1)=sigma(1,:)
!        sigmatrans(:,2)=sigma(2,:)
!        sigmatrans(:,3)=sigma(3,:)
!        term2=matmul(sigmatrans,vlat)
!        volvel=term2(1,1)+term2(2,2)+term2(3,3)
!        term2=-2.d0*volvel/vol*vlat
!
!!Third term**********************************************
!!This term was taken from the moldyn manual. It is slightly different in the cleveland paper!!!
!!        term3=matmul(sigmatrans,sigma)
!!        term3=matmul(term3,latdottrans)
!!        term3=matmul(term3,vlat)
!!From cleveland
!        term3=matmul(vlat,sigmatrans)
!        term3=matmul(term3,sigma)
!        term3=matmul(term3,latdottrans)
!        trace3=term3(1,1)+term3(2,2)+term3(3,3)
!        term3=1.d0/vol**2*trace3*latvec 
!
!!Fourth term*********************************************
!        term4=matmul(vlat,sigmatrans)
!        term4=matmul(term4,vlat)
!        term4=1.d0/vol*term4
! 
!!Fifth term**********************************************
!        term5_1=matmul(vlat,sigmatrans)
!        term5_1=matmul(term5_1,sigma)
!        term5_1=matmul(term5_1,latdottrans)
!        term5_2=matmul(sigma,latdottrans)
!        term5_2=matmul(term5_2,vlat)
!        term5_2=matmul(term5_2,sigmatrans)
!        term5=1.d0/vol**2*matmul((term5_1-term5_2),latvec)
!!Sum
!        acclat=term1+term2+term3+term4+term5
!elseif(md_type==3) then
!!Compute the acceleration of the cell
!        acclat=velmat/vol-str_matrix
!!Here the pressure is applied
!        acclat=acclat-pressure
!!Scale with lattice mass
!        acclat=acclat/latmass*vol**(4.d0/3.d0)
!!Multiply with sigma from left
!        acclat=matmul(acclat,sigma)
!        acclat=matmul(acclat,f0inv)
!else
!stop "Wrong option in MD"
!endif
!end subroutine
!


subroutine fpos_flat(pressure,fpos,flat,strten,fcart,latvec,md_type) 
!Computes the pure generalized forces on atom and cell (no contributions from velocities)
use mod_interface
use global, only: nat
implicit none
integer:: iat,i,j,md_type
real(8),dimension(3,nat):: fcart,fpos
real(8),dimension(3,3)  :: latvec,tmplat,pressure,a,velmat,sigma,flat,str_matrix
real(8),dimension(3,3)  :: term1,term2,term3,term4,term5,term5_1,term5_2,sigmatrans,f0inv
real(8):: amass(nat),latmass,crossp(3),strten(6),vol,vpostmp(3),volvel,trace3,vol_1_3
           flat=0.d0
!Get volume
           a=latvec
           vol= a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)+&
                a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)
!Get full stress matrix
           str_matrix(1,1)=strten(1)
           str_matrix(2,2)=strten(2)
           str_matrix(3,3)=strten(3)
           str_matrix(1,2)=strten(6)
           str_matrix(2,1)=strten(6)
           str_matrix(1,3)=strten(5)
           str_matrix(3,1)=strten(5)
           str_matrix(2,3)=strten(4)
           str_matrix(3,2)=strten(4)
!Compute sigma
        call cross_product(latvec(:,2),latvec(:,3),crossp); sigma(:,1)=crossp
        call cross_product(latvec(:,3),latvec(:,1),crossp); sigma(:,2)=crossp
        call cross_product(latvec(:,1),latvec(:,2),crossp); sigma(:,3)=crossp
!Compute the acceleration of the cell
        flat=-str_matrix
!Here the pressure is applied
        flat=flat-pressure
if(md_type.ne.4) then
!Multiply with sigma from left
        flat=matmul(flat,sigma)
!Convert cartesian forces to reduced forces on atoms
        call invertmat(latvec,tmplat,3)
        do iat=1,nat
          fpos(:,iat)=matmul(tmplat,fcart(:,iat))
        enddo
else
!Andersen
!Be careful: the positions are not the real positions but scaled with vol**(1/3)
          vol_1_3=vol**(1.d0/3.d0)
          do iat=1,nat
              fpos(:,iat)=fcart(:,iat)/vol_1_3
          enddo
!Andersen MD
!Compute the hydrostatic pressure stuff
        flat(1,1)=(flat(1,1)+flat(2,2)+flat(3,3))/3.d0
endif
end subroutine


subroutine ekin_at_lat(amass,latmass,latvec,vpos,vlat,ekinat,ekinlat,f0,md_type,nat)
use mod_interface
implicit none
integer:: iat,i,md_type,nat
real(8):: latvec(3,3),vpos(3,nat),vlat(3,3),ekinat,ekinlat,rkin,vposcurtmp(3),crossp(3),f0(3,3),vol
real(8):: latmass,amass(nat),lattrans(3,3),latdottrans(3,3),ekintrace(3,3),sigma(3,3),sigmatrans(3,3)
!Compute the kinetic energies:
  rkin=0.d0
  do iat=1,nat
     vposcurtmp=matmul(latvec,vpos(:,iat))
     rkin=rkin+amass(iat)*(vposcurtmp(1)**2+vposcurtmp(2)**2+vposcurtmp(3)**2)
  enddo
  ekinat=0.5d0*rkin
  rkin=0.d0
!  do i=1,3
!     rkin=rkin+(vlat(1,i)**2+vlat(2,i)**2+vlat(3,i)**2)*latmass
!  enddo
!  ekinlat=0.5d0*rkin
!For this we first need to calculate the matrix g^-1*dg/dt=ginv*gdot=gtot
          do i=1,3
            lattrans(:,i)=latvec(i,:)
            latdottrans(:,i)=vlat(i,:)
          enddo
if(md_type==1) then
  ekintrace=matmul(latdottrans,vlat)
  rkin=ekintrace(1,1)+ekintrace(2,2)+ekintrace(3,3)
  ekinlat=rkin*latmass*0.5d0
elseif(md_type==2) then
!sigma 
        call cross_product(latvec(:,2),latvec(:,3),crossp); sigma(:,1)=crossp
        call cross_product(latvec(:,3),latvec(:,1),crossp); sigma(:,2)=crossp
        call cross_product(latvec(:,1),latvec(:,2),crossp); sigma(:,3)=crossp
        sigmatrans(:,1)=sigma(1,:)
        sigmatrans(:,2)=sigma(2,:)
        sigmatrans(:,3)=sigma(3,:)
        ekintrace=matmul(vlat,sigmatrans)
        ekintrace=matmul(ekintrace,sigma)
        ekintrace=matmul(ekintrace,latdottrans)
  rkin=ekintrace(1,1)+ekintrace(2,2)+ekintrace(3,3)
  ekinlat=rkin*latmass*0.5d0
elseif(md_type==3) then
  call getvol(latvec,vol)
  ekintrace=matmul(latdottrans,f0)  
  ekintrace=matmul(ekintrace,vlat)  
  rkin=ekintrace(1,1)+ekintrace(2,2)+ekintrace(3,3)
  ekinlat=rkin*latmass*0.5d0/vol**(4.d0/3.d0)  
endif
!  write(*,*) "Ekinetic",ekinlat,rkin
end subroutine 


!**********************************************************************************************
subroutine ekin_at_lat_andersen(amass,latmass,latvec,vpos,vlat,vvol,ekinat,ekinlat,f0,md_type,nat)
use mod_interface
implicit none
integer:: iat,i,md_type,nat
real(8):: latvec(3,3),vpos(3,nat),vlat(3,3),ekinat,ekinlat,rkin,vposcurtmp(3),crossp(3),f0(3,3),vol,vvol
real(8):: latmass,amass(nat),lattrans(3,3),latdottrans(3,3),ekintrace(3,3),sigma(3,3),sigmatrans(3,3),vol_1_3
rkin=0.d0
if(md_type.lt.4) then
!Compute the kinetic energies:
  do iat=1,nat
     vposcurtmp=matmul(latvec,vpos(:,iat))
     rkin=rkin+amass(iat)*(vposcurtmp(1)**2+vposcurtmp(2)**2+vposcurtmp(3)**2)
  enddo
!  do i=1,3
!     rkin=rkin+(vlat(1,i)**2+vlat(2,i)**2+vlat(3,i)**2)*latmass
!  enddo
!  ekinlat=0.5d0*rkin
!For this we first need to calculate the matrix g^-1*dg/dt=ginv*gdot=gtot
          do i=1,3
            lattrans(:,i)=latvec(i,:)
            latdottrans(:,i)=vlat(i,:)
          enddo
else
  call getvol(latvec,vol)
  vol_1_3=vol**(1.d0/3.d0)
  do iat=1,nat
     vposcurtmp=vpos(:,iat)*vol_1_3
     rkin=rkin+amass(iat)*(vposcurtmp(1)**2+vposcurtmp(2)**2+vposcurtmp(3)**2)
  enddo
endif
  ekinat=0.5d0*rkin
  rkin=0.d0

if(md_type==1) then
  ekintrace=matmul(latdottrans,vlat)
  rkin=ekintrace(1,1)+ekintrace(2,2)+ekintrace(3,3)
  ekinlat=rkin*latmass*0.5d0
elseif(md_type==2) then
!sigma 
        call cross_product(latvec(:,2),latvec(:,3),crossp); sigma(:,1)=crossp
        call cross_product(latvec(:,3),latvec(:,1),crossp); sigma(:,2)=crossp
        call cross_product(latvec(:,1),latvec(:,2),crossp); sigma(:,3)=crossp
        sigmatrans(:,1)=sigma(1,:)
        sigmatrans(:,2)=sigma(2,:)
        sigmatrans(:,3)=sigma(3,:)
        ekintrace=matmul(vlat,sigmatrans)
        ekintrace=matmul(ekintrace,sigma)
        ekintrace=matmul(ekintrace,latdottrans)
  rkin=ekintrace(1,1)+ekintrace(2,2)+ekintrace(3,3)
  ekinlat=rkin*latmass*0.5d0
elseif(md_type==3) then
  call getvol(latvec,vol)
  ekintrace=matmul(latdottrans,f0)  
  ekintrace=matmul(ekintrace,vlat)  
  rkin=ekintrace(1,1)+ekintrace(2,2)+ekintrace(3,3)
  ekinlat=rkin*latmass*0.5d0/vol**(4.d0/3.d0)  
elseif(md_type==4) then
  rkin=vvol**2*latmass
  ekinlat=0.5d0*rkin
endif
!  write(*,*) "Ekinetic",ekinlat,rkin
end subroutine 


!**********************************************************************************************
subroutine MD_ANDERSEN_MHM     (parini,parres,latvec_in,xred_in,fcart_in,strten_in,vel_in,vel_lat_in,vvol_in,etot_in,iprec,counter,folder)
 use mod_interface
 use global, only: target_pressure_habohr,target_pressure_gpa,nat,ntypat,znucl,amu,amutmp,typat
 use global, only: char_type,units,usewf_md
 use global, only: fixat,fixlat
 use defs_basis
 use interface_code
 use mod_parini, only: typ_parini
!Implements constant cell shape MD from Andersen,
!I havent checked the equations of motion myself, taken at the moment from
!http://www.grs-sim.de/cms/upload/Carloni/Presentations/Strodel4.pdf
!The coordinates of interest are the xred of atoms and the cell volume vol
implicit none
 type(typ_parini), intent(in):: parini
 type(typ_parini), intent(inout):: parres
 real(8) :: latvec_in(3,3),xred_in(3,nat),vel_in(3,nat),fcart_in(3,nat),etot_in,strten_in(6),vel_lat_in,vvol_in
!*********************************************************************
!Variables for my MD part
! real(8), parameter :: Ha_eV=27.21138386_dp ! 1 Hartree, in eV
! real(8), parameter :: kb_HaK=8.617343d-5/Ha_eV ! Boltzmann constant in Ha/K
! real(8), parameter :: amu_emass=1.660538782d-27/9.10938215d-31 ! 1 atomic mass unit, in electronic mass
 real(8):: amass(nat)
 real(8):: pressure_ener
 real(8):: pressure_md(3,3)
 real(8):: pressure
 real(8):: vol
 real(8):: vol0
 real(8):: volprev
 real(8),dimension(3,nat):: xcart
 real(8),dimension(3,nat):: fposcur
 real(8),dimension(3,nat):: accposcur
 real(8),dimension(3,nat):: accpospred
 real(8),dimension(3,nat):: accposprev
 real(8),dimension(3,nat):: fpospred
 real(8),dimension(3,nat):: vpospred
 real(8),dimension(3,nat):: poscur
 real(8),dimension(3,nat):: vxyz
 real(8),dimension(3,nat):: vposcur
 real(8),dimension(3,nat):: pospred
 real(8),dimension(3,nat):: dxred
 real(8),dimension(3,3):: latvec
 real(8),dimension(3,3):: latcur
 real(8),dimension(3,3):: latvec0
 real(8),dimension(3,3):: latinv
 real(8),dimension(3,3):: unitmat
 real(8),dimension(3,3):: elmatcur
 real(8),dimension(3,3):: acclatcur
 real(8),dimension(3,3):: acclatpred
 real(8),dimension(3,3):: acclatprev
 real(8),dimension(3,3):: dlatvec
 real(8),dimension(3,3):: flatcur
 real(8),dimension(3,3):: flatpred
 real(8):: accvolcur
 real(8):: accvolpred
 real(8):: accvolprev
 real(8):: fvolcur
 real(8):: fvolpred
 real(8),dimension(3,3):: velmatcur
 real(8),dimension(3,3):: tmplat
 real(8),dimension(3,3):: tmplatt
 real(8),dimension(3,3):: tpred
 real(8),dimension(3,3):: f0
 real(8),dimension(3,3):: f0inv
 real(8),dimension(3,3):: vlatcur
 real(8),dimension(3,3):: vlatpred
 real(8):: volcur
 real(8):: volpred
 real(8):: vvolcur
 real(8),dimension(3,3):: lattrans
 real(8),dimension(3,3):: latdottrans
 real(8),dimension(3,3):: a 
 real(8),dimension(3,3):: g
 real(8),dimension(3,3):: gdot
 real(8),dimension(3,3):: ginv
 real(8),dimension(3,3):: gtot
 real(8),dimension(3,3):: sigma
 real(8),dimension(3,3):: sigmatrans
 real(8),dimension(3,3):: latpred
 real(8),dimension(3,3):: tcur
 real(8):: vvolpred
 real(8):: volpred_1_3
 real(8):: vol_1_3_inv
 real(8):: vol_1_3
 real(8),dimension(3,3):: velmatpred
 real(8),dimension(3,3):: vlat
 real(8),dimension(3,3):: str_matrix
 real(8):: vposcurtmp(3)
 real(8):: vpospredtmp(3)
 real(8):: crossp(3)
 real(8):: velcm(3)
 real(8):: latmass
 real(8):: latmassinv
 real(8):: latmass0
 real(8):: ekinatom
 real(8):: ekinatom_prev
 real(8):: ekinlat
 real(8):: ekinlat_prev
 real(8):: rkin
 real(8):: enthalpy
 real(8):: enmin1
 real(8):: enmin2
 real(8):: ent_pos_0
 real(8):: en0000
 real(8):: e_rxyz
 real(8):: econs_max
 real(8):: econs_min
 real(8):: torquenrm 
 real(8):: ecut_tmp
 real(8):: toldff_tmp
 real(8):: counter
 real(8):: dt
 real(8):: dt_ratio
 integer:: i
 integer:: j
 integer:: iat
 integer:: itime
 integer:: nummax
 integer:: nummin
 integer:: iprec 
 integer:: options 
 integer:: md_type
 logical:: getwfk

 character(40)::filename,folder
 character(4) ::fn4
 write(*,*) "# PARAMETERS: Ha_eV,kb_HaK,amu_emass -->",Ha_eV,kb_HaK,amu_emass
!MD-type: 1 for PR and 2 for Cleveland and 3 for Wentzcovitch
  md_type=parres%md_algo


!Assign masses to each atom (for MD)
 do iat=1,nat
   amass(iat)=amu_emass*amu(typat(iat))
   write(*,'(a,i5,2(1x,es15.7))') " # MD: iat, AMU, EM: ", iat, amu(typat(iat)),amass(iat)
 enddo
!******************************************************************
!Here we split the routine in the Abinit native part and my new implentation
 write(*,'(a)') ' # Entering standalone ANDERSEN MD  '
!The "reduced" coordinates in Andersen are quite different from the ones in PR
!Set temporary variables, initially
  vxyz(:,:)=vel_in(:,:)
  pressure=target_pressure_habohr
  unitmat=0.d0
  do i=1,3
    unitmat(i,i)=1.d0
  enddo
  latmass0=parini%bmass*amu_emass !This means we use the barostat mass as the lattice mass (in ELECTRON MASS)
  vlat=vel_lat_in  !The initial cell volume velocity

  itime=0
  dt=parres%dtion_md

!Set options=1 for Velocity Verlet of cell dynamics
!Set options=2 for Normal Verlet of cell dynamics
!Set options=3 for Beeman integration scheme, corrector-predictor
  options=parres%md_integrator
  if(options.lt.1.or.options.gt.3) stop "Wrong algo option"


  write(*,'(a,i3)') " # MD Algorithm: ANDERSEN, MD Integrator: ",options

!Now we run my implementation of ANDERSEN MD
pressure_md=0.d0
pressure_md(1,1)=1.d0;pressure_md(2,2)=1.d0;pressure_md(3,3)=1.d0
pressure_ener=0.d0;pressure_md=pressure_md*pressure  !Here the pressure is not passed to the energyandforces, so we move on the ENERGY surface

!Transform the initial velocities given in cartesian coordinates into velocities of the intenal coordinates
       latvec=latvec_in
       call getvol(latvec,vol)
       vol_1_3=vol**(1.d0/3.d0) 
       vol_1_3_inv=1.d0/vol_1_3 
       do iat=1,nat
          vposcur(:,iat)=vxyz(:,iat)*vol_1_3_inv
       enddo

!Keep track of volume
!In ANDERSEN, latvec0 is normed to have volume 0, so computing latvec0*vol^(1/3) will always give the current lattice vector
        vol0=vol
        latvec0=latvec/vol**(1.d0/3.d0)
        volprev=vol

!!This is for the original volume-mass factor
        latmass=latmass0
        latmassinv=1.d0/latmass
        write(*,*) "Latmasse",latmass

!Initialize internal coordinates and lattice vectors
        call rxyz_int2cart(latvec,xred_in,xcart,nat)
        poscur=xcart*vol_1_3_inv           !call rxyz_cart2int(latvec,poscur,rxyz,nat)
        volcur=vol
        vvolcur=vvol_in

!MHM: initialize variable to track the minima/maxima***********
    write(*,*) '# MINHOP start MD'
    nummax=0
    nummin=0
    enmin1=0.d0
    en0000=0.d0
    econs_max=-1.d100
    econs_min=1.d100
!**************************************************************

!Compute initial kinetic energies
  latcur=latvec0*vol_1_3
  rkin=0.d0
  do iat=1,nat
     vposcurtmp=vposcur(:,iat)*vol_1_3
     rkin=rkin+amass(iat)*(vposcurtmp(1)**2+vposcurtmp(2)**2+vposcurtmp(3)**2)
  enddo
  ekinatom=0.5d0*rkin
  rkin=vvolcur**2*latmass
  ekinlat=0.5d0*rkin


!MHM: Write output to file in every step***********************************
!INITIAL STEP, STILL THE SAME STRUCTURE AS INPUT
       write(*,*) "Pressure, Energy",pressure,etot_in
       call get_enthalpy(latvec_in,etot_in,pressure,enthalpy)
       ent_pos_0=enthalpy
       en0000=enthalpy-ent_pos_0
       write(*,'(a,i5,1x,4(1x,1pe17.10),3(1x,i2))') ' MD: it,enthalpy,pV,ekinat,ekinlat,nmax,nmin,mdmin ',&
             &itime,enthalpy,pressure*vol,ekinatom,ekinlat,nummax,nummin,parres%mdmin
       write(fn4,'(i4.4)') itime
       filename=trim(folder)//"posmd."//fn4//".ascii"
       units=units
       write(*,*) "# Writing the positions in MD:",filename
       call write_atomic_file_ascii(parini,filename,nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &char_type(1:ntypat),ntypat,typat,fixat,fixlat,etot_in,pressure,enthalpy,en0000)
!*********************************************************************




call acceleration(pressure_md,accposcur,acclatcur,accvolcur,vposcur,&
     &vlatcur,vvolcur,strten_in,fcart_in,latcur,amass,latmass,f0inv,md_type,nat) 
        accposprev=accposcur
!        acclatprev=acclatcur
        accvolprev=accvolcur

!Eliminate the unnecessary atmic and cell components
call elim_fixed_at(nat,accposcur);call elim_fixed_at(nat,accposprev)!;call elim_fixed_lat(latcur,acclatcur);call elim_fixed_lat(latcur,acclatprev)
call elim_fixed_at(nat,vposcur)!;call elim_fixed_lat(latcur,vlatcur)


        do itime=1,parini%nmd_dynamics

!          if(itime.ne.1) e_rxyz=enthalpy  !e_rxyz=e_rxyz+pressure_md(1,1)*vol   
!Check the torque on the cell for rotation
          call torque_cell(latcur,vlatcur,torquenrm)
!Now perform the step
          velcm=0.d0
          rkin=0.d0
          if(options==1) then
!For velocity verlet
!             latpred(:,:)=latcur(:,:) + dt*vlatcur(:,:) + (.5d0*dt*dt*acclatcur)
             dlatvec(:,:)=dt*vlatcur(:,:) + (.5d0*dt*dt*acclatcur)
          do iat=1,nat
              dxred(:,iat)=dt*vposcur(:,iat) + 0.5d0*dt*dt*accposcur(:,iat)
!             pospred(:,iat)=poscur(:,iat) + dt*vposcur(:,iat) + dt*dt*accposcur(:,iat)
          enddo
          call propagate(nat,poscur,latcur,dxred,dlatvec,pospred,latpred)
          do iat=1,nat
             vpospred(:,iat)=(pospred(:,iat)-poscur(:,iat))/dt   !Update dummy variable vposcur for next iteration
          enddo
            if(md_type==4) then
               volpred=volcur+dt*vvolcur+0.5d0*dt*dt*accvolcur       
            endif
          elseif(options==2) then
!Instead if normal verlet is used:
!              latpred(:,:)=latcur(:,:) + dt*vlatcur + dt*dt*acclatcur
              dlatvec(:,:)=dt*vlatcur + dt*dt*acclatcur
          do iat=1,nat
              dxred(:,iat)=dt*vposcur(:,iat) + dt*dt*accposcur(:,iat)
!             pospred(:,iat)=poscur(:,iat) + dt*vposcur(:,iat) + dt*dt*accposcur(:,iat)
          enddo
          call propagate(nat,poscur,latcur,dxred,dlatvec,pospred,latpred)
          vlatpred(:,:)=(latpred-latcur)/dt
          do iat=1,nat
              vpospred(:,iat)=(pospred(:,iat)-poscur(:,iat))/dt   !Update dummy variable vposcur for next iteration
          enddo
            if(md_type==4) then
               volpred=volcur+dt*vvolcur+dt*dt*accvolcur       
               vvolpred=(volpred-volcur)/dt
            endif
          elseif(options==3) then
!Predictor part of Beeman for cell
!              latpred(:,:)=latcur(:,:) + dt*vlatcur + dt*dt/6.d0*(4.d0*acclatcur(:,:)-acclatprev(:,:))
              dlatvec(:,:)=dt*vlatcur + dt*dt/6.d0*(4.d0*acclatcur(:,:)-acclatprev(:,:))
          do iat=1,nat
!Predictor part of Beeman
!             pospred(:,iat)=poscur(:,iat) + dt*vposcur(:,iat) + dt*dt/6.d0*(4.d0*accposcur(:,iat)-accposprev(:,iat))             
             dxred(:,iat)=dt*vposcur(:,iat) + dt*dt/6.d0*(4.d0*accposcur(:,iat)-accposprev(:,iat))             
          enddo
          call propagate(nat,poscur,latcur,dxred,dlatvec,pospred,latpred)
          vlatpred(:,:)=vlatcur(:,:) + 0.5d0*(3.d0*acclatcur(:,:)-acclatprev(:,:))  
          do iat=1,nat
             vpospred(:,iat)=vposcur(:,iat)+0.5d0*(3.d0*accposcur(:,iat)-accposprev(:,iat))
          enddo
            if(md_type==4) then
               volpred=volcur+dt*vvolcur+dt*dt/6.d0*(4.d0*accvolcur-accvolprev)       
               vvolpred=volcur+0.5d0*(3.d0*accvolcur-accvolprev)
            endif
          endif

if(md_type==4) then
  do iat=1,nat
     vposcurtmp=vposcur(:,iat)*vol_1_3
     rkin=rkin+amass(iat)*(vposcurtmp(1)**2+vposcurtmp(2)**2+vposcurtmp(3)**2)
     velcm=velcm+vposcurtmp*amass(iat)
  enddo
else
  do iat=1,nat
             vposcurtmp=matmul(latcur,vposcur(:,iat))
             rkin=rkin+amass(iat)*(vposcurtmp(1)**2+vposcurtmp(2)**2+vposcurtmp(3)**2)
             velcm=velcm+vposcurtmp*amass(iat)
  enddo
endif


!!!!Now calculate steps for the lattice according to the velocity verlet algorithm
!!!          if(options==1) then
!!!          elseif(options==2) then
!!!          elseif(options==3) then
!!!          else
!!!              stop "Wrong option for the Cell part in MD"
!!!          endif



!Kinetic energy of atoms
          ekinatom=0.5d0*rkin
          if(parini%verb.gt.0) write(*,'(a,es15.7)')"Velocity of CM: ", sqrt(velcm(1)**2+velcm(2)**2+velcm(3)**2)
!Kinetic energy according to Parrinello Rahman
          rkin=0.d0
if(md_type==4) then
  rkin=vvolcur**2*latmass
  ekinlat=0.5d0*rkin
else
          do i=1,3
             rkin=rkin+(vlatcur(1,i)**2+vlatcur(2,i)**2+vlatcur(3,i)**2)*latmass
          enddo
          ekinlat=0.5d0*rkin
endif
          rkin=ekinlat+ekinatom
if(parini%verb.gt.0) write(*,'(a,3(1x,es15.7))') " # Torquenrm, Ekin, Enthalpy: ",torquenrm, rkin,enthalpy
!Update counter
          enmin2=enmin1
          enmin1=en0000

!Here we perform the force call
!****************************************************************************************************************        
!****************************************************************************************************************        
!Get the current acell_in and rprim_in, and also the atomic positions
if(md_type==4) then
       volpred_1_3=volpred**(1.d0/3.d0)
       latvec_in=latvec0*volpred_1_3
       latpred=latvec_in
       xcart=pospred*volpred_1_3
       call rxyz_cart2int(latvec_in,xred_in,xcart,nat)
else
       xred_in=pospred 
       vel_in=0.d0
       latvec_in=latpred
endif
       if(itime==1)  then
           getwfk=.false.
       elseif(itime.ne.1.and.usewf_md) then
           getwfk=.true.
       else
           getwfk=.false.
       endif
       call get_energyandforces_single(parini,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,getwfk)
!Single point calculation finished. Now returning to MD part. All output variables are now in *_in
!****************************************************************************************************************        

!For velocity verlet of cell
        if(options==1) then
           call acceleration(pressure_md,accpospred,acclatpred,accvolpred,vpospred,&
                &vlatcur,vvolcur,strten_in,fcart_in,latpred,amass,latmass,f0inv,md_type,nat) 
           vlatpred(:,:)=vlatcur(:,:)+0.5d0*dt*(acclatpred+acclatcur)
           if(md_type==4) vvolpred=vvolcur+0.5d0*dt*(accvolpred+accvolcur)
        elseif(options==3) then
!Corrector part of Beeman. Note that a fixed number of iterations are used!!!
           call ekin_at_lat_andersen(amass,latmass,latpred,vpospred,vlatpred,vvolpred,ekinatom_prev,ekinlat_prev,f0,md_type,nat)
           do i=1,5
             call acceleration(pressure_md,accpospred,acclatpred,accvolpred,vpospred,&
                  &vlatpred,vvolpred,strten_in,fcart_in,latpred,amass,latmass,f0inv,md_type,nat) 
             do iat=1,nat
             vpospred(:,iat)=vposcur(:,iat)+dt/6.d0*(2.d0*accpospred(:,iat)+5.d0*accposcur(:,iat)-accposprev(:,iat))
             enddo
             vlatpred=vlatcur+dt/6.d0*(2.d0*acclatpred+5.d0*acclatcur-acclatprev)
             if(md_type==4) vvolpred=vvolcur+dt/6.d0*(2.d0*accvolpred+5.d0*accvolcur-accvolprev)
             call ekin_at_lat_andersen(amass,latmass,latpred,vpospred,vlatpred,vvolpred,ekinatom,ekinlat,f0,md_type,nat)
if(parini%verb.gt.1)write(* ,'(a,i5,1x,es15.5,es15.5)') " # Beeman: Corrector relative convergence: it, atoms, lattice: ",&
                  & i,abs(ekinatom-ekinatom_prev)/ekinatom,abs(ekinlat-ekinlat_prev)/max(ekinlat,1.d-100)
        !     write(67,'(a,i5,1x,es15.5,es15.5)') " # Beeman: Corrector relative convergence: it, atoms, lattice: ",&
        !          & i,abs(ekinatom-ekinatom_prev)/ekinatom,abs(ekinlat-ekinlat_prev)/max(ekinlat,1.d-100)
             ekinlat_prev=ekinlat
             ekinatom_prev=ekinatom
           enddo 
        endif  
call acceleration(pressure_md,accpospred,acclatpred,accvolpred,vpospred,&
     &vlatpred,vvolpred,strten_in,fcart_in,latpred,amass,latmass,f0inv,md_type,nat) 

!Eliminate the unnecessary atmic and cell components
call elim_fixed_at(nat,accpospred);call elim_fixed_at(nat,vpospred)!;call elim_fixed_lat(latcur,acclatpred);call elim_fixed_lat(latcur,vlatpred)

!Compute the "predicted" kinetic energies:
call ekin_at_lat_andersen(amass,latmass,latpred,vpospred,vlatpred,vvolpred,ekinatom,ekinlat,f0,md_type,nat)


!MHM: Write output to file in every step***********************************
       counter=real(itime,8)
       call get_enthalpy(latvec_in,etot_in,pressure,enthalpy)
       call getvol(latvec_in,vol)
       en0000=enthalpy-ent_pos_0
       if (itime >= 3 .and. enmin1 > enmin2 .and. enmin1 > en0000)  nummax=nummax+1
       if (itime >= 3 .and. enmin1 < enmin2 .and. enmin1 < en0000)  nummin=nummin+1
!       write(67,'(a,i5,1x,4(1x,1pe17.10),3(1x,i2))') ' MD: it,enthalpy,pV,ekinat,ekinlat,nmax,nmin,mdmin ',&
!             &itime,enthalpy,pressure*vol,ekinatom,ekinlat,nummax,nummin,parres%mdmin
       write(*,'(a,i5,1x,4(1x,1pe17.10),3(1x,i2))') ' # MD: it,enthalpy,pV,ekinat,ekinlat,nmax,nmin,mdmin ',&
             &itime,enthalpy,pressure*vol,ekinatom,ekinlat,nummax,nummin,parres%mdmin
       write(fn4,'(i4.4)') itime
       filename=trim(folder)//"posmd."//fn4//".ascii"
       units=units
       write(*,*) "# Writing the positions in MD: ",filename
       call write_atomic_file_ascii(parini,filename,nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &char_type(1:ntypat),ntypat,typat,fixat,fixlat,etot_in,pressure,enthalpy,en0000)
!*********************************************************************

!        econs_max=max(econs_max,rkin+e_rxyz)
!        econs_min=min(econs_min,rkin+e_rxyz)
        econs_max=max(econs_max,rkin+enthalpy)
        econs_min=min(econs_min,rkin+enthalpy)
!        write(6,*) istep,e_rxyz-e_pos,nummax,nummin
       if (nummin.ge.parres%mdmin) then
          if (nummax.ne.nummin) &
               write(*,*) '# WARNING: nummin,nummax',nummin,nummax

!          write(67,*) " MD finished: exiting!"
          write(*,*) " MD finished: exiting!"
          exit
       endif


!Update the variables for next iteration
        accposprev=accposcur
        accposcur=accpospred
        acclatprev=acclatcur
        acclatcur=acclatpred
        velmatcur=velmatpred
        flatcur=flatpred
        fposcur=fpospred
        poscur=pospred
        vlatcur=vlatpred
        vposcur=vpospred
        latcur=latpred

        accvolprev=accvolcur
        accvolcur=accvolpred
        volcur=volpred
        vvolcur=vvolpred
 
!        e_rxyz=e_rxyz_pred
     enddo 
!Adjust MD stepsize
!Minimum number of steps per crossed minimum is 15, average should be parini%nit_per_min
     if(parini%auto_dtion_md) then
!       dt_ratio=real(itime,8)/real(parres%mdmin,8) !old version version
       dt_ratio=real(itime,8)/real(nummin,8) 
       if(dt_ratio.lt.real(parini%nit_per_min,8)) then
         parres%dtion_md=parres%dtion_md*1.d0/1.1d0
       else
         parres%dtion_md=parres%dtion_md*1.1d0 
       endif
     parres%dtion_md=min(parres%dtion_md,dt*dt_ratio/15.d0)
     write(*,'(3(a,es10.2))') " # MD: steps per minium: ",dt_ratio,&
           &", parres%dtion_md set to: ",parres%dtion_md,", upper boundary: ",dt*dt_ratio/15.d0 
     endif
   

!MD stopped, now do relaxation
     write(*,'(a,i5,es15.7,es15.7)') ' # EXIT MD ',itime,enthalpy,etot_in

end subroutine
!**********************************************************************************************

!**********************************************************************************************
subroutine MD_PR_MHM_OLD    (parini,parres,latvec_in,xred_in,fcart_in,strten_in,vel_in,vel_lat_in,etot_in,iprec,counter,folder)
 use mod_interface
 use global, only: target_pressure_habohr,target_pressure_gpa,nat,ntypat,znucl,amu,amutmp,typat
 use global, only: char_type,units,usewf_md,fixat,fixlat
 use defs_basis
 use interface_code
 use mod_parini, only: typ_parini
implicit none
 type(typ_parini), intent(in):: parini
 type(typ_parini), intent(inout):: parres
 real(8) :: latvec_in(3,3),xred_in(3,nat),vel_in(3,nat),fcart_in(3,nat),etot_in,strten_in(6),vel_lat_in(3,3)
!*********************************************************************
!Variables for my MD part
! real(8), parameter :: Ha_eV=27.21138386_dp ! 1 Hartree, in eV
! real(8), parameter :: kb_HaK=8.617343d-5/Ha_eV ! Boltzmann constant in Ha/K
! real(8), parameter :: amu_emass=1.660538782d-27/9.10938215d-31 ! 1 atomic mass unit, in electronic mass
 real(8):: amass(nat)
 real(8):: pressure_ener
 real(8):: pressure_md(3,3)
 real(8):: pressure
 real(8):: vol
 real(8):: vol0
 real(8):: volprev
 real(8),dimension(3,nat):: fposcur
 real(8),dimension(3,nat):: fpospred
 real(8),dimension(3,nat):: vpospred
 real(8),dimension(3,nat):: poscur
 real(8),dimension(3,nat):: vxyz
 real(8),dimension(3,nat):: vposcur
 real(8),dimension(3,nat):: pospred
 real(8),dimension(3,3):: latvec
 real(8),dimension(3,3):: latinv
 real(8),dimension(3,3):: unitmat
 real(8),dimension(3,3):: elmatcur
 real(8),dimension(3,3):: flatcur
 real(8),dimension(3,3):: flatpred
 real(8),dimension(3,3):: velmatcur
 real(8),dimension(3,3):: tmplat
 real(8),dimension(3,3):: tmplatt
 real(8),dimension(3,3):: tpred
 real(8),dimension(3,3):: f0
 real(8),dimension(3,3):: f0inv
 real(8),dimension(3,3):: latcur
 real(8),dimension(3,3):: vlatcur
 real(8),dimension(3,3):: vlat
 real(8),dimension(3,3):: lattrans
 real(8),dimension(3,3):: latdottrans
 real(8),dimension(3,3):: a 
 real(8),dimension(3,3):: g
 real(8),dimension(3,3):: gdot
 real(8),dimension(3,3):: ginv
 real(8),dimension(3,3):: gtot
 real(8),dimension(3,3):: sigma
 real(8),dimension(3,3):: latpred
 real(8),dimension(3,3):: tcur
 real(8),dimension(3,3):: vlatpred
 real(8),dimension(3,3):: velmatpred
 real(8),dimension(3,3):: str_matrix
 real(8):: vposcurtmp(3)
 real(8):: vpospredtmp(3)
 real(8):: crossp(3)
 real(8):: velcm(3)
 real(8):: latmass
 real(8):: latmassinv
 real(8):: latmass0
 real(8):: ekinatom
 real(8):: ekinlat
 real(8):: rkin
 real(8):: enthalpy
 real(8):: enmin1
 real(8):: enmin2
 real(8):: ent_pos_0
 real(8):: en0000
 real(8):: e_rxyz
 real(8):: econs_max
 real(8):: econs_min
 real(8):: torquenrm 
 real(8):: ecut_tmp
 real(8):: toldff_tmp
 real(8):: counter
 real(8):: dt
 integer:: i
 integer:: j
 integer:: iat
 integer:: itime
 integer:: nummax
 integer:: nummin
 integer:: options
 integer:: iprec 
 logical:: getwfk

 character(40)::filename,folder
 character(4) ::fn4
 write(*,*) "# PARAMETERS: Ha_eV,kb_HaK,amu_emass -->",Ha_eV,kb_HaK,amu_emass


!Assign masses to each atom (for MD)
 do iat=1,nat
   amass(iat)=amu_emass*amu(typat(iat))
   write(*,'(a,i5,2(1x,es15.7))') " # MD: iat, AMU, EM: ", iat, amu(typat(iat)),amass(iat)
 enddo
!******************************************************************
!NEW VERISON OF MD: Directly implemented, simplest Parrinello-Rahman MD

!Here we split the routine in the Abinit native part and my new implentation
 write(*,'(a)') ' # Entering standalone PR MD '
!Set temporary variables, initially
  vxyz(:,:)=vel_in(:,:)
  pressure=target_pressure_habohr
  unitmat=0.d0
  do i=1,3
    unitmat(i,i)=1.d0
  enddo
  latmass0=parini%bmass*amu_emass !This means we use the barostat mass as the lattice mass (in ELECTRON MASS)
  vlat=vel_lat_in  !The initial cell velocity
  itime=0
  dt=parres%dtion_md

!Set options=1 for Velocity Verlet of cell dynamics
!Set options=2 for Normal Verlet of cell dynamics
  options=2  


!Now we run my implementation of PR MD
pressure_md=0.d0
pressure_md(1,1)=1.d0;pressure_md(2,2)=1.d0;pressure_md(3,3)=1.d0
pressure_ener=0.d0;pressure_md=pressure_md*pressure  !Here the pressure is not passed to the energyandforces, so we move on the ENERGY surface

!Transform the initial velocities given in cartesian coordinates into velocities of the intenal coordinates
!       call acell_rprim2latvec(latvec,acell_in,rprim_in)  !First we need to get latvec
       latvec=latvec_in
       call invertmat(latvec,latinv,3)
       do iat=1,nat
          vposcur(:,iat)=matmul(latinv,vxyz(:,iat))
       enddo

!Creating the initial velocity matrix as part of the stress tensor       
        velmatcur=0.d0
        do iat=1,nat
           do i=1,3
              do j=1,3
              velmatcur(i,j)=velmatcur(i,j)+vxyz(i,iat)*vxyz(j,iat)
              enddo
           enddo
        enddo
        a=latvec
        vol= a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)+&
             a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)

!ATTENTION: THIS IS FOR THE ASSUMPTION THAT THE STRUCTURE IS TOTALLY RELAXED INITIALLY!!!
!           flatcur=velmatcur/vol

!ATTENTION: THIS IS WITH TAKING INTO ACCOUNT THE STRESS TENSOR GIVEN IN INPUT
           str_matrix(1,1)=strten_in(1)
           str_matrix(2,2)=strten_in(2)
           str_matrix(3,3)=strten_in(3)
           str_matrix(1,2)=strten_in(6)
           str_matrix(2,1)=strten_in(6)
           str_matrix(1,3)=strten_in(5)
           str_matrix(3,1)=strten_in(5)
           str_matrix(2,3)=strten_in(4)
           str_matrix(3,2)=strten_in(4)
           flatcur=velmatcur/vol-str_matrix
           flatcur=flatcur-pressure_md

!ATTENTION: THIS IS FOR THE ASSUMPTION THAT THE STRUCTURE IS TOTALLY RELAXED INITIALLY!!!
!Initial force on atoms
!           fposcur=0.d0

!ATTENTION: THIS IS WITH TAKING INTO ACCOUNT THE STRESS TENSOR GIVEN IN INPUT
!Convert cartesian initial input forces to reduced forces on atoms
          call invertmat(latvec,tmplat,3)
          do iat=1,nat
            fposcur(:,iat)=matmul(tmplat,fcart_in(:,iat))
          enddo

!Keep track of volume
        vol0=vol
        volprev=vol

!Initialize f0
        call cross_product(latvec(:,2),latvec(:,3),crossp); tmplat(:,1)=crossp
        call cross_product(latvec(:,3),latvec(:,1),crossp); tmplat(:,2)=crossp
        call cross_product(latvec(:,1),latvec(:,2),crossp); tmplat(:,3)=crossp
        do i=1,3
        tmplatt(:,i)=tmplat(i,:)
        enddo

!This is for the original Rahman Parrinello. You also need to change the volume-mass factor
        f0=unitmat
        f0inv=unitmat
        latmass=latmass0
        latmassinv=1.d0/latmass

!Initialize internal coordinates and lattice vectors
        poscur=xred_in           !call rxyz_cart2int(latvec,poscur,rxyz,nat)
        latcur=latvec
        vlatcur=vlat


!!Write every step
!        call rxyz_int2cart(latcur,poscur,rxyz,nat)
!        call wtpos_inter(nat,rxyz,latcur,0)

!MHM: initialize variable to track the minima/maxima***********
    write(*,*) '# MINHOP start MD'
    nummax=0
    nummin=0
    enmin1=0.d0
    en0000=0.d0
    econs_max=-1.d100
    econs_min=1.d100
!**************************************************************

!Compute initial kinetic energies
  rkin=0.d0
  do iat=1,nat
     vposcurtmp=matmul(latcur,vposcur(:,iat))
     rkin=rkin+amass(iat)*(vposcurtmp(1)**2+vposcurtmp(2)**2+vposcurtmp(3)**2)
  enddo
  ekinatom=0.5d0*rkin
  rkin=0.d0
  do i=1,3
     rkin=rkin+(vlatcur(1,i)**2+vlatcur(2,i)**2+vlatcur(3,i)**2)*latmass
  enddo
  ekinlat=0.5d0*rkin



!MHM: Write output to file in every step***********************************
!INITIAL STEP, STILL THE SAME STRUCTURE AS INPUT
       write(*,*) "Pressure, Energy",pressure,etot_in
       call get_enthalpy(latvec_in,etot_in,pressure,enthalpy)
       ent_pos_0=enthalpy
       en0000=enthalpy-ent_pos_0
!       write(67,'(a,i5,1x,4(1x,1pe17.10),3(1x,i2))') ' MD: it,enthalpy,pV,ekinat,ekinlat,nmax,nmin,mdmin ',&
!             &itime,enthalpy,pressure*vol,ekinatom,ekinlat,nummax,nummin,parres%mdmin
       write(*,'(a,i5,1x,4(1x,1pe17.10),3(1x,i2))') ' # MD: it,enthalpy,pV,ekinat,ekinlat,nmax,nmin,parres%mdmin ',&
             &itime,enthalpy,pressure*vol,ekinatom,ekinlat,nummax,nummin,parres%mdmin
       write(fn4,'(i4.4)') itime
       filename=trim(folder)//"posmd."//fn4//".ascii"
       units=units
       write(*,*) "# Writing the positions in MD:",filename
       call write_atomic_file_ascii(parini,filename,nat,units,poscur,latvec_in,fcart_in,strten_in,&
            &char_type(1:ntypat),ntypat,typat,fixat,fixlat,etot_in,pressure,enthalpy,en0000)
!*********************************************************************

        do itime=1,parini%nmd_dynamics

!          if(itime.ne.1) e_rxyz=enthalpy  !e_rxyz=e_rxyz+pressure_md(1,1)*vol   
!Check the torque on the cell for rotation
          call torque_cell(latcur,vlatcur,torquenrm)
!First calculate the steps for the internal coordinate, according to simple verlet algorithm
!For this we first need to calculate the matrix g^-1*dg/dt=ginv*gdot=gtot
          do i=1,3
            lattrans(:,i)=latcur(i,:)
            latdottrans(:,i)=vlatcur(i,:)
          enddo
          gdot=matmul(latdottrans,latcur)+matmul(lattrans,vlatcur)
          g=matmul(lattrans,latcur)
          call invertmat(g,ginv,3)
          gtot=matmul(ginv,gdot)

!Now perform the step
          velcm=0.d0
          rkin=0.d0
          do iat=1,nat
             !pospred(:,iat)=poscur(:,iat) + dt*vposcur(:,iat) + atmassinv(Kinds(iat))*dt*dt*fposcur(:,iat)&
             !               & - atmassinv(Kinds(iat))*dt*dt*matmul(gtot,vposcur(:,iat))
             pospred(:,iat)=poscur(:,iat) + dt*vposcur(:,iat) + 1.d0/amass(iat)*dt*dt*fposcur(:,iat)&
                            & - dt*dt*matmul(gtot,vposcur(:,iat))
!                            & - 1.d0/amass(iat)*dt*dt*matmul(gtot,vposcur(:,iat))
             vposcurtmp=matmul(latcur,vposcur(:,iat))
             rkin=rkin+amass(iat)*(vposcurtmp(1)**2+vposcurtmp(2)**2+vposcurtmp(3)**2)
             vpospred(:,iat)=(pospred(:,iat)-poscur(:,iat))/dt   !Update dummy variable vposcur for next iteration
             velcm=velcm+vposcurtmp*amass(iat)
          enddo
          ekinatom=0.5d0*rkin
          if(parini%verb.gt.0) write(*,'(a,es15.7)')"Velocity of CM: ", sqrt(velcm(1)**2+velcm(2)**2+velcm(3)**2)

!Now calculate steps for the lattice according to the velocity verlet algorithm
!For simplicity, just use the stress tensor at the moment for the force on the lattice
          call cross_product(latcur(:,2),latcur(:,3),crossp); sigma(:,1)=crossp
          call cross_product(latcur(:,3),latcur(:,1),crossp); sigma(:,2)=crossp
          call cross_product(latcur(:,1),latcur(:,2),crossp); sigma(:,3)=crossp

!Compute the F
          tcur=matmul(sigma,f0inv)
 
          if(options==1) then
!For velocity verlet
              latpred(:,:)=latcur(:,:) + dt*vlatcur(:,:) + (.5d0*dt*dt*latmassinv)*matmul(flatcur(:,:),tcur(:,:))
          elseif(options==2) then
!Instead if normal verlet is used:
              latpred(:,:)=latcur(:,:) + dt*vlatcur + dt*dt*latmassinv*matmul(flatcur(:,:),tcur(:,:))
              vlatpred(:,:)=(latpred-latcur)/dt
          else
              stop "Wrong option for the Cell part in MD"
          endif

!Kinetic energy according to Parrinello Rahman
          rkin=0.d0
          do i=1,3
             rkin=rkin+(vlatcur(1,i)**2+vlatcur(2,i)**2+vlatcur(3,i)**2)*latmass
          enddo
          ekinlat=0.5d0*rkin
          rkin=ekinlat+ekinatom
          write(*,'(a,3(1x,es15.7))') " # Torquenrm, Ekin, Enthalpy: ",torquenrm, rkin,enthalpy
!Update counter
          enmin2=enmin1
          enmin1=en0000

!Update velocity part of stress tensor
        velmatpred=0.d0
        do iat=1,nat
           vpospredtmp=matmul(latcur,vpospred(:,iat))
           do i=1,3
              do j=1,3
              velmatpred(i,j)=velmatpred(i,j)+vpospredtmp(i)*vpospredtmp(j)
              enddo
           enddo
        enddo

!Here we perform the force call
!****************************************************************************************************************        
!****************************************************************************************************************        
!Get the current acell_in and rprim_in, and also the atomic positions
       xred_in=pospred 
       vel_in=0.d0
       latvec_in=latpred
       if(itime==1)  then
           getwfk=.false.
       elseif(itime.ne.1.and.usewf_md) then
           getwfk=.true.
       else
           getwfk=.false.
       endif
       call get_energyandforces_single(parini,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,getwfk)
!Single point calculation finished. Now returning to MD part. All output variables are now in *_in
!****************************************************************************************************************        

!Convert cartesian forces to reduced forces on atoms
        call invertmat(latpred,tmplat,3)
        do iat=1,nat
          fpospred(:,iat)=matmul(tmplat,fcart_in(:,iat))
        enddo

!Transform the forces on the lattice vectors into the stress tensore according to the paper Tomas Bucko and Jurg Hafner
!This procedure is no longer necessary, since we have the correct stress tensor now!
!        do i=1,3
!           tmplat(:,i)=latpred(i,:)
!        enddo
!        flatpred=matmul(flatpred,tmplat)
           str_matrix(1,1)=strten_in(1) 
           str_matrix(2,2)=strten_in(2) 
           str_matrix(3,3)=strten_in(3) 
           str_matrix(1,2)=strten_in(6) 
           str_matrix(2,1)=strten_in(6) 
           str_matrix(1,3)=strten_in(5) 
           str_matrix(3,1)=strten_in(5) 
           str_matrix(2,3)=strten_in(4) 
           str_matrix(3,2)=strten_in(4) 
        a=latpred
        vol= a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)+&
             a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)

        flatpred=velmatpred/vol-str_matrix
!Here the pressure is applied
        flatpred=flatpred-pressure_md
        volprev=vol

        call cross_product(latpred(:,2),latpred(:,3),crossp); sigma(:,1)=crossp
        call cross_product(latpred(:,3),latpred(:,1),crossp); sigma(:,2)=crossp
        call cross_product(latpred(:,1),latpred(:,2),crossp); sigma(:,3)=crossp
        tpred=matmul(sigma,f0inv)

!For velocity verlet of cell
        if (options==1)  vlatpred(:,:)=vlatcur(:,:)+0.5d0*dt*latmassinv*(matmul(flatpred,tpred)+matmul(flatcur(:,:),tcur(:,:)))
!Compute the "predicted" kinetic energies:
  rkin=0.d0
  do iat=1,nat
     vposcurtmp=matmul(latpred,vpospred(:,iat))
     rkin=rkin+amass(iat)*(vposcurtmp(1)**2+vposcurtmp(2)**2+vposcurtmp(3)**2)
  enddo
  ekinatom=0.5d0*rkin
  rkin=0.d0
  do i=1,3
     rkin=rkin+(vlatpred(1,i)**2+vlatpred(2,i)**2+vlatpred(3,i)**2)*latmass
  enddo
  ekinlat=0.5d0*rkin


!MHM: Write output to file in every step***********************************
       counter=real(itime,8)
       call get_enthalpy(latvec_in,etot_in,pressure,enthalpy)
       en0000=enthalpy-ent_pos_0
       if (itime >= 3 .and. enmin1 > enmin2 .and. enmin1 > en0000)  nummax=nummax+1
       if (itime >= 3 .and. enmin1 < enmin2 .and. enmin1 < en0000)  nummin=nummin+1
      ! write(67,'(a,i5,1x,4(1x,1pe17.10),3(1x,i2))') ' MD: it,enthalpy,pV,ekinat,ekinlat,nmax,nmin,mdmin ',&
      !       &itime,enthalpy,pressure*vol,ekinatom,ekinlat,nummax,nummin,parres%mdmin
       write(*,'(a,i5,1x,4(1x,1pe17.10),3(1x,i2))') ' # MD: it,enthalpy,pV,ekinat,ekinlat,nmax,nmin,mdmin ',&
             &itime,enthalpy,pressure*vol,ekinatom,ekinlat,nummax,nummin,parres%mdmin
       write(fn4,'(i4.4)') itime
       filename=trim(folder)//"posmd."//fn4//".ascii"
       units=units
       write(*,*) "# Writing the positions in MD: ",filename
       call write_atomic_file_ascii(parini,filename,nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &char_type(1:ntypat),ntypat,typat,fixat,fixlat,etot_in,pressure,enthalpy,en0000)
!*********************************************************************

!        econs_max=max(econs_max,rkin+e_rxyz)
!        econs_min=min(econs_min,rkin+e_rxyz)
        econs_max=max(econs_max,rkin+enthalpy)
        econs_min=min(econs_min,rkin+enthalpy)
!        write(6,*) istep,e_rxyz-e_pos,nummax,nummin
       if (nummin.ge.parres%mdmin) then
          if (nummax.ne.nummin) &
               write(*,*) '# WARNING: nummin,nummax',nummin,nummax

!          write(67,*) " MD finished: exiting!"
          write(*,*) " MD finished: exiting!"
          exit
       endif


!Update the variables for next iteration
        velmatcur=velmatpred
        flatcur=flatpred
        fposcur=fpospred
        poscur=pospred
        vlatcur=vlatpred
        vposcur=vpospred
        latcur=latpred
!        e_rxyz=e_rxyz_pred
     enddo 
!MD stopped, now do relaxation
     write(*,'(a,i5,es15.7,es15.7)') ' # EXIT MD ',itime,enthalpy,etot_in
stop
end subroutine 

!**********************************************************************************************

subroutine GEOPT_FIRE_MHM(parini,latvec_in,xred_in,fcart_in,strten_in,vel_in,vel_lat_in,vvol_in,etot_in,iprec,counter,folder)
 use mod_interface
 use global, only: target_pressure_habohr,target_pressure_gpa,nat,ntypat,znucl,amu,amutmp,typat
 use global, only: char_type
 use global, only: units,usewf_geopt,max_kpt,fixat,fixlat,correctalg,ka1,kb1,kc1,confine
 use defs_basis
 use mod_fire
 use interface_code
 use modsocket, only: sock_extra_string
 use mod_parini, only: typ_parini
implicit none
 type(typ_parini), intent(in):: parini
 real(8) :: latvec_in(3,3),xred_in(3,nat),vel_in(3,nat),fcart_in(3,nat),etot_in,strten_in(6),vel_lat_in(3,3),vvol_in
!*********************************************************************
!Variables for my MD part
! real(8), parameter :: Ha_eV=27.21138386_dp ! 1 Hartree, in eV
! real(8), parameter :: kb_HaK=8.617343d-5/Ha_eV ! Boltzmann constant in Ha/K
! real(8), parameter :: amu_emass=1.660538782d-27/9.10938215d-31 ! 1 atomic mass unit, in electronic mass
 real(8):: amass(nat)
 real(8):: pressure_ener
 real(8):: pressure_md(3,3)
 real(8):: pressure
 real(8):: vol
 real(8):: vol0
 real(8):: volprev
 real(8),dimension(3,nat):: xcart
 real(8),dimension(3,nat):: fposcur
 real(8),dimension(3,nat):: fpospred
 real(8),dimension(3,nat):: vpospred
 real(8),dimension(3,nat):: poscur
 real(8),dimension(3,nat):: vxyz
 real(8),dimension(3,nat):: vposcur
 real(8),dimension(3,nat):: pospred
 real(8),dimension(3,nat):: accposcur
 real(8),dimension(3,nat):: accpospred
 real(8),dimension(3,nat):: accposprev
 real(8),dimension(3,nat):: dxred
 real(8),dimension(3,3):: dlatvec
 real(8),dimension(3,3):: latvec
 real(8),dimension(3,3):: latvec0
 real(8),dimension(3,3):: latinv
 real(8),dimension(3,3):: unitmat
 real(8),dimension(3,3):: elmatcur
 real(8),dimension(3,3):: flatcur
 real(8),dimension(3,3):: flatpred
 real(8),dimension(3,3):: velmatcur
 real(8),dimension(3,3):: tmplat
 real(8),dimension(3,3):: tmplatt
 real(8),dimension(3,3):: tpred
 real(8),dimension(3,3):: f0
 real(8),dimension(3,3):: f0inv
 real(8),dimension(3,3):: latcur
 real(8),dimension(3,3):: vlatcur
 real(8),dimension(3,3):: vlat
 real(8),dimension(3,3):: lattrans
 real(8),dimension(3,3):: latdottrans
 real(8),dimension(3,3):: a 
 real(8),dimension(3,3):: g
 real(8),dimension(3,3):: gdot
 real(8),dimension(3,3):: ginv
 real(8),dimension(3,3):: gtot
 real(8),dimension(3,3):: sigma
 real(8),dimension(3,3):: sigmatrans
 real(8),dimension(3,3):: latpred
 real(8),dimension(3,3):: tcur
 real(8),dimension(3,3):: vlatpred
 real(8),dimension(3,3):: velmatpred
 real(8),dimension(3,3):: str_matrix
 real(8),dimension(3,3):: acclatcur
 real(8),dimension(3,3):: acclatpred
 real(8),dimension(3,3):: acclatprev
 real(8):: accvolcur
 real(8):: accvolpred
 real(8):: accvolprev
 real(8):: fvolcur
 real(8):: fvolpred
 real(8):: volcur
 real(8):: volpred
 real(8):: vvolcur
 real(8):: vvolpred
 real(8):: volpred_1_3
 real(8):: vol_1_3_inv
 real(8):: vol_1_3
 real(8):: vposcurtmp(3)
 real(8):: vpospredtmp(3)
 real(8):: crossp(3)
 real(8):: velcm(3)
 real(8):: latmass
 real(8):: latmassinv
 real(8):: latmass0
 real(8):: ekinatom
 real(8):: ekinlat
 real(8):: ekinlat_prev
 real(8):: ekinatom_prev
 real(8):: rkin
 real(8):: enthalpy
 real(8):: enthalpy_min
 real(8):: enmin1
 real(8):: enmin2
 real(8):: ent_pos_0
 real(8):: en0000
 real(8):: e_rxyz
 real(8):: econs_max
 real(8):: econs_min
 real(8):: torquenrm 
 real(8):: dt 
 real(8):: ecut_tmp
 real(8):: toldff_tmp
 real(8):: dstr(6)
 real(8):: strtarget(6)
 real(8):: counter 
 integer:: i
 integer:: j
 integer:: iat
 integer:: itime
 integer:: nummax
 integer:: nummin
 integer:: options
 integer:: iexit
 integer:: iprec
 integer:: md_type
 logical:: getwfk
 character(40)::filename,folder
 character(4) ::fn4

! FIRE VARIABLES
 real(8):: alpha,P,P_at,P_lat,fmax,fmax_at,fmax_lat,fall(3,nat+3),fallnorm,vall(3,nat+3),vallnorm
 integer:: nstep,istr
 logical:: multiprec
 integer:: iprec_cur
 real(8):: tolmxf_switch,alpha_lat,f_lat
 logical:: cellfix_done
 real(8):: cellfix_switch
!Latvec_variable
 integer:: latvec_io
latvec_io=0
!multiprec is hardcoded and, if true, starts a geopt with iprec==2, and then switches 
!to iprec==1 when the fmax==tolmxf_switch. The switch only occurs once
 max_kpt=.false.
 multiprec=.true.
 tolmxf_switch=1.d-3
 enthalpy_min=1.d10
 cellfix_switch=1.d-3
 cellfix_done=.false.


!Lattice alpha
alpha_lat=20.d0

!Assign masses to each atom (for FIRE, all atoms have the same mass)
 do iat=1,nat
!   amass(iat)=amu_emass*1.d0
   amass(iat)=amu_emass*amu(typat(iat))
   if(parini%verb.gt.0) write(*,'(a,i5,2(1x,es15.7))') " # FIRE: iat, AMU, EM: ", iat, amass(iat)/amu_emass,amass(iat)
 enddo

!Set FIRE parameters
alpha=alphastart
!if(fixlat(7)) alpha=1.d0
nstep=0
latmass0=latmass_rel_fire*sum(amass(:))/real(nat,8)
dt=parini%paropt_geopt%dt_start
vel_in=0.d0
vel_lat_in=0.d0
P=0.d0
P_at=0.d0
P_lat=0.d0

!MD-type: 1 for PR and 2 for Cleveland and 3 for Wentzcovitch 4 Andersen
  if(.not.fixlat(7)) then 
      md_type=1
  else
      md_type=4
      latmass0=0.005d0*latmass_rel_fire*sum(amass(:))/real(nat,8)
  endif

   write(*,'(a,i5,2(1x,es15.7))') " # FIRE: latmass, AMU, EM: ", iat, latmass0/amu_emass,latmass0

!******************************************************************
!NEW VERISON OF MD: Directly implemented, simplest Parrinello-Rahman MD
!Can be accessed if the ionmov of idtset2 is equal to -13

!Here we split the routine in the Abinit native part and my new implentation
write(*,'(a)') ' # Entering FIRE based on standalone PR MD '
!Set temporary variables, initially
  vxyz(:,:)=vel_in(:,:)
  pressure=target_pressure_habohr
  unitmat=0.d0
  do i=1,3
    unitmat(i,i)=1.d0
  enddo
  vlat=vel_lat_in  !The initial cell velocity
  itime=0

!Set options=1 for Velocity Verlet of cell dynamics
!Set options=2 for Normal Verlet of cell dynamics
!Set options=3 for Beeman integration scheme, corrector-predictor
  options=3
  if(options.lt.1.or.options.gt.3) stop "Wrong integrator option"


  if(md_type.lt.1.or.md_type.gt.4) stop "Wrong algo option"

if(parini%verb.gt.0) write(*,'(a,i3,a,i3)') " # GEOPT Algorithm: ",md_type, ", MD Integrator: ",options

!Now we run my implementation of PR MD
pressure_md=0.d0
pressure_md(1,1)=1.d0;pressure_md(2,2)=1.d0;pressure_md(3,3)=1.d0
pressure_ener=0.d0;pressure_md=pressure_md*pressure  !Here the pressure is not passed to the energyandforces, so we move on the ENERGY surface

!Transform the initial velocities given in cartesian coordinates into velocities of the intenal coordinates
       latvec=latvec_in
       call getvol(latvec,vol)
if(md_type==4) then
       vol_1_3=vol**(1.d0/3.d0)
       vol_1_3_inv=1.d0/vol_1_3
       do iat=1,nat
          vposcur(:,iat)=vxyz(:,iat)*vol_1_3_inv
       enddo
else
       call invertmat(latvec,latinv,3)
       do iat=1,nat
          vposcur(:,iat)=matmul(latinv,vxyz(:,iat))
       enddo
endif

!Compute f0
        call cross_product(latvec(:,2),latvec(:,3),crossp); sigma(:,1)=crossp
        call cross_product(latvec(:,3),latvec(:,1),crossp); sigma(:,2)=crossp
        call cross_product(latvec(:,1),latvec(:,2),crossp); sigma(:,3)=crossp
        sigmatrans(:,1)=sigma(1,:)
        sigmatrans(:,2)=sigma(2,:)
        sigmatrans(:,3)=sigma(3,:)
        f0=matmul(sigmatrans,sigma)
        call invertmat(f0,f0inv,3)


!Here we perform the initial force call
!****************************************************************************************************************        
!****************************************************************************************************************        
       vel_in=0.d0
       getwfk=.false.
       write(fn4,'(i4.4)') 0
       sock_extra_string="FIRE"//trim(fn4)
       call get_energyandforces_single(parini,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,getwfk)
!Single point calculation finished. Now returning to MD part. All output variables are now in *_in
!****************************************************************************************************************        
!FIRE: check for convergence
call convcheck(nat,latvec_in,fcart_in,strten_in,target_pressure_habohr,parini%paropt_geopt%strfact,fmax,fmax_at,fmax_lat,parini%paropt_geopt%fmaxtol,iexit)
!!!!!Compute maximal component of forces, EXCLUDING any fixed components
!!! fmax=0.0d0
!!! do iat=1,nat
!!!   do i=1,3
!!!!     if (dtsets(1)%iatfix(i,iat) /= 1) then
!!!       if( abs(fcart_in(i,iat)) >= fmax ) fmax=abs(fcart_in(i,iat))
!!!!     end if
!!!   end do
!!! end do
!!! strtarget=0.d0
!!! strtarget(1:3)=-target_pressure_habohr
!!! dstr(:)=strten_in(:)-strtarget(:)
!!!!Eventually take into account the stress
!!! do istr=1,6
!!!     if(abs(dstr(istr))*strfact >= fmax ) fmax=abs(dstr(istr))*strfact
!!! end do
!!! iexit=0
!!! if(fmax.lt.tolmxf) iexit=1
!!!
!Initial iprec after running the first force call
 if(multiprec) iprec=2

 fall=0.d0
 call fpos_flat(pressure_md,fall(:,1:nat),fall(:,nat+1:nat+3),strten_in,fcart_in,latvec_in,md_type)
 if(md_type==4) f_lat=fall(1,nat+1)

!Keep track of volume
        vol0=vol
        latvec0=latvec/vol**(1.d0/3.d0)
        volprev=vol

        if(md_type==1) then
!!This is for the original Rahman Parrinello. You also need to change the volume-mass factor
        latmass=latmass0
        latmassinv=1.d0/latmass
        elseif(md_type==2) then
!This is for cleveland
        latmass=latmass0/vol**(4.d0/3.d0)
        latmassinv=1.d0/latmass
        elseif(md_type==3) then
!This is for Wentzcovitch
        latmass=latmass0
        latmassinv=1.d0/latmass
        elseif(md_type==4) then
!This is for Andersen
        latmass=latmass0
        latmassinv=1.d0/latmass
        endif

!Initialize internal coordinates and lattice vectors
if(md_type==4) then
        call rxyz_int2cart(latvec,xred_in,xcart,nat)
        poscur=xcart*vol_1_3_inv           !call rxyz_cart2int(latvec,poscur,rxyz,nat)
        volcur=vol
        vvolcur=vvol_in
        latcur=latvec
else
        poscur=xred_in           !call rxyz_cart2int(latvec,poscur,rxyz,nat)
        latcur=latvec
        vlatcur=vlat
endif
!!Write every step
!        call rxyz_int2cart(latcur,poscur,rxyz,nat)
!        call wtpos_inter(nat,rxyz,latcur,0)

!MHM: initialize variable to track the minima/maxima***********
    write(*,*) '# MINHOP start FIRE'
    nummax=0
    nummin=0
    enmin1=0.d0
    en0000=0.d0
    econs_max=-1.d100
    econs_min=1.d100
!**************************************************************

!MHM: Write output to file in every step***********************************
!INITIAL STEP, STILL THE SAME STRUCTURE AS INPUT
       write(*,*) "Pressure, Energy",pressure,etot_in
       call get_enthalpy(latvec_in,etot_in,pressure,enthalpy)
       ent_pos_0=enthalpy
       en0000=enthalpy-ent_pos_0
if(parini%verb.gt.0) then
       write(fn4,'(i4.4)') itime
       filename=trim(folder)//"posgeopt."//fn4//".ascii"
       units=units
       write(*,*) "# Writing the positions in FIRE:",filename
       call write_atomic_file_ascii(parini,filename,nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &char_type(1:ntypat),ntypat,typat,fixat,fixlat,etot_in,pressure,enthalpy,en0000)
       if(parini%verb.ge.3) then
       filename=trim(folder)//"posgeopt."//fn4//".vasp"
       call write_atomic_file_poscar(parini,filename,nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &char_type(1:ntypat),ntypat,typat,fixat,fixlat,etot_in,pressure,enthalpy,en0000)
       endif
endif
!*********************************************************************

!Exit if the forces are already converged at the beginning
        if(iexit==1) then
          write(*,'(a,i4,2(1x,es25.15))') " # FIRE converged before entering iterations", itime,enthalpy,fmax
          max_kpt=.false.
          return 
        endif

call acceleration_fire(pressure_md,accposcur,acclatcur,accvolcur,vposcur,&
     &vlatcur,vvolcur,strten_in,fcart_in,latcur,amass,latmass,f0inv,md_type,nat)
        accposprev=accposcur
        acclatprev=acclatcur
        accvolprev=accvolcur

!Eliminate the unnecessary atmic and cell components
call elim_fixed_at(nat,accposcur);call elim_fixed_at(nat,accposprev)
if(md_type.ne.4) call elim_fixed_lat(latcur,acclatcur)
if(md_type.ne.4) call elim_fixed_lat(latcur,acclatprev)
call elim_fixed_at(nat,vposcur)
if(md_type.ne.4) call elim_fixed_lat(latcur,vlatcur)


!FIRE cycles
        do itime=1,parini%paropt_geopt%nit
!          if(itime.ne.1) e_rxyz=enthalpy  !e_rxyz=e_rxyz+pressure_md(1,1)*vol   
!Check the torque on the cell for rotation
          call torque_cell(latcur,vlatcur,torquenrm)
!Now perform the step
          velcm=0.d0
          rkin=0.d0
          if(options==1) then
!For velocity verlet
!             latpred(:,:)=latcur(:,:) + dt*vlatcur(:,:) + (.5d0*dt*dt*acclatcur)
             dlatvec(:,:)=dt*vlatcur(:,:) + (.5d0*dt*dt*acclatcur)
          do iat=1,nat
              dxred(:,iat)=dt*vposcur(:,iat) + dt*dt*accposcur(:,iat)
!             pospred(:,iat)=poscur(:,iat) + dt*vposcur(:,iat) + dt*dt*accposcur(:,iat)
          enddo
          call propagate(nat,poscur,latcur,dxred,dlatvec,pospred,latpred)
          do iat=1,nat
             vpospred(:,iat)=(pospred(:,iat)-poscur(:,iat))/dt   !Update dummy variable vposcur for next iteration
          enddo
            if(md_type==4) then
!               volpred=volcur+dt*vvolcur+0.5d0*dt*dt*accvolcur
               volpred=volcur+alpha_lat*fall(1,nat+1)
            endif
          elseif(options==2) then
!Instead if normal verlet is used:
!              latpred(:,:)=latcur(:,:) + dt*vlatcur + dt*dt*acclatcur
              dlatvec(:,:)=dt*vlatcur + dt*dt*acclatcur
          do iat=1,nat
              dxred(:,iat)=dt*vposcur(:,iat) + dt*dt*accposcur(:,iat)
!             pospred(:,iat)=poscur(:,iat) + dt*vposcur(:,iat) + dt*dt*accposcur(:,iat)
          enddo
          call propagate(nat,poscur,latcur,dxred,dlatvec,pospred,latpred)
          vlatpred(:,:)=(latpred-latcur)/dt
          do iat=1,nat
              vpospred(:,iat)=(pospred(:,iat)-poscur(:,iat))/dt   !Update dummy variable vposcur for next iteration
          enddo
            if(md_type==4) then
!               volpred=volcur+dt*vvolcur+dt*dt*accvolcur
!               vvolpred=(volpred-volcur)/dt
               volpred=volcur+alpha_lat*fall(1,nat+1)
            endif
          elseif(options==3) then
!Predictor part of Beeman for cell
!              latpred(:,:)=latcur(:,:) + dt*vlatcur + dt*dt/6.d0*(4.d0*acclatcur(:,:)-acclatprev(:,:))
              dlatvec(:,:)=dt*vlatcur + dt*dt/6.d0*(4.d0*acclatcur(:,:)-acclatprev(:,:))
          do iat=1,nat
!Predictor part of Beeman
!             pospred(:,iat)=poscur(:,iat) + dt*vposcur(:,iat) + dt*dt/6.d0*(4.d0*accposcur(:,iat)-accposprev(:,iat))             
             dxred(:,iat)=dt*vposcur(:,iat) + dt*dt/6.d0*(4.d0*accposcur(:,iat)-accposprev(:,iat))
          enddo
          call propagate(nat,poscur,latcur,dxred,dlatvec,pospred,latpred)
          vlatpred(:,:)=vlatcur(:,:) + 0.5d0*(3.d0*acclatcur(:,:)-acclatprev(:,:))
          do iat=1,nat
             vpospred(:,iat)=vposcur(:,iat)+0.5d0*(3.d0*accposcur(:,iat)-accposprev(:,iat))
          enddo
            if(md_type==4) then
!               volpred=volcur+dt*vvolcur+dt*dt/6.d0*(4.d0*accvolcur-accvolprev)
!               vvolpred=volcur+0.5d0*(3.d0*accvolcur-accvolprev)
               volpred=volcur+alpha_lat*fall(1,nat+1)
            endif
          endif

if(md_type==4) then
  do iat=1,nat
     vposcurtmp=vposcur(:,iat)*vol_1_3
     rkin=rkin+amass(iat)*(vposcurtmp(1)**2+vposcurtmp(2)**2+vposcurtmp(3)**2)
     velcm=velcm+vposcurtmp*amass(iat)
  enddo
else
  do iat=1,nat
             vposcurtmp=matmul(latcur,vposcur(:,iat))
             rkin=rkin+amass(iat)*(vposcurtmp(1)**2+vposcurtmp(2)**2+vposcurtmp(3)**2)
             velcm=velcm+vposcurtmp*amass(iat)
  enddo
endif



!!!
!!!
!!!
!!!!Now perform the step
!!!          velcm=0.d0
!!!          rkin=0.d0
!!!          do iat=1,nat
!!!          if(options==1.or.options==2) then
!!!             pospred(:,iat)=poscur(:,iat) + dt*vposcur(:,iat) + dt*dt*accposcur(:,iat)
!!!             vpospred(:,iat)=(pospred(:,iat)-poscur(:,iat))/dt   !Update dummy variable vposcur for next iteration
!!!          elseif(options==3) then
!!!!Predictor part of Beeman
!!!             pospred(:,iat)=poscur(:,iat) + dt*vposcur(:,iat) + dt*dt/6.d0*(4.d0*accposcur(:,iat)-accposprev(:,iat))
!!!             vpospred(:,iat)=vposcur(:,iat)+0.5d0*(3.d0*accposcur(:,iat)-accposprev(:,iat))
!!!          endif
!!!             vposcurtmp=matmul(latcur,vposcur(:,iat))
!!!             rkin=rkin+amass(iat)*(vposcurtmp(1)**2+vposcurtmp(2)**2+vposcurtmp(3)**2)
!!!             velcm=velcm+vposcurtmp*amass(iat)
!!!          enddo
!!!
!!!!Now calculate steps for the lattice according to the velocity verlet algorithm
!!!          if(options==1) then
!!!!For velocity verlet
!!!              latpred(:,:)=latcur(:,:) + dt*vlatcur(:,:) + (.5d0*dt*dt*acclatcur)
!!!          elseif(options==2) then
!!!!Instead if normal verlet is used:
!!!              latpred(:,:)=latcur(:,:) + dt*vlatcur + dt*dt*acclatcur
!!!              vlatpred(:,:)=(latpred-latcur)/dt
!!!          elseif(options==3) then
!!!!Predictor part of Beeman for cell
!!!              latpred(:,:)=latcur(:,:) + dt*vlatcur + dt*dt/6.d0*(4.d0*acclatcur(:,:)-acclatprev(:,:))
!!!              vlatpred(:,:)=vlatcur(:,:) + 0.5d0*(3.d0*acclatcur(:,:)-acclatprev(:,:))
!!!          else
!!!              stop "Wrong option for the Cell part in MD"
!!!          endif

!Kinetic energy of atoms
          ekinatom=0.5d0*rkin
          if(parini%verb.gt.0) write(*,'(a,es15.7)')"Velocity of CM: ", sqrt(velcm(1)**2+velcm(2)**2+velcm(3)**2)
!Kinetic energy according to Parrinello Rahman
          rkin=0.d0
if(md_type==4) then
  rkin=vvolcur**2*latmass
  ekinlat=0.5d0*rkin
else
          do i=1,3
             rkin=rkin+(vlatcur(1,i)**2+vlatcur(2,i)**2+vlatcur(3,i)**2)*latmass
          enddo
          ekinlat=0.5d0*rkin
          rkin=ekinlat+ekinatom
endif
!          call ekin_at_lat(amass,latmass,latcur,vposcur,vlatcur,ekinatom,ekinlat,f0,md_type,nat)
          rkin=ekinlat+ekinatom
if(parini%verb.gt.0) write(*,'(a,4(1x,es15.7))') " # Torquenrm, Ekin, Enthalpy: ",torquenrm, ekinlat,ekinatom,enthalpy
!Update counter
          enmin2=enmin1
          enmin1=en0000

!Here we perform the force call
!****************************************************************************************************************        
!****************************************************************************************************************        
!Get the current acell_in and rprim_in, and also the atomic positions
if(md_type==4) then
       volpred_1_3=volpred**(1.d0/3.d0)
       latvec_in=latvec0*volpred_1_3
       latpred=latvec_in
       xcart=pospred*volpred_1_3
       call rxyz_cart2int(latvec_in,xred_in,xcart,nat)
else
       xred_in=pospred
       vel_in=0.d0
       latvec_in=latpred
endif
!       getwfk=.true.
       if(usewf_geopt) then
           getwfk=.true.
       else
           getwfk=.false.
       endif
       if(iprec==2.and.multiprec.and.fmax.lt.tolmxf_switch) then
           getwfk=.false.
           iprec=1
           enthalpy_min=1.d10
       endif
       write(fn4,'(i4.4)') itime
       sock_extra_string="FIRE"//trim(fn4)
       call get_energyandforces_single(parini,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,getwfk)
!Single point calculation finished. Now returning to MD part. All output variables are now in *_in
!****************************************************************************************************************        

!FIRE: check for convergence
call convcheck(nat,latvec_in,fcart_in,strten_in,target_pressure_habohr,parini%paropt_geopt%strfact,fmax,fmax_at,fmax_lat,parini%paropt_geopt%fmaxtol,iexit)
!!!!Compute maximal component of forces, EXCLUDING any fixed components
!!! fmax=0.0d0
!!! do iat=1,nat
!!!   do i=1,3
!!!!     if (dtsets(1)%iatfix(i,iat) /= 1) then
!!!       if( abs(fcart_in(i,iat)) >= fmax ) fmax=abs(fcart_in(i,iat))
!!!!     end if
!!!   end do
!!! end do
!!! strtarget=0.d0
!!! strtarget(1:3)=-target_pressure_habohr
!!! dstr(:)=strten_in(:)-strtarget(:)
!!!!Eventually take into account the stress
!!! do istr=1,6
!!!     if(abs(dstr(istr))*strfact >= fmax ) fmax=abs(dstr(istr))*strfact
!!! end do
!!! iexit=0
!!! if(fmax.lt.tolmxf) iexit=1


!For velocity verlet of cell
        if(options==1) then
           call acceleration_fire(pressure_md,accpospred,acclatpred,accvolpred,vpospred,&
                &vlatcur,vvolcur,strten_in,fcart_in,latpred,amass,latmass,f0inv,md_type,nat)
           vlatpred(:,:)=vlatcur(:,:)+0.5d0*dt*(acclatpred+acclatcur)
           if(md_type==4) vvolpred=0.d0! vvolpred=vvolcur+0.5d0*dt*(accvolpred+accvolcur)
        elseif(options==3) then
!Corrector part of Beeman. Note that a fixed number of iterations are used!!!
           call ekin_at_lat_andersen(amass,latmass,latpred,vpospred,vlatpred,vvolpred,ekinatom_prev,ekinlat_prev,f0,md_type,nat)
           do i=1,5
             call acceleration_fire(pressure_md,accpospred,acclatpred,accvolpred,vpospred,&
                  &vlatpred,vvolpred,strten_in,fcart_in,latpred,amass,latmass,f0inv,md_type,nat)
             do iat=1,nat
             vpospred(:,iat)=vposcur(:,iat)+dt/6.d0*(2.d0*accpospred(:,iat)+5.d0*accposcur(:,iat)-accposprev(:,iat))
             enddo
             vlatpred=vlatcur+dt/6.d0*(2.d0*acclatpred+5.d0*acclatcur-acclatprev)
             if(md_type==4) vvolpred=0.d0! vvolpred=vvolcur+dt/6.d0*(2.d0*accvolpred+5.d0*accvolcur-accvolprev)
             call ekin_at_lat_andersen(amass,latmass,latpred,vpospred,vlatpred,vvolpred,ekinatom,ekinlat,f0,md_type,nat)
if(parini%verb.gt.1)write(* ,'(a,i5,1x,es15.5,es15.5)') " # Beeman: Corrector relative convergence: it, atoms, lattice: ",&
                  & i,abs(ekinatom-ekinatom_prev)/ekinatom,abs(ekinlat-ekinlat_prev)/max(ekinlat,1.d-100)
!             write(67,'(a,i5,1x,es15.5,es15.5)') " # Beeman: Corrector relative convergence: it, atoms, lattice: ",&
!                  & i,abs(ekinatom-ekinatom_prev)/ekinatom,abs(ekinlat-ekinlat_prev)/max(ekinlat,1.d-100)
             ekinlat_prev=ekinlat
             ekinatom_prev=ekinatom
           enddo
        endif
!call acceleration(pressure_md,accpospred,acclatpred,accvolpred,vpospred,vlatpred,vvolpred,strten_in,fcart_in,latcur,amass,latmass,f0inv,md_type,nat)!BUG?
call acceleration_fire(pressure_md,accpospred,acclatpred,accvolpred,vpospred,&
     &vlatpred,vvolpred,strten_in,fcart_in,latpred,amass,latmass,f0inv,md_type,nat)

!Eliminate the unnecessary atmic and cell components
call elim_fixed_at(nat,accpospred);call elim_fixed_at(nat,vpospred)
if(md_type.ne.4) call elim_fixed_lat(latcur,acclatpred)
if(md_type.ne.4) call elim_fixed_lat(latcur,vlatpred)
!Compute the "predicted" kinetic energies:
call ekin_at_lat_andersen(amass,latmass,latpred,vpospred,vlatpred,vvolpred,ekinatom,ekinlat,f0,md_type,nat)


!MHM: Write output to file in every step***********************************
       counter=real(itime,8)
       call get_enthalpy(latvec_in,etot_in,pressure,enthalpy)
       call getvol(latvec_in,vol)
       en0000=enthalpy-ent_pos_0
       if (itime >= 3 .and. enmin1 > enmin2 .and. enmin1 > en0000)  nummax=nummax+1
       if (itime >= 3 .and. enmin1 < enmin2 .and. enmin1 < en0000)  nummin=nummin+1
!       if (mpi_enreg%me == 0) write(67,'(a,i5,1x,1pe17.10)') 'FIRE ' ,itime,enthalpy
!       if (mpi_enreg%me == 0) write(*,'(a,i5,1x,1pe17.10)') ' #FIRE ',itime,enthalpy
if(parini%verb.gt.0) then
       write(fn4,'(i4.4)') itime
       filename=trim(folder)//"posgeopt."//fn4//".ascii"
       units=units
       write(*,*) "# Writing the positions in FIRE: ",filename
       call write_atomic_file_ascii(parini,filename,nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &char_type(1:ntypat),ntypat,typat,fixat,fixlat,etot_in,pressure,enthalpy,en0000)
       if(parini%verb.ge.3) then
       filename=trim(folder)//"posgeopt."//fn4//".vasp"
       call write_atomic_file_poscar(parini,filename,nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &char_type(1:ntypat),ntypat,typat,fixat,fixlat,etot_in,pressure,enthalpy,en0000)
       endif

       write(*,'(a,i4,6(1x,es17.8),3(1x,es9.2),1x,i4,1x,i4)') " # GEOPT FIRE    ",&
              &itime,enthalpy, fmax, fmax_at,fmax_lat,rkin,P,P_at,P_lat,dt,nstep,iprec
endif
        if(iexit==1) then
          write(*,'(a,i4,2(1x,es25.15))') " # FIRE converged", itime,enthalpy,fmax
          exit
        endif
!*********************************************************************

!        econs_max=max(econs_max,rkin+e_rxyz)
!        econs_min=min(econs_min,rkin+e_rxyz)
        econs_max=max(econs_max,rkin+enthalpy)
        econs_min=min(econs_min,rkin+enthalpy)
!        write(6,*) istep,e_rxyz-e_pos,nummax,nummin

!Update the variables for next iteration
        accposprev=accposcur
        accposcur=accpospred
        acclatprev=acclatcur
        acclatcur=acclatpred
        velmatcur=velmatpred
        flatcur=flatpred
        fposcur=fpospred
        poscur=pospred
        latcur=latpred

        accvolprev=accvolcur
        accvolcur=accvolpred
if(md_type==4)   volcur=volpred

!        vvolcur=vvolpred
!        vlatcur=vlatpred
!        vposcur=vpospred
!FIRE velocity update ****************
        call fpos_flat(pressure_md,fall(:,1:nat),fall(:,nat+1:nat+3),strten_in,fcart_in,latpred,md_type)
!        do iat=1,nat
!          fall(:,iat)=accpospred(:,iat)*amass(iat)
!        enddo
!        fall(:,nat+1:nat+3)=acclatpred*latmass 
        vall=0.d0
        vall(:,1:nat)=vpospred(:,1:nat)         

!        if(.not.fixlat(7)) then
        vall(:,nat+1:nat+3)=vlatpred(:,:)     
        if(md_type==4) then
            vall(:,nat+1:nat+3)=0.d0
            vall(1,nat+1)=vvolpred
        endif
!           else
!           vall(:,nat+1:nat+3)=0.d0
!        endif    
        P=0.d0
        do i=1,nat+3
         P=P+fall(1,i)*vall(1,i)+fall(2,i)*vall(2,i)+fall(3,i)*vall(3,i) 
         if(i==nat) then
            P_at=P
!            write(16,'(a,e15.7,$)') " P_at=",P_at
         elseif(i==nat+3) then
            P_lat=P-P_at
!            write(16,'(a,e15.7,$)') " P_lat=",P_lat
         endif       
        enddo
!        write(16,'(a,e15.7)') " P=",P

!Slight Modification of FIRE: Here we consider P as two different contributions, one from the atoms and one from the lattice
!        P=min(P_at,P_lat)
        fallnorm=0.d0
        do i=1,nat+3
         fallnorm=fallnorm+fall(1,i)**2+fall(2,i)**2+fall(3,i)**2
        enddo       
        fall=fall/sqrt(fallnorm)
        vallnorm=0.d0
        do i=1,nat+3
         vallnorm=vallnorm+vall(1,i)**2+vall(2,i)**2+vall(3,i)**2
        enddo       
        vallnorm=sqrt(vallnorm)
        vposcur=(1.d0-alpha)*vpospred+alpha*fall(:,1:nat)*vallnorm
        vlatcur=(1.d0-alpha)*vlatpred+alpha*fall(:,nat+1:nat+3)*vallnorm
        vvolcur=(1.d0-alpha)*vvolpred+alpha*fall(1,nat+1)*vallnorm
!Feedback on energy if gradient fails
        if(enthalpy.gt.enthalpy_min+1.d-3) then
             P=-1.d0
             enthalpy_min=enthalpy
        endif


         if(fmax.lt.cellfix_switch.and..not.cellfix_done.and.(.not.(any(fixlat).or.any(fixat).or.confine.ge.1))) then
!Only perform the cell correction once, presumably close to the end of the optimization run
             cellfix_done=.true.
             call correct_latvec(latcur,poscur,nat,correctalg,latvec_io)
             if(latvec_io.ne.0) then
               max_kpt=.false.
               ka1=0;kb1=0;kc1=0
             endif
         endif 

!FIRE timestep update ****************
!        if(P.le.0.d0.or.(fmax.lt.cellfix_switch.and..not.cellfix_done.and.(.not.(any(fixlat).or.any(fixat).or.confine.ge.1)))) then
        if(P.le.0.d0.or.latvec_io.ne.0) then
         latvec_io=0
         nstep=0
         dt=max(dt*fdec,dtmin)
!         accposcur=0.d0
!         acclatcur=0.d0
!         accvolcur=0.d0
         vposcur=0.d0
         vlatcur=0.d0
         vvolcur=0.d0
         accposprev=0.d0
         acclatprev=0.d0
         accvolprev=0.d0
         alpha=alphastart
         if((multiprec.and.itime.ge.parini%paropt_geopt%nit/2).or.&
          &(fmax.lt.1.0d0*tolmxf_switch)) max_kpt=.true.
         elseif(P.gt.0.d0 .and. nstep.gt.Nmin) then
           dt=min(dt*finc,dtmax)
!           alpha=max(alpha*falpha,0.1d0)!alpha*falpha
           alpha=alpha*falpha
        endif 
        nstep=nstep+1
!*************************************
!Feedback on cell
        if(md_type==4) then
          if(f_lat*fall(1,nat+1).gt.0) then
             alpha_lat=alpha_lat*1.05d0
          else
             alpha_lat=alpha_lat*0.7d0
          endif
        f_lat=fall(1,nat+1)
        endif


        call acceleration_fire(pressure_md,accposcur,acclatcur,accvolcur,vposcur,&
             &vlatcur,vvolcur,strten_in,fcart_in,latcur,amass,latmass,f0inv,md_type,nat)
!Eliminate the unnecessary atmic and cell components
        call elim_fixed_at(nat,accposcur);call elim_fixed_at(nat,vposcur)
if(md_type.ne.4) call elim_fixed_lat(latcur,acclatcur)
if(md_type.ne.4) call elim_fixed_lat(latcur,vlatcur)


     enthalpy_min=min(enthalpy,enthalpy_min)
     enddo 
!FIRE
     write(*,'(a,i5,es15.7,es15.7)') ' # EXIT FIRE ',itime,enthalpy,etot_in
     max_kpt=.false.

end subroutine GEOPT_FIRE_MHM 

!************************************************************************************

subroutine get_char_type(filename,nat,char_type,typat,ntypat)
use mod_interface
implicit none
integer:: nat,natin,iat,ntypat,nfound,typat(ntypat),ierror
character(40):: filename
real(8):: pos(3,nat),dproj(6)
character(2):: char_type(ntypat)
character(2):: char_single

char_type="NA" 
nfound=0

open(unit=46,file=trim(filename),iostat=ierror)
 if (ierror /= 0) then
    write(*,*) ' COULD not read file ',filename
    stop
 end if

read(46,*)natin
if(natin.ne.nat) stop "Number of atoms not consistent with abinit input"
read(46,*) dproj(1:3)
read(46,*) dproj(4:6)

do iat=1,nat
    read(46,*) pos(:,iat),char_single
    if(char_type(typat(iat))=="NA") then
!       write(*,'(a,a,i3,i3)') "New atom type character found ",trim(char_single),typat(iat),iat
       nfound=nfound+1
       if (nfound.gt.ntypat) stop "Too many different atom characters found!" 
       char_type(typat(iat))=trim(char_single)
    elseif(trim(char_type(typat(iat))).ne.trim(char_single)) then
       write(*,'(a,a,a,i3,i3)') "Already atom type character assigned:",char_type(typat(iat)),char_single,typat(iat),iat
    endif    
!    write(*,'(a,i3,a,a)') "Type character of atom ",iat,": ",char_type(typat(iat))
enddo
close(46)
!write(*,*) "Finished reading atom characters"

end subroutine

!************************************************************************************

 subroutine dproj2latvec(dproj,latvec)
 use mod_interface
 !This subroutine will convert the distance and projective representation of 
 !a periodic cell (dxx,dyx,dyy,dzx,dzy,dzz) into a 
 !lattice vektor format (vec1(:,1),vec2(:,2),vec3(:,3)) with dxx oriented into x direction
 implicit none
 real*8:: dproj(6),latvec(3,3)

 latvec(:,:)=0.d0
 latvec(1,1)=dproj(1)
 latvec(1,2)=dproj(2)
 latvec(2,2)=dproj(3)
 latvec(1,3)=dproj(4)
 latvec(2,3)=dproj(5)
 latvec(3,3)=dproj(6)
 return
 end subroutine

!************************************************************************************

 subroutine latvec2dproj(dproj,latvec,rotmat,rxyz,nat)
 use mod_interface, except_this_one=>norm
 !This subroutine will convert the lattice vector representation of thei
 !periodic cell (vec1,vec2,vec3) into the projective representation (dxx,dyx,dyy,dzx,dzy,dzz)
 !The cell will thus be rotated. The rotational matrix is stored in rotmat as an operator rotmat
 !and the atomic position rxyz are transformed into the new coordination sizstem as well
 implicit none
 integer,intent(in)  :: nat
 integer             :: iat,i
 real*8,intent(inout):: dproj(6),latvec(3,3),rotmat(3,3),rxyz(3,nat)
 real*8  :: tempvec(3),rotmat1(3,3),rotmat2(3,3),crossp(3),alpha,latvect(3,3)
 real*8  :: eps,axe(3),norm,rxyzt(3),rotmatt(3,3)
 eps=1.d-6
 !Calculating dxx
 dproj(1)=sqrt(latvec(1,1)*latvec(1,1)+latvec(2,1)*latvec(2,1)+latvec(3,1)*latvec(3,1))

 !Calculate the first rotation to align the first axis in x direction
 rotmat1(:,:)=0.d0
 do i=1,3
    rotmat1(i,i)=1.d0
 enddo
 tempvec(:)=0.d0
 tempvec(1)=1.d0
 !tempvec is the x-unit vector
 call cross_product(latvec(:,1),tempvec,crossp)
! if (abs(crossp(1)).lt.eps*1.d-1 .and. abs(crossp(2)).lt.eps*1.d-1 .and. abs(crossp(3)).lt.eps*1.d-1) goto 1001 !no rotation needed
 if (abs(crossp(1)).lt.eps*1.d-1 .and. abs(crossp(2)).lt.eps*1.d-1 .and. abs(crossp(3)).lt.eps*1.d-1) then
!Check if the latvec points along the positive x-rection
    alpha=dot_product(tempvec(:),latvec(:,1))/sqrt(dot_product(latvec(:,1),latvec(:,1)))
    if(alpha.gt.0.d0) then
        goto 1001 !no rotation needed
    else
!Rotate along the y-axis
    axe(:)=0.d0
    axe(2)=1.d0
    alpha=dacos(-1.d0)
    call rotation(rotmat1,alpha,axe)
    latvec(:,:)=matmul(rotmat1(:,:),latvec(:,:))
        goto 1001
    endif
 endif
 norm=sqrt(crossp(1)*crossp(1)+crossp(2)*crossp(2)+crossp(3)*crossp(3))
 axe(:)=crossp(:)/norm
 alpha=dacos(dot_product(tempvec,latvec(:,1))/dproj(1))
 call rotation(rotmat1,alpha,axe)
 latvec(:,:)=matmul(rotmat1(:,:),latvec(:,:))
! call DGEMM('N','N',3,3,3,1.d0,rotmat,3,latvec,3,0.d0,latvect,3)
! latvec=latvect

 1001 continue
 if(latvec(2,1).gt.eps .or. latvec(3,1).gt.eps) then
 write(*,*) "Error in 1. rotation",latvec(2,1),latvec(3,1)
 stop
 endif
 !Calculate the second rotation to align the second axis in xy plane 
 rotmat2(:,:)=0.d0
 do i=1,3
    rotmat2(i,i)=1.d0
 enddo
! axe(:)=0.d0
! axe(1)=1.d0
 tempvec(:)=latvec(:,2)
 tempvec(1)=0.d0
 call cross_product(tempvec,(/0.d0,1.d0,0.d0/),axe)

! if (abs(axe(1)).lt.eps*1.d-1 .and. abs(axe(2)).lt.eps*1.d-1 .and. abs(axe(3)).lt.eps*1.d-1 .and. tempvec(2).gt.0.d0) goto 1002 !no rotation needed
 if (abs(axe(1)).lt.eps*1.d-1 .and. abs(axe(2)).lt.eps*1.d-1 .and. abs(axe(3)).lt.eps*1.d-1) then
 if (tempvec(2).gt.0.d0) then  !     goto 1002 !no rotation needed
!Check if the latvec points along the positive x-rection
    goto 1002 !no rotation needed
 else
!Rotate along the y-axis
    axe(:)=0.d0
    axe(1)=1.d0
    alpha=dacos(-1.d0)
    call rotation(rotmat2,alpha,axe)
    latvec(:,:)=matmul(rotmat2(:,:),latvec(:,:))
        goto 1002
 endif
 endif

 norm=sqrt(axe(1)*axe(1)+axe(2)*axe(2)+axe(3)*axe(3))
 axe=axe/norm
 call cross_product(axe,latvec(:,2),crossp)

 if (abs(crossp(1)).lt.eps*1.d-1 .and. abs(crossp(2)).lt.eps*1.d-1 .and. crossp(3).gt.0.d0) goto 1002
 norm=sqrt(axe(1)*axe(1)+axe(2)*axe(2)+axe(3)*axe(3))
 axe=axe/norm
 call cross_product(axe,latvec(:,2),crossp)
! if (abs(crossp(1)).lt.eps*1.d-1 .and. abs(crossp(2)).lt.eps*1.d-1 .and. crossp(3).gt.0.d0) goto 1002 !no rotation needed
 if (abs(axe(1)).lt.eps*1.d-1 .and. abs(axe(2)).lt.eps*1.d-1 .and. abs(axe(3)).lt.eps*1.d-1 .and. tempvec(2).gt.0.d0) goto 1002 !no rotation needed
 norm=sqrt(axe(1)*axe(1)+axe(2)*axe(2)+axe(3)*axe(3))
 axe=axe/norm
 call cross_product(axe,latvec(:,2),crossp)
! if (abs(crossp(1)).lt.eps*1.d-1 .and. abs(crossp(2)).lt.eps*1.d-1 .and. crossp(3).gt.0.d0) goto 1002 !no rotation needed

 if (abs(crossp(1)).lt.eps*1.d-1 .and. abs(crossp(2)).lt.eps*1.d-1 .and. crossp(3).gt.0.d0) then
!Check if the latvec points along the positive y-rection
    alpha=dot_product(axe(:),latvec(:,2))/sqrt(dot_product(latvec(:,2),latvec(:,2)))
    if(alpha.gt.0.d0) then
        goto 1002 !no rotation needed
    else
!Rotate along the y-axis
    axe(:)=0.d0
    axe(1)=1.d0
    alpha=dacos(-1.d0)
    call rotation(rotmat2,alpha,axe)
    latvec(:,:)=matmul(rotmat2(:,:),latvec(:,:))
        goto 1002
    endif
 endif
 norm=sqrt(tempvec(2)*tempvec(2)+tempvec(3)*tempvec(3))
 alpha=dot_product((/0.d0,1.d0,0.d0/),tempvec(:))/norm
 alpha=dacos(alpha)
 call rotation(rotmat2,alpha,axe)
 latvec(:,:)=matmul(rotmat2(:,:),latvec(:,:))
! call DGEMM('N','N',3,3,3,1.d0,rotmat,3,latvec,3,0.d0,latvect,3)
! latvec=latvect
 1002 continue
 if(latvec(3,2).gt.eps) then! stop "Error in 2. rotation" 
! write(*,*) latvec(:,1)
! write(*,*) latvec(:,2)
 write(*,*) "Error in 2. rotation"
 stop
 endif
 if(latvec(3,3).lt.0.d0) stop "Error in orientation of the cell"

 !The total rotational matrix:
 rotmat=matmul(rotmat2(:,:),rotmat1(:,:))
! call DGEMM('N','N',3,3,3,1.d0,rotmat2,3,rotmat1,3,0.d0,rotmatt,3)
! rotmat=rotmatt 

 !Apply rotation on all atoms
 do iat=1,nat
     rxyz(:,iat)=matmul(rotmat,rxyz(:,iat))
!     call DGEMM('N','N',3,3,3,1.d0,rotmat,3,latvec,3,0.d0,latvect,3)
!     rxyz(:,iat)=rxyzt(:)
 enddo

 !Calculate all other elements of dproj
 dproj(2)=latvec(1,2)
 dproj(3)=latvec(2,2)
 dproj(4)=latvec(1,3)
 dproj(5)=latvec(2,3)
 dproj(6)=latvec(3,3)
 end subroutine

!************************************************************************************

 subroutine cross_product(a,b,crossp)
 use mod_interface
 !a very simple implementation of the cross product
 implicit none
 real(8)::a(3),b(3)
 real(8)::crossp(3)
 crossp(1)=a(2)*b(3)-a(3)*b(2)
 crossp(2)=a(3)*b(1)-a(1)*b(3)
 crossp(3)=a(1)*b(2)-a(2)*b(1)
 return
 end subroutine

 subroutine dot_p(a,b,dotp)
 use mod_interface
 !a very simple implementation of the dot product
 implicit none
 real(8)::a(3),b(3)
 real(8)::dotp(3)
 integer::i
 dotp=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
 end subroutine

!************************************************************************************

 subroutine rotation(rotmat,angle,axe)
 use mod_interface
 !This subroutine will calculate the rotational matrix rotmat for a
 !3-dim vector around an axis 'axe' by the angle 'angle'.
 implicit none
 real(8),INTENT(IN) :: angle
 real(8),INTENT(IN) :: axe(3)
 real(8):: rotator(3,3)
 real(8):: rotmat(3,3),cosang,sinang
 cosang=dcos(angle)
 sinang=dsin(angle)


 !Define Rotation Matrix
 rotator(1,1)=cosang+(axe(1)**2)*(1.d0-cosang)
 rotator(1,2)=axe(1)*axe(2)*(1.d0-cosang)-axe(3)*sinang
 rotator(1,3)=axe(1)*axe(3)*(1.d0-cosang)+axe(2)*sinang
                                                
 rotator(2,1)=axe(2)*axe(1)*(1.d0-cosang)+axe(3)*sinang
 rotator(2,2)=cosang+(axe(2)**2)*(1.d0-cosang) 
 rotator(2,3)=axe(2)*axe(3)*(1.d0-cosang)-axe(1)*sinang
                                              
 rotator(3,1)=axe(3)*axe(1)*(1.d0-cosang)-axe(2)*sinang
 rotator(3,2)=axe(3)*axe(2)*(1.d0-cosang)+axe(1)*sinang
 rotator(3,3)=cosang+(axe(3)**2)*(1.d0-cosang)
 rotmat(:,:)=rotator(:,:)

 !do i=1,3
 !   vector2(i)=rotator(i,1)*vector(1)+rotator(i,2)*vector(2)+rotator(i,3)*vector(3)
 !enddo
 !vector(:)=vector2(:)
 end subroutine rotation

!************************************************************************************

 subroutine rxyz_int2cart(latvec,rxyzint,rxyzcart,nat)
 use mod_interface
 !This subrouine will convert the internal coordinates into cartesian coordinates
 implicit none
 real(8):: rxyzint(3,nat), rxyzcart(3,nat),latvec(3,3)
 integer:: nat,iat
 do iat=1,nat
  rxyzcart(:,iat)=matmul(latvec,rxyzint(:,iat))
 enddo
 end subroutine rxyz_int2cart 

!************************************************************************************

 subroutine rxyz_cart2int(latvec,rxyzint,rxyzcart,nat)
 use mod_interface
 !This subrouine will convert the internal coordinates into cartesian coordinates
 implicit none
 real(8):: rxyzint(3,nat), rxyzcart(3,nat),latvec(3,3),latvecinv(3,3)
 integer:: nat,iat
 call invertmat(latvec,latvecinv,3)
 do iat=1,nat
  rxyzint(:,iat)=matmul(latvecinv,rxyzcart(:,iat))
 enddo
 end subroutine rxyz_cart2int

!************************************************************************************

 subroutine invertmat(mat,matinv,n)
 use mod_interface
 implicit none
 real(8),intent(in) :: mat(n,n)
 integer               :: n
 real(8),allocatable   :: WORK(:)
 real(8)               :: matinv(n,n),det(3),a(n,n),div
 integer               :: IPIV(n), INFO
 integer               :: LDWORK
 !Here only for a 3*3 matrix
 if (n==3) then
 a=mat
 div=(a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)+a(1,2)*a(2,3)*a(3,1)+&
 &a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)) 
 div=1.d0/div
      matinv(1,1) = (a(2,2)*a(3,3)-a(2,3)*a(3,2))*div
      matinv(1,2) =-(a(1,2)*a(3,3)-a(1,3)*a(3,2))*div
      matinv(1,3) = (a(1,2)*a(2,3)-a(1,3)*a(2,2))*div
      matinv(2,1) =-(a(2,1)*a(3,3)-a(2,3)*a(3,1))*div
      matinv(2,2) = (a(1,1)*a(3,3)-a(1,3)*a(3,1))*div
      matinv(2,3) =-(a(1,1)*a(2,3)-a(1,3)*a(2,1))*div
      matinv(3,1) = (a(2,1)*a(3,2)-a(2,2)*a(3,1))*div
      matinv(3,2) =-(a(1,1)*a(3,2)-a(1,2)*a(3,1))*div
      matinv(3,3) = (a(1,1)*a(2,2)-a(1,2)*a(2,1))*div
 else
 !General n*n matrix 
 matinv=mat
 allocate(WORK(n))
 call  DGETRF( n, n, matinv, n, IPIV, INFO )
 if (info.ne.0) stop "Error in DGETRF"
 LDWORK=-1
 call  DGETRI( n, matinv, n, IPIV, WORK,LDWORK , INFO )
 LDWORK=WORK(1)
 deallocate(WORK)
 allocate(WORK(LDWORK))
 call  DGETRI( n, matinv, n, IPIV, WORK,LDWORK , INFO )
 if (info.ne.0) stop "Error in DGETRI"
 endif
 end subroutine

!************************************************************************************

 subroutine latvec2acell_rprim(latvec,acell,rprim)
 use mod_interface
!This routine will split up the latvec into two parts, 
!acell and rprim, where rprim is normalized to 1
 implicit none
 real(8):: latvec(3,3), rprim(3,3), acell(3)
 integer:: i
 do i=1,3
   acell(i)=sqrt(latvec(1,i)**2+latvec(2,i)**2+latvec(3,i)**2)
   rprim(:,i)=latvec(:,i)/acell(i)
 enddo
 end subroutine latvec2acell_rprim

!************************************************************************************

 subroutine acell_rprim2latvec(latvec,acell,rprim)
 use mod_interface
!This routine will combine 
!acell and rprim to latvec
 implicit none
 real(8):: latvec(3,3), rprim(3,3), acell(3)
 integer:: i
 do i=1,3
   latvec(:,i)=acell(i)*rprim(:,i)
 enddo
 end subroutine acell_rprim2latvec

!************************************************************************************

!subroutine get_enthalpy(acell,rprim,energy,pressure,enthalpy)
subroutine get_enthalpy(latvec,energy,pressure,enthalpy)
use mod_interface
!This routine will compute the enthalpy within the given units
implicit none
integer:: nat,natin,iat
character(40):: filename,units
real(8):: acell(3),v(3,3),ucvol,pressure,latvec(3,3),energy,enthalpy

!latvec(:,1)=acell(1)*rprim(:,1)
!latvec(:,2)=acell(2)*rprim(:,2)
!latvec(:,3)=acell(3)*rprim(:,3)

!Compute cell volume
 v=latvec
 ucvol= v(1,1)*v(2,2)*v(3,3)-v(1,1)*v(2,3)*v(3,2)-v(1,2)*v(2,1)*v(3,3)+&
        v(1,2)*v(2,3)*v(3,1)+v(1,3)*v(2,1)*v(3,2)-v(1,3)*v(2,2)*v(3,1)
enthalpy=energy+pressure*ucvol
end subroutine

!************************************************************************************

function round(enerd,accur)
  implicit none
  real*8 enerd,accur,round
  integer*8 ii
  ii=enerd/accur
  round=ii*accur
  !           write(*,'(a,1pe24.17,1x,i17,1x,1pe24.17)') 'enerd,ii,round',enerd,ii,round
  return
end function round

!************************************************************************************

subroutine wtioput(ediff,ekinetic,ekinetic_max,nsoften)
  use mod_interface
  implicit none
  integer:: nsoften
  real(8):: ediff, ekinetic,ekinetic_max
  open(unit=11,file='ioput',status='unknown')
  write(11,'(3(1x,1pe24.17)1x,a)') ediff,ekinetic,ekinetic_max,' ediff, temperature, maximal temperature'
  close(11)
END SUBROUTINE wtioput

!************************************************************************************

subroutine hunt(xx,n,x,jlo)
  use mod_interface
  implicit none
  !C x is in interval [xx(jlo),xx(jlow+1)[ ; xx(0)=-Infinity ; xx(n+1) = Infinity
  !Arguments
  integer :: jlo,n
  real(kind=8) :: x,xx(n)
  !Local variables
  integer :: inc,jhi,jm
  logical :: ascnd
  if (n.le.0) stop 'hunt'
  if (n == 1) then
     if (x.ge.xx(1)) then
        jlo=1
     else
        jlo=0
     endif
     return
  endif
  ascnd=xx(n).ge.xx(1)
  if(jlo.le.0.or.jlo.gt.n)then
     jlo=0
     jhi=n+1
     goto 3
  endif
  inc=1
  if(x.ge.xx(jlo).eqv.ascnd)then
1    continue
     jhi=jlo+inc
     if(jhi.gt.n)then
        jhi=n+1
     else if(x.ge.xx(jhi).eqv.ascnd)then
        jlo=jhi
        inc=inc+inc
        goto 1
     endif
  else
     jhi=jlo
2    continue
     jlo=jhi-inc
     if(jlo.lt.1)then
        jlo=0
     else if(x.lt.xx(jlo).eqv.ascnd)then
        jhi=jlo
        inc=inc+inc
        goto 2
     endif
  endif
3 continue
  if(jhi-jlo == 1)then
     if(x == xx(n))jlo=n
     if(x == xx(1))jlo=1
     return
  endif
  jm=(jhi+jlo)/2
  if(x.ge.xx(jm).eqv.ascnd)then
     jlo=jm
  else
     jhi=jm
  endif
  goto 3
END SUBROUTINE hunt

!************************************************************************************
subroutine insert(nlminx,nlmin,fp_len,nat,k_e_wpos,e_wpos,ent_wpos,fp_wpos,wpos_red,&
  &wpos_latvec,wpos_fcart,wpos_strten,spg_wpos,spgtol_wpos,fdos_wpos,&
  &e_arr,ent_arr,fp_arr,pl_arr,lat_arr,f_arr,str_arr,spg_arr,spgtol_arr,dos_arr,ct_arr)
  ! inserts the energy re_wpos at position k_e_wpos and shifts up all other energies
!  implicit real*8 (a-h,o-z)
  use mod_interface
  implicit none
  integer:: fp_len,ct_arr(nlminx),spg_arr(nlminx),nat,iat,spg_wpos
  integer:: k, nlmin, k_e_wpos, nlminx,i
  real(8):: e_wpos, ent_wpos, wpos_red(3,nat),wpos_latvec(3,3),spgtol_wpos,fdos_wpos,fp_wpos(fp_len)
  real(8):: e_arr(nlminx),ent_arr(nlminx),fp_arr(fp_len,nlminx),pl_arr(3,nat,nlminx),f_arr(3,nat,nlminx)
  real(8):: lat_arr(3,3,nlminx),spgtol_arr(nlminx),dos_arr(nlminx),fp1(fp_len),fp2(fp_len),str_arr(6,nlminx)
  real(8):: wpos_fcart(3,nat),wpos_strten(6)
              do k=nlmin-1,k_e_wpos+1,-1
                 e_arr(k+1)=e_arr(k)
                 ent_arr(k+1)=ent_arr(k)
                 do i=1,fp_len
                 fp_arr(i,k+1)=fp_arr(i,k)
                 enddo
                 do iat=1,nat
                 pl_arr(1,iat,k+1)=pl_arr(1,iat,k)
                 pl_arr(2,iat,k+1)=pl_arr(2,iat,k)
                 pl_arr(3,iat,k+1)=pl_arr(3,iat,k)
                 f_arr(1,iat,k+1)=f_arr(1,iat,k)
                 f_arr(2,iat,k+1)=f_arr(2,iat,k)
                 f_arr(3,iat,k+1)=f_arr(3,iat,k)
                 enddo
                 str_arr(:,k+1)=str_arr(:,k)
                 do i=1,3
                 lat_arr(1,i,k+1)=lat_arr(1,i,k)
                 lat_arr(2,i,k+1)=lat_arr(2,i,k)
                 lat_arr(3,i,k+1)=lat_arr(3,i,k)
                 enddo
                 spg_arr(k+1)=spg_arr(k)
                 spgtol_arr(k+1)=spgtol_arr(k)
                 dos_arr(k+1)=dos_arr(k)
                 ct_arr(k+1)=ct_arr(k)
              enddo
              e_arr(k_e_wpos+1)=e_wpos
              ent_arr(k_e_wpos+1)=ent_wpos
              ct_arr(k_e_wpos+1)=1
                 do i=1,fp_len
                 fp_arr(i,k_e_wpos+1)=fp_wpos(i)
                 enddo
                 do iat=1,nat
                 pl_arr(1,iat,k_e_wpos+1)=wpos_red(1,iat)
                 pl_arr(2,iat,k_e_wpos+1)=wpos_red(2,iat)
                 pl_arr(3,iat,k_e_wpos+1)=wpos_red(3,iat)
                 f_arr(1,iat,k_e_wpos+1)=wpos_fcart(1,iat)
                 f_arr(2,iat,k_e_wpos+1)=wpos_fcart(2,iat)
                 f_arr(3,iat,k_e_wpos+1)=wpos_fcart(3,iat)
                 enddo
                 str_arr(:,k_e_wpos+1)=wpos_strten(:)
                 do i=1,3
                 lat_arr(1,i,k_e_wpos+1)=wpos_latvec(1,i)
                 lat_arr(2,i,k_e_wpos+1)=wpos_latvec(2,i)
                 lat_arr(3,i,k_e_wpos+1)=wpos_latvec(3,i)
                 enddo
                 spg_arr(k_e_wpos+1)=spg_wpos
                 spgtol_arr(k_e_wpos+1)=spgtol_wpos
                 dos_arr(k_e_wpos+1)=fdos_wpos
       write(*,*) '  -----   INSERT -----------'
       return
END SUBROUTINE insert

!!subroutine insert(nlminx,nlmin,k_e_wpos,rent_wpos,re_wpos,earr)
!!  ! inserts the energy re_wpos at position k_e_wpos and shifts up all other energies
!!!  implicit real*8 (a-h,o-z)
!!  implicit none
!!  real(8):: earr(0:nlminx,3), re_wpos, rent_wpos
!!  integer:: k, nlmin, k_e_wpos, nlminx
!!  do k=nlmin-1,k_e_wpos+1,-1
!!     earr(k+1,1)=earr(k,1)
!!     earr(k+1,2)=earr(k,2)
!!     earr(k+1,3)=earr(k,3)
!!  enddo
!!  earr(k_e_wpos+1,1)=rent_wpos
!!  earr(k_e_wpos+1,2)=re_wpos
!!  earr(k_e_wpos+1,3)=1.d0
!!  return
!!END SUBROUTINE insert

!************************************************************************************

subroutine save_low_conf(nat,npmin,npminx,ent_wpos,e_wpos,pos,latvec,spg,spgtol,fdos,elocmin,poslocmin,latlocmin)
  !C save configuration if it is among the lowest ones in energy
  use mod_interface
!  implicit real*8 (a-h,o-z)
  implicit none
  integer:: iat,nat, npmin, npminx, kmax, k 
  real(8):: e_wpos, ent_wpos, emax,spg,spgtol,fdos
  real(8):: elocmin(npminx,5)
  real(8):: pos(3,nat),latvec(3,3),poslocmin(3,nat,npminx),latlocmin(3,3,npminx)

  if (npmin.le.npminx) then
     kmax=npmin
     elocmin(kmax,1)=ent_wpos
     elocmin(kmax,2)=e_wpos
     elocmin(kmax,3)=spg
     elocmin(kmax,4)=spgtol
     elocmin(kmax,5)=fdos
     do iat=1,nat
        poslocmin(1,iat,kmax)=pos(1,iat)
        poslocmin(2,iat,kmax)=pos(2,iat)
        poslocmin(3,iat,kmax)=pos(3,iat)
     enddo
     latlocmin(:,:,kmax)=latvec(:,:)
  else
     ! find configuration kmax that is highest in energy
     emax=-1.d100
     do k=1,npminx
        if (elocmin(k,1).gt.emax) then
           emax=elocmin(k,1)
           kmax=k
        endif
     enddo
     if (ent_wpos.lt.elocmin(kmax,1)) then
        elocmin(kmax,1)=ent_wpos
        elocmin(kmax,2)=e_wpos
        elocmin(kmax,3)=spg
        elocmin(kmax,4)=spgtol
        elocmin(kmax,5)=fdos
        do iat=1,nat
           poslocmin(1,iat,kmax)=pos(1,iat)
           poslocmin(2,iat,kmax)=pos(2,iat)
           poslocmin(3,iat,kmax)=pos(3,iat)
        enddo
        latlocmin(:,:,kmax)=latvec(:,:)
     endif
  endif
  return
END SUBROUTINE save_low_conf

!************************************************************************************

subroutine winter(parini,nat,units,ent_pos,e_pos,pos_red,pos_latvec,pos_fcart,pos_strten,nlminx,nlmin,npminx,& 
   &ent_arr,e_arr,ct_arr,spg_arr,spgtol_arr,dos_arr,pl_arr,lat_arr,f_arr,str_arr,fp_arr,fp_len,ent_delta,fp_delta,& 
   &eref,ediff,ekinetic,ekinetic_max,dt,nsoften,char_type,ntypat,typat,fixat,fixlat,pressure)
  use mod_interface
  use mod_parini, only: typ_parini
  implicit none
  type(typ_parini), intent(in):: parini
  integer, intent(in) :: nlminx,nlmin,nsoften,nat,npminx,fp_len
  real(8), intent(in) :: eref,ediff,ekinetic,dt,e_pos,ent_pos,ekinetic_max,ent_delta,fp_delta
  real(8), intent(in) :: pos_latvec(3,3) 
  real(8), intent(in) :: pos_strten(6) 
  real(8), dimension(nlminx),      intent(in) :: ent_arr,e_arr,spgtol_arr,dos_arr
  real(8), dimension(3,3,nlminx),  intent(in) :: lat_arr
  real(8), dimension(6,nlminx),    intent(in) :: str_arr
  real(8), dimension(3,nat,nlminx),intent(in) :: pl_arr,f_arr
  real(8), dimension(fp_len,nlminx),intent(in):: fp_arr
  integer, dimension(nlminx),      intent(in) :: ct_arr,spg_arr
  character(2), intent(in):: char_type(ntypat) 
  integer, intent(in):: ntypat 
  integer, intent(in):: typat(nat) 
  real(8), intent(in):: pressure 
  real(8), intent(in):: pos_red(3,nat) 
  real(8), intent(in):: pos_fcart(3,nat) 
  logical :: fixat(nat),fixlat(7)

  
!local variables
  character(len=5) :: fn
  character(len=40) :: filename
  character(len=40) :: units 
  integer :: mm,k,n_arr
!  real(8):: acell(3), rprim(3,3)

!     call latvec2acell_rprim(pos_latvec,acell,rprim)
     filename="poscur.ascii"
     call write_atomic_file_ascii(parini,filename,nat,units,pos_red,pos_latvec,pos_fcart,pos_strten,char_type,&
          &ntypat,typat,fixat,fixlat,e_pos,pressure,ent_pos,e_pos)
!     call write_atomic_file('poscur',re_pos,pos,at,'')
     write(*,*) ' wrote poscur.ascii for RESTART',ent_pos
     

     
     open(unit=12,file='earr.dat',status='unknown')
     mm=min(nlmin,nlminx)
     write(12,'(2(i10),a)') mm,nlminx+1,&
          '          # No. of minima already found, no. of minima to be found in consecutive run'
!!     write(12,'(e24.17,1x,a)') eref,'   eref'
!!     write(12,'(e24.17,1x,a)') accur,'   accur'
       write(12,'(2(e14.6),1x,a)') ent_delta,fp_delta," # delta_enthalpy, delta_fingerprint"
     do k=1,mm
!        write(12,'(e24.17,1x,e24.17,1x,1pe17.10)') earr(k,1),earr(k,2),earr(k,3)
        write(12,'(i5.5,1x,e24.17,1x,e24.17,1x,i5,1x,i5,1x,e15.7,e15.7)') &
                 &k,ent_arr(k),e_arr(k),ct_arr(k),spg_arr(k),spgtol_arr(k),dos_arr(k)
        if(k.le.npminx) then
! write poslow files
           write(fn,'(i5.5)') k
!generate filename and open files
           filename="poslow"//fn//".ascii"
           call write_atomic_file_ascii(parini,filename,nat,units,pl_arr(1:3,1:nat,k),lat_arr(:,:,k),f_arr(1:3,1:nat,k),str_arr(:,k),&
           &char_type,ntypat,typat,fixat,fixlat,e_arr(k),pressure,ent_arr(k),e_arr(k))
        endif
     enddo
     write(*,*) ' wrote poslow files',nlmin
     write(*,*) ' wrote earr.dat for  RESTART'
     close(12)
     
     call wtioput(ediff,ekinetic,ekinetic_max,nsoften)
     write(*,*) ' wrote ioput for  RESTART'

!Write binaries
     n_arr=3*nat*mm
     filename="poslow.bin"
     call bin_write(filename,pl_arr,n_arr)
     filename="fcart.bin"
     call bin_write(filename,f_arr,n_arr)
     n_arr=6*mm
     filename="strten.bin"
     call bin_write(filename,str_arr,n_arr)
     n_arr=3*3*mm
     filename="latvec.bin"
     call bin_write(filename,lat_arr,n_arr)
     n_arr=fp_len*mm
     filename="fp.bin"
     call bin_write(filename,fp_arr,n_arr)
     write(*,*) ' wrote binary files poslow.bin, fcart.bin, strten.bin, latvec.bin and fp.bin for RESTART'


!!     call wtpos(npminx,nlminx,nlmin,npmin,poslocmin,latlocmin,earr,elocmin,char_type,&
!!          &ntypat,typat,fixat,fixlat,pressure,units,nat,eref)

END SUBROUTINE winter

!************************************************************************************

!!!subroutine wtpos(npminx,nlminx,nlmin,npmin,poslocmin,latlocmin,earr,elocmin,char_type,&
!!!           &ntypat,typat,fixat,fixlat,pressure,units,nat,eref)
!!!  implicit none
!!!  integer, intent(in) :: npminx,nlminx,nlmin,npmin
!!!  real(8), dimension(npminx,5), intent(in) :: elocmin
!!!  real(8), dimension(0:nlminx,3), intent(in) :: earr
!!!  real(8), dimension(3,nat,npminx), intent(in) :: poslocmin
!!!  real(8), dimension(3,3,npminx), intent(in) :: latlocmin
!!!  character(2), intent(in):: char_type(ntypat) 
!!!  integer, intent(in):: ntypat,nat
!!!  integer, intent(in):: typat(nat) 
!!!  real(8), intent(in):: pressure 
!!!  logical :: fixat(nat),fixlat(7)
!!!
!!!  !local variables
!!!  character(len=5) :: fn
!!!  character(len=40) :: filename
!!!  character(len=40) :: units 
!!!  integer :: k,kk,i
!!!  real(8):: acell(3), eref
!!!  real(8),allocatable:: elocmin_sorted(:,:)
!!!  real(8):: tmp1(5)
!!!  allocate(elocmin_sorted(min(npminx,npmin),5))
!!!       write(*,*) 'nlmin,nlminx,npmin,npminx',nlmin,nlminx,npmin,npminx
!!!       do i=1,min(40,nlmin,nlminx)
!!!         write(*,'(i4,e24.17)') i,earr(i,1)
!!!       enddo
!!!  do k=1,min(npmin,npminx)
!!!             write(*,'(a,i4,e24.17)') 'k,elocmin(k,1)',k,elocmin(k,1)
!!!
!!!
!!!     !C Classify the configuration in the global ranking
!!!     kk=0
!!!     find_kk : do
!!!        kk=kk+1
!!!        if (kk > min(nlmin,nlminx)) then 
!!!           write(*,*) 'ranking error for',k
!!!           stop 
!!!        endif
!!!!        if (earr(kk,1) == elocmin(k,1)) exit find_kk
!!!        if (abs(earr(kk,1) - elocmin(k,1)) .lt. 1.d-12 ) then 
!!!             write(*,*) 'match ',abs(earr(kk,1) - elocmin(k,1)),kk,k
!!!             exit find_kk
!!!        endif
!!!     end do find_kk
!!!
!!!     if (kk <= npminx) then
!!!
!!!                write(*,'(a,i4,i4,1x,1pe21.14)') 'k,kk,elocmin(k,1)',k,kk,elocmin(k,1)
!!!
!!!        !C generate filename and open files
!!!           write(fn,'(i5.5)') kk
!!!           filename="poslow"//fn//".ascii"
!!!           call write_atomic_file_ascii(filename,nat,units,poslocmin(1:3,1:nat,k),latlocmin(:,:,k),char_type,&
!!!           &ntypat,typat,fixat,fixlat,elocmin(k,2)+eref,pressure,elocmin(k,1),elocmin(k,2))
!!!           elocmin_sorted(kk,:)=elocmin(k,:)
!!!     endif
!!!
!!!  end do
!!!
!!!!Write also the file strlist.dat 
!!!  open(unit=787,file="strlist.dat")
!!!  do k=1,min(npminx,npmin)
!!!           write(787,'(i5.5,1x,es25.15,es25.15,1x,i5,es15.7,es15.7)')  k, elocmin_sorted(k,1:2),&
!!!           &int(elocmin_sorted(k,3)),elocmin_sorted(k,4:5)   
!!!  enddo
!!!  close(787)
!!!  deallocate(elocmin_sorted)
!!!
!!!END SUBROUTINE wtpos

!************************************************************************************

subroutine torque_cell(latvec0,vlat,torquenrm)
use mod_interface
implicit none
real(8), intent(in)    :: latvec0(3,3)
real(8), intent(inout) :: vlat(3,3),torquenrm
real(8) :: torque(3),crossp(3),sx,sy,sz,tmax,cx,cy,cz
integer :: i,ii,it,itmax
real(8) :: unitmat(3,3),rotmat(3,3),rotmatall(3,3),xaxis(3),axis(3),latvec(3,3),tnorm,angle,axisnorm
       latvec=latvec0
!Initialize

       torque=0.d0
       do i=1,3
       call cross_product(latvec(:,i),vlat(:,i),crossp)
       torque=torque+crossp
       enddo
       torquenrm=sqrt(torque(1)**2+torque(2)**2+torque(3)**2)
end subroutine torque_cell

!************************************************************************************

!Various methods to initialize the velocities for the MD part of Minima Hopping
!GAUSSIAN DISTRIBUTION**********************************************************
      subroutine gausdist(nat,vxyz,amass)
      use mod_interface
!generates 3*nat random numbers distributed according to  exp(-.5*vxyz**2)
      implicit none!real*8 (a-h,o-z)
      real:: s1,s2
      real(8):: t1,t2,tt,amass(nat)
! On Intel the random_number can take on the values 0. and 1.. To prevent overflow introduce eps
      real(8),parameter:: eps=1.d-8
      real(8),dimension(3*nat)::  vxyz
      integer:: nat,i
      do i=1,3*nat-1,2
        call random_number(s1)
        t1=eps+(1.d0-2.d0*eps)*dble(s1)
        call random_number(s2)
        t2=dble(s2)
        tt=sqrt(-2.d0*log(t1))
        vxyz(i)=tt*cos(6.28318530717958648d0*t2)
        vxyz(i+1)=tt*sin(6.28318530717958648d0*t2)
      enddo
        call random_number(s1)
        t1=eps+(1.d0-2.d0*eps)*dble(s1)
        call random_number(s2)
        t2=dble(s2)
        tt=sqrt(-2.d0*log(t1))
        vxyz(3*nat)=tt*cos(6.28318530717958648d0*t2)

        call elim_moment(nat,vxyz,amass)
      return
      end subroutine

!************************************************************************************

!GAUSSIAN DISTRIBUTION FOR THE CELL VECTORS***************************************
      subroutine gausdist_cell(latvec,vlat)
      use mod_interface
! generates 3*3 random numbers distributed according to  exp(-.5*vxyz**2) for the cell vectors
      implicit none
      integer:: i
      real:: s1,s2
      real(8) :: t1,t2,tt
! On Intel the random_number can take on the values 0. and 1.. To prevent overflow introduce eps
      real(8),parameter:: eps=1.d-8
      real(8)::  vlat(9),latvec(9)

      do i=1,3*3-1,2
        call random_number(s1)
        t1=eps+(1.d0-2.d0*eps)*dble(s1)
        call random_number(s2)
        t2=dble(s2)
        tt=sqrt(-2.d0*log(t1))
        vlat(i)=tt*cos(6.28318530717958648d0*t2)
        vlat(i+1)=tt*sin(6.28318530717958648d0*t2)
      enddo
        call random_number(s1)
        t1=eps+(1.d0-2.d0*eps)*dble(s1)
        call random_number(s2)
        t2=dble(s2)
        tt=sqrt(-2.d0*log(t1))
        vlat(3*3)=tt*cos(6.28318530717958648d0*t2)

        call elim_torque_cell(latvec,vlat)
      return
      end subroutine

!************************************************************************************

        subroutine elim_moment(nat,vxyz,atmass)
        use mod_interface
        implicit none
        real(8):: vxyz(3,nat),sx,sz,sy,atmass(nat)
        integer:: iat,nat       
 
        sx=0.d0 ; sy=0.d0 ; sz=0.d0
        do iat=1,nat
        sx=sx+vxyz(1,iat)*atmass(iat)
        sy=sy+vxyz(2,iat)*atmass(iat)
        sz=sz+vxyz(3,iat)*atmass(iat)
        enddo
        sx=sx/nat ; sy=sy/nat ; sz=sz/nat
        do iat=1,nat
        vxyz(1,iat)=vxyz(1,iat)-sx/atmass(iat)
        vxyz(2,iat)=vxyz(2,iat)-sy/atmass(iat)
        vxyz(3,iat)=vxyz(3,iat)-sz/atmass(iat)
        enddo
        return
        end subroutine

!************************************************************************************

subroutine elim_torque_cell(latvec0,vlat)
use mod_interface
implicit none
real(8), intent(in)    :: latvec0(3,3)
real(8), intent(inout) :: vlat(3,3)
real(8) :: torque(3),crossp(3),sx,sy,sz,tmax,cx,cy,cz
integer :: i,ii,it,itmax
real(8) :: unitmat(3,3),rotmat(3,3),rotmatall(3,3),xaxis(3),axis(3),latvec(3,3),tnorm,angle,axisnorm
       latvec=latvec0
       itmax=5000
!Initialize
       unitmat=0.d0
       do i=1,3
       unitmat(i,i)=1.d0
       enddo
       rotmatall=unitmat
       xaxis=(/1.d0,0.d0,0.d0/)

it=0
do
it=it+1
       torque=0.d0
       do i=1,3
       call cross_product(latvec(:,i),vlat(:,i),crossp)
       torque=torque+crossp
       enddo
       !write(*,'(3(e11.3),i6)')torque,loop
       if (it.ge.itmax) goto 1001
       if (torque(1)**2+torque(2)**2+torque(3)**2.lt.1.d-22) goto 1000


        tnorm=sqrt(torque(1)**2+torque(2)**2+torque(3)**2)
        angle=-acos(dot_product(torque,xaxis)/tnorm)
        call cross_product(xaxis,torque,axis)
        axisnorm=sqrt(axis(1)**2+axis(2)**2+axis(3)**2)
        rotmat=unitmat
        if(axisnorm.gt.0.d0) then
        axis=axis/axisnorm
        call rotation(rotmat,angle,axis)
        endif
        rotmatall=matmul(rotmat,rotmatall)
        sy=0.d0 ; sz=0.d0
        do i=1,3
          latvec(:,i)=matmul(rotmat,latvec(:,i))
          vlat(:,i)=matmul(rotmat,vlat(:,i))
          sy=sy+latvec(2,i)**2
          sz=sz+latvec(3,i)**2
        enddo
         cx=torque(1)/(sz+sy)*0.9d0
         do i=1,3
         vlat(2,i)=vlat(2,i)+cx*latvec(3,i)
         vlat(3,i)=vlat(3,i)-cx*latvec(2,i)
         enddo
enddo
1001  write(100,'(a,3(e11.3))') 'WARNING REMAINING TORQUE',torque

1000  call invertmat(rotmatall,rotmat,3)
      do i=1,3
      vlat(:,i)=matmul(rotmat,vlat(:,i))
      enddo
      return
end subroutine elim_torque_cell


!************************************************************************************

subroutine init_vel(parini,vel,vel_lat,vel_vol,latvec,pos_red,latmass,temp,nsoften,folder)
 use mod_interface
 use global, only: target_pressure_habohr,target_pressure_gpa,nat,ntypat,znucl
 use global, only: amu,amutmp,typat,char_type,fixat,fixlat,mol_soften,bc
 use defs_basis
 use mod_parini, only: typ_parini
implicit none
 type(typ_parini), intent(in):: parini
 real(8) :: acell_in(3),xred_in(3,nat),vel_in(3,nat),fcart_in(3,nat),etot_in,strten_in(6),vel_lat_in(3,3)
!*********************************************************************************************
 real(8):: vel(3,nat),temp,pos_red(3,nat),vcm(3),vel_vol
 integer:: i,iat,idim,nsoften
 real(8):: amass(nat),s1,s2,v2gauss,vtest,rescale_vel,vel_lat(3,3),latvec(3,3),latmass
! real(8), parameter :: Ha_eV=27.21138386d0 ! 1 Hartree, in eV
! real(8), parameter :: kb_HaK=8.617343d-5/Ha_eV ! Boltzmann constant in Ha/K
! real(8), parameter :: amu_emass=1.660538782d-27/9.10938215d-31 ! 1 atomic mass unit, in electronic mass
 real(8), parameter :: temp_fac_lat=1.d-1 !This percentage of the temperature that should be given to the lattice 
 real(8):: pressure,curv0,curv,res,count_soft
 character(40):: folder

  vel=0.d0
  vel_vol=0.d0
  vel_lat=0.d0

  write(*,*) "# Initializing velocities"
  write(*,*) "# PARAMETERS: Ha_eV,kb_HaK,amu_emass -->",Ha_eV,kb_HaK,amu_emass
  write(*,*) "# Temp",temp

!Assign masses to each atom (for MD)
 do iat=1,nat
 amass(iat)=amu_emass*amu(typat(iat))
  if(parini%verb.gt.0) write(*,'(a,i5,1x,es15.7,1x,es15.7)') " # iat, AMU, EM: ",iat,amu(typat(iat)),amass(iat)
 end do


pressure=target_pressure_habohr

if(mol_soften) then
         if(any(fixat).or.any(fixlat)) stop "Fixed atoms or cell not yet implemented for molecular softening"
!Init atomic velocities
         call init_rotvels(nat,pos_red,latvec,temp,amass,vel)
         if(parini%verb.gt.0) then
         do iat=1,nat
            write(*,'(a,i5,3(1x,es15.7))') " # iat, VEL: ",iat ,vel(:,iat)
         end do
         endif
!Now get some cell velocity
         call gausdist_cell(latvec,vel_lat)
!Soften the velocities of lattice
         if(.not.fixlat(7)) then
         call soften_lat(parini,latvec,pos_red,vel_lat,curv0,curv,res,pressure,count_soft,amass,nsoften,folder)
         call elim_fixed_lat(latvec,vel_lat)
         endif
else
!Get random Gaussian distributed atomic velocities
         if(.not.(all(fixat))) then
           call gausdist(nat,vel,amass)
!Soften the velocities of atoms
           call soften_pos(parini,latvec,pos_red,vel,curv0,curv,res,pressure,count_soft,amass,nsoften,folder)
!Get rid of center-of-mass velocity, taking also into account fixed atoms (This has already been done in gausdist for all free atoms, but lets do it again...)
           if(.not.any(fixat(:))) then
             s1=sum(amass(:))
             do idim=1,3
               s2=sum(amass(:)*vel(idim,:))
               vel(idim,:)=vel(idim,:)-s2/s1
             end do
           else
             if(trim(parini%potential_potential)=="lenosky_tb_lj") then
                write(*,'(a)') " Eliminating LJ atom velocities"  
                do iat=1,nat
                  if(int(znucl(typat(iat))).gt.200) vel(:,iat)=0.d0
                enddo
             endif
             
             write(*,'(a)') " Eliminating fixed atom velocities"  
             s1=0.d0
             vcm=0.d0
             do iat=1,nat
                if(.not.fixat(iat)) then
                  s1=s1+amass(iat)
                  vcm(:)=vcm(:)+amass(iat)*vel(:,iat)
                endif
             enddo
             if (.not. (nat.gt.1 .and. (nat-count(fixat))==1)) then !Dont eliminate center of mass if there is only one atom to move
             do iat=1,nat
                vel(:,iat)=vel(:,iat)-vcm(:)/s1
             enddo
             endif
             call elim_fixed_at(nat,vel)
           endif
!Recompute v2gauss
           v2gauss=0.d0
           vtest=0.d0
           do iat=1,nat
             do idim=1,3
               v2gauss=v2gauss+vel(idim,iat)*vel(idim,iat)*amass(iat)
               vtest=vtest+vel(idim,iat)/(3.d0*nat)
             end do
           end do
!Now rescale the velocities to give the exact temperature
         rescale_vel=sqrt(3.d0*nat*kb_HaK*temp/v2gauss)
         vel(:,:)=vel(:,:)*rescale_vel
         if(parini%verb.gt.0) then
           do iat=1,nat
              write(*,'(a,i5,3(1x,es15.7))') " # iat, VEL: ",iat ,vel(:,iat)
           end do
         endif
         endif

!Now get some cell velocity
         if(.not.(all(fixlat(1:6))).and.(.not.fixlat(7)).and.bc.ne.2) then
           call gausdist_cell(latvec,vel_lat)
!Soften the velocities of lattice
           call soften_lat(parini,latvec,pos_red,vel_lat,curv0,curv,res,pressure,count_soft,amass,nsoften,folder)
           call elim_fixed_lat(latvec,vel_lat)
         endif
endif



!Rescale to get the correct "temperature" for the cell. This is chosen by a factor "temp_fac_lat"
if(fixlat(7)) then !Andersen MD
         call random_number(vel_vol)
         vel_vol=(vel_vol-0.5d0)*2.d0
         vel_vol=vel_vol*sqrt(kb_HaK*temp*temp_fac_lat)
elseif(.not.all(fixlat(1:6)).and.bc.ne.2) then
!        Recompute v2gauss
         v2gauss=0.d0
         vtest=0.d0
         do i=1,3
           do idim=1,3
             v2gauss=v2gauss+vel_lat(idim,i)*vel_lat(idim,i)*latmass
             vtest=vtest+vel_lat(idim,i)/(3.d0*3.d0)
           end do
         end do
!        Now rescale the velocities to give the exact temperature*temp_fac_lat
         rescale_vel=sqrt(3.d0*3.d0*kb_HaK*temp*temp_fac_lat/v2gauss)
         vel_lat(:,:)=vel_lat(:,:)*rescale_vel
         do i=1,3
           write(*,'(a,i5,3(1x,es15.7))') " # lat, VEL: ",i,vel_lat(:,i)
         end do
endif
end subroutine init_vel

!************************************************************************************

        subroutine soften_pos(parini,latvec,pos_red0,ddcart,curv0,curv,res,pressure,count_soft,amass,nsoft,folder)
 use mod_interface, except_this_one=>norm
 use global, only: target_pressure_habohr,target_pressure_gpa,nat,ntypat,znucl,amu,amutmp,typat
 use global, only: char_type,alpha_at,units,usewf_soften,auto_soft,fixat,fixlat
 use defs_basis
 use interface_code
 use modsocket, only: sock_extra_string
 use mod_parini, only: typ_parini
implicit none
 type(typ_parini), intent(in):: parini
 real(8) :: acell_in(3),xred_in(3,nat),vel_in(3,nat),etot_in,strten_in(6),vel_lat_in(3,3)
!*******************************************************************
        integer:: nsoft,i,it,nit,iprec,iat
        real(8):: curv0,curv,res,pressure,count_soft,alpha
        real(8):: rxyz(3*nat)
        real(8):: latvec(9),latvec_in(9)
        real(8):: ddcart(3*nat)
        real(8):: rxyzcart(3*nat)
        real(8):: flat(9)
        real(8):: pos_red0(3*nat)
        real(8):: pos_red_in(3*nat)
        real(8):: amass(nat)
        real(8):: wlat(9),wlatold(9),fxyzcart(3*nat)
        real(8), allocatable :: wpos(:),fxyz(:)
        real(8):: eps_dd
        real(8):: etot
        real(8):: etot0
        real(8):: fd2
        real(8):: res1 
        real(8):: res2
        real(8):: sdd
        real(8):: sdf
        real(8):: tt
        logical:: getwfk
        character(40):: filename,folder
        character(4):: fn4         
        real(8):: pos_prev(3*nat),dir_prev(3*nat),dir(3*nat),angle,norm
        logical:: decrease
        write(*,'(a,i5)')" # Entering SOFTENING routine for ATOMS, nsoft= ",nsoft 
        if(auto_soft) write(*,'(a)')" # Automatic softening activated" 
        decrease=.false.
!        rxyzcart=rxyz        
!First transform the atomic positions from internal to external coordinates
        call rxyz_int2cart(latvec,pos_red0,rxyz,nat)

        nit=nsoft
        eps_dd=5.d-1
!        alpha=3.d0    ! step size for  Si
!        alpha=1.d0    ! step size for  C 
        alpha=alpha_at
        allocate(wpos(3*nat),fxyz(3*nat))

!         call  rxyz_int2cart(latvec,rxyz,rxyzcart,nat)
!         call  energyandforces(nat,latvec,rxyzcart,fxyzcart,flat,pressure,etot0,count_soft)                                    !
call rxyz_cart2int(latvec,pos_red_in,rxyz,nat)
latvec_in=latvec
iprec=1;getwfk=.false.
       write(fn4,'(i4.4)') 0
       sock_extra_string="SOFTAT"//trim(fn4)
call get_energyandforces_single(parini,latvec_in,pos_red_in,fxyz,strten_in,etot_in,iprec,getwfk)
call get_enthalpy(latvec_in,etot_in,pressure,etot0)

!         call  fxyz_cart2int(nat,fxyzcart,fxyz,latvec)

! normalize initial guess
        sdd=0.d0
        call elim_fixed_at(nat,ddcart)
        do i=1,3*nat
        sdd=sdd+ddcart(i)**2
        enddo
        sdd=eps_dd/sqrt(sdd)
        do i=1,3*nat
        ddcart(i)=sdd*ddcart(i)
        enddo

 do it=1,nit
        do i=1,3*nat
        wpos(i)=rxyz(i)+ddcart(i)
        enddo

!Here we check for direction variation during softening
        dir_prev=dir
        dir=wpos-pos_prev
        if(it.ge.3) then
           do iat=1,nat  
              angle=dot_product(dir((iat-1)*3+1:(iat)*3),dir_prev((iat-1)*3+1:(iat)*3))
              norm=sqrt(dot_product(dir((iat-1)*3+1:(iat)*3),dir((iat-1)*3+1:(iat)*3)))+&
                   sqrt(dot_product(dir_prev((iat-1)*3+1:(iat)*3),dir_prev((iat-1)*3+1:(iat)*3)))
              angle=angle/norm
              if(angle.lt.0.d0) decrease=.true.
           enddo
        endif

        pos_prev=wpos


call rxyz_cart2int(latvec,pos_red_in,wpos,nat)
latvec_in=latvec
iprec=1
if(usewf_soften) then
    getwfk=.true.
else
    getwfk=.false.
endif
       write(fn4,'(i4.4)') it
       sock_extra_string="SOFTAT"//trim(fn4)
call get_energyandforces_single(parini,latvec_in,pos_red_in,fxyz,strten_in,etot_in,iprec,getwfk)
call get_enthalpy(latvec_in,etot_in,pressure,etot)

        fd2=2.d0*(etot-etot0)/eps_dd**2

        sdf=0.d0
        sdd=0.d0
        do i=1,3*nat
        sdf=sdf+ddcart(i)*fxyz(i)
        sdd=sdd+ddcart(i)*ddcart(i)
        enddo

        curv=-sdf/sdd
        if (it.eq.1) curv0=curv
        res=0.d0
        tt=0.d0

        do i=1,3*nat
        tt=tt+fxyz(i)**2
        fxyz(i)=fxyz(i)+curv*ddcart(i)
        res=res+fxyz(i)**2
        enddo
        
if(parini%verb.gt.0) write(*,'(a,(e13.5),i5)') ' # SOFTEN: ',res,it
        res=sqrt(res)
        tt=sqrt(tt)
if(it==1) write(*,'(a,i5,4(e13.5),e18.10)')' # SOFTEN: init atomic  it,fnrm,res,curv,fd2,etot ',it,tt,res,curv,fd2,etot-etot0
if(parini%verb.gt.0) then
        write(*,'(a,i5,4(e13.5),e18.10)') ' # SOFTEN: it,fnrm,res,curv,fd2,etot ',it,tt,res,curv,fd2,etot-etot0
        write(fn4,'(i4.4)') it 
        filename=trim(folder)//'possoft_at'//fn4//'.ascii'
        call write_atomic_file_ascii(parini,filename,nat,units,pos_red_in,latvec_in,fxyz,strten_in,char_type,&
             &ntypat,typat,fixat,fixlat,etot_in,pressure,curv,res)
endif
        call elim_fixed_at(nat,fxyz)
        do i=1,3*nat
        wpos(i)=wpos(i)+alpha*fxyz(i)
        enddo

        do i=1,3*nat
        ddcart(i)=wpos(i)-rxyz(i)
        enddo
         
        call elim_moment(nat,ddcart(1:3*nat),amass)
        call elim_fixed_at(nat,ddcart)
                         
        sdd=0.d0
        do i=1,3*nat
        sdd=sdd+ddcart(i)*ddcart(i)
        enddo

!        if (res.le.curv*eps_dd*5.d-1) goto 1000
        sdd=eps_dd/sqrt(sdd)
        do i=1,3*nat
        ddcart(i)=ddcart(i)*sdd
        enddo
      enddo

!       write(*,*) '# No convergence in low_cur_dir',res
1000   continue
      
write(*,'(a,i5,4(e13.5),e18.10)')' # SOFTEN: final atomic it,fnrm,res,curv,fd2,etot ',it,tt,res,curv,fd2,etot-etot0
!Decrease the stepsize if necessary
       if(auto_soft.and.nsoft.lt.3) write(*,'(a,e18.10)') ' # SOFTEN: increase nsoft for auto adjustment'
       if(auto_soft.and.nsoft.ge.3) then
          if(decrease) then 
             alpha_at=alpha_at/1.1d0
          else
             alpha_at=alpha_at*1.1d0
          endif
       write(*,'(a,e18.10)') ' # SOFTEN: new alpha_at :  ',alpha_at
       endif


       deallocate(wpos,fxyz)
    end subroutine

!************************************************************************************

        subroutine soften_lat(parini,latvec,pos_red0,ddlat,curv0,curv,res,pressure,count_soft,amass,nsoft,folder)
 use mod_interface, except_this_one=>norm
 use global, only: target_pressure_habohr,target_pressure_gpa,nat,ntypat,znucl,amu,amutmp,typat
 use global, only: char_type,alpha_at,alpha_lat,units,usewf_soften,auto_soft,fixat,fixlat
 use defs_basis
 use interface_code
 use modsocket, only: sock_extra_string
 use mod_parini, only: typ_parini
 implicit none
 type(typ_parini), intent(in):: parini
 real(8) :: acell_in(3),xred_in(3,nat),vel_in(3,nat),fcart_in(3,nat),etot_in,strten_in(6),vel_lat_in(3,3)
!*******************************************************************
        integer:: nsoft,i,it,nit,iprec,iat
        real(8):: curv0,curv,res,pressure,count_soft,alpha,alphalat
        real(8):: rxyz(3*nat)
        real(8):: latvec(9),latvec_in(9)
        real(8):: dd(3*nat)
        real(8):: ddcart(3*nat)
        real(8):: rxyzcart(3*nat)
        real(8):: ddlat(9)
        real(8):: ddall(3*nat+9)
        real(8):: flat(9)
        real(8):: pos_red0(3*nat)
        real(8):: amass(nat)
        real(8):: wlat(9),wlatold(9),fxyzcart(3*nat)
        real(8), allocatable :: wpos(:),fxyz(:)
        real(8):: eps_dd
        real(8):: etot
        real(8):: etot0
        real(8):: fd2
        real(8):: res1 
        real(8):: res2
        real(8):: sdd
        real(8):: sdf
        real(8):: tt
        real(8):: vol
        character(40):: filename,folder
        character(4):: fn4         
        logical:: getwfk
        real(8):: lat_prev(3*3),dir_prev(3*3),dir(3*3),angle,norm
        logical:: decrease
        decrease=.false.

        write(*,'(a,i5)')" # Entering SOFTENING routine for LATTICE, nsoft= ",nsoft
        if(auto_soft) write(*,'(a)')" # Automatic softening activated"

!        ddcart=dd
        rxyz=pos_red0
!        rxyzcart=rxyz        

        nit=nsoft
!        eps_dd=1.d0
        eps_dd=1.d0
!        alphalat=1.d0 !step size for Lattice
        alphalat=alpha_lat

latvec_in=latvec
iprec=1;getwfk=.false.
       write(fn4,'(i4.4)') 0
       sock_extra_string="SOFTLAT"//trim(fn4)
call get_energyandforces_single(parini,latvec_in,rxyz,fcart_in,strten_in,etot_in,iprec,getwfk)
call strten2flat(strten_in,flat,latvec,pressure)
call get_enthalpy(latvec_in,etot_in,pressure,etot0)
!Multiply the force with the unit cell volume to get the correct force
!call getvol(latvec,vol)
!flat=flat*vol

!         call  fxyz_cart2int(nat,fxyzcart,fxyz,latvec)

! normalize initial guess
9909 continue
        sdd=0.d0
        do i=1,9
        sdd=sdd+ddlat(i)**2
        enddo
        sdd=eps_dd/sqrt(sdd)
        do i=1,9
        ddlat(i)=sdd*ddlat(i)
        enddo
!Check if this eps_dd will invert the cell
        do i=1,9
        wlat(i)=latvec(i)+ddlat(i)
        enddo
        call getvol(wlat,vol)
        if(vol.le.0.d0) then
              eps_dd=eps_dd*0.5d0
              write(*,*) "Epsdd set to ",eps_dd,vol
              call getvol(latvec,vol)
              
              goto 9909
        endif


 do it=1,nit
        do i=1,9
        wlat(i)=latvec(i)+ddlat(i)
        enddo

!Here we check for direction variation during softening
        dir_prev=dir
        dir=wlat-lat_prev
        if(it.ge.3) then
           do iat=1,3  
              angle=dot_product(dir((iat-1)*3+1:(iat)*3),dir_prev((iat-1)*3+1:(iat)*3))
              norm=sqrt(dot_product(dir((iat-1)*3+1:(iat)*3),dir((iat-1)*3+1:(iat)*3)))+&
                   sqrt(dot_product(dir_prev((iat-1)*3+1:(iat)*3),dir_prev((iat-1)*3+1:(iat)*3)))
              angle=angle/norm
              if(angle.lt.0.d0) decrease=.true.
           enddo
        endif

        lat_prev=wlat

latvec_in=wlat
iprec=1
if(usewf_soften) then
    getwfk=.true.
else
    getwfk=.false.
endif
       write(fn4,'(i4.4)') it
       sock_extra_string="SOFTLAT"//trim(fn4)
call get_energyandforces_single(parini,latvec_in,rxyz,fcart_in,strten_in,etot_in,iprec,getwfk)
call strten2flat(strten_in,flat,wlat,pressure)
call get_enthalpy(latvec_in,etot_in,pressure,etot)
!Multiply the force with the unit cell volume to get the correct force
!call getvol(wlat,vol)
!flat=flat*vol

        fd2=2.d0*(etot-etot0)/eps_dd**2

        sdf=0.d0
        sdd=0.d0
        do i=1,9
        sdf=sdf+ddlat(i)*flat(i)
        sdd=sdd+ddlat(i)*ddlat(i)
        enddo

        curv=-sdf/sdd
        if (it.eq.1) curv0=curv
        res=0.d0
        res1=0.d0
        res2=0.d0
        tt=0.d0

        do i=1,9
        tt=tt+flat(i)**2
        flat(i)=flat(i)+curv*ddlat(i)
        res=res+flat(i)**2
        enddo

if(parini%verb.gt.0) write(*,'(a,(e13.5),i5)') ' # SOFTEN: ',res,it
        res=sqrt(res)
        tt=sqrt(tt)
if(it==1) write(*,'(a,i5,4(e13.5),e18.10)')' # SOFTEN: init lattice  it,fnrm,res,curv,fd2,etot ',it,tt,res,curv,fd2,etot-etot0
if(parini%verb.gt.0) then
        write(*,'(a,i5,4(e13.5),e18.10)') ' # SOFTEN: it,fnrm,res,curv,fd2,etot ',it,tt,res,curv,fd2,etot-etot0
        write(fn4,'(i4.4)') it 
        filename=trim(folder)//'possoft_lat'//fn4//'.ascii'
        call write_atomic_file_ascii(parini,filename,nat,units,rxyz,latvec_in,fcart_in,strten_in,char_type,ntypat,&
        &typat,fixat,fixlat,etot_in,pressure,curv,res)
endif

        wlatold=wlat
        do i=1,9
        wlat(i)=wlat(i)+alphalat*flat(i)
!        write(*,*) " # FLAT",flat(i)
        enddo

        do i=1,9
        ddlat(i)=wlat(i)-latvec(i)
        enddo
         
        call elim_torque_cell(latvec,ddlat)
                         
        sdd=0.d0
        do i=1,9
        sdd=sdd+ddlat(i)*ddlat(i)
        enddo

!        if (res.le.curv*eps_dd*1.d-2) goto 1000
        sdd=eps_dd/sqrt(sdd)
        do i=1,9
        ddlat(i)=ddlat(i)*sdd
        enddo
      enddo

!       write(*,*) '# No convergence in low_cur_dir',res
1000   continue

write(*,'(a,i5,4(e13.5),e18.10)')' # SOFTEN: final lattice it,fnrm,res,curv,fd2,etot ',it,tt,res,curv,fd2,etot-etot0
!Decrease the stepsize if necessary
       if(auto_soft.and.nsoft.lt.3) write(*,'(a,e18.10)') ' # SOFTEN: increase nsoft for auto adjustment'
       if(auto_soft.and.nsoft.ge.3) then
          if(decrease) then
             alpha_lat=alpha_lat/1.1d0
          else
             alpha_lat=alpha_lat*1.1d0
          endif
       write(*,'(a,es18.10,1x,es18.10)') ' # SOFTEN: new alpha_at, alpha_lat : ',alpha_at,alpha_lat
       endif


    end subroutine

!************************************************************************************

 subroutine fxyz_cart2int(nat,fxyz_cart,fxyz_int,latvec)
 use mod_interface
 !This subrtouine will transform theforces initially in the cartesian system into the internal coordinates with respect to the
 !cell vectors provided in the latvec
 implicit none
 real(8):: fxyz_cart(3,nat),fxyz_int(3,nat),latvec(3,3),transmat(3,3)
 integer:: nat,iat
 transmat(1,:)=latvec(:,1)
 transmat(2,:)=latvec(:,2)
 transmat(3,:)=latvec(:,3)
 do iat=1,nat
 fxyz_int(:,iat)=matmul(transmat,fxyz_cart(:,iat))
 enddo
 end subroutine fxyz_cart2int 

!************************************************************************************

 subroutine strten2flat(strten,flat,latvec,press)
 use mod_interface
 !flat is the force on the lettice vector per unit cell volume
 implicit none
 real(8):: strten(6),flat(3,3),latvec(3,3),press,pressmat(3,3),str_matrix(3,3),latvect(3,3),latvectinv(3,3),vol
 pressmat=0.d0
 pressmat(1,1)=press
 pressmat(2,2)=press
 pressmat(3,3)=press
           str_matrix(1,1)=strten(1)
           str_matrix(2,2)=strten(2)
           str_matrix(3,3)=strten(3)
           str_matrix(1,2)=strten(6)
           str_matrix(2,1)=strten(6)
           str_matrix(1,3)=strten(5)
           str_matrix(3,1)=strten(5)
           str_matrix(2,3)=strten(4)
           str_matrix(3,2)=strten(4)
           flat=str_matrix+pressmat
           flat=-flat
 latvect(:,1)=latvec(1,:)
 latvect(:,2)=latvec(2,:)
 latvect(:,3)=latvec(3,:)
 call invertmat(latvect,latvectinv,3)
 flat=matmul(flat,latvectinv)
! I guess i forgot to divide through volume
 call getvol(latvec,vol)
 flat=flat*vol
 end subroutine

!*************************************************************************

 subroutine backtocell(nat,latvec,rxyz_red)
 use mod_interface
 !This subroutine will transform back all atoms into the periodic cell
 !defined by the 3 lattice vectors in latvec=[v1.v2.v3]
 implicit none
 integer:: nat,i,iat,j
 real(8) :: latvec(3,3), rxyz_red(3,nat), rxyz(3,nat), crossp(3),a(3),b(3), nvec(3,3), dist(6),eps,count
 real(8) :: v(3,3),vol
 logical:: neccesary
!First check if the volume is positive
 call getvol(latvec,vol)
 if(vol.le.0.d0) stop "Negative volume during backtocell"
 do iat=1,nat
        rxyz_red(1,iat)=modulo(modulo(rxyz_red(1,iat),1.d0),1.d0)
        rxyz_red(2,iat)=modulo(modulo(rxyz_red(2,iat),1.d0),1.d0)
        rxyz_red(3,iat)=modulo(modulo(rxyz_red(3,iat),1.d0),1.d0)
 enddo
! v=latvec
! vol= v(1,1)*v(2,2)*v(3,3)-v(1,1)*v(2,3)*v(3,2)-v(1,2)*v(2,1)*v(3,3)+&
!      v(1,2)*v(2,3)*v(3,1)+v(1,3)*v(2,1)*v(3,2)-v(1,3)*v(2,2)*v(3,1)
! if(vol.le.0.d0) stop "Negative volume during backtocell"
! call rxyz_int2cart(latvec,rxyz_red,rxyz,nat)
!
! !To really be on the safe side, the translation vector can be shortened by  a factor eps in order
! !to get the atom into the cell. 
!!  eps=1.d-10
!  eps=1.d-15
! ! eps=1.d0-eps
! count=0.d0
! neccesary=.true.
! do while(neccesary)
! neccesary=.false.
! count=count+1.d0
! !generate 3 normal vectors of the 3 planes
! call nveclatvec(latvec,nvec)
! do iat=1,nat
! !3 planes through origin (xy,yz,zx)
! do i=1,3
! dist(i)=DOT_PRODUCT(rxyz(:,iat),nvec(:,i))
!! if(dist(i).lt.0.d0) then
! if(dist(i).lt.-abs(dist(i))*eps) then
!! write(*,*) "unten 1",i,iat
! rxyz(:,iat)=rxyz(:,iat)+latvec(:,mod(i+1,3)+1)!*eps
! neccesary=.true.
! endif
!
! !3 planes on top/side/back (xy,yz,zx)
! dist(i+3)=DOT_PRODUCT(rxyz(:,iat)-latvec(:,mod(i+1,3)+1),nvec(:,i))
!! if(dist(i+3).gt.0.d0) then
! if(dist(i+3).gt.abs(dist(i+3))*eps) then
!! write(*,*) "unten 1",i,iat
!! if(dist(i+3).gt.eps) then
!! write(*,*) "oben 2",i,iat
! rxyz(:,iat)=rxyz(:,iat)-latvec(:,mod(i+1,3)+1)!*eps
! neccesary=.true.
! endif
! enddo
! enddo
! if(count.gt.1.d6) stop "Too many iterations in back-to-cell"
! enddo
! 
! call rxyz_cart2int(latvec,rxyz_red,rxyz,nat)
 end subroutine

!************************************************************************************

 subroutine nveclatvec(latvec,nvec)
 !Will calculate the normalized normal vector to the 3 planes of the cell
 use mod_interface, except_this_one=>norm
 implicit none
 real*8, intent(in) :: latvec(3,3)
 real*8, intent(out):: nvec(3,3)
 real*8             :: a(3),b(3),crossp(3),norm
 integer:: i
 do i=1,3
 a=latvec(:,i)
 b=latvec(:,mod(i,3)+1)
 call cross_product(a,b,crossp)
 norm=dsqrt(crossp(1)*crossp(1)+crossp(2)*crossp(2)+crossp(3)*crossp(3))
 nvec(:,i)=crossp(:)/norm
 enddo
 end subroutine

!************************************************************************************

subroutine getvol(latvec,vol)
use mod_interface
implicit none
real(8):: latvec(3,3),v(3,3),vol
 v=latvec
 vol= v(1,1)*v(2,2)*v(3,3)-v(1,1)*v(2,3)*v(3,2)-v(1,2)*v(2,1)*v(3,3)+&
      v(1,2)*v(2,3)*v(3,1)+v(1,3)*v(2,1)*v(3,2)-v(1,3)*v(2,2)*v(3,1)
end subroutine

!************************************************************************************

subroutine correct_latvec(latvec,pos_red,nat,correctalg,iout)
use mod_interface
implicit none
integer:: correctalg,nat,iproc,iout
real(8):: latvec(3,3),pos_red(3,nat),latvec0(3,3),diff(9)
latvec0=latvec
select case(correctalg)
 case(1)
  !The oganov method
  iproc=1
  call correct_latvec_oganov(latvec,pos_red,nat,iproc)
 case(2)
  !Niggli reduction
  call fixcell_niggli(nat,latvec,pos_red)
 case default
  stop "Currently not implemented method for cell correction!"
end select 
diff(1:3)=latvec(:,1)-latvec0(:,1)
diff(4:6)=latvec(:,2)-latvec0(:,2)
diff(7:9)=latvec(:,3)-latvec0(:,3)
if(dot_product(diff,diff).lt.1.d-14) then
 write(*,*) "Cell not changed"
 write(*,*) latvec
 write(*,*) latvec0
 iout=0
else
 iout=1
 write(*,*) "Cell changed"
 write(*,*) latvec
 write(*,*) latvec0
endif


end subroutine

!************************************************************************************
 
 subroutine correct_latvec_oganov(latvec,pos_red,nat,iproc)
 use mod_interface, except_this_one=>norm
 !use cell_utils
 !This subroutine will use the algorithm proposed by oganov and glass (J.Phys,Cond.Mat 20,2008) to perform a transformation of the lattice vectors into an equivalent
 !system where the length of all cell vectors are similar (no nasty angles).
 implicit none
 real(8)              :: latvec(3,3),rxyz(3,nat),pos_red(3,nat)  
 logical              :: correct 
 real(8)              :: val_inter,norm_half,norm !The val_... are the absolute value of the projections
 real(8)              :: a(3),b(3),c(3) ! the three latticevectors
 real(8)              :: tempvec(3),v(3,3),vol,sign_inter
 integer              :: i,nat,counter,iproc

 call backtocell(nat,latvec,pos_red)
 call rxyz_int2cart(latvec,pos_red,rxyz,nat)
 a=latvec(:,1);b=latvec(:,2);c=latvec(:,3)
 !check volume
 v=latvec
 vol= v(1,1)*v(2,2)*v(3,3)-v(1,1)*v(2,3)*v(3,2)-v(1,2)*v(2,1)*v(3,3)+&
      v(1,2)*v(2,3)*v(3,1)+v(1,3)*v(2,1)*v(3,2)-v(1,3)*v(2,2)*v(3,1)
! write(*,*) "initial vol",vol
 !Check the 6 criteria ab,ba,ac,ca,bc,cb
 correct=.true.
 counter=0
 do while(correct) 
 counter=counter+1
 correct=.false.
 !ab
 val_inter=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
 sign_inter=sign(1.d0,val_inter)
 norm=b(1)*b(1)+b(2)*b(2)+b(3)*b(3) 
 norm=dsqrt(norm)
 norm_half=norm*0.5d0
 val_inter=abs(val_inter/norm)
 if(val_inter.gt.norm_half) then
! write(*,*) "Transform ab"
 correct=.true.
 val_inter=val_inter/norm
 a(:)=a(:)-ceiling(val_inter)*sign_inter*b(:)
 latvec(:,1)=a(:)
 endif


 !ba
 val_inter=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
 sign_inter=sign(1.d0,val_inter)
 norm=a(1)*a(1)+a(2)*a(2)+a(3)*a(3) 
 norm=dsqrt(norm)
 norm_half=norm*0.5d0
 val_inter=abs(val_inter/norm)
 if(val_inter.gt.norm_half) then
! write(*,*) "Transform ba"
 correct=.true.
 val_inter=val_inter/norm
 b(:)=b(:)-ceiling(val_inter)*sign_inter*a(:)
 latvec(:,2)=b(:)
 endif


 !ac
 val_inter=a(1)*c(1)+a(2)*c(2)+a(3)*c(3)
 sign_inter=sign(1.d0,val_inter)
 norm=c(1)*c(1)+c(2)*c(2)+c(3)*c(3)
 norm=dsqrt(norm)
 norm_half=norm*0.5d0
 val_inter=abs(val_inter/norm)
 if(val_inter.gt.norm_half) then
! write(*,*) "Transform ac"
 correct=.true.
 val_inter=val_inter/norm
 a(:)=a(:)-ceiling(val_inter)*sign_inter*c(:)
 latvec(:,1)=a(:)
 endif


 !ca
 val_inter=a(1)*c(1)+a(2)*c(2)+a(3)*c(3)
 sign_inter=sign(1.d0,val_inter)
 norm=a(1)*a(1)+a(2)*a(2)+a(3)*a(3)
 norm=dsqrt(norm)
 norm_half=norm*0.5d0
 val_inter=abs(val_inter/norm)
 if(val_inter.gt.norm_half)then
! write(*,*) "Transform ca"
 correct=.true.
 val_inter=val_inter/norm
 c(:)=c(:)-ceiling(val_inter)*sign_inter*a(:)
 latvec(:,3)=c(:)
 endif


 !bc
 val_inter=b(1)*c(1)+b(2)*c(2)+b(3)*c(3)
 sign_inter=sign(1.d0,val_inter)
 norm=c(1)*c(1)+c(2)*c(2)+c(3)*c(3)
 norm=dsqrt(norm)
 norm_half=norm*0.5d0
 val_inter=abs(val_inter/norm)
 if(val_inter.gt.norm_half)then
! write(*,*) "Transform bc"
 correct=.true.
 val_inter=val_inter/norm
 b(:)=b(:)-ceiling(val_inter)*sign_inter*c(:)
 latvec(:,2)=b(:)
 endif


 !cb
 val_inter=b(1)*c(1)+b(2)*c(2)+b(3)*c(3)
 sign_inter=sign(1.d0,val_inter)
 norm=b(1)*b(1)+b(2)*b(2)+b(3)*b(3)
 norm=dsqrt(norm)
 norm_half=norm*0.5d0
 val_inter=abs(val_inter/norm)
 if(val_inter.gt.norm_half)then
! write(*,*) "Transform cb"
 correct=.true.
 val_inter=val_inter/norm
 c(:)=c(:)-ceiling(val_inter)*sign_inter*b(:)
 latvec(:,3)=c(:)
 endif
 
 enddo

 if(iproc==0) write(*,*)"# Number of correct cycles" , counter
 v=latvec
 vol= v(1,1)*v(2,2)*v(3,3)-v(1,1)*v(2,3)*v(3,2)-v(1,2)*v(2,1)*v(3,3)+&
      v(1,2)*v(2,3)*v(3,1)+v(1,3)*v(2,1)*v(3,2)-v(1,3)*v(2,2)*v(3,1)
! write(100,*) "final vol",vol
 call backtocell_cart(nat,latvec,rxyz)
 call rxyz_cart2int(latvec,pos_red,rxyz,nat)
 end subroutine

!************************************************************************************

 subroutine backtocell_cart(nat,latvec,rxyz)
 !This subroutine will transform back all atoms into the periodic cell
 !defined by the 3 lattice vectors in latvec=[v1.v2.v3]
 use mod_interface
 implicit none
 integer:: nat,i,iat,j
 real(8) :: latvec(3,3), rxyz(3,nat), crossp(3),a(3),b(3), nvec(3,3), dist(6),eps,count
 real(8) :: v(3,3),vol,rxyz_red(3,nat)
 logical:: neccesary
 call getvol(latvec,vol)
   !First check if the volume is positive
   if(vol.le.0.d0) stop "Negative volume during backtocell"
   call rxyz_cart2int(latvec,rxyz_red,rxyz,nat)
   do iat=1,nat
       rxyz_red(1,iat)=modulo(modulo(rxyz_red(1,iat),1.d0),1.d0)
       rxyz_red(2,iat)=modulo(modulo(rxyz_red(2,iat),1.d0),1.d0)
       rxyz_red(3,iat)=modulo(modulo(rxyz_red(3,iat),1.d0),1.d0)
   enddo
   call rxyz_int2cart(latvec,rxyz_red,rxyz,nat)
!!First check if the volume is positive
! v=latvec
! vol= v(1,1)*v(2,2)*v(3,3)-v(1,1)*v(2,3)*v(3,2)-v(1,2)*v(2,1)*v(3,3)+&
!      v(1,2)*v(2,3)*v(3,1)+v(1,3)*v(2,1)*v(3,2)-v(1,3)*v(2,2)*v(3,1)
! if(vol.le.0.d0) stop "Negative volume during backtocell"
!
! !To really be on the safe side, the translation vector can be shortened by  a factor eps in order
! !to get the atom into the cell. 
!!  eps=1.d-10
!  eps=1.d-15
! ! eps=1.d0-eps
! count=0.d0
! neccesary=.true.
! do while(neccesary)
! neccesary=.false.
! count=count+1.d0
! !generate 3 normal vectors of the 3 planes
! call nveclatvec(latvec,nvec)
! do iat=1,nat
! !3 planes through origin (xy,yz,zx)
! do i=1,3
! dist(i)=DOT_PRODUCT(rxyz(:,iat),nvec(:,i))
!! if(dist(i).lt.0.d0) then
! if(dist(i).lt.-abs(dist(i))*eps) then
!! write(*,*) "unten 1",i,iat
! rxyz(:,iat)=rxyz(:,iat)+latvec(:,mod(i+1,3)+1)!*eps
! neccesary=.true.
! endif
!
! !3 planes on top/side/back (xy,yz,zx)
! dist(i+3)=DOT_PRODUCT(rxyz(:,iat)-latvec(:,mod(i+1,3)+1),nvec(:,i))
!! if(dist(i+3).gt.0.d0) then
! if(dist(i+3).gt.abs(dist(i+3))*eps) then
!! write(*,*) "unten 1",i,iat
!! if(dist(i+3).gt.eps) then
!! write(*,*) "oben 2",i,iat
! rxyz(:,iat)=rxyz(:,iat)-latvec(:,mod(i+1,3)+1)!*eps
! neccesary=.true.
! endif
! enddo
! enddo
! if(count.gt.1.d6) stop "Too many iterations in back-to-cell"
! enddo
 end subroutine

!************************************************************************************

!!  subroutine read_params()
!!  use mod_interface
!!  use defs_basis
!!  use mod_fire,   only:dtmin, dtmax
!!  use minpar, only:parmin_bfgs
!!  use global, only: target_pressure_habohr,target_pressure_gpa,nat,ntypat,znucl,amu,amutmp,typat,ntime_md,char_type
!!  use global, only: nsoften,alpha_at,alpha_lat,ntime_geopt,parres%mdmin,dtion_fire,dtion_md,tolmxf,strfact,dtion_fire_min
!!  use global, only: dtion_fire_max,ka,kb,kc,dkpt1,dkpt2,usewf_geopt,usewf_soften,usewf_md,geopt_method,alphax_at
!!  use global, only: alphax_lat,findsym,finddos,auto_soft,mdmin_max,mdmin_min,auto_mdmin,parini%md_algo,parini%md_integrator,auto_dtion_md
!!  use global, only: parini%nit_per_min,fixat,fixlat,rcov,mol_soften,fragarr,code,auto_kpt
!!  implicit none
!!  integer:: itype,n
!!  real(8):: tmp_val
!!  character(250):: all_line
!!  character(4):: tmp_ch
!!  !The convention within abinit are followed:
!!  !First atomic data, then MD parameters, and then MD parameters
!!  !Example:
!!  ! 120      target_pressure_gpa
!!  ! 4        nat     
!!  ! 2        ntypat    
!!  ! 14    1  znucl     
!!  ! 28    1  amu     
!!  ! 1 2 2 2  typat  
!!  ! 200      ntime_md
!!  ! 200      ntime_geopt
!!  ! 2.d0     bmass
!!  ! Auto 6 12       parres%mdmin or mdmin_min mdmin_max
!!  ! 8 3.0 2.0 F nsoften, alpha_at, alpha_lat,mol_soften
!!  ! FIRE     !Method of geopt: FIRE or BFGS
!!  ! 20.d0 1.d0 50.d0    dtion_fire, min, max .or. alphax_at, alphax_lat 
!!  ! 1.d-4    tolmxf   
!!  ! 100      strfact 
!!  ! .true. .false. .false. usewf_geopt,usewf_soften,usewf_md
!!  ! .true. .false. !Determines if symmetry (isotropy)/FDOS calculation should be used
!!  ! Auto dfine dcourse or ka kb kc    kpt mesh 
!!   open(unit=12,file="params.in") 
!!   read(12,*) target_pressure_gpa   !Target pressure in GPA
!!   target_pressure_habohr=target_pressure_gpa/HaBohr3_GPA
!!   write(*,'(a,2(1x,es15.7))') " # Target pressure for the run in GPa and Ha/Bohr**3: ",target_pressure_gpa,target_pressure_habohr
!!   read(12,*) nat                   !Number of atoms
!!   read(12,*) ntypat                !Number of atom types
!!   if(.not.allocated(znucl)) allocate(znucl(ntypat))
!!   if(.not.allocated(char_type)) allocate(char_type(ntypat))
!!   if(.not.allocated(amu)) allocate(amu(ntypat))
!!   if(.not.allocated(amutmp)) allocate(amutmp(ntypat))
!!   if(.not.allocated(rcov)) allocate(rcov(ntypat))
!!   read(12,*) znucl(1:ntypat)       !Nuclei charge 
!!   read(12,*) amutmp(1:ntypat)      !Atomic mass used for MD and Fire, if all 0, then automatic determination
!!   if(.not.allocated(typat)) allocate(typat(nat))          
!!   if(.not.allocated(fixat))then
!!      allocate(fixat(nat))
!!      fixat(:)=.false.
!!   endif
!!   if(.not.allocated(fragarr))then
!!      allocate(fragarr(nat))
!!      fragarr(:)=-1
!!   endif
!!   read(12,*) typat(1:nat)          !Types of atoms  
!!   read(12,*) ntime_md, parini%md_algo, parini%md_integrator          !Maximum number of iterations during MD, parini%md_algo,
!parini%md_integrator
!!   if(parini%md_algo==4) fixlat(7)=.true.;if(fixlat(7)) parini%md_algo=4
!!   read(12,*) ntime_geopt           !Maximum number of iterations during GEOPT
!!   read(12,*) bmass                 !Cell mass during MD and FIRE
!!  ! read(12,*) parres%mdmin                 !Number of enthalpy minima crossed unit stop MD
!!  !Block parres%mdmin****************
!!   read(12,'(a250)') all_line
!!   n = len_trim(all_line)
!!   read(all_line(1:n),*) tmp_ch
!!    if(trim(tmp_ch)=="Auto") then
!!  !   if(.not.parini%auto_mdmin) then
!!       read(all_line(1:n),*) tmp_ch,mdmin_min,mdmin_max  !Number of enthalpy minima crossed unit stop MD
!!       if(.not.parini%auto_mdmin) parres%mdmin=mdmin_min
!!  !   endif
!!     parini%auto_mdmin=.true.
!!    else
!!     read(all_line(1:n),*)  parres%mdmin  !Number of enthalpy minima crossed unit stop MD
!!     parini%auto_mdmin=.false.
!!    endif
!!  !Block parres%mdmin****************
!!  !Block soften****************
!!   read(12,'(a250)') all_line
!!   n = len_trim(all_line)
!!   read(all_line(1:n),*) tmp_ch
!!    if(trim(tmp_ch)=="Auto") then
!!     if(.not.auto_soft.and.alpha_at.lt.0.d0.and.alpha_lat.lt.0.d0) then
!!       read(all_line(1:n),*) tmp_ch,nsoften, alpha_at, alpha_lat,mol_soften  !Number of softening steps, softening stepsize for atoms and lattice
!!     else
!!       read(all_line(1:n),*) tmp_ch,nsoften,tmp_val,tmp_val,mol_soften
!!     endif
!!     auto_soft=.true.
!!    else
!!     read(all_line(1:n),*) nsoften, alpha_at, alpha_lat,mol_soften   !Number of softening steps, softening stepsize for atoms and lattice
!!     auto_soft=.false.
!!    endif
!!  !Block soften****************
!!  !Block MD timestep***********
!!   read(12,'(a250)') all_line
!!   n = len_trim(all_line)
!!   read(all_line(1:n),*) tmp_ch
!!    if(trim(tmp_ch)=="Auto") then
!!     if(.not.parini%auto_dtion_md) then
!!       read(all_line(1:n),*) tmp_ch, dtion_md,parini%nit_per_min    !(Auto) MD timestep, target number of iterations per minimum
!!     else
!!       read(all_line(1:n),*) tmp_ch, tmp_val ,parini%nit_per_min    !(Auto) MD timestep, target number of iterations per minimum
!!     endif
!!     parini%auto_dtion_md=.true.
!!    else
!!     read(all_line(1:n),*) dtion_md              !(Auto) MD timestep
!!     parini%auto_dtion_md=.false.
!!    endif
!!  !Block MD timestep***********
!!   read(12,*) geopt_method          !Either FIRE or BFGS
!!   if    (trim(geopt_method)=="FIRE") then 
!!     read(12,*) dtion_fire,dtion_fire_min,dtion_fire_max        !Initial timestep for FIRE, Min, Max
!!   elseif(trim(geopt_method)=="MBFGS".or.trim(geopt_method)=="RBFGS") then
!!     read(12,*) alphax_at,alphax_lat    !Stepsize for atoms and for lattice in BFGS
!!   else
!!     stop "Wrong options for geometry optimizer, only FIRE or RBFGS or MBFGS accepted"
!!   endif
!!  !Copy parameters of fire to the fire module
!!      dtmin=dtion_fire_min
!!      dtmax=dtion_fire_max
!!  !Copy parameters of bfgs to the bfgs module
!!      parmin_bfgs%betax=alphax_at
!!      parmin_bfgs%betax_lat=alphax_lat
!!   read(12,*) tolmxf                !Force tolerance for GEOPT convergance 
!!   read(12,*) strfact               !Factor to multiply stress 
!!   read(12,*) usewf_geopt,usewf_soften,usewf_md  !Determines if the previous wavefunction should be read for GEOPT, SOFTENING and MD
!!   read(12,*) findsym,finddos       !Determines if symmetry (isotropy)/FDOS calculation should be used
!!  ! read(12,*) units                !Either angstroem or bohr
!!  !Block KPT****************
!!   read(12,'(a250)') all_line
!!   n = len_trim(all_line)
!!   read(all_line(1:n),*) tmp_ch
!!    if(trim(tmp_ch)=="Auto") then
!!      read(all_line(1:n),*) tmp_ch,dkpt1,dkpt2
!!      ka=0;kb=0;kc=0
!!      auto_kpt=.true.
!!    else
!!     read(all_line(1:n),*) ka,kb,kc
!!     dkpt1=0.d0
!!     dkpt2=0.d0
!!     auto_kpt=.false.
!!    endif
!!  !Block KPT****************
!!  
!!  
!!  
!!  
!!  
!!  
!!  ! read(12,*) ka,kb,kc              !For fixed kpoint mesh
!!  !!To generate automatic kpoin mesh
!!  ! if(ka==0.and.kb==0.and.kc==0)then
!!  ! read(12,*) dkpt1,dkpt2
!!  ! else
!!  ! dkpt1=0.d0
!!  ! dkpt2=0.d0
!!  ! endif
!!   close(12) 
!!  
!!  !Initiallize LJ parameter if required
!!   if(trim(code)=="blj") call blj_init_parameter()
!!  
!!  !Get the correct atomic masses and atomic character
!!   do itype=1,ntypat
!!     call atmdata(amu(itype),rcov(itype),char_type(itype),znucl(itype))
!!   enddo
!!  
!!  !Put replace amu from file if desired
!!   if(amutmp(1).ne.0.d0) amu=amutmp
!!  end subroutine

!************************************************************************************

subroutine pathintegral(parini,latvec,xred)
 use mod_interface
 use global, only: nat,ntypat,znucl,amu,typat,char_type,units,target_pressure_habohr,fixat,fixlat
 use defs_basis
 use interface_code
! Main program to test potential subroutines
!       use parameters 
 use mod_parini, only: typ_parini
       implicit none
       type(typ_parini), intent(in):: parini
       integer count1,count2,count_rate,count_max,lwork,info,nint,i,iat,idispl,irep,iprec,istr,jstr,ilat,jlat
! nat: number of atoms (e.g. vacancy)
!       parameter(nint=2*32*1+1)
       real(8):: rxyz0(3,nat),fxyz(3,nat),displ(3,nat)
       real(8),allocatable:: work(:),simpson(:)
       real(8):: evals(3),s2(3,3),dmat(3,3),dproj(6),rotmat(3,3),xred(3,nat)
       real(8):: stepsize_at,stepsize_lat,t1,t2,t3,path,xred_in(3,nat),latvec_in(3,3),strten_in(6)
       real(8):: str_matrix(3,3),transformed(3,3),transformed_inv(3,3),fcart_in(3,nat)
       real(8):: ener,vol,sumx,sumy,sumz,ener0,etot_in,tstressall
       real(8):: dlat(6),latvec(3,3),latvecinv(3,3),stress(3,3),displat(3,3),tstress(3,3)
       real(8):: stressmod(3,3),trafo(3,3),latvectrans(3,3),s2t(3,3),metric(3,3)
       real(8):: angmom(3,3),latvecold(3,3),ener_comp1,time,dener,ener_comp2,count
       character(2)::atom
       character(40):: filename
       character(4):: fn4
       logical:: getwfk,atoms,cell
       open(unit=88,file="pathparams.in")
       read(88,*) atoms,cell !Read which degrees of freedom should be tested, usually only one of them
       read(88,*) stepsize_at,stepsize_lat !Stepsize of atomic and cell degrees of freedom
       read(88,*) nint !Number of steps along the path 
       if(.not.atoms) stepsize_at=0.d0
       if(.not.cell) stepsize_lat=0.d0
       if(atoms.and.abs(stepsize_at).le.1.d-12) stop "Increase stepsize for atoms"
       if(cell.and.abs(stepsize_lat).le.1.d-12) stop "Increase stepsize for lattice"
       allocate(simpson(nint)) 

       write(*,*) '# Testing ',trim(parini%potential_potential),' potential'
       if (mod(nint,2).ne.1) stop '# nint has to be odd'
       simpson(1)=1.d0/3.d0
       simpson(2)=4.d0/3.d0
       do i=3,nint-2,2
       simpson(i)=2.d0/3.d0
       simpson(i+1)=4.d0/3.d0
       enddo
       simpson(nint)=1.d0/3.d0

       count=0.d0
       call rxyz_int2cart(latvec,xred,rxyz0,nat)
! create random displacements (use sin() instead of rand())
!       stepsize=1.d-3
       call random_number(displ)
       displ=(displ-5.d-1)*2.d0
       displ=displ*stepsize_at
!       do iat=1,nat
!       displ(1,iat)=stepsize*sin(iat+.2d0)
!       displ(2,iat)=stepsize*sin(iat+.4d0)
!       displ(3,iat)=stepsize*sin(iat+.7d0)
!!       displ(1,iat)=stepsize*abs(sin(iat+.2d0))
!!       displ(2,iat)=stepsize*abs(sin(iat+.4d0))
!!       displ(3,iat)=stepsize*abs(sin(iat+.7d0))
!       enddo
       call random_number(displat)
       displat=(displat-5.d-1)*2.d0
       displat=displat*stepsize_lat   


!      displat=0.d0   
!      displ=0.d0

! calculate energy at equilibrium geometry (diamond structure)
! and at two additional points along the displacements
      call cpu_time(t1)
      call system_clock(count1,count_rate,count_max)
       path=0.d0
       do 1000,idispl=1,nint
       do 10,irep=1,1
       call rxyz_cart2int(latvec,xred_in,rxyz0,nat)
       latvec_in=latvec
       iprec=1; getwfk=.true.; if(idispl==1) getwfk=.false.
       call get_energyandforces_single(parini,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,getwfk);count=count+1
       fxyz=fcart_in
       ener=etot_in
!Convert strten into "stress", which is actually the forces on the cell vectors
       str_matrix(1,1)=strten_in(1)
       str_matrix(2,2)=strten_in(2)
       str_matrix(3,3)=strten_in(3)
       str_matrix(1,2)=strten_in(6)
       str_matrix(2,1)=strten_in(6)
       str_matrix(1,3)=strten_in(5)
       str_matrix(3,1)=strten_in(5)
       str_matrix(2,3)=strten_in(4)
       str_matrix(3,2)=strten_in(4)
       call getvol(latvec,vol)
       transformed(:,1)=latvec(1,:) 
       transformed(:,2)=latvec(2,:) 
       transformed(:,3)=latvec(3,:) 
       call invertmat(transformed,transformed_inv,3)
       stress=(-vol*matmul(str_matrix,transformed_inv))
!       call energyandforces(nat,latvec,rxyz0,fxyz,stress,0.d0,ener,count)
!****************************************************************************************************************        
     write(fn4,'(i4.4)') idispl
     filename='pospath_'//fn4//'.ascii'
     call write_atomic_file_ascii(parini,filename,nat,units,xred_in,latvec_in,fcart_in,strten_in,&
          &char_type,ntypat,typat,fixat,fixlat,etot_in,target_pressure_habohr,etot_in,0.d0)


10      continue
       if (idispl.eq.1) ener0=ener

!       if (idispl.eq.1 .or. idispl.eq.nint) then
! check whether total force vanishes
       sumx=0.d0
       sumy=0.d0
       sumz=0.d0
       do 8483,iat=1,nat
           sumx=sumx+fxyz(1,iat)
           sumy=sumy+fxyz(2,iat)
8483        sumz=sumz+fxyz(3,iat)

       write(*,'(a,i5,2(x,e19.12))') ' # idispl,ener,ener/nat',idispl,ener,ener/nat
       write(*,'(a,3(x,e10.3))') &
         ' # Sum of x, y and z component of forces:',sumx,sumy,sumz
!       endif


! integrate force*displacement
       t1=0.d0
       t2=0.d0
       t3=0.d0
       tstressall=0.d0
       
       do iat=1,nat
       t1=t1-fxyz(1,iat)*displ(1,iat)
       t2=t2-fxyz(2,iat)*displ(2,iat)
       t3=t3-fxyz(3,iat)*displ(3,iat)
       enddo
       
       do istr=1,3
       do jstr=1,3
       tstressall=tstressall-stress(istr,jstr)*displat(istr,jstr)
       enddo
       enddo
      


       path=path+simpson(idispl)*(tstressall+t1+t2+t3)


       trafo=0.d0
      
       latvecold=latvec
       call invertmat(latvec,latvecinv,3) 
       do ilat=1,3
       do jlat=1,3
       latvec(ilat,jlat)=latvec(ilat,jlat)+displat(ilat,jlat)
       enddo
       enddo


       call updaterxyz(latvecold,latvec,rxyz0,nat)


! next positions along path
       do iat=1,nat
       rxyz0(1,iat)=rxyz0(1,iat)+displ(1,iat)
       rxyz0(2,iat)=rxyz0(2,iat)+displ(2,iat)
       rxyz0(3,iat)=rxyz0(3,iat)+displ(3,iat)
       enddo


1000        continue
      call cpu_time(t2)
      call system_clock(count2,count_rate,count_max)

      call latvec2dproj(dproj,latvec,rotmat,rxyz0,nat)

!Calculate first energz with initial cell
      call rxyz_cart2int(latvec,xred_in,rxyz0,nat)
      latvec_in=latvec
      iprec=1; getwfk=.false. 
      call get_energyandforces_single(parini,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,getwfk)
      fxyz=fcart_in
      ener_comp1=etot_in
!Convert strten into "stress", which is actually the forces on the cell vectors
       str_matrix(1,1)=strten_in(1)
       str_matrix(2,2)=strten_in(2)
       str_matrix(3,3)=strten_in(3)
       str_matrix(1,2)=strten_in(6)
       str_matrix(2,1)=strten_in(6)
       str_matrix(1,3)=strten_in(5)
       str_matrix(3,1)=strten_in(5)
       str_matrix(2,3)=strten_in(4)
       str_matrix(3,2)=strten_in(4)
       call getvol(latvec,vol)
       transformed(:,1)=latvec(1,:)
       transformed(:,2)=latvec(2,:) 
       transformed(:,3)=latvec(3,:)
       call invertmat(transformed,transformed_inv,3)
       stress=(-vol*matmul(str_matrix,transformed_inv))
!      call energyandforces(nat,latvec,rxyz0,fxyz,stress,0.d0,ener_comp1,count1)
!****************************************************************************************************************   
      

      write(*,*) '# Final forces on the transformed cell'
      do iat=1,3
        write(*,'(a,3(es25.15))') " #  ",stress(:,iat)
      enddo
!Calculate angular momentum
      call cross_product(latvec(:,1),stress(:,1),angmom(:,1))
      call cross_product(latvec(:,2),stress(:,2),angmom(:,2))
      call cross_product(latvec(:,3),stress(:,3),angmom(:,3))
      write(*,*) '# Angular momentum on the cell'
      write(*,'(a,3(es25.15))') " #  ", angmom(:,1)+angmom(:,2)+angmom(:,3) 


      call latvec2dproj(dproj,latvec,rotmat,rxyz0,nat)

! compare energy difference with  force*displacement to check correctness of forces
       dener=ener-ener0
       write(*,*) '# Check correctness of forces'
       write(*,*) '# Difference of total energies ',dener
       write(*,*) '# Integral force*displacement  ',path
       write(*,*) '# Difference ',path-dener
       write(*,*) '# number of force evaluations (count)',int(count)
!       write(*,*) ' #CPU time ', t2-t1
!       time=(count2-count1)/float(count_rate)
       time=(t2-t1)
       write(*,*) '# elapsed time ', time
       write(*,*) '# time/count, time/(count*nat)',time/count, time/(count*nat)
!       write(*,*) 'Energy difference with diagonalized latvec:',ener_comp1-ener_comp2
       end subroutine

!****************************************************************************************************************   

subroutine plot_fp_grid(parini,nlminx,nlmin,nat,fp_len,fp_arr,lat_arr,pl_arr)
use mod_parini, only: typ_parini
use mod_interface
implicit none
type(typ_parini), intent(in):: parini
integer:: nlminx,nlmin,fp_len,i,kk,nat
real(8):: fp_arr(fp_len,nlminx),fp_dist
real(8):: tmp_acell(3),tmp_real,tmp_rprim(3,3),lat_arr(3,3,nlminx),pl_arr(3,nat,nlminx),randpos(3)
character(5):: fn,fn2
open(unit=11,file="fingerprint_grid.plot")
!!open(unit=11,file="oganov_grid.plot")
!!open(unit=12,file="calypso_grid.plot")
!!open(unit=14,file="malypso_grid.plot")
!!open(unit=13,file="e_grid.plot")
!!open(unit=15,file="list.plot")
!!!RANDOM ROTATIONS
!do kk=1,nlmin
!   tmp_rprim=0.d0
!   tmp_rprim(1,1)=1.d0
!   tmp_rprim(2,2)=1.d0
!   tmp_rprim(3,3)=1.d0
!!!!   call random_number(tmp_acell)
!!!!   tmp_acell=tmp_acell-0.5d0*3.d0
!!!!   call random_number(tmp_real)
!!!!   tmp_real=tmp_real-0.5d0*3.d0
!!!!   tmp_acell=tmp_acell/sqrt(tmp_acell(1)**2+tmp_acell(2)**2+tmp_acell(3)**2)
!!!!   call rotation(tmp_rprim,tmp_real,tmp_acell)
!   lat_arr(:,1,kk)=matmul(tmp_rprim,lat_arr(:,1,1))
!   lat_arr(:,2,kk)=matmul(tmp_rprim,lat_arr(:,2,1))
!   lat_arr(:,3,kk)=matmul(tmp_rprim,lat_arr(:,3,1))
!!!!   pl_arr(:,:,kk)=pl_arr(:,:,1)
!   call random_number(randpos)
!   do i=1,nat
!   pl_arr(:,i,kk)=pl_arr(:,i,1)+randpos(:)
!       pl_arr(1,i,kk)=modulo(modulo(pl_arr(1,i,kk),1.d0),1.d0)
!       pl_arr(2,i,kk)=modulo(modulo(pl_arr(2,i,kk),1.d0),1.d0)
!       pl_arr(3,i,kk)=modulo(modulo(pl_arr(3,i,kk),1.d0),1.d0)
!   enddo 
!enddo
!!!Get the oganov fp
!!!fp_method=11 !11: Oganov FP, 12: CALYPSO FP, 13: Modified CALYPSO 21: molecular gaussian overlap, 22: molecular sprint
!!  deallocate(fp_arr)
!!  call init_fp(fp_len)
!!  allocate(fp_arr(fp_len,nlminx))
!do kk=1,nlmin
!   call get_fp(fp_len,pl_arr(:,:,kk),lat_arr(:,:,kk),fp_arr(:,kk))
!enddo

!Get the fp distances and write to file
do kk=1,nlmin
write(fn,'(i5.5)') kk
open(unit=22,file="fingerprint_"//fn)
do i=1,fp_len
  write(22,*) fp_arr(i,kk)
enddo
close(22)
do i=1,nlmin
   write(fn2,'(i5.5)') i
!   open(unit=24,file="fingerprint_assign"//fn//"_"//fn2)
   call get_fp_distance(parini,fp_len,fp_arr(:,kk),fp_arr(:,i),fp_dist)
!   close(24)
   write(11,*) kk,i,fp_dist
!!   write(13,*) kk,i,abs(e_arr(i)-e_arr(kk))
!!   write(15,*) 11,fp_dist
enddo
write(11,*) " "
!!write(13,*) " "
enddo
close(11)
!!stop

!!!Get the calypso fp
!!fp_method=12 !11: Oganov FP, 12: CALYPSO FP, 21: molecular gaussian overlap, 22: molecular sprint
!!  deallocate(fp_arr)
!!  call init_fp(fp_len)
!!  allocate(fp_arr(fp_len,nlminx))
!!do kk=1,nlmin
!!   call get_fp(fp_len,pl_arr(:,:,kk),lat_arr(:,:,kk),fp_arr(:,kk))
!!enddo
!!
!!!Get the fp distances and write to file
!!do kk=1,nlmin
!!do i=1,nlmin
!!   call get_fp_distance(fp_len,fp_arr(:,kk),fp_arr(:,i),fp_dist)
!!   write(12,*) kk,i,fp_dist
!!   write(15,*) 12,fp_dist
!!enddo
!!write(12,*) " "
!!enddo
!!
!!!Get the malypso fp
!!fp_method=13 
!!  deallocate(fp_arr)
!!  call init_fp(fp_len)
!!  allocate(fp_arr(fp_len,nlminx))
!!do kk=1,nlmin
!!   call get_fp(fp_len,pl_arr(:,:,kk),lat_arr(:,:,kk),fp_arr(:,kk))
!!enddo
!!
!!!Get the fp distances and write to file
!!do kk=1,nlmin
!!do i=1,nlmin
!!   call get_fp_distance(fp_len,fp_arr(:,kk),fp_arr(:,i),fp_dist)
!!   write(14,*) kk,i,fp_dist
!!   write(15,*) 13,fp_dist
!!enddo
!!write(14,*) " "
!!enddo
!!
!!
!!!Here we write the plot script for gnuplot
!!open(unit=29,file="Howtoplot_ed")
!!write(29,'(a)') "set pm3d map"
!!write(29,'(a,i,a)') "set xrange[0:",nlmin+1,"]"
!!write(29,'(a,i,a)') "set yrange[0:",nlmin+1,"]"
!!write(29,'(a)') "set size square"
!!write(29,'(a)') "plot 'e_grid.plot' with image"
!!close(29)
!!open(unit=29,file="Howtoplot_oganov")
!!write(29,'(a)') "set pm3d map"
!!write(29,'(a,i,a)') "set xrange[0:",nlmin+1,"]"
!!write(29,'(a,i,a)') "set yrange[0:",nlmin+1,"]"
!!write(29,'(a)') "set size square"
!!write(29,'(a)') "plot 'oganov_grid.plot' with image"
!!close(29)
!!open(unit=29,file="Howtoplot_calypso")
!!write(29,'(a)') "set pm3d map"
!!write(29,'(a,i,a)') "set xrange[0:",nlmin+1,"]"
!!write(29,'(a,i,a)') "set yrange[0:",nlmin+1,"]"
!!write(29,'(a)') "set size square"
!!write(29,'(a)') "plot 'calypso_grid.plot' with image"
!!close(29)
!!open(unit=29,file="Howtoplot_malypso")
!!write(29,'(a)') "set pm3d map"
!!write(29,'(a,i,a)') "set xrange[0:",nlmin+1,"]"
!!write(29,'(a,i,a)') "set yrange[0:",nlmin+1,"]"
!!write(29,'(a)') "set size square"
!!write(29,'(a)') "plot 'malypso_grid.plot' with image"
!!close(29)
!!stop
!Here we write the plot script for gnuplot
open(unit=29,file="Howtoplot_fingerprint")
write(29,'(a)') "set pm3d map"
write(29,'(a,i5,a)') "set xrange[0:",nlmin+1,"]"
write(29,'(a,i5,a)') "set yrange[0:",nlmin+1,"]"
write(29,'(a)') "set size square"
write(29,'(a)') "plot 'fingerprint_grid.plot' with image"
close(29)

end subroutine

!****************************************************************************************************************   

subroutine rotate_like_crazy(parini,latvec,xred,tolmin,tolmax,ntol)
 use mod_interface
 use global, only: nat,ntypat,znucl,amu,typat,char_type,units,target_pressure_habohr,target_pressure_gpa
 use global, only: fixat,fixlat,findsym,fragarr,ka,kb,kc
 use defs_basis
 use interface_code
! Main program to test potential subroutines
!       use parameters 
 use mod_parini, only: typ_parini
       implicit none
       type(typ_parini), intent(in):: parini
       integer::  iprec,nstruct,i,ntol,spgint
! nat: number of atoms (e.g. vacancy)
!       parameter(nint=2*32*1+1)
       real(8):: latvec(3,3),xred(3,nat),pinit,pfinal,psteps,latvec0(3,3),xred0(3,nat),vel_vol_in
       real(8):: vel_in(3,nat),vel_lat_in(3,3),fcart(3,nat),strten(6),axis(3),rotmat(3,3),angle
       real(8):: counter,count_geopt,enthalpy,energy,vol,ext_press,tolmin,tolmax,spgtol_pos,spg_pos
       character(2)::atom
       character(40):: filename,filenameout,folder
       character(400):: command
       character(7):: fn
       character(5):: fn5
       logical:: getwfk,reuse_str,file_exists,readfix,readfrag
!Define current directory
  folder=""
!Read from the rotation file
  open(unit=45,file="rotate.in")
  read(45,*) nstruct
  close(45)

!Now run over all poslows found in the folder
do i=1,nstruct

call random_number(axis(:))
call random_number(angle)
        axis=axis-0.5d0
        angle=(angle-0.5d0)*10.d0
        axis=axis/sqrt(sum(axis(:)**2))
        call rotation(rotmat,angle,axis)
        latvec=matmul(rotmat,latvec)
!Compute energy
        call get_energyandforces_single(parini,latvec,xred,fcart,strten,energy,iprec,getwfk)
        call get_enthalpy(latvec,energy,target_pressure_habohr,enthalpy)
        call getvol(latvec,vol)
        ext_press=strten(1)+strten(2)+strten(3)
        ext_press=-ext_press/3.d0*HaBohr3_GPa
        write(fn5,'(i5.5)') i
        filename = 'poslow'//fn5
        filenameout=trim(filename)//".rotated.ascii"
        call write_atomic_file_ascii(parini,filenameout,nat,units,xred,latvec,fcart,strten,&
             &char_type,ntypat,typat,fixat,fixlat,energy,target_pressure_habohr,ext_press,enthalpy)
        open(unit=2,file="Enthalpies_rotated",access="append")
        write(2,'(a,5(1x,es25.15),i5)')trim(filename),target_pressure_gpa,&
             &enthalpy*Ha_ev,energy*Ha_ev,ext_press,vol*Bohr_Ang**3,spgint
        close(2)
enddo
write(*,'(a)') " # Finished rotating structures. Output in _DIR folders and in the file 'Enthalpies'"
write(*,'(a)') " # Have a nice day!"

end subroutine

!****************************************************************************************************************   
!****************************************************************************************************************   

subroutine poslowrelax(parini,latvec,xred,tolmin,tolmax,ntol)
 use mod_interface
 use global, only: nat,ntypat,znucl,amu,typat,char_type,units,target_pressure_habohr,target_pressure_gpa
 use global, only: fixat,fixlat,findsym,fragarr,ka,kb,kc
 use defs_basis
 use interface_code
 use mod_parini, only: typ_parini
! Main program to test potential subroutines
       implicit none
       type(typ_parini), intent(in):: parini
!       use parameters 
       integer::  iprec,nstruct,i,ntol,spgint
! nat: number of atoms (e.g. vacancy)
!       parameter(nint=2*32*1+1)
       real(8):: latvec(3,3),xred(3,nat),pinit,pfinal,psteps,latvec0(3,3),xred0(3,nat),vel_vol_in
       real(8):: vel_in(3,nat),vel_lat_in(3,3),fcart(3,nat),strten(6)
       real(8):: counter,count_geopt,enthalpy,energy,vol,ext_press,tolmin,tolmax,spgtol_pos,spg_pos
       character(2)::atom
       character(40):: filename,filenameout,folder
       character(400):: command
       character(7):: fn
       character(5):: fn5
       logical:: getwfk,reuse_str,file_exists,readfix,readfrag
!Define current directory
  folder=""

open(unit=4,file="poslowrelax.in")
read(4,*) nstruct
close(4)

count_geopt=0.d0
open(unit=2,file="Enthalpies_poslow")
write(2,'(a)')" # poslow.index     TrgtPressure(GPa)        Enthalpy(eV)              Energy(eV)&
              &                 ExtPressure(GPa)          Volume(ang^3)    spg"
close(2)
!Now run over all poslows found in the folder
do i=1,nstruct
        write(fn5,'(i5.5)') i
        filename = 'poslow'//fn5//'.ascii'
        INQUIRE(FILE=trim(filename), EXIST=file_exists)  
        readfix=.false.
        readfrag=.false.
if(file_exists) then
   call read_atomic_file_ascii(filename,nat,units,xred,latvec,fcart,strten,fixat,fixlat,readfix,fragarr,readfrag,enthalpy,energy)
   write(*,'(a,es15.7)') " # Now running "//trim(filename)//" at pressure GPa ",target_pressure_gpa
   !First remove all posgeopt files
   call system('rm -f posgeopt.*.ascii')
   !Call geometry optimizer 
     vel_in=0.d0
     vel_vol_in=0.d0
     vel_lat_in=0.d0
      iprec=1
      if(parini%geopt_ext) then
        call geopt_external(parini,latvec,xred,fcart,strten,energy,iprec,ka,kb,kc,counter)
      else
        if(parini%paropt_geopt%approach=="RBFGS")  call GEOPT_RBFGS_MHM(parini,latvec,xred,fcart,strten,energy,iprec,counter,folder)
        if(parini%paropt_geopt%approach=="MBFGS")  call GEOPT_MBFGS_MHM(parini,latvec,xred,fcart,strten,energy,iprec,counter,folder)
        if(parini%paropt_geopt%approach=="FIRE")   call GEOPT_FIRE_MHM(parini,latvec,xred,fcart,strten,vel_in,vel_lat_in,vel_vol_in,energy,iprec,counter,folder)
        if(parini%paropt_geopt%approach=="SQNM")   call     GEOPT_SQNM(parini,latvec,xred,fcart,strten,energy,iprec,counter,folder)
        if(parini%paropt_geopt%approach=="QBFGS")  call    GEOPT_QBFGS(parini,latvec,xred,fcart,strten,energy,iprec,counter,folder)
        if(parini%paropt_geopt%approach=="SD")     call     GEOPT_SD  (parini,latvec,xred,fcart,strten,energy,iprec,counter,folder)
      endif
!Check for symmetry
   if(findsym) then
      call find_symmetry(parini,nat,xred,latvec,typat,tolmin,tolmax,ntol,spgtol_pos,spgint)
      spg_pos=real(spgint,8)
   else
      spg_pos=0.d0
      spgtol_pos=0.d0
   endif
!Update GEOPT counter
       count_geopt=count_geopt+counter     
       write(*,'(a,i7)') " # Counter of GEOPT updated: ", int(count_geopt)
       call get_enthalpy(latvec,energy,target_pressure_habohr,enthalpy)
       call getvol(latvec,vol)
       ext_press=strten(1)+strten(2)+strten(3)
       ext_press=-ext_press/3.d0*HaBohr3_GPa
       filenameout=trim(filename)//"relaxed.ascii"
       call write_atomic_file_ascii(parini,filenameout,nat,units,xred,latvec,fcart,strten,&
             &char_type,ntypat,typat,fixat,fixlat,energy,target_pressure_habohr,ext_press,enthalpy)
       open(unit=2,file="Enthalpies_poslow",access="append")
       write(2,'(a,5(1x,es25.15),i5)')trim(filename),target_pressure_gpa,&
             &enthalpy*Ha_ev,energy*Ha_ev,ext_press,vol*Bohr_Ang**3,spgint
       close(2)
       command="mkdir "//trim(filename)//"_DIR"
       call system(trim(command)) 
       command="cp posgeopt.*.ascii "//trim(filename)//"_DIR"
       call system(trim(command)) 
endif
enddo
write(*,'(a)') " # Finished relaxing structures. Output in _DIR folders and in the file 'Enthalpies'"
write(*,'(a)') " # Have a nice day!"


end subroutine

!****************************************************************************************************************   

!****************************************************************************************************************   

subroutine enthalpyrelax(parini,latvec,xred,tolmin,tolmax,ntol,findsym)
 use mod_interface
 use global, only: nat,ntypat,znucl,amu,typat,char_type,units,target_pressure_habohr,target_pressure_gpa
 use global, only: fixat,fixlat
 use defs_basis
 use interface_code
 use mod_parini, only: typ_parini
! Main program to test potential subroutines
       implicit none
       type(typ_parini), intent(in):: parini
       integer::  ka,kb,kc,iprec
! nat: number of atoms (e.g. vacancy)
!       parameter(nint=2*32*1+1)
       real(8):: latvec(3,3),xred(3,nat),pinit,pfinal,psteps,pcur,latvec0(3,3),xred0(3,nat),vel_vol_in
       real(8):: vel_in(3,nat),vel_lat_in(3,3),fcart(3,nat),strten(6)
       real(8):: counter,count_geopt,enthalpy,energy,vol,ext_press
       character(2)::atom
       character(40):: filename,folder
       character(400):: command
       character(7):: fn
       logical:: getwfk,reuse_str
       real(8):: tolmin,tolmax,spgtol
       integer:: ntol,spgint
       logical:: findsym
       open(unit=88,file="enthalpyparams.in")
       read(88,*) pinit,pfinal,psteps,reuse_str !Read which degrees of freedom should be tested, usually only one of them
!Define current directory
  folder=""


latvec0=latvec
xred0=xred
count_geopt=0.d0
open(unit=2,file="Enthalpies")
write(2,'(a)')" #   TrgtPressure(GPa)        Enthalpy(eV)              Energy(eV)&
                &                 ExtPressure(GPa)          Volume(ang^3)             SPG"

close(2)
!Now run over all pressures
pcur=pinit
do while (pcur.le.pfinal)
!pcur=pinit,pfinal,psteps

write(*,'(a,es15.7)') " # Now running at pressure GPa ",pcur
if(.not.reuse_str) latvec=latvec0
if(.not.reuse_str) xred=xred0
!First remove all posgeopt files
call system('rm -f posgeopt.*.ascii')
 target_pressure_gpa=pcur                                   !Target pressure in GPA
 target_pressure_habohr=target_pressure_gpa/HaBohr3_GPA
!Call geometry optimizer 
  vel_in=0.d0
  vel_lat_in=0.d0
   iprec=1
   if(parini%geopt_ext) then
     call geopt_external(parini,latvec,xred,fcart,strten,energy,iprec,ka,kb,kc,counter)
   else
     if(parini%paropt_geopt%approach=="RBFGS")  call GEOPT_RBFGS_MHM(parini,latvec,xred,fcart,strten,energy,iprec,counter,folder)
     if(parini%paropt_geopt%approach=="MBFGS")  call GEOPT_MBFGS_MHM(parini,latvec,xred,fcart,strten,energy,iprec,counter,folder)
     if(parini%paropt_geopt%approach=="FIRE")   call GEOPT_FIRE_MHM(parini,latvec,xred,fcart,strten,vel_in,vel_lat_in,vel_vol_in,energy,iprec,counter,folder)
     if(parini%paropt_geopt%approach=="SQNM")   call     GEOPT_SQNM(parini,latvec,xred,fcart,strten,energy,iprec,counter,folder)
     if(parini%paropt_geopt%approach=="QBFGS")  call    GEOPT_QBFGS(parini,latvec,xred,fcart,strten,energy,iprec,counter,folder)
     if(parini%paropt_geopt%approach=="SD")     call     GEOPT_SD  (parini,latvec,xred,fcart,strten,energy,iprec,counter,folder)
   endif
!Check for symmetry
   if(findsym) then
      call find_symmetry(parini,nat,xred,latvec,typat,tolmin,tolmax,ntol,spgtol,spgint)
   else
      spgint=0
      spgtol=0.d0
   endif
!Update GEOPT counter
    count_geopt=count_geopt+counter     
    write(*,'(a,i7)') " # Counter of GEOPT updated: ", int(count_geopt)
    call get_enthalpy(latvec,energy,target_pressure_habohr,enthalpy)
    call getvol(latvec,vol)
    ext_press=strten(1)+strten(2)+strten(3)
    ext_press=-ext_press/3.d0*HaBohr3_GPa
    write(fn,'(f7.2)') pcur
    fn = repeat( '0', 7-len_trim(adjustl(fn))) // adjustl(fn)
    write(*,*) fn
    filename="posenth."//fn//".ascii"
    call write_atomic_file_ascii(parini,filename,nat,units,xred,latvec,fcart,strten,&
          &char_type,ntypat,typat,fixat,fixlat,energy,target_pressure_habohr,ext_press,enthalpy)
    open(unit=2,file="Enthalpies",access="append")
    write(2,'(5(1x,es25.15),3x,i5)') target_pressure_gpa,enthalpy*Ha_ev,energy*Ha_ev,ext_press,vol*Bohr_Ang**3,spgint
    close(2)
    command="mkdir GPa"//fn
    call system(trim(command)) 
    command="cp posgeopt.*.ascii GPa"//fn
    call system(trim(command)) 
pcur=pcur+psteps
enddo
write(*,'(a)') " # Finished relaxing structures. Output in GPa folders and in the file 'Enthalpies'"
write(*,'(a)') " # Have a nice day!"


end subroutine

!****************************************************************************************************************   

subroutine varvol(parini,latvec,xred,tolmin,tolmax,ntol,findsym)
 use mod_interface
 use global, only: nat,ntypat,znucl,amu,typat,char_type,units,target_pressure_habohr,target_pressure_gpa
 use global, only: fixat,fixlat
 use defs_basis
 use interface_code
 use mod_parini, only: typ_parini
! Main program to test potential subroutines
       implicit none
       type(typ_parini), intent(in):: parini
       integer::  ka,kb,kc,iprec
! nat: number of atoms (e.g. vacancy)
!       parameter(nint=2*32*1+1)
       real(8):: latvec(3,3),xred(3,nat),latvec0(3,3),xred0(3,nat)
       real(8):: vel_in(3,nat),vel_lat_in(3,3),fcart(3,nat),strten(6)
       real(8):: counter,count_geopt,enthalpy,energy,vol,ext_press
       character(2)::atom
       character(40):: filename,folder
       character(400):: command
       character(12):: fn
       logical:: getwfk,reuse_str
       real(8):: tolmin,tolmax,spgtol
       integer:: ntol,spgint,itime
       logical:: findsym,is_percentage
       real(8):: vfmin,vfmax,vfsteps,vmin,vmax,vsteps,vcur,vol0,vinit,vfinal,convert,alpha
       open(unit=88,file="varvol.in")
       read(88,*) vfmin,vfmax,vfsteps,is_percentage 
!Convert units: usually the unis are given as in the poscur. Assume everything is in volume/atom
if(trim(units)=="angstroem") then
  convert=Bohr_Ang
else
  convert =1.d0
endif

!Define current directory
folder=""
latvec0=latvec
xred0=xred
call getvol(latvec0,vol0)
  if(.not.is_percentage) then
     vfmin=nat*vfmin/vol0/convert**3
     vfmax=nat*vfmax/vol0/convert**3
     vfsteps=nat*vfsteps/vol0/convert**3
  endif
   

open(unit=2,file="VolEnergies")
write(2,'(a)')" #   TrgtPressure(GPa)        Enthalpy(eV)              Energy(eV)&
                &                 Energy(eV/atom)           ExtPressure(GPa)          Volume(ang^3)             Volume(ang^3/atom)        Volume(%)                 SPG"

close(2)
!Now run over all volumes
vinit=vfmin*vol0
vfinal=vfmax*vol0
vsteps=vfsteps*vol0
vcur=vinit
itime=1
do while (vcur.le.vfinal+1.d-10)
   write(*,'(a,es15.7)') " # Now running at volume per atom ",vcur*convert**3/nat
!Translate vcur to a lattice vector
   alpha=vcur/vol0
   alpha=alpha**(1.d0/3.d0)
   latvec=latvec0*alpha
!First remove all posgeopt files
   call system('rm -f posgeopt.*.ascii')
       if(itime==1)  then
           getwfk=.false.
       else
           getwfk=.true.
       endif
   call get_energyandforces_single(parini,latvec,xred0,fcart,strten,energy,iprec,getwfk)
!Check for symmetry
   if(findsym) then
      call find_symmetry(parini,nat,xred,latvec,typat,tolmin,tolmax,ntol,spgtol,spgint)
   else
      spgint=0
      spgtol=0.d0
   endif
   call get_enthalpy(latvec,energy,target_pressure_habohr,enthalpy)
   call getvol(latvec,vol)
   ext_press=strten(1)+strten(2)+strten(3)
   ext_press=-ext_press/3.d0*HaBohr3_GPa
   if(is_percentage) then
      write(fn,'(f12.3)') vcur/vol0
   else
      write(fn,'(f12.3)') vcur*convert**3/nat
   endif
   fn = repeat( '0', 12-len_trim(adjustl(fn))) // adjustl(fn)
   filename="posvol."//fn//".ascii"
   call write_atomic_file_ascii(parini,filename,nat,units,xred,latvec,fcart,strten,&
         &char_type,ntypat,typat,fixat,fixlat,energy,target_pressure_habohr,ext_press,enthalpy)
   open(unit=2,file="VolEnergies",access="append")
   write(2,'(7(1x,es25.15),1x,es25.8,3x,i5)') target_pressure_gpa,enthalpy*Ha_ev,energy*Ha_ev,energy*Ha_ev/nat,ext_press,vol*Bohr_Ang**3,vol*Bohr_Ang**3/nat,vol/vol0,spgint
   close(2)
   vcur=vcur+vsteps
   itime=itime+1
enddo
write(*,'(a)') " # Finished relaxing structures. Output in GPa folders and in the file 'VolEnergies'"
write(*,'(a)') " # Have a nice day!"


end subroutine

!****************************************************************************************************************   

 subroutine updaterxyz(latvecold,latvecnew,rxyz,nat)
 use mod_interface
 !This subroutine will update the atomic positions in the cartesian coordinates after the cell shape has been changed according
 !to the change in the lattice vectors thus keeping the relative coordinates of all atoms constant
 implicit none
 real(8), intent(in)   :: latvecold(3,3), latvecnew(3,3)
 real(8), intent(inout):: rxyz(3,nat)
 integer, intent(in)   :: nat
 real(8)               :: latvecold_inv(3,3),trafo(3,3)
 integer               :: iat
 call invertmat(latvecold,latvecold_inv,3)
 trafo=matmul(latvecnew,latvecold_inv) !Transform atomic positions according to cell
 do iat=1,nat
 rxyz(:,iat)=matmul(trafo,rxyz(:,iat))
 enddo
 end subroutine

!****************************************************************************************************************   


subroutine k_expansion(latvec,xred,ka,kb,kc,k_latvec,k_xcart)
use mod_interface
!This routine expands the cell defined by latevec to a supercell
!of dimension ka,kb,kc. The atomic positions in real space
!will be returned in k_xcart. k_nat will then be the number of
!all atoms in the supercell and is ka*kb*kc*nat
use global, only: nat
implicit none
real(8):: latvec(3,3),k_latvec(3,3),k_xcart(3,nat,ka,kb,kc),xred(3,nat) 
integer:: iat,k,l,m,ka,kb,kc
do k=1,ka
do l=1,kb
do m=1,kc
do iat=1,nat
   k_xcart(:,iat,k,l,m)=matmul(latvec,xred(:,iat))+&
   &real(k-1,8)*latvec(:,1)+real(l-1,8)*latvec(:,2)+real(m-1,8)*latvec(:,3)
enddo
enddo
enddo
enddo
k_latvec(:,1)=real(ka,8)*latvec(:,1)
k_latvec(:,2)=real(kb,8)*latvec(:,2)
k_latvec(:,3)=real(kc,8)*latvec(:,3)
end subroutine

!************************************************************************************

subroutine elim_fixed_at(nat,x)
use mod_interface
use global, only: fixat
implicit none
integer:: iat,nat
real(8):: x(3,nat)
!write(*,*) "# Eliminiate"
do iat=1,nat
  if(fixat(iat)) then
     x(:,iat)=0.d0

  endif
!  write(*,*) x(:,iat)
enddo
end subroutine

!************************************************************************************

subroutine elim_fixed_lat(latvec,x)
use mod_interface
use global, only: fixlat,bc
implicit none
integer:: i,k
real(8):: x(3,3),latvec(3,3),lenlat,tmpvec(3)
real(8):: len1,len2,len3,tmp1(3),tmp2(3),tmp3(3),tmp1len,tmp2len,tmp3len
!The fixlat has 7 components:
!a,b,c,alha,beta,gamma, and for fixed cell shape (volume fluctuation)

!We should perform the projection out of the forces self-consistently
do k=1,10
if(fixlat(7)) then
!Here we implement cell fluctuation
!There are only forces along the cell vectors
   call diagcomp(latvec,x)
   return
elseif(all(fixlat(1:6))) then
!When the whole cell is fixed
   x=0.d0
   return
elseif(bc==2) then
!If the boundary condition is free
   x=0.d0
   return
else
!Treat the case where a, b, or c are fixed
do i=1,3
   if(fixlat(i)) then
!Project out the component of x onto latvec
   lenlat=dot_product(latvec(:,i),latvec(:,i))
   tmpvec=dot_product(x(:,i),latvec(:,i))*latvec(:,i)/lenlat
   x(:,i)=x(:,i)-tmpvec(:)
   endif   
enddo

!Now eliminate the compoments of x which would change the angle between the lattice vectors
!The correct way:
do i=1,3
   if(fixlat(i+3)) then
   len1=dot_product(latvec(:,modulo(i,3)+1),latvec(:,modulo(i,3)+1))
   len2=dot_product(latvec(:,modulo(i+1,3)+1),latvec(:,modulo(i+1,3)+1))
   call cross_product(latvec(:,modulo(i,3)+1),latvec(:,modulo(i+1,3)+1),tmpvec)
   call cross_product(tmpvec,latvec(:,modulo(i,3)+1),tmp1)
   call cross_product(tmpvec,latvec(:,modulo(i+1,3)+1),tmp2)
   tmp1len=dot_product(tmp1,tmp1)  
   tmp2len=dot_product(tmp2,tmp2)  
!Project out these components
   tmpvec=dot_product(x(:,modulo(i,3)+1),tmp1)*tmp1/tmp1len
   x(:,modulo(i,3)+1)=x(:,modulo(i,3)+1)-tmpvec
   tmpvec=dot_product(x(:,modulo(i+1,3)+1),tmp2)*tmp2/tmp2len
   x(:,modulo(i+1,3)+1)=x(:,modulo(i+1,3)+1)-tmpvec
   endif
enddo
endif
!write(*,*) "elim_fixed_lat norm", sqrt(sum(x*x))
enddo
!!The simple way:
!!This means only keep the components along the lattice vectors...
!do i=1,3
!   if(fixlat(i+3)) then
!!   write(*,*) i+3,modulo(i,3)+1,modulo(i+1,3)+1
!   lenlat=dot_product(latvec(:,modulo(i,3)+1),latvec(:,modulo(i,3)+1))
!   tmpvec=dot_product(x(:,modulo(i,3)+1),latvec(:,modulo(i,3)+1))*latvec(:,modulo(i,3)+1)/lenlat
!!   write(*,*) lenlat
!   x(:,modulo(i,3)+1)=tmpvec(:)
!   lenlat=dot_product(latvec(:,modulo(i+1,3)+1),latvec(:,modulo(i+1,3)+1))
!   tmpvec=dot_product(x(:,modulo(i+1,3)+1),latvec(:,modulo(i+1,3)+1))*latvec(:,modulo(i+1,3)+1)/lenlat
!!   write(*,*) lenlat
!   x(:,modulo(i+1,3)+1)=tmpvec(:)
!   endif   
!enddo

!!write(*,*) x(:,1)
!!write(*,*) x(:,2)
!!write(*,*) x(:,3)
end subroutine

!************************************************************************************

subroutine diagcomp(latvec,x)
use mod_interface
implicit none
real(8):: latvec(3,3),x(3,3),xnrm,latvect(3,3),latvecinv(3,3),sigma(3,3)
real(8):: len1,len2,len3,tmp1,tmp2,tmp3,tmpvec(3),vol
!There are only forces along the cell vectors
!Convert the target stress to the target gradient, assuming cell volume to be 1
   len1=sqrt(dot_product(latvec(:,1),latvec(:,1)))
   len2=sqrt(dot_product(latvec(:,2),latvec(:,2)))
   len3=sqrt(dot_product(latvec(:,3),latvec(:,3)))
   tmp1=dot_product(x(:,1),latvec(:,1))/len1
   tmp2=dot_product(x(:,2),latvec(:,2))/len2
   tmp3=dot_product(x(:,3),latvec(:,3))/len3
!!!!Convert the target stress to the target gradient, assuming cell volume to be 1
   latvect(:,1)=latvec(1,:)
   latvect(:,2)=latvec(2,:)
   latvect(:,3)=latvec(3,:)
   sigma=matmul(x,latvect)
   xnrm=sigma(1,1)+sigma(2,2)+sigma(3,3)
   xnrm=xnrm/3.d0
   sigma=0.d0
   sigma(1,1)=xnrm
   sigma(2,2)=xnrm
   sigma(3,3)=xnrm
   call invertmat(latvect,latvecinv,3)
   x=matmul(sigma,latvecinv)

   call getvol(latvec,vol)
   vol=vol**(1.d0/3.d0)
!   xnrm=xnrm/vol
!   xnrm=(tmp1/len1+tmp2/len2+tmp3/len3)/3.d0*vol
   xnrm=(x(1,1)+x(2,2)+x(3,3))/3.d0
   vol=1.d0/vol
   x(:,1)=latvec(:,1)*vol*xnrm
   x(:,2)=latvec(:,2)*vol*xnrm
   x(:,3)=latvec(:,3)*vol*xnrm
end subroutine


!************************************************************************************

subroutine slab_stress(flat,fix_z)
use mod_interface
!This routine will eliminate all z-components of the first two cell forces,
!and all x-y-components of the last cell force 
implicit none
integer:: i,j
real(8):: flat(3,3),ekin1,ekin2
logical:: fix_z
ekin1=0.d0
do i=1,3
 do j=1,3
 ekin1=ekin1+flat(i,j)**2
 enddo
enddo

flat(3,1)=0.d0
flat(3,2)=0.d0
flat(1,3)=0.d0
flat(2,3)=0.d0
if(fix_z) flat(3,3)=0.d0

ekin2=0.d0
do i=1,3
 do j=1,3
 ekin2=ekin2+flat(i,j)**2
 enddo
enddo
write(*,'(a,3(1x,es20.10))') "# Befor/After eliminating restricted components: ",ekin1,ekin2,ekin1-ekin2
end subroutine

!************************************************************************************

subroutine propagate(nat,xred,latvec0,dxred,dlatvec,xredout,latvecout)
use mod_interface
!This subroutine will propagate the coordinates of the atoms and the cells according to the
!value of fixed or free degrees of freedom according to
!xred=xred+dxred,latvec=latvec+dlatvec
use global, only: fixat, fixlat
implicit none
integer::nat,i,iat,j
real(8):: xred(3,nat),latvec(3,3),dxred(3,nat),dlatvec(3,3),xredout(3,nat),latvecout(3,3),len1,len2
real(8):: orig_angle(3),new_angle(3),axis(3),rotmat(3,3),center(3),latvec0(3,3)

latvec=latvec0
!write(*,*) "In propagate", fixlat
!Eliminate components not to be changed
if(any(fixlat)) call elim_fixed_lat(latvec,dlatvec)
if(any(fixat))  call elim_fixed_at(nat,dxred)

!Propagate
xredout=xred+dxred
latvecout=latvec+dlatvec
!if(fixlat(7)) write(*,*) "We have fixed cell shape"


!We need to fix the angles since we propagated along a linear direction, not rotational
if(any(fixlat(4:6)).and.(.not.all(fixlat(1:6)))) then
do j=1,10 !Do iteratively
!  write(*,*) "We have fixed cell angles",fixlat(4:6)
  do i=1,3
  if(fixlat(3+i)) then
!First get the original angle
  if(j==1)   orig_angle(i)=acos(dot_product(latvec(:,modulo(i,3)+1),latvec(:,modulo(i+1,3)+1))/&
             &sqrt(dot_product(latvec(:,modulo(i,3)+1),latvec(:,modulo(i,3)+1))*&
             &dot_product(latvec(:,modulo(i+1,3)+1),latvec(:,modulo(i+1,3)+1))))
             new_angle(i) =acos(dot_product(latvecout(:,modulo(i,3)+1),latvecout(:,modulo(i+1,3)+1))/&
             &sqrt(dot_product(latvecout(:,modulo(i,3)+1),latvecout(:,modulo(i,3)+1))*&
             &dot_product(latvecout(:,modulo(i+1,3)+1),latvecout(:,modulo(i+1,3)+1))))
   write(*,'(a,2(es25.15),1x,i5)') " # Old and New angle: ",orig_angle(i),new_angle(i),j
   call cross_product(latvecout(:,modulo(i,3)+1),latvecout(:,modulo(i+1,3)+1),axis)
   axis=axis/sqrt(dot_product(axis,axis))

   len1=sqrt(dot_product(latvecout(:,modulo(i,3)+1),latvecout(:,modulo(i,3)+1)))
   len2=sqrt(dot_product(latvecout(:,modulo(i+1,3)+1),latvecout(:,modulo(i+1,3)+1)))

   center=0.5d0*(latvecout(:,modulo(i,3)+1)/len1+latvecout(:,modulo(i+1,3)+1)/len2)
   center=center/sqrt(dot_product(center,center))

   call rotation(rotmat,-0.5d0*orig_angle(i),axis)
   latvecout(:,modulo(i,3)+1)=matmul(rotmat,center)*len1
   call rotation(rotmat,0.5d0*orig_angle(i),axis)
   latvecout(:,modulo(i+1,3)+1)=matmul(rotmat,center)*len2
   endif
   enddo
   latvec=latvecout
enddo
endif

!If the cell length were to be fixed, we have to rescale again
if(any(fixlat(1:3))) then
!  write(*,*) "We have fixed cell length",fixlat(1:3)
  do i=1,3
  if(fixlat(i)) then
   len1=dot_product(latvec0(:,i),latvec0(:,i))
   len2=dot_product(latvecout(:,i),latvecout(:,i))
   latvecout(:,i)=latvecout(:,i)*sqrt(len1/len2)
  endif
  enddo
endif
end subroutine propagate

!************************************************************************************

subroutine convcheck(nat,latvec_in,fcart_in,strten_in,target_pressure_habohr,strfact,fmax,fmax_at,fmax_lat,tolmxf,iexit)
use mod_interface
use global, only: fixat, fixlat
use defs_basis
implicit none
integer:: nat, iexit,iat,istr,i
real(8):: latvec_in(3,3),fcart_in(3,nat),strten_in(6),target_pressure_habohr,fmax,dstr(6)
real(8):: tolmxf,strtarget(6),strfact,fmax_at,fmax_lat
real(8):: gradtarget(3,3),strmattarget(3,3),latvect(3,3),latvectinv(3,3),flattarget(3,3)
real(8):: strmat(3,3),flat(3,3),dflat(3,3),dstrmat(3,3),dist_ang(6)
!Compute maximal component of forces, EXCLUDING any fixed components
fmax_at=0.0d0
fmax_lat=0.0d0
if(.not.any(fixat)) then
 do iat=1,nat
   do i=1,3
       if(abs(fcart_in(i,iat)) >= fmax_at ) fmax_at=abs(fcart_in(i,iat))
   end do
 end do
else
 do iat=1,nat
   do i=1,3
     if (.not.fixat(iat)) then
       if( abs(fcart_in(i,iat)) >= fmax_at ) fmax_at=abs(fcart_in(i,iat))
     end if
   end do
 end do
endif


strtarget=0.d0
strtarget(1:3)=-target_pressure_habohr
if(.not.any(fixlat)) then
 dstr(:)=strten_in(:)-strtarget(:)
!Evaluate the convergence
 do istr=1,6
     if(abs(dstr(istr))*strfact >= fmax_lat ) fmax_lat=abs(dstr(istr))*strfact
 end do
elseif(fixlat(7)) then
 dstr(:)=strten_in(:)-strtarget(:)
 fmax_lat=strfact*abs(((dstr(1)+dstr(2)+dstr(3))/3.d0))
! write(*,'(a,5(es25.15))') "Pressure", HaBohr3_GPa*(strten_in(1)+strten_in(2)+strten_in(3))/3.d0,HaBohr3_GPa*strten_in(1),HaBohr3_GPa*strten_in(2),HaBohr3_GPa*strten_in(3),abs(((dstr(1)+dstr(2)+dstr(3))/3.d0))*HaBohr3_GPa
! call dist_latvec2ang(dist_ang,latvec_in,pi)
! write(*,'(a,6(es25.15))') "angdeg", dist_ang
elseif(all(fixlat(1:6))) then
 fmax_lat=0.d0
else
!Convert sigma target to the gradient target
!Get full stress matrix of target
 strmattarget(1,1)=strtarget(1)
 strmattarget(2,2)=strtarget(2)
 strmattarget(3,3)=strtarget(3)
 strmattarget(1,2)=strtarget(6)
 strmattarget(2,1)=strtarget(6)
 strmattarget(1,3)=strtarget(5)
 strmattarget(3,1)=strtarget(5)
 strmattarget(2,3)=strtarget(4)
 strmattarget(3,2)=strtarget(4)
!Convert the target stress to the target gradient, assuming cell volume to be 1
 latvect(:,1)=latvec_in(1,:)
 latvect(:,2)=latvec_in(2,:)
 latvect(:,3)=latvec_in(3,:)
 call invertmat(latvect,latvectinv,3)
 flattarget=matmul(strmattarget,latvectinv)
!Get full stress matrix
 strmat(1,1)=strten_in(1)
 strmat(2,2)=strten_in(2)
 strmat(3,3)=strten_in(3)
 strmat(1,2)=strten_in(6)
 strmat(2,1)=strten_in(6)
 strmat(1,3)=strten_in(5)
 strmat(3,1)=strten_in(5)
 strmat(2,3)=strten_in(4)
 strmat(3,2)=strten_in(4)
!Convert the target stress to the target gradient, assuming cell volume to be 1
 flat=matmul(strmat,latvectinv)
!Eliminate fixed lattice components from the lattice gradient difference
 dflat=flat-flattarget
 call elim_fixed_lat(latvec_in,dflat)
!Transform back
 dstrmat=-matmul(dflat,latvect)
!Extract components
 dstr(1)=dstrmat(1,1)
 dstr(2)=dstrmat(2,2)
 dstr(3)=dstrmat(3,3)
 dstr(4)=dstrmat(2,3)
 dstr(5)=dstrmat(1,3)
 dstr(6)=dstrmat(1,2)
!Evaluate the convergence
 do istr=1,6
     if(abs(dstr(istr))*strfact >= fmax_lat ) fmax_lat=abs(dstr(istr))*strfact
 end do
endif
fmax=max(fmax_at,fmax_lat)
 iexit=0
 if(fmax.lt.tolmxf) iexit=1
!write(*,*) "FLAT;FAT", fmax_lat,fmax_at
end subroutine

!************************************************************************************

subroutine dist_ang2latvec(dist_ang,latvec,pi)
use mod_interface
!This subroutine will generate the lattice vector representation of the cell
!from the length/angle representation
implicit none
real(8):: dist_ang(6),latvec(3,3),pi,convang
convang=pi/180.d0

latvec(1,1)=dist_ang(1)
latvec(2,1)=0.d0
latvec(3,1)=0.d0
latvec(1,2)=dist_ang(2)*cos(convang*dist_ang(6))
latvec(2,2)=dist_ang(2)*sin(convang*dist_ang(6))
latvec(3,2)=0.d0
latvec(1,3)=dist_ang(3)*cos(convang*dist_ang(5))
latvec(2,3)=(dist_ang(2)*dist_ang(3)*cos(convang*dist_ang(4))-latvec(1,2)*latvec(1,3))/latvec(2,2)
latvec(3,3)=sqrt(dist_ang(3)**2-latvec(1,3)**2-latvec(2,3)**2)
end subroutine

!************************************************************************************

subroutine dist_latvec2ang(dist_ang,latvec,pi)
use mod_interface
!This subroutine will generate the angdeg represenation of the cell from the lattice vectors
implicit none
real(8):: dist_ang(6),latvec(3,3),pi,convang
convang=180.d0/pi
dist_ang(1)=sqrt(dot_product(latvec(:,1),latvec(:,1)))
dist_ang(2)=sqrt(dot_product(latvec(:,2),latvec(:,2)))
dist_ang(3)=sqrt(dot_product(latvec(:,3),latvec(:,3)))
dist_ang(4)=acos(dot_product(latvec(:,2),latvec(:,3))/(dist_ang(2)*dist_ang(3)))*convang
dist_ang(5)=acos(dot_product(latvec(:,3),latvec(:,1))/(dist_ang(3)*dist_ang(1)))*convang
dist_ang(6)=acos(dot_product(latvec(:,1),latvec(:,2))/(dist_ang(1)*dist_ang(2)))*convang
end subroutine

!************************************************************************************

subroutine fragments(latvec,xred,nfrag,xcart,fragarr,fragsize)
use mod_interface
use global, only: nat,ntypat,znucl,rcov,typat,char_type
implicit none
real(8),dimension(3,nat), INTENT(IN) :: xred
real(8):: latvec(3,3),rotmat(3,3),dproj(6)
real(8),allocatable:: pos(:,:)
integer :: nfrag, nfragold
logical :: occured,niter
real(8)::  dist, mindist, angle, vec(3), cmass(3), velcm(3), bondlength, bfactor,rnrmi,scpr
real(8):: ekin,vcm1,vcm2,vcm3,ekin0,scale,xcart(3,nat)
integer::iat, jat, nmax(1), imin(2),ifrag
integer, dimension(nat):: fragarr,fragsize(nat)

bfactor=1.1d0
fragarr(:)=0                  !Array, which atom belongs to which fragment
nfrag=0                       !Number of fragments
fragsize=0                    !Size of each fragment

!Calculate number of fragments and fragmentlist of the atoms
do
  nfragold=nfrag
  do iat=1,nat                !Check the first atom that isn't part of a cluster yet
    if(fragarr(iat)==0) then
      nfrag=nfrag+1
      fragarr(iat)=nfrag
      fragsize(nfrag)=fragsize(nfrag)+1
      xcart(:,iat)=matmul(latvec,xred(:,iat))
      exit
    endif
  enddo
  if (nfragold==nfrag) exit
7000 niter=.false.
  do iat=1,nat                !Check if all the other atoms are part of the current cluster
    do jat=1,nat
    bondlength=rcov(typat(iat))+rcov(typat(jat))
    if(nfrag==fragarr(iat) .AND. jat.ne.iat .AND. fragarr(jat)==0) then
!         call pbc_distance1(latvec,xred(:,iat),xred(:,jat),dist)
         call pbc_distance2(latvec,xred(:,iat),xcart(:,iat),xred(:,jat),xcart(:,jat),dist)
         if(dist<(bfactor*bondlength)**2) then
         fragarr(jat)=nfrag
         fragsize(nfrag)=fragsize(nfrag)+1
         niter=.true.
         endif
    endif
    enddo
  enddo
  if(niter) then
  goto 7000
  endif
enddo


open(unit=46,file="pos_fragment.ascii")
write(46,*) "Fragmentation in cartesian coordinates"
allocate(pos(3,nat))
pos=xcart
call latvec2dproj(dproj,latvec,rotmat,pos,nat)
write(46,*) dproj(1:3)
write(46,*) dproj(4:6)
do iat=1,nat
  write(46,'(3(1x,es25.15),2x,a2,1x,a1,i4)')       pos(:,iat),trim(char_type(typat(iat))),"#",fragarr(iat)
enddo
deallocate(pos)
close(46)


end subroutine fragments

!************************************************************************************

subroutine pbc_distance0(latvec,xred_1,xred_2,distance2,dxyz)
use mod_interface
!This routine computes the distance of the reduced coordinates xred_1 and xred_2 and applies 
!periodic boundary conditions to them and computes the squared distance
implicit none
integer:: i
real(8):: xred_1(3),xred_2(3),diff(3),distance2,latvec(3,3),dxyz(3)
diff=xred_2-xred_1
do i=1,3
    if(diff(i).le.-0.5d0) then
        diff(i)=diff(i)+1.d0
    elseif(diff(i).gt.0.5d0) then
        diff(i)=diff(i)-1.d0
    endif
enddo
dxyz=matmul(latvec,diff)
distance2=dot_product(dxyz,dxyz)
end subroutine

!************************************************************************************

subroutine pbc_distance1(latvec,xred_1,xred_2,distance2)
use mod_interface
!This routine computes the distance of the reduced coordinates xred_1 and xred_2 and applies 
!periodic boundary conditions to them and computes the squared distance. Minimal image convention
implicit none
integer:: i
real(8):: xred_1(3),xred_2(3),diff(3),distance2,latvec(3,3)
diff=xred_2-xred_1
do i=1,3
    if(diff(i).le.-0.5d0) then
        diff(i)=diff(i)+1.d0
    elseif(diff(i).gt.0.5d0) then
        diff(i)=diff(i)-1.d0
    endif
enddo
distance2=dot_product(matmul(latvec,diff),matmul(latvec,diff))
end subroutine

!************************************************************************************

subroutine pbc_distance2(latvec,xred_1,xcart_1,xred_2,xcart_2,distance2)
use mod_interface
!This routine computes the distance of the reduced coordinates xred_1 and xred_2 and applies 
!periodic boundary conditions to them and computes the squared distance
implicit none
integer:: i,k,l,m
real(8):: xred_1(3),xred_2(3),diff(3),distance2,latvec(3,3),xcart_1(3),xcart_2(3),xcart_tmp(3),xcart_20(3),xcart_10(3)
real(8):: xred_10(3),xred_20(3),distance2_tmp
xred_10(1)=modulo(modulo(xred_1(1),1.d0),1.d0)
xred_10(2)=modulo(modulo(xred_1(2),1.d0),1.d0)
xred_10(3)=modulo(modulo(xred_1(3),1.d0),1.d0)
xcart_10=matmul(latvec,xred_10)
xred_20(1)=modulo(modulo(xred_2(1),1.d0),1.d0)
xred_20(2)=modulo(modulo(xred_2(2),1.d0),1.d0)
xred_20(3)=modulo(modulo(xred_2(3),1.d0),1.d0)
xcart_20=matmul(latvec,xred_20)
distance2=1.d5
do k=-1,1
do l=-1,1
do m=-1,1
  xcart_tmp=xcart_20+real(k,8)*latvec(:,1)+real(l,8)*latvec(:,2)+real(m,8)*latvec(:,3)
  diff=xcart_10-xcart_tmp
  distance2_tmp=dot_product(diff,diff)
  if(distance2_tmp.lt.distance2) then
    distance2=distance2_tmp
    xcart_2=xcart_tmp+xcart_1-xcart_10
  endif
enddo
enddo
enddo
end subroutine


!************************************************************************************

!!subroutine rand_sphere(radxyz,pi)
!!!random spherical distribution
!!!This routine will compute a random uniformli distributed point on a sphere of radius 1
!!!Uniform distribution of points on a sphere:
!!!http://mathworld.wolfram.com/SpherePointPicking.html
!!implicit none
!!real(8), intent(in):: pi
!!real(8):: u,v,theta,phi,radxyz(3),sinphi,costheta,cosphi,sintheta
!!
!!
!!  call RANDOM_NUMBER(v)
!!  call RANDOM_NUMBER(u)
!!
!!  theta=u*2.d0*pi
!!  phi=acos(2.d0*v-1.d0)
!!
!!  sinphi=sin(phi)
!!  costheta=cos(theta)
!!  cosphi=cos(phi)
!!  sintheta=sin(theta)
!!
!!  radxyz(1)=costheta*sinphi
!!  radxyz(2)=sintheta*sinphi
!!  radxyz(3)=cosphi
!!end subroutine


!************************************************************************************

subroutine inertia_tensor(nat,xcart,cmass,amass,intens)
use mod_interface
!This routine computes the inertia tensor with respect to the center of mass of a system with nat atoms
implicit none
integer:: nat,iat,i,j
real(8):: xcart(3,nat),amass(nat),intens(3,3),cmass(3),xtmp(3),dist2

intens=0.d0
do iat=1,nat
xtmp=xcart(:,iat)-cmass(:)
dist2=dot_product(xtmp,xtmp)
  do i=1,3
  do j=1,i
     intens(i,j)=intens(i,j)+amass(iat)*(dist2*delta_kronecker(i,j)-xtmp(i)*xtmp(j))
  enddo
  enddo
enddo

intens(1,2)=intens(2,1)
intens(1,3)=intens(3,1)
intens(2,3)=intens(3,2)
end subroutine

!************************************************************************************

subroutine rot_ener(omega,intens,erot)
use mod_interface
!This routine computes the rotational kinetic energy of a system with known tensor of intertia and a given 
!angular velocity omega, based on its orientation and magnitude
implicit none
real(8):: omega(3),intens(3,3),erot
erot= 0.5d0*dot_product(omega,matmul(intens,omega))
end subroutine


!************************************************************************************

subroutine init_rotvels(nat,xred,latvec,temp,amass,vel)
use mod_interface
!This routine will first find the correct partitioning of the system into molecules, then assign
!rotational and translational velocities to these molecules according to the tempereature temp
use global, only: char_type,typat,units
use defs_basis, only: Bohr_Ang,pi,kb_HaK
implicit none
integer,intent(in):: nat
real(8),intent(in):: xred(3,nat),latvec(3,3),temp,amass(nat)
real(8),intent(out):: vel(3,nat)
integer:: iat,nfrag,ifrag,LWORK,INFO,idim,omega_opt
real(8):: xcart(3,nat),ekin_rot,ekin_trans,rotmat(3,3),dproj(6),angbohr,erot_tmp,ekin_tot,v2gauss,vtest,tmp(3)
real(8):: etrans_tmp,latvec_tmp(3,3),diag_inert(3,3),dir(3),weight(3),frac_rot(2)
real(8),allocatable:: cmass(:,:),vel_cmass(:,:),intens(:,:,:),fragxcart(:,:,:),fragmass(:,:),fragvel(:,:,:),masstot(:)
real(8),allocatable:: WORK(:),eval(:),frag_ekin_rot(:),omega(:,:)
integer,allocatable:: counter(:)
integer, dimension(nat):: fragarr,fragsize
!Exact total kinetic energy at given temperature
ekin_tot=3.d0*real(nat,8)*kb_HaK*temp

!Direction to bias the rotational component:
!1 Random
!2 Along lowest intertia tensor direction
!3 Weighted according to the principle inertia
!4 A weighted random distribution of the rotational direction: principle axis weighted with 1, 1/2, 1/4
!5 A weighted random distribution of the rotational direction: principle axis weighted with 1, 1/3, 1/6
!6 A weighted random distribution of the rotational direction: principle axis weighted with 1, 1/4, 1/8
!7 A weighted random distribution of the rotational direction: principle axis weighted with 1, 1/2, 1/4
!8 A weighted random distribution of the rotational direction: principle axis weighted with 1, 1/3, 1/9
!9 A weighted random distribution of the rotational direction: principle axis weighted with 1, 1/4, 1/16
omega_opt=3

weight=1.d0
SELECT CASE (omega_opt)
   CASE (2)
      weight=(/1.d0, 0.d0, 0.d0/)
   CASE (3)
      weight=(/1.d0, 1.d0, 1.d0/)
   CASE (4)
      weight=(/1.d0, 1.d0/2.d0, 1/4.d0/)
   CASE (5)
      weight=(/1.d0, 1.d0/3.d0, 1/6.d0/)
   CASE (6)
      weight=(/1.d0, 1.d0/4.d0, 1/8.d0/)
   CASE (7)
      weight=(/1.d0, 1.d0/2.d0, 1/4.d0/)
   CASE (8)
      weight=(/1.d0, 1.d0/3.d0, 1/9.d0/)
   CASE (9)
      weight=(/1.d0, 1.d0/4.d0, 1/16.d0/)
  CASE DEFAULT
      weight=1.d0
END SELECT


!Partition the kinetic energy to go into the translational and the rotational energy
frac_rot(1)=0.3d0;frac_rot(2)=0.7d0 !Lower and upper fractional boundary for rotational energy
call random_number(ekin_rot)
ekin_rot=ekin_rot*(frac_rot(2)-frac_rot(1))+frac_rot(1)
ekin_trans=(1.d0-ekin_rot)*ekin_tot
ekin_rot=ekin_rot*ekin_tot
vel=0.d0

!First identify the fragments or molecules in the cell
call fragments(latvec,xred,nfrag,xcart,fragarr,fragsize)

!Allocate the arrays corresponding to molecular quantities
allocate(cmass(3,nfrag),vel_cmass(3,nfrag),intens(3,3,nfrag))
allocate(fragxcart(3,maxval(fragsize),nfrag),fragmass(maxval(fragsize),nfrag),counter(nfrag))
allocate(frag_ekin_rot(nfrag),omega(3,nfrag),fragvel(3,maxval(fragsize),nfrag),masstot(nfrag))

!Setup the xcart array for all nfrag clusters
counter=0
do iat=1,nat
   counter(fragarr(iat))=counter(fragarr(iat))+1
   fragxcart(:,counter(fragarr(iat)),fragarr(iat))=xcart(:,iat)
   fragmass(counter(fragarr(iat)),fragarr(iat))=amass(iat)
enddo

!Write the connected molecules
open(unit=2,file="xcart.molecule.ascii")
angbohr=1.d0
if(trim(units)=="angstroem") then
  angbohr=Bohr_Ang
endif
write(2,*) nat
latvec_tmp=latvec
call latvec2dproj(dproj,latvec_tmp,rotmat,xcart,nat)
write(2,*) dproj(1:3)*angbohr
write(2,*) dproj(4:6)*angbohr
  do iat=1,nat
       write(2,'(3(1x,es25.15),2x,a2)') angbohr*xcart(:,iat),trim(char_type(typat(iat)))
  enddo
close(2)

!Compute the inertia tensor stuff
cmass=0.d0
!Now we will compute for each of the fragments the inertia tensor
do ifrag=1,nfrag
  masstot(ifrag)=0.d0
  !Compute the center of mass of every fragment
  do iat=1,fragsize(ifrag)
      masstot(ifrag)=masstot(ifrag)+fragmass(iat,ifrag)
      cmass(:,ifrag)=cmass(:,ifrag)+fragmass(iat,ifrag)*fragxcart(:,iat,ifrag)
  enddo
  cmass(:,ifrag)=cmass(:,ifrag)/masstot(ifrag)
  !Compute intertia tensor
  call inertia_tensor(fragsize(ifrag),fragxcart(:,1:fragsize(ifrag),ifrag),cmass(:,ifrag),&
       &fragmass(1:fragsize(ifrag),ifrag),intens(:,:,ifrag))
!  write(*,*) "ifrag,mass",ifrag,masstot(ifrag)
!  write(*,*) "CM", cmass(:,ifrag)
!  write(*,*) intens(:,1,ifrag)
!  write(*,*) intens(:,2,ifrag)
!  write(*,*) intens(:,3,ifrag)
!LWORK=-1
!allocate(WORK(1),eval(3))
!        call DSYEV('V','L',3,intens(:,:,ifrag),3,eval,WORK,LWORK,INFO)
!        if (info.ne.0) stop 'DSYEV in correct_hessin'
!LWORK=WORK(1)
!deallocate(WORK)
!allocate(WORK(LWORK))
!        call DSYEV('V','L',3,intens(:,:,ifrag),3,eval,WORK,LWORK,INFO)
!        if (info.ne.0) stop 'DSYEV in correct_hessin'
!        write(*,*) '---   App. eigenvalues in a.u. -------------'
!        do iat=1,3
!         write(*,'(1x,es25.15)') eval(iat)
!        enddo
!deallocate(work,eval)
end do

!Assign the target kinetic rotational energies of each fragment
!And aso the rotational exis, not scaled yet omega
fragvel=0.d0
call random_number(frag_ekin_rot)
!write(*,*) "RANDOM",frag_ekin_rot
!Scale to 1
frag_ekin_rot=frag_ekin_rot/sum(frag_ekin_rot)*ekin_rot
do ifrag=1,nfrag
!Only run if molecule is larger than 1 atom
if(fragsize(ifrag).gt.2) then
!Set omega
   if(omega_opt==1) then
!Either random direction
   call rand_sphere(omega(:,ifrag),pi)
!... or along the smallest inertia tensor  
   elseif(omega_opt.le.9) then
         !Diagonalize the tensor
          diag_inert(:,:)=intens(:,:,ifrag)
          allocate(WORK(1),eval(3))
          LWORK=-1
          call DSYEV('V','L',3,diag_inert,3,eval,WORK,LWORK,INFO)
          if (info.ne.0) stop 'DSYEV in init_rotvels'
          LWORK=WORK(1)
          deallocate(WORK)
          allocate(WORK(LWORK))
          call DSYEV('V','L',3,diag_inert,3,eval,WORK,LWORK,INFO)
          if (info.ne.0) stop 'DSYEV in correct_hessin'
!Determine direction of rotation, meaning the sign of the vector, and the weight
              write(*,'(a,i5,3(es25.15))') "Weight e  : ",ifrag,eval
              write(*,'(a,i5,3(es25.15))') "Weight w  : ",ifrag,weight
              call random_number(dir)
              dir=(dir-0.5d0)*2.d0
              write(*,'(a,i5,3(es25.15))') "Weight r  : ",ifrag,dir
              dir=dir*weight
              write(*,'(a,i5,3(es25.15))') "Weight rw : ",ifrag,dir
              dir(2)=dir(2)*eval(1)/eval(2);dir(3)=dir(3)*eval(1)/eval(3)
              write(*,'(a,i5,3(es25.15))') "Weight rwe: ",ifrag,dir
              omega(:,ifrag)=diag_inert(:,1)*dir(1)+diag_inert(:,2)*dir(2)+diag_inert(:,3)*dir(3)
          deallocate(work,eval)
          else
              stop "Wrong option for omega_opt"
    endif
       
       
          
!Set the correct rotational velocity omega
   call rot_ener(omega(:,ifrag),intens(:,:,ifrag),erot_tmp)
!   write(*,*) "Init",frag_ekin_rot(ifrag),erot_tmp
   omega(:,ifrag)=omega(:,ifrag)*sqrt(frag_ekin_rot(ifrag)/erot_tmp)
   call rot_ener(omega(:,ifrag),intens(:,:,ifrag),erot_tmp)
!   write(*,*) "Scaled",frag_ekin_rot(ifrag),erot_tmp
!Assign rotational velocities to the atoms in the fragments
   call assign_vel(fragsize(ifrag),fragxcart(:,1:fragsize(ifrag),ifrag),&
        &cmass(:,ifrag),omega(:,ifrag),fragvel(:,1:fragsize(ifrag),ifrag))
!Recompute v2gauss
         v2gauss=0.d0
         vtest=0.d0
         do iat=1,fragsize(ifrag)
           do idim=1,3
             v2gauss=v2gauss+fragvel(idim,iat,ifrag)**2*amass(iat)
             vtest=vtest+fragvel(idim,iat,ifrag)/(3.d0*real(nat,8))
           end do
         end do
!   write(*,*) "TEST", 0.5d0*v2gauss,vtest
endif
enddo

!Assign translational kinetic energy
!Get random Gaussian distributed atomic velocities
call gausdist(nfrag,vel_cmass,masstot)
call elim_moment(nfrag,vel_cmass,masstot)
!Set the correct translational velocity
etrans_tmp=0.d0
do ifrag=1,nfrag
  etrans_tmp=etrans_tmp+0.5d0*dot_product(vel_cmass(:,ifrag),vel_cmass(:,ifrag))*masstot(ifrag)
enddo
!   write(*,*) "Init",ekin_trans,etrans_tmp
   vel_cmass(:,:)=vel_cmass(:,:)*sqrt(ekin_trans/etrans_tmp)
etrans_tmp=0.d0
do ifrag=1,nfrag
  etrans_tmp=etrans_tmp+0.5d0*dot_product(vel_cmass(:,ifrag),vel_cmass(:,ifrag))*masstot(ifrag)
enddo
!   write(*,*) "Scaled",ekin_trans,etrans_tmp
   !Recompute v2gauss
         v2gauss=0.d0
         vtest=0.d0
         do iat=1,nfrag
           do idim=1,3
             v2gauss=v2gauss+vel_cmass(idim,iat)**2*masstot(iat)
             vtest=vtest+vel_cmass(idim,iat)/(3.d0*real(nfrag,8))
           end do
         end do
!   write(*,*) "TEST", 0.5d0*v2gauss,vtest

!Distribute the CM velocities to each cluster
do ifrag=1,nfrag
   do iat=1,fragsize(ifrag)
   fragvel(:,iat,ifrag)=fragvel(:,iat,ifrag)+vel_cmass(:,ifrag)
   enddo
end do

!Transform back to the current velocity structure
counter=0
do iat=1,nat
   counter(fragarr(iat))=counter(fragarr(iat))+1
   vel(:,iat)=fragvel(:,counter(fragarr(iat)),fragarr(iat))
enddo

!Check total energy
etrans_tmp=0.d0
do iat=1,nat
   etrans_tmp=etrans_tmp+0.5d0*dot_product(vel(:,iat),vel(:,iat))*amass(iat)
enddo

!write(*,*) "TOTAL", etrans_tmp,ekin_tot

deallocate(cmass,vel_cmass,intens,fragxcart,fragmass,counter)
deallocate(frag_ekin_rot,omega,fragvel,masstot)

end subroutine

!************************************************************************************

subroutine assign_vel(nat,xcart,cmass,omega,vel)
use mod_interface
!This routine will take the cartesian coordinates and the center of mass and compute the 
!velocities of every atom in the system with nat atoms
implicit none
integer:: nat,iat
real(8):: xcart(3,nat),cmass(3),xtmp(3),vel(3,nat),omega(3)
do iat=1,nat
 xtmp=xcart(:,iat)-cmass(:)
 call cross_product(omega,xtmp,vel(:,iat))
enddo
end subroutine


!************************************************************************************

subroutine refragment(fragarr,nat)
use mod_interface
!This routine will rearragne the integer array fragarr such that the fragment indexes are well 
!assigned in ascending order, and new array indexes are assigned to atoms with 
!values <0
implicit none
integer:: fragarr(nat),nat,iat,jat,cnt,find,fragarr_tmp(nat)
cnt=0
fragarr_tmp=0
if(all(fragarr.le.0)) then
 do iat=1,nat
    fragarr(iat)=iat
 enddo
else
!first assign new values to existing fragments
do iat=1,nat
   if(fragarr(iat).gt.0.and.fragarr_tmp(iat)==0) then
      cnt=cnt+1
      do jat=1,nat 
        if(fragarr(iat)==fragarr(jat)) fragarr_tmp(jat)=cnt
      enddo
   endif
enddo
!Now all unassigned clusters
do iat=1,nat
  if(fragarr(iat).le.0) then
    cnt=cnt+1
    fragarr_tmp(iat)=cnt
  endif 
enddo
fragarr=fragarr_tmp
endif
end subroutine

!************************************************************************************

subroutine make_linked_list(fragarr,fragsize,lhead,llist,nat,nmol)
use mod_interface
!This subroutine will create a linked list to represent molecules.
!lhead is an array of length nat, of which only nmol will be significant. Their entries point to a
!position in array lhead, the first atom in the molecule, 
!llist_at is an array with entries linking the atoms in the molecule. Termination flag is 0
implicit none
integer:: fragarr(nat),nat,iat,nmol,ifrag
integer:: lhead(nmol),llist(nat),fragsize(nmol)
llist=0
lhead=0
do iat=1,nat
   llist(iat)=lhead(fragarr(iat)) 
   lhead(fragarr(iat))=iat
enddo
call get_fragsize(fragsize,lhead,llist,nat,nmol)
!do ifrag=1,nmol
!   iat=lhead(ifrag)
!   do while (iat.ne.0)
!     write(*,*) ifrag, iat
!     iat=llist(iat)
!   enddo
!enddo
end subroutine

!************************************************************************************

subroutine get_fcm_torque(fcm,torque,fcart,quat,xcart_mol,lhead,llist,nat,nmol)
use mod_interface
!Computes the total force on molecules and the torques,assuming xcart_mol has molecules with CM at origin
implicit none
integer:: nat,nmol,iat,ifrag,llist(nat),lhead(nmol)
real(8),intent(in):: fcart(3,nat),quat(4,nmol),xcart_mol(3,nat)
real(8):: fcm(3,nmol),torque(3,nmol),crossp(3),xtmp(3),rotmat(3,3)
do ifrag=1,nmol
   call quat2rotmat(rotmat,quat(:,ifrag))
   iat=lhead(ifrag)
   fcm(:,ifrag)=0.d0
   torque(:,ifrag)=0.d0
   do while (iat.ne.0)
     fcm(:,ifrag)=fcm(:,ifrag)+fcart(:,iat)   
     xtmp=matmul(rotmat,xcart_mol(:,iat))
     call cross_product(xtmp,fcart(:,iat),crossp)
     torque(:,ifrag)=torque(:,ifrag)+crossp
     iat=llist(iat)
   enddo
enddo
end subroutine

!************************************************************************************

subroutine expand_rigid(latvec,xred_cm,quat,xcart_mol,lhead,llist,nat,nmol,xred_in)
use mod_interface
!This routine will produce the desired, expanded form of xred_in which can be directly fed into 
!the single point energy routine. Input are quat and the xred_cm, the center of masses of each molecule
!in the cell
implicit none
real(8),intent(in):: latvec(3,3),xred_cm(3,nmol),quat(4,nmol),xcart_mol(3,nat)
real(8),intent(out):: xred_in(3,nat)
real(8):: rotmat(3,3),xcart_tmp(3,nat),cmass(3,nmol)
integer:: nat,nmol,iat,imol,llist(nat),lhead(nmol)
call rxyz_int2cart(latvec,xred_cm,cmass,nmol)
do imol=1,nmol
   call quat2rotmat(rotmat,quat(:,imol))
   iat=lhead(imol)
   do while (iat.ne.0)
     xcart_tmp(:,iat)=matmul(rotmat,xcart_mol(:,iat))
     xcart_tmp(:,iat)=xcart_tmp(:,iat)+cmass(:,imol)
     iat=llist(iat)
   enddo
enddo
!Produce reduced coordinates
call rxyz_cart2int(latvec,xred_in,xcart_tmp,nat)
end subroutine expand_rigid

!************************************************************************************

subroutine init_cm_mol(latvec,xred,xcart_shifted,xred_cm,quat,amass,masstot,intens,inprin,inaxis,lhead,llist,nat,nmol)
use mod_interface
!This routine will get the cm and shift all molecular units into xcart_shifted
!and write, for each molecule, an xyz file containing the shifted molecular unit
!We will also express the center of masses in reduced coordinates (xred_cm), and of
!course intiallize the orientation vectors quat to 0
use defs_basis, only: Ha_eV,Bohr_Ang
use global, only: char_type,ntypat,typat,units
implicit none
real(8),intent(in):: latvec(3,3),xred(3,nat)
integer:: nat,nmol,iat,llist(nat),lhead(nmol),fragsize(nmol),imol,jmol,kmol
real(8):: xcart_in(3,nat),xcart_shifted(3,nat),cmass(3,nmol),amass(nat)
real(8):: masstot(nmol),angbohr,quat(4,nmol),xred_cm(3,nmol),xcart_tmp(3,nat)
real(8):: circular(3,3),tol,rot_c(3,3),rot_all(3,3),inprin(3,nmol),intens(3,3,nmol)
real(8):: inaxis(3,3,nmol),ident(3,3),tmp(3,nmol),quat_tmp(4),tmp_real(4),tmp_mat(3,3)
character(5):: fn
logical:: symtop(nmol),tmp_logical
tol=1.d-10
!Convert to cartesian coordinates
call rxyz_int2cart(latvec,xred,xcart_in,nat)
!Get the center of masses and save it in xred_cm
call get_cmass(cmass,masstot,xcart_in,amass,lhead,llist,nat,nmol)
call rxyz_cart2int(latvec,xred_cm,cmass,nmol)
!Compute the fragment sizes
call get_fragsize(fragsize,lhead,llist,nat,nmol)
!Get inertia tensor
call get_inertia_tensor(intens,inprin,inaxis,cmass,xcart_in,amass,lhead,llist,nat,nmol)
!Circular matrix that will rotate xyz to zxy
circular(1,:)=(/0.d0,0.d0,1.d0/)
circular(2,:)=(/1.d0,0.d0,0.d0/)
circular(3,:)=(/0.d0,1.d0,0.d0/)
!Identity matrix
ident=0.d0
ident(1,1)=1.d0
ident(2,2)=1.d0
ident(3,3)=1.d0
!Initiallize the orientation vector
quat=0.d0
!Compute symmetric tops and rotate accordingly
symtop=.false.
rot_c=ident
   write(*,*) inprin
do imol=1,nmol
 if(fragsize(imol).ge.2) then
   if(abs(inprin(3,imol)-inprin(1,imol)).lt.tol) then
     rot_c=circular
     symtop(imol)=.true.
   elseif(abs(inprin(2,imol)-inprin(3,imol)).lt.tol)then  
     rot_c=matmul(circular,circular)
     symtop(imol)=.true.
   elseif(abs(inprin(1,imol)-inprin(2,imol)).lt.tol)then
     symtop(imol)=.true.
   endif
   rot_all=matmul(rot_c,transpose(inaxis(:,:,imol)))
   call rotmat2quat(rot_all,quat_tmp) 
   call qinv(quat_tmp,quat(:,imol))
  !Rotate tha shit
  !do imol=1,nmol
   inaxis(:,:,imol)=matmul(rot_c,inaxis(:,:,imol))
   inaxis(:,:,imol)=matmul(inaxis(:,:,imol),transpose(inaxis(:,:,imol)))
   inprin(:,imol)=matmul(rot_c,inprin(:,imol))
   if(symtop(imol)) then !Average the principle axis
       inprin(1,imol)=0.5d0*(inprin(1,imol)+inprin(2,imol))
       inprin(2,imol)=inprin(1,imol)
   endif
 else
   symtop(imol)=.true.
 endif
   iat=lhead(imol)
   do while (iat.ne.0)
     xcart_shifted(:,iat)=xcart_in(:,iat)-cmass(:,imol)
     xcart_shifted(:,iat)=matmul(rot_all,xcart_shifted(:,iat)) 
     iat=llist(iat)
   enddo
!enddo
enddo
!Rearrange the indexes of lhead, such that the first molecules are arbitrary, then such with symmetric tops, then such with 2 atoms, then with 1 atom
do imol=1,nmol
do jmol=2,nmol
 if(((fragsize(jmol).gt.fragsize(jmol-1)).and.(symtop(jmol-1).and..not.symtop(jmol))).or.&      !To move the smallest molecules down the list
   &((fragsize(jmol).gt.fragsize(jmol-1)).and.(symtop(jmol-1).and.symtop(jmol)))) then !For the case that we have other symmetric molecules
   iat=              lhead(jmol-1);      lhead(jmol-1)=     lhead(jmol);      lhead(jmol)=iat
   tmp_real=        quat(:,jmol-1);     quat(:,jmol-1)=    quat(:,jmol);     quat(:,jmol)=tmp_real
   iat=           fragsize(jmol-1);   fragsize(jmol-1)=  fragsize(jmol);   fragsize(jmol)=iat
   tmp_logical=     symtop(jmol-1);     symtop(jmol-1)=    symtop(jmol);     symtop(jmol)=tmp_logical
   tmp_real(1:3)=  cmass(:,jmol-1);    cmass(:,jmol-1)=   cmass(:,jmol);    cmass(:,jmol)=tmp_real(1:3)
   tmp_mat=intens     (:,:,jmol-1); intens(:,:,jmol-1)=intens(:,:,jmol); intens(:,:,jmol)=tmp_mat
   tmp_mat=inaxis     (:,:,jmol-1); inaxis(:,:,jmol-1)=inaxis(:,:,jmol); inaxis(:,:,jmol)=tmp_mat
   tmp_real(1:3)= inprin(:,jmol-1);   inprin(:,jmol-1)=  inprin(:,jmol);   inprin(:,jmol)=tmp_real(1:3)  
   tmp_real(1)=    masstot(jmol-1);    masstot(jmol-1)=   masstot(jmol);    masstot(jmol)=tmp_real(1)
   tmp_real(1:3)=xred_cm(:,jmol-1);  xred_cm(:,jmol-1)= xred_cm(:,jmol);  xred_cm(:,jmol)=tmp_real(1:3)
 endif  
enddo
enddo
write(*,*) symtop
!Get inertia tensor
xcart_tmp=0.d0
!call get_inertia_tensor(intens,inprin,inaxis,xcart_tmp,xcart_shifted,amass,lhead,llist,nat,nmol)
       do iat=1,nmol
          write(122,*) "Fragsize",fragsize(iat)
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

if(trim(units)=="angstroem") then
  angbohr=Bohr_Ang
elseif(trim(units)=="bohr") then
  angbohr=1.d0
else
  stop "Wrong file unit format"
endif

!!do imol=1,nmol
!!iat=lhead(imol)
!!   call quat2rotmat(rot_all,quat(:,imol)) 
!!   do while (iat.ne.0)
!!     write(222,'(a2,2x,3(es25.15))') trim(char_type(typat(iat))),(cmass(:,imol)+matmul(rot_all,xcart_shifted(:,iat)))*angbohr
!!     iat=llist(iat)
!!   enddo
!!enddo

do imol=1,nmol
   write(fn,'(i5.5)') imol
   open(unit=22,file="posfrag"//fn//".xyz")
   write(22,*)  fragsize(imol), trim(units)
   write(22,*)  
   iat=lhead(imol)
   do while (iat.ne.0)
     write(22,'(a2,2x,3(es25.15))') trim(char_type(typat(iat))),xcart_shifted(:,iat)*angbohr
     iat=llist(iat)
   enddo
   close(22)
enddo
end subroutine
!contains

!************************************************************************************

subroutine get_cmass(cmass,masstot,xcart,amass,lhead,llist,nat,nmol)
!This subroutine will compute the center of mass and the total mass of the fragment
!given the cartesian coordinates xcart and the linked cell lists and atomic masses
implicit none
integer:: nat,nmol,iat,ifrag,lhead(nmol),llist(nat)
real(8):: xcart(3,nat),amass(nat),masstot(nmol),cmass(3,nmol)
cmass=0.d0
do ifrag=1,nmol
  masstot(ifrag)=0.d0
  !Compute the center of mass of every fragment
   iat=lhead(ifrag)
   do while (iat.ne.0)
        masstot(ifrag)=masstot(ifrag)+amass(iat)
        cmass(:,ifrag)=cmass(:,ifrag)+amass(iat)*xcart(:,iat)
        iat=llist(iat)
    enddo
  cmass(:,ifrag)=cmass(:,ifrag)/masstot(ifrag)
enddo
end subroutine

!************************************************************************************

subroutine get_inertia_tensor(intens,inprin,inaxis,cmass,xcart,amass,lhead,llist,nat,nmol)
use mod_interface
!This routine will compute the intertia tensors of all molecules involved
!given the center of mass and atomic masses, the principle  moments of inertia inprin,
!and the axis of the inertia tensor for each molecule
use global, only: fragsize
implicit none
integer:: nat,nmol,iat,ifrag,i,j,llist(nat),lhead(nmol),LWORK,info
real(8):: xcart(3,nat),amass(nat),cmass(3,nmol),intens(3,3,nmol),dist2,xtmp(3)
real(8):: inprin(3,nmol),inaxis(3,3,nmol),diag_inert(3,3),tmp_vec(3),tmp_val
real(8),allocatable:: work(:),eval(:)
!Compute the tensor of intertia
do ifrag=1,nmol
   iat=lhead(ifrag)
   intens(:,:,ifrag)=0.d0
   do while (iat.ne.0)
       xtmp=xcart(:,iat)-cmass(:,ifrag)
       dist2=dot_product(xtmp,xtmp)
         do i=1,3
         do j=1,i
            intens(i,j,ifrag)=intens(i,j,ifrag)+amass(iat)*(dist2*delta_kronecker(i,j)-xtmp(i)*xtmp(j))
         enddo
         enddo
       iat=llist(iat)
   enddo
   intens(1,2,ifrag)=intens(2,1,ifrag)
   intens(1,3,ifrag)=intens(3,1,ifrag)
   intens(2,3,ifrag)=intens(3,2,ifrag)
enddo

!Compute the priciple axis and the principle moments of inertia
!Get the correct array size
diag_inert(:,:)=0.d0
allocate(WORK(1),eval(3))
LWORK=-1
call DSYEV('V','L',3,diag_inert,3,eval,WORK,LWORK,INFO)
if (info.ne.0) stop 'DSYEV in get_inertia_tensor'
LWORK=WORK(1)
deallocate(WORK)
allocate(WORK(LWORK))

do ifrag=1,nmol
if(fragsize(ifrag)==1) then !We have a single atom
  inprin(:,ifrag)=0.d0
  inaxis(:,1,ifrag)=(/1.d0,0.d0,0.d0/)
  inaxis(:,2,ifrag)=(/0.d0,1.d0,0.d0/)
  inaxis(:,3,ifrag)=(/0.d0,0.d0,1.d0/)
else
!Diagonalize the tensor
          diag_inert(:,:)=intens(:,:,ifrag)
          call DSYEV('V','L',3,diag_inert,3,eval,WORK,LWORK,INFO)
          if (info.ne.0) stop 'DSYEV in get_inertia_tensor'
          if(eval(1).lt.-1.d-10) stop "Negative diagonal element in inertia tensor!!!"
          inaxis(:,1,ifrag)=sign(1.d0,eval(1))*diag_inert(:,1)
          inaxis(:,2,ifrag)=sign(1.d0,eval(2))*diag_inert(:,2)
          inaxis(:,3,ifrag)=sign(1.d0,eval(3))*diag_inert(:,3)
          inprin(1,ifrag)=abs(eval(1))
          inprin(2,ifrag)=abs(eval(2))
          inprin(3,ifrag)=abs(eval(3))
endif
enddo
deallocate(work,eval)
end subroutine


!************************************************************************************

subroutine get_fragsize(fragsize,lhead,llist,nat,nmol)
use mod_interface
implicit none
integer:: nat,nmol,iat,ifrag,llist(nat),lhead(nmol),fragsize(nmol)
fragsize=0
do ifrag=1,nmol
   iat=lhead(ifrag)
   do while (iat.ne.0)
     fragsize(ifrag)=fragsize(ifrag)+1
     iat=llist(iat)
   enddo
enddo
end subroutine

!************************************************************************************

subroutine rbmd_symasym_s1(T_t,L_t,dt,L_til_t)
use mod_interface
!J.Chem.Phys. 110,7,3291,Matubayasi
!Symmetric top: we assume Ix=Iy
!First step: equation 15
!We construct the vector of L^tilde at T=t: L_til_t
!Input is the L_t and T_t, the angular momentum and
!the torque at time T=t
implicit none
real(8),intent(in) :: T_t(3), L_t(3),dt
real(8),intent(out):: L_til_t(3)
!integer:: nmol
L_til_t=L_t+0.5d0*dt*T_t
end subroutine

!************************************************************************************

subroutine rbmd_sym_s23(Inprin,L_til_t,dt,L_til_t5,L_til_t10)
!J.Chem.Phys. 110,7,3291,Matubayasi
!Symmetric top: we assume Ix=Iy
!Second and third steps: equation 17 and 18
!We construct the vector of L^tilde at T=t+0.5dt: L_til_t5
!Then, We construct the vector of L^tilde at T=t+dt: L_til_t10
!Input is the principal inertia tensor Inprin, symmetric Ix=Iz, and L^tilde at T=t, L_til_t
use mod_interface
implicit none
real(8),intent(in) :: Inprin(3), L_til_t(3),dt
real(8),intent(out):: L_til_t10(3),L_til_t5(3)
real(8):: Ixz,Iangle,ca,sa
!integer:: nmol,imol
!do imol=1,nmol
  Ixz=1.d0/Inprin(1)-1.d0/Inprin(2)
  Iangle=0.5*dt*Ixz*L_til_t(3)
  ca=cos(Iangle)
  sa=sin(Iangle)

!Eq.17
  L_til_t5(1)=ca*L_til_t(1)-sa*L_til_t(2)
  L_til_t5(2)=sa*L_til_t(1)-sa*L_til_t(2)
  L_til_t5(3)=   L_til_t(3)
!Eq.18
  L_til_t10(1)=ca*L_til_t5(1)-sa*L_til_t5(2)
  L_til_t10(2)=sa*L_til_t5(1)-sa*L_til_t5(2)
  L_til_t10(3)=    L_til_t(3)
!enddo
end subroutine

!************************************************************************************

subroutine rbmd_symasym_s4(Inprin,L_til_t5,quat_t,dt,quat_t10)
!J.Chem.Phys. 110,7,3291,Matubayasi
!Symmetric top: we assume Ix=Iy
!Fourth step: equations 3 and 9/10
!We construct the new quaternion at time T=t+dt: quat_t10
!We thus need to first obtain omega_tilde at T=t+0.5*dt by using L_tilde at T=t+0.5*dt
use mod_interface, except_this_one=>norm
implicit none
real(8),intent(in) :: Inprin(3),L_til_t5(3),dt,quat_t(4)
real(8),intent(out):: quat_t10(4)
real(8):: omega_til_t5(3),expmat(4,4),norm
real(8):: ca,sa,angle
!real(8),dimension(4,4):: A_omega
!integer:: nmol,imol

!do imol=1,nmol
!Construct omega according to equation 3
  omega_til_t5(1)=L_til_t5(1)/Inprin(1)
  omega_til_t5(2)=L_til_t5(2)/Inprin(2)
  omega_til_t5(3)=L_til_t5(3)/Inprin(3)
!Do the crazy stuff at equation 10
  expmat=0.d0
  norm=sqrt(omega_til_t5(1)+omega_til_t5(2)+omega_til_t5(3))  
  angle=0.5d0*dt*norm
  ca=cos(angle)
  sa=sin(angle)
  expmat(1,1)=ca
  expmat(2,2)=ca
  expmat(3,3)=ca
  expmat(4,4)=ca
  expmat=expmat+sa/norm*A_omega(omega_til_t5(1:3))
  quat_t10(:)=matmul(expmat,quat_t)
!enddo
end subroutine rbmd_symasym_s4
!contains
!end subroutine

function A_omega(omega)
real(8),dimension(4,4):: A_omega
real(8)::omega(3)
A_omega=0.d0
A_omega(2,1)=omega(1)
A_omega(3,1)=omega(2)
A_omega(4,1)=omega(3)
A_omega(3,2)=-omega(3)
A_omega(4,2)=omega(2)
A_omega(4,3)=-omega(1)

A_omega(1,2)=-A_omega(2,1)
A_omega(1,3)=-A_omega(3,1)
A_omega(1,4)=-A_omega(4,1)
A_omega(2,3)=-A_omega(3,2)
A_omega(2,4)=-A_omega(4,2)
A_omega(3,4)=-A_omega(4,3)
end function A_omega

!************************************************************************************

subroutine rbmd_symasym_s5(T_t10,L_til_t10,dt,L_t10)
!J.Chem.Phys. 110,7,3291,Matubayasi
!Symmetric top: we assume Ix=Iy
!Fifth  step, fourth step according to the manuscript: equation 19
!We construct the vector of L at T=t+dt: L_t10
!Input is the L_til_t10 (L^tilde at T=t+dt) and the torque at T=t+dt: T_t10
use mod_interface
implicit none
real(8),intent(in) :: T_t10(3), L_til_t10(3),dt
real(8),intent(out):: L_t10(3)
!integer:: nmol
L_t10=L_til_t10+0.5d0*dt*T_t10
end subroutine

!************************************************************************************

subroutine rbmd_asym_s23(Inprin,L_til_t,dt,L_til_t5,L_til_t10)
!J.Chem.Phys. 110,7,3291,Matubayasi
!Asymmetric top: Ix.ne.Iy.ne.Iz
!Second and third steps: equation 21 and 22
!We construct the vector of L^tilde at T=t+0.5dt: L_til_t5
!Then, We construct the vector of L^tilde at T=t+dt: L_til_t10
!Input is the principal inertia tensor Inprin, symmetric Ix=Iz, and L^tilde at T=t, L_til_t
implicit none
real(8),intent(in) :: Inprin(3), L_til_t(3),dt
real(8),intent(out):: L_til_t10(3),L_til_t5(3)
real(8):: Izy,Izx,Ianglezy,Ianglezx,cazy,sazy,cazx,sazx,Lz1,Lz2
!integer:: nmol,imol
!do imol=1,nmol
  Izy=1.d0/Inprin(3)-1.d0/Inprin(2)
  Izx=1.d0/Inprin(3)-1.d0/Inprin(1)
  Ianglezy=0.5*dt*Izy*L_til_t(2)
  cazy=cos(Ianglezy)
  sazy=sin(Ianglezy)
!Eq.21
  L_til_t5(1)= cazy*L_til_t(1)-sazy*L_til_t(3)
  Lz1        =-sazy*L_til_t(1)+cazy*L_til_t(3)
  Ianglezx=0.5*dt*Izx*L_til_t5(1)
  cazx=cos(Ianglezx)
  sazx=sin(Ianglezx)
  L_til_t5(2)= cazx*L_til_t(2)-sazx*Lz1
  L_til_t5(3)= sazx*L_til_t(2)+cazx*Lz1
!Eq.22
  L_til_t10(2)=   cazx*L_til_t5(2)-sazx*L_til_t5(3)
  Lz2         =  -sazx*L_til_t5(2)+cazx*L_til_t5(3)
  Ianglezy=0.5*dt*Izy*L_til_t10(2)
  cazy=cos(Ianglezy)
  sazy=sin(Ianglezy)
  L_til_t10(1)= cazy*L_til_t5(1)+sazy*Lz2
  L_til_t10(3)=-sazy*L_til_t5(1)+cazy*Lz2
!enddo
end subroutine

!************************************************************************************

subroutine rbmd_driver(quat_t,T_t,L_t,quat_t10,T_t10,L_t10,dt,inprin,&
           &fragsize,symtop,nmol)
!Currently implemented only to propagate the assymtric tops with fragsize.ge.3
use mod_interface
implicit none
real(8),intent(in) :: inprin(3,nmol),L_t(3,nmol),T_t(3,nmol),quat_t(4,nmol),dt,T_t10(3,nmol)
integer,intent(in) :: nmol,fragsize(nmol)
logical,intent(in) :: symtop(nmol)
real(8),intent(out):: L_t10(3,nmol),quat_t10(4,nmol)
real(8) :: L_til_t(3,nmol),L_til_t5(3,nmol),L_til_t10(3,nmol)
integer :: imol
do imol=1,nmol
    if(.not.symtop(imol).and.fragsize(imol).ge.3) then
    call rbmd_symasym_s1(T_t(:,imol),L_t(:,imol),dt,L_til_t(:,imol))
    call rbmd_asym_s23(Inprin(:,imol),L_til_t(:,imol),dt,L_til_t5(:,imol),L_til_t10(:,imol))
    call rbmd_symasym_s4(Inprin(:,imol),L_til_t5(:,imol),quat_t(:,imol),dt,quat_t10(:,imol))
    call rbmd_symasym_s5(T_t10(:,imol),L_til_t10(:,imol),dt,L_t10(:,imol))
    endif
enddo
end subroutine


!end program MHMoF
!****************************************************************************************************************   

 subroutine find_kpt(k1, k2, k3, lat, gridden)
! This code will define the KPT mesh based on the desired grid density
   use mod_interface
   implicit none
   integer, intent(out) :: k1,k2,k3
   real(8), intent(in)  :: lat(3,3), gridden
   
   integer :: i, j
   real(8) :: lat1, lat2, lat3
   real(8) :: angles(3), cos_arr(3)
   real(8) :: glat(3,3), crossp(3), vol, glen(3), a(3,3)
   real(8) :: pi
   
   pi = acos(-1.d0)
   
   lat1 = sqrt(lat(1,1)**2 + lat(2,1)**2 + lat(3,1)**2)
   lat2 = sqrt(lat(1,2)**2 + lat(2,2)**2 + lat(3,2)**2)
   lat3 = sqrt(lat(1,3)**2 + lat(2,3)**2 + lat(3,3)**2)
   
   cos_arr(1) = (lat(1,2)*lat(1,3) + lat(2,2)*lat(2,3) + lat(3,2)*lat(3,3))/(lat2*lat3)
   cos_arr(2) = (lat(1,1)*lat(1,3) + lat(2,1)*lat(2,3) + lat(3,1)*lat(3,3))/(lat1*lat3)
   cos_arr(3) = (lat(1,1)*lat(1,2) + lat(2,1)*lat(2,2) + lat(3,1)*lat(3,2))/(lat1*lat2)
   
   angles(1) = (acos(cos_arr(1))/pi)*180.d0
   angles(2) = (acos(cos_arr(2))/pi)*180.d0
   angles(3) = (acos(cos_arr(3))/pi)*180.d0
   
   call getvol(lat,vol)
   
   call cross_product(lat(:,2), lat(:,3), crossp(:))
   glat(:,1) = 2.d0*pi*crossp(:)/vol
   call cross_product(lat(:,3), lat(:,1), crossp(:))
   glat(:,2) = 2.d0*pi*crossp(:)/vol
   call cross_product(lat(:,1), lat(:,2), crossp(:))
   glat(:,3) = 2.d0*pi*crossp(:)/vol
   
   !Compute the correct kpts
   glen(1) = sqrt(glat(1,1)**2 + glat(2,1)**2 + glat(3,1)**2)
   glen(2) = sqrt(glat(1,2)**2 + glat(2,2)**2 + glat(3,2)**2)
   glen(3) = sqrt(glat(1,3)**2 + glat(2,3)**2 + glat(3,3)**2)
   
   call track_kpt(gridden, glen(1), k1)
   call track_kpt(gridden, glen(2), k2)
   call track_kpt(gridden, glen(3), k3)
   
 end subroutine find_kpt
 !contains
   
   subroutine track_kpt(gridden, glen, kpt)
     use mod_interface
     implicit none
     real(8), intent(in) :: gridden, glen
     integer :: kpt,j
     real(8) :: d_test
   real(8) :: pi
   
   pi = acos(-1.d0)
     
     kpt = int(glen/(gridden*2.d0*pi))
     if (kpt == 0) kpt = 1
     d_test=glen/(kpt*2.d0*pi)
     if (d_test.ge.gridden) then
       do j = 1, 25
         kpt = kpt + j
         d_test = glen/(kpt*2.d0*pi)
         if (d_test.le.gridden) exit
       enddo
     endif
   end subroutine track_kpt
   

!***

!************************************************************************************
subroutine MD_MHM_ROT(parini,parres,latvec_in,xred_in,xred_cm_in,xcart_mol,quat_in,fcart_in,strten_in,&
                      &vel_in,vel_cm_in,vel_lat_in,l_in,vvol_in,etot_in,&
                      &masstot,intens,inprin,inaxis,lhead,llist,nmol,iprec,counter,folder)
 use mod_interface
 use global, only: target_pressure_habohr,target_pressure_gpa,nat,ntypat,znucl,amu,amutmp,typat
 use global, only: char_type,units,usewf_md
 use global, only: fixat,fixlat
 use defs_basis
 use interface_code
 use mod_parini, only: typ_parini
implicit none
 type(typ_parini), intent(in):: parini
 type(typ_parini), intent(inout):: parres
 real(8) ::latvec_in(3,3),xred_in(3,nat),vel_in(3,nat),fcart_in(3,nat),etot_in,strten_in(6),vel_lat_in(3,3),vvol_in
 real(8):: amass(nmol)
 real(8):: pressure_ener
 real(8):: pressure_md(3,3)
 real(8):: pressure
 real(8):: vol
 real(8):: vol0
 real(8):: volprev
!Rigid stuff
 real(8),dimension(3,nmol):: xred_cm_in
 real(8),dimension(3,nmol):: fcart_cm
 real(8),dimension(3,nmol):: torque
 real(8),dimension(3,nat) :: xcart_mol
 real(8),dimension(3,nmol):: l_in
 real(8),dimension(3,nmol):: vel_cm_in
 real(8),dimension(4,nmol):: quat_in
 real(8),dimension(4,nmol):: quatcur
 real(8),dimension(4,nmol):: quatpred
 real(8),dimension(3,3,nmol):: intens
 real(8),dimension(3,nmol):: inprin
 real(8),dimension(nmol):: masstot
 real(8),dimension(3,3,nmol):: inaxis
 integer,dimension(nmol):: lhead
 integer,dimension(nat):: llist
 integer:: nmol
!Rigid stuff
 real(8),dimension(3,nmol):: xcart
 real(8),dimension(3,nmol):: fposcur
 real(8),dimension(3,nmol):: accposcur
 real(8),dimension(3,nmol):: accpospred
 real(8),dimension(3,nmol):: accposprev
 real(8),dimension(3,nmol):: fpospred
 real(8),dimension(3,nmol):: vpospred
 real(8),dimension(3,nmol):: poscur
 real(8),dimension(3,nmol):: vxyz
 real(8),dimension(3,nmol):: vposcur
 real(8),dimension(3,nat):: pospred
 real(8),dimension(3,nat):: dxred
 real(8),dimension(3,3):: dlatvec
 real(8),dimension(3,3):: latvec
 real(8),dimension(3,3):: latvec0
 real(8),dimension(3,3):: latinv
 real(8),dimension(3,3):: unitmat
 real(8),dimension(3,3):: elmatcur
 real(8),dimension(3,3):: acclatcur
 real(8),dimension(3,3):: acclatpred
 real(8),dimension(3,3):: acclatprev
 real(8),dimension(3,3):: flatcur
 real(8),dimension(3,3):: flatpred
 real(8),dimension(3,3):: velmatcur
 real(8),dimension(3,3):: tmplat
 real(8),dimension(3,3):: tmplatt
 real(8),dimension(3,3):: tpred
 real(8),dimension(3,3):: f0
 real(8),dimension(3,3):: f0inv
 real(8),dimension(3,3):: latcur
 real(8),dimension(3,3):: vlatcur
 real(8),dimension(3,3):: vlat
 real(8),dimension(3,3):: lattrans
 real(8),dimension(3,3):: latdottrans
 real(8),dimension(3,3):: a 
 real(8),dimension(3,3):: g
 real(8),dimension(3,3):: gdot
 real(8),dimension(3,3):: ginv
 real(8),dimension(3,3):: gtot
 real(8),dimension(3,3):: sigma
 real(8),dimension(3,3):: sigmatrans
 real(8),dimension(3,3):: latpred
 real(8),dimension(3,3):: tcur
 real(8),dimension(3,3):: vlatpred
 real(8),dimension(3,3):: velmatpred
 real(8),dimension(3,3):: str_matrix
 real(8):: accvolcur
 real(8):: accvolpred
 real(8):: accvolprev
 real(8):: fvolcur
 real(8):: fvolpred
 real(8):: volcur
 real(8):: volpred
 real(8):: vvolcur
 real(8):: vvolpred
 real(8):: volpred_1_3
 real(8):: vol_1_3_inv
 real(8):: vol_1_3

 real(8):: vposcurtmp(3)
 real(8):: vpospredtmp(3)
 real(8):: crossp(3)
 real(8):: velcm(3)
 real(8):: latmass
 real(8):: latmassinv
 real(8):: latmass0
 real(8):: ekinatom
 real(8):: ekinatom_prev
 real(8):: ekinlat
 real(8):: ekinlat_prev
 real(8):: rkin
 real(8):: enthalpy
 real(8):: enmin1
 real(8):: enmin2
 real(8):: ent_pos_0
 real(8):: en0000
 real(8):: e_rxyz
 real(8):: econs_max
 real(8):: econs_min
 real(8):: torquenrm 
 real(8):: ecut_tmp
 real(8):: toldff_tmp
 real(8):: counter
 real(8):: dt
 integer:: i
 integer:: j
 integer:: iat
 integer:: itime
 integer:: nummax
 integer:: nummin
 integer:: iprec 
 integer:: options 
 integer:: md_type
 logical:: getwfk

 character(40)::filename,folder
 character(4) ::fn4
 if(any(fixat).and.nmol.ne.nat) stop "MD with fixed atoms not implemented yet!"
 write(*,*) "# PARAMETERS: Ha_eV,kb_HaK,amu_emass -->",Ha_eV,kb_HaK,amu_emass


!Assign masses to each atom (for MD)
 do iat=1,nmol
   amass(iat)=masstot(iat)
   write(*,'(a,i5,1(1x,es15.7))') " # MD: imol, AMU, EM: ", iat,amass(iat)
 enddo
!******************************************************************
!NEW VERISON OF MD: Directly implemented, simplest Parrinello-Rahman and other MD

!Here we split the routine in the Abinit native part and my new implentation
md_type=parres%md_algo
if(md_type==1) then
 write(*,'(a)') ' # Entering standalone Parrinello Rahman MD '
elseif(md_type==2) then
 write(*,'(a)') ' # Entering standalone Cleveland MD '
elseif(md_type==3) then
 write(*,'(a)') ' # Entering standalone Wentzcovitch MD '
elseif(md_type==4) then
 write(*,'(a)') ' # Entering standalone Andersen MD '
endif

!The "reduced" coordinates in Andersen are quite different from the ones in PR
!Set temporary variables, initially
  vxyz(:,:)=vel_cm_in(:,:)
  pressure=target_pressure_habohr
  unitmat=0.d0
  do i=1,3
    unitmat(i,i)=1.d0
  enddo
  latmass0=parini%bmass*amu_emass !This means we use the barostat mass as the lattice mass (in ELECTRON MASS)
  vlat=vel_lat_in  !The initial cell velocity
  itime=0
  dt=parres%dtion_md

!Set options=1 for Velocity Verlet of cell dynamics
!Set options=2 for Normal Verlet of cell dynamics
!Set options=3 for Beeman integration scheme, corrector-predictor
  options=parres%md_integrator
  if(options.lt.1.or.options.gt.3) stop "Wrong algo option"

!MD-type: 1 for PR and 2 for Cleveland and 3 for Wentzcovitch and 4 for Andersen
  md_type=parres%md_algo
  if(md_type.lt.1.or.md_type.gt.4) stop "Wrong integrator option"

  write(*,'(a,i3,a,i3)') " # MD Algorithm: ",md_type, ", MD Integrator: ",options

!Now we run my implementation of MD
pressure_md=0.d0
pressure_md(1,1)=1.d0;pressure_md(2,2)=1.d0;pressure_md(3,3)=1.d0
pressure_ener=0.d0;pressure_md=pressure_md*pressure  !Here the pressure is not passed to the energyandforces, so we move on the ENERGY surface

!Transform the initial velocities given in cartesian coordinates into velocities of the intenal coordinates
       latvec=latvec_in
       call getvol(latvec,vol)
if(md_type==4) then
       vol_1_3=vol**(1.d0/3.d0)
       vol_1_3_inv=1.d0/vol_1_3
       do iat=1,nmol
          vposcur(:,iat)=vxyz(:,iat)*vol_1_3_inv
       enddo
else
       call invertmat(latvec,latinv,3)
       do iat=1,nmol
          vposcur(:,iat)=matmul(latinv,vxyz(:,iat))
       enddo
endif

!Compute f0
        call cross_product(latvec(:,2),latvec(:,3),crossp); sigma(:,1)=crossp
        call cross_product(latvec(:,3),latvec(:,1),crossp); sigma(:,2)=crossp
        call cross_product(latvec(:,1),latvec(:,2),crossp); sigma(:,3)=crossp
        sigmatrans(:,1)=sigma(1,:)
        sigmatrans(:,2)=sigma(2,:)
        sigmatrans(:,3)=sigma(3,:)
        f0=matmul(sigmatrans,sigma)
        call invertmat(f0,f0inv,3)


!!ATTENTION: THIS IS WITH TAKING INTO ACCOUNT THE STRESS TENSOR GIVEN IN INPUT

!Keep track of volume
        vol0=vol
        latvec0=latvec/vol**(1.d0/3.d0)
        volprev=vol

        if(md_type==1) then
!!This is for the original Rahman Parrinello. You also need to change the volume-mass factor
        latmass=latmass0
        latmassinv=1.d0/latmass
        elseif(md_type==2) then
!This is for cleveland
        latmass=latmass0/vol**(4.d0/3.d0)
        latmassinv=1.d0/latmass
        elseif(md_type==3) then
!This is for Wentzcovitch
        latmass=latmass0
        latmassinv=1.d0/latmass
        elseif(md_type==4) then
!This is for Andersen
        latmass=latmass0
        latmassinv=1.d0/latmass
        endif

!Initialize internal coordinates and lattice vectors
if(md_type==4) then
        call rxyz_int2cart(latvec,xred_cm_in,xcart,nmol)
        poscur=xcart*vol_1_3_inv           !call rxyz_cart2int(latvec,poscur,rxyz,nat)
        volcur=vol
        vvolcur=vvol_in
else
        poscur=xred_cm_in           !call rxyz_cart2int(latvec,poscur,rxyz,nat)
        latcur=latvec
        vlatcur=vlat
endif

!Initiallize quat
quatcur=quat_in
quatpred=quat_in

!MHM: initialize variable to track the minima/maxima***********
    write(*,*) '# MINHOP start MD'
    nummax=0
    nummin=0
    enmin1=0.d0
    en0000=0.d0
    econs_max=-1.d100
    econs_min=1.d100
!**************************************************************

!Compute initial kinetic energies
!Here we need to add rotational energies!!!
if(md_type==4) then
  latcur=latvec0*vol_1_3
  rkin=0.d0
  do iat=1,nmol
     vposcurtmp=vposcur(:,iat)*vol_1_3
     rkin=rkin+amass(iat)*(vposcurtmp(1)**2+vposcurtmp(2)**2+vposcurtmp(3)**2)
  enddo
  ekinatom=0.5d0*rkin
  rkin=vvolcur**2*latmass
  ekinlat=0.5d0*rkin
else
  rkin=0.d0
  do iat=1,nmol
     vposcurtmp=matmul(latcur,vposcur(:,iat))
     rkin=rkin+amass(iat)*(vposcurtmp(1)**2+vposcurtmp(2)**2+vposcurtmp(3)**2)
  enddo
  ekinatom=0.5d0*rkin
  rkin=0.d0
  do i=1,3
     rkin=rkin+(vlatcur(1,i)**2+vlatcur(2,i)**2+vlatcur(3,i)**2)*latmass
  enddo
  ekinlat=0.5d0*rkin
endif


!MHM: Write output to file in every step***********************************
!INITIAL STEP, STILL THE SAME STRUCTURE AS INPUT
       call expand_rigid(latvec_in,xred_cm_in,quat_in,xcart_mol,lhead,llist,nat,nmol,xred_in)
       write(*,*) "Pressure, Energy",pressure,etot_in
       call get_enthalpy(latvec_in,etot_in,pressure,enthalpy)
       ent_pos_0=enthalpy
       en0000=enthalpy-ent_pos_0
       write(*,'(a,i5,1x,4(1x,1pe17.10),3(1x,i2))') ' # MD: it,enthalpy,pV,ekinat,ekinlat,nmax,nmin,mdmin ',&
             &itime,enthalpy,pressure*vol,ekinatom,ekinlat,nummax,nummin,parres%mdmin
       write(fn4,'(i4.4)') itime
       filename=trim(folder)//"posmd."//fn4//".ascii"
       units=units
       write(*,*) "# Writing the positions in MD:",filename
       call write_atomic_file_ascii(parres,filename,nat,units,xred_in,latvec_in,fcart_in,strten_in,&
       &char_type(1:ntypat),ntypat,typat,fixat,fixlat,etot_in,pressure,enthalpy,en0000)
!*********************************************************************


call get_fcm_torque(fcart_cm,torque,fcart_in,quat_in,xcart_mol,lhead,llist,nat,nmol)
call acceleration(pressure_md,accposcur,acclatcur,accvolcur,vposcur,vlatcur,vvolcur,&
                  strten_in,fcart_cm,latcur,amass,latmass,f0inv,md_type,nmol)
        accposprev=accposcur
        acclatprev=acclatcur
        accvolprev=accvolcur

!Eliminate the unnecessary atmic and cell components
call elim_fixed_at(nmol,accposcur);call elim_fixed_at(nmol,accposprev)
if(md_type.ne.4) call elim_fixed_lat(latcur,acclatcur)
if(md_type.ne.4) call elim_fixed_lat(latcur,acclatprev)
call elim_fixed_at(nmol,vposcur)
if(md_type.ne.4) call elim_fixed_lat(latcur,vlatcur)


        do itime=1,parres%nmd_dynamics

!          if(itime.ne.1) e_rxyz=enthalpy  !e_rxyz=e_rxyz+pressure_md(1,1)*vol   
!Check the torque on the cell for rotation
          call torque_cell(latcur,vlatcur,torquenrm)
!Now perform the step
          velcm=0.d0
          rkin=0.d0
          if(options==1) then
!For velocity verlet
!             latpred(:,:)=latcur(:,:) + dt*vlatcur(:,:) + (.5d0*dt*dt*acclatcur)
             dlatvec(:,:)=dt*vlatcur(:,:) + (.5d0*dt*dt*acclatcur)
          do iat=1,nmol
              dxred(:,iat)=dt*vposcur(:,iat) + 0.5d0*dt*dt*accposcur(:,iat)  !0.5 was missing before
!             pospred(:,iat)=poscur(:,iat) + dt*vposcur(:,iat) + dt*dt*accposcur(:,iat)
          enddo
          call propagate(nmol,poscur,latcur,dxred,dlatvec,pospred,latpred)
          do iat=1,nmol
             vpospred(:,iat)=(pospred(:,iat)-poscur(:,iat))/dt   !Update dummy variable vposcur for next iteration
          enddo
            if(md_type==4) then
               volpred=volcur+dt*vvolcur+0.5d0*dt*dt*accvolcur
            endif
          elseif(options==2) then
!Instead if normal verlet is used:
              dlatvec(:,:)=dt*vlatcur + dt*dt*acclatcur
          do iat=1,nmol
              dxred(:,iat)=dt*vposcur(:,iat) + dt*dt*accposcur(:,iat)
          enddo
          call propagate(nmol,poscur,latcur,dxred,dlatvec,pospred,latpred)
          vlatpred(:,:)=(latpred-latcur)/dt
          do iat=1,nmol
              vpospred(:,iat)=(pospred(:,iat)-poscur(:,iat))/dt   !Update dummy variable vposcur for next iteration
          enddo
            if(md_type==4) then
               volpred=volcur+dt*vvolcur+dt*dt*accvolcur
               vvolpred=(volpred-volcur)/dt
            endif
          elseif(options==3) then
!Predictor part of Beeman for cell
              dlatvec(:,:)=dt*vlatcur + dt*dt/6.d0*(4.d0*acclatcur(:,:)-acclatprev(:,:))
          do iat=1,nmol
!Predictor part of Beeman
             dxred(:,iat)=dt*vposcur(:,iat) + dt*dt/6.d0*(4.d0*accposcur(:,iat)-accposprev(:,iat))             
          enddo
          call propagate(nmol,poscur,latcur,dxred,dlatvec,pospred,latpred)
          vlatpred(:,:)=vlatcur(:,:) + 0.5d0*(3.d0*acclatcur(:,:)-acclatprev(:,:))  
          do iat=1,nmol
             vpospred(:,iat)=vposcur(:,iat)+0.5d0*(3.d0*accposcur(:,iat)-accposprev(:,iat))
          enddo
            if(md_type==4) then
               volpred=volcur+dt*vvolcur+dt*dt/6.d0*(4.d0*accvolcur-accvolprev)
               vvolpred=volcur+0.5d0*(3.d0*accvolcur-accvolprev)
            endif
          endif


!Rotational part to predict the correct quatpred
!Implement the case for assymetric top
!call rbmd_driver






if(md_type==4) then
  do iat=1,nmol
     vposcurtmp=vposcur(:,iat)*vol_1_3
     rkin=rkin+amass(iat)*(vposcurtmp(1)**2+vposcurtmp(2)**2+vposcurtmp(3)**2)
     velcm=velcm+vposcurtmp*amass(iat)
  enddo
else
  do iat=1,nmol
             vposcurtmp=matmul(latcur,vposcur(:,iat))
             rkin=rkin+amass(iat)*(vposcurtmp(1)**2+vposcurtmp(2)**2+vposcurtmp(3)**2)
             velcm=velcm+vposcurtmp*amass(iat)
  enddo
endif

!Kinetic energy of atoms
          ekinatom=0.5d0*rkin
          if(parres%verb.gt.0) write(*,'(a,es15.7)')"Velocity of CM: ", sqrt(velcm(1)**2+velcm(2)**2+velcm(3)**2)
!Kinetic energy according to Parrinello Rahman
          rkin=0.d0
if(md_type==4) then
  rkin=vvolcur**2*latmass
  ekinlat=0.5d0*rkin
else
          do i=1,3
             rkin=rkin+(vlatcur(1,i)**2+vlatcur(2,i)**2+vlatcur(3,i)**2)*latmass
          enddo
          ekinlat=0.5d0*rkin
endif
          rkin=ekinlat+ekinatom
if(parres%verb.gt.0) write(*,'(a,3(1x,es15.7))') " # Torquenrm, Ekin, Enthalpy: ",torquenrm, rkin,enthalpy
!Update counter
          enmin2=enmin1
          enmin1=en0000

!Here we perform the force call
!****************************************************************************************************************        
!Get the current acell_in and rprim_in, and also the atomic positions
if(md_type==4) then
       volpred_1_3=volpred**(1.d0/3.d0)
       latvec_in=latvec0*volpred_1_3
       latpred=latvec_in
       xcart=pospred*volpred_1_3
       call rxyz_cart2int(latvec_in,xred_cm_in,xcart,nmol)
       call expand_rigid(latvec_in,xred_cm_in,quatpred,xcart_mol,lhead,llist,nat,nmol,xred_in)
else
       latvec_in=latpred
       xred_cm_in=pospred
       call expand_rigid(latvec_in,pospred,quatpred,xcart_mol,lhead,llist,nat,nmol,xred_in)
endif
       if(itime==1)  then
           getwfk=.false.
       elseif(itime.ne.1.and.usewf_md) then
           getwfk=.true.
       else
           getwfk=.false.
       endif
       call get_energyandforces_single(parres,latvec_in,xred_cm_in,fcart_in,strten_in,etot_in,iprec,getwfk)
       call get_fcm_torque(fcart_cm,torque,fcart_in,quatpred,xcart_mol,lhead,llist,nat,nmol)
!Single point calculation finished. Now returning to MD part. All output variables are now in *_in
!****************************************************************************************************************        

!For velocity verlet of cell
        if(options==1) then
           call acceleration(pressure_md,accpospred,acclatpred,accvolpred,vpospred,&
                &vlatcur,vvolcur,strten_in,fcart_cm,latpred,amass,latmass,f0inv,md_type,nmol) 
           vlatpred(:,:)=vlatcur(:,:)+0.5d0*dt*(acclatpred+acclatcur)
           if(md_type==4) vvolpred=vvolcur+0.5d0*dt*(accvolpred+accvolcur)
        elseif(options==3) then
!Corrector part of Beeman. Note that a fixed number of iterations are used!!!
           call ekin_at_lat_andersen(amass,latmass,latpred,vpospred,vlatpred,vvolpred,ekinatom_prev,ekinlat_prev,f0,md_type,nmol)
           do i=1,5
             call acceleration(pressure_md,accpospred,acclatpred,accvolpred,vpospred,&
                  &vlatpred,vvolpred,strten_in,fcart_cm,latpred,amass,latmass,f0inv,md_type,nmol) 
             do iat=1,nmol
             vpospred(:,iat)=vposcur(:,iat)+dt/6.d0*(2.d0*accpospred(:,iat)+5.d0*accposcur(:,iat)-accposprev(:,iat))
             enddo
             vlatpred=vlatcur+dt/6.d0*(2.d0*acclatpred+5.d0*acclatcur-acclatprev)
             if(md_type==4) vvolpred=vvolcur+dt/6.d0*(2.d0*accvolpred+5.d0*accvolcur-accvolprev)
             call ekin_at_lat_andersen(amass,latmass,latpred,vpospred,vlatpred,vvolpred,ekinatom,ekinlat,f0,md_type,nmol)
if(parres%verb.gt.1)write(* ,'(a,i5,1x,es15.5,es15.5)') " # Beeman: Corrector relative convergence: it, atoms, lattice: ",&
                  & i,abs(ekinatom-ekinatom_prev)/ekinatom,abs(ekinlat-ekinlat_prev)/max(ekinlat,1.d-100)
!             write(67,'(a,i5,1x,es15.5,es15.5)') " # Beeman: Corrector relative convergence: it, atoms, lattice: ",&
!                  & i,abs(ekinatom-ekinatom_prev)/ekinatom,abs(ekinlat-ekinlat_prev)/max(ekinlat,1.d-100)
             ekinlat_prev=ekinlat
             ekinatom_prev=ekinatom
           enddo 
        endif  
call acceleration(pressure_md,accpospred,acclatpred,accvolpred,vpospred,&
     &vlatpred,vvolpred,strten_in,fcart_cm,latpred,amass,latmass,f0inv,md_type,nmol)

!Eliminate the unnecessary atmic and cell components
call elim_fixed_at(nmol,accpospred);call elim_fixed_at(nmol,vpospred)
if(md_type.ne.4) call elim_fixed_lat(latcur,acclatpred)
if(md_type.ne.4) call elim_fixed_lat(latcur,vlatpred)

!Compute the "predicted" kinetic energies:
call ekin_at_lat_andersen(amass,latmass,latpred,vpospred,vlatpred,vvolpred,ekinatom,ekinlat,f0,md_type,nmol)


!MHM: Write output to file in every step***********************************
       counter=real(itime,8)
       call get_enthalpy(latvec_in,etot_in,pressure,enthalpy)
       call getvol(latvec_in,vol)
       volpred=vol
       en0000=enthalpy-ent_pos_0
       if (itime >= 3 .and. enmin1 > enmin2 .and. enmin1 > en0000)  nummax=nummax+1
       if (itime >= 3 .and. enmin1 < enmin2 .and. enmin1 < en0000)  nummin=nummin+1
!       write(67,'(a,i5,1x,4(1x,1pe17.10),3(1x,i2))') ' MD: it,enthalpy,pV,ekinat,ekinlat,nmax,nmin,mdmin ',&
!             &itime,enthalpy,pressure*vol,ekinatom,ekinlat,nummax,nummin,parres%mdmin
       write(*,'(a,i5,1x,4(1x,1pe17.10),3(1x,i2))') ' # MD: it,enthalpy,pV,ekinat,ekinlat,nmax,nmin,mdmin ',&
             &itime,enthalpy,pressure*vol,ekinatom,ekinlat,nummax,nummin,parres%mdmin
       write(fn4,'(i4.4)') itime
       filename=trim(folder)//"posmd."//fn4//".ascii"
       units=units
       write(*,*) "# Writing the positions in MD: ",filename
       call write_atomic_file_ascii(parres,filename,nat,units,xred_in,latvec_in,fcart_in,strten_in,&
       &char_type(1:ntypat),ntypat,typat,fixat,fixlat,etot_in,pressure,enthalpy,en0000)
!*********************************************************************
        econs_max=max(econs_max,rkin+enthalpy)
        econs_min=min(econs_min,rkin+enthalpy)
       if (nummin.ge.parres%mdmin) then
          if (nummax.ne.nummin) &
               write(*,*) '# WARNING: nummin,nummax',nummin,nummax

!          write(67,*) " MD finished: exiting!"
          write(*,*) " MD finished: exiting!"
          exit
       endif
!Update the variables for next iteration
        accposprev=accposcur
        accposcur=accpospred
        acclatprev=acclatcur
        acclatcur=acclatpred
        velmatcur=velmatpred
        flatcur=flatpred
        fposcur=fpospred
        poscur=pospred
        vlatcur=vlatpred
        vposcur=vpospred
        latcur=latpred

        accvolprev=accvolcur
        accvolcur=accvolpred
        volcur=volpred
        vvolcur=vvolpred

     enddo 
!Adjust MD stepsize
!Minimum number of steps per crossed minimum is 15, average should be parres%nit_per_min
     if(parres%auto_dtion_md) then
       if(real(itime,8)/real(parres%mdmin,8).lt.real(parres%nit_per_min,8)) then
         parres%dtion_md=parres%dtion_md*1.d0/1.1d0
       else
         parres%dtion_md=parres%dtion_md*1.1d0 
       endif
     parres%dtion_md=min(parres%dtion_md,real(itime,8)*dt/(real(parres%mdmin,8)*15.d0))
     write(*,'(3(a,es10.2))') " # MD: steps per minium: ",real(itime,8)/real(parres%mdmin,8),&
           &", parres%dtion_md set to: ",parres%dtion_md,", upper boundary: ",real(itime,8)*dt/(real(parres%mdmin,8)*15.d0) 
     endif
!MD stopped, now do relaxation
     write(*,'(a,i5,es15.7,es15.7)') ' # EXIT MD ',itime,enthalpy,etot_in

end subroutine

!**********************************************************************************************

subroutine init_fp(fp_len,latvec)
!This routine will initiallize the parameters for the fingerprinting
!For 10<fp_method<20: fully periodic systems
!For 20<fp_method<30: molecular systems
use mod_interface
use fingerprint
use global, only: ntypat,nat,typat,units
use defs_basis, only: Bohr_Ang,pi
implicit none
integer:: fp_len,iat,nmax
real(8):: convert,latvec(3,3),vol
!Convert cutoff, sigma and dbin into atomic units
convert=1.d0
if(units=="angstroem") then
  convert=1.d0/Bohr_Ang
endif

!Get recomended size of the continuous oganov FP
if(fp_method==15.or.fp_method==16) then
   call getvol(latvec,vol)
   call estimate_nmax_per_atom(vol,nat,ntypat,fp_rcut*convert,pi,nmax)
   write(*,*) vol,fp_rcut*convert,nmax,nat,pi,ntypat
   write(*,'(a,i10)') " # Estimated fingerprint size for COGANOV and CAOGANOV: ",nmax
   if(nmax.gt.fp_at_nmax) write(*,'(a,i10,i10)') " # WARNING: FPATNMAX too small!", fp_at_nmax, nmax
endif

select case(fp_method)
  case(11)!Oganov method
!     read(56,*) fp_11_rcut,fp_11_sigma,fp_11_dbin
     fp_11_rcut=fp_rcut*convert
     fp_11_sigma=fp_sigma*convert
     fp_11_dbin=fp_dbin*convert
     fp_11_fp_size=ceiling(fp_11_rcut/fp_11_dbin)
     fp_11_fp_dim=ntypat*(ntypat+1)/2
     fp_len=fp_11_fp_size*fp_11_fp_dim
     if(.not.allocated(fp_11_nkinds_sum)) allocate(fp_11_nkinds_sum(ntypat))
     fp_11_nkinds_sum=0
     do iat=1,nat
       fp_11_nkinds_sum(typat(iat))=fp_11_nkinds_sum(typat(iat))+1
     enddo
  case(12)!Calypso method
!TEMPORARY READ FROM STDINPUT
!     read(*,*) tmpvar,fp_12_nl
!     read(56,*) tmpvar,fp_12_nl
     fp_12_nl=fp_nl
     fp_12_fp_dim=ntypat*(ntypat+1)/2
     fp_len=fp_12_fp_dim*fp_12_nl
     if(.not.allocated(fp_12_r_cut)) allocate(fp_12_r_cut(fp_12_fp_dim))
     fp_12_r_cut=fp_rcut*convert
  case(13)!Modified Calypso method
!TEMPORARY READ FROM STDINPUT
!     read(*,*) tmpvar,fp_13_nl
!     read(56,*) tmpvar,fp_13_nl
     fp_12_nl=fp_nl
     fp_13_fp_dim=ntypat!*(ntypat+1)/2
     fp_len=fp_13_nl*nat*fp_13_fp_dim
     if(.not.allocated(fp_13_r_cut)) allocate(fp_13_r_cut(fp_13_fp_dim))
     fp_13_r_cut=fp_rcut*convert
  case(14)!xyz2sm
     fp_len=3*fp_14_m*nat
  case(15)!Continuous Oganov method
     fp_15_rcut=fp_rcut*convert
     fp_15_sigma=fp_sigma*convert
     fp_15_fp_size=fp_at_nmax
     fp_15_fp_dim=ntypat*(ntypat+1)/2
     fp_len=3*fp_15_fp_size*fp_15_fp_dim
!Careful: the FP has 3 entries for the gaussians
     if(.not.allocated(fp_15_nkinds_sum)) allocate(fp_15_nkinds_sum(ntypat))
     fp_15_nkinds_sum=0
     do iat=1,nat
       fp_15_nkinds_sum(typat(iat))=fp_15_nkinds_sum(typat(iat))+1
     enddo
  case(16)!Continuous Atomic Oganov method
     fp_16_rcut=fp_rcut*convert
     fp_16_sigma=fp_sigma*convert
     fp_16_fp_size=fp_at_nmax
     fp_16_fp_dim=ntypat*(ntypat+1)/2
     fp_len=3*fp_16_fp_size*ntypat*nat
!Careful: the FP has 3 entries for the gaussians, nmax entries for all possible neighbors, ndim possible AB interaction, nat atomic lists of gaussians
     if(.not.allocated(fp_16_nkinds_sum)) allocate(fp_16_nkinds_sum(ntypat))
     fp_16_nkinds_sum=0
     do iat=1,nat
       fp_16_nkinds_sum(typat(iat))=fp_16_nkinds_sum(typat(iat))+1
     enddo
  case(17)
     fp_len=fp_17_lseg*(ntypat+1)*nat
  case(18) !Molecular gaussian orbital fingerprint
     fp_len=fp_18_lseg*fp_18_molecules_sphere*fp_18_principleev*fp_18_molecules
  case(21)!Gaussian molecular overlap
!The method only has a FP of length nat
     fp_len=nat  !If we only have stype orbitals, alse fp_len=4*nat
  case default
     stop "Wrong choice for FP"
end select
close(56)
end subroutine

!**********************************************************************************************

subroutine get_fp(fp_len,pos_red,latvec,fp)
!This routine will initiallize the parameters for the fingerprinting
!For 10<fp_method<20: fully periodic systems
!For 20<fp_method<30: molecular systems
use mod_interface
use fingerprint, only: fp_15_fp_size, fp_method, fp_11_rcut, fp_11_sigma, fp_11_dbin
use fingerprint, only: fp_11_fp_size, fp_11_nkinds_sum, fp_11_fp_dim, fp_17_natx_sphere
use fingerprint, only: fp_12_r_cut, fp_12_fp_dim, fp_16_fp_size, fp_12_nl, fp_13_nl
use fingerprint, only: fp_13_r_cut, fp_16_fp_dim
use fingerprint
use global, only: ntypat,nat,typat,rcov,char_type
implicit none
integer:: fp_len,iat,natmol
real(8):: fp(fp_len),pos_red(3,nat),latvec(3,3),rxyz(3,nat),vol,rcov_arr(nat),fp_coganov_atomic(3,fp_15_fp_size,ntypat,nat)
real(8):: rvan(nat) !nat*molecules)
character(len=2):: finalchar(nat) ! dimension(nat*molecules)

select case(fp_method)
  case(11)!Oganov method
     call getvol(latvec,vol)
     write(*,'(a,es15.7)') " # Suggested minimal cutoff radius for Oganov FP: ", vol**(1.d0/3.d0)*2.d0
     call rxyz_int2cart(latvec,pos_red,rxyz,nat)
     call get_fp_oganov(nat,rxyz,latvec,fp_11_rcut,fp_11_sigma,fp_11_dbin,&
          &typat,ntypat,fp_11_nkinds_sum,fp_11_fp_size,fp_11_fp_dim,fp)
  case(12)!Calypso method
     call rxyz_int2cart(latvec,pos_red,rxyz,nat)
     call get_fp_calypso(nat,rxyz,latvec,fp_12_r_cut,typat,ntypat,fp_12_fp_dim,fp_12_nl,fp)
  case(13)!Modified Calypso method
     call rxyz_int2cart(latvec,pos_red,rxyz,nat)
     call get_fp_malypso(nat,rxyz,rcov,latvec,fp_13_r_cut,typat,ntypat,fp_13_fp_dim,fp_13_nl,fp)
  case(14)!xyz2sm fingerprint
     do iat=1,nat
        rcov_arr(iat)=rcov(typat(iat))
     enddo
     call rxyz_int2cart(latvec,pos_red,rxyz,nat)
     call xyz2sm(nat,latvec,rxyz,rcov_arr,fp_14_w1,fp_14_w2,fp_14_m,fp)
  case(15)!C-Oganov fingerprint
     call rxyz_int2cart(latvec,pos_red,rxyz,nat)
     call get_fp_coganov(nat,rxyz,latvec,fp_15_rcut,fp_15_sigma,rcov,&
          &typat,ntypat,fp_15_nkinds_sum,fp_15_fp_size,fp_15_fp_dim,fp)
  case(16)!C-Atomic-Oganov fingerprint
     call rxyz_int2cart(latvec,pos_red,rxyz,nat)
     call get_fp_coganov_atomic(nat,rxyz,latvec,fp_16_rcut,fp_16_sigma,rcov,&
          &typat,ntypat,fp_16_nkinds_sum,fp_16_fp_size,fp_16_fp_dim,fp)
  case(17)!GOM
     do iat = 1, nat
        rcov_arr(iat) = rcov(typat(iat))
     end do
     call rxyz_int2cart(latvec,pos_red,rxyz,nat)
     call get_fp_gauss(nat, ntypat, fp_17_natx_sphere, typat, fp_17_lseg, fp_17_width_cutoff,&
          & fp_17_nex_cutoff, latvec, rxyz, rcov_arr, fp)
  case(18)!MOLGOM
!This fingerprint wants to have the number of atoms per molecule
     natmol=nat/fp_18_molecules     
     if(natmol*fp_18_molecules.ne.nat) stop "Something wrong with the number of molecules"
     call findmolecule(rxyz,latvec,finalchar,pos_red,char_type,typat,ntypat,natmol)
!
! Assign the Van-der-Waals radii
!
     do iat = 1, nat
         call sym2rvan(finalchar(iat), rvan(iat))
     end do
!
!from this system a fingerprint is taken
!
     call periodic_fingerprint(rxyz,latvec,finalchar,rvan,fp,natmol)

  case(21)!Gaussian molecular overlap
     do iat=1,nat
        rcov_arr(iat)=rcov(typat(iat))
     enddo
     call rxyz_int2cart(latvec,pos_red,rxyz,nat)
     call fingerprint_gaussmol(nat,fp_len,rxyz,rcov_arr,fp)
  case default
     stop "Wrong choice for FP"
end select
end subroutine

!**********************************************************************************************

subroutine get_fp_distance(parini,fp_len,fp1,fp2,fp_dist)
use mod_parini, only: typ_parini
!This routine will initiallize the parameters for the fingerprinting
!For 10<fp_method<20: fully periodic systems
!For 20<fp_method<30: molecular systems
use mod_interface
use fingerprint
use global, only: ntypat,nat,typat
use defs_basis, only: pi
implicit none
type(typ_parini), intent(in):: parini
integer:: fp_len
real(8):: fp(fp_len),pos_red(3,nat),latvec(3,3),rxyz(3,nat),fp1(fp_len),fp2(fp_len),fp_dist
select case(fp_method)
  case(11)!Oganov method
        call get_cosinedistance(parini,fp1,fp2,fp_11_fp_size,fp_11_fp_dim,ntypat,fp_11_nkinds_sum,fp_dist)
  case(12)!Calypso method
        call get_distance_calypso(fp1,fp2,fp_12_fp_dim,fp_12_nl,fp_dist)
  case(13)!Modified Calypso method
        call get_distance_malypso(fp1,fp2,fp_13_fp_dim,nat,typat,fp_13_nl,fp_dist)
  case(14)!xyz2sm
        call get_distance_xyz2sm(nat,typat,fp_len/nat,fp1,fp2,fp_dist)
  case(15)!Continuous Oganov
        call get_cosinedistance_coganov(fp1,fp2,fp_15_fp_size,fp_15_fp_dim,ntypat,fp_15_nkinds_sum,fp_15_rcut,pi,fp_dist)
  case(16)!Continuous Atomic Oganov
        call get_cosinedistance_coganov_atomic(fp1,fp2,nat,fp_16_fp_size,fp_16_fp_dim,typat,ntypat,fp_16_nkinds_sum,fp_16_rcut,pi,fp_dist)
  case(17)!GOM
        call get_distance_gauss(fp1, fp2, fp_17_lseg, nat, ntypat, typat, fp_dist)
  case(18)!MOLGOM
        call get_distance_molgom(fp1,fp2,fp_dist,fp_18_lseg,fp_18_molecules,fp_18_molecules_sphere,fp_18_principleev)
  case(21)!Gaussian molecular overlap
        call fpdistance_gaussmol(fp_len,fp1,fp2,fp_dist)
  case default
     stop "Wrong choice for FP"
end select
end subroutine

!**********************************************************************************************

subroutine identical(parini,nlminx,nlmin,fp_method,fp_len,ent_wpos,fp_wpos,ent_arr,fp_arr,&
           &ent_delta,fp_delta,newmin,kid,fp_dist_min,k_e_wpos,n_unique,n_nonuni,lid,nid)
use mod_parini, only: typ_parini
use mod_interface
implicit none
type(typ_parini), intent(in):: parini
integer:: nlminx,nlmin,fp_len,kid,k_e_wpos,n_unique,n_nonuni
integer:: i,l,klow,k,khigh,fp_method,lid(nlminx),nid
real(8):: fp_arr(fp_len,nlminx),fp_wpos(fp_len),ent_arr(nlminx),ent_wpos,fp_delta,ent_delta,fp_dist_min,fp_dist
logical newmin,inrange

inrange=.false.

!C  check whether new minimum
call hunt(ent_arr,min(nlmin,nlminx),ent_wpos,k_e_wpos)
newmin=.true.
!!do i=1,nlmin
!!write(*,'(a,i3,5(e24.17))') ' # check identical: enarr ',i,ent_arr(i),(fp_arr(l,i),l=1,2)
!!enddo
write(*,'(a,e24.17,i3,5(e24.17))') ' # ID: Check identical: ent_wpos,k_e_wpos ',ent_wpos,k_e_wpos!,(wfp(l),l=1,2)

! find lowest configuration that might be identical
klow=k_e_wpos
do k=k_e_wpos,1,-1
if (ent_wpos-ent_arr(k).lt.0.d0) stop 'zeroA'
if (ent_wpos-ent_arr(k).lt.ent_delta) inrange=.true.
if (ent_wpos-ent_arr(k).gt.ent_delta) exit
klow=k
enddo

! find highest  configuration that might be identical
khigh=k_e_wpos+1
do k=k_e_wpos+1,nlmin
if (ent_arr(k)-ent_wpos.lt.0.d0) stop 'zeroB'
if (ent_arr(k)-ent_wpos.lt.ent_delta) inrange=.true.
if (ent_arr(k)-ent_wpos.gt.ent_delta) exit
khigh=k
enddo

nid=0
if(.not.inrange) then
  write(*,'(a)') ' # ID: there is no minimum in the given enthalpy range'
  write(*,'(a,L4)') ' # ID: wpos is a new minimum: ',newmin
  return  
endif
write(*,'(a,i5,i5)') ' # ID: Check k bounds: ',max(1,klow),min(nlmin,khigh)
fp_dist_min=1.d100
do k=max(1,klow),min(nlmin,khigh)
call get_fp_distance(parini,fp_len,fp_wpos,fp_arr(:,k),fp_dist)
write(*,'(a,i5,es25.15)') ' # ID: Check fp_dist: ',k,fp_dist
!!write(*,'(a,20(e10.3))') '(MH) fp_wpos', (fp_wpos(i),i=1,fp_len)
!!write(*,'(a,20(e10.3))') '(MH) fp_arr ', (fp_arr(i,k),i=1,fp_len)
if (fp_dist.lt.fp_delta) then
    write(*,'(a,i5)') ' # ID: wpos is identical to ',k
    newmin=.false.
    nid=nid+1
    lid(nid)=k
    if (fp_dist.lt.fp_dist_min) then
       fp_dist_min=fp_dist
       kid=k
       k_e_wpos=kid
    endif
endif
enddo

   write(*,'(a,L4)') ' # ID: wpos is a new minimum: ',newmin
   if (nid.gt.1) write(*,*) '# WARNING: more than one identical configuration found'
!          if (nsm.gt.1) write(100+iproc,*) 'WARNING: more than one identical configuration found'
if (nid.eq.1) n_unique=n_unique+1
if (nid.gt.1) n_nonuni=n_nonuni+1

return
end
                     
!**********************************************************************************************

subroutine replace(nlminx,nlmin,fp_len,nat,kid,e_wpos,ent_wpos,fp_wpos,wpos_red,&
  &wpos_latvec,spg_wpos,spgtol_wpos,fdos_wpos,&
  &e_arr,ent_arr,fp_arr,pl_arr,lat_arr,spg_arr,spgtol_arr,dos_arr,ct_arr,findsym)
!Replace the structure kid with wpos, only if the symmetry index is higher in wpos
  use mod_interface
  implicit none
  integer:: fp_len,ct_arr(nlminx),spg_arr(nlminx),nat,iat,spg_wpos
  integer:: k, nlmin,nlminx,i,kid
  real(8):: e_wpos, ent_wpos, wpos_red(3,nat),wpos_latvec(3,3),spgtol_wpos,fdos_wpos,fp_wpos(fp_len)
  real(8):: e_arr(nlminx),ent_arr(nlminx),fp_arr(fp_len,nlminx),pl_arr(3,nat,nlminx)
  real(8):: lat_arr(3,3,nlminx),spgtol_arr(nlminx),dos_arr(nlminx)
  logical:: findsym
  write(*,*) "KKKID",kid,nlmin,nlminx,fp_len,nat
  write(*,*) spg_wpos,spg_arr(kid),spg_wpos,spg_arr(kid),ent_wpos,ent_arr(kid)
  if((findsym.and.((spg_wpos.gt.spg_arr(kid)).or.(spg_wpos.eq.spg_arr(kid).and.ent_wpos.lt.ent_arr(kid)))).or.&
   &(.not.findsym.and.ent_wpos.lt.ent_arr(kid))) then
   write(*,'(a,i4,L3,i4,i4,es15.7,es15.7)') " # Replace array elements: ID,Findsym,SPG_ARR,SPG_NEW,ENT_ARR,ENT_NEW  ",&
   &kid,findsym,spg_arr(kid),spg_wpos,ent_arr(kid),ent_wpos
   e_arr(kid)=e_wpos   
   ent_arr(kid)=ent_wpos
   fp_arr(:,kid)=fp_wpos(:)
   pl_arr(:,:,kid)=wpos_red(:,:)   
   lat_arr(:,:,kid)=wpos_latvec(:,:)
   spgtol_arr(kid)=spgtol_wpos
   dos_arr(kid)=fdos_wpos   
  endif
end subroutine

!**********************************************************************************************

 subroutine dist2plane(point,nvec,ppoint,dist)
 !This subroutine will calculate  the distance between a plane and a point in space
 !The point is 'point', the normalized normal vector of the plane is 'nvec', 'ppoint' is an arbitrary point on the plane
 !and the output is the distance 'dist'  
 use mod_interface
 implicit none
 real(8), intent(in) :: point(3),nvec(3),ppoint(3)
 real(8), intent(out):: dist
 integer             :: i
 real(8)             :: p,nvectmp(3)
 nvectmp(:)=nvec(:)!/sqrt(nvec(1)*nvec(1)+nvec(2)*nvec(2)+nvec(3)*nvec(3))
 p=DOT_PRODUCT(nvectmp,ppoint) 
 p=-p
 dist=DOT_PRODUCT(nvectmp,point)+p
 end subroutine

!**********************************************************************************************

 subroutine dist2line(point,ppoint1,ppoint2,dist)
 !This subroutine will calculate a the distance between a plane and a line in space
 !The point is 'point', 'ppoint1' and 'ppoint2' are arbitrary points on the line
 !and the output is the distance 'dist'  
 use mod_interface
 implicit none
 real(8), intent(in) :: point(3),ppoint1(3),ppoint2(3)
 real(8), intent(out):: dist
 integer             :: i
 real(8)             :: v1(3),v2(3),vnrm,crossp(3)
 v1(:)=ppoint2(:)-ppoint1(:)
 vnrm=sqrt(v1(1)*v1(1)+v1(2)*v1(2)+v1(3)*v1(3))
 v1(:)=v1(:)/vnrm
 v2(:)=point(:)-ppoint1(:)
 call cross_product(v1,v2,crossp)
 dist=sqrt(crossp(1)*crossp(1)+crossp(2)*crossp(2)+crossp(3)*crossp(3))
 end subroutine

!**********************************************************************************************

subroutine compare_lammps(parini)
use mod_interface
use global
use interface_code
use defs_basis
use mod_parini, only: typ_parini
! Main program to test potential subroutines
implicit none
type(typ_parini), intent(in):: parini
integer:: k,l,n,iat,iprec
real(8):: latvec(3,3),xred(3,nat),xcart(3,nat),f_lammps(3,nat),f(3,nat),e_lammps,e,tmp_r,tmp_i,tilts(6),latvec_in(3,3),strten(6),latvec_box(3,3)
character(400):: line_log,line_dump
logical:: getwfk
open(unit=6,file="lammps.log")
open(unit=7,file="dump.atom")
open(unit=8,file="tilts")
read(8,*) tilts
latvec=0.d0
latvec(1,1)=tilts(1)
latvec(2,2)=tilts(2)
latvec(3,3)=tilts(3)
latvec(1,2)=tilts(4)
latvec(1,3)=tilts(5)
latvec(2,3)=tilts(6)
write(*,*) "LATVEC FROM TILTS"
write(*,*) latvec

latvec_box=0.d0
latvec_box(1,1)=tilts(1)+MAX(0.0d0,tilts(4),tilts(5),tilts(4)+tilts(5))-MIN(0.0d0,tilts(4),tilts(5),tilts(4)+tilts(5))
latvec_box(2,2)=tilts(2)+MAX(0.0d0,tilts(6))-MIN(0.0d0,tilts(6))
latvec_box(3,3)=tilts(3)


!!xlo_bound = xlo + MIN(0.0,tilts(4),tilts(5),tilts(4)+tilts(5))
!!xhi_bound = xhi + MAX(0.0,tilts(4),tilts(5),tilts(4)+tilts(5))
!!ylo_bound = ylo + MIN(0.0,tilts(6))
!!yhi_bound = yhi + MAX(0.0,tilts(6))
!!zlo_bound = zlo
!!zhi_bound = zhi 


do while(.true.)
   read(6,'(a400)',end=99) line_log
   n = len_trim(line_log)
   k = index(line_log(1:n),"Step PotEng Press Fmax Fnorm")
   if(k.ne.0) then
     do while(.true.)
       read(6,'(a400)',end=99) line_log
       n = len_trim(line_log)
       read(line_log,*) tmp_i,e_lammps
       write(*,*) tmp_i,e_lammps
       do while(.true.)
           read(7,'(a400)',end=99) line_dump
           n = len_trim(line_dump)
           l = index(line_dump(1:n),"ITEM: TIMESTEP")
           if(l.ne.0) then
              read(7,*) iat
              if(iat.ne.tmp_i) stop "Wrong index"
              cycle
           endif

           l = index(line_dump(1:n),"ITEM: ATOMS")
           if(l.ne.0) then
             do iat=1,nat
               read(7,*) tmp_i,tmp_i,xcart(:,iat),f_lammps(:,iat) 
!               read(7,*) tmp_i,tmp_i,xred(:,iat),f_lammps(:,iat) 
             enddo
             exit
           endif
       enddo
!!       do iat=1,nat
!!          xcart(:,iat)=matmul(latvec_box,xred(:,iat))
!!       enddo
       call rxyz_cart2int(latvec,xred,xcart,nat)
       latvec_in=latvec/Bohr_ang
       call get_energyandforces_single(parini,latvec_in,xred,f,strten,e,iprec,getwfk)
       f=f*Ha_eV/Bohr_Ang
       e=e*Ha_eV
       write(*,*) tmp_i,e_lammps
       if(any(abs(f-f_lammps).gt.1.d-10)) then
          write(*,*) "ENERGIES"
          write(*,'(2es25.15)') e,e_lammps
          write(*,*) "ERROR in forces"
          do iat=1,nat
            write(*,'(6es25.15)') f(:,iat),f_lammps(:,iat)
          enddo
          stop
       endif
       if(abs(e-e_lammps).gt.1.d-10) then
          write(*,*) "ERROR in energies"
          write(*,'(2es25.15)') e,e_lammps
          stop
       endif

    enddo
  endif
enddo
99 continue
close(7)
close(8)
end subroutine

!**********************************************************************************************

subroutine bin_write(filename,array,n)
use mod_interface
implicit none
integer:: n
real(8):: array(n)
character(40):: filename
open(unit = 120, status='replace',file=trim(filename),form='unformatted')  ! create a new file, or overwrite an existing on
write(120) array
close(120) ! close the file
!!open(unit = 120, status='replace',file=trim(filename),form='formatted')  ! create a new file, or overwrite an existing on
!!write(120,*) array
!!close(120) ! close the file
end subroutine

!**********************************************************************************************

subroutine rotmat_fcart_stress(latvec_init,latvec_trans,rotmat)
!This subroutine will compute a rotation matrix, which transforms
!fcart_trans into the original orientation forces fcart by fcart=matmul(rotmat,fcart_trans)
!stress_trans into the original orientation stress by stress=rotmat*stress_trans*rotnat^T
use mod_interface
implicit none
real(8):: latvec_init(3,3),latvec_trans(3,3),latvec_trans_inv(3,3),rotmat(3,3)
call invertmat(latvec_trans,latvec_trans_inv,3)
rotmat=matmul(latvec_init,latvec_trans_inv)
end subroutine

!**********************************************************************************************

subroutine rotate_stresstensor(strten,rotmat)
!This subroutine will rotate the stress tensor by rotmat according to rotmat*stress*rotmat^T
use mod_interface
implicit none
real(8):: strten(6),rotmat(3,3),stress(3,3)
        stress(1,1) =  strten(1) 
        stress(2,2) =  strten(2) 
        stress(3,3) =  strten(3) 
        stress(2,1) =  strten(6) 
        stress(3,1) =  strten(5) 
        stress(3,2) =  strten(4) 
        stress(1,2) =  stress(2,1)
        stress(1,3) =  stress(3,1)
        stress(2,3) =  stress(3,2)
           stress=matmul(rotmat,matmul(stress,transpose(rotmat)))
        strten(1) =  stress(1,1)
        strten(2) =  stress(2,2)
        strten(3) =  stress(3,3)
        strten(6) =  stress(2,1)
        strten(5) =  stress(3,1)
        strten(4) =  stress(3,2)
end subroutine


!**********************************************************************************************

subroutine bin_read(filename,array,n)
use mod_interface
implicit none
integer:: n
real(8):: array(n)
character(40):: filename
open(unit = 120, status='old',file=trim(filename),form='unformatted')  ! open an existing file
read(120) array ! read the data into array x, of the appropriate data type
close(120) ! close the file
end subroutine

!**********************************************************************************************
subroutine print_logo()
use mod_interface
implicit none  
!!!!     write(*,'(a,23x,a)')' #',' ____           _           _ _        __  __ _       _   _             '
!!!!     write(*,'(a,23x,a)')' #','|  _ \ ___ _ __(_) ___   __| (_) ___  |  \/  (_)_ __ | | | | ___  _ __  '
!!!!     write(*,'(a,23x,a)')' #','| |_) / _ \ /__| |/ _ \ / _  | |/ __| | |\/| | | /_ \| |_| |/ _ \| /_ \ '
!!!!     write(*,'(a,23x,a)')' #','|  __/  __/ |  | | (_) | (_| | | (__  | |  | | | | | |  _  | (_) | |_) |'
!!!!     write(*,'(a,23x,a)')' #','|_|   \___|_|  |_|\___/ \__,_|_|\___| |_|  |_|_|_| |_|_| |_|\___/| .__/ '
!!!!     write(*,'(a,23x,a)')' #','                                                                 |_|    '
!!!!     write(*,'(a,23x,a)')' #','based          _______ __             '
!!!!     write(*,'(a,23x,a)')' #','  ____  ____  / ____(_) /___  _____   '
!!!!     write(*,'(a,23x,a)')' #',' / __ \/ __ \/ /_  / / // _ \/ ___/   '
!!!!     write(*,'(a,23x,a)')' #','/ /_/ / / / / __/ / / //  __(__  )    '
!!!!     write(*,'(a,23x,a)')' #','\____/_/ /_/_/   /_/_/ \___/____/     '
!!!!     write(*,'(a,23x,a)')' #','                                      '
!!!!     write(*,'(a,23x,a)')' #',''
!!!!     write(*,'(a,23x,a)')' #','Redesigned at Basel University,        '
!!!!     write(*,'(a,23x,a)')' #','M. Amsler, S. Goedecker, Computational Physics Group'
!!!!     write(*,'(a,23x,a)')' #',''
!!!!     write(*,'(a,23x,a)')' #','----> you can grep the standard output #'
!!!!     write(*,'(a,23x,a)')' #',' '








write(*,'(a,23x,a)')' #',"  "
write(*,'(a,23x,a)')' #',"       ---_ ......._-_--.  "
write(*,'(a,23x,a)')' #',"      (|\ /      / /| \  \  "
write(*,'(a,23x,a)')' #',"      /  /     .'  -=-'   `.  "
write(*,'(a,23x,a)')' #',"     /  /    .'             )  "
write(*,'(a,23x,a)')' #',"   _/  /   .'        _.)   /  "
write(*,'(a,23x,a)')' #',"  / o   o        _.-' /  .'  "
write(*,'(a,23x,a)')' #',"  \          _.-'    / .'*|  "
write(*,'(a,23x,a)')' #',"   \______.-'//    .'.' \*|  "
write(*,'(a,23x,a)')' #',"    \|  \ | //   .'.' _ |*|  "
write(*,'(a,23x,a)')' #',"     `   \|//  .'.'_ _ _|*|  "
write(*,'(a,23x,a)')' #',"      .  .// .'.' | _ _ \*|      __  __ _       _                             "
write(*,'(a,23x,a)')' #',"      \`-|\_/ /    \ _ _ \*\    |  \/  (_)_ __ | |__   ___   ___ __ _  ___    "
write(*,'(a,23x,a)')' #',"       `/'\__/      \ _ _ \*\   | |\/| | | '_ \| '_ \ / _ \ / __/ _` |/ _ \   "
write(*,'(a,23x,a)')' #',"      /^|            \ _ _ \*   | |  | | | | | | | | | (_) | (__ (_| | (_) |  "
write(*,'(a,23x,a)')' #',"     '  `             \ _ _ \   |_|  |_|_|_| |_|_| |_|\___/ \___\__,_|\___/      "
write(*,'(a,23x,a)')' #',"   "
write(*,'(a,23x,a)')' #',"   "
write(*,'(a,23x,a)')' #',"  Minima Hopping for Crystal Optimization, Version 3.0  "
write(*,'(a,23x,a)')' #',"  "
write(*,'(a,3x,a)')' #',"  Maximilian Amsler  "
write(*,'(a,3x,a)')' #',"  Computational Physics Group  "
write(*,'(a,3x,a)')' #',"  University of Basel  "
write(*,'(a,3x,a)')' #',"    "
write(*,'(a,3x,a)')' #',"  Please cite the following references:  "
write(*,'(a,3x,a)')' #',"  M. Amsler and S. Goedecker: Crystal structure prediction using the minima hopping method,"
write(*,'(a,3x,a)')' #',"  J. Chem. Phys. 133, 224104 (2010)   "
write(*,'(a,3x,a)')' #',"  S. Goedecker: Minima hopping: An efficient search method for the global minimum of the potential energy surface of complex molecular systems,"
write(*,'(a,3x,a)')' #',"  J. Chem. Phys. 120, 9911 (2004)  "
write(*,'(a,3x,a)')' #',"  "
write(*,'(a,3x,a)')' #',"  ----> you can grep the standard output for #'  "
write(*,'(a,3x,a)')' #',"  "

end subroutine
!**********************************************************************************************
