!!MODULE String_Utility 
!!   IMPLICIT NONE 
!!   PRIVATE 
!!   PUBLIC :: StrUpCase 
!!   PUBLIC :: StrLowCase 
!!   CHARACTER( * ), PRIVATE, PARAMETER :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz' 
!!   CHARACTER( * ), PRIVATE, PARAMETER :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' 
!!CONTAINS 
!!   FUNCTION StrUpCase ( Input_String ) RESULT ( Output_String ) 
!!     ! -- Argument and result 
!!     CHARACTER( * ), INTENT( IN )     :: Input_String 
!!     CHARACTER( LEN( Input_String ) ) :: Output_String 
!!     ! -- Local variables 
!!     INTEGER :: i, n 
!!     ! -- Copy input string 
!!     Output_String = Input_String 
!!     ! -- Loop over string elements 
!!     DO i = 1, LEN( Output_String ) 
!!       ! -- Find location of letter in lower case constant string 
!!       n = INDEX( LOWER_CASE, Output_String( i:i ) ) 
!!       ! -- If current substring is a lower case letter, make it upper case 
!!       IF ( n /= 0 ) Output_String( i:i ) = UPPER_CASE( n:n ) 
!!     END DO 
!!   END FUNCTION StrUpCase 
!!   FUNCTION StrLowCase ( Input_String ) RESULT ( Output_String ) 
!!     ! -- Argument and result 
!!     CHARACTER( * ), INTENT( IN )     :: Input_String 
!!     CHARACTER( LEN( Input_String ) ) :: Output_String 
!!     ! -- Local variables 
!!     INTEGER :: i, n 
!!     ! -- Copy input string 
!!     Output_String = Input_String 
!!     ! -- Loop over string elements 
!!     DO i = 1, LEN( Output_String ) 
!!       ! -- Find location of letter in upper case constant string 
!!       n = INDEX( UPPER_CASE, Output_String( i:i ) ) 
!!       ! -- If current substring is an upper case letter, make it lower case 
!!       IF ( n /= 0 ) Output_String( i:i ) = LOWER_CASE( n:n ) 
!!     END DO 
!!   END FUNCTION StrLowCase 
!!END MODULE String_Utility 

subroutine params_read(parini)
use mod_parini, only: typ_parini
use String_Utility
use defs_basis
use interface_ipi
use interface_msock
use mod_fire,   only:dtmin, dtmax
use minpar, only:parmin_bfgs
              
use steepest_descent, only: sd_beta_lat,sd_beta_at
use fingerprint, only: & 
   fp_method,&!All
   fp_12_nl,&                            !CALYPSO parameters
   fp_13_nl,&                            !Modified CALYPSO parameters
   fp_18_large_vanradius



implicit none
type(typ_parini), intent(inout):: parini
integer:: calls=0
integer:: itype,n,k1,k2,l,m,mdmin_in,n_lj,iat,n2,kpt_abc(3),i,j
real(8):: tmp_val,alpha_lat_in,alpha_at_in,dtion_md_in,dkpt_12(2)
character(450):: all_line
character(4):: ch_tmp
character(15):: find_string
character(1):: fn1
character(2):: fn2
character(40):: filename
logical:: found,nat_found,ntypat_found,znucl_found,read_poscur,file_exists
!Init found variables
found=.false.
nat_found=.false.
ntypat_found=.false.
znucl_found=.false.
read_poscur=.false.


!The parser will start reading some parameters first that are needed to allocate some arrays
!These variables are: nat, ntypat, nconfine
parini%nconfine=1
file_exists=.false.
filename="params_new.in"
INQUIRE(FILE=trim(filename), EXIST=file_exists)
if(file_exists) then
  open(unit=12,file="params_new.in")
else
  stop "File params_new.in must be provided!!!"
endif
  do while(.true.)
   read(12,'(a450)',end=99)all_line
   n = len_trim(all_line)
!NAT
   if(.not.nat_found) then
   call parsescalar_int("NAT",3,all_line(1:n),n,parini%nat,found)
   if(found) nat_found=.true.
   if(found) cycle
   endif
!NTYPE
   if(.not.ntypat_found) then
   call parsescalar_int("NTYPE",5,all_line(1:n),n,parini%ntypat_global,found)
   if(found) ntypat_found=.true.
   if(found) cycle
   endif
!ZNUCL
   if(.not.znucl_found) then
   call parsescalar_real("ZNUCL",5,all_line(1:n),n,tmp_val,found)
   if(found) znucl_found=.true.
   if(found) cycle
   endif
!CONFNCONF
   call parsescalar_int("CONFNCONF",9,all_line(1:n),n,parini%nconfine,found)
   if(found) cycle
   enddo
99 continue
close(12)



if(.not.nat_found.or..not.ntypat_found.or..not.znucl_found) then
   write(*,*) "NAT or NTYPE or ZNUCL not found in params_new.in"
!Try to read it from poscur.ascii
   filename='poscur.ascii'
   INQUIRE(FILE=trim(filename), EXIST=file_exists)
   if(file_exists) then
     write(*,*) "Trying to read them from ",filename !Will also read ZNUCL here
     call ascii_getsystem(parini,filename)
     read_poscur=.true.
   endif
   if(.not.read_poscur) then
!Try to read the vasp input file for the system 
   file_exists=.false.
   filename="poscur.vasp"
   INQUIRE(FILE=trim(filename), EXIST=file_exists)
   if(file_exists) then
     write(*,*) "Trying to read them from ",filename !Will also read ZNUCL here
     call poscar_getsystem(parini,filename)
     read_poscur=.true.
   endif
   endif
endif

deallocate(parini%znucl)
deallocate(parini%char_type)
deallocate(parini%amu)
deallocate(parini%rcov)
deallocate(parini%typat_global)
!Allocate the arrays
 if(.not.allocated(parini%znucl))       then;   allocate(parini%znucl(parini%ntypat_global))                       ; parini%znucl=0                 ; endif
 if(.not.allocated(parini%char_type))   then;   allocate(parini%char_type(parini%ntypat_global))                   ; parini%char_type="  "          ; endif
 if(.not.allocated(parini%amu))         then;   allocate(parini%amu(parini%ntypat_global))                         ; parini%amu=0                   ; endif
 if(.not.allocated(parini%rcov))        then;   allocate(parini%rcov(parini%ntypat_global))                        ; parini%rcov=0                  ; endif
 if(.not.allocated(parini%typat_global))       then;   allocate(parini%typat_global(parini%nat))                          ;
     parini%typat_global=0                 ; endif
 if(.not.allocated(parini%fixat))       then;   allocate(parini%fixat(parini%nat))                          ; parini%fixat=.false.           ; endif
 if(.not.allocated(parini%fragarr))     then;   allocate(parini%fragarr(parini%nat))                        ; parini%fragarr=0               ; endif
 if(.not.allocated(parini%conf_dim))    then;   allocate(parini%conf_dim     (parini%nconfine))             ; parini%conf_dim=0              ; endif
 if(.not.allocated(parini%conf_av))     then;   allocate(parini%conf_av      (parini%nconfine))             ; parini%conf_av=0               ; endif
 if(.not.allocated(parini%conf_exp))    then;   allocate(parini%conf_exp     (parini%nconfine))             ; parini%conf_exp=0              ; endif
 if(.not.allocated(parini%conf_prefac)) then;   allocate(parini%conf_prefac  (parini%nconfine))             ; parini%conf_prefac=0           ; endif
 if(.not.allocated(parini%conf_cut))    then;   allocate(parini%conf_cut     (parini%nconfine))             ; parini%conf_cut=0              ; endif
 if(.not.allocated(parini%conf_eq))     then;   allocate(parini%conf_eq      (parini%nconfine))             ; parini%conf_eq=0               ; endif
 if(.not.allocated(parini%conf_list))   then;   allocate(parini%conf_list    (parini%nat,parini%nconfine))         ; parini%conf_list=0             ; endif
 if(.not.allocated(parini%conf_nat))    then;   allocate(parini%conf_nat     (parini%nconfine))             ; parini%conf_nat=0              ; endif
 if(.not.allocated(parini%conf_cartred))then;   allocate(parini%conf_cartred (parini%nconfine))             ; parini%conf_cartred="C"        ; endif


!Feed arrays with standard parameters
call params_defaults(parini,mdmin_in,dtion_md_in,alpha_lat_in,alpha_at_in,read_poscur)

!Read znucl
if(.not.read_poscur) then
   open(unit=12,file="params_new.in")
   found=.false.
   do while(.true.)
      read(12,'(a450)',end=98)all_line
      n = len_trim(all_line)
   !Znucl
      call parsearray_real("ZNUCL",5,all_line(1:n),n,parini%znucl(1:parini%ntypat_global),parini%ntypat_global,found)
      if(found) exit
   enddo
   98 continue
   close(12)
   if(.not.found) then
     stop "ZNUCL not found in params_new.in/poscur.ascii/poscur.vasp"
   endif
endif

!Get the correct atomic masses and atomic character
 do itype=1,parini%ntypat_global
   call atmdata(parini%amu(itype),parini%rcov(itype),parini%char_type(itype),parini%znucl(itype))
 enddo

!Read the other variables
open(unit=12,file="params_new.in")
  do while(.true.)
   read(12,'(a450)',end=97)all_line
   n = len_trim(all_line)
!Press
   call parsescalar_real("PRESS",5,all_line(1:n),n,parini%target_pressure_gpa,found)
   if(found) parini%target_pressure_habohr=parini%target_pressure_gpa/HaBohr3_GPA
   if(found) cycle
!Amu
   call parsearray_real("AMU",3,all_line(1:n),n,parini%amu(1:parini%ntypat_global),parini%ntypat_global,found)
   if(found) cycle
!TYPAT
   call parsearray_int("TYPAT",5,all_line(1:n),n,parini%typat_global(1:parini%nat),parini%nat,found)
   if(found) cycle
!MDNIT
   call parsescalar_int("MDNIT",5,all_line(1:n),n,parini%nmd_dynamics,found)
   if(found) cycle
!MDALGO
   call parsescalar_int("MDALGO",6,all_line(1:n),n,parini%md_algo,found)
   if(found) cycle
!MDINTEGRATOR
   call parsescalar_int("MDINT",5,all_line(1:n),n,parini%md_integrator,found)
   if(found) cycle
!MDPRESSCOMP
   call parsescalar_real("MDPRESSCOMP",11,all_line(1:n),n,parini%md_presscomp,found)
   if(found) cycle
!GEONIT
   call parsescalar_int("GEONIT",6,all_line(1:n),n,parini%paropt_geopt%nit,found)
   if(found) cycle
!CELLMASS
   call parsescalar_real("CELLMASS",8,all_line(1:n),n,parini%bmass,found)
   if(found) cycle
!VOIDS
   call parse_logical("VOIDS",5,all_line(1:n),n,parini%voids,found)
   if(found) cycle
!COREREP
   call parse_logical("COREREP",7,all_line(1:n),n,parini%core_rep,found)
   if(found) cycle
!Block mdmin****************
!AUTO_MDMIN
   call parse_logical("AUTO_MDMIN",10,all_line(1:n),n,parini%auto_mdmin,found)
   if(found) cycle
!MDMININIT
   call parsescalar_int("MDMININIT",9,all_line(1:n),n,mdmin_in,found)
   if(found) cycle
!MDMINMIN
   call parsescalar_int("MDMINMIN",8,all_line(1:n),n,parini%mdmin_min,found)
   if(found) cycle
!MDMINMAX
   call parsescalar_int("MDMINMAX",8,all_line(1:n),n,parini%mdmin_max,found)
   if(found) cycle
!Block mdmin****************
!Block soften****************
!AUTO_SOFT
   call parse_logical("AUTO_SOFT",9,all_line(1:n),n,parini%auto_soft,found)
   if(found) cycle
!SOFTLAT
   call parsescalar_real("SOFTLAT",7,all_line(1:n),n,alpha_lat_in,found)
   if(found) cycle
!SOFTAT
   call parsescalar_real("SOFTAT",6,all_line(1:n),n,alpha_at_in,found)
   if(found) cycle
!MOLSOFT
   call parse_logical("MOLSOFT",7,all_line(1:n),n,parini%mol_soften,found)
   if(found) cycle
!SOFTNIT
   call parsescalar_int("SOFTNIT",7,all_line(1:n),n,parini%nsoften_minhopp,found)
   if(found) cycle
!Block soften****************
!Block MD timestep***********
!AUTO_MDDT
   call parse_logical("AUTO_MDDT",9,all_line(1:n),n,parini%auto_dtion_md,found)
   if(found) cycle
!MDDTINIT
   call parsescalar_real("MDDTINIT",8,all_line(1:n),n,dtion_md_in,found)
   if(found) cycle
!MDDTIPM
   call parsescalar_int("MDDTIPM",7,all_line(1:n),n,parini%nit_per_min,found)
   if(found) cycle
!MDENCON
   call parse_logical("MDENCON",7,all_line(1:n),n,parini%energy_conservation,found)
   if(found) cycle
!Block MD timestep***********
!Block GEOPT*****************
!GEOALGO
   call parsescalar_string("GEOALGO",7,all_line(1:n),n,parini%paropt_geopt%approach,5,found)
   if(found) parini%paropt_geopt%approach=StrUpCase(parini%paropt_geopt%approach)
   if(found) cycle
!GEOFIREDTINIT
   call parsescalar_real("GEOFIREDTINIT",13,all_line(1:n),n,parini%paropt_geopt%dt_start,found)
   if(found) cycle
!GEOFIREDTMIN
   call parsescalar_real("GEOFIREDTMIN",12,all_line(1:n),n,parini%paropt_geopt%dtmin,found)
   if(found) cycle
!GEOFIREDTMAX
   call parsescalar_real("GEOFIREDTMAX",12,all_line(1:n),n,parini%paropt_geopt%dtmax,found)
   if(found) cycle
!GEOHESSLAT
   call parsescalar_real("GEOHESSLAT",10,all_line(1:n),n,parini%alphax_lat,found)
   if(found) cycle
!GEOHESSAT
   call parsescalar_real("GEOHESSAT",9,all_line(1:n),n,parini%alphax_at,found)
   if(found) cycle
!GEOTOLMXF
   call parsescalar_real("GEOTOLMXF",9,all_line(1:n),n,parini%paropt_geopt%fmaxtol,found)
   if(found) cycle
!GEOSTRFACT
   call parsescalar_real("STRFACT",7,all_line(1:n),n,parini%paropt_geopt%strfact,found)
   if(found) cycle
!GEOEXT
   call parse_logical("GEOEXT",6,all_line(1:n),n,parini%geopt_ext,found)
   if(found) cycle
!GEOSQNMNHIS
   call parsescalar_int ("GEOSQNMNHIST",12,all_line(1:n),n,parini%paropt_geopt%nhist,found)
   if(found) cycle
!GEOSQNMMAXRISE
   call parsescalar_real("GEOSQNMMAXRISE",14,all_line(1:n),n,parini%paropt_geopt%maxrise,found)
   if(found) cycle
!GEOSQNMCUTOFF
   call parsescalar_real("GEOSQNMCUTOFF",13,all_line(1:n),n,parini%paropt_geopt%cutoffRatio,found)
   if(found) cycle
!GEOSQNMSTEEP
   call parsescalar_real("GEOSQNMSTEEP",12,all_line(1:n),n,parini%paropt_geopt%steepthresh,found)
   if(found) cycle
!GEOSQNMTRUSTR
   call parsescalar_real("GEOSQNMTRUSTR",13,all_line(1:n),n,parini%paropt_geopt%trustr,found)
   if(found) cycle
!GEOQBFGSNDIM
   call parsescalar_int ("GEOQBFGSNDIM",12,all_line(1:n),n,parini%qbfgs_bfgs_ndim,found)
   if(found) cycle
!GEOQBFGSTRTI
   call parsescalar_real("GEOQBFGSTRI",11,all_line(1:n),n,parini%qbfgs_trust_radius_ini,found)
   if(found) cycle
!GEOQBFGSTRTMIN
   call parsescalar_real("GEOQBFGSTRMIN",13,all_line(1:n),n,parini%qbfgs_trust_radius_min,found)
   if(found) cycle
!GEOQBFGSTRTMAX
   call parsescalar_real("GEOQBFGSTRMAX",13,all_line(1:n),n,parini%qbfgs_trust_radius_max,found)
   if(found) cycle
!GEOQBFGSW1
   call parsescalar_real("GEOQBFGSW1",10,all_line(1:n),n,parini%qbfgs_w_1,found)
   if(found) cycle
!GEOQBFGSW2
   call parsescalar_real("GEOQBFGSW1",10,all_line(1:n),n,parini%qbfgs_w_2,found)
   if(found) cycle
!Block GEOPT*****************
!USEWFGEO
   call parse_logical("USEWFGEO",8,all_line(1:n),n,parini%usewf_geopt,found)
   if(found) cycle
!USEWFSOFT
   call parse_logical("USEWFSOFT",9,all_line(1:n),n,parini%usewf_soften,found)
   if(found) cycle
!USEWFMD
   call parse_logical("USEWFMD",7,all_line(1:n),n,parini%usewf_md,found)
   if(found) cycle
!FINDSYM
   call parse_logical("FINDSYM",7,all_line(1:n),n,parini%findsym,found)
   if(found) cycle
!FINDDOS
   call parse_logical("FINDDOS",7,all_line(1:n),n,parini%finddos,found)
   if(found) cycle
!Block KPT****************
!AUTO_KPT
   call parse_logical("AUTO_KPT",8,all_line(1:n),n,parini%auto_kpt,found)
   if(found) cycle
!KPTMESH
   call parsearray_int("KPTMESH",7,all_line(1:n),n,kpt_abc(1:3),3,found)
   if(found) then
     parini%ka=kpt_abc(1)
     parini%kb=kpt_abc(2)
     parini%kc=kpt_abc(3)
   endif
   if(found) cycle
!KPTDEN
   call parsearray_real("KPTDEN",6,all_line(1:n),n,dkpt_12(1:2),2,found)
   if(found) then
     parini%dkpt1=dkpt_12(1)
     parini%dkpt2=dkpt_12(2)
   endif
   if(found) cycle
!Block KPT****************
!Block FINGERPRINT****************
!FPMETHOD
   call parsescalar_string("FPMETHOD",8,all_line(1:n),n,parini%fp_method_ch,20,found)
   if(found) parini%fp_method_ch=StrUpCase(parini%fp_method_ch)
   if(found) cycle
!FPCUT
   call parsescalar_real("FPCUT",5,all_line(1:n),n,parini%fp_rcut,found)
   if(found) cycle
!FPDBIN
   call parsescalar_real("FPDBIN",6,all_line(1:n),n,parini%fp_dbin,found)
   if(found) cycle
!FPSIGMA
   call parsescalar_real("FPSIGMA",7,all_line(1:n),n,parini%fp_sigma,found)
   if(found) cycle
!FPNL
   call parsescalar_int("FPNL",4,all_line(1:n),n,parini%fp_nl,found)
   if(found) cycle
!FPPOWER
   call parsescalar_int("FPPOWER",7,all_line(1:n),n,parini%fp_14_m,found)
   if(found) cycle
!FPGAUSSFAC1
   call parsescalar_real("FPGAUSSFAC1",11,all_line(1:n),n,parini%fp_14_w1,found)
   if(found) cycle
!FPGAUSSFAC2
   call parsescalar_real("FPGAUSSFAC2",11,all_line(1:n),n,parini%fp_14_w2,found)
   if(found) cycle
!FPNMAX
   call parsescalar_int("FPATNMAX",8,all_line(1:n),n,parini%fp_at_nmax,found)
   if(found) cycle

!FPWIDTHCUT
!   call parsescalar_real("FPWIDTHCUT",10,all_line(1:n),n,fp_17_width_cutoff,found)
!FPNATX
   call parsescalar_int("FPNATX",6,all_line(1:n),n,parini%fp_17_natx_sphere,found)
   if(found) cycle
!FPORBITAL
   call parsescalar_string("FPORBITAL",9,all_line(1:n),n,parini%fp_17_orbital,2,found)
   parini%fp_18_orbital=parini%fp_17_orbital
   if(found) cycle
!FPNEXCUT
   call parsescalar_real("FPNEXCUT",8,all_line(1:n),n,parini%fp_17_nex_cutoff,found)
   parini%fp_18_nex_cutoff=int(parini%fp_17_nex_cutoff)
   if(found) cycle
!FPPRINCIPLEEV
   call parsescalar_int("FPPRINCIPLEEV",13,all_line(1:n),n,parini%fp_18_principleev,found)
   if(found) cycle
!FPMOLECULES
   call parsescalar_int("FPMOLECULES",11,all_line(1:n),n,parini%fp_18_molecules,found)
   if(found) cycle
!FPEXPA
   call parsescalar_int("FPEXPA",6,all_line(1:n),n,parini%fp_18_expaparameter,found)
   if(found) cycle
!FPMOLSPHERE
   call parsescalar_int("FPMOLSPHERE",11,all_line(1:n),n,parini%fp_18_molecules_sphere,found)
   if(found) cycle
!FPWIDTHCUT
   call parsescalar_real("FPWIDTHCUT",10,all_line(1:n),n,parini%fp_18_width_cutoff,found)
   if(found) cycle
!FPWIDTHOVER
   call parsescalar_real("FPWIDTHOVER",11,all_line(1:n),n,parini%fp_18_width_overlap,found)
   if(found) cycle
!Block FINGERPRINT****************


!Block LAYER-CONFINEMENT****************
!CONFINEMENT
   call parse_logical("CONFINEMENT",11,all_line(1:n),n,parini%use_confine,found)
   if(found) cycle
!CONFCARTRED
   call parsearray_string("CONFCARTRED",11,all_line(1:n),n,parini%conf_cartred,1,parini%nconfine,found)
   if(found) cycle
!CONFDIM
   call parsearray_int("CONFDIM",7,all_line(1:n),n,parini%conf_dim(1:parini%nconfine),parini%nconfine,found)
   if(found) cycle
!CONFEXP
   call parsearray_int("CONFEXP",7,all_line(1:n),n,parini%conf_exp(1:parini%nconfine),parini%nconfine,found)
   if(found) cycle
!CONFPREFAC
   call parsearray_real("CONFPREFAC",10,all_line(1:n),n,parini%conf_prefac(1:parini%nconfine),parini%nconfine,found)
   if(found) cycle
!CONFCUT
   call parsearray_real("CONFCUT",7,all_line(1:n),n,parini%conf_cut(1:parini%nconfine),parini%nconfine,found)
   if(found) cycle
!CONFAV
   call parsearray_int("CONFAV",6,all_line(1:n),n,parini%conf_av(1:parini%nconfine),parini%nconfine,found)
   if(found) cycle
!CONFEQ
   call parsearray_real("CONFEQ",6,all_line(1:n),n,parini%conf_eq(1:parini%nconfine),parini%nconfine,found)
   if(found) cycle
!CONFNAT
   call parsearray_int("CONFNAT",7,all_line(1:n),n,parini%conf_nat(1:parini%nconfine),parini%nconfine,found)
   if(found) cycle
!CONFLIST#
!Go through all number of confinements
   do i=1,parini%nconfine
      if(i.lt.10) then
        write(fn1,'(i1.1)') i
        write(find_string,'(a,a)') "CONFLIST"//fn1 
        call exist_string(trim(find_string),9,all_line(1:n),n,found)
        if(found.and.parini%conf_nat(i)==parini%nat) then
          do j=1,parini%nat
            parini%conf_list(j,i)=j
          enddo
        else
          call parsearray_int(trim(find_string),9,all_line(1:n),n,parini%conf_list(1:parini%conf_nat(i),i),parini%conf_nat(i),found)
        endif
      else
        stop "Loop not implemented for nconfine greater than 9"
      endif
   enddo
!Block LAYER-CONFINEMENT****************

!Block IPI_SOCKET****************
!IPIINET
   call parsescalar_int("IPIINET",7,all_line(1:n),n,parini%sock_inet,found)
   if(found) cycle
!IPIPORT
   call parsescalar_int("IPIPORT",7,all_line(1:n),n,parini%sock_port,found)
   if(found) cycle
!IPIHOST
   call parsescalar_string("IPIHOST",7,all_line(1:n),n,parini%sock_host,1024,found)
   if(found) cycle
!IPIECUTWF
   call parsearray_real("IPIECUTWF",9,all_line(1:n),n,parini%sock_ecutwf,2,found)
   if(found) cycle
!Block IPI_SOCKET****************

!Block MSOCK****************
!MSOCKINET
   call parsescalar_int("SOCKINET",8,all_line(1:n),n,parini%sock_inet,found)
   if(found) cycle
!MSOCKPORT
   call parsescalar_int("SOCKPORT",8,all_line(1:n),n,parini%sock_port,found)
   if(found) cycle
!MSOCKHOST
   call parsescalar_string("SOCKHOST",8,all_line(1:n),n,parini%sock_host,1024,found)
   if(found) cycle
!MSOCKECUTWF
   call parsearray_real("SOCKECUTWF",10,all_line(1:n),n,parini%sock_ecutwf,2,found)
   if(found) cycle

!Block MSOCK****************

!VERBOSE
   call parsescalar_int("VERBOSE",7,all_line(1:n),n,parini%verb,found)
   if(found) cycle
!BOUNDARY
   call parsescalar_int("BOUNDARY",8,all_line(1:n),n,parini%bc,found)
   if(found) cycle
!CODE
   call parsescalar_string("CODE",4,all_line(1:n),n,parini%potential_potential,20,found)
   if(found) parini%potential_potential = StrLowCase( parini%potential_potential )
   if(found) cycle
  enddo
97 continue
close(12)


!Post Processing
!MDMIN
  if(calls==0) then
      if(.not.parini%auto_mdmin) then
           parini%mdmin=mdmin_in
      else
           parini%mdmin=max(mdmin_in,parini%mdmin_min)
      endif
  endif
!SOFTEN
  if(calls==0.or..not.parini%auto_soft) parini%alpha_lat=alpha_lat_in
  if(calls==0.or..not.parini%auto_soft) parini%alpha_at=alpha_at_in
!MDTIMESTE=P
  if(calls==0.or..not.parini%auto_dtion_md) parini%dtion_md=dtion_md_in
!KPT
  if(parini%auto_kpt) then
    parini%ka=0;parini%kb=0;parini%kc=0
  else
    parini%dkpt1=0.d0
    parini%dkpt2=0.d0
  endif

!Initiallize confinement
if(parini%use_confine) call  init_confinement_parser(parini)

!Initiallize LJ parameter if required
if(trim(parini%potential_potential)=="blj".and.calls==0) call blj_init_parameter(parini)

!Initiallize LJ parameter if required
if(trim(parini%potential_potential)=="mlj") call mlj_init_parameter(parini)

!Initiallize TB-LJ parameter if required
if(trim(parini%potential_potential)=="lenosky_tb_lj".and.calls==0) then
  call check_lenosky_tb_lj(parini)
  n_lj=0
  do iat=1,parini%nat
     if(parini%znucl(parini%typat_global(iat)).gt.200) n_lj=n_lj+1
  enddo
  call lenosky_tb_lj_init_parameter()
endif

!Initiallize voids
if(parini%voids.and.calls==0) then
  call check_voids(parini) 
  call voids_init_parameter()
endif

!Initiallize tersoff
if(trim(parini%potential_potential)=="tersoff".and.calls==0) then
  call init_tersoff(parini)
endif

!Initiallize edip
if(trim(parini%potential_potential)=="edip".and.calls==0) then
  call init_edip(parini)
endif

!Initiallize ipi
if(trim(parini%potential_potential)=="ipi".and.calls==0) then
  call init_ipi(parini,parini%nat)
endif

!Initiallize msock
if(trim(parini%potential_potential)=="msock".and.calls==0) then
  call init_msock(parini)
endif



!Assign correct parameters for fingerprint, not allocating any arrays!
call fp_assign(parini)

!Increase calls to the routine
!if(calls==0) call params_echo(parini)
if(calls==0) call params_check(parini)
calls=calls+1

!Check range of parameters
call params_check(parini)

!Copy parameters of fire to the fire module
    dtmin=parini%paropt_geopt%dtmin
    dtmax=parini%paropt_geopt%dtmax
!Copy parameters of bfgs to the bfgs module
    parmin_bfgs%betax=parini%alphax_at
    parmin_bfgs%betax_lat=parini%alphax_lat
!Copy parameters to sqnm module
    parini%paropt_geopt%beta_lat=parini%alphax_lat
    parini%paropt_geopt%beta_at=parini%alphax_at
!Copy parameters to sd module
    sd_beta_lat=parini%alphax_lat
    sd_beta_at=parini%alphax_at
end subroutine
!************************************************************************************
subroutine params_read_for_yaml(parini)
    use mod_parini, only: typ_parini
    use mod_fire,   only: dtmin, dtmax
    use minpar, only: parmin_bfgs
    use steepest_descent, only: sd_beta_lat,sd_beta_at
    use interface_ipi
    use interface_msock
    use String_Utility
    implicit none
    type(typ_parini), intent(inout):: parini
    integer, save:: calls=0
    integer:: iat, n_lj

    !Upcase important parameter:
    parini%paropt_geopt%approach = StrUpCase ( parini%paropt_geopt%approach )
    parini%fp_method_ch = StrUpCase ( parini%fp_method_ch )
    parini%fp_17_orbital = StrUpCase ( parini%fp_17_orbital )
    parini%fp_18_orbital = StrUpCase ( parini%fp_18_orbital )

!Post Processing
!KPT
  if(parini%auto_kpt) then
    parini%ka=0;parini%kb=0;parini%kc=0
  else
    parini%dkpt1=0.d0
    parini%dkpt2=0.d0
  endif

!Initiallize confinement
if(parini%use_confine) call  init_confinement_parser(parini)

!Initiallize LJ parameter if required
if(trim(parini%potential_potential)=="blj".and.calls==0) call blj_init_parameter(parini)

!Initiallize LJ parameter if required
if(trim(parini%potential_potential)=="mlj") call mlj_init_parameter(parini)

!Initiallize TB-LJ parameter if required
if(trim(parini%potential_potential)=="lenosky_tb_lj".and.calls==0) then
  call check_lenosky_tb_lj(parini)
  n_lj=0
  do iat=1,parini%nat
     if(parini%znucl(parini%typat_global(iat)).gt.200) n_lj=n_lj+1
  enddo
  call lenosky_tb_lj_init_parameter()
endif

!Initiallize voids
if(parini%voids.and.calls==0) then
  call check_voids(parini) 
  call voids_init_parameter()
endif

!Initiallize tersoff
if(trim(parini%potential_potential)=="tersoff".and.calls==0) then
  call init_tersoff(parini)
endif

!Initiallize edip
if(trim(parini%potential_potential)=="edip".and.calls==0) then
  call init_edip(parini)
endif

!Initiallize ipi
if(trim(parini%potential_potential)=="ipi".and.calls==0) then
  call init_ipi(parini,parini%nat)
endif

!Initiallize msock
if(trim(parini%potential_potential)=="msock".and.calls==0) then
  call init_msock(parini)
endif



!Assign correct parameters for fingerprint, not allocating any arrays!
call fp_assign(parini)

!Increase calls to the routine
!if(calls==0) call params_echo(parini)
!if(calls==0) call params_check(parini)
calls=calls+1

!Copy parameters of fire to the fire module
    dtmin=parini%paropt_geopt%dtmin
    dtmax=parini%paropt_geopt%dtmax
!Copy parameters of bfgs to the bfgs module
    parmin_bfgs%betax=parini%alphax_at
    parmin_bfgs%betax_lat=parini%alphax_lat
!Copy parameters to sqnm module
    parini%paropt_geopt%beta_lat=parini%alphax_lat
    parini%paropt_geopt%beta_at=parini%alphax_at
!Copy parameters to sd module
    sd_beta_lat=parini%alphax_lat
    sd_beta_at=parini%alphax_at

    !call params_echo(parini)
    !Check range of parameters


    call params_check(parini)
end subroutine params_read_for_yaml
!************************************************************************************
subroutine params_defaults(parini,mdmin_in,dtion_md_in,alpha_lat_in,alpha_at_in,read_poscur)
use defs_basis
use mod_fire,   only:dtmin, dtmax
use minpar, only:parmin_bfgs
use fingerprint, only: & 
   fp_method,&!All
   fp_12_nl,&                            !CALYPSO parameters
   fp_13_nl,&                            !Modified CALYPSO parameters
   fp_18_large_vanradius
   
use mod_parini, only: typ_parini
implicit none
type(typ_parini), intent(inout):: parini
integer:: mdmin_in,itype,i,j
real(8):: dtion_md_in,alpha_lat_in,alpha_at_in
logical:: read_poscur
!These are the default variables
parini%target_pressure_gpa=0.d0
parini%target_pressure_habohr=0.d0
!Get the correct atomic masses and atomic character
 do itype=1,parini%ntypat_global
   call atmdata(parini%amu(itype),parini%rcov(itype),parini%char_type(itype),parini%znucl(itype))
 enddo
parini%voids=.false.
parini%core_rep=.false.
if(.not.read_poscur) parini%typat_global(1:parini%nat)=1
parini%nmd_dynamics=300
parini%md_algo=1
parini%md_integrator=3
parini%md_presscomp=-0.d0
parini%paropt_geopt%nit=300
parini%bmass=1.d0
parini%auto_mdmin=.false.
mdmin_in=1
parini%mdmin_min=2
parini%mdmin_max=2
parini%energy_conservation=.false.
parini%auto_soft=.false.
alpha_lat_in=1.d0
alpha_at_in=1.d0
parini%mol_soften=.false.
parini%nsoften_minhopp=10
parini%auto_dtion_md=.false.
dtion_md_in=20.d0
parini%nit_per_min=25.d0
parini%paropt_geopt%approach="FIRE"
parini%paropt_geopt%dt_start=10.d0
parini%paropt_geopt%dtmin=1.d0
parini%paropt_geopt%dtmax=80.d0
parini%alphax_lat=1.d0
parini%alphax_at=1.d0
parini%paropt_geopt%fmaxtol=2.d-4
parini%paropt_geopt%strfact=100.d0
parini%usewf_geopt=.false.
parini%usewf_soften=.false.
parini%usewf_md=.false.
parini%findsym=.false.
parini%finddos=.false.
parini%auto_kpt=.true.
parini%ka=1;parini%kb=1;parini%kc=1
parini%dkpt1=0.04d0
parini%dkpt2=0.06d0
parini%bc=1
parini%verb=3
parini%potential_potential="vasp"
!Define if the external optimizer should be used. Only available for:
parini%geopt_ext=.false.

parini%fp_rcut=15.d0
fp_method=11
parini%fp_method_ch="OGANOV"
parini%fp_nl=6
parini%fp_sigma=0.02d0
parini%fp_dbin= 0.05d0
fp_12_nl=6
fp_13_nl=6
parini%fp_14_m=3
parini%fp_14_w1=1.d0
parini%fp_14_w2=1.5d0
parini%fp_at_nmax=10000
parini%fp_17_nex_cutoff=3
parini%fp_17_width_cutoff=parini%fp_rcut/sqrt(2.d0*parini%fp_17_nex_cutoff)
parini%fp_17_orbital='S'
parini%fp_17_lseg=1
parini%fp_17_natx_sphere=75

parini%fp_18_orbital='S'
parini%fp_18_principleev = 6
parini%fp_18_lseg=1
parini%fp_18_molecules=1
parini%fp_18_expaparameter = 4
parini%fp_18_nex_cutoff = 3
parini%fp_18_molecules_sphere = 50
parini%fp_18_width_cutoff = 1.d0
parini%fp_18_width_overlap = 1.d0
fp_18_large_vanradius = 1.7d0/0.52917720859d0

!SQNM
parini%paropt_geopt%beta_lat=1.d0
parini%paropt_geopt%beta_at=1.d0
parini%paropt_geopt%nhist=10
parini%paropt_geopt%maxrise=1.d-6
parini%paropt_geopt%cutoffRatio=1.d-4
parini%paropt_geopt%steepthresh=1.d0
parini%paropt_geopt%trustr=0.1d0

!QBFGS
parini%qbfgs_bfgs_ndim=1
parini%qbfgs_trust_radius_max=0.5d0
parini%qbfgs_trust_radius_min=1.d-3
parini%qbfgs_trust_radius_ini=0.5D0
parini%qbfgs_w_1=0.01D0
parini%qbfgs_w_2=0.5D0

parini%use_confine=.false.
parini%conf_cartred="C"
parini%conf_dim=1
parini%conf_exp=4
parini%conf_prefac=1.d-2
parini%conf_cut=1.d0
parini%conf_av=2
parini%conf_eq=0
parini%conf_nat=parini%nat
   do i=1,parini%nconfine
          do j=1,parini%nat
            parini%conf_list(j,i)=j
          enddo
   enddo
!Block LAYER-CONFINEMENT****************
!!Block IPI_SOCKET****************
!      ipi_inet=0 !0 for unix socket, 1 for tcp
!      ipi_port=3141
!      ipi_host="mh-driver"
!      ipi_ecutwf=1.d0
!Block IPI_SOCKET****************
      parini%sock_inet=0 !0 for unix socket, 1 for tcp
      parini%sock_port=3141
      parini%sock_host="mh-driver"
      parini%sock_ecutwf=1.d0
end subroutine

!************************************************************************************

subroutine params_check(parini)
use defs_basis
use mod_fire,   only:dtmin, dtmax
use minpar, only:parmin_bfgs
use fingerprint, only: & 
   fp_method,&!All
   fp_12_nl,&                            !CALYPSO parameters
   fp_13_nl,&                            !Modified CALYPSO parameters
   fp_18_large_vanradius
use mod_parini, only: typ_parini
implicit none
type(typ_parini), intent(in):: parini
integer:: i,j
!This routine will check if all values are valid...
if(any(parini%typat_global(:).lt.1)) stop "Error in typat"
if(any(parini%amu(:).le.0.d0)) stop "Error in amu"
if(any(parini%rcov(:).le.0.d0)) stop "Error in rcov"
if(any(parini%znucl(:).le.0)) stop "Error in znucl"
if(parini%nmd_dynamics.lt.1) stop "Error in parini%nmd_dynamics"
if(parini%md_algo.lt.1.or.parini%md_algo.gt.4) stop "Error in parini%md_algo"
if(parini%md_integrator.lt.1.or.parini%md_integrator.gt.3) stop "Error in parini%md_integrator"
if(parini%paropt_geopt%nit.lt.0) stop "Error in parini%paropt_geopt%nit"
if(parini%bmass.le.0.d0) stop "Error in bmass"
if(parini%mdmin_min.lt.0) stop "Error in parini%mdmin_min"
if(parini%mdmin_max.lt.parini%mdmin_min) stop "Error in parini%mdmin_max"
if(parini%alpha_lat.le.0.d0) stop "Error in alpha_lat"
if(parini%alpha_at.le.0.d0) stop "Error in alpha_at"
if(parini%nsoften_minhopp.lt.1) stop "Error in nsoften"
if(parini%dtion_md.le.0.d0) stop "Error in parini%dtion_md"
if(parini%nit_per_min.le.0) stop "Error in parini%nit_per_min"
if(trim(parini%paropt_geopt%approach).ne."FIRE".and.&
   trim(parini%paropt_geopt%approach).ne."MBFGS".and.&
   trim(parini%paropt_geopt%approach).ne."RBFGS".and.&
   trim(parini%paropt_geopt%approach).ne."SQNM".and.&
   trim(parini%paropt_geopt%approach).ne."QBFGS".and.&
   trim(parini%paropt_geopt%approach).ne."SD") &
   stop "Error in parini%paropt_geopt%approach"
if(parini%paropt_geopt%dt_start.lt.parini%paropt_geopt%dtmin.or.parini%paropt_geopt%dt_start.gt.parini%paropt_geopt%dtmax) stop "Error in parini%paropt_geopt%dt_start"
if(parini%paropt_geopt%dtmin.le.0.d0) stop "Error in parini%paropt_geopt%dtmin"
if(parini%paropt_geopt%dtmax.lt.parini%paropt_geopt%dtmin) stop "Error in parini%paropt_geopt%dtmax"
if(parini%alphax_lat.le.0.d0) stop "Error in alphax_lat"
if(parini%alphax_at.le.0.d0) stop "Error in alphax_at"
if(parini%paropt_geopt%fmaxtol.le.0.d0) stop "Error in parini%paropt_geopt%fmaxtol"
if(parini%paropt_geopt%strfact.le.0.d0) stop "Error in parini%paropt_geopt%strfact"
if(parini%ka.lt.0) stop "Error in ka"
if(parini%kb.lt.0) stop "Error in kb"
if(parini%kc.lt.0) stop "Error in kc"
if(parini%dkpt1.lt.0.d0) stop "Error in dkpt1"
if(parini%dkpt2.lt.0.d0) stop "Error in dkpt2"
if(parini%bc.lt.1.or.parini%bc.gt.3) stop "Error in bc"
if(parini%verb.lt.0.or.parini%verb.gt.3) stop "Error in verb"
if(trim(parini%fp_method_ch).ne."OGANOV".and.trim(parini%fp_method_ch).ne."BCM".and.trim(parini%fp_method_ch).ne."ATORB".and.&
  &trim(parini%fp_method_ch).ne."XYZ2SM".and.trim(parini%fp_method_ch).ne."GAUSS".and.trim(parini%fp_method_ch).ne."COGANOV".and.&
  &trim(parini%fp_method_ch).ne."CAOGANOV".and.trim(parini%fp_method_ch).ne."GOM".and.trim(parini%fp_method_ch).ne."MOLGOM") stop "Error in fp_method_ch"
if(parini%fp_rcut.le.0.d0) stop "Error in fp_rcut"
if(parini%fp_dbin.le.0.d0) stop "Error in fp_dbin"
if(parini%fp_sigma.le.0.d0) stop "Error in fp_sigma"
if(parini%fp_nl.le.0) stop "Error in fp_nl"
if(parini%fp_14_m.lt.1) stop "Error in fp_14_m"
if(parini%fp_14_w1.lt.0.d0) stop "Error in p_14_w1"
if(parini%fp_14_w2.lt.parini%fp_14_w1) stop "Error in p_14_w2"
if(parini%fp_at_nmax.lt.0) stop "Error in fp_at_nmax"
if(trim(parini%fp_17_orbital).ne.'S'.and.trim(parini%fp_17_orbital).ne.'SP') stop "Error in fp_17_orbital"
if(trim(parini%fp_18_orbital).ne.'S'.and.trim(parini%fp_18_orbital).ne.'SP') stop "Error in fp_17_orbital"
if(parini%fp_18_principleev.lt.0) stop "Error in fp_18_principleev"
if(parini%fp_18_molecules.lt.1) stop "Error in fp_18_molecules"
if(parini%fp_18_expaparameter.lt.1) stop "Error in fp_18_expaparameter"
if(parini%fp_18_nex_cutoff.lt.1) stop "Error in fp_18_nex_cutoff"
if(parini%fp_18_molecules_sphere.lt.0) stop "Error in fp_18_molecules_sphere"
if(parini%fp_18_width_cutoff.lt.0.d0) stop "Error in fp_18_width_cutoff"
if(parini%fp_18_width_overlap.lt.0.d0) stop "Error in fp_18_width_overlap"
do i=1,parini%nconfine
  if(.not.(parini%conf_cartred(i).eq."C".or.parini%conf_cartred(i).eq."c".or.&
          &parini%conf_cartred(i).eq."K".or.parini%conf_cartred(i).eq."k".or.&
          &parini%conf_cartred(i).eq."R".or.parini%conf_cartred(i).eq."r".or.&
          &parini%conf_cartred(i).eq."D".or.parini%conf_cartred(i).eq."d")) stop "Error in conf_cartred"
  if(parini%conf_dim(i).lt.1.or.parini%conf_dim(i).gt.3) stop "Error in conf_dim"
  if(parini%conf_exp(i).lt.1) stop "Error in conf_exp"
  if(parini%conf_prefac(i).lt.0.d0) stop "Error in conf_prefac"
  if(parini%conf_av(i).lt.1.or.parini%conf_av(i).gt.2) stop "Error in parini%conf_av"
          do j=1,parini%nat
            if(parini%conf_list(j,i).lt.1.or.parini%conf_list(j,i).gt.parini%nat) stop "Error in conf_list"
          enddo
enddo
!if(ipi_inet.lt.0 .or. ipi_inet.gt.1) stop "Error in ipi_inet: must be 0 for unix socket, 1 for tcp"
!if(ipi_port.lt.1) stop "Error in ipi_port"
!if(ipi_ecutwf(1).lt.0.d0) stop "Error in ipi_ecutwfc"
!if(ipi_ecutwf(2).lt.0.d0) stop "Error in ipi_ecutwfc"

if(parini%sock_inet.lt.0 .or. parini%sock_inet.gt.1) stop "Error in sock_inet: must be 0 for unix socket, 1 for tcp"
if(parini%sock_port.lt.1) stop "Error in sock_port"
if(parini%sock_ecutwf(1).lt.0.d0) stop "Error in sock_ecutwfc"
if(parini%sock_ecutwf(2).lt.0.d0) stop "Error in sock_ecutwfc"
!SQNM
if(parini%paropt_geopt%beta_lat.le.0d0) stop "Error in sqnm_beta_lat"
if(parini%paropt_geopt%beta_at.le.0.d0) stop "Error in parini%paropt_geopt%beta_at"
if(parini%paropt_geopt%nhist.lt.1)      stop "Error in sqnm_nhist"
if(parini%paropt_geopt%maxrise.le.0.d0) stop "Error in sqnm_maxrise"
if(parini%paropt_geopt%cutoffRatio.le.0.d0) stop "Error in sqnm_cutoffRatio"
if(parini%paropt_geopt%steepthresh.le.0.d0) stop "Error in sqnm_steepthresh"
if(parini%paropt_geopt%trustr.le.0.d0)  stop "Error in sqnm_trustr"


end subroutine

!************************************************************************************

subroutine params_echo(parini)
use defs_basis
use String_Utility 
use mod_fire,   only:dtmin, dtmax
use minpar, only:parmin_bfgs
use fingerprint, only: & 
   fp_method,&!All
   fp_12_nl,&                            !CALYPSO parameters
   fp_13_nl,&                            !Modified CALYPSO parameters
   fp_18_large_vanradius
use mod_parini, only: typ_parini
implicit none
type(typ_parini), intent(in):: parini
integer:: i,j
character(1):: fn1
character(60):: formatting,string
!Echo the variables
write(*,'(a)')             " ################################ Echo params_new.in #############################"
write(*,'(a,i5)')          " # VERBOSITY     ", parini%verb
write(*,'(a)')             " # SYSTEM parameters *************************************************************"
write(*,'(a,i5)')          " # BOUNDARY      ", parini%bc
write(*,'(a,es15.7)')      " # PRESS         ", parini%target_pressure_gpa
write(*,'(a,L3)')          " # VOIDS         ", parini%voids
write(*,'(a,L3)')          " # COREREP       ", parini%core_rep
write(*,'(a)')             " # COMPUTE parameters ************************************************************"
write(*,'(a,a)')           " # CODE          ", trim(parini%potential_potential)
write(*,'(a,L3)')          " # FINDSYM       ", parini%findsym
write(*,'(a,L3)')          " # FINDDOS       ", parini%finddos
write(*,'(a,L3)')          " # USEWFGEO      ", parini%usewf_geopt
write(*,'(a,L3)')          " # USEWFMD       ", parini%usewf_md
write(*,'(a,L3)')          " # USEWFSOFT     ", parini%usewf_soften
write(*,'(a)')             " # Atomic parameters *************************************************************"
write(*,'(a,i5)')          " # NAT           ", parini%nat
write(*,'(a,i5)')          " # NTYPE         ", parini%ntypat_global
write(formatting,'(a,i5,a)') '(a,',parini%ntypat_global,'i4)'
write(*,trim(formatting))  " # ZNUCL         ", int(parini%znucl(:))
write(formatting,'(a,i5,a)') '(a,',parini%nat,'i3)'
write(*,trim(formatting))  " # TYPAT         ", parini%typat_global(:)
write(formatting,'(a,i5,a)') '(a,',parini%ntypat_global,'i4)'
write(*,trim(formatting))  " # AMU           ", int(parini%amu(:))
write(*,'(a)')             " # MD parameters *****************************************************************"
write(*,'(a,i5)')          " # MDALGO        ", parini%md_algo
write(*,'(a,es15.7)')      " # MDPRESSCOMP   ", parini%md_presscomp
write(*,'(a,i5)')          " # MDINT         ", parini%md_integrator
write(*,'(a,i5)')          " # MDNIT         ", parini%nmd_dynamics
write(*,'(a,L3)')          " # AUTO_MDDT     ", parini%auto_dtion_md
write(*,'(a,L3)')          " # MDENCON       ", parini%energy_conservation
write(*,'(a,es15.7)')      " # MDDTINIT      ", parini%dtion_md
write(*,'(a,i5)')          " # MDDTIPM       ", parini%nit_per_min
write(*,'(a,L3)')          " # AUTO_MDMIN    ", parini%auto_mdmin
write(*,'(a,i5)')          " # MDMININIT     ", parini%mdmin
write(*,'(a,i5)')          " # MDMINMIN      ", parini%mdmin_min
write(*,'(a,i5)')          " # MDMINMAX      ", parini%mdmin_max
write(*,'(a,es15.7)')      " # CELLMASS      ", parini%bmass
write(*,'(a)')             " # GEOPT parameters **************************************************************"
write(*,'(a,L3)')          " # GEOEXT        ", parini%geopt_ext
write(*,'(a,i5)')          " # GEONIT        ", parini%paropt_geopt%nit
write(*,'(a,es15.7)')      " # GEOTOLMXF     ", parini%paropt_geopt%fmaxtol
write(*,'(a,es15.7)')      " # STRFACT       ", parini%paropt_geopt%strfact
if(.not.parini%geopt_ext) then
write(*,'(a,a)')           " # GEOALGO       ", trim(parini%paropt_geopt%approach) 
if(trim(parini%paropt_geopt%approach)=="FIRE") then
write(*,'(a,es15.7)')      " # GEOFIREDTINIT ", parini%paropt_geopt%dt_start
write(*,'(a,es15.7)')      " # GEOFIREDTMIN  ", parini%paropt_geopt%dtmin
write(*,'(a,es15.7)')      " # GEOFIREDTMAX  ", parini%paropt_geopt%dtmax
elseif(trim(parini%paropt_geopt%approach)=="MBFGS".or.&
       trim(parini%paropt_geopt%approach)=="RBFGS".or.&
       trim(parini%paropt_geopt%approach)=="SQNM".or. &
       trim(parini%paropt_geopt%approach)=="SD") then
write(*,'(a,es15.7)')      " # GEOHESSAT     ", parini%alphax_at
write(*,'(a,es15.7)')      " # GEOHESSLAT    ", parini%alphax_lat
endif
if(trim(parini%paropt_geopt%approach)=="QBFGS") then
write(*,'(a,i5)'    )      " # GEOQBFGSNDIM  ", parini%qbfgs_bfgs_ndim
write(*,'(a,es15.7)')      " # GEOQBFGSTRI   ", parini%qbfgs_trust_radius_ini
write(*,'(a,es15.7)')      " # GEOQBFGSTRMIN ", parini%qbfgs_trust_radius_min
write(*,'(a,es15.7)')      " # GEOQBFGSTRMAX ", parini%qbfgs_trust_radius_max
write(*,'(a,es15.7)')      " # GEOQBFGSW1    ", parini%qbfgs_w_1
write(*,'(a,es15.7)')      " # GEOQBFGSW2    ", parini%qbfgs_w_2
endif
if(trim(parini%paropt_geopt%approach)=="SQNM") then
write(*,'(a,i5)')          " # GEOSQNMNHIST  ", parini%paropt_geopt%nhist
write(*,'(a,es15.7)')      " # GEOSQNMMAXRISE", parini%paropt_geopt%maxrise
write(*,'(a,es15.7)')      " # GEOSQNMCUTOFF ", parini%paropt_geopt%cutoffRatio
write(*,'(a,es15.7)')      " # GEOSQNMSTEEP  ", parini%paropt_geopt%steepthresh
write(*,'(a,es15.7)')      " # GEOSQNMTRUSTR ", parini%paropt_geopt%trustr
endif
endif 
write(*,'(a)')             " # SOFTEN parameters *************************************************************"
write(*,'(a,L3)')          " # AUTO_SOFT     ", parini%auto_soft
write(*,'(a,L3)')          " # MOLSOFT       ", parini%mol_soften
write(*,'(a,es15.7)')      " # SOFTAT        ", parini%alpha_at
write(*,'(a,es15.7)')      " # SOFTLAT       ", parini%alpha_lat
write(*,'(a,i5)')          " # SOFTNIT       ", parini%nsoften_minhopp
write(*,'(a)')             " # KPOINTS parameters ************************************************************"
write(*,'(a,L3)')          " # AUTO_KPT     ", parini%auto_kpt
if(.not.parini%auto_kpt) then
write(*,'(a,3i5)')         " # KPTMESH      ", parini%ka,parini%kb,parini%kc
else
write(*,'(a,2es15.7)')     " # KPTDEN       ", parini%dkpt1,parini%dkpt2
endif
write(*,'(a)')             " # FINGERPRINT parameters ********************************************************"
write(*,'(a,a)')           " # FPMETHOD      ", trim(parini%fp_method_ch)
if(parini%bc==1.or.parini%bc==3) then
write(*,'(a,es15.7)')      " # FPCUT         ", parini%fp_rcut
endif
if(trim(parini%fp_method_ch)=="OGANOV") then
write(*,'(a,es15.7)')      " # FPDBIN        ", parini%fp_dbin
write(*,'(a,es15.7)')      " # FPSIGMA       ", parini%fp_sigma
elseif(trim(parini%fp_method_ch)=="BCM".or.trim(parini%fp_method_ch)=="ATORB") then
write(*,'(a,i5)')          " # FPNL          ", parini%fp_nl
elseif(trim(parini%fp_method_ch)=="XYZ2SM") then
write(*,'(a,i5)')          " # FPPOWER       ", parini%fp_14_m
write(*,'(a,es15.7)')      " # FPGAUSSFAC1   ", parini%fp_14_w1
write(*,'(a,es15.7)')      " # FPGAUSSFAC2   ", parini%fp_14_w2
elseif(trim(parini%fp_method_ch)=="COGANOV") then
write(*,'(a,es15.7)')      " # FPSIGMA       ", parini%fp_sigma
write(*,'(a,i5)')          " # FPATNMAX      ", parini%fp_at_nmax
elseif(trim(parini%fp_method_ch)=="CAOGANOV") then
write(*,'(a,es15.7)')      " # FPSIGMA       ", parini%fp_sigma
write(*,'(a,i5)')          " # FPATNMAX        ", parini%fp_at_nmax
elseif(trim(parini%fp_method_ch)=="GOM") then
write(*,'(a,i5)')          " # FPNATX        ", parini%fp_17_natx_sphere
write(*,'(a,i5)')          " # FPLSEG        ", parini%fp_17_lseg
write(*,'(a,a)')           " # FPORBITAL     ", parini%fp_17_orbital
write(*,'(a,es15.7)')      " # FPNEXCUT      ", parini%fp_17_nex_cutoff
elseif(trim(parini%fp_method_ch)=="MOLGOM") then
write(*,'(a,a)')           " # FPORBITAL     ", parini%fp_18_orbital
write(*,'(a,i5)')          " # FPNEXCUT      ", parini%fp_18_nex_cutoff
write(*,'(a,i5)')          " # FPPRINCIPLEEV ", parini%fp_18_principleev
write(*,'(a,i5)')          " # FPMOLECULES   ", parini%fp_18_molecules
write(*,'(a,i5)')          " # FPEXPA        ", parini%fp_18_expaparameter
write(*,'(a,i5)')          " # FPMOLSPHERE   ", parini%fp_18_molecules_sphere
write(*,'(a,es15.7)')      " # FPWIDTHCUT    ", parini%fp_18_width_cutoff
write(*,'(a,es15.7)')      " # FPWIDTHOVER   ", parini%fp_18_width_overlap
endif
write(*,'(a)')             " # CONFINEMENT parameters ********************************************************"
write(*,'(a,L3)')          " # CONFINEMENT   ",parini%use_confine
if(parini%use_confine) then
write(*,'(a,i5)')          " # CONFNCONF     ",parini%nconfine
write(formatting,'(a,i5,a)') '(a,',parini%nconfine,'A4)'
write(*,trim(formatting))  " # CONFCARTRED   ",parini%conf_cartred
write(formatting,'(a,i5,a)') '(a,',parini%nconfine,'i4)'
write(*,trim(formatting))  " # CONFDIM       ",parini%conf_dim
write(formatting,'(a,i5,a)') '(a,',parini%nconfine,'i4)'
write(*,trim(formatting))  " # CONFEXP       ",parini%conf_exp
write(formatting,'(a,i5,a)') '(a,',parini%nconfine,'es9.2)'
write(*,trim(formatting))  " # CONFPREFAC    ",parini%conf_prefac
write(formatting,'(a,i5,a)') '(a,',parini%nconfine,'es9.2)'
write(*,trim(formatting))  " # CONFCUT       ",parini%conf_cut
write(formatting,'(a,i5,a)') '(a,',parini%nconfine,'i4)'
write(*,trim(formatting))  " # CONFAV        ",parini%conf_av
write(formatting,'(a,i5,a)') '(a,',parini%nconfine,'es9.2)'
write(*,trim(formatting))  " # CONFEQ        ",parini%conf_eq
write(formatting,'(a,i5,a)') '(a,',parini%nconfine,'i4)'
write(*,trim(formatting))  " # CONFNAT       ",parini%conf_nat
!Go through all number of confinements
   do i=1,parini%nconfine
      if(i.lt.10) then
        write(fn1,'(i1.1)') i
        write(string,'(a)') "CONFLIST"//fn1
        write(formatting,'(a,i5,a)') '(a,',parini%conf_nat(i),'i4)'
        write(*,trim(formatting))  " # "//trim(string)//"     ",parini%conf_list(1:parini%conf_nat(i),i)
      else
        stop "Loop not implemented for nconfine greater than 9"
      endif
   enddo
endif
if(StrLowCase(trim(adjustl(parini%potential_potential)))=="ipi") then
write(*,'(a)')             " # IPI parameters ****************************************************************"
write(formatting,'(a)') '(a,i4)'
write(*,trim(formatting))  " # IPIINET       ",parini%sock_inet
write(formatting,'(a)') '(a,i8)'
write(*,trim(formatting))  " # IPIPORT       ",parini%sock_port
write(formatting,'(a)') '(a,a)'
write(*,trim(formatting))  " # IPIHOST       ",trim(adjustl(parini%sock_host))
write(formatting,'(a)') '(a,2f10.4)'
write(*,trim(formatting))  " # IPIECUTWF     ",parini%sock_ecutwf(:)
endif
if(StrLowCase(trim(adjustl(parini%potential_potential)))=="msock") then
write(*,'(a)')             " # MSOCK parameters **************************************************************"
write(formatting,'(a)') '(a,i4)'
write(*,trim(formatting))  " # SOCKINET      ",parini%sock_inet
write(formatting,'(a)') '(a,i8)'
write(*,trim(formatting))  " # SOCKPORT      ",parini%sock_port
write(formatting,'(a)') '(a,a)'
write(*,trim(formatting))  " # SOCKHOST      ",trim(adjustl(parini%sock_host))
write(formatting,'(a)') '(a,2f10.4)'
write(*,trim(formatting))  " # SOCKECUTWF    ",parini%sock_ecutwf(:)
endif
write(*,'(a)')             " ############################ END Echo params_new.in #############################"
end subroutine

!************************************************************************************
subroutine fp_assign(parini)
use mod_parini, only: typ_parini
use fingerprint 
implicit none
type(typ_parini), intent(inout):: parini
if(parini%bc==2) then !Molecular systems
  fp_method=21
  parini%fp_method_ch="GAUSS"
elseif(parini%bc==1) then !Crystals
  if(trim(parini%fp_method_ch)=="OGANOV") then
    fp_method=11
  elseif(trim(parini%fp_method_ch)=="BCM") then
    fp_method=12
  elseif(trim(parini%fp_method_ch)=="ATORB") then
    fp_method=13
  elseif(trim(parini%fp_method_ch)=="XYZ2SM") then
    fp_method=14
  elseif(trim(parini%fp_method_ch)=="COGANOV") then
    fp_method=15
  elseif(trim(parini%fp_method_ch)=="CAOGANOV") then
    fp_method=16
  elseif(trim(parini%fp_method_ch)=="GOM") then
    fp_method=17
    if(trim(parini%fp_17_orbital)=="S")then
    parini%fp_17_lseg=1
    elseif(trim(parini%fp_17_orbital)=="SP")then
    parini%fp_17_lseg=4 
    endif
    parini%fp_17_width_cutoff=parini%fp_rcut/sqrt(2.d0*parini%fp_17_nex_cutoff)
  elseif(trim(parini%fp_method_ch)=="MOLGOM") then
    fp_method=18
    if(trim(parini%fp_18_orbital)=="S")then
      parini%fp_18_lseg=1
    elseif(trim(parini%fp_18_orbital)=="SP")then
      parini%fp_18_lseg=4 
    endif
!    parini%fp_18_width_cutoff=fp_rcut/sqrt(2.d0*parini%fp_18_nex_cutoff)
  endif
endif    
end subroutine

