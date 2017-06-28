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
use global, only: target_pressure_habohr,target_pressure_gpa,nat,ntypat,znucl,amu,amutmp,typat,ntime_md,char_type,&
                &nsoften,alpha_at,alpha_lat,ntime_geopt,bmass,mdmin,dtion_fire,dtion_md,tolmxf,strfact,dtion_fire_min,&
                &dtion_fire_max,ka,kb,kc,dkpt1,dkpt2,usewf_geopt,usewf_soften,usewf_md,geopt_method,alphax_at,&
                &alphax_lat,findsym,finddos,auto_soft,mdmin_max,mdmin_min,auto_mdmin,md_algo,md_integrator,auto_dtion_md,&
                &nit_per_min,fixat,fixlat,rcov,mol_soften,fragarr,code,auto_kpt,bc,geopt_ext,energy_conservation,use_confine,&
                &voids,core_rep,md_presscomp
use sqnm,   only: sqnm_beta_lat,sqnm_beta_at,sqnm_nhist,sqnm_maxrise,sqnm_cutoffRatio,sqnm_steepthresh,sqnm_trustr
use qbfgs,  only: qbfgs_bfgs_ndim,qbfgs_trust_radius_max,qbfgs_trust_radius_min,qbfgs_trust_radius_ini,qbfgs_w_1,qbfgs_w_2
use steepest_descent, only: sd_beta_lat,sd_beta_at
use modsocket, only:sock_inet,sock_port,sock_host,sock_ecutwf
use fingerprint, only: & 
   fp_rcut,fp_method,fp_method_ch,fp_nl,&!All
   fp_sigma,fp_dbin,&              !Oganov parameters
   fp_12_nl,&                            !CALYPSO parameters
   fp_13_nl,&                            !Modified CALYPSO parameters
   fp_14_m,fp_14_w1,fp_14_w2,&          !xyz2sm parameters
   fp_at_nmax,&
   fp_17_width_cutoff,fp_17_nex_cutoff,fp_17_natx_sphere,fp_17_lseg,fp_17_orbital,&
   fp_18_orbital,fp_18_principleev,fp_18_lseg,fp_18_molecules,&
   fp_18_expaparameter,fp_18_nex_cutoff,fp_18_molecules_sphere,fp_18_width_cutoff,&
   fp_18_width_overlap,fp_18_large_vanradius


use confinement

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
nconfine=1
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
   call parsescalar_int("NAT",3,all_line(1:n),n,nat,found)
   if(found) nat_found=.true.
   if(found) cycle
   endif
!NTYPE
   if(.not.ntypat_found) then
   call parsescalar_int("NTYPE",5,all_line(1:n),n,ntypat,found)
   if(found) ntypat_found=.true.
   if(found) cycle
   endif
!ZNUCL
   if(.not.znucl_found) then
   call parsescalar_real("ZNUCL",5,all_line(1:n),n,tmp_val,found)
   if(found) znucl_found=.true.
   if(found) exit
   endif
!CONFNCONF
   call parsescalar_int("CONFNCONF",9,all_line(1:n),n,nconfine,found)
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
     call ascii_getsystem(filename)
     read_poscur=.true.
   endif
   if(.not.read_poscur) then
!Try to read the vasp input file for the system 
   file_exists=.false.
   filename="poscur.vasp"
   INQUIRE(FILE=trim(filename), EXIST=file_exists)
   if(file_exists) then
     write(*,*) "Trying to read them from ",filename !Will also read ZNUCL here
     call poscar_getsystem(filename)
     read_poscur=.true.
   endif
   endif
endif

!Allocate the arrays
 if(.not.allocated(znucl))       then;   allocate(znucl(ntypat))                       ; znucl=0                 ; endif
 if(.not.allocated(char_type))   then;   allocate(char_type(ntypat))                   ; char_type="  "          ; endif
 if(.not.allocated(amu))         then;   allocate(amu(ntypat))                         ; amu=0                   ; endif
 if(.not.allocated(amutmp))      then;   allocate(amutmp(ntypat))                      ; amutmp=0                ; endif
 if(.not.allocated(rcov))        then;   allocate(rcov(ntypat))                        ; rcov=0                  ; endif
 if(.not.allocated(typat))       then;   allocate(typat(nat))                          ; typat=0                 ; endif
 if(.not.allocated(fixat))       then;   allocate(fixat(nat))                          ; fixat=.false.           ; endif
 if(.not.allocated(fragarr))     then;   allocate(fragarr(nat))                        ; fragarr=0               ; endif
 if(.not.allocated(conf_dim))    then;   allocate(conf_dim     (nconfine))             ; conf_dim=0              ; endif
 if(.not.allocated(conf_av))     then;   allocate(conf_av      (nconfine))             ; conf_av=0               ; endif
 if(.not.allocated(conf_exp))    then;   allocate(conf_exp     (nconfine))             ; conf_exp=0              ; endif
 if(.not.allocated(conf_prefac)) then;   allocate(conf_prefac  (nconfine))             ; conf_prefac=0           ; endif
 if(.not.allocated(conf_cut))    then;   allocate(conf_cut     (nconfine))             ; conf_cut=0              ; endif
 if(.not.allocated(conf_eq))     then;   allocate(conf_eq      (nconfine))             ; conf_eq=0               ; endif
 if(.not.allocated(conf_list))   then;   allocate(conf_list    (nat,nconfine))         ; conf_list=0             ; endif
 if(.not.allocated(conf_nat))    then;   allocate(conf_nat     (nconfine))             ; conf_nat=0              ; endif
 if(.not.allocated(conf_cartred))then;   allocate(conf_cartred (nconfine))             ; conf_cartred="C"        ; endif


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
      call parsearray_real("ZNUCL",5,all_line(1:n),n,znucl(1:ntypat),ntypat,found)
      if(found) exit
   enddo
   98 continue
   close(12)
   if(.not.found) then
     stop "ZNUCL not found in params_new.in/poscur.ascii/poscur.vasp"
   endif
endif

!Get the correct atomic masses and atomic character
 do itype=1,ntypat
   call atmdata(amu(itype),rcov(itype),char_type(itype),znucl(itype))
 enddo

!Read the other variables
open(unit=12,file="params_new.in")
  do while(.true.)
   read(12,'(a450)',end=97)all_line
   n = len_trim(all_line)
!Press
   call parsescalar_real("PRESS",5,all_line(1:n),n,target_pressure_gpa,found)
   if(found) target_pressure_habohr=target_pressure_gpa/HaBohr3_GPA
   if(found) cycle
!Amu
   call parsearray_real("AMU",3,all_line(1:n),n,amu(1:ntypat),ntypat,found)
   if(found) cycle
!TYPAT
   call parsearray_int("TYPAT",5,all_line(1:n),n,typat(1:nat),nat,found)
   if(found) cycle
!MDNIT
   call parsescalar_int("MDNIT",5,all_line(1:n),n,ntime_md,found)
   if(found) cycle
!MDALGO
   call parsescalar_int("MDALGO",6,all_line(1:n),n,md_algo,found)
   if(found) cycle
!MDINTEGRATOR
   call parsescalar_int("MDINT",5,all_line(1:n),n,md_integrator,found)
   if(found) cycle
!MDPRESSCOMP
   call parsescalar_real("MDPRESSCOMP",11,all_line(1:n),n,md_presscomp,found)
   if(found) cycle
!GEONIT
   call parsescalar_int("GEONIT",6,all_line(1:n),n,ntime_geopt,found)
   if(found) cycle
!CELLMASS
   call parsescalar_real("CELLMASS",8,all_line(1:n),n,bmass,found)
   if(found) cycle
!VOIDS
   call parse_logical("VOIDS",5,all_line(1:n),n,voids,found)
   if(found) cycle
!COREREP
   call parse_logical("COREREP",7,all_line(1:n),n,core_rep,found)
   if(found) cycle
!Block mdmin****************
!AUTO_MDMIN
   call parse_logical("AUTO_MDMIN",10,all_line(1:n),n,auto_mdmin,found)
   if(found) cycle
!MDMININIT
   call parsescalar_int("MDMININIT",9,all_line(1:n),n,mdmin_in,found)
   if(found) cycle
!MDMINMIN
   call parsescalar_int("MDMINMIN",8,all_line(1:n),n,mdmin_min,found)
   if(found) cycle
!MDMINMAX
   call parsescalar_int("MDMINMAX",8,all_line(1:n),n,mdmin_max,found)
   if(found) cycle
!Block mdmin****************
!Block soften****************
!AUTO_SOFT
   call parse_logical("AUTO_SOFT",9,all_line(1:n),n,auto_soft,found)
   if(found) cycle
!SOFTLAT
   call parsescalar_real("SOFTLAT",7,all_line(1:n),n,alpha_lat_in,found)
   if(found) cycle
!SOFTAT
   call parsescalar_real("SOFTAT",6,all_line(1:n),n,alpha_at_in,found)
   if(found) cycle
!MOLSOFT
   call parse_logical("MOLSOFT",7,all_line(1:n),n,mol_soften,found)
   if(found) cycle
!SOFTNIT
   call parsescalar_int("SOFTNIT",7,all_line(1:n),n,nsoften,found)
   if(found) cycle
!Block soften****************
!Block MD timestep***********
!AUTO_MDDT
   call parse_logical("AUTO_MDDT",9,all_line(1:n),n,auto_dtion_md,found)
   if(found) cycle
!MDDTINIT
   call parsescalar_real("MDDTINIT",8,all_line(1:n),n,dtion_md_in,found)
   if(found) cycle
!MDDTIPM
   call parsescalar_int("MDDTIPM",7,all_line(1:n),n,nit_per_min,found)
   if(found) cycle
!MDENCON
   call parse_logical("MDENCON",7,all_line(1:n),n,energy_conservation,found)
   if(found) cycle
!Block MD timestep***********
!Block GEOPT*****************
!GEOALGO
   call parsescalar_string("GEOALGO",7,all_line(1:n),n,geopt_method,5,found)
   if(found) geopt_method=StrUpCase(geopt_method)
   if(found) cycle
!GEOFIREDTINIT
   call parsescalar_real("GEOFIREDTINIT",13,all_line(1:n),n,dtion_fire,found)
   if(found) cycle
!GEOFIREDTMIN
   call parsescalar_real("GEOFIREDTMIN",12,all_line(1:n),n,dtion_fire_min,found)
   if(found) cycle
!GEOFIREDTMAX
   call parsescalar_real("GEOFIREDTMAX",12,all_line(1:n),n,dtion_fire_max,found)
   if(found) cycle
!GEOHESSLAT
   call parsescalar_real("GEOHESSLAT",10,all_line(1:n),n,alphax_lat,found)
   if(found) cycle
!GEOHESSAT
   call parsescalar_real("GEOHESSAT",9,all_line(1:n),n,alphax_at,found)
   if(found) cycle
!GEOTOLMXF
   call parsescalar_real("GEOTOLMXF",9,all_line(1:n),n,tolmxf,found)
   if(found) cycle
!GEOSTRFACT
   call parsescalar_real("STRFACT",7,all_line(1:n),n,strfact,found)
   if(found) cycle
!GEOEXT
   call parse_logical("GEOEXT",6,all_line(1:n),n,geopt_ext,found)
   if(found) cycle
!GEOSQNMNHIS
   call parsescalar_int ("GEOSQNMNHIST",12,all_line(1:n),n,sqnm_nhist,found)
   if(found) cycle
!GEOSQNMMAXRISE
   call parsescalar_real("GEOSQNMMAXRISE",14,all_line(1:n),n,sqnm_maxrise,found)
   if(found) cycle
!GEOSQNMCUTOFF
   call parsescalar_real("GEOSQNMCUTOFF",13,all_line(1:n),n,sqnm_cutoffRatio,found)
   if(found) cycle
!GEOSQNMSTEEP
   call parsescalar_real("GEOSQNMSTEEP",12,all_line(1:n),n,sqnm_steepthresh,found)
   if(found) cycle
!GEOSQNMTRUSTR
   call parsescalar_real("GEOSQNMTRUSTR",13,all_line(1:n),n,sqnm_trustr,found)
   if(found) cycle
!GEOQBFGSNDIM
   call parsescalar_int ("GEOQBFGSNDIM",12,all_line(1:n),n,qbfgs_bfgs_ndim,found)
   if(found) cycle
!GEOQBFGSTRTI
   call parsescalar_real("GEOQBFGSTRI",11,all_line(1:n),n,qbfgs_trust_radius_ini,found)
   if(found) cycle
!GEOQBFGSTRTMIN
   call parsescalar_real("GEOQBFGSTRMIN",13,all_line(1:n),n,qbfgs_trust_radius_min,found)
   if(found) cycle
!GEOQBFGSTRTMAX
   call parsescalar_real("GEOQBFGSTRMAX",13,all_line(1:n),n,qbfgs_trust_radius_max,found)
   if(found) cycle
!GEOQBFGSW1
   call parsescalar_real("GEOQBFGSW1",10,all_line(1:n),n,qbfgs_w_1,found)
   if(found) cycle
!GEOQBFGSW2
   call parsescalar_real("GEOQBFGSW1",10,all_line(1:n),n,qbfgs_w_2,found)
   if(found) cycle
!Block GEOPT*****************
!USEWFGEO
   call parse_logical("USEWFGEO",8,all_line(1:n),n,usewf_geopt,found)
   if(found) cycle
!USEWFSOFT
   call parse_logical("USEWFSOFT",9,all_line(1:n),n,usewf_soften,found)
   if(found) cycle
!USEWFMD
   call parse_logical("USEWFMD",7,all_line(1:n),n,usewf_md,found)
   if(found) cycle
!FINDSYM
   call parse_logical("FINDSYM",7,all_line(1:n),n,findsym,found)
   if(found) cycle
!FINDDOS
   call parse_logical("FINDDOS",7,all_line(1:n),n,finddos,found)
   if(found) cycle
!Block KPT****************
!AUTO_KPT
   call parse_logical("AUTO_KPT",8,all_line(1:n),n,auto_kpt,found)
   if(found) cycle
!KPTMESH
   call parsearray_int("KPTMESH",7,all_line(1:n),n,kpt_abc(1:3),3,found)
   if(found) then
     ka=kpt_abc(1)
     kb=kpt_abc(2)
     kc=kpt_abc(3)
   endif
   if(found) cycle
!KPTDEN
   call parsearray_real("KPTDEN",6,all_line(1:n),n,dkpt_12(1:2),2,found)
   if(found) then
     dkpt1=dkpt_12(1)
     dkpt2=dkpt_12(2)
   endif
   if(found) cycle
!Block KPT****************
!Block FINGERPRINT****************
!FPMETHOD
   call parsescalar_string("FPMETHOD",8,all_line(1:n),n,fp_method_ch,20,found)
   if(found) fp_method_ch=StrUpCase(fp_method_ch)
   if(found) cycle
!FPCUT
   call parsescalar_real("FPCUT",5,all_line(1:n),n,fp_rcut,found)
   if(found) cycle
!FPDBIN
   call parsescalar_real("FPDBIN",6,all_line(1:n),n,fp_dbin,found)
   if(found) cycle
!FPSIGMA
   call parsescalar_real("FPSIGMA",7,all_line(1:n),n,fp_sigma,found)
   if(found) cycle
!FPNL
   call parsescalar_int("FPNL",4,all_line(1:n),n,fp_nl,found)
   if(found) cycle
!FPPOWER
   call parsescalar_int("FPPOWER",7,all_line(1:n),n,fp_14_m,found)
   if(found) cycle
!FPGAUSSFAC1
   call parsescalar_real("FPGAUSSFAC1",11,all_line(1:n),n,fp_14_w1,found)
   if(found) cycle
!FPGAUSSFAC2
   call parsescalar_real("FPGAUSSFAC2",11,all_line(1:n),n,fp_14_w2,found)
   if(found) cycle
!FPNMAX
   call parsescalar_int("FPATNMAX",8,all_line(1:n),n,fp_at_nmax,found)
   if(found) cycle

!FPWIDTHCUT
!   call parsescalar_real("FPWIDTHCUT",10,all_line(1:n),n,fp_17_width_cutoff,found)
!FPNATX
   call parsescalar_int("FPNATX",6,all_line(1:n),n,fp_17_natx_sphere,found)
   if(found) cycle
!FPORBITAL
   call parsescalar_string("FPORBITAL",9,all_line(1:n),n,fp_17_orbital,2,found)
   fp_18_orbital=fp_17_orbital
   if(found) cycle
!FPNEXCUT
   call parsescalar_real("FPNEXCUT",8,all_line(1:n),n,fp_17_nex_cutoff,found)
   fp_18_nex_cutoff=int(fp_17_nex_cutoff)
   if(found) cycle
!FPPRINCIPLEEV
   call parsescalar_int("FPPRINCIPLEEV",13,all_line(1:n),n,fp_18_principleev,found)
   if(found) cycle
!FPMOLECULES
   call parsescalar_int("FPMOLECULES",11,all_line(1:n),n,fp_18_molecules,found)
   if(found) cycle
!FPEXPA
   call parsescalar_int("FPEXPA",6,all_line(1:n),n,fp_18_expaparameter,found)
   if(found) cycle
!FPMOLSPHERE
   call parsescalar_int("FPMOLSPHERE",11,all_line(1:n),n,fp_18_molecules_sphere,found)
   if(found) cycle
!FPWIDTHCUT
   call parsescalar_real("FPWIDTHCUT",10,all_line(1:n),n,fp_18_width_cutoff,found)
   if(found) cycle
!FPWIDTHOVER
   call parsescalar_real("FPWIDTHOVER",11,all_line(1:n),n,fp_18_width_overlap,found)
   if(found) cycle
!Block FINGERPRINT****************


!Block LAYER-CONFINEMENT****************
!CONFINEMENT
   call parse_logical("CONFINEMENT",11,all_line(1:n),n,use_confine,found)
   if(found) cycle
!CONFCARTRED
   call parsearray_string("CONFCARTRED",11,all_line(1:n),n,conf_cartred,1,nconfine,found)
   if(found) cycle
!CONFDIM
   call parsearray_int("CONFDIM",7,all_line(1:n),n,conf_dim(1:nconfine),nconfine,found)
   if(found) cycle
!CONFEXP
   call parsearray_int("CONFEXP",7,all_line(1:n),n,conf_exp(1:nconfine),nconfine,found)
   if(found) cycle
!CONFPREFAC
   call parsearray_real("CONFPREFAC",10,all_line(1:n),n,conf_prefac(1:nconfine),nconfine,found)
   if(found) cycle
!CONFCUT
   call parsearray_real("CONFCUT",7,all_line(1:n),n,conf_cut(1:nconfine),nconfine,found)
   if(found) cycle
!CONFAV
   call parsearray_int("CONFAV",6,all_line(1:n),n,conf_av(1:nconfine),nconfine,found)
   if(found) cycle
!CONFEQ
   call parsearray_real("CONFEQ",6,all_line(1:n),n,conf_eq(1:nconfine),nconfine,found)
   if(found) cycle
!CONFNAT
   call parsearray_int("CONFNAT",7,all_line(1:n),n,conf_nat(1:nconfine),nconfine,found)
   if(found) cycle
!CONFLIST#
!Go through all number of confinements
   do i=1,nconfine
      if(i.lt.10) then
        write(fn1,'(i1.1)') i
        write(find_string,'(a,a)') "CONFLIST"//fn1 
        call exist_string(trim(find_string),9,all_line(1:n),n,found)
        if(found.and.conf_nat(i)==nat) then
          do j=1,nat
            conf_list(j,i)=j
          enddo
        else
          call parsearray_int(trim(find_string),9,all_line(1:n),n,conf_list(1:conf_nat(i),i),conf_nat(i),found)
        endif
      else
        stop "Loop not implemented for nconfine greater than 9"
      endif
   enddo
!Block LAYER-CONFINEMENT****************

!Block IPI_SOCKET****************
!IPIINET
   call parsescalar_int("IPIINET",7,all_line(1:n),n,sock_inet,found)
   if(found) cycle
!IPIPORT
   call parsescalar_int("IPIPORT",7,all_line(1:n),n,sock_port,found)
   if(found) cycle
!IPIHOST
   call parsescalar_string("IPIHOST",7,all_line(1:n),n,sock_host,1024,found)
   if(found) cycle
!IPIECUTWF
   call parsearray_real("IPIECUTWF",9,all_line(1:n),n,sock_ecutwf,2,found)
   if(found) cycle
!Block IPI_SOCKET****************

!Block MSOCK****************
!MSOCKINET
   call parsescalar_int("SOCKINET",8,all_line(1:n),n,sock_inet,found)
   if(found) cycle
!MSOCKPORT
   call parsescalar_int("SOCKPORT",8,all_line(1:n),n,sock_port,found)
   if(found) cycle
!MSOCKHOST
   call parsescalar_string("SOCKHOST",8,all_line(1:n),n,sock_host,1024,found)
   if(found) cycle
!MSOCKECUTWF
   call parsearray_real("SOCKECUTWF",10,all_line(1:n),n,sock_ecutwf,2,found)
   if(found) cycle

!Block MSOCK****************

!VERBOSE
   call parsescalar_int("VERBOSE",7,all_line(1:n),n,parini%verb,found)
   if(found) cycle
!BOUNDARY
   call parsescalar_int("BOUNDARY",8,all_line(1:n),n,bc,found)
   if(found) cycle
!CODE
   call parsescalar_string("CODE",4,all_line(1:n),n,code,20,found)
   if(found) code = StrLowCase( code )
   if(found) cycle
  enddo
97 continue
close(12)


!Post Processing
!MDMIN
  if(calls==0) then
      if(.not.auto_mdmin) then
           mdmin=mdmin_in
      else
           mdmin=max(mdmin_in,mdmin_min)
      endif
  endif
!SOFTEN
  if(calls==0.or..not.auto_soft) alpha_lat=alpha_lat_in
  if(calls==0.or..not.auto_soft) alpha_at=alpha_at_in
!MDTIMESTE=P
  if(calls==0.or..not.auto_dtion_md) dtion_md=dtion_md_in
!KPT
  if(auto_kpt) then
    ka=0;kb=0;kc=0
  else
    dkpt1=0.d0
    dkpt2=0.d0
  endif

!Initiallize confinement
if(use_confine) call  init_confinement_parser()

!Initiallize LJ parameter if required
if(trim(code)=="blj".and.calls==0) call blj_init_parameter()

!Initiallize LJ parameter if required
if(trim(code)=="mlj") call mlj_init_parameter()

!Initiallize TB-LJ parameter if required
if(trim(code)=="lenosky_tb_lj".and.calls==0) then
  call check_lenosky_tb_lj()
  n_lj=0
  do iat=1,nat
     if(znucl(typat(iat)).gt.200) n_lj=n_lj+1
  enddo
  call lenosky_tb_lj_init_parameter()
endif

!Initiallize voids
if(voids.and.calls==0) then
  call check_voids() 
  call voids_init_parameter()
endif

!Initiallize tersoff
if(trim(code)=="tersoff".and.calls==0) then
  call init_tersoff()
endif

!Initiallize edip
if(trim(code)=="edip".and.calls==0) then
  call init_edip()
endif

!Initiallize ipi
if(trim(code)=="ipi".and.calls==0) then
  call init_ipi()
endif

!Initiallize msock
if(trim(code)=="msock".and.calls==0) then
  call init_msock()
endif



!Assign correct parameters for fingerprint, not allocating any arrays!
call fp_assign()

!Increase calls to the routine
if(calls==0) call params_echo(parini)
if(calls==0) call params_check(parini)
calls=calls+1

!Check range of parameters
call params_check(parini)

!Copy parameters of fire to the fire module
    dtmin=dtion_fire_min
    dtmax=dtion_fire_max
!Copy parameters of bfgs to the bfgs module
    parmin_bfgs%betax=alphax_at
    parmin_bfgs%betax_lat=alphax_lat
!Copy parameters to sqnm module
    sqnm_beta_lat=alphax_lat
    sqnm_beta_at=alphax_at
!Copy parameters to sd module
    sd_beta_lat=alphax_lat
    sd_beta_at=alphax_at
end subroutine

!************************************************************************************
subroutine params_defaults(parini,mdmin_in,dtion_md_in,alpha_lat_in,alpha_at_in,read_poscur)
use defs_basis
use mod_fire,   only:dtmin, dtmax
use minpar, only:parmin_bfgs
use global, only: target_pressure_habohr,target_pressure_gpa,nat,ntypat,znucl,amu,amutmp,typat,ntime_md,char_type,&
                &nsoften,alpha_at,alpha_lat,ntime_geopt,bmass,mdmin,dtion_fire,dtion_md,tolmxf,strfact,dtion_fire_min,&
                &dtion_fire_max,ka,kb,kc,dkpt1,dkpt2,usewf_geopt,usewf_soften,usewf_md,geopt_method,alphax_at,&
                &alphax_lat,findsym,finddos,auto_soft,mdmin_max,mdmin_min,auto_mdmin,md_algo,md_integrator,auto_dtion_md,&
                &nit_per_min,fixat,fixlat,rcov,mol_soften,fragarr,code,auto_kpt,bc,geopt_ext,energy_conservation,use_confine,&
                &voids,core_rep,md_presscomp
use sqnm,   only: sqnm_beta_lat,sqnm_beta_at,sqnm_nhist,sqnm_maxrise,sqnm_cutoffRatio,sqnm_steepthresh,sqnm_trustr
use qbfgs,  only: qbfgs_bfgs_ndim,qbfgs_trust_radius_max,qbfgs_trust_radius_min,qbfgs_trust_radius_ini,qbfgs_w_1,qbfgs_w_2
use modsocket, only:sock_inet,sock_port,sock_host,sock_ecutwf
use fingerprint, only: & 
   fp_rcut,fp_method,fp_method_ch,fp_nl,&!All
   fp_sigma,fp_dbin,&              !Oganov parameters
   fp_12_nl,&                            !CALYPSO parameters
   fp_13_nl,&                            !Modified CALYPSO parameters
   fp_14_m,fp_14_w1,fp_14_w2,&          !xyz2sm parameters
   fp_at_nmax,&
   fp_17_width_cutoff,fp_17_nex_cutoff,fp_17_natx_sphere,fp_17_lseg,fp_17_orbital,&
   fp_18_orbital,fp_18_principleev,fp_18_lseg,fp_18_molecules,&
   fp_18_expaparameter,fp_18_nex_cutoff,fp_18_molecules_sphere,fp_18_width_cutoff,&
   fp_18_width_overlap,fp_18_large_vanradius
   
use confinement
use mod_parini, only: typ_parini
implicit none
type(typ_parini), intent(inout):: parini
integer:: mdmin_in,itype,i,j
real(8):: dtion_md_in,alpha_lat_in,alpha_at_in
logical:: read_poscur
!These are the default variables
target_pressure_gpa=0.d0
target_pressure_habohr=0.d0
!Get the correct atomic masses and atomic character
 do itype=1,ntypat
   call atmdata(amu(itype),rcov(itype),char_type(itype),znucl(itype))
 enddo
voids=.false.
core_rep=.false.
if(.not.read_poscur) typat(1:nat)=1
ntime_md=300
md_algo=1
md_integrator=3
md_presscomp=-0.d0
ntime_geopt=300
bmass=1.d0
auto_mdmin=.false.
mdmin_in=1
mdmin_min=2
mdmin_max=2
energy_conservation=.false.
auto_soft=.false.
alpha_lat_in=1.d0
alpha_at_in=1.d0
mol_soften=.false.
nsoften=10
auto_dtion_md=.false.
dtion_md_in=20.d0
nit_per_min=25.d0
geopt_method="FIRE"
dtion_fire=10.d0
dtion_fire_min=1.d0
dtion_fire_max=80.d0
alphax_lat=1.d0
alphax_at=1.d0
tolmxf=2.d-4
strfact=100.d0
usewf_geopt=.false.
usewf_soften=.false.
usewf_md=.false.
findsym=.false.
finddos=.false.
auto_kpt=.true.
ka=1;kb=1;kc=1
dkpt1=0.04d0
dkpt2=0.06d0
bc=1
parini%verb=3
code="vasp"
!Define if the external optimizer should be used. Only available for:
geopt_ext=.false.

fp_rcut=15.d0
fp_method=11
fp_method_ch="OGANOV"
fp_nl=6
fp_sigma=0.02d0
fp_dbin= 0.05d0
fp_12_nl=6
fp_13_nl=6
fp_14_m=3
fp_14_w1=1.d0
fp_14_w2=1.5d0
fp_at_nmax=10000
fp_17_nex_cutoff=3
fp_17_width_cutoff=fp_rcut/sqrt(2.d0*fp_17_nex_cutoff)
fp_17_orbital='S'
fp_17_lseg=1
fp_17_natx_sphere=75

fp_18_orbital='S'
fp_18_principleev = 6
fp_18_lseg=1
fp_18_molecules=1
fp_18_expaparameter = 4
fp_18_nex_cutoff = 3
fp_18_molecules_sphere = 50
fp_18_width_cutoff = 1.d0
fp_18_width_overlap = 1.d0
fp_18_large_vanradius = 1.7d0/0.52917720859d0

!SQNM
sqnm_beta_lat=1.d0
sqnm_beta_at=1.d0
sqnm_nhist=10
sqnm_maxrise=1.d-6
sqnm_cutoffRatio=1.d-4
sqnm_steepthresh=1.d0
sqnm_trustr=0.1d0

!QBFGS
qbfgs_bfgs_ndim=1
qbfgs_trust_radius_max=0.5d0
qbfgs_trust_radius_min=1.d-3
qbfgs_trust_radius_ini=0.5D0
qbfgs_w_1=0.01D0
qbfgs_w_2=0.5D0

use_confine=.false.
conf_cartred="C"
conf_dim=1
conf_exp=4
conf_prefac=1.d-2
conf_cut=1.d0
conf_av=2
conf_eq=0
conf_nat=nat
   do i=1,nconfine
          do j=1,nat
            conf_list(j,i)=j
          enddo
   enddo
!Block LAYER-CONFINEMENT****************
!!Block IPI_SOCKET****************
!      ipi_inet=0 !0 for unix socket, 1 for tcp
!      ipi_port=3141
!      ipi_host="mh-driver"
!      ipi_ecutwf=1.d0
!Block IPI_SOCKET****************
      sock_inet=0 !0 for unix socket, 1 for tcp
      sock_port=3141
      sock_host="mh-driver"
      sock_ecutwf=1.d0
end subroutine

!************************************************************************************

subroutine params_check(parini)
use defs_basis
use mod_fire,   only:dtmin, dtmax
use minpar, only:parmin_bfgs
use global, only: target_pressure_habohr,target_pressure_gpa,nat,ntypat,znucl,amu,amutmp,typat,ntime_md,char_type,&
                &nsoften,alpha_at,alpha_lat,ntime_geopt,bmass,mdmin,dtion_fire,dtion_md,tolmxf,strfact,dtion_fire_min,&
                &dtion_fire_max,ka,kb,kc,dkpt1,dkpt2,usewf_geopt,usewf_soften,usewf_md,geopt_method,alphax_at,&
                &alphax_lat,findsym,finddos,auto_soft,mdmin_max,mdmin_min,auto_mdmin,md_algo,md_integrator,auto_dtion_md,&
                &nit_per_min,fixat,fixlat,rcov,mol_soften,fragarr,code,auto_kpt,bc,voids,core_rep,md_presscomp
use sqnm,   only: sqnm_beta_lat,sqnm_beta_at,sqnm_nhist,sqnm_maxrise,sqnm_cutoffRatio,sqnm_steepthresh,sqnm_trustr
use modsocket, only:sock_inet,sock_port,sock_host,sock_ecutwf
use fingerprint, only: & 
   fp_rcut,fp_method,fp_method_ch,fp_nl,&!All
   fp_sigma,fp_dbin,&              !Oganov parameters
   fp_12_nl,&                            !CALYPSO parameters
   fp_13_nl,&                            !Modified CALYPSO parameters
   fp_14_m,fp_14_w1,fp_14_w2,&          !xyz2sm parameters
   fp_at_nmax,&
   fp_17_width_cutoff,fp_17_nex_cutoff,fp_17_natx_sphere,fp_17_lseg,fp_17_orbital,&
   fp_18_orbital,fp_18_principleev,fp_18_lseg,fp_18_molecules,&
   fp_18_expaparameter,fp_18_nex_cutoff,fp_18_molecules_sphere,fp_18_width_cutoff,&
   fp_18_width_overlap,fp_18_large_vanradius
use confinement
use mod_parini, only: typ_parini
implicit none
type(typ_parini), intent(in):: parini
integer:: i,j
!This routine will check if all values are valid...
if(any(typat(:).lt.1)) stop "Error in typat"
if(any(amu(:).le.0.d0)) stop "Error in amu"
if(any(rcov(:).le.0.d0)) stop "Error in rcov"
if(any(znucl(:).le.0)) stop "Error in znucl"
if(ntime_md.lt.1) stop "Error in ntime_md"
if(md_algo.lt.1.or.md_algo.gt.4) stop "Error in md_algo"
if(md_integrator.lt.1.or.md_integrator.gt.3) stop "Error in md_integrator"
if(ntime_geopt.lt.0) stop "Error in ntime_geopt"
if(bmass.le.0.d0) stop "Error in bmass"
if(mdmin_min.lt.0) stop "Error in mdmin_min"
if(mdmin_max.lt.mdmin_min) stop "Error in mdmin_max"
if(alpha_lat.le.0.d0) stop "Error in alpha_lat"
if(alpha_at.le.0.d0) stop "Error in alpha_at"
if(nsoften.lt.1) stop "Error in nsoften"
if(dtion_md.le.0.d0) stop "Error in dtion_md"
if(nit_per_min.le.0) stop "Error in nit_per_min"
if(trim(geopt_method).ne."FIRE".and.&
   trim(geopt_method).ne."MBFGS".and.&
   trim(geopt_method).ne."RBFGS".and.&
   trim(geopt_method).ne."SQNM".and.&
   trim(geopt_method).ne."QBFGS".and.&
   trim(geopt_method).ne."SD") &
   stop "Error in geopt_method"
if(dtion_fire.lt.dtion_fire_min.or.dtion_fire.gt.dtion_fire_max) stop "Error in dtion_fire"
if(dtion_fire_min.le.0.d0) stop "Error in dtion_fire_min"
if(dtion_fire_max.lt.dtion_fire_min) stop "Error in dtion_fire_max"
if(alphax_lat.le.0.d0) stop "Error in alphax_lat"
if(alphax_at.le.0.d0) stop "Error in alphax_at"
if(tolmxf.le.0.d0) stop "Error in tolmxf"
if(strfact.le.0.d0) stop "Error in strfact"
if(ka.lt.0) stop "Error in ka"
if(kb.lt.0) stop "Error in kb"
if(kc.lt.0) stop "Error in kc"
if(dkpt1.lt.0.d0) stop "Error in dkpt1"
if(dkpt2.lt.0.d0) stop "Error in dkpt2"
if(bc.lt.1.or.bc.gt.3) stop "Error in bc"
if(parini%verb.lt.0.or.parini%verb.gt.3) stop "Error in verb"
if(trim(fp_method_ch).ne."OGANOV".and.trim(fp_method_ch).ne."BCM".and.trim(fp_method_ch).ne."ATORB".and.&
  &trim(fp_method_ch).ne."XYZ2SM".and.trim(fp_method_ch).ne."GAUSS".and.trim(fp_method_ch).ne."COGANOV".and.&
  &trim(fp_method_ch).ne."CAOGANOV".and.trim(fp_method_ch).ne."GOM".and.trim(fp_method_ch).ne."MOLGOM") stop "Error in fp_method_ch"
if(fp_rcut.le.0.d0) stop "Error in fp_rcut"
if(fp_dbin.le.0.d0) stop "Error in fp_dbin"
if(fp_sigma.le.0.d0) stop "Error in fp_sigma"
if(fp_nl.le.0) stop "Error in fp_nl"
if(fp_14_m.lt.1) stop "Error in fp_14_m"
if(fp_14_w1.lt.0.d0) stop "Error in p_14_w1"
if(fp_14_w2.lt.fp_14_w1) stop "Error in p_14_w2"
if(fp_at_nmax.lt.0) stop "Error in fp_at_nmax"
if(trim(fp_17_orbital).ne.'S'.and.trim(fp_17_orbital).ne.'SP') stop "Error in fp_17_orbital"
if(trim(fp_18_orbital).ne.'S'.and.trim(fp_18_orbital).ne.'SP') stop "Error in fp_17_orbital"
if(fp_18_principleev.lt.0) stop "Error in fp_18_principleev"
if(fp_18_molecules.lt.1) stop "Error in fp_18_molecules"
if(fp_18_expaparameter.lt.1) stop "Error in fp_18_expaparameter"
if(fp_18_nex_cutoff.lt.1) stop "Error in fp_18_nex_cutoff"
if(fp_18_molecules_sphere.lt.0) stop "Error in fp_18_molecules_sphere"
if(fp_18_width_cutoff.lt.0.d0) stop "Error in fp_18_width_cutoff"
if(fp_18_width_overlap.lt.0.d0) stop "Error in fp_18_width_overlap"
do i=1,nconfine
  if(.not.(conf_cartred(i).eq."C".or.conf_cartred(i).eq."c".or.&
          &conf_cartred(i).eq."K".or.conf_cartred(i).eq."k".or.&
          &conf_cartred(i).eq."R".or.conf_cartred(i).eq."r".or.&
          &conf_cartred(i).eq."D".or.conf_cartred(i).eq."d")) stop "Error in conf_cartred"
  if(conf_dim(i).lt.1.or.conf_dim(i).gt.3) stop "Error in conf_dim"
  if(conf_exp(i).lt.1) stop "Error in conf_exp"
  if(conf_prefac(i).lt.0.d0) stop "Error in conf_prefac"
  if(conf_av(i).lt.1.or.conf_av(i).gt.2) stop "Error in conf_av"
          do j=1,nat
            if(conf_list(j,i).lt.1.or.conf_list(j,i).gt.nat) stop "Error in conf_list"
          enddo
enddo
!if(ipi_inet.lt.0 .or. ipi_inet.gt.1) stop "Error in ipi_inet: must be 0 for unix socket, 1 for tcp"
!if(ipi_port.lt.1) stop "Error in ipi_port"
!if(ipi_ecutwf(1).lt.0.d0) stop "Error in ipi_ecutwfc"
!if(ipi_ecutwf(2).lt.0.d0) stop "Error in ipi_ecutwfc"

if(sock_inet.lt.0 .or. sock_inet.gt.1) stop "Error in sock_inet: must be 0 for unix socket, 1 for tcp"
if(sock_port.lt.1) stop "Error in sock_port"
if(sock_ecutwf(1).lt.0.d0) stop "Error in sock_ecutwfc"
if(sock_ecutwf(2).lt.0.d0) stop "Error in sock_ecutwfc"
!SQNM
if(sqnm_beta_lat.le.0d0) stop "Error in sqnm_beta_lat"
if(sqnm_beta_at.le.0.d0) stop "Error in sqnm_beta_at"
if(sqnm_nhist.lt.1)      stop "Error in sqnm_nhist"
if(sqnm_maxrise.le.0.d0) stop "Error in sqnm_maxrise"
if(sqnm_cutoffRatio.le.0.d0) stop "Error in sqnm_cutoffRatio"
if(sqnm_steepthresh.le.0.d0) stop "Error in sqnm_steepthresh"
if(sqnm_trustr.le.0.d0)  stop "Error in sqnm_trustr"


end subroutine

!************************************************************************************

subroutine params_echo(parini)
use defs_basis
use String_Utility 
use mod_fire,   only:dtmin, dtmax
use minpar, only:parmin_bfgs
use global, only: target_pressure_habohr,target_pressure_gpa,nat,ntypat,znucl,amu,amutmp,typat,ntime_md,char_type,&
                &nsoften,alpha_at,alpha_lat,ntime_geopt,bmass,mdmin,dtion_fire,dtion_md,tolmxf,strfact,dtion_fire_min,&
                &dtion_fire_max,ka,kb,kc,dkpt1,dkpt2,usewf_geopt,usewf_soften,usewf_md,geopt_method,alphax_at,&
                &alphax_lat,findsym,finddos,auto_soft,mdmin_max,mdmin_min,auto_mdmin,md_algo,md_integrator,auto_dtion_md,&
                &nit_per_min,fixat,fixlat,rcov,mol_soften,fragarr,code,auto_kpt,bc,geopt_ext,energy_conservation,use_confine,&
                &voids,core_rep,md_presscomp
use sqnm,   only: sqnm_beta_lat,sqnm_beta_at,sqnm_nhist,sqnm_maxrise,sqnm_cutoffRatio,sqnm_steepthresh,sqnm_trustr
use qbfgs,  only: qbfgs_bfgs_ndim,qbfgs_trust_radius_max,qbfgs_trust_radius_min,qbfgs_trust_radius_ini,qbfgs_w_1,qbfgs_w_2
use modsocket, only:sock_inet,sock_port,sock_host,sock_ecutwf
use fingerprint, only: & 
   fp_rcut,fp_method,fp_method_ch,fp_nl,&!All
   fp_sigma,fp_dbin,&              !Oganov parameters
   fp_12_nl,&                            !CALYPSO parameters
   fp_13_nl,&                            !Modified CALYPSO parameters
   fp_14_m,fp_14_w1,fp_14_w2,&          !xyz2sm parameters
   fp_at_nmax,&
   fp_17_width_cutoff,fp_17_nex_cutoff,fp_17_natx_sphere,fp_17_lseg,fp_17_orbital,&
   fp_18_orbital,fp_18_principleev,fp_18_lseg,fp_18_molecules,&
   fp_18_expaparameter,fp_18_nex_cutoff,fp_18_molecules_sphere,fp_18_width_cutoff,&
   fp_18_width_overlap,fp_18_large_vanradius
use confinement
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
write(*,'(a,i5)')          " # BOUNDARY      ", bc
write(*,'(a,es15.7)')      " # PRESS         ", target_pressure_gpa
write(*,'(a,L3)')          " # VOIDS         ", voids
write(*,'(a,L3)')          " # COREREP       ", core_rep
write(*,'(a)')             " # COMPUTE parameters ************************************************************"
write(*,'(a,a)')           " # CODE          ", trim(code)
write(*,'(a,L3)')          " # FINDSYM       ", findsym
write(*,'(a,L3)')          " # FINDDOS       ", finddos
write(*,'(a,L3)')          " # USEWFGEO      ", usewf_geopt
write(*,'(a,L3)')          " # USEWFMD       ", usewf_md
write(*,'(a,L3)')          " # USEWFSOFT     ", usewf_soften
write(*,'(a)')             " # Atomic parameters *************************************************************"
write(*,'(a,i5)')          " # NAT           ", nat
write(*,'(a,i5)')          " # NTYPE         ", ntypat
write(formatting,'(a,i5,a)') '(a,',ntypat,'i4)'
write(*,trim(formatting))  " # ZNUCL         ", int(znucl(:))
write(formatting,'(a,i5,a)') '(a,',nat,'i3)'
write(*,trim(formatting))  " # TYPAT         ", typat(:)
write(formatting,'(a,i5,a)') '(a,',ntypat,'i4)'
write(*,trim(formatting))  " # AMU           ", int(amu(:))
write(*,'(a)')             " # MD parameters *****************************************************************"
write(*,'(a,i5)')          " # MDALGO        ", md_algo
write(*,'(a,es15.7)')      " # MDPRESSCOMP   ", md_presscomp
write(*,'(a,i5)')          " # MDINT         ", md_integrator
write(*,'(a,i5)')          " # MDNIT         ", ntime_md
write(*,'(a,L3)')          " # AUTO_MDDT     ", auto_dtion_md
write(*,'(a,L3)')          " # MDENCON       ", energy_conservation
write(*,'(a,es15.7)')      " # MDDTINIT      ", dtion_md
write(*,'(a,i5)')          " # MDDTIPM       ", nit_per_min
write(*,'(a,L3)')          " # AUTO_MDMIN    ", auto_mdmin
write(*,'(a,i5)')          " # MDMININIT     ", mdmin
write(*,'(a,i5)')          " # MDMINMIN      ", mdmin_min
write(*,'(a,i5)')          " # MDMINMAX      ", mdmin_max
write(*,'(a,es15.7)')      " # CELLMASS      ", bmass
write(*,'(a)')             " # GEOPT parameters **************************************************************"
write(*,'(a,L3)')          " # GEOEXT        ", geopt_ext
write(*,'(a,i5)')          " # GEONIT        ", ntime_geopt
write(*,'(a,es15.7)')      " # GEOTOLMXF     ", tolmxf
write(*,'(a,es15.7)')      " # STRFACT       ", strfact
if(.not.geopt_ext) then
write(*,'(a,a)')           " # GEOALGO       ", trim(geopt_method) 
if(trim(geopt_method)=="FIRE") then
write(*,'(a,es15.7)')      " # GEOFIREDTINIT ", dtion_fire
write(*,'(a,es15.7)')      " # GEOFIREDTMIN  ", dtion_fire_min
write(*,'(a,es15.7)')      " # GEOFIREDTMAX  ", dtion_fire_max
elseif(trim(geopt_method)=="MBFGS".or.&
       trim(geopt_method)=="RBFGS".or.&
       trim(geopt_method)=="SQNM".or. &
       trim(geopt_method)=="SD") then
write(*,'(a,es15.7)')      " # GEOHESSAT     ", alphax_at
write(*,'(a,es15.7)')      " # GEOHESSLAT    ", alphax_lat
endif
if(trim(geopt_method)=="QBFGS") then
write(*,'(a,i5)'    )      " # GEOQBFGSNDIM  ",qbfgs_bfgs_ndim
write(*,'(a,es15.7)')      " # GEOQBFGSTRI   ",qbfgs_trust_radius_ini
write(*,'(a,es15.7)')      " # GEOQBFGSTRMIN ",qbfgs_trust_radius_min
write(*,'(a,es15.7)')      " # GEOQBFGSTRMAX ",qbfgs_trust_radius_max
write(*,'(a,es15.7)')      " # GEOQBFGSW1    ",qbfgs_w_1
write(*,'(a,es15.7)')      " # GEOQBFGSW2    ",qbfgs_w_2
endif
if(trim(geopt_method)=="SQNM") then
write(*,'(a,i5)')          " # GEOSQNMNHIST  ",sqnm_nhist
write(*,'(a,es15.7)')      " # GEOSQNMMAXRISE",sqnm_maxrise
write(*,'(a,es15.7)')      " # GEOSQNMCUTOFF ",sqnm_cutoffRatio
write(*,'(a,es15.7)')      " # GEOSQNMSTEEP  ",sqnm_steepthresh
write(*,'(a,es15.7)')      " # GEOSQNMTRUSTR ",sqnm_trustr
endif
endif 
write(*,'(a)')             " # SOFTEN parameters *************************************************************"
write(*,'(a,L3)')          " # AUTO_SOFT     ", auto_soft
write(*,'(a,L3)')          " # MOLSOFT       ", mol_soften
write(*,'(a,es15.7)')      " # SOFTAT        ", alpha_at
write(*,'(a,es15.7)')      " # SOFTLAT       ", alpha_lat
write(*,'(a,i5)')          " # SOFTNIT       ", nsoften
write(*,'(a)')             " # KPOINTS parameters ************************************************************"
write(*,'(a,L3)')          " # AUTO_KPT     ", auto_kpt
if(.not.auto_kpt) then
write(*,'(a,3i5)')         " # KPTMESH      ", ka,kb,kc
else
write(*,'(a,2es15.7)')     " # KPTDEN       ", dkpt1,dkpt2
endif
write(*,'(a)')             " # FINGERPRINT parameters ********************************************************"
write(*,'(a,a)')           " # FPMETHOD      ", trim(fp_method_ch)
if(bc==1.or.bc==3) then
write(*,'(a,es15.7)')      " # FPCUT         ", fp_rcut
endif
if(trim(fp_method_ch)=="OGANOV") then
write(*,'(a,es15.7)')      " # FPDBIN        ", fp_dbin
write(*,'(a,es15.7)')      " # FPSIGMA       ", fp_sigma
elseif(trim(fp_method_ch)=="BCM".or.trim(fp_method_ch)=="ATORB") then
write(*,'(a,i5)')          " # FPNL          ", fp_nl
elseif(trim(fp_method_ch)=="XYZ2SM") then
write(*,'(a,i5)')          " # FPPOWER       ", fp_14_m
write(*,'(a,es15.7)')      " # FPGAUSSFAC1   ", fp_14_w1
write(*,'(a,es15.7)')      " # FPGAUSSFAC2   ", fp_14_w2
elseif(trim(fp_method_ch)=="COGANOV") then
write(*,'(a,es15.7)')      " # FPSIGMA       ", fp_sigma
write(*,'(a,i5)')          " # FPATNMAX      ", fp_at_nmax
elseif(trim(fp_method_ch)=="CAOGANOV") then
write(*,'(a,es15.7)')      " # FPSIGMA       ", fp_sigma
write(*,'(a,i5)')          " # FPATNMAX        ", fp_at_nmax
elseif(trim(fp_method_ch)=="GOM") then
write(*,'(a,i5)')          " # FPNATX        ", fp_17_natx_sphere
write(*,'(a,i5)')          " # FPLSEG        ", fp_17_lseg
write(*,'(a,a)')           " # FPORBITAL     ", fp_17_orbital
write(*,'(a,es15.7)')      " # FPNEXCUT      ", fp_17_nex_cutoff
elseif(trim(fp_method_ch)=="MOLGOM") then
write(*,'(a,a)')           " # FPORBITAL     ", fp_18_orbital
write(*,'(a,i5)')          " # FPNEXCUT      ", fp_18_nex_cutoff
write(*,'(a,i5)')          " # FPPRINCIPLEEV ", fp_18_principleev
write(*,'(a,i5)')          " # FPMOLECULES   ", fp_18_molecules
write(*,'(a,i5)')          " # FPEXPA        ", fp_18_expaparameter
write(*,'(a,i5)')          " # FPMOLSPHERE   ", fp_18_molecules_sphere
write(*,'(a,es15.7)')      " # FPWIDTHCUT    ", fp_18_width_cutoff
write(*,'(a,es15.7)')      " # FPWIDTHOVER   ", fp_18_width_overlap
endif
write(*,'(a)')             " # CONFINEMENT parameters ********************************************************"
write(*,'(a,L3)')          " # CONFINEMENT   ",use_confine
if(use_confine) then
write(*,'(a,i5)')          " # CONFNCONF     ",nconfine
write(formatting,'(a,i5,a)') '(a,',nconfine,'A4)'
write(*,trim(formatting))  " # CONFCARTRED   ",conf_cartred
write(formatting,'(a,i5,a)') '(a,',nconfine,'i4)'
write(*,trim(formatting))  " # CONFDIM       ",conf_dim
write(formatting,'(a,i5,a)') '(a,',nconfine,'i4)'
write(*,trim(formatting))  " # CONFEXP       ",conf_exp
write(formatting,'(a,i5,a)') '(a,',nconfine,'es9.2)'
write(*,trim(formatting))  " # CONFPREFAC    ",conf_prefac
write(formatting,'(a,i5,a)') '(a,',nconfine,'es9.2)'
write(*,trim(formatting))  " # CONFCUT       ",conf_cut
write(formatting,'(a,i5,a)') '(a,',nconfine,'i4)'
write(*,trim(formatting))  " # CONFAV        ",conf_av
write(formatting,'(a,i5,a)') '(a,',nconfine,'es9.2)'
write(*,trim(formatting))  " # CONFEQ        ",conf_eq
write(formatting,'(a,i5,a)') '(a,',nconfine,'i4)'
write(*,trim(formatting))  " # CONFNAT       ",conf_nat
!Go through all number of confinements
   do i=1,nconfine
      if(i.lt.10) then
        write(fn1,'(i1.1)') i
        write(string,'(a)') "CONFLIST"//fn1
        write(formatting,'(a,i5,a)') '(a,',conf_nat(i),'i4)'
        write(*,trim(formatting))  " # "//trim(string)//"     ",conf_list(1:conf_nat(i),i)
      else
        stop "Loop not implemented for nconfine greater than 9"
      endif
   enddo
endif
if(StrLowCase(trim(adjustl(code)))=="ipi") then
write(*,'(a)')             " # IPI parameters ****************************************************************"
write(formatting,'(a)') '(a,i4)'
write(*,trim(formatting))  " # IPIINET       ",sock_inet
write(formatting,'(a)') '(a,i8)'
write(*,trim(formatting))  " # IPIPORT       ",sock_port
write(formatting,'(a)') '(a,a)'
write(*,trim(formatting))  " # IPIHOST       ",trim(adjustl(sock_host))
write(formatting,'(a)') '(a,2f10.4)'
write(*,trim(formatting))  " # IPIECUTWF     ",sock_ecutwf(:)
endif
if(StrLowCase(trim(adjustl(code)))=="msock") then
write(*,'(a)')             " # MSOCK parameters **************************************************************"
write(formatting,'(a)') '(a,i4)'
write(*,trim(formatting))  " # SOCKINET      ",sock_inet
write(formatting,'(a)') '(a,i8)'
write(*,trim(formatting))  " # SOCKPORT      ",sock_port
write(formatting,'(a)') '(a,a)'
write(*,trim(formatting))  " # SOCKHOST      ",trim(adjustl(sock_host))
write(formatting,'(a)') '(a,2f10.4)'
write(*,trim(formatting))  " # SOCKECUTWF    ",sock_ecutwf(:)
endif
write(*,'(a)')             " ############################ END Echo params_new.in #############################"
end subroutine

!************************************************************************************
subroutine fp_assign()
use fingerprint 
use global, only: bc
implicit none
if(bc==2) then !Molecular systems
  fp_method=21
  fp_method_ch="GAUSS"
elseif(bc==1) then !Crystals
  if(trim(fp_method_ch)=="OGANOV") then
    fp_method=11
  elseif(trim(fp_method_ch)=="BCM") then
    fp_method=12
  elseif(trim(fp_method_ch)=="ATORB") then
    fp_method=13
  elseif(trim(fp_method_ch)=="XYZ2SM") then
    fp_method=14
  elseif(trim(fp_method_ch)=="COGANOV") then
    fp_method=15
  elseif(trim(fp_method_ch)=="CAOGANOV") then
    fp_method=16
  elseif(trim(fp_method_ch)=="GOM") then
    fp_method=17
    if(trim(fp_17_orbital)=="S")then
    fp_17_lseg=1
    elseif(trim(fp_17_orbital)=="SP")then
    fp_17_lseg=4 
    endif
    fp_17_width_cutoff=fp_rcut/sqrt(2.d0*fp_17_nex_cutoff)
  elseif(trim(fp_method_ch)=="MOLGOM") then
    fp_method=18
    if(trim(fp_18_orbital)=="S")then
      fp_18_lseg=1
    elseif(trim(fp_18_orbital)=="SP")then
      fp_18_lseg=4 
    endif
!    fp_18_width_cutoff=fp_rcut/sqrt(2.d0*fp_18_nex_cutoff)
  endif
endif    
end subroutine

