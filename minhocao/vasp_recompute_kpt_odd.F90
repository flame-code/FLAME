program recompute_kpt
implicit none
real(8):: tmp,dkpt1,dkpt2,latvec(3,3),angbohr,scaling,dkpt_12(2)
integer:: ka,kb,kc,vasp_kpt_mode,n,kpt_abc(3)
logical:: tmpl,auto_kpt,found
character(8):: tmp_ch 
character(250):: all_line
angbohr=1.d0/0.52917720859d0

open(unit=12,file="params_new.in")
  do while(.true.)
   read(12,'(a250)',end=97)all_line
   n = len_trim(all_line)
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
enddo

!To generate automatic kpoin mesh
! if(ka==0.and.kb==0.and.kc==0)then
! read(12,*) dkpt1,dkpt2
! else
! dkpt1=0.d0
! dkpt2=0.d0
! endif
 97 continue
  if(AUTO_KPT) then
    ka=0;kb=0;kc=0
  else
    dkpt1=0.d0
    dkpt2=0.d0
  endif

 close(12)
 vasp_kpt_mode=2

 open(unit=12,file="CONTCAR")
 read(12,*) tmp_ch
 read(12,*) scaling
 read(12,*) latvec(:,1)
 read(12,*) latvec(:,2)
 read(12,*) latvec(:,3)
 latvec=latvec*scaling
 latvec=latvec*angbohr
call make_kpt_vasp_odd(latvec,dkpt1,ka,kb,kc,vasp_kpt_mode)
end 



!!!program recompute_kpt
!!!implicit none
!!!real(8):: tmp,dkpt1,dkpt2,latvec(3,3),angbohr,scaling
!!!integer:: ka,kb,kc,vasp_kpt_mode,n
!!!logical:: tmpl
!!!character(8):: tmp_ch 
!!!character(250):: all_line
!!!angbohr=1.d0/0.52917720859d0
!!!
!!! open(unit=12,file="params.in")
!!!  read(12,'(a250)') all_line             !Target pressure in GPA
!!!  read(12,'(a250)') all_line             !Number of atoms
!!!  read(12,'(a250)') all_line             !Number of atom types
!!!  read(12,'(a250)') all_line             !Nuclei charge 
!!!  read(12,'(a250)') all_line             !Atomic mass used for MD and Fire, if all 0, then automatic determination
!!!  read(12,'(a250)') all_line             !Types of atoms  
!!!  read(12,'(a250)') all_line             !Maximum number of iterations during MD
!!!  read(12,'(a250)') all_line             !Maximum number of iterations during GEOPT
!!!  read(12,'(a250)') all_line             !Cell mass during MD and FIRE
!!!  read(12,'(a250)') all_line             !Number of enthalpy minima crossed unit stop MD
!!!  read(12,'(a250)') all_line             !nsoften, alpha_at, alpha_lat
!!!  read(12,'(a250)') all_line             !Initial timestep for MD
!!!  read(12,'(a250)') all_line             !Method of geopt: FIRE or BFGS
!!!  read(12,'(a250)') all_line             !Initial timestep for FIRE, Min, Max
!!!  read(12,'(a250)') all_line             !Force tolerance for GEOPT convergance 
!!!  read(12,'(a250)') all_line             !Factor to multiply stress 
!!!  read(12,'(a250)') all_line             !Determines if the previous wavefunction should be read for GEOPT, SOFTENING and MD
!!!  read(12,'(a250)') all_line             !Determines if symmetry (isotropy)/FDOS calculation should be used
!!!! read(12,*) ka,kb,kc              !For fixed kpoint mesh
!!!!Block KPT****************
!!! read(12,'(a250)') all_line
!!! n = len_trim(all_line)
!!! read(all_line(1:n),*) tmp_ch
!!!  if(trim(tmp_ch)=="Auto") then
!!!    read(all_line(1:n),*) tmp_ch,dkpt1,dkpt2
!!!    ka=0;kb=0;kc=0
!!!  else
!!!   read(all_line(1:n),*) ka,kb,kc
!!!    dkpt1=0.d0
!!!    dkpt2=0.d0
!!!  endif
!!!!Block KPT****************
!!!
!!!
!!!!To generate automatic kpoin mesh
!!!! if(ka==0.and.kb==0.and.kc==0)then
!!!! read(12,*) dkpt1,dkpt2
!!!! else
!!!! dkpt1=0.d0
!!!! dkpt2=0.d0
!!!! endif
!!! close(12)
!!! vasp_kpt_mode=2
!!!
!!! open(unit=12,file="CONTCAR")
!!! read(12,*) tmp_ch
!!! read(12,*) scaling
!!! read(12,*) latvec(:,1)
!!! read(12,*) latvec(:,2)
!!! read(12,*) latvec(:,3)
!!! latvec=latvec*scaling
!!! latvec=latvec*angbohr
!!!call make_kpt_vasp(latvec,dkpt1,ka,kb,kc,vasp_kpt_mode)
!!!end 

!************************************************************************************
subroutine make_kpt_vasp_odd(latvec,dkpt,ka,kb,kc,vasp_kpt_mode)
!This routine will append some informations to a file already containing some informations about the abininit runs
!The informations appended are:
!-The atomic informations
!A file KPOINTS is generated. There are two options available:
!Automated generation if vasp_kpt_mode==1
!Monkhorst pack if vasp_kpt_mode==2 
!ATTENTION:
!The meaning of dkpt1 and dkpt2 is different depending on vasp_kpt_mode:
!accuracy is given by the integer length of dkpt for vasp_kpt_mode==1 (10 for insulators, 100 for metals)
!accuracy is 2pi/bohr*dkpt for vasp_kpt_mode==2 
implicit none
real(8):: dproj(6),acell(3),rprim(3,3),latvec(3,3),dkpt,angbohr
integer:: iat,iprec,ka,kb,kc,itype,vasp_kpt_mode
logical:: getwfk
character(1):: fn
character(150):: command
angbohr=1.d0/0.52917720859d0 

!Kpoint mesh
open(unit=87,file="KPOINTS")
write(87,'(a,i5)') "# Definition of the k-point mesh ",vasp_kpt_mode
write(87,'(i5)') 0
if(dkpt==0.d0) then
write(87,'(a)') "Gamma"!"Monkhorst Pack"
write(87,'(3(1x,i5),a)') ka,kb,kc,"  # Number of gridpoints in each dimension"
write(87,'(3(1x,i5),a)') 0,0,0,"  # Shifts"
elseif(vasp_kpt_mode==2) then
call find_kpt(ka,kb,kc,latvec,dkpt)
if(MOD(ka,2)==0) ka=ka+1
if(MOD(kb,2)==0) kb=kb+1
if(MOD(kc,2)==0) kc=kc+1
write(87,'(a)') "Gamma"!"Monkhorst Pack"
write(87,'(3(1x,i5),a)') ka,kb,kc,"  # Number of gridpoints in each dimension"
write(87,'(3(1x,i5),a)') 0,0,0,"  # Shifts"
elseif(vasp_kpt_mode==1) then
write(87,'(a)') "Auto"
write(87,'(i5,a)') dkpt," # K-mesh length"
else
stop "Wrong kpt option"
endif
write(*,'(a,3(1x,i5))') " # KPT mesh re-set as follows: ",ka,kb,kc
close(87)
end

subroutine find_kpt(k1,k2,k3,lat,gridden)
!This code will define the KPT mesh based on the desired grid density
  implicit none
  integer       :: k1,k2,k3,i,j
  real(8)       :: lat1, lat2, lat3
  real(8)       :: angles(3),cos_arr(3)
  real(8)       :: lat(3,3),glat(3,3),crossp(3),vol,gridden,glen(3),a(3,3)
  real(8)       :: pi
  pi=acos(-1.d0)

  lat1=sqrt(lat(1,1)**2+lat(2,1)**2+lat(3,1)**2)
  lat2=sqrt(lat(1,2)**2+lat(2,2)**2+lat(3,2)**2)
  lat3=sqrt(lat(1,3)**2+lat(2,3)**2+lat(3,3)**2)
  cos_arr(1)=(lat(1,2)*lat(1,3)+lat(2,2)*lat(2,3)+lat(3,2)*lat(3,3))/lat2/lat3
  cos_arr(2)=(lat(1,1)*lat(1,3)+lat(2,1)*lat(2,3)+lat(3,1)*lat(3,3))/lat1/lat3
  cos_arr(3)=(lat(1,1)*lat(1,2)+lat(2,1)*lat(2,2)+lat(3,1)*lat(3,2))/lat1/lat2
  angles(1)=(acos(cos_arr(1))/pi)*180.d0
  angles(2)=(acos(cos_arr(2))/pi)*180.d0
  angles(3)=(acos(cos_arr(3))/pi)*180.d0

  a=lat
  vol= a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)+&
       a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)
  call cross_product(lat(:,2),lat(:,3),crossp(:))
  glat(:,1)=2.d0*pi*crossp(:)/vol
  call cross_product(lat(:,3),lat(:,1),crossp(:))
  glat(:,2)=2.d0*pi*crossp(:)/vol
  call cross_product(lat(:,1),lat(:,2),crossp(:))
  glat(:,3)=2.d0*pi*crossp(:)/vol
!Compute the correct kpts
  glen(1)=sqrt(glat(1,1)**2+glat(2,1)**2+glat(3,1)**2)
  glen(2)=sqrt(glat(1,2)**2+glat(2,2)**2+glat(3,2)**2)
  glen(3)=sqrt(glat(1,3)**2+glat(2,3)**2+glat(3,3)**2)
  call track_kpt(gridden,glen(1),k1)
  call track_kpt(gridden,glen(2),k2)
  call track_kpt(gridden,glen(3),k3)
end subroutine find_kpt

subroutine track_kpt(gridden,glen,kpt)
  implicit none
  integer :: kpt,j
  real(8) :: glen,gridden,d_test,pi 
  pi=acos(-1.d0)
  kpt=int(glen/(gridden*2.d0*pi))
  if (kpt == 0) kpt = 1
  d_test=glen/(kpt*2.d0*pi)
  if (d_test.ge.gridden) then
     do j=1,25
        kpt=kpt+j
        d_test=glen/(kpt*2.d0*pi)
        if (d_test.le.gridden) exit
     enddo
  endif
end subroutine track_kpt
!************************************************************************************

 subroutine cross_product(a,b,crossp)
 !a very simple implementation of the cross product
 implicit none
 real(8)::a(3),b(3)
 real(8)::crossp(3)
 crossp(1)=a(2)*b(3)-a(3)*b(2)
 crossp(2)=a(3)*b(1)-a(1)*b(3)
 crossp(3)=a(1)*b(2)-a(2)*b(1)
 return
 end subroutine

