module global_expand
implicit none
logical:: reduced
integer:: nat
integer,allocatable:: fragarr(:),fragarr_in(:)
integer:: verb                  !0: very little output, 1: normal output, 2: folders for geopt and md, 3: output stress and forces
end module

program main
use global_expand
!this routine will take all the poslow files in the current folder and expand
!them into a multiple in each dimension and overwrite the current files.
implicit none
integer:: i, nat_in,expansions(3),ntypat,iat,k,l,m
real(8),allocatable::xred_in(:,:),xred(:,:),xcart(:,:),fcart_in(:,:),fcart(:,:),posnoise(:,:)
real(8)::displ,latvec_in(3,3),latvec(3,3),strten_in(6),strten(6),energy,enthalpy,ext_press
logical,allocatable:: fixat_in(:),fixat(:)
logical:: fixlat(7)             !Contains the information of the cell constraints: a,b,c,alpha,beta,gamma,cellshape
logical:: file_exists
logical:: readfrag, readfix
character(5):: fn5
character(40):: units,filename
character(2),allocatable:: char_typat(:),char_typat_in(:)
real(8):: target_pressure_habohr
!Set units
units="angstroem"

!Get the number of atoms from some poslow file
do i=1,10000000
        write(fn5,'(i5.5)') i
        filename = 'poslow'//fn5//'.ascii'
        INQUIRE(FILE=trim(filename), EXIST=file_exists)  
        if(file_exists) then
          open(unit=2,file=trim(filename))
          read(2,*) nat_in
          close(2)
          exit
        endif
enddo

!Read the file containing the expansions and decide if random displacements need
!to be made
displ=0.d0
open(unit=2,file="expansion.in")
read(2,*) expansions(:)
read(2,*) displ
close(2)


!Allocate all necessary arrays
nat=expansions(1)*expansions(2)*expansions(3)*nat_in
allocate(xred_in(3,nat_in),xred(3,nat),fcart_in(3,nat_in),fcart(3,nat),fragarr_in(nat_in),fragarr(nat),char_typat_in(nat_in),char_typat(nat),fixat_in(nat_in),fixat(nat),posnoise(3,nat),&
        &xcart(3,nat))




!Now run over all poslows found in the folder
do i=1,10000000
        write(fn5,'(i5.5)') i
        filename = 'poslow'//fn5//'.ascii'
        INQUIRE(FILE=trim(filename), EXIST=file_exists)  
        if(.not.file_exists) then
                write(*,*) trim(filename)//" not found! Exiting"
                stop
        endif
        readfix=.false.
        readfrag=.false.
   if(file_exists) then
   call read_atomic_file_ascii(filename,nat_in,units,xred_in,latvec_in,fcart_in,strten_in,char_typat_in,fixat_in,fixlat,readfix,fragarr,readfrag,enthalpy,energy)
   write(*,'(a)') " # Now running "//trim(filename)
   call k_expansion(nat_in,latvec_in,xred_in,expansions(1),expansions(2),expansions(3),latvec,xcart)
   do k=1,nat
   char_typat(k)=char_typat_in(modulo(k-1,nat_in)+1)
   fixat(k)=fixat_in(modulo(k-1,nat_in)+1)
   fragarr(k)=k
   enddo
   !Add noise to the atoms
   call random_number(posnoise)
   posnoise=displ*(posnoise-0.5d0)
   xcart=xcart+posnoise
   call rxyz_cart2int(latvec,xred,xcart,nat)
   !Write the new atomic file
!       filename = trim(filename)//".out"
       !WARNING: the following is probably broken since target_pressure_habohr is a local variable.
       call write_atomic_file_ascii(filename,nat,units,xred,latvec,fcart,strten,&
             &char_typat,fixat,fixlat,energy,target_pressure_habohr,ext_press,enthalpy)
endif
enddo




end program
subroutine k_expansion(nat,latvec,xred,ka,kb,kc,k_latvec,k_xcart)
!This routine expands the cell defined by latevec to a supercell
!of dimension ka,kb,kc. The atomic positions in real space
!will be returned in k_xcart. k_nat will then be the number of
!all atoms in the supercell and is ka*kb*kc*nat
implicit none
real(8):: latvec(3,3),k_latvec(3,3),k_xcart(3,nat,ka,kb,kc),xred(3,nat)
integer:: iat,k,l,m,ka,kb,kc,nat
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
      
!**********************************************************************************************

subroutine read_atomic_file_ascii(filename,nat,units,xred,latvec,fcart,strten,char_typat,&
           &fixat,fixlat,readfix,fragarr,readfrag,printval1,printval2)
!subroutine read_atomic_file_ascii(filename,nat,units,xred,latvec,acell,rprim,printval1,printval2)
!This routine will read the file "filename" ins ascii file format
!The output unit will always be in atomic units (bohr)
!So if units==angstroem, the file will be converted
!   if units==bohr, the positions will not be changed
use defs_basis, only: Bohr_Ang,Ha_eV
use global_expand, only: reduced
implicit none
integer:: nat,natin,iat,ierror,io,n,k,fragarr(nat),fragarr_tmp,lhead(nat),llist(nat),nmol
logical:: fixat(nat),fixlat(7),readfix,reduced_tmp,readfrag
character(1) :: isfixed,ch_tmp
character(2) :: element,char_typat(nat)
character(40):: filename,units
real(8):: pos(3,nat),xred(3,nat),latvec(3,3),dproj(6),strten(6),fcart(3,nat)
real(8):: angbohr,evhartree,enthalpy_at,printval1,printval2
character(200):: line
reduced_tmp=.false.
if(readfrag) fragarr=-1
if(trim(units)=="angstroem") then
!  angbohr=1.889725989d0
  angbohr=1.d0/Bohr_Ang
  evhartree=1.d0/Ha_eV
elseif(trim(units)=="bohr") then
  angbohr=1.d0
  evhartree=1.d0
else
  stop "Wrong file unit format"
endif

open(unit=46,file=trim(filename),iostat=ierror)
 if (ierror /= 0) then
    write(*,*) ' COULD not read file ',filename
    stop
 end if

!If the input file is called "poscur.ascii", then dont read the input enthalpy and energy, since they might not be provided
if(trim(filename)=="poscur.ascii") then
read(46,*)natin
else
read(46,*)natin, printval1,printval2
endif
if(natin.ne.nat) stop "Number of atoms not consistent with params.in input"
read(46,*) dproj(1:3)
read(46,*) dproj(4:6)
do iat=1,nat
1010  read(46,'(a200)') line

  if(iat==1) then  
    n = len_trim(line)
    k = index(line(1:n),"reduced")
    if(k.ne.0) then 
         reduced_tmp=.true.
         goto 1010
    endif
  endif

  if(iat==1) then  
    n = len_trim(line)
    k = index(line(1:n),"fixlat",back=.true.)
    if(k.ne.0) then
      k = scan(line(1:n),"t",.true.)
      if(readfix) then
         read(line(k+1:n),*,iostat=io) fixlat(:)
         if(io.lt.0) stop "Could not read lattice constraints"
      endif 
      goto 1010
    endif
  endif

  if(readfix) then
     read(line,*,iostat=io) pos(:,iat),char_typat(iat),isfixed
     if(io.lt.0)  isfixed="a"
     if(isfixed=="f" .or. isfixed=="F")  then
         fixat(iat)=.true.
     else
         fixat(iat)=.false.
     endif
     read(line,*,iostat=io) pos(:,iat),char_typat(iat),fragarr_tmp,isfixed
     if(io.lt.0)  isfixed="a"
     if(isfixed=="f" .or. isfixed=="F")  then
         fixat(iat)=.true.
     else
         fixat(iat)=.false.
     endif
  endif
  if(readfrag) then
     read(line,*,iostat=io) pos(:,iat),char_typat(iat),fragarr_tmp
     if(io.lt.0) then 
         fragarr_tmp=-1 
     else
         fragarr(iat)=fragarr_tmp
     endif
     read(line,*,iostat=io) pos(:,iat),char_typat(iat),isfixed,fragarr_tmp
     if(io.lt.0) then 
         fragarr_tmp=-1 
     else
         fragarr(iat)=fragarr_tmp
     endif
  endif
  read(line,*,iostat=io) pos(:,iat),char_typat(iat)
      if(io.lt.0) then 
         write(*,*) "File read error in",trim(filename)
         stop
      endif 
  if(any(fixat).and.fragarr_tmp.ge.1) stop "Fragmented structure not supported with fixed atoms"
enddo
strten=0.d0
fcart=0.d0
if(trim(filename).ne."poscur.ascii") then
  read(46,*,end=99) ch_tmp,strten(:)
  do iat=1,nat
     read(46,*,end=99) ch_tmp,fcart(:,iat)
  enddo
  fcart=fcart*evhartree/angbohr
endif
99 close(46)


dproj=dproj*angbohr
call dproj2latvec(dproj,latvec)
if(reduced_tmp) then
xred=pos
else
pos=pos*angbohr
call rxyz_cart2int(latvec,xred,pos,nat)
endif

if(trim(filename)=="poscur.ascii") reduced=reduced_tmp
!!if(readfrag) then
!!call refragment(fragarr,nat)
!!call make_linked_list(fragarr,lhead,llist,nat,nmol)
!!endif
!stop
end subroutine


!************************************************************************************

subroutine write_atomic_file_ascii(filename,nat,units,xred,latvec0,fcart,strten,char_typat,&
           &fixat,fixlat,energy,pressure,printval1,printval2)
!This routine will write the file "filename" in ascii file format
!The input unit will always be in atomic units (bohr, hartree), but the output can be specified by the vaule in "units"
!So if units==angstroem, the file will be converted to angstroem
!   if units==bohr, the positions will not be changed
use defs_basis, only: Ha_eV,Bohr_Ang,HaBohr3_GPa
use global_expand, only: reduced,fragarr,verb

implicit none
integer:: nat,natin,iat
character(40):: filename,units
real(8):: pos(3,nat),xred(3,nat),latvec(3,3),latvec0(3,3),dproj(6),rotmat(3,3),v(3,3),ucvol,fcart(3,nat),strten(6)
real(8):: angbohr,hartree2ev,in_GPA,in_ang3,int_press
real(8):: energy, etotal, enthalpy, enthalpy_at,pressure,printval1,printval2
character(2):: char_typat(nat)
logical:: fixat(nat),fixlat(7)
if(trim(units)=="angstroem") then
!  angbohr=1.d0/1.889725989d0
!  hartree2ev=27.211396132d0
!  in_GPA=29421.033d0
!  in_ang3=0.148184743d0
!  int_press=160.217646200d0
  angbohr=Bohr_Ang
  hartree2ev=Ha_eV
  in_GPA=HaBohr3_GPa
  in_ang3=Bohr_Ang**3
  int_press=HaBohr3_GPa/Ha_eV*Bohr_Ang**3
elseif(trim(units)=="bohr") then
  angbohr=1.d0
  hartree2ev=1.d0
  in_GPA=1.d0
  in_ang3=1.d0
  int_press=1.d0
else
  stop "Wrong file unit format"
endif

latvec=latvec0

do iat=1,nat
   pos(:,iat)=matmul(latvec,xred(:,iat))
enddo

call latvec2dproj(dproj,latvec,rotmat,pos,nat)

!Compute cell volume
 v=latvec
 ucvol= v(1,1)*v(2,2)*v(3,3)-v(1,1)*v(2,3)*v(3,2)-v(1,2)*v(2,1)*v(3,3)+&
      v(1,2)*v(2,3)*v(3,1)+v(1,3)*v(2,1)*v(3,2)-v(1,3)*v(2,2)*v(3,1)
!Convert units
!!!The second entry in the first line is the enthalpy in atomic units
!!enthalpy_at=energy+pressure*ucvol
dproj=dproj*angbohr
pos=pos*angbohr
if(reduced) pos=xred
etotal=energy*hartree2ev
enthalpy=energy*hartree2ev+pressure*in_GPA/int_press*ucvol*in_ang3
ucvol=ucvol*in_ang3

open(unit=46,file=trim(filename))
if(trim(units)=="angstroem") then
  write(46,'(i4,1x,es25.15,1x,es25.15,1x,a20,es15.7,a14,es15.7,a13,es15.7)')&
  & nat,printval1,printval2,"         energy(eV)=",etotal," enthalpy(eV)=",enthalpy,&
  &" ucvol(ang3)=",ucvol
elseif(trim(units)=="bohr") then
  write(46,'(i4,1x,es25.15,1x,es25.15,1x,a20,es15.7,a14,es15.7,a13,es15.7)')&
  & nat,printval1,printval2,"         energy(Ha)=",etotal," enthalpy(Ha)=",enthalpy,&
  &" ucvol(bhr3)=",ucvol
endif
write(46,*) dproj(1:3)
write(46,*) dproj(4:6)
if(reduced) write(46,'(a)') "#keyword:reduced"
if(any(fixlat)) write(46,'(a,7(1x,L))') "#keyword:fixlat ",fixlat(:)
do iat=1,nat
  if(fixat(iat)) then
      write(46,'(3(1x,es25.15),2x,a2,2x,a1)') pos(:,iat),trim(char_typat(iat)),"F"
  elseif(all(fragarr.gt.0)) then 
      write(46,'(3(1x,es25.15),2x,a2,2x,i5)')       pos(:,iat),trim(char_typat(iat)),fragarr(iat)
  else 
      write(46,'(3(1x,es25.15),2x,a2)')       pos(:,iat),trim(char_typat(iat))
  endif
enddo
if(verb.ge.3) then
write(46,'(a,6(1x,es25.15),a)')  "# ",strten(:)," Stress tensor in GPa"
write(46,'(a,3(1x,es25.15),a)')  "# ",matmul(rotmat,fcart(:,1))*hartree2ev/angbohr," Forces in file units"
do iat=2,nat
   write(46,'(a,3(1x,es25.15),a)')  "# ",matmul(rotmat,fcart(:,iat))*hartree2ev/angbohr
enddo
endif
close(46)
end subroutine

!************************************************************************************

subroutine get_char_type(filename,nat,char_type,typat,ntypat)
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
!************************************************************************************

 subroutine rxyz_int2cart(latvec,rxyzint,rxyzcart,nat)
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
!************************************************************************************

 subroutine dproj2latvec(dproj,latvec)
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
 !a very simple implementation of the cross product
 implicit none
 real(8)::a(3),b(3)
 real(8)::crossp(3)
 crossp(1)=a(2)*b(3)-a(3)*b(2)
 crossp(2)=a(3)*b(1)-a(1)*b(3)
 crossp(3)=a(1)*b(2)-a(2)*b(1)
 return
 end subroutine
!************************************************************************************

 subroutine rotation(rotmat,angle,axe)
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
