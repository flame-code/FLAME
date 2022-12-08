program main
implicit none
integer:: nat,nsteps,istep,iat,ntype,i,n,j,k
real(8),allocatable:: xred1(:,:)
real(8)::latvec1(3,3),dproj1(6)
real(8),allocatable:: xred2(:,:),xcart(:,:)
real(8)::latvec2(3,3),dproj2(6)
real(8),allocatable::xred(:,:)
real(8)::latvec_step(3,3),latvec(3,3),scaling
real(8):: printval1,printval2,energy,pressure
character(140):: filename1,filename2,filename,units,tmp_ch
character(2),allocatable::typchar(:),typ_at(:)
character(4):: fn
integer, allocatable:: itype(:),nitype(:)
character(150):: all_line,all_line_tmp
logical:: vasp_5,red
units="bohr"
nsteps=40
!Get the two filenames
write(*,*) "Name of file:"
read(*,*) filename1
write(*,*) "Is it VASP5 format? (T for true, F for false)"
read(*,*) vasp_5

!Get nat and allocate some arrays
open(unit=2,file=trim(filename1))
read(2,*) tmp_ch

!Read the first structure to initialize stuff
!read(all_line(1:n),*) typchar(:)
read(2,*,end=99) scaling
read(2,*,end=99) latvec(:,1)
read(2,*,end=99) latvec(:,2)
read(2,*,end=99) latvec(:,3)
if(vasp_5) read(2,'(a150)',end=99) all_line_tmp
read(2,'(a150)',end=99) all_line
all_line=" "//trim(all_line)
!write(*,*) all_line
ntype=0
do i=1,150-1
if(all_line(i:i)==" ".and.all_line(i+1:i+1).ne." ") then
ntype=ntype+1
endif
enddo
!write(*,*) ntype
allocate(typchar(ntype),nitype(ntype))

read(all_line,*) nitype(:)
if(.not.vasp_5) then
write(*,*) "Enter elements symbols "
read(*,*)typchar(:)
endif
if(vasp_5) read(all_line_tmp,*) typchar(:)
!write(*,*) typchar(:)

read(2,'(a150)',end=99)all_line
!write(*,*) all_line
n = len_trim(all_line)
k = index(all_line(1:n),"S") !In selective dynamics
  if(k.ne.0) then  !ignore this line
  read(2,'(a150)',end=99)all_line
  endif 
all_line=adjustl(all_line)
all_line=trim(all_line)
all_line=all_line(1:1)
if(trim(all_line).eq."D".or.trim(all_line).eq."d") then
  red=.true.
elseif(trim(all_line).eq."C".or.trim(all_line).eq."c") then
  red=.false.
else
  stop "Coordinates must be either Direct or Cartesian"
endif
!Allocate stuff
nat=0
do i=1,ntype
nat=nat+nitype(i)
enddo
!write(*,*) nat
allocate(xred(3,nat),itype(nat),typ_at(nat))
iat=0
do i=1,ntype
 do j=1,nitype(i)
 iat=iat+1
 itype(iat)=i
 typ_at(iat)=typchar(i)
 enddo
enddo
!Read initial positions
do iat=1,nat
   read(2,*,end=99) xred(:,iat)
enddo
istep=0
  latvec=latvec*scaling
  if(.not.red) then
   allocate(xcart(3,nat))
   xcart=xred
   call rxyz_cart2int(latvec,xred,xcart,nat)
  endif
  energy=0.d0
  pressure=0.d0
  printval1=0.d0
  printval2=0.d0
  filename=trim(filename1)//".ascii"
  call write_atomic_file_ascii(filename,nat,units,xred,latvec,energy,pressure,printval1,printval2,typ_at)   


99 continue


endprogram




!**********************************************************************************************

subroutine read_atomic_file_ascii(filename,nat,units,xred,latvec,printval1,printval2,typ_at)
!subroutine read_atomic_file_ascii(filename,nat,units,xred,latvec,acell,rprim,printval1,printval2)
!This routine will read the file "filename" ins ascii file format
!The output unit will always be in atomic units (bohr)
!So if units==angstroem, the file will be converted
!   if units==bohr, the positions will not be changed
implicit none
integer:: nat,natin,iat,ierror
character(140):: filename,units
real(8):: pos(3,nat),xred(3,nat),latvec(3,3),dproj(6)
real(8):: angbohr, enthalpy_at,printval1,printval2
character(2)::typ_at(nat)

if(trim(units)=="angstroem") then
  angbohr=1.889725989d0
elseif(trim(units)=="bohr") then
  angbohr=1.d0
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
if(natin.ne.nat) stop "Number of atoms not consistent with abinit input"
read(46,*) dproj(1:3)
read(46,*) dproj(4:6)
do iat=1,nat
  read(46,*) pos(:,iat),typ_at(iat)
enddo
close(46)
dproj=dproj*angbohr
pos=pos*angbohr
call dproj2latvec(dproj,latvec)

!acell(1)=dproj(1)
!acell(2)=sqrt(dproj(2)**2+dproj(3)**2)
!acell(3)=sqrt(dproj(4)**2+dproj(5)**2+dproj(6)**2)
!
!rprim(1,1)=1.d0
!rprim(2,1)=0.d0
!rprim(3,1)=0.d0
!
!rprim(1,2)=dproj(2)/acell(2)
!rprim(2,2)=dproj(3)/acell(2)
!rprim(3,2)=0.d0
!
!rprim(1,3)=dproj(4)/acell(3)
!rprim(2,3)=dproj(5)/acell(3)
!rprim(3,3)=dproj(6)/acell(3)


call rxyz_cart2int(latvec,xred,pos,nat)

end subroutine


!************************************************************************************

subroutine write_atomic_file_ascii(filename,nat,units,xred,latvec0,energy,pressure,printval1,printval2,typ_at)
!This routine will write the file "filename" in ascii file format
!The input unit will always be in atomic units (bohr, hartree), but the output can be specified by the vaule in "units"
!So if units==angstroem, the file will be converted to angstroem
!   if units==bohr, the positions will not be changed
implicit none
integer:: nat,natin,iat
character(140):: filename,units
real(8):: pos(3,nat),xred(3,nat),latvec(3,3),latvec0(3,3),dproj(6),rotmat(3,3),v(3,3),ucvol
real(8):: angbohr,hartree2ev,in_GPA,in_ang3,int_press
real(8):: energy, etotal, enthalpy, enthalpy_at,pressure,printval1,printval2
character(2)::typ_at(nat)

if(trim(units)=="angstroem") then
  angbohr=1.d0/1.889725989d0
  hartree2ev=27.211396132d0
  in_GPA=29421.033d0
  in_ang3=0.148184743d0
  int_press=160.217646200d0
elseif(trim(units)=="bohr") then
  angbohr=1.d0
  hartree2ev=1.d0
  in_GPA=1.d0
  in_ang3=1.d0
  int_press=1.d0
else
  stop "Wrong file unit format"
endif

!latvec(:,1)=acell(1)*rprim(:,1)
!latvec(:,2)=acell(2)*rprim(:,2)
!latvec(:,3)=acell(3)*rprim(:,3)

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
do iat=1,nat
  write(46,'(3(1x,es25.15),2x,a2)') pos(:,iat),typ_at(iat)
enddo
close(46)
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
 if(latvec(3,3).lt.0.d0)then
 write(*,*) latvec(:,1)
 write(*,*) latvec(:,2)
 write(*,*) latvec(:,3)
 stop "Error in orientation of the cell"
 endif
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

 subroutine dot_p(a,b,dotp)
 !a very simple implementation of the dot product
 implicit none
 real(8)::a(3),b(3)
 real(8)::dotp(3)
 integer::i
 dotp=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
 end subroutine

!************************************************************************************

 subroutine rotation(rotmat,angle,axe)
 !This subroutine will calculate the rotational matrix rotmat for a
 !3-dim vector around an axis 'axe' by the angle 'angle'.
 implicit none
 real(8),INTENT(IN) :: angle
 real(8),INTENT(IN) :: axe(3)
 real(8):: rotator(3,3)
 real(8):: rotmat(3,3)

 !Define Rotation Matrix
 rotator(1,1)=dcos(angle)+(axe(1)**2)*(1.d0-dcos(angle))
 rotator(1,2)=axe(1)*axe(2)*(1.d0-dcos(angle))-axe(3)*dsin(angle)
 rotator(1,3)=axe(1)*axe(3)*(1.d0-dcos(angle))+axe(2)*dsin(angle)

 rotator(2,1)=axe(2)*axe(1)*(1.d0-dcos(angle))+axe(3)*dsin(angle)
 rotator(2,2)=dcos(angle)+(axe(2)**2)*(1.d0-dcos(angle))
 rotator(2,3)=axe(2)*axe(3)*(1.d0-dcos(angle))-axe(1)*dsin(angle)

 rotator(3,1)=axe(3)*axe(1)*(1.d0-dcos(angle))-axe(2)*dsin(angle)
 rotator(3,2)=axe(3)*axe(2)*(1.d0-dcos(angle))+axe(1)*dsin(angle)
 rotator(3,3)=dcos(angle)+(axe(3)**2)*(1.d0-dcos(angle))
 rotmat(:,:)=rotator(:,:)

 !do i=1,3
 !   vector2(i)=rotator(i,1)*vector(1)+rotator(i,2)*vector(2)+rotator(i,3)*vector(3)
 !enddo
 !vector(:)=vector2(:)
 end subroutine rotation

!************************************************************************************
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
 a=mat
 div=(a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)+&
 &a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)) 
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
! !General n*n matrix 
! matinv=mat
! allocate(WORK(n))
! call  DGETRF( n, n, matinv, n, IPIV, INFO )
! if (info.ne.0) stop "Error in DGETRF"
! LDWORK=-1
! call  DGETRI( n, matinv, n, IPIV, WORK,LDWORK , INFO )
! LDWORK=WORK(1)
! deallocate(WORK)
! allocate(WORK(LDWORK))
! call  DGETRI( n, matinv, n, IPIV, WORK,LDWORK , INFO )
! if (info.ne.0) stop "Error in DGETRI"
 end subroutine

!************************************************************************************
