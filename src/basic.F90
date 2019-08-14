!*****************************************************************************************
subroutine elim_white_space(string)
    implicit none
    character(256), intent(inout):: string
    !local variables 
    integer:: i, j
    do i=1,256
        if(string(1:1)==' ') then
            do j=1,256-i
                string(j:j)=string(j+1:j+1)
            enddo
            string(256-i+1:256-i+1)=''
        else
            exit
        endif
    enddo
end subroutine elim_white_space
!*****************************************************************************************
function delta_kronecker(i,j) result(delta)
    implicit none
    integer, intent(in):: i, j
    !local variables
    real(8):: delta
    if(i==j) then
        delta=1.d0
    else
        delta=0.d0
    endif
end function delta_kronecker
!*****************************************************************************************
!This subroutine will transform back all atoms into the periodic cell
!defined by the 3 lattice vectors in latvec=[v1.v2.v3]
subroutine backtocell_alborz(nat,latvec,rxyz_red)
    implicit none
    integer, intent(in):: nat
    real(8), intent(in):: latvec(3,3)
    real(8), intent(inout):: rxyz_red(3,nat)
    !local variables 
    real(8):: vol
    integer:: iat
    !First check if the volume is positive
    call getvol_alborz(latvec,vol)
    if(vol.le.0.d0) stop "Negative volume during backtocell_alborz"
    do iat=1,nat
        rxyz_red(1,iat)=modulo(modulo(rxyz_red(1,iat),1.d0),1.d0)
        rxyz_red(2,iat)=modulo(modulo(rxyz_red(2,iat),1.d0),1.d0)
        rxyz_red(3,iat)=modulo(modulo(rxyz_red(3,iat),1.d0),1.d0)
    enddo
end subroutine backtocell_alborz
!*****************************************************************************************
subroutine getvol_alborz(cellvec,vol)
    implicit none
    real(8), intent(in):: cellvec(3,3)
    real(8), intent(out):: vol
    vol=cellvec(1,1)*cellvec(2,2)*cellvec(3,3)-cellvec(1,1)*cellvec(2,3)*cellvec(3,2)- &
        cellvec(1,2)*cellvec(2,1)*cellvec(3,3)+cellvec(1,2)*cellvec(2,3)*cellvec(3,1)+ &
        cellvec(1,3)*cellvec(2,1)*cellvec(3,2)-cellvec(1,3)*cellvec(2,2)*cellvec(3,1)
    if(vol<0.d0) stop 'ERROR: Negative volume!'
end subroutine getvol_alborz
!*****************************************************************************************
!This routine computes the distance of the reduced coordinates xred_1 and xred_2 and applies 
!periodic boundary conditions to them and computes the squared distance
subroutine pbc_distance1_alborz(cellvec,xred_1,xred_2,distance2,dxyz)
    implicit none
    real(8), intent(in):: cellvec(3,3), xred_1(3), xred_2(3)
    real(8), intent(inout):: distance2, dxyz(3)
    !local variables 
    integer:: i
    real(8):: diff(3)
    diff(1:3)=xred_2(1:3)-xred_1(1:3)
    do i=1,3
        if(.not. diff(i)>-0.5d0) then
            diff(i)=diff(i)+1.d0
        elseif(diff(i)>0.5d0) then
            diff(i)=diff(i)-1.d0
        endif
    enddo
    !dxyz=matmul(cellvec,diff)
    dxyz(1)=cellvec(1,1)*diff(1)+cellvec(1,2)*diff(2)+cellvec(1,3)*diff(3)
    dxyz(2)=cellvec(2,1)*diff(1)+cellvec(2,2)*diff(2)+cellvec(2,3)*diff(3)
    dxyz(3)=cellvec(3,1)*diff(1)+cellvec(3,2)*diff(2)+cellvec(3,3)*diff(3)
    distance2=dxyz(1)*dxyz(1)+dxyz(2)*dxyz(2)+dxyz(3)*dxyz(3)
end subroutine pbc_distance1_alborz
!*****************************************************************************************
!This subroutine will return how many periodic expansions for each lattice vector 
!direction are necessary for the periodic boundary conditions
!with for the given rcut. nec1,nec2,nec3 for latvec(:,1),latvec(:,2),latvec(:,3)
subroutine n_rep_dim_alborz(cellvec,rcut,nec1,nec2,nec3)
    implicit none
    real(8), intent(in):: cellvec(3,3), rcut
    integer, intent(out):: nec1, nec2, nec3
    !local variables 
    real(8):: zero(3), dist(3)
    real(8):: vn(3,3) !normalized normal vector of the planes form by cell vectors
    integer:: i
    nec1=0 ; nec2=0 ; nec3=0
    call nveclatvec_alborz(cellvec,vn)
    zero=(/0.d0,0.d0,0.d0/)
    do i=1,3
        call dist2plane_alborz(cellvec(1,mod(i+1,3)+1),vn(1,i),zero,dist(i))
        !write(*,*) "rcut",i,rcut, dist
    enddo
    nec1=int(rcut/dist(2))+1
    nec2=int(rcut/dist(3))+1
    nec3=int(rcut/dist(1))+1
end subroutine n_rep_dim_alborz
!*****************************************************************************************
!Will calculate the normalized normal vector to the 3 planes of the cell
subroutine nveclatvec_alborz(cellvec,vn)
    implicit none
    real(8), intent(in) :: cellvec(3,3)
    real(8), intent(out):: vn(3,3)
    !local variables 
    real(8):: v(3), vnorm
    integer:: i
    do i=1,3
        call cross_product_alborz(cellvec(1,i),cellvec(1,mod(i,3)+1),v)
        vnorm=sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
        vn(1,i)=v(1)/vnorm
        vn(2,i)=v(2)/vnorm
        vn(3,i)=v(3)/vnorm
    enddo
end subroutine nveclatvec_alborz
!*****************************************************************************************
!This subroutine will calculate  the distance between a plane and a point in space
!The point is 'r1', the normalized normal vector of the plane is 'nvec', 
!'r0' is an arbitrary point on the plane !and the output is the distance 'dist'  
subroutine dist2plane_alborz(r1,vn,r0,dist)
    implicit none
    real(8), intent(in):: r1(3), vn(3), r0(3)
    real(8), intent(out):: dist
    !local variables 
    real(8):: rd(3)
    rd(1)=r1(1)-r0(1)
    rd(2)=r1(2)-r0(2)
    rd(3)=r1(3)-r0(3)
    dist=abs(rd(1)*vn(1)+rd(2)*vn(2)+rd(3)*vn(3))
end subroutine dist2plane_alborz
!*****************************************************************************************
!This routine will write the file "filename" in ascii file format
!The input unit will always be in atomic units (bohr, hartree), but the output can be specified by the vaule in "units"
!So if units==angstroem, the file will be converted to angstroem
!   if units==bohr, the positions will not be changed
subroutine write_atomic_file_ascii_alborz(filename,nat,xred,latvec0,energy,pressure,printval1,printval2,kinds)
implicit none
integer:: nat,natin,iat
character(40):: filename,units
real(8):: pos(3,nat),xred(3,nat),latvec(3,3),latvec0(3,3),dproj(6),rotmat(3,3),v(3,3),ucvol
real(8):: angbohr,hartree2ev,in_GPA,in_ang3,int_press
real(8):: energy, etotal, enthalpy, enthalpy_at,pressure,printval1,printval2
integer:: Kinds(nat)
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
  angbohr=1.d0
  hartree2ev=1.d0
  in_GPA=1.d0
  in_ang3=1.d0
  int_press=1.d0
endif

!latvec(:,1)=acell(1)*rprim(:,1)
!latvec(:,2)=acell(2)*rprim(:,2)
!latvec(:,3)=acell(3)*rprim(:,3)

latvec=latvec0

do iat=1,nat
   pos(:,iat)=matmul(latvec,xred(:,iat))
enddo

call latvec2dproj_alborz(dproj,latvec,rotmat,pos,nat)

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
  write(46,'(i4,1x,es25.15,1x,es25.15,1x,a20,es15.7,a14,es15.7,a13,es15.7)')&
  & nat,printval1,printval2,"         energy(Ha)=",etotal," enthalpy(Ha)=",enthalpy,&
  &" ucvol(bhr3)=",ucvol
write(46,*) dproj(1:3)
write(46,*) dproj(4:6)
do iat=1,nat
if(kinds(iat)==1)  write(46,'(3(1x,es25.15),2x,a2)') pos(:,iat),"A "
if(kinds(iat)==2)  write(46,'(3(1x,es25.15),2x,a2)') pos(:,iat),"B "
enddo
close(46)
end subroutine write_atomic_file_ascii_alborz
!*****************************************************************************************
!This subroutine will convert the distance and projective representation of 
!a periodic cell (dxx,dyx,dyy,dzx,dzy,dzz) into a 
!lattice vektor format (vec1(:,1),vec2(:,2),vec3(:,3)) with dxx oriented into x direction
subroutine dproj2latvec_alborz(dproj,cellvec)
    implicit none
    real(8), intent(in):: dproj(6)
    real(8), intent(out):: cellvec(3,3)
    !local variables
    cellvec(1,1)=dproj(1) ; cellvec(1,2)=dproj(2) ; cellvec(1,3)=dproj(4)
    cellvec(2,1)=0.d0     ; cellvec(2,2)=dproj(3) ; cellvec(2,3)=dproj(5)
    cellvec(3,1)=0.d0     ; cellvec(3,2)=0.d0     ; cellvec(3,3)=dproj(6) 
end subroutine dproj2latvec_alborz
!*****************************************************************************************
!This subroutine will convert the lattice vector representation of thei
!periodic cell (vec1,vec2,vec3) into the projective representation (dxx,dyx,dyy,dzx,dzy,dzz)
!The cell will thus be rotated. The rotational matrix is stored in rotmat as an operator rotmat
!and the atomic position rxyz are transformed into the new coordination sizstem as well
subroutine latvec2dproj_alborz(dproj,latvec,rotmat,rxyz,nat)
    implicit none
    integer, intent(in):: nat
    real(8),intent(inout):: dproj(6), latvec(3,3), rotmat(3,3), rxyz(3,nat)
    !local variables
    real(8):: tempvec(3), rotmat1(3,3), rotmat2(3,3), crossp(3), alpha, latvect(3,3)
    real(8):: eps,axe(3), anrm,rxyzt(3), rotmatt(3,3)
    integer:: iat,i
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
    call cross_product_alborz(latvec(:,1),tempvec,crossp)
    if (abs(crossp(1)).lt.eps*1.d-1 .and. abs(crossp(2)).lt.eps*1.d-1 .and. abs(crossp(3)).lt.eps*1.d-1) goto 1001 !no rotation needed
    anrm=sqrt(crossp(1)*crossp(1)+crossp(2)*crossp(2)+crossp(3)*crossp(3))
    axe(:)=crossp(:)/anrm
    alpha=dacos(dot_product(tempvec,latvec(:,1))/dproj(1))
    call rotation_alborz(alpha,axe,rotmat1)
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
    call cross_product_alborz(tempvec,(/0.d0,1.d0,0.d0/),axe)
    if (abs(axe(1)).lt.eps*1.d-1 .and. abs(axe(2)).lt.eps*1.d-1 .and. abs(axe(3)).lt.eps*1.d-1 .and. tempvec(2).gt.0.d0) goto 1002 !no rotation needed
    anrm=sqrt(axe(1)*axe(1)+axe(2)*axe(2)+axe(3)*axe(3))
    axe=axe/anrm
    call cross_product_alborz(axe,latvec(:,2),crossp)
    if (abs(crossp(1)).lt.eps*1.d-1 .and. abs(crossp(2)).lt.eps*1.d-1 .and. crossp(3).gt.0.d0) goto 1002 !no rotation needed
    anrm=sqrt(tempvec(2)*tempvec(2)+tempvec(3)*tempvec(3))
    alpha=dot_product((/0.d0,1.d0,0.d0/),tempvec(:))/anrm
    alpha=dacos(alpha)
    call rotation_alborz(alpha,axe,rotmat2)
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
end subroutine latvec2dproj_alborz
!*****************************************************************************************
!a very simple implementation of the cross product
subroutine cross_product_alborz(a,b,c)
    implicit none
    real(8), intent(in):: a(3), b(3)
    real(8), intent(out):: c(3)
    c(1)=a(2)*b(3)-a(3)*b(2)
    c(2)=a(3)*b(1)-a(1)*b(3)
    c(3)=a(1)*b(2)-a(2)*b(1)
end subroutine cross_product_alborz
!*****************************************************************************************
!This subroutine will calculate the rotational matrix rotmat for a
!3-dim vector around an axis 'axe' by the angle 'angle'.
subroutine rotation_alborz(angle,axe,rotmat)
    implicit none
    real(8), intent(in):: angle
    real(8), intent(in):: axe(3)
    real(8), intent(out):: rotmat(3,3)
    !local variables

    !Define Rotation Matrix
    rotmat(1,1)=cos(angle)+(axe(1)**2)*(1.d0-cos(angle))
    rotmat(1,2)=axe(1)*axe(2)*(1.d0-cos(angle))-axe(3)*dsin(angle)
    rotmat(1,3)=axe(1)*axe(3)*(1.d0-cos(angle))+axe(2)*dsin(angle)

    rotmat(2,1)=axe(2)*axe(1)*(1.d0-cos(angle))+axe(3)*dsin(angle)
    rotmat(2,2)=cos(angle)+(axe(2)**2)*(1.d0-cos(angle))
    rotmat(2,3)=axe(2)*axe(3)*(1.d0-cos(angle))-axe(1)*dsin(angle)

    rotmat(3,1)=axe(3)*axe(1)*(1.d0-cos(angle))-axe(2)*dsin(angle)
    rotmat(3,2)=axe(3)*axe(2)*(1.d0-cos(angle))+axe(1)*dsin(angle)
    rotmat(3,3)=cos(angle)+(axe(3)**2)*(1.d0-cos(angle))

    !do i=1,3
    !   vector2(i)=rotator(i,1)*vector(1)+rotator(i,2)*vector(2)+rotator(i,3)*vector(3)
    !enddo
    !vector(:)=vector2(:)
end subroutine rotation_alborz
!*****************************************************************************************
!This subrtouine will transform theforces initially in the cartesian system into the 
!internal coordinates with respect to the cell vectors provided in the cv
subroutine fxyz_cart2int_alborz(nat,v_cart,cv,v_int)
    implicit none
    integer, intent(in):: nat
    real(8), intent(in):: v_cart(3,nat), cv(3,3)
    real(8), intent(out):: v_int(3,nat)
    !local variables 
    integer:: iat
    do iat=1,nat
        v_int(1,iat)=cv(1,1)*v_cart(1,iat)+cv(2,1)*v_cart(2,iat)+cv(3,1)*v_cart(3,iat)
        v_int(2,iat)=cv(1,2)*v_cart(1,iat)+cv(2,2)*v_cart(2,iat)+cv(3,2)*v_cart(3,iat)
        v_int(3,iat)=cv(1,3)*v_cart(1,iat)+cv(2,3)*v_cart(2,iat)+cv(3,3)*v_cart(3,iat)
    enddo
end subroutine fxyz_cart2int_alborz
!*****************************************************************************************
subroutine fxyz_red2cart(nat,fint,cv,fcart)
    implicit none
    integer, intent(in):: nat
    real(8), intent(in):: fint(3,nat), cv(3,3)
    real(8), intent(out):: fcart(3,nat)
    !local variables 
    integer:: iat
    real(8):: cvinv(3,3)
    call invertmat_alborz(cv,cvinv)
    do iat=1,nat
        fcart(1,iat)=cvinv(1,1)*fint(1,iat)+cvinv(2,1)*fint(2,iat)+cvinv(3,1)*fint(3,iat)
        fcart(2,iat)=cvinv(1,2)*fint(1,iat)+cvinv(2,2)*fint(2,iat)+cvinv(3,2)*fint(3,iat)
        fcart(3,iat)=cvinv(1,3)*fint(1,iat)+cvinv(2,3)*fint(2,iat)+cvinv(3,3)*fint(3,iat)
    enddo
end subroutine fxyz_red2cart
!*****************************************************************************************
subroutine count_words(str,n)
    implicit none
    character(*), intent(in):: str
    integer, intent(out):: n
    !local variables 
    integer:: ind, k
    character(50):: word
    character(256):: str2
    str2=str
    n=0
    word=''
    do
        str2=adjustl(str2)
        read(str2,*,iostat=k) word
        !write(*,'(i3,1x,a,i5)') n,trim(word),len_trim(word)
        if(len_trim(word)==0) return
        ind=index(str2,trim(word))
        !write(*,'(2i3)') ind,ind+len_trim(word)-1
        str2(ind:ind+len_trim(word)-1)=' '
        !write(*,*) 'updated str2',trim(str2)
        n=n+1
        word=''
        !write(*,*)
    enddo
end subroutine count_words
!*****************************************************************************************
subroutine count_substring(str1,str2,n)
    implicit none
    character(*), intent(in) :: str1, str2
    integer, intent(out):: n
    !local variables 
    integer:: m, ind
    stop 'ERROR: this subroutine is not tested'
    n=0
    if(len(str2)==0) return
    m=1
    do 
        ind=index(str1(m:),str2)
        if(ind==0) return
        n=n+1
        m=m+ind+len(str2)
    enddo
end subroutine count_substring
!*****************************************************************************************
