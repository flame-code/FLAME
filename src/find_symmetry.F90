subroutine find_symmetry(parini,nat,xred,latvec,typat,tolmin,tolmax,ntol,tolcur,spgcur)
use mod_parini, only: typ_parini
use tb_lj_params
use void_lj_params, only: nat_atoms
use yaml_output
implicit none
type(typ_parini), intent(in):: parini
real(8):: xred(3,nat)
real(8):: latvec(3,3)
real(8):: tol,tolfact,tolmax,tolmin,tolcur,tolcur2
integer:: itol,ntol,spgcur,spg,nat
integer:: spgcur2,spg2
integer:: typat(nat),nat_sym


spgcur=0
spgcur2=0
tolcur2=0.d0
tolfact=(tolmax/tolmin)**(1.d0/real(ntol,8))

if(trim(parini%potential_potential)=="lenosky_tb_lj") then
  nat_sym=n_silicon+n_h
elseif(parini%voids) then
  nat_sym=nat_atoms
else
  nat_sym=nat
endif

#if defined(SPGLIB)
call yaml_sequence_open('SPGLIB space group iterations')
do itol=0,ntol
tol=tolmin*tolfact**itol
!call check_symmetry(nat,xred,latvec,typat,spg,tol)
!if(spg.gt.spgcur) then
!   spgcur=spg
!   tolcur=tol
!endif
!write(*,'(a,i5,es15.7)') " # FINDSYM SPG, TOL      : ",spg,tol

call get_spg(nat_sym,xred(:,1:nat_sym),latvec,typat(1:nat_sym),tol,spg2)

if(spg2.gt.spgcur2) then
   spgcur2=spg2
   tolcur2=tol
endif
!write(*,'(a,i5,es15.7)') " # SPGLIB  SPG, TOL      : ",spg2,tol
enddo
call yaml_sequence_close()
#endif

!write(*,'(a,i5,es15.7)') " # FINAL FINDSYM SPG, TOL: ",spgcur,tolcur
call yaml_mapping_open('FINAL SPGLIB',flow=.true.)
call yaml_map('SPG',spgcur2,fmt='(i8)')
call yaml_map('TOL',tolcur2,fmt='(es15.7)')
call yaml_mapping_close()
!write(*,'(a,i5,es15.7)') " # FINAL SPGLIB  SPG, TOL: ",spgcur2,tolcur2
spgcur=spgcur2
tolcur=tolcur2
end subroutine find_symmetry



subroutine check_symmetry(parini,nat,xred,latvec,typat,spg,tol)
use mod_parini, only: typ_parini
use tb_lj_params
implicit none
type(typ_parini), intent(in):: parini
real(8):: xred(3,nat)
real(8):: latvec(3,3)
real(8):: angbohr,tol
integer:: nat,nat_sym,iat, jat,count,n,k
integer:: typat(nat),spg
character(40):: type_format,tmp_ch
character(150):: all_line
parameter(angbohr=1.889725989d0)

spg=0

if(trim(parini%potential_potential)=="lenosky_tb_lj") then
  nat_sym=n_silicon+n_h
else
  nat_sym=nat
endif

!Prepare input file
open(unit=45,file="findsym.in")
write(45,'(a)') "Input for FINDSYM, from filename "
write(45,'(es15.7,a)') tol," Tolerance"
write(45,'(a)') "1    form of lattice parameters: to be entered as vectors"
write(45,'(3(1x,f15.8),a)') latvec(:,1), " Lattice vector 1" 
write(45,'(3(1x,f15.8),a)') latvec(:,2), " Lattice vector 2"
write(45,'(3(1x,f15.8),a)') latvec(:,3), " Lattice vector 3"
write(45,'(a)') "2    form of primitive lattice vectors"
write(45,'(a)') "P    unknown centering"
write(45,'(i5,a)') nat_sym, "    number of atoms"
!call latvec2dist_ang(latvec,angdeg)
!write(*,*) "a,b,c,alpha,beta,gamma"
!write(*,'(6(1x,es25.16))') angdeg



write(type_format,'(a,i0,a)') "(",nat_sym,"(1x,i0),a)"
write(45,type_format) typat(1:nat_sym)  ,"  Atom kinds"

do iat=1,nat_sym
write(45,'(3(1x,es25.16),a,i4)') xred(:,iat), "  Reduced Coordinates ", iat
enddo 
close(45)
!Run findsym, be sure to store the output in "findsym.out"
 call system("rm -f findsym.out")
 call system("./run_findsym.sh")
! call system("sleep 1")

!Now start reading the results
!First check if there was an error during the run
open(unit=45,file="findsym.out")
do while(.true.)
 read(45,'(a150)',end=99) all_line
!!write(*,*) all_line
 n = len_trim(all_line)

 k = index(all_line(1:n),"rror")
 if(k==0) k = index(all_line(1:n),"bombed")
  if(k.ne.0) then
    spg=-1
    close(45)
    return
  endif
enddo
99 continue
close(45)

!If i get here, this means the file is fine
open(unit=45,file="findsym.out")
do while(.true.)
 read(45,'(a150)',end=88) all_line
!!write(*,*) all_line
 n = len_trim(all_line)

 k = index(all_line(1:n),"_symmetry_Int_Tables_number")
  if(k.ne.0) then
    read(all_line,*) tmp_ch, spg
    goto 88
  endif
enddo
88 continue
close(45)

if(spg.le.0) stop "Something wrong with symmetry finder"
end subroutine





