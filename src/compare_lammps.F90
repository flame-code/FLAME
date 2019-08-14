subroutine compare_lammps(parini,parres)
use global
use interface_code
use defs_basis
use mod_parini, only: typ_parini
! Main program to test potential subroutines
implicit none
type(typ_parini), intent(in):: parini
type(typ_parini), intent(inout):: parres
integer:: k,l,n,iat,iprec
real(8):: latvec(3,3),xred(3,parini%nat),xcart(3,parini%nat),f_lammps(3,parini%nat),f(3,parini%nat),e_lammps,e,tmp_r,tmp_i,tilts(6),latvec_in(3,3),strten(6),latvec_box(3,3)
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
             do iat=1,parini%nat
               read(7,*) tmp_i,tmp_i,xcart(:,iat),f_lammps(:,iat) 
!               read(7,*) tmp_i,tmp_i,xred(:,iat),f_lammps(:,iat) 
             enddo
             exit
           endif
       enddo
!!       do iat=1,nat
!!          xcart(:,iat)=matmul(latvec_box,xred(:,iat))
!!       enddo
       call rxyz_cart2int(latvec,xred,xcart,parini%nat)
       latvec_in=latvec/Bohr_ang
       call get_energyandforces_single(parini,parres,latvec_in,xred,f,strten,e,iprec,getwfk)
       f=f*Ha_eV/Bohr_Ang
       e=e*Ha_eV
       write(*,*) tmp_i,e_lammps
       if(any(abs(f-f_lammps).gt.1.d-10)) then
          write(*,*) "ENERGIES"
          write(*,'(2es25.15)') e,e_lammps
          write(*,*) "ERROR in forces"
          do iat=1,parini%nat
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

