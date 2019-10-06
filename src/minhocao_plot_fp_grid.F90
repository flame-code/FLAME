subroutine plot_fp_grid(parini,nlminx,nlmin,nat,fp_len,fp_arr,lat_arr,pl_arr)
use mod_parini, only: typ_parini
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
!!  call init_fp(parini,fp_len)
!!  allocate(fp_arr(fp_len,nlminx))
!do kk=1,nlmin
!   call get_fp(parini,fp_len,pl_arr(:,:,kk),lat_arr(:,:,kk),fp_arr(:,kk))
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
!!  call init_fp(parini,fp_len)
!!  allocate(fp_arr(fp_len,nlminx))
!!do kk=1,nlmin
!!   call get_fp(parini,fp_len,pl_arr(:,:,kk),lat_arr(:,:,kk),fp_arr(:,kk))
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
!!  call init_fp(parini,fp_len)
!!  allocate(fp_arr(fp_len,nlminx))
!!do kk=1,nlmin
!!   call get_fp(parini,fp_len,pl_arr(:,:,kk),lat_arr(:,:,kk),fp_arr(:,kk))
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
