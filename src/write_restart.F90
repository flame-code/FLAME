subroutine winter(parini,nat,units,ent_pos,e_pos,pos_red,pos_latvec,pos_fcart,pos_strten,nlminx,nlmin,npminx,& 
   &ent_arr,e_arr,ct_arr,spg_arr,spgtol_arr,dos_arr,pl_arr,lat_arr,f_arr,str_arr,fp_arr,fp_len,ent_delta,fp_delta,& 
   &eref,ediff,ekinetic,ekinetic_max,dt,nsoften,char_type,ntypat,typat,fixat,fixlat,pressure)
  use mod_parini, only: typ_parini
  use yaml_output
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
    call yaml_mapping_open('wrote poscur.ascii for RESTART',flow=.true.)
    call yaml_map('ent_pos',ent_pos)
    call yaml_mapping_close()
     !write(*,*) ' wrote poscur.ascii for RESTART',ent_pos
     

     
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
     call yaml_map('wrote poslow files',nlmin)
     !write(*,*) ' wrote poslow files',nlmin
     !write(*,*) ' wrote earr.dat for  RESTART'
     call yaml_comment('wrote earr.dat for  RESTART',hfill='~')
     close(12)
     
     call wtioput(ediff,ekinetic,ekinetic_max,parini%nsoften_minhopp)
     !write(*,*) ' wrote ioput for  RESTART'
     call yaml_comment('wrote ioput for  RESTART',hfill='~')

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
     call yaml_comment('Binary files {poslow,fcart,strten,latvec,fp}.bin written for RESTART.')
     !write(*,*) ' wrote binary files poslow.bin, fcart.bin, strten.bin, latvec.bin and fp.bin for RESTART'


!!     call wtpos(npminx,nlminx,nlmin,npmin,poslocmin,latlocmin,earr,elocmin,char_type,&
!!          &ntypat,typat,fixat,fixlat,pressure,units,nat,eref)

END SUBROUTINE winter

!**********************************************************************************************

subroutine wtioput(ediff,ekinetic,ekinetic_max,nsoften)
  implicit none
  integer:: nsoften
  real(8):: ediff, ekinetic,ekinetic_max
  open(unit=11,file='ioput',status='unknown')
  write(11,'(3(1x,1pe24.17)1x,a)') ediff,ekinetic,ekinetic_max,' ediff, temperature, maximal temperature'
  close(11)
END SUBROUTINE wtioput
