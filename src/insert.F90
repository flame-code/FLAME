subroutine insert(nlminx,nlmin,fp_len,nat,k_e_wpos,e_wpos,ent_wpos,fp_wpos,wpos_red,&
  &wpos_latvec,wpos_fcart,wpos_strten,spg_wpos,spgtol_wpos,fdos_wpos,&
  &e_arr,ent_arr,fp_arr,pl_arr,lat_arr,f_arr,str_arr,spg_arr,spgtol_arr,dos_arr,ct_arr)
  ! inserts the energy re_wpos at position k_e_wpos and shifts up all other energies
!  implicit real*8 (a-h,o-z)
  use yaml_output
  implicit none
  integer:: fp_len,ct_arr(nlminx),spg_arr(nlminx),nat,iat,spg_wpos
  integer:: k, nlmin, k_e_wpos, nlminx,i
  real(8):: e_wpos, ent_wpos, wpos_red(3,nat),wpos_latvec(3,3),spgtol_wpos,fdos_wpos,fp_wpos(fp_len)
  real(8):: e_arr(nlminx),ent_arr(nlminx),fp_arr(fp_len,nlminx),pl_arr(3,nat,nlminx),f_arr(3,nat,nlminx)
  real(8):: lat_arr(3,3,nlminx),spgtol_arr(nlminx),dos_arr(nlminx),fp1(fp_len),fp2(fp_len),str_arr(6,nlminx)
  real(8):: wpos_fcart(3,nat),wpos_strten(6)
              do k=nlmin-1,k_e_wpos+1,-1
                 e_arr(k+1)=e_arr(k)
                 ent_arr(k+1)=ent_arr(k)
                 do i=1,fp_len
                 fp_arr(i,k+1)=fp_arr(i,k)
                 enddo
                 do iat=1,nat
                 pl_arr(1,iat,k+1)=pl_arr(1,iat,k)
                 pl_arr(2,iat,k+1)=pl_arr(2,iat,k)
                 pl_arr(3,iat,k+1)=pl_arr(3,iat,k)
                 f_arr(1,iat,k+1)=f_arr(1,iat,k)
                 f_arr(2,iat,k+1)=f_arr(2,iat,k)
                 f_arr(3,iat,k+1)=f_arr(3,iat,k)
                 enddo
                 str_arr(:,k+1)=str_arr(:,k)
                 do i=1,3
                 lat_arr(1,i,k+1)=lat_arr(1,i,k)
                 lat_arr(2,i,k+1)=lat_arr(2,i,k)
                 lat_arr(3,i,k+1)=lat_arr(3,i,k)
                 enddo
                 spg_arr(k+1)=spg_arr(k)
                 spgtol_arr(k+1)=spgtol_arr(k)
                 dos_arr(k+1)=dos_arr(k)
                 ct_arr(k+1)=ct_arr(k)
              enddo
              e_arr(k_e_wpos+1)=e_wpos
              ent_arr(k_e_wpos+1)=ent_wpos
              ct_arr(k_e_wpos+1)=1
                 do i=1,fp_len
                 fp_arr(i,k_e_wpos+1)=fp_wpos(i)
                 enddo
                 do iat=1,nat
                 pl_arr(1,iat,k_e_wpos+1)=wpos_red(1,iat)
                 pl_arr(2,iat,k_e_wpos+1)=wpos_red(2,iat)
                 pl_arr(3,iat,k_e_wpos+1)=wpos_red(3,iat)
                 f_arr(1,iat,k_e_wpos+1)=wpos_fcart(1,iat)
                 f_arr(2,iat,k_e_wpos+1)=wpos_fcart(2,iat)
                 f_arr(3,iat,k_e_wpos+1)=wpos_fcart(3,iat)
                 enddo
                 str_arr(:,k_e_wpos+1)=wpos_strten(:)
                 do i=1,3
                 lat_arr(1,i,k_e_wpos+1)=wpos_latvec(1,i)
                 lat_arr(2,i,k_e_wpos+1)=wpos_latvec(2,i)
                 lat_arr(3,i,k_e_wpos+1)=wpos_latvec(3,i)
                 enddo
                 spg_arr(k_e_wpos+1)=spg_wpos
                 spgtol_arr(k_e_wpos+1)=spgtol_wpos
                 dos_arr(k_e_wpos+1)=fdos_wpos
       call yaml_comment('INSERT')
       !write(*,*) '  -----   INSERT -----------'
       return
END SUBROUTINE insert
