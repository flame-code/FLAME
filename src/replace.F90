subroutine replace(nlminx,nlmin,fp_len,nat,kid,e_wpos,ent_wpos,fp_wpos,wpos_red,&
  &wpos_latvec,spg_wpos,spgtol_wpos,fdos_wpos,&
  &e_arr,ent_arr,fp_arr,pl_arr,lat_arr,spg_arr,spgtol_arr,dos_arr,ct_arr,findsym)
!Replace the structure kid with wpos, only if the symmetry index is higher in wpos
  implicit none
  integer:: fp_len,ct_arr(nlminx),spg_arr(nlminx),nat,iat,spg_wpos
  integer:: k, nlmin,nlminx,i,kid
  real(8):: e_wpos, ent_wpos, wpos_red(3,nat),wpos_latvec(3,3),spgtol_wpos,fdos_wpos,fp_wpos(fp_len)
  real(8):: e_arr(nlminx),ent_arr(nlminx),fp_arr(fp_len,nlminx),pl_arr(3,nat,nlminx)
  real(8):: lat_arr(3,3,nlminx),spgtol_arr(nlminx),dos_arr(nlminx)
  logical:: findsym
  write(*,*) "KKKID",kid,nlmin,nlminx,fp_len,nat
  write(*,*) spg_wpos,spg_arr(kid),spg_wpos,spg_arr(kid),ent_wpos,ent_arr(kid)
  if((findsym.and.((spg_wpos.gt.spg_arr(kid)).or.(spg_wpos.eq.spg_arr(kid).and.ent_wpos.lt.ent_arr(kid)))).or.&
   &(.not.findsym.and.ent_wpos.lt.ent_arr(kid))) then
   write(*,'(a,i4,L3,i4,i4,es15.7,es15.7)') " # Replace array elements: ID,Findsym,SPG_ARR,SPG_NEW,ENT_ARR,ENT_NEW  ",&
   &kid,findsym,spg_arr(kid),spg_wpos,ent_arr(kid),ent_wpos
   e_arr(kid)=e_wpos   
   ent_arr(kid)=ent_wpos
   fp_arr(:,kid)=fp_wpos(:)
   pl_arr(:,:,kid)=wpos_red(:,:)   
   lat_arr(:,:,kid)=wpos_latvec(:,:)
   spgtol_arr(kid)=spgtol_wpos
   dos_arr(kid)=fdos_wpos   
  endif
end subroutine

