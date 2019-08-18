subroutine identical(parini,nlminx,nlmin,fp_method,fp_len,ent_wpos,fp_wpos,ent_arr,fp_arr,&
           &ent_delta,fp_delta,newmin,kid,fp_dist_min,k_e_wpos,n_unique,n_nonuni,lid,nid)
use mod_parini, only: typ_parini
use yaml_output
implicit none
type(typ_parini), intent(in):: parini
integer:: nlminx,nlmin,fp_len,kid,k_e_wpos,n_unique,n_nonuni
integer:: i,l,klow,k,khigh,fp_method,lid(nlminx),nid
real(8):: fp_arr(fp_len,nlminx),fp_wpos(fp_len),ent_arr(nlminx),ent_wpos,fp_delta,ent_delta,fp_dist_min,fp_dist
logical newmin,inrange

inrange=.false.

!check whether new minimum
call hunt(ent_arr,min(nlmin,nlminx),ent_wpos,k_e_wpos)
newmin=.true.
!!do i=1,nlmin
!!write(*,'(a,i3,5(e24.17))') ' # check identical: enarr ',i,ent_arr(i),(fp_arr(l,i),l=1,2)
!!enddo
call yaml_mapping_open('Check identical',flow=.true.)
call yaml_map('ent_wpos',ent_wpos,fmt='(e24.17)')
call yaml_map('k_e_wpos',k_e_wpos,fmt='(i8)')
call yaml_mapping_close()
!write(*,'(a,e24.17,i3,5(e24.17))') ' # ID: Check identical: ent_wpos,k_e_wpos ',ent_wpos,k_e_wpos!,(wfp(l),l=1,2)

! find lowest configuration that might be identical
klow=k_e_wpos
do k=k_e_wpos,1,-1
if (ent_wpos-ent_arr(k).lt.0.d0) stop 'zeroA'
if (ent_wpos-ent_arr(k).lt.ent_delta) inrange=.true.
if (ent_wpos-ent_arr(k).gt.ent_delta) exit
klow=k
enddo

! find highest  configuration that might be identical
khigh=k_e_wpos+1
do k=k_e_wpos+1,nlmin
if (ent_arr(k)-ent_wpos.lt.0.d0) stop 'zeroB'
if (ent_arr(k)-ent_wpos.lt.ent_delta) inrange=.true.
if (ent_arr(k)-ent_wpos.gt.ent_delta) exit
khigh=k
enddo

nid=0
if(.not.inrange) then
    call yaml_comment('ID: there is no minimum in the given enthalpy range')
    call yaml_map('wpos is a new minimum',newmin)
  !write(*,'(a)') ' # ID: there is no minimum in the given enthalpy range'
  !write(*,'(a,L4)') ' # ID: wpos is a new minimum: ',newmin
  return  
endif
write(*,'(a,i5,i5)') ' # ID: Check k bounds: ',max(1,klow),min(nlmin,khigh)
fp_dist_min=1.d100
do k=max(1,klow),min(nlmin,khigh)
call get_fp_distance(parini,fp_len,fp_wpos,fp_arr(:,k),fp_dist)
write(*,'(a,i5,es25.15)') ' # ID: Check fp_dist: ',k,fp_dist
!!write(*,'(a,20(e10.3))') '(MH) fp_wpos', (fp_wpos(i),i=1,fp_len)
!!write(*,'(a,20(e10.3))') '(MH) fp_arr ', (fp_arr(i,k),i=1,fp_len)
if (fp_dist.lt.fp_delta) then
    write(*,'(a,i5)') ' # ID: wpos is identical to ',k
    newmin=.false.
    nid=nid+1
    lid(nid)=k
    if (fp_dist.lt.fp_dist_min) then
       fp_dist_min=fp_dist
       kid=k
       k_e_wpos=kid
    endif
endif
enddo

   call yaml_map('wpos is a new minimum',newmin)
   !write(*,'(a,L4)') ' # ID: wpos is a new minimum: ',newmin
   if (nid.gt.1) write(*,*) '# WARNING: more than one identical configuration found'
!          if (nsm.gt.1) write(100+iproc,*) 'WARNING: more than one identical configuration found'
if (nid.eq.1) n_unique=n_unique+1
if (nid.gt.1) n_nonuni=n_nonuni+1

return
end subroutine
