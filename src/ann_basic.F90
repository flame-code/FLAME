!*****************************************************************************************
subroutine ann_allocate(natmax,ann_arr)
    use mod_interface
    use mod_ann, only: typ_ann_arr
    implicit none
    integer, intent(in):: natmax
    type(typ_ann_arr), intent(inout):: ann_arr
    !local variables
    integer:: istat, ng
    ng=ann_arr%ann(1)%nn(0)
    !allocate(ann_arr%yall(ng,natmax),stat=istat)
    !if(istat/=0) stop 'ERROR: unable to allocate array ann_arr%yall.'
    !allocate(ann_arr%y0d(ng,3,natmax,natmax),stat=istat)
    !if(istat/=0) stop 'ERROR: unable to allocate array ann_arr%y0d.'
    !allocate(ann_arr%y0dr(ng,9,natmax,natmax),stat=istat)
    !if(istat/=0) stop 'ERROR: unable to allocate array ann_arr%y0dr.'
end subroutine ann_allocate
!*****************************************************************************************
subroutine ann_deallocate(ann_arr)
    use mod_interface
    use mod_ann, only: typ_ann_arr
    implicit none
    type(typ_ann_arr), intent(inout):: ann_arr
    !local variables
    integer:: istat
    !deallocate(ann_arr%yall,stat=istat)
    !if(istat/=0) stop 'ERROR: unable to deallocate array ann_arr%yall.'
    !deallocate(ann_arr%y0d,stat=istat)
    !if(istat/=0) stop 'ERROR: unable to deallocate array ann_arr%y0d.'
    !deallocate(ann_arr%y0dr,stat=istat)
    !if(istat/=0) stop 'ERROR: unable to deallocate array ann_arr%y0dr.'
end subroutine ann_deallocate
!*****************************************************************************************
