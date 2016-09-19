!*****************************************************************************************
subroutine ann_allocate(ekf,ann_arr)
    use mod_interface
    use mod_ann, only: typ_ann_arr, typ_ekf
    use dynamic_memory
    implicit none
    type(typ_ekf), intent(in):: ekf
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
    ann_arr%fat_chi=f_malloc0([1.to.3,1.to.ann_arr%natmax],id='fat_chi')
    ann_arr%chi_i=f_malloc0([1.to.ann_arr%natmax],id='ann_arr%chi_i')
    ann_arr%chi_o=f_malloc0([1.to.ann_arr%natmax],id='ann_arr%chi_o')
    ann_arr%chi_d=f_malloc0([1.to.ann_arr%natmax],id='ann_arr%chi_d')
    ann_arr%a=f_malloc0([1.to.(ann_arr%natmax+1)*(ann_arr%natmax+1)],id='a: aq=-chi')
    ann_arr%g_per_atom=f_malloc([1.to.ekf%num(1),1.to.ann_arr%natmax],id='g_per_atom') !HERE
    !symfunc%linked_lists%maxbound_rad is assumed 10000
    ann_arr%fatpq=f_malloc([1.to.3,1.to.10000],id='fatpq')
    ann_arr%stresspq=f_malloc([1.to.3,1.to.3,1.to.10000],id='stresspq')
end subroutine ann_allocate
!*****************************************************************************************
subroutine ann_deallocate(ann_arr)
    use mod_interface
    use mod_ann, only: typ_ann_arr
    use dynamic_memory
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
    call f_free(ann_arr%chi_i)
    call f_free(ann_arr%chi_o)
    call f_free(ann_arr%chi_d)
    call f_free(ann_arr%a)
    call f_free(ann_arr%fat_chi)
    call f_free(ann_arr%g_per_atom)
    call f_free(ann_arr%fatpq)
    call f_free(ann_arr%stresspq)
end subroutine ann_deallocate
!*****************************************************************************************
