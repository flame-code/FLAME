!*****************************************************************************************
subroutine ann_allocate(opt_ann,ann_arr)
    use mod_interface
    use mod_ann, only: typ_ann_arr
    use mod_opt_ann, only: typ_opt_ann
    use dynamic_memory
    implicit none
    type(typ_opt_ann), intent(in):: opt_ann
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
    allocate(ann_arr%fat_chi(1:3,1:ann_arr%natmax))
    allocate(ann_arr%chi_i(1:ann_arr%natmax))
    allocate(ann_arr%chi_o(1:ann_arr%natmax))
    allocate(ann_arr%chi_d(1:ann_arr%natmax))
    allocate(ann_arr%a(1:(ann_arr%natmax+1)*(ann_arr%natmax+1)))
    ann_arr%fat_chi=0.d0
    ann_arr%chi_i=0.d0
    ann_arr%chi_o=0.d0
    ann_arr%chi_d=0.d0
    ann_arr%a=0.d0
    allocate(ann_arr%g_per_atom(1:opt_ann%num(1),1:ann_arr%natmax))
    !symfunc%linked_lists%maxbound_rad is assumed 10000
    allocate(ann_arr%fatpq(1:3,1:10000))
    allocate(ann_arr%stresspq(1:3,1:3,1:10000))
    allocate(ann_arr%ipiv(1:ann_arr%natmax+1))
    allocate(ann_arr%qq(1:ann_arr%natmax+1))
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
    deallocate(ann_arr%chi_i)
    deallocate(ann_arr%chi_o)
    deallocate(ann_arr%chi_d)
    deallocate(ann_arr%a)
    deallocate(ann_arr%fat_chi)
    deallocate(ann_arr%g_per_atom)
    deallocate(ann_arr%fatpq)
    deallocate(ann_arr%stresspq)
    deallocate(ann_arr%ipiv)
    deallocate(ann_arr%qq)
end subroutine ann_deallocate
!*****************************************************************************************
