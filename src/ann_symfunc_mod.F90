!*****************************************************************************************
module mod_symfunc
    use mod_linked_lists, only: typ_linked_lists
    use mod_symfunc_data, only: typ_symfunc_data
    implicit none
    private
    public:: typ_symfunc_data, typ_symfunc
    type, extends(typ_symfunc_data):: typ_symfunc
        contains
        procedure, public, pass(self):: get_symfunc
    end type typ_symfunc
    type, public:: typ_symfunc_arr
        integer:: nconf=-1
        type(typ_symfunc), allocatable:: symfunc(:)
    end type typ_symfunc_arr
contains
!*****************************************************************************************
subroutine get_symfunc(self,parini,ann_arr,atoms,apply_gbounds)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    use mod_symfunc_bond, only: symmetry_functions_driver_bond
    use time_profiling
    !use mod_timing , only: TCAT_SYMFUNC_COMPUT
    use dynamic_memory
    implicit none
    class(typ_symfunc), intent(inout):: self
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    logical, intent(in):: apply_gbounds
    !local variables
    integer:: i0, iat, jat, isat, i, ib, istat, ng, ig
    integer:: iats, iate, mat, mproc, ibs, ibe
    real(8):: gleft
    call f_routine(id='get_symfunc')
    !call f_timing(TCAT_SYMFUNC_COMPUT,'ON')
    !-----------------------------------------------------------------
    !first index of "y0d" is for number of symmetry function
    !second index of "y0d" is for x,y,z
    !3th index of "y0d" is for number of atoms
    !-----------------------------------------------------------------
    bondbased: if(parini%bondbased_ann) then
        call symmetry_functions_driver_bond(parini,ann_arr,atoms,self%typ_symfunc_data)
            do ib=1,self%linked_lists%maxbound_rad
                if(apply_gbounds) then
                    !normalization of y and y0d
                    do i0=1,ann_arr%ann(1)%nn(0)
                        gleft=ann_arr%ann(1)%gbounds(1,i0)
                        self%y(i0,ib)=(self%y(i0,ib)-gleft)*ann_arr%ann(1)%two_over_gdiff(i0)-1.d0
                        self%y0d_bond(i0,ib)=self%y0d_bond(i0,ib)*ann_arr%ann(1)%two_over_gdiff(i0)
                    enddo
                endif
                if(ann_arr%ann(1)%nn(0)/=ann_arr%ann(1)%ng1+ann_arr%ann(1)%ng2+ann_arr%ann(1)%ng3+ann_arr%ann(1)%ng4) then
                    write(*,'(a74,2i7)') 'ERROR: inconsistency between # of input nodes and # of symmetry functions:', &
                    ann_arr%ann(1)%nn(0),ann_arr%ann(1)%ng1+ann_arr%ann(1)%ng2+ann_arr%ann(1)%ng3+ann_arr%ann(1)%ng4
                    stop
                endif
            enddo
    else bondbased
        if (parini%symfunc_type_ann=='behler') then 
            call symmetry_functions_driver(parini,ann_arr,atoms,self%typ_symfunc_data)
        else
            call symmetry_functions_driver_stefan(parini,ann_arr,atoms,self%typ_symfunc_data)
        endif
        if(apply_gbounds) then
            do iat=1,atoms%nat
                isat=atoms%itypat(iat)
                !normalization of y and y0d
                do i0=1,ann_arr%ann(isat)%nn(0)
                    gleft=ann_arr%ann(isat)%gbounds(1,i0)
                    self%y(i0,iat)=(self%y(i0,iat)-gleft)*ann_arr%ann(isat)%two_over_gdiff(i0)-1.d0
                enddo
            enddo
            if(parini%mpi_env%nproc>1) then
                mat=atoms%nat/parini%mpi_env%nproc
                iats=parini%mpi_env%iproc*mat+1
                mproc=mod(atoms%nat,parini%mpi_env%nproc)
                iats=iats+max(0,parini%mpi_env%iproc-parini%mpi_env%nproc+mproc)
                if(parini%mpi_env%iproc>parini%mpi_env%nproc-mproc-1) mat=mat+1
                iate=iats+mat-1
            else
                iats=1
                iate=atoms%nat
                !mat=atoms%nat
            endif
            ibs=self%linked_lists%prime_bound(iats)
            ibe=self%linked_lists%prime_bound(iate+1)-1
            do ib=ibs,ibe !1,self%linked_lists%maxbound_rad
                iat=self%linked_lists%bound_rad(1,ib)
                isat=atoms%itypat(iat)
                do i0=1,ann_arr%ann(isat)%nn(0)
                    gleft=ann_arr%ann(isat)%gbounds(1,i0)
                    self%y0d(i0,1,ib)=self%y0d(i0,1,ib)*ann_arr%ann(isat)%two_over_gdiff(i0)
                    self%y0d(i0,2,ib)=self%y0d(i0,2,ib)*ann_arr%ann(isat)%two_over_gdiff(i0)
                    self%y0d(i0,3,ib)=self%y0d(i0,3,ib)*ann_arr%ann(isat)%two_over_gdiff(i0)
                    self%y0dr(i0,1,ib)=self%y0dr(i0,1,ib)*ann_arr%ann(isat)%two_over_gdiff(i0)
                    self%y0dr(i0,2,ib)=self%y0dr(i0,2,ib)*ann_arr%ann(isat)%two_over_gdiff(i0)
                    self%y0dr(i0,3,ib)=self%y0dr(i0,3,ib)*ann_arr%ann(isat)%two_over_gdiff(i0)
                    self%y0dr(i0,4,ib)=self%y0dr(i0,4,ib)*ann_arr%ann(isat)%two_over_gdiff(i0)
                    self%y0dr(i0,5,ib)=self%y0dr(i0,5,ib)*ann_arr%ann(isat)%two_over_gdiff(i0)
                    self%y0dr(i0,6,ib)=self%y0dr(i0,6,ib)*ann_arr%ann(isat)%two_over_gdiff(i0)
                    self%y0dr(i0,7,ib)=self%y0dr(i0,7,ib)*ann_arr%ann(isat)%two_over_gdiff(i0)
                    self%y0dr(i0,8,ib)=self%y0dr(i0,8,ib)*ann_arr%ann(isat)%two_over_gdiff(i0)
                    self%y0dr(i0,9,ib)=self%y0dr(i0,9,ib)*ann_arr%ann(isat)%two_over_gdiff(i0)
                enddo
            enddo
        endif
        !---------------------------------------------
        do iat=1,atoms%nat
            isat=atoms%itypat(iat)
            ng=ann_arr%ann(isat)%nn(0)
            if(ng/=ann_arr%ann(isat)%ng1+ann_arr%ann(isat)%ng2+ann_arr%ann(isat)%ng3+ &
                ann_arr%ann(isat)%ng4+ann_arr%ann(isat)%ng5+ann_arr%ann(isat)%ng6) then
                write(*,'(a74,2i7)') &
                    'ERROR: inconsistency between # of input nodes and # of symmetry functions:', &
                    ann_arr%ann(isat)%nn(0),ann_arr%ann(isat)%ng1+ann_arr%ann(isat)%ng2+ &
                    ann_arr%ann(isat)%ng3+ann_arr%ann(isat)%ng4+ann_arr%ann(isat)%ng5+ann_arr%ann(isat)%ng6
                stop
            endif
        enddo
    endif bondbased
    !call f_timing(TCAT_SYMFUNC_COMPUT,'OF')
    call f_release_routine()
end subroutine get_symfunc
!*****************************************************************************************
end module mod_symfunc
!*****************************************************************************************
