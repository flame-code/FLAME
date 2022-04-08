!*****************************************************************************************
module mod_ann
    use dictionaries
    use mod_linked_lists, only: typ_linked_lists
    use mod_electrostatics, only: typ_poisson
    !use mod_radpots_cent2, only: typ_radpots_cent2
    implicit none
    private
    public:: ann_arr_deallocate
    public:: convert_x_ann, convert_ann_x, convert_x_ann_arr
    public:: convert_ann_epotd
    public:: typ_ann_arr
    type, public:: typ_chi_arr
        real(8), allocatable:: chis(:)
    end type typ_chi_arr
    type, public:: typ_ann
        type(dictionary), pointer :: dict=>null()
        integer:: nl !number of hidden layer plus one
        integer:: nn(0:10)
        !integer:: n0=-1
        !integer:: n1=-1
        !integer:: n2=-1
        !integer:: n_all=n0*n1+n1*n2+n2+n1+n2+1
        integer:: ng1=-1
        integer:: ng2=-1
        integer:: ng3=-1
        integer:: ng4=-1
        integer:: ng5=-1
        integer:: ng6=-1
        real(8), allocatable:: a(:,:,:)
        real(8), allocatable:: b(:,:)
        real(8), allocatable:: x(:,:)
        real(8), allocatable:: y(:,:)
        real(8), allocatable:: yd(:,:)
        real(8), allocatable:: ad(:,:)
        real(8), allocatable:: bd(:,:)
        real(8), allocatable:: d(:)
        real(8), allocatable:: gbounds(:,:)
        real(8), allocatable:: two_over_gdiff(:)
        real(8):: gausswidth
        real(8):: gausswidth_ion
        real(8):: chi0
        real(8):: hardness
        real(8):: spring_const
        real(8):: qtarget
        real(8):: zion
        real(8):: qinit
        real(8):: rionic
        real(8):: ener_ref
        real(8):: ampl_chi=-1.d0
        real(8):: prefactor_chi=-1.d0
        character(20):: method

        !The 1st type of symmetry functions introduced by Behler
        real(8), allocatable:: g1eta(:)
        real(8), allocatable:: g1rs(:)

        !The 2nd type of symmetry functions introduced by Behler
        real(8), allocatable:: g2eta(:)
        real(8), allocatable:: g2rs(:)
        integer, allocatable:: g2i(:)

        !The 3rd type of symmetry functions introduced by Behler
        real(8), allocatable:: g3kappa(:)

        !The 4th type of symmetry functions introduced by Behler
        real(8), allocatable:: g4eta(:)
        real(8), allocatable:: g4zeta(:)
        real(8), allocatable:: g4lambda(:)

        !The 5th type of symmetry functions introduced by Behler
        real(8), allocatable:: g5eta(:)
        real(8), allocatable:: g5zeta(:)
        real(8), allocatable:: g5lambda(:)
        integer, allocatable:: g5i(:,:)

        !some other variables
        real(8), allocatable:: his(:,:)

        character(256):: hlines(10)
        contains
        procedure, public, pass(self):: ann_allocate
        
    end type typ_ann
    
    type:: typ_ann_arr
        logical:: exists_yaml_file = .false.
        integer:: iunit
        integer:: nann=-1
        integer:: natmax=1000
        integer:: nweight_max=-1
        logical:: compute_symfunc=.true.
        logical:: cal_force=.true.
        logical:: linear_rho_pot_initiated=.false. 
        logical:: trial_energy_required=.false.
        character(30):: event='unknown'
        character(50):: approach='unknown'
        character(50):: syslinsolver='unknown'
        real(8):: rcut=-1.d0
        real(8):: ener_ref
        real(8):: epot_es
        real(8):: fchi_angle
        real(8):: fchi_norm
        real(8):: dpm_err 
        real(8):: dpm_rmse
        !real(8), allocatable:: yall(:,:)
        !real(8), allocatable:: y0d(:,:,:)
        integer:: natsum(10)
        !real(8):: repfac(10,10)
        real(8):: reprcut(10,10)
        real(8):: qmax(10)
        real(8):: qmin(10)
        real(8):: qsum(10)
        real(8):: chi_max(10)
        real(8):: chi_min(10)
        real(8):: chi_sum(10)
        real(8):: chi_delta(10)
        real(8), allocatable:: yall_bond(:,:,:)
        real(8), allocatable:: y0d_bond(:,:,:,:)
        !real(8), allocatable:: y0dr(:,:,:)
        integer, allocatable:: loc(:)
        integer, allocatable, public:: num(:)
        real(8), allocatable:: a(:)
        real(8), allocatable:: chi_i(:)
        real(8), allocatable:: chi_o(:)
        real(8), allocatable:: chi_d(:)
        real(8), allocatable:: fat_chi(:,:)
        real(8), allocatable:: dqat_weights(:,:)
        real(8), allocatable:: g_per_atom(:,:)
        real(8), allocatable:: g_per_bond(:,:,:)
        real(8), allocatable:: fatpq(:,:)
        real(8), allocatable:: stresspq(:,:,:)
        integer, allocatable:: ipiv(:)
        real(8), allocatable:: qq(:)
        type(typ_ann), allocatable:: ann(:)
        !type(typ_radpots_cent2):: radpots_cent2
        type(typ_chi_arr), allocatable:: chi_tmp(:)
        type(typ_chi_arr), allocatable:: chi_ref_train(:)
        type(typ_chi_arr), allocatable:: chi_ref_valid(:)
        contains
        procedure, public, pass(self):: init_ann_arr
        procedure, public, pass(self):: fini_ann_arr
        procedure, public, pass(self):: set_number_of_ann
    end type typ_ann_arr
    type, public:: typ_cent
        real(8), allocatable:: gwi(:)
        real(8), allocatable:: gwe(:)
        real(8), allocatable:: gwit(:)
        real(8), allocatable:: rel(:,:)
        real(8), allocatable:: qgrad(:)
        real(8), allocatable:: rgrad(:,:)
        type(typ_poisson):: poisson
    end type typ_cent
contains
!*****************************************************************************************
subroutine set_number_of_ann(self,nann)
    implicit none
    class(typ_ann_arr), intent(inout):: self
    integer, intent(in):: nann
    !local variables
    self%nann=nann
end subroutine set_number_of_ann
!*****************************************************************************************
subroutine init_ann_arr(self)
    !use mod_opt_ann, only: typ_opt_ann
    use yaml_output
    use dynamic_memory
    implicit none
    class(typ_ann_arr), intent(inout):: self
    !local variables
    integer:: nw, ialpha, iann, nwtot
    do iann=1,self%nann
        nw=0
        do ialpha=1,self%ann(iann)%nl
            nw=nw+(self%ann(iann)%nn(ialpha-1)+1)*self%ann(iann)%nn(ialpha)
        enddo
        self%nweight_max=max(self%nweight_max,nw)
    enddo
    call ann_arr_allocate(self)
    self%num=f_malloc0([1.to.self%nann],id='self%num')
    self%loc=f_malloc0([1.to.self%nann],id='self%loc')
    nwtot=0
    call yaml_sequence_open('EKF') !,flow=.true.)
    do iann=1,self%nann
        do ialpha=1,self%ann(iann)%nl
            self%num(iann)=self%num(iann)+(self%ann(iann)%nn(ialpha-1)+1)*self%ann(iann)%nn(ialpha)
        enddo
        self%loc(iann)=nwtot+1
        nwtot=nwtot+self%num(iann)
        call yaml_sequence(advance='no')
        call yaml_map('iann',iann)
        call yaml_map('loc',self%loc(iann))
        call yaml_map('num',self%num(iann))
        call yaml_map('n',nwtot)
        !write(*,'(a,3i5)') 'EKF: ',self%loc(iann),self%num(iann),nwtot
    enddo
    call yaml_sequence_close()
end subroutine init_ann_arr
!*****************************************************************************************
subroutine fini_ann_arr(self)
    !use mod_opt_ann, only: typ_opt_ann
    use dynamic_memory
    implicit none
    class(typ_ann_arr), intent(inout):: self
    !local variables
    call ann_arr_deallocate(self)
    call f_free(self%num)
    call f_free(self%loc)
    !call self%radpots_cent2%fini_radpots_cent2()
end subroutine fini_ann_arr
!*****************************************************************************************
subroutine ann_allocate(self,nnmax,nl)
    implicit none
    class(typ_ann), intent(inout):: self
    integer, intent(in):: nnmax, nl
    !local variables
    !integer::
    allocate(self%a(nnmax,nnmax,nl))
    allocate(self%b(nnmax,nl))
    allocate(self%x(nnmax,nl))
    allocate(self%y(nnmax,0:nl))
    allocate(self%yd(nnmax,nl))
    allocate(self%ad(nnmax*nnmax,nl))
    allocate(self%bd(nnmax,nl))
    allocate(self%d(nnmax))
    allocate(self%gbounds(2,nnmax))
    allocate(self%two_over_gdiff(nnmax))
    allocate(self%g1eta(nnmax))
    allocate(self%g1rs(nnmax))
    allocate(self%g2eta(nnmax))
    allocate(self%g2rs(nnmax))
    allocate(self%g2i(nnmax))
    allocate(self%g3kappa(nnmax))
    allocate(self%g4eta(nnmax))
    allocate(self%g4zeta(nnmax))
    allocate(self%g4lambda(nnmax))
    allocate(self%g5eta(nnmax))
    allocate(self%g5zeta(nnmax))
    allocate(self%g5lambda(nnmax))
    allocate(self%g5i(2,nnmax))
    allocate(self%his(1000,nnmax))
end subroutine ann_allocate
!*****************************************************************************************
subroutine ann_arr_allocate(ann_arr)
    !use mod_opt_ann, only: typ_opt_ann
    use dynamic_memory
    implicit none
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
    if(.not. allocated(ann_arr%a)) then
        allocate(ann_arr%a(1:(ann_arr%natmax+1)*(ann_arr%natmax+1)))
    end if
    ann_arr%fat_chi=0.d0
    ann_arr%chi_i=0.d0
    ann_arr%chi_o=0.d0
    ann_arr%chi_d=0.d0
    if(allocated(ann_arr%a)) then
        ann_arr%a=0.d0
    end if
    allocate(ann_arr%dqat_weights(ann_arr%nweight_max,ann_arr%natmax))
    allocate(ann_arr%g_per_atom(ann_arr%nweight_max,ann_arr%natmax))
    !symfunc%linked_lists%maxbound_rad is assumed 10000
    allocate(ann_arr%fatpq(1:3,1:20000))
    allocate(ann_arr%stresspq(1:3,1:3,1:20000))
    allocate(ann_arr%ipiv(1:ann_arr%natmax+1))
    allocate(ann_arr%qq(1:ann_arr%natmax+1))
    allocate(ann_arr%yall_bond(10,10,10))
    allocate(ann_arr%y0d_bond(10,3,10,10))
end subroutine ann_arr_allocate
!*****************************************************************************************
subroutine ann_arr_deallocate(ann_arr)
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
    if(allocated(ann_arr%a)) then
        deallocate(ann_arr%a)
    endif
    deallocate(ann_arr%fat_chi)
    deallocate(ann_arr%dqat_weights)
    deallocate(ann_arr%g_per_atom)
    deallocate(ann_arr%fatpq)
    deallocate(ann_arr%stresspq)
    if (allocated(ann_arr%ipiv)) deallocate(ann_arr%ipiv)
    deallocate(ann_arr%qq)
    deallocate(ann_arr%yall_bond)
    deallocate(ann_arr%y0d_bond)
end subroutine ann_arr_deallocate
!*****************************************************************************************
subroutine convert_x_ann(n,x,ann)
    implicit none
    integer, intent(in):: n
    real(8), intent(in):: x(n)
    type(typ_ann), intent(inout):: ann
    !local variables
    integer:: i, j, l, ialpha
    l=0
    do ialpha=1,ann%nl
        do j=1,ann%nn(ialpha)
            do i=1,ann%nn(ialpha-1)
                l=l+1
                ann%a(i,j,ialpha)=x(l)
            enddo
        enddo
        do i=1,ann%nn(ialpha)
            l=l+1
            ann%b(i,ialpha)=x(l)
        enddo
    enddo
    if(l/=n) stop 'ERROR: l/=n'
end subroutine convert_x_ann
!*****************************************************************************************
subroutine convert_ann_x(n,x,ann)
    implicit none
    integer, intent(in):: n
    real(8), intent(inout):: x(n)
    type(typ_ann), intent(in):: ann
    !local variables
    integer:: i, j, l, ialpha
    l=0
    do ialpha=1,ann%nl
        do j=1,ann%nn(ialpha)
            do i=1,ann%nn(ialpha-1)
                l=l+1
                x(l)=ann%a(i,j,ialpha)
            enddo
        enddo
        do i=1,ann%nn(ialpha)
            l=l+1
            x(l)=ann%b(i,ialpha)
        enddo
    enddo
    if(l/=n) stop 'ERROR: l/=n'
end subroutine convert_ann_x
!*****************************************************************************************
subroutine convert_x_ann_arr(nwtot,x,ann_arr)
    use dynamic_memory
    implicit none
    integer, intent(in):: nwtot
    real(8), intent(in):: x(nwtot)
    type(typ_ann_arr), intent(inout):: ann_arr
    !local variables
    integer:: ia, n
    n=sum(ann_arr%num(1:ann_arr%nann))
    if(n/=nwtot) then
        write(*,*) 'ERROR: inconsistency in total number of weights in convert_x_ann_arr'
        stop
    endif
    do ia=1,ann_arr%nann
        call convert_x_ann(ann_arr%num(ia),x(ann_arr%loc(ia)),ann_arr%ann(ia))
    enddo
end subroutine convert_x_ann_arr
!*****************************************************************************************
subroutine convert_ann_epotd(ann,n,epotd)
    implicit none
    type(typ_ann), intent(in):: ann
    integer, intent(in):: n
    real(8), intent(inout):: epotd(n)
    !local variables
    integer:: i, ij, l, ialpha
    l=0
    do ialpha=1,ann%nl
        do ij=1,ann%nn(ialpha)*ann%nn(ialpha-1)
            l=l+1
            epotd(l)=ann%ad(ij,ialpha)
        enddo
        do i=1,ann%nn(ialpha)
            l=l+1
            epotd(l)=ann%bd(i,ialpha)
        enddo
    enddo
    if(l/=n) stop 'ERROR: l/=n'
end subroutine convert_ann_epotd
!*****************************************************************************************
end module mod_ann
!*****************************************************************************************
!module data_point
!    implicit none
!    integer:: m
!    integer:: natmax=-1
!    integer, allocatable:: natarr(:)
!    real(8), allocatable:: pnt(:), f_p(:)
!    real(8), allocatable:: ratall(:,:,:), epotall(:)
!end module data_point
!*****************************************************************************************
module mod_parlm
    implicit none
    type typ_parlm
        real(8):: ftol=1.d-8
        real(8):: xtol=1.d-8
        real(8):: gtol=1.d-8
        real(8):: factor=100.d0
        integer:: maxfev=1000
        integer:: nprint=1
        integer:: n=0
        integer:: mode, info, nfev, njev
        integer:: iter
        integer:: icontinue
        integer:: iflag
        logical:: finish
        real(8):: epsmch
        real(8):: fnorm
        !real(8):: fnorm1
        real(8):: xnorm
        real(8):: gnorm
        real(8):: pnorm
        real(8):: par
        real(8):: delta
        real(8):: actred
        real(8):: prered
        real(8):: ratio
        real(8), allocatable:: wa1(:), wa2(:), wa3(:), wa4(:), qtf(:)
        real(8), allocatable:: x(:)
        real(8), allocatable:: fvec(:)
        real(8), allocatable:: fjac(:,:)
        real(8), allocatable:: diag(:)
        integer, allocatable:: ipvt(:)
    end type typ_parlm
end module mod_parlm
!*****************************************************************************************
