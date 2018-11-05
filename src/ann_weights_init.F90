!*****************************************************************************************
subroutine set_annweights(parini,opt_ann)
    use mod_interface
    use mod_parini, only: typ_parini
    !use mod_ann, only: 
    use mod_opt_ann, only: typ_opt_ann
    use mod_processors, only: iproc, mpi_comm_abz
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_opt_ann), intent(inout):: opt_ann
    !local variables
    integer:: ierr, i
    real(8):: tt
    character(50):: approach
#if defined(MPI)
    include 'mpif.h'
#endif
    !approach='uniform'
    approach='pure_electrostatic'
    !approach='type_dependent'
    if(iproc==0) then
        if(trim(approach)=='uniform') then
            call random_number(opt_ann%x)
            opt_ann%x(1:opt_ann%n)=(opt_ann%x(1:opt_ann%n)-0.5d0)*2.d0*parini%ampl_rand
        elseif(trim(approach)=='pure_electrostatic') then
            do i=1,opt_ann%n
                if(i<=opt_ann%n/2) then
                    call random_number(opt_ann%x(i))
                    tt=-2.d0*parini%ampl_rand
                    opt_ann%x(i)=(opt_ann%x(i)-0.5d0)*tt
                else
                    opt_ann%x(i)=-opt_ann%x(i-opt_ann%n/2)
                endif
            enddo
            !opt_ann%x(opt_ann%n/2)=0.d0
            !opt_ann%x(opt_ann%n)=0.d0
        elseif(trim(approach)=='type_dependent') then
            do i=1,opt_ann%n
                if(i<=opt_ann%n/4) then
                    call random_number(opt_ann%x(i))
                    tt=2.d0*parini%ampl_rand
                    opt_ann%x(i)=(opt_ann%x(i)-0.5d0)*tt
                elseif(i<=opt_ann%n/2) then
                    opt_ann%x(i)=-opt_ann%x(i-opt_ann%n/4)
                elseif(i<=3*opt_ann%n/4) then
                    call random_number(opt_ann%x(i))
                    tt=2.d0*parini%ampl_rand
                    opt_ann%x(i)=-(opt_ann%x(i)-0.5d0)*tt
                else
                    opt_ann%x(i)=-opt_ann%x(i-opt_ann%n/4)
                endif
            enddo
        else
            stop 'ERROR: unknown approach in set_annweights'
        endif

    endif
#if defined(MPI)
    call MPI_BCAST(opt_ann%x,opt_ann%n,MPI_DOUBLE_PRECISION,0,mpi_comm_abz,ierr)
#endif
end subroutine set_annweights
!*****************************************************************************************
