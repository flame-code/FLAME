!*****************************************************************************************
subroutine set_annweights(parini,ekf)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ekf
    use mod_processors, only: iproc, mpi_comm_abz
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ekf), intent(inout):: ekf
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
            call random_number(ekf%x)
            ekf%x(1:ekf%n)=(ekf%x(1:ekf%n)-0.5d0)*2.d0*parini%ampl_rand
        elseif(trim(approach)=='pure_electrostatic') then
            do i=1,ekf%n
                if(i<=ekf%n/2) then
                    call random_number(ekf%x(i))
                    tt=-2.d0*parini%ampl_rand
                    ekf%x(i)=(ekf%x(i)-0.5d0)*tt
                else
                    ekf%x(i)=-ekf%x(i-ekf%n/2)
                endif
            enddo
            !ekf%x(ekf%n/2)=0.d0
            !ekf%x(ekf%n)=0.d0
        elseif(trim(approach)=='type_dependent') then
            do i=1,ekf%n
                if(i<=ekf%n/4) then
                    call random_number(ekf%x(i))
                    tt=2.d0*parini%ampl_rand
                    ekf%x(i)=(ekf%x(i)-0.5d0)*tt
                elseif(i<=ekf%n/2) then
                    ekf%x(i)=-ekf%x(i-ekf%n/4)
                elseif(i<=3*ekf%n/4) then
                    call random_number(ekf%x(i))
                    tt=2.d0*parini%ampl_rand
                    ekf%x(i)=-(ekf%x(i)-0.5d0)*tt
                else
                    ekf%x(i)=-ekf%x(i-ekf%n/4)
                endif
            enddo
        else
            stop 'ERROR: unknown approach in set_annweights'
        endif

    endif
#if defined(MPI)
    call MPI_BCAST(ekf%x,ekf%n,MPI_DOUBLE_PRECISION,0,mpi_comm_abz,ierr)
#endif
end subroutine set_annweights
!*****************************************************************************************
