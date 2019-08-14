!*****************************************************************************************
subroutine cal_matvec_mpi(n,p,g,v1)
    use mod_processors, only: iproc, nproc, mpi_comm_abz
    integer, intent(in):: n
    real(8), intent(in):: p(n,n), g(n)
    real(8), intent(out):: v1(n)
    !local variables
    integer:: i, j, m, mm, ierr
    real(8):: tt
    real(8), allocatable:: v2(:)
    integer, save:: narr(2,0:99), ni(2), nlarr(0:99), nnarr(0:99)
    integer, save:: icall=0
    !real(8):: time1, time2, time3
    !real(8), save:: dtime1=0.d0, dtime2=0.d0
#if defined(MPI)
    include 'mpif.h'
    !integer:: status_mpi(MPI_STATUS_SIZE)
    associate(MPI_DP=>MPI_DOUBLE_PRECISION)
    icall=icall+1
    if(icall==1) then
    ni(2)=n/nproc
    mm=mod(n,nproc)
    ni(1)=ni(2)*iproc+1
    if(mm/=0) then
        ni(1)=ni(1)+iproc
        if(iproc<mm) then
            ni(2)=ni(2)+1
        endif
    endif
    call MPI_ALLGATHER(ni,2,MPI_INTEGER,narr,2,MPI_INTEGER,mpi_comm_abz,ierr)
    nlarr(0:nproc-1)=narr(1,0:nproc-1)-1
    nnarr(0:nproc-1)=narr(2,0:nproc-1)
    endif
    allocate(v2(n))
    !now each MPI process calculates ni(2) elements
    !call cpu_time(time1)
    m=mod(n,5)
    do j=ni(1),ni(1)+ni(2)-1
        tt=0.d0
        do i=1,m
            tt=tt+p(i,j)*g(i)
        enddo
        do i=m+1,n,5
            tt=tt+p(i+0,j)*g(i+0) &
                 +p(i+1,j)*g(i+1) &
                 +p(i+2,j)*g(i+2) &
                 +p(i+3,j)*g(i+3) &
                 +p(i+4,j)*g(i+4)
        enddo
        v2(j)=tt
    enddo
    !call cpu_time(time2)
    call MPI_ALLGATHERV(v2(ni(1)),ni(2),MPI_DP,v1,nnarr,nlarr,MPI_DP,mpi_comm_abz,ierr)
    !call cpu_time(time3)
    !dtime1=dtime1+time2-time1
    !dtime2=dtime2+time3-time2
    !if(icall==1962) then
    !    write(*,'(a,i3,2f10.2)') 'cal_matvec_mpi ',iproc,dtime1,dtime2
    !endif
    end associate
    deallocate(v2)
!!!  #else
#endif
end subroutine cal_matvec_mpi
!*****************************************************************************************
