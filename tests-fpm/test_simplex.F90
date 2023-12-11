!*****************************************************************************************
subroutine test_simplex()
    use iso_fortran_env, only: error_unit, output_unit
    use mod_colors, only: green_passed, red_failed
    implicit none
    !local variables
    real(8):: errmax, step, ftol !, time1, time2
    real(8), allocatable:: vertices(:,:), fval(:)
    integer:: ndim, iter
    external:: cal_cost_for_simplex
    ftol=1.d-10
    step=0.d0
    ndim=2
    allocate(vertices(ndim,ndim+1),fval(ndim+1))
    vertices(1,1)=1.0d0
    vertices(2,1)=1.5d0
    vertices(1,2)=2.d0
    vertices(2,2)=2.d0
    vertices(1,3)=3.d0
    vertices(2,3)=3.d0
    !call cpu_time(time1)
    call simplex(vertices,fval,step,ndim,ftol,cal_cost_for_simplex,iter)
    !write(*,'(a,3f12.8)') 'V1: a,b= ',vertices(1,1),vertices(2,1),fval(1)
    !write(*,'(a,3f12.8)') 'V2: a,b= ',vertices(1,2),vertices(2,2),fval(2)
    !write(*,'(a,3f12.8)') 'V3: a,b= ',vertices(1,3),vertices(2,3),fval(3)
    !call cpu_time(time2)
    !write(*,'(a,f8.3)') 'time= ',time2-time1
    errmax=0.d0
    errmax=max(errmax,abs(vertices(1,1)-2.d0))
    errmax=max(errmax,abs(vertices(2,1)-3.d0))
    errmax=max(errmax,abs(vertices(1,2)-2.d0))
    errmax=max(errmax,abs(vertices(2,2)-3.d0))
    errmax=max(errmax,abs(vertices(1,3)-2.d0))
    errmax=max(errmax,abs(vertices(2,3)-3.d0))
    deallocate(vertices,fval)
    if(errmax<1.d-6) then
        write(output_unit,'(2a)') green_passed,' in test_simplex: errmax'
    else
        write(error_unit,'(2a,es14.5)') red_failed,' in test_simplex: errmax=  ',errmax
        call exit(1)
    end if
end subroutine test_simplex
!*****************************************************************************************
subroutine cal_cost_for_simplex(ndim,vertex,cost)
    implicit none
    integer, intent(in) :: ndim
    real(8), intent(in) :: vertex(ndim)
    real(8), intent(out) :: cost
    !local variables
    real(8):: rmse, tt, x, a, b
    integer:: n, i
    real(8), allocatable:: points(:,:)
    n=100
    allocate(points(2,n))
    a=2.d0
    b=3.d0
    do i=1,n
        x=1.d-2*real(i,kind=8)
        points(1,i)=x
        points(2,i)=a*x+b+real(2*mod(i,2)-1,kind=8)*1.d-5*x**2
        write(21,*) points(1,i),points(2,i)
    enddo
    a=vertex(1)
    b=vertex(2)
    !write(31,'(a,2f8.3)') 'a,b= ',a,b
    cost=0.d0
    do i=1,n
        x=points(1,i)
        tt=a*x+b
        cost=cost+(tt-points(2,i))**2
    enddo
    cost=sqrt(cost)
    !write(*,*) 'COST ',cost
    deallocate(points)
end subroutine cal_cost_for_simplex
!*****************************************************************************************
