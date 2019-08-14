!*****************************************************************************************
subroutine cal_architecture(ann,epot)
    use mod_ann, only: typ_ann
    implicit none
    type(typ_ann), intent(inout):: ann
    real(8), intent(inout):: epot
    !local variables
    if(ann%nl==3) then
        call cal_architecture_2hiddenlayer(ann,epot)
    else
        write(*,'(a,i7,a)') 'ERROR: ANN with ',ann%nl,' hidden layers not implemented'
    endif
end subroutine cal_architecture
!*****************************************************************************************
!This routine lack implementation of forces so cannot be called yet.
subroutine cal_architecture_1hiddenlayer(ann,epot)
    use mod_ann, only: typ_ann
    implicit none
    type(typ_ann), intent(inout):: ann
    real(8), intent(inout):: epot
    !local variables
    integer:: i, j
    real(8):: tt
    stop 'ERROR: This routine lack implementation of forces so cannot be called yet.'
    if(ann%nl/=2) then
        write(*,'(a,i3)') 'ERROR: this routine works only for ann%nl=2, while ann%nl= ',ann%nl
        stop
    endif
    do j=1,ann%nn(1)
        tt=0.d0
        do i=1,ann%nn(0)
            tt=tt+ann%a(i,j,1)*ann%y(i,0)
        enddo
        ann%y(j,1)=tanh(ann%b(j,1)+tt)
    enddo
    tt=0.d0
    do j=1,ann%nn(1)
        tt=tt+ann%a(j,1,2)*ann%y(j,1)
    enddo
    epot=ann%b(1,2)+tt
    !-------------------------------------------------------
end subroutine cal_architecture_1hiddenlayer
!*****************************************************************************************
subroutine cal_architecture_2hiddenlayer(ann,epot)
    use mod_ann, only: typ_ann
    implicit none
    type(typ_ann), intent(inout):: ann
    real(8), intent(inout):: epot
    !local variables
    integer:: i, j, k, l, lp
    real(8):: tt, tt2
    real(8):: c(100), ff(1000,1000)
    if(ann%nl/=3) then
        write(*,'(a,i3)') 'ERROR: this routine works only for ann%nl=3, while ann%nl= ',ann%nl
        stop
    endif
    !-------------------------------------------------------
    do j=1,ann%nn(1)
        tt=0.d0
        do i=1,ann%nn(0)
            tt=tt+ann%a(i,j,1)*ann%y(i,0)
        enddo
        ann%x(j,1)=ann%b(j,1)+tt
        tt2=tanh(ann%x(j,1))
        ann%y(j,1)=tt2
        ann%yd(j,1)=1.d0-tt2**2
    enddo
    !-------------------------------------------------------
    do k=1,ann%nn(2)
        tt=0.d0
        do j=1,ann%nn(1)
            tt=tt+ann%a(j,k,2)*ann%y(j,1)
        enddo
        ann%x(k,2)=ann%b(k,2)+tt
        tt2=tanh(ann%x(k,2))
        ann%y(k,2)=tt2
        ann%yd(k,2)=1.d0-tt2**2
    enddo
    !-------------------------------------------------------
    tt=0.d0
    do k=1,ann%nn(2)
        tt=tt+ann%a(k,1,3)*ann%y(k,2)
    enddo
    ann%x(1,3)=ann%b(1,3)+tt
    ann%y(1,3)=ann%x(1,3)
    ann%yd(1,3)=1.d0
    !-------------------------------------------------------
    epot=ann%y(1,3)
    !c(1:100)=0.d0
    !here "l" = "k" in my notes(page 8) and here "lp" = "k'"  in my notes.
    do l=1,ann%nn(1)
        tt=0.d0
        do lp=1,ann%nn(2)
            tt=tt+ann%a(lp,1,3)*ann%yd(lp,2)*ann%a(l,lp,2)
        enddo
        c(l)=tt
    enddo
    do j=1,ann%nn(0)
        tt=0.d0
        do l=1,ann%nn(1)
            tt=tt+ann%yd(l,1)*ann%a(j,l,1)*c(l)
        enddo
        ann%d(j)=tt
    enddo
end subroutine cal_architecture_2hiddenlayer
!*****************************************************************************************
subroutine cal_architecture_der(ann,epot)
    use mod_ann, only: typ_ann
    implicit none
    type(typ_ann), intent(inout):: ann
    real(8), intent(inout):: epot
    !local variables
    if(ann%nl==3) then
        call cal_architecture_der_2hiddenlayer(ann,epot)
    else
        write(*,'(a,i7,a)') 'ERROR: ANN with ',ann%nl,' hidden layers not implemented'
    endif
end subroutine cal_architecture_der
!****************************************************************************************
subroutine cal_architecture_der_1hiddenlayer(ann,epot)
    use mod_ann, only: typ_ann
    implicit none
    type(typ_ann), intent(inout):: ann
    real(8), intent(inout):: epot
    !local variables
    integer:: i, j, k, l, lp
    real(8):: tt
    if(ann%nl/=3) then
        write(*,'(a,i3)') 'ERROR: this routine works only for ann%nl=3, while ann%nl= ',ann%nl
        stop
    endif
    !-------------------------------------------------------
    do j=1,ann%nn(1)
        tt=0.d0
        do i=1,ann%nn(0)
            tt=tt+ann%a(i,j,1)*ann%y(i,0)
        enddo
        ann%x(j,1)=ann%b(j,1)+tt
        ann%y(j,1)=tanh(ann%x(j,1))
    enddo
    !-------------------------------------------------------
    tt=0.d0
    do j=1,ann%nn(1)
        tt=tt+ann%a(j,1,2)*ann%y(j,1)
    enddo
    ann%x(1,2)=ann%b(1,2)+tt
    epot=ann%x(1,2)
    !------------------------------------------------------
    do j=1,ann%nn(1)
         tt=ann%a(j,1,2)*(1.d0-ann%y(j,1)**2)
        ann%bd(j,1)=tt
        do lp=1,ann%nn(0)
            ann%ad(lp+(j-1)*ann%nn(0),2)=tt*ann%y(lp,0)
        enddo
    enddo    
    !-------------------------------------------------------
    tt=1.d0
    ann%bd(1,2)=tt
    do l=1,ann%nn(1)      
        ann%ad(l,2)=tt*ann%y(l,0)
    enddo
    !-------------------------------------------------------
end subroutine cal_architecture_der_1hiddenlayer
!*****************************************************************************************
subroutine cal_architecture_der_2hiddenlayer(ann,epot)
    use mod_ann, only: typ_ann
    implicit none
    type(typ_ann), intent(inout):: ann
    real(8), intent(inout):: epot
    !local variables
    integer:: i, j, k, l, lp
    real(8):: tt
    if(ann%nl/=3) then
        write(*,'(a,i3)') 'ERROR: this routine works only for ann%nl=3, while ann%nl= ',ann%nl
        stop
    endif
    !-------------------------------------------------------
    do j=1,ann%nn(1)
        tt=0.d0
        do i=1,ann%nn(0)
            tt=tt+ann%a(i,j,1)*ann%y(i,0)
        enddo
        ann%x(j,1)=ann%b(j,1)+tt
        ann%y(j,1)=tanh(ann%x(j,1))
    enddo
    !-------------------------------------------------------
    do k=1,ann%nn(2)
        tt=0.d0
        do j=1,ann%nn(1)
            tt=tt+ann%a(j,k,2)*ann%y(j,1)
        enddo
        ann%x(k,2)=ann%b(k,2)+tt
        ann%y(k,2)=tanh(ann%x(k,2))
    enddo
    !-------------------------------------------------------
    tt=0.d0
    do k=1,ann%nn(2)
        tt=tt+ann%a(k,1,3)*ann%y(k,2)
    enddo
    ann%x(1,3)=ann%b(1,3)+tt
    epot=ann%x(1,3)
    !-------------------------------------------------------
    do l=1,ann%nn(1)
        tt=0.d0
        do k=1,ann%nn(2)
            tt=tt+ann%a(k,1,3)*ann%a(l,k,2)*(1.d0-ann%y(k,2)**2)
        enddo
        tt=tt*(1.d0-ann%y(l,1)**2)
        ann%bd(l,1)=tt
        do lp=1,ann%nn(0)
            !ann%ad(lp,l,1)=tt*ann%y(lp,0)
            ann%ad(lp+(l-1)*ann%nn(0),1)=tt*ann%y(lp,0)
        enddo
    enddo
    !-------------------------------------------------------
    do l=1,ann%nn(2)
        tt=ann%a(l,1,3)*(1.d0-ann%y(l,2)**2)
        ann%bd(l,2)=tt
        do lp=1,ann%nn(1)
            !ann%ad(lp,l,2)=tt*ann%y(lp,1)
            ann%ad(lp+(l-1)*ann%nn(1),2)=tt*ann%y(lp,1)
        enddo
    enddo
    !-------------------------------------------------------
    tt=1.d0
    ann%bd(1,3)=tt
    do l=1,ann%nn(2)
        !ann%ad(l,1,3)=tt*ann%y(l,2)
        ann%ad(l,3)=tt*ann%y(l,2)
    enddo
    !-------------------------------------------------------
end subroutine cal_architecture_der_2hiddenlayer
!*****************************************************************************************
