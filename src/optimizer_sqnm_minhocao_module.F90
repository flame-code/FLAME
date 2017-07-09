!! @file
!! @author Stefan Goedecker and Bastian Schaefer
!! @section LICENCE
!!    Copyright (C) 2014 UNIBAS
!!    This file is not freely distributed.
!!    A licence is necessary from UNIBAS

module module_sqn
    implicit none
    private

    public :: modify_gradient_minhocao
    public :: getSubSpaceEvecEval_minhocao
!    public :: findbonds
!    public :: projectbond

contains

subroutine modify_gradient_minhocao(nat,ndim,rrr,eval,res,fxyz,alpha,alphalat_scale,dd)
!    use module_base
    implicit none
    !parameters
    integer, intent(in) :: nat
    integer, intent(in) :: ndim
    real(8), intent(out):: dd(3,nat+3)
    real(8), intent(in) :: fxyz(3,nat+3)
    real(8), intent(in) :: rrr(3,nat+3,ndim)
    real(8), intent(in) :: eval(ndim)
    real(8), intent(in) :: res(ndim)
    real(8), intent(in) :: alpha
    !internal
    integer :: iat,i,l
    real(8) :: scpr(ndim)
    real(8) :: tt
    real(8) :: alpha_lat,alphalat_scale
    ! decompose gradient

    do iat=1,nat+3
        do l=1,3
            dd(l,iat)=-fxyz(l,iat)
        enddo
    enddo
    
    do i=1,ndim
        scpr(i)=0.0d0
        do iat=1,nat!+3
            do l=1,3
                scpr(i)=scpr(i)-fxyz(l,iat)*rrr(l,iat,i)
            enddo
        enddo
        do iat=1,nat!+3
            do l=1,3
                dd(l,iat)=dd(l,iat)-scpr(i)*rrr(l,iat,i)
            enddo
        enddo
    enddo

    !simple sd in space orthogonal to relevant subspace
    alpha_lat=alpha*alphalat_scale
    do iat=1,nat
        do l=1,3
            dd(l,iat)=dd(l,iat)*alpha
        enddo
    enddo
    do iat=nat+1,nat+3
        do l=1,3
            dd(l,iat)=dd(l,iat)*alpha_lat
        enddo
    enddo

    do i=1,ndim
    !quasi newton in relevant subspace
        tt=scpr(i)/sqrt(eval(i)**2+res(i)**2)
        do iat=1,nat+3
            do l=1,3
                dd(l,iat)=dd(l,iat)+tt*rrr(l,iat,i)
            enddo
        enddo
    enddo
end subroutine

subroutine getSubSpaceEvecEval_minhocao(label,verbosity,nat,nhist,nhistx&
                              ,ndim,cutoffratio,lwork,work,rxyz,&
                              fxyz,aa,rr,ff,rrr,fff,eval,res,success)
!    use module_base
!    use yaml_output
    !hard-coded parameters:
    !threshold for linear dependency:
    !if (eval(idim)/eval(nhist).gt.1.d-4) then
    implicit none
    !parameters
    integer, intent(in) :: verbosity,nat,nhist,nhistx,lwork
    character(len=*), intent(in) :: label
    integer, intent(out) :: ndim
    real(8), intent(in) :: rxyz(3,nat+3,0:nhistx),fxyz(3,nat+3,0:nhistx)
    real(8), intent(out) :: aa(nhistx,nhistx),eval(nhistx)
    real(8), intent(out) :: work(lwork)
    real(8), intent(out) :: rr(3,nat+3,0:nhistx), ff(3,nat+3,0:nhistx)
    real(8), intent(out) :: rrr(3,nat+3,0:nhistx), fff(3,nat+3,0:nhistx)
    real(8), intent(out) :: res(nhistx)
    real(8), intent(in) :: cutoffratio
    logical, intent(out) :: success
    !internal
    integer :: i,j,l,iat,info,idim,jdim
    real(8) :: tt
    real(8) :: rnorm(nhistx)
    logical :: debug
    debug=.false.

    success=.false.

    ! calculate norms
    do i=1,nhist
        rnorm(i)=0.0d0
         do iat=1,nat!+3
             do l=1,3
                rnorm(i)=rnorm(i) + (rxyz(l,iat,i)-rxyz(l,iat,i-1))**2
             enddo
         enddo
         rnorm(i)=1.0d0/sqrt(rnorm(i))
    enddo

    !find linear dependencies via diagonalization of overlap matrix   
    !build overlap matrix:
    do i=1,nhist
        do j=1,nhist
            aa(i,j)=0.0d0
            do iat=1,nat!+3
                do l=1,3
                aa(i,j)=aa(i,j) + (rxyz(l,iat,i)-rxyz(l,iat,i-1))&
                &*(rxyz(l,iat,j)-rxyz(l,iat,j-1))
                enddo
            enddo
            aa(i,j)=aa(i,j)*rnorm(i)*rnorm(j)
         enddo
    enddo

    call dsyev('V',"L",nhist,aa,nhistx,eval,work,lwork,info)
!    if (info.ne.0) then
    if (debug) then
!        call yaml_warning(trim(adjustl(label))//' 1st DSYEV '//&
!        '(Overlapmatrix) in getSupSpaceEvecEval failed with info: '//&
!        trim(yaml_toa(info))//', iproc: '//trim(yaml_toa(iproc)))
        write(*,*) trim(adjustl(label)),' 1st DSYEV (Overlapmatrix) in getSupSpaceEvecEval failed with info: ',&
        info
        return
!        stop 'info'
    endif
!    if(iproc==0 .and. verbosity>=3)then
    if(debug)then
        do i=1,nhist
!            call yaml_scalar(trim(adjustl(label))//' Overlap '//&
!            'eigenvalues: '//trim(yaml_toa(i))//' '//&
!            trim(yaml_toa(eval(i))))
            write(*,*)  trim(adjustl(label)),' Overlap eigenvalues: ',&
            i,eval(i)
        enddo
    endif

    do idim=1,nhist
        do iat=1,nat!+3
            do l=1,3
                rr(l,iat,idim)=0.0d0
                ff(l,iat,idim)=0.0d0
            enddo
        enddo
    enddo

    ndim=0
    do idim=1,nhist
        !remove linear dependencies by using
        !the overlap-matrix eigenvalues:
        if (eval(idim)/eval(nhist).gt.cutoffratio) then    ! HERE
            ndim=ndim+1

            do jdim=1,nhist
                do iat=1,nat!+3
                    do l=1,3
                         rr(l,iat,ndim)=rr(l,iat,ndim)+&
                                 aa(jdim,idim)*rnorm(jdim)*&
                                 (rxyz(l,iat,jdim)-rxyz(l,iat,jdim-1))
                         ff(l,iat,ndim)=ff(l,iat,ndim)-&
                                 aa(jdim,idim)*rnorm(jdim)*&
                                 (fxyz(l,iat,jdim)-fxyz(l,iat,jdim-1))

                    enddo
                enddo
            enddo

            do iat=1,nat!+3
                do l=1,3
                    rr(l,iat,ndim)=rr(l,iat,ndim)/&
                                    sqrt(abs(eval(idim)))
                    ff(l,iat,ndim)=ff(l,iat,ndim)/&
                                    sqrt(abs(eval(idim)))
                enddo
            enddo
        endif
    enddo

    ! Hessian matrix in significant orthogonal subspace
    do i=1,ndim
        do j=1,ndim
            aa(i,j)=0.0d0
            do iat=1,nat!+3
                do l=1,3
                    aa(i,j)=aa(i,j) + .5d0*(rr(l,iat,i)*ff(l,iat,j)+&
                                      &rr(l,iat,j)*ff(l,iat,i))
                enddo
            enddo
        enddo
    enddo

    call dsyev('V',"L",ndim,aa,nhistx,eval,work,lwork,info)
!    if (info.ne.0) then
    if (debug) then
!        call yaml_warning(trim(adjustl(label))//' 2nd DSYEV '//&
!        '(subpsace hessian) in getSupSpaceEvecEval failed with info: '&
!         //trim(yaml_toa(info))//', iproc:'//trim(yaml_toa(iproc)))
        write(*,*) trim(adjustl(label)),' 2nd DSYEV (subpsace hessian) in getSupSpaceEvecEval failed with info: ',&
         info
        return
!        stop 'info'
    endif

    ! calculate vectors in full 3*nat-dim space
    do i=1,ndim
        do iat=1,nat!+3
            do l=1,3
                rrr(l,iat,i)=0.0d0
                fff(l,iat,i)=0.0d0
            enddo
        enddo
    enddo

    do i=1,ndim
        tt=0.0d0
        do j=1,ndim
            do iat=1,nat!+3
                do l=1,3
                    rrr(l,iat,i)=rrr(l,iat,i) + aa(j,i)*rr(l,iat,j)
                    fff(l,iat,i)=fff(l,iat,i) + aa(j,i)*ff(l,iat,j)
                enddo
            enddo
        enddo
        do iat=1,nat!+3
            do l=1,3
                tt=tt+(fff(l,iat,i)-eval(i)*rrr(l,iat,i))**2
            enddo
        enddo
        !residuue according to Weinstein criterion
        res(i)=sqrt(tt)
!        if(iproc==0 .and. verbosity>=3)&
!        if(verbosity>=2)&
        if(debug)&
!            call yaml_scalar(trim(adjustl(label))//' i, '//&
!            'eigenvalue, residue: '//trim(yaml_toa(i))//' '//&
!            trim(yaml_toa(eval(i)))//' '//trim(yaml_toa(res(i))))
            write(*,*) trim(adjustl(label)),' i, eigenvalue, residue: ',i,&
            eval(i),res(i)
    enddo
    success=.true.
end subroutine

!!subroutine findbonds(label,verbosity,nat,rcov,pos,nbond,&
!!                    iconnect)
!!!has to be called before findsad (if operating in biomolecule mode)
!!!    use module_base
!!!    use yaml_output
!!    implicit none
!!    !parameters
!!    integer, intent(in) :: verbosity,nat
!!    character(len=*), intent(in) :: label
!!    real(8), intent(in) :: rcov(nat)
!!    real(8), intent(in) :: pos(3,nat)
!!    integer, intent(out) :: nbond
!!    integer, intent(out) :: iconnect(2,1000)
!!    !internal
!!    integer :: iat,jat
!!    real(8) :: dist2
!!    nbond=0
!!    do iat=1,nat
!!        do jat=1,iat-1
!!            dist2=(pos(1,iat)-pos(1,jat))**2+&
!!                  (pos(2,iat)-pos(2,jat))**2+&
!!                  (pos(3,iat)-pos(3,jat))**2
!!            if (dist2.le.(1.2d0*(rcov(iat)+rcov(jat)))**2) then
!!                nbond=nbond+1
!!                if (nbond.gt.1000) stop &
!!                     'nbond>1000, increase size of iconnect in '//&
!!                     'routine which calls subroutine findbonds'
!!                iconnect(1,nbond)=iat
!!                iconnect(2,nbond)=jat
!!            endif
!!        enddo
!!    enddo
!!!    if(iproc==0.and.verbosity>=2)&
!!    if(verbosity>=2)&
!!!        call yaml_scalar(trim(adjustl(label))//&
!!!        ' Found'//trim(yaml_toa(nbond))//' bonds.')
!!        write(*,*) trim(adjustl(label)),&
!!        ' Found',nbond,' bonds.'
!!end subroutine
!!
!!subroutine projectbond(nat,nbond,rat,fat,fstretch,iconnect,wold,&
!!                        alpha_stretch0,alpha_stretch)
!!!    use module_base, only: 8
!!    implicit none
!!    integer, intent(in) :: nat
!!    integer, intent(in) :: nbond
!!    real(8), intent(in) :: rat(3,nat)
!!    real(8), intent(inout) :: fat(3,nat)
!!    real(8), intent(inout) :: fstretch(3,nat)
!!    integer, intent(in) :: iconnect(2,nbond)
!!    real(8), intent(inout) :: wold(nbond)
!!    real(8), intent(in) :: alpha_stretch0
!!    real(8), intent(inout) :: alpha_stretch
!!    !internal
!!    integer :: iat,jat,ibond,jbond,l,nsame,info
!!    real(8) :: ss(nbond,nbond),w(nbond),vv(3,nat,nbond)
!!    real(8) :: per
!!    !functions
!!    real(8) :: ddot
!!    
!!    
!!    fstretch=0.0d0
!!
!!    !|v_i> := |rat_k>-|rat_l>
!!    !|F>=sum_i c_i* |v_i>
!!    !<v_j|F> = sum_i c_i <v_j|v_i>
!!
!!    ! set up positional overlap matrix
!!    vv=0.0d0
!!    do ibond=1,nbond
!!        iat=iconnect(1,ibond)
!!        jat=iconnect(2,ibond)
!!        do l=1,3
!!            vv(l,iat,ibond)=rat(l,jat)-rat(l,iat)
!!            vv(l,jat,ibond)=rat(l,iat)-rat(l,jat)
!!        enddo
!!    enddo
!!    
!!    ss=0.0d0
!!    w=0.0d0
!!    do ibond=1,nbond
!!        do jbond=1,nbond
!!            ss(ibond,jbond)=&
!!                       ddot(3*nat,vv(1,1,ibond),1,vv(1,1,jbond),1)
!!        enddo
!!        w(ibond)=ddot(3*nat,vv(1,1,ibond),1,fat(1,1),1)
!!    enddo
!!    
!!    nsame=0
!!    do ibond=1,nbond
!!        if ( wold(ibond)*w(ibond).gt.0.0d0) nsame=nsame+1
!!        wold(ibond)=w(ibond)
!!    enddo
!!    !determine feedback on streching components of force
!!    per=real(nsame,8)/nbond
!!    if (per.gt. .66d0) then
!!        alpha_stretch=alpha_stretch*1.10d0
!!    else
!!        alpha_stretch=max(1.d-2*alpha_stretch0,&
!!                           alpha_stretch/1.10d0)
!!    endif
!!
!!    call DPOSV('L', nbond, 1, ss, nbond, w, nbond, info )
!!    if (info.ne.0) then
!!        write(*,*)'info',info
!!        stop 'info DPOSV in minenergyforces'
!!    endif
!!
!!    ! calculate projected force
!!    fstretch=0.0d0
!!    do ibond=1,nbond
!!        do iat=1,nat
!!            do l=1,3
!!                fstretch(l,iat)=fstretch(l,iat)+w(ibond)*&
!!                                vv(l,iat,ibond)
!!            enddo
!!        enddo
!!    enddo
!!
!!    !
!!    do iat=1,nat
!!        do l=1,3
!!            fat(l,iat)=fat(l,iat)-fstretch(l,iat)
!!        enddo
!!    enddo
!!end subroutine
!!
end module
