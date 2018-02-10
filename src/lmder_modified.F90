!*****************************************************************************************
subroutine init_lmder_modified(parlm,m,ldfjac)
    !use mod_interface
    use mod_parlm, only: typ_parlm
    implicit none
    type(typ_parlm), intent(inout):: parlm
    integer, intent(in):: m, ldfjac
    if(parlm%n<1) stop 'ERROR: n<1'
    if(m<parlm%n) stop 'ERROR: m<n'
    if(ldfjac<m) stop 'ERROR: ldfjac<m'
    allocate(parlm%x(parlm%n))
    allocate(parlm%fvec(m))
    allocate(parlm%fjac(ldfjac,parlm%n))
    allocate(parlm%diag(parlm%n))
    allocate(parlm%ipvt(parlm%n))
    allocate(parlm%qtf(parlm%n))
    allocate(parlm%wa1(parlm%n))
    allocate(parlm%wa2(parlm%n))
    allocate(parlm%wa3(parlm%n))
    allocate(parlm%wa4(m))
    parlm%mode=1
    parlm%icontinue=0
    parlm%finish=.false.
end subroutine init_lmder_modified
!*****************************************************************************************
subroutine final_lmder_modified(parlm)
    !use mod_interface
    use mod_parlm, only: typ_parlm
    implicit none
    type(typ_parlm), intent(inout):: parlm
    deallocate(parlm%x)
    deallocate(parlm%fvec)
    deallocate(parlm%fjac)
    deallocate(parlm%diag)
    deallocate(parlm%ipvt)
    deallocate(parlm%qtf)
    deallocate(parlm%wa1)
    deallocate(parlm%wa2)
    deallocate(parlm%wa3)
    deallocate(parlm%wa4)
end subroutine final_lmder_modified
!*****************************************************************************************
subroutine lmder_modified(parlm,m,ldfjac)
    !use mod_interface
    use mod_parlm, only: typ_parlm
    implicit none
      type(typ_parlm), intent(inout):: parlm
      integer:: m,ldfjac
!     **********
!
!     subroutine lmder
!
!     the purpose of lmder is to minimize the sum of the squares of
!     m nonlinear functions in n variables by a modification of
!     the levenberg-marquardt algorithm. the user must provide a
!     subroutine which calculates the functions and the jacobian.
!
!     the subroutine statement is
!
!       subroutine lmder(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,
!                        maxfev,diag,mode,factor,nprint,info,nfev,
!                        njev,ipvt,qtf,wa1,wa2,wa3,wa4)
!
!     where
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions and the jacobian. fcn must
!         be declared in an external statement in the user
!         calling program, and should be written as follows.
!
!         subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
!         integer m,n,ldfjac,iflag
!         double precision x(n),fvec(m),fjac(ldfjac,n)
!         ----------
!         if iflag = 1 calculate the functions at x and
!         return this vector in fvec. do not alter fjac.
!         if iflag = 2 calculate the jacobian at x and
!         return this matrix in fjac. do not alter fvec.
!         ----------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of lmder.
!         in this case set iflag to a negative integer.
!
!       m is a positive integer input variable set to the number
!         of functions.
!
!       n is a positive integer input variable set to the number
!         of variables. n must not exceed m.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length m which contains
!         the functions evaluated at the output x.
!
!       fjac is an output m by n array. the upper n by n submatrix
!         of fjac contains an upper triangular matrix r with
!         diagonal elements of nonincreasing magnitude such that
!
!                t     t           t
!               p *(jac *jac)*p = r *r,
!
!         where p is a permutation matrix and jac is the final
!         calculated jacobian. column j of p is column ipvt(j)
!         (see below) of the identity matrix. the lower trapezoidal
!         part of fjac contains information generated during
!         the computation of r.
!
!       ldfjac is a positive integer input variable not less than m
!         which specifies the leading dimension of the array fjac.
!
!       ftol is a nonnegative input variable. termination
!         occurs when both the actual and predicted relative
!         reductions in the sum of squares are at most ftol.
!         therefore, ftol measures the relative error desired
!         in the sum of squares.
!
!       xtol is a nonnegative input variable. termination
!         occurs when the relative error between two consecutive
!         iterates is at most xtol. therefore, xtol measures the
!         relative error desired in the approximate solution.
!
!       gtol is a nonnegative input variable. termination
!         occurs when the cosine of the angle between fvec and
!         any column of the jacobian is at most gtol in absolute
!         value. therefore, gtol measures the orthogonality
!         desired between the function vector and the columns
!         of the jacobian.
!
!       maxfev is a positive integer input variable. termination
!         occurs when the number of calls to fcn with iflag = 1
!         has reached maxfev.
!
!       diag is an array of length n. if mode = 1 (see
!         below), diag is internally set. if mode = 2, diag
!         must contain positive entries that serve as
!         multiplicative scale factors for the variables.
!
!       mode is an integer input variable. if mode = 1, the
!         variables will be scaled internally. if mode = 2,
!         the scaling is specified by the input diag. other
!         values of mode are equivalent to mode = 1.
!
!       factor is a positive input variable used in determining the
!         initial step bound. this bound is set to the product of
!         factor and the euclidean norm of diag*x if nonzero, or else
!         to factor itself. in most cases factor should lie in the
!         interval (.1,100.).100. is a generally recommended value.
!
!       nprint is an integer input variable that enables controlled
!         printing of iterates if it is positive. in this case,
!         fcn is called with iflag = 0 at the beginning of the first
!         iteration and every nprint iterations thereafter and
!         immediately prior to return, with x, fvec, and fjac
!         available for printing. fvec and fjac should not be
!         altered. if nprint is not positive, no special calls
!         of fcn with iflag = 0 are made.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0  improper input parameters.
!
!         info = 1  both actual and predicted relative reductions
!                   in the sum of squares are at most ftol.
!
!         info = 2  relative error between two consecutive iterates
!                   is at most xtol.
!
!         info = 3  conditions for info = 1 and info = 2 both hold.
!
!         info = 4  the cosine of the angle between fvec and any
!                   column of the jacobian is at most gtol in
!                   absolute value.
!
!         info = 5  number of calls to fcn with iflag = 1 has
!                   reached maxfev.
!
!         info = 6  ftol is too small. no further reduction in
!                   the sum of squares is possible.
!
!         info = 7  xtol is too small. no further improvement in
!                   the approximate solution x is possible.
!
!         info = 8  gtol is too small. fvec is orthogonal to the
!                   columns of the jacobian to machine precision.
!
!       nfev is an integer output variable set to the number of
!         calls to fcn with iflag = 1.
!
!       njev is an integer output variable set to the number of
!         calls to fcn with iflag = 2.
!
!       ipvt is an integer output array of length n. ipvt
!         defines a permutation matrix p such that jac*p = q*r,
!         where jac is the final calculated jacobian, q is
!         orthogonal (not stored), and r is upper triangular
!         with diagonal elements of nonincreasing magnitude.
!         column j of p is column ipvt(j) of the identity matrix.
!
!       qtf is an output array of length n which contains
!         the first n elements of the vector (q transpose)*fvec.
!
!       wa1, wa2, and wa3 are work arrays of length n.
!
!       wa4 is a work array of length m.
!
!     subprograms called
!
!       user-supplied ...... fcn
!
!       minpack-supplied ... dpmpar,enorm,lmpar,qrfac
!
!       fortran-supplied ... dabs,dmax1,dmin1,dsqrt,mod
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      integer i,j,l
      double precision dirder,fnorm1, &
                       one,p1,p5,p25,p75,p0001, &
                       sum,temp,temp1,temp2,zero
      double precision dpmpar,enorm
      data one,p1,p5,p25,p75,p0001,zero &
           /1.0d0,1.0d-1,5.0d-1,2.5d-1,7.5d-1,1.0d-4,0.0d0/
      if(parlm%icontinue==400) then
        parlm%icontinue=0
        goto 400
      endif
      if(parlm%icontinue==500) then
        parlm%icontinue=0
        goto 500
      endif
      if(parlm%icontinue==600) then
        parlm%icontinue=0
        goto 600
      endif
      if(parlm%icontinue==700) then
        parlm%icontinue=0
        goto 700
      endif
      if(parlm%icontinue==800) then
        parlm%icontinue=0
        goto 800
      endif
!
!     epsmch is the machine precision.
!
      parlm%epsmch = dpmpar(1)
!
      parlm%info = 0
      parlm%iflag = 0
      parlm%nfev = 0
      parlm%njev = 0
!
!     check the input parameters for errors.
!
      if (parlm%n .le. 0 .or. m .lt. parlm%n .or. ldfjac .lt. m &
          .or. parlm%ftol .lt. zero .or. parlm%xtol .lt. zero .or. parlm%gtol .lt. zero &
          .or. parlm%maxfev .le. 0 .or. parlm%factor .le. zero) go to 300
      if (parlm%mode==2) then
      do j = 1, parlm%n
         if (parlm%diag(j) .le. zero) go to 300
      enddo
      endif
!
!     evaluate the function at the starting point
!     and calculate its norm.
!
      parlm%iflag = 1
      if(parlm%icontinue/=400) then
        parlm%icontinue=400
        goto 900
      endif
  400 continue
      parlm%nfev = 1
      if (parlm%iflag .lt. 0) go to 300
      parlm%fnorm = enorm(m,parlm%fvec)
!
!     initialize levenberg-marquardt parameter and iteration counter.
!
      parlm%par = zero
      parlm%iter = 1
!
!     beginning of the outer loop.
!
   30 continue
!
!        calculate the jacobian matrix.
!
         parlm%iflag = 2
      if(parlm%icontinue/=500) then
        parlm%icontinue=500
        goto 900
      endif
  500 continue
         parlm%njev = parlm%njev + 1
         if (parlm%iflag .lt. 0) go to 300
!
!        if requested, call fcn to enable printing of iterates.
!
         if (parlm%nprint .le. 0) go to 40
         parlm%iflag = 0
         if (mod(parlm%iter-1,parlm%nprint) .eq. 0) then
      if(parlm%icontinue/=600) then
        parlm%icontinue=600
        goto 900
      endif
         endif
  600 continue
         if (parlm%iflag .lt. 0) go to 300
   40    continue
!
!        compute the qr factorization of the jacobian.
!
         call qrfac(m,parlm%n,parlm%fjac,ldfjac,.true.,parlm%ipvt,parlm%n,parlm%wa1,parlm%wa2,parlm%wa3)
!
!        on the first iteration and if mode is 1, scale according
!        to the norms of the columns of the initial jacobian.
!
         if (parlm%iter==1) then
         if (parlm%mode/=2) then
         do j = 1, parlm%n
            parlm%diag(j) = parlm%wa2(j)
            if (parlm%wa2(j) .eq. zero) parlm%diag(j) = one
         enddo
         endif
!
!        on the first iteration, calculate the norm of the scaled x
!        and initialize the step bound delta.
!
         do j = 1, parlm%n
            parlm%wa3(j) = parlm%diag(j)*parlm%x(j)
         enddo
         parlm%xnorm = enorm(parlm%n,parlm%wa3)
         parlm%delta = parlm%factor*parlm%xnorm
         if (parlm%delta .eq. zero) parlm%delta = parlm%factor
         endif
!
!        form (q transpose)*fvec and store the first n components in
!        qtf.
!
         do i = 1, m
            parlm%wa4(i) = parlm%fvec(i)
         enddo
         do j = 1, parlm%n
            if (parlm%fjac(j,j)/=zero) then
            sum = zero
            do i = j, m
               sum = sum + parlm%fjac(i,j)*parlm%wa4(i)
            enddo
            temp = -sum/parlm%fjac(j,j)
            do i = j, m
               parlm%wa4(i) = parlm%wa4(i) + parlm%fjac(i,j)*temp
            enddo
            endif
            parlm%fjac(j,j) = parlm%wa1(j)
            parlm%qtf(j) = parlm%wa4(j)
         enddo
!
!        compute the norm of the scaled gradient.
!
         parlm%gnorm = zero
         if (parlm%fnorm/=zero) then
         do j = 1, parlm%n
            l = parlm%ipvt(j)
            if (parlm%wa2(l)/=zero) then
            sum = zero
            do i = 1, j
               sum = sum + parlm%fjac(i,j)*(parlm%qtf(i)/parlm%fnorm)
            enddo
            parlm%gnorm = dmax1(parlm%gnorm,dabs(sum/parlm%wa2(l)))
            endif
         enddo
         endif
!
!        test for convergence of the gradient norm.
!
         if (parlm%gnorm .le. parlm%gtol) parlm%info = 4
         if (parlm%info .ne. 0) go to 300
!
!        rescale if necessary.
!
         if (parlm%mode/=2) then
         do j = 1, parlm%n
            parlm%diag(j) = dmax1(parlm%diag(j),parlm%wa2(j))
         enddo
         endif
!
!        beginning of the inner loop.
!
  200    continue
!
!           determine the levenberg-marquardt parameter.
!
            call lmpar(parlm%n,parlm%fjac,ldfjac,parlm%ipvt,parlm%diag,parlm%qtf,parlm%delta,parlm%par,parlm%wa1,parlm%wa2, &
                       parlm%wa3,parlm%wa4)
!
!           store the direction p and x + p. calculate the norm of p.
!
            do j = 1, parlm%n
               parlm%wa1(j) = -parlm%wa1(j)
               parlm%wa2(j) = parlm%x(j) + parlm%wa1(j)
               parlm%wa3(j) = parlm%diag(j)*parlm%wa1(j)
            enddo
            parlm%pnorm = enorm(parlm%n,parlm%wa3)
!
!           on the first iteration, adjust the initial step bound.
!
            if (parlm%iter .eq. 1) parlm%delta = dmin1(parlm%delta,parlm%pnorm)
!
!           evaluate the function at x + p and calculate its norm.
!
            parlm%iflag = 1
      if(parlm%icontinue/=700) then
        parlm%icontinue=700
        goto 900
      endif
  700 continue
            parlm%nfev = parlm%nfev + 1
            if (parlm%iflag .lt. 0) go to 300
            fnorm1 = enorm(m,parlm%wa4)
!
!           compute the scaled actual reduction.
!
            parlm%actred = -one
            if (p1*fnorm1 .lt. parlm%fnorm) parlm%actred = one - (fnorm1/parlm%fnorm)**2
!
!           compute the scaled predicted reduction and
!           the scaled directional derivative.
!
            do j = 1, parlm%n
               parlm%wa3(j) = zero
               l = parlm%ipvt(j)
               temp = parlm%wa1(l)
               do i = 1, j
                  parlm%wa3(i) = parlm%wa3(i) + parlm%fjac(i,j)*temp
            enddo
            enddo
            temp1 = enorm(parlm%n,parlm%wa3)/parlm%fnorm
            temp2 = (dsqrt(parlm%par)*parlm%pnorm)/parlm%fnorm
            parlm%prered = temp1**2 + temp2**2/p5
            dirder = -(temp1**2 + temp2**2)
!
!           compute the ratio of the actual to the predicted
!           reduction.
!
            parlm%ratio = zero
            if (parlm%prered .ne. zero) parlm%ratio = parlm%actred/parlm%prered
!
!           update the step bound.
!
            if (.not. (parlm%ratio .gt. p25)) then
               if (parlm%actred .ge. zero) temp = p5
               if (parlm%actred .lt. zero) &
                  temp = p5*dirder/(dirder + p5*parlm%actred)
               if (p1*fnorm1 .ge. parlm%fnorm .or. temp .lt. p1) temp = p1
               parlm%delta = temp*dmin1(parlm%delta,parlm%pnorm/p1)
               parlm%par = parlm%par/temp
            elseif (.not. (parlm%par .ne. zero .and. parlm%ratio .lt. p75)) then
               parlm%delta = parlm%pnorm/p5
               parlm%par = p5*parlm%par
               endif
!
!           test for successful iteration.
!
            if (.not. (parlm%ratio<p0001)) then
!
!           successful iteration. update x, fvec, and their norms.
!
            do j = 1, parlm%n
               parlm%x(j) = parlm%wa2(j)
               parlm%wa2(j) = parlm%diag(j)*parlm%x(j)
            enddo
            do i = 1, m
               parlm%fvec(i) = parlm%wa4(i)
            enddo
            parlm%xnorm = enorm(parlm%n,parlm%wa2)
            parlm%fnorm = fnorm1
            parlm%iter = parlm%iter + 1
            endif
!
!           tests for convergence.
!
            if (dabs(parlm%actred) .le. parlm%ftol .and. parlm%prered .le. parlm%ftol &
                .and. p5*parlm%ratio .le. one) parlm%info = 1
            if (parlm%delta .le. parlm%xtol*parlm%xnorm) parlm%info = 2
            if (dabs(parlm%actred) .le. parlm%ftol .and. parlm%prered .le. parlm%ftol &
                .and. p5*parlm%ratio .le. one .and. parlm%info .eq. 2) parlm%info = 3
            if (parlm%info .ne. 0) go to 300
!
!           tests for termination and stringent tolerances.
!
            if (parlm%nfev .ge. parlm%maxfev) parlm%info = 5
            if (dabs(parlm%actred) .le. parlm%epsmch .and. parlm%prered .le. parlm%epsmch &
                .and. p5*parlm%ratio .le. one) parlm%info = 6
            if (parlm%delta .le. parlm%epsmch*parlm%xnorm) parlm%info = 7
            if (parlm%gnorm .le. parlm%epsmch) parlm%info = 8
            if (parlm%info .ne. 0) go to 300
!
!           end of the inner loop. repeat if iteration unsuccessful.
!
            if (parlm%ratio .lt. p0001) go to 200
!
!        end of the outer loop.
!
         go to 30
  300 continue
!
!     termination, either normal or user imposed.
!
      if (parlm%iflag .lt. 0) parlm%info = parlm%iflag
      parlm%iflag = 0
      if (parlm%nprint .gt. 0) then
      if(parlm%icontinue/=800) then
        parlm%icontinue=800
        goto 900
      endif
      endif
  800 continue
!
!     last card of subroutine lmder.
!
      parlm%finish=.true.
  900 continue
      end
