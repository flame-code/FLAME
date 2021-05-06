!SUBROUTINE APC(N,A,F,Z)
! Modified by Ali Sadeghi to get real*8 matrix A(N,N) and converted to F90 
!
! SOLUTION OF THE LINEAR MIN-SUM ASSIGNMENT PROBLEM.
! HUNGARIAN METHOD. COMPLEXITY O(N**3).
!
! MEANING OF THE INPUT PARAMETERS:
! N      = NUMBER OF ROWS AND COLUMNS OF THE COST MATRIX.
! A(I,J) = COST OF THE ASSIGNMENT OF ROW  I  TO COLUMN  J .
! ON RETURN, THE INPUT PARAMETERS ARE UNCHANGED.
!
! MEANING OF THE OUTPUT PARAMETERS:
! F(I) = COLUMN ASSIGNED TO ROW  I .
! Z    = COST OF THE OPTIMAL ASSIGNMENT =
!      = A(1,F(1)) + A(2,F(2)) + ... + A(N,F(N)) .
!
!
! THE CODE IS BASED ON THE HUNGARIAN METHOD AS DESCRIBED BY
! LAWLER (COMBINATORIAL OPTIMIZATION : NETWORKS AND
! MATROIDS, HOLT, RINEHART AND WINSTON, NEW YORK, 1976).
! THE ALGORITHMIC PASCAL-LIKE DESCRIPTION OF THE CODE IS
! GIVEN IN G.CARPANETO, S.MARTELLO AND P.TOTH, ALGORITHMS AND
! CODES FOR THE ASSIGNMENT PROBLEM, ANNALS OF OPERATIONS
! RESEARCH 7, 1988.
!
! SUBROUTINE APC DETERMINES THE INITIAL DUAL AND PARTIAL
! PRIMAL SOLUTIONS AND THEN SEARCHES FOR AUGMENTING PATHS
! UNTIL ALL ROWS AND COLUMNS ARE ASSIGNED.
!
! MEANING OF THE MAIN INTERNAL VARIABLES:
! FB(J) = ROW ASSIGNED TO COLUMN  J .
! M     = NUMBER OF INITIAL ASSIGNMENTS.
! U(I)  = DUAL VARIABLE ASSOCIATED WITH ROW  I .
! V(J)  = DUAL VARIABLE ASSOCIATED WITH COLUMN  J .
!
! APC NEEDS THE FOLLOWING SUBROUTINES: INCR_inalborz
!                                      INIT_inalborz
!                                      PATH_inalborz
!
! THIS WORK WAS SUPPORTED BY  C.N.R. , ITALY.
subroutine hung(N,A,F,Z)
    implicit none
      integer:: n
      real(8)::  A(n,n),Z,U(n),V(n)
      integer F(N),FB(n), RC(n)
      integer:: M,I,J
! SEARCH FOR THE INITIAL DUAL AND PARTIAL PRIMAL SOLUTIONS.
      CALL INIT_inalborz(N,A,F,M,U,V,FB,RC)
      IF ( M .NE. N ) then 
! SOLUTION OF THE REDUCED PROBLEM.
      DO  I=1,N
        IF ( F(I) == 0 ) THEN 
! DETERMINATION OF AN AUGMENTING PATH STARTING FROM ROW  I .
        CALL PATH_inalborz(N,A,I,F,J,U,V,FB,RC)
! ASSIGNMENT OF ROW  I  AND COLUMN  J .
        CALL INCR_inalborz(n,F,J,FB,RC)
        ENDIF
      ENDDO    
      ENDIF
! COMPUTATION OF THE SOLUTION COST  Z .
      Z = sum(u(1:N)) + sum(V(1:N))
end subroutinE hung
!********************************************************************************
subroutine INCR_inalborz(n,F,J,FB,RC)
implicit none
! ASSIGNMENT OF COLUMN  J .
      integer:: n,I,J,JJ,  F(n),FB(n),RC(n)
   10 I = RC(J)
      FB(J) = I
      JJ = F(I)
      F(I) = J
      J = JJ
      IF ( J > 0 ) GO TO 10
      RETURN
end subroutine INCR_inalborz
!********************************************************************************
subroutine INIT_inalborz(N,A,F,M,U,V,FB,P)
! SEARCH FOR THE INITIAL DUAL AND PARTIAL PRIMAL SOLUTIONS.
! P(I) = FIRST UNSCANNED COLUMN OF ROW  I .
    implicit none
      integer:: n,m, F(n),FB(n),P(n)
      real(8) A(n,n) , U(n),V(n)
      real(8), parameter :: INF = 1.d9
      real(8) min, IA
      integer i,j, k,R, JMIN, KK
! PHASE 1 .
      M = 0
      F(1:N)=0
      FB(1:N)=0
! SCANNING OF THE COLUMNS ( INITIALIZATION OF  V(J) ).
      DO 40 J=1,N
        MIN = INF
        DO 30 I=1,N
          IA = A(I,J)
          IF ( IA .GT. MIN ) GO TO 30
          IF ( IA .LT. MIN ) GO TO 20
          IF ( F(I) .NE. 0 ) GO TO 30
   20     MIN = IA
          R = I
   30   CONTINUE
        V(J) = MIN
        IF ( F(R) .NE. 0 ) GO TO 40
! ASSIGNMENT OF COLUMN  J  TO ROW  R .
        M = M + 1
        FB(J) = R
        F(R) = J
        U(R) = 0.d0
        P(R) = J + 1
   40 CONTINUE
! PHASE 2 .
! SCANNING OF THE UNASSIGNED ROWS ( UPDATING OF  U(I) ).
      DO 110 I=1,N
        IF ( F(I) .NE. 0 ) GO TO 110
        MIN = INF
        DO 60 K=1,N
          IA = A(I,K) - V(K)
          IF ( IA .GT. MIN )  GO TO 60
          IF ( IA .LT. MIN )  GO TO 50
          IF ( FB(K) .NE. 0 ) GO TO 60
          IF ( FB(J) .EQ. 0 ) GO TO 60
   50     MIN = IA
          J = K
   60   CONTINUE
        U(I) = MIN
        JMIN = J
        IF ( FB(J) .EQ. 0 ) GO TO 100
        DO 80 J=JMIN,N
          IF ( A(I,J) - V(J) .GT. MIN ) GO TO 80
          R = FB(J)
          KK = P(R)
          IF ( KK .GT. N ) GO TO 80
          DO 70 K=KK,N
            IF ( FB(K) .GT. 0 ) GO TO 70
            IF ( A(R,K) - U(R) - V(K) .EQ. 0.d0 ) GO TO 90
   70     CONTINUE
          P(R) = N + 1
   80   CONTINUE
        GO TO 110
! REASSIGNMENT OF ROW  R  AND COLUMN  K .
   90   F(R) = K
        FB(K) = R
        P(R) = K + 1
! ASSIGNMENT OF COLUMN  J  TO ROW  I .
  100   M = M + 1
        F(I) = J
        FB(J)= I
        P(I) = J + 1
  110 CONTINUE
      RETURN
end subroutine INIT_inalborz
!********************************************************************************
subroutine PATH_inalborz(N,A,II,F,JJ,U,V,FB,RC)
! DETERMINATION OF AN AUGMENTING PATH STARTING FROM
! UNASSIGNED ROW  II  AND TERMINATING AT UNASSIGNED COLUMN
! JJ , WITH UPDATING OF DUAL VARIABLES  U(I)  AND  V(J) .
!
! MEANING OF THE MAIN INTERNAL VARIABLES:
! LR(L) = L-TH LABELLED ROW ( L=1,NLR ).
! PI(J) = MIN ( A(I,J) - U(I) - V(J) , SUCH THAT ROW  I  IS
!         LABELLED AND NOT EQUAL TO  FB(J) ).
! RC(J) = ROW PRECEDING COLUMN  J  IN THE CURRENT
!         ALTERNATING PATH.
! UC(L) = L-TH UNLABELLED COLUMN ( L=1,NUC ).
      implicit none
      integer:: N 
      real(8)::  A(n,n),U(n),V(N),PI(n), IA, MIN
      integer:: F(N),LR(n),UC(n)
      integer:: FB(n),RC(n)
      real(8), parameter :: INF = 1.d9
      integer::  i,j,k,L,ii,jj,NUC,NLR,R
! INITIALIZATION.
      LR(1) = II
      DO 10 K=1,N
        PI(K) = A(II,K) - U(II) - V(K)
        RC(K) = II
        UC(K) = K
   10 CONTINUE
      NUC = N
      NLR = 1
      GO TO 40
! SCANNING OF THE LABELLED ROWS.
   20 R = LR(NLR)
      DO 30 L=1,NUC
        J = UC(L)
        IA = A(R,J) - U(R) - V(J)
        IF ( IA .GE. PI(J) ) GO TO 30
        PI(J) = IA
        RC(J) = R
   30 CONTINUE
! SEARCH FOR A ZERO ELEMENT IN AN UNLABELLED COLUMN.
   40 DO 50 L=1,NUC
        J = UC(L)
        IF ( PI(J) .EQ. 0.d0 ) GO TO 100
   50 CONTINUE
! UPDATING OF THE DUAL VARIABLES  U(I)  AND  V(J) .
      MIN = INF
      DO 60 L=1,NUC
        J = UC(L)
        IF ( MIN .GT. PI(J) ) MIN = PI(J)
   60 CONTINUE
      DO 70 L=1,NLR
        R = LR(L)
        U(R) = U(R) + MIN
   70 CONTINUE
      DO 90 J=1,N
        IF ( PI(J) .EQ. 0.d0 ) GO TO 80
        PI(J) = PI(J) - MIN
        GO TO 90
   80   V(J) = V(J) - MIN
   90 CONTINUE
      GO TO 40
  100 IF ( FB(J) .EQ. 0 ) GO TO 110
! LABELLING OF ROW  FB(J)  AND REMOVAL OF THE LABEL  OF COLUMN  J .
      NLR = NLR + 1
      LR(NLR) = FB(J)
      UC(L) = UC(NUC)
      NUC = NUC - 1
      GO TO 20
! DETERMINATION OF THE UNASSIGNED COLUMN  J .
  110 JJ = J
      RETURN
end subroutine PATH_inalborz
!******************************************************************************** 
