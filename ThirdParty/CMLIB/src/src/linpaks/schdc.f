 
      SUBROUTINE SCHDC(A,LDA,P,WORK,JPVT,JOB,INFO)
C***BEGIN PROLOGUE  SCHDC
C***DATE WRITTEN   790319   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D2B1B
C***KEYWORDS  CHOLESKY DECOMPOSITION,LINEAR ALGEBRA,LINPACK,MATRIX,
C             POSITIVE DEFINITE
C***AUTHOR  DONGARRA, J., (ANL)
C           STEWART, G. W., (U. OF MARYLAND)
C***PURPOSE  Computes the Cholesky decomposition of A POSITIVE DEFINITE
C            matrix.  A pivoting option allows the user to estimate the
C            condition of a POSITIVE DEFINITE matrix or determine the
C            rank of a POSITIVE SEMIDEFINITE matrix.
C***DESCRIPTION
C
C     SCHDC computes the Cholesky decomposition of a positive definite
C     matrix.  A pivoting option allows the user to estimate the
C     condition of a positive definite matrix or determine the rank
C     of a positive semidefinite matrix.
C
C     On Entry
C
C         A      REAL(LDA,P).
C                A contains the matrix whose decomposition is to
C                be computed.  Only the upper half of A need be stored.
C                The lower part of the array A is not referenced.
C
C         LDA    INTEGER.
C                LDA is the leading dimension of the array A.
C
C         P      INTEGER.
C                P is the order of the matrix.
C
C         WORK   REAL.
C                WORK is a work array.
C
C         JPVT   INTEGER(P).
C                JPVT contains integers that control the selection
C                of the pivot elements, if pivoting has been requested.
C                Each diagonal element A(K,K)
C                is placed in one of three classes according to the
C                value of JPVT(K).
C
C                   If JPVT(K) .GT. 0, then X(K) is an initial
C                                      element.
C
C                   If JPVT(K) .EQ. 0, then X(K) is a free element.
C
C                   If JPVT(K) .LT. 0, then X(K) is a final element.
C
C                Before the decomposition is computed, initial elements
C                are moved by symmetric row and column interchanges to
C                the beginning of the array A and final
C                elements to the end.  Both initial and final elements
C                are frozen in place during the computation and only
C                free elements are moved.  At the K-th stage of the
C                reduction, if A(K,K) is occupied by a free element
C                it is interchanged with the largest free element
C                A(L,L) with L .GE. K.  JPVT is not referenced if
C                JOB .EQ. 0.
C
C        JOB     INTEGER.
C                JOB is an integer that initiates column pivoting.
C                If JOB .EQ. 0, no pivoting is done.
C                If JOB .NE. 0, pivoting is done.
C
C     On Return
C
C         A      A contains in its upper half the Cholesky factor
C                of the matrix A as it has been permuted by pivoting.
C
C         JPVT   JPVT(J) contains the index of the diagonal element
C                of a that was moved into the J-th position,
C                provided pivoting was requested.
C
C         INFO   contains the index of the last positive diagonal
C                element of the Cholesky factor.
C
C     For positive definite matrices INFO = P is the normal return.
C     For pivoting with positive semidefinite matrices INFO will
C     in general be less than P.  However, INFO may be greater than
C     the rank of A, since rounding error can cause an otherwise zero
C     element to be positive.  Indefinite systems will always cause
C     INFO to be less than P.
C
C     LINPACK.  This version dated 03/19/79 .
C     J. J. Dongarra and G. W. Stewart, Argonne National Laboratory and
C     University of Maryland.
C
C
C     BLAS SAXPY,SSWAP
C     Fortran SQRT
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  SAXPY,SSWAP
C***END PROLOGUE  SCHDC
      INTEGER LDA,P,JPVT(*),JOB,INFO
      REAL A(LDA,*),WORK(*)
C
      INTEGER PU,PL,PLP1,I,J,JP,JT,K,KB,KM1,KP1,L,MAXL
      REAL TEMP
      REAL MAXDIA
      LOGICAL SWAPK,NEGK
C***FIRST EXECUTABLE STATEMENT  SCHDC
      PL = 1
      PU = 0
      INFO = P
      IF (JOB .EQ. 0) GO TO 160
C
C        PIVOTING HAS BEEN REQUESTED. REARRANGE THE
C        THE ELEMENTS ACCORDING TO JPVT.
C
         DO 70 K = 1, P
            SWAPK = JPVT(K) .GT. 0
            NEGK = JPVT(K) .LT. 0
            JPVT(K) = K
            IF (NEGK) JPVT(K) = -JPVT(K)
            IF (.NOT.SWAPK) GO TO 60
               IF (K .EQ. PL) GO TO 50
                  CALL SSWAP(PL-1,A(1,K),1,A(1,PL),1)
                  TEMP = A(K,K)
                  A(K,K) = A(PL,PL)
                  A(PL,PL) = TEMP
                  PLP1 = PL + 1
                  IF (P .LT. PLP1) GO TO 40
                  DO 30 J = PLP1, P
                     IF (J .GE. K) GO TO 10
                        TEMP = A(PL,J)
                        A(PL,J) = A(J,K)
                        A(J,K) = TEMP
                     GO TO 20
   10                CONTINUE
                     IF (J .EQ. K) GO TO 20
                        TEMP = A(K,J)
                        A(K,J) = A(PL,J)
                        A(PL,J) = TEMP
   20                CONTINUE
   30             CONTINUE
   40             CONTINUE
                  JPVT(K) = JPVT(PL)
                  JPVT(PL) = K
   50          CONTINUE
               PL = PL + 1
   60       CONTINUE
   70    CONTINUE
         PU = P
         IF (P .LT. PL) GO TO 150
         DO 140 KB = PL, P
            K = P - KB + PL
            IF (JPVT(K) .GE. 0) GO TO 130
               JPVT(K) = -JPVT(K)
               IF (PU .EQ. K) GO TO 120
                  CALL SSWAP(K-1,A(1,K),1,A(1,PU),1)
                  TEMP = A(K,K)
                  A(K,K) = A(PU,PU)
                  A(PU,PU) = TEMP
                  KP1 = K + 1
                  IF (P .LT. KP1) GO TO 110
                  DO 100 J = KP1, P
                     IF (J .GE. PU) GO TO 80
                        TEMP = A(K,J)
                        A(K,J) = A(J,PU)
                        A(J,PU) = TEMP
                     GO TO 90
   80                CONTINUE
                     IF (J .EQ. PU) GO TO 90
                        TEMP = A(K,J)
                        A(K,J) = A(PU,J)
                        A(PU,J) = TEMP
   90                CONTINUE
  100             CONTINUE
  110             CONTINUE
                  JT = JPVT(K)
                  JPVT(K) = JPVT(PU)
                  JPVT(PU) = JT
  120          CONTINUE
               PU = PU - 1
  130       CONTINUE
  140    CONTINUE
  150    CONTINUE
  160 CONTINUE
      DO 270 K = 1, P
C
C        REDUCTION LOOP.
C
         MAXDIA = A(K,K)
         KP1 = K + 1
         MAXL = K
C
C        DETERMINE THE PIVOT ELEMENT.
C
         IF (K .LT. PL .OR. K .GE. PU) GO TO 190
            DO 180 L = KP1, PU
               IF (A(L,L) .LE. MAXDIA) GO TO 170
                  MAXDIA = A(L,L)
                  MAXL = L
  170          CONTINUE
  180       CONTINUE
  190    CONTINUE
C
C        QUIT IF THE PIVOT ELEMENT IS NOT POSITIVE.
C
         IF (MAXDIA .GT. 0.0E0) GO TO 200
            INFO = K - 1
C     ......EXIT
            GO TO 280
  200    CONTINUE
         IF (K .EQ. MAXL) GO TO 210
C
C           START THE PIVOTING AND UPDATE JPVT.
C
            KM1 = K - 1
            CALL SSWAP(KM1,A(1,K),1,A(1,MAXL),1)
            A(MAXL,MAXL) = A(K,K)
            A(K,K) = MAXDIA
            JP = JPVT(MAXL)
            JPVT(MAXL) = JPVT(K)
            JPVT(K) = JP
  210    CONTINUE
C
C        REDUCTION STEP. PIVOTING IS CONTAINED ACROSS THE ROWS.
C
         WORK(K) = SQRT(A(K,K))
         A(K,K) = WORK(K)
         IF (P .LT. KP1) GO TO 260
         DO 250 J = KP1, P
            IF (K .EQ. MAXL) GO TO 240
               IF (J .GE. MAXL) GO TO 220
                  TEMP = A(K,J)
                  A(K,J) = A(J,MAXL)
                  A(J,MAXL) = TEMP
               GO TO 230
  220          CONTINUE
               IF (J .EQ. MAXL) GO TO 230
                  TEMP = A(K,J)
                  A(K,J) = A(MAXL,J)
                  A(MAXL,J) = TEMP
  230          CONTINUE
  240       CONTINUE
            A(K,J) = A(K,J)/WORK(K)
            WORK(J) = A(K,J)
            TEMP = -A(K,J)
            CALL SAXPY(J-K,TEMP,WORK(KP1),1,A(KP1,J),1)
  250    CONTINUE
  260    CONTINUE
  270 CONTINUE
  280 CONTINUE
      RETURN
      END
