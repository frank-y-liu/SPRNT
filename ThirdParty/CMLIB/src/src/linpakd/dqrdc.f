      SUBROUTINE DQRDC(X,LDX,N,P,QRAUX,JPVT,WORK,JOB)
C***BEGIN PROLOGUE  DQRDC
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D5
C***KEYWORDS  DECOMPOSITION,DOUBLE PRECISION,LINEAR ALGEBRA,LINPACK,
C             MATRIX,ORTHOGONAL TRIANGULAR
C***AUTHOR  STEWART, G. W., (U. OF MARYLAND)
C***PURPOSE  Uses Householder transformations to compute the Qr factori-
C            zation of N by P matrix X.  Column pivoting is optional.
C***DESCRIPTION
C
C     DQRDC uses Householder transformations to compute the QR
C     factorization of an N by P matrix X.  Column pivoting
C     based on the 2-norms of the reduced columns may be
C     performed at the user's option.
C
C     On Entry
C
C        X       DOUBLE PRECISION(LDX,P), where LDX .GE. N.
C                X contains the matrix whose decomposition is to be
C                computed.
C
C        LDX     INTEGER.
C                LDX is the leading dimension of the array X.
C
C        N       INTEGER.
C                N is the number of rows of the matrix X.
C
C        P       INTEGER.
C                P is the number of columns of the matrix X.
C
C        JPVT    INTEGER(P).
C                JPVT contains integers that control the selection
C                of the pivot columns.  The K-th column X(K) of X
C                is placed in one of three classes according to the
C                value of JPVT(K).
C
C                   If JPVT(K) .GT. 0, then X(K) is an initial
C                                      column.
C
C                   If JPVT(K) .EQ. 0, then X(K) is a free column.
C
C                   If JPVT(K) .LT. 0, then X(K) is a final column.
C
C                Before the decomposition is computed, initial columns
C                are moved to the beginning of the array X and final
C                columns to the end.  Both initial and final columns
C                are frozen in place during the computation and only
C                free columns are moved.  At the K-th stage of the
C                reduction, if X(K) is occupied by a free column
C                it is interchanged with the free column of largest
C                reduced norm.  JPVT is not referenced if
C                JOB .EQ. 0.
C
C        WORK    DOUBLE PRECISION(P).
C                WORK is a work array.  WORK is not referenced if
C                JOB .EQ. 0.
C
C        JOB     INTEGER.
C                JOB is an integer that initiates column pivoting.
C                If JOB .EQ. 0, no pivoting is done.
C                If JOB .NE. 0, pivoting is done.
C
C     On Return
C
C        X       X contains in its upper triangle the upper
C                triangular matrix R of the QR factorization.
C                Below its diagonal X contains information from
C                which the orthogonal part of the decomposition
C                can be recovered.  Note that if pivoting has
C                been requested, the decomposition is not that
C                of the original matrix X but that of X
C                with its columns permuted as described by JPVT.
C
C        QRAUX   DOUBLE PRECISION(P).
C                QRAUX contains further information required to recover
C                the orthogonal part of the decomposition.
C
C        JPVT    JPVT(K) contains the index of the column of the
C                original matrix that has been interchanged into
C                the K-th column, if pivoting was requested.
C
C     LINPACK.  This version dated 08/14/78 .
C     G. W. Stewart, University of Maryland, Argonne National Lab.
C
C     DQRDC uses the following functions and subprograms.
C
C     BLAS DAXPY,DDOT,DSCAL,DSWAP,DNRM2
C     Fortran DABS,DMAX1,MIN0,DSQRT
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  DAXPY,DDOT,DNRM2,DSCAL,DSWAP
C***END PROLOGUE  DQRDC
      INTEGER LDX,N,P,JOB
      INTEGER JPVT(*)
      DOUBLE PRECISION X(LDX,*),QRAUX(*),WORK(*)
C
      INTEGER J,JP,L,LP1,LUP,MAXJ,PL,PU
      DOUBLE PRECISION MAXNRM,DNRM2,TT
      DOUBLE PRECISION DDOT,NRMXL,T
      LOGICAL NEGJ,SWAPJ
C
C***FIRST EXECUTABLE STATEMENT  DQRDC
      PL = 1
      PU = 0
      IF (JOB .EQ. 0) GO TO 60
C
C        PIVOTING HAS BEEN REQUESTED.  REARRANGE THE COLUMNS
C        ACCORDING TO JPVT.
C
         DO 20 J = 1, P
            SWAPJ = JPVT(J) .GT. 0
            NEGJ = JPVT(J) .LT. 0
            JPVT(J) = J
            IF (NEGJ) JPVT(J) = -J
            IF (.NOT.SWAPJ) GO TO 10
               IF (J .NE. PL) CALL DSWAP(N,X(1,PL),1,X(1,J),1)
               JPVT(J) = JPVT(PL)
               JPVT(PL) = J
               PL = PL + 1
   10       CONTINUE
   20    CONTINUE
         PU = P
         DO 50 JJ = 1, P
            J = P - JJ + 1
            IF (JPVT(J) .GE. 0) GO TO 40
               JPVT(J) = -JPVT(J)
               IF (J .EQ. PU) GO TO 30
                  CALL DSWAP(N,X(1,PU),1,X(1,J),1)
                  JP = JPVT(PU)
                  JPVT(PU) = JPVT(J)
                  JPVT(J) = JP
   30          CONTINUE
               PU = PU - 1
   40       CONTINUE
   50    CONTINUE
   60 CONTINUE
C
C     COMPUTE THE NORMS OF THE FREE COLUMNS.
C
      IF (PU .LT. PL) GO TO 80
      DO 70 J = PL, PU
         QRAUX(J) = DNRM2(N,X(1,J),1)
         WORK(J) = QRAUX(J)
   70 CONTINUE
   80 CONTINUE
C
C     PERFORM THE HOUSEHOLDER REDUCTION OF X.
C
      LUP = MIN0(N,P)
      DO 200 L = 1, LUP
         IF (L .LT. PL .OR. L .GE. PU) GO TO 120
C
C           LOCATE THE COLUMN OF LARGEST NORM AND BRING IT
C           INTO THE PIVOT POSITION.
C
            MAXNRM = 0.0D0
            MAXJ = L
            DO 100 J = L, PU
               IF (QRAUX(J) .LE. MAXNRM) GO TO 90
                  MAXNRM = QRAUX(J)
                  MAXJ = J
   90          CONTINUE
  100       CONTINUE
            IF (MAXJ .EQ. L) GO TO 110
               CALL DSWAP(N,X(1,L),1,X(1,MAXJ),1)
               QRAUX(MAXJ) = QRAUX(L)
               WORK(MAXJ) = WORK(L)
               JP = JPVT(MAXJ)
               JPVT(MAXJ) = JPVT(L)
               JPVT(L) = JP
  110       CONTINUE
  120    CONTINUE
         QRAUX(L) = 0.0D0
         IF (L .EQ. N) GO TO 190
C
C           COMPUTE THE HOUSEHOLDER TRANSFORMATION FOR COLUMN L.
C
            NRMXL = DNRM2(N-L+1,X(L,L),1)
            IF (NRMXL .EQ. 0.0D0) GO TO 180
               IF (X(L,L) .NE. 0.0D0) NRMXL = DSIGN(NRMXL,X(L,L))
               CALL DSCAL(N-L+1,1.0D0/NRMXL,X(L,L),1)
               X(L,L) = 1.0D0 + X(L,L)
C
C              APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS,
C              UPDATING THE NORMS.
C
               LP1 = L + 1
               IF (P .LT. LP1) GO TO 170
               DO 160 J = LP1, P
                  T = -DDOT(N-L+1,X(L,L),1,X(L,J),1)/X(L,L)
                  CALL DAXPY(N-L+1,T,X(L,L),1,X(L,J),1)
                  IF (J .LT. PL .OR. J .GT. PU) GO TO 150
                  IF (QRAUX(J) .EQ. 0.0D0) GO TO 150
                     TT = 1.0D0 - (DABS(X(L,J))/QRAUX(J))**2
                     TT = DMAX1(TT,0.0D0)
                     T = TT
                     TT = 1.0D0 + 0.05D0*TT*(QRAUX(J)/WORK(J))**2
                     IF (TT .EQ. 1.0D0) GO TO 130
                        QRAUX(J) = QRAUX(J)*DSQRT(T)
                     GO TO 140
  130                CONTINUE
                        QRAUX(J) = DNRM2(N-L,X(L+1,J),1)
                        WORK(J) = QRAUX(J)
  140                CONTINUE
  150             CONTINUE
  160          CONTINUE
  170          CONTINUE
C
C              SAVE THE TRANSFORMATION.
C
               QRAUX(L) = X(L,L)
               X(L,L) = -NRMXL
  180       CONTINUE
  190    CONTINUE
  200 CONTINUE
      RETURN
      END
