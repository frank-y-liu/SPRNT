      SUBROUTINE DPOFA(A,LDA,N,INFO)
C***BEGIN PROLOGUE  DPOFA
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D2B1B
C***KEYWORDS  DOUBLE PRECISION,FACTOR,LINEAR ALGEBRA,LINPACK,MATRIX,
C             POSITIVE DEFINITE
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Factors a double precision SYMMETRIC POSITIVE DEFINITE
C            matrix.
C***DESCRIPTION
C
C     DPOFA factors a double precision symmetric positive definite
C     matrix.
C
C     DPOFA is usually called by DPOCO, but it can be called
C     directly with a saving in time if  RCOND  is not needed.
C     (time for DPOCO) = (1 + 18/N)*(time for DPOFA) .
C
C     On Entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the symmetric matrix to be factored.  Only the
C                diagonal and upper triangle are used.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C     On Return
C
C        A       an upper triangular matrix  R  so that  A = TRANS(R)*R
C                where  TRANS(R)  is the transpose.
C                The strict lower triangle is unaltered.
C                If  INFO .NE. 0 , the factorization is not complete.
C
C        INFO    INTEGER
C                = 0  for normal return.
C                = K  signals an error condition.  The leading minor
C                     of order  K  is not positive definite.
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C
C     Subroutines and Functions
C
C     BLAS DDOT
C     Fortran DSQRT
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  DDOT
C***END PROLOGUE  DPOFA
      INTEGER LDA,N,INFO
      DOUBLE PRECISION A(LDA,*)
C
      DOUBLE PRECISION DDOT,T
      DOUBLE PRECISION S
      INTEGER J,JM1,K
C     BEGIN BLOCK WITH ...EXITS TO 40
C
C
C***FIRST EXECUTABLE STATEMENT  DPOFA
         DO 30 J = 1, N
            INFO = J
            S = 0.0D0
            JM1 = J - 1
            IF (JM1 .LT. 1) GO TO 20
            DO 10 K = 1, JM1
               T = A(K,J) - DDOT(K-1,A(1,K),1,A(1,J),1)
               T = T/A(K,K)
               A(K,J) = T
               S = S + T*T
   10       CONTINUE
   20       CONTINUE
            S = A(J,J) - S
C     ......EXIT
            IF (S .LE. 0.0D0) GO TO 40
            A(J,J) = DSQRT(S)
   30    CONTINUE
         INFO = 0
   40 CONTINUE
      RETURN
      END
