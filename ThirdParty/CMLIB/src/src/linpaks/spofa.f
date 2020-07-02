      SUBROUTINE SPOFA(A,LDA,N,INFO)
C***BEGIN PROLOGUE  SPOFA
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D2B1B
C***KEYWORDS  FACTOR,LINEAR ALGEBRA,LINPACK,MATRIX,POSITIVE DEFINITE
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Factors a real SYMMETRIC POSITIVE DEFINITE matrix.
C***DESCRIPTION
C
C     SPOFA factors a real symmetric positive definite matrix.
C
C     SPOFA is usually called by SPOCO, but it can be called
C     directly with a saving in time if  RCOND  is not needed.
C     (Time for SPOCO) = (1 + 18/N)*(Time for SPOFA) .
C
C     On Entry
C
C        A       REAL(LDA, N)
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
C     BLAS SDOT
C     Fortran SQRT
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  SDOT
C***END PROLOGUE  SPOFA
      INTEGER LDA,N,INFO
      REAL A(LDA,*)
C
      REAL SDOT,T
      REAL S
      INTEGER J,JM1,K
C     BEGIN BLOCK WITH ...EXITS TO 40
C
C
C***FIRST EXECUTABLE STATEMENT  SPOFA
         DO 30 J = 1, N
            INFO = J
            S = 0.0E0
            JM1 = J - 1
            IF (JM1 .LT. 1) GO TO 20
            DO 10 K = 1, JM1
               T = A(K,J) - SDOT(K-1,A(1,K),1,A(1,J),1)
               T = T/A(K,K)
               A(K,J) = T
               S = S + T*T
   10       CONTINUE
   20       CONTINUE
            S = A(J,J) - S
C     ......EXIT
            IF (S .LE. 0.0E0) GO TO 40
            A(J,J) = SQRT(S)
   30    CONTINUE
         INFO = 0
   40 CONTINUE
      RETURN
      END
