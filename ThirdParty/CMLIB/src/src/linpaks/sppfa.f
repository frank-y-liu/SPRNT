      SUBROUTINE SPPFA(AP,N,INFO)
C***BEGIN PROLOGUE  SPPFA
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D2B1B
C***KEYWORDS  FACTOR,LINEAR ALGEBRA,LINPACK,MATRIX,PACKED,
C             POSITIVE DEFINITE
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  SPPFA factors a real SYMMETRIC POSITIVE DEFINITE matrix
C            (packed form)
C***DESCRIPTION
C
C     SPPFA factors a real symmetric positive definite matrix
C     stored in packed form.
C
C     SPPFA is usually called by SPPCO, but it can be called
C     directly with a saving in time if  RCOND  is not needed.
C     (Time for SPPCO) = (1 + 18/N)*(Time for SPPFA) .
C
C     On Entry
C
C        AP      REAL (N*(N+1)/2)
C                the packed form of a symmetric matrix  A .  The
C                columns of the upper triangle are stored sequentially
C                in a one-dimensional array of length  N*(N+1)/2 .
C                See comments below for details.
C
C        N       INTEGER
C                the order of the matrix  A .
C
C     On Return
C
C        AP      an upper triangular matrix  R , stored in packed
C                form, so that  A = TRANS(R)*R .
C
C        INFO    INTEGER
C                = 0  for normal return.
C                = K  if the leading minor of order  K  is not
C                     positive definite.
C
C
C     Packed Storage
C
C          The following program segment will pack the upper
C          triangle of a symmetric matrix.
C
C                K = 0
C                DO 20 J = 1, N
C                   DO 10 I = 1, J
C                      K = K + 1
C                      AP(K) = A(I,J)
C             10    CONTINUE
C             20 CONTINUE
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
C***END PROLOGUE  SPPFA
      INTEGER N,INFO
      REAL AP(*)
C
      REAL SDOT,T
      REAL S
      INTEGER J,JJ,JM1,K,KJ,KK
C     BEGIN BLOCK WITH ...EXITS TO 40
C
C
C***FIRST EXECUTABLE STATEMENT  SPPFA
         JJ = 0
         DO 30 J = 1, N
            INFO = J
            S = 0.0E0
            JM1 = J - 1
            KJ = JJ
            KK = 0
            IF (JM1 .LT. 1) GO TO 20
            DO 10 K = 1, JM1
               KJ = KJ + 1
               T = AP(KJ) - SDOT(K-1,AP(KK+1),1,AP(JJ+1),1)
               KK = KK + K
               T = T/AP(KK)
               AP(KJ) = T
               S = S + T*T
   10       CONTINUE
   20       CONTINUE
            JJ = JJ + J
            S = AP(JJ) - S
C     ......EXIT
            IF (S .LE. 0.0E0) GO TO 40
            AP(JJ) = SQRT(S)
   30    CONTINUE
         INFO = 0
   40 CONTINUE
      RETURN
      END
