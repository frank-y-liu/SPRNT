      SUBROUTINE CPPFA(AP,N,INFO)
C***BEGIN PROLOGUE  CPPFA
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D2D1B
C***KEYWORDS  COMPLEX,FACTOR,LINEAR ALGEBRA,LINPACK,MATRIX,PACKED,
C             POSITIVE DEFINITE
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Factors a COMPLEX HERMITIAN POSITIVE DEFINITE matrix
C            (packed form)
C***DESCRIPTION
C
C     CPPFA factors a complex Hermitian positive definite matrix
C     stored in packed form.
C
C     CPPFA is usually called by CPPCO, but it can be called
C     directly with a saving in time if  RCOND  is not needed.
C     (Time for CPPCO) = (1 + 18/N)*(Time for CPPFA) .
C
C     On Entry
C
C        AP      COMPLEX (N*(N+1)/2)
C                the packed form of a Hermitian matrix  A .  The
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
C                form, so that  A = CTRANS(R)*R .
C
C        INFO    INTEGER
C                = 0  for normal return.
C                = K  If the leading minor of order  K  is not
C                     positive definite.
C
C
C     Packed Storage
C
C          The following program segment will pack the upper
C          triangle of a Hermitian matrix.
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
C     BLAS CDOTC
C     Fortran AIMAG,CMPLX,CONJG,REAL,SQRT
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  CDOTC
C***END PROLOGUE  CPPFA
      INTEGER N,INFO
      COMPLEX AP(*)
C
      COMPLEX CDOTC,T
      REAL S
      INTEGER J,JJ,JM1,K,KJ,KK
C     BEGIN BLOCK WITH ...EXITS TO 40
C
C
C***FIRST EXECUTABLE STATEMENT  CPPFA
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
               T = AP(KJ) - CDOTC(K-1,AP(KK+1),1,AP(JJ+1),1)
               KK = KK + K
               T = T/AP(KK)
               AP(KJ) = T
               S = S + REAL(T*CONJG(T))
   10       CONTINUE
   20       CONTINUE
            JJ = JJ + J
            S = REAL(AP(JJ)) - S
C     ......EXIT
            IF (S .LE. 0.0E0 .OR. AIMAG(AP(JJ)) .NE. 0.0E0) GO TO 40
            AP(JJ) = CMPLX(SQRT(S),0.0E0)
   30    CONTINUE
         INFO = 0
   40 CONTINUE
      RETURN
      END
