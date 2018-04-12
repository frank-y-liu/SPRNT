      SUBROUTINE CPOFA(A,LDA,N,INFO)
C***BEGIN PROLOGUE  CPOFA
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D2D1B
C***KEYWORDS  COMPLEX,FACTOR,LINEAR ALGEBRA,LINPACK,MATRIX,
C             POSITIVE DEFINITE
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Factors a COMPLEX HERMITIAN POSITIVE DEFINITE matrix.
C***DESCRIPTION
C
C     CPOFA factors a complex Hermitian positive definite matrix.
C
C     CPOFA is usually called by CPOCO, but it can be called
C     directly with a saving in time if  RCOND  is not needed.
C     (Time for CPOCO) = (1 + 18/N)*(Time for CPOFA) .
C
C     On Entry
C
C        A       COMPLEX(LDA, N)
C                the Hermitian matrix to be factored.  Only the
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
C        A       an upper triangular matrix  R  so that  A =
C                CTRANS(R)*R where  CTRANS(R)  is the conjugate
C                transpose.  The strict lower triangle is unaltered.
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
C     BLAS CDOTC
C     Fortran AIMAG,CMPLX,CONJG,REAL,SQRT
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  CDOTC
C***END PROLOGUE  CPOFA
      INTEGER LDA,N,INFO
      COMPLEX A(LDA,*)
C
      COMPLEX CDOTC,T
      REAL S
      INTEGER J,JM1,K
C     BEGIN BLOCK WITH ...EXITS TO 40
C
C
C***FIRST EXECUTABLE STATEMENT  CPOFA
         DO 30 J = 1, N
            INFO = J
            S = 0.0E0
            JM1 = J - 1
            IF (JM1 .LT. 1) GO TO 20
            DO 10 K = 1, JM1
               T = A(K,J) - CDOTC(K-1,A(1,K),1,A(1,J),1)
               T = T/A(K,K)
               A(K,J) = T
               S = S + REAL(T*CONJG(T))
   10       CONTINUE
   20       CONTINUE
            S = REAL(A(J,J)) - S
C     ......EXIT
            IF (S .LE. 0.0E0 .OR. AIMAG(A(J,J)) .NE. 0.0E0) GO TO 40
            A(J,J) = CMPLX(SQRT(S),0.0E0)
   30    CONTINUE
         INFO = 0
   40 CONTINUE
      RETURN
      END
