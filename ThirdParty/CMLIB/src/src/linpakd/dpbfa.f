      SUBROUTINE DPBFA(ABD,LDA,N,M,INFO)
C***BEGIN PROLOGUE  DPBFA
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D2B2
C***KEYWORDS  BANDED,DOUBLE PRECISION,FACTOR,LINEAR ALGEBRA,LINPACK,
C             MATRIX,POSITIVE DEFINITE
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Factors a d.p. SYMMETRIC POSITIVE DEFINITE matrix stored
C            in band form
C***DESCRIPTION
C
C     DPBFA factors a double precision symmetric positive definite
C     matrix stored in band form.
C
C     DPBFA is usually called by DPBCO, but it can be called
C     directly with a saving in time if  RCOND  is not needed.
C
C     On Entry
C
C        ABD     DOUBLE PRECISION(LDA, N)
C                the matrix to be factored.  The columns of the upper
C                triangle are stored in the columns of ABD and the
C                diagonals of the upper triangle are stored in the
C                rows of ABD .  See the comments below for details.
C
C        LDA     INTEGER
C                the leading dimension of the array  ABD .
C                LDA must be .GE. M + 1 .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        M       INTEGER
C                the number of diagonals above the main diagonal.
C                0 .LE. M .LT. N .
C
C     On Return
C
C        ABD     an upper triangular matrix  R , stored in band
C                form, so that  A = TRANS(R)*R .
C
C        INFO    INTEGER
C                = 0  for normal return.
C                = K  if the leading minor of order  K  is not
C                     positive definite.
C
C     Band Storage
C
C           If  A  is a symmetric positive definite band matrix,
C           the following program segment will set up the input.
C
C                   M = (band width above diagonal)
C                   DO 20 J = 1, N
C                      I1 = MAX0(1, J-M)
C                      DO 10 I = I1, J
C                         K = I-J+M+1
C                         ABD(K,J) = A(I,J)
C                10    CONTINUE
C                20 CONTINUE
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C
C     Subroutines and Functions
C
C     BLAS DDOT
C     Fortran MAX0,DSQRT
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  DDOT
C***END PROLOGUE  DPBFA
      INTEGER LDA,N,M,INFO
      DOUBLE PRECISION ABD(LDA,*)
C
      DOUBLE PRECISION DDOT,T
      DOUBLE PRECISION S
      INTEGER IK,J,JK,K,MU
C     BEGIN BLOCK WITH ...EXITS TO 40
C
C
C***FIRST EXECUTABLE STATEMENT  DPBFA
         DO 30 J = 1, N
            INFO = J
            S = 0.0D0
            IK = M + 1
            JK = MAX0(J-M,1)
            MU = MAX0(M+2-J,1)
            IF (M .LT. MU) GO TO 20
            DO 10 K = MU, M
               T = ABD(K,J) - DDOT(K-MU,ABD(IK,JK),1,ABD(MU,J),1)
               T = T/ABD(M+1,JK)
               ABD(K,J) = T
               S = S + T*T
               IK = IK - 1
               JK = JK + 1
   10       CONTINUE
   20       CONTINUE
            S = ABD(M+1,J) - S
C     ......EXIT
            IF (S .LE. 0.0D0) GO TO 40
            ABD(M+1,J) = DSQRT(S)
   30    CONTINUE
         INFO = 0
   40 CONTINUE
      RETURN
      END
