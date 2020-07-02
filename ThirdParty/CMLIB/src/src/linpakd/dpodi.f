      SUBROUTINE DPODI(A,LDA,N,DET,JOB)
C***BEGIN PROLOGUE  DPODI
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D2B1B,D3B1B
C***KEYWORDS  DETERMINANT,DOUBLE PRECISION,FACTOR,INVERSE,
C             LINEAR ALGEBRA,LINPACK,MATRIX,POSITIVE DEFINITE
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Computes the determinant and inverse of a certain double
C            precision SYMMETRIC POSITIVE DEFINITE matrix (see abstract)
C            using the factors computed by DPOCO, DPOFA or DQRDC.
C***DESCRIPTION
C
C     DPODI computes the determinant and inverse of a certain
C     double precision symmetric positive definite matrix (see below)
C     using the factors computed by DPOCO, DPOFA or DQRDC.
C
C     On Entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the output  A  from DPOCO or DPOFA
C                or the output  X  from DQRDC.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        JOB     INTEGER
C                = 11   both determinant and inverse.
C                = 01   inverse only.
C                = 10   determinant only.
C
C     On Return
C
C        A       If DPOCO or DPOFA was used to factor  A , then
C                DPODI produces the upper half of INVERSE(A) .
C                If DQRDC was used to decompose  X , then
C                DPODI produces the upper half of inverse(TRANS(X)*X)
C                where TRANS(X) is the transpose.
C                Elements of  A  below the diagonal are unchanged.
C                If the units digit of JOB is zero,  A  is unchanged.
C
C        DET     DOUBLE PRECISION(2)
C                determinant of  A  or of  TRANS(X)*X  if requested.
C                Otherwise not referenced.
C                Determinant = DET(1) * 10.0**DET(2)
C                with  1.0 .LE. DET(1) .LT. 10.0
C                or  DET(1) .EQ. 0.0 .
C
C     Error Condition
C
C        A division by zero will occur if the input factor contains
C        a zero on the diagonal and the inverse is requested.
C        It will not occur if the subroutines are called correctly
C        and if DPOCO or DPOFA has set INFO .EQ. 0 .
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C
C     Subroutines and Functions
C
C     BLAS DAXPY,DSCAL
C     Fortran MOD
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  DAXPY,DSCAL
C***END PROLOGUE  DPODI
      INTEGER LDA,N,JOB
      DOUBLE PRECISION A(LDA,*)
      DOUBLE PRECISION DET(2)
C
      DOUBLE PRECISION T
      DOUBLE PRECISION S
      INTEGER I,J,JM1,K,KP1
C
C     COMPUTE DETERMINANT
C
C***FIRST EXECUTABLE STATEMENT  DPODI
      IF (JOB/10 .EQ. 0) GO TO 70
         DET(1) = 1.0D0
         DET(2) = 0.0D0
         S = 10.0D0
         DO 50 I = 1, N
            DET(1) = A(I,I)**2*DET(1)
C        ...EXIT
            IF (DET(1) .EQ. 0.0D0) GO TO 60
   10       IF (DET(1) .GE. 1.0D0) GO TO 20
               DET(1) = S*DET(1)
               DET(2) = DET(2) - 1.0D0
            GO TO 10
   20       CONTINUE
   30       IF (DET(1) .LT. S) GO TO 40
               DET(1) = DET(1)/S
               DET(2) = DET(2) + 1.0D0
            GO TO 30
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
C
C     COMPUTE INVERSE(R)
C
      IF (MOD(JOB,10) .EQ. 0) GO TO 140
         DO 100 K = 1, N
            A(K,K) = 1.0D0/A(K,K)
            T = -A(K,K)
            CALL DSCAL(K-1,T,A(1,K),1)
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 90
            DO 80 J = KP1, N
               T = A(K,J)
               A(K,J) = 0.0D0
               CALL DAXPY(K,T,A(1,K),1,A(1,J),1)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
C
C        FORM  INVERSE(R) * TRANS(INVERSE(R))
C
         DO 130 J = 1, N
            JM1 = J - 1
            IF (JM1 .LT. 1) GO TO 120
            DO 110 K = 1, JM1
               T = A(K,J)
               CALL DAXPY(K,T,A(1,J),1,A(1,K),1)
  110       CONTINUE
  120       CONTINUE
            T = A(J,J)
            CALL DSCAL(J,T,A(1,J),1)
  130    CONTINUE
  140 CONTINUE
      RETURN
      END
