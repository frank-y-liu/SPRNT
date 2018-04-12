      SUBROUTINE SPTSL(N,D,E,B)
C***BEGIN PROLOGUE  SPTSL
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D2B2A
C***KEYWORDS  LINEAR ALGEBRA,LINPACK,MATRIX,POSITIVE DEFINITE,SOLVE,
C             TRIDIAGONAL
C***AUTHOR  DONGARRA, J., (ANL)
C***PURPOSE  Solves the system A*X=B where A is POSITIVE DEFINITE
C            and TRIDIAGONAL
C***DESCRIPTION
C
C     SPTSL given a positive definite tridiagonal matrix and a right
C     hand side will find the solution.
C
C     On Entry
C
C        N        INTEGER
C                 is the order of the tridiagonal matrix.
C
C        D        REAL(N)
C                 is the diagonal of the tridiagonal matrix.
C                 On output, D is destroyed.
C
C        E        REAL(N)
C                 is the offdiagonal of the tridiagonal matrix.
C                 E(1) through E(N-1) should contain the
C                 offdiagonal.
C
C        B        REAL(N)
C                 is the right hand side vector.
C
C     On Return
C
C        B        contains the solution.
C
C     LINPACK.  This version dated 08/14/78 .
C     Jack Dongarra, Argonne National Laboratory.
C
C     No externals
C     Fortran MOD
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  SPTSL
      INTEGER N
      REAL D(*),E(*),B(*)
C
      INTEGER K,KBM1,KE,KF,KP1,NM1,NM1D2
      REAL T1,T2
C
C     CHECK FOR 1 X 1 CASE
C
C***FIRST EXECUTABLE STATEMENT  SPTSL
      IF (N .NE. 1) GO TO 10
         B(1) = B(1)/D(1)
      GO TO 70
   10 CONTINUE
         NM1 = N - 1
         NM1D2 = NM1/2
         IF (N .EQ. 2) GO TO 30
            KBM1 = N - 1
C
C           ZERO TOP HALF OF SUBDIAGONAL AND BOTTOM HALF OF
C           SUPERDIAGONAL
C
            DO 20 K = 1, NM1D2
               T1 = E(K)/D(K)
               D(K+1) = D(K+1) - T1*E(K)
               B(K+1) = B(K+1) - T1*B(K)
               T2 = E(KBM1)/D(KBM1+1)
               D(KBM1) = D(KBM1) - T2*E(KBM1)
               B(KBM1) = B(KBM1) - T2*B(KBM1+1)
               KBM1 = KBM1 - 1
   20       CONTINUE
   30    CONTINUE
         KP1 = NM1D2 + 1
C
C        CLEAN UP FOR POSSIBLE 2 X 2 BLOCK AT CENTER
C
         IF (MOD(N,2) .NE. 0) GO TO 40
            T1 = E(KP1)/D(KP1)
            D(KP1+1) = D(KP1+1) - T1*E(KP1)
            B(KP1+1) = B(KP1+1) - T1*B(KP1)
            KP1 = KP1 + 1
   40    CONTINUE
C
C        BACK SOLVE STARTING AT THE CENTER, GOING TOWARDS THE TOP
C        AND BOTTOM
C
         B(KP1) = B(KP1)/D(KP1)
         IF (N .EQ. 2) GO TO 60
            K = KP1 - 1
            KE = KP1 + NM1D2 - 1
            DO 50 KF = KP1, KE
               B(K) = (B(K) - E(K)*B(K+1))/D(K)
               B(KF+1) = (B(KF+1) - E(KF)*B(KF))/D(KF+1)
               K = K - 1
   50       CONTINUE
   60    CONTINUE
         IF (MOD(N,2) .EQ. 0) B(1) = (B(1) - E(1)*B(2))/D(1)
   70 CONTINUE
      RETURN
      END
