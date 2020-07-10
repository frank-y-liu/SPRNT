      SUBROUTINE DSIDI(A,LDA,N,KPVT,DET,INERT,WORK,JOB)
C***BEGIN PROLOGUE  DSIDI
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D2B1A,D3B1A
C***KEYWORDS  DETERMINANT,DOUBLE PRECISION,FACTOR,INVERSE,
C             LINEAR ALGEBRA,LINPACK,MATRIX,SYMMETRIC
C***AUTHOR  BUNCH, J., (UCSD)
C***PURPOSE  Computes the determinant, inertia and inverse of a d.p.
C            SYMMETRIC matrix using the factors from DSIFA.
C***DESCRIPTION
C
C     DSIDI computes the determinant, inertia and inverse
C     of a double precision symmetric matrix using the factors from
C     DSIFA.
C
C     On Entry
C
C        A       DOUBLE PRECISION(LDA,N)
C                the output from DSIFA.
C
C        LDA     INTEGER
C                the leading dimension of the array A.
C
C        N       INTEGER
C                the order of the matrix A.
C
C        KPVT    INTEGER(N)
C                the pivot vector from DSIFA.
C
C        WORK    DOUBLE PRECISION(N)
C                work vector.  Contents destroyed.
C
C        JOB     INTEGER
C                JOB has the decimal expansion  ABC  where
C                   if  C .NE. 0, the inverse is computed,
C                   if  B .NE. 0, the determinant is computed,
C                   if  A .NE. 0, the inertia is computed.
C
C                For example, JOB = 111  gives all three.
C
C     On Return
C
C        Variables not requested by JOB are not used.
C
C        A      contains the upper triangle of the inverse of
C               the original matrix.  The strict lower triangle
C               is never referenced.
C
C        DET    DOUBLE PRECISION(2)
C               determinant of original matrix.
C               DETERMINANT = DET(1) * 10.0**DET(2)
C               with 1.0 .LE. DABS(DET(1)) .LT. 10.0
C               or DET(1) = 0.0.
C
C        INERT  INTEGER(3)
C               the inertia of the original matrix.
C               INERT(1)  =  number of positive eigenvalues.
C               INERT(2)  =  number of negative eigenvalues.
C               INERT(3)  =  number of zero eigenvalues.
C
C     Error Condition
C
C        A division by zero may occur if the inverse is requested
C        and  DSICO  has set RCOND .EQ. 0.0
C        or  DSIFA  has set  INFO .NE. 0 .
C
C     LINPACK.  This version dated 08/14/78 .
C     James Bunch, Univ. Calif. San Diego, Argonne Nat. Lab
C
C     Subroutines and Functions
C
C     BLAS DAXPY,DCOPY,DDOT,DSWAP
C     Fortran DABS,IABS,MOD
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  DAXPY,DCOPY,DDOT,DSWAP
C***END PROLOGUE  DSIDI
      INTEGER LDA,N,JOB
      DOUBLE PRECISION A(LDA,*),WORK(*)
      DOUBLE PRECISION DET(2)
      INTEGER KPVT(*),INERT(3)
C
      DOUBLE PRECISION AKKP1,DDOT,TEMP
      DOUBLE PRECISION TEN,D,T,AK,AKP1
      INTEGER J,JB,K,KM1,KS,KSTEP
      LOGICAL NOINV,NODET,NOERT
C***FIRST EXECUTABLE STATEMENT  DSIDI
      NOINV = MOD(JOB,10) .EQ. 0
      NODET = MOD(JOB,100)/10 .EQ. 0
      NOERT = MOD(JOB,1000)/100 .EQ. 0
C
      IF (NODET .AND. NOERT) GO TO 140
         IF (NOERT) GO TO 10
            INERT(1) = 0
            INERT(2) = 0
            INERT(3) = 0
   10    CONTINUE
         IF (NODET) GO TO 20
            DET(1) = 1.0D0
            DET(2) = 0.0D0
            TEN = 10.0D0
   20    CONTINUE
         T = 0.0D0
         DO 130 K = 1, N
            D = A(K,K)
C
C           CHECK IF 1 BY 1
C
            IF (KPVT(K) .GT. 0) GO TO 50
C
C              2 BY 2 BLOCK
C              USE DET (D  S)  =  (D/T * C - T) * T  ,  T = DABS(S)
C                      (S  C)
C              TO AVOID UNDERFLOW/OVERFLOW TROUBLES.
C              TAKE TWO PASSES THROUGH SCALING.  USE  T  FOR FLAG.
C
               IF (T .NE. 0.0D0) GO TO 30
                  T = DABS(A(K,K+1))
                  D = (D/T)*A(K+1,K+1) - T
               GO TO 40
   30          CONTINUE
                  D = T
                  T = 0.0D0
   40          CONTINUE
   50       CONTINUE
C
            IF (NOERT) GO TO 60
               IF (D .GT. 0.0D0) INERT(1) = INERT(1) + 1
               IF (D .LT. 0.0D0) INERT(2) = INERT(2) + 1
               IF (D .EQ. 0.0D0) INERT(3) = INERT(3) + 1
   60       CONTINUE
C
            IF (NODET) GO TO 120
               DET(1) = D*DET(1)
               IF (DET(1) .EQ. 0.0D0) GO TO 110
   70             IF (DABS(DET(1)) .GE. 1.0D0) GO TO 80
                     DET(1) = TEN*DET(1)
                     DET(2) = DET(2) - 1.0D0
                  GO TO 70
   80             CONTINUE
   90             IF (DABS(DET(1)) .LT. TEN) GO TO 100
                     DET(1) = DET(1)/TEN
                     DET(2) = DET(2) + 1.0D0
                  GO TO 90
  100             CONTINUE
  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
  140 CONTINUE
C
C     COMPUTE INVERSE(A)
C
      IF (NOINV) GO TO 270
         K = 1
  150    IF (K .GT. N) GO TO 260
            KM1 = K - 1
            IF (KPVT(K) .LT. 0) GO TO 180
C
C              1 BY 1
C
               A(K,K) = 1.0D0/A(K,K)
               IF (KM1 .LT. 1) GO TO 170
                  CALL DCOPY(KM1,A(1,K),1,WORK,1)
                  DO 160 J = 1, KM1
                     A(J,K) = DDOT(J,A(1,J),1,WORK,1)
                     CALL DAXPY(J-1,WORK(J),A(1,J),1,A(1,K),1)
  160             CONTINUE
                  A(K,K) = A(K,K) + DDOT(KM1,WORK,1,A(1,K),1)
  170          CONTINUE
               KSTEP = 1
            GO TO 220
  180       CONTINUE
C
C              2 BY 2
C
               T = DABS(A(K,K+1))
               AK = A(K,K)/T
               AKP1 = A(K+1,K+1)/T
               AKKP1 = A(K,K+1)/T
               D = T*(AK*AKP1 - 1.0D0)
               A(K,K) = AKP1/D
               A(K+1,K+1) = AK/D
               A(K,K+1) = -AKKP1/D
               IF (KM1 .LT. 1) GO TO 210
                  CALL DCOPY(KM1,A(1,K+1),1,WORK,1)
                  DO 190 J = 1, KM1
                     A(J,K+1) = DDOT(J,A(1,J),1,WORK,1)
                     CALL DAXPY(J-1,WORK(J),A(1,J),1,A(1,K+1),1)
  190             CONTINUE
                  A(K+1,K+1) = A(K+1,K+1) + DDOT(KM1,WORK,1,A(1,K+1),1)
                  A(K,K+1) = A(K,K+1) + DDOT(KM1,A(1,K),1,A(1,K+1),1)
                  CALL DCOPY(KM1,A(1,K),1,WORK,1)
                  DO 200 J = 1, KM1
                     A(J,K) = DDOT(J,A(1,J),1,WORK,1)
                     CALL DAXPY(J-1,WORK(J),A(1,J),1,A(1,K),1)
  200             CONTINUE
                  A(K,K) = A(K,K) + DDOT(KM1,WORK,1,A(1,K),1)
  210          CONTINUE
               KSTEP = 2
  220       CONTINUE
C
C           SWAP
C
            KS = IABS(KPVT(K))
            IF (KS .EQ. K) GO TO 250
               CALL DSWAP(KS,A(1,KS),1,A(1,K),1)
               DO 230 JB = KS, K
                  J = K + KS - JB
                  TEMP = A(J,K)
                  A(J,K) = A(KS,J)
                  A(KS,J) = TEMP
  230          CONTINUE
               IF (KSTEP .EQ. 1) GO TO 240
                  TEMP = A(KS,K+1)
                  A(KS,K+1) = A(K,K+1)
                  A(K,K+1) = TEMP
  240          CONTINUE
  250       CONTINUE
            K = K + KSTEP
         GO TO 150
  260    CONTINUE
  270 CONTINUE
      RETURN
      END
