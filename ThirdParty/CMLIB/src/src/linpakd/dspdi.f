      SUBROUTINE DSPDI(AP,N,KPVT,DET,INERT,WORK,JOB)
C***BEGIN PROLOGUE  DSPDI
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D2B1A,D3B1A
C***KEYWORDS  DETERMINANT,DOUBLE PRECISION,FACTOR,INVERSE,
C             LINEAR ALGEBRA,LINPACK,MATRIX,PACKED,SYMMETRIC
C***AUTHOR  BUNCH, J., (UCSD)
C***PURPOSE  Computes the determinant, inertia and inverse of a
C            double precision SYMMETRIC matrix using the factors from
C            DSPFA, where the matrix is stored in packed form.
C***DESCRIPTION
C
C     DSPDI computes the determinant, inertia and inverse
C     of a double precision symmetric matrix using the factors from
C     DSPFA, where the matrix is stored in packed form.
C
C     On Entry
C
C        AP      DOUBLE PRECISION (N*(N+1)/2)
C                the output from DSPFA.
C
C        N       INTEGER
C                the order of the matrix A.
C
C        KPVT    INTEGER(N)
C                the pivot vector from DSPFA.
C
C        WORK    DOUBLE PRECISION(N)
C                work vector.  Contents ignored.
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
C        AP     contains the upper triangle of the inverse of
C               the original matrix, stored in packed form.
C               The columns of the upper triangle are stored
C               sequentially in a one-dimensional array.
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
C        A division by zero will occur if the inverse is requested
C        and  DSPCO  has set RCOND .EQ. 0.0
C        or  DSPFA  has set  INFO .NE. 0 .
C
C     LINPACK.  This version dated 08/14/78 .
C     James Bunch, Univ. Calif. San Diego, Argonne Nat. Lab.
C
C     Subroutines and Functions
C
C     BLAS DAXPY,DCOPY,DDOT,DSWAP
C     Fortran DABS,IABS,MOD
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  DAXPY,DCOPY,DDOT,DSWAP
C***END PROLOGUE  DSPDI
      INTEGER N,JOB
      DOUBLE PRECISION AP(*),WORK(*)
      DOUBLE PRECISION DET(2)
      INTEGER KPVT(*),INERT(3)
C
      DOUBLE PRECISION AKKP1,DDOT,TEMP
      DOUBLE PRECISION TEN,D,T,AK,AKP1
      INTEGER IJ,IK,IKP1,IKS,J,JB,JK,JKP1
      INTEGER K,KK,KKP1,KM1,KS,KSJ,KSKP1,KSTEP
      LOGICAL NOINV,NODET,NOERT
C***FIRST EXECUTABLE STATEMENT  DSPDI
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
         IK = 0
         DO 130 K = 1, N
            KK = IK + K
            D = AP(KK)
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
                  IKP1 = IK + K
                  KKP1 = IKP1 + K
                  T = DABS(AP(KKP1))
                  D = (D/T)*AP(KKP1+1) - T
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
            IK = IK + K
  130    CONTINUE
  140 CONTINUE
C
C     COMPUTE INVERSE(A)
C
      IF (NOINV) GO TO 270
         K = 1
         IK = 0
  150    IF (K .GT. N) GO TO 260
            KM1 = K - 1
            KK = IK + K
            IKP1 = IK + K
            KKP1 = IKP1 + K
            IF (KPVT(K) .LT. 0) GO TO 180
C
C              1 BY 1
C
               AP(KK) = 1.0D0/AP(KK)
               IF (KM1 .LT. 1) GO TO 170
                  CALL DCOPY(KM1,AP(IK+1),1,WORK,1)
                  IJ = 0
                  DO 160 J = 1, KM1
                     JK = IK + J
                     AP(JK) = DDOT(J,AP(IJ+1),1,WORK,1)
                     CALL DAXPY(J-1,WORK(J),AP(IJ+1),1,AP(IK+1),1)
                     IJ = IJ + J
  160             CONTINUE
                  AP(KK) = AP(KK) + DDOT(KM1,WORK,1,AP(IK+1),1)
  170          CONTINUE
               KSTEP = 1
            GO TO 220
  180       CONTINUE
C
C              2 BY 2
C
               T = DABS(AP(KKP1))
               AK = AP(KK)/T
               AKP1 = AP(KKP1+1)/T
               AKKP1 = AP(KKP1)/T
               D = T*(AK*AKP1 - 1.0D0)
               AP(KK) = AKP1/D
               AP(KKP1+1) = AK/D
               AP(KKP1) = -AKKP1/D
               IF (KM1 .LT. 1) GO TO 210
                  CALL DCOPY(KM1,AP(IKP1+1),1,WORK,1)
                  IJ = 0
                  DO 190 J = 1, KM1
                     JKP1 = IKP1 + J
                     AP(JKP1) = DDOT(J,AP(IJ+1),1,WORK,1)
                     CALL DAXPY(J-1,WORK(J),AP(IJ+1),1,AP(IKP1+1),1)
                     IJ = IJ + J
  190             CONTINUE
                  AP(KKP1+1) = AP(KKP1+1)
     1                         + DDOT(KM1,WORK,1,AP(IKP1+1),1)
                  AP(KKP1) = AP(KKP1)
     1                       + DDOT(KM1,AP(IK+1),1,AP(IKP1+1),1)
                  CALL DCOPY(KM1,AP(IK+1),1,WORK,1)
                  IJ = 0
                  DO 200 J = 1, KM1
                     JK = IK + J
                     AP(JK) = DDOT(J,AP(IJ+1),1,WORK,1)
                     CALL DAXPY(J-1,WORK(J),AP(IJ+1),1,AP(IK+1),1)
                     IJ = IJ + J
  200             CONTINUE
                  AP(KK) = AP(KK) + DDOT(KM1,WORK,1,AP(IK+1),1)
  210          CONTINUE
               KSTEP = 2
  220       CONTINUE
C
C           SWAP
C
            KS = IABS(KPVT(K))
            IF (KS .EQ. K) GO TO 250
               IKS = (KS*(KS - 1))/2
               CALL DSWAP(KS,AP(IKS+1),1,AP(IK+1),1)
               KSJ = IK + KS
               DO 230 JB = KS, K
                  J = K + KS - JB
                  JK = IK + J
                  TEMP = AP(JK)
                  AP(JK) = AP(KSJ)
                  AP(KSJ) = TEMP
                  KSJ = KSJ - (J - 1)
  230          CONTINUE
               IF (KSTEP .EQ. 1) GO TO 240
                  KSKP1 = IKP1 + KS
                  TEMP = AP(KSKP1)
                  AP(KSKP1) = AP(KKP1)
                  AP(KKP1) = TEMP
  240          CONTINUE
  250       CONTINUE
            IK = IK + K
            IF (KSTEP .EQ. 2) IK = IK + K + 1
            K = K + KSTEP
         GO TO 150
  260    CONTINUE
  270 CONTINUE
      RETURN
      END
