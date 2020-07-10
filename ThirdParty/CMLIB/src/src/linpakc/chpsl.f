      SUBROUTINE CHPSL(AP,N,KPVT,B)
C***BEGIN PROLOGUE  CHPSL
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D2D1A
C***KEYWORDS  COMPLEX,HERMITIAN,LINEAR ALGEBRA,LINPACK,MATRIX,PACKED,
C             SOLVE
C***AUTHOR  BUNCH, J., (UCSD)
C***PURPOSE  Solves the COMPLEX HERMITIAN system  A*X=B  using factors
C            of CHPFA.
C***DESCRIPTION
C
C     CHISL solves the complex Hermitian system
C     A * X = B
C     using the factors computed by CHPFA.
C
C     On Entry
C
C        AP      COMPLEX(N*(N+1)/2)
C                the output from CHPFA.
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        KVPT    INTEGER(N)
C                the pivot vector from CHPFA.
C
C        B       COMPLEX(N)
C                the right hand side vector.
C
C     On Return
C
C        B       the solution vector  X .
C
C     Error Condition
C
C        A division by zero may occur if  CHPCO  has set RCOND .EQ. 0.0
C        or  CHPFA  has set INFO .NE. 0  .
C
C     To compute  INVERSE(A) * C  where  C  is a matrix
C     with  P  columns
C           CALL CHPFA(AP,N,KVPT,INFO)
C           IF (INFO .NE. 0) GO TO ...
C           DO 10 J = 1, P
C              CALL CHPSL(AP,N,KVPT,C(1,J))
C        10 CONTINUE
C
C     LINPACK.  This version dated 08/14/78 .
C     James Bunch, Univ. Calif. San Diego, Argonne Nat. Lab.
C
C     Subroutines and Functions
C
C     BLAS CAXPY,CDOTC
C     Fortran CONJG,IABS
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  CAXPY,CDOTC
C***END PROLOGUE  CHPSL
      INTEGER N,KPVT(*)
      COMPLEX AP(*),B(*)
C
      COMPLEX AK,AKM1,BK,BKM1,CDOTC,DENOM,TEMP
      INTEGER IK,IKM1,IKP1,K,KK,KM1K,KM1KM1,KP
C
C     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND
C     D INVERSE TO B.
C
C***FIRST EXECUTABLE STATEMENT  CHPSL
      K = N
      IK = (N*(N - 1))/2
   10 IF (K .EQ. 0) GO TO 80
         KK = IK + K
         IF (KPVT(K) .LT. 0) GO TO 40
C
C           1 X 1 PIVOT BLOCK.
C
            IF (K .EQ. 1) GO TO 30
               KP = KPVT(K)
               IF (KP .EQ. K) GO TO 20
C
C                 INTERCHANGE.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
   20          CONTINUE
C
C              APPLY THE TRANSFORMATION.
C
               CALL CAXPY(K-1,B(K),AP(IK+1),1,B(1),1)
   30       CONTINUE
C
C           APPLY D INVERSE.
C
            B(K) = B(K)/AP(KK)
            K = K - 1
            IK = IK - K
         GO TO 70
   40    CONTINUE
C
C           2 X 2 PIVOT BLOCK.
C
            IKM1 = IK - (K - 1)
            IF (K .EQ. 2) GO TO 60
               KP = IABS(KPVT(K))
               IF (KP .EQ. K - 1) GO TO 50
C
C                 INTERCHANGE.
C
                  TEMP = B(K-1)
                  B(K-1) = B(KP)
                  B(KP) = TEMP
   50          CONTINUE
C
C              APPLY THE TRANSFORMATION.
C
               CALL CAXPY(K-2,B(K),AP(IK+1),1,B(1),1)
               CALL CAXPY(K-2,B(K-1),AP(IKM1+1),1,B(1),1)
   60       CONTINUE
C
C           APPLY D INVERSE.
C
            KM1K = IK + K - 1
            KK = IK + K
            AK = AP(KK)/CONJG(AP(KM1K))
            KM1KM1 = IKM1 + K - 1
            AKM1 = AP(KM1KM1)/AP(KM1K)
            BK = B(K)/CONJG(AP(KM1K))
            BKM1 = B(K-1)/AP(KM1K)
            DENOM = AK*AKM1 - 1.0E0
            B(K) = (AKM1*BK - BKM1)/DENOM
            B(K-1) = (AK*BKM1 - BK)/DENOM
            K = K - 2
            IK = IK - (K + 1) - K
   70    CONTINUE
      GO TO 10
   80 CONTINUE
C
C     LOOP FORWARD APPLYING THE TRANSFORMATIONS.
C
      K = 1
      IK = 0
   90 IF (K .GT. N) GO TO 160
         IF (KPVT(K) .LT. 0) GO TO 120
C
C           1 X 1 PIVOT BLOCK.
C
            IF (K .EQ. 1) GO TO 110
C
C              APPLY THE TRANSFORMATION.
C
               B(K) = B(K) + CDOTC(K-1,AP(IK+1),1,B(1),1)
               KP = KPVT(K)
               IF (KP .EQ. K) GO TO 100
C
C                 INTERCHANGE.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
  100          CONTINUE
  110       CONTINUE
            IK = IK + K
            K = K + 1
         GO TO 150
  120    CONTINUE
C
C           2 X 2 PIVOT BLOCK.
C
            IF (K .EQ. 1) GO TO 140
C
C              APPLY THE TRANSFORMATION.
C
               B(K) = B(K) + CDOTC(K-1,AP(IK+1),1,B(1),1)
               IKP1 = IK + K
               B(K+1) = B(K+1) + CDOTC(K-1,AP(IKP1+1),1,B(1),1)
               KP = IABS(KPVT(K))
               IF (KP .EQ. K) GO TO 130
C
C                 INTERCHANGE.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
  130          CONTINUE
  140       CONTINUE
            IK = IK + K + K + 1
            K = K + 2
  150    CONTINUE
      GO TO 90
  160 CONTINUE
      RETURN
      END
