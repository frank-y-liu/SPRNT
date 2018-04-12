      SUBROUTINE DSPFA(AP,N,KPVT,INFO)
C***BEGIN PROLOGUE  DSPFA
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D2B1A
C***KEYWORDS  DOUBLE PRECISION,FACTOR,LINEAR ALGEBRA,LINPACK,MATRIX,
C             PACKED,SYMMETRIC
C***AUTHOR  BUNCH, J., (UCSD)
C***PURPOSE  Factors a double precision SYMMETRIC matrix stored in
C            packed form by elimination with symmetric pivoting.
C***DESCRIPTION
C
C     DSPFA factors a double precision symmetric matrix stored in
C     packed form by elimination with symmetric pivoting.
C
C     To solve  A*X = B , follow DSPFA by DSPSL.
C     To compute  INVERSE(A)*C , follow DSPFA by DSPSL.
C     To compute  DETERMINANT(A) , follow DSPFA by DSPDI.
C     To compute  INERTIA(A) , follow DSPFA by DSPDI.
C     To compute  INVERSE(A) , follow DSPFA by DSPDI.
C
C     On Entry
C
C        AP      DOUBLE PRECISION (N*(N+1)/2)
C                the packed form of a symmetric matrix  A .  The
C                columns of the upper triangle are stored sequentially
C                in a one-dimensional array of length  N*(N+1)/2 .
C                See comments below for details.
C
C        N       INTEGER
C                the order of the matrix  A .
C
C     Output
C
C        AP      a block diagonal matrix and the multipliers which
C                were used to obtain it stored in packed form.
C                The factorization can be written  A = U*D*TRANS(U)
C                where  U  is a product of permutation and unit
C                upper triangular matrices, TRANS(U) is the
C                transpose of  U , and  D  is block diagonal
C                with 1 by 1 and 2 by 2 blocks.
C
C        KPVT    INTEGER(N)
C                an integer vector of pivot indices.
C
C        INFO    INTEGER
C                = 0  normal value.
C                = K  if the K-th pivot block is singular.  This is
C                     not an error condition for this subroutine,
C                     but it does indicate that DSPSL or DSPDI may
C                     divide by zero if called.
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
C                      AP(K)  = A(I,J)
C             10    CONTINUE
C             20 CONTINUE
C
C     LINPACK.  This version dated 08/14/78 .
C     James Bunch, Univ. Calif. San Diego, Argonne Nat. Lab.
C
C     Subroutines and Functions
C
C     BLAS DAXPY,DSWAP,IDAMAX
C     Fortran DABS,DMAX1,DSQRT
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  DAXPY,DSWAP,IDAMAX
C***END PROLOGUE  DSPFA
      INTEGER N,KPVT(*),INFO
      DOUBLE PRECISION AP(*)
C
      DOUBLE PRECISION AK,AKM1,BK,BKM1,DENOM,MULK,MULKM1,T
      DOUBLE PRECISION ABSAKK,ALPHA,COLMAX,ROWMAX
      INTEGER IDAMAX,IJ,IJJ,IK,IKM1,IM,IMAX,IMAXP1,IMIM,IMJ,IMK
      INTEGER J,JJ,JK,JKM1,JMAX,JMIM,K,KK,KM1,KM1K,KM1KM1,KM2,KSTEP
      LOGICAL SWAP
C
C     INITIALIZE
C
C     ALPHA IS USED IN CHOOSING PIVOT BLOCK SIZE.
C***FIRST EXECUTABLE STATEMENT  DSPFA
      ALPHA = (1.0D0 + DSQRT(17.0D0))/8.0D0
C
      INFO = 0
C
C     MAIN LOOP ON K, WHICH GOES FROM N TO 1.
C
      K = N
      IK = (N*(N - 1))/2
   10 CONTINUE
C
C        LEAVE THE LOOP IF K=0 OR K=1.
C
C     ...EXIT
         IF (K .EQ. 0) GO TO 200
         IF (K .GT. 1) GO TO 20
            KPVT(1) = 1
            IF (AP(1) .EQ. 0.0D0) INFO = 1
C     ......EXIT
            GO TO 200
   20    CONTINUE
C
C        THIS SECTION OF CODE DETERMINES THE KIND OF
C        ELIMINATION TO BE PERFORMED.  WHEN IT IS COMPLETED,
C        KSTEP WILL BE SET TO THE SIZE OF THE PIVOT BLOCK, AND
C        SWAP WILL BE SET TO .TRUE. IF AN INTERCHANGE IS
C        REQUIRED.
C
         KM1 = K - 1
         KK = IK + K
         ABSAKK = DABS(AP(KK))
C
C        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
C        COLUMN K.
C
         IMAX = IDAMAX(K-1,AP(IK+1),1)
         IMK = IK + IMAX
         COLMAX = DABS(AP(IMK))
         IF (ABSAKK .LT. ALPHA*COLMAX) GO TO 30
            KSTEP = 1
            SWAP = .FALSE.
         GO TO 90
   30    CONTINUE
C
C           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
C           ROW IMAX.
C
            ROWMAX = 0.0D0
            IMAXP1 = IMAX + 1
            IM = IMAX*(IMAX - 1)/2
            IMJ = IM + 2*IMAX
            DO 40 J = IMAXP1, K
               ROWMAX = DMAX1(ROWMAX,DABS(AP(IMJ)))
               IMJ = IMJ + J
   40       CONTINUE
            IF (IMAX .EQ. 1) GO TO 50
               JMAX = IDAMAX(IMAX-1,AP(IM+1),1)
               JMIM = JMAX + IM
               ROWMAX = DMAX1(ROWMAX,DABS(AP(JMIM)))
   50       CONTINUE
            IMIM = IMAX + IM
            IF (DABS(AP(IMIM)) .LT. ALPHA*ROWMAX) GO TO 60
               KSTEP = 1
               SWAP = .TRUE.
            GO TO 80
   60       CONTINUE
            IF (ABSAKK .LT. ALPHA*COLMAX*(COLMAX/ROWMAX)) GO TO 70
               KSTEP = 1
               SWAP = .FALSE.
            GO TO 80
   70       CONTINUE
               KSTEP = 2
               SWAP = IMAX .NE. KM1
   80       CONTINUE
   90    CONTINUE
         IF (DMAX1(ABSAKK,COLMAX) .NE. 0.0D0) GO TO 100
C
C           COLUMN K IS ZERO.  SET INFO AND ITERATE THE LOOP.
C
            KPVT(K) = K
            INFO = K
         GO TO 190
  100    CONTINUE
         IF (KSTEP .EQ. 2) GO TO 140
C
C           1 X 1 PIVOT BLOCK.
C
            IF (.NOT.SWAP) GO TO 120
C
C              PERFORM AN INTERCHANGE.
C
               CALL DSWAP(IMAX,AP(IM+1),1,AP(IK+1),1)
               IMJ = IK + IMAX
               DO 110 JJ = IMAX, K
                  J = K + IMAX - JJ
                  JK = IK + J
                  T = AP(JK)
                  AP(JK) = AP(IMJ)
                  AP(IMJ) = T
                  IMJ = IMJ - (J - 1)
  110          CONTINUE
  120       CONTINUE
C
C           PERFORM THE ELIMINATION.
C
            IJ = IK - (K - 1)
            DO 130 JJ = 1, KM1
               J = K - JJ
               JK = IK + J
               MULK = -AP(JK)/AP(KK)
               T = MULK
               CALL DAXPY(J,T,AP(IK+1),1,AP(IJ+1),1)
               IJJ = IJ + J
               AP(JK) = MULK
               IJ = IJ - (J - 1)
  130       CONTINUE
C
C           SET THE PIVOT ARRAY.
C
            KPVT(K) = K
            IF (SWAP) KPVT(K) = IMAX
         GO TO 190
  140    CONTINUE
C
C           2 X 2 PIVOT BLOCK.
C
            KM1K = IK + K - 1
            IKM1 = IK - (K - 1)
            IF (.NOT.SWAP) GO TO 160
C
C              PERFORM AN INTERCHANGE.
C
               CALL DSWAP(IMAX,AP(IM+1),1,AP(IKM1+1),1)
               IMJ = IKM1 + IMAX
               DO 150 JJ = IMAX, KM1
                  J = KM1 + IMAX - JJ
                  JKM1 = IKM1 + J
                  T = AP(JKM1)
                  AP(JKM1) = AP(IMJ)
                  AP(IMJ) = T
                  IMJ = IMJ - (J - 1)
  150          CONTINUE
               T = AP(KM1K)
               AP(KM1K) = AP(IMK)
               AP(IMK) = T
  160       CONTINUE
C
C           PERFORM THE ELIMINATION.
C
            KM2 = K - 2
            IF (KM2 .EQ. 0) GO TO 180
               AK = AP(KK)/AP(KM1K)
               KM1KM1 = IKM1 + K - 1
               AKM1 = AP(KM1KM1)/AP(KM1K)
               DENOM = 1.0D0 - AK*AKM1
               IJ = IK - (K - 1) - (K - 2)
               DO 170 JJ = 1, KM2
                  J = KM1 - JJ
                  JK = IK + J
                  BK = AP(JK)/AP(KM1K)
                  JKM1 = IKM1 + J
                  BKM1 = AP(JKM1)/AP(KM1K)
                  MULK = (AKM1*BK - BKM1)/DENOM
                  MULKM1 = (AK*BKM1 - BK)/DENOM
                  T = MULK
                  CALL DAXPY(J,T,AP(IK+1),1,AP(IJ+1),1)
                  T = MULKM1
                  CALL DAXPY(J,T,AP(IKM1+1),1,AP(IJ+1),1)
                  AP(JK) = MULK
                  AP(JKM1) = MULKM1
                  IJJ = IJ + J
                  IJ = IJ - (J - 1)
  170          CONTINUE
  180       CONTINUE
C
C           SET THE PIVOT ARRAY.
C
            KPVT(K) = 1 - K
            IF (SWAP) KPVT(K) = -IMAX
            KPVT(K-1) = KPVT(K)
  190    CONTINUE
         IK = IK - (K - 1)
         IF (KSTEP .EQ. 2) IK = IK - (K - 2)
         K = K - KSTEP
      GO TO 10
  200 CONTINUE
      RETURN
      END
