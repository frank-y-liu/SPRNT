      SUBROUTINE DPPCO(AP,N,RCOND,Z,INFO)
C***BEGIN PROLOGUE  DPPCO
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D2B1B
C***KEYWORDS  CONDITION,DOUBLE PRECISION,FACTOR,LINEAR ALGEBRA,LINPACK,
C             MATRIX,PACKED,POSITIVE DEFINITE
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Factors a double precision SYMMETRIC POSITIVE DEFINITE
C            matrix stored in packed form and estimates the condition of
C            the matrix.
C***DESCRIPTION
C
C     DPPCO factors a double precision symmetric positive definite
C     matrix stored in packed form
C     and estimates the condition of the matrix.
C
C     If  RCOND  is not needed, DPPFA is slightly faster.
C     To solve  A*X = B , follow DPPCO by DPPSL.
C     To compute  INVERSE(A)*C , follow DPPCO by DPPSL.
C     To compute  DETERMINANT(A) , follow DPPCO by DPPDI.
C     To compute  INVERSE(A) , follow DPPCO by DPPDI.
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
C     On Return
C
C        AP      an upper triangular matrix  R , stored in packed
C                form, so that  A = TRANS(R)*R .
C                If  INFO .NE. 0 , the factorization is not complete.
C
C        RCOND   DOUBLE PRECISION
C                an estimate of the reciprocal condition of  A .
C                For the system  A*X = B , relative perturbations
C                in  A  and  B  of size  EPSILON  may cause
C                relative perturbations in  X  of size  EPSILON/RCOND .
C                If  RCOND  is so small that the logical expression
C                           1.0 + RCOND .EQ. 1.0
C                is true, then  A  may be singular to working
C                precision.  In particular,  RCOND  is zero  if
C                exact singularity is detected or the estimate
C                underflows.  If INFO .NE. 0 , RCOND is unchanged.
C
C        Z       DOUBLE PRECISION(N)
C                a work vector whose contents are usually unimportant.
C                If  A  is singular to working precision, then  Z  is
C                an approximate null vector in the sense that
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C                If  INFO .NE. 0 , Z  is unchanged.
C
C        INFO    INTEGER
C                = 0  for normal return.
C                = K  signals an error condition.  The leading minor
C                     of order  K  is not positive definite.
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
C     LINPACK DPPFA
C     BLAS DAXPY,DDOT,DSCAL,DASUM
C     Fortran DABS,DMAX1,DREAL,DSIGN
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  DASUM,DAXPY,DDOT,DPPFA,DSCAL
C***END PROLOGUE  DPPCO
      INTEGER N,INFO
      DOUBLE PRECISION AP(*),Z(*)
      DOUBLE PRECISION RCOND
C
      DOUBLE PRECISION DDOT,EK,T,WK,WKM
      DOUBLE PRECISION ANORM,S,DASUM,SM,YNORM
      INTEGER I,IJ,J,JM1,J1,K,KB,KJ,KK,KP1
C
C     FIND NORM OF A
C
C***FIRST EXECUTABLE STATEMENT  DPPCO
      J1 = 1
      DO 30 J = 1, N
         Z(J) = DASUM(J,AP(J1),1)
         IJ = J1
         J1 = J1 + J
         JM1 = J - 1
         IF (JM1 .LT. 1) GO TO 20
         DO 10 I = 1, JM1
            Z(I) = Z(I) + DABS(AP(IJ))
            IJ = IJ + 1
   10    CONTINUE
   20    CONTINUE
   30 CONTINUE
      ANORM = 0.0D0
      DO 40 J = 1, N
         ANORM = DMAX1(ANORM,Z(J))
   40 CONTINUE
C
C     FACTOR
C
      CALL DPPFA(AP,N,INFO)
      IF (INFO .NE. 0) GO TO 180
C
C        RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C        ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
C        THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
C        GROWTH IN THE ELEMENTS OF W  WHERE  TRANS(R)*W = E .
C        THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
C
C        SOLVE TRANS(R)*W = E
C
         EK = 1.0D0
         DO 50 J = 1, N
            Z(J) = 0.0D0
   50    CONTINUE
         KK = 0
         DO 110 K = 1, N
            KK = KK + K
            IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))
            IF (DABS(EK-Z(K)) .LE. AP(KK)) GO TO 60
               S = AP(KK)/DABS(EK-Z(K))
               CALL DSCAL(N,S,Z,1)
               EK = S*EK
   60       CONTINUE
            WK = EK - Z(K)
            WKM = -EK - Z(K)
            S = DABS(WK)
            SM = DABS(WKM)
            WK = WK/AP(KK)
            WKM = WKM/AP(KK)
            KP1 = K + 1
            KJ = KK + K
            IF (KP1 .GT. N) GO TO 100
               DO 70 J = KP1, N
                  SM = SM + DABS(Z(J)+WKM*AP(KJ))
                  Z(J) = Z(J) + WK*AP(KJ)
                  S = S + DABS(Z(J))
                  KJ = KJ + J
   70          CONTINUE
               IF (S .GE. SM) GO TO 90
                  T = WKM - WK
                  WK = WKM
                  KJ = KK + K
                  DO 80 J = KP1, N
                     Z(J) = Z(J) + T*AP(KJ)
                     KJ = KJ + J
   80             CONTINUE
   90          CONTINUE
  100       CONTINUE
            Z(K) = WK
  110    CONTINUE
         S = 1.0D0/DASUM(N,Z,1)
         CALL DSCAL(N,S,Z,1)
C
C        SOLVE R*Y = W
C
         DO 130 KB = 1, N
            K = N + 1 - KB
            IF (DABS(Z(K)) .LE. AP(KK)) GO TO 120
               S = AP(KK)/DABS(Z(K))
               CALL DSCAL(N,S,Z,1)
  120       CONTINUE
            Z(K) = Z(K)/AP(KK)
            KK = KK - K
            T = -Z(K)
            CALL DAXPY(K-1,T,AP(KK+1),1,Z(1),1)
  130    CONTINUE
         S = 1.0D0/DASUM(N,Z,1)
         CALL DSCAL(N,S,Z,1)
C
         YNORM = 1.0D0
C
C        SOLVE TRANS(R)*V = Y
C
         DO 150 K = 1, N
            Z(K) = Z(K) - DDOT(K-1,AP(KK+1),1,Z(1),1)
            KK = KK + K
            IF (DABS(Z(K)) .LE. AP(KK)) GO TO 140
               S = AP(KK)/DABS(Z(K))
               CALL DSCAL(N,S,Z,1)
               YNORM = S*YNORM
  140       CONTINUE
            Z(K) = Z(K)/AP(KK)
  150    CONTINUE
         S = 1.0D0/DASUM(N,Z,1)
         CALL DSCAL(N,S,Z,1)
         YNORM = S*YNORM
C
C        SOLVE R*Z = V
C
         DO 170 KB = 1, N
            K = N + 1 - KB
            IF (DABS(Z(K)) .LE. AP(KK)) GO TO 160
               S = AP(KK)/DABS(Z(K))
               CALL DSCAL(N,S,Z,1)
               YNORM = S*YNORM
  160       CONTINUE
            Z(K) = Z(K)/AP(KK)
            KK = KK - K
            T = -Z(K)
            CALL DAXPY(K-1,T,AP(KK+1),1,Z(1),1)
  170    CONTINUE
C        MAKE ZNORM = 1.0
         S = 1.0D0/DASUM(N,Z,1)
         CALL DSCAL(N,S,Z,1)
         YNORM = S*YNORM
C
         IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM
         IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0
  180 CONTINUE
      RETURN
      END
