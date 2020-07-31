      SUBROUTINE STRCO(T,LDT,N,RCOND,Z,JOB)
C***BEGIN PROLOGUE  STRCO
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D2A3
C***KEYWORDS  CONDITION,FACTOR,LINEAR ALGEBRA,LINPACK,MATRIX,TRIANGULAR
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Estimates the condition of a real TRIANGULAR matrix.
C***DESCRIPTION
C
C     STRCO estimates the condition of a real triangular matrix.
C
C     On Entry
C
C        T       REAL(LDT,N)
C                T contains the triangular matrix.  The zero
C                elements of the matrix are not referenced, and
C                the corresponding elements of the array can be
C                used to store other information.
C
C        LDT     INTEGER
C                LDT is the leading dimension of the array T.
C
C        N       INTEGER
C                N is the order of the system.
C
C        JOB     INTEGER
C                = 0         T  is lower triangular.
C                = nonzero   T  is upper triangular.
C
C     On Return
C
C        RCOND   REAL
C                an estimate of the reciprocal condition of  T .
C                For the system  T*X = B , relative perturbations
C                in  T  and  B  of size  EPSILON  may cause
C                relative perturbations in  X  of size  EPSILON/RCOND .
C                If  RCOND  is so small that the logical expression
C                           1.0 + RCOND .EQ. 1.0
C                is true, then  T  may be singular to working
C                precision.  In particular,  RCOND  is zero  if
C                exact singularity is detected or the estimate
C                underflows.
C
C        Z       REAL(N)
C                a work vector whose contents are usually unimportant.
C                If  T  is close to a singular matrix, then  Z  is
C                an approximate null vector in the sense that
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C
C     Subroutines and Functions
C
C     BLAS SAXPY,SSCAL,SASUM
C     Fortran ABS,AMAX1,SIGN
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  SASUM,SAXPY,SSCAL
C***END PROLOGUE  STRCO
      INTEGER LDT,N,JOB
      REAL T(LDT,*),Z(*)
      REAL RCOND
C
      REAL W,WK,WKM,EK
      REAL TNORM,YNORM,S,SM,SASUM
      INTEGER I1,J,J1,J2,K,KK,L
      LOGICAL LOWER
C***FIRST EXECUTABLE STATEMENT  STRCO
      LOWER = JOB .EQ. 0
C
C     COMPUTE 1-NORM OF T
C
      TNORM = 0.0E0
      DO 10 J = 1, N
         L = J
         IF (LOWER) L = N + 1 - J
         I1 = 1
         IF (LOWER) I1 = J
         TNORM = AMAX1(TNORM,SASUM(L,T(I1,J),1))
   10 CONTINUE
C
C     RCOND = 1/(NORM(T)*(ESTIMATE OF NORM(INVERSE(T)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  T*Z = Y  AND  TRANS(T)*Y = E .
C     TRANS(T)  IS THE TRANSPOSE OF T .
C     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
C     GROWTH IN THE ELEMENTS OF Y .
C     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
C
C     SOLVE TRANS(T)*Y = E
C
      EK = 1.0E0
      DO 20 J = 1, N
         Z(J) = 0.0E0
   20 CONTINUE
      DO 100 KK = 1, N
         K = KK
         IF (LOWER) K = N + 1 - KK
         IF (Z(K) .NE. 0.0E0) EK = SIGN(EK,-Z(K))
         IF (ABS(EK-Z(K)) .LE. ABS(T(K,K))) GO TO 30
            S = ABS(T(K,K))/ABS(EK-Z(K))
            CALL SSCAL(N,S,Z,1)
            EK = S*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = ABS(WK)
         SM = ABS(WKM)
         IF (T(K,K) .EQ. 0.0E0) GO TO 40
            WK = WK/T(K,K)
            WKM = WKM/T(K,K)
         GO TO 50
   40    CONTINUE
            WK = 1.0E0
            WKM = 1.0E0
   50    CONTINUE
         IF (KK .EQ. N) GO TO 90
            J1 = K + 1
            IF (LOWER) J1 = 1
            J2 = N
            IF (LOWER) J2 = K - 1
            DO 60 J = J1, J2
               SM = SM + ABS(Z(J)+WKM*T(K,J))
               Z(J) = Z(J) + WK*T(K,J)
               S = S + ABS(Z(J))
   60       CONTINUE
            IF (S .GE. SM) GO TO 80
               W = WKM - WK
               WK = WKM
               DO 70 J = J1, J2
                  Z(J) = Z(J) + W*T(K,J)
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         Z(K) = WK
  100 CONTINUE
      S = 1.0E0/SASUM(N,Z,1)
      CALL SSCAL(N,S,Z,1)
C
      YNORM = 1.0E0
C
C     SOLVE T*Z = Y
C
      DO 130 KK = 1, N
         K = N + 1 - KK
         IF (LOWER) K = KK
         IF (ABS(Z(K)) .LE. ABS(T(K,K))) GO TO 110
            S = ABS(T(K,K))/ABS(Z(K))
            CALL SSCAL(N,S,Z,1)
            YNORM = S*YNORM
  110    CONTINUE
         IF (T(K,K) .NE. 0.0E0) Z(K) = Z(K)/T(K,K)
         IF (T(K,K) .EQ. 0.0E0) Z(K) = 1.0E0
         I1 = 1
         IF (LOWER) I1 = K + 1
         IF (KK .GE. N) GO TO 120
            W = -Z(K)
            CALL SAXPY(N-KK,W,T(I1,K),1,Z(I1),1)
  120    CONTINUE
  130 CONTINUE
C     MAKE ZNORM = 1.0
      S = 1.0E0/SASUM(N,Z,1)
      CALL SSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
      IF (TNORM .NE. 0.0E0) RCOND = YNORM/TNORM
      IF (TNORM .EQ. 0.0E0) RCOND = 0.0E0
      RETURN
      END
