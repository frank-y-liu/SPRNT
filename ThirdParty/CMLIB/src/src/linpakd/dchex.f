      SUBROUTINE DCHEX(R,LDR,P,K,L,Z,LDZ,NZ,C,S,JOB)
C***BEGIN PROLOGUE  DCHEX
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D7B
C***KEYWORDS  CHOLESKY DECOMPOSITION,DOUBLE PRECISION,EXCHANGE,
C             LINEAR ALGEBRA,LINPACK,MATRIX,POSITIVE DEFINITE
C***AUTHOR  STEWART, G. W., (U. OF MARYLAND)
C***PURPOSE  Updates the Cholesky factorization  A=TRANS(R)*R  of a
C            POSITIVE DEFINITE matrix A of order P under diagonal
C            permutations of the form  TRANS(E)*A*E  where E is a
C            permutation matrix.
C***DESCRIPTION
C
C     DCHEX updates the Cholesky factorization
C
C                   A = TRANS(R)*R
C
C     of a positive definite matrix A of order P under diagonal
C     permutations of the form
C
C                   TRANS(E)*A*E
C
C     where E is a permutation matrix.  Specifically, given
C     an upper triangular matrix R and a permutation matrix
C     E (which is specified by K, L, and JOB), DCHEX determines
C     an orthogonal matrix U such that
C
C                           U*R*E = RR,
C
C     where RR is upper triangular.  At the users option, the
C     transformation U will be multiplied into the array Z.
C     If A = TRANS(X)*X, so that R is the triangular part of the
C     QR factorization of X, then RR is the triangular part of the
C     QR factorization of X*E, i.e. X with its columns permuted.
C     For a less terse description of what DCHEX does and how
C     it may be applied, see the LINPACK guide.
C
C     The matrix Q is determined as the product U(L-K)*...*U(1)
C     of plane rotations of the form
C
C                           (    C(I)       S(I) )
C                           (                    ) ,
C                           (    -S(I)      C(I) )
C
C     where C(I) is double precision.  The rows these rotations operate
C     on are described below.
C
C     There are two types of permutations, which are determined
C     by the value of JOB.
C
C     1. Right circular shift (JOB = 1).
C
C         The columns are rearranged in the following order.
C
C                1,...,K-1,L,K,K+1,...,L-1,L+1,...,P.
C
C         U is the product of L-K rotations U(I), where U(I)
C         acts in the (L-I,L-I+1)-plane.
C
C     2. Left circular shift (JOB = 2).
C         The columns are rearranged in the following order
C
C                1,...,K-1,K+1,K+2,...,L,K,L+1,...,P.
C
C         U is the product of L-K rotations U(I), where U(I)
C         acts in the (K+I-1,K+I)-plane.
C
C     On Entry
C
C         R      DOUBLE PRECISION(LDR,P), where LDR .GE. P.
C                R contains the upper triangular factor
C                that is to be updated.  Elements of R
C                below the diagonal are not referenced.
C
C         LDR    INTEGER.
C                LDR is the leading dimension of the array R.
C
C         P      INTEGER.
C                P is the order of the matrix R.
C
C         K      INTEGER.
C                K is the first column to be permuted.
C
C         L      INTEGER.
C                L is the last column to be permuted.
C                L must be strictly greater than K.
C
C         Z      DOUBLE PRECISION(LDZ,N)Z), where LDZ .GE. P.
C                Z is an array of NZ P-vectors into which the
C                transformation U is multiplied.  Z is
C                not referenced if NZ = 0.
C
C         LDZ    INTEGER.
C                LDZ is the leading dimension of the array Z.
C
C         NZ     INTEGER.
C                NZ is the number of columns of the matrix Z.
C
C         JOB    INTEGER.
C                JOB determines the type of permutation.
C                       JOB = 1  right circular shift.
C                       JOB = 2  left circular shift.
C
C     On Return
C
C         R      contains the updated factor.
C
C         Z      contains the updated matrix Z.
C
C         C      DOUBLE PRECISION(P).
C                C contains the cosines of the transforming rotations.
C
C         S      DOUBLE PRECISION(P).
C                S contains the sines of the transforming rotations.
C
C     LINPACK.  This version dated 08/14/78 .
C     G. W. Stewart, University of Maryland, Argonne National Lab.
C
C     DCHEX uses the following functions and subroutines.
C
C     BLAS DROTG
C     Fortran MIN0
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  DROTG
C***END PROLOGUE  DCHEX
      INTEGER LDR,P,K,L,LDZ,NZ,JOB
      DOUBLE PRECISION R(LDR,*),Z(LDZ,*),S(*)
      DOUBLE PRECISION C(*)
C
      INTEGER I,II,IL,IU,J,JJ,KM1,KP1,LMK,LM1
      DOUBLE PRECISION RJP1J,T
C
C     INITIALIZE
C
C***FIRST EXECUTABLE STATEMENT  DCHEX
      KM1 = K - 1
      KP1 = K + 1
      LMK = L - K
      LM1 = L - 1
C
C     PERFORM THE APPROPRIATE TASK.
C
      GO TO (10,130), JOB
C
C     RIGHT CIRCULAR SHIFT.
C
   10 CONTINUE
C
C        REORDER THE COLUMNS.
C
         DO 20 I = 1, L
            II = L - I + 1
            S(I) = R(II,L)
   20    CONTINUE
         DO 40 JJ = K, LM1
            J = LM1 - JJ + K
            DO 30 I = 1, J
               R(I,J+1) = R(I,J)
   30       CONTINUE
            R(J+1,J+1) = 0.0D0
   40    CONTINUE
         IF (K .EQ. 1) GO TO 60
            DO 50 I = 1, KM1
               II = L - I + 1
               R(I,K) = S(II)
   50       CONTINUE
   60    CONTINUE
C
C        CALCULATE THE ROTATIONS.
C
         T = S(1)
         DO 70 I = 1, LMK
            CALL DROTG(S(I+1),T,C(I),S(I))
            T = S(I+1)
   70    CONTINUE
         R(K,K) = T
         DO 90 J = KP1, P
            IL = MAX0(1,L-J+1)
            DO 80 II = IL, LMK
               I = L - II
               T = C(II)*R(I,J) + S(II)*R(I+1,J)
               R(I+1,J) = C(II)*R(I+1,J) - S(II)*R(I,J)
               R(I,J) = T
   80       CONTINUE
   90    CONTINUE
C
C        IF REQUIRED, APPLY THE TRANSFORMATIONS TO Z.
C
         IF (NZ .LT. 1) GO TO 120
         DO 110 J = 1, NZ
            DO 100 II = 1, LMK
               I = L - II
               T = C(II)*Z(I,J) + S(II)*Z(I+1,J)
               Z(I+1,J) = C(II)*Z(I+1,J) - S(II)*Z(I,J)
               Z(I,J) = T
  100       CONTINUE
  110    CONTINUE
  120    CONTINUE
      GO TO 260
C
C     LEFT CIRCULAR SHIFT
C
  130 CONTINUE
C
C        REORDER THE COLUMNS
C
         DO 140 I = 1, K
            II = LMK + I
            S(II) = R(I,K)
  140    CONTINUE
         DO 160 J = K, LM1
            DO 150 I = 1, J
               R(I,J) = R(I,J+1)
  150       CONTINUE
            JJ = J - KM1
            S(JJ) = R(J+1,J+1)
  160    CONTINUE
         DO 170 I = 1, K
            II = LMK + I
            R(I,L) = S(II)
  170    CONTINUE
         DO 180 I = KP1, L
            R(I,L) = 0.0D0
  180    CONTINUE
C
C        REDUCTION LOOP.
C
         DO 220 J = K, P
            IF (J .EQ. K) GO TO 200
C
C              APPLY THE ROTATIONS.
C
               IU = MIN0(J-1,L-1)
               DO 190 I = K, IU
                  II = I - K + 1
                  T = C(II)*R(I,J) + S(II)*R(I+1,J)
                  R(I+1,J) = C(II)*R(I+1,J) - S(II)*R(I,J)
                  R(I,J) = T
  190          CONTINUE
  200       CONTINUE
            IF (J .GE. L) GO TO 210
               JJ = J - K + 1
               T = S(JJ)
               CALL DROTG(R(J,J),T,C(JJ),S(JJ))
  210       CONTINUE
  220    CONTINUE
C
C        APPLY THE ROTATIONS TO Z.
C
         IF (NZ .LT. 1) GO TO 250
         DO 240 J = 1, NZ
            DO 230 I = K, LM1
               II = I - KM1
               T = C(II)*Z(I,J) + S(II)*Z(I+1,J)
               Z(I+1,J) = C(II)*Z(I+1,J) - S(II)*Z(I,J)
               Z(I,J) = T
  230       CONTINUE
  240    CONTINUE
  250    CONTINUE
  260 CONTINUE
      RETURN
      END
