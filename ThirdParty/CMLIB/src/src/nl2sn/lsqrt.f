      SUBROUTINE LSQRT(N1, N, L, A, IRC)
C
C  ***  COMPUTE ROWS N1 THROUGH N OF THE CHOLESKY FACTOR  L  OF
C  ***  A = L*(L**T),  WHERE  L  AND THE LOWER TRIANGLE OF  A  ARE BOTH
C  ***  STORED COMPACTLY BY ROWS (AND MAY OCCUPY THE SAME STORAGE).
C  ***  IRC = 0 MEANS ALL WENT WELL.  IRC = J MEANS THE LEADING
C  ***  PRINCIPAL  J X J  SUBMATRIX OF  A  IS NOT POSITIVE DEFINITE --
C  ***  AND  L(J*(J+1)/2)  CONTAINS THE (NONPOS.) REDUCED J-TH DIAGONAL.
C
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C  ***  PARAMETERS  ***
C
      INTEGER N1, N, IRC
      REAL L(*), A(*)
C     DIMENSION L(N*(N+1)/2), A(N*(N+1)/2)
C
C  ***  LOCAL VARIABLES  ***
C
      INTEGER I, IJ, IK, IM1, I0, J, JK, JM1, J0, K
      REAL T, TD, ZERO
C
C  ***  INTRINSIC FUNCTIONS  ***
C/+
      REAL  SQRT
C/
C/6
      DATA ZERO/0.E+0/
C/7
C     PARAMETER (ZERO=0.E+0)
C/
C
C  ***  BODY  ***
C
      I0 = N1 * (N1 - 1) / 2
      DO 50 I = N1, N
         TD = ZERO
         IF (I .EQ. 1) GO TO 40
         J0 = 0
         IM1 = I - 1
         DO 30 J = 1, IM1
              T = ZERO
              IF (J .EQ. 1) GO TO 20
              JM1 = J - 1
              DO 10 K = 1, JM1
                   IK = I0 + K
                   JK = J0 + K
                   T = T + L(IK)*L(JK)
 10                CONTINUE
 20           IJ = I0 + J
              J0 = J0 + J
              T = (A(IJ) - T) / L(J0)
              L(IJ) = T
              TD = TD + T*T
 30           CONTINUE
 40      I0 = I0 + I
         T = A(I0) - TD
         IF (T .LE. ZERO) GO TO 60
         L(I0) =  SQRT(T)
 50      CONTINUE
C
      IRC = 0
      GO TO 999
C
 60   L(I0) = T
      IRC = I
C
 999  RETURN
C
C  ***  LAST CARD OF LSQRT  ***
      END
