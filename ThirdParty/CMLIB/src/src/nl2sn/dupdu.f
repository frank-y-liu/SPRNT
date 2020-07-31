      SUBROUTINE DUPDU(D, HDIAG, IV, LIV, LV, N, V)
C
C  ***  UPDATE SCALE VECTOR D FOR HUMSL  ***
C
C  ***  PARAMETER DECLARATIONS  ***
C
      INTEGER LIV, LV, N
      INTEGER IV(LIV)
      REAL D(N), HDIAG(N), V(LV)
C
C  ***  LOCAL VARIABLES  ***
C
      INTEGER DTOLI, D0I, I
      REAL T, VDFAC
C
C  ***  INTRINSIC FUNCTIONS  ***
C/+
      REAL  ABS, AMAX1,  SQRT
C/
C  ***  SUBSCRIPTS FOR IV AND V  ***
C
      INTEGER DFAC, DTOL, DTYPE, NITER
C/6
      DATA DFAC/41/, DTOL/59/, DTYPE/16/, NITER/31/
C/7
C     PARAMETER (DFAC=41, DTOL=59, DTYPE=16, NITER=31)
C/
C
C-------------------------------  BODY  --------------------------------
C
      I = IV(DTYPE)
      IF (I .EQ. 1) GO TO 10
         IF (IV(NITER) .GT. 0) GO TO 999
C
 10   DTOLI = IV(DTOL)
      D0I = DTOLI + N
      VDFAC = V(DFAC)
      DO 20 I = 1, N
         T = AMAX1( SQRT( ABS(HDIAG(I))), VDFAC*D(I))
         IF (T .LT. V(DTOLI)) T = AMAX1(V(DTOLI), V(D0I))
         D(I) = T
         DTOLI = DTOLI + 1
         D0I = D0I + 1
 20      CONTINUE
C
 999  RETURN
C  ***  LAST CARD OF DUPDU FOLLOWS  ***
      END
