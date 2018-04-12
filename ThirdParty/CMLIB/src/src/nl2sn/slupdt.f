      SUBROUTINE SLUPDT(A, COSMIN, P, SIZE, STEP, U, W, WCHMTD, WSCALE,
     1                  Y)
C
C  ***  UPDATE SYMMETRIC  A  SO THAT  A * STEP = Y  ***
C  ***  (LOWER TRIANGLE OF  A  STORED ROWWISE       ***
C
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C  ***  PARAMETER DECLARATIONS  ***
C
      INTEGER P
      REAL A(*), COSMIN, SIZE, STEP(P), U(P), W(P),
     1                 WCHMTD(P), WSCALE, Y(P)
C     DIMENSION A(P*(P+1)/2)
C
C  ***  LOCAL VARIABLES  ***
C
      INTEGER I, J, K
      REAL DENMIN, SDOTWM, T, UI, WI
C
C     ***  CONSTANTS  ***
      REAL HALF, ONE, ZERO
C
C  ***  INTRINSIC FUNCTIONS  ***
C/+
      REAL  ABS, AMIN1
C/
C  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  ***
C
      EXTERNAL  SLVMUL
C
C/6
C     DATA HALF/0.5E+0/, ONE/1.E+0/, ZERO/0.E+0/
C/7
      PARAMETER (HALF=0.5E+0, ONE=1.E+0, ZERO=0.E+0)
C/
C
C-----------------------------------------------------------------------
C
      SDOTWM = SDOT(P,STEP,1,WCHMTD,1)
      DENMIN = COSMIN * SNRM2(P,STEP,1) * SNRM2(P,WCHMTD,1)
      WSCALE = ONE
      IF (DENMIN .NE. ZERO) WSCALE = AMIN1(ONE,  ABS(SDOTWM/DENMIN))
      T = ZERO
      IF (SDOTWM .NE. ZERO) T = WSCALE / SDOTWM
      DO 10 I = 1, P
 10      W(I) = T * WCHMTD(I)
      CALL SLVMUL(P, U, A, STEP)
      T = HALF * (SIZE * SDOT(P,STEP,1,U,1)  -  SDOT(P,STEP,1,Y,1))
      DO 20 I = 1, P
 20      U(I) = T*W(I) + Y(I) - SIZE*U(I)
C
C  ***  SET  A = A + U*(W**T) + W*(U**T)  ***
C
      K = 1
      DO 40 I = 1, P
         UI = U(I)
         WI = W(I)
         DO 30 J = 1, I
              A(K) = SIZE*A(K) + UI*W(J) + WI*U(J)
              K = K + 1
 30           CONTINUE
 40      CONTINUE
C
 999  RETURN
C  ***  LAST CARD OF SLUPDT FOLLOWS  ***
      END
