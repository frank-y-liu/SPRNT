      SUBROUTINE CROTG(CA,CB,C,S)
C***BEGIN PROLOGUE  CROTG
C***DATE WRITTEN   790101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D1B10
C***KEYWORDS  BLAS,COMPLEX,COMPLEX GIVENS TRANSFORMATION
C***AUTHOR  (NONE)
C***PURPOSE  Construct a complex Givens transformation.
C***DESCRIPTION
C
C     Complex Givens transformation
C
C      Construct the Givens transformation
C
C             (C    S)
C       G  =  (      ),  C**2 + CABS(S)**2 =1,
C             (-S   C)
C
C      which zeros the second entry of the complex 2-vector (CA,CB)**T
C
C    The quantity CA/CABS(CA)*NORM(CA,CB) overwrites CA in storage.
C
C    INPUT:
C        CA (Complex)
C        CB (Complex)
C
C    Output:
C        CA (Complex)      CA/CABS(CA)*NORM(CA,CB)
C        C  (Real)
C        S  (Complex)
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  CROTG
      COMPLEX CA,CB,S
      REAL C
      REAL NORM,SCALE
      COMPLEX ALPHA
C***FIRST EXECUTABLE STATEMENT  CROTG
      IF (CABS(CA) .NE. 0.) GO TO 10
         C = 0.
         S = (1.,0.)
         CA = CB
         GO TO 20
   10 CONTINUE
         SCALE = CABS(CA) + CABS(CB)
         NORM = SCALE * SQRT((CABS(CA/SCALE))**2 + (CABS(CB/SCALE))**2)
         ALPHA = CA /CABS(CA)
         C = CABS(CA) / NORM
         S = ALPHA * CONJG(CB) / NORM
         CA = ALPHA * NORM
   20 CONTINUE
      RETURN
      END
