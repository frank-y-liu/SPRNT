      FUNCTION R9PAK(Y,N)
C***BEGIN PROLOGUE  R9PAK
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  A6B
C***KEYWORDS  FNLIB,PACK
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Pack a base 2 exponent into a floating point number.
C***DESCRIPTION
C
C Pack a base 2 exponent into floating point number Y.  This
C routine is almost the inverse of R9UPAK.  It is not exactly
C the inverse, because ABS(X) need not be between 0.5 and
C 1.0.  If both R9PAK and 2.0**N were known to be in range, we
C could compute
C       R9PAK = Y * 2.0**N .
C***REFERENCES  (NONE)
C***ROUTINES CALLED  I1MACH,R1MACH,R9UPAK,XERROR
C***END PROLOGUE  R9PAK
       DATA NMIN, NMAX / 2*0 /
       DATA A1N210 / 3.321928094 887362 E0/
C***FIRST EXECUTABLE STATEMENT  R9PAK
       IF(NMIN.NE.0) GO TO 10
       A1N2B = 1.0
       IF (I1MACH(7).NE.2) A1N2B = R1MACH(5)*A1N210
       NMIN =A1N2B*FLOAT(I1MACH(12))
       NMAX = A1N2B*FLOAT(I1MACH(13))
C
 10    CALL R9UPAK(Y,R9PAK,NY)
C
       NSUM = N + NY
       IF (NSUM.LT.NMIN) GO TO 40
       IF (NSUM.GT.NMAX) CALL XERROR ( 'R9PAK  PACKED NUMBER OVERFLOWS',
     1 30, 2, 2)
C
       IF (NSUM.EQ.0) RETURN
       IF (NSUM.GT.0) GO TO 30
C
 20    R9PAK = 0.5*R9PAK
       NSUM = NSUM + 1
       IF(NSUM.NE.0) GO TO 20
       RETURN
C
30     R9PAK = 2.0*R9PAK
       NSUM = NSUM - 1
       IF(NSUM.NE.0) GO TO 30
       RETURN
C
40     CALL XERROR( 'R9PAK  PACKED NUMBER UNDERFLOWS', 31, 1, 1)
       R9PAK = 0.0
       RETURN
C
      END
