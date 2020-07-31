      DOUBLE PRECISION FUNCTION D9PAK(Y,N)
C***BEGIN PROLOGUE  D9PAK
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  A6B
C***KEYWORDS  FNLIB,PACK
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Pack a base 2 exponent into a double precision floating
C            point number.
C***DESCRIPTION
C
C Pack a base 2 exponent into floating point number X.  This routine is
C almost the inverse of D9UPAK.  It is not exactly the inverse, because
C DABS(X) need not be between 0.5 and 1.0.  If both D9PAK and 2.d0**N
C were known to be in range we could compute
C               D9PAK = X *2.0d0**N
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH,D9UPAK,I1MACH,XERROR
C***END PROLOGUE  D9PAK
        DOUBLE PRECISION Y, A1N2B,A1N210,D1MACH
        DATA NMIN,NMAX / 2*0 /
        DATA A1N210 / 3.321928094 8873623478 7031942948 9 D0 /
C***FIRST EXECUTABLE STATEMENT  D9PAK
        IF(NMIN.NE.0) GO TO 10
        A1N2B = 1.0D0
        IF(I1MACH(7).NE.2) A1N2B=D1MACH(5)*A1N210
        NMIN = A1N2B*DBLE(FLOAT(I1MACH(15)))
        NMAX = A1N2B*DBLE(FLOAT(I1MACH(16)))
C
 10     CALL D9UPAK(Y,D9PAK,NY)
C
        NSUM=N+NY
        IF(NSUM.LT.NMIN)GO TO 40
        IF(NSUM.GT.NMAX)CALL XERROR ( 'D9PAK   PACKED NUMBER OVERFLOWS',
     1 31,1, 2)
C
        IF (NSUM.EQ.0) RETURN
        IF(NSUM.GT.0) GO TO 30
C
 20     D9PAK = 0.5D0*D9PAK
        NSUM=NSUM+1
        IF(NSUM.NE.0) GO TO 20
        RETURN
C
 30     D9PAK = 2.0D0*D9PAK
        NSUM=NSUM - 1
        IF (NSUM.NE.0) GO TO 30
        RETURN
C
 40     CALL XERROR ( 'D9PAK   PACKED NUMBER UNDERFLOWS', 32, 1, 1)
        D9PACK = 0.0D0
        RETURN
C
      END
