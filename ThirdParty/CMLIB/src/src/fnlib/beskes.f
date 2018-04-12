      SUBROUTINE BESKES(XNU,X,NIN,BKE)
C***BEGIN PROLOGUE  BESKES
C***DATE WRITTEN   770601   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  C10B3
C***KEYWORDS  BESSEL FUNCTION,EXPONENTIALLY SCALED,FRACTIONAL ORDER,
C             MODIFIED BESSEL FUNCTION,SEQUENCE,THIRD KIND
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes a sequence of exponentially scaled modified Bessel
C            functions of the third kind of fractional order.
C***DESCRIPTION
C
C BESKES computes a sequence of exponentially scaled
C (i.e., multipled by EXP(X)) modified Bessel
C functions of the third kind of order XNU + I at X, where X .GT. 0,
C XNU lies in (-1,1), and I = 0, 1, ... , NIN - 1, if NIN is positive
C and I = 0, -1, ... , NIN + 1, if NIN is negative.  On return, the
C vector BKE(.) contains the results at X for order starting at XNU.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  R1MACH,R9KNUS,XERROR
C***END PROLOGUE  BESKES
      DIMENSION BKE(*)
      DATA ALNBIG / 0. /
C***FIRST EXECUTABLE STATEMENT  BESKES
      IF (ALNBIG.EQ.0.) ALNBIG = ALOG (R1MACH(2))
C
      V = ABS(XNU)
      N = IABS(NIN)
C
      IF (V.GE.1.) CALL XERROR ( 'BESKES  ABS(XNU) MUST BE LT 1', 29,
     1  2, 2)
      IF (X.LE.0.) CALL XERROR ( 'BESKES  X IS LE 0', 17, 3, 2)
      IF (N.EQ.0) CALL XERROR ( 'BESKES  N THE NUMBER IN THE SEQUENCE IS
     1 0', 41, 4, 2)
C
      CALL R9KNUS (V, X, BKE(1), BKNU1, ISWTCH)
      IF (N.EQ.1) RETURN
C
      VINCR = SIGN (1.0, FLOAT(NIN))
      DIRECT = VINCR
      IF (XNU.NE.0.) DIRECT = VINCR*SIGN(1.0,XNU)
      IF (ISWTCH.EQ.1 .AND. DIRECT.GT.0.) CALL XERROR ( 'BESKES  X SO SM
     1ALL BESSEL K-SUB-XNU+1 OVERFLOWS', 47, 5, 2)
      BKE(2) = BKNU1
C
      IF (DIRECT.LT.0.) CALL R9KNUS (ABS(XNU+VINCR), X, BKE(2), BKNU1,
     1  ISWTCH)
      IF (N.EQ.2) RETURN
C
      VEND = ABS(XNU+FLOAT(NIN)) - 1.0
      IF ((VEND-0.5)*ALOG(VEND)+ 0.27 - VEND*(ALOG(X)-.694).GT.ALNBIG)
     1  CALL XERROR (  'BESKES  X SO SMALL OR ABS(NU) SO BIG THAT BESSEL
     2 K-SUB-NU OVERFLOWS', 67, 5, 2)
C
      V = XNU
      DO 10 I=3,N
        V = V + VINCR
        BKE(I) = 2.0*V*BKE(I-1)/X + BKE(I-2)
 10   CONTINUE
C
      RETURN
      END
