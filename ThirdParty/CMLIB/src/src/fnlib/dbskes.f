      SUBROUTINE DBSKES(XNU,X,NIN,BKE)
C***BEGIN PROLOGUE  DBSKES
C***DATE WRITTEN   770601   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  C10B3
C***KEYWORDS  BESSEL FUNCTION,DOUBLE PRECISION,EXPONENTIALLY SCALED,
C             FRACTIONAL ORDER,MODIFIED BESSEL FUNCTION,SEQUENCE,
C             SPECIAL FUNCTION,THIRD KIND
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes a d.p. sequence of exponentially scaled modified
C            Bessel functions of the third kind of fractional order.
C***DESCRIPTION
C
C DBSKES(XNU,X,NIN,BKE) computes a double precision sequence
C of exponentially scaled modified Bessel functions
C of the third kind of order XNU + I at X, where X .GT. 0,
C XNU lies in (-1,1), and I = 0, 1, ... , NIN - 1, if NIN is positive
C and I = 0, -1, ... , NIN + 1, if NIN is negative.  On return, the
C vector BKE(.) contains the results at X for order starting at XNU.
C XNU, X, and BKE are double precison.  NIN is integer.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH,D9KNUS,XERROR
C***END PROLOGUE  DBSKES
      DOUBLE PRECISION XNU, X, BKE(*), BKNU1, V, VINCR, VEND, ALNBIG,
     1  D1MACH
      DATA ALNBIG / 0.D0 /
C***FIRST EXECUTABLE STATEMENT  DBSKES
      IF (ALNBIG.EQ.0.D0) ALNBIG = DLOG (D1MACH(2))
C
      V = DABS(XNU)
      N = IABS(NIN)
C
      IF (V.GE.1.D0) CALL XERROR ( 'DBSKES  ABS(XNU) MUST BE LT 1', 29,
     1  2, 2)
      IF (X.LE.0.D0) CALL XERROR ( 'DBSKES  X IS LE 0', 17, 3, 2)
      IF (N.EQ.0) CALL XERROR ( 'DBSKES  N THE NUMBER IN THE SEQUENCE IS
     1 0', 41, 4, 2)
C
      CALL D9KNUS (V, X, BKE(1), BKNU1, ISWTCH)
      IF (N.EQ.1) RETURN
C
      VINCR = SIGN (1.0, FLOAT(NIN))
      DIRECT = VINCR
      IF (XNU.NE.0.D0) DIRECT = VINCR*DSIGN(1.D0, XNU)
      IF (ISWTCH.EQ.1 .AND. DIRECT.GT.0.) CALL XERROR ( 'DBSKES  X SO SM
     1ALL BESSEL K-SUB-XNU+1 OVERFLOWS', 47, 5, 2)
      BKE(2) = BKNU1
C
      IF (DIRECT.LT.0.) CALL D9KNUS (DABS(XNU+VINCR), X, BKE(2), BKNU1,
     1  ISWTCH)
      IF (N.EQ.2) RETURN
C
      VEND = DABS (XNU+DBLE(FLOAT(NIN))) - 1.0D0
      IF ((VEND-.5D0)*DLOG(VEND)+0.27D0-VEND*(DLOG(X)-.694D0).GT.ALNBIG)
     1  CALL XERROR (  'DBSKES  X SO SMALL OR ABS(NU) SO BIG THAT BESSEL
     2 K-SUB-NU OVERFLOWS', 67, 5, 2)
C
      V = XNU
      DO 10 I=3,N
        V = V + VINCR
        BKE(I) = 2.0D0*V*BKE(I-1)/X + BKE(I-2)
 10   CONTINUE
C
      RETURN
      END
