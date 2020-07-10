      FUNCTION INITS(OS,NOS,ETA)
C***BEGIN PROLOGUE  INITS
C***DATE WRITTEN   770401   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C3A2
C***KEYWORDS  INITIALIZE,ORTHOGONAL SERIES,SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Initializes an orthogonal series so that it defines the
C            number of terms to carry in the series to meet a specified
C            error.
C***DESCRIPTION
C
C Initialize the orthogonal series so that INITS is the number of terms
C needed to insure the error is no larger than ETA.  Ordinarily, ETA
C will be chosen to be one-tenth machine precision.
C
C             Input Arguments --
C OS     array of NOS coefficients in an orthogonal series.
C NOS    number of coefficients in OS.
C ETA    requested accuracy of series.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  XERROR
C***END PROLOGUE  INITS
      DIMENSION OS(NOS)
C***FIRST EXECUTABLE STATEMENT  INITS
      IF (NOS.LT.1) CALL XERROR ( 'INITS   NUMBER OF COEFFICIENTS LT 1',
     1 35, 2, 2)
C
      ERR = 0.
      DO 10 II=1,NOS
        I = NOS + 1 - II
        ERR = ERR + ABS(OS(I))
        IF (ERR.GT.ETA) GO TO 20
 10   CONTINUE
C
C     ... Extra clause added to the check below.
C     ... Machine eps is just smaller than the CS coefficients.
C
 20   CONTINUE
      IF ((I.EQ.NOS) .AND. (ERR .GT. 10.0E0*ETA)) 
     +   CALL XERROR ( 'INITS   ETA MAY BE TOO SMALL', 28,
     1  1, 2)
      INITS = I
C
      RETURN
      END
