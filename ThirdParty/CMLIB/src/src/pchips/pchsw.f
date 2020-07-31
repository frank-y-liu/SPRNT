      SUBROUTINE PCHSW(DFMAX,IEXTRM,D1,D2,H,SLOPE,IERR)
C***BEGIN PROLOGUE  PCHSW
C***REFER TO  PCHCS
C***ROUTINES CALLED  R1MACH,XERROR
C***REVISION DATE  870707   (YYMMDD)
C***DESCRIPTION
C
C         PCHSW:  PCHCS Switch Excursion Limiter.
C
C     Called by  PCHCS  to adjust D1 and D2 if necessary to insure that
C     the extremum on this interval is not further than DFMAX from the
C     extreme data value.
C
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C        INTEGER  IEXTRM, IERR
C        REAL  DFMAX, D1, D2, H, SLOPE
C
C        CALL  PCHSW (DFMAX, IEXTRM, D1, D2, H, SLOPE, IERR)
C
C   Parameters:
C
C     DFMAX -- (input) maximum allowed difference between F(IEXTRM) and
C           the cubic determined by derivative values D1,D2.  (assumes
C           DFMAX.GT.0.)
C
C     IEXTRM -- (input) index of the extreme data value.  (assumes
C           IEXTRM = 1 or 2 .  Any value .NE.1 is treated as 2.)
C
C     D1,D2 -- (input) derivative values at the ends of the interval.
C           (Assumes D1*D2 .LE. 0.)
C          (output) may be modified if necessary to meet the restriction
C           imposed by DFMAX.
C
C     H -- (input) interval length.  (Assumes  H.GT.0.)
C
C     SLOPE -- (input) data SLOPE on the interval.
C
C     IERR -- (output) error flag.  should be zero.
C           If IERR=-1, assumption on D1 and D2 is not satisfied.
C           If IERR=-2, quadratic equation locating extremum has
C                       negative descriminant (should never occur).
C
C    -------
C    WARNING:  This routine does no validity-checking of arguments.
C    -------
C
C
C  Fortran intrinsics used:  ABS, SIGN, SQRT.
C
C***END PROLOGUE  PCHSW
C
C ----------------------------------------------------------------------
C
C  Programmed by:  Fred N. Fritsch,  FTS 532-4275, (415) 422-4275,
C                  Mathematics and Statistics Division,
C                  Lawrence Livermore National Laboratory.
C
C  Change record:
C     82-08-05   Converted to SLATEC library version.
C     87-07-07   Replaced DATA statement for SMALL with a use of R1MACH.
C
C ----------------------------------------------------------------------
C
C  Programming notes:
C
C     To produce a double precision version, simply:
C        a. Change PCHSW to DPCHSW wherever it occurs,
C        b. Change PCHCS to DPCHCS wherever it occurs,
C        c. Change R1MACH to D1MACH wherever it occurs,
C        d. Change all references to the Fortran intrinsics to their
C           double precision equivalents,
C        e. Change the real declarations to double precision, and
C        f. Change constants ZERO, ONE, TWO, ... to double precision.
C
C  DECLARE ARGUMENTS.
C
      INTEGER  IEXTRM, IERR
      REAL  DFMAX, D1, D2, H, SLOPE
C
C  DECLARE LOCAL VARIABLES.
C
      REAL  CP, DMAX, FACT, LAMBDA, NU, ONE, PHI, RADCAL, RHO, SIGMA,
     *      SMALL, THAT, THIRD, THREE, TWO, ZERO
      DATA  ZERO /0./,  ONE /1./,  TWO /2./,  THREE /3./, FACT /100./
C        THIRD SHOULD BE SLIGHTLY LESS THAN 1/3.
      DATA  THIRD /0.33333/
C
C  NOTATION AND GENERAL REMARKS.
C
C     RHO IS THE RATIO OF THE DATA SLOPE TO THE DERIVATIVE BEING TESTED.
C     LAMBDA IS THE RATIO OF D2 TO D1.
C     THAT = T-HAT(RHO) IS THE NORMALIZED LOCATION OF THE EXTREMUM.
C     PHI IS THE NORMALIZED VALUE OF P(X)-F1 AT X = XHAT = X-HAT(RHO),
C           WHERE  THAT = (XHAT - X1)/H .
C        THAT IS, P(XHAT)-F1 = D*H*PHI,  WHERE D=D1 OR D2.
C     SIMILARLY,  P(XHAT)-F2 = D*H*(PHI-RHO) .
C
C      SMALL SHOULD BE A FEW ORDERS OF MAGNITUDE GREATER THAN MACHEPS.
C***FIRST EXECUTABLE STATEMENT  PCHSW
      SMALL = FACT*R1MACH(4)
C
C  DO MAIN CALCULATION.
C
      IF (D1 .EQ. ZERO)  THEN
C
C        SPECIAL CASE -- D1.EQ.ZERO .
C
C          IF D2 IS ALSO ZERO, THIS ROUTINE SHOULD NOT HAVE BEEN CALLED.
         IF (D2 .EQ. ZERO)  GO TO 5001
C
         RHO = SLOPE/D2
C          EXTREMUM IS OUTSIDE INTERVAL WHEN RHO .GE. 1/3 .
         IF (RHO .GE. THIRD)  GO TO 5000
         THAT = (TWO*(THREE*RHO-ONE)) / (THREE*(TWO*RHO-ONE))
         PHI = THAT**2 * ((THREE*RHO-ONE)/THREE)
C
C          CONVERT TO DISTANCE FROM F2 IF IEXTRM.NE.1 .
         IF (IEXTRM .NE. 1)  PHI = PHI - RHO
C
C          TEST FOR EXCEEDING LIMIT, AND ADJUST ACCORDINGLY.
         DMAX = DFMAX / (H*ABS(PHI))
         IF (ABS(D2) .GT. DMAX)  D2 = SIGN (DMAX, D2)
      ELSE
C
         RHO = SLOPE/D1
         LAMBDA = -D2/D1
         IF (D2 .EQ. ZERO)  THEN
C
C           SPECIAL CASE -- D2.EQ.ZERO .
C
C             EXTREMUM IS OUTSIDE INTERVAL WHEN RHO .GE. 1/3 .
            IF (RHO .GE. THIRD)  GO TO 5000
            CP = TWO - THREE*RHO
            NU = ONE - TWO*RHO
            THAT = ONE / (THREE*NU)
         ELSE
            IF (LAMBDA .LE. ZERO)  GO TO 5001
C
C           NORMAL CASE -- D1 AND D2 BOTH NONZERO, OPPOSITE SIGNS.
C
            NU = ONE - LAMBDA - TWO*RHO
            SIGMA = ONE - RHO
            CP = NU + SIGMA
            IF (ABS(NU) .GT. SMALL)  THEN
               RADCAL = (NU - (TWO*RHO+ONE))*NU + SIGMA**2
               IF (RADCAL .LT. ZERO)  GO TO 5002
               THAT = (CP - SQRT(RADCAL)) / (THREE*NU)
            ELSE
               THAT = ONE/(TWO*SIGMA)
            ENDIF
         ENDIF
         PHI = THAT*((NU*THAT - CP)*THAT + ONE)
C
C          CONVERT TO DISTANCE FROM F2 IF IEXTRM.NE.1 .
         IF (IEXTRM .NE. 1)  PHI = PHI - RHO
C
C          TEST FOR EXCEEDING LIMIT, AND ADJUST ACCORDINGLY.
         DMAX = DFMAX / (H*ABS(PHI))
         IF (ABS(D1) .GT. DMAX)  THEN
            D1 = SIGN (DMAX, D1)
            D2 = -LAMBDA*D1
         ENDIF
      ENDIF
C
C  NORMAL RETURN.
C
 5000 CONTINUE
      IERR = 0
      RETURN
C
C  ERROR RETURNS.
C
 5001 CONTINUE
C     D1 AND D2 BOTH ZERO, OR BOTH NONZERO AND SAME SIGN.
      IERR = -1
      CALL XERROR ('PCHSW -- D1 AND/OR D2 INVALID'
     *           , 29, IERR, 1)
      RETURN
C
 5002 CONTINUE
C     NEGATIVE VALUE OF RADICAL (SHOULD NEVER OCCUR).
      IERR = -2
      CALL XERROR ('PCHSW -- NEGATIVE RADICAL'
     *           , 25, IERR, 1)
      RETURN
C------------- LAST LINE OF PCHSW FOLLOWS ------------------------------
      END
