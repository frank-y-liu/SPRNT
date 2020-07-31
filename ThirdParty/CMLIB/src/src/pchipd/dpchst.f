      DOUBLE PRECISION FUNCTION DPCHST(ARG1,ARG2)
C***BEGIN PROLOGUE  DPCHST
C***REFER TO  DPCHCE,DPCHCI,DPCHCS,DPCHIM
C***ROUTINES CALLED  (NONE)
C***REVISION DATE  870707   (YYMMDD)
C***DESCRIPTION
C
C         DPCHST:  DPCHIP Sign-Testing Routine.
C
C
C     Returns:
C        -1. if ARG1 and ARG2 are of opposite sign.
C         0. if either argument is zero.
C        +1. if ARG1 and ARG2 are of the same sign.
C
C     The object is to do this without multiplying ARG1*ARG2, to avoid
C     possible over/underflow problems.
C
C  Fortran intrinsics used:  SIGN.
C
C***END PROLOGUE  DPCHST
C
C ----------------------------------------------------------------------
C
C  Programmed by:  Fred N. Fritsch,  FTS 532-4275, (415) 422-4275,
C                  Mathematics and Statistics Division,
C                  Lawrence Livermore National Laboratory.
C
C  Change record:
C     82-08-05   Converted to SLATEC library version.
C
C ----------------------------------------------------------------------
C
C  Programming notes:
C
C     To produce a single precision version, simply:
C        a. Change DPCHST to PCHST wherever it occurs,
C        b. Change all references to the Fortran intrinsics to their
C           single presision equivalents,
C        c. Change the double precision declarations to real, and
C        d. Change the constants  ZERO  and  ONE  to single precision.
C
C  DECLARE ARGUMENTS.
C
      DOUBLE PRECISION  ARG1, ARG2
C
C  DECLARE LOCAL VARIABLES.
C
      DOUBLE PRECISION  ONE, ZERO
      DATA  ZERO /0.D0/,  ONE/1.D0/
C
C  PERFORM THE TEST.
C
C***FIRST EXECUTABLE STATEMENT  DPCHST
      DPCHST = DSIGN(ONE,ARG1) * DSIGN(ONE,ARG2)
      IF ((ARG1.EQ.ZERO) .OR. (ARG2.EQ.ZERO))  DPCHST = ZERO
C
      RETURN
C------------- LAST LINE OF DPCHST FOLLOWS -----------------------------
      END
