      REAL FUNCTION PCHST(ARG1,ARG2)
C***BEGIN PROLOGUE  PCHST
C***REFER TO  PCHCE,PCHCI,PCHCS,PCHIM
C***ROUTINES CALLED  (NONE)
C***REVISION DATE  870707   (YYMMDD)
C***DESCRIPTION
C
C         PCHST:  PCHIP Sign-Testing Routine.
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
C***END PROLOGUE  PCHST
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
C     To produce a double precision version, simply:
C        a. Change PCHST to DPCHST wherever it occurs,
C        b. Change all references to the Fortran intrinsics to their
C           double presision equivalents,
C        c. Change the real declarations to double precision, and
C        d. Change the constants  ZERO  and  ONE  to double precision.
C
C  DECLARE ARGUMENTS.
C
      REAL  ARG1, ARG2
C
C  DECLARE LOCAL VARIABLES.
C
      REAL  ONE, ZERO
      DATA  ZERO /0./,  ONE /1./
C
C  PERFORM THE TEST.
C
C***FIRST EXECUTABLE STATEMENT  PCHST
      PCHST = SIGN(ONE,ARG1) * SIGN(ONE,ARG2)
      IF ((ARG1.EQ.ZERO) .OR. (ARG2.EQ.ZERO))  PCHST = ZERO
C
      RETURN
C------------- LAST LINE OF PCHST FOLLOWS ------------------------------
      END
