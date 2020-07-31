      SUBROUTINE DPCHMC(N,X,F,D,INCFD,SKIP,ISMON,IERR)
C***BEGIN PROLOGUE  DPCHMC
C***DATE WRITTEN   820518   (YYMMDD)
C***REVISION DATE  870707   (YYMMDD)
C***CATEGORY NO.  E1B
C***KEYWORDS  LIBRARY=SLATEC(PCHIP),
C             TYPE=DOUBLE PRECISION(PCHMC-S DPCHMC-D),
C             CUBIC HERMITE INTERPOLATION,MONOTONE INTERPOLATION,
C             PIECEWISE CUBIC INTERPOLATION
C***AUTHOR  FRITSCH, F. N., (LLNL)
C             MATHEMATICS AND STATISTICS DIVISION
C             LAWRENCE LIVERMORE NATIONAL LABORATORY
C             P.O. BOX 808  (L-316)
C             LIVERMORE, CA  94550
C             FTS 532-4275, (415) 422-4275
C***PURPOSE  Check a cubic Hermite function for monotonicity.
C***DESCRIPTION
C
C       **** Double Precision version of PCHMC ****
C
C          DPCHMC:  Piecewise Cubic Hermite Monotonicity Checker.
C
C     Checks the cubic Hermite function defined by  N, X, F, D  for
C     monotonicity.
C
C     To provide compatibility with DPCHIM and DPCHIC, includes an
C     increment between successive values of the F- and D-arrays.
C
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C        PARAMETER  (INCFD = ...)
C        INTEGER  N, ISMON(N), IERR
C        DOUBLE PRECISION  X(N), F(INCFD,N), D(INCFD,N)
C        LOGICAL  SKIP
C
C        CALL  DPCHMC (N, X, F, D, INCFD, SKIP, ISMON, IERR)
C
C   Parameters:
C
C     N -- (input) number of data points.  (Error return if N.LT.2 .)
C
C     X -- (input) real*8 array of independent variable values.  The
C           elements of X must be strictly increasing:
C                X(I-1) .LT. X(I),  I = 2(1)N.
C           (Error return if not.)
C
C     F -- (input) real*8 array of function values.  F(1+(I-1)*INCFD) is
C           the value corresponding to X(I).
C
C     D -- (input) real*8 array of derivative values.  D(1+(I-1)*INCFD)
C           is the value corresponding to X(I).
C
C     INCFD -- (input) increment between successive values in F and D.
C           (Error return if  INCFD.LT.1 .)
C
C     SKIP -- (input/output) logical variable which should be set to
C           .TRUE. if the user wishes to skip checks for validity of
C           preceding parameters, or to .FALSE. otherwise.
C           This will save time in case these checks have already
C           been performed.
C           SKIP will be set to .TRUE. on normal return.
C
C     ISMON -- (output) integer array indicating on which intervals the
C           PCH function defined by  N, X, F, D  is monotonic.
C           For data interval [X(I),X(I+1)],
C             ISMON(I) = -1  if function is strictly decreasing;
C             ISMON(I) =  0  if function is constant;
C             ISMON(I) =  1  if function is strictly increasing;
C             ISMON(I) =  2  if function is non-monotonic;
C             ISMON(I) =  3  if unable to determine.  (This means that
C                            the D-values are near the boundary of the
C                            monotonicity region.  A small increase pro-
C                            duces non-monotonicity; decrease, strict
C                            monotonicity.)
C           The above applies to I=1(1)N-1.  ISMON(N) indicates whether
C              the entire function is monotonic on [X(1),X(N)].
C
C     IERR -- (output) error flag.
C           Normal return:
C              IERR = 0  (no errors).
C           "Recoverable" errors:
C              IERR = -1  if N.LT.2 .
C              IERR = -2  if INCFD.LT.1 .
C              IERR = -3  if the X-array is not strictly increasing.
C          (The ISMON-array has not been changed in any of these cases.)
C               NOTE:  The above errors are checked in the order listed,
C                   and following arguments have **NOT** been validated.
C
C***REFERENCES  F.N.FRITSCH AND R.E.CARLSON, 'MONOTONE PIECEWISE CUBIC
C                 INTERPOLATION,' SIAM J.NUMER.ANAL. 17, 2 (APRIL 1980),
C                 238-246.
C***ROUTINES CALLED  DCHFMC,XERROR
C***END PROLOGUE  DPCHMC
C
C ----------------------------------------------------------------------
C
C  Change record:
C     82-08-04   Converted to SLATEC library version.
C     83-12-01   Reversed order of subscripts of F and D, so that the
C                routine will work properly when INCFD.GT.1  .
C     87-07-07   Corrected XERROR calls for d.p. name(s).
C
C ----------------------------------------------------------------------
C
C  Programming notes:
C
C     To produce a single precision version, simply:
C        a. Change DPCHMC to PCHMC wherever it occurs,
C        b. Change DCHFMC to CHFMC wherever it occurs, and
C        c. Change the double precision declarations to real.
C
C  DECLARE ARGUMENTS.
C
      INTEGER N, INCFD, ISMON(N), IERR
      DOUBLE PRECISION  X(N), F(INCFD,N), D(INCFD,N)
      LOGICAL  SKIP
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER I, NSEG
      DOUBLE PRECISION  DELTA
      INTEGER DCHFMC
C
C  VALIDITY-CHECK ARGUMENTS.
C
C***FIRST EXECUTABLE STATEMENT  DPCHMC
      IF (SKIP)  GO TO 5
C
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  I = 2, N
         IF ( X(I).LE.X(I-1) )  GO TO 5003
    1 CONTINUE
      SKIP = .TRUE.
C
C  FUNCTION DEFINITION IS OK -- GO ON.
C
    5 CONTINUE
      NSEG = N - 1
      DO 90  I = 1, NSEG
         DELTA = (F(1,I+1)-F(1,I))/(X(I+1)-X(I))
C                   -------------------------------
         ISMON(I) = DCHFMC (D(1,I), D(1,I+1), DELTA)
C                   -------------------------------
         IF (I .EQ. 1)  THEN
            ISMON(N) = ISMON(1)
         ELSE
C           NEED TO FIGURE OUT CUMULATIVE MONOTONICITY FROM FOLLOWING
C           'MULTIPLICATION TABLE'--
C
C                    *      I S M O N (I)
C                     *  -1   0   1   2   3
C               I      *--------------------*
C               S   -1 I -1  -1   2   2   3 I
C               M    0 I -1   0   1   2   3 I
C               O    1 I  2   1   1   2   3 I
C               N    2 I  2   2   2   2   2 I
C              (N)   3 I  3   3   3   2   3 I
C                      *--------------------*
C
C           IF EQUAL OR ALREADY DECLARED NONMONOTONIC, NO CHANGE NEEDED.
            IF ((ISMON(I).NE.ISMON(N)) .AND. (ISMON(N).NE.2))  THEN
               IF ( MAX0(ISMON(I), ISMON(N)) .GT. 1)  THEN
C                 AT LEAST ONE IS EITHER 'NO' OR 'MAYBE'.
                  IF (ISMON(I) .EQ. 2)  THEN
                     ISMON(N) = 2
                  ELSE
                     ISMON(N) = 3
                  ENDIF
               ELSE IF (ISMON(I)*ISMON(N) .LT. 0)  THEN
C                 BOTH MONOTONIC, BUT IN OPPOSITE SENSES.
                  ISMON(N) = 2
               ELSE
C                 AT THIS POINT, ONE IS ZERO, THE OTHER IS +-1.
                  ISMON(N) = ISMON(N) + ISMON(I)
               ENDIF
            ENDIF
         ENDIF
   90 CONTINUE
C
C  NORMAL RETURN.
C
      IERR = 0
      RETURN
C
C  ERROR RETURNS.
C
 5001 CONTINUE
C     N.LT.2 RETURN.
      IERR = -1
      CALL XERROR ('DPCHMC -- NUMBER OF DATA POINTS LESS THAN TWO'
     *           , 0, IERR, 1)
      RETURN
C
 5002 CONTINUE
C     INCFD.LT.1 RETURN.
      IERR = -2
      CALL XERROR ('DPCHMC -- INCREMENT LESS THAN ONE'
     *           , 0, IERR, 1)
      RETURN
C
 5003 CONTINUE
C     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      CALL XERROR ('DPCHMC -- X-ARRAY NOT STRICTLY INCREASING'
     *           , 0, IERR, 1)
      RETURN
C------------- LAST LINE OF DPCHMC FOLLOWS -----------------------------
      END
