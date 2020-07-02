      REAL FUNCTION PCHID(N,X,F,D,INCFD,SKIP,IA,IB,IERR)
C***BEGIN PROLOGUE  PCHID
C***DATE WRITTEN   820723   (YYMMDD)
C***REVISION DATE  870707   (YYMMDD)
C***CATEGORY NO.  E1B,H2A2
C***KEYWORDS  LIBRARY=SLATEC(PCHIP),
C             TYPE=SINGLE PRECISION(PCHID-S DPCHID-D),
C             CUBIC HERMITE INTERPOLATION,NUMERICAL INTEGRATION,
C             QUADRATURE
C***AUTHOR  FRITSCH, F. N., (LLNL)
C             MATHEMATICS AND STATISTICS DIVISION
C             LAWRENCE LIVERMORE NATIONAL LABORATORY
C             P.O. BOX 808  (L-316)
C             LIVERMORE, CA  94550
C             FTS 532-4275, (415) 422-4275
C***PURPOSE  Evaluate the definite integral of a piecewise cubic
C            Hermite function over an interval whose endpoints are
C            data points.
C***DESCRIPTION
C
C          PCHID:  Piecewise Cubic Hermite Integrator, Data limits
C
C     Evaluates the definite integral of the cubic Hermite function
C     defined by  N, X, F, D  over the interval [X(IA), X(IB)].
C
C     To provide compatibility with PCHIM and PCHIC, includes an
C     increment between successive values of the F- and D-arrays.
C
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C        PARAMETER  (INCFD = ...)
C        INTEGER  N, IA, IB, IERR
C        REAL  X(N), F(INCFD,N), D(INCFD,N)
C        LOGICAL  SKIP
C
C        VALUE = PCHID (N, X, F, D, INCFD, SKIP, IA, IB, IERR)
C
C   Parameters:
C
C     VALUE -- (output) value of the requested integral.
C
C     N -- (input) number of data points.  (Error return if N.LT.2 .)
C
C     X -- (input) real array of independent variable values.  The
C           elements of X must be strictly increasing:
C                X(I-1) .LT. X(I),  I = 2(1)N.
C           (Error return if not.)
C
C     F -- (input) real array of function values.  F(1+(I-1)*INCFD) is
C           the value corresponding to X(I).
C
C     D -- (input) real array of derivative values.  D(1+(I-1)*INCFD) is
C           the value corresponding to X(I).
C
C     INCFD -- (input) increment between successive values in F and D.
C           (Error return if  INCFD.LT.1 .)
C
C     SKIP -- (input/output) logical variable which should be set to
C           .TRUE. if the user wishes to skip checks for validity of
C           preceding parameters, or to .FALSE. otherwise.
C           This will save time in case these checks have already
C           been performed (say, in PCHIM or PCHIC).
C           SKIP will be set to .TRUE. on return with IERR = 0 or -4.
C
C     IA,IB -- (input) indices in X-array for the limits of integration.
C           both must be in the range [1,N].  (Error return if not.)
C           No restrictions on their relative values.
C
C     IERR -- (output) error flag.
C           Normal return:
C              IERR = 0  (no errors).
C           "Recoverable" errors:
C              IERR = -1  if N.LT.2 .
C              IERR = -2  if INCFD.LT.1 .
C              IERR = -3  if the X-array is not strictly increasing.
C              IERR = -4  if IA or IB is out of range.
C                (Value has not been computed in any of these cases.)
C               NOTE:  The above errors are checked in the order listed,
C                   and following arguments have **NOT** been validated.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  XERROR
C***END PROLOGUE  PCHID
C
C ----------------------------------------------------------------------
C
C  Change record:
C     82-08-04   Converted to SLATEC library version.
C
C ----------------------------------------------------------------------
C
C  Programming notes:
C
C     To produce a double precision version, simply:
C        a. Change PCHID to DPCHID wherever it occurs,
C        b. Change the real declarations to double precision,  and
C        c. Change the constants ZERO, HALF, SIX to double precision.
C
C  DECLARE ARGUMENTS.
C
      INTEGER  N, INCFD, IA, IB, IERR
      REAL  X(N), F(INCFD,N), D(INCFD,N)
      LOGICAL  SKIP
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER  I, IUP, LOW
      REAL  H, HALF, SIX, SUM, VALUE, ZERO
C
C  INITIALIZE.
C
      DATA  ZERO /0./,  HALF /0.5/,  SIX /6./
C
C  VALIDITY-CHECK ARGUMENTS.
C
C***FIRST EXECUTABLE STATEMENT  PCHID
      IF (SKIP)  GO TO 5
C
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  I = 2, N
         IF ( X(I).LE.X(I-1) )  GO TO 5003
    1 CONTINUE
C
C  FUNCTION DEFINITION IS OK, GO ON.
C
    5 CONTINUE
      SKIP = .TRUE.
      IF ((IA.LT.1) .OR. (IA.GT.N))  GO TO 5004
      IF ((IB.LT.1) .OR. (IB.GT.N))  GO TO 5004
      IERR = 0
C
C  COMPUTE INTEGRAL VALUE.
C
      IF (IA .EQ. IB)  THEN
         VALUE = ZERO
      ELSE
         LOW = MIN0(IA, IB)
         IUP = MAX0(IA, IB) - 1
         SUM = ZERO
         DO 10  I = LOW, IUP
            H = X(I+1) - X(I)
            SUM = SUM + H*( (F(1,I) + F(1,I+1)) +
     *                      (D(1,I) - D(1,I+1))*(H/SIX) )
   10    CONTINUE
         VALUE = HALF * SUM
         IF (IA .GT. IB)  VALUE = -VALUE
      ENDIF
C
C  NORMAL RETURN.
C
      PCHID = VALUE
      RETURN
C
C  ERROR RETURNS.
C
 5001 CONTINUE
C     N.LT.2 RETURN.
      IERR = -1
      CALL XERROR ('PCHID -- NUMBER OF DATA POINTS LESS THAN TWO'
     *           , 44, IERR, 1)
      RETURN
C
 5002 CONTINUE
C     INCFD.LT.1 RETURN.
      IERR = -2
      CALL XERROR ('PCHID -- INCREMENT LESS THAN ONE'
     *           , 32, IERR, 1)
      RETURN
C
 5003 CONTINUE
C     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      CALL XERROR ('PCHID -- X-ARRAY NOT STRICTLY INCREASING'
     *           , 40, IERR, 1)
      RETURN
C
 5004 CONTINUE
C     IA OR IB OUT OF RANGE RETURN.
      IERR = -4
      CALL XERROR ('PCHID -- IA OR IB OUT OF RANGE'
     *           , 30, IERR, 1)
      RETURN
C------------- LAST LINE OF PCHID FOLLOWS ------------------------------
      END
