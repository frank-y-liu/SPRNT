      DOUBLE PRECISION FUNCTION DPCHIA(N,X,F,D,INCFD,SKIP,A,B,IERR)
C***BEGIN PROLOGUE  DPCHIA
C***DATE WRITTEN   820730   (YYMMDD)
C***REVISION DATE  870707   (YYMMDD)
C***CATEGORY NO.  E3,H2A2
C***KEYWORDS  LIBRARY=SLATEC(PCHIP),
C             TYPE=DOUBLE PRECISION(PCHIA-S DPCHIA-D),
C             CUBIC HERMITE INTERPOLATION,NUMERICAL INTEGRATION,
C             QUADRATURE
C***AUTHOR  FRITSCH, F. N., (LLNL)
C             MATHEMATICS AND STATISTICS DIVISION
C             LAWRENCE LIVERMORE NATIONAL LABORATORY
C             P.O. BOX 808  (L-316)
C             LIVERMORE, CA  94550
C             FTS 532-4275, (415) 422-4275
C***PURPOSE  Evaluate the definite integral of a piecewise cubic
C            Hermite function over an arbitrary interval.
C***DESCRIPTION
C
C       **** Double Precision version of PCHIA ****
C
C          DPCHIA:  Piecewise Cubic Hermite Integrator, Arbitrary limits
C
C     Evaluates the definite integral of the cubic Hermite function
C     defined by  N, X, F, D  over the interval [A, B].
C
C     To provide compatibility with DPCHIM and DPCHIC, includes an
C     increment between successive values of the F- and D-arrays.
C
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C        PARAMETER  (INCFD = ...)
C        INTEGER  N, IERR
C        DOUBLE PRECISION  X(N), F(INCFD,N), D(INCFD,N), A, B
C        DOUBLE PRECISION  VALUE, DPCHIA
C        LOGICAL  SKIP
C
C        VALUE = DPCHIA (N, X, F, D, INCFD, SKIP, A, B, IERR)
C
C   Parameters:
C
C     VALUE -- (output) value of the requested integral.
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
C           been performed (say, in DPCHIM or DPCHIC).
C           SKIP will be set to .TRUE. on return with IERR.GE.0 .
C
C     A,B -- (input) the limits of integration.
C           NOTE:  There is no requirement that [A,B] be contained in
C                  [X(1),X(N)].  However, the resulting integral value
C                  will be highly suspect, if not.
C
C     IERR -- (output) error flag.
C           Normal return:
C              IERR = 0  (no errors).
C           Warning errors:
C              IERR = 1  if  A  is outside the interval [X(1),X(N)].
C              IERR = 2  if  B  is outside the interval [X(1),X(N)].
C              IERR = 3  if both of the above are true.  (Note that this
C                        means that either [A,B] contains data interval
C                        or the intervals do not intersect at all.)
C           "Recoverable" errors:
C              IERR = -1  if N.LT.2 .
C              IERR = -2  if INCFD.LT.1 .
C              IERR = -3  if the X-array is not strictly increasing.
C                (Value has not been computed in any of these cases.)
C               NOTE:  The above errors are checked in the order listed,
C                   and following arguments have **NOT** been validated.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DCHFIV,DPCHID,XERROR
C***END PROLOGUE  DPCHIA
C
C ----------------------------------------------------------------------
C
C  Change record:
C     82-08-04   Converted to SLATEC library version.
C     87-07-07   Corrected conversion to double precision.
C     87-07-07   Corrected XERROR calls for d.p. name(s).
C
C ----------------------------------------------------------------------
C
C  Programming notes:
C
C     To produce a single precision version, simply:
C        a. Change DPCHIA to PCHIA wherever it occurs,
C        b. Change DPCHID to PCHID wherever it occurs,
C        c. Change DCHFIV to CHFIV wherever it occurs,
C        d. Change the double precision declarations to real,  and
C        e. Change the constant  ZERO  to single precision.
C
C  DECLARE ARGUMENTS.
C
      INTEGER N, INCFD, IERR
      DOUBLE PRECISION  X(N), F(INCFD,N), D(INCFD,N), A, B
      LOGICAL  SKIP
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER I, IA, IB, IERD, IERV, IL, IR
      DOUBLE PRECISION  VALUE, XA, XB, ZERO
      DOUBLE PRECISION  DCHFIV, DPCHID
C
C  INITIALIZE.
C
      DATA  ZERO /0.D0/
C
C  VALIDITY-CHECK ARGUMENTS.
C
C***FIRST EXECUTABLE STATEMENT  DPCHIA
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
      IERR = 0
      IF ( (A.LT.X(1)) .OR. (A.GT.X(N)) )  IERR = IERR + 1
      IF ( (B.LT.X(1)) .OR. (B.GT.X(N)) )  IERR = IERR + 2
C
C  COMPUTE INTEGRAL VALUE.
C
      IF (A .EQ. B)  THEN
         VALUE = ZERO
      ELSE
         XA = DMIN1 (A, B)
         XB = DMAX1 (A, B)
         IF (XB .LE. X(2))  THEN
C           INTERVAL IS TO LEFT OF X(2), SO USE FIRST CUBIC.
C                   --------------------------------------------
            VALUE = DCHFIV (X(1),X(2), F(1,1),F(1,2),
     *                                D(1,1),D(1,2), A, B, IERV)
C                   --------------------------------------------
            IF (IERV .LT. 0)  GO TO 5004
         ELSE IF (XA .GE. X(N-1))  THEN
C           INTERVAL IS TO RIGHT OF X(N-1), SO USE LAST CUBIC.
C                   -----------------------------------------------
            VALUE = DCHFIV(X(N-1),X(N), F(1,N-1),F(1,N),
     *                                 D(1,N-1),D(1,N), A, B, IERV)
C                   -----------------------------------------------
            IF (IERV .LT. 0)  GO TO 5004
         ELSE
C           'NORMAL' CASE -- XA.LT.XB, XA.LT.X(N-1), XB.GT.X(2).
C      ......LOCATE IA AND IB SUCH THAT
C               X(IA-1).LT.XA.LE.X(IA).LE.X(IB).LE.XB.LE.X(IB+1)
            IA = 1
            DO 10  I = 1, N-1
               IF (XA .GT. X(I))  IA = I + 1
   10       CONTINUE
C             IA = 1 IMPLIES XA.LT.X(1) .  OTHERWISE,
C             IA IS LARGEST INDEX SUCH THAT X(IA-1).LT.XA,.
C
            IB = N
            DO 20  I = N, IA, -1
               IF (XB .LT. X(I))  IB = I - 1
   20       CONTINUE
C             IB = N IMPLIES XB.GT.X(N) .  OTHERWISE,
C             IB IS SMALLEST INDEX SUCH THAT XB.LT.X(IB+1) .
C
C     ......COMPUTE THE INTEGRAL.
            IERV = 0
            IF (IB .LT. IA)  THEN
C              THIS MEANS IB = IA-1 AND
C                 (A,B) IS A SUBSET OF (X(IB),X(IA)).
C                      ------------------------------------------------
               VALUE = DCHFIV (X(IB),X(IA), F(1,IB),F(1,IA),
     *                                     D(1,IB),D(1,IA), A, B, IERV)
C                      ------------------------------------------------
               IF (IERV .LT. 0)  GO TO 5004
            ELSE
C
C              FIRST COMPUTE INTEGRAL OVER (X(IA),X(IB)).
               IF (IB .EQ. IA)  THEN
                  VALUE = ZERO
               ELSE
C                         ---------------------------------------------
                  VALUE = DPCHID (N, X, F, D, INCFD, SKIP, IA, IB, IERD)
C                         ---------------------------------------------
                  IF (IERD .LT. 0)  GO TO 5005
               ENDIF
C
C              THEN ADD ON INTEGRAL OVER (XA,X(IA)).
               IF (XA .LT. X(IA))  THEN
                  IL = MAX0 (1, IA-1)
                  IR = IL + 1
C                                 -------------------------------------
                  VALUE = VALUE + DCHFIV (X(IL),X(IR), F(1,IL),F(1,IR),
     *                                D(1,IL),D(1,IR), XA, X(IA), IERV)
C                                 -------------------------------------
                  IF (IERV .LT. 0)  GO TO 5004
               ENDIF
C
C              THEN ADD ON INTEGRAL OVER (X(IB),XB).
               IF (XB .GT. X(IB))  THEN
                  IR = MIN0 (IB+1, N)
                  IL = IR - 1
C                                 -------------------------------------
                  VALUE = VALUE + DCHFIV (X(IL),X(IR), F(1,IL),F(1,IR),
     *                                D(1,IL),D(1,IR), X(IB), XB, IERV)
C                                 -------------------------------------
                  IF (IERV .LT. 0)  GO TO 5004
               ENDIF
C
C              FINALLY, ADJUST SIGN IF NECESSARY.
               IF (A .GT. B)  VALUE = -VALUE
            ENDIF
         ENDIF
      ENDIF
C
C  NORMAL RETURN.
C
      DPCHIA = VALUE
      RETURN
C
C  ERROR RETURNS.
C
 5001 CONTINUE
C     N.LT.2 RETURN.
      IERR = -1
      CALL XERROR ('DPCHIA -- NUMBER OF DATA POINTS LESS THAN TWO'
     *           , 0, IERR, 1)
      RETURN
C
 5002 CONTINUE
C     INCFD.LT.1 RETURN.
      IERR = -2
      CALL XERROR ('DPCHIA -- INCREMENT LESS THAN ONE'
     *           , 0, IERR, 1)
      RETURN
C
 5003 CONTINUE
C     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      CALL XERROR ('DPCHIA -- X-ARRAY NOT STRICTLY INCREASING'
     *           , 0, IERR, 1)
      RETURN
C
 5004 CONTINUE
C     TROUBLE IN DCHFIV.  (SHOULD NEVER OCCUR.)
      IERR = -4
      CALL XERROR ('DPCHIA -- TROUBLE IN DCHFIV'
     *           , 0, IERR, 1)
      RETURN
C
 5005 CONTINUE
C     TROUBLE IN DPCHID.  (SHOULD NEVER OCCUR.)
      IERR = -5
      CALL XERROR ('DPCHIA -- TROUBLE IN DPCHID'
     *           , 0, IERR, 1)
      RETURN
C------------- LAST LINE OF DPCHIA FOLLOWS -----------------------------
      END
