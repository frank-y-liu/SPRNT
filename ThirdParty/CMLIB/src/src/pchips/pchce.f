      SUBROUTINE PCHCE(IC,VC,N,X,H,SLOPE,D,INCFD,IERR)
C***BEGIN PROLOGUE  PCHCE
C***REFER TO  PCHIC
C***ROUTINES CALLED  PCHDF,PCHST,XERROR
C***REVISION DATE  870707   (YYMMDD)
C***DESCRIPTION
C
C          PCHCE:  PCHIC End Derivative Setter.
C
C    Called by PCHIC to set end derivatives as requested by the user.
C    It must be called after interior derivative values have been set.
C                      -----
C
C    To facilitate two-dimensional applications, includes an increment
C    between successive values of the D-array.
C
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C        PARAMETER  (INCFD = ...)
C        INTEGER  IC(2), N, IERR
C        REAL  VC(2), X(N), H(N), SLOPE(N), D(INCFD,N)
C
C        CALL  PCHCE (IC, VC, N, X, H, SLOPE, D, INCFD, IERR)
C
C   Parameters:
C
C     IC -- (input) integer array of length 2 specifying desired
C           boundary conditions:
C           IC(1) = IBEG, desired condition at beginning of data.
C           IC(2) = IEND, desired condition at end of data.
C           ( see prologue to PCHIC for details. )
C
C     VC -- (input) real array of length 2 specifying desired boundary
C           values.  VC(1) need be set only if IC(1) = 2 or 3 .
C                    VC(2) need be set only if IC(2) = 2 or 3 .
C
C     N -- (input) number of data points.  (assumes N.GE.2)
C
C     X -- (input) real array of independent variable values.  (the
C           elements of X are assumed to be strictly increasing.)
C
C     H -- (input) real array of interval lengths.
C     SLOPE -- (input) real array of data slopes.
C           If the data are (X(I),Y(I)), I=1(1)N, then these inputs are:
C                  H(I) =  X(I+1)-X(I),
C              SLOPE(I) = (Y(I+1)-Y(I))/H(I),  I=1(1)N-1.
C
C     D -- (input) real array of derivative values at the data points.
C           The value corresponding to X(I) must be stored in
C                D(1+(I-1)*INCFD),  I=1(1)N.
C          (output) the value of D at X(1) and/or X(N) is changed, if
C           necessary, to produce the requested boundary conditions.
C           no other entries in D are changed.
C
C     INCFD -- (input) increment between successive values in D.
C           This argument is provided primarily for 2-D applications.
C
C     IERR -- (output) error flag.
C           Normal return:
C              IERR = 0  (no errors).
C           Warning errors:
C              IERR = 1  if IBEG.LT.0 and D(1) had to be adjusted for
C                        monotonicity.
C              IERR = 2  if IEND.LT.0 and D(1+(N-1)*INCFD) had to be
C                        adjusted for monotonicity.
C              IERR = 3  if both of the above are true.
C
C    -------
C    WARNING:  This routine does no validity-checking of arguments.
C    -------
C
C  Fortran intrinsics used:  ABS, IABS.
C
C***END PROLOGUE  PCHCE
C
C ----------------------------------------------------------------------
C
C  Programmed by:  Fred N. Fritsch,  FTS 532-4275, (415) 422-4275,
C                  Mathematics and Statistics Division,
C                  Lawrence Livermore National Laboratory.
C
C  Change record:
C     82-08-05   Converted to SLATEC library version.
C     87-07-07   Minor corrections made to prologue..
C
C ----------------------------------------------------------------------
C
C  Programming notes:
C
C     1. The function  PCHST(ARG1,ARG2)  is assumed to return zero if
C        either argument is zero, +1 if they are of the same sign, and
C        -1 if they are of opposite sign.
C     2. To produce a double precision version, simply:
C        a. Change PCHCE to DPCHCE wherever it occurs,
C        b. Change PCHDF to DPCHCE wherever it occurs,
C        c. Change PCHST to DPCHST wherever it occurs,
C        d. Change all references to the Fortran intrinsics to their
C           double presision equivalents,
C        e. Change the real declarations to double precision, and
C        f. Change the constants ZERO, HALF, ... to double precision.
C     3. One could reduce the number of arguments and amount of local
C        storage, at the expense of reduced code clarity, by passing in
C        the array WK (rather than splitting it into H and SLOPE) and
C        increasing its length enough to incorporate STEMP and XTEMP.
C     4. The two monotonicity checks only use the sufficient conditions.
C        Thus, it is possible (but unlikely) for a boundary condition to
C        be changed, even though the original interpolant was monotonic.
C        (At least the result is a continuous function of the data.)
C
C  DECLARE ARGUMENTS.
C
      INTEGER  IC(2), N, INCFD, IERR
      REAL  VC(2), X(N), H(N), SLOPE(N), D(INCFD,N)
C
C  DELCARE LOCAL VARIABLES.
C
      INTEGER  IBEG, IEND, IERF, INDEX, J, K
      REAL  HALF, STEMP(3), THREE, TWO, XTEMP(4), ZERO
      REAL  PCHDF, PCHST
C
C  INITIALIZE.
C
      DATA  ZERO /0./,  HALF /0.5/,  TWO /2./,  THREE /3./
C
C***FIRST EXECUTABLE STATEMENT  PCHCE
      IBEG = IC(1)
      IEND = IC(2)
      IERR = 0
C
C  SET TO DEFAULT BOUNDARY CONDITIONS IF N IS TOO SMALL.
C
      IF ( IABS(IBEG).GT.N )  IBEG = 0
      IF ( IABS(IEND).GT.N )  IEND = 0
C
C  TREAT BEGINNING BOUNDARY CONDITION.
C
      IF (IBEG .EQ. 0)  GO TO 2000
      K = IABS(IBEG)
      IF (K .EQ. 1)  THEN
C        BOUNDARY VALUE PROVIDED.
         D(1,1) = VC(1)
      ELSE IF (K .EQ. 2)  THEN
C        BOUNDARY SECOND DERIVATIVE PROVIDED.
         D(1,1) = HALF*( (THREE*SLOPE(1) - D(1,2)) - HALF*VC(1)*H(1) )
      ELSE IF (K .LT. 5)  THEN
C        USE K-POINT DERIVATIVE FORMULA.
C        PICK UP FIRST K POINTS, IN REVERSE ORDER.
         DO 10  J = 1, K
            INDEX = K-J+1
C           INDEX RUNS FROM K DOWN TO 1.
            XTEMP(J) = X(INDEX)
            IF (J .LT. K)  STEMP(J) = SLOPE(INDEX-1)
   10    CONTINUE
C                 -----------------------------
         D(1,1) = PCHDF (K, XTEMP, STEMP, IERF)
C                 -----------------------------
         IF (IERF .NE. 0)  GO TO 5001
      ELSE
C        USE 'NOT A KNOT' CONDITION.
         D(1,1) = ( THREE*(H(1)*SLOPE(2) + H(2)*SLOPE(1))
     *             - TWO*(H(1)+H(2))*D(1,2) - H(1)*D(1,3) ) / H(2)
      ENDIF
C
      IF (IBEG .GT. 0)  GO TO 2000
C
C  CHECK D(1,1) FOR COMPATIBILITY WITH MONOTONICITY.
C
      IF (SLOPE(1) .EQ. ZERO)  THEN
         IF (D(1,1) .NE. ZERO)  THEN
            D(1,1) = ZERO
            IERR = IERR + 1
         ENDIF
      ELSE IF ( PCHST(D(1,1),SLOPE(1)) .LT. ZERO)  THEN
         D(1,1) = ZERO
         IERR = IERR + 1
      ELSE IF ( ABS(D(1,1)) .GT. THREE*ABS(SLOPE(1)) )  THEN
         D(1,1) = THREE*SLOPE(1)
         IERR = IERR + 1
      ENDIF
C
C  TREAT END BOUNDARY CONDITION.
C
 2000 CONTINUE
      IF (IEND .EQ. 0)  GO TO 5000
      K = IABS(IEND)
      IF (K .EQ. 1)  THEN
C        BOUNDARY VALUE PROVIDED.
         D(1,N) = VC(2)
      ELSE IF (K .EQ. 2)  THEN
C        BOUNDARY SECOND DERIVATIVE PROVIDED.
         D(1,N) = HALF*( (THREE*SLOPE(N-1) - D(1,N-1)) +
     *                                           HALF*VC(2)*H(N-1) )
      ELSE IF (K .LT. 5)  THEN
C        USE K-POINT DERIVATIVE FORMULA.
C        PICK UP LAST K POINTS.
         DO 2010  J = 1, K
            INDEX = N-K+J
C           INDEX RUNS FROM N+1-K UP TO N.
            XTEMP(J) = X(INDEX)
            IF (J .LT. K)  STEMP(J) = SLOPE(INDEX)
 2010    CONTINUE
C                 -----------------------------
         D(1,N) = PCHDF (K, XTEMP, STEMP, IERF)
C                 -----------------------------
         IF (IERF .NE. 0)  GO TO 5001
      ELSE
C        USE 'NOT A KNOT' CONDITION.
         D(1,N) = ( THREE*(H(N-1)*SLOPE(N-2) + H(N-2)*SLOPE(N-1))
     *             - TWO*(H(N-1)+H(N-2))*D(1,N-1) - H(N-1)*D(1,N-2) )
     *                                                         / H(N-2)
      ENDIF
C
      IF (IEND .GT. 0)  GO TO 5000
C
C  CHECK D(1,N) FOR COMPATIBILITY WITH MONOTONICITY.
C
      IF (SLOPE(N-1) .EQ. ZERO)  THEN
         IF (D(1,N) .NE. ZERO)  THEN
            D(1,N) = ZERO
            IERR = IERR + 2
         ENDIF
      ELSE IF ( PCHST(D(1,N),SLOPE(N-1)) .LT. ZERO)  THEN
         D(1,N) = ZERO
         IERR = IERR + 2
      ELSE IF ( ABS(D(1,N)) .GT. THREE*ABS(SLOPE(N-1)) )  THEN
         D(1,N) = THREE*SLOPE(N-1)
         IERR = IERR + 2
      ENDIF
C
C  NORMAL RETURN.
C
 5000 CONTINUE
      RETURN
C
C  ERROR RETURN.
C
 5001 CONTINUE
C     ERROR RETURN FROM PCHDF.
C   *** THIS CASE SHOULD NEVER OCCUR ***
      IERR = -1
      CALL XERROR ('PCHCE -- ERROR RETURN FROM PCHDF'
     *           , 32, IERR, 1)
      RETURN
C------------- LAST LINE OF PCHCE FOLLOWS ------------------------------
      END
