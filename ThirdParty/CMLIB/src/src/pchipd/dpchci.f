      SUBROUTINE DPCHCI(N,H,SLOPE,D,INCFD)
C***BEGIN PROLOGUE  DPCHCI
C***REFER TO  DPCHIC
C***ROUTINES CALLED  DPCHST
C***REVISION DATE  870707   (YYMMDD)
C***DESCRIPTION
C
C          DPCHCI:  DPCHIC Initial Derivative Setter.
C
C    Called by DPCHIC to set derivatives needed to determine a monotone
C    piecewise cubic Hermite interpolant to the data.
C
C    Default boundary conditions are provided which are compatible
C    with monotonicity.  If the data are only piecewise monotonic, the
C    interpolant will have an extremum at each point where monotonicity
C    switches direction.
C
C    To facilitate two-dimensional applications, includes an increment
C    between successive values of the D-array.
C
C    The resulting piecewise cubic Hermite function should be identical
C    (within roundoff error) to that produced by DPCHIM.
C
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C        PARAMETER  (INCFD = ...)
C        INTEGER  N
C        DOUBLE PRECISION  H(N), SLOPE(N), D(INCFD,N)
C
C        CALL  DPCHCI (N, H, SLOPE, D, INCFD)
C
C   Parameters:
C
C     N -- (input) number of data points.
C           If N=2, simply does linear interpolation.
C
C     H -- (input) real*8 array of interval lengths.
C     SLOPE -- (input) real*8 array of data slopes.
C           If the data are (X(I),Y(I)), I=1(1)N, then these inputs are:
C                  H(I) =  X(I+1)-X(I),
C              SLOPE(I) = (Y(I+1)-Y(I))/H(I),  I=1(1)N-1.
C
C     D -- (output) real*8 array of derivative values at data points.
C           If the data are monotonic, these values will determine a
C           a monotone cubic Hermite function.
C           The value corresponding to X(I) is stored in
C                D(1+(I-1)*INCFD),  I=1(1)N.
C           No other entries in D are changed.
C
C     INCFD -- (input) increment between successive values in D.
C           This argument is provided primarily for 2-D applications.
C
C    -------
C    WARNING:  This routine does no validity-checking of arguments.
C    -------
C
C  Fortran intrinsics used:  DABS, DMAX1, DMIN1.
C
C***END PROLOGUE  DPCHCI
C
C ----------------------------------------------------------------------
C
C  Programmed by:  Fred N. Fritsch,  FTS 532-4275, (415) 422-4275,
C                  Mathematics and Statistics Division,
C                  Lawrence Livermore National Laboratory.
C
C  Change record:
C     82-06-01   Modified end conditions to be continuous functions of
C                data when monotonicity switches in next interval.
C     82-06-02   1. Modified formulas so end conditions are less prone
C                   to over/underflow problems.
C                2. Minor modification to HSUM calculation.
C     82-08-05   Converted to SLATEC library version.
C
C ----------------------------------------------------------------------
C
C  Programming notes:
C
C     1. The function  DPCHST(ARG1,ARG2)  is assumed to return zero if
C        either argument is zero, +1 if they are of the same sign, and
C        -1 if they are of opposite sign.
C     2. To produce a single precision version, simply:
C        a. Change DPCHCI to PCHCI wherever it occurs,
C        b. Change DPCHST to PCHST wherever it occurs,
C        c. Change all references to the Fortran intrinsics to their
C           single presision equivalents,
C        d. Change the double precision declarations to real, and
C        e. Change the constants ZERO and THREE to single precision.
C
C  DECLARE ARGUMENTS.
C
      INTEGER N, INCFD
      DOUBLE PRECISION  H(N), SLOPE(N), D(INCFD,N)
C
C  DECLARE LOCAL VARIBLES.
C
      INTEGER I, NLESS1
      DOUBLE PRECISION  DEL1, DEL2, DMAX, DMIN, DRAT1, DRAT2, HSUM,
     *      HSUMT3, THREE, W1, W2, ZERO
      DOUBLE PRECISION  DPCHST
C
C  INITIALIZE.
C
      DATA  ZERO /0.D0/, THREE/3.D0/
C***FIRST EXECUTABLE STATEMENT  DPCHCI
      NLESS1 = N - 1
      DEL1 = SLOPE(1)
C
C  SPECIAL CASE N=2 -- USE LINEAR INTERPOLATION.
C
      IF (NLESS1 .GT. 1)  GO TO 10
      D(1,1) = DEL1
      D(1,N) = DEL1
      GO TO 5000
C
C  NORMAL CASE  (N .GE. 3).
C
   10 CONTINUE
      DEL2 = SLOPE(2)
C
C  SET D(1) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
C     SHAPE-PRESERVING.
C
      HSUM = H(1) + H(2)
      W1 = (H(1) + HSUM)/HSUM
      W2 = -H(1)/HSUM
      D(1,1) = W1*DEL1 + W2*DEL2
      IF ( DPCHST(D(1,1),DEL1) .LE. ZERO)  THEN
         D(1,1) = ZERO
      ELSE IF ( DPCHST(DEL1,DEL2) .LT. ZERO)  THEN
C        NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
         DMAX = THREE*DEL1
         IF (DABS(D(1,1)) .GT. DABS(DMAX))  D(1,1) = DMAX
      ENDIF
C
C  LOOP THROUGH INTERIOR POINTS.
C
      DO 50  I = 2, NLESS1
         IF (I .EQ. 2)  GO TO 40
C
         HSUM = H(I-1) + H(I)
         DEL1 = DEL2
         DEL2 = SLOPE(I)
   40    CONTINUE
C
C        SET D(I)=0 UNLESS DATA ARE STRICTLY MONOTONIC.
C
         D(1,I) = ZERO
         IF ( DPCHST(DEL1,DEL2) .LE. ZERO)  GO TO 50
C
C        USE BRODLIE MODIFICATION OF BUTLAND FORMULA.
C
         HSUMT3 = HSUM+HSUM+HSUM
         W1 = (HSUM + H(I-1))/HSUMT3
         W2 = (HSUM + H(I)  )/HSUMT3
         DMAX = DMAX1( DABS(DEL1), DABS(DEL2) )
         DMIN = DMIN1( DABS(DEL1), DABS(DEL2) )
         DRAT1 = DEL1/DMAX
         DRAT2 = DEL2/DMAX
         D(1,I) = DMIN/(W1*DRAT1 + W2*DRAT2)
C
   50 CONTINUE
C
C  SET D(N) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
C     SHAPE-PRESERVING.
C
      W1 = -H(N-1)/HSUM
      W2 = (H(N-1) + HSUM)/HSUM
      D(1,N) = W1*DEL1 + W2*DEL2
      IF ( DPCHST(D(1,N),DEL2) .LE. ZERO)  THEN
         D(1,N) = ZERO
      ELSE IF ( DPCHST(DEL1,DEL2) .LT. ZERO)  THEN
C        NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
         DMAX = THREE*DEL2
         IF (DABS(D(1,N)) .GT. DABS(DMAX))  D(1,N) = DMAX
      ENDIF
C
C  NORMAL RETURN.
C
 5000 CONTINUE
      RETURN
C------------- LAST LINE OF DPCHCI FOLLOWS -----------------------------
      END
