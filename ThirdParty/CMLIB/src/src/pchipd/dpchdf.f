      DOUBLE PRECISION FUNCTION DPCHDF(K,X,S,IERR)
C***BEGIN PROLOGUE  DPCHDF
C***REFER TO  DPCHCE,DPCHSP
C***ROUTINES CALLED  XERROR
C***REVISION DATE  870707   (YYMMDD)
C***DESCRIPTION
C
C          DPCHDF:   DPCHIP Finite Difference Formula
C
C     Uses a divided difference formulation to compute a K-point approx-
C     imation to the derivative at X(K) based on the data in X and S.
C
C     Called by  DPCHCE  and  DPCHSP  to compute 3- and 4-point boundary
C     derivative approximations.
C
C ----------------------------------------------------------------------
C
C     On input:
C        K      is the order of the desired derivative approximation.
C               K must be at least 3 (error return if not).
C        X      contains the K values of the independent variable.
C               X need not be ordered, but the values **MUST** be
C               distinct.  (Not checked here.)
C        S      contains the associated slope values:
C                  S(I) = (F(I+1)-F(I))/(X(I+1)-X(I)), I=1(1)K-1.
C               (Note that S need only be of length K-1.)
C
C     On return:
C        S      will be destroyed.
C        IERR   will be set to -1 if K.LT.2 .
C        DPCHDF  will be set to the desired derivative approximation if
C               IERR=0 or to zero if IERR=-1.
C
C ----------------------------------------------------------------------
C
C  Reference:  Carl de Boor, A Practical Guide to Splines, Springer-
C              Verlag (New York, 1978), pp. 10-16.
C
C***END PROLOGUE  DPCHDF
C
C ----------------------------------------------------------------------
C
C  Programmed by:  Fred N. Fritsch,  FTS 532-4275, (415) 422-4275,
C                  Mathematics and Statistics Division,
C                  Lawrence Livermore National Laboratory.
C
C  Change record:
C     82-08-05   Converted to SLATEC library version.
C     87-07-07   Corrected XERROR calls for d.p. name(s).
C
C ----------------------------------------------------------------------
C
C  Programming notes:
C
C     To produce a single precision version, simply:
C        a. Change DPCHDF to PCHDF wherever it occurs,
C        b. Change the double precision declarations to real, and
C        c. Change the constant ZERO to single precision.
C
C  DECLARE ARGUMENTS.
C
      INTEGER K, IERR
      DOUBLE PRECISION  X(K), S(K)
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER I, J
      DOUBLE PRECISION  VALUE, ZERO
      DATA  ZERO /0.D0/
C
C  CHECK FOR LEGAL VALUE OF K.
C
C***FIRST EXECUTABLE STATEMENT  DPCHDF
      IF (K .LT. 3)  GO TO 5001
C
C  COMPUTE COEFFICIENTS OF INTERPOLATING POLYNOMIAL.
C
      DO 10  J = 2, K-1
         DO 9  I = 1, K-J
            S(I) = (S(I+1)-S(I))/(X(I+J)-X(I))
    9    CONTINUE
   10 CONTINUE
C
C  EVALUATE DERIVATIVE AT X(K).
C
      VALUE = S(1)
      DO 20  I = 2, K-1
         VALUE = S(I) + VALUE*(X(K)-X(I))
   20 CONTINUE
C
C  NORMAL RETURN.
C
      IERR = 0
      DPCHDF = VALUE
      RETURN
C
C  ERROR RETURN.
C
 5001 CONTINUE
C     K.LT.3 RETURN.
      IERR = -1
      CALL XERROR ('DPCHDF -- K LESS THAN THREE'
     *           , 0, IERR, 1)
      DPCHDF = ZERO
      RETURN
C------------- LAST LINE OF DPCHDF FOLLOWS -----------------------------
      END
