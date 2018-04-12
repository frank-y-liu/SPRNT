      SUBROUTINE DBKNOT(X,N,K,T)
C***BEGIN PROLOGUE  DBKNOT
C***REFER TO  DB2INK,DB3INK
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***END PROLOGUE  DBKNOT
C
C  --------------------------------------------------------------------
C  DBKNOT CHOOSES A KNOT SEQUENCE FOR INTERPOLATION OF ORDER K AT THE
C  DATA POINTS X(I), I=1,..,N.  THE N+K KNOTS ARE PLACED IN THE ARRAY
C  T.  K KNOTS ARE PLACED AT EACH ENDPOINT AND NOT-A-KNOT END
C  CONDITIONS ARE USED.  THE REMAINING KNOTS ARE PLACED AT DATA POINTS
C  IF N IS EVEN AND BETWEEN DATA POINTS IF N IS ODD.  THE RIGHTMOST
C  KNOT IS SHIFTED SLIGHTLY TO THE RIGHT TO INSURE PROPER INTERPOLATION
C  AT X(N) (SEE PAGE 350 OF THE REFERENCE).
C  DOUBLE PRECISION VERSION OF BKNOT.
C  --------------------------------------------------------------------
C
C  ------------
C  DECLARATIONS
C  ------------
C
C  PARAMETERS
C
      INTEGER
     *        N, K
      DOUBLE PRECISION
     *     X(N), T(*)
C
C  LOCAL VARIABLES
C
      INTEGER
     *        I, J, IPJ, NPJ, IP1
      DOUBLE PRECISION
     *     RNOT
C
C
C  ----------------------------
C  PUT K KNOTS AT EACH ENDPOINT
C  ----------------------------
C
C     (SHIFT RIGHT ENPOINTS SLIGHTLY -- SEE PG 350 OF REFERENCE)
      RNOT = X(N) + 0.10D0*( X(N)-X(N-1) )
      DO 110 J=1,K
         T(J) = X(1)
         NPJ = N + J
         T(NPJ) = RNOT
  110 CONTINUE
C
C  --------------------------
C  DISTRIBUTE REMAINING KNOTS
C  --------------------------
C
      IF (MOD(K,2) .EQ. 1)  GO TO 150
C
C     CASE OF EVEN K --  KNOTS AT DATA POINTS
C
      I = (K/2) - K
      JSTRT = K+1
      DO 120 J=JSTRT,N
         IPJ = I + J
         T(J) = X(IPJ)
  120 CONTINUE
      GO TO 200
C
C     CASE OF ODD K --  KNOTS BETWEEN DATA POINTS
C
  150 CONTINUE
      I = (K-1)/2 - K
      IP1 = I + 1
      JSTRT = K + 1
      DO 160 J=JSTRT,N
         IPJ = I + J
         T(J) = 0.50D0*( X(IPJ) + X(IPJ+1) )
  160 CONTINUE
  200 CONTINUE
C
      RETURN
      END
