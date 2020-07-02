      SUBROUTINE BVDER(X,Y,YP,G,IPAR)
C***BEGIN PROLOGUE  BVDER
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***REFER TO  BVSUP
C
C **********************************************************************
C     NFC = Number of base solution vectors
C
C     NCOMP = Number of components per solution vector
C
C              1 -- Nonzero particular solution
C     INHOMO =
C              2 OR 3 -- Zero particular solution
C
C             0 -- Inhomogeneous vector term G(X) identically zero
C     IGOFX =
C             1 -- Inhomogeneous vector term G(X) not identically zero
C
C     G = Inhomogeneous vector term G(X)
C
C     XSAV = Previous value of X
C
C     C = Normalization factor for the particular solution
C
C              0   ( IF  NEQIVP = 0 )
C     IVP =
C           Number of differential equations integrated due to
C           the original boundary value problem   ( if  NEQIVP .GT. 0 )
C
C     NOFST - For problems with auxiliary initial value equations,
C             NOFST communicates to the routine FMAT how to access
C             the dependent variables corresponding to this initial
C             value problem.  For example, during any call to FMAT,
C             the first dependent variable for the initial value
C             problem is in position  Y(NOFST + 1).
C             See example in SAND77-1328.
C **********************************************************************
C***ROUTINES CALLED  FMAT,GVEC,UIVP,UVEC
C***COMMON BLOCKS    ML8SZ,MLIVP
C***END PROLOGUE  BVDER
C
      DIMENSION Y(*),YP(*),G(*)
C
C **********************************************************************
C
      COMMON /ML8SZ/ NFC,NCOMP,INHOMO,IGOFX,XSAV,C,IVP
C
C **********************************************************************
C     THE COMMON BLOCK BELOW IS USED TO COMMUNICATE WITH THE USER
C     SUPPLIED SUBROUTINE FMAT.  THE USER SHOULD NOT ALTER THIS
C     COMMON BLOCK.
C
      COMMON /MLIVP/ NOFST
C **********************************************************************
C
C***FIRST EXECUTABLE STATEMENT  BVDER
      IF (IVP .GT. 0)  CALL UIVP(X,Y(IVP+1),YP(IVP+1))
      NOFST = IVP
      NA = 1
      DO 1 K = 1,NFC
      CALL FMAT(X,Y(NA),YP(NA))
      NOFST = NOFST - NCOMP
    1 NA = NA + NCOMP
C
      IF (INHOMO .NE. 1)  RETURN
      CALL FMAT(X,Y(NA),YP(NA))
C
      IF (IGOFX .EQ. 0)  RETURN
      IF (X .EQ. XSAV)  GO TO 2
      IF (IVP .EQ. 0)  CALL GVEC(X,G)
      IF (IVP .GT. 0)  CALL UVEC(X,Y(IVP+1),G)
      XSAV = X
C
C     IF THE USER HAS CHOSEN NOT TO NORMALIZE THE PARTICULAR
C     SOLUTION, THEN C IS DEFINED IN BVPOR TO BE 1.0
C
    2 DO 3 J = 1,NCOMP
      L = NA + J - 1
    3 YP(L) = YP(L)  +  G(J) / C
C
      RETURN
      END
