      REAL FUNCTION CV(XVAL,NDATA,NCONST,NORD,NBKPT,BKPT,W)
C***BEGIN PROLOGUE  CV
C***REFER TO  FC
C
C     This is a companion function subprogram for FC( ).
C     The documentation for FC( ) has more complete
C     usage instructions.
C
C      Evaluate the variance function of the curve obtained
C      by the constrained B-spline fitting subprogram, FC( ).
C
C     Revised Feb. 14, 1979.
C***ROUTINES CALLED  BSPLVN,SDOT
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***END PROLOGUE  CV
      DIMENSION BKPT(NBKPT), W(*), V(40)
C***FIRST EXECUTABLE STATEMENT  CV
      ZERO = 0.
      MDG = NBKPT - NORD + 3
      MDW = NBKPT - NORD + 1 + NCONST
      IS = MDG*(NORD+1) + 2*MAX0(NDATA,NBKPT) + NBKPT + NORD**2
      LAST = NBKPT - NORD + 1
      ILEFT = NORD
   10 IF (.NOT.(XVAL.GE.BKPT(ILEFT+1) .AND. ILEFT.LT.LAST-1)) GO TO 20
      ILEFT = ILEFT + 1
      GO TO 10
   20 CALL BSPLVN(BKPT, NORD, 1, XVAL, ILEFT, V(NORD+1))
      ILEFT = ILEFT - NORD + 1
      IP = MDW*(ILEFT-1) + ILEFT + IS
      N = NBKPT - NORD
      DO 30 I=1,NORD
         V(I) = SDOT(NORD,W(IP),1,V(NORD+1),1)
         IP = IP + MDW
   30 CONTINUE
      CV = AMAX1(SDOT(NORD,V,1,V(NORD+1),1),ZERO)
C
C     SCALE THE VARIANCE SO IT IS AN UNBIASED ESTIMATE.
      CV = CV/FLOAT(MAX0(NDATA-N,1))
      RETURN
      END
