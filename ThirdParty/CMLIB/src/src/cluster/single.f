      SUBROUTINE SINGLE(X, COUNT, AVE, SD, XMIN, XMAX, SSQ)
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
C   PURPOSE
C   -------
C
C      INCORPORATES A NEW VALUE INTO THE SUMMARY STATISTICS
C
C   INPUT PARAMETERS
C   ----------------
C
C   SEE SUBROUTINE BUILD FOR PARAMETER DESCRIPTIONS.
C
C   REFERENCE
C   ---------
C
C     HARTIGAN, J. A. (1975).  CLUSTERING ALGORITHMS, JOHN WILEY &
C        SONS, INC., NEW YORK.  PAGE 109.
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
      IF(COUNT.EQ.0.) THEN
         AVE=0.
         SD=0.
         XMIN=R1MACH(2)
         XMAX=-R1MACH(2)
         SSQ=0.
      ENDIF
      COUNT=COUNT+1.
      AVE=AVE+(X-AVE)/COUNT
      IF(COUNT.NE.1.) SSQ=SSQ+COUNT*(X-AVE)**2/(COUNT-1.)
      SD=(SSQ/COUNT)**0.5
      IF(XMIN.GT.X) XMIN=X
      IF(XMAX.LT.X) XMAX=X
      RETURN
      END
