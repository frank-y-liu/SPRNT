      SUBROUTINE INTERV(XT,LXT,X,ILEFT,MFLAG)
C***BEGIN PROLOGUE  INTERV
C***REFER TO  FC
C
C Computes largest ILEFT in (1,LXT) such that XT(ILEFT) .LE. X
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  INTERV
      DIMENSION XT(LXT)
      DATA ILO /1/
C***FIRST EXECUTABLE STATEMENT  INTERV
      IHI = ILO + 1
      IF (IHI .LT. LXT)                GO TO 20
         IF (X .GE. XT(LXT))           GO TO 110
         IF (LXT .LE. 1)               GO TO 90
         ILO = LXT - 1
                                       GO TO 21
   20 IF (X .GE. XT(IHI))              GO TO 40
   21 IF (X .GE. XT(ILO))              GO TO 100
C *** NOW X .LT. XT(IHI) . FIND LOWER BOUND
   30 ISTEP = 1
   31 IHI = ILO
      ILO = IHI - ISTEP
      IF (ILO .LE. 1)                  GO TO 35
      IF (X .GE. XT(ILO))              GO TO 50
      ISTEP = ISTEP*2
                                       GO TO 31
   35 ILO = 1
      IF (X .LT. XT(1))                GO TO 90
                                       GO TO 50
C *** NOW X .GE. XT(ILO) . FIND UPPER BOUND
   40 ISTEP = 1
   41 ILO = IHI
      IHI = ILO + ISTEP
      IF (IHI .GE. LXT)                GO TO 45
      IF (X .LT. XT(IHI))              GO TO 50
      ISTEP = ISTEP*2
                                       GO TO 41
   45 IF (X .GE. XT(LXT))              GO TO 110
      IHI = LXT
C *** NOW XT(ILO) .LE. X .LT. XT(IHI) . NARROW THE INTERVAL
   50 MIDDLE = (ILO + IHI)/2
      IF (MIDDLE .EQ. ILO)             GO TO 100
C     NOTE. IT IS ASSUMED THAT MIDDLE = ILO IN CASE IHI = ILO+1
      IF (X .LT. XT(MIDDLE))           GO TO 53
         ILO = MIDDLE
                                       GO TO 50
   53    IHI = MIDDLE
                                       GO TO 50
C *** SET OUTPUT AND RETURN
   90 MFLAG = -1
      ILEFT = 1
                                       RETURN
  100 MFLAG = 0
      ILEFT = ILO
                                       RETURN
  110 MFLAG = 1
      ILEFT = LXT
                                       RETURN
      END
