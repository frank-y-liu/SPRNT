      SUBROUTINE JACOBN(N,T,Y,DFDY,MATDIM,ML,MU)
C***BEGIN PROLOGUE JACOBN
C***REFER TO SDRIV3
C
C   THIS IS A USER-PROVIDED ROUTINE, WHICH IS PROVIDED
C   HERE TO AVOID LOADER ERRORS ONLY.
C
C***ROUTINES CALLED NONE
C***END PROLOGUE JACOBN
      REAL Y(*),DFDY(MATDIM,*)
      CALL XERROR('SDRIV  -- JACOBN IS A USER-SUPPLED ROUTINE',43,1,2)
      RETURN
      END