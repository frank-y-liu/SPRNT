      SUBROUTINE SDASLV(NEQ,DELTA,WM,IWM)
C
C***BEGIN PROLOGUE  SDASLV
C***REFER TO  SDASSL
C***ROUTINES CALLED SGESL,SGBSL
C***COMMON BLOCKS    SDA001
C***DATE WRITTEN   830315   (YYMMDD)
C***REVISION DATE  830315   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***END PROLOGUE  SDASLV
C-----------------------------------------------------------------------
C     this routine manages the solution of the linear
C     system arising in the newton iteration.
C     matrices and real temporary storage and
C     real information are stored in the array wm.
C     integer matrix information is stored in
C     the array iwm.
C     for a dense matrix, the linpack routine
C     dgesl is called.
C     for a banded matrix,the linpack routine
C     dgbsl is called.
C-----------------------------------------------------------------------
C
      DIMENSION DELTA(*),WM(*),IWM(*)
      COMMON/SDA001/NPD,NTEMP,LML,LMU,
     *  LMXORD,LMTYPE,
     *  LNST,LNRE,LNJE,LETF,LCTF,LIPVT
C
      MTYPE=IWM(LMTYPE)
      GO TO(100,100,300,400,400),MTYPE
C
C     dense matrix
100   CALL SGESL(WM(NPD),NEQ,NEQ,IWM(LIPVT),DELTA,0)
      RETURN
C
C     dummy section for mtype=3
300   CONTINUE
      RETURN
C
C     banded matrix
400   MEBAND=2*IWM(LML)+IWM(LMU)+1
      CALL SGBSL(WM(NPD),MEBAND,NEQ,IWM(LML),
     *  IWM(LMU),IWM(LIPVT),DELTA,0)
      RETURN
C------end of subroutine sdaslv------
      END
