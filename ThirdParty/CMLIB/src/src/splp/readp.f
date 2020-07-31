      SUBROUTINE READP(IPAGE,LIST,RLIST,LPAGE,IREC)
C***BEGIN PROLOGUE  READP
C***REFER TO  SPLP
C     READ RECORD NUMBER IRECN, OF LENGTH LPG, FROM UNIT
C     NUMBER IPAGEF INTO THE STORAGE ARRAY, LIST(*).
C     READ RECORD  IRECN+1, OF LENGTH LPG, FROM UNIT NUMBER
C     IPAGEF INTO THE STORAGE ARRAY RLIST(*).
C
C     TO CONVERT THIS PROGRAM UNIT TO DOUBLE PRECISION CHANGE
C     /REAL (12 BLANKS)/ TO /DOUBLE PRECISION/.
C***ROUTINES CALLED  XERRWV
C***END PROLOGUE  READP
      INTEGER LIST(*)
      REAL RDUM
      REAL             RLIST(*)
C***FIRST EXECUTABLE STATEMENT  READP
      IPAGEF=IPAGE
      LPG   =LPAGE
      IRECN=IREC
      READ(IPAGEF,REC=IRECN,ERR=100)(LIST(I),I=1,LPG)
      READ(IPAGEF,REC=IRECN+1,ERR=100)(RLIST(I),I=1,LPG)
      RETURN
 100  CONTINUE
      IOPT=1
      NERR=100
      CALL XERRWV('SPLP( ). ERROR IN READP().',26,NERR,IOPT,2,LPG,IRECN,
     * 0,RDUM,RDUM)
      RETURN
      END