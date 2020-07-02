      SUBROUTINE SPLPO(IPAGE,LPAGE)
C
C     1. OPEN UNIT NUMBER IPAGEF AS A RANDOM ACCESS FILE.
C
C     2. THE RECORD LENGTH IS CONSTANT=LPG.
C
C     THIS IS FORTRAN 77 STANDARD FOR THE 'OPEN' STATEMENT.
      IPAGEF=IPAGE
      LPG   =LPAGE
      OPEN(UNIT=IPAGEF,IOSTAT=IOS,ERR=100,STATUS='UNKNOWN',
     *ACCESS='DIRECT',FORM='UNFORMATTED',RECL=LPG)
      RETURN
 100  CONTINUE
      IOPT=1
      NERR=100
      CALL XERRWV('SPLP( ). OPEN() HAS ERROR FLAG (I1).',36,NERR,IOPT,
     *1,IOS,IDUM,0,RDUM,RDUM)
      RETURN
      END
