*SODPER
      SUBROUTINE SODPER
     +   (INFO,LUNERR,SHORT,
     +   N,NP,M,
     +   LDSCLD,LDRHO,
     +   LWKMN,LIWKMN,
     +   SCLD,SCLB,W,RHO,
     +   MSGB,ISODR,MSGX,
     +   XPLUSD,LDXPD,NROW,NETA,NTOL)
C***BEGIN PROLOGUE  SODPER
C***REFER TO  SODR,SODRC
C***ROUTINES CALLED  SODPE1,SODPE2,SODPHD
C***DATE WRITTEN   860529   (YYMMDD)
C***REVISION DATE  870204   (YYMMDD)
C***CATEGORY NO.  G2E,I1B1
C***KEYWORDS  ORTHOGONAL DISTANCE REGRESSION,
C             NONLINEAR LEAST SQUARES,
C             ERRORS IN VARIABLES
C***AUTHOR  BOGGS, PAUL T.
C             OPTIMIZATION GROUP/SCIENTIFIC COMPUTING DIVISION
C             NATIONAL BUREAU OF STANDARDS, GAITHERSBURG, MD 20899
C           BYRD, RICHARD H.
C             DEPARTMENT OF COMPUTER SCIENCE
C             UNIVERSITY OF COLORADO, BOULDER, CO 80309
C           DONALDSON, JANET R.
C             OPTIMIZATION GROUP/SCIENTIFIC COMPUTING DIVISION
C             NATIONAL BUREAU OF STANDARDS, BOULDER, CO 80303-3328
C           SCHNABEL, ROBERT B.
C             DEPARTMENT OF COMPUTER SCIENCE
C             UNIVERSITY OF COLORADO, BOULDER, CO 80309
C             AND
C             OPTIMIZATION GROUP/SCIENTIFIC COMPUTING DIVISION
C             NATIONAL BUREAU OF STANDARDS, BOULDER, CO 80303-3328
C***PURPOSE  CONTROLLING ROUTINE FOR PRINTING ERROR REPORTS.
C***END PROLOGUE  SODPER
C
C  VARIABLE DECLARATIONS (ALPHABETICALLY)
C
      INTEGER D1
C        THE FIRST DIGIT OF INFO.
      INTEGER D2
C        THE SECOND DIGIT OF INFO.
      INTEGER D3
C        THE THIRD DIGIT OF INFO.
      INTEGER D4
C        THE FOURTH DIGIT OF INFO.
      INTEGER D5
C        THE FIFTH DIGIT OF INFO.
      LOGICAL HEAD
C        THE INDICATOR VARIABLE USED TO DESIGNATE WHETHER THE PACKAGE
C        HEADING IS TO BE PRINTED (HEAD=.TRUE.) OR NOT (HEAD=.FALSE.).
      INTEGER INFO
C        AN INDICATOR VARIABLE, USED PRIMARILY TO DESIGNATE WHY THE
C        COMPUTATIONS WERE STOPPED.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      LOGICAL ISODR
C        THE INDICATOR VARIABLE USED TO DESIGNATE WHETHER THE SOLUTION
C        IS TO BE FOUND BY ODR (ISODR=.TRUE.) OR BY OLS (ISODR=.FALSE.).
      INTEGER LDRHO
C        THE LEADING DIMENSION OF ARRAY RHO.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      INTEGER LDSCLD
C        THE LEADING DIMENSION OF ARRAY SCLD.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      INTEGER LDXPD
C        THE LEADING DIMENSION OF ARRAY XPLUSD.
      INTEGER LIWKMN
C        THE MINIMUM ACCEPTABLE LENGTH OF ARRAY IWORK.
      INTEGER LUNERR
C        THE LOGICAL UNIT NUMBER USED FOR ERROR MESSAGES.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      INTEGER LWKMN
C        THE MINIMUM ACCEPTABLE LENGTH OF ARRAY WORK.
      INTEGER M
C        THE NUMBER OF COLUMNS OF DATA IN THE INDEPENDENT VARIABLE.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      INTEGER MSGB(NP+1)
C        THE ERROR CHECKING RESULTS FOR THE JACOBIAN WRT BETA.
      INTEGER MSGX(M+1)
C        THE ERROR CHECKING RESULTS FOR THE JACOBIAN WRT X.
      INTEGER N
C        THE NUMBER OF OBSERVATIONS.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      INTEGER NETA
C        THE NUMBER OF RELIABLE DIGITS IN THE MODEL.
      REAL SCLB(NP)
C        THE SCALE OF EACH BETA.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      REAL SCLD(LDSCLD,M)
C        THE SCALE OF EACH DELTA.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      INTEGER NP
C        THE NUMBER OF FUNCTION PARAMETERS.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      INTEGER NROW
C        THE NUMBER OF THE ROW OF THE INDEPENDENT VARIABLE ARRAY AT
C        WHICH THE DERIVATIVE IS TO BE CHECKED.
      INTEGER NTOL
C        THE NUMBER OF DIGITS OF AGREEMENT REQUIRED BETWEEN THE
C        FINITE DIFFERENCE AND THE USER-SUPPLIED DERIVATIVES.
      REAL RHO(LDRHO,M)
C        THE DELTA WEIGHTS.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      LOGICAL SHORT
C        THE INDICATOR VARIABLE USED TO DESIGNATE WHETHER THE USER HAS
C        INVOKED ODRPACK BY THE SHORT-CALL (SHORT=.TRUE.) OR THE LONG-
C        CALL (SHORT=.FALSE.).
      INTEGER UNIT
C        THE LOGICAL UNIT NUMBER USED FOR ERROR MESSAGES.
      REAL W(N)
C        THE OBSERVATIONAL ERROR WEIGHTS.
C        (FOR DETAILS, SEE ODRPACK REFERENCE GUIDE.)
      REAL XPLUSD(LDXPD,M)
C        THE ARRAY X + DELTA.
C
C
C***FIRST EXECUTABLE STATEMENT  SODPER
C
C
C  SET LOGICAL UNIT NUMBER FOR ERROR REPORT
C
      IF (LUNERR.EQ.0) THEN
         RETURN
      ELSE IF (LUNERR.LT.0) THEN
         UNIT = 6
      ELSE
         UNIT = LUNERR
      END IF
C
C  PRINT HEADING
C
      HEAD = .TRUE.
      CALL SODPHD(HEAD,UNIT)
C
      IF (INFO.LT.0) THEN
C
C  PRINT APPROPRIATE ERROR MESSAGE FOR USER INVOKED STOP FROM FUN OR JAC
C
         WRITE (UNIT,1000) INFO
      ELSE
C
C  PRINT APPROPRIATE ERROR MESSAGES FOR ODRPACK INVOKED STOP
C
         D1 = MOD(INFO,100000)/10000
         D2 = MOD(INFO,10000)/1000
         D3 = MOD(INFO,1000)/100
         D4 = MOD(INFO,100)/10
         D5 = MOD(INFO,10)
C
         IF (D1.GE.1 .AND. D1.LE.3) THEN
C
C  PRINT APPROPRIATE MESSAGES FOR ERRORS IN
C     PROBLEM SPECIFICATION PARAMETERS
C     DIMENSION SPECIFICATION PARAMETERS
C     NUMBER OF GOOD DIGITS IN X
C     OBSERVATIONAL ERROR WEIGHTS
C     DELTA WEIGHTS
C
            CALL SODPE1
     +         (UNIT,D1,D2,D3,D4,D5,
     +         N,M,
     +         LDSCLD,LDRHO,
     +         LWKMN,LIWKMN,
     +         SCLD,SCLB,NP,W,RHO)
         ELSE
            CALL SODPE2
     +         (UNIT,
     +         NP,M,
     +         MSGB,ISODR,MSGX,
     +         XPLUSD,LDXPD,NROW,NETA,NTOL)
         END IF
      END IF
C
C  PRINT CORRECT FORM OF CALL STATEMENT
C
      IF (INFO.LT.0 .OR.
     +    (D1.GE.1 .AND. D1.LE.3) .OR.
     +    (D1.EQ.4 .AND. (D2.EQ.2 .OR. D3.EQ.2))) THEN
         IF (SHORT) THEN
            WRITE (UNIT,2100)
         ELSE
            WRITE (UNIT,2200)
         END IF
      END IF
C
      RETURN
C
C  FORMAT STATEMENTS
C
 1000 FORMAT
     +   (//' THE USER HAS REQUESTED THAT THE COMPUTATIONS BE STOPPED'/
     +      ' IN USER-SUPPLIED SUBROUTINE FUN OR JAC.  THE USER-SET'/
     +      ' VALUE OF IFLAG IS ', I5)
 2100 FORMAT
     +   (//' THE CORRECT FORM OF THE CALL STATEMENT IS '//
     +      '       CALL SODR'/
     +      '      +     (FUN,JAC,'/
     +      '      +     N,M,NP,'/
     +      '      +     X,LDX,'/
     +      '      +     Y,'/
     +      '      +     BETA,'/
     +      '      +     RHO,LDRHO,'/
     +      '      +     JOB,'/
     +      '      +     WORK,LWORK,IWORK,LIWORK,'/
     +      '      +     INFO)')
 2200 FORMAT
     +   (//' THE CORRECT FORM OF THE CALL STATEMENT IS '//
     +      '       CALL SODRC'/
     +      '      +     (FUN,JAC,'/
     +      '      +     N,M,NP,'/
     +      '      +     X,LDX,IFIXX,LDIFX,SCLD,LDSCLD,'/
     +      '      +     Y,'/
     +      '      +     BETA,IFIXB,SCLB,'/
     +      '      +     RHO,LDRHO,W,'/
     +      '      +     JOB,TAUFAC,'/
     +      '      +     SSTOL,PARTOL,MAXIT,'/
     +      '      +     IPRINT,LUNERR,LUNRPT,'/
     +      '      +     WORK,LWORK,IWORK,LIWORK,'/
     +      '      +     INFO)')
C
      END
