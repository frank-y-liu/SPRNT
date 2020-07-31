*SPACK
      SUBROUTINE SPACK
     +   (N2,N1,V1,V2,IFIX)
C***BEGIN PROLOGUE  SPACK
C***REFER TO  SODR,SODRC
C***ROUTINES CALLED  SCOPY
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
C***PURPOSE  SELECT THE UNFIXED ELEMENTS OF V2 AND RETURN THEM IN V1
C***END PROLOGUE  SPACK
C
C  VARIABLE DECLARATIONS (ALPHABETICALLY)
C
      INTEGER I
C        AN INDEXING VARIABLE.
      INTEGER IFIX(N2)
C        THE INDICATOR VALUES USED TO DESIGNATE WHETHER THE INDIVIDUAL
C        ELEMENTS OF V2 ARE FIXED AT THEIR INPUT VALUES OR NOT.
C        (FOR DETAILS, SEE DISCUSSION OF IFIXB AND IFIXX IN PROLOGUE OF
C        SUBROUTINE SODR OR SODRC.)
      INTEGER N1
C        THE NUMBER OF ITEMS IN V1.
      INTEGER N2
C        THE NUMBER OF ITEMS IN V2.
      REAL V1(N2)
C        THE VECTOR OF THE UNFIXED ITEMS FROM V2.
      REAL V2(N2)
C        THE VECTOR OF THE FIXED AND UNFIXED ITEMS FROM WHICH THE
C        UNFIXED ELEMENTS ARE TO BE EXTRACTED.
C
C
C***FIRST EXECUTABLE STATEMENT  SPACK
C
C
      N1 = 0
      IF (IFIX(1).GE.0) THEN
         DO 10 I=1,N2
            IF (IFIX(I).NE.0) THEN
               N1 = N1+1
               V1(N1) = V2(I)
            END IF
   10    CONTINUE
      ELSE
         N1 = N2
         CALL SCOPY(N2,V2,1,V1,1)
      END IF
C
      RETURN
      END