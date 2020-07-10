      SUBROUTINE BOXES(MM, M, N, A, CLAB, RLAB, TITLE, XMIN, XMAX, NC,
     *                 IVAR, IERR, OUNIT)
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
C   PURPOSE
C   -------
C
C      PRODUCES A PLOT IN THE FORM OF A 3-DIMENSIONAL BOX OF THREE
C      OR MORE VARIABLES FOR EACH CASE.
C
C   DESCRIPTION
C   -----------
C
C   1.  THE ORIGIN OF THE COORDINATE SYSTEM IS THE LOWER RIGHT CORNER
C       OF EACH BOX.  THE VARIABLES APPEAR IN THE ORDER THEY ARE PLACED
C       IN THE IVAR ARRAY.  FOR EXAMPLE, THE FIRST VARIABLE IN THE
C       X-DIRECTION WILL BE THE (1,1) ELEMENT OF THE ARRAY IVAR AND
C       WILL BE REPRESENTED ON THE BOX AS STARTING AT THE ORIGIN MOVING
C       TO THE LEFT UNTIL THE FIRST DIVIDER.  THE SECOND VARIABLE WILL
C       BE THE (1,2) ELEMENT OF THE ARRAY IVAR AND WILL BE REPRESENTED
C       ON THE BOX AS STARTING AT THE FIRST DIVIDER AND MOVE LEFT UNTIL
C       THE NEXT DIVIDER, ETC.  THE FIRST VARIABLE IN THE Y-DIRECTION
C       WILL BE THE (2,1) ELEMENT OF IVAR AND WILL BE REPRESENTED ON
C       THE BOX AS STARTING AT THE ORIGIN AND MOVING VERTICALLY UNTIL
C       THE FIRST DIVISION.  THE FIRST VARIABLE IN THE Z-DIRECTION IS
C       THE (3,1) ELEMENT OF IVAR AND WILL BE REPRESENTED ON THE BOX AS
C       STARTING AT THE ORIGIN AND MOVING IN A DIAGONAL DIRECTION UNTIL
C       THE FIRST DIVISION.  THE BOXES ARE DRAWN ON FORTRAN UNIT OUNIT.
C
C       FOR EXAMPLE, IF VARIABLES 1 AND 2 ARE TO BE PLOTTED ON THE
C       X-AXIS, VARIABLES 3 AND 4 ON THE Y-AXIS, AND VARIABLE 5 ON THE
C       Z-AXIS, SET
C
C             IVAR(1,1) = 1      IVAR(1,2) = 2
C             IVAR(2,1) = 3      IVAR(2,2) = 4
C             IVAR(3,1) = 5
C
C       AN OUTPUT BOX COULD LOOK LIKE THIS:
C
C                 -------------------------------
C                /                   /         /!
C               /                   /         / !
C              /                   /         /  !
C             /                   /         /   !
C            /                   /         /    !
C           /                   /         /     !
C          /                   /         /     /!
C         /-------------------/---------/     / !
C         !                   !         !    /  !
C         !                   !         !   /   !
C        E!                   !         !  /    !
C         !                   !         ! /    /
C         !                   !         !/    /
C         !-------------------!---------!    /F
C         !                   !         !   /
C         !                   !         !  /
C        D!                   !         ! /
C         !                   !         !/
C         !-------------------!---------!A
C                   C              B
C
C       WHERE    A  IS THE ORIGIN OF THE COORDINATE SYSTEM
C                B  IS THE FIRST VARIABLE
C                C  IS THE SECOND VARIABLE
C                D  IS THE THIRD VARIABLE
C                E  IS THE FOURTH VARIABLE
C                F  IS THE FIFTH VARIABLE
C
C   INPUT PARAMETERS
C   ----------------
C
C   MM    INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         THE LEADING DIMENSION OF THE MATRIX A.  MUST BE AT LEAST M.
C
C   M     INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         THE NUMBER OF CASES.
C
C   N     INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         THE NUMBER OF VARIABLES.
C
C   A     REAL MATRIX WHOSE FIRST DIMENSION MUST BE MM AND SECOND
C            DIMENSION MUST BE AT LEAST N (UNCHANGED ON OUTPUT).
C         THE DATA MATRIX.
C
C         A(I,J) IS THE VALUE FOR THE J-TH VARIABLE FOR THE I-TH CASE.
C
C   CLAB  VECTOR OF 4-CHARACTER VARIABLES DIMENSIONED AT LEAST N
C            (UNCHANGED ON OUTPUT).
C         THE LABELS OF THE VARIABLES.
C
C   RLAB  VECTOR OF 4-CHARACTER VARIABLES DIMENSIONED AT LEAST M
C            (UNCHANGED ON OUTPUT).
C         THE LABELS OF THE CASES.
C
C   TITLE 10-CHARACTER VARIABLE (UNCHANGED ON OUTPUT).
C         THE TITLE OF THE DATA SET.
C
C   XMIN  INTEGER VECTOR DIMENSIONED AT LEAST N (UNCHANGED ON OUTPUT).
C         XMIN(I) HOLDS THE MINIMUM VALUE TO BE PLOTTED FOR VARIABLE I.
C
C   XMAX  INTEGER VECTOR DIMENSIONED AT LEAST N (UNCHANGED ON OUTPUT).
C         XMAX(I) HOLDS THE MAXIMUM VALUE TO BE PLOTTED FOR VARIABLE I.
C
C         IF XMIN(I) .GE. XMAX(I) ON INPUT, THEIR VALUES WILL BE
C            DETERMINED BY THE ROUTINE.
C
C   NC    INTEGER VECTOR DIMENSIONED AT LEAST N (UNCHANGED ON OUTPUT).
C         NC(I) HOLDS THE MAXIMUM NUMBER OF DIVISIONS FOR VARIABLE I.
C
C         THE TOTAL NUMBER OF SIDE DIVISIONS FOR THE Y AND Z AXES
C            TOGETHER MUST BE NO GREATER THAN 20.  AN ERROR WILL
C            OCCUR OTHERWISE.
C
C   IVAR  INTEGER MATRIX WHOSE FIRST DIMENSION MUST BE 3 AND WHOSE
C            SECOND DIMENSION MUST BE AT LEAST N (UNCHANGED ON OUTPUT).
C         THE FIRST ROW OF IVAR HOLDS THE VARIABLES TO BE PLOTTED ON THE
C            X-AXIS.  THE SECOND AND THIRD ROWS HOLDS THE VARIABLES FOR
C            THE Y-AXIS AND THE Z-AXIS, RESPECTIVELY.
C
C   OUNIT INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         UNIT NUMBER FOR OUTPUT.
C
C   OUTPUT PARAMETER
C   ----------------
C
C   IERR  INTEGER SCALAR.
C         ERROR FLAG.
C
C         IERR = 0, NO ERROR WAS DETECTED DURING EXECUTION
C
C         IERR = 1, THE NUMBER OF SIDE DIVISIONS FOR THE Y- AND
C                   Z-COORDINATES IS GREATER THAN 20.  REDUCE THE
C                   APPROPRIATE VALUES OF NC.
C
C   REFERENCES
C   ----------
C
C     HARTIGAN, J. A. (1975).  CLUSTERING ALGORITHMS, JOHN WILEY &
C        SONS, INC., NEW YORK.  PAGES 39, 56-57.
C
C     HARTIGAN, J. A. (1975) PRINTER GRAPHICS FOR CLUSTERING. JOURNAL OF
C        STATISTICAL COMPUTATION AND SIMULATION. VOLUME 4,PAGES 187-213.
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
      DIMENSION A(MM,*),XMIN(*),XMAX(*),NC(*),IVAR(3,*)
      CHARACTER*4 CLAB(*), RLAB(*)
      CHARACTER*10 TITLE
      CHARACTER*1 B, P(21,132)
      INTEGER NCS(3), NW(3), OUNIT
      DATA B/' '/
C
C     SET UP DIMENSIONS OF BOX
C
      IERR = 0
      DO 10 J=1,N
   10    IF (XMIN(J) .GE. XMAX(J)) CALL RANGE(M,A(1,J),XMIN(J),XMAX(J))
      DO 20 K=1,3
         NW(K)=0
         DO 20 I=1,N
            JL=IVAR(K,I)
   20    IF (JL.NE.0) NW(K)=NW(K)+NC(JL)
C
C     WRITE RANGE AND BOX ASSIGNMENTS
C
      IF (OUNIT .LE. 0) RETURN
      WRITE(OUNIT,1)
    1 FORMAT('1 VARIABLE ASSIGNMENTS ')
      WRITE(OUNIT,2)
    2 FORMAT('     Z'/'  X-/'/'    !'/'    Y')
      DO 40 K=1,3
         IF(K.EQ.1) WRITE(OUNIT,3)
    3    FORMAT(//' X VARIABLE RANGES'/' ',17('-'))
         IF(K.EQ.2) WRITE(OUNIT,4)
    4    FORMAT(//' Y VARIABLE RANGES'/' ',17('-'))
         IF(K.EQ.3) WRITE(OUNIT,5)
    5    FORMAT(//' Z VARIABLE RANGES'/' ',17('-'))
         DO 30 J=1,N
            IF(IVAR(K,J).EQ.0) GO TO 40
            JL=IVAR(K,J)
            WRITE(OUNIT,6) CLAB(JL), XMIN(JL),XMAX(JL),NC(JL)
   30    CONTINUE
   40 CONTINUE
    6 FORMAT(' VARIABLE : ',A4,'  MINIMUM : ',F8.2,'    MAXIMUM : ',
     * F8.3,'   NUMBER OF INTERVALS : ',I4)
      WRITE(OUNIT,7)
    7 FORMAT('1')
      NDD=NW(3)+NW(2)+1
      IF(NDD.GT.21) THEN
         IERR = 1
         WRITE(OUNIT,8) NDD
    8    FORMAT(1X, 'THE TOTAL NUMBER OF DIVISIONS FOR THE 2ND ',/,
     *          ' AND 3RD SIDES IS ', I4, ' WHICH IS GREATER THAN 20')
         RETURN
      ENDIF
      NR=132/(NW(1)+NW(3)+1)
      IF(NR.EQ.0) RETURN
C
C     LOOP THROUGH CASES
C
      NLIN=0
      NSKIP=(NW(1)+NW(3)+1)/4-1
      IF(NSKIP.LT.1) NSKIP=1
      DO 150 II=1,M,NR
         IL=II
         IU=II+NR-1
         IF(IU.GT.M) IU=M
         DO 50 I=1,21
            DO 50 J=1,132
   50          P(I,J)=B
         NDD=0
         DO 110 I=IL,IU
            JS=(I-IL+1)*(NW(1)+NW(3)+1)-NW(3)
            IS=NW(3)+1
            DO 60 K=1,3
               NCS(K)=0
               DO 60 J=1,N
                  IF(IVAR(K,J).NE.0) THEN
                     JN=IVAR(K,J)
                     NCH=KC(A(I,JN),XMIN(JN),XMAX(JN),NC(JN))
                     NCS(K)=NCS(K)+NCH
                  ENDIF
   60          CONTINUE
            NNI = NCS(2)+NCS(3)+2
            IF(NNI.GE.NDD) NDD=NNI
            DO 70 K=1,3
   70          CALL LBOX(IS,JS,NCS(K),P,K)
C
C     INSERT BOX FOR CASE I
C
            DO 100 KB=1,3
               K=4-KB
               ISS=IS
               JSS=JS
               DO 90 J=1,N
                  IF(IVAR(K,J).EQ.0) GO TO 100
                  JN=IVAR(K,J)
                  NCH=KC(A(I,JN),XMIN(JN),XMAX(JN),NC(JN))
                  IF(K.EQ.1) JSS=JSS-NCH
                  IF(K.EQ.2) ISS=ISS+NCH
                  IF(K.EQ.3) THEN
                     JSS=JSS+NCH
                     ISS=ISS-NCH
                  ENDIF
                  DO 80 KK=1,3
   80                IF(KK.NE.K) CALL LBOX(ISS,JSS,NCS(KK),P,KK)
   90          CONTINUE
  100       IF(K.EQ.3) P(ISS,JSS)=B
  110    CONTINUE
         NLIN=NLIN+NDD
C
C     FINALLY, PRINT OUT ARRAY
C
         DO 140 I=1,21
            DO 120 J=1,132
  120          IF(P(I,J).NE.B) GOTO 130
            GOTO 140
  130       WRITE(OUNIT,9)(P(I,J),J=1,132)
  140    CONTINUE
    9    FORMAT(1X,132A1)
         WRITE(OUNIT,11) ((B,J=1,NSKIP),RLAB(I),I=IL,IU)
   11    FORMAT(33A4)
         IF(NLIN+NW(2)+NW(3).GT.100) THEN
            NLIN=0
            WRITE(OUNIT,7)
         ENDIF
  150 CONTINUE
      RETURN
      END
