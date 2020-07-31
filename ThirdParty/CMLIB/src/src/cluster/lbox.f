      SUBROUTINE LBOX(I, J, NL, A, IS)
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
C   PURPOSE
C   -------
C
C      DRAWS A LINE IN THE MATRIX OF BOXES
C
C   INPUT PARAMETERS
C   ----------------
C
C   I,J   INTEGER SCALARS (UNCHANGED ON OUTPUT).
C         THE ROW AND COLUMN OF THE STARTING POSITION FOR THE LINE.
C
C   NL    INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         THE LENGTH OF THE LINE.
C
C   A     MATRIX OF 1-CHARACTER VARIABLES DIMENSIONED 21 BY 132
C            (CHANGED ON OUTPUT).
C         THE MATRIX OF BOXES.
C
C   IS    INTEGER SCALAR (UNCHANGED ON OUTPUT).
C         THE FLAG THAT DETERMINES THE DIRECTION OF THE LINE.
C
C         IS = 1   HORIZONTAL LINE.
C         IS = 2   VERTICAL LINE.
C         IS = 3   SLANTED LINE.
C
C   REFERENCES
C   ----------
C
C     HARTIGAN, J. A. (1975).  CLUSTERING ALGORITHMS, JOHN WILEY &
C        SONS, INC., NEW YORK.  PAGE 55.
C
C     HARTIGAN, J. A. (1975) PRINTER GRAPHICS FOR CLUSTERING. JOURNAL OF
C        STATISTICAL COMPUTATION AND SIMULATION. VOLUME 4,PAGES 187-213.
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
      CHARACTER*1 D,V,H,S,A(21,*)
      DATA D,V,H,S/'.','!','_','/'/
C
      IF(I.LT.1.OR.I.GT.132) RETURN
      IF(J.LT.1.OR.J.GT.132) RETURN
      IF(A(I,J).NE.V) A(I,J)=S
      IF (IS.EQ.1)THEN
         DO 10 K=1,NL
            II=I
            JJ=J-K
   10       IF(JJ.GT.0) A(II,JJ)=H
      ELSE IF (IS.EQ.2) THEN
         DO 20 K=1,NL
            JJ=J
            II=I+K
   20       IF(II.GT.0) A(II,JJ)=V
      ELSE IF (IS.EQ.3) THEN
         DO 30 K=1,NL
            II=I-K
            JJ=J+K
   30       IF(K.NE.NL) A(II,JJ)=S
      ENDIF
      RETURN
      END
