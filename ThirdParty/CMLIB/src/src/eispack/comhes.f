      SUBROUTINE COMHES(NM,N,LOW,IGH,AR,AI,INT)
C***BEGIN PROLOGUE  COMHES
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D4C1B2
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Reduces complex general matrix to complex upper Hess
C            form using stabilized elementary similarity
C            transformations.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure COMHES,
C     NUM. MATH. 12, 349-368(1968) by Martin and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
C
C     Given a COMPLEX GENERAL matrix, this subroutine
C     reduces a submatrix situated in rows and columns
C     LOW through IGH to upper Hessenberg form by
C     stabilized elementary similarity transformations.
C
C     On INPUT
C
C        NM must be set to the row dimension of two-dimensional
C          array parameters as declared in the calling program
C          dimension statement.
C
C        N is the order of the matrix.
C
C        LOW and IGH ARE integers determined by the balancing
C          subroutine  CBAL.  If  CBAL  has not been used,
C          set LOW=1, IGH=N.
C
C        AR and AI contain the real and imaginary parts,
C          respectively, of the complex input matrix.
C
C     On OUTPUT
C
C        AR and AI contain the real and imaginary parts,
C          respectively, of the Hessenberg matrix.  The
C          multipliers which were used in the reduction
C          are stored in the remaining triangles under the
C          Hessenberg matrix.
C
C        INT contains information on the rows and columns
C          interchanged in the reduction.
C          Only elements LOW through IGH are used.
C
C     Calls CDIV for complex division.
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  CDIV
C***END PROLOGUE  COMHES
C
      INTEGER I,J,M,N,LA,NM,IGH,KP1,LOW,MM1,MP1
      REAL AR(NM,N),AI(NM,N)
      REAL XR,XI,YR,YI
      INTEGER INT(IGH)
C
C***FIRST EXECUTABLE STATEMENT  COMHES
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200
C
      DO 180 M = KP1, LA
         MM1 = M - 1
         XR = 0.0E0
         XI = 0.0E0
         I = M
C
         DO 100 J = M, IGH
            IF (ABS(AR(J,MM1)) + ABS(AI(J,MM1))
     1         .LE. ABS(XR) + ABS(XI)) GO TO 100
            XR = AR(J,MM1)
            XI = AI(J,MM1)
            I = J
  100    CONTINUE
C
         INT(M) = I
         IF (I .EQ. M) GO TO 130
C     .......... INTERCHANGE ROWS AND COLUMNS OF AR AND AI ..........
         DO 110 J = MM1, N
            YR = AR(I,J)
            AR(I,J) = AR(M,J)
            AR(M,J) = YR
            YI = AI(I,J)
            AI(I,J) = AI(M,J)
            AI(M,J) = YI
  110    CONTINUE
C
         DO 120 J = 1, IGH
            YR = AR(J,I)
            AR(J,I) = AR(J,M)
            AR(J,M) = YR
            YI = AI(J,I)
            AI(J,I) = AI(J,M)
            AI(J,M) = YI
  120    CONTINUE
C     .......... END INTERCHANGE ..........
  130    IF (XR .EQ. 0.0E0 .AND. XI .EQ. 0.0E0) GO TO 180
         MP1 = M + 1
C
         DO 160 I = MP1, IGH
            YR = AR(I,MM1)
            YI = AI(I,MM1)
            IF (YR .EQ. 0.0E0 .AND. YI .EQ. 0.0E0) GO TO 160
            CALL CDIV(YR,YI,XR,XI,YR,YI)
            AR(I,MM1) = YR
            AI(I,MM1) = YI
C
            DO 140 J = M, N
               AR(I,J) = AR(I,J) - YR * AR(M,J) + YI * AI(M,J)
               AI(I,J) = AI(I,J) - YR * AI(M,J) - YI * AR(M,J)
  140       CONTINUE
C
            DO 150 J = 1, IGH
               AR(J,M) = AR(J,M) + YR * AR(J,I) - YI * AI(J,I)
               AI(J,M) = AI(J,M) + YR * AI(J,I) + YI * AR(J,I)
  150       CONTINUE
C
  160    CONTINUE
C
  180 CONTINUE
C
  200 RETURN
      END
