      SUBROUTINE STOR1(U,YH,V,YP,NTEMP,NDISK,NTAPE)
C***BEGIN PROLOGUE  STOR1
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***REFER TO  BVSUP
C
C **********************************************************************
C             0 -- Storage at output points.
C     NTEMP =
C             1 -- Temporary storage
C **********************************************************************
C***ROUTINES CALLED  (NONE)
C***COMMON BLOCKS    ML8SZ
C***END PROLOGUE  STOR1
C
      DIMENSION U(*),YH(*),V(*),YP(*)
C
C **********************************************************************
C
      COMMON /ML8SZ/ NFC,NCOMP,INHOMO,IGOFX,XSAV,C,IVP
C
C **********************************************************************
C
C***FIRST EXECUTABLE STATEMENT  STOR1
      NCTNF = NCOMP * NFC
      DO 10 J = 1,NCTNF
   10 U(J) = YH(J)
      IF (INHOMO .EQ. 1)  GO TO 30
C
C   ZERO PARTICULAR SOLUTION
C
      IF (NTEMP .EQ. 1)  RETURN
      DO 20 J = 1,NCOMP
   20 V(J) = 0.
      GO TO 70
C
C   NONZERO PARTICULAR SOLUTION
C
   30 IF (NTEMP .EQ. 0)  GO TO 50
C
      DO 40 J = 1,NCOMP
   40 V(J) = YP(J)
      RETURN
C
   50 DO 60 J = 1,NCOMP
   60 V(J) = C * YP(J)
C
C  IS OUTPUT INFORMATION TO BE WRITTEN TO DISK
C
   70 IF (NDISK .EQ. 1)  WRITE (NTAPE) (V(J),J=1,NCOMP),(U(J),J=1,NCTNF)
C
      RETURN
      END
