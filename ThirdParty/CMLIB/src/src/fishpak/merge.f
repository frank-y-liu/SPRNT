      SUBROUTINE MERGE (TCOS, I1, M1, I2, M2, I3)
C***BEGIN PROLOGUE  MERGE
C***SUBSIDIARY
C***PURPOSE  Subsidiary to GENBUN
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (MERGE-S, CMPMRG-C)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     This subroutine merges two ascending strings of numbers in the
C     array TCOS.  The first string is of length M1 and starts at
C     TCOS(I1+1).  The second string is of length M2 and starts at
C     TCOS(I2+1).  The merged string goes into TCOS(I3+1).
C
C***SEE ALSO  GENBUN
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  MERGE
      DIMENSION       TCOS(*)
C
C
C***FIRST EXECUTABLE STATEMENT  MERGE
      J1 = 1
      J2 = 1
      J = I3
      IF (M1 .EQ. 0) GO TO 107
      IF (M2 .EQ. 0) GO TO 104
  101 J = J+1
      L = J1+I1
      X = TCOS(L)
      L = J2+I2
      Y = TCOS(L)
      IF (X-Y) 102,102,103
  102 TCOS(J) = X
      J1 = J1+1
      IF (J1 .GT. M1) GO TO 106
      GO TO 101
  103 TCOS(J) = Y
      J2 = J2+1
      IF (J2 .LE. M2) GO TO 101
      IF (J1 .GT. M1) GO TO 109
  104 K = J-J1+1
      DO 105 J=J1,M1
         M = K+J
         L = J+I1
         TCOS(M) = TCOS(L)
  105 CONTINUE
      GO TO 109
  106 CONTINUE
      IF (J2 .GT. M2) GO TO 109
  107 K = J-J2+1
      DO 108 J=J2,M2
         M = K+J
         L = J+I2
         TCOS(M) = TCOS(L)
  108 CONTINUE
  109 CONTINUE
      RETURN
      END
